import numpy as np
from .params import InitParams, InputParams
from .recode_header import ReCoDeHeader
from .em_reader import MRCReader, SEQReader, emfile
from .recoder import reduce_L4
from .misc import rc_cfg as rc, get_dtype_string
from .fileutils import read_file
import os
import sys
from pathlib import Path
from datetime import datetime
import time
import argparse
import math
from numba import jit
import c_recode
import warnings
import scipy.ndimage as nd
from skimage.measure import regionprops

import zlib
import blosc
import zstandard as zstd
import snappy
import lz4.frame
import lzma
import bz2

class ReCoDeWriter():

    def __init__(self, 
        image_filename, dark_data=None, dark_filename='', 
        output_directory='', input_params=None, params_filename='', mode='batch', 
        validation_frame_gap=-1, log_filename='recode.log', run_name='run', 
        verbosity=0, use_C=False, max_count=-1, chunk_time_in_sec=0, node_id=0):
        """
        Validates and holds user specified parameters for initializing ReCoDe.

        Parameters
        ----------
        image_filename : string
            file to be processed
        dark_data : numpy array
            calibration_data (required if dark_filename is None), 
            if both dark_data and dark_filename are given, dark_data will be used and dark_filename will be ignored
        dark_filename : string
            file containing calibration data (required if dark_data is None)
        output_directory : string
            location where processed data will be written to
        input_params : InputParams or dict
            object of type pyrecode.params.InputParams or a dictionary containing input parameters (required if params_file is None),
            if both input_params and params_filename are given, input_params will be used and params_filename will be ignored
        params_filename : string
            file containing input parameters (required if input_params is None)
        mode: string
            can be either 'batch' or 'stream', indicating offline and online processing modes respectively
        validation_frame_gap : int
            number of frames to skip before saving validation frames
        log_filename : string
            the name of the log file
        run_name : string
            the name used to identify this run in the log file
        verbosity : int
            0, 1 or 2
        use_C : boolean
            indicates whether the optimized C implementation should be used
        max_count : int
            maximum number of data chunks (files) to process when mode='stream', ignored when mode='batch'
        chunk_time_in_sec : int
            acquisition time of data chunks (files) when mode='stream', ignored when mode='batch'
        node_id : int
            index of the processing node / thread in a distributed / multiprocess call
        """

        # parse and validate initialization params
        self._init_params = InitParams( mode, 
                    image_filename, output_directory, dark_filename, 
                    params_filename, validation_frame_gap, log_filename,
                    run_name, verbosity, use_C)

        # parse and validate input params
        if input_params is None:
            self._input_params = InputParams()
            self._input_params.load(Path(self._init_params._params_filename))
        else:
            self._input_params = input_params

        if not self._input_params.validate():
            raise ValueError('Invalid input params')

        # check if initialization and input params are consistent
        if self._init_params.use_C:
            if self._input_params.source_dtype != np.uint16 or self._input_params.target_dtype != np.uint16:
                raise ValueError ('use_C=True can only be used if source and target dtypes are both unsigned 16-bit')

        # create ReCoDe header
        self._header = ReCoDeHeader()
        self._header.create(self._init_params, self._input_params)
        self._header.print()
        if not self._header.validate():
            raise ValueError('Invalid ReCoDe header created')

        # load and validate calibration frame
        if dark_data is None:
            if self._input_params.dark_file_type == rc.FILE_TYPE_MRC:
                t = MRCReader(str(self._init_params.dark_filename))
            elif self._input_params.dark_file_type == rc.FILE_TYPE_SEQ:
                t = SEQReader(str(self._init_params.dark_filename))
            elif self._input_params.dark_file_type == rc.FILE_TYPE_BINARY:
                t = read_file(str(self._init_params.dark_filename), self._header._rc_header['ny'], self._header._rc_header['nx'], self._input_params.source_dtype)
            else:
                raise NotImplementedError("No implementation available for loading calibration file of type 'Other'")
        else:
            t = dark_data
        
        if self._header._rc_header['ny'] != t.shape[0] or self._header._rc_header['nx'] != t.shape[1]:
            raise RuntimeError('Data and Calibration frames have different shapes')

        self._calibration_frame = t
        self._calibration_frame_p_threshold = self._calibration_frame + self._input_params.dark_threshold_epsilon
        if self._input_params.dark_file_type in [rc.FILE_TYPE_MRC, rc.FILE_TYPE_SEQ]:
            t.close()

        # ensure calibration frame has the same data type as source
        self._src_dtype = np.dtype(get_dtype_string(self._header._rc_header["source_dtype"]))
        if self._calibration_frame.dtype != self._src_dtype:
            warnings.warn('Calibration data type not same as source. Attempting to cast.')
            self._calibration_frame = self._calibration_frame.astype(self._src_dtype)
            self._calibration_frame_p_threshold = self._calibration_frame_p_threshold.astype(self._src_dtype)

        self._node_id = node_id

        self._intermediate_file_name = None
        self._intermediate_file = None
        self._validation_file_name = None
        self._validation_file = None

        self._buffer_sz = None
        self._rct_buffer = None
        self._rct_buffer_fill_position = -1
        self._available_buffer_space = self._buffer_sz
        self._frame_sz = None
        self._frame_buffer = None
        self._chunk_offset = None
        self._num_frames_in_part = None

        # variables used for counting on validation frames
        self._vc_struct = nd.generate_binary_structure(2,2)
        self._vc_dose_rate = 0.0
        self._vc_roi = {'x_start':None, 'y_start': None, 'nx': None, 'ny': None}
        self._vc_n_pixels = None


    def start(self):
        '''
        Prepare for processing based on available input parameters. Assume source data is not available at this moment.
        Create output files and internal buffers
        '''

        # create part-file
        base_filename = self._init_params.image_filename.stem
        self._intermediate_file_name = os.path.join(self._init_params._output_directory, base_filename + '.rc' + str(self._input_params.reduction_level) + '_part' + '{0:03d}'.format(self._node_id))
        self._intermediate_file = open(self._intermediate_file_name, 'wb')

        # serialize ReCoDe header
        self._header.serialize_to (self._intermediate_file)

        # serialize source header
        # self._source.serialize_header(str(self._intermediate_file_name))

        # create validation file
        if self._init_params.validation_frame_gap > 0:
            self._validation_file_name = os.path.join(self._init_params._output_directory, base_filename + '_part' + '{0:03d}'.format(self._node_id) + '_validation_frames.bin')
            self._validation_file = open(self._validation_file_name, 'wb')

        # create buffer to hold reduced_compressed data
        self._buffer_sz = 1000            # best to ensure buffer size is large enough to hold the expected amount of data to be processed by this thread for a single chunk
        self._rct_buffer = bytearray(self._buffer_sz)
        self._rct_buffer_fill_position = -1
        self._available_buffer_space = self._buffer_sz

        # self._bytes_per_pixel = np.dtype(get_dtype_string(self._header._rc_header["source_dtype"])).itemsize
        self._bytes_per_pixel = self._src_dtype.itemsize
        self._n_pixels_in_frame = self._header._rc_header['ny'] * self._header._rc_header['nx']
        self._frame_sz = np.uint64(self._n_pixels_in_frame)*self._bytes_per_pixel
        self._frame_buffer = bytearray(self._buffer_sz)
        self._n_bytes_in_binary_image = math.ceil(self._n_pixels_in_frame/8)

        if self._init_params.use_C:
            self._c_reader = c_recode.Reader()
            _max_sz = int(math.ceil((self._n_pixels_in_frame*self._input_params.source_bit_depth*1.0)/8.0))
            self._pixvals = memoryview(bytearray(self._n_pixels_in_frame*self._bytes_per_pixel))
            self._packed_pixvals = memoryview(bytearray(_max_sz))

        self._chunk_offset = 0
        self._num_frames_in_part = 0

        # initialize validation counting parameters
        self._vc_roi['nx'] = min(self._header._rc_header['nx'], 128)
        self._vc_roi['ny'] = min(self._header._rc_header['ny'], 128)
        self._vc_roi['x_start'] = math.floor((self._header._rc_header['nx'] - self._vc_roi['nx']) / 2.0)
        self._vc_roi['y_start'] = math.floor((self._header._rc_header['ny'] - self._vc_roi['ny']) / 2.0)
        self._vc_n_pixels = self._vc_roi['nx']*self._vc_roi['ny']

        if self._input_params.compression_scheme == 1:      #zstd
            self._compressor_context = zstd.ZstdCompressor(level=self._input_params.compression_level, write_content_size=False)

    def _do_sanity_checks(self, data=None):
        '''
        Source data is now available. Validate if input params agree with data. 
        This function is called separately for each chunk by each thread.
        Load source header here.
        '''
        # get source file headers
        if data is None:
            if self._input_params.source_file_type == rc.FILE_TYPE_MRC:
                self._source = MRCReader(str(self._init_params.image_filename))
                self._source_shape = self._source.shape

            elif self._input_params.source_file_type == rc.FILE_TYPE_SEQ:
                self._source = SEQReader(str(self._init_params.image_filename))
                self._source_shape = self._source.shape

            elif self._input_params.source_file_type == rc.FILE_TYPE_BINARY:
                self._source_shape = tuple([self._header._rc_header['nz'], self._header._rc_header['ny'], self._header._rc_header['nx']])

            else:
                raise NotImplementedError("No implementation available for loading calibration file of type 'Other'")
        else:
            self._source = data
            self._source_shape = self._source.shape
        

        # Validate that user given header values agree with source (MRC or SEQ) Header values of source file
        if self._source_shape[1] != self._header._rc_header['ny']:
            raise RuntimeError('Expected height does not match height in source file')

        if self._source_shape[2] != self._header._rc_header['nx']:
            raise RuntimeError('Expected width does not match width in source file')
        
        if self._input_params.num_frames == -1:
            self._header._rc_header['nz'] = self._source_shape[0]
        else:
            if self._input_params.num_frames > self._source_shape[0]:
                raise RuntimeError('Number of frames requested in config file is larger than available in source file')
            else:
                self._header._rc_header['nz'] = self._input_params.num_frames
        
        # close source file to reduce read overhead
        if self._input_params.dark_file_type in [rc.FILE_TYPE_MRC, rc.FILE_TYPE_SEQ]:
            self._source.close()

    def run(self, data=None):
        '''
        Source data is now available. This function is called separately for each chunk by each thread.
        Each thread will independently read and process data
        '''
        run_metrics = {}

        # do sanity checks
        self._do_sanity_checks(data)

        # determine the number of available frames for this thread
        if self._init_params.mode == 'batch':
            n_frames_in_chunk = self._input_params.nz
        elif self._init_params.mode == 'stream':
            n_frames_in_chunk = self._source_shape[0]
        n_frames_per_thread = math.ceil((n_frames_in_chunk*1.0) / (self._input_params.num_threads*1.0))
        frame_offset = self._node_id * n_frames_per_thread
        available_frames = min(n_frames_per_thread, n_frames_in_chunk - frame_offset)

        print(n_frames_in_chunk, n_frames_per_thread, frame_offset, available_frames, self._chunk_offset)

        # read the thread-specific data from chunk into memory
        stt=datetime.now()
        if data is None:
            with emfile(str(self._init_params.image_filename), self._input_params.source_file_type, mode="r") as f:
                data = f[frame_offset:frame_offset+available_frames]
        if data.dtype != self._src_dtype:
            warnings.warn('Source data type not as specified. Attempting to cast.')
            data = data.astype(self._src_dtype)
        run_metrics['run_data_read_time'] = datetime.now()-stt

        # process frames
        run_start=datetime.now()
        for count, frame in enumerate(data):
            
            absolute_frame_index = self._chunk_offset + frame_offset + count
            compressed_frame_length, _metrics, binary_frame = self._reduce_compress(frame, absolute_frame_index)
            
            # if buffer doesn't have enough space to hold new data offload buffer data to file
            if self._available_buffer_space < compressed_frame_length:
                self._offload_buffer()
                self._rct_buffer_fill_position = 0
                self._available_buffer_space = self._buffer_sz
                
            # copy the last received frame that was never written to buffer
            self._rct_buffer[self._rct_buffer_fill_position:self._rct_buffer_fill_position+compressed_frame_length] = self._frame_buffer[:compressed_frame_length]
            self._rct_buffer_fill_position += compressed_frame_length
            self._available_buffer_space -= compressed_frame_length

            # serialize validation data
            if self._init_params.validation_frame_gap > 0:
                if absolute_frame_index % self._init_params.validation_frame_gap == 0:
                    # serialize validation frame
                    self._validation_file.write(frame.tobytes())
                    # count for dose rate estimation
                    _vc_frame = binary_frame[self._vc_roi['y_start']:self._vc_roi['y_start']+self._vc_roi['ny'], self._vc_roi['x_start']:self._vc_roi['x_start']+self._vc_roi['nx']]
                    labeled_foreground, num_features = nd.measurements.label(_vc_frame, structure=self._vc_struct)
                    self._vc_dose_rate = num_features/self._vc_n_pixels
                    if 'run_dose_rates' in run_metrics:
                        run_metrics['run_dose_rates'].append(self._vc_dose_rate)
                    else:
                        run_metrics['run_dose_rates'] = [self._vc_dose_rate]

            for key in _metrics:
                if key in run_metrics:
                    run_metrics[key] += _metrics[key]
                else:
                    run_metrics[key] = _metrics[key]

        self._chunk_offset += available_frames
        self._num_frames_in_part += available_frames

        run_metrics['run_time'] = datetime.now()-run_start
        run_metrics['run_frames'] = n_frames_in_chunk
        return run_metrics
        
    def _reduce_compress(self, frame, absolute_frame_index):

        run_metrics = {}

        start=datetime.now()

        stt=datetime.now()
        binary_frame = frame > self._calibration_frame_p_threshold  # ensure dtypes match
        pixel_intensities = frame[binary_frame] - self._calibration_frame_p_threshold[binary_frame]
        # cast pixel intensities to appropriate type
        print(type(pixel_intensities[0]))
        run_metrics['frame_thresholding_time'] = datetime.now()-stt

        # packed_binary_frame = self._pack_binary_frame(binary_frame)
        stt=datetime.now()
        packed_binary_frame = bytearray(_pack_binary_frame(binary_frame, self._n_bytes_in_binary_image))
        run_metrics['frame_binary_image_packing_time'] = datetime.now()-stt

        stt=datetime.now()
        if self._input_params.source_bit_depth % 8 == 0:
            packed_pixel_intensities = pixel_intensities.tobytes()
        else:
            if self._init_params.use_C:
                n_pixels = len(pixel_intensities)
                n_packed = int(math.ceil((n_pixels*self._input_params.source_bit_depth*1.0)/8.0))
                b = pixel_intensities.tobytes()
                self._pixvals[:len(b)] = b
                self._c_reader._bit_pack_pixel_intensities(n_packed, n_pixels, self._input_params.source_bit_depth, self._pixvals, self._packed_pixvals)
                packed_pixel_intensities = np.frombuffer(self._packed_pixvals, dtype=np.uint8, count=n_packed)
            else:
                packed_pixel_intensities = _bit_pack(pixel_intensities, self._input_params.source_bit_depth)
        _n_bytes_in_packed_pixvals = len(packed_pixel_intensities)
        run_metrics['frame_pixel_intensity_packing_time'] = datetime.now()-stt

        if self._input_params.rc_operation_mode == 0:   # Reduce Only

            # get compresed sizes
            compressed_frame_length =  8 + self._n_bytes_in_binary_image + _n_bytes_in_packed_pixvals   # the additional 8 is for the frame_id and _n_bytes_in_packed_pixvals
            if compressed_frame_length > self._frame_sz:
                raise ValueError('Buffer size smaller than compressed data size')

            # copy all elements of d into self._compressed_frame_data
            for index, val in enumerate([absolute_frame_index, _n_bytes_in_packed_pixvals]):
                t = np.zeros([1], dtype=np.uint32)
                t[0] = val
                self._frame_buffer[index*4:(index+1)*4] = t.tobytes()
            self._frame_buffer[8:8+self._n_bytes_in_binary_image] = packed_binary_frame
            self._frame_buffer[8+self._n_bytes_in_binary_image:8+self._n_bytes_in_binary_image+_n_bytes_in_packed_pixvals] = packed_pixel_intensities

        elif self._input_params.rc_operation_mode == 1: # Reduce-Compress

            # compress
            stt=datetime.now()
            #  compressed_binary_frame = zlib.compress(packed_binary_frame, self._input_params.compression_level)
            compressed_binary_frame = self._compress(packed_binary_frame)
            run_metrics['frame_binary_image_compression_time'] = datetime.now()-stt
            stt=datetime.now()
            # compressed_packed_pixel_intensities = zlib.compress(packed_pixel_intensities, self._input_params.compression_level)
            compressed_packed_pixel_intensities = self._compress(packed_pixel_intensities)
            run_metrics['frame_pixel_intensity_compression_time'] = datetime.now()-stt

            # get compresed sizes
            _n_compressed_binary_frame = len(compressed_binary_frame)
            _n_compressed_packed_pixel_intensities = len(compressed_packed_pixel_intensities)
            compressed_frame_length =  16 + _n_compressed_binary_frame + _n_compressed_packed_pixel_intensities   # the additional 16 comes from the frame_id and the three sizes
            if compressed_frame_length > self._frame_sz:
                raise ValueError('Buffer size smaller than compressed data size')

            # copy all elements of d into self._compressed_frame_data
            for index, val in enumerate([absolute_frame_index, _n_compressed_binary_frame, _n_compressed_packed_pixel_intensities]):
                t = np.zeros([1], dtype=np.uint32)
                t[0] = val
                self._frame_buffer[index*4:(index+1)*4] = t.tobytes()
            self._frame_buffer[16:16+_n_compressed_binary_frame] = compressed_binary_frame
            self._frame_buffer[16+_n_compressed_binary_frame:16+_n_compressed_binary_frame+_n_compressed_packed_pixel_intensities] = compressed_packed_pixel_intensities

            print('No. of foreground pixels', len(pixel_intensities), 
                  '_n_bytes_in_packed_pixvals', _n_bytes_in_packed_pixvals, 
                  '_n_compressed_binary_frame', _n_compressed_binary_frame, 
                  '_n_compressed_packed_pixel_intensities', _n_compressed_packed_pixel_intensities)

            '''
            for i in range(10):
                print(i, 'packed_pixel_intensities', packed_pixel_intensities[i], 'pixel_intensities', pixel_intensities[i])
            '''
        else:
            raise ValueError('Unknown RC Operation Mode')

        run_metrics['frame_time'] = datetime.now()-start

        return compressed_frame_length, run_metrics, binary_frame

    def _pack_binary_frame(self, binary_frame):
        x = bytearray(self._n_bytes_in_binary_image)
        count = 0
        index = 0
        for b in binary_frame.flatten():
            if b == 1:
                x[count] += pow(2,index)*b
            index += 1
            if index == 8:
                count += 1
                index = 0
        return x

    def close(self):
        # clear buffer (send remaining data to file)
        self._offload_buffer()

        # serialize the true number of frames and the process id ReCoDe header: self._num_frames_in_part
        self._header.update('nz', self._num_frames_in_part)
        self._intermediate_file.seek(0)
        self._header.serialize_to (self._intermediate_file)

        # close part file
        self._intermediate_file.close()

        # close validation file
        if self._init_params.validation_frame_gap > 0:
            self._validation_file.close()


    def _offload_buffer(self):
        self._intermediate_file.write(self._rct_buffer[:self._rct_buffer_fill_position])
        self._intermediate_file.flush()

    def _compress(self, packed_binary_frame):

        if self._input_params.compression_scheme == 0:        #zlib
            return zlib.compress(packed_binary_frame, self._input_params.compression_level)

        elif self._input_params.compression_scheme == 1:      #zstd
            return self._compressor_context.compress(packed_binary_frame)

        elif self._input_params.compression_scheme == 2:      #lz4
            return lz4.frame.compress(packed_binary_frame, compression_level=self._input_params.compression_level, store_size=False)

        elif self._input_params.compression_scheme == 3:      #snappy
            return snappy.compress(packed_binary_frame)

        elif self._input_params.compression_scheme == 4:      #bzip
            return bz2.compress(packed_binary_frame, compresslevel=self._input_params.compression_level)

        elif self._input_params.compression_scheme == 5:      #lzma
            return lzma.compress(packed_binary_frame, preset=self._input_params.compression_level)

        elif self._input_params.compression_scheme == 6:      #blosc_zlib
            return blosc.compress(packed_binary_frame, clevel=self._input_params.compression_level, cname='zlib', shuffle=blosc.BITSHUFFLE)

        elif self._input_params.compression_scheme == 7:      #blosc_zstd
            return blosc.compress(packed_binary_frame, clevel=self._input_params.compression_level, cname='zstd', shuffle=blosc.BITSHUFFLE)

        elif self._input_params.compression_scheme == 8:      #blosc_lz4
            return blosc.compress(packed_binary_frame, clevel=self._input_params.compression_level, cname='lz4', shuffle=blosc.BITSHUFFLE)

        elif self._input_params.compression_scheme == 9:      #blosc_snappy
            return blosc.compress(packed_binary_frame, clevel=self._input_params.compression_level, cname='snappy', shuffle=blosc.BITSHUFFLE)

        elif self._input_params.compression_scheme == 10:      #blosclz
            return blosc.compress(packed_binary_frame, clevel=self._input_params.compression_level, cname='blosclz', shuffle=blosc.BITSHUFFLE)

        elif self._input_params.compression_scheme == 11:      #blosc_lz4hc
            return blosc.compress(packed_binary_frame, clevel=self._input_params.compression_level, cname='lz4hc', shuffle=blosc.BITSHUFFLE)

        else:
            raise NotImplementedError('compression scheme not implemented')

def print_run_metrics(run_metrics):
    for key in run_metrics:
        if key.startswith('frame_'):
            print(key, "\t", run_metrics[key]/run_metrics['run_frames'], "\t", run_metrics[key]/run_metrics['frame_time'])
        else:
            if key == 'run_dose_rates':
                print(key, "\t", run_metrics[key], "\t", 'Avg.=', np.mean(run_metrics[key]))
            else:
                print(key, "\t", run_metrics[key])


@jit(nopython=True)
def _pack_binary_frame(binary_frame, n_bytes_in_binary_image):
    x = np.zeros((n_bytes_in_binary_image), dtype=np.uint8)
    count = 0
    index = 0
    for b in binary_frame.flatten():
        if b == 1:
            x[count] |= (1 << index)
        index += 1
        if index == 8:
            count += 1
            index = 0
    return x


@jit(nopython=True)
def _bit_pack(pixel_intensities, source_bit_depth):
    n_pixels = len(pixel_intensities)
    n_packed = int(math.ceil((n_pixels*source_bit_depth*1.0)/8.0))
    packed = np.zeros((n_packed), dtype=np.uint8)
    bp = 0          # bit position in jth byte of target byte array
    j = 0           # byte number in target byte array
    for pixval in pixel_intensities:
        for i in range(source_bit_depth):
            if pixval & (1 << i):
                packed[j] = packed[j] | (1 << bp)
            bp += 1
            if bp == 8:
                j += 1
                bp = 0
    return packed


if __name__== "__main__":

    parser = argparse.ArgumentParser(description='ReCoDe Queue Manager')
    parser.add_argument('--image_filename', dest='image_filename', action='store', default='', help='path of folder containing data (typically inside RAM disk for on-the-fly)')
    parser.add_argument('--calibration_file', dest='calibration_file', action='store', default='', help='path to calibration file')
    parser.add_argument('--out_dir', dest='out_dir', action='store', default='', help='output directory')
    parser.add_argument('--params_file', dest='params_file', action='store', default='', help='path to params file')
    parser.add_argument('--mode', dest='mode', action='store', default='batch', help='batch or stream')
    parser.add_argument('--validation_frame_gap', dest='validation_frame_gap', action='store', type=int, default=-1, help='validation frame gap')
    parser.add_argument('--log_file', dest='log_file', action='store', default='', help='path to log file')
    parser.add_argument('--run_name', dest='run_name', action='store', default='run_1', help='run name')
    parser.add_argument('--verbosity', dest='verbosity', action='store', type=int, default=0, help='verbosity level')
    parser.add_argument('--use_c', dest='use_c', action='store_true', help='')
    parser.add_argument('--max_count', dest='max_count', action='store', type=int, default=1, help='the number of chunks to process')
    parser.add_argument('--chunk_time_in_sec', dest='chunk_time_in_sec', action='store', type=int, default=1, help='seconds of data contained in each chunk')
    

    args = parser.parse_args()

    writer = ReCoDeWriter (
        args.image_filename, 
        args.calibration_file, 
        output_directory=args.out_dir, 
        input_params=None, 
        params_filename=args.params_file, 
        mode='batch', 
        validation_frame_gap=-1, 
        log_filename='recode.log', 
        run_name='run', 
        verbosity=0, use_C=False, max_count=-1, chunk_time_in_sec=0, node_id=0)

    writer.start()
    run_metrics = writer.run()
    writer.close()
    print(run_metrics)


    '''
    print(self._input_params.source_bit_depth)
    a = np.random.randint(0, high=4096, size=10)
    pa = _bit_pack(a, 12)
    print(pa)
    '''