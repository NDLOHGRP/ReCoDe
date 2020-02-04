import os
import numpy as np
import c_recode
from .recode_header import ReCoDeHeader
from scipy.sparse import coo_matrix
from pathlib import Path

class ReCoDeReader():
    def __init__(self, file, is_intermediate=False):
        self._source_filename = file
        self._current_frame_index = 0
        self._c_reader = c_recode.Reader()
        if is_intermediate:
            self._is_intermediate = 1
        else:
            self._is_intermediate = 0
        self._header = None
        self._frame_metadata = None
        self._seek_table = None
        self._rc_header = None               # initialized in _load_header
        self._frame_data_start_position = 0  # the byte where frame data for frame 0 starts; initialized in _load_header
        self._frame_sz = 0                   # size of frames in bytes; initialized in _load_header
        self._frame_buffer = None            # initialized in _load_header
        self._n_frame_metadata_elems = 3
        self._sz_frame_metadata_elemes = 4
        self._sz_frame_metadata = self._n_frame_metadata_elems*self._sz_frame_metadata_elemes

    def open(self):
        self._open()
        self._load_header()
        self._create_read_buffers()
        self._load_frame_metadata()

    def _open(self):
        self._c_reader._open_file(self._source_filename)
        
    def _load_header(self):
        self._rc_header = ReCoDeHeader()
        self._rc_header.load(self._source_filename)
        self._rc_header.print()
        self._header = self._rc_header._rc_header
        
        if self._is_intermediate:
            self._frame_data_start_position = self._rc_header._rc_header_length
        else:
            self._frame_data_start_position = self._rc_header._rc_header_length + self._header['nz']*self._sz_frame_metadata
        self._frame_sz = np.uint64(self._header["ny"])*np.uint64(self._header["nx"])*np.uint64(3)*np.uint64(2)
        self._frame_buffer = memoryview(bytes(self._frame_sz))
        return self._rc_header._rc_header

    def _create_read_buffers(self):
        assert 'nz' in self._header, print('Attempting to set persistent variables before reading header')
        self._c_reader._create_read_buffers(self._header["ny"], self._header["nx"], self._header["bit_depth"])
    
    def _load_frame_metadata(self):
        self._frame_metadata = np.zeros((self._header['nz'],3), dtype=np.uint32)
        if not self._is_intermediate:
            assert 'nz' in self._header, print('Attempting to read seek table before reading header')
            # n_compressed_bytes_in_binary_image, n_compressed_bytes_in_pixvals, n_bytes_in_packed_pixvals
            with open(self._source_filename, "rb") as f:
                self._rc_header.skip_header(f)
                for _frame_index in range(self._header['nz']):
                    p = f.read(self._sz_frame_metadata)
                    if not p:
                        return 0
                    self._frame_metadata[_frame_index] = np.frombuffer(p, dtype=np.uint32)

            self._seek_table = np.zeros((self._header['nz'],2), dtype=np.uint32)
            self._seek_table[0,:] = 0
            self._seek_table[1:,0] = np.sum(self._frame_metadata[:-1,0:2], axis=1)
            self._seek_table[:,1] = np.cumsum(self._seek_table[:,0])

            # print(self._frame_metadata)
            # print(self._seek_table)

        return 0

    def get_true_shape(self):
        return tuple([self._header["nz"], self._header["ny"], self._header["nx"]])

    def _get_shape(self):
        return tuple([self._header["nz"], self._header["ny"], self._header["nx"]])

    def _get_dtype(self):
        return (self._stack.dtype)

    def _get_sub_volume(self, slice_z, slice_y, slice_x):
        raise NotImplementedError

    def _get_frame(self, z):
        assert not self._is_intermediate, print("Random acceess is not available for intermediate files")
        assert z < self._header['nz'], print('Cannot get frame' + str(z) + ' of dataset with' + str(self._header['nz']) + 'frames')
        s = self._frame_data_start_position
        self._c_reader._fseek(s + self._seek_table[z,1], int(0))
        N = self._c_reader._get_frame_sparse_L1(self._frame_metadata[z,0], self._frame_metadata[z,1], self._frame_metadata[z,2], self._is_intermediate, self._frame_buffer)
        sparse_d = self._make_coo_frame(self._frame_buffer, N)
        if N != 0:
            self._current_frame_index = z + 1
        return {z: sparse_d}

    '''
    def _get_next_frame_depricated (self):
        print("Using a deprecated version of this function")
        z = self._current_frame_index
        if z == 0:
            self._c_reader._fseek(self._frame_data_start_position, int(0))
        assert z < self._header['nz'], print('Cannot get frame' + str(z) + ' of dataset with' + str(self._header['nz']) + 'frames')
        # self._c_reader._fseek(self._seek_table[z,0], int(1))
        # print('_sz =', self._frame_sz)
        frame_data = memoryview(bytes(self._frame_sz))
        N = self._c_reader._get_frame_sparse_L1(self._frame_metadata[z,0], self._frame_metadata[z,1], self._frame_metadata[z,2], self._is_intermediate, frame_data)
        sparse_d = self._make_coo_frame(frame_data, N)
        # print(sparse_d)
        # print('Num elements =', N)
        if N != 0:
            self._current_frame_index += 1
        return {z: sparse_d}
    '''

    def _get_next_frame (self):
        z = self._current_frame_index
        if z == 0:
            self._c_reader._fseek(self._frame_data_start_position, int(0))
        assert z < self._header['nz'], print('Cannot get frame' + str(z) + ' of dataset with' + str(self._header['nz']) + 'frames')
        
        if self._is_intermediate:
            frame_id = self._c_reader._get_frame_sparse_L1(self._frame_metadata[z,0], self._frame_metadata[z,1], self._frame_metadata[z,2], self._is_intermediate, 1, self._frame_buffer)
        else:
            frame_id = z
        N = self._c_reader._get_frame_sparse_L1(self._frame_metadata[z,0], self._frame_metadata[z,1], self._frame_metadata[z,2], self._is_intermediate, 0, self._frame_buffer)
        if N >= 0:
            a = np.frombuffer(self._frame_buffer, dtype=np.uint16, count=N*3)
            sparse_d = self._make_coo_frame(a, N)
            self._current_frame_index += 1
            return {frame_id: sparse_d}
        else:
            print('Reached EoF after ' + str(z) + ' frames: Quitting')
            self._header['nz'] = self._current_frame_index
            return None
            
        return {z: sparse_d}

    def _make_coo_frame(self, frame_data, N):
        d = np.asarray(frame_data, dtype=np.uint16)
        d = np.transpose(np.reshape(d[:N*3], [N, 3]))
        sparse_d = coo_matrix((d[2], (d[0],d[1])), shape=(self._header['ny'], self._header['nx']), dtype=np.uint16)
        return sparse_d

    def close(self):
        self._c_reader._close_file()

'''
To be done: 
0. saving the frames in a dictionary makes the call to np.save too slow due to slow for loop. instead save the frame ids as a separate list and the COO frames as a separate list. Save the two arrays in a dict with two keys {frame_ids: [], frame_data: []}
1. add option to merge into a rc1 file (not a npy dict). the rc1 file can be read sequentially which is necessary for large datasets. npy dicts need to be loaded all at once into memory.
2. add optional parameter 'validate_part_files' that ensures that all part files are indeed parts of the same file (based on header). Implement and use header.is_equal(other_header) function in recode_header.py.
3. add multithreading support
4. add option to automatically find part files, based on name
'''
def merge_parts(folder_path, base_filename, num_parts, with_duplicate_check=False):
    _folder_path = Path(folder_path)
    frames = {}
    lens = {}
    for index in range(num_parts):
        intermediate_file_name = os.path.join(_folder_path, base_filename + '_part' + '{0:03d}'.format(index))
        reader = ReCoDeReader(intermediate_file_name, is_intermediate=True)
        reader.open()
        print('Processing file ' + str(index) + ' of ' + str(num_parts))
        frames_i = {}
        for i in range(reader._get_shape()[0]):
            d = reader._get_next_frame()
            if d is not None:
                frames_i.update(d)
            else:
                break
        reader.close()
        if with_duplicate_check:
            for frame_id in frames_i:
                if frame_id in frames:
                    print('Found duplicate frame id ' + str(frame_id) + ' in ' + intermediate_file_name)
        frames.update(frames_i)
        lens[index] = len(frames_i)
    print('Extracted ' + str(len(frames)) + ' frames: ')
    print(lens)
    return frames


if __name__== "__main__":

    file_name = 'D:/cbis/GitHub/ReCoDe/scratch/400fps_dose_43.rc1'
    intermediate_file_name = 'D:/cbis/GitHub/ReCoDe/scratch/400fps_dose_43.rc1_part000'

    reader = ReCoDeReader(file_name, is_intermediate=False)
    reader.open()
    frame_data = reader._get_frame(5)
    for i in range(3):
        frame_data = reader._get_next_frame()    
    print(frame_data)
    reader._close()

