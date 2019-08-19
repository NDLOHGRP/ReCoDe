from .params import InitParams, InputParams
from .recode_header import ReCoDeHeader
from .em_reader import MRCReader, SEQReader
from .recoder import reduce_L4

class ReCoDeWriter():
    
    def __init__(self, image_filename, dark_filename, output_directory='', input_params=None, params_filename='', validation_frame_gap=-1, log_filename='recode.log', run_name='run', verbosity=0, use_C=False):
        """
        Validates and holds user specified parameters for initializing ReCoDe.

        Parameters
        ----------
        image_filename : string
            file to be processed
        dark_filename : string
            file containing calibration data (required ir mode is 'rc')
        output_directory : string
            location where processed data will be written to
        params_filename : string
            file containing input parameters (required ir mode is 'rc')
        validation_frame_gap : int
            number of frames to skip before saving validation frames
        log_filename : string
            the name of the log file
        run_name : string
            the name used to identify this run in the log file
        verbosity : int
            0, 1 or 2
        use_C : boolean
            indicates if the optimized C implementation will be used
        """
        self._image_filename = image_filename
        self._dark_filename = dark_filename
        self._output_directory = output_directory
        self._params_filename = params_filename
        self._validation_frame_gap = validation_frame_gap
        self._log_filename = log_filename
        self._run_name = run_name
        self._verbosity = verbosity
        self._use_C = use_C

        assert self._dark_filename != '', 'Dark filename cannot be empty'
        assert self._image_filename != '', 'Image filename cannot be empty'

        self._init_params = InitParams('rc', self._image_filename, self._output_directory, self._dark_filename, self._params_filename, self._validation_frame_gap, self._log_filename, self._run_name, self._verbosity, self._use_C)

        if input_params is None:
            self._input_params = InputParams(self._init_params._params_filename)
        else:
            self._input_params = input_params
        
        # EMWriterBase.__init__(self, self._init_params._source_file_name, 'rc')

    def create(self):

        self._header = ReCoDeHeader(self._init_params, self._input_params)

        # Validate that user given header values agree with MRC and SEQ Header values
        if self._input_params.source_file_type == 'mrc' or self._input_params.source_file_type == 'mrcs':
            
            self._source = MRCReader(self._init_params.image_filename, self._input_params.source_file_type)
            self._source_header = self._source.header
            self._header._nx = self._source_header['nx']
            self._header._ny = self._source_header['ny']
            if self._input_params.num_frames == -1:
                self._header._nz = self._source_header['nz']
            self._header._bit_depth = self._source_header['dtype']
        
        elif self._input_params.source_file_type == 'seq':
            self._source = SEQReader(self._init_params.image_filename, self._input_params.source_file_type)
            self._source_header = self._source.header
            if self._input_params.num_frames == -1:
                self._header._nz = self._source_header['nz']

        self._header.validate()

        if self._init_params.use_C:
            assert self._header._bit_depth == 15, 'The optimized C version of ReCoDe only supports unsigned 16-bit data'

        if self._input_params.dark_file_type == 'mrc' or self._input_params.dark_file_type == 'mrcs':
            t = MRCReader(self._init_params.dark_filename, self._input_params.dark_file_type)
        elif self._input_params.dark_file_type == 'seq':
            t = SEQReader(self._init_params.dark_filename, self._input_params.dark_file_type)

        assert self._source.shape[1] == t.shape[1] and self._source.shape[2] == t.shape[2], 'Data and Calibration frames have different shapes'
        self._calibration_frame = t[0]

        self._file_handles = []

    def write(self):
        counted_frames = reduce_L4(self._source, self._calibration_frame, self._input_params.dark_threshold_epsilon, offset=0, max_frames=10)
        
    def flush(self):
        self._file_handle.flush()
        
    def close(self):
        
        if self._has_new_data:
            self.write()

        self.flush()
        self._file_handle.close()

    def validate(self):
        #To Do: Add checks for shape and dtype
        return True

