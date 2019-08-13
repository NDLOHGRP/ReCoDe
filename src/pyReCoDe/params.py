from abc import ABC, abstractmethod
import numpy as np

class InitParams():

    def __init__(self, mode, image_filename, output_directory, dark_filename='', params_filename='', validation_frame_gap=-1, log_filename='recode.log', run_name='run', verbosity=0, use_C=False):
        """
        Validates and holds user specified parameters for initializing ReCoDe.

        Parameters
        ----------
        mode : string
            the run mode 'rc' for reduction compression and 'de' for decompression expansion
        verbosity : int
            0, 1 or 2
        validation_frame_gap : int
            number of frames to skip before saving validation frames
        image_filename : string
            file to be processed
        dark_filename : string
            file containing calibration data (required ir mode is 'rc')
        params_filename : string
            file containing input parameters (required ir mode is 'rc')
        output_directory : string
            location where processed data will be written to
        log_filename : string
            the name of the log file
        run_name : string
            the name used to identify this run in the log file
        use_C : boolean
            indicates if the optimized C implementation will be used
        """
        self._mode = mode
        self._verbosity = verbosity
        self._validation_frame_gap = validation_frame_gap
        self._image_filename = image_filename
        self._dark_filename = dark_filename
        self._params_filename = params_filename
        self._output_directory = output_directory
        self._log_filename = log_filename
        self._run_name = run_name
        self._use_C = use_C

        assert self._validate_input_params() == True, self._show_usage()


    def _validate_input_params(self):

        if self._output_directory == '':
            print ('Output Directory cannot be empty')
            return False

        if self._image_filename == '':
            print ('Image filename cannot be empty')
            return False

        if self._mode == 'rc' and self._dark_filename == '':
            print ('Dark filename cannot be empty when mode is rc')
            return False

        if self._mode == 'rc' and self._params_filename == '':
            print ('Parameter filename cannot be empty when mode is rc')
            return False

        if self._verbosity > 2:
            self._verbosity = 2

        if self._verbosity < 0:
            self._verbosity = 0

        return True

    @property
    def mode(self):
        """Returns the run mode: rc or de
        """
        return self._mode

    @property
    def verbosity(self):
        """Returns the verbosity level: 0-2
        """
        return self._verbosity

    @property
    def validation_frame_gap(self):
        """Returns the gap between validation frames
        """
        return self._validation_frame_gap

    @property
    def image_filename(self):
        """Returns the path to the file to be processed
        """
        return self._image_filename

    @property
    def dark_filename(self):
        """Returns the path to the file containing calibration data
        """
        return self._dark_filename

    @property
    def params_filename(self):
        """Returns the path to input params file
        """
        return self._params_filename

    @property
    def output_directory(self):
        """Returns the output directory
        """
        return self._output_directory

    @property
    def log_filename(self):
        """Returns the name of the log file
        """
        return self._log_filename

    @property
    def run_name(self):
        """Returns the name used to identify this run in the log file
        """
        return self._run_name

    @property
    def use_C(self):
        """Returns indicator showing if the optimized C implementation will be used
        """
        return self._use_C

    def _show_usage(self):
        print("Usage:\n")
        print("ReCoDe -rc -i ARG -d ARG -p ARG -o ARG [-v ARG] [-?] [--help]")
        print("ReCoDe -de -i ARG -o ARG [-v ARG] [-?] [--help]")
        print("")
        print("-rc:    Perform Reduction-Compression (Either -rc or -de must be specified)")
        print("-de:    Perform Decompression-Expansion (Either -rc or -de must be specified)")
        print("-i:     (Required) Image file to be compressed when using -rc and ReCoDe file to be decompressed when using -de")
        print("-o:     (Required) Output directory")
        print("-d:     Dark file (Required when using -rc)")
        print("-p:     Params file (Required when using -rc)")
        print("-v:     Verbosity level (0 or 1, Optional)")
        print("-l:     Log file name (Optional)")
        print("-n:     Run name (Optional). Used while logging.")
        print("-vf:    Gap between validation frames. (Optional). If not specified no validation frames are saved.")
        print("-h:     Displays this help (Optional)")
        print("-help:  Displays this help (Optional)")
        print("--help: Displays this help (Optional)")
        print("-?:     Displays this help (Optional)")


class InputParams():

    def __init__(self, params_filename):
        self._params_filename = params_filename

    def load(self):
        # make all non filename strings lowercase with .lower()
        self._reduction_level = 0
        self._rc_operation_mode = 0
        self._dark_threshold_epsilon = 0
        self._bit_depth = 0
        self._source_bit_depth = 0
        self._num_cols = 0
        self._num_rows = 0
        self._num_frames = 0
        self._frame_offset = 0
        self._num_dark_frames = 0
        self._dark_frame_offset = 0
        self._keep_part_files = 0
        self._num_threads = 0
        self._L2_statistics = 0
        self._L4_centroiding = 0
        self._compression_scheme = 0
        self._compression_level = 0
        self._keep_dark_data = 0
        self._source_file_type = 0
        self._source_header_length = 0
        self._dark_file_type = 0
        self._dark_header_length = 0