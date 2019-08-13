import numpy as np

from emio import emfile
from .params import InputParams
from .recode_utils import _read_data_from_file, _read_bytearray_from_stream

class ReCoDeHeader():

    def __init__(self, init_params, input_params):
        self._uid = 158966344846346
        self._version_major = 0
        self._version_minor = 1
        self._reduction_level = input_params.reduction_level
        self._recode_operation_mode = input_params.recode_operation_mode
        self._bit_depth = input_params.bit_depth
        self._nx = input_params.nx
        self._ny = input_params.ny
        self._nz = input_params.nz
        self._L2_statistics = input_params.L2_statistics
        self._L4_centroiding = input_params.L4_centroiding
        self._compression_scheme = input_params.compression_scheme
        self._compression_level = input_params.compression_level
        self._source_file_type = input_params.source_file_type
        self._source_header_length = input_params.source_header_length
        self._source_header_position = 0
        self._source_file_name = init_params.image_filename
        self._dark_file_name = init_params.dark_filename
        self._dark_threshold_epsilon = input_params.dark_threshold_epsilon
        self._has_dark_data = input_params.has_dark_data
        self._frame_offset = input_params.frame_offset
        self._dark_frame_offset = input_params.dark_frame_offset
        self._num_dark_frames = input_params.num_dark_frames
        self._source_bit_depth = input_params.source_bit_depth
        self._futures = np.array((255), dtype=np.uint8)

    def validate():
        # check that no fields are left empty
        return True