import numpy as np
import mrcfile
import pims

def _read_data_from_file(init_params, input_params):

    if input_params._source_file_type == 'mrc' or input_params._source_file_type == 'mrcs':
        mrc = mrcfile.open(init_params._source_file_name, mode='r')
        return mrc.header, mrc.data

    self._init_params._source_file_name, self._input_params._source_file_type, self._input_params._source_bit_depth
    num_bytes = 0
    return np.array((num_bytes)), num_bytes

def _read_bytearray_from_stream(_iostream, number_of_bytes):
    """Read a 'bytearray' from the stream.
    This default implementation relies on the stream implementing the
    'io.BufferedIOBase.readinto' method to avoid copying the new
    array while creating the mutable 'bytearray'. Subclasses
    should override this if their stream does not support
    'io.BufferedIOBase.readinto'
    Returns:
        A 2-tuple of the 'bytearray' and the number of bytes that
        were read from the stream.
    """
    result_array = bytearray(number_of_bytes)
    bytes_read = _iostream.readinto(result_array)
    return result_array, bytes_read