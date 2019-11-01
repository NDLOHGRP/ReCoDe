import numpy as np
import PyReCoDe
from .recode_header import ReCoDeHeader

class ReCoDeReader():
    def __init__(self, file):
        self._source_filename = file
        self._header = {}
        self._seek_table = {}
        self._rc_header = None
        self._c_reader = PyReCoDe.Reader()

    def _open(self):
        self._c_reader._open_file(self._source_filename)
        
    def _load_header(self):
        self._rc_header = ReCoDeHeader()
        self._rc_header.load(self._source_filename)
        self._rc_header.print()
        self._header = self._rc_header._rc_header
        return self._rc_header._rc_header

    def _load_seek_table(self):
        assert 'nz' in self._header, print('Attempting to read seek table before reading header')
        _nz = int(self._header['nz'])
        self._seek_table = np.zeros((_nz,3), dtype=np.uint32)
        # n_compressed_bytes_in_binary_image, n_compressed_bytes_in_pixvals, n_bytes_in_packed_pixvals
        with open(self._source_filename, "rb") as f:
            self._rc_header.skip_header(f)
            for _frame_index in range(_nz):
                p = f.read(12)
                if not p:
                    return
                self._seek_table[_frame_index] = np.frombuffer(p, dtype=np.uint32)
        a = np.sum(self._seek_table, axis=1)
        b = np.expand_dims(np.cumsum(a), axis=1)
        print(self._seek_table.shape, b.shape)
        self._seek_table = np.append(self._seek_table, b, axis=1)
        print(self._seek_table)

    def get_true_shape(self):
        return tuple([self._header["nz"], self._header["ny"], self._header["nx"]])

    def _get_shape(self):
        return tuple([self._header["nz"], self._header["ny"], self._header["nx"]])

    def _get_dtype(self):
        return (self._stack.dtype)

    def _get_sub_volume(self, slice_z, slice_y, slice_x):
        raise NotImplementedError

    def _get_frame(self, z_index):
        self._c_reader._fseek(int(self._rc_header._rc_header_length + self._seek_table[z_index,3]), int(0))
        byte_array = bytearray(self._header["ny"]*self._header["nx"]*3*2)
        frame_data = memoryview(byte_array)
        self._c_reader._get_frame_sparse_L1(self._seek_table[z_index,0], self._seek_table[z_index,1], self._seek_table[z_index,2], frame_data)
        d = np.asarray(frame_data)
        print(d)
        return d

    def _close(self):
        self._c_reader._close_file()


if __name__== "__main__":

    reader = ReCoDeReader('D:/cbis/GitHub/ReCoDe/scratch/400fps_dose_43.rc1')
    reader._open()
    reader._load_header()
    reader._load_seek_table()
    reader._get_frame(3)
    reader._close()

