import numpy as np

class ReCoDeReader():
    def __init__(self, file, source_type):
        try:
            import mrcfile
            self._mrcfile = mrcfile
        except ImportError:
            print("Reading MRC files requires mrcfile to be installed")
            return
        # EMReaderBase.__init__(self, file, source_type, False)

    def _open(self):
        try:
            self._file_handle = self._mrcfile.open(self._source_filename, mode="r")
        except ValueError:
            self._file_handle = self._mrcfile.open(self._source_filename, mode="r", permissive=True)
        self._stack = self._file_handle.data
        
    def _load_header(self):
        record = self._file_handle.header
        names = record.dtype.names
        return dict(zip(names, [record[item] for item in names]))

    def get_true_shape(self):
        return self._stack.shape

    def _get_shape(self):
        return tuple([self._header["nz"], self._header["ny"], self._header["nx"]])

    def _get_dtype(self):
        return (self._stack.dtype)

    def _get_sub_volume(self, slice_z, slice_y, slice_x):
        if self._file_handle.is_image_stack():
            container = self._stack[slice_z, slice_y, slice_x]
        elif self._file_handle.is_single_image():
            container = self._stack[np.newaxis,slice_y,slice_x]
        else:
            raise NotImplementedError
        return io.emcreate(container)

    def _get_frame(self, z_index):
        if self._file_handle.is_single_image():
            container = self._stack[np.newaxis,:,:]
        elif self._file_handle.is_image_stack():
            container = self._stack[z_index][np.newaxis,:,:]
        else:
            raise NotImplementedError
        return io.emcreate(container)

    def close(self):
        self._file_handle.close()

class MRCReaderCorrupted(MRCReader):
    def __init__(self, *args, nz=1, ny=1024, nx=1024, dtype='float32', byteorder='>'):
        self.nz = nz
        self.ny = ny 
        self.nx = nx
        self._dtype = np.dtype(dtype)
        self._shape = (self.nz, self.ny, self.nx)
        self.byteorder = byteorder
        MRCReader.__init__(self, *args)

    def _open(self):
        self._file_handle = self._mrcfile.open(self._source_filename, permissive=True, mode="r")
        dt = self.dtype
        dt = dt.newbyteorder(self.byteorder)
        nbytes = dt.itemsize*np.prod(self.shape)
        data_bytes = self._file_handle._iostream.read(nbytes)
        self._stack = np.frombuffer(data_bytes, dtype=dt).reshape(self.shape)
    
    def _get_shape(self):
        return (self.nz, self.ny, self.nx)
