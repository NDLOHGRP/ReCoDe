import numpy as np
from scipy.sparse import csr_matrix
from ctypes import *

class ReCoDe():

    def __init__(self, network_model_dir='', network_weights_dir=''):
        self.recode = cdll.LoadLibrary("D:\cbis\GitHub\ReCoDe\win32\Release\ReCoDe.dll")
        print(self.recode)

    def reduce_compress(self, image_filename, dark_filename, params_filename, output_directory='/', verbosity=0):
        a = self.recode.py_rc(image_filename, dark_filename, params_filename, output_directory, verbosity)
        print(a)
        
    def decompress_expand(self, recode_filename, output_directory='/', sparse=True, verbosity=0):
        self.recode.py_de.restype = c_void_p
        a = self.recode.py_de(recode_filename, output_directory, sparse, verbosity)
        ptr = cast(a, POINTER(c_ushort))

        print("C returned: ")
        row = np.zeros((ptr[0]), dtype=np.uint16)
        col = np.zeros((ptr[0]), dtype=np.uint16)
        image_nx = ptr[1]
        image_ny = ptr[2]
        for i in range(ptr[0]):
            row[i] = ptr[i*2+3]
            col[i] = ptr[i*2+4]
        
        data = np.ones((ptr[0]), dtype=np.uint16)
        sparse_frame = csr_matrix((data, (row, col)), shape=(image_nx, image_ny)).toarray()
        
        print (sparse_frame)
        return sparse_frame
            
            
if __name__ == "__main__":
    
    recoder = ReCoDe()
    # recoder.reduce_compress("","","")
    recoder.decompress_expand("")