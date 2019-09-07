'''
Usage:
"D:\Program Files\Python\Python37\python_d" pyrecode_c
'''

import sys
import numpy as np
# sys.path.append('../../win32/x64/Debug')
# sys.path.append('../../win32/x64/Release')
sys.path.append('D:/cbis/GitHub/ReCoDe/win32/x64/Release')

for p in sys.path:
    print(p)

import PyReCoDe
from PyReCoDe import _compress_stream, _decompress_stream

print( dir(PyReCoDe) )
print(PyReCoDe.__doc__)

for r in range(10):
    byte_array = bytearray(1024)
    for i in range(1024):
        byte_array[i] = i%(r+2)
    data = memoryview(byte_array)
    compressed_data = memoryview(bytes(1024))
    decompressed_data = memoryview(bytes(1024))
    n_compressed_bytes = _compress_stream(0,1,data,1024,compressed_data)
    print("n_compressed_bytes =", n_compressed_bytes)
    _decompress_stream(0,1,compressed_data,decompressed_data,n_compressed_bytes,1024)
    a1 = np.asarray(data)
    a2 = np.asarray(decompressed_data)
    state = np.array_equal(a1, a2)
    print(r, state)

'''
for d in decompressed_data:
    print(d)
'''

# data.release()
# compressed_data.release()
# decompressed_data.release()
