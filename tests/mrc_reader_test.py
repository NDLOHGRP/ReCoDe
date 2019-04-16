import numpy as np
import mrcfile

t_data = np.loadtxt('D:\cbis\GitHub\ReCoDe\win32\Release\dark_frame.txt', dtype=np.uint16, delimiter=' ')
nx = np.shape(t_data)[0]
ny = np.shape(t_data)[1]

data =np.zeros((nx,ny), dtype=np.uint16)

count = 0
mrcfilename = 'Z:/Benedikt_Data/20181126_testing/20181126_06240_dark_RawImages.mrc'
with mrcfile.open(mrcfilename, permissive=True) as mrc:
    data = mrc.data[0,:,:]

print(np.shape(data))
diff = t_data - data
print(np.shape(diff))
print(np.array_equal(data, t_data))
np.savetxt('../scratch/diff.txt', diff, delimiter='\t', fmt='%d')