import numpy as np
import os
import sys
from os.path import isfile, join
import tifffile as tiff
import scipy.ndimage as nd
import pims
from skimage.measure import regionprops
from scipy.sparse import coo_matrix
import mrcfile
import subprocess
from sklearn.neighbors import NearestNeighbors
from datetime import datetime
import time

sys.path.append('../src')
from pyrecode.recoder import reduce_L1, reduce_L2, reduce_L3, reduce_L4
from pyrecode.params import InputParams

class BinaryDataReader:
    def __init__(self, shape, bytes_per_pixel, dtype):
        self.shape = shape
        self.bytes_per_pixel = bytes_per_pixel
        self.dtype = dtype
        self.nBytesPerFrame = self.shape[0]*self.shape[1]*self.bytes_per_pixel

    def extract_binary_frame(self, binary_filename, frame_index):
        s1 = frame_index*self.nBytesPerFrame
        f = open(binary_filename, "rb")
        f.seek(s1)
        image_bytes = f.read(self.nBytesPerFrame)
        f.close()
        flattened_frame = np.frombuffer(image_bytes, dtype=self.dtype)
        frame = np.zeros((1, self.shape[0], self.shape[1]))
        frame[0] = flattened_frame.reshape((self.shape[0], self.shape[1]))	
        return frame

    def extract_binary_frames (self, binary_filename, frame_start_index, num_frames):
        s1 = frame_start_index*self.nBytesPerFrame
        f = open(binary_filename, "rb")
        f.seek(s1)
        frames = np.zeros((num_frames, self.shape[0], self.shape[1]))
        for i in range(num_frames):
            image_bytes = f.read(self.nBytesPerFrame)
            flattened_frame = np.frombuffer(image_bytes, dtype=self.dtype)
            frames[i] = flattened_frame.reshape((self.shape[0], self.shape[1]))
        f.close()
        return frames

def run(cmd_str):
    cmd_str = str(os.path.normpath(cmd_str))
    process = subprocess.Popen(cmd_str.split(), stdout=subprocess.PIPE)
    result = []
    for line in process.stdout:
        result.append(line)
    errcode = process.returncode
    for line in result:
        print(line)
    if errcode is not None:
        return -1
    else:
        return 0

def _match_frames_L1 (a,b):
    if np.array_equal (a,b):
        return 0
    else:
        return -1

def _match_frames_L4 (a,b):
    (y1,x1) = np.nonzero(a[10:-10,10:-10])
    (y2,x2) = np.nonzero(b[10:-10,10:-10])
    nbrs = NearestNeighbors(n_neighbors=1, algorithm='auto').fit(np.stack((y1,x1),1))
    distances, indices = nbrs.kneighbors(np.stack((y2,x2),1))
    cnt = np.sum(distances > np.sqrt(2))
    print('mismatch count =', str(cnt))
    if cnt == 0:
        return 0
    else:
        return -1

def _match (_decoded_frames, py_reduced_frames, reduction_level):
    for a,b in zip(_decoded_frames, py_reduced_frames):
        if reduction_level == 1:
            status = _match_frames_L1 (a,b)
        elif reduction_level == 4:
            status = _match_frames_L4 (a,b)
        if status != 0:
            return -1
        else:
            continue
    return 0

def _test (py_reduced_frames, reduction_level, test_dir):
    # make params file for reduction_level
    ip = InputParams()
    ip.load('../config/recode_params_3.txt')
    ip.reduction_level = reduction_level
    ip.serialize('../config/recode_params_test.txt')

    status = run('../win32/x64/Release/ReCoDe -rc -i ../data/Frames_0_3_12-13-59.436.bin -d ../data/Dark_Frame_12-23-00.232.bin -p ../config/recode_params_test.txt -o ../scratch/' + test_dir)
    if status == 0:
        if reduction_level == 1:
            status = run('../win32/x64/Release/ReCoDe -de -i ../scratch/' + test_dir + '/Frames_0_3_12-13-59.436.bin.rc1 -o ../scratch/' + test_dir)
        elif reduction_level == 4:
            status = run('../win32/x64/Release/ReCoDe -de -i ../scratch/' + test_dir + '/Frames_0_3_12-13-59.436.bin.rc4 -o ../scratch/' + test_dir)
        if status == 0:
            reader = BinaryDataReader([512,4096], 2, np.uint16)
            _decoded_frames = reader.extract_binary_frames('../scratch/' + test_dir + '/Recoded_Frames_0_3_12-13-59.436.bin', 0, 4)
            np.savetxt('../scratch/' + test_dir + '/Recoded_Frames_0_3_12-13-59.436.txt', _decoded_frames[0], fmt='%d')
            np.savetxt('../scratch/' + test_dir + '/Recoded_Frames_1_3_12-13-59.436.txt', _decoded_frames[1], fmt='%d')
            np.savetxt('../scratch/' + test_dir + '/Recoded_Frames_2_3_12-13-59.436.txt', _decoded_frames[2], fmt='%d')
            np.savetxt('../scratch/' + test_dir + '/Recoded_Frames_3_3_12-13-59.436.txt', _decoded_frames[3], fmt='%d')
            return _match (_decoded_frames, py_reduced_frames, reduction_level)
        else:
            return -1
    else:
        return -1

def _on_the_fly_test (py_reduced_frames, reduction_level, test_dir):
    # make params file for reduction_level
    ip = InputParams()
    ip.load('../config/recode_params_3.txt')
    ip.reduction_level = reduction_level
    ip.serialize('../config/recode_params_test.txt')

    # status = run('../win32/x64/Release/OnTheFlyReCoDe -test -i ../data/Frames_0_3_12-13-59.436.bin -d ../data/Dark_Frame_12-23-00.232.bin -p ../config/recode_params_test.txt -o ../scratch -l recode.log -n run1 -v 1 -vf -1')
    status = run('../win32/x64/Release/OnTheFlyReCoDe -test -i ../data/large_files/17-21-33.138_Dark_Ref_3.seq -d ../scratch/temp/003.seq -o ../scratch/'+test_dir+' -p ../config/recode_params.txt -l recode.log -n _run1 -v 1 -vf 100')
    time.sleep(2)
    return -1
    if status == 0:
        status = run('../win32/x64/Release/ReCoDe -merge -i ../scratch/'+test_dir+'/_run1.rc1_part000 -o ../scratch/'+test_dir+' -p ../config/recode_params.txt -n _run1')
        if status == 0:
            if reduction_level == 1:
                status = run('../win32/x64/Release/ReCoDe -de -i ../scratch/'+test_dir+'/17-21-33.138_Dark_Ref_3.seq.rc1 -o ../scratch/'+test_dir)
            elif reduction_level == 4:
                status = run('../win32/x64/Release/ReCoDe -de -i ../scratch/Frames_0_3_12-13-59.436.bin.rc4 -o ../scratch')
            if status == 0:
                reader = BinaryDataReader([512,4096], 2, np.uint16)
                _decoded_frames = reader.extract_binary_frames('../scratch/Recoded_17-21-33.138_Dark_Ref_3.bin', 0, 4)
                np.savetxt('../scratch/'+test_dir+'/Recoded_Frames_0_3_12-13-59.436.txt', _decoded_frames[0], fmt='%d')
                np.savetxt('../scratch/'+test_dir+'/Recoded_Frames_1_3_12-13-59.436.txt', _decoded_frames[1], fmt='%d')
                np.savetxt('../scratch/'+test_dir+'/Recoded_Frames_2_3_12-13-59.436.txt', _decoded_frames[2], fmt='%d')
                np.savetxt('../scratch/'+test_dir+'/Recoded_Frames_3_3_12-13-59.436.txt', _decoded_frames[3], fmt='%d')
                return _match (_decoded_frames, py_reduced_frames, reduction_level)
            else:
                return -1
        else:
            print ('Merge failed')
            return -1
    else:
        print ('Simulated on-the-fly failed')
        return -1
    
def _test_1 ():

    now = datetime.now()
    # timestamp = str(datetime.timestamp(now))
    timestamp = now.strftime("%m_%d_%Y-%H_%M_%S")

    _test_dir = os.path.join("../scratch", timestamp)
    if os.path.exists(_test_dir):
        for files in os.listdir(_test_dir):
            os.remove(files)
    else:
        os.mkdir(_test_dir)

    reader = BinaryDataReader([512,4096], 2, np.uint16)
    data = reader.extract_binary_frames('../data/Frames_0_3_12-13-59.436.bin', 0, 4)
    dark = reader.extract_binary_frame('../data/Dark_Frame_12-23-00.232.bin', 0)

    for reduction_level in [1]:
        
        if reduction_level == 1:
            py_reduced_frames, _ = reduce_L1 (data, dark[0], 0, max_frames=-1)
        elif reduction_level == 2:
            py_reduced_frames, _ = reduce_L2 (data, dark[0], 0, max_frames=-1)
        elif reduction_level == 3:
            py_reduced_frames, _ = reduce_L3 (data, dark[0], 0, max_frames=-1)
        elif reduction_level == 4:
            py_reduced_frames, _ = reduce_L4 (data, dark[0], 0, max_frames=-1)

        np.savetxt('../scratch/Counted_Frame_3.txt', py_reduced_frames[3], fmt='%d')

        status = _test (py_reduced_frames, reduction_level, timestamp)
        if status == 0:
            print("Passed")
        else:
            print("Failed")
    

if __name__== "__main__":

    # _test_1 ()
    # exit(0)

    now = datetime.now()
    # timestamp = str(datetime.timestamp(now))
    timestamp = now.strftime("%m_%d_%Y-%H_%M_%S")

    _test_dir = os.path.join("../scratch", timestamp)
    if os.path.exists(_test_dir):
        for files in os.listdir(_test_dir):
            os.remove(files)
    else:
        os.mkdir(_test_dir)
    test_dir = _test_dir

    dataset = pims.open('../data/large_files/17-21-33.138_Dark_Ref_3.seq')
    darkset = pims.open('../data/large_files/17-21-33.138_Dark_Ref_3.seq')
    data = np.asarray(dataset[0:4])
    dark = np.asarray(darkset[0:1])

    py_reduced_frames, _ = reduce_L1 (data, dark[0], 20, max_frames=-1)
    np.savetxt('../scratch/'+test_dir+'/Counted_Frame_3.txt', py_reduced_frames[3], fmt='%d')

    # status = _on_the_fly_test (py_reduced_frames, reduction_level, timestamp)

    reader = BinaryDataReader([512,4096], 2, np.uint16)
    _decoded_frames = reader.extract_binary_frames('../scratch/Recoded_17-21-33.138_Dark_Ref_3.bin', 0, 4)
    np.savetxt('../scratch/'+test_dir+'/Recoded_Frames_0_3_12-13-59.436.txt', _decoded_frames[0], fmt='%d')
    np.savetxt('../scratch/'+test_dir+'/Recoded_Frames_1_3_12-13-59.436.txt', _decoded_frames[1], fmt='%d')
    np.savetxt('../scratch/'+test_dir+'/Recoded_Frames_2_3_12-13-59.436.txt', _decoded_frames[2], fmt='%d')
    np.savetxt('../scratch/'+test_dir+'/Recoded_Frames_3_3_12-13-59.436.txt', _decoded_frames[3], fmt='%d')
    status = _match (_decoded_frames, py_reduced_frames, 1)

    if status == 0:
        print("Passed")
    else:
        print("Failed")