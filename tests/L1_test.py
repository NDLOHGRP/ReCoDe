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

def match_frames(a,b):
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

def match_L4(L4_decoded, L4_counted_frames):
    for a,b in zip(L4_decoded, L4_counted_frames):
        status = match_frames(a,b)
        if status != 0:
            return -1
        else:
            continue
    return 0

def test_L4(L4_counted_frames):
    status = run('../win32/Debug/ReCoDe -rc -i ../data/Frames_0_3_12-13-59.436.bin -d ../data/Dark_Frame_12-23-00.232.bin -p ../config/recode_params_bd.txt -o ../scratch')
    if status == 0:
        status = run('../win32/Debug/ReCoDe -de -i ../scratch/Frames_0_3_12-13-59.436.bin.rc4 -o ../scratch')
        if status == 0:
            reader = BinaryDataReader([512,4096], 2, np.uint16)
            L4_decoded = reader.extract_binary_frames('../scratch/Recoded_Frames_0_3_12-13-59.436.bin', 0, 4)
            np.savetxt('../scratch/Recoded_Frames_0_3_12-13-59.436.txt', L4_decoded[0], fmt='%d')
            np.savetxt('../scratch/Recoded_Frames_1_3_12-13-59.436.txt', L4_decoded[1], fmt='%d')
            np.savetxt('../scratch/Recoded_Frames_2_3_12-13-59.436.txt', L4_decoded[2], fmt='%d')
            np.savetxt('../scratch/Recoded_Frames_3_3_12-13-59.436.txt', L4_decoded[3], fmt='%d')
            return match_L4 (L4_decoded, L4_counted_frames)
        else:
            return -1
    else:
        return -1

def reduce_L4(data, dark, threshold, max_frames=10):
    """
    Extracts features of secondary electron puddles
    Returns puddle size, frame index, centroid y and centroid x of detected secondary electron puddles.
    Parameters
    ----------
    data : numpy array of shape [frames, rows, columns]
    dark : numpy array of shape [rows, columns]
    threshold : used to separate signal from noise
        
    Returns
    -------
    feature_dict : dictionary with frame indices as keys and list of N 1x4 arrays as value. The N 1x4 arrays are [centroid y, centroid x, puddle size, mean intensity] for the N secondary electron puddles
    """
    print(data.shape)
    n_Frames = data.shape[0]
    if max_frames < 0 or max_frames > n_Frames:
        max_frames = n_Frames

    data = data.astype(np.int32)
    dark = dark.astype(np.int32)

    feature_dict = {}
    s = nd.generate_binary_structure(2,2)
    counted_frames = np.zeros((n_Frames, data.shape[1], data.shape[2]))
    for frame_index in range(max_frames):
        dark_subtracted = data[frame_index] - dark
        dark_subtracted_binary = dark_subtracted > threshold
        labeled_foreground, num_features = nd.measurements.label(dark_subtracted_binary, structure=s)
        print('num_features =', num_features)
        regions = regionprops(labeled_foreground, dark_subtracted)
        centroids = [reg['weighted_centroid'] for reg in regions]
        centroids = (np.round(centroids)).astype(np.uint16)
        C = np.swapaxes(centroids,0,1)
        indices = np.ravel_multi_index(C, [512, 4096])
        np.put(counted_frames[frame_index], indices, 1)
        sm = coo_matrix( ( np.ones(centroids.shape[0]), (centroids[:,0],centroids[:,1]) ), shape=(512, 4096))
        feature_dict[frame_index] = sm
    return counted_frames, feature_dict
    
if __name__== "__main__":

    reader = BinaryDataReader([512,4096], 2, np.uint16)
    data = reader.extract_binary_frames('../data/Frames_0_3_12-13-59.436.bin', 0, 4)
    # data = reader.extract_binary_frame('../data/Frame_0_12-13-59.436.bin', 0)
    dark = reader.extract_binary_frame('../data/Dark_Frame_12-23-00.232.bin', 0)
    L4_counted_frames, _ = reduce_L4(data, dark[0], 0, max_frames=-1)
    np.savetxt('../scratch/Counted_Frame_2.txt', L4_counted_frames[2], fmt='%d')
    status = test_L4(L4_counted_frames)
    if status == 0:
        print("Passed")
    else:
        print("Failed")
    