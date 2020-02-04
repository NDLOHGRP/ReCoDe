import numpy as np
import os
import sys
from scipy.sparse import coo_matrix
import scipy.ndimage as nd
from skimage.measure import regionprops
from numba import jit
from datetime import datetime

'''
To Do: 
1. centroiding scheme option for _get_centroids_2D and implementing 'max' centroiding scheme. For 'unweighted' centroiding scheme set frame = b_frame
2. multithreading support
3. support for RecodeReader object in addition to dictionary of frames
'''
def L1_to_L4 (L1_frames, calibration_frame=None, n_Frames=-1):

    if n_Frames < 1:
        n_Frames = len(L1_frames)
    
    fid = list(L1_frames.keys())[0]
    sz = L1_frames[fid].shape
    _n_pixels_in_frame = sz[0]*sz[1]*1.
    t = np.ones(int(_n_pixels_in_frame), dtype=np.bool)
    
    if calibration_frame is None:
        calibration_frame = np.zeros((sz[0], sz[1]), dtype=np.uint16)
    
    s = nd.generate_binary_structure(2,2)
    cd = {}

    start=datetime.now()
    for key in range(n_Frames):
        if key in L1_frames:
            
            s1=datetime.now()
            frame = L1_frames[key].todense()
            dark_subtracted_binary = frame > calibration_frame
            labeled_foreground, num_features = nd.measurements.label(dark_subtracted_binary, structure=s)
            sc=datetime.now()
            # centroids = _get_centroids_2D (labeled_foreground, dark_subtracted_binary, frame)
            centroids = _get_centroids_2D_nb (labeled_foreground, dark_subtracted_binary, frame)
            tc = datetime.now()-sc
            centroids = (np.round(centroids)).astype(np.uint16)
            cd[key] = coo_matrix((t[:len(centroids)], (centroids[:,1], centroids[:,0])), shape=(sz[0], sz[1]), dtype=np.bool)
            t1 = datetime.now()-s1

            print(key)
            print('Dose Rate =', num_features/(_n_pixels_in_frame))
            print('Processing time: ' + str(t1) + ', ' + str(tc))
    print ('Total processing time: ' + str(datetime.now()-start))
    return cd


def _get_centroids_2D(labelled_image, b_frame, frame):
    
    n_cols = b_frame.shape[1]
    _pixels = np.argwhere(b_frame)

    centroids = {}
    for p in _pixels:
        L = labelled_image[p[0], p[1]]
        v = (frame[p[0], p[1]])*1.
        if L in centroids:
            centroids[L]['x'] += v*p[0]
            centroids[L]['y'] += v*p[1]
            centroids[L]['w'] += v
        else:
            centroids[L] = {'x':v*p[0], 'y':v*p[1], 'w':v}
            
    centroid_arr = np.zeros((len(centroids),2))
    for i,c in enumerate(centroids):
        centroid_arr[i,0] = centroids[c]['x'] / centroids[c]['w']
        centroid_arr[i,1] = centroids[c]['y'] / centroids[c]['w']

    return centroid_arr


@jit(nopython=True)
def _get_centroids_2D_nb (labelled_image, b_frame, frame):
    
    centroids = dict()
    centroids[0] = np.zeros(3, dtype=np.float32)

    n_rows = b_frame.shape[0]
    n_cols = b_frame.shape[1]

    for r in range(n_rows):
        for c in range(n_cols):
            if b_frame[r][c]:
                label = labelled_image[r][c]
                v = frame[r][c]*1.
                if label in centroids:
                    centroids[label][0] += v*r
                    centroids[label][1] += v*c
                    centroids[label][2] += v
                else:
                    centroids[label] = np.zeros(3, dtype=np.float32)
                    centroids[label][0] = v*r
                    centroids[label][1] = v*c
                    centroids[label][2] = v
            
    centroid_arr = np.zeros((len(centroids)-1,2))
    for i,label in enumerate(centroids):
        if label > 0:
            centroid_arr[i-1,0] = centroids[label][0] / centroids[label][2]
            centroid_arr[i-1,1] = centroids[label][1] / centroids[label][2]

    return centroid_arr



if __name__== "__main__":

    test_frame = [
    [0,0,0,0,0,0,0,1,1],
    [0,1,1,0,0,0,0,0,0],
    [0,1,1,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0],
    [1,0,0,1,0,0,0,0,0],
    [0,1,0,0,1,0,0,0,0],
    [0,0,1,0,0,0,0,0,0],
    [0,1,1,0,0,0,0,0,1],
    [1,0,0,0,0,0,0,1,1]]

    test_b_frame = [
    [0,0,0,0,0,0,0,1,1],
    [0,1,1,0,0,0,0,0,0],
    [0,1,1,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0],
    [1,0,0,1,0,0,0,0,0],
    [0,1,0,0,1,0,0,0,0],
    [0,0,1,0,0,0,0,0,0],
    [0,1,1,0,0,0,0,0,1],
    [1,0,0,0,0,0,0,1,1]]

    labelled_image = [
    [0,0,0,0,0,0,0,1,1],
    [0,2,2,0,0,0,0,0,0],
    [0,2,2,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0],
    [3,0,0,4,0,0,0,0,0],
    [0,3,0,0,4,0,0,0,0],
    [0,0,3,0,0,0,0,0,0],
    [0,3,3,0,0,0,0,0,5],
    [3,0,0,0,0,0,0,5,5]]

    c = _get_centroids_2D_nb (
        np.array(labelled_image, dtype=np.uint16), 
        np.array(test_b_frame, dtype=np.bool), 
        np.array(test_frame, dtype=np.uint16)
    )
    start = time.time()
    for i in range(500):
        c = _get_centroids_2D_nb (
            np.array(labelled_image, dtype=np.uint16), 
            np.array(test_b_frame, dtype=np.bool), 
            np.array(test_frame, dtype=np.uint16)
        )
    end = time.time()
    print("Elapsed (after compilation) = %s" % (end - start))
    print(c)