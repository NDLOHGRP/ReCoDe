import numpy as np
import os
import sys
from scipy.sparse import coo_matrix
import tifffile as tiff
import scipy.ndimage as nd
from skimage.measure import regionprops
from datetime import datetime
import numba
from numba import jit
from numba import njit
from numba import types
from numba.typed import Dict

def ff_cca (binary_frame):
    bf = np.pad(binary_frame, ((1,1),(1,1)), mode='constant')
    n_cols = bf.shape[1]
    states = np.array(bf).flatten()
    fg_pixels = np.argwhere(bf)
    CCs = {}
    neighbors = {}
    for fgp in fg_pixels:
        ni = {}
        for ij in [[-1,-1],[-1,0],[-1,1],[0,1],[1,1],[1,0],[1,-1],[0,-1]]:
            x = fgp[0]-ij[0]
            y = fgp[1]-ij[1]
            li = x*n_cols+y
            if states[li]:
                ni[li] = [x,y]
        if len(ni) > 0:
            neighbors[fgp[0]*n_cols + fgp[1]] = ni
    # print(neighbors)

    for fgp in fg_pixels:
        li = fgp[0]*n_cols + fgp[1]
        if states[li]:
            cc = {}
            cc[li] = 1
            states[li] = 0
            if li in neighbors:
                for n in neighbors[li]:
                    if states[n]:
                        cc[n] = 1
                        states[n] = 0
            CCs[li] = cc

    '''
    for key in CCs:
        if len(CCs[key]) > 1:
            print(CCs[key])
    '''


def compute_gain_L1 (L1_frames, n_Frames=-1):

    '''
    ToDo: 
    1. custom frame shape
    2. histogram binning for different datatypes
        a. For unsigned integers max no. of bins is 4096, if data has range greater than that, then truncate, else use the smaller range
        b. Special handling is required for signed integers
        c. Floats will always be binned to a specific range, say 4096
    '''

    if n_Frames < 1:
        n_Frames = len(L1_frames)
    print (n_Frames)
    frame = np.zeros((512, 4096), dtype=np.uint16)
    gain_thresholds = np.zeros((512, 4096), dtype=np.uint16)
    s = nd.generate_binary_structure(2,2)
    dose_rate = 0
    start=datetime.now()
    nf = np.min([n_Frames, 100])
    for key in range(nf):
        if key in L1_frames:
            print(key)
            rcv = L1_frames[key]
            frame.fill(0)
            frame[rcv[:,0], rcv[:,1]] =  rcv[:,2]
            dark_subtracted_binary = frame > 0
            labeled_foreground, num_features = nd.measurements.label(dark_subtracted_binary, structure=s)
            dose_rate += num_features/(512.*4096.)
    dose_rate /= (nf*1.)
    n_events = int(np.floor(n_Frames*dose_rate))
    print("dose_rate =", dose_rate, ", n_events =", n_events)
    print("Time to calculate avg. dose rate:", datetime.now()-start)

    # get pixel histograms
    start=datetime.now()
    hists = np.zeros((512*4096, 4096), dtype=np.uint16)
    for key in range(n_Frames):
        if key in L1_frames:
            # print(key)
            rcv = L1_frames[key]
            pixel_indices = np.ravel_multi_index([rcv[:,0], rcv[:,1]], [512, 4096])
            hists[pixel_indices, 4095-rcv[:,2]] += 1
    print("Time to get histograms:", datetime.now()-start)

    start=datetime.now()
    c_hists = np.cumsum(hists, axis=1)
    print(len(c_hists), len(c_hists[0]))
    print(c_hists[500000])
    print("Time to get cumulative histograms:", datetime.now()-start)

    start=datetime.now()
    t = np.argmax(c_hists > n_events, axis=1)
    print(len(t))
    print(c_hists[500000,t[500000]-1], c_hists[500000,t[500000]], c_hists[500000,t[500000]+1])
    print("Time to get thresholds:", datetime.now()-start)

    threshs = 4096 - t
    '''
    for i,g in enumerate(threshs):
        if g != 4096:
            print(i, g, c_hists[i,:], c_hists[i,g-1], c_hists[i,g], c_hists[i,g+1])
    '''
    calibration_frame = np.reshape(threshs, [512,4096])
    print(calibration_frame)
    return calibration_frame

def L1_to_L4_old (L1_frames, calibration_frame, n_Frames=-1):
    if n_Frames < 1:
        n_Frames = len(L1_frames)
    frame = np.zeros((512, 4096), dtype=np.uint16)
    s = nd.generate_binary_structure(2,2)
    cd = {}
    start=datetime.now()
    # for key in range(100):
    for key in range(n_Frames):
        if key in L1_frames:
            print(key)
            rcv = L1_frames[key]
            frame.fill(0)
            frame[rcv[:,0], rcv[:,1]] =  rcv[:,2]
            dark_subtracted_binary = frame > calibration_frame
            labeled_foreground, num_features = nd.measurements.label(dark_subtracted_binary, structure=s)
            regions = regionprops(labeled_foreground, frame)
            centroids = [reg['weighted_centroid'] for reg in regions]
            centroids = (np.round(centroids)).astype(np.uint16)
            cd[key] = centroids
            print('Dose Rate =', num_features/(512.*4096.))
    print (datetime.now()-start)
    return cd

def L1_to_L4 (L1_frames, calibration_frame, n_Frames=-1):
    if n_Frames < 1:
        n_Frames = len(L1_frames)
    frame = np.zeros((512, 4096), dtype=np.uint16)
    s = nd.generate_binary_structure(2,2)
    cd = {}
    start=datetime.now()
    # for key in range(100):
    for key in range(n_Frames):
        if key in L1_frames:
            print(key)
            rcv = L1_frames[key]
            
            s1=datetime.now()
            frame.fill(0)
            t1 = datetime.now()-s1

            s2=datetime.now()
            frame[rcv[:,0], rcv[:,1]] =  rcv[:,2]
            t2 = datetime.now()-s2

            s3=datetime.now()
            dark_subtracted_binary = frame > calibration_frame
            t3 = datetime.now()-s3

            s4=datetime.now()
            labeled_foreground, num_features = nd.measurements.label(dark_subtracted_binary, structure=s)
            t4 = datetime.now()-s4

            s5=datetime.now()
            # regions = regionprops(labeled_foreground, frame)
            # centroids = [reg['weighted_centroid'] for reg in regions]
            centroids = get_centroids_2D (labeled_foreground, dark_subtracted_binary, frame)
            # centroids = (np.round(centroids)).astype(np.uint16)
            t5 = datetime.now()-s5

            cd[key] = centroids
            print('Dose Rate =', num_features/(512.*4096.))
            print(t1, t2, t3, t4, t5)
    print (datetime.now()-start)
    return cd

@jit(nopython=True, parallel=True)
def _median(a):
    print(a.shape)
    b = np.zeros(a.shape[0])
    for i in numba.prange(a.shape[0]):
        b[i] = np.median(a[i])
    return b

@jit(nopython=True, parallel=True)
def _hists(a):
    h = np.zeros((a.shape[0], 4096))
    for p in numba.prange(a.shape[0]):
        for f in range(a.shape[1]):
            h[p,a[p,f]] += 1
    return h

@jit(nopython=True, parallel=True)
def _hists_np(a):
    print(a.shape)
    _bins = np.arange(0,256)
    h = np.zeros((a.shape[0],255))
    for i in numba.prange(a.shape[0]):
        h[i] = np.histogram(a[i], bins=_bins)[0]
    return h


def get_centroids_2D(labelled_image, b_frame, frame):
    
    n_cols = b_frame.shape[1]
    _pixels = np.argwhere(b_frame)

    centroids = {}
    for p in _pixels:
        L = labelled_image[p[0], p[1]]
        v = frame[p[0], p[1]]
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


if __name__== "__main__":

    '''
    a = np.random.randint(0, high=4096, size=(512*4096, 1200), dtype=np.uint16)
    b = _median(a[:5,:5])
    start=datetime.now()
    b = _median(a)
    print (datetime.now()-start)
    sys.exit(0)
    '''

    '''
    a = np.random.randint(0, high=4096, size=(512*4096, 1200), dtype=np.uint16)
    h = _hists(a[:5,:5])
    start=datetime.now()
    h = _hists_np(a)
    print (datetime.now()-start)
    sys.exit(0)
    '''
    # foo()

    '''
    frame = np.random.randint(0, high=100, size=[512,4096])
    b_frame = frame > 98

    s = nd.generate_binary_structure(2,2)

    start=datetime.now()
    labeled_foreground, num_features = nd.measurements.label(b_frame, structure=s)
    props = regionprops(labeled_foreground, frame)
    centroids_1 = [prop.weighted_centroid for prop in props]
    centroids_1 = (np.round(centroids_1)).astype(np.uint16)
    print (datetime.now()-start)

    start=datetime.now()
    labeled_foreground, _ = nd.measurements.label(b_frame, structure=s)
    # centroids = get_centroids(labeled_foreground.flatten(), b_frame, frame.flatten())
    centroids_2 = get_centroids_2D (labeled_foreground, b_frame, frame)
    centroids_2 = (np.round(centroids_2)).astype(np.uint16)
    t = datetime.now()-start
    # print(centroids)
    print (t)
    print('Results Identical =', np.array_equal(centroids_1, centroids_2))

    sys.exit(0)
    '''

    calibration_available = True

    if calibration_available:
        _calibration_frame = np.load('/scratch/loh/abhik/30-Aug-2019/400fps_dose_39_40_gain_ref/400fps_dose_39_40_calibration_frame.npy')
    else:
        f = np.load('/scratch/loh/abhik/30-Aug-2019/400fps_dose_39_40_gain_ref/400fps_dose_39_40.npy')
        _L1_frames_gain_ref = f.item()
        _calibration_frame = compute_gain_L1(_L1_frames_gain_ref, n_Frames=-1)
        np.save('/scratch/loh/abhik/30-Aug-2019/400fps_dose_39_40_gain_ref/400fps_dose_39_40_calibration_frame.npy', _calibration_frame)

    f1 = np.load('/scratch/loh/abhik/30-Aug-2019/400fps_dose_43/400fps_dose_43.npy')
    _L1_frames = f1.item()
    _counted_dictionary = L1_to_L4(_L1_frames, _calibration_frame, n_Frames=100)
    # np.save('/scratch/loh/abhik/30-Aug-2019/400fps_dose_43/counted_2_400fps_dose_39_40.npy', _counted_dictionary)
    sys.exit(0)

    d = np.load('/scratch/loh/abhik/30-Aug-2019/400fps_dose_43/counted_2_400fps_dose_39_40.npy')
    _counted_decoded_frames = d.item()
    nFrames = 1200
    frame = np.zeros((512, 4096), dtype=np.uint16)
    sum_frame = np.zeros((512, 4096), dtype=np.uint16)
    for i in range(nFrames):
        rc = _counted_decoded_frames[i]
        print(np.shape(rc))
        print (np.max(rc[:,0]), np.max(rc[:,1]))
        frame.fill(0)
        frame[rc[:,0], rc[:,1]] = 1
        sum_frame = np.add(sum_frame,frame)
        print(i, np.sum(frame))
        # tiff.imsave(os.path.join('/home/abhik/code/ReCoDe/scratch','frame_' + str(i) + '.tiff'), frame)
    tiff.imsave(os.path.join('/home/abhik/code/ReCoDe/scratch','sum_image_2.tiff'), sum_frame)