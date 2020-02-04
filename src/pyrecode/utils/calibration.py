import numpy as np
from datetime import datetime
import pims
import os
import argparse
import scipy.ndimage as nd
from skimage.measure import regionprops
from datetime import datetime

def _save_file(arr, filename):
    arr.astype('uint16').tofile(filename)

def _count_events(frame, t):
    _binary_frame = frame > t
    s = nd.generate_binary_structure(2,2)
    labeled_foreground, num_features = nd.measurements.label(_binary_frame, structure=s)
    return num_features, np.sum(_binary_frame)

def make_calibration_frames (filepath, nFrames, n_stats_frames, n_sigmas, savepath='', filename_prefix=''):
    
    if not filename_prefix.endswith('_'):
        filename_prefix += '_'
    
    start=datetime.now()
    fp = pims.open(filepath)
    d = np.zeros((nFrames, fp[0].shape[0], fp[0].shape[1]), dtype=np.uint16)
    for frame_index in range(nFrames):
        d[frame_index] = fp[frame_index]

    _m = np.median(d, axis=0)
    _std = np.std(d, axis=0)
    print("Calibration time:", datetime.now()-start)

    n_pixels_in_frame = d[0].shape[0]*d[0].shape[1]    
    for i in range(n_sigmas):
        t = np.floor(_m + _std*i).astype(np.uint16)
        _save_file(t, os.path.join(savepath, filename_prefix + "dark_ref_" + str(i) + ".bin"))
        n_events = 0
        p_foreground_pixels = 0
        for f in range(nFrames-n_stats_frames, nFrames):
            n_e, n_fp = _count_events(d[f].astype(np.uint16), t)
            n_events += n_e
            p_foreground_pixels += (n_fp/n_pixels_in_frame)
        avg_n_events = n_events/n_stats_frames
        avg_p_foreground_pixels = p_foreground_pixels/n_stats_frames
        print("Avg. prop. foreground pixels for sigma=" + str(i) + " is: " + str(avg_p_foreground_pixels))
        print("Avg. electron count for sigma=" + str(i) + " is: " + str(avg_n_events))
        print("Avg. dose rate for sigma=" + str(i) + " is: " + str(avg_n_events/n_pixels_in_frame))
        print("")

if __name__== "__main__":

    parser = argparse.ArgumentParser(description='ReCoDe Queue Manager')
    parser.add_argument('--flatfield_filepath', dest='filepath', action='store', default='', help='path to the flat-field illuminated file to be used for calibration')
    parser.add_argument('--n_frames', dest='n_frames', action='store', type=int, default=100, help='number of frames to use for calibration')
    parser.add_argument('--n_stats_frames', dest='n_stats_frames', action='store', type=int, default=10, help='number of frames on which to estimate dose rate')
    parser.add_argument('--n_sigmas', dest='n_sigmas', action='store', type=int, default=4, help='the number of sigmas to try')
    parser.add_argument('--savepath', dest='savepath', action='store', default='', help='path to folder where calibration frames are to be stored')
    parser.add_argument('--save_prefix', dest='filename_prefix', action='store', default='', help='prefix for calibration filename')
    args = parser.parse_args()
    make_calibration_frames (args.filepath, args.n_frames, args.n_stats_frames, args.n_sigmas, args.savepath, args.filename_prefix)
    