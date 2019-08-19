import numpy as np
import scipy.ndimage as nd
from skimage.measure import regionprops
from scipy.sparse import coo_matrix

def reduce_L1 (data, dark, threshold, max_frames=10):
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

    _thresh_frame = dark + threshold

    feature_dict = {}
    s = nd.generate_binary_structure(2,2)
    reduced_frames = np.zeros((n_Frames, data.shape[1], data.shape[2]))
    for frame_index in range(max_frames):
        dark_subtracted = data[frame_index] - dark
        _binary_ind = dark_subtracted > 0
        (y1,x1) = np.nonzero(_binary_ind)
        pixvals = dark_subtracted[_binary_ind]
        reduced_frames[frame_index][_binary_ind] = pixvals
        sm = coo_matrix((pixvals, (y1,x1)), shape=(512, 4096))
        print("Num Features:", len(y1))
        feature_dict[frame_index] = sm
    return reduced_frames, feature_dict

def reduce_L2 (data, dark, threshold, max_frames=10):
    return True

def reduce_L3 (data, dark, threshold, max_frames=10):
    return True

# def reduce_L4 (data, dark, threshold, output_directory, id, offset=0, max_frames=10):
def reduce_L4 (data, dark, threshold, max_frames=10):
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