import numpy as np
import os
import sys
from os.path import isfile, join
import tifffile as tiff
import scipy.ndimage as nd
import pims
from skimage.measure import regionprops

dirpath = '/scratch/loh/abhik/10-April-2019/DE16'

# load dark noise
dark_data = pims.open(os.path.join(dirpath,'16-11-02.184_Dark_Ref_Pre_30mins_400fps_4096_512_OffsetY=1024.seq'))
darks = np.max(dark_data[0:4000], axis=0)
print(darks.shape)

'''
ref_data = pims.open(os.path.join(dirpath,'16-54-49.041_1fps_Beam_Blanker_Location.seq'), dtype='uint16')
for i in range(7):
    rf_i = ref_data[i]
    tiff.imsave('Ref_Image_' + str(i) + '.tiff', (ref_data[i] > darks).astype(np.int16))
sys.exit(0)
'''
d = []
for chunk in range(5):

    validation_data = np.fromfile(os.path.join(dirpath,'validation_frames_400fps_15min_4096_512_OffsetY=1024/_run1_part00' + str(chunk) + '_validation_frames.bin'), dtype='uint16')

    frame_sz = 4096* 512
    frame_indices = [validation_data[i] for i in range(0, len(validation_data), frame_sz+1)]
    print(frame_indices)

    i = 0

    for i in range(len(frame_indices)):

        frame_index = frame_indices[i]
        frame_start = i*(frame_sz+1)+1
        frame_end = frame_start + frame_sz
        frame_i = validation_data[frame_start:frame_end]
        image = frame_i.reshape(512, 4096)

        dark_corrected = (image > darks).astype(np.int16)

        # print(image)
        # print(validation_data[0:100])
        # tiff.imsave('Frame_' + str(frame_index) +'.tiff', image)
        # tiff.imsave('Dark_Corrected_Frame_' + str(frame_index) +'.tiff', dark_corrected)

        s = nd.generate_binary_structure(2,2)
        labeled_foreground, num_features = nd.measurements.label(dark_corrected, structure=s)
        regions = regionprops(labeled_foreground, dark_corrected)
        print(num_features)

        roi_1 = {'x':1750, 'y':50, 'w':150, 'h':400}
        roi_2 = {'x':2350, 'y':50, 'w':150, 'h':400}

        count_1 = 0
        count_2 = 0
        for region in regions:
            c = region.centroid
            if (c[0] > roi_1['y'] and c[0] < roi_1['y']+roi_1['h'] and
                c[1] > roi_1['x'] and c[1] < roi_1['x']+roi_1['w']):
                count_1 += 1
                
            if (c[0] > roi_2['y'] and c[0] < roi_2['y']+roi_2['h'] and
                c[1] > roi_2['x'] and c[1] < roi_2['x']+roi_2['w']):
                count_2 += 1
                
        d_i = [chunk, i, frame_indices[i], count_1, count_2]
        d.append(d_i)
        print(d_i)
        
np.savetxt('counts_over_time.txt', d, fmt="%d", delimiter='\t')