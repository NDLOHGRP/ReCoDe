from pyrecode.recode_reader import ReCoDeReader, merge_parts
from pyrecode.utils.converters import L1_to_L4, recalibrate_L1, read_dark_ref
import numpy as np
import os
import sys


_data_folder = '/scratch/loh/abhik/28-Feb-2010'
# _datasets = ['500Kmag_300fps_NoCorrection']
_datasets = [
    '300Kmag_300fps_NoCorrection','300Kmag_300fps_1',
    '200Kmag_300fps_NoCorrection','120Kmag_300fps_NoCorrection','120Kmag_300fps_LowerDose',
    '120Kmag_300fps_HigherDose','200Kmag_300fps','300Kmag_300fps_6_DeepInside',
    '300Kmag_300fps_5_HigherDose','300Kmag_300fps_4_Edge','300Kmag_300fps_3','300Kmag_300fps_2']

_ref_2 = np.transpose(read_dark_ref(os.path.join(_data_folder, 'DE16', '_dark_ref_2.bin'), (1024,1024), np.uint16))
_ref_4A = np.transpose(read_dark_ref(os.path.join(_data_folder, 'temp', '__dark_ref_4A.bin'), (1024,1024), np.uint16))
# _ref_2 = read_dark_ref(os.path.join(_data_folder, 'DE16', '_dark_ref_2.bin'), (1024,1024), np.uint16)
# _ref_4A = read_dark_ref(os.path.join(_data_folder, 'temp', '__dark_ref_4A.bin'), (1024,1024), np.uint16)

for i in range(len(_datasets)):
    
    print("Processing File: " + _datasets[i])
    print("=========================================")
    
    fname = os.path.join(_data_folder, _datasets[i], _datasets[i] + '.rc1.npy')
    x = np.load(fname)
    merged_frames = x.item()

    new_L1 = recalibrate_L1(merged_frames, np.uint16, original_calibration_frame = _ref_2, new_calibration_frame=_ref_4A, epsilon=5.0)
    np.save(os.path.join(_data_folder, _datasets[i], _datasets[i] + '_4A_5.rc1.npy'), new_L1)
    
    counted_d = L1_to_L4(new_L1, area_threshold=1)
    np.save(os.path.join(_data_folder, _datasets[i], _datasets[i] + '_4A_5.rc4.npy'), counted_d)

sys.exit(0)

_data_folder = '/scratch/loh/abhik/28-Feb-2010'
_datasets = [
    '300Kmag_300fps_NoCorrection_Lattice','500Kmag_300fps_NoCorrection','300Kmag_300fps_NoCorrection','300Kmag_300fps_1',
    '200Kmag_300fps_NoCorrection','120Kmag_300fps_NoCorrection','120Kmag_300fps_LowerDose',
    '120Kmag_300fps_HigherDose','200Kmag_300fps','300Kmag_300fps_6_DeepInside',
    '300Kmag_300fps_5_HigherDose','300Kmag_300fps_4_Edge','300Kmag_300fps_3','300Kmag_300fps_2']

for i in range(len(_datasets)):
    print("Processing File: " + _datasets[i])
    print("=========================================")
    for filename in os.listdir(os.path.join(_data_folder, _datasets[i])):
        if filename.endswith(".rc1_part000"):
            _tag = os.path.splitext(filename)[0]
            print(_tag)
    merged_frames = merge_parts(os.path.join(_data_folder, _datasets[i]), _tag + '.rc1', 20)
    fname = os.path.join(_data_folder, _datasets[i], _datasets[i] + '.rc1.npy')
    np.save(fname, merged_frames)

    for k in range(len(merged_frames)):
        if k not in merged_frames:
            print('Frame ' + str(k) + ' missing in ' + fname)

    counted_d = L1_to_L4(merged_frames)
    np.save(os.path.join(_data_folder, _datasets[i], _datasets[i] + '.rc4.npy'), counted_d)

'''
_data_folder = '/scratch/loh/abhik/12_Feb_2020'
_datasets = ['L1_stacking_ref1_t1_60k','L1_stacking_ref1_t2_60k','L1_stacking_ref1_t3_60k', 
             'L1_stacking_ref2_t1_60k','L1_stacking_ref2_t2_60k','L1_stacking_ref2_t3_60k', 
             'L1_stacking_ref3_t1_60k','L1_stacking_ref3_t2_60k','L1_stacking_ref3_t3_60k']
for i in range(len(_datasets)):
    for filename in os.listdir(os.path.join(_data_folder, _datasets[i])):
        if filename.endswith(".rc1_part000"):
            _tag = os.path.splitext(filename)[0]
            print(_tag)
    merged_frames = merge_parts(os.path.join(_data_folder, _datasets[i]), _tag + '.rc1', 12)
    fname = os.path.join(_data_folder, _datasets[i], _datasets[i] + '.rc1.npy')
    np.save(fname, merged_frames)
'''

'''
_data_folder = '/scratch/loh/abhik/23-Jan-2020'
_datasets = ['gold_nano_1k_1k_dr=0.03','gold_nano_1k_1k_dr=0.06','nano_cluster_low_dose_rate_threshold_2', 'nano_cluster_low_dose_rate_threshold_1', 'carbon_edge', 'nano_cluster']
_tags = ['gold_nano_1k_1k', 'gold_nano_1k_1k', 'gold_nano_1k_1k_2', 'gold_nano_1k_1k', 'gold_nano_1k_1k', 'gold_nano_1k_1k']
for i in range(len(_datasets)):
    merged_frames = merge_parts(os.path.join(_data_folder, _datasets[i]), _tags[i] + '.rc1', 12)
    fname = os.path.join(_data_folder, _datasets[i], _tags[i] + '.rc1.npy')
    np.save(fname, merged_frames)
    x = np.load(fname)
    d = x.item()
    counted_d = L1_to_L4(d)
    np.save(os.path.join(_data_folder, _datasets[i], _tags[i] + '.rc4.npy'), counted_d)
'''

'''
_data_folder = '/scratch/loh/abhik/23-Jan-2020'
_dataset = 'gold_nano_1k_1k_dr=0.06'
_tag = 'gold_nano_1k_1k'
fname = os.path.join(_data_folder, _dataset, _tag + '.rc1.npy')
x = np.load(fname)
d = x.item()
counted_d = L1_to_L4(d, n_Frames=600)
sys.exit(0)
'''