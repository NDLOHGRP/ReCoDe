from pyrecode.recode_reader import ReCoDeReader, merge_parts
from pyrecode.utils.converters import L1_to_L4
import numpy as np
import os
import sys


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
_data_folder = '/scratch/loh/abhik/23-Jan-2020'
_dataset = 'gold_nano_1k_1k_dr=0.06'
_tag = 'gold_nano_1k_1k'
fname = os.path.join(_data_folder, _dataset, _tag + '.rc1.npy')
x = np.load(fname)
d = x.item()
counted_d = L1_to_L4(d, n_Frames=600)
sys.exit(0)
'''