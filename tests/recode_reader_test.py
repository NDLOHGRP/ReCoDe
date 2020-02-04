from pyrecode.recode_reader import ReCoDeReader, merge_parts
from pyrecode.utils.converters import L1_to_L4
import numpy as np
import os
import sys
from scipy.sparse import csr_matrix

'''
file_name = 'D:/cbis/GitHub/ReCoDe/scratch/400fps_dose_43.rc1'
reader = ReCoDeReader(file_name, is_intermediate=False)
reader.open()
frame_data = reader._get_frame(0)
for i in range(2):
    frame_data = reader._get_next_frame()    
    print(frame_data)
reader.close()
'''

'''
intermediate_file_name = 'D:/cbis/GitHub/ReCoDe/scratch/400fps_dose_43.rc1_part000'
reader = ReCoDeReader(intermediate_file_name, is_intermediate=True)
reader.open()
# frame_data = reader._get_frame(5)
for i in range(3):
    frame_data = reader._get_next_frame()
    print(frame_data)
reader.close()
'''

'''
intermediate_file_name = 'D:/cbis/GitHub/ReCoDe/scratch/400fps_dose_43.rc1_part000'
reader = ReCoDeReader(intermediate_file_name, is_intermediate=True)
reader.open()
for i in range(0,reader._get_shape()[0]):
    frame_data = reader._get_next_frame()
    print(frame_data)
    if frame_data is None:
        break
reader.close()
'''

'''
_data_folder = '/scratch/loh/abhik/23-Jan-2020/nano_cluster'
_tag = 'gold_nano_1k_1k'
frames = {}
fid_summary = {}
num_parts = 12
for index in range(num_parts):
    intermediate_file_name = os.path.join(_data_folder, _tag + '.rc1_part' + '{0:03d}'.format(index))
    reader = ReCoDeReader(intermediate_file_name, is_intermediate=True)
    reader.open()
    print('Processing file ' + str(index) + ' of ' + str(num_parts))
    frames_i = {}
    for i in range(reader._get_shape()[0]):
        d = reader._get_next_frame()
        fid = list(d.keys())[0]
        print(fid)
        if d is not None:
            frames_i[fid] = d[fid]
        else:
            break
    reader.close()
    fids = list(frames_i.keys())
    fids = list(map(int, fids))
    fid_summary[index] = {'Start': np.min(fids), 'End': np.max(fids), 'Count': len(fids)}
    print(fids)
print('Extracted ' + str(len(frames)) + ' frames: ')
print(fid_summary)
sys.exit(0)
'''

_data_folder = '/scratch/loh/abhik/23-Jan-2020/gold_nano_1k_1k_dr=0.06'
_tag = 'gold_nano_1k_1k'
fname = os.path.join(_data_folder, _tag + '.rc1.npy')
x = np.load(fname)
d = x.item()
counted_d = L1_to_L4(d)
np.save(os.path.join(_data_folder, _tag + '.rc4.npy'), counted_d)
sys.exit(0)


_data_folder = '/scratch/loh/abhik/23-Jan-2020/gold_nano_1k_1k_dr=0.06'
_tag = 'gold_nano_1k_1k'
intermediate_file_name = os.path.join(_data_folder, _tag + '.rc1_part000')
reader = ReCoDeReader(intermediate_file_name, is_intermediate=True)
reader.open()
frames = {}
nz = reader._get_shape()[0]
nz = 5
for i in range(nz):
    d = reader._get_next_frame()
    print(d)
    if d is not None:
        frames.update(d)
    else:
        break
reader.close()


'''
Merge part files
'''
merged_frames = merge_parts(_data_folder, _tag + '.rc1', 12)
np.save(os.path.join(_data_folder, _tag + '.rc1'), merged_frames)


'''
Sum dictionary of 2-d frames in COO format to a single dense frame
'''
_out_folder = '/home/abhik/code/ReCoDe/scratch'
sz = frames[0].shape
view = np.zeros(sz)
for frame_id in frames:
    a = frames[frame_id].toarray()
    view += frames[frame_id].toarray()
np.save(os.path.join(_out_folder, 'view_' + _tag + '.npy'), view)

'''
Convert dictionary of 2-d frames in COO format to a 3-d matrix in CSR format
tot_nnz = 0
for frame in frames:
    tot_nnz += frame.nnz
data = csr_matrix((3, 4), dtype=np.int8)
'''
