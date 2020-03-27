import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pyrecode.recode_reader import ReCoDeReader
import pims
import imageio
from mpl_toolkits.axes_grid1 import make_axes_locatable
from multiprocessing import Process, Queue
import sys

'''
ToDo: make_animation, make_gif
'''


'''
ToDo: Add support for getting frames from intermediate file and recode file in addition to dictionaries
      Should automatically check for type of d and read accordingly
'''
def plot_frames(
    d, src='sparse_dict', frame_rate=-1, pixel_size=-1, fractionation=1, 
    start=0, reps=1, title=None, vmax_ptile=90, 
    fft=False, fft_type='fft_of_sum_frames', nx_fft=256, ny_fft=256,
    show_plot=True, figsize=(20,10), savepath=None, dpi=600):
    
    assert src in ['seq', 'sparse_dict'], 'Unknown src'
    
    if src == 'seq':
        a = d[0]
    elif src == 'sparse_dict':
        a = d[0].toarray()
    
    df = fractionation
    nx, ny = a.shape
    Nf = reps
    
    if fft:
        assert fft_type in ['fft_of_sum_frames','sum_of_ffts'], 'Unknown param value for fft_tpye'
        ry = int(ny_fft/2)
        rx = int(nx_fft/2)
        
    assert show_plot or savepath is not None, 'Either show_plot or savepath must be selected'
    
    for i in range(start,start+(df*Nf),df):
        
        view = np.zeros((nx,ny))
        part_view = np.zeros((nx,ny))
        if fft:
            sum_ffts = np.zeros((ny_fft,nx_fft))
            
        for index in range(i,i+df):
            
            if src == 'seq':
                a = d[index]
            elif src == 'sparse_dict':
                a = d[index].toarray()
            
            view = np.add(view, a)
            part_view = np.add(part_view, a)
            
            if fft and fft_type == 'sum_of_ffts':
                if ny_fft < ny:
                    # cy = np.random.randint(ry,ny-ry)
                    cy = int(ny_fft/2)
                else:
                    cy = int(ny_fft/2)
                
                if nx_fft < nx:
                    # cx = np.random.randint(rx,nx-rx)
                    cx = int(nx_fft/2)
                else:
                    cx = int(nx_fft/2)
                    
                if index%30 == 0:
                    part_view = np.zeros((nx,ny))
                    
                sum_ffts = np.add(sum_ffts, 
                             np.abs(np.fft.fftshift(np.fft.fft2(np.fft.fftshift(
                                 part_view[cy-ry:cy+ry,cx-rx:cx+rx])))))
                
        if i == start:
            _vmax = np.percentile(view,vmax_ptile)

        if fft:
            fig, axs = plt.subplots(nrows=1, ncols=2, figsize=figsize)
            ax = axs[0]
        else:
            fig, ax = plt.subplots(figsize=figsize)
            
        im = ax.imshow(view, vmin=np.min(view), vmax=_vmax)
        
        _title_str = ''
        if title is not None:
            _title_str += title
        if frame_rate > -1:
            timestamp = 'Time = {0:.3f} to {1:.3f} seconds ({2:d} frames)'.format(1.*i/frame_rate, 1.*(i+df)/frame_rate, df)
            _title_str += '\n' + timestamp
        
        dose_rate = np.sum(view)/(nx*1.0*ny*1.0*df)
        dr = 'Dose Rate for Dataset = {0:0.4f} e/p/f'.format(dose_rate)
        fd = 'Dose in Frame = {0:0.4f} e/p/f'.format(dose_rate*fractionation)
        td = 'Total Dose Suffered at the Start of Frame = {0:.3f} e/p'.format(i*dose_rate)
        _title_str += '\n' + dr + '. ' + fd + '\n' + td
        
        ax.set_title(_title_str)

        if pixel_size > -1:
            _x_ticks =  [j for j in range(0,nx+1,200)]
            _y_ticks =  [j for j in range(0,ny+1,200)]
            _x_labels = ['{0:0.2f}'.format(pixel_size*j) for j in range(0,nx+1,200)]
            _y_labels = ['{0:0.2f}'.format(pixel_size*j) for j in range(0,ny+1,200)]

            ax.set_xticks(_x_ticks)
            ax.set_xticklabels(_x_labels)
            ax.set_xlabel('nm')

            ax.set_yticks(_y_ticks)
            ax.set_yticklabels(_y_labels)
            ax.set_ylabel('nm')
        else:
            ax.set_xlabel('pixels')
            ax.set_ylabel('pixels')
            
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        cbar = fig.colorbar(im, cax=cax)

        if fft:
            if fft_type == 'fft_of_sum_frames':
                for bootstrap in range(100):
                    if ny_fft < ny:
                        cy = np.random.randint(ry,ny-ry)
                    else:
                        cy = int(ny_fft/2)

                    if nx_fft < nx:
                        cx = np.random.randint(rx,nx-rx)
                    else:
                        cx = int(nx_fft/2)
                    img = np.abs(np.fft.fftshift(np.fft.fft2(np.fft.fftshift(view[cy-ry:cy+ry,cx-rx:cx+rx]))))
                    sum_ffts = np.add(sum_ffts,img)
            
            img_1 = np.log(sum_ffts)
            img_1[np.isinf(img_1)] = np.max(img_1[np.logical_not(np.isinf(img_1))])
            
            im1 = axs[1].imshow(img_1, cmap='gray')
            
            if fft_type == 'fft_of_sum_frames':
                axs[1].set_title('FFT of Summed Frames')
            elif fft_type == 'sum_of_ffts':
                axs[1].set_title('Sum of FFTs')
            
            resolutions = (nx_fft/np.arange(-int(nx_fft/2),int(nx_fft/2)))*(pixel_size*10)
            _x_ticks_fft =  [j for j in range(0,ny_fft+1,int(ny_fft/10))]
            _y_ticks_fft =  [j for j in range(0,nx_fft+1,int(nx_fft/10))]
            _x_labels_fft = ['{0:0.2f}'.format(resolutions[j]) for j in range(0,nx_fft+1,int(nx_fft/10))]
            _y_labels_fft = ['{0:0.2f}'.format(resolutions[j]) for j in range(0,ny_fft+1,int(ny_fft/10))]

            axs[1].set_xticks(_x_ticks_fft)
            axs[1].set_xticklabels(_x_labels_fft)
            axs[1].set_xlabel('Resolution ($\AA$)')

            axs[1].set_yticks(_y_ticks_fft)
            axs[1].set_yticklabels(_y_labels_fft)
            axs[1].set_ylabel('Resolution ($\AA$)')
            
            divider = make_axes_locatable(axs[1])
            cax = divider.append_axes("right", size="5%", pad=0.1)
            cbar = fig.colorbar(im1, cax=cax)
            
        if savepath is not None:
            fig.savefig(savepath, dpi=dpi)
            
        if show_plot:
            plt.show()

        plt.close()

def _plot_frames(data, args_list):
    for args in args_list:
        plot_frames (
            data, 
            src='sparse_dict',
            frame_rate=args['frame_rate'], 
            pixel_size=args['pixel_size'], 
            fractionation=args['fractionation'], 
            start=args['start'], 
            reps=args['reps'], 
            title=args['title'],
            fft=args['fft'],
            fft_type=args['fft_type'],
            nx_fft=args['nx_fft'], 
            ny_fft=args['ny_fft'],
            show_plot=args['show_plot'],
            savepath=args['savepath']
        )

def plot_frames_mt (data, args, nt=1):

    if nt < 2:
        return _plot_frames (data, args)

    q = Queue()
    procs = []
    indices = np.array_split(np.arange(len(args)), nt)
    print(indices)

    for i in range(nt):
        args_list = [args[j] for j in indices[i]]
        p = Process(target=_plot_frames, args=(data, args_list))
        procs.append(p)
        p.start()

    for proc in procs:
        proc.join()


if __name__== "__main__":

    _data_folder = '/scratch/loh/abhik/28-Feb-2010'
    _datasets = ['300Kmag_300fps_4_Edge','']

    '''
    _datasets = [
        '300Kmag_300fps_NoCorrection_Lattice',
        '500Kmag_300fps_NoCorrection',
        '300Kmag_300fps_NoCorrection',
        '300Kmag_300fps_1',
        '200Kmag_300fps_NoCorrection',
        '120Kmag_300fps_NoCorrection',
        '120Kmag_300fps_LowerDose',
        '120Kmag_300fps_HigherDose',
        '200Kmag_300fps',
        '300Kmag_300fps_6_DeepInside',
        '300Kmag_300fps_5_HigherDose',
        '300Kmag_300fps_3',
        '300Kmag_300fps_2']
    '''
    
    _magnifications = {}
    _magnifications['300Kmag_300fps_NoCorrection_Lattice'] = 300
    _magnifications['500Kmag_300fps_NoCorrection'] = 500
    _magnifications['300Kmag_300fps_NoCorrection'] = 300
    _magnifications['300Kmag_300fps_1'] = 300
    _magnifications['200Kmag_300fps_NoCorrection'] = 200
    _magnifications['120Kmag_300fps_NoCorrection'] = 120
    _magnifications['120Kmag_300fps_LowerDose'] = 120
    _magnifications['120Kmag_300fps_HigherDose'] = 120
    _magnifications['200Kmag_300fps'] = 200
    _magnifications['300Kmag_300fps_6_DeepInside'] = 300
    _magnifications['300Kmag_300fps_5_HigherDose'] = 300
    _magnifications['300Kmag_300fps_4_Edge'] = 300
    _magnifications['300Kmag_300fps_3'] = 300
    _magnifications['300Kmag_300fps_2'] = 300

    for dataset_index in range(len(_datasets)):
        for recode_level in [4,1]:
            for correction_type in ['','_4A_5']:

                fname = os.path.join(_data_folder, _datasets[dataset_index], _datasets[dataset_index] + correction_type + '.rc' + str(recode_level) + '.npy')
                x = np.load(fname)
                data = x.item()

                magnification = _magnifications[_datasets[dataset_index]]
                _pixel_size = 0.102/(magnification/60)

                args = []
                for _start in [900]:
                    for _fractionation in [300, 600, 900, 1800, 3600]:
                        for _fft_type in ['sum_of_ffts', 'fft_of_sum_frames']:

                            args_i = {
                                'fname': fname, 
                                'frame_rate':300, 
                                'pixel_size':_pixel_size, 
                                'fractionation':_fractionation, 
                                'start':_start, 
                                'reps':1, 
                                'title':_datasets[dataset_index],
                                'fft':True,
                                'fft_type':_fft_type,
                                'nx_fft':256, 
                                'ny_fft':256,
                                'show_plot':False,
                                'savepath':os.path.join('/scratch/loh/abhik/28-Feb-2010/views_3', _datasets[dataset_index] + 
                                                                                                '_start=' + str(_start) + 
                                                                                                '_frac=' + str(_fractionation) + 
                                                                                                '_' + _fft_type +
                                                                                                '_rc=' + str(recode_level) + 
                                                                                                '_' + correction_type +
                                                                                                '.png')
                            }
                            args.append(args_i)

                plot_frames_mt (data, args, nt=20)
