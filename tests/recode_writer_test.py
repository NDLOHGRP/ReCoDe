from pyrecode.recode_writer import ReCoDeWriter, print_run_metrics, _bit_pack
import argparse
import numpy as np

if __name__== "__main__":

    a = np.zeros((2), dtype=np.uint16)
    a[0] = 78
    a[1] = 145
    pa=_bit_pack(a, 12)
    print(pa)

    parser = argparse.ArgumentParser(description='ReCoDe Queue Manager')
    parser.add_argument('--image_filename', dest='image_filename', action='store', default='', help='path of folder containing data (typically inside RAM disk for on-the-fly)')
    parser.add_argument('--calibration_file', dest='calibration_file', action='store', default='', help='path to calibration file')
    parser.add_argument('--out_dir', dest='out_dir', action='store', default='', help='output directory')
    parser.add_argument('--params_file', dest='params_file', action='store', default='', help='path to params file')
    parser.add_argument('--mode', dest='mode', action='store', default='batch', help='batch or stream')
    parser.add_argument('--validation_frame_gap', dest='validation_frame_gap', action='store', type=int, default=-1, help='validation frame gap')
    parser.add_argument('--log_file', dest='log_file', action='store', default='', help='path to log file')
    parser.add_argument('--run_name', dest='run_name', action='store', default='run_1', help='run name')
    parser.add_argument('--verbosity', dest='verbosity', action='store', type=int, default=0, help='verbosity level')
    parser.add_argument('--use_c', dest='use_c', action='store_true', help='')
    parser.add_argument('--max_count', dest='max_count', action='store', type=int, default=1, help='the number of chunks to process')
    parser.add_argument('--chunk_time_in_sec', dest='chunk_time_in_sec', action='store', type=int, default=1, help='seconds of data contained in each chunk')
    

    args = parser.parse_args()

    writer = ReCoDeWriter (
        args.image_filename, 
        args.calibration_file, 
        output_directory=args.out_dir, 
        input_params=None, 
        params_filename=args.params_file, 
        mode='batch', 
        validation_frame_gap=50, 
        log_filename='recode.log', 
        run_name='run', 
        verbosity=0, use_C=False, max_count=-1, chunk_time_in_sec=0, node_id=0)

    writer.start()
    run_metrics = writer.run()
    writer.close()
    print_run_metrics(run_metrics)