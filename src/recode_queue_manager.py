import os
import sys
from os.path import isfile, join
import time
import zmq
import queue
import argparse
import subprocess
from pathlib import Path

'''
#=========Params==========
# path = "R:\\"
path = "R:\DE16"
max_count = 12
chunk_time_in_sec = 5
#=========Params==========
'''

'''
D:\StreamPix\Abhik\ReCoDe\win32\x64\Release\OnTheFlyReCoDe 
-d "D:\StreamPix\Abhik\23-Aug-2019\DE16\16-41-57.509_Dark_Ref2_700fps.seq" 
-o "D:\StreamPix\Abhik\ReCoDe\scratch" 
-p "D:\StreamPix\Abhik\ReCoDe\config\recode_params.txt" 
-i abc 
-l recode.log 
-v 1 
-vf 2000 
-n 720fps_dose2_nt12
'''

parser = argparse.ArgumentParser(description='ReCoDe Queue Manager')
parser.add_argument('--path', dest='path', action='store', default='', help='path of folder containing data (typically inside RAM disk for on-the-fly)')
parser.add_argument('--max_count', dest='max_count', action='store', type=int, default=1, help='the number of chunks to process')
parser.add_argument('--chunk_time_in_sec', dest='chunk_time_in_sec', action='store', type=int, default=1, help='seconds of data contained in each chunk')
parser.add_argument('--calibration_file', dest='calibration_file', action='store', default='', help='path to calibration file')
parser.add_argument('--out_dir', dest='out_dir', action='store', default='', help='output directory')
parser.add_argument('--params_file', dest='params_file', action='store', default='', help='path to params file')
parser.add_argument('--out_file_name', dest='out_file_name', action='store', default='_', help='name to be prepended to part files')
parser.add_argument('--log_file', dest='log_file', action='store', default='', help='path to log file')
parser.add_argument('--verbosity', dest='verbosity', action='store', type=int, default=0, help='verbosity level')
parser.add_argument('--validation_frame_gap', dest='validation_frame_gap', action='store', type=int, default=-1, help='validation frame gap')
parser.add_argument('--run_name', dest='run_name', action='store', default='run_1', help='run name')

args = parser.parse_args()
print(args.path)
print(args.max_count)
print(args.chunk_time_in_sec)

path = args.path
max_count = args.max_count
chunk_time_in_sec = args.chunk_time_in_sec

_recode_server_exe = Path('D:/cbis/GitHub/ReCoDe/win32/x64/Release/OnTheFlyReCoDe')
# _recode_server_exe = Path('D:/StreamPix/Abhik/ReCoDe/win32/x64/Release/OnTheFlyReCoDe')
command = _recode_server_exe + ' -d ' + Path(args.calibration_file) + ' -o ' + Path(args.out_dir) + ' -p ' + Path(args.params_file) + ' -i ' + Path(args.out_file_name) + ' -l ' + Path(args.log_file) + ' -v ' + str(args.verbosity) + ' -vf ' + str(args.validation_frame_gap) + ' -n ' + args.run_name
print(command)
# os.system("start /wait cmd /c " + command)
subprocess.Popen(command, creationflags=subprocess.CREATE_NEW_CONSOLE)

assert os.path.isdir(path), 'Directory ' + path + ' not found.'

for f in os.listdir(path):
    os.remove(os.path.join(path,f))

with open(os.path.join('..','temp','interrupt_state'), 'w') as f:
    f.write('0\n')

context = zmq.Context()
print ("Connecting to server...")
socket = context.socket(zmq.REQ)
socket.connect ("tcp://127.0.0.1:18534")
print ("Ready")

q = queue.Queue()
queued_files = {}
isFirst = True

count = 0

has_interrupt = False
while(not has_interrupt and count < max_count-1):

    f_list = [f for f in os.listdir(path) if isfile(join(path, f)) and f.endswith(".seq")]
    
    has_new_file = False
    for f in f_list:
        if not f in queued_files:
            q.put(f)
            queued_files[f] = 1
            has_new_file = True
            
    if (has_new_file and q.qsize() > 1):
        
        fname = q.get()
        
        if isFirst:
            isFirst = False
            continue
        
        # Rename
        new_fname = join(path, "Next_Stream.seq")
        os.rename(join(path, fname), new_fname)
        
        # Read and Process
        print("Monitor Thread: Processing " + fname)
        start_time = time.time()
        socket.send_string(join(path, "Next_Stream.seq"))
        message = socket.recv()
        print(message)
        elapsed_time = time.time() - start_time
        print("Monitor Thread: Processed in " + str(elapsed_time) + " seconds.")
        count += 1
        
        # Delete
        print("Monitor Thread: Clearing " + path)
        start_time = time.time()
        os.remove(join(path, "Next_Stream.seq"))
        elapsed_time = time.time() - start_time
        print("Monitor Thread: Cleared in " + str(elapsed_time) + " seconds.\n")
        
        # Check for interrupts
        with open(os.path.join('..','temp','interrupt_state'), 'r') as f:
            content = f.readline()
            content = content.strip()
            has_interrupt = int(content) == 1
            print("Interrupted")
    
# print(f_list)
# print('{0:d} files in Queue.'.format(q.qsize()))

# Process the last chunk
fname = q.get()

# Wait for the last file to be written (necessary when processing is faster than acquisition)
time.sleep(chunk_time_in_sec+1)
    
# Rename
new_fname = join(path, "Next_Stream.seq")
os.rename(join(path, fname), new_fname)
        
# Read and Process
print("Monitor Thread: Processing " + fname)
start_time = time.time()
socket.send_string(join(path, "Next_Stream.seq"))
message = socket.recv()
print(message)
elapsed_time = time.time() - start_time
print("Monitor Thread: Processed in " + str(elapsed_time) + " seconds.")

print("Monitor Thread: shutting down")
socket.send_string("shutdown")
message = socket.recv()
print("ReCoDe Server: " + message.decode("utf-8"))
sys.exit(0)