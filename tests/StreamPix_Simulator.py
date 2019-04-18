import numpy as np
import pims
import os
import sys
from os.path import isfile, join
import time
import zmq
import queue

def load_file(fpath):
    with open(fpath, "rb") as binary_file:
        data = bytearray(binary_file.read())
    return data
    
def write_file(fpath, data):
    with open(fpath, "wb") as write_file:
        write_file.write(data)


ITERATIONS = 4
THROUGHPUT_IN_SECS_PER_CHUNK = 5

source_directory = "D:\\cbis\\GitHub\\ReCoDe\\scratch\\temp"
target_directory = "R:\\"

# copy all files from source directory to temp directory
f_list = [f for f in os.listdir(source_directory) if isfile(join(source_directory, f)) and f.endswith(".seq")]
seq_data = []
for f in f_list:
    seq_data.append(load_file(join(source_directory, f)))
        
print (len(seq_data))

# start serving chunks at [THROUGHPUT_IN_SECS_PER_CHUNK]
index = 0
for m in range(ITERATIONS):
    for data in seq_data:
        start_time = time.time()
        fname = join(target_directory, str(index) + '.seq')
        write_file(fname, data)
        elapsed_time = time.time() - start_time
        print("Iteration " + str(m) + ": File '" + fname + "' written in " + str(elapsed_time) + " seconds.")
        time.sleep(THROUGHPUT_IN_SECS_PER_CHUNK)
        index += 1