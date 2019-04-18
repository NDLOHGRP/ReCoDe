import os
import sys
from os.path import isfile, join
import time
import zmq
import queue

# path = "R:\\"
path = "R:\DE16"

context = zmq.Context()
print ("Connecting to server...")
socket = context.socket(zmq.REQ)
socket.connect ("tcp://127.0.0.1:18534")
print ("Ready")

q = queue.Queue()
queued_files = {}
isFirst = True

count = 0
max_count = 3

while(count < max_count-1):

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
    
# print(f_list)
# print('{0:d} files in Queue.'.format(q.qsize()))

# Process the last chunk
fname = q.get()

# Wait for the last file to be written (necessary when processing is faster than acquisition)
time.sleep(10)
    
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