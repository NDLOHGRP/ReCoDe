import os
from os.path import isfile, join
import time
import zmq

f1_path = "R:/buff1"
f2_path = "R:/buff2"

curr_target = 2
clear_target = False

context = zmq.Context()
print ("Connecting to server...")
socket = context.socket(zmq.REQ)
socket.connect ("tcp://127.0.0.1:18534")


while(True):

    if (curr_target == 1):
        path = f1_path
    elif (curr_target == 2):
        path = f2_path
        
    f_list = [f for f in os.listdir(path) if isfile(join(path, f)) and f.endswith(".seq")]
    
    if (clear_target):
        # Rename
        os.rename(join(path, f_list[0]), join(path, "Next_Stream.seq"))
        # Read and Process
        print("Processing " + path)
        socket.send_string(str(curr_target))
        message = socket.recv()
        print(message)
        # Delete
        print("Clearing " + path)
        start_time = time.time()
        os.remove(join(path, "Next_Stream.seq"))
        elapsed_time = time.time() - start_time
        print("Cleared in " + str(elapsed_time) + " seconds.")
        
        clear_target = False
        f_list = [f for f in os.listdir(path) if isfile(join(path, f))]
    
    
    if (len(f_list) > 0):
        # flip target and request to clear
        print("Data arrived in " + str(curr_target) + ", Flipping.")
        if (curr_target == 1):
            curr_target = 2
        else:
            curr_target = 1
            
        clear_target = True