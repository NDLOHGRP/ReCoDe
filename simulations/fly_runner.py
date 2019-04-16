import subprocess

name = "ReCoDeServerTest"
buffer_size = 1024*1024
nImageFrames = 700
min_threads = 3
max_threads = 9
threads_step = 3
rc_operation_mode = 1   # 1 = reduce and compress; 0 = reduce only
replicates = 2

darkFile = "D:/Current Downloads/Dark_Frame_12-23-00.232.bin"
destination_directory = "D:/Current Downloads/on_the_fly_output/"

# dose_rates = [0.8, 1.6, 3.2, 6.4]
dose_rates = [0.8]

for dose_rate in dose_rates:

    if (dose_rate == 0.8):
        imageFile = "D:/Current Downloads/12-06-23.598_Part.bin"
    elif (dose_rate == 1.6):
        imageFile = "D:/Current Downloads/12-06-23.598_Part.bin"
    elif (dose_rate == 3.2):
        imageFile = "D:/Current Downloads/12-06-23.598_Part.bin"
    elif (dose_rate == 6.4):
        imageFile = "D:/Current Downloads/12-06-23.598_Part.bin"

    subprocess.run(["D:/cbis/GitHub/ReCoDe/win32/x64/Release/OnTheFlyReCoDe.exe", name, str(buffer_size), str(nImageFrames), str(min_threads), str(max_threads), str(threads_step), str(dose_rate), str(rc_operation_mode), str(replicates), "D:/Current Downloads/12-06-23.598_Part.bin", "D:/Current Downloads/Dark_Frame_12-23-00.232.bin", "D:/Current Downloads/on_the_fly_output/"])