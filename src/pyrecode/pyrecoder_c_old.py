import sys
import pyrecoder

print( dir(pyrecoder) )
print(pyrecoder.__doc__)

recoder = pyrecoder.Recoder()
recoder.open_file("Dark_Frame_12-23-00.232.bin")

mv = memoryview(bytes(1024))
state = recoder.get_recode_header(mv)
print('state=',state)
for i in range(10):
    print(mv[i])
# mv.release()

mv = memoryview(bytes(10))
state = recoder.get_frame(1, mv)
print('state=',state)
for i in range(10):
    print(mv[i])


sys.exit(0)

    
mv = recoder.get_frame(1)
for d in mv:
    print(d)
mv = recoder.get_next_frame()
mv = recoder.get_frames(1,10)
recoder.close_file()




f = pyrecoder.open()
print(f)

'''
mv = pyrecoder.read_frame()
for d in mv:
    print(d)
mv.release()
'''

with pyrecoder.read_frame() as mv:
    for d in mv:
        print(d)
        
        
        
# [213, 1, 158, 1, 197, 1, 156, 1, 189, 1]