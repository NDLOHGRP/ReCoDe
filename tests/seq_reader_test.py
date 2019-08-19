import numpy as np
import pims

# image = pims.open('D://cbis//images//20161207//14-37-50.811.seq')
# image = pims.open('D:/cbis/GitHub/ReCoDe/scratch/temp/007.seq')
# image = pims.open('R://003.seq')
image = pims.open('D:/cbis/images/SequenceBlocks/17-21-33.138_Dark_Ref_3.seq')

for index,frame in enumerate(image):
    print('Frame', index, ': 110,99=', frame[110,99], '220,1111=', frame[220,1111], '330,2222=', frame[330,2222], '440,3333=', frame[440,3333], '500,4000=', frame[500,4000])

# np.savetxt("D:/cbis/GitHub/ReCoDe/scratch/frame_0_PIMS.txt", image[300], fmt="%d", delimiter=" ")
# print(np.shape(image))