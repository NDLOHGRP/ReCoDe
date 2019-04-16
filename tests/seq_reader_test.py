import numpy as np
import pims

# image = pims.open('D://cbis//images//20161207//14-37-50.811.seq')
image = pims.open('D:/cbis/GitHub/ReCoDe/scratch/temp/007.seq')

# np.savetxt("D:/cbis/GitHub/ReCoDe/scratch/frame_0_PIMS.txt", image[300], fmt="%d", delimiter=" ")

print(np.shape(image))