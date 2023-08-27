import openpiv
from openpiv import tools, process, validation, filters, scaling, pyprocess
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import display

base="/home/xbao/LAB_test/PLUME-3_51V-B25_2-11032008/ANALYSIS/_DEBUG/DEBUG/"
frame_a=np.loadtxt(base+"f1_1602118472.46.txt")
frame_b=np.loadtxt(base+"f2_1602118472.55.txt")

#Need to validate the result(get rid of spurious vectors) multiple times
import time

import openpiv.gpu_process
reload(openpiv.gpu_process)

# gpu code parametes
min_window_size = 32
overlap_ratio = 0.5
coarse_factor = 1
nb_iter_max = 2
dt = 1 # sec

a=time.time()
# First time is slow as the GPU modules need to compile. Once they are compiled, they stay compiled.
#Every time you run this after the first time it will be fast.
x, y, u, v, mask = openpiv.gpu_process.WiDIM( frame_a.astype(np.int32), frame_b.astype(np.int32), np.ones(frame_a.shape, dtype=np.int32),
                                                     min_window_size, 
                                                     overlap_ratio,
                                                     coarse_factor,
                                                     dt,
                                                     nb_iter_max = nb_iter_max,
                                                     trust_1st_iter = 0,
                                                     validation_iter = 3)
a=time.time()-a
print "GPU: time is",a

f1=frame_a
plt.figure(figsize=(f1.shape[1]/50,f1.shape[0]/50))
plt.quiver(x[0],y[::-1,0],u,v,color='r')
# plt.gca().invert_yaxis()
plt.imshow(frame_a,cmap="gray")
plt.show()
