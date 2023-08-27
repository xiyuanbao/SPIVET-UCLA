# WHN 10/2/09
# XBao 04/16/2022
# FTLE script.  Creates a full progression.

from spivet import pivlib, flolib
from spivet.flolib import floftle
from numpy import *

import os
# set tmp dir to a large space
#https://stackoverflow.com/questions/11697214/how-to-set-the-tmpdir-environment-variable-to-another-directory
tdir = "/data/xbao/tmp_python"
os.environ["TMPDIR"]=tdir
try:
    os.mkdir(tdir)
except:
    pass

# >>>>> SETUP <<<<<
#pdifname  = "../PIVDATA-WFV.ex2" 
pdifname  = "PLNRVEL.ex2"    # Input file for PIVData object.

velkey = "U-MF-SN-ZD-GS"            # Velocity key for vort, divrg, trace.
tssdiv = 100               # Time steps per Epoch for tracing.

# Passive tracer setup for mpsvtrace().
epslc    = slice(3,7)                  # Epoch slice for processing.
direc    = -1
eplim    = [0,6] 
ntrpc    = [1,4,4]                      # Number of tracers per cell per call.
ncalls   = 1                            # Number of calls. >1 not implemented in FTLE core lib.
irspath  = "PLANE-FTLETRACE-BWD-NOBACKUP"     # Intermidate results output path.
pdofname = "PLANE-FTLEFIELD-BWD.ex2"          # Output file for PIVData object.
trace    = True                         # Advect tracers.
interp   = ['C','C']                    # Interpolation scheme.
max_proc = 10
# >>>>> END USER MODIFIABLE CODE <<<<< 

# pd = pivlib.loadpivdata(pdifname)
times = []
params = (pdifname,velkey,tssdiv,epslc,direc,eplim,ntrpc,ncalls,irspath,pdofname,
            trace,interp)
from multiprocessing import Pool

def single_FTLE(e,params):
    pdifname,velkey,tssdiv,epslc,direc,eplim,ntrpc,ncalls,irspath,pdofname,\
            trace,interp = params
    pd = pivlib.loadpivdata(pdifname)
    ecnt = e-epslc.start
    if (direc == 1):
        pepslc = slice(e,e+1,1)
        eisz   = eplim[1] -e
    else:
        pepslc = slice(e,e-1,-1)
        eisz   = e -eplim[0]

    ftledict = floftle.ftleinit(pd,pepslc,eisz,ntrpc,ncalls,tssdiv)

    if ( trace ):
        floftle.ftletrace(pd,velkey,irspath,ftledict,interp=interp)


    tfpd = floftle.ftlecomp(irspath,ftledict)    

    tfpd.save(str(ecnt)+"_"+pdofname)
    return ecnt
#     return [ecnt,tfpd]#PIVDATA error!
    
myPool = Pool(min(epslc.stop-epslc.start,max_proc))
procs=[]
ecnt  = 0
for e in xrange(epslc.start,epslc.stop):
    p=myPool.apply_async(single_FTLE,(e,params))
    procs.append(p)
    
    
myPool.close()
myPool.join()
output = [p.get() for p in procs]

# # sort pds
# s=[]
# for i in range(len(output)):
#     s.append(output[i][0])
# sort_index = argsort(s)    #numpy
# times = []
# for i in sort_index:
#     if output[i][0]==0:
#         fpd = output[i][1]
#         times.append(fpd[0].time)
#     else:
#         tfpd = output[i][1]
#         fpd.addVars(output[i][0],tfpd[0].getVars())
#         times.append(tfpd[0].time)       
# fpd.setTimes(times)
# fpd.save(pdofname)

times = []
for ecnt in range(len(output)):
    tfpd = pivlib.loadpivdata(str(ecnt)+"_"+pdofname)
    if ecnt == 0:
        fpd = tfpd
        times.append(fpd[0].time)
    else:
        fpd.addVars(ecnt,tfpd[0].getVars())
        times.append(tfpd[0].time)
    os.remove(str(ecnt)+"_"+pdofname)
fpd.setTimes(times)
fpd.save(pdofname)

os.system('ls -l '+tdir)
#     if (direc == 1):
#         pepslc = slice(e,e+1,1)
#         eisz   = eplim[1] -e
#     else:
#         pepslc = slice(e,e-1,-1)
#         eisz   = e -eplim[0]

#     ftledict = floftle.ftleinit(pd,pepslc,eisz,ntrpc,ncalls,tssdiv)

#     if ( trace ):
#         floftle.ftletrace(pd,velkey,irspath,ftledict,interp=interp
# 				,pdifname=pdifname
# 			)

#     if ( ecnt == 0 ):
#         fpd = floftle.ftlecomp(irspath,ftledict)    
#         times.append(fpd[0].time)
#     else:
#         tfpd = floftle.ftlecomp(irspath,ftledict)    
#         fpd.addVars(ecnt,tfpd[0].getVars())
#         times.append(tfpd[0].time)        

#     ecnt = ecnt +1

# fpd.setTimes(times)
# fpd.save(pdofname)
