# WHN 10/3/09

# Extract planar velocity field.

from spivet import pivlib
from numpy import *

# >>>>> SETUP <<<<<                                                             
pdifname  = "../PIVDATA-WFV.ex2"    # Input file for PIVData object.
pdofname  = "SIDE-PLNRVEL.ex2"           # Output file name for PIVData object.
velkey    = "U-MF-SN-ZD-GS"         # Velocity key for vort, divrg, trace.

zplane = 67

# >>>>> END USER MODIFIABLE CODE <<<<< 

pd     = pivlib.loadpivdata(pdifname)
oog    = array(pd.origin)
oog[0] = oog[0] +zplane*pd.cellsz[0]
opd    = pivlib.PIVData(pd.cellsz,oog)

ncells = pd[0].eshape
for e in xrange(len(pd)):
    vel  = pd[e][velkey]
    pvel = vel[:,zplane,...].reshape(3,1,ncells[1],ncells[2])
    pvel = pivlib.cpivvar(pvel,vel.name,vel.units)
    pvel[0,...] = 0

    opd.addVars(e,pvel)

opd.setTimes(pd.getTimes())
opd.save(pdofname)
