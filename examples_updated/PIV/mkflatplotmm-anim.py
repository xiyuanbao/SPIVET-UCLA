# WHN 4/27/09.
#
# Constructs two lateral composition plots from forward tracing.  
# Unlike mkflatplot-anim, this script annotates the y-axis in mm
# starting at the bottom of the tank.
#

from spivet import pivlib, flolib
from numpy import *
import pylab, os

# >>>>> SETUP <<<<<

cdatap = "FLAT-CDATA"          # File containing composition parameters.
tdpath = "FTRACE"             # Input data path.
odpath = "FTRACEMM-ANIM"        # Output data path.
viewz  = True                  # Set to false to look along x-axis.
ptsz   = 0.1                   # Plot point size [pts].
dpi    = 600                   # Plot DPI.

# >>>>> END USER MODIFIABLE CODE <<<<<

def parsefn(fn):
    fn  = fn[0:-4]
    fnc = fn.rsplit("-",2)

    fnp = {"CALL":int( fnc[1].strip("CALL") ),
           "TS":int( fnc[2].strip("TS") )}

    return fnp

# Initialization.
if ( not os.path.exists(odpath) ):
    os.mkdir(odpath)

cdata  = pivlib.pklload(cdatap)
cntrc  = cdata['cntrc']    # In node-centered coordinates.
blkext = cdata['blkext']
ncells = cdata['ncells']
cellsz = cdata['cellsz']

# Setup the rbndx array and other plot parameters.
if ( viewz ):
    rbndx   = [ [cntrc[0]-0.5,cntrc[0]+0.5], 
                [-0.5,iinfo(int).max], 
                [-0.5,iinfo(int).max] ]  # Needs to be cell-centered.
    tcsx    = 2   # Tracer coordinate scatter x-component.
    tcsy    = 1   # Tracer coordinate scatter y-component.
    sxcntr  = cntrc[1]*cellsz[2] 
    xlim    = cellsz[2]*array([blkext[1,0],blkext[1,1]]) -sxcntr
    aspect  = 1.
else:
    rbndx   = [ [-0.5,iinfo(int).max], 
                [-0.5,iinfo(int).max], 
                [cntrc[1]-0.5,cntrc[1]+0.5] ]  # Needs to be cell-centered.
    tcsx    = 0   # Tracer coordinate scatter x-component.
    tcsy    = 1   # Tracer coordinate scatter y-component.
    sxcntr  = -cntrc[0]*cellsz[0]
    xlim    = abs(cellsz[0])*array([blkext[0,0],blkext[0,1]]) -sxcntr
    aspect  = 1.

# Get list of tracer files.
dl   = os.listdir(tdpath)
tsfp = []
ifp  = []
cnt  = 0
for f in xrange(len(dl)):
    if ( dl[cnt][-4::] != '.ex2' ):
        dl.pop(cnt)
    else:
        dl[cnt] = "%s/%s" % (tdpath,dl[cnt])
        cnt = cnt +1

# Create y-ticks.
ymx    = (ncells[1] -1. +0.5)*cellsz[1]

# Create the plots.
ts = 0
while ( len(dl) > 0 ):
    print "Processing TS%i ..." % ts

    fl  = []
    cnt = 0 
    for f in xrange(len(dl)):
        fn  = dl[cnt]
        fnc = parsefn(fn)

        if ( fnc['TS'] == ts ):
            fl.append( dl.pop(cnt) )
        else:
            cnt = cnt +1
    print fl
    [tcrd,trid] = flolib.mptassmblr(fl,rbndx)
    
    msk = trid[:,0] < 32767

    fig = pylab.figure()
    pylab.scatter(abs(cellsz[tcsx])*tcrd[msk,tcsx]-sxcntr,
                  ymx -cellsz[1]*tcrd[msk,tcsy],
                  s=ptsz,c=trid[msk,0],
                  edgecolors='none',cmap=pivlib.cm.wcontrast,vmin=0,vmax=32767)
    ax = pylab.gca()
    ax.set_aspect(aspect)
    ax.set_xlim(xlim)
    ax.set_ylim([0,ymx+cellsz[1]/2.])
    pylab.ylabel("Height [mm]")
    pylab.xlabel("Width [mm]")

    pylab.savefig("%s/FLATCOMP-TS%04i.png" % (odpath,ts),dpi=dpi)

    pylab.close(fig)

    ts = ts +1
