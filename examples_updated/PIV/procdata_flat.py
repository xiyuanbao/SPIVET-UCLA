"""
Filename:  procdata_flat.py
Copyright (C) 2007 William Newsome
 
This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details, published at 
http://www.gnu.org/copyleft/gpl.html

/////////////////////////////////////////////////////////////////
Description:
    Post process PIV results.  Composition is varied laterally
    and tracers are forward advected.
"""

from spivet import pivlib
from spivet import flolib
from numpy import *

# >>>>> SETUP <<<<<
pdifname  = "PIVDATA-WFV.ex2"       # Input file for PIVData object.

velkey = "U-MF-SN-ZD-GS"            # Velocity key for vort, divrg, trace.
tssdiv = 3                         # Time steps per Epoch for tracing.

# Passive tracer setup for mpsvtrace().
m_epslc   = slice(0,37)               # Epoch slice for processing.
m_ntrpc   = 20                      # Number of tracers per cell per call.
m_ncalls  = 30                      # Number of calls.
m_irspath = "FTRACE"               # Intermidate results output path.
m_adjid   = False                    # Force tracers to source ID.
m_hist    = 1                        # History level.
m_ctcsp   = True                     # Coerce tracer coords to single precision.

def tcfun(tcrd,ntrpc,pd):
    """
    Tracer coordinate adjustment callback function.
    """
    # >>>> USER SPECS <<<<    
    nlyrs = 6
    lyrh  = 9.  # [mm]
    # >>>> END USER SPECS <<<<

    ncells = pd[0].eshape
    cellsz = pd.cellsz

    tdal = nlyrs*lyrh/cellsz[1]  # Total depth of all layers [cells].

    msk = tcrd[:,1] >= ncells[1] -1. +0.5 -tdal

    return tcrd[msk,:]


def icomp(tcrd,pd):
    """
    Tracer composition callback function.
    """
    # >>>> USER SPECS <<<<
    nlyrs  = 6
    
    # 21.75 evens the stripes, while 21.5 centers over the heater.
    cntrc = [21.25,29.25]   # Center tracer coordinates [z,x].

    lyrh = 9.     # Block height [mm].  Following Farnetani:2002, should be 9.
    bw   = 23.    # Block width [mm].   Following Farnetani:2002.
    # >>>> END USER SPECS <<<<

    ncells = pd[0].eshape
    cellsz = pd.cellsz

    # Compute max permissible radius along x-axis.
    rmx   = min( ncells[2] -cntrc[1] -0.5, cntrc[1] +0.5 )  # To cell faces.
    rmx   = rmx*cellsz[2]
    rmxsq = rmx**2

    # Compute radii and cell coordinates.  Need to work in mm since cells are
    # not cubes.
    zcen = cntrc[0]*cellsz[0]
    xcen = cntrc[1]*cellsz[2]

    rztcrd = cellsz[0]*tcrd[:,0] -zcen
    rytcrd = cellsz[1]*tcrd[:,1]
    rxtcrd = cellsz[2]*tcrd[:,2] -xcen
    rsq    = rztcrd**2 +rxtcrd**2
    del rztcrd, rxtcrd

    # Setup trid.
    nbpl = int( round( (rmx -bw/2.)/bw +0.5 ) ) # Blocks/layer.    
    tnb  = nbpl*nlyrs                           # Total number of blocks.
    cmul = 32767/tnb

    brad      = bw*arange(nbpl +1,dtype=float)
    brad[1::] = brad[1::] -bw/2.
    brsq      = brad**2

    trid    = empty(tcrd.shape[0],dtype=int16)
    trid[:] = 32767

    cnt = 0
    for l in range( nlyrs ):
        print "LAYER %i:" % (l),

        ymn = cellsz[1]*(ncells[1] -1. +0.5) -(l+1)*lyrh
        ymx = ymn +lyrh
 
        lmsk = ( rytcrd >  ymn ) \
              *( rytcrd <= ymx*(1. +finfo(float).eps) )

        for b in range( nbpl ):
            cval = cnt*cmul
            print " %5i" % (cval),

            bmsk = lmsk*( rsq >= brsq[b] ) \
                       *( rsq <= brsq[b+1]*(1. +finfo(float).eps) )

            trid[bmsk] = cval
            
            cnt = cnt +1

        print ""

    # Save some info for later plotting.  Save everything in [cells] except
    # cellsz and origin.
    ext = array( [[(zcen +brad[-1])/cellsz[0],(zcen -brad[-1])/cellsz[0]],
                  [(xcen -brad[-1])/cellsz[2],(xcen +brad[-1])/cellsz[2]]] )
    cdata = {'ncells':ncells,
             'cellsz':cellsz,
             'origin':pd.origin,
             'cntrc':array(cntrc),
             'blkext':ext}

    pivlib.pkldump(cdata,"FLAT-CDATA")

    return trid


def initsrc( src, csdiv, cellsz ):
    # src will be passed in filled with -1.0.  cellsz is valid for the src
    # array.
    return

    ncells = src.shape

    src[:,0,:]  = 32767
    src[0,...]  = 32767
    src[-1,...] = 32767
    src[:,:,0]  = 32767
    src[:,:,-1] = 32767

# >>>>> END USER MODIFIABLE CODE <<<<<

pd = pivlib.loadpivdata(pdifname)

# Advect tracers.
print "Advecting tracers with mpsvtrace()."
ncells = pd[0][velkey].shape[1:4]

src        = empty(ncells,dtype=int16)
src[:,...] = -1
initsrc(src,array((1,1,1)),pd.cellsz)

flolib.mpsvtrace(pd[m_epslc],
                 velkey,
                 icomp,
                 m_ntrpc,
                 m_ncalls,
                 tssdiv,
                 m_irspath,
                 src=src,
                 adjid=m_adjid,
                 hist=m_hist,
                 ctcsp=m_ctcsp,
                 tcfun=tcfun)

