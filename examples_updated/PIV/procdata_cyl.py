"""
Filename:  procdata_cyl.py
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
    Post process PIV results.
"""

from spivet import pivlib
from spivet import flolib
from numpy import *

# >>>>> SETUP <<<<<
pdifname  = "PIVDATA-WFV.ex2"       # Input file for PIVData object.

velkey = "U-MF-SN-ZD-GS"            # Velocity key for vort, divrg, trace.
tssdiv = 3                          # Time steps per Epoch for tracing.

# Passive tracer setup for mpsvtrace() and rpsvtrace().  mpsvtrace()-specific
# parameters are prefixed with m_, while those limitied to rpsvtrace()
# are prefixed by r_.
m_trace = False                    # Advect passive tracers with mpsvtrace.
m_epslc = slice(0,32)               # Epoch slice for processing.
m_ntrpc = 10                      # Number of tracers per cell per call.
m_ncalls = 13                      # Number of calls.
m_adjid = False                    # Force tracers to source ID.
m_tcofname  = "TCRD"                # Output file for tracer coordinates.
m_tidofname = "TRID"                # Output file for tracer ID.

r_trace  = True                    # Advect passive tracers with rpsvtrace.
r_epslc  = slice(0,51)               # Epoch slice for processing.
r_csdiv  = (2,2,2)                  # Composition subdivision factor (z,y,x).
r_cinterp = True                     # Composition interpolation flag.
r_cisep  = None                    # Continuous integration starting Epoch.
r_pdofname=pdifname                # Output file for PIVData object. 
r_interp  = ['L','L']              # Interpolation scheme.

def initcomp( comp, csdiv, cellsz ):
    # comp will be passed in as all zeros.  cellsz is valid for the comp
    # array.
    ncells = comp.shape[1:4]
    sf     = 1 #csdiv[1]
    nlyrs  = ncells[1]/sf
    cmul   = 32767/(nlyrs -1)

    comp[0,:,-sf::,:] = 0.
    for i in range(1,nlyrs):
        sndx = (-i -1)*sf
        endx = sndx +sf
        yslc = slice(sndx,endx)

        xslc = slice(0,ncells[2])

        comp[0,:,yslc,xslc] = i*cmul

def initsrc( src, csdiv, cellsz ):
    # src will be passed in filled with -1.0.  cellsz is valid for the src
    # array.
    sf     = 1 #csdiv[1]
    ncells = src.shape
    nlyrs  = ncells[1]/sf
    cmul   = 32767/(nlyrs -1)

    zpo = 0

    src[:,0,...] = 32767
    src[:,-sf::,:] = 0.

    trad = cellsz[2]*(ncells[2] -1.)/2.
    zcen = cellsz[0]*(ncells[0] -1.)/2.
    xcen = trad
    ndxm = indices((ncells[0],ncells[2]),dtype=float)
    zcrd = cellsz[0]*ndxm[0,...] -zcen
    xcrd = cellsz[2]*ndxm[1,...] -trad

    zcrd = zcrd.reshape(zcrd.size)
    xcrd = xcrd.reshape(xcrd.size)
    rad  = sqrt(zcrd*zcrd +xcrd*xcrd)

    msk  = rad > ( trad -abs(cellsz[0]) )
    zcrd = (compress(msk,zcrd) +zcen)/cellsz[0]
    xcrd = (compress(msk,xcrd) +xcen)/cellsz[2]
    zcrd = zcrd.round().astype(int)
    xcrd = xcrd.round().astype(int)

    for i in range(1,nlyrs):
        sndx = (-i -1)*sf
        endx = sndx +sf
        yslc = slice(sndx,endx)

        """
        src[:, yslc, 0:3] \
            = src[ :, yslc, -3::] \
            = src[ zpo:(zpo+3), yslc,  :] \
            = src[(ncells[0]-zpo-3):(ncells[0]-zpo), yslc,  :] \
            = i*cmul
        """

        src[zcrd,yslc,xcrd] = i*cmul

def initrbndx( ncells ):
    # Will be passed the number of cells in the velocity array.  Must
    # return a 3x2 integer array specifying the rbndx for rpsvtrace().
    # To use the full dataset, return None.
    rbndx = [ [0,ncells[0]],
              [0,ncells[1]],
              [0,ncells[2]] ]

    return rbndx

# Pathlines setup.
ptrace = False                       # Advect passive tracers for pathlines.
epslc  = slice(0,12)                 # Epoch slice for processing.
ids    = [0,1,2]                     # Integer ID's for pathline groups.
pthofname = "PATHLINES-MID"          # Output file name.    
pthdesc   = "PATHLINES"              # File description.

def inittcrd( id, ncells ):
    # Must return an lx3 array containing initial tracer coordinates,
    # with l equal to the number of tracers desired.
    spc = 1

    """
    # For 2D array.
    ndxmat = indices((ncells[0]/spc,ncells[2]/spc),dtype=float)
    ndxmat = ndxmat*spc
    ntrcrs = ndxmat[0,...].size

    zcrd = ndxmat[0,...].reshape(ntrcrs)
    ycrd = empty(ntrcrs,dtype=float)
    xcrd = ndxmat[1,...].reshape(ntrcrs)
    """
    # For 1D array.
    ztrcrs = ncells[0]/spc
    xtrcrs = 0 #ncells[2]/spc
    ntrcrs = ztrcrs +xtrcrs

    zcrd                = empty(ntrcrs,dtype=float)
    zcrd[0:ztrcrs]      = spc*arange(ztrcrs) 
    zcrd[ztrcrs:ntrcrs] = ncells[0]/2

    xcrd                = empty(ntrcrs,dtype=float)
    xcrd[0:ztrcrs]      = ncells[2]/2
    xcrd[ztrcrs:ntrcrs] = spc*arange(xtrcrs)

    ycrd = empty(ntrcrs,dtype=float)

    # Common.
    if ( id == 1 ):
        ycrd[:] = ncells[1] -1
    if ( id == 2 ):
        ycrd[:] = ncells[1] -5        
    if ( id == 3 ):
        ycrd[:] = ncells[1] -10
    if ( id == 4 ):
        ycrd[:] = ncells[1] -30

    tcrd = array([zcrd,ycrd,xcrd])
    return tcrd.transpose()        

# >>>>> END USER MODIFIABLE CODE <<<<<

pd = pivlib.loadpivdata(pdifname)
pd.addsQA("procdata.py")

# Advect tracers.
if ( m_trace ):
    print "Advecting tracers with mpsvtrace()."
    ncells = pd[0][velkey].shape[1:4]

    icomp = pivlib.PIVVar((1,ncells[0],ncells[1],ncells[2]),
                          name="COMPOSITION",
                          units="NA",
                          vtype=pd[0][velkey].vtype,
                          dtype=int16)
    initcomp(icomp,array((1,1,1)),pd.cellsz)

    src        = empty(ncells,dtype=int16)
    src[:,...] = -1
    initsrc(src,array((1,1,1)),pd.cellsz)

    [tcrd,trid] = flolib.mpsvtrace(pd[m_epslc],
                                   velkey,
                                   icomp,
                                   m_ntrpc,
                                   m_ncalls,
                                   tssdiv,
                                   src=src,
                                   adjid=m_adjid)

    fh = open(m_tcofname,'wb')
    tcrd.tofile(fh)
    fh.close()

    fh = open(m_tidofname,'wb')
    trid.tofile(fh)
    fh.close()

if ( r_trace ):
    print "Advecting tracers with rpsvtrace()."
    ncells = pd[0][velkey].shape[1:4]
    cellsz = pd.cellsz/r_csdiv
    rbndx  = array(initrbndx(ncells))
    rsize  = rbndx[:,1] -rbndx[:,0]
    xrsize = r_csdiv*rsize

    icomp = pivlib.PIVVar((1,xrsize[0],xrsize[1],xrsize[2]),
                          name="COMPOSITION",
                          units="NA",
                          vtype=pd[0][velkey].vtype,
                          dtype=int16)
    initcomp(icomp,r_csdiv,cellsz)

    src        = empty(xrsize,dtype=int16)
    src[:,...] = -1
    initsrc(src,r_csdiv,cellsz)

    tpd = flolib.rpsvtrace(pd[r_epslc],
                           velkey,
                           icomp,
                           tssdiv,
                           r_csdiv,
                           src=src,
                           rbndx=rbndx,
                           cinterp=r_cinterp,
                           cisep=r_cisep,
                           interp=r_interp)

    mpdofname = r_pdofname
    if ( r_pdofname == pdifname ):
        tst = (array(r_csdiv) == 1).all()
        if ( tst and ( len(tpd) == len(pd) ) and ( r_cisep == None ) ):
            for e in range( len(pd) ):
                pd[e].addVars( tpd[e]['COMPOSITION'] )
            pd.save(r_pdofname,4)
        else:
            mpdofname = "COMP-%s" % r_pdofname
            tpd.save(mpdofname,4)
    else:
        tpd.save(r_pdofname,4)

    print "Composition saved as: %s" % mpdofname

# Compute pathlines.
if ( ptrace ):
    print "Computing pathlines."
    ncells = pd[0][velkey].shape[1:4]

    xthlst = []
    for id in ids:
        itcrd = inittcrd(id,ncells)
        thlst = flolib.pthtrace(pd[epslc],velkey,itcrd,tssdiv)

        xthlst.append([id,thlst])

    print "Saving pathlines."
    flolib.pthwrite(xthlst,pthofname,pthdesc)

