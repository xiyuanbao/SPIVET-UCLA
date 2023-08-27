from numpy import *
from scipy import interpolate
from spivet import pivlib, tlclib, spivetconf
from spivet.steputil import _parsefn, _fretrieve
import sys, os, urlparse
from spivet import compat


def assemble_epoch(accobpath):
    """
    Assembles the individual PIVDATA files that have been accumulated
    in accobpath into one PIVData object.  Each PIVDATA file is
    expected to hole a single Epoch.  Returns the object.
    """
    # Find out how many Epochs need to be processed.
    dl = os.listdir(accobpath)
    dl.sort()
    cp = 0
    for f in range( len(dl) ):
        if ( not dl[cp].startswith('PIVDATA') ):
            dl.pop(cp)
        else:
            cp = cp +1

    # Expects file names to have the form
    #     PIVDATA-E%i.ex2
    _cmprule = lambda x,y: cmp( int(x[0:-4].split('-')[1][1::]),
                                int(y[0:-4].split('-')[1][1::]) )
    dl.sort(_cmprule)

    # Assemble the results.
    for i in range(len(dl)):
        pd = pivlib.loadpivdata( "%s/%s" % (accobpath,dl[i]) )

        if ( i == 0 ):
            epd = pivlib.PIVData(pd.cellsz,pd.origin,pd.desc)

        epd.append(pd[0])
    epd.save("PIVDATA.ex2")
    return epd

accobpath = "./OUTPUT/FG0"
assemble_epoch(accobpath)
