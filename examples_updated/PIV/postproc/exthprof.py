# WHN 5/26/09.
#
# Extract heater profile from a time-series dataset.

from spivet import pivlib
from numpy import *
import pylab

pylab.ion()

# ----- USER SETUP -----

ifpath = "DAQANALYSIS/HTRTEMP"
ofpath = "HTRPROFILE"

sdt   = 0.005   # Signal sampling period [s]. 
pdt   = 0.5     # Period for profile [s].
olapf = 0.5     # Overlap fraction (percent/100).

# ----- END USER SETUP -----

oset  = int( round( pdt/sdt ) )
ext   = int( oset*olapf/(2.*(1. -olapf)) )
bcntr = int( ceil( oset/2. ) ) +ext
bsize = oset +2*ext

data = pivlib.pklload(ifpath)
npts = data.size

nblks = ( npts -2*ext )/oset
nnpts = nblks*oset +2*ext

print "NBLKS %i" % nblks

# Trim the data.
data = data[0:nnpts]

# Build the profile.
prof = data[0:(nnpts-2*ext):oset]

for i in range( 1, bsize ):
    prof = prof +data[i:(nnpts-2*ext+i):oset]

prof = prof/bsize

# Store the results.
pdict = {"prof":prof,"dt":pdt}
pivlib.pkldump(pdict,ofpath)

# Plot the results.
dta = arange( prof.size )*pdt
fig = pylab.figure()
pylab.plot(dta,prof)
pylab.title("Heater Profile")
pylab.ylabel(r"Temperature [$^{\circ}$C]")
pylab.xlabel("Time [s]")
pylab.savefig("%s-FULL.png" % ofpath)
pylab.close(fig)

fig = pylab.figure()
ndx = int(1200./pdt)
pylab.plot(dta[0:(ndx+1)],prof[0:(ndx+1)])
pylab.title("Heater Profile")
pylab.ylabel(r"Temperature [$^{\circ}$C]")
pylab.xlabel("Time [s]")
pylab.grid()
pylab.savefig("%s-INIT.png" % ofpath)
pylab.close(fig)
