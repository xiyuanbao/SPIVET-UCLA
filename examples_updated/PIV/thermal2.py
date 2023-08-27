# WHN 5/18/09.
#
# Very simple finite volume model to reconstruct the temperature field.
# The method is that of Jameson, Schmidt, and Turkel (AIAA, 1981).  
# The general approach is to use the PIV velocity field and only solve
# the continuity and energy equations.  Continuity must be solved to correct
# for noise in the velocity field which can produce a non-zero velocity
# divergence.  The non-zero divergence results in too much mass, and hence 
# rho*T, being advected into or out of cells.  By solving the continuity
# equation, we can compensate for most (but not all) of the temperature error 
# that would otherwise result from the noise in the velocity field.  This trick
# is implemented as follows:
#    1) Continuity and energy are advanced a full Epoch time step.  At this
#       point each cell contains too much or too little energy due to the
#       spurious rho fluctuation that in turn is caused by velocity field
#       noise.
#    2) The cell energy is divided by the rho from continuity.  This provides
#       the corrected cell temperature.
#    3) A corrected rho is computed using the temperature from step 2 and
#       the syrup's thermal expansivity.
#
# Equations solved (in non-dimensional form):
#    p/pt' integral(rho' dv') = - integral(rho' u'.ds')
#
#    Cp' p/pt' integral(rho' T' dv') = integral(kappa' grad'(T').ds') 
#                                     -integral(Cp' rho' T' u'.ds')
#
# using
#   t     = x0**2/kappa0 t'
#   u     = kappa0/x0    u'
#   x     = x0           x'
#   rho   = rho0         rho'
#   T     = T0           T'
#   Cp    = Cp0          Cp'
#   kappa = kappa0       kappa'

from spivet import pivlib, flolib
from scipy import ndimage
from numpy import *

# ----- USER SETUP -----

# Results will be stored in ifpath provided epslc is the size of the PIVData
# object, otherwise the results will be stored in a file with THRML appended
# to the name given by ifpath.
ifpath  = "PIVDATA-SYNC.ex2"
mofpath = "NTMODEL.ex2"
vname   = "U-MF-SN-ZD-GS"
velsf   = 0.001  # Velocity scalefactor ( mm/s -> m/s ).
cszsf   = 0.001  # pd.cellsz scalefactor ( mm -> m ).

# Set output variable parameters for temperature and rho.  These variables
# will be stored in ifpath using the parameters specified.  Parameters are
# [oset,sf,units], such that the output variable is
#    ovar = sf*(ncv*nvar+oset)
# where ncv is the corresponding entry in ncval below, and nvar is the 
# numerical variable.
ovp = [[-273.,1.,"degC"],[0.,0.001,"NA"]]
#Xiyuan
epcohlast=39
epslc = slice(0,epcohlast+1)  # Epoch slice to process.  Must start at 0.

# Nondimensional constants [rho0,T0,k0,Cp0,x0,u0,t0], where rho is the density,
# T temperature, k thermal conductivity, Cp specific heat capacity, x a length.
# These should be expressed in kg-m-s units.
ncval = {'rho0':1441.,    # Syrup density at 25.2 degC [kg/m^3]
         'T0':25.2 +273., # Bath temp [K].
         'k0':0.35,       # Thermal conductivity of syrup [W/m-K]
         'Cp0':2278.,     # Specific heat capacity of syrup [J/kg-K].
         'x0':0.265}      # Height of the tank [m].

# Bath temperature.  
btemp = ( 25.2 +273. )/ncval['T0']  

# Heater temperature or profile data.  If a temperature is specified, it
# should be non-dimensionalized.  If a path to a dictionary is provided,
# the data in the dictionary should have units of degC for temp and sec for
# time.  Two entries must exist in the dictionary, 'prof' and 'dt'.  'prof',
# a list or array, gives the heater temperature versus time, and 'dt' 
# specifies the period between samples.  hti, the time origin for the 'prof' 
# data should also be specified here (in units of sec).
#
# A few remarks on hti.  
# 1) Epoch 0 as stored in the PIVDATA file captures the quiescent state of 
# the flow and as such is valid until the point in time at which the heater 
# is turned on.  The time value for Epoch 0, as computed by SPIVET after 
# running the synchronize step, is set to the timestamp of Frame 0 of the 
# last Plane (this is t=0).  Note that the SPIV system still has to collect 
# Frame 1, but t=0 nevertheless corresponds to the end of image collection 
# for Frame 0 of Epoch 0.
#
# 2) This code's 'clock' runs off of Epoch time (ie, Epoch 0 is taken as 
# t = 0) as described in 1).  
#
# 3) The data acquisition system is generally activated before the SPIV 
# system, and the heater is turned on shortly after the SPIV system collects 
# the last Epoch 0 image (from Frame 1).  hti, then, represents the time value
# in the heater profile trace that corresponds to t=0 as described above. 
htemp = "HTRPROFILE"  # pivlib.pkldump file containing the profile dictionary.
hti   = 223. +40.     # ***** UPDATE EACH DATASET.  Set at end of E0.

# Need to build the tank bottom since it has a big impact on how heat is
# transfered to (from) the syrup.  The thickness is specified by the equivalent
# number of cells that it would span using pivdata.cellsz.  This cell size
# is called the coarse cell size by the code.
#
# Also need to specify the heater radius and the cell which roughly corresponds
# to the heater's center.  Note, the center cell is most easily determined 
# using a tracer field or perhaps the velocity field in Paraview or Visit.
#
# The heater in the UM tank doesn't span the full thickness (along y) of the
# tank bottom.  chtrsy permits specification of which y-cells are actually
# occupied by the heater.
# 
# The SPIVET coordinate system is used here.  z points away from the camera,
# y points vertically down.
cnbtmcy = 11                # Coarse cell count for tank bottom along y-axis.
chtrsy  = slice(-cnbtmcy,
                -cnbtmcy+2) # Coarse heater cell slice along y-axis.
#Xiyuan: heater position may be problematic
chtrcc  = [21,29]           # Heater center cell in coarse mesh [k,i].
htrrad  = 0.01/ncval['x0']  # Heater radius (non-dimensionalized).

# Specify cell and timestep subdivision factors.  
csdiv  = [2,1,1]           # Cell subdivision factors [k,j,i] (for fine mesh).
tssdiv = 10                # Time step subdivision.

# These are the artificial dissipation constants.  They are heuristically
# determined by experimentation and should not need to be modified unless
# oscillations are apparent.  The first value controls the dissipation
# from second order differences, and the second term controls the dissipation
# from fourth order differences.  See the Jameson:1981 paper mentioned
# above for more details.
adc = [1./8.,1./256.]       # Artificial disspation constants.

# Material properties.
def matprop(mesh,htrc,nbtmcy,ncval):
    """
    ----
    
    mesh            # Mesh dictionary.  See bldmesh() for details.
    htrc            # lx3 list of heater cell indices.  See gethtrc().
    nbtmcy          # Number of bottom cells (fine) for tank bottom.
    ncval           # Non-dimensionalization constant dictionary.

    ----
    
    Set material properties as needed.  

    This code was not designed to handle strong variations in material
    properties (it can't handle copper and plexiglass at the same time,
    for example).  Consequently, sharp discontinuities in materials
    may require increased artificial dissipation to decrease oscillations
    or may cause the code to become unstable.

    All cells must have rhoc*, k*, and Cp*.  rhoc* defines a linear density
    versus temperature mapping.  rhoc* = [c0,c1] such that rho* = c0 +c1*(T*).
    k* is the non-dimensionalized thermal conductivity, and Cp* is the
    non-dimensionalized specific heat capacity.

    Returns [rhoc*,k*,Cp*].
    """
    # Initialization.
    mshape = mesh['shape']
    pshape = array([1,mshape[0],mshape[1],mshape[2]])

    rho0 = ncval['rho0']
    T0   = ncval['T0']
    k0   = ncval['k0']
    Cp0  = ncval['Cp0']

    # k* = k/k0, Cp* = Cp/Cp0.
    kstar  = empty(pshape,dtype=float)    
    cpstar = empty(pshape,dtype=float)

    # rhoc must be non-dimensionalized.  rho = rhoc[0] +rhoc[1]*(T*)
    # where T* is the non-dimensional absolute temperature.
    rhoc   = []
    rhoc.append( empty(pshape,dtype=float) )
    rhoc.append( empty(pshape,dtype=float) )

    # Syrup.  
    kstar[   0, :, 0:(mshape[1] -nbtmcy), :] = 1.
    cpstar[  0, :, 0:(mshape[1] -nbtmcy), :] = 1.
    rhoc[0][ 0, :, 0:(mshape[1] -nbtmcy), :] = (1452. +273.*0.447)/rho0
    rhoc[1][ 0, :, 0:(mshape[1] -nbtmcy), :] = -0.447*T0/rho0    

    # Acrylic.  Thermal conductivity of syrup is almost identical to that
    # of acrylic.  So we set all acrylic properties equal to those of
    # syrup at 25 degC.
    kstar[   0, :, (mshape[1]-nbtmcy):mshape[1], :] = 1.
    cpstar[  0, :, (mshape[1]-nbtmcy):mshape[1], :] = 1.
    rhoc[0][ 0, :, (mshape[1]-nbtmcy):mshape[1], :] = 1440./rho0
    rhoc[1][ 0, :, (mshape[1]-nbtmcy):mshape[1], :] = 0.    

    # Annular gap between heater and false bottom.  Provides an air gap
    # between the heater and the false bottom to quasi-match the lab setup.
    ancy = htrc[:,1].max() -htrc[:,1].min() +1 # cells thick along y.

    chcells = htrc[ htrc[:,1] == mshape[1] -nbtmcy ]

    andx = indices((mshape[0],ancy,mshape[2]))  
    andx[1,...] = andx[1,...] +mshape[1] -nbtmcy 

    cssec = zeros((mshape[0],mshape[2]),dtype=int)
    cssec[chcells[:,0],chcells[:,2]] = 1

    tmplt = zeros((3,3),dtype=int)
    tmplt[:,1] = 1
    tmplt[1,:] = 1

    ci  = ndimage.convolve(cssec,tmplt)

    msk = ci >= 1
    msk[chcells[:,0],chcells[:,2]] = False
    msk = msk.reshape(msk.shape[0],1,msk.shape[1])
    msk = msk.repeat(ancy,1)

    nacells = msk.sum()

    annc = []
    for i in xrange(3):
        annc.append( andx[i,...][msk].reshape(nacells) )
    annc = array(annc)
    annc = annc.transpose()
    
    kstar[   0, annc[:,0], annc[:,1], annc[:,2] ] = 0.05
    cpstar[  0, annc[:,0], annc[:,1], annc[:,2] ] = 1.
    rhoc[0][ 0, annc[:,0], annc[:,1], annc[:,2] ] = 1440./rho0
    rhoc[1][ 0, annc[:,0], annc[:,1], annc[:,2] ] = 0.

    # The UM heater has a low thermal conductivity plug that mounts the heater
    # to the tank.  So heat transfer from the bottom of the heater to the
    # tank bottom is small.  Code below simulates this plug below the heater.
    ymx = htrc[:,1].max()

    kstar[   0, htrc[:,0], (ymx+1)::, htrc[:,2] ] = 0.05 
    cpstar[  0, htrc[:,0], (ymx+1)::, htrc[:,2] ] = 1. 
    rhoc[0][ 0, htrc[:,0], (ymx+1)::, htrc[:,2] ] = 1440./rho0 
    rhoc[1][ 0, htrc[:,0], (ymx+1)::, htrc[:,2] ] = 0.

    # Copper.  As mentioned above, the code can't handle copper and syrup
    # (thermal diffusion in copper is much too fast).  We set the heater
    # properties to those of syrup and force these cells to have the heater
    # temperature.
    kstar[   0, htrc[:,0], htrc[:,1], htrc[:,2] ] = 1. 
    cpstar[  0, htrc[:,0], htrc[:,1], htrc[:,2] ] = 1. 
    rhoc[0][ 0, htrc[:,0], htrc[:,1], htrc[:,2] ] = 1440./rho0 
    rhoc[1][ 0, htrc[:,0], htrc[:,1], htrc[:,2] ] = 0.

    return [rhoc,kstar,cpstar]

# Boundary conditions.
def getbcs(mesh,htrc,nbtmcy,ncval):
    """
    ----
    
    mesh            # Mesh dictionary.  See bldmesh() for details.
    htrc            # lx3 list of heater cell indices.  See gethtrc().
    nbtmcy          # Number of bottom cells (fine) for tank bottom.
    ncval           # Non-dimensionalization constant dictionary.

    ----

    Set boundary conditions as needed.  Two types of boundary conditions
    are supported: constant temperature and constant velocity.

    Boundary conditions must be of the type FixedValueBC.  See documentation
    below on that class for details.

    Returns [tempbc,htempbc,velbc], where tempbc is for general temperature,
    htempbc is for the heater, and velbc is velocity BC's.
    """
    # Initialization.
    mshape  = mesh['shape']

    # Cell indices.
    ndxm    = indices((mshape[0],nbtmcy,mshape[2]),dtype=int)
    ndxm[1] = ndxm[1] +mshape[1] -nbtmcy

    # Grab all cells that are not heater cells.
    rhndx      = htrc.copy()  # Points into ndxm array.
    rhndx[:,1] = rhndx[:,1] -mshape[1] +nbtmcy

    nhmsk = ones(ndxm[0].shape,dtype=bool)
    nhmsk[rhndx[:,0],rhndx[:,1],rhndx[:,2]] = False

    ndxm  = ndxm.reshape((3,ndxm.size/3))
    ndxm  = ndxm.transpose()
    nhmsk = nhmsk.reshape(ndxm.shape[0])
    
    # Non-heater bottom cells.
    bcells = ndxm[nhmsk,:]

    # Get tank bottom edge cell indices along y-axis.
    bemsk = bcells[:,1] == ( mshape[1] -1 )
    bedge = bcells[bemsk,:]
    
    # Get a ring of cells around the top of the domain (ie, tank).  The 
    # artificial dissipation scheme doesn't reach into this ring (because of
    # the simplistic way the derivatives are implemented), so oscillations can 
    # occur.  We'll set a constant temp BC on this ring.
    tring = indices((mshape[0],2,mshape[2]),dtype=int)
    tring = tring.reshape(3,tring[0].size).transpose()
    msk   = ( ( tring[:,0] < 2 ) + ( tring[:,0] >= mshape[0] -2 ) ) \
           +( ( tring[:,2] < 2 ) + ( tring[:,2] >= mshape[2] -2 ) )
    tring = tring[msk,:]

    # Get corner columns.  As with the upper ring, artifical disspipation 
    # doesn't reach here.
    cndx = indices((mshape[0],mshape[1],mshape[2]),dtype=int)
    cndx = cndx.reshape(3,cndx[0].size).transpose()

    zmsk = ( cndx[:,0] < 2 ) + ( cndx[:,0] >= mshape[0] -2 )
    xmsk = ( cndx[:,2] < 2 ) + ( cndx[:,2] >= mshape[2] -2 )
    cmsk = zmsk*xmsk

    cndx = cndx[cmsk,:]

    # Setup tempbc's.
    tempbc = []
    tempbc.append( FixedValueBC(bedge,0,btemp) )
    tempbc.append( FixedValueBC(tring,0,btemp) )
    tempbc.append( FixedValueBC(cndx,0,btemp) )

    # Get heater profile.
    htempbc = []
    if ( isinstance(htemp,str) ):
        pdict = pivlib.pklload(htemp)
        prof  = (pdict['prof'] +273.)/ncval['T0']  
        dt    = pdict['dt']/ncval['t0']
        ti    = hti/ncval['t0']

        htempbc.append( FixedValueBC(htrc,0,prof,ti,dt) )
    else:
        tempbc.append( FixedValueBC(htrc,0,htemp) )

    # Setup velbc's.
    velbc = []
    velbc.append( FixedValueBC(bcells,slice(0,3),0.) )
    velbc.append( FixedValueBC(htrc,slice(0,3),0.) )

    return [tempbc,htempbc,velbc]

# ----- END USER SETUP -----

class FixedValueBC:
    """
    General class for fixed-value boundary conditions.
    """
    def __init__(self,cndx,cmp,value,ti=0.,dt=1.):
        """
        ----

        cndx        # lx3 array of cell indices.
        cmp         # Variable component affected by BC.
        value       # Value to set for BC.
        ti=0.       # Initial time for time-varying value series.
        dt=1.       # Time increment between values in series.

        ----

        The variable component affected by the BC can be either a scalar
        value or a python slice object.

        Value can be a 1D list of values representing a time varying series.
        In this case, set ti to the time origin for the series and dt to
        the period between series values.  Both ti and dt should be 
        non-dimensionalized.
        """
        self.cndx  = cndx
        self.cmp   = cmp
        self.value = value
        self.ti    = ti
        self.dt    = dt

    def apply(self,var):
        """
        Apply the boundary condition.  All variables must be shaped like a
        PIVVar.  Variable is modified in place.
        """
        # Initialization.
        cndx  = self.cndx
        cmp   = self.cmp
        value = self.value

        # Application
        var[cmp,cndx[:,0],cndx[:,1],cndx[:,2]] = value

    def tapply(self,var,t):
        """
        Apply the boundary condition at time t.  All variables must be shaped
        like a PIVVar.  Variable is modified in place.
        """
        # Initialization.
        cndx  = self.cndx
        cmp   = self.cmp
        ndx   = int( round( (t +self.ti)/self.dt ) )
        value = self.value[ndx]

        # Application.
        var[cmp,cndx[:,0],cndx[:,1],cndx[:,2]] = value


def gethtrc(htrcc,htrrad,htrsy,mesh):
    """
    ----

    htrcc           # Cell location of heater center (fine) [k,i].
    htrrad          # Non-dimensionalized heater radius.
    htrsy           # Heater cell y-slice (fine).
    mesh            # Mesh dictionary.  See bldmesh() for details.

    ----

    Returns an lx3 list specifying the k,j,i indices of a circular heater
    located at htrcc.  
    """
    # Initialization.
    mshape = mesh['shape']
    cellsz = mesh['cellsz']

    shtrcc = array(htrcc)*array([cellsz[0],cellsz[2]])

    # Compute radii.
    ndxm  = indices(mshape,dtype=int)
    sndxz = cellsz[0]*ndxm[0] -shtrcc[0]
    sndxx = cellsz[2]*ndxm[2] -shtrcc[1]

    rad = sqrt( sndxz**2 +sndxx**2 )

    hmsk            = zeros(mshape,dtype=bool)
    hmsk[:,htrsy,:] = True
    hmsk            = hmsk*( rad -htrrad <= 0. )

    nhtrc = hmsk.sum()
    print "Heater cells: %i" % nhtrc

    # Compute indices.
    htrc = []
    for i in xrange(3):
        htrc.append( ndxm[i,...][hmsk].reshape(nhtrc) )
    htrc = array(htrc)
    htrc = htrc.transpose()

    return htrc

def bldmesh(pd,cszsf,csdiv,nbtmcy,ncval):
    """
    ----

    pd              # PIVData object containing input data.
    cszsf           # Cell size scale factor.
    csdiv           # Cell subdivision factor.
    nbtmcy          # Number of bottom cells (fine) for tank bottom.
    ncval           # Non-dimensionalization constant dictionary.    

    ----
    Construct a mesh dictionary.

    Returns mesh, a dictionary.
    """
    # Initialization.
    csdiv = array(csdiv)
    eshape = array(pd[0].eshape)
    cellsz = cszsf*array(pd.cellsz)/(ncval['x0']*csdiv)

    mshape  = csdiv*eshape +array([0,nbtmcy,0])
    tncells = mshape.prod()

    print "mshape (k,j,i): (%i,%i,%i)" % tuple(mshape)

    # Face areas and index orientation.  If the index orientation is
    # positive, then the positive axis and the direction of increasing
    # index are aligned.  For our setup, the positive z-axis points in the
    # direction of decreasing k-index.
    farea = array([ cellsz[1]*cellsz[2],    # k-face
                    cellsz[0]*cellsz[2],    # j-face
                    cellsz[0]*cellsz[1] ])  # i-face
    farea = abs(farea)

    ndxor = sign(cellsz)

    # Cell volume.
    cvol = abs( cellsz.prod() )

    mesh = {'shape':mshape,
            'cellsz':cellsz,
            'farea':farea,
            'ndxor':ndxor,
            'cvol':cvol}

    return mesh

def sdivvel(pd,vname,velsf,u0,csdiv,mesh):
    """
    ----

    pd              # PIVData object containing input data.
    vname           # Name of the velocity variable in pd.
    velsf           # Velocity scale factor.
    u0              # Velocity non-dimensionalization constant.
    csdiv           # Cell subdivision factor.
    mesh            # Mesh dictionary.  See bldmesh() for details.

    ----

    Subdivide the velocity field and copy results into an array of correct
    shape.  Also normalize velocity by u0.

    Returns vel, a list of all velocities that will be used in the
    calculations.
    """
    # Initialization.
    csdiv = array(csdiv)

    mshape   = mesh['shape']
    tnmcells = mshape.prod()

    # Get coordinates for the cells.
    crd = indices(mshape,dtype=float)
    for i in xrange(3):
        crd[i,...] = crd[i,...]/csdiv[i] -(csdiv[i] -1.)/(2.*csdiv[i])
    
    crd = crd.reshape(3,tnmcells).transpose()

    # Interpolate the velocity field.
    vel = []
    for e in pd:
        var  = e[vname]
        ivar = flolib.svinterp(var,crd)
        ivar = ivar.transpose().reshape([3,mshape[0],mshape[1],mshape[2]])

        vel.append(velsf*ivar/u0)

    return vel

def initT(btemp,mesh):
    """
    ----

    btemp           # Non-dimensionalized bath temperature.
    mesh            # Mesh dictionary.  See bldmesh() for details.

    ----

    Initialize the temperature field.

    Returns temp, a properly inialized temperature field.
    """
    mshape = mesh['shape']

    temp        = empty((1,mshape[0],mshape[1],mshape[2]),dtype=float)
    temp[:,...] = btemp

    return temp

def rhofun(temp,rhoc):
    """
    ----

    temp            # Temperature field.
    rhoc            # Density function constants.

    ----

    Compute density based on temperature using:
        rho = rhoc[0] -rhoc[1]*temp

    Returns rho.
    """
    return rhoc[0] +rhoc[1]*temp

def getfv(var,axis):
    """
    ---

    var             # Variable to compute at face centers.
    axis            # Spatial axis along which face values are computed.

    ----

    Computes values of var at face centers.  The value at the i+1/2 face will 
    be set to 
        fv[i+1/2] = ( var[i] +var[i+1] )/2.

    The outermost cells faces (ie, those on the domain perimeter) will be 
    assiged a face value equal to that of the parent cell.  This is equivalent 
    to a Neumann BC where the gradient is set to zero.

    Returns fvar, the values at face centers along the specified axis.  fvar
    will have the same shape as var except along axis, where the size of fvar
    will be one larger than var.
    """
    # Initialization.
    vshape          = array( var.shape )
    xvshape         = vshape.copy()
    xvshape[axis+1] = xvshape[axis+1] +1

    fvar = empty(xvshape,dtype=float)

    # Get slices for perimeter cells on the 0 face and copy over the values.
    slc = []
    for i in xrange(3):
        slc.append( slice(0, vshape[i+1]) )

    slc[axis] = slice( 0, 1 )

    fvar[:,slc[0],slc[1],slc[2]] = var[:,slc[0],slc[1],slc[2]]

    # Get slices for perimeter cells on the max face and copy over the values.
    slc  = []
    fslc = []
    for i in xrange(3):
        slc.append(  slice(0, vshape[i+1]) )
        fslc.append( slice(0, vshape[i+1]) )

    slc[axis]  = slice( vshape[axis+1] -1,  vshape[axis+1]  )
    fslc[axis] = slice( xvshape[axis+1] -1, xvshape[axis+1] )

    fvar[:,fslc[0],fslc[1],fslc[2]] = var[:,slc[0],slc[1],slc[2]]

    # Get slices for interior cells.  Slices are ordered as minus, plus for the
    # axis.
    slc  = []
    fslc = []
    for i in xrange(3):
        slc.append( slice(0, vshape[i+1]) )
        slc.append( slice(0, vshape[i+1]) )

        fslc.append( slice(0, vshape[i+1]) )

    slc[2*axis]   = slice( 0, vshape[axis+1] -1  )
    slc[2*axis+1] = slice( 1, vshape[axis+1] )

    fslc[axis] = slice( 1, xvshape[axis+1] -1 )

    # Computation.
    fvar[:,fslc[0],fslc[1],fslc[2]] = 0.5*( var[:,slc[0],slc[2],slc[4]] \
                                           +var[:,slc[1],slc[3],slc[5]] )

    return fvar

def fluxcomp(var,mesh):
    """
    ----

    var             # Flux variable.
    mesh            # Mesh dictionary.  See bldmesh() for details.

    ----

    Computes a flux integral of the form int( var . ds ) for all cells
    in the domain.

    Returns flux, the residual flux into or out of each cell.
    """
    # Initialization.
    mshape = mesh['shape']
    farea  = mesh['farea']
    ndxor  = mesh['ndxor']

    vshape = var.shape

    # Get values at face centers.  fv will be 1 entry larger than the domain
    # along each axis.
    fv = []
    for i in xrange(3):
        fv.append( farea[i]*getfv(var,i)[i,...] )

    # Compute the cell fluxes along each axis and then combine.  The outer
    # loop is over each axis.  The inner loop builds the slices.  For each
    # axis, a special set of slices must be created to run over the cells
    # along that axis.  
    flux = []
    for i in xrange(3):
        slc = []
        for s in xrange(3):
            slc.append( slice(0, vshape[s+1]) )
            slc.append( slice(0, vshape[s+1]) )

        # Replace the axis slices with a special, staggered set.
        slc[2*i]   = slice( 0, vshape[i+1]    )
        slc[2*i+1] = slice( 1, vshape[i+1] +1 )
        
        flux.append( fv[i][slc[1],slc[3],slc[5]] -fv[i][slc[0],slc[2],slc[4]] )

    flux = ndxor[0]*flux[0] +ndxor[1]*flux[1] +ndxor[2]*flux[2]

    return flux.reshape([1,vshape[1],vshape[2],vshape[3]])

def extrhoT(flovar,rhoc):
    """
    ----

    flovar          # The flow variable array.
    rhoc            # Density function constants.

    ----

    flovar contains the rho computed using the continuity equation and
    rhoT determined via the energy equation.  extrhoT() computes the
    corrected temperature and then uses that temperature to compute a 
    corrected rho.

    Returns [rho, T], where rho is the correct rho based on temperature of
    fluid.
    """
    # Extract T and compute rho.
    temp = flovar[1,...]/flovar[0,...]
    rho  = rhofun(temp,rhoc)

    return [rho,temp]

def artdiss(rhoT,adc,mesh,dt):
    """
    ----

    rhoT            # rhoT component of flovar.
    adc             # Artificial dissipation constants.
    mesh            # Mesh dictionary.  See bldmesh() for details.
    dt              # Time step.

    ----

    Applies artificial dissipation using the JST scheme (Jameson:1981).

    dt should be set to st (step time) from the RK4 scheme.

    Returns diss, the dissipation for rhoT that can be added directly to
    the RHS after the RHS has been normalized by cell volume and Cp and diss
    has been divided by cell volume.
    """
    # Initialization.
    rtshape  = array( rhoT.shape )
    xrtshape = rtshape.copy()

    nfac = mesh['cvol']/dt

    xrtshape[1::] = xrtshape[1::] +1

    switch = zeros([3,1,rtshape[1],rtshape[2],rtshape[3]],dtype=float)
    diss   = [ zeros([1,xrtshape[1], rtshape[2], rtshape[3]],dtype=float),
               zeros([1, rtshape[1],xrtshape[2], rtshape[3]],dtype=float),
               zeros([1, rtshape[1], rtshape[2],xrtshape[3]],dtype=float) ]

    # All slices for a given axis will be ordered from max to min.
    
    # Build the switch.  The switch is only needed for cells 1:-1 (see below).
    # The switch for cell i is given by
    #     switch = abs( rhoT[i+1] -2*rhoT[i] +rhoT[i-1] )
    #              --------------------------------------
    #                   rhoT[i+1] +2*rhoT[i] +rhoT[i-1]
    for i in xrange(3):
        slc  = [[],[],[]]
        for ax in xrange(3):
            for s in xrange(3):
                slc[ax].append( slice(1,rtshape[ax+1]-1) )

        for s in xrange(3):
            slc[i][s] = slice(2-s,rtshape[i+1]-s)

        switch[i,0,slc[0][1],slc[1][1],slc[2][1]] = \
            abs (     rhoT[0,slc[0][0],slc[1][0],slc[2][0]]
                  -2.*rhoT[0,slc[0][1],slc[1][1],slc[2][1]]
                     +rhoT[0,slc[0][2],slc[1][2],slc[2][2]] ) \
               /(     rhoT[0,slc[0][0],slc[1][0],slc[2][0]]
                  +2.*rhoT[0,slc[0][1],slc[1][1],slc[2][1]]
                     +rhoT[0,slc[0][2],slc[1][2],slc[2][2]] ) \

    # Compute the dissipation values at the faces.  The dissipation will only
    # be computed for cells 2:-2 (just for ease in handling the boundaries).
    # The disspation for the remaining cells will be set to zero.
    for i in xrange(3):
        eslc = [[],[],[]]
        rslc = [[],[],[]]
        for ax in xrange(3):
            for s in xrange(2):
                eslc[ax].append( slice(2,rtshape[ax+1]-2) )
            for s in xrange(4):
                rslc[ax].append( slice(2,rtshape[ax+1]-2) )

        for s in xrange(2):
            eslc[i][s] = slice(2-s,rtshape[i+1]-1-s)

        # This set of slices corresponds to the terms (in order):
        #   rhoT[i+2], rhoT[i+1], rhoT[i], rhoT[i-1]
        for s in xrange(4):
            rslc[i][s] = slice(3-s,rtshape[i+1]-s)

        eps2 = adc[0]*maximum(switch[i,0,eslc[0][0],eslc[1][0],eslc[2][0]],
                              switch[i,0,eslc[0][1],eslc[1][1],eslc[2][1]])
        eps4 = maximum(0.,adc[1] -eps2)

        cmp2 = eps2*( rhoT[0,rslc[0][1],rslc[1][1],rslc[2][1]] \
                     -rhoT[0,rslc[0][2],rslc[1][2],rslc[2][2]] )

        cmp4 = eps4*(    rhoT[0,rslc[0][0],rslc[1][0],rslc[2][0]] \
                     -3.*rhoT[0,rslc[0][1],rslc[1][1],rslc[2][1]] \
                     +3.*rhoT[0,rslc[0][2],rslc[1][2],rslc[2][2]] \
                        -rhoT[0,rslc[0][3],rslc[1][3],rslc[2][3]] )

        diss[i][0,2:-2,2:-2,2:-2] = cmp2 -cmp4

    # Compute the residual.
    diss[0] = diss[0][0,1::,:,:] -diss[0][0,0:-1,:,:]
    diss[1] = diss[1][0,:,1::,:] -diss[1][0,:,0:-1,:]
    diss[2] = diss[2][0,:,:,1::] -diss[2][0,:,:,0:-1]

    return nfac*(diss[0] +diss[1] +diss[2])

def rk4fun(flovar,st,fargs):
    """
    RK4 RHS function for rhoT.

    fargs[0] ---- U at time t(i).
    fargs[1] ---- U at time t(i+1).
    fargs[2] ---- t(i)
    fargs[3] ---- rt, where rt = t -t(i).
    fargs[4] ---- dt, where dt = t(i+1) -t(i).
    fargs[5] ---- [rhoc,kstar,cpstar,ncval]
    fargs[6] ---- mesh
    fargs[7] ---- [tempbc,htempbc]
    fargs[8] ---- adc
    """
    # Initialization.
    rhoc   = fargs[5][0]
    kstar  = fargs[5][1]
    cpstar = fargs[5][2]
    ncval  = fargs[5][3]
    mesh   = fargs[6]

    cvol   = mesh['cvol']
    cellsz = mesh['cellsz']
    mshape = mesh['shape']

    rhs = empty([2,1,mshape[0],mshape[1],mshape[2]],dtype=float)

    # Interpolate the velocity field in time.
    srt  = fargs[3] +st
    asrt = abs(srt)

    if ( asrt > 0. ):
        ft   = (srt)/fargs[4] 
        ivel = (1. -ft)*fargs[0] +ft*fargs[1]
    else:
        ivel = fargs[0]

    # Use the rho computed by continuity for RK4 steps.  Will compute a
    # temperature based rho later.
    temp = flovar[1,...]/flovar[0,...]
    rho  = flovar[0,...]

    # ----- CONTINUITY -----
    # Compute int( rho u . ds ).
    rhs[0,...] = -1./cvol*fluxcomp(rho*ivel,mesh)    
    if ( st > finfo(float).eps and adc != None ):    
        diss       = artdiss(flovar[0,...],adc,mesh,st)
        rhs[0,...] = rhs[0,...] +1./cvol*diss

    # ----- ENERGY -----
    # Compute grad(T).
    gradT = []
    for i in xrange(3):
        gradT.append( flolib.d_di(temp,i,cellsz[i])[0,...] )

    gradT = array( gradT )

    # Compute integral( k gradT . ds ).
    gtflux = fluxcomp(kstar*gradT,mesh)

    # Compute integral( Cp rho T u . ds ).
    rtu     = flovar[1,...]*ivel
    rtuflux = fluxcomp(cpstar*rtu,mesh)

    # Compute the RHS.
    cvcpr      = 1./(cvol*cpstar)
    rhs[1,...] = cvcpr*(gtflux -rtuflux)

    if ( st > finfo(float).eps and adc != None ):
        diss  = cpstar*artdiss(flovar[1,...],adc,mesh,st)
        dmin  = diss.min()
        dmax  = diss.max()
        dmean = abs(diss).mean()

        rhs[1,...] = rhs[1,...] +cvcpr*diss
    else:
        dmin  = 0.
        dmax  = 0.
        dmean = 0.

    fmt = "| %10.3e | %10.3e %10.3e %10.3e | %10.3e %10.3e %10.3e | " \
         +"%10.3e %10.3e %10.3e | %10.3e %10.3e %10.3e"
    print fmt % (st,
                 gtflux.min(),gtflux.max(),abs(gtflux).mean(),
                 rtuflux.min(),rtuflux.max(),abs(rtuflux).mean(),
                 dmin,dmax,dmean,
                 rhs.min(),rhs.max(),abs(rhs).mean())

    return rhs

def desdiv(var,csdiv):
    """
    ----
    
    var             # Variable to de-subdivide.
    csdiv           # Cell subdivision factors.

    ----
    
    De-subdivides a variable.

    Returns dvar.
    """
    # Initialization.
    csdiv = array(csdiv)
    
    # De-subdivide.
    dvar = var[:,0::csdiv[0],:,:]
    for i in xrange( 1, csdiv[0] ):
        dvar = dvar +var[:,i::csdiv[0],:,:]
    var = dvar

    dvar = var[:,:,0::csdiv[1],:]
    for i in xrange( 1, csdiv[1] ):
        dvar = dvar +var[:,:,i::csdiv[1],:]
    var = dvar

    dvar = var[:,:,:,0::csdiv[2]]
    for i in xrange( 1, csdiv[2] ):
        dvar = dvar +var[:,:,:,i::csdiv[2]]

    dvar = dvar/csdiv.prod()

    return dvar

# Initialization.
pd     = pivlib.loadpivdata(ifpath)
useipd = True
#Xiyuan
print epslc.stop,epslc.start,len(pd)
if ( epslc.stop -epslc.start == len(pd) ):
    opd = pd
else:
    opd = pivlib.PIVData(pd.cellsz,pd.origin,"RECONSTRUCTED TEMPERATURE")

    useipd = False

pd = pd[epslc]

csdiv = array(csdiv)
oos   = pd.cellsz*(csdiv -1.)/(2.*csdiv)
npd   = pivlib.PIVData(pd.cellsz/csdiv,pd.origin-oos,"NUMERICAL MODEL")

nbtmcy   = csdiv[1]*cnbtmcy
htrsy    = slice(csdiv[1]*chtrsy.start,csdiv[1]*chtrsy.stop)
htrcc    = array(chtrcc)
htrcc[0] = csdiv[0]*htrcc[0]
htrcc[1] = csdiv[2]*htrcc[1]  

ncval['t0'] = ncval['rho0']*ncval['Cp0']*ncval['x0']**2/ncval['k0']
ncval['u0'] = ncval['x0']/ncval['t0']
#Xiyuan
print "u0=",ncval['u0']
# Setup the mesh dictionary.
mesh   = bldmesh(pd,cszsf,csdiv,nbtmcy,ncval)
mshape = mesh['shape']
print mesh

# Get heater cells.
htrc = gethtrc(htrcc,htrrad,htrsy,mesh)

# Get material properties and boundary conditions.
[rhoc,kstar,cpstar]    = matprop(mesh,htrc,nbtmcy,ncval)
[tempbc,htempbc,velbc] = getbcs(mesh,htrc,nbtmcy,ncval)

# Get subdivided velocities and apply velocity BC's.  This only needs to be
# done once.
vel = sdivvel(pd,vname,velsf,ncval['u0'],csdiv,mesh)
for e in xrange( len(vel) ):
    for b in velbc:
        b.apply(vel[e])

# Initialize T and compute rho.
temp = initT(btemp,mesh)
for b in tempbc:
    b.apply(temp)

for b in htempbc:
    b.tapply(temp,0.)

rho  = rhofun(temp,rhoc)

# Flow variables will be stored as a large array ordered as [frho, rhoT],
# where frho is the fictitous rho generated by a non-zero velocity divergence.
flovar = empty([2,1,mshape[0],mshape[1],mshape[2]],dtype=float)

flovar[0,...] = rho
flovar[1,...] = rho*temp

npd.addVars( 0, [pivlib.cpivvar(ncval['T0']*temp,"T-NU"),
                 pivlib.cpivvar(ncval['rho0']*rho,"RHO-NU"),
                 pivlib.cpivvar(ncval['u0']*vel[0],"U-NU")] )

lr_temp = temp[:,:,0:(mshape[1]-nbtmcy),:]
lr_temp = ovp[0][1]*(ncval['T0']*desdiv(lr_temp,csdiv) +ovp[0][0])

lr_rho  = rho[:,:,0:(mshape[1]-nbtmcy),:]
lr_rho  = ovp[1][1]*(ncval['rho0']*desdiv(lr_rho,csdiv) +ovp[1][0])

opd.addVars(0, [pivlib.cpivvar(lr_temp,"T-NU",ovp[0][2]),
                pivlib.cpivvar(lr_rho,"RHO-NU",ovp[1][2])] )

sstr = "********************************************************************"\
      +"********************************************************************"\
      +"****"

# Main loop.  Note: To visualize how the code is executed, at each Epoch,
# et, the code is stepping FROM the previous Epoch (et -1) to the current
# (et).
print "Computing ..."
fargs    = range(9)
fargs[5] = [rhoc,kstar,cpstar,ncval]
fargs[6] = mesh
fargs[7] = [tempbc,htempbc]
fargs[8] = adc
for et in xrange( 1, len(pd) ):
    print "| Epoch %3i %s" % (et,sstr)
    edt = (pd[et].time -pd[et -1].time)/ncval['t0']
    idt = edt/tssdiv

    fargs[0] = vel[et -1]
    fargs[1] = vel[et]
    fargs[2] = (pd[et -1].time -pd[0].time)/ncval['t0'] # t(et -1)
    fargs[4] = edt

    for it in xrange( tssdiv ):
        fargs[3] = it*idt
        ittime   = fargs[2] +fargs[3]

        fmt = "| IT %3i --- "\
             +"| DFUSN -------------------------- "\
             +"| ADVEC -------------------------- "\
             +"| ADISS -------------------------- "\
             +"| RESID --------------------------" 

        print fmt % it
        print "| ST         " \
             +"| MIN        MAX        MEAN (ABS) " \
             +"| MIN        MAX        MEAN (ABS) " \
             +"| MIN        MAX        MEAN (ABS) " \
             +"| MIN        MAX        MEAN (ABS)" \

        flovar = flolib.floutil.rk4(flovar,idt,rk4fun,fargs)
    
        # Extract rho and T.
        [rho,temp] = extrhoT(flovar,rhoc)

        for b in tempbc:
            b.apply(temp)

        for b in htempbc:
            b.tapply(temp,ittime)

        # Update flovar.
        flovar[0,...] = rho
        flovar[1,...] = rho*temp

    # Store the results.
    npd.addVars( et, [pivlib.cpivvar(ncval['T0']*temp,"T-NU"),
                      pivlib.cpivvar(ncval['rho0']*rho,"RHO-NU"),
                      pivlib.cpivvar(ncval['u0']*vel[et],"U-NU")] )

    lr_temp = temp[:,:,0:(mshape[1]-nbtmcy),:]
    lr_temp = ovp[0][1]*(ncval['T0']*desdiv(lr_temp,csdiv) +ovp[0][0])

    lr_rho  = rho[:,:,0:(mshape[1]-nbtmcy),:]
    lr_rho  = ovp[1][1]*(ncval['rho0']*desdiv(lr_rho,csdiv) +ovp[1][0])

    opd.addVars( et, [pivlib.cpivvar(lr_temp,"T-NU",ovp[0][2]),
                      pivlib.cpivvar(lr_rho,"RHO-NU",ovp[1][2])] )

npd.save(mofpath)

if ( useipd ):
    opdpath = ifpath
else:
    opd.setTimes( pd.getTimes() )
    opdpath = "%s-THRML.ex2" % (ifpath[0:-4])
   
opd.save(opdpath)
