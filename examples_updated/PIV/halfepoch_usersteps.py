"""
Filename:  usersteps.py
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
  SPIVET user steps.

Contents:
  class coloradj
  class dltae0
  class vorticity
  class divergence
  class buoyancy
  class buoyflux
  class strainrate
  class viscstress
  class viscosity
  class spgravity
  class pvadj
  class eqtime
  class synchronize
  class prune
  class zdpvadj

  _splrep()
  _splev()

"""

from spivet import pivlib, flolib
from spivet import steps
from scipy import interpolate

from numpy import *


#################################################################
#
class coloradj(steps.spivetstep):
    """
    STEP: Adjust color

    Adjusts color of RGB images using the X-Rite ColorChecker 
    calibration.

    ---- carriage inputs ----
    imgchnls - (M) List of RGB image channels for each image (ie, a list 
               of rank 2).
               
    ---- carriage ouptuts ----
    imgchnls - Color corrected image channels.  Images will be adjusted 
               in place.

    ---- config dictionary contents ----
    adjmat --- (M) 3x3 color adjustment array from the calibration.
    
               Can be a URL to the pickled array.
    lbpath --- (O) Path to store COLORADJ array if necessary.  Only 
               needed if this object is retrieved.  Defaults to
               'COLORADJ'. 
    hsmthd --- (O) Hue separation method.  Defaults to 0.  Set to
               None to return RGB channels.  Otherwise, HSI channels
               will be returned.
    """
    def __init__(self,carriage={},config={}):
        steps.spivetstep.__init__(self,carriage,config)
        
    def execute(self):
        # Initialization.
        carriage = self.m_carriage
        imgchnls = carriage['imgchnls']
        
        adjmat = self.m_config['adjmat']
            
        if ( self.m_config.has_key('lbpath') ):
            lbpath = self.m_config['lbpath']
        else:
            lbpath = 'COLORADJ'
            
        if ( isinstance(adjmat,basestring) ):
            pth    = steps.url2local(adjmat,lbpath)
            adjmat = pivlib.pklload(pth[0])

            self.m_config['adjmat'] = adjmat
            
        if ( self.m_config.has_key('hsmthd') ):
            hsmthd = self.m_config['hsmthd']
        else:
            hsmthd = 0        

        # Correct the images and convert to HSI.
        for img in imgchnls:
            r    = img[0]
            g    = img[1]
            b    = img[2]
            argb = []
            for c in range(3):
                ac = adjmat[c,0]*r +adjmat[c,1]*g +adjmat[c,2]*b
                argb.append( ac.clip(0.,1.) )

            if ( hsmthd != None ):
                argb = pivlib.rgb2hsi(argb,hsmthd)

            r[:,...] = argb[0] 
            g[:,...] = argb[1]
            b[:,...] = argb[2]
            
            
#################################################################
#
class dltae0(steps.spivetstep):
    """
    STEP: Compute variable anomaly from Epoch 0.  All Epochs will
    be processed.

    ---- carriage inputs ----
    *pivdata - (M) Input PIVData object containing the variable.
               
    ---- carriage ouptuts ----
    *pivdata - Input PIVData object with DLTAE0-VAR variable added.
    
    ---- config dictionary contents ----
    varnm ---- (M) Name of variable or list of names to process.  
    pdname --- (O) Name of the PIVData object in the carriage to convert.
               Defaults to 'dpivdata'.
    """
    def __init__(self,carriage={},config={}):
        steps.spivetstep.__init__(self,carriage,config)
        
    def execute(self):
        # Initialization.
        carriage = self.m_carriage
        
        varnm = self.m_config['varnm']
        if ( not isinstance(varnm,list) ):
            varnm = [varnm]
        
        if ( self.m_config.has_key('pdname') ):
            pdname = self.m_config['pdname']
        else:
            pdname = 'dpivdata'
        
        pd = carriage[pdname]

        # Compute the delta.
        for vn in varnm:        
            dvn   = "DLTAE0-%s" % vn
            var0  = pd[0][vn]
            units = var0.units
            vtype = var0.vtype
            for e in range(0,len(pd)):
                var  = pd[e][vn]

                dvar = var -var0
                dvar.setAttr(dvn,units,vtype)
                
                pd[e].addVars(dvar)
                
                
#################################################################
#
class vorticity(steps.spivetstep):
    """
    STEP: Compute vorticity.  All Epochs will be processed.

    ---- carriage inputs ----
    *pivdata - (M) Input PIVData object containing the velocity 
               variable.
               
    ---- carriage ouptuts ----
    *pivdata - Input PIVData object with VORT variable added.
    
    ---- config dictionary contents ----
    velnm ---- (O) Name of variable containing the velocity field.
               Defaults to 'U'  
    pdname --- (O) Name of the PIVData object in the carriage to convert.
               Defaults to 'dpivdata'.
    """
    def __init__(self,carriage={},config={}):
        steps.spivetstep.__init__(self,carriage,config)
        
    def execute(self):
        # Initialization.
        carriage = self.m_carriage
        config   = self.m_config

        if ( config.has_key('velnm') ):
            velnm = config['velnm']
        else:
            velnm = 'U'
            
        if ( config.has_key('pdname') ):
            pdname = config['pdname']
        else:
            pdname = 'dpivdata'
        
        pd = carriage[pdname]
        
        # Compute vorticity.
        for e in pd:
            e.addVars( flolib.vorticity(e[velnm],pd.cellsz) )
            
            
#################################################################
#
class divergence(steps.spivetstep):
    """
    STEP: Compute divergence.  All Epochs will be processed.

    ---- carriage inputs ----
    *pivdata - (M) Input PIVData object containing the variable.
               
    ---- carriage ouptuts ----
    *pivdata - Input PIVData object with the divergence variable added.
               Variable name will be prepended with 'DEL.'.
    
    ---- config dictionary contents ----
    varnm ---- (O) Name of variable containing the vector field.
               Defaults to 'U'  
    pdname --- (O) Name of the PIVData object in the carriage to convert.
               Defaults to 'dpivdata'.
    """
    def __init__(self,carriage={},config={}):
        steps.spivetstep.__init__(self,carriage,config)
        
    def execute(self):
        # Initialization.
        carriage = self.m_carriage
        config   = self.m_config

        if ( config.has_key('varnm') ):
            varnm = config['varnm']
        else:
            varnm = 'U'
            
        if ( config.has_key('pdname') ):
            pdname = config['pdname']
        else:
            pdname = 'dpivdata'
        
        pd     = carriage[pdname]
        cellsz = pd.cellsz
        
        # Compute divergence.
        for e in pd:
            var = e[varnm]
            
            dwdz = flolib.d_di(var[0,...],0,cellsz[0])
            dvdy = flolib.d_di(var[1,...],1,cellsz[1])
            dudx = flolib.d_di(var[2,...],2,cellsz[2])
    
            div = dwdz +dvdy +dudx
            div.setAttr("DEL.%s" % varnm,"%s_MM" % (var.units))        
        
            e.addVars( div )
            
            
#################################################################
#
class buoyancy(steps.spivetstep):
    """
    STEP: Compute instantaneous buoyancy.  All Epochs will be 
    processed.
    
    Buoyancy here is defined as
        B = -drho g 
    where drho is the difference between density at the current
    Epoch and a reference density, and g is the gravitational 
    acceleration.

    The output buoyancy will have units of mN/mm^3.

    ---- carriage inputs ----
    *pivdata - (M) Input PIVData object containing the variable rho.
               
    ---- carriage ouptuts ----
    *pivdata - Input PIVData object with the buoyancy variable added.
               Variable name will be 'BUOYANCY'.
    
    ---- config dictionary contents ----
    varnm ---- (O) Name of variable containing the specific gravity or
               density.  Defaults to 'RHO-NU'  
    varsf ---- (O) Multiplicative scale factor to convert rho into 
               units of kg/mm^3.  Defaults to 1.
    gaccel --- (O) Gravitational acceleration constant.  Defaults to
               9800 mm/s^2.    
    pdname --- (O) Name of the PIVData object in the carriage to convert.
               Defaults to 'dpivdata'.
    rhoref --- (O) Reference density against which buoyancy is calculated.
               If specified, the reference density must have the same
               units as the variable provided by varnm.  If set to None, 
               the density field of of Epoch 0 will be used instead.  
               Defaults to None.
    """
    def __init__(self,carriage={},config={}):
        steps.spivetstep.__init__(self,carriage,config)
        
    def execute(self):
        # Initialization.
        carriage = self.m_carriage
        config   = self.m_config

        if ( config.has_key('varnm') ):
            varnm = config['varnm']
        else:
            varnm = 'RHO-NU'
        
        if ( config.has_key('varsf') ):
            varsf = config['varsf']
        else:
            varsf = 1.
            
        if ( config.has_key('gaccel') ):
            gaccel = config['gaccel']
        else:
            gaccel = 9800.
        
        if ( config.has_key('pdname') ):
            pdname = config['pdname']
        else:
            pdname = 'dpivdata'

        if ( config.has_key('rhoref') ):
            rhoref = config['rhoref']
        else:
            rhoref = None
            
        pd     = carriage[pdname]
        cellsz = pd.cellsz

        # Compute buoyancy.
        rho0 = pd[0][varnm].copy()
        if ( rhoref != None ):
            rho0[:,...] = rhoref
        for e in range(len(pd)):
            var = pd[e][varnm]
            
            buoy = gaccel*varsf*(rho0 -var)
            buoy.setAttr("BUOYANCY","mN_MM3",vtype=var.vtype)
            
            pd[e].addVars( buoy )            
            
            
#################################################################
#
class buoyflux(steps.spivetstep):
    """
    STEP: Compute buoyancy flux.  All Epochs will be processed.
    
    Buoyancy flux is defined as the amount of buoyancy transported 
    per unit area per unit time
        bflux = -drho g u
    where drho is the difference between density at the current
    Epoch and a reference density, and g is the gravitational 
    acceleration, and u is fluid velocity.  Note that the flux is a 
    vector and will be evaluated at cell centers.

    The output of buoyflux will have units of kPa/s.

    ---- carriage inputs ----
    *pivdata - (M) Input PIVData object containing the variables rho
               and u.
               
    ---- carriage ouptuts ----
    *pivdata - Input PIVData object with the buoyancy flux variable 
               added.  Variable name will be 'BUOYFLUX'.
    
    ---- config dictionary contents ----
    varnm ---- (O) List of variable names containing the density (or
               specific gravity) and velocity.  Defaults to 
               ['RHO-NU', 'U'].  
    varsf ---- (O) Multiplicative scale factor to convert rho into 
               units of kg/mm^3.  Defaults to 1.
    gaccel --- (O) Gravitational acceleration constant.  Defaults to
               9800 mm/s^2.    
    pdname --- (O) Name of the PIVData object in the carriage to convert.
               Defaults to 'dpivdata'.
    rhoref --- (O) Reference density against which buoyancy is calculated.
               If specified, the reference density must have the same
               units as the variable provide by varnm.  If set to None, 
               the density field of of Epoch 0 will be used instead.  
               Defaults to None.               
    """
    def __init__(self,carriage={},config={}):
        steps.spivetstep.__init__(self,carriage,config)
        
    def execute(self):
        # Initialization.
        carriage = self.m_carriage
        config   = self.m_config

        if ( config.has_key('varnm') ):
            varnm = config['varnm']
        else:
            varnm = ['RHO-NU','U']
        
        if ( config.has_key('varsf') ):
            varsf = config['varsf']
        else:
            varsf = 1.
        
        if ( config.has_key('gaccel') ):
            gaccel = config['gaccel']
        else:
            gaccel = 9800.
        
        if ( config.has_key('pdname') ):
            pdname = config['pdname']
        else:
            pdname = 'dpivdata'
            
        if ( config.has_key('rhoref') ):
            rhoref = config['rhoref']
        else:
            rhoref = None            
        
        pd     = carriage[pdname]
        cellsz = pd.cellsz
        
        # Compute buoyancy flux.
        rho0 = pd[0][varnm[0]].copy()
        if ( rhoref != None ):
            rho0[:,...] = rhoref

        bfshape    = array(rho0.shape)
        bfshape[0] = 3
        
        for e in range(len(pd)):
            rho = pd[e][varnm[0]]
            vel = pd[e][varnm[1]]
            
            bflux = varsf*gaccel*(rho0 -rho)*vel
            bflux.setAttr("BUOYFLUX","kPA_S",vtype=rho.vtype)
            
            pd[e].addVars( bflux )            
            

#################################################################
#
class strainrate(steps.spivetstep):
    """
    STEP: Compute strain rate tensor.  All Epochs will be processed.

    The strain rate tensor is given by
        s[l,m]= 1/2*( du_l/dx_m +du_m/dx_l )
        
    Recall that u and x will be processed in the order [z,y,x], so
    the full tensor will be
          / 2*dw/dz       dw/dy +dv/dz  dw/dx +du/dz \
      0.5*| dw/dy +dv/dz  2*dv/dy       dv/dx +du/dy |
          \ dw/dx +du/dz  dv/dx +du/dy  2*du/dx      /
    where w,v,u are the velocity components along the z,y,x-axes
    respectively.
    
    Derivatives will be computed using second-order central differences
    on all cells except those that line the edge of the domain.  For
    domain-edge cells, the derivatives will be computed using first
    order differences.
    
    The output of strainrate will have units of 1/s.

    ---- carriage inputs ----
    *pivdata - (M) Input PIVData object containing the velocity
               variable.
               
    ---- carriage ouptuts ----
    *pivdata - Input PIVData object with the strain rate variable 
               added.  Variable name will be 'STRNRATE'.
    
    ---- config dictionary contents ----
    varnm ---- (O) Variable name containing the velocity field.
               Defaults to 'U'.
    pdname --- (O) Name of the PIVData object in the carriage to convert.
               Defaults to 'pivdata'.
    """
    def __init__(self,carriage={},config={}):
        steps.spivetstep.__init__(self,carriage,config)
        
    def execute(self):
        # Initialization.
        carriage = self.m_carriage
        config   = self.m_config

        if ( config.has_key('varnm') ):
            varnm = config['varnm']
        else:
            varnm = 'U'
        
        if ( config.has_key('pdname') ):
            pdname = config['pdname']
        else:
            pdname = 'pivdata'
        
        pd     = carriage[pdname]
        cellsz = pd.cellsz
        ncells = pd[0].eshape
        
        # Compute strain rate tensor.
        for e in pd:
            vel  = e[varnm]

            # Get the gradient of the velocity field.
            der = []
            for vc in xrange(3):
                for d in xrange(3):
                    der.append( flolib.d_di(vel[vc,...],d,cellsz[d])  )
            
            # Build the tensor.
            s11 = der[0]                 # dw/dz
            s22 = der[4]                 # dv/dy
            s33 = der[8]                 # du/dx
            s12 = 0.5*(der[1] +der[3])   # 0.5*(dw/dy +dv/dz)
            s13 = 0.5*(der[2] +der[6])   # 0.5*(dw/dx +du/dz)
            s23 = 0.5*(der[5] +der[7])   # 0.5*(dv/dx +du/dy)
            
            stnsr = array([s11,s12,s13,s12,s22,s23,s13,s23,s33])
            stnsr = stnsr.reshape([9,ncells[0],ncells[1],ncells[2]])
            stnsr = pivlib.cpivvar(stnsr,"STRNRATE","1_S",vtype=vel.vtype)
        
            e.addVars(stnsr)
        

#################################################################
#
class viscstress(steps.spivetstep):
    """
    STEP: Compute viscous stress tensor.  All Epochs will be processed.

    The stress tensor is given by

        tau[l,m] = mu*(du_l/dx_m +du_m/dx_l) -2*k_lm*mu*dot(del,u)/3

    where mu is the dynamic viscosity, u is the velocity, k is the 
    Kronecker delta, dot() represents the dot product, and del is 
    the gradient operator.
    
    User can enforce dot(del,u) = 0 by setting the corresponding
    variable name to None in the config dictionary.  The viscous stress
    tensor will then be computed ignoring the divergence term.

    Recall that u and x will be processed in the order [z,y,x], so
    the full tensor will be
        / tau_zz tau_zy tau_zx \
        | tau_yz tau_yy tau_yx |
        \ tau_xz tau_xy tau_xx /
    
    The output of viscstress will have units of kPa.

    ---- carriage inputs ----
    *pivdata - (M) Input PIVData object containing the viscosity, 
               strain rate, and velocity divergence variables.
               
    ---- carriage ouptuts ----
    *pivdata - Input PIVData object with the viscous stress tensor 
               variable added.  Variable name will be 'VSTRESS'.
    
    ---- config dictionary contents ----
    varnm ---- (O) List containing the variable names for mu (viscosity),
               strain rate, and the divergence of the velocity field.
               Pass None for the third element of the list to set
               dot(del,u) = 0 (NOTE: Pass the Python object None, not the
               string 'None').  Defaults to ['MU','STRNRATE','DEL.U'].
    pdname --- (O) Name of the PIVData object in the carriage to convert.
               Defaults to 'pivdata'.
    """
    def __init__(self,carriage={},config={}):
        steps.spivetstep.__init__(self,carriage,config)
        
    def execute(self):
        # Initialization.
        carriage = self.m_carriage
        config   = self.m_config

        if ( config.has_key('varnm') ):
            varnm = config['varnm']
        else:
            varnm = ['MU','STRNRATE','DEL.U']
        
        if ( config.has_key('pdname') ):
            pdname = config['pdname']
        else:
            pdname = 'pivdata'
        
        pd     = carriage[pdname]
        cellsz = pd.cellsz
        ncells = pd[0].eshape
        
        # Compute stress tensor.
        for e in pd:
            mu = e[varnm[0]]
            sr = e[varnm[1]]

            stnsr = 2.*mu*sr

            if ( varnm[2] != None ):
                dv  = e[varnm[2]]
                ddu = -2./3.*mu*dv            
                for c in xrange(3):
                    stnsr[c*4,...] = stnsr[c*4,...] +ddu
           
            stnsr.setAttr("VSTRESS","kPA")
                    
            e.addVars(stnsr)        
        

#################################################################
#
class viscosity(steps.spivetstep):
    """
    STEP: Compute viscosity.  All Epochs will be processed.

    ---- carriage inputs ----
    *pivdata - (M) Input PIVData object containing the temperature
               field variable.
               
    ---- carriage ouptuts ----
    *pivdata - Input PIVData object with MU variable added.
    
    ---- config dictionary contents ----
    cal ------ (M) Viscosity calibration coefficients for function of the
               form:
                   visc = a0*exp(a1 +a2*T +a3*T^2)
               The a* are the coefficients.  These coefficients should be
               such that the viscosity has units of kPa-s.
    tmpnm ---- (O) Name of variable containing the temperature field.
               Defaults to 'T'  
    pdname --- (O) Name of the PIVData object in the carriage to convert.
               Defaults to 'dpivdata'.
    """
    def __init__(self,carriage={},config={}):
        steps.spivetstep.__init__(self,carriage,config)
        
    def execute(self):
        # Initialization.
        carriage = self.m_carriage
        config   = self.m_config

        cal = config['cal']
        print "Viscosity coefficients: %g, %g, %g, %g" % \
            (cal[0],cal[1],cal[2],cal[3])
        
        if ( config.has_key('tmpnm') ):
            tmpnm = config['tmpnm']
        else:
            tmpnm = 'T'
            
        if ( config.has_key('pdname') ):
            pdname = config['pdname']
        else:
            pdname = 'dpivdata'
        
        pd     = carriage[pdname]
        cellsz = pd.cellsz
                
        # Compute viscosity.
        for e in pd:
            tmp = e[tmpnm]
            
            visc = cal[0]*exp(cal[1] +cal[2]*tmp +cal[3]*tmp**2)
            visc.setAttr('MU','kPA-S')
            
            e.addVars(visc)
            
            
#################################################################
#
class spgravity(steps.spivetstep):
    """
    STEP: Compute specific gravity.  All Epochs will be processed.

    ---- carriage inputs ----
    *pivdata - (M) Input PIVData object containing the temperature
               field variable.
               
    ---- carriage ouptuts ----
    *pivdata - Input PIVData object with SPGRAVITY variable added.
    
    ---- config dictionary contents ----
    cal ------ (M) Specific gravity calibration coefficients for function 
               of the form:
                   spgravity = a0*T +a1
               The a* are the coefficients.
    tmpnm ---- (O) Name of variable containing the temperature field.
               Defaults to 'T'  
    pdname --- (O) Name of the PIVData object in the carriage to convert.
               Defaults to 'dpivdata'.
    """
    def __init__(self,carriage={},config={}):
        steps.spivetstep.__init__(self,carriage,config)
        
    def execute(self):
        # Initialization.
        carriage = self.m_carriage
        config   = self.m_config

        cal = config['cal']
        print "Specific gravity coefficients: %g, %g" % \
            (cal[0],cal[1])
        
        if ( config.has_key('tmpnm') ):
            tmpnm = config['tmpnm']
        else:
            tmpnm = 'T'
            
        if ( config.has_key('pdname') ):
            pdname = config['pdname']
        else:
            pdname = 'dpivdata'
        
        pd     = carriage[pdname]
        cellsz = pd.cellsz
                
        # Compute specific gravity.
        for e in pd:
            tmp = e[tmpnm]
            
            spg = cal[0]*tmp +cal[1]
            spg.setAttr('SPGRAVITY','NA')
            
            e.addVars(spg)
            
            
#################################################################
#
class pvadj(steps.spivetstep):
    """
    STEP: Adjusts a variable by subtracting a constant from the cells
    of each plane.  The primary motivation for this step is to adjust 
    the velocity to compensate for known linear slide variation.  All 
    Epochs in the PIVData object will be processed.

    If pvadj is used to correct velocity values for linear slide
    positional errors, then pvadja should be an array of displacement
    errors and tnrm set to True.  

    ---- carriage inputs ----
    *pivdata - (M) Input PIVData object containing the variable.
               NOTE: If tnrm = True, then the PIVData object must
               also contain the PLNR-DLTATIME variable created with
               the recordtime step.
               
    ---- carriage ouptuts ----
    *pivdata - Input PIVData object with adjusted variable appended.
               Adjusted variable name will be appended with '-PA'.
    
    ---- config dictionary contents ----
    varnm ---- (M) Name of variable.
    varcmp --- (M) Integer index into the PIVVar specifying which variable 
               component (eg, 0 for the z-component) to adjust.  
    pvadja --- (M) 1D array of planar 'errors.'  There should be one entry 
               for each plane in the dataset.  The pvadja value for
               a given plane will be subtracted from the varcmp-component 
               of all cells in the plane.  If tnrm = True, then the
               pvadja values will be normalized by the corresponding 
               PLNR-DLTATIME value for the plane prior to subtracting
               from the varcmp-component. 
               
               pvadja can also be a URL to a file created using 
               pivlib.pkldump().
    tnrm ----- (O) Boolean flag specifying whether the values of pvadja
               should be normalized by the PLNR-DLTATIME values prior
               to subtracting from the varcmp-component.  Defaults to 
               False.
    pdname --- (O) Name of the PIVData object in the carriage to convert.
               Defaults to 'dpivdata'.
    lbpath --- (O) Path to store pvadja array if necessary.  Only
               needed if pvadja objects are retrieved.  Defaults to
               'CALIBRATION'.               
    """
    def __init__(self,carriage={},config={}):
        steps.spivetstep.__init__(self,carriage,config)
        
    def execute(self):
        # Initialization.
        carriage = self.m_carriage
        
        varnm  = self.m_config['varnm']
        varcmp = self.m_config['varcmp']

        if ( self.m_config.has_key('lbpath') ):
            lbpath = self.m_config['lbpath']
        else:
            lbpath = 'CALIBRATION'
        
        pvadja = self.m_config['pvadja']
        if ( isinstance(pvadja,basestring) ):
            pth    = steps.url2local(pvadja,lbpath)
            pvadja = pivlib.pklload(pth[0])
        
        if ( self.m_config.has_key('tnrm') ):
            tnrm = self.m_config['tnrm']
        else:
            tnrm = False
        
        if ( self.m_config.has_key('pdname') ):
            pdname = self.m_config['pdname']
        else:
            pdname = 'dpivdata'
        
        pd = carriage[pdname]

        # Compute the delta.
        dvn   = "%s-PA" % varnm
        units = pd[0][varnm].units
        vtype = pd[0][varnm].vtype        
        for e in pd:
            var  = e[varnm]
            avar = var.copy()

            if ( tnrm ):
                pdtv    = e['PLNR-DLTATIME']
                npvadja = pvadja/pdtv[0,:,0,0]
            else:
                npvadja = pvadja
            
            for s in range(avar.shape[1]):
                 avar[varcmp,s,...] = avar[varcmp,s,...] -npvadja[s]

            avar.setAttr(dvn,units,vtype)
            e.addVars(avar)

            
#################################################################
#
def _splrep(x,y,s):
    """ 
    ----
    
    x              # splrep() x argument.
    y              # splrep() y argument.
    s              # splrep() s argument (the smoothness). 
    
    ----
    
    Helper function for synchronize step.  _splrep() is a wrapper
    around Scipy's splrep() function that the synchronize step can
    call using map().
    
    Returns spl, the spline object. 
    """
    spl = interpolate.splrep(x,y,s=s)
    return spl


#################################################################
#
def _splev(x,spl):
    """ 
    ----
    
    x              # splev() x argument.
    spl            # splev() tck argument (the spline object).
    
    ----
    
    Helper function for synchronize step.  _splev() is a wrapper
    around Scipy's splev() function that the synchronize step can
    call using map().
    
    Returns val, the interpolated value. 
    """
    val = interpolate.splev(x,spl)
    return val


#################################################################
#
class eqtime(steps.spivetstep):
    """
    STEP: Equipartitions total elapsed time between Epochs.  Each 
    variable in varnm is temporally interpolated using cubic splines 
    such that the Epochs within the PIVDATA object are equally spaced 
    in time.  The primary motivation for this step is to correct 
    for jitter in elapsed time between Epochs. 

    The eqtime step should only be run after all Epochs are
    available in the PIVData object since data will be splined across
    those Epochs.

    ---- carriage inputs ----
    *pivdata - (M) Input PIVData object containing the variable(s) to
               be temporally adjusted.  
                                 
    ---- carriage ouptuts ----
    *pivdata - A new PIVData object containing only the adjusted 
               variables.  Adjusted variable names will be the
               same as input.  The new PIVData object on the carriage
               will be named 
                   "%s-eq" % pdname
    
    ---- config dictionary contents ----
    varnm ---- (M) Name of variable or list of variables to adjust.    
    pdname --- (O) Name of the PIVData object in the carriage containing
               the variables.  Defaults to 'dpivdata'.
    """
    def __init__(self,carriage={},config={}):
        steps.spivetstep.__init__(self,carriage,config)
        
    def execute(self):
        # Initialization.
        carriage = self.m_carriage
        
        varnm = self.m_config['varnm']
        if ( not isinstance(varnm,list) ):
            varnm = [ varnm ]

        if ( self.m_config.has_key('pdname') ):
            pdname = self.m_config['pdname']
        else:
            pdname = 'dpivdata'
        
        pd   = carriage[pdname]        
        nepc = len(pd)
        if ( nepc < 4 ):
            raise ValueError("PIVData must contain at least 4 Epochs")

        # Setup the new PIVData object.
        npdname = "%s-eq" % pdname
        npd     = pivlib.PIVData(pd.cellsz,pd.origin)

        # Get Epoch times and equipartition.
        times  = array(pd.getTimes())
        dt     = (times[-1] -times[0])/(nepc -1)
        qtimes = dt*arange(nepc,dtype=float) +times[0]
        
        ncells  = array(pd[0].eshape)
        tncells = ncells.prod()

        atimes = times.reshape([1,nepc])
        atimes = atimes.repeat(tncells,0)
        
        aqtimes = qtimes.reshape([1,nepc])
        aqtimes = aqtimes.repeat(tncells,0)
        
        # Create the smoothing vector.  Set svec to zero for all cells
        # so that splrep does interpolation.
        svec = zeros(tncells)

        # Adjust.
        for vn in varnm:
            var = []
            for e in pd:
                var.append(e[vn])
            var  = array(var)
            avar = empty(var.shape,dtype=float)
            
            # var has shape [nepochs,ncmp,nz,ny,nx]
            ncmp = var.shape[1]
            for c in xrange(ncmp):
                cmp = var[:,c,:,...]
                cmp = cmp.reshape((nepc,tncells)).transpose()
                
                spl = map(_splrep,atimes,cmp,svec)
                for e in xrange(nepc):
                    acmp = array( map(_splev,aqtimes[:,e],spl) )
                    avar[e,c,...] = acmp.reshape((ncells[0],ncells[1],ncells[2]))
            
            for e in xrange(nepc):
                pivvar = pivlib.PIVVar((ncmp,ncells[0],ncells[1],ncells[2]),
                                       vn,
                                       pd[e][vn].units,
                                       dtype=float,vtype=pd[e][vn].vtype,
                                       data=avar[e,...])
                npd.addVars(e,pivvar)
                
        # Store the results.
        npd.setTimes(qtimes)
        carriage[npdname] = npd


#################################################################
#
class synchronize(steps.spivetstep):
    """
    STEP: Interpolates a variable or list of variables in time using 
    splines and adjusts the variables as though the planes were all
    taken simultaneously.  This step is primarily aimed at 'degrouping'
    datasets where images from some planes are taken in groups.

    The variables will be synchronized by evolving planar data forward
    in time.  That is, planes acquired earlier will be evolved as though
    they were taken later.  This approach is meant to leverage existing
    data instead of extrapolating backward in time at the start of the
    experiment.  As a result, however, the final PIVData object will
    contain one less Epoch (the last Epoch will be discared).
    
    Since we interpolate to 0.5 epoch here, all epoch will be kept in final result.

    The synchronize step should only be run after all Epochs are
    available.

    ---- carriage inputs ----
    *pivdata - (M) Input PIVData object containing the variable(s) to
               be synchronized.  PIVData object must contain the PLNR-TIME
               variable (created with recordtime step).  
               
               NOTE: All variables being synchronized must have the same 
               vtype as PLNR-TIME.  No checks are performed.
               
    ---- carriage ouptuts ----
    *pivdata - A new PIVData object with synchronized variables appended.
               Synchronized variable names will be appended with '-SN'.
               The new PIVData object will replace the old pivdata
               object named pdname on the carriage.
    
    ---- config dictionary contents ----
    varnm ---- (M) Name of variable or list of variables to synchronize.    
    pdname --- (O) Name of the PIVData object in the carriage containing
               the variables.  Defaults to 'dpivdata'.
    """
    def __init__(self,carriage={},config={}):
        steps.spivetstep.__init__(self,carriage,config)
        
    def execute(self):
        # Initialization.
        carriage = self.m_carriage
        
        varnm = self.m_config['varnm']
        if ( not isinstance(varnm,list) ):
            varnm = [ varnm ]

        if ( self.m_config.has_key('pdname') ):
            pdname = self.m_config['pdname']
        else:
            pdname = 'dpivdata'
        
        pd   = carriage[pdname]        
        nepc = len(pd)
        if ( nepc < 4 ):
            raise ValueError("PIVData must contain at least 4 Epochs")

        # Extract the PLNR-TIME variable and adjust it.
        plnrt = []
        for e in pd:
            plnrt.append(e['PLNR-TIME'])
        plnrt = array(plnrt)
        
        ncells   = plnrt.shape[2::]
        tnpcells = ncells[1]*ncells[2]   # Planar cells.
        tncells  = ncells[0]*tnpcells    # All cells.
        
        aplnrt = empty(plnrt.shape,dtype=float)
        for e in range(nepc):
            aplnrt[e,...] = (plnrt[e,...].min()*0.5 + plnrt[e,...].max()*0.5)

        plnrt  = plnrt.reshape((nepc,tncells)).transpose()
        aplnrt = aplnrt.reshape((nepc,tncells)).transpose()

        # Create the smoothing vector.  Set svec to zero for all cells
        # so that splrep does interpolation.
        svec = zeros(tncells)

        # Synchronize.
        for vn in varnm:
            var = []
            for e in pd:
                var.append(e[vn])
            var  = array(var)
            avar = empty(var.shape,dtype=float)
            
            # var has shape [nepochs,ncmp,nz,ny,nx]
            ncmp = var.shape[1]
            for c in range(ncmp):
                cmp = var[:,c,:,...]
                cmp = cmp.reshape((nepc,tncells)).transpose()
                
                spl = map(_splrep,plnrt,cmp,svec)
                for e in range(nepc):
                    acmp = array( map(_splev,aplnrt[:,e],spl) )
                    avar[e,c,...] = acmp.reshape((ncells[0],ncells[1],ncells[2]))
            
            for e in range(nepc):
                pivvar = pivlib.PIVVar((ncmp,ncells[0],ncells[1],ncells[2]),
                                       "%s-SN" % vn,
                                       pd[e][vn].units,
                                       dtype=float,vtype=pd[e][vn].vtype,
                                       data=avar[e,...])
                pd[e].addVars(pivvar)
                
        # Store the results.
        pd.setTimes(aplnrt[0])
        pd = pd[0:(nepc)]
        carriage[pdname] = pd
        
        
#################################################################
#
class prune(steps.spivetstep):
    """
    STEP: Prunes a variable, or list of variables, from a PIVData
    object.

    ---- carriage inputs ----
    *pivdata - (M) PIVData object for variable pruning.

    ---- carriage outputs ----
    *pivdata - The input PIVData object with variables pruned.

    ---- config dictionary contents ----
    varnm ---- (M) Name of variable or list of names to delete from the
               PIVData object.  
    pdname --- (O) Name of the PIVData object in the carriage to convert.
               Eg, 'epivdata'.  Defaults to 'pivdata'.
    """
    def __init__(self,carriage={},config={}):
        steps.spivetstep.__init__(self,carriage,config)

    def execute(self):
        # Initialization.
        carriage = self.m_carriage
        if ( self.m_config.has_key('pdname') ):
            pdname = self.m_config['pdname']
        else:
            pdname = 'pivdata'

        pd = carriage[pdname]

        varnm = self.m_config['varnm']
        if ( not isinstance(varnm,list) ):
            varnm = [ varnm ]
        
        # Prune.
        for e in range( len(pd) ):
            for vn in varnm:
                pd[e].pop(vn)

                
#################################################################
#                
class zdpvadj(steps.spivetstep):
    """
    STEP:  Like pvadj, the zdpvadj step adjusts a variable by 
    subtracting a constant from the cells of each plane.  And also
    like pvadj, zdpvadj is motivated by the need to compensate for
    linear slide positional variability.
    
    The pvadj step applies a 'known' correction to the computed 
    velocity field.  This known correction is determined by
    imaging a quiescent tank several times and forming a statistical
    estimate of the linear slide positional error.  The resulting
    correction should compensate for repeatable position errors that 
    arise from lash, manufacturing concentricity variation, etc.  
    pvadj does not, however, compensate for positional errors that 
    arise due to belt tension fluctuations resulting from ambient
    temperature variation, etc.  Furthermore, an analysis of the 
    statistical distribution of the positional errors used to 
    compute the pvadj adjustment values shows a large variance.
    Hence, using only pvadj to correct the velocity field will likely
    not remove all of the error.
    
    zdpvadj works on a completely different principle.  zdpvadj assumes
    a constant density fluid and computes the velocity flux through
    a control volume.  The control volume is aligned with the rectilinear
    coordinates of the dataset, with one of the k-faces (the normal to a
    CV k-face is parallel to the z-axis) fixed at 
        k = kmaxindex -1/2 
    where kmaxindex is the largest cell index along the z-axis.  The 
    velocity at cell k = kmaxindex is known from the no-slip condition.  
    The other k-face, initially positioned at k = 1/2 where the no-slip 
    condition also holds, is then iteratively walked toward the 
    k = kmaxindex -1/2 k-face (which is always fixed).  Given that the 
    in-plane displacements (ie, x,y-components of displacement) are more 
    accurate than the out of plane (the degree to which this is true 
    depends on several factors, like the angle between the cameras), any 
    residual flux through this control volume must be primarily due to an 
    erroneous z-velocity on the mobile k-face.  zdpvadj assumes that the 
    residual flux over the CV is due to a plane-constant error in the 
    z-velocity at the mobile face.  A z-velocity correction is then 
    formulated and applied to k-face cells just inside the CV.  The CV 
    face is subsequently moved to the next k-plane, and process is 
    repeated. 
 
    The control volume is a rectangle with 6 'large' faces.  Each 
    perimeter cell of the CV contributes one subface to the larger CV 
    face.
    
    All Epochs in the PIVData object will be processed.
    
    ---- carriage inputs ----
    *pivdata - (M) Input PIVData object containing the variable to be
               adjusted.
               
    ---- carriage ouptuts ----
    *pivdata - Input PIVData object with adjusted variable appended.
               Adjusted variable name will be appended with '-ZD'.
    
    ---- config dictionary contents ----
    varnm ---- (M) Name of variable.  The z-component (ie, component 0)
               of this variable will be adjusted.
    cvbndx --- (M) 2x2 array specifying the [[jmn,jmx],[imn,imx]] bounds
               of the control volume along the y and x axes respectively.
               In standard Python slice convention, jmn and imn will be
               the indices of the first cell inside the control volume,
               while jmx and imx will be the indices of the cells just
               outsize the CV.
    kfhcw ---- (O) CV k-face halo cell average z-velocity values over
               the k-face.  'Halo' cells are the perimeter of cells just 
               outside the CV.  The k=0 and k=ncells[0] -1 halo cells (ie,
               at the rear and front of the tank) are the cells that are
               assumed to have a known average z-velocity, w (w = 0.0 if 
               the no-slip condition holds).  kfhcw is a two element array 
               specifying this known mean z-velocity at k=[0,ncells[0]-1].  
               
               If kfhcw is not specified or set to None, kfhcw will
               be estimated using the following procedure.  A CV enclosing
               all cells in z, but otherwise using cvbndx, will be
               constructed.  The u,v volumetric flux (ie, through the i,j
               faces) will be computed.  By continuity for constant
               density fluids, any residual volumetric flux remaining 
               must flow through the k-faces. The zdpvadj step will divide
               this expected k-face flux evenly between the two k-faces.
               Then for each k-face, a velocity correction will be applied
               based on the difference between the expected k-face flux
               and the k-face flux computed from the unadjusted variable.
    rlxf ----- (O) Relaxation factor for z-velocity corrections.  The
               formulation as described so far assumes that all velocity
               components are known exactly.  In reality, the computed
               velocities contain noise from experimental and computational
               error.  As a result, applying the full velocity correction 
               estimated using the above methods can cause the technique to
               oscillate (one k-plane is under-corrected, the next over-
               corrected, etc).  rlxf is used to scale the correction
               and damp these oscillations.  That is
                   azvcor = rlxf*zvcor
               where azvcor is the applied correction, and zcvor is
               the raw correction using the CV flux.  Defaults to 0.6.
    nits ----- (O) Number of iterations to run the adjustment routine.
               Given that rlxf defaults to a number less than 1.0, the
               z-velocity could be undercorrected if only one iteration
               of the technique is run.  nits allows the full process
               described above to be run iteratively.  Defaults to 3.
    pdname --- (O) Name of the PIVData object in the carriage to convert.
               Defaults to 'dpivdata'.
    """
    def __init__(self,carriage={},config={}):
        steps.spivetstep.__init__(self,carriage,config)

    def execute(self):    
        print "STARTING: zdpvadj"
        # Initialization.
        carriage = self.m_carriage
        
        varnm  = self.m_config['varnm']
        cvbndx = self.m_config['cvbndx']
        cvbndx = array(cvbndx)

        if ( self.m_config.has_key('kfhcw') ):
            kfhcw  = self.m_config['kfhcw']
            
            if ( kfhcw == None ):
                kfhcw  = [0.,0.]
                hkfhcw = False
            else:
                hkfhcw = True
        else:
            kfhcw  = [0.,0.]
            hkfhcw = False
            
        if ( self.m_config.has_key('rlxf') ):
            rlxf = self.m_config['rlxf']
        else:
            rlxf = 0.6
            
        if ( self.m_config.has_key('nits') ):
            nits = self.m_config['nits']
        else:
            nits = 3
        
        if ( self.m_config.has_key('pdname') ):
            pdname = self.m_config['pdname']
        else:
            pdname = 'dpivdata'
        
        pd = carriage[pdname]
                
        cellsz = pd.cellsz
        ncells = pd[0][varnm].shape[1::]
    
        scellsz = sign(cellsz)
    
        # Compute subface areas.
        ifa = abs( cellsz[0]*cellsz[1] )
        jfa = abs( cellsz[0]*cellsz[2] )
        kfa = abs( cellsz[1]*cellsz[2] )
    
        # In what follows, the CV is taken as the perimeter faces that enclose
        # the cells of cvbndx.  So if cell[k,j,i] is a corner cell just inside the
        # CV, then centers for three subfaces of the CV lie at [k-1/2,j,i], 
        # [k,j-1/2,i], and [k,j,i-1/2].  
        #
        # The k index is oriented with z, j with y, and i with x.
    
        # Python slice boundaries for CV-enclosed cells.
        kmn = 1
        kmx = ncells[0] -1
        jmn = cvbndx[0,0]
        jmx = cvbndx[0,1]
        imn = cvbndx[1,0]
        imx = cvbndx[1,1]
    
        # Indices for halo cells.  These halo cells are the cells that 
        # surround the CV (ie, outside of CV).  If the CV has the same 
        # shape as var, then set the halo cell value equal to the 
        # perimeter cell values of var (ie, use edge extension).
        jmnn = max(jmn -1,0)
        jmxn = min(jmx,ncells[1]-1)
        imnn = max(imn -1,0)
        imxn = min(imx,ncells[2]-1)
    
        # Compute number of cells and area for k-face of CV.
        ncvkfcells = (jmx -jmn)*(imx -imn)  # Number of CV cells in one CV k-face.
        cvkfa      = ncvkfcells*kfa         # Area of CV k-face.
    
        # Adjust the variable.
        cnt = 0
        for e in pd:
            print "Epoch %i" % cnt
            var  = e[varnm]
            avar = var.copy()
            avar.setAttr("%s-ZD" % var.name,var.units,var.vtype)
    
            for it in range(nits):    
                # Build velocity values at subface centers.  
                ujm  = 0.5*jfa*( avar[1, :, jmnn,    imn:imx] 
                                +avar[1, :, jmn,     imn:imx] )  # -u.da on jmn-1/2
                ujp  = 0.5*jfa*( avar[1, :, jmx-1,   imn:imx] 
                                +avar[1, :, jmxn,    imn:imx] )  # u.da on jmx-1/2
                uim  = 0.5*ifa*( avar[2, :, jmn:jmx, imnn] 
                                +avar[2, :, jmn:jmx, imn] )      # -u.da on imn-1/2
                uip  = 0.5*ifa*( avar[2, :, jmn:jmx, imx-1] 
                                +avar[2, :, jmn:jmx, imxn] )     # u.da on imx-1/2
        
                if ( not hkfhcw ):
                    # Compute the total flux though the j,i faces.  
                    # Construct a temporary CV that includes all z-cells
                    # and use edge extension.  Any residual flux in or out
                    # of the CV will be evenly split between the two
                    # k-faces.
                    flux = scellsz[1]*( ujp.sum() 
                                       -ujm.sum() ) \
                          +scellsz[2]*( uip.sum() 
                                       -uim.sum() )
                  
                    kfhcw[0] = -scellsz[0]*0.5*flux/cvkfa
                    kfhcw[1] = -kfhcw[0]
                      
                # u.da on kmx-1/2
                ukp  = 0.5*kfa*( avar[0, kmx-1, jmn:jmx, imn:imx] +kfhcw[1] ) 
                ukps = ukp.sum()
        
                # Adjust the halo cells for the k-axis. 
                epaw          = avar[0,0,jmn:jmx,imn:imx].mean()
                avar[0,0,...] = avar[0,0,...] -rlxf*(epaw -kfhcw[0])
        
                print " | w_err it %i, Plane %i: %f" % (it,0,epaw -kfhcw[0])
        
                # Compute the erroneous z-component on the kmn-1/2 CV face.
                for k in range(kmn,kmx):
                    # -u.da on kmn-1/2.
                    ukm  = 0.5*kfa*( avar[0, k-1, jmn:jmx, imn:imx]  \
                                    +avar[0, k,   jmn:jmx, imn:imx] ) 
                    
                    # flux = integral( u.da ) over CV.  scellsz takes care
                    # of face orientation with respect to world coordinates.
                    # For our particular setup, the kmx-1/2 face is in 
                    # negative z.
                    flux = scellsz[0]*( ukps -ukm.sum() ) \
                          +scellsz[1]*( ujp[k:kmx,...].sum() 
                                       -ujm[k:kmx,...].sum() ) \
                          +scellsz[2]*( uip[k:kmx,...].sum() 
                                       -uim[k:kmx,...].sum() )
        
                    # Allocate the flux to the kmn - 1/2 face and correct avar.
                    # If the velocities of the cells on the k and k-1 planes were
                    # known exactly, then the 'correct' w flux across the CV kmn -1/2
                    # face would be
                    #
                    #   fc = -scellsz[0]*kfa*(w[k,...].sum() +w[k-1,...].sum() )/2.
                    #
                    # where w is the z-component of velocity.  But the observed flux 
                    # across the CV face is contaminated by an erroneous w-component 
                    # caused by slide positional error. 
                    #
                    #   fo = -scellsz[0]*kfa*( w[k,...].sum() 
                    #                         +ncvkfcells*w_err 
                    #                         +w[k-1,...].sum() )/2.
                    #
                    # If we assume that the k = kmn - 1 and k = kmx velocities (ie,
                    # the kfhcw) are well-known and 'error-free', then the 
                    # w-velocities for the cells of each plane can be updated 
                    # incrementally (so the k-1 plane will always be corrected).  
                    # This arrangement leaves the full w_err to the k plane.  Note 
                    # that flux = fo -fc.
                    w_err         = -2.*scellsz[0]*flux/cvkfa         
                    avar[0,k,...] = avar[0,k,...] -rlxf*w_err
        
                    print " | w_err it %i, Plane %i: %f" % (it,k,w_err)
        
                # Adjust the halo cells for the k-axis.
                epaw            = avar[0,kmx,jmn:jmx,imn:imx].mean()
                avar[0,kmx,...] = avar[0,kmx,...] -rlxf*(epaw -kfhcw[1])
        
                print " | w_err it %i, Plane %i: %f" % (it,kmx,epaw -kfhcw[1])
                print " | "
    
            e.addVars(avar)
            cnt = cnt +1   
        
        print " | EXITING: zdpvadj."


#################################################################
#
class pressure(steps.spivetstep):
    """
    STEP: Compute dynamic pressure for thermal convection.  All 
    Epochs will be processed.

    The momentum equation for Stokes flow is taken to be

        del(p') = dot(del,tau) +rho'*f

    where P' is the dynamic pressure, tau is the viscous stress
    tensor, rho' is variation in density, and f is the body force
    vector.  In integral form, this becomes

        integral(p'*da)_S = integral( (dot(del,tau) +rho'*f)dv )_V

    where S and V are the cell surface and volume, respectively. 
    And discretizing provides
    
        p'[l+1] = V*( dot(del,dot(tau,c)) +rho'*dot(f,c) )[l]/A[l+1] 
                 +p'[l]*A[l]/A[l+1]

    where l indicates a cell face index parallel to the c direction
    (with c varying from 0 to 2), and A[l] is the area of face l.
    
    Using the discrete form above, pressure is computed on cell faces
    by marching away from the (i,j,k) = (0,0,0) cell, where the
    user has specified a reference pressure value.  Neumann 
    boundary conditions are used for the three sides of the outer 
    perimeter of the flow domain of which the outer faces of cell 
    (i,j,k) = (0,0,0) are a part.  So if p'[1/2] represents the
    pressure at the cell center for cell (0,0,0), then p'[0] on the
    outer face will be p'[0] = p'[1/2] (ie, the gradient is zero).
    This approach is illustrated by the 2D graphic below.  Cell
    (0,0,0) in this six-cell domain is denoted by a 0, and the outer 
    perimeter of the flow domain subject to the Neumann BC is marked 
    with a '>' or 'v'.

           vvvvvvvvvvv    
         > |---|---|---|
         > | 0 |   |   |  
         > |---|---|---|
         > |   |   |   |
         > |---|---|---|
    
    
    After computing pressure values at all 6 cell faces, the 
    pressure at the cell center will be set to the face average.  

    NOTE: A valid density field and viscous stress tensor must be 
    available.  The user does not need to compute rho' as it will
    be computed as
        rho' = rho[current Epoch] -rho[Epoch 0]
        
    Units of density and viscous stresses should be kg/mm^3 and kPa, 
    respectively.  
    
    Output pressure will have units of kPa.

    ---- carriage inputs ----
    *pivdata - (M) Input PIVData object containing the density and
               viscous stress tensor.
               
    ---- carriage ouptuts ----
    *pivdata - Input PIVData object with pressure variable added.
               Variable name will be 'P'.
    
    ---- config dictionary contents ----
    varnm ---- (O) List of variable names for density and viscous
               stress tensor.  Defaults to ['RHO-NU','VSTRESS'].
    rhosf ---- (O) Multiplicative scale factor to convert rho into 
               units of kg/mm^3.  Defaults to 1.               
    c0p ------ (O) Pressure for cell (i,j,k) = (0,0,0).  Should be in
               units of kPa.  Defaults to 0.
    bforce --- (O) Body force vector [z,y,x] in units of mm/s^2.  
               Specifying gravity in laboratory coordinates would be 
               [0,9800,0].  Defaults to [0,0,0]. 
    pdname --- (O) Name of the PIVData object in the carriage to convert.
               Defaults to 'dpivdata'.
    """
    def __init__(self,carriage={},config={}):
        steps.spivetstep.__init__(self,carriage,config)
        
    def execute(self):
        # Initialization.
        carriage = self.m_carriage
        config   = self.m_config

        if ( config.has_key('varnm') ):
            varnm = config['varnm']
        else:
            varnm = ['RHO-NU','VSTRESS']
            
        if ( config.has_key('rhosf') ):
            rhosf = config['rhosf']
        else:
            rhosf = 1.
            
        if ( config.has_key('c0p') ):
            c0p = config['c0p']
        else:
            c0p = 0.
            
        if ( config.has_key('bforce') ):
            bforce = array( config['bforce'] )
        else:
            bforce = array([0.,0.,0.])
            
        if ( config.has_key('pdname') ):
            pdname = config['pdname']
        else:
            pdname = 'dpivdata'
        
        pd = carriage[pdname]
        
        nrm = 1./6.
        
        # Allocate storage for the face pressures.
        ncells = array( pd[0].eshape )
        kfp    = empty([ncells[0]+1,ncells[1],  ncells[2]  ],dtype=float)
        jfp    = empty([ncells[0],  ncells[1]+1,ncells[2]  ],dtype=float)
        ifp    = empty([ncells[0],  ncells[1],  ncells[2]+1],dtype=float)

        # Allocate storage for V*( dot(del,dot(tau,c)) +rho'*dot(f,c) ).
        rhs = empty([3,ncells[0],ncells[1],ncells[2]],dtype=float)

        # Compute face areas and cell volume.
        cellsz = pd.cellsz
        kfa    = sign(cellsz[0])*abs(cellsz[1]*cellsz[2])
        jfa    = sign(cellsz[1])*abs(cellsz[0]*cellsz[2])
        ifa    = sign(cellsz[2])*abs(cellsz[0]*cellsz[1]) 
               
        vol = abs( cellsz[0]*cellsz[1]*cellsz[2] )
               
        # Get rho0. 
        rho0 = pd[0][varnm[0]]
                
        for e in pd:
            rho = rhosf*(e[varnm[0]] -rho0)
            tau = e[varnm[1]]
            
            pres = zeros([ncells[0],ncells[1],ncells[2]],dtype=float)
            
            # Compute dot(del,tau).
            ddt = []
            for i in xrange(3):
                for j in xrange(3):
                    ddt.append( flolib.d_di(tau[i*3+j,...],i,cellsz[i]) )
            ddt = [ ddt[0] +ddt[3] +ddt[6],
                    ddt[1] +ddt[4] +ddt[7],
                    ddt[2] +ddt[5] +ddt[8] ]
            
            # Compute V*( dot(del,dot(tau,c)) +rho'*dot(f,c) )
            for i in xrange(3):
                rhs[i,...] = vol*(ddt[i] +rho*bforce[i])
            
            # Compute face pressures for cells subject to the Neumann
            # boundary condition.
            kfp[0,0,0] = 0
            jfp[0,0,0] = 0
            ifp[0,0,0] = 0

            # Extend 0 cell along i-axis.            
            for i in xrange(1,ncells[2]+1):
                ifp[0,0,i] = rhs[2,0,0,i-1]/ifa +ifp[0,0,i-1]
                
            pres[0,0,:] = 0.5*(ifp[0,0,1::] +ifp[0,0,0:-1])
            kfp[0,0,:]  = pres[0,0,:]
            jfp[0,0,:]  = pres[0,0,:]
            
            # Now extend the just computed i-row down the j-axis.
            for j in xrange(1,ncells[1]+1):
                jfp[0,j,:] = rhs[1,0,j-1,:]/jfa +jfp[0,j-1,:]
                
            pres[0,:,:] = 0.5*(jfp[0,1::,:] +jfp[0,0:-1,:])
            kfp[0,:,:]  = pres[0,:,:]
            
            # Advance the k=0 plane forward.
            for k in xrange(1,ncells[0]+1):
                kfp[k,:,:] = rhs[0,k-1,:,:]/kfa +kfp[k-1,:,:]
                
            # Compute face pressure on remaining Neumann boundaries
            # (j=0, i=0)
            pres[:,0,:] = 0.5*(kfp[1::,0,:] +kfp[0:-1,0,:])
            jfp[:,0,:]  = pres[:,0,:]
            
            pres[:,:,0] = 0.5*(kfp[1::,:,0] +kfp[0:-1,:,0])
            ifp[:,:,0]  = pres[:,:,0]
            
            # Compute the pressure on the remaining j and i faces.
            for j in xrange(1,ncells[1]+1):
                jfp[:,j,:] = rhs[1,:,j-1,:]/jfa +jfp[:,j-1,:]
                 
            for i in xrange(1,ncells[2]+1):
                ifp[:,:,i] = rhs[2,:,:,i-1]/ifa +ifp[:,:,i-1]
                
            # Now compute the cell-centered pressure.
            pres[:,:,:] = kfp[1::,:,:] +kfp[0:-1,:,:] \
                         +jfp[:,1::,:] +jfp[:,0:-1,:] \
                         +ifp[:,:,1::] +ifp[:,:,0:-1]
            
            pres = nrm*pres
            
            e.addVars( pivlib.cpivvar(pres,'P','kPA') )
