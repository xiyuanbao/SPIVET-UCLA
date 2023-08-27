#
# WHN 3/20/09.
# Recipe for synchronizing planes.
#


from spivet import pivlib, steps
from numpy import *

steps.enable_user_steps()
import usersteps as usteps

# Set the reps for z-velocity filtering.  Update the denominator with the
# inter-frame delay.
zreps = 0.1*0.174/14.4

# CONFIG: conf_pzmedfltr 
conf_pzmedfltr = {'varnm':'U',
                  'planar':True,
                  'rthsf':1.7,
                  'reps':zreps,
                  'mfits':3,
                  'mfdim':5,
                  'cndx':0,
                  'pdname':'pivdata'}

# CONFIG: synchronize
conf_synchronize = {'varnm':'U-MF',#'T-MF'
                    'pdname':'pivdata'}

# CONFIG: zdpvadj
conf_zdpvadj = {'varnm':'U-MF-SN',
                'cvbndx':[[0,71],[0,98]],#should be x y range
                'rlxf':0.6,
                'nits':3,
                'kfhcw':[0.,0.],
                'pdname':'pivdata'}

"""
# CONFIG: conf_zmedfltr 
conf_zmedfltr = {'varnm':'U-MF-SN-ZD',
                 'planar':False,
                 'rthsf':1.7,
                 'reps':zreps,
                 'mfits':3,
                 'mfdim':5,
                 'cndx':0,
                 'pdname':'pivdata'}
"""

# CONFIG: gsmooth
conf_gsmooth = {'pdname':'pivdata',
                'varnm':'U-MF-SN-ZD',#'T-MF-SN'
                'gbsd':0.7}

# Execute.
pd = pivlib.loadpivdata('PIVDATA.ex2')
carriage = {'pivdata':pd}

t = steps.medfltr()
t.setConfig(conf_pzmedfltr)
t.setCarriage(carriage)
t.execute()

t = usteps.synchronize()
t.setConfig(conf_synchronize)
t.setCarriage(carriage)
t.execute()

t = usteps.zdpvadj()
t.setConfig(conf_zdpvadj)
t.setCarriage(carriage)
t.execute()

"""
t = steps.medfltr()
t.setConfig(conf_zmedfltr)
t.setCarriage(carriage)
t.execute()
"""

t = steps.gsmooth()
t.setConfig(conf_gsmooth)
t.setCarriage(carriage)
t.execute()

pd = carriage['pivdata']
pd.save("PIVDATA-SYNC.ex2")

