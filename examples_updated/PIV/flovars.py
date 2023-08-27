#
# WHN 3/20/09.
# Recipe for additional flow variable construction.
#

from spivet import steps, pivlib
from numpy import *

steps.enable_user_steps()
import usersteps as usteps

# CONFIG: dltae0
conf_dltae0 = {'varnm':'T-NU',
               'pdname':'pivdata'}

# CONFIG: vorticity
conf_vorticity = {'pdname':'pivdata',
                  'velnm':'U-MF-SN-ZD-GS'}

# CONFIG: divergence
conf_divergence = {'pdname':'pivdata',
                   'varnm':'U-MF-SN-ZD-GS'}

# CONFIG: viscosity
conf_viscosity = {'cal':[0.001,13.77,-0.1468,5.355e-4],
                  'tmpnm':'T-NU',
                  'pdname':'pivdata'}

# CONFIG: spgravity
conf_spgravity = {'cal':[-4.466e-4,1.452],
                    'tmpnm':'T-NU',
                    'pdname':'pivdata'}

# Execute.
pd = pivlib.loadpivdata('PIVDATA-WFV.ex2')
carriage = {'pivdata':pd}

t = usteps.dltae0()
t.setConfig(conf_dltae0)
t.setCarriage(carriage)
t.execute()

t = usteps.vorticity()
t.setConfig(conf_vorticity)
t.setCarriage(carriage)
t.execute()

t = usteps.divergence()
t.setConfig(conf_divergence)
t.setCarriage(carriage)
t.execute()

t = usteps.viscosity()
t.setConfig(conf_viscosity)
t.setCarriage(carriage)
t.execute()

t = usteps.spgravity()
t.setConfig(conf_spgravity)
t.setCarriage(carriage)
t.execute()

# Store results.
pd.save('PIVDATA-WFV.ex2')
