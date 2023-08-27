#
# WHN 3/20/09.
# After synchronizing, the PIVDATA object contains a lot of files that we
# don't need.  This script is meant to be run after flovars, and will 
# remove the unneeded historical variables from the PIVDATA-WFV file.
# These historical variables will remain in the PIVDATA-SYNC file should
# they be needed again.
#

from spivet import pivlib, steps

steps.enable_user_steps()
import usersteps as usteps

# CONFIG: prune
conf_prune = {'varnm':[#'T','T-MF','T-MFFLG','T-MF-SN',
                       'U','U-MF','U-MF-SN','U-MF-SN-ZD',
                       'U-MFFLG',
                       #'PLNR-TIME',
                       'PLNR-DLTATIME',
                       'R0','R0INAC','R0-MF','R0-MFFLG',
                       'R0-MF-MF','R0-MF-MFFLG','R0CMAX',
                       'R0CMAX-TR',
                       'R1','R1INAC','R1-MF','R1-MFFLG',
                       'R1-MF-MF','R1-MF-MFFLG','R1CMAX',
                       'R1CMAX-TR'
                       ],
              'pdname':'pivdata'}


# Execute.
pd = pivlib.loadpivdata('PIVDATA-SYNC.ex2')
carriage = {'pivdata':pd}

t = usteps.prune()
t.setConfig(conf_prune)
t.setCarriage(carriage)
t.execute()

pd = carriage['pivdata']
pd.save("PIVDATA-WFV.ex2")
