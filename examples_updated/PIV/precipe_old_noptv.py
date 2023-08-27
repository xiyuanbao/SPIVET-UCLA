#
# WHN 3/20/09.
# Recipe for datasets collected using planar groups.
#

from spivet import steps
from numpy import *

steps.enable_user_steps()
import usersteps as usteps
root_dir='file:///home/xbao/LAB_test/'
# PIVDICT
wicsp = [root_dir+'SYRUPCAL-06132008/CALIBRATION/WICSP_CAM0',root_dir+'SYRUPCAL-06132008/CALIBRATION/WICSP_CAM1']
camcal = [root_dir+'SYRUPCAL-06132008/CALIBRATION/CAMCAL_CAM0',
root_dir+'SYRUPCAL-06132008/CALIBRATION/CAMCAL_CAM1']
tlccal = root_dir+'TLCCAL-09282008/TLCCALIBRATION/TLCCAL_CAM1'

padpix = 32
imdim  = array([1357, 1030]) +2*array([padpix,padpix])
drbndx = array([[54+padpix,-23-padpix],
                [35+padpix,-35-padpix]])
rbndx  = array([[0,imdim[0]],[0,imdim[1]]]) +drbndx

pivdict={
    'gp_rbndx':rbndx,
    'gp_bsize':(32,32),
    'gp_bolap':(0,0),
    'gp_bsdiv':2,
    'ir_eps':0.003,
    'ir_maxits':100, 
    'ir_mineig':0.05,
    'ir_imthd':'C',
    'ir_iedge':0,
    'ir_tps_csrp':0.1,
    'ir_tps_ithp':90,
    'ir_tps_wsize':[5,5],
    'ir_tps_sdmyf':0.1,
    'ir_tps_alpha':0.8,
    'ir_tps_beta':1.5,
    'ir_tps_csize':15,
    'ir_tps_nits':30,
    'ir_tps_scit':5,
    'ir_tps_annl':0.98,
    'of_maxdisp':(34,34), 
    'of_rmaxdisp':(5,5), 
    'of_hrsc':False,
    'of_nccth':-999,
    'of_highp':False,
    'of_bsdcmx':0.,
    'pg_wicsp':wicsp,
    'pg_camcal':camcal,
    'tc_tlccal':tlccal,
    'tc_tlccam':1,
    'tc_ilvec':[0.,0.,-1.],
    'tc_interp':False
}

# CONFIG: loadimg
conf_loadimg = {'rgb':True}

# CONFIG: coloradj
adjmat = root_dir+'COLORADJ/COLORADJ'
conf_coloradj = {'adjmat':adjmat}

# CONFIG: dewarpimg
conf_dewarpimg = {'wicsp':wicsp,'camcal':camcal,'padpix':padpix}

# CONFIG: oflow2d
conf_oflow2d = {'pivdict':pivdict}

# CONFIG: pdmedfltr
conf_pdmedfltr3 = {'varnm':['R0','R1'],
                  'planar':True,
                  'rthsf':3.,
                  'reps':0.1,
                  'mfits':3,
                  'mfdim':3,
                  'pdname':'pivdata'}

# CONFIG: pdmedfltr
conf_pdmedfltr5 = {'varnm':['R0-MF','R1-MF'],
                  'planar':True,
                  'rthsf':4.,
                  'reps':0.1,
                  'mfits':3,
                  'mfdim':5,
                  'pdname':'pivdata'}

# CONFIG: pdgsmooth
conf_pdgsmooth = {'varnm':['R0-MF-MF','R1-MF-MF'],
                  'planar':True,
                  'gbsd':1.,
                  'pdname':'pivdata'}

# CONFIG: refine2dof
mf3d = {'filter':'medfltr','fdim':3,'rthsf':4.,'reps':0.1,'planar':True,'nit':3,
        'cndx':None}
mf5d = {'filter':'medfltr','fdim':5,'rthsf':4.,'reps':0.1,'planar':True,'nit':3,
        'cndx':None}
gsd = {'filter':'gsmooth','gbsd':1.,'planar':True,'nit':1,'cndx':None}

rpivdict = pivdict.copy()

ifurl = root_dir+"R2DIMGFLTR/R2DIMGFLTR"

##### CHECK.  The second term is for the extra two cells along y.
crbndx = pivdict['gp_bsdiv']*array([[27,39],[8,30-8]]) +array([[1,1],[0,0]])

conf_refine2dof = {'pivdict':rpivdict,
                   'varnm':['R0-MF-MF-GS','R1-MF-MF-GS'],
                   'crbndx':crbndx,
                   'rfactor':1.2,
                   'its':20,
                   'eps':0.01,
                   'planes':[],
                   'cfltrprm':[mf3d,mf5d,gsd],
                   'ifltr':[ifurl,ifurl]}

# CONFIG: oflow3d
conf_oflow3d = {'pivdict':pivdict,
                'varnms':['R0-MF-MF-GS-TR','R1-MF-MF-GS-TR','R0INAC','R1INAC'],
                'zcellsz':-5.0}

# CONFIG: disp2vel
conf_disp2vel = {}

# CONFIG: recordtime
conf_recordtime = {}

# CONFIG: bimedfltr
conf_bimedfltr = {'rbndx':pivdict['gp_rbndx'],
                  'bsize':pivdict['gp_bsize'],
                  'prccam':[False,True],
                  'prcchnl':[True,False,False]}

# CONFIG: tlcmask
conf_tlcmask = {'pivdict':pivdict,
                'prccam':[False,True]}

# CONFIG: sctfcomp
conf_sctfcomp = {'pivdict':pivdict,
                 'zcellsz':-5.0,
                 'tlcp0swz':0.}

# CONFIG: loop_plane
conf_loop_plane = {'steps':[['loadimg','spivet.steps',conf_loadimg],
                            ['coloradj','usersteps',conf_coloradj],
                            ['dewarpimg','spivet.steps',conf_dewarpimg],
                            ['oflow2d','spivet.steps',conf_oflow2d],
                            ['medfltr','spivet.steps',conf_pdmedfltr3],
                            ['medfltr','spivet.steps',conf_pdmedfltr5],
                            ['gsmooth','spivet.steps',conf_pdgsmooth],
                            ['refine2dof','spivet.steps',conf_refine2dof],
                            ['oflow3d','spivet.steps',conf_oflow3d],
                            ['disp2vel','spivet.steps',conf_disp2vel],
                            ['recordtime','spivet.steps',conf_recordtime],
                            ['bimedfltr','spivet.steps',conf_bimedfltr],
                            ['tlcmask','spivet.steps',conf_tlcmask],
                            ['sctfcomp','spivet.steps',conf_sctfcomp]]}


# CONFIG: ptmedfltr
conf_ptmedfltr = {'varnm':['T'],
                  'planar':True,
                  'rthsf':2.,
                  'reps':0.15,
                  'mfits':1,
                  'mfdim':9,
                  'pdname':'epivdata'}

# CONFIG: loop_epoch
conf_loop_epoch = {'steps':[['loop_plane','spivet.steps',conf_loop_plane],
                            ['medfltr','spivet.steps',conf_ptmedfltr]],
                   'parallel':True}

# Execute.
carriage = {'lbfpath':'EFRMLSTS'}

t = steps.loop_epoch()
t.setConfig(conf_loop_epoch)
t.setCarriage(carriage)
t.execute()

# Store results.
pd = carriage['dpivdata'][0]
pd.save('PIVDATA.ex2')
