#!/bin/bash
# XB 4/9/2022
# A reference flow for PIV and postprocessing
# better run each script by hand
  
# generated file list for each epoch in folder EFRMLSTS/
python bldsched.py

# choose one of the recipe for piv (*precipe*py)
# here hybrid CCPIV and optical flow 
# using modified Brox 2004 method
# with parallel GPU/CPU multilprocessing
python cc_opyf_precipe.py

# optional, if an error occured during the piv step,
# we end up with a series seperate exodus file for each epoch
# once we rerun piv for failed epochs, run assemble.py to
# generate a combined file
#python assemble.py

# interpolate to a selected plane
# modify sync.py and usersteps.py
# to choose different plane other than the last
# e.g. 0.5 epoch 
python sync.py

# to remove unwanted entry in ex2 file
python prunewfv.py

#optional below

# additional flow variable construction
python flovars.py

# trace flow in a cylinder or plane
python procdata_cyl.py
python procdata_flat.py

# temperature reconstruction
# assuming full space velocity available
# and advect temperature from IC/BC
python thermal2.py
