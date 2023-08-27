from spivet import steps
import os

# CONFIG: bldsched
root_dir='file:///home/xbao/LAB/'
bfileurl = root_dir+"04132022/DATA"
conf_bldsched = {'bfileurl':bfileurl}

# Get files.
files = os.listdir('../DATA')
c = 0
for i in range(len(files)):
    f = files[c]
    if ( f[-4::] != ".tif" ):
        files.pop(c)
    else:
        c = c +1

# Run the step.
carriage = {'files':files}

t = steps.bldsched()
t.setConfig(conf_bldsched)
t.setCarriage(carriage)
t.execute()
