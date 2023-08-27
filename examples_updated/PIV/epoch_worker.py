# XB 4/9/2022
# library used to run parallel and hybrid GPU/CPU PIV 
from numpy import *
from scipy import interpolate
from spivet import pivlib, tlclib, spivetconf
from spivet.steputil import _parsefn, _fretrieve
import sys, os, urlparse

from spivet import compat

from spivet.steps import global_accumulator
import pickle
import socket, traceback

import subprocess as sp
import time

from random import seed
from random import random
def get_gpu_memory():
    _output_to_list = lambda x: x.decode('ascii').split('\n')[:-1]

    ACCEPTABLE_AVAILABLE_MEMORY = 1024
    COMMAND = "nvidia-smi --query-gpu=memory.free --format=csv"
    memory_free_info = _output_to_list(sp.check_output(COMMAND.split()))[1:]
    memory_free_values = [int(x.split()[0]) for i, x in enumerate(memory_free_info)]
#     print(memory_free_values)
    return memory_free_values

def wait_gpu(interval, mem_need):
	while True:
                mem = get_gpu_memory()[0]
                if mem < mem_need:
                        time.sleep(interval)
                else:
                        break


def _loop_epoch_worker_gpu(carriage,steps,cnt,obpath):
    """
    Performs the loop work. sflg will be set if an exception is thrown.
    If an exception is thrown, then sepdvn or epivdata will be set
    to the traceback text.
    
    If running in parallel, stores a PIVDATA object file directly
    to the NWS.  Returns [sflg,sepdvn,stdout,stderr] where
        sflg ------ Success flag.  0 if successfull, non-zero otherwise.
        sepdvn ---- The variable name on the NWS where the PIVDATA
                    is stored.
        stdout ---- The stdout transcript.
        stderr ---- The stderr transcript.
        
    If running locally, returns [sflg,epivdata,stdout,stderr] where
        sflg ------ Success flag.  0 if successfull, non-zero otherwise.
        epivdata -- A valid PIVData object.
        stdout ---- The stdout transcript.
        stderr ---- The stderr transcript.        
    """
    import socket, traceback
    #import openpiv.gpu_process
    # Initialization.
    epdfn = "PIVDATA.ex2"

    lsofn      = "stdout"+str(cnt)
    sys.stdout = open(lsofn,'w', 0)

    lsefn      = "stderr"+str(cnt)
    sys.stderr = open(lsefn,'w', 0)

    rhost = socket.gethostname()
    print "Running on: %s" % rhost

    #steps = globals()[wrkrstepsvn]
    try:
        rnk = SleighRank
        parallel = True
        print "Running remotely."
    except:
        print "Running locally."
        parallel = False

    try:
        if ( parallel ):
            svnext = "%i-%i" % (SleighRank,random.randint(0,1000000))

        # Create the step objects and store the config dictionaries.
        steplst = []
        for step in steps:
            cname   = step[0]
            modname = step[1]
            config  = step[2]

            __import__(modname)
            module = sys.modules[modname]

            cobj = module.__dict__[cname]()
            cobj.setConfig(config)

            steplst.append(cobj)

        # Execute the steps.
        for step in steplst:
            step.setCarriage(carriage)
            step.execute()

        epd = carriage['epivdata']
        if ( parallel ):
            # Write the PIVDATA object to a file and then load in directly
            # into the network space.  This is necessary to prevent Python 
            # ASCII serialization of the PIVDATA object (which can cause the 
            # process to run out of memory).
            epd.save(epdfn)

            sepdvn = "PIVDATA-%s-%s" % (rhost,svnext)
            SleighUserNws.declare(sepdvn,'single')
            fh = open(epdfn,'rb')
            SleighUserNws.storeFile(sepdvn,fh)
            fh.close()

        sflg = 0
        if ( parallel ):
            pdo = sepdvn
        else:
            pdo = epd

    except KeyboardInterrupt:
        # Need to handle Ctrl-C explicitly and die.
        raise

    except:
        sflg = 1
        traceback.print_exc(file=sys.stderr)
        pdo  = traceback.format_exc()

    # Close out the local logs.
    sys.stdout.close()
    sys.stdout = sys.__stdout__

    sys.stderr.close()
    sys.stderr = sys.__stderr__

    # Load the log files.
    fh      = open(lsofn,'r')
    stdoutf = fh.read()
    fh.close()

    fh      = open(lsefn,'r')
    stderrf = fh.read()
    fh.close()

    assemble_flag=global_accumulator(obpath,[[cnt,[sflg,pdo,stdoutf,stderrf]]])
    return assemble_flag

  
if __name__ =="__main__":
	sys.argv[1]
	with open(sys.argv[1], "rb") as fp:   # Unpickling
		crg = pickle.load(fp) 
	with open(sys.argv[2], "rb") as fp:   # Unpickling
		steps = pickle.load(fp)
	cnt = int(sys.argv[3])
	obpath = sys.argv[4]
	try:
		wait_time = int(cnt/20.0)*600+cnt
		time.sleep(wait_time)
		
		interval = 10
		mem_need = 3000#Actually 389MB to 1200MB
		#wait until gpu memory is available
		while True:
			
			wait_gpu(interval, mem_need)
			seed(cnt)
			for i in range(cnt):
				#wait_time = 10*random()
				wait_time =2* cnt
				wait_gpu(wait_time, 3500)
			print "start cnt:",cnt
			try:
				flag = _loop_epoch_worker_gpu(crg, steps, cnt,  obpath)
				print "epoch finished, flag:",str(flag)
				break
			except:
				print "except here!"
				continue
		print flag	
	except KeyboardInterrupt:
        	# Need to handle Ctrl-C explicitly and die.
        	raise	
	#print _loop_epoch_worker_gpu(crg, steps, cnt,  obpath)
