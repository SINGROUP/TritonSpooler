import numpy
from spooler import Spooler
from psi4utils import StatusID
import psi4utils as utils
import pickle, os, glob
import time
from collections import Counter


spl = Spooler()

#cidlist = numpy.arange(23200,23400) # these got interrupted by the scratch FS fail! fuck
# CID 34694, 37174 seems to not exist!
#cidlist = numpy.arange(36500, 37000) # DONE - triton crashed during this probably!
cidlist = numpy.arange(47000, 48000)

# determine which molecules need downloading
cidTodo = []
for cid in cidlist:
	jobID = "{}_0".format(cid)
	if jobID in spl.jobs.keys():
		continue
	else:
		cidTodo.append(cid)

print("CIDs to download:",cidTodo)


# from here we have to get molecule xyzs and create jobs for the spooler
newjobs = {}

if len(cidTodo) > 0:
        jobs = utils.XYZ_pubchem_fetch_all(cidTodo)
        #print(jobs)
        newjobs = {**jobs, **newjobs}
        jobIDs = newjobs.keys()
        jobSts = [newjobs[k].state for k in jobIDs]
        print('SUMMARY:')
        print(Counter(jobSts))

        fbin = open("spooler.resupply.data", "wb")
        pickle.dump(newjobs, fbin)
        fbin.close()

        print("added new jobs to the spooler")
        time.sleep(1)
        os.system("touch spooler.resupply.lock")

