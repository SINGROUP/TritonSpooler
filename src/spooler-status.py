import numpy
from spooler import Spooler
from psi4utils import StatusID
import pickle, os, glob


spl = Spooler()

print("active tasks:")
keys = list(spl.tasks.keys())
keys.sort()
for t in keys: 
        print(t,spl.tasks[t], spl.jobs[spl.tasks[t]].state)


print("Active jobs:",len(spl.tasks))
#lst = [k for k in spl.jobs.keys() if spl.jobs[k].state != StatusID.STATUS_COMPLETED]
#lst.sort()
#print(lst)
print("completed:",len(spl.GetJobsByState(StatusID.STATUS_COMPLETED)))


print("showing state specifics...")
states = [StatusID.STATUS_GEOOPT_SCF_FAIL, StatusID.STATUS_GEOOPT_FAIL, StatusID.STATUS_GEOOPT_FAIL_SEGFAULT, StatusID.STATUS_GEOOPT_FAIL_NOCHEMSPIDER, StatusID.STATUS_CCSD_FAIL, StatusID.STATUS_CCSD_FAIL_UNKO, StatusID.STATUS_CCSD_FAIL_TIME]
for s in states:
	lst = spl.GetJobsByState(s)
	print("State: ", s)
	for a in lst:
		print(a)


# MANUAL FIX FOR MANUALLY FIXED RUNS
'''
#jobs = ["282_0", "608_0", "584_0", "435_0"]
for j in spl.jobs.keys():
	#spl.jobs[j].state = StatusID.STATUS_COMPLETED
	spl.jobs[j].geoopt_failstate = 0

spl.SaveStatus()
'''
