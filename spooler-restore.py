import numpy
from spooler import Spooler
from psi4utils import StatusID
import psi4utils as utils
import pickle, os, glob, sys
import time



if os.path.exists("spooler.restore.lock"):
	print("there is already a restore queue! try later")
	sys.exit()

# provide a list of jobIDs and the state they should end up with
restorees = {}
#restorees["21447_0"] = StatusID.STATUS_GEOOPT_FINISHED
restorees["47922_0"] = StatusID.STATUS_GEOOPT_FINISHED


'''
fin = open("autorestore.txt", "r")
lines = fin.readlines()
fin.close()
for l in lines:
        restorees[l.strip()] = StatusID.STATUS_GEOOPT_SETUP
print(restorees)
'''

fbin = open("spooler.restore.data", "wb")
pickle.dump(restorees, fbin)
fbin.close()

print("added new restore jobs to the spooler")
time.sleep(1)
os.system("touch spooler.restore.lock")

