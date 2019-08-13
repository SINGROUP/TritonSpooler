import numpy
from spooler import Spooler
import psi4utils as utils
import pickle, os, glob




spl = Spooler()
utils.SetSpooler(spl)

'''
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
if len(cidTodo) > 0:
	jobs = utils.XYZ_pubchem_fetch_all(cidTodo)
	#print(jobs)
	spl.jobs = {**jobs, **spl.jobs}
	spl.SaveStatus()
'''

'''
# remove irreparable fails? --------------------------------------------
toremove = []
for j in spl.jobs.keys():
	if spl.jobs[j].state == utils.StatusID.STATUS_XYZ_FAIL_PUBCHEM:
		toremove.append(j)
for j in toremove:
	print("PRESPOOLER: removing job {}".format(spl.jobs[j]))
	del spl.jobs[j]
# --- ------------------------------------------------------------------
'''
'''
# recovery feature
# list of folders to reintriduce to the spooler since they failed
toadd = glob.glob("./recovery/*/molecule_*")

for folder in toadd:
	
	w = folder.split('/')
	newstate = utils.StatusID[w[2]]
	jobID = w[3].split('ecule_')[1]
	spl.jobs[jobID].state = newstate
	
	print("PRESPOOLER: recovring job {}".format(spl.jobs[jobID]))
	os.system("mv {} .".format(folder))

	' ''
	# this is allows us to change some of the job parameters (memory, time, ...) on recovery
	recoveryinfo = "{}/recovery.data".format(folder)
	if os.path.exists(recoveryinfo):
		fbin = open(recoveryinfo, "rb")
		newprops = pickle.load(fbin)
		fbin.close()
		for key in newprops.keys():
			
	' ''

'''
print(spl.jobs)

spl.Spool()




'''
electrons		time ccsd [min]
110			ok
136			ok
'''
