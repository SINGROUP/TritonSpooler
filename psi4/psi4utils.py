import subprocess
import numpy
from enum import Enum
import json
import os
import glob
from spooler import SpoolerJob
import math

from mendeleev import element as Element
from chemspipy import ChemSpider
import pubchempy
cs = ChemSpider("xONbOGbUKNBxBwtz9jdNfq6UxHGkVr1m")




MAX_HEAVY_ATOMS = 24
MAX_ELECTRONS = MAX_HEAVY_ATOMS * 6

class StatusID(Enum):
	
	STATUS_XYZ_FAIL_PUBCHEM = 0
	STATUS_XYZ_FAIL_ATOMS = 1
	STATUS_XYZ_FAIL_ELECT = 2
	STATUS_XYZ_FAIL_ROHF = -1
	STATUS_XYZ_FAIL_ION = -2
	STATUS_XYZ_FAIL_BASIS = -3
	STATUS_STARTUP = 3					# xyz was downloaded but no more!
	STATUS_GEOOPT_SETUP = 4			# molecule folder is there (with the geoopt input files ready to roll)
	STATUS_GEOOPT_RUNNING = 5		# geoopt has be submitted to the queue
	STATUS_GEOOPT_FINISHED = 6	# geoopt is not running anymore
	
	STATUS_GEOOPT_SCF_FAIL = -110
	STATUS_GEOOPT_FAIL = -100		# geoopt has failed miserably and unrecoverably (unlikely)
	STATUS_GEOOPT_FAIL_BASIS = -120
	STATUS_GEOOPT_FAIL_SEGFAULT = -140
	STATUS_GEOOPT_FAIL_NOCHEMSPIDER = -150
	STATUS_GEOOPT_FAIL_MEM = -160
	
	
	STATUS_GEOOPT_RELOAD = 101	# need to load the last checkpoint and reoptimize
	STATUS_GEOOPT_DONE = 100		# geoopt is good and done
	

	STATUS_CCSD_SETUP = 200			# ccsdrun has all needed input files, ready to start running
	STATUS_CCSD_RUNNING = 201		# ccsd is running
	STATUS_CCSD_FINISHED = 202		# ccsd is no longer running
	
	STATUS_CCSD_FAIL = -300
	STATUS_CCSD_TIME = 301			# ccsd needs more time...
	STATUS_CCSD_DONE =  300
	
	STATUS_CCSD_FAIL_MEM = -400 # exceeded memory limit
	STATUS_CCSD_FAIL_TIME= -450 # exceeded memory limit
	STATUS_CCSD_FAIL_UNKO= -470 # exceeded memory limit
	
	
	STATUS_CCSD_FAIL_TOTAL = -1000 # ccsd failed completely, nothing more to do
	
	STATUS_COMPLETED = 1000			# itz all done!


class Simulation(SpoolerJob):
	
	def __init__(self, jobID, types=None, xyz=None):
		
		super().__init__(jobID, None)
		
		self.inchikey = ""
		self.nthreads = 10
		self.queue = "batch"
		self.arch = "[skl]"
		self.time = 120
		self.cachelevel = 2
		self.memory = 1000
		self.timefactor = 1
		
		self.types = types
		self.xyz = xyz
		
		self.geoopt_failstate = 0
		
		self.directory =  "molecule_{}".format(self.ID)
		
		self.input_go =  "{}/geoopt.py".format(self.directory)
		self.job_go = "{}/geoopt.job".format(self.directory)
		self.output_go = "{}/geoopt.out".format(self.directory)
		
		self.input_cc =  "{}/ccsdrun.py".format(self.directory)
		self.job_cc = "{}/ccsdrun.job".format(self.directory)
		self.output_cc = "{}/ccsdrun.out".format(self.directory)
		
		# compute useful properties
		self.nheavy = 0
		self.nelect = 0
		if not(types is None):
			for e in types:
				if e != 1: self.nheavy += 1
				self.nelect += e
		# ---

	
	def LoadLastXYZ(self):
		
		natm = self.xyz.shape[0]
		
		# find all occurrences of: "Writing optimization data to binary file."
		marker = "Writing optimization data to binary file."
		fin = open(self.output_go, "r")
		lines = fin.readlines()
		fin.close()
		
		indexes = []
		for i in range(len(lines)):
			if (marker in lines[i]) and len(lines)-i > natm+4:
				if "Structure for next" in lines[i+1]:
					indexes.append(i)
		
		# get the last good one
		if len(indexes) > 0: indexes = indexes[-1]

		lines = lines[indexes+3 : ]
		lines = lines[0:natm]
		print(lines)
		
		for i in range(natm):
			w = lines[i].strip().split()
			print(w)
			for c in range(3):
				self.xyz[i,c] = float(w[c+1])
	
	
	def XYZ2String(self):
		
		strout = "FLAGXYZ\n"
		for i in range(len(self.types)):
			strout = strout + "{}\t{} {} {}\n".format(self.types[i], self.xyz[i,0], self.xyz[i,1], self.xyz[i,2])
		
		return strout
		
	
	def GetTimeString(self):

		t = int(self.time * self.timefactor)
		t = min(t, 5*24*60)
		dd = math.floor(t/(24*60))
		t -= dd * 24 * 60
		hh = math.floor(t/60)
		t -= hh * 60
		mm = t

		t = "{:1d}-{:02d}:{:02d}:00".format(dd,hh,mm)
		return t

	def WriteInput_geoopt(self):
		
		# copy the file
		if self.geoopt_failstate != 2:
			os.system("cp template-geoopt.py " + self.input_go)
		else:
			os.system("cp template-geoopt-mp2.py " + self.input_go)
		
		# write xyz
		os.system("perl -i -pe 's/FLAGXYZ\n/{}/g' {} ".format(self.XYZ2String(), self.input_go))
		
		# write parameters
		os.system("sed -i -e 's/CCCACHELEVEL/{}/g' {}".format(self.cachelevel, self.input_go))
		os.system("sed -i -e 's/NUMTHREADS/{}/g' {}".format(self.nthreads, self.input_go))

		# copy the job file if it doesnt already exist
		os.system("cp template-geoopt.job " + self.job_go)
		
		# write parameters
		os.system("sed -i -e 's/NTHR/{}/g' {}".format(self.nthreads, self.job_go))
		os.system("sed -i -e 's/CPUTYPE/{}/g' {}".format(self.arch, self.job_go))
		os.system("sed -i -e 's/QUEUE/{}/g' {}".format(self.queue, self.job_go))
		os.system("sed -i -e 's/TIMELIM/{}/g' {}".format(self.GetTimeString(), self.job_go))
		os.system("sed -i -e 's/MEMLIM/{}/g' {}".format(self.memory, self.job_go))
		
		
	def WriteInput_ccsd(self):
		
		# copy the file
		os.system("cp template-ccsd.py " + self.input_cc)
		
		# write xyz
		os.system("perl -i -pe 's/FLAGXYZ\n/{}/g' {} ".format(self.XYZ2String(), self.input_cc))
		
		# write parameters
		os.system("sed -i -e 's/CCCACHELEVEL/{}/g' {}".format(self.cachelevel, self.input_cc))
		os.system("sed -i -e 's/NUMTHREADS/{}/g' {}".format(self.nthreads, self.input_cc))

		# copy the job file if it doesnt already exist
		os.system("cp template-ccsd.job " + self.job_cc)
		
		# write parameters
		os.system("sed -i -e 's/NTHR/{}/g' {}".format(self.nthreads, self.job_cc))
		os.system("sed -i -e 's/CPUTYPE/{}/g' {}".format(self.arch, self.job_cc))
		os.system("sed -i -e 's/QUEUE/{}/g' {}".format(self.queue, self.job_cc))
		os.system("sed -i -e 's/TIMELIM/{}/g' {}".format(self.GetTimeString(), self.job_cc))
		os.system("sed -i -e 's/MEMORY/{}/g' {}".format(self.memory, self.job_cc))


def TimeString2Mins(timestring):
        
	w = timestring.split('-')
	dd = int(w[0])
	others = w[1]
	w = others.split(':')
	hh = int(w[0])
	mm = int(w[1])

	totmins = mm + hh*60 + dd * 24 * 60
	return totmins

def Mins2TimeString(mins):
        
	t = mins
	dd = math.floor(t/(24*60))
	t -= dd * 24 * 60
	hh = math.floor(t/60)
	t -= hh * 60
	mm = t

	t = "{:1d}-{:02d}:{:02d}:00".format(dd,hh,mm)

	return t

PSI4JobInfo = {}

#
# returns a dictionary of jobID: status
def XYZ_pubchem_fetch(cid):
	
	#print(cid)
	comp = pubchempy.Compound.from_cid(int(cid))
	inchikey = comp.inchikey
	charge = comp.to_dict()['charge']
	print(cid, inchikey, charge)
	
	cmd = 'curl -s "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/JSON?record_type=3d"'.format(cid)
	req = subprocess.run([cmd], shell=True, stdout=subprocess.PIPE)
	info = json.loads(req.stdout.decode('utf-8'))

	
	# key = jobID as string
	# value = Status code
	jobs = {}
	jobID = "{}_0".format(cid) # temporary jobID
	jobobj = Simulation(jobID, None, None)
	
	
	# could not get molecule, failed job
	if "Fault" in info.keys():
		print("Cannot get molecule CID {} from pubchem!".format(cid))
		jobobj.state = StatusID.STATUS_XYZ_FAIL_PUBCHEM
		jobs[jobobj.ID] = jobobj
		return jobs
	
	if charge != 0:
		print("Molecule CID {} is charged!".format(cid))
		jobobj.state = StatusID.STATUS_XYZ_FAIL_ION
		jobs[jobobj.ID] = jobobj
		return jobs
		
	# code here => we got a workable pubchem output
	info = info["PC_Compounds"][0]
	atoms = info['atoms']
	elements = numpy.asarray(atoms['element'], dtype=numpy.int32)
	if not ((elements < 37).all()):
		print("Molecule CID {} has basisless atoms!".format(cid))
		jobobj.state = StatusID.STATUS_XYZ_FAIL_BASIS
		jobs[jobobj.ID] = jobobj
		return jobs
	
	jobs[jobID] = Simulation(jobID, types=elements, xyz=None)
	
	# ignore the molecule if it is too big
	if jobs[jobID].nheavy > MAX_HEAVY_ATOMS:
		print('too big')
		jobs[jobID].state = StatusID.STATUS_XYZ_FAIL_ATOMS
		return jobs
	
	if jobs[jobID].nelect > MAX_ELECTRONS:
		print('too many e-')
		jobs[jobID].state = StatusID.STATUS_XYZ_FAIL_ELECT
		return jobs
	
	if jobs[jobID].nelect % 2 == 1:
		print('too openshell')
		jobs[jobID].state = StatusID.STATUS_XYZ_FAIL_ROHF
		return jobs
	
	# code here means that the molecule is not too big

	conformerInfo = info['coords'][0]['conformers']
	jobs = {} # clean the dictionary
	
	for i in range(len(conformerInfo)):

		conformer = conformerInfo[i]
		cx = numpy.asarray(conformer['x'])
		cy = numpy.asarray(conformer['y'])
		cz = numpy.asarray(conformer['z'])
	
		jobID = "{}_{}".format(cid, i)

		
		natm = elements.shape[0]
		xyz = numpy.zeros((natm, 3))
		for a in range(natm):
			xyz[a,0] = cx[a]
			xyz[a,1] = cy[a]
			xyz[a,2] = cz[a]
		

		# check if the structure is sorted from heavier to lighter atom
		idxs = numpy.argsort(-elements)
		elements = elements[idxs]
		xyz = xyz[idxs]
		# ------------------------------------------------------------------
		
		# put the coordinates in the input file
		jobs[jobID] = Simulation(jobID, elements, xyz)
		jobs[jobID].state = StatusID.STATUS_STARTUP
		jobs[jobID].inchikey = inchikey
		# ------------------------------------------------------------------
		# JOB IS NOW GOOD TO GO!
		
	
	return jobs
	
def XYZ_pubchem_fetch_all(cids):
	
	jobs = {}
	for cid in cids:
		js = XYZ_pubchem_fetch(cid)
		jobs = {**jobs, **js}

	return jobs



def SetSpooler(spl):
	
	spl.actions[StatusID.STATUS_STARTUP] = Action_startup
	spl.actions[StatusID.STATUS_GEOOPT_SETUP] = Action_geoopt_setup
	spl.actions[StatusID.STATUS_GEOOPT_RUNNING] = Action_geoopt_running
	spl.actions[StatusID.STATUS_GEOOPT_FINISHED] = Action_geoopt_finished
	#spl.actions[StatusID.STATUS_GEOOPT_RELOAD] = Action_geoopt_reload
	
	#spl.actions[StatusID.STATUS_GEOOPT_FAIL_SEGFAULT] = Action_geoopt_fail_segfault
	#spl.actions[StatusID.STATUS_GEOOPT_FAIL_MEM] = Action_geoopt_fail_mem
	
	
	spl.actions[StatusID.STATUS_GEOOPT_DONE] = Action_geoopt_done
	spl.actions[StatusID.STATUS_CCSD_SETUP] = Action_ccsd_setup
	spl.actions[StatusID.STATUS_CCSD_RUNNING] = Action_ccsd_running
	spl.actions[StatusID.STATUS_CCSD_FINISHED] = Action_ccsd_finished
	spl.actions[StatusID.STATUS_CCSD_DONE] = Action_ccsd_done
	spl.actions[StatusID.STATUS_CCSD_FAIL] = Action_ccsd_fail
	#spl.actions[StatusID.STATUS_CCSD_FAIL_MEM] = Action_ccsd_fail_mem
	
	
	return


## Creates the compute folder with the initial input files.
## Run for jobs in status STATUS_STARTUP
# \param taskID: the taskID provided by the spooler
def Action_startup(job, taskID):
	
	print("calling Action_startup for task {} - {}".format(taskID, job))
	
	# make folder and copy template input
	os.system("mkdir -p " + job.directory)
	
	
	# some flags depends on how big the run is
	# set the cache level according to natoms
	if job.nheavy <= 6: # small jobs
		job.time = 60
		nthr = 10
	elif job.nheavy <= 11: # medium jobs
		timelim = 60
		nthr = 20
	else: # larger jobs
		timelim = 180
		nthr = 40

	job.WriteInput_geoopt()
	
	
	# the task is finished, and this is the new status
	job.state = StatusID.STATUS_GEOOPT_SETUP
	return True
	

## The point of this function is to set the correct taskID in the SLURM job script
## so that we can later read the taskID from squeue output
##
def Action_geoopt_setup(job, taskID):
	
	print("calling Action_geoopt_setup for task {} - {}".format(taskID, job))
	
	# tweak the job file
	os.system("sed -i 's/SBATCH -J .*/SBATCH -J go_{}/g' {}".format(taskID, job.job_go))
	
	slurms = glob.glob("{}/slurm-*".format(job.directory))
	
	# send to queue
	os.chdir(job.directory)
	if len(slurms) > 0:
		os.system("rm slurm-*")
	os.system("sbatch geoopt.job")
	os.chdir("../")

	job.state = StatusID.STATUS_GEOOPT_RUNNING
	
	print("ACTION GEOOPT_SETUP: GEOOPT {} submitted".format(job.ID))
	
	return False
	
## This is executed when the job is in the queue, to check if it is still there
def Action_geoopt_running(job, taskID):
	
	print("calling Action_geoopt_running for task {} - {}".format(taskID, job))
	
	# check if this calculation is still running
	req = subprocess.run(["squeue -u federif1 | grep 'go_{}'".format(taskID)], shell=True, stdout=subprocess.PIPE)
	req = req.stdout.decode('utf-8')
	
	if len(req) > 0:

		print("TASK {} JOB {} still geoopt_running...".format(taskID, job.ID))
		return False


	job.state = StatusID.STATUS_GEOOPT_FINISHED
	return True

## Called when the job has left the slurm queue
## checks if the simulation finished or there was some problem
def Action_geoopt_finished(job, taskID):
	
	print("calling Action_geoopt_finished for task {} - {}".format(taskID, job))

	# check if the run finished for real
	xyzfile = job.directory + "/GEOM-B3LYP.xyz"
	
	slurmfile = glob.glob("{}/slurm-*".format(job.directory))
	slurmfile.sort()
	if len(slurmfile) > 0: slurmfile = slurmfile[-1]
	
	
	# check if it was a segfault
	req = subprocess.run(['grep "Exceeded job memory limit" {}'.format(slurmfile)], shell=True, stdout=subprocess.PIPE)
	req = req.stdout.decode('utf-8')
	if len(req) > 5:
		print("GEOOPT {}: it went out of memory!".format(job.ID))
		job.state = StatusID.STATUS_GEOOPT_FAIL_MEM
		return True
	
	
	# check if it was out of memory
	req = subprocess.run(['grep "Segmentation fault" {}'.format(slurmfile)], shell=True, stdout=subprocess.PIPE)
	req = req.stdout.decode('utf-8')
	if len(req) > 5:
		print("GEOOPT {}: there was a segfault!".format(job.ID))
		job.state = StatusID.STATUS_GEOOPT_FAIL_SEGFAULT
		return True
	
	
	# check for bad basis set/element
	req = subprocess.run(['grep "BasisSetNotFound" {}'.format(slurmfile)], shell=True, stdout=subprocess.PIPE)
	req = req.stdout.decode('utf-8')
	if len(req) > 5:
		print("GEOOPT {}: no basis for an element in the system!".format(job.ID))
		job.state = StatusID.STATUS_GEOOPT_FAIL_BASIS
		return True


	# no output file... maybe triton crashed? retry (might get stuck in infinite loop!
	if not os.path.exists(job.output_go):
		print("GEOOPT {}: output file not found?! retrying".format(job.ID))
		job.state = StatusID.STATUS_GEOOPT_SETUP
		return True
	
	
	# check if there was an SCF failure
	req = subprocess.run(['grep "PsiException: Could not converge SCF iterations in" {}'.format(job.output_go)], shell=True, stdout=subprocess.PIPE)
	req = req.stdout.decode('utf-8')
	if len(req) > 5:
		print("GEOOPT {}: SCF took too long!".format(job.ID))
		job.state = StatusID.STATUS_GEOOPT_SCF_FAIL
		return True



	# check if optking failed in max steps... reload and resubmit!
	req = subprocess.run(['grep "**** Optimization has failed! (in 200 steps)" {}'.format(job.output_go)], shell=True, stdout=subprocess.PIPE)
	req = req.stdout.decode('utf-8')
	if len(req) > 5:
		print("GEOOPT {}: ran out of iterations...".format(job.ID))
		job.LoadLastXYZ()
		job.WriteInput_geoopt()
		job.state = StatusID.STATUS_GEOOPT_SETUP
		return True


	# check if there was some other optking failure
	req = subprocess.run(['grep "Optimizer: Optimization failed!" {} slurm*'.format(job.output_go)], shell=True, stdout=subprocess.PIPE)
	req = req.stdout.decode('utf-8')
	if len(req) > 5:
		print("GEOOPT {}: SCF took too long!".format(job.ID))
		job.state = StatusID.STATUS_GEOOPT_FAIL
		return True

	# check how many valid configurations can be seen in the output file
	fin = open(job.output_go, "r")
	lines = fin.readlines()
	fin.close()
	marker = "Writing optimization data to binary file."
	found = False
	for i in range(len(lines)):
		if marker in lines[i] and i+len(job.types)+4 < len(lines):
			found = True
			break
	if not found:
		print("GEOOPT {}: no valid rstart structures found in the output!".format(job.ID))
		job.state = StatusID.STATUS_GEOOPT_FAIL
		return True



	# code here => output file was found!
	# parse for good end flag
	req = subprocess.run(['grep "**** Optimization is complete!" {}'.format(job.output_go)], shell=True, stdout=subprocess.PIPE)
	req = req.stdout.decode('utf-8')
	if len(req) > 5:
		
		if os.path.exists(xyzfile):
			print("GEOOPT {}: completed!".format(job.ID))
			# load the xyz from the dft output
			xyz = numpy.loadtxt(xyzfile)
			xyz = xyz[:,1:]
			job.xyz = xyz
			job.state = StatusID.STATUS_GEOOPT_DONE
			return True
		else:
			print("GEOOPT {}: completed but no xyzfile... reload and restart".format(job.ID))
			job.LoadLastXYZ()
			job.state = StatusID.STATUS_STARTUP
			return True
	
	else:
		# **** Optimization is complete! was not in the output file.... bad!
		print("GEOOPT {}: run did not really finish".format(job.ID))
		job.LoadLastXYZ()
		job.state = StatusID.STATUS_STARTUP
		return True
	
	
	
	# code should never get here!
	job.state = StatusID.STATUS_GEOOPT_FAIL
	return True
	

def Action_geoopt_fail_segfault(job, taskID):
	
	print("calling Action_geoopt_fail_segfault for task {} - {}".format(taskID, job))
	
	# the GEOOPT failed due to segfault
	# 1: try again with a geometry from chemspider
	if job.geoopt_failstate == 0:
		
		query = cs.filter_inchikey(job.inchikey)
		results = cs.filter_results(query);
		
		job.geoopt_failstate = 1 # mark it so we know we already passed by here!
		
		# if chemspider doesnt have it... we cant do much! just go to the next fail level... MP2
		if len(results) == 0:
			job.state = StatusID.STATUS_GEOOPT_FAIL_SEGFAULT
			return True
		
		# code means we have some results on pubchem
		# take the first one (maybe not the best practice!)
		chemspiID = results[0]
		
		try:
			result = cs.get_details(chemspiID, fields=['Mol3D'])['mol3D']
			lines = result.split('\n')
			lines = [l.strip() for l in lines]
			lines = lines[4: 4+job.elements.shape[0]]
		except:
			print("chemspider failed too...")
			job.state = StatusID.STATUS_GEOOPT_FAIL_NOCHEMSPIDER
			return True
		
		xyz = numpy.zeros((elements.shape[0],3))
		types = numpy.zeros(elements.shape[0], dtype=numpy.int32)
		
		for i in range(len(lines)):
			l = lines[i]
			w = l.split()
			for c in range(3): xyz[i,c] = float(w[c])
			types[i] = Element(w[3]).atomic_number
		
		idxs = numpy.argsort(types)
		xyz = xyz[idxs]
		
		# rewrite te input file with new xyz
		job.xyz = xyz
		job.WriteInput_geoopt()
		job.state = StatusID.STATUS_GEOOPT_SETUP
		return True
		
	
	# 2: if we are already using the chemspider geometry, try with MP2
	if job.geoopt_failstate == 1:
		# we already tried the pubchem trick... didnt work
		job.geoopt_failstate = 2
		job.timefactor *= 4
		# rewrite the input but using the MP2 relaxation
		job.WriteInput_geoopt()
	
	
	
	if job.geoopt_failstate == 2:
		# we also tried the MP2 but still segfault!
		# give up here!
		job.state = StatusID.STATUS_GEOOPT_FAIL
		return True
	
	
	return True
	
	
def Action_geoopt_fail_mem(job, taskID):
	
	print("calling Action_geoopt_fail_mem for task {} - {}".format(taskID, job))
	
	# if we already tried... itz cocked!
	if job.memory == 4000:
		job.state = StatusID.STATUS_GEOOPT_FAIL
		return True
	
	
	job.arch = "hsw"
	job.nthreads = 12
	job.memory = 4000
	job.WriteInput_geoopt()
	
	job.state = StatusID.STATUS_GEOOPT_SETUP
	return True






# DEPRECATED
def Action_geoopt_reload(job, taskID):
	
	# check if the run finished for real
	GeometryReload(jobID)
	
	# mark for resubmission
	job.state = StatusID.STATUS_STARTUP
	return True
	


## Prepare the CCSD calculation
#
def Action_geoopt_done(job, taskID):
	
	print("calling Action_geoopt_done for task {} - {}".format(taskID, job))

	if job.arch != "hsw":

		# some flags depends on how big the run is
		# set the cache level according to natoms
		if job.nelect <= 48: # small jobs
			job.time = 60
			job.nthreads = 20
			job.arch = "[skl]"
			job.memory = 2000

		elif job.nelect <= 100: # medium jobs
			job.cachelevel = 1
			job.arch = "[skl]"
			job.time = 60*2
			job.nthreads = 20
			job.memory = 2000

		elif job.nelect <= 130:
			job.cachelevel = 0
			job.arch = "[hsw]"
			job.time = 12*60
			job.nthreads = 24
			job.memory = 5000

		else: # larger jobs
			job.cachelevel = 0
			job.arch = "[hsw]"
			job.time = 24*60
			job.nthreads = 24
			job.memory = 5000


	job.WriteInput_ccsd()
	job.state = StatusID.STATUS_CCSD_SETUP

	# the task is finished, and this is the new status
	return True


def Action_ccsd_setup(job, taskID):
	
	print("calling Action_ccsd_setup for task {} - {}".format(taskID, job))
	
	# tweak the job file
	os.system("sed -i 's/SBATCH -J .*/SBATCH -J cc_{}/g' {}".format(taskID, job.job_cc))
	
	req = subprocess.run(['grep "QUEUE" {}'.format(job.job_cc)], shell=True, stdout=subprocess.PIPE)
	req = req.stdout.decode('utf-8')
	if len(req) > 5:
		job.state = StatusID.STATUS_GEOOPT_DONE
		return True


	# send to queue
	os.chdir(job.directory)
	os.system("rm slurm-*")
	os.system("sbatch ccsdrun.job")
	os.chdir("../")

	print("ACTION CCSD_SETUP: CCSD {} submitted".format(job.ID))
	job.state = StatusID.STATUS_CCSD_RUNNING
	return False


def Action_ccsd_running(job, taskID):
	
	print("calling Action_ccsd_running for task {} - {}".format(taskID, job))
	
	# check if this calculation is still running
	req = subprocess.run(["squeue -u federif1 | grep 'cc_{} '".format(taskID)], shell=True, stdout=subprocess.PIPE)
	req = req.stdout.decode('utf-8')
	
	if len(req) > 0:
		print("JOB {} still ccsd_running on TASK {}...".format(job.ID, taskID))
		return False

	job.state = StatusID.STATUS_CCSD_FINISHED
	return True


def Action_ccsd_finished(job, taskID):
	
	print("calling Action_ccsd_finished for task {} - {}".format(taskID, job))
	
	# check if the run finished for real
	slurmfile = glob.glob("{}/slurm-*.out".format(job.directory))
	
	# if there is no slurm output at all, try to rerun - should not happen!
	if len(slurmfile) == 0:
		print("CCSD {}: slurm output not found, something is wrong here".format(job.ID))
		job.state = StatusID.STATUS_CCSD_SETUP
		return True
	slurmfile = slurmfile[0]
	
	
	# no output file... maybe triton crashed? retry (might get stuck in infinite loop!
	# it is possible we got bad constraints?
	if not os.path.exists(job.output_cc):
		print("CCSD {}: output file not found?! retrying".format(job.ID))
		job.state = StatusID.STATUS_CCSD_SETUP
		return True
	
	# check if the final output is all there
	req = subprocess.run(['grep "Properties computed using the CC density matrix" {}'.format(job.output_cc)], shell=True, stdout=subprocess.PIPE)
	req = req.stdout.decode('utf-8')
	if len(req) > 5:
		print("CCSD {} completed.".format(job.ID))
		job.state = StatusID.STATUS_CCSD_DONE
		return True


	# code here => ccsd did not finish... find out why
	req = subprocess.run(['grep "error: Exceeded job memory limit" {}'.format(slurmfile)], shell=True, stdout=subprocess.PIPE)
	req = req = req.stdout.decode('utf-8')
	if len(req) > 5:
		print("CCSD {} ran out of memory.".format(job.ID))
		job.state = StatusID.STATUS_CCSD_FAIL_MEM
		return True
	
	
	# could have failed because of time limit?
	req = subprocess.run(['grep "DUE TO TIME LIMIT" {}'.format(slurmfile)], shell=True, stdout=subprocess.PIPE)
	req = req = req.stdout.decode('utf-8')
	if len(req) > 5:
		print("CCSD {} not enough time.".format(job.ID))
		if job.time * job.timefactor >= 5*24*60:
			print("CCSD {} already ran at max time".format(job.ID))
			job.state = StatusID.STATUS_CCSD_FAIL_TIME
			return True
		else:
			job.timefactor += 0.5
			job.state = StatusID.STATUS_GEOOPT_DONE
			return True
	
	job.state = StatusID.STATUS_CCSD_FAIL
	return True
	


def Action_ccsd_fail(job, taskID):
	
	slurmfile = glob.glob("{}/slurm-*.out".format(job.directory))
	slurmfile.sort()
	
	if len(slurmfile) == 0:
		print("CCSD {}: slurm output not found, something is wrong here".format(job.ID))
		job.state = StatusID.STATUS_CCSD_SETUP
		return True
	slurmfile = slurmfile[-1]
	
	req = subprocess.run(['grep "(error writing to file)" {}'.format(slurmfile)], shell=True, stdout=subprocess.PIPE)
	req = req = req.stdout.decode('utf-8')
	if len(req) > 5:
		if hasattr(job, "redo"):
			print("CCSD {}: job failed with weird failed to write message... and it was already retried! FAIL")
			job.state = StatusID.STATUS_CCSD_FAIL_UNKO
			return True
		else:
			print("CCSD {}: job failed with weird failed to write message... retry")
			job.redo = 1
			job.state = StatusID.STATUS_CCSD_SETUP
			return True
	
	
	req = subprocess.run(['grep "PSIO_ERROR: Failed to write" {}'.format(slurmfile)], shell=True, stdout=subprocess.PIPE)
	req = req = req.stdout.decode('utf-8')
	if len(req) > 5:
		print("CCSD {}: job failed with weird failed to write message... retry")
		job.state = StatusID.STATUS_CCSD_SETUP
		return True
	
	job.state = StatusID.STATUS_CCSD_FAIL_UNKO
	return True


# deprecated for now
def Action_ccsd_fail_mem(jobID, taskID):
	
	print("calling Action_ccsd_fail_mem for task {} job {}".format(taskID, jobID))
	
	moldir		= "molecule_{}".format(jobID)
	xyzfile 	= "{}/GEOM-B3LYP.xyz".format(moldir)
	molinput	= "{}/ccsdrun.py".format(moldir)
	outfile		= "{}/ccsdrun.out".format(moldir)
	jobfile		= "{}/ccsdrun.job".format(moldir)
	
	# if the job was in the hugemem queue, no point retrying
	req = subprocess.run(['grep "hugemem" {}'.format(jobfile)], shell=True, stdout=subprocess.PIPE)
	req = req = req.stdout.decode('utf-8')
	if len(req) > 2:
		print("CCSD {} out of memory on HUGEMEM, quitting.".format(jobID))
		return True, StatusID.STATUS_CCSD_FAIL
	
	# else, try to send it to the hugemem queue
	# ...
	os.system("sed -i 's/SBATCH -p .*/SBATCH -p hugemem/g' {}".format(jobfile))
	os.system("sed -i 's/#SBATCH --constraint.*/### no constraint ###/g' {}".format(jobfile))
	
	return True, StatusID.STATUS_CCSD_SETUP
	


def Action_ccsd_done(job, taskID):
	
	print("calling Action_ccsd_done for task {} - {}".format(taskID, job))
	
	# copy the file
	os.system("mv {} completed/.".format(job.directory))
	job.state = StatusID.STATUS_COMPLETED
	return True
	

