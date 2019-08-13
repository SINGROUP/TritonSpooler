import os
from shutil import copyfile
import math
import random
import numpy
import fnmatch
import sys
import glob
import time
import json
import subprocess
import pickle


'''


each job in the spooler has a status code (Enum, user given)

job state can be associated to ACTION or CHECK functions

When a job state has an associated ACTION, the spooler will take a free taskID, attach the action to it,
and perform the user defined action.
The action function then returns the new state of the job, if it changed, and a suggested waiting time.

When a job state has an associated CHECK, ...? what is the difference

for each possible state the user has to supply an action that the spooler will execute
if there is None action for a state, the spooler will do nothing and consider that job concluded somehow

at each iteration the spooler looks at the jobs and performs an action

(job, status) -> action



'''


## Class that represents a general spooler job.
#
# This is meant to be abstract since it provides nothing more than the basics.
# 
#
# Derive from this class to implement custom behaviour.
#
class SpoolerJob(object):
	

	## Default constructor needs to know the jobID and the state
	def __init__(self, jobID, state):
	
		## Unique identifier for this job.
		# There cannot be two jobs with same ID in the spooler.		
		self.ID = jobID

		## State of the job. 
		# Can be anything equatable value type (number or string).
		self.state = state
		
	
	## Convenient printing function.
	def __repr__(self):
		return "<job {} - state: {}>".format(self.ID, self.state)



### Spooler object.
class Spooler(object):
	
	## Constructor.
	def __init__(self):
		
		## Maximum amount of concurrent active job tasks.
		self.maxTasks = 40

		## Spooler job check cycle period.
		self.cycleTime = 30

		## Index of the current checkpoint file.
		self.fileindex = 0


		## Dictionary of job objects.
		# key: jobID (string)
		# val: job object (SpoolerJob)
		self.jobs = {}
		
		## Dictionary of currently active tasks.
		# key: taskID (int)
		# val: jobID of tasks sent to the queue
		self.tasks = {}
		
		## Available taskIDs.
		self.taskIDFree = set(range(self.maxTasks))
		
		
		## Dictionary of actions to perform for any job state
		# User provided!
		# key: job status descriptor (up to the user)
		# val: function(jobobj, taskID)
		self.actions = {}
		
		
		## Load the status of the spooler from spooler.status file.
		self.LoadStatus()
		
		
		
		
	## Loads the status of the spooler from spooler.status
	# If the file does not exist, it will create an empty one.
	#
	def LoadStatus(self):
		
		if os.path.exists("spooler.status"):
			
			print("SPOOLER: loading status data.")
			
			fbin = open("spooler.status", "rb")
			tmp = pickle.load(fbin)
			self.jobs = tmp[0]
			self.tasks = tmp[1]
			fbin.close()

		else:
			print("SPOOLER: status data missing, making new.")
			self.SaveStatus()
			
		
		return



	## Saves the current spooler status in a checkpoint file.
	# Each time this executes, the ĺatest spooler.status is copied to spooler.X.status,
	# to preserve some history in case of corruption due to Triton failures.
	#
	def SaveStatus(self):

		print("SPOOLER: saving status data.")
		# make a copy of the status file
		copyfile("spooler.status", "spooler.{}.status".format(self.fileindex))
		self.fileindex += 1
		if self.fileindex == 5:
			self.fileindex = 0

		fbin = open("spooler.status", "wb")
		pickle.dump([self.jobs, self.tasks], fbin)
		fbin.close()

		return
	
	
	## Loads the resupply file if one is ready.
	# The file is deleted afterwards.
	#
	def Resupply(self):
		
		if os.path.exists("spooler.resupply.lock"):
			fbin = open("spooler.resupply.data", "rb")
			newjobs = pickle.load(fbin)
			fbin.close()
			
			print('resupply new jobs: ', newjobs)
			
			self.jobs = {**newjobs, **self.jobs}
			os.system("rm spooler.resupply.*")
		
		return

	
	## Loads the restore file if one is ready.
	# The file is deleted afterwards.
	#
	def Restore(self):
		
		if os.path.exists("spooler.restore.lock"):
			fbin = open("spooler.restore.data", "rb")
			newjobs = pickle.load(fbin)
			fbin.close()

			for k in newjobs.keys():
				if k in self.jobs.keys():
					self.jobs[k].state = newjobs[k]
			
			os.system("rm spooler.restore.data spooler.restore.lock")

		return

	
	## Get a list of jobs that need to be done.
	# 
	# 
	# 	
	def GetJobsTodo(self):
		
		todos = []
		for jobID in self.jobs.keys():
			job = self.jobs[jobID]
			hasAction = job.state in self.actions.keys()
			#hasCheck = state in self.checks.keys()
			hasRunningTask = jobID in list(self.tasks.values())
			
			if (hasAction) and (not hasRunningTask):
				todos.append(jobID)
		
		#print("GJTD",todos)
		return list(set(todos))
	


	def GetJobsByState(self, state):

		result = []
		for k in self.jobs.keys():
			if self.jobs[k].state == state:
				result.append(self.jobs[k])

		return result

	
	
	## Main spooler cycle.
	#
	def Spool(self, **kwargs):
		
		
		
		# give an initial check to the jobs that are supposed to be running...
		for taskID in self.tasks.keys():
			self.taskIDFree.remove(taskID)
		
		# resupply just in case
		self.Resupply(); self.Restore()
		
		looper = True
		depletion = False

		if kwargs is not None:
			for k in kwargs.keys():
				v = kwargs[k]
				if k == "depletion" and v == True:
					looper = False
					depletion = True

		print("spooler starting: jobs: ", self.jobs)

		while True:
			
			todo = self.GetJobsTodo() # ??? how is this done?
			print("SPOOLER: iteration info")
			print("\tJobs to do [{}]: {}".format(len(todo),todo))
			print("\tFree slots: {}".format(len(self.taskIDFree)))
			print("\tActive tasks: [{}] {}".format(len(self.tasks.keys()), str(self.tasks)))
			
			if len(self.tasks) == 0 and len(todo) == 0:
				print("SPOOLER: all jobs finished, nothing more to do.")
				break

			if depletion and len(todo) < self.maxTasks:
				print("SPOOLER: pool depleting... stopping the loop.")
				break



			# job to do and free task slot
			if len(todo) > 0 and len(self.taskIDFree) > 0:
				print("SPOOLER: time to start a new job...")
				lst = todo
				jobID = lst.pop()
				job = self.jobs[jobID]
				state = job.state
				
				if state in self.actions.keys():
					
					taskID = self.taskIDFree.pop()
					self.tasks[taskID] = jobID
					print("SPOOLER: chosen job {} <= task {}".format(jobID, taskID))
					finished = self.actions[state](job, taskID)
					
					if finished:
					
						print("SPOOLER: task {} finished immediately.".format(taskID))
						del self.tasks[taskID]
						self.taskIDFree.add(taskID)
					
					else:
						# job did not end immediately... put it on a checklist - maybe not needed
						print("SPOOLER: task {} is on the checklist now.".format(taskID))
						#self.tasksToCheck[taskID] = jobID
					
					print("SPOOLER: post-action {}".format(job, self.jobs[jobID]))
					self.SaveStatus()
					continue
						
			
			# check the currently running tasks and perform their actions
			removers = []
			for taskID in self.tasks.keys():
				jobID = self.tasks[taskID]
				job = self.jobs[jobID]
				finished = self.actions[job.state](job, taskID)
				
				if finished:
					print("SPOOLER: running task {} ({}) ended".format(taskID, job))
					removers.append(taskID)
					
			
			for taskID in removers:
				del self.tasks[taskID]
				self.taskIDFree.add(taskID)
				self.SaveStatus()
				
			if len(removers) > 0:
				continue

			time.sleep(self.cycleTime)
			self.Resupply()
			self.Restore()
			
		
		return






