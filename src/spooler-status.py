import numpy
from spooler import Spooler
import pickle, os, glob

# Create a spooler to load the status file
spl = Spooler()

print("active tasks:")
keys = list(spl.tasks.keys())
keys.sort()
for t in keys: 
        print(t,spl.tasks[t], spl.jobs[spl.tasks[t]].state)


print("Active jobs:",len(spl.tasks))
