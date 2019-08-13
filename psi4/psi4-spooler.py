import numpy
from spooler import Spooler
import psi4utils as utils
import pickle, os, glob



# create a spooler
spl = Spooler()

# configure it
utils.SetSpooler(spl)

# run the spooler cycle
spl.Spool()

