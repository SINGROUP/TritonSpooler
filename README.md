# TritonSpooler
A simple and flexible job spooler for Triton (or any other cluster/supercomputer)


The spooler is a python object that manages a set of jobs.
A job is a python object with an ID and a state, that represents a set of operations to be completed. For example, in the Psi4 application a job consists of optimizing a molecule geometry and then run CCSD. In reality the job is split into many more small operations.
The spooler must be provided with a dictionary of state-function pairs in order to know what task to perform for each job, depending on their states. If a state does not appear in the dictionary, the spooler will consider jobs in that state completed.
When all jobs are in such states, the spooler will stop.


## Folder content
. asd
. lol
. asd
