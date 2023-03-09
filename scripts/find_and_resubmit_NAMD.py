import numpy as np
import os
import subprocess as sp

# WARNING: Do not use this script if some jobs are still running
# This script assumes TRAJ/traj-0 finished correctly.

NTRAJ         = 201 # Look in the first {NTRAJ} directories
RESUBMIT_FLAG = True # Whether to resubmit the unfinished jobs or not

#############################################
DIRS       = [f"TRAJ/traj-{traj}/" for traj in range(NTRAJ)]
NSTEPS_REF = 0 # TRAJ/traj-0 used as reference
NSTEPS     = 0
count_good = 0

def resubmit():
    if ( RESUBMIT_FLAG == True ): 
        print(f"\tResubmitted: TRAJ/traj-{traj}")
        sp.call("sbatch submit.SQD",shell=True)
    else:
        print(f"\tFailed: TRAJ/traj-{traj}")

for traj, path in enumerate(DIRS):
    os.chdir(path)
    #print("/".join(os.getcwd().split("/")[-2:]))
    if ( os.path.exists("MD_OUTPUT/") ):
        if ( os.path.isfile("MD_OUTPUT/Population.dat") ):
            POP    = np.loadtxt("MD_OUTPUT/Population.dat")
            NSTEPS = len(POP)
            if ( traj == 0 ):
                NSTEPS_REF  = NSTEPS
                count_good += 1
                os.chdir("../../")
                continue
            else:
                if ( NSTEPS == NSTEPS_REF ):
                    count_good += 1
                    os.chdir("../../")
                    continue
                else:
                    resubmit()
                    os.chdir("../../")
        else:
            resubmit()
            os.chdir("../../")
    else:
        resubmit()
        os.chdir("../../")

print(f"There were {count_good} good trajectories out of {NTRAJ}.")