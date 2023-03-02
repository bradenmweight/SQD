import numpy as np
from matplotlib import pyplot as plt
import subprocess as sp
import os

NTRAJ   = 201  # USER INPUT
dtI     = 0.1 # USER INPUT
NSTEPS  = 1500 # USER INPUT

####### KEEP FALSE UNLESS ABSOLUTELY SURE #######
resubmit_bad = False
##################################################

TRAJ_DIRS = [ f"TRAJ/traj-{j}/MD_OUTPUT/" for j in range(NTRAJ) ]

NSTATES    = 0
TIME       = []
POP_COEFF  = []
POP_AS     = []
good_traj  = 0
for ind, dir in enumerate(TRAJ_DIRS):
    try:
        tmp = np.loadtxt(f"{dir}/Population.dat")[:NSTEPS]
    except OSError:
        print(f"No good: {dir}")
        if ( resubmit_bad == True ):
            os.chdir(f"{'/'.join(dir.split('/')[:-2])}")
            sp.call("sbatch submit.SQD",shell=True)
            os.chdir("../../")
        continue
    if ( ind == 0 ): # Let's assume the first trajectory is finished and okay.
        NSTATES   = len(tmp[0])-3 # Not time-step or total population or AS
        POP_COEFF = np.zeros( (NSTEPS,NSTATES) )
        POP_AS    = np.zeros( (NSTEPS,NSTATES) )
        TIME[:]   = tmp[:,0] * dtI
    if ( len(tmp) < NSTEPS ):
        continue
    POP_COEFF[:,:] += tmp[:,2:-1]
    AS_t = tmp[:,1].astype(int)
    for step in range( NSTEPS ):
        POP_AS[step,AS_t[step]] += 1
    good_traj += 1

POP_COEFF /= good_traj
POP_AS  /= good_traj
print("Number of good (or finished) trajectories:", good_traj, "of", NTRAJ)

# Save the average population to a file
OUTPOP = np.zeros(( NSTEPS, NSTATES+3 ))
OUTPOP[:,0] = TIME
for state in range( NSTATES ):
    OUTPOP[:,state+1] = POP_COEFF[:,state]
OUTPOP[:,-1] = np.sum(POP_COEFF[:,:],axis=-1)
np.savetxt( f"Population_average-{good_traj}.dat", OUTPOP, fmt="%2.5f" )

# Make simple plot

#for state in range( NSTATES ):
#    plt.plot( TIME, POP_COEFF[:,state], lw=3, label=f"S{state}" )
#plt.plot( TIME, np.sum(POP_COEFF[:,:],axis=-1), c="black", alpha=0.5, lw=3, label=f"Total" )
plt.plot( TIME, POP_COEFF[:,2], c='orange', lw=3, label="<COEFF>" )
plt.plot( TIME, POP_AS[:,2], c='orange', lw=3, label=f"<AS>" )

plt.legend()
plt.xlim(0,TIME[-1])
plt.ylim(0,1.01 )
plt.savefig(f"Population_average-{good_traj}.jpg",dpi=600)


