import numpy as np
from matplotlib import pyplot as plt
import subprocess as sp

NTRAJ   = 101  # USER INPUT
dtI     = 1.0  # USER INPUT
NSTEPS  = 1000  # USER INPUT


#### DO NOT MODIFY BELOW THIS POINT ####
TRAJ_DIRS = [ f"TRAJ/traj-{j}/MD_OUTPUT/" for j in range(NTRAJ) ]
DATA_DIR = "PLOTS_DATA"
sp.call("mkdir -p PLOTS_DATA",shell=True)

NSTATES = 0
#NSTEPS  = 0
TIME = []
POP  = []
good_traj = 0
for ind, dir in enumerate(TRAJ_DIRS):
    try:
        tmp = np.loadtxt(f"{dir}/Population.dat")[:NSTEPS]
    except OSError:
        print(f"No good: {dir}")
        continue
    if ( ind == 0 ): # Let's assume the first trajectory is finished and okay.
        #NSTEPS  = len(tmp)
        NSTATES = len(tmp[0])-2 # Not time-step or total population
        POP = np.zeros( (NSTEPS,NSTATES) )
        TIME[:] = tmp[:,0] * dtI
    if ( len(tmp) < NSTEPS ):
        continue
    POP[:,:] += tmp[:,1:-1]
    good_traj += 1

POP /= good_traj
print("Number of good (or finished) trajectories:", good_traj, "of", NTRAJ)

# Save the average population to a file
OUTPOP = np.zeros(( NSTEPS, NSTATES+2 ))
OUTPOP[:,0] = TIME
for state in range( NSTATES ):
    OUTPOP[:,state+1] = POP[:,state]
OUTPOP[:,-1] = np.sum(POP[:,:],axis=-1)
np.savetxt( f"{DATA_DIR}/Population_average-{good_traj}.dat", OUTPOP, fmt="%2.5f" )

# Make simple plot

for state in range( NSTATES ):
    plt.plot( TIME, POP[:,state], lw=3, label=f"S{state}" )
plt.plot( TIME, np.sum(POP[:,:],axis=-1), c="black", alpha=0.5, lw=3, label=f"Total" )

plt.legend()
plt.xlim(0,TIME[-1])
plt.ylim(0,1.005)
plt.xlabel("Time (fs)", fontsize=15)
plt.ylabel("Population", fontsize=15)
plt.title(f"# of Trajectories: {good_traj}", fontsize=15)
plt.savefig(f"{DATA_DIR}/Population_average-{good_traj}.jpg",dpi=600)


