import numpy as np
from matplotlib import pyplot as plt

NTRAJ   = 400  # USER INPUT
dtI     = 0.25 # USER INPUT # TODO Get this automatically


TRAJ_DIRS = [ f"TRAJ/traj-{j}/MD_OUTPUT/" for j in range(NTRAJ) ]

NSTATES = 0
NSTEPS  = 0
TIME = []
POP  = []
good_traj = 0
for ind, dir in enumerate(TRAJ_DIRS):
    tmp = np.loadtxt(f"{dir}/Population.dat")
    if ( ind == 0 ): # Let's assume the first trajectory is finished and okay.
        NSTEPS  = len(tmp)
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
np.savetxt( f"Population_average-{good_traj}.dat", OUTPOP, fmt="%2.5f" )

# Make simple plot

#for state in range( NSTATES ):
#    plt.plot( TIME, POP[:,state], lw=3, label=f"S{state}" )
#plt.plot( TIME, np.sum(POP[:,:],axis=-1), c="black", alpha=0.5, lw=3, label=f"Total" )

# Plot only S2
plt.plot( TIME, POP[:,2], lw=3, label=f"S2" )

plt.legend()
plt.xlim(0,TIME[-1])
plt.ylim(0,1)
plt.savefig(f"Population_average-{good_traj}.jpg",dpi=600)


