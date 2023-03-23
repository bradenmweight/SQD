import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
import subprocess as sp

# This script will search for unfinished trajectories.
# WARNING: Only run this script when all jobs have 
#          finished/failed in this directory.
# To change submission details for failed jobs, 
#          change ./submit.SQD

# SYNTAX: python3 track_angle_allTRAJ.py

# USER INPUT
TRACK_ATOMS = [31,22,23] # 1,2,3,4,...
NTRAJ       = 201
dtI         = 1.0 # fs


#### DO NOT MODIFY BELOW THIS POINT ####
DIRS = [f"TRAJ/traj-{traj}/MD_OUTPUT/" for traj in range(NTRAJ)]
DATA_DIR = "PLOTS_DATA"
sp.call("mkdir -p PLOTS_DATA",shell=True)

# Read all coordinates
NATOMS   = 0
NSTEPS   = 0
GEOMS    = 0
TYPES    = []
ANGLE      = []
bad_traj = 0

for traj, path in enumerate(DIRS):
    lines = open(f"{path}/trajectory.xyz","r").readlines()
    if ( traj >= 1 and len(lines) < NSTEPS*(NATOMS+2) ):
        ANGLE = np.delete(DIH,traj,axis=0)
        #print( traj, len(lines), NSTEPS*(NATOMS+2), ANGLE.shape )
        print(f"Traj = {traj} did not finish properly.")
        bad_traj += 1
        continue

    step  = 0
    for count, line in enumerate( lines ):
        if ( count == 0 and traj == 0 ):
            NATOMS = int(line)
            NSTEPS = len(lines) // (NATOMS+2)
            GEOMS  = np.zeros(( NTRAJ, NSTEPS, NATOMS, 3 ))
            ANGLE    = np.zeros(( NTRAJ, NSTEPS ))
        if ( (count) % (NATOMS+2) == 0 ):
            for at in range( NATOMS ):
                if ( step == 0 and traj == 0 ):
                    TYPES.append( lines[count+2+at].split()[0] )
                GEOMS[traj-bad_traj,step,at,:] = np.array( lines[count+2+at].split()[1:] ).astype(float)
            step += 1

    if ( traj == 0 ):
        # Say which atom types (and labels) we are trying to track
        print("\nFinding time-dependent dihedral angle for the following atom types (labels):")
        print(f"\
{TYPES[TRACK_ATOMS[0]-1]}({TRACK_ATOMS[0]}) \
{TYPES[TRACK_ATOMS[1]-1]}({TRACK_ATOMS[1]}) \
{TYPES[TRACK_ATOMS[2]-1]}({TRACK_ATOMS[2]})")

    # Extract important atoms for simplicity
    A_ATOMS        = np.zeros(( NSTEPS, 4, 3 ))
    A_ATOMS[:,0,:] =  GEOMS[traj-bad_traj,:,TRACK_ATOMS[0]-1,:]
    A_ATOMS[:,1,:] =  GEOMS[traj-bad_traj,:,TRACK_ATOMS[1]-1,:]
    A_ATOMS[:,2,:] =  GEOMS[traj-bad_traj,:,TRACK_ATOMS[2]-1,:]

    # Calculate the dihedral given by the four atoms
    for step in range( len(ANGLE[0,:]) ):
        b10 = A_ATOMS[step,1,:] - A_ATOMS[step,0,:]
        b21 = A_ATOMS[step,2,:] - A_ATOMS[step,1,:]

        norm10 = np.linalg.norm(b10)
        norm21 = np.linalg.norm(b21)

        cosPHI = np.dot(b10,b21) / norm10 / norm21

        ANGLE[traj,step] = 180 - np.degrees(np.arccos(cosPHI))

# Save average (and error in) dihedral to file
ANG_AVE = np.average(ANGLE[:,:],axis=0)
ANG_STD = np.std(ANGLE[:,:],axis=0)
np.savetxt(f"{DATA_DIR}/ANGLE_AVE_STD_{TRACK_ATOMS[0]}_{TRACK_ATOMS[1]}_{TRACK_ATOMS[2]}.dat", np.c_[ANG_AVE, ANG_STD] )

# FIND SPECIFIC TRAJECTORIES WITH MAX AND MIN ANGLES FOR VISUALIZATION
traj_min = [-1,1000]
traj_max = [-1,-1000]
for traj in range(NTRAJ-bad_traj):
    tmp1 = np.min(ANGLE[traj,:])
    tmp2 = np.max(ANGLE[traj,:])
    if ( tmp1 < traj_min[1] ):
        traj_min[0] = traj
        traj_min[1] = tmp1
    if ( tmp2 > traj_max[1] ):
        traj_max[0] = traj
        traj_max[1] = tmp2

print( f"Minimum dihedral angle trajectory label (DIH) = {traj_min[0]+1} ({round(traj_min[1],2)})" )
print( f"Maximum dihedral angle trajectory label (DIH) = {traj_max[0]+1} ({round(traj_max[1],2)})" )

# Make a plot of the average dihedral for all trajectories
for traj in range( NTRAJ-bad_traj ):
    plt.plot( np.arange(NSTEPS)*dtI, ANGLE[traj,:], c='black', lw=2, alpha=0.2 )
plt.plot( np.arange(NSTEPS)*dtI, ANG_AVE, c='red', lw=3, alpha=0.8 )
plt.xlim(0,(NSTEPS-1)*dtI)
plt.xlabel("Time (fs)",fontsize=15)
plt.ylabel(f"\
{TYPES[TRACK_ATOMS[0]-1]}-\
{TYPES[TRACK_ATOMS[1]-1]}-\
{TYPES[TRACK_ATOMS[2]-1]} Angle ($^o$)",fontsize=15)
plt.savefig(f"{DATA_DIR}/ANGLE_allTRAJ_{TRACK_ATOMS[0]}_{TRACK_ATOMS[1]}_{TRACK_ATOMS[2]}.jpg",dpi=600)
plt.clf()


# Make a density plot of the dihedral
SIG  = 5 # Gaussian width (degrees)
MIN  = 0 # Degrees
MAX  = 180 # Degrees
NPTS = 100
cmap = "binary" # "hot_r"
PLOT_MAX = 1.0

AGRID = np.linspace( MIN,MAX,NPTS  )
SPEC  = np.zeros(( NPTS, NSTEPS  ))

for pt in range(NPTS):
    SPEC[pt,:] += np.sum( np.exp( -( AGRID[pt] - ANGLE[:,:] )**2 / 2 / SIG**2 ), axis=0 )

# Normalize to MAX(SPEC) = 1.0
SPEC /= np.max(SPEC) # /= NTRAJ-bad_traj

# Save 2D data
np.savetxt(f"{DATA_DIR}/ANGLE_density_{TRACK_ATOMS[0]}_{TRACK_ATOMS[1]}_{TRACK_ATOMS[2]}.dat", SPEC)

# Remove small values for better coloring
#MIN_THRESHOLD = 0.1
#SPEC[ SPEC < MIN_THRESHOLD ] = 0

plt.contourf( np.arange(NSTEPS)*dtI, AGRID, SPEC[:,:], cmap=cmap, levels=100, vmax=PLOT_MAX )

norm = mpl.colors.Normalize(vmin=0, vmax=PLOT_MAX)
plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),pad=0.01)

plt.xlabel("Time (fs)",fontsize=15)
plt.ylabel(f"\
{TYPES[TRACK_ATOMS[0]-1]}-\
{TYPES[TRACK_ATOMS[1]-1]}-\
{TYPES[TRACK_ATOMS[2]-1]} Angle ($^o$)",fontsize=15)
plt.savefig(f"{DATA_DIR}/ANGLE_density_{TRACK_ATOMS[0]}_{TRACK_ATOMS[1]}_{TRACK_ATOMS[2]}.jpg",dpi=600)
plt.clf()

