import numpy as np
from matplotlib import pyplot as plt

# USER INPUT
TRACK_ATOMS = [3,7,23,25] # 1,2,3,4,...



#### DO NOT MODIFY BELOW THIS POINT ####
# Read all coordinates
lines = open("trajectory.xyz","r").readlines()
NATOMS = 0
NSTEPS = 0
GEOMS  = 0
TYPES  = []
step   = 0
for count, line in enumerate( lines ):
    if ( count == 0 ):
        NATOMS = int(line)
        NSTEPS = len(lines) // (NATOMS+2)
        GEOMS  = np.zeros(( NSTEPS, NATOMS, 3 ))
    if ( (count) % (NATOMS+2) == 0 ):
        for at in range( NATOMS ):
            if ( step == 0 ):
                TYPES.append( lines[count+2+at].split()[0] )
            GEOMS[step,at,:] = np.array( lines[count+2+at].split()[1:] ).astype(float)
        step += 1

# Say which atom types (and labels) we are trying to track
print("\nFinding time-dependent dihedral angle for the following atom types (labels):")
print(f"{TYPES[TRACK_ATOMS[0]]}({TRACK_ATOMS[0]}) \
        {TYPES[TRACK_ATOMS[1]]}({TRACK_ATOMS[1]}) \
        {TYPES[TRACK_ATOMS[2]]}({TRACK_ATOMS[2]}) \
        {TYPES[TRACK_ATOMS[3]]}({TRACK_ATOMS[3]})")

# Extract important atoms for simplicity
D_ATOMS        = np.zeros(( NSTEPS, 4, 3 ))
D_ATOMS[:,0,:] =  GEOMS[:,TRACK_ATOMS[0]-1,:]
D_ATOMS[:,1,:] =  GEOMS[:,TRACK_ATOMS[1]-1,:]
D_ATOMS[:,2,:] =  GEOMS[:,TRACK_ATOMS[2]-1,:]
D_ATOMS[:,3,:] =  GEOMS[:,TRACK_ATOMS[3]-1,:]

# Calculate the dihedral given by the four atoms
DIH = np.zeros(( NSTEPS ))
for step in range( NSTEPS ):
    b0 = D_ATOMS[step,1,:] - D_ATOMS[step,0,:]
    b1 = D_ATOMS[step,2,:] - D_ATOMS[step,1,:]
    b2 = D_ATOMS[step,3,:] - D_ATOMS[step,2,:]

    b0crossb1 = np.cross(b0,b1)
    b1crossb2 = np.cross(b1,b2)

    norm01 = np.linalg.norm(b0crossb1)
    norm12 = np.linalg.norm(b1crossb2)

    cosPHI = np.dot(b0crossb1,b1crossb2) / norm01 / norm12

    DIH[step] = np.degrees(np.arccos(cosPHI))

# Save data to file
np.savetxt( f"DIH_{TRACK_ATOMS[0]}-{TRACK_ATOMS[1]}-{TRACK_ATOMS[2]}-{TRACK_ATOMS[3]}.dat", DIH , header=f"" )

# Make a plot of this dehedral
dtI = 1.0 # fs
plt.plot( np.arange(NSTEPS)*dtI, DIH )
plt.xlim(0,(NSTEPS-1)*dtI)
plt.xlabel("Time (fs)",fontsize=15)
plt.ylabel(f"{TYPES[TRACK_ATOMS[0]]}-\
{TYPES[TRACK_ATOMS[1]]}-\
{TYPES[TRACK_ATOMS[2]]}-\
{TYPES[TRACK_ATOMS[3]]} Dihedral Angle ($^o$)",fontsize=15)
plt.savefig(f"DIH_{TRACK_ATOMS[0]}-{TRACK_ATOMS[1]}-{TRACK_ATOMS[2]}-{TRACK_ATOMS[3]}.jpg",dpi=600)


