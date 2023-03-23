import numpy as np
from matplotlib import pyplot as plt

# USER INPUT
TRACK_ATOMS = [31,22,23] # 1,2,3,4,...



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
print(f"{TYPES[TRACK_ATOMS[0]-1]}({TRACK_ATOMS[0]}) \
        {TYPES[TRACK_ATOMS[1]-1]}({TRACK_ATOMS[1]}) \
        {TYPES[TRACK_ATOMS[2]-1]}({TRACK_ATOMS[2]})")

# Extract important atoms for simplicity
A_ATOMS        = np.zeros(( NSTEPS, 3, 3 ))
A_ATOMS[:,0,:] =  GEOMS[:,TRACK_ATOMS[0]-1,:]
A_ATOMS[:,1,:] =  GEOMS[:,TRACK_ATOMS[1]-1,:]
A_ATOMS[:,2,:] =  GEOMS[:,TRACK_ATOMS[2]-1,:]

# Calculate the angle given by the four atoms
ANGLE = np.zeros(( NSTEPS ))
for step in range( NSTEPS ):
    b10 = A_ATOMS[step,1,:] - A_ATOMS[step,0,:]
    b21 = A_ATOMS[step,2,:] - A_ATOMS[step,1,:]

    norm10 = np.linalg.norm(b10)
    norm21 = np.linalg.norm(b21)

    cosPHI = np.dot(b10,b21) / norm10 / norm21

    ANGLE[step] = 180 - np.degrees(np.arccos(cosPHI))

# Save data to file
np.savetxt( f"ANGLE_{TRACK_ATOMS[0]}-{TRACK_ATOMS[1]}-{TRACK_ATOMS[2]}.dat", ANGLE , header=f"" )

# Make a plot of this dehedral
dtI = 1.0 # fs
plt.plot( np.arange(NSTEPS)*dtI, ANGLE )
plt.xlim(0,(NSTEPS-1)*dtI)
plt.xlabel("Time (fs)",fontsize=15)
plt.ylabel(f"{TYPES[TRACK_ATOMS[0]-1]}-\
{TYPES[TRACK_ATOMS[1]-1]}-\
{TYPES[TRACK_ATOMS[2]-1]} Angle ($^o$)",fontsize=15)
plt.savefig(f"ANGLE_{TRACK_ATOMS[0]}-{TRACK_ATOMS[1]}-{TRACK_ATOMS[2]}.jpg",dpi=600)


