import numpy as np
import subprocess as sp
import os

##### Braden M. Weight, January 13, 2023 #####

# Run this in the MD_OUTPUT folder of an NVT job.
# In the get_Globals() funciton, define some parameters.
# Then one should move the new TRAJ folder to another location 
#     for better organization.

def get_Globals():
    global N_EQUIL, N_SAVE_EVERY
    N_EQUIL = 40 # Skip the first {N_EQUIL} step for equilibration
    N_SAVE_EVERY = 2 # Save every {N_SAVE_EVERY} step after {N_EQUIL} steps

def read_XYZ(file01):
    XYZ_File = open(file01,"r").readlines()
    NAtoms = int(XYZ_File[0])
    NSteps = len(XYZ_File) // (NAtoms+2)
    print("NSteps =", NSteps)
    XYZ_DATA   = np.zeros(( NSteps, NAtoms, 3 ))
    XYZ_LABELS = []
    for step in range( NSteps ):
        for at in range( NAtoms+2 ):
            line_ind = step * (NAtoms+2) + at
            #print(line_ind)
            t = XYZ_File[line_ind].split()
            if ( at not in [0,1] and len(t) == 4 ):
                if (step == 0): XYZ_LABELS.append( t[0] )
                XYZ_DATA[step,at-2,:] = np.array([ float(t[1]), float(t[2]), float(t[3]) ]) # Ang
                


    return XYZ_LABELS, XYZ_DATA

def write_XYZ(file01,XYZ_LABELS,XYZ_DATA):
    NAtoms = len(XYZ_LABELS)
    f = open(file01,"w")
    f.write(f"{NAtoms}\n")
    f.write(f"_____\n")
    for at in range( NAtoms ):
        f.write("%s  %2.6f  %2.6f  %2.6f\n" % (XYZ_LABELS[at],XYZ_DATA[at,0],XYZ_DATA[at,1],XYZ_DATA[at,2]))
    f.close()

def make_initital_conditions(XYZ_LABELS,GEOMS,VELOC):
    sp.call("mkdir TRAJ",shell=True)
    NSteps = len(GEOMS)
    NAtoms = len(GEOMS[0])
    NTraj = (NSteps - N_EQUIL - 1) / N_SAVE_EVERY
    print(f"There will be {NTraj} trajectories, rounded down.")
    
    GEOMS_SAVE = np.zeros(( int(np.floor(NTraj)), NAtoms, 3 ))    
    VELOC_SAVE = np.zeros(( int(np.floor(NTraj)), NAtoms, 3 ))    

    count = 0
    for step in range( N_EQUIL+1, NSteps ):
        if ( step % N_SAVE_EVERY == 0 ):
            if ( os.path.exists(f"TRAJ/traj-{count}") ):
                sp.call(f"rm -r TRAJ/traj-{count}",shell=True)
            sp.call(f"mkdir TRAJ/traj-{count}",shell=True)
            write_XYZ(f"TRAJ/traj-{count}/geometry_input.xyz", XYZ_LABELS, GEOMS[step])
            write_XYZ(f"TRAJ/traj-{count}/velocity_input.xyz", XYZ_LABELS, VELOC[step])
            count += 1

def main():
    get_Globals()
    XYZ_LABELS, GEOMS = read_XYZ("trajectory.xyz")
    XYZ_LABELS, VELOC = read_XYZ("velocity.xyz")

    make_initital_conditions(XYZ_LABELS,GEOMS,VELOC)




if ( __name__ == "__main__" ):
    main()