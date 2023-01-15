import numpy as np
from random import gauss
import properties

import Eh
import spinLSC

def get_Force(DYN_PROPERTIES):
    if ( DYN_PROPERTIES["NAMD_METHOD"] == "EH" ):
        return Eh.get_Force(DYN_PROPERTIES)
    elif ( DYN_PROPERTIES["NAMD_METHOD"] == "SPINLSC" ):
        return spinLSC.get_Force(DYN_PROPERTIES)
    else:
        print("No NAMD_METHOD found.")
        exit()


# def get_random_XStep( DYN_PROPERTIES ):
#     LANGEVIN_LAMBDA = DYN_PROPERTIES["LANGEVIN_LAMBDA"] / 1000 / 27.2114 # meV --> a.u. # TODO Change units in read_input.py
#     TEMP  = DYN_PROPERTIES["TEMP"] * (0.025 / 300) / 27.2114 # K -> KT (a.u.) # TODO Change units in read_input.py
#     NAtoms  = DYN_PROPERTIES["NAtoms"]
#     masses = np.array([ np.array([m,m,m]) for m in DYN_PROPERTIES["MASSES"] ]) # TODO Change the shape in read_input.py
#     dtI     = DYN_PROPERTIES["dtI"]


    
#     return DYN_PROPERTIES, R_RAND

# def get_damped_XStep( DYN_PROPERTIES ):
#     V       = DYN_PROPERTIES["Atom_velocs_new"]
#     LANGEVIN_LAMBDA = DYN_PROPERTIES["LANGEVIN_LAMBDA"] / 1000 / 27.2114 # meV --> a.u. # TODO Change units in read_input.py
#     masses = np.array([ np.array([m,m,m]) for m in DYN_PROPERTIES["MASSES"] ])
#     F_DAMP          = -1.0 * LANGEVIN_LAMBDA * DYN_PROPERTIES["Atom_velocs_new"] * masses
#     cL = 2*M*dtI/(2*M+dtI*LANGEVIN_LAMBDA) # Langevin constant
#     R_DAMP += cL * (1-dtI/cL) * V # updated nuclear position
#     return DYN_PROPERTIES, R_DAMP

# def get_damped_VStep( DYN_PROPERTIES ):
#     V       = DYN_PROPERTIES["Atom_velocs_new"]
#     LANGEVIN_LAMBDA = DYN_PROPERTIES["LANGEVIN_LAMBDA"] / 1000 / 27.2114 # meV --> a.u. # TODO Change units in read_input.py
#     masses = np.array([ np.array([m,m,m]) for m in DYN_PROPERTIES["MASSES"] ])
#     dtI     = DYN_PROPERTIES["dtI"]
#     aL = (2*masses-dtI*LANGEVIN_LAMBDA)/(2*masses+dtI*LANGEVIN_LAMBDA) # Langevin constant
#     v += aL * ((1 - 1/aL) * (v - 0.5 * (F1 + F2) * dtN / M) + (1 - 1/aL) * 0.5 * F1 * dtI / M + (bL + bL/aL) * η)


def Nuclear_X_Step(DYN_PROPERTIES):

    masses  = DYN_PROPERTIES["MASSES"]
    dtI     = DYN_PROPERTIES["dtI"]
    NAtoms  = DYN_PROPERTIES["NAtoms"]

    masses = np.array([ np.array([m,m,m]) for m in masses ])

    # Save previous step
    DYN_PROPERTIES["Atom_coords_old"] = DYN_PROPERTIES["Atom_coords_new"] * 1.0

    # Propagate nuclear coordinates
    DYN_PROPERTIES["FORCE_NEW"] = get_Force(DYN_PROPERTIES)
    a = DYN_PROPERTIES["FORCE_NEW"] / masses

    DYN_PROPERTIES["Atom_coords_new"] += DYN_PROPERTIES["Atom_velocs_new"] * dtI + 0.5000000 * a[:,:] * dtI*dtI

    # FUNCTIONALITY FOR LANGEVIN DYNAMICS
    if ( DYN_PROPERTIES["MD_ENSEMBLE"] == "NVT" ):
        if ( DYN_PROPERTIES["NVT_TYPE"] == "LANGEVIN" ):

            LANGEVIN_LAMBDA = DYN_PROPERTIES["LANGEVIN_LAMBDA"] / 1000 / 27.2114 # meV --> a.u. # TODO Change units in read_input.py
            TEMP  = DYN_PROPERTIES["TEMP"] * (0.025 / 300) / 27.2114 # K -> KT (a.u.) # TODO Change units in read_input.py

            DYN_PROPERTIES["G_RAND"] = gauss(0,1) # Gaussian random number
            SIGMA = np.sqrt(LANGEVIN_LAMBDA*dtI*TEMP/(2*masses**2)) # Langevin constant
            cL = 2*masses*dtI/(2*masses+dtI*LANGEVIN_LAMBDA) # Langevin constant
            R_RAND = (1-dtI/cL) * 0.5 * DYN_PROPERTIES["FORCE_NEW"] * dtI / masses + SIGMA * DYN_PROPERTIES["G_RAND"]
            R_DAMP = cL * (1-dtI/cL) * DYN_PROPERTIES["Atom_velocs_new"]

            DYN_PROPERTIES["Atom_coords_new"] += R_DAMP + R_RAND

        #elif ( DYN_PROPERTIES["NVT_TYPE"] == "RESCALE" ):
            # Do nothing for atomic positions

    return DYN_PROPERTIES

def Nuclear_V_Step(DYN_PROPERTIES):

    masses  = DYN_PROPERTIES["MASSES"]
    dtI     = DYN_PROPERTIES["dtI"]
    NAtoms  = DYN_PROPERTIES["NAtoms"]

    masses = np.array([ np.array([m,m,m]) for m in masses ])

    # Save previous step
    DYN_PROPERTIES["Atom_velocs_old"] = DYN_PROPERTIES["Atom_velocs_new"] * 1.0
    DYN_PROPERTIES["FORCE_OLD"] = DYN_PROPERTIES["FORCE_NEW"] * 1.0 # Store old force
    DYN_PROPERTIES["Atom_velocs_old"] += DYN_PROPERTIES["Atom_velocs_new"]

    # Get new force
    DYN_PROPERTIES["FORCE_NEW"] = get_Force(DYN_PROPERTIES)

    # Compute accelerations
    anew = DYN_PROPERTIES["FORCE_NEW"] / masses 
    aold = DYN_PROPERTIES["FORCE_OLD"] / masses

    DYN_PROPERTIES["Atom_velocs_new"] += 0.5000000 * (aold[:,:] + anew[:,:]) * dtI

    # FUNCTIONALITY FOR LANGEVIN DYNAMICS
    if ( DYN_PROPERTIES["MD_ENSEMBLE"] == "NVT" ):
        if ( DYN_PROPERTIES["NVT_TYPE"] == "LANGEVIN" ):

            LANGEVIN_LAMBDA = DYN_PROPERTIES["LANGEVIN_LAMBDA"] / 1000 / 27.2114 # meV --> a.u. # TODO Change units in read_input.py
            TEMP  = DYN_PROPERTIES["TEMP"] * (0.025 / 300) / 27.2114 # K -> KT (a.u.) # TODO Change units in read_input.py

            SIGMA = np.sqrt(LANGEVIN_LAMBDA*dtI*TEMP/(2*masses**2)) # Langevin constant
            F1 = DYN_PROPERTIES["FORCE_NEW"]
            F2 = DYN_PROPERTIES["FORCE_NEW"]
            aL = (2*masses-dtI*LANGEVIN_LAMBDA)/(2*masses+dtI*LANGEVIN_LAMBDA) # Langevin constant

            V_DAMP = aL * (1 - 1/aL) * (DYN_PROPERTIES["Atom_velocs_new"] - 0.5 * (F1 + F2) * dtI / masses)
            V_RAND = (1 - 1/aL) * 0.5 * F1 * dtI / masses + (SIGMA + SIGMA/aL) * DYN_PROPERTIES["G_RAND"]

            DYN_PROPERTIES["Atom_velocs_new"] += V_DAMP + V_RAND

        if ( DYN_PROPERTIES["NVT_TYPE"] == "RESCALE" ):
            if ( DYN_PROPERTIES["MD_STEP"] % DYN_PROPERTIES["RESCALE_FREQ"] == 0 ):
                TEMP_GOAL      = DYN_PROPERTIES["TEMP"] * (0.025 / 300) / 27.2114 # K -> KT (a.u.) # TODO Change units in read_input.py
                TEMP_NOW       = properties.computer_Temperature(DYN_PROPERTIES)* (0.025 / 300) / 27.2114  # K -> KT (a.u.)
                scale_factor   = np.sqrt( TEMP_GOAL / TEMP_NOW )
                DYN_PROPERTIES["Atom_velocs_new"] *= scale_factor

    return DYN_PROPERTIES