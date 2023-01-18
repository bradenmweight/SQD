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

def do_Langevin_XStep(DYN_PROPERTIES):
    """
    Mark Tuckerman -- Stat. Mech.: Theory and Mol. Simulation
    Chapter 15.5 Page 594 Eq. 15.5.18
    """
    masses  = np.array([ np.array([m,m,m]) for m in DYN_PROPERTIES["MASSES"] ])
    dtI     = DYN_PROPERTIES["dtI"]
    NAtoms  = DYN_PROPERTIES["NAtoms"]
    LANGEVIN_LAMBDA = DYN_PROPERTIES["LANGEVIN_LAMBDA"] / 1000 / 27.2114 # meV --> a.u. # TODO Change units in read_input.py
    TEMP  = DYN_PROPERTIES["TEMP"] * (0.025 / 300) / 27.2114 # K -> KT (a.u.) # TODO Change units in read_input.py

    DYN_PROPERTIES["G_RAND_EPSILON"] = np.array([gauss(0,1) for dof in range(3*NAtoms)], dtype=float).reshape((NAtoms,3)) # Gaussian random number
    DYN_PROPERTIES["G_RAND_THETA"]   = np.array([gauss(0,1) for dof in range(3*NAtoms)], dtype=float).reshape((NAtoms,3)) # Gaussian random number

    # Difference in acceleration and damped velocity
    a_ORIG = DYN_PROPERTIES["FORCE_NEW"]/masses  # Original acceleration
    a_DAMP = LANGEVIN_LAMBDA * DYN_PROPERTIES["Atom_velocs_new"] # Acceleration due to damping
    SIGMA = np.sqrt(2 * TEMP * LANGEVIN_LAMBDA / masses) # Gaussian Width
    RANDOM_FAC = 0.5 * DYN_PROPERTIES["G_RAND_EPSILON"] + 1/(2*np.sqrt(3)) * DYN_PROPERTIES["G_RAND_THETA"]

    # A(t) has units of position
    # Store this for VStep without changing.
    DYN_PROPERTIES["LANGEVIN_A"] = 0.5 * dtI**2 * (a_ORIG - a_DAMP) + SIGMA * dtI**(3/2) * RANDOM_FAC
    
    DYN_PROPERTIES["Atom_coords_new"] += dtI * DYN_PROPERTIES["Atom_velocs_new"] + DYN_PROPERTIES["LANGEVIN_A"]

    return DYN_PROPERTIES



def do_Langevin_VStep(DYN_PROPERTIES):
    """
    Mark Tuckerman -- Stat. Mech.: Theory and Mol. Simulation
    Chapter 15.5 Page 594 Eq. 15.5.18
    """
    masses  = np.array([ np.array([m,m,m]) for m in DYN_PROPERTIES["MASSES"] ])
    dtI     = DYN_PROPERTIES["dtI"]
    NAtoms  = DYN_PROPERTIES["NAtoms"]
    LANGEVIN_LAMBDA = DYN_PROPERTIES["LANGEVIN_LAMBDA"] / 1000 / 27.2114 # meV --> a.u. # TODO Change units in read_input.py
    TEMP  = DYN_PROPERTIES["TEMP"] * (0.025 / 300) / 27.2114 # K -> KT (a.u.) # TODO Change units in read_input.py

    SIGMA = np.sqrt(2 * TEMP * LANGEVIN_LAMBDA / masses) # Gaussian Width

    F0 = DYN_PROPERTIES["FORCE_OLD"]
    F1 = DYN_PROPERTIES["FORCE_NEW"]
    LANGEVIN_A = DYN_PROPERTIES["LANGEVIN_A"]

    DYN_PROPERTIES["Atom_velocs_new"] += 0.5000000 * dtI * ( F0 + F1 ) / masses - \
                                         dtI * LANGEVIN_LAMBDA * DYN_PROPERTIES["Atom_velocs_new"] + \
                                         SIGMA * np.sqrt(dtI) * DYN_PROPERTIES["G_RAND_EPSILON"] - \
                                         LANGEVIN_LAMBDA * LANGEVIN_A

    return DYN_PROPERTIES


def do_VELOC_RESCALE(DYN_PROPERTIES):
    TEMP_GOAL      = DYN_PROPERTIES["TEMP"] * (0.025 / 300) / 27.2114 # K -> KT (a.u.) # TODO Change units in read_input.py
    TEMP_NOW       = properties.compute_Temperature(DYN_PROPERTIES)* (0.025 / 300) / 27.2114  # K -> KT (a.u.)
    scale_factor   = np.sqrt( TEMP_GOAL / TEMP_NOW )
    DYN_PROPERTIES["Atom_velocs_new"] *= scale_factor
    return DYN_PROPERTIES


def Nuclear_X_Step(DYN_PROPERTIES):

    masses  = DYN_PROPERTIES["MASSES"]
    dtI     = DYN_PROPERTIES["dtI"]
    NAtoms  = DYN_PROPERTIES["NAtoms"]

    masses = np.array([ np.array([m,m,m]) for m in masses ])

    # Save previous step
    DYN_PROPERTIES["Atom_coords_old"] = DYN_PROPERTIES["Atom_coords_new"] * 1.0

    # Get the new force
    DYN_PROPERTIES["FORCE_NEW"] = get_Force(DYN_PROPERTIES)
    a = DYN_PROPERTIES["FORCE_NEW"] / masses

    # FUNCTIONALITY FOR CANONICAL ENSEMBLE
    if ( DYN_PROPERTIES["MD_ENSEMBLE"] == "NVT" ):
        if ( DYN_PROPERTIES["NVT_TYPE"] == "LANGEVIN" ):
            DYN_PROPERTIES = do_Langevin_XStep(DYN_PROPERTIES)
        
        # IF DYN_PROPERTIES["NVT_TYPE"] == "RESCALE", POSITION IS UNAFFECTED
    
    elif ( DYN_PROPERTIES["MD_ENSEMBLE"] == "NVE" ):
        DYN_PROPERTIES["Atom_coords_new"] += DYN_PROPERTIES["Atom_velocs_new"] * dtI + 0.5000000 * a[:,:] * dtI*dtI

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

    # FUNCTIONALITY FOR CANONICAL ENSEMBLE
    if ( DYN_PROPERTIES["MD_ENSEMBLE"] == "NVT" ):
        if ( DYN_PROPERTIES["NVT_TYPE"] == "LANGEVIN" ):
            DYN_PROPERTIES = do_Langevin_VStep(DYN_PROPERTIES)

        elif ( DYN_PROPERTIES["NVT_TYPE"] == "RESCALE" ):
            if ( DYN_PROPERTIES["MD_STEP"] % DYN_PROPERTIES["RESCALE_FREQ"] == 0 ):
                DYN_PROPERTIES = do_VELOC_RESCALE(DYN_PROPERTIES)

    elif ( DYN_PROPERTIES["MD_ENSEMBLE"] == "NVE" ):
        DYN_PROPERTIES["Atom_velocs_new"] += 0.5000000 * (aold[:,:] + anew[:,:]) * dtI

    return DYN_PROPERTIES