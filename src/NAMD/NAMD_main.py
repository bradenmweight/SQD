import numpy as np
import sys
import subprocess as sp
from time import time

import read_input
import nuclear_propagation
import output
import rotation

### NAMD METHODS ###
import Eh
import spinLSC

# HOW TO HANDLE THESE PATHS BETTER ?
sys.path.append("/scratch/bweight/software/many_molecule_many_mode_NAMD/src/ELECTRONIC_STRUCTURE_CONTROL/")
sys.path.append("/scratch/bweight/software/many_molecule_many_mode_NAMD/src/WFN_OVERLAP/PYTHON/")

import G16_TD


# This code with be the main control code for the NAMD.
# We will make calls to electronic structure and
#   wavefunction overlaps here, which are hanfled elsewhere.
# Additionally, the mixed-quantum classical or semi-classical
#   dynamics will be handled elsewhere.

def initialize_mapping(DYN_PROPERTIES):
    """
    Wrapper for initializing mapping variables
    """
    if ( DYN_PROPERTIES["NAMD_METHOD"] == "EH" ):
        return Eh.initialize_mapping(DYN_PROPERTIES)
    elif ( DYN_PROPERTIES["NAMD_METHOD"] == "SPINLSC" ):
        return spinLSC.initialize_mapping(DYN_PROPERTIES)

def propagage_Mapping(DYN_PROPERTIES):
    """
    Wrapper for electronic propagation
    """
    if ( DYN_PROPERTIES["NAMD_METHOD"] == "EH" ):
        return Eh.propagage_Mapping(DYN_PROPERTIES)
    elif ( DYN_PROPERTIES["NAMD_METHOD"] == "SPINLSC" ):
        return spinLSC.propagage_Mapping(DYN_PROPERTIES)

def rotate_Mapping(DYN_PROPERTIES):
    """
    Wrapper for the transformation of mapping variables
    """
    if ( DYN_PROPERTIES["NAMD_METHOD"] == "EH" ):
        return Eh.rotate_Mapping(DYN_PROPERTIES)
    if ( DYN_PROPERTIES["NAMD_METHOD"] == "SPINLSC" ):
        return spinLSC.rotate_Mapping(DYN_PROPERTIES)




def main( ):
    DYN_PROPERTIES = read_input.read()
    DYN_PROPERTIES = read_input.initialize_MD_variables(DYN_PROPERTIES)

    # Remove COM motion and angular velocity
    # Do we need to do this at every step. Probably should at least remove COM.
    # With Wigner-sampled geometries and velocities, this is already done.
    if ( DYN_PROPERTIES["REMOVE_COM_MOTION"] == True ):
        DYN_PROPERTIES = rotation.shift_COM(DYN_PROPERTIES)
    if ( DYN_PROPERTIES["REMOVE_ANGULAR_VELOCITY"] == True ):
        DYN_PROPERTIES = rotation.remove_rotations(DYN_PROPERTIES)


    # Initialize electronic DOFs
    DYN_PROPERTIES = initialize_mapping(DYN_PROPERTIES)

    # Perform first electronic structure calculation
        # Get diagonal energies and gradients
    DYN_PROPERTIES = G16_TD.main(DYN_PROPERTIES)
    DYN_PROPERTIES["FORCE_NEW"] = nuclear_propagation.get_Force(DYN_PROPERTIES)
    output.save_data(DYN_PROPERTIES)

    # Start main MD loop
    for step in range( DYN_PROPERTIES["NSteps"] ):
        T_STEP_START = time()
        #print(f"Working on step {step} of { DYN_PROPERTIES['NSteps'] }")

        # Propagate nuclear coordinates
        DYN_PROPERTIES = nuclear_propagation.Nuclear_X_Step(DYN_PROPERTIES)

        # Perform jth electronic structure calculation
            # Get diagonal energies and grad
            # Get overlap and NACT
        DYN_PROPERTIES["MD_STEP"] += 1 # This needs to be exactly here for technical reasons.
        T0 = time()
        DYN_PROPERTIES = G16_TD.main(DYN_PROPERTIES)
        print( "Total QM took %2.2f s." % (time() - T0) )

        if ( DYN_PROPERTIES["NStates"] >= 2 ):
            # Propagate electronic DOFs (Interpolated Hamiltonian)
            T0 = time()
            DYN_PROPERTIES = propagage_Mapping(DYN_PROPERTIES)
            print( "Electronic propagation took %2.2f s." % (time() - T0) )

            # Rotate electronic DOFs from t0 basis to t1 basis
            DYN_PROPERTIES = rotate_Mapping(DYN_PROPERTIES)

        # Propagate nuclear momenta
        DYN_PROPERTIES = nuclear_propagation.Nuclear_V_Step(DYN_PROPERTIES)

        # Remove COM motion and angular velocity
        # Do we need to do this at every step. Probably should at least remove COM.
        # With Wigner-sampled geometries and velocities, this is already done.
        if ( DYN_PROPERTIES["REMOVE_COM_MOTION"] == True ):
            DYN_PROPERTIES = rotation.shift_COM(DYN_PROPERTIES)
        if ( DYN_PROPERTIES["REMOVE_ANGULAR_VELOCITY"] == True ):
            DYN_PROPERTIES = rotation.remove_rotations(DYN_PROPERTIES)

        output.save_data(DYN_PROPERTIES)

        print( "Total MD Step took %2.2f s." % (time() - T_STEP_START) )

if ( __name__ == "__main__" ):
    main()