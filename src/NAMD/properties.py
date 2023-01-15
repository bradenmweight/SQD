import numpy as np

### NAMD METHODS ###
import Eh
import spinLSC

def get_density_matrix( DYN_PROPERTIES ):
    if ( DYN_PROPERTIES["NAMD_METHOD"] == "EH" ):
        return Eh.get_density_matrix(DYN_PROPERTIES)
    elif ( DYN_PROPERTIES["NAMD_METHOD"] == "SPINLSC" ):
        return spinLSC.get_density_matrix(DYN_PROPERTIES)
    else:
        print("NAMD_METHOD not recognized. Quitting.")
        exit()


def compute_KE(DYN_PROPERTIES):
    KE = 0.0
    for at in range( DYN_PROPERTIES["NAtoms"] ):
        KE += 0.50000000 * DYN_PROPERTIES["MASSES"][at] * np.linalg.norm(DYN_PROPERTIES["Atom_velocs_new"][at,:])**2
    DYN_PROPERTIES["KE"] = KE
    return DYN_PROPERTIES

def compute_PE(DYN_PROPERTIES):
    PE = 0.0
    RHO = get_density_matrix(DYN_PROPERTIES)
    if ( DYN_PROPERTIES["NStates"] >= 2 ):
        for state in range( DYN_PROPERTIES["NStates"] ):
            PE += RHO[state,state].real * DYN_PROPERTIES["DIAG_ENERGIES_NEW"][state]
        DYN_PROPERTIES["PE"] = PE
    else:
        DYN_PROPERTIES["PE"] = DYN_PROPERTIES["DIAG_ENERGIES_NEW"]

    return DYN_PROPERTIES

def computer_Temperature(DYN_PROPERTIES):
    DYN_PROPERTIES = compute_KE(DYN_PROPERTIES)
    NAtoms = DYN_PROPERTIES['NAtoms']

    KE = DYN_PROPERTIES["KE"] * 27.2114 # a.u. --> eV

    T = (2/3) * KE / NAtoms * (300 / 0.025) # eV --> K

    return T