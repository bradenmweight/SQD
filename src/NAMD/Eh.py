import numpy as np
import random

import properties


def initialize_mapping(DYN_PROPERTIES):

    NStates = DYN_PROPERTIES["NStates"]
    ISTATE  = DYN_PROPERTIES["ISTATE"]
    #DYN_PROPERTIES["ZPE"] = 0.0 # Ehrenfest has no ZPE 

    ### MMST Style Initialization ###
    z = np.zeros(( NStates ), dtype=complex)
    z[ISTATE] = 1.0 + 0.0j # Ehrenfest has no electronic sampling

    DYN_PROPERTIES["MAPPING_VARS"] = z

    # Check initial density matrix
    #RHO = get_density_matrix(DYN_PROPERTIES)
    #print("Initial Density Matrix:")
    #print( RHO )

    return DYN_PROPERTIES

def get_Force(DYN_PROPERTIES):
    
    dEad    = DYN_PROPERTIES["DIAG_GRADIENTS"] # NStates x NAtoms x 3 (a.u.)
    NAtoms  = DYN_PROPERTIES["NAtoms"]
    NStates = DYN_PROPERTIES["NStates"]
    MD_STEP = DYN_PROPERTIES["MD_STEP"]
    Ead     = DYN_PROPERTIES["DIAG_ENERGIES_NEW"]
    z = DYN_PROPERTIES["MAPPING_VARS"]


    F = np.zeros(( NAtoms, 3 ))
    rho = np.real( properties.get_density_matrix(DYN_PROPERTIES) )
    if ( DYN_PROPERTIES["CPA"] == True ):
        print("Using CPA forces. F = F(G.S.)")
        F[:,:] = -dEad[0,:,:] # G.S. Forces Only -- For Classical Path Approximation
    else:
        for j in range( NStates ):
            for k in range( j, NStates ):
                if ( j == k ):
                    F[:,:] -= dEad[j,:,:] * rho[j,j]
                else:
                    if ( MD_STEP >= 1 ):
                        NACR    = DYN_PROPERTIES["NACR_APPROX_NEW"] # NStates x NStates x NAtoms x 3 (a.u.)
                        Ejk = Ead[j] - Ead[k]
                        F[:,:] -= 2 * rho[j,k] * NACR[j,k,:,:] * Ejk  # Double count upper triangle

    return F

def rotate_t0_to_t1(S, A): # Recall, we perform TD-DFT with one additional state. Already removed from overlap by this point.
    if ( len(A.shape) == 1 ):
        return S.T @ A
    elif( len(A.shape) == 2 ):
        return S.T @ A @ S
    else:
        print("Shape of rotating object not correct." )
        print(f"Needs to be either 1D or 2D numpy array. Received {len(A.shape)}D array.")

def rotate_t1_to_t0(S, A): # Recall, we perform TD-DFT with one additional state. Already removed from overlap by this point.
    if ( len(A.shape) == 1 ):
        return S @ A
    elif( len(A.shape) == 2 ):
        return S @ A @ S.T
    else:
        print("Shape of rotating object not correct." )
        print(f"Needs to be either 1D or 2D numpy array. Received {len(A.shape)}D array.")

def propagage_Mapping(DYN_PROPERTIES):
    NStates = DYN_PROPERTIES["NStates"]
    z       = DYN_PROPERTIES["MAPPING_VARS"]
    dtI     = DYN_PROPERTIES["dtI"]

    Ead_old = DYN_PROPERTIES["DIAG_ENERGIES_OLD"] # diag in t0 basis
    Ead_new = DYN_PROPERTIES["DIAG_ENERGIES_NEW"] # diag in t1 basis
    OVERLAP = DYN_PROPERTIES["OVERLAP_NEW"]       # <t0|t1>

    # print("Propagating in diagonal basis.")
    z = rotate_t0_to_t1( OVERLAP, z ) # Transform to t1 basis
    z = np.exp( -1j * Ead_new * dtI ) * z # Diagonal propagation in t1 basis

    DYN_PROPERTIES["MAPPING_VARS"] = z

    #print("Propagation Norm (1):")
    #POP = np.sum((0.500000 * np.outer( np.conjugate(z), z ))[np.diag_indices(len(z))])
    #print(np.real(np.round(POP,8)))

    return DYN_PROPERTIES

def rotate_Mapping(DYN_PROPERTIES):
    z = DYN_PROPERTIES["MAPPING_VARS"]
    S = DYN_PROPERTIES["OVERLAP_NEW"]

    #print("Rotation Norm (0):")
    #POP = np.sum((0.500000 * np.outer( np.conjugate(z), z ))[np.diag_indices(len(z))])
    #print(np.real(np.round(POP,8)))
    z = rotate_t0_to_t1( S, z )
    #print("Rotation Norm (1):")
    #POP = np.sum((0.500000 * np.outer( np.conjugate(z), z ))[np.diag_indices(len(z))])
    #print(np.real(np.round(POP,8)))

    DYN_PROPERTIES["MAPPING_VARS"] = z
    
    DYN_PROPERTIES = check_Mapping_Normalization(DYN_PROPERTIES)

    return DYN_PROPERTIES

def get_density_matrix( DYN_PROPERTIES ):
    z = DYN_PROPERTIES["MAPPING_VARS"]
    return np.outer( np.conjugate(z), z )

def check_Mapping_Normalization(DYN_PROPERTIES):
    z = DYN_PROPERTIES["MAPPING_VARS"]
    POP = np.real(properties.get_density_matrix( DYN_PROPERTIES )[np.diag_indices(DYN_PROPERTIES["NStates"])])
    norm = np.sum( POP )
    #print(f"Electronic Norm.: {np.round(norm,4)}")
    if ( abs(1.0 - norm) > 1e-5 and abs(1.0 - norm) < 1e-2 ):
        print(f"Mapping norm is wrong: {np.round(norm,4)} != 1.00000")
    if( abs(1.0 - norm) > 1e-2 ):
        print(f"Electronic Norm.: {np.round(norm,4)}")
        print(f"ERROR: Mapping norm is VERY wrong: {np.round(norm,4)} != 1.00")
        print("ERROR: Check if we should renormalize. If not, this trajectory may be trash.")
        if ( DYN_PROPERTIES["FORCE_MAP_NORM"] == True ):
            DYN_PROPERTIES["MAPPING_VARS"] /= np.sqrt(norm)
        else:
            print("\tUser chose not to renormalize.")
    
    return DYN_PROPERTIES