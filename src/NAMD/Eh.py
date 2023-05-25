import numpy as np
import random

import properties


def initialize_mapping(DYN_PROPERTIES):

    NStates = DYN_PROPERTIES["NStates"]
    ISTATE  = DYN_PROPERTIES["ISTATE"]
    #DYN_PROPERTIES["ZPE"] = 0.0 # Ehrenfest has no ZPE 

    """
    ### MMST Style Initialization ###
    #z = np.zeros(( NStates ), dtype=complex)
    #z[ISTATE] = 1.0 + 0.0j # Ehrenfest has no electronic sampling
    """

    ### Spin-mapping Style Initialization ###
    Rw = 2*np.sqrt(NStates+1) # Radius of W Sphere

    # Initialize mapping radii
    r = np.zeros(( NStates )) # np.ones(( NStates )) * np.sqrt(gw)
    r[ISTATE] = np.sqrt( 2 )

    # Set mapping variables
    z = np.zeros(( NStates ),dtype=complex)
    for i in range(NStates):
        phi = random.random() * 2 * np.pi # Azimuthal Angle -- Always Random
        z[i] = r[i] * ( np.cos( phi ) + 1j * np.sin( phi ) )
  
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

    Zreal = np.real(z) * 1.0
    Zimag = np.imag(z) * 1.0

    OVERLAP  = (DYN_PROPERTIES["OVERLAP_NEW"])

    Hamt0 = np.zeros(( NStates, NStates )) # t0 basis
    Hamt1 = np.zeros(( NStates, NStates )) # t1 basis

    #### t0 Ham ####
    Ead_old    = DYN_PROPERTIES["DIAG_ENERGIES_OLD"]
    E_GS_t0    = Ead_old[0] * 1.0
    Hamt0[:,:] = np.diag(Ead_old) 
    Hamt0     -= np.identity(NStates) * E_GS_t0

    #### t1 Ham in t0 basis ####
    Ead_new    = DYN_PROPERTIES["DIAG_ENERGIES_NEW"]
    Hamt1[:,:] = np.diag(Ead_new)
        
    Hamt1 = rotate_t1_to_t0( DYN_PROPERTIES["OVERLAP_NEW"] , Hamt1 ) # Rotate to t0 basis

    Hamt1 -= np.identity(NStates) * E_GS_t0 # Subtract t0 reference energy 'after' rotation to t0 basis

    dtE    = DYN_PROPERTIES["dtE"]
    ESTEPS = DYN_PROPERTIES["ESTEPS"]

    if ( DYN_PROPERTIES["EL_PROP"] == "VV" ):
        """
        Propagate with second-order symplectic (Velocity-Verlet-like)
        """
        #print("Propagation Norm (0):")
        POP = np.sum((0.500000 * np.outer( np.conjugate(z), z ))[np.diag_indices(len(z))])
        #print(np.real(np.round(POP,8)))
        for step in range( ESTEPS ):
            """
            ARK: Do we need linear interpolation ? 
            # If we ignore it, we can analytically evolve the MVs in the diagonal basis.
            # BMW, ~ time-saved is probably too small to implement this. 
            #      ~ Although, it would be more accurate overall.
            """

            # Linear interpolation betwen t0 and t1
            if ( DYN_PROPERTIES["EL_INTERPOLATION"] ):
                H = Hamt0 + (step)/(ESTEPS) * ( Hamt1 - Hamt0 )
            else:
                H = Hamt1
            # Propagate Imaginary first by dt/2
            Zimag -= 0.5000000 * H @ Zreal * dtE

            # Propagate Real by full dt
            Zreal += H @ Zimag * dtE
            
            # Propagate Imaginary final by dt/2
            Zimag -= 0.5000000 * H @ Zreal * dtE

    elif ( DYN_PROPERTIES["EL_PROP"] == "RK" ):
        """
        Propagate with explicit 4th-order Runge-Kutta
        """

        def get_H( step, dt ):
            # Linear interpolation betwen t0 and t1
            if ( DYN_PROPERTIES["EL_INTERPOLATION"] ):
                H = Hamt0 + (step)/(ESTEPS) * ( Hamt1 - Hamt0 )
            else:
                H = Hamt1
            return H

        def f( y, H ):
            return -1j * H @ y

        yt = z.copy()

        for step in range( ESTEPS ):

            k1 = f(yt, get_H( step, 0 ))
            k2 = f(yt + k1*dtE/2, get_H( step, k1*dtE/2 ))
            k3 = f(yt + k2*dtE/2, get_H( step, k2*dtE/2 ))
            k4 = f(yt + k3*dtE, get_H( step, k3*dtE ))

            yt += 1/6 * ( k1 + 2*k2 + 2*k3 + k4 ) * dtE

        z = yt

    DYN_PROPERTIES["MAPPING_VARS"] = z

    #print("Propagation Norm (1):")
    POP = np.sum((0.500000 * np.outer( np.conjugate(z), z ))[np.diag_indices(len(z))])
    #print(np.real(np.round(POP,8)))

    return DYN_PROPERTIES

def rotate_Mapping(DYN_PROPERTIES):
    z = DYN_PROPERTIES["MAPPING_VARS"]
    S = DYN_PROPERTIES["OVERLAP_NEW"]

    #print("Rotation Norm (0):")
    POP = np.sum((0.500000 * np.outer( np.conjugate(z), z ))[np.diag_indices(len(z))])
    #print(np.real(np.round(POP,8)))
    z = rotate_t0_to_t1( S, z )
    #print("Rotation Norm (1):")
    POP = np.sum((0.500000 * np.outer( np.conjugate(z), z ))[np.diag_indices(len(z))])
    #print(np.real(np.round(POP,8)))

    DYN_PROPERTIES["MAPPING_VARS"] = z
    
    DYN_PROPERTIES = check_Mapping_Normalization(DYN_PROPERTIES)

    return DYN_PROPERTIES

def get_density_matrix( DYN_PROPERTIES ):
    z = DYN_PROPERTIES["MAPPING_VARS"]
    return 0.500000 * np.outer( np.conjugate(z), z )

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