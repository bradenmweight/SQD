import numpy as np
import random

import properties


def initialize_mapping(DYN_PROPERTIES):

    NStates = DYN_PROPERTIES["NStates"]
    ISTATE  = DYN_PROPERTIES["ISTATE"]

    ### Spin-mapping Style Initialization ###
    Rw = 2*np.sqrt(NStates+1) # Radius of W Sphere
    DYN_PROPERTIES["ZPE"] = (2/NStates) * (np.sqrt(NStates + 1) - 1) # We choose the W-sphere
    gw = DYN_PROPERTIES["ZPE"]

    # Initialize mapping radii
    r = np.ones(( NStates )) * np.sqrt(gw)
    r[ISTATE] = np.sqrt( 2 + gw )

    # Set real mapping variables
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
    Ead     = DYN_PROPERTIES["DIAG_ENERGIES_NEW"]

    if ( DYN_PROPERTIES["MD_STEP"] >= 1 ):
        NACR    = DYN_PROPERTIES["NACR_APPROX_NEW"] # NStates x NStates x NAtoms x 3 (a.u.)
    else:
        NACR    = np.zeros(( NStates, NStates, NAtoms, 3 )) # NStates x NStates x NAtoms x 3 (a.u.)

    z = DYN_PROPERTIES["MAPPING_VARS"]

    # Do we need to find the state-independent force ?
    # It will already be included in the dEad terms at this point...
    #F0 = np.einsum("jad->ad", dEad[:,:,:]) / NStates

    F = np.zeros(( NAtoms, 3 ))
    rho = np.real( properties.get_density_matrix(DYN_PROPERTIES) )
    for j in range( NStates ):
        for k in range( j, NStates ):
            if ( j == k ):
                F[:,:] -= dEad[j,:,:] * rho[j,j]
            else:
                Ejk = Ead[j] - Ead[k]
                F[:,:] -= 2 * rho[j,k] * NACR[j,k,:,:] * Ejk  # Double count upper triangle


    return F

def rotate_t0_to_t1(S, A): # Recall, we perform TD-DFT with one additional state. Already removed from overlap.
    if ( len(A.shape) == 1 ):
        return S.T @ A
        #return S @ A
    elif( len(A.shape) == 2 ):
        return S.T @ A @ S
        #return S @ A @ S.T
    else:
        print("Shape of rotating object not correct." )
        print(f"Needs to be either 1D or 2D numpy array. Received {len(A.shape)}D array.")

def rotate_t1_to_t0(S, A): # Recall, we perform TD-DFT with one additional state. Already removed from overlap at this point.
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

    OVERLAP  = (DYN_PROPERTIES["OVERLAP_NEW"]) # Recall, we perform TD-DFT with one additional state. Already removed.

    Hamt0 = np.zeros(( NStates, NStates )) # t0 basis
    Hamt1 = np.zeros(( NStates, NStates )) # t1 basis

    #### t0 Ham ####
    # ADD NACT = dR/dt.NACR TO OFF-DIAGONALS
    if ( DYN_PROPERTIES["MD_STEP"] >= 2 ): 
        NACR_OLD  = DYN_PROPERTIES["NACR_APPROX_OLD"]
        VELOC_OLD = DYN_PROPERTIES["Atom_velocs_old"]
    Ead_old    = DYN_PROPERTIES["DIAG_ENERGIES_OLD"]
    E_GS_t0    = Ead_old[0] * 1.0
    Hamt0      = np.diag(Ead_old) 
    Hamt0     -= np.identity(NStates) * E_GS_t0

    #### t1 Ham in t0 basis ####
    NACR_NEW   = DYN_PROPERTIES["NACR_APPROX_NEW"]
    VELOC_NEW  = DYN_PROPERTIES["Atom_velocs_new"]
    Ead_new    = DYN_PROPERTIES["DIAG_ENERGIES_NEW"]
    Hamt1[:,:] = np.diag(Ead_new) #- np.identity(NStates) * E_GS_t0 # Is it okay to subtract this here ? Identity will also rotate, no ? Maybe bad.
        
    # TODO Check the direction of this rotation
    Hamt1 = rotate_t1_to_t0( DYN_PROPERTIES["OVERLAP_NEW"] , Hamt1 ) # Rotate to t0 basis

    # Shift by reference energy for mapping evolution
    Hamt1 -= np.identity(NStates) * E_GS_t0 # Subtract t0 reference 'after' rotation to t0 basis

    dtE    = DYN_PROPERTIES["dtE"]
    ESTEPS = DYN_PROPERTIES["ESTEPS"]

    if ( DYN_PROPERTIES["EL_PROP"] == "VV" ):
        """
        Propagate with second-order symplectic (Velocity-Verlet-like)
        """

        for step in range( ESTEPS ):
            """
            ARK: Do we need linear interpolation ? 
            # If we ignore it, we can analytically evolve the MVs in the diagonal basis.
            # BMW, ~ time-saved is probably too small to implement this. 
            #      ~ Although, it would be more accurate overeall.
            """

            # Linear interpolation betwen t0 and t1
            H = Hamt0 + (step)/(ESTEPS) * ( Hamt1 - Hamt0 )

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
            H = Hamt0 + (step + dt/dtE)/(ESTEPS) * ( Hamt1 - Hamt0 )
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

    return DYN_PROPERTIES

def rotate_Mapping(DYN_PROPERTIES):
    z = DYN_PROPERTIES["MAPPING_VARS"]
    S = DYN_PROPERTIES["OVERLAP_NEW"]

    z = rotate_t0_to_t1( S, z )

    DYN_PROPERTIES["MAPPING_VARS"] = z
    
    DYN_PROPERTIES = check_Mapping_Normalization(DYN_PROPERTIES)

    return DYN_PROPERTIES

def get_density_matrix( DYN_PROPERTIES ):
    z = DYN_PROPERTIES["MAPPING_VARS"]
    gw = DYN_PROPERTIES["ZPE"]
    NStates = DYN_PROPERTIES["NStates"]
    return 0.500000 * (np.outer( np.conjugate(z), z ) - gw * np.identity(NStates) )

def check_Mapping_Normalization(DYN_PROPERTIES):
    z = DYN_PROPERTIES["MAPPING_VARS"]
    POP = np.real(properties.get_density_matrix( DYN_PROPERTIES )[np.diag_indices(DYN_PROPERTIES["NStates"])])
    norm = np.sum( POP )
    if ( abs(1.0 - norm) > 1e-5 and abs(1.0 - norm) < 1e-2 ):
        print(f"Mapping norm is wrong: {np.round(norm,4)} != 1.00000")
    if( abs(1.0 - norm) > 1e-2 ):
        print(f"Electronic Norm.: {np.round(norm,4)}")
        print(f"ERROR: Mapping norm is VERY wrong: {np.round(norm,4)} != 1.00")
        print("ERROR: Check if we should renormalize. If not, this trajectory may be trash.")
        if ( DYN_PROPERTIES["FORCE_MAP_NORM"] == True ):
            print("For spin-LSC, the populations can be negative. Is this renormalization okay ?")
            DYN_PROPERTIES["MAPPING_VARS"] /= np.sqrt(norm)
        else:
            print("\tUser chose not to renormalize.")
    
    return DYN_PROPERTIES