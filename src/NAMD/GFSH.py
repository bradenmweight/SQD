import numpy as np
import random

import properties


def initialize_mapping(DYN_PROPERTIES):

    NStates = DYN_PROPERTIES["NStates"]
    ISTATE  = DYN_PROPERTIES["ISTATE"]

    # Set mapping variables
    z = np.zeros(( NStates ),dtype=complex)
    z[ISTATE] = 1.0 + 0.0j

    DYN_PROPERTIES["MAPPING_VARS"] = z
    DYN_PROPERTIES["ACTIVE_STATE"] = ISTATE

    return DYN_PROPERTIES


def get_Force(DYN_PROPERTIES):
    
    dEad    = DYN_PROPERTIES["DIAG_GRADIENTS"] # NStates x NAtoms x 3 (a.u.)
    NAtoms  = DYN_PROPERTIES["NAtoms"]
    NStates = DYN_PROPERTIES["NStates"]
    MD_STEP = DYN_PROPERTIES["MD_STEP"]
    Ead     = DYN_PROPERTIES["DIAG_ENERGIES_NEW"]
    z = DYN_PROPERTIES["MAPPING_VARS"]
    AS = DYN_PROPERTIES["ACTIVE_STATE"]

    F = np.zeros(( NAtoms, 3 ))
    if ( DYN_PROPERTIES["CPA"] == True ):
        print("Using CPA forces. F = F(G.S.)")
        F[:,:] = -dEad[0,:,:] # G.S. Forces Only -- For Classical Path Approximation
    else:
        # Surface hopping forces only come from currecnt active state (AS)
        # These are only diagonal forces. Never from off-diagonal non-adiabatic NACR.
        print(f"Using active state forces. F = F(S{AS})")
        F[:,:] = -dEad[AS,:,:]

    return F


def rotate_t0_to_t1(S, A): # Recall, we perform TD-DFT with one additional state. Already removed from overlap by this point.
    if ( len(A.shape) == 1 ):
        return S.T @ A
    elif( len(A.shape) == 2 ):
        return S.T @ A @ S
    else:
        print("Shape of rotating object not correct." )
        print(f"Needs to be either 1D or 2D numpy array. Received {len(A.shape)}D array.")

def rotate_t1_to_t0(S, A): # Recall, we perform TD-DFT with one additional state. Already removed from overlap.
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
    MD_STEP = DYN_PROPERTIES["MD_STEP"]
    DYN_PROPERTIES["MAPPING_VARS_OLD"] = z * 1.0

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
        #POP = np.sum( np.real( properties.get_density_matrix(DYN_PROPERTIES) )[np.diag_indices(NStates)])
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

    #print("Propagation Norm (1):")
    #POP = np.sum( np.real( properties.get_density_matrix(DYN_PROPERTIES) )[np.diag_indices(NStates)])
    #print(np.real(np.round(POP,8)))

    return DYN_PROPERTIES



def rotate_Mapping(DYN_PROPERTIES):
    z       = DYN_PROPERTIES["MAPPING_VARS"]
    S       = DYN_PROPERTIES["OVERLAP_NEW"]
    NStates = DYN_PROPERTIES["NStates"]

    #print("Rotation Norm (0):")
    #POP = np.sum( np.real( properties.get_density_matrix(DYN_PROPERTIES) )[np.diag_indices(NStates)])
    #print(np.real(np.round(POP,8)))
    z = rotate_t0_to_t1( S, z )
    #print("Rotation Norm (1):")
    #POP = np.sum( np.real( properties.get_density_matrix(DYN_PROPERTIES) )[np.diag_indices(NStates)])
    #print(np.real(np.round(POP,8)))

    DYN_PROPERTIES["MAPPING_VARS"] = z
    
    DYN_PROPERTIES = check_Mapping_Normalization(DYN_PROPERTIES)

    DYN_PROPERTIES = get_hop(DYN_PROPERTIES)

    return DYN_PROPERTIES

def get_density_matrix( DYN_PROPERTIES ):
    z = DYN_PROPERTIES["MAPPING_VARS"]
    return np.outer( np.conjugate(z), z )


def check_Mapping_Normalization(DYN_PROPERTIES):
    z = DYN_PROPERTIES["MAPPING_VARS"]
    POP = np.real(properties.get_density_matrix( DYN_PROPERTIES )[np.diag_indices(DYN_PROPERTIES["NStates"])])
    norm = np.sum( POP )
    print(f"Electronic Norm.: {np.round(norm,4)}")
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





##### SH-Specific Functions #####

def get_hop_prob(DYN_PROPERTIES):

    AS      = DYN_PROPERTIES["ACTIVE_STATE"]
    Z_OLD   = DYN_PROPERTIES["MAPPING_VARS_OLD"]
    Z_NEW   = DYN_PROPERTIES["MAPPING_VARS"]
    dtI     = DYN_PROPERTIES["dtI"]
    NStates = DYN_PROPERTIES["NStates"]
    
    RHO_OLD = np.real( np.outer( np.conjugate(Z_OLD), Z_OLD ) )
    RHO_NEW = np.real( properties.get_density_matrix(DYN_PROPERTIES) )

    POP_OLD = np.real( RHO_OLD[np.diag_indices(NStates)] )
    POP_NEW = np.real( RHO_NEW[np.diag_indices(NStates)] )

    # If active state increased in population, GFSH says not to hop
    if ( POP_NEW[AS] > POP_OLD[AS] ):
        print(f"Population in ACTIVE_STATE (AS = S{AS}) increased. Skipping hop check.")
        return None



    #### THIS IS CENTRAL QUANTITY ####
    HOP_PROB = np.zeros(( NStates ))
    ##################################


    ### THIS IS GLOBAL FLUX ROUTINE ###
    SUM_POP_DECREASE = 0.0
    for j in range( NStates ):
        if ( POP_OLD[j] > POP_NEW[j] ):
            SUM_POP_DECREASE += POP_OLD[j] - POP_NEW[j] # Positive quantity
    print("SUM_POP_DECREASE =", SUM_POP_DECREASE)
    for j in range( NStates ):
        if ( j != AS ):
            dP_j  =     POP_NEW[j]  - POP_OLD[j]
            dP_AS = -1*(POP_NEW[AS] - POP_OLD[AS]) # This one is backwards
            A = dP_j  / POP_OLD[AS]
            B = dP_AS / SUM_POP_DECREASE
            HOP_PROB[j] = A * B
            print( "j, dP_j, dP_AS, Prob_j:", j, dP_j, dP_AS, HOP_PROB[j] )
            # TODO -- Need to write CPA probability, too
            # P_CPA ~ MIN[ 0, p_j * exp(-dE_{AS-->j}/kT) ]
            if ( DYN_PROPERTIES["CPA"] == True ):
                Ead = DYN_PROPERTIES["DIAG_ENERGIES_NEW"]
                dE  = Ead[j] - Ead[AS] # NEW - OLD
                kT  = DYN_PROPERTIES["TEMP"] * (0.025 / 300) / 27.2114 # TODO CHANGE TEMP IN 'read_input.py'
                print( "dE, kT, exp[-dE/kT]", dE, kT, np.exp(-dE/kT) )
                HOP_PROB[j] = np.max([ HOP_PROB[j] * np.exp(-dE/kT), 0.0 ])
            else:
                HOP_PROB[j] = np.max([ HOP_PROB[j], 0.0 ])

    print( "HOP PROB:", HOP_PROB )
    return HOP_PROB


def get_hop(DYN_PROPERTIES):
    AS      = DYN_PROPERTIES["ACTIVE_STATE"] * 1 # "1", not "1.0"
    Z_OLD   = DYN_PROPERTIES["MAPPING_VARS_OLD"] * 1.0
    Z_NEW   = DYN_PROPERTIES["MAPPING_VARS"] * 1.0
    NStates = DYN_PROPERTIES["NStates"]
    masses  = DYN_PROPERTIES["MASSES"]
    VELOC   = DYN_PROPERTIES["Atom_velocs_new"] * 1.0
    Ead     = DYN_PROPERTIES["DIAG_ENERGIES_NEW"] * 1.0
    
    
    HOP_PROB = get_hop_prob(DYN_PROPERTIES)
    if ( type(HOP_PROB) == type(None) ): # Skip hop
        return DYN_PROPERTIES

    PROB_SUM = np.zeros(( NStates ))
    for j in range( NStates ):
        for k in range( j ):
            PROB_SUM[j] += HOP_PROB[k]

    AS_OLD = AS * 1
    AS_NEW = AS * 1

    RANDOM = random.random()
    if ( 0.0 <= RANDOM and RANDOM < PROB_SUM[0] ):
        AS_NEW = 0
    else:
        for j in range( 1, NStates ):
            print("RANDOM:", PROB_SUM[j-1], RANDOM, PROB_SUM[j] )
            if ( PROB_SUM[j-1] <= RANDOM and RANDOM < PROB_SUM[j] ):
                print(f"SATISFIED PROBABILITY CRITERION {AS} --> {j-1}")
                AS_NEW = j-1
    
    if ( AS_NEW != AS_OLD ):
        print(f"TRY HOPPING FROM {AS} --> {j}")
        if ( DYN_PROPERTIES["CPA"] == True ):
            PDIFF = Ead[AS_NEW] - Ead[AS_OLD]
            KIN = properties.compute_KE(DYN_PROPERTIES)
            if ( AS_NEW > AS_OLD and KIN < np.abs(PDIFF) ):
                print(f"REJECTING HOP ({AS_OLD}-->{AS_NEW}) BASED ON LACK OF K.E.:")
                DYN_PROPERTIES["HOP_FLAG"] = False
                return DYN_PROPERTIES
            else:
                # Calculate the uniform scaling factor based on energy criterion
                print("Rescaling velocities uniformly with CPA routine.")
                scale_factor = np.sqrt( 1 - PDIFF / KIN )
                VELOC[:,:] *= scale_factor
        else:
            NACR    = DYN_PROPERTIES["NACR_APPROX_NEW"]
            a = 0.50000 * np.einsum("ad,a->",NACR[AS_OLD,AS_NEW,:,:]**2, 1/masses[:] )
            b =           np.einsum("ad,ad->",NACR[AS_OLD,AS_NEW,:,:]**2, VELOC[:,:] )
            dE = Ead[AS_NEW] - Ead[AS_OLD]
            discriminant = b**2 - 4 * a * dE
            if ( discriminant < 0.0 ):
                print(f"Rejecting hop due to lack of K.E. (discriminant = {round(discriminant,6)}).")
                DYN_PROPERTIES["HOP_FLAG"] = False
                return DYN_PROPERTIES
            elif ( b < 0.0 ):
                gamma = (b + np.sqrt(discriminant)) / 2 / a
            else:
                gamma = (b - np.sqrt(discriminant)) / 2 / a
            VELOC[:,:] -= gamma * np.einsum("ad,a->", NACR[AS_OLD,AS_NEW,:,:], 1/masses[:] )
    
        print(f" ##### Hopping from S{AS_OLD}-->S{AS_NEW} #####")
    
        DYN_PROPERTIES["HOP_FLAG"] = True
        DYN_PROPERTIES["ACTIVE_STATE"]    = AS_NEW
        DYN_PROPERTIES["Atom_velocs_new"] = VELOC


    # TODO -- Add decoherence corrections here

    return DYN_PROPERTIES
