import numpy as np

def reorder_all_properties( FROM_STATE,TO_STATE,DYN_PROPERTIES ):

    # ENERGY
    E = DYN_PROPERTIES["DIAG_ENERGIES_NEW"]
    print( "ENERGY\n",E )
    E[TO_STATE], E[FROM_STATE] = E[FROM_STATE], E[TO_STATE]
    DYN_PROPERTIES["DIAG_ENERGIES_NEW"] = E
    print( "ENERGY\n",E )

    # OVERLAP
    OVLP = DYN_PROPERTIES["OVERLAP_NEW"]
    print( "OVLP\n",OVLP )
    OVLP[[TO_STATE,FROM_STATE],:] = OVLP[[FROM_STATE,TO_STATE],:]
    OVLP[:,[TO_STATE,FROM_STATE]] = OVLP[:,[FROM_STATE,TO_STATE]]
    #OVLP[TO_STATE,FROM_STATE], OVLP[FROM_STATE,TO_STATE] = OVLP[FROM_STATE,TO_STATE], OVLP[TO_STATE,FROM_STATE] 
    DYN_PROPERTIES["OVERLAP_NEW"] = OVLP
    print( "OVLP\n",OVLP )

    # NACT
    NACT = DYN_PROPERTIES["NACT_NEW"]
    print( "NACT\n",NACT )
    NACT[[TO_STATE,FROM_STATE],:] = NACT[[FROM_STATE,TO_STATE],:]
    NACT[:,[TO_STATE,FROM_STATE]] = NACT[:,[FROM_STATE,TO_STATE]]
    #NACT[TO_STATE,FROM_STATE], NACT[FROM_STATE,TO_STATE] = NACT[FROM_STATE,TO_STATE], NACT[TO_STATE,FROM_STATE] 
    DYN_PROPERTIES["NACT_NEW"] = NACT
    print( "NACT\n",NACT )

    # NACR
    NACR = DYN_PROPERTIES["NACR_APPROX_NEW"]
    #print( "NACR\n",NACR )
    NACR[[TO_STATE,FROM_STATE],:] = NACR[[FROM_STATE,TO_STATE],:]
    NACR[:,[TO_STATE,FROM_STATE]] = NACR[:,[FROM_STATE,TO_STATE]]
    #NACR[TO_STATE,FROM_STATE], NACR[FROM_STATE,TO_STATE] = NACR[FROM_STATE,TO_STATE], NACR[TO_STATE,FROM_STATE] 
    DYN_PROPERTIES["NACR_APPROX_NEW"] = NACR
    #print( "NACR\n",NACR )

    # Diagonal Gradients
    GRAD = DYN_PROPERTIES["DIAG_GRADIENTS"]
    #print( "GRAD\n",GRAD )
    GRAD[[TO_STATE,FROM_STATE],:,:] = GRAD[[FROM_STATE,TO_STATE],:,:]
    DYN_PROPERTIES["DIAG_GRADIENTS"] = GRAD
    #print( "GRAD\n",GRAD)

    # Mapping Variables
    MAP = DYN_PROPERTIES["MAPPING_VARS"]
    print( "MAP\n",MAP )
    MAP[TO_STATE], MAP[FROM_STATE] = MAP[FROM_STATE], MAP[TO_STATE]
    DYN_PROPERTIES["MAPPING_VARS"] = MAP
    print( "MAP\n",MAP )

    return DYN_PROPERTIES

def check_state_labels(DYN_PROPERTIES):
    NStates   = DYN_PROPERTIES['NStates']
    ENERGY    = DYN_PROPERTIES["DIAG_ENERGIES_NEW"] * 1
    ORDER_OLD = [ j for j in range(NStates) ]
    ORDER_NEW = [ j for j in range(NStates) ]
    for j in range( NStates ):
        for k in range( j+1, NStates ):
            if ( ENERGY[j] > ENERGY[k] ): # THIS WOULD BE WRONG, IF TRUE, SINCE j < k
                print("\n\nCORRECTED STATE ORDERING BY ENERGY:")
                ORDER_NEW[j] = k
                ORDER_NEW[k] = j
                TMP          = ENERGY[j]
                ENERGY[j]    = ENERGY[k]
                ENERGY[k]    = TMP
                print("\tOLD:", ORDER_OLD)
                print("\tNEW:", ORDER_NEW)
                print("\tReordering all properties\n\n")
                DYN_PROPERTIES = reorder_all_properties(j,k,DYN_PROPERTIES)


    return DYN_PROPERTIES






"""
def reorder_all_properties(FROM_STATE,TO_STATE,DYN_PROPERTIES):

    # ENERGY
    E = DYN_PROPERTIES["DIAG_ENERGIES_NEW"]
    print( "ENERGY\n",E )
    E[TO_STATE], E[FROM_STATE] = E[FROM_STATE], E[TO_STATE]
    DYN_PROPERTIES["DIAG_ENERGIES_NEW"] = E
    print( "ENERGY\n",E )


    # OVERLAP
    OVLP = DYN_PROPERTIES["OVERLAP_NEW"]
    print( "OVLP\n",OVLP )
    OVLP[TO_STATE,:], OVLP[FROM_STATE,:] = OVLP[FROM_STATE,:], OVLP[TO_STATE,:]
    OVLP[:,TO_STATE], OVLP[:,FROM_STATE] = OVLP[:,FROM_STATE], OVLP[:,TO_STATE]
    DYN_PROPERTIES["OVERLAP_NEW"] = OVLP
    print( "OVLP\n",OVLP )

    # NACT
    NACT = DYN_PROPERTIES["NACT_NEW"]
    print( "NACT\n",NACT )
    NACT[TO_STATE,:], NACT[FROM_STATE,:] = NACT[FROM_STATE,:], NACT[TO_STATE,:]
    NACT[:,TO_STATE], NACT[:,FROM_STATE] = NACT[:,FROM_STATE], NACT[:,TO_STATE]
    DYN_PROPERTIES["NACT_NEW"] = NACT
    print( "NACT\n",NACT )

    # NACR
    NACR = DYN_PROPERTIES["NACR_APPROX_NEW"]
    print( "NACR\n",NACR )
    NACR[TO_STATE,:,:,:], NACR[FROM_STATE,:,:,:] = NACR[FROM_STATE,:,:,:], NACR[TO_STATE,:,:,:]
    NACR[:,TO_STATE,:,:], NACR[:,FROM_STATE,:,:] = NACR[:,FROM_STATE,:,:], NACR[:,TO_STATE,:,:]
    DYN_PROPERTIES["NACR_APPROX_NEW"] = NACR
    print( "NACR\n",NACR )

    # Diagonal Gradients
    GRAD = DYN_PROPERTIES["DIAG_GRADIENTS"]
    print( "GRAD\n",GRAD )
    GRAD[TO_STATE,:,:], GRAD[FROM_STATE,:,:] = GRAD[FROM_STATE,:,:], GRAD[TO_STATE,:]
    DYN_PROPERTIES["DIAG_GRADIENTS"] = GRAD
    print( "GRAD\n",GRAD)

    # Mapping Variables
    MAP = DYN_PROPERTIES["MAPPING_VARS"]
    print( "MAP\n",MAP )
    MAP[TO_STATE], MAP[FROM_STATE] = MAP[FROM_STATE], MAP[TO_STATE]
    DYN_PROPERTIES["MAPPING_VARS"] = MAP
    print( "MAP\n",MAP )

    return DYN_PROPERTIES
"""


# This is an ad-hoc solution. Do not use this one.
"""
def check_for_trivial_crossing(OVERLAP_CORR,DYN_PROPERTIES):

    def flip_row_column( j,k,NStates,MAT ):
        if ( len(MAT.shape) >= 2 and len(MAT) == NStates and len(MAT[0]) == NStates ):
            tmp1 = MAT[:,j] * 1.0
            tmp2 = MAT[j,:] * 1.0
            tmp3 = MAT[:,k] * 1.0
            tmp4 = MAT[k,:] * 1.0
            MAT[:,j] = tmp3 * 1.0
            MAT[:,k] = tmp1 * 1.0
            MAT[j,:] = tmp4 * 1.0
            MAT[k,:] = tmp2 * 1.0
            return MAT
        elif( len(MAT.shape) == 1 ): # Mapping variables
            tmp    = MAT[j] * 1.0
            MAT[j] = MAT[k] * 1.0
            MAT[k] = tmp * 1.0
            return MAT
        
        elif ( len(MAT.shape) >= 2 and len(MAT) == NStates and len(MAT[0]) != NStates ): # Gets energies and gradients
            tmp    = MAT[j] * 1.0
            MAT[j] = MAT[k] * 1.0
            MAT[k] = tmp * 1.0
            return MAT

    NStates = DYN_PROPERTIES["NStates"]
    OVERLAP_TMP = OVERLAP_CORR * 1.0
    for j in range( NStates ):
        for k in range( j+1, NStates ):
            if ( abs(OVERLAP_CORR[j,k]) > 0.9 ): # This will happen across trivial crossings
                SIGN_JK = OVERLAP_CORR[j,k] / abs(OVERLAP_CORR[j,k])
                OVERLAP_CORR = flip_row_column( j,k,NStates,OVERLAP_CORR ) # Switch the identity of the two states
                SIGN_DIAG   = OVERLAP_CORR[j,j] / abs(OVERLAP_CORR[j,j])
                SIGN_JK_NEW = OVERLAP_CORR[j,k] / abs(OVERLAP_CORR[j,k])
                # Make diagonals positive.
                if ( SIGN_DIAG < 0 ):
                    OVERLAP_CORR[j,j] *= -1
                elif ( SIGN_DIAG > 0 ):
                    OVERLAP_CORR[k,k] *= -1
                # Keep original sign of the off-diagonal since we already phase-corrected.
                if ( SIGN_JK_NEW != SIGN_JK ):
                    OVERLAP_CORR[j,k] *= -1

                for PROPERTY in ["MAPPING_VARS", "DIAG_ENERGIES_NEW", "DIAG_GRADIENTS"]:
                    DYN_PROPERTIES[PROPERTY] = flip_row_column( j,k,NStates,DYN_PROPERTIES[PROPERTY] )
                for key, value in DYN_PROPERTIES.items():
                    print (key)

    return OVERLAP_CORR, DYN_PROPERTIES
"""