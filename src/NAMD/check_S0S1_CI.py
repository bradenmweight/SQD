from properties import get_density_matrix
from wrappers import initialize_mapping



def main( DYN_PROPERTIES ):

    if ( DYN_PROPERTIES["S0S1_CI_FLAG"] == True ): # Check if we are near S0/S1 CI

        E       = DYN_PROPERTIES["DIAG_ENERGIES_NEW"] * 27.2114
        NSTATES = DYN_PROPERTIES["NStates"]
        
        for state in range( 1, NSTATES ):
            if ( abs( E[state] - E[0] ) < 0.1 ): # Arbitrary threshold (eV)
                # Change settings to BOMD in S0 (GS)
                DYN_PROPERTIES["BOMD"]            = True
                DYN_PROPERTIES["ISTATE"]          = 0
                DYN_PROPERTIES                    = initialize_mapping(DYN_PROPERTIES)
                ###DYN_PROPERTIES["MAPPING_VARS"]    = np.zeros( (NSTATES), dtype=complex )
                ###DYN_PROPERTIES["MAPPING_VARS"][0] = 1 + 0j
                # Keep the same number of excited states in the calculation
                print("WARNING:\n\tFound S0/S1 conical intersection. Starting BOMD in S0.")
                break
    return DYN_PROPERTIES