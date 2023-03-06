import numpy as np
import os




def read_Gradients(DIAG_GRADIENTS,DYN_PROPERTIES):

    NStates = DYN_PROPERTIES["NStates"]
    BOMD    = DYN_PROPERTIES["BOMD"]
    ISTATE  = DYN_PROPERTIES["ISTATE"]

    def read_FCHK_GRAD(NAtoms):
        lines = open("geometry.fchk","r").readlines()
        grads = []
        for count, line in enumerate(lines):
            t = line.split()
            if ( len(t) == 5 and t[:2] == "Cartesian Gradient".split() ):
                counter = 1
                while (True):
                    s = lines[count+counter].split()
                    columns = len(s)
                    for col in range(columns):
                        grads.append( float(s[col]) )
                    counter += 1
                    if ( len(grads) == 3*NAtoms ):
                        break
        return np.array(grads,dtype=float).reshape((NAtoms,3))

    if ( DYN_PROPERTIES["CPA"] == False ):

        if ( BOMD == True ):
            if ( ISTATE == 0 ):
                os.chdir("GS_NEW/")
                DIAG_GRADIENTS[0,:,:] = read_FCHK_GRAD( DYN_PROPERTIES["NAtoms"] )
                os.chdir("../")
            else:
                os.chdir(f"TD_NEW_S{ISTATE}/")
                DIAG_GRADIENTS[ISTATE,:,:] = read_FCHK_GRAD( DYN_PROPERTIES["NAtoms"] )
                os.chdir("../")
        else:

            os.chdir("GS_NEW/")
            DIAG_GRADIENTS[0,:,:] = read_FCHK_GRAD( DYN_PROPERTIES["NAtoms"] )
            os.chdir("../")

            if ( NStates >= 2 ):
                for state in range( 1, NStates ):
                    os.chdir(f"TD_NEW_S{state}/")
                    DIAG_GRADIENTS[state,:,:] = read_FCHK_GRAD( DYN_PROPERTIES["NAtoms"] )
                    os.chdir("../")
    
    return DIAG_GRADIENTS

def main(DYN_PROPERTIES):

    NStates = DYN_PROPERTIES["NStates"]
    NAtoms  = DYN_PROPERTIES["NAtoms"]


    DIAG_GRADIENTS = np.zeros(( NStates, NAtoms, 3)) # Diagonal gradients
    read_Gradients(DIAG_GRADIENTS,DYN_PROPERTIES)

    #for state in range( NStates ):
    #    np.savetxt(f"DIAG_GRADIENT_S{state}.dat", DIAG_GRADIENTS[state,:,:], header=f"Diagonal Gradients (S{state}) (NAtoms x 3)", fmt="%2.8f" )

    DYN_PROPERTIES["DIAG_GRADIENTS"] = DIAG_GRADIENTS

    return DYN_PROPERTIES


def read_XYZ():
    XYZ_File = open("geometry_new.xyz","r").readlines()
    NAtoms = int(XYZ_File[0])
    Atom_labels = []
    Atom_coords_new = np.zeros(( NAtoms, 3 ))
    for count, line in enumerate(XYZ_File[2:]):
        t = line.split()
        Atom_labels.append( t[0] )
        Atom_coords_new[count,:] = np.array([ float(t[1]), float(t[2]), float(t[3]) ])

    return Atom_labels, Atom_coords_new

if ( __name__ == "__main__" ):
    
    Atom_labels, Atom_coords_new = read_XYZ()

    DYN_PROPERTIES = {"Atom_labels":Atom_labels, "Atom_coords_new":Atom_coords_new }
    DYN_PROPERTIES["Atom_coords_old"] = DYN_PROPERTIES["Atom_coords_new"] + 0.1
    DYN_PROPERTIES["NStates"] = 4
    DYN_PROPERTIES["NAtoms"] = len(DYN_PROPERTIES["Atom_labels"])
    DYN_PROPERTIES["FUNCTIONAL"] = "BLYP"
    DYN_PROPERTIES["CHARGE"] = 0
    DYN_PROPERTIES["MULTIPLICITY"] = 1
    DYN_PROPERTIES["BASIS_SET"] = "sto-3g"
    DYN_PROPERTIES["MEMORY"] = 5
    DYN_PROPERTIES["NCPUS"] = 1
    DYN_PROPERTIES["MD_STEP"] = 1
    DYN_PROPERTIES["RUN_ELEC_STRUC"] = "USE_CURRENT_NODE" # "SUBMIT_SBATCH", "USE_CURRENT_NODE", "TEST"
    DYN_PROPERTIES["SBATCH_G16"] = "~/submit_scripts/submit.gaussian" # For "SUBMIT_SBATCH" in previous only
    DYN_PROPERTIES = main(DYN_PROPERTIES)