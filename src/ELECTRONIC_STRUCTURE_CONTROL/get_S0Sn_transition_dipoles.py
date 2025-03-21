import numpy as np
import subprocess as sp
import os

def read_Dipoles(TRANS_DIPOLES,NStates,BOMD,ISTATE):

    if ( BOMD == True ):
        if ( ISTATE == 0 and NStates == 1 ):
            # Read GS
            os.chdir(f"GS_NEW/")
            sp.call( "grep 'Dipole moment (field-independent basis, Debye):' geometry.out -A 1 | tail -n 1 | awk '{print $2, $4, $6}' > DIP_TMP.dat" ,shell=True)
            TRANS_DIPOLES[0,:] = np.loadtxt("DIP_TMP.dat") / 2.5 # Debye to a.u. # (1,3)
            os.chdir("../")

        elif ( ISTATE == 0 and NStates >= 2 ):
            # Read GS
            os.chdir(f"GS_NEW/")
            sp.call( "grep 'Dipole moment (field-independent basis, Debye):' geometry.out -A 1 | tail -n 1 | awk '{print $2, $4, $6}' > DIP_TMP.dat" ,shell=True)
            TRANS_DIPOLES[0,:] = np.loadtxt("DIP_TMP.dat") / 2.5 # Debye to a.u. # (1,3)
            os.chdir("../")

            # Read Excited States
            os.chdir(f"TD_NEW_S1/")
            sp.call( "grep 'Ground to excited state transition electric dipole moments' geometry.out -A %d | tail -n %d | awk '{print $2, $3, $4}' > DIP_TMP.dat" % (NStates,NStates-1), shell=True)
            TRANS_DIPOLES[1:,:] = np.loadtxt("DIP_TMP.dat") / 2.5 # Debye to a.u. # (N-1,3)
            os.chdir("../")
        
        elif ( ISTATE != 0 and NStates >= 2 ):
            # Read Excited States
            os.chdir(f"TD_NEW_S{ISTATE}/")
            sp.call( "grep 'Ground to excited state transition electric dipole moments' geometry.out -A %d | tail -n %d | awk '{print $2, $3, $4}' > DIP_TMP.dat" % (NStates,NStates-1), shell=True)
            TRANS_DIPOLES[1:,:] = np.loadtxt("DIP_TMP.dat") / 2.5 # Debye to a.u. # (N-1,3)
            os.chdir("../")
    
    else:
        os.chdir("GS_NEW/")
        sp.call( "grep 'Dipole moment (field-independent basis, Debye):' geometry.out -A 1 | tail -n 1 | awk '{print $2, $4, $6}' > DIP_TMP.dat" ,shell=True)
        TRANS_DIPOLES[0,:] = np.loadtxt("DIP_TMP.dat") / 2.5 # Debye to a.u. # (1,3)
        os.chdir("../")
        os.chdir("TD_NEW_S1/")
        sp.call( "grep 'Ground to excited state transition electric dipole moments' geometry.out -A %d | tail -n %d | awk '{print $2, $3, $4}' > DIP_TMP.dat" % (NStates,NStates-1), shell=True)
        TRANS_DIPOLES[1:,:] = np.loadtxt("DIP_TMP.dat") / 2.5 # Debye to a.u. # (N-1,3)
        os.chdir("../")
    

    return TRANS_DIPOLES


def main(DYN_PROPERTIES):
    
    NStates = DYN_PROPERTIES["NStates"]
    BOMD    = DYN_PROPERTIES["BOMD"]
    ISTATE  = DYN_PROPERTIES["ISTATE"]

    TRANS_DIPOLES = np.zeros(( NStates, 3 )) # Diagonal gradients
    TRANS_DIPOLES = read_Dipoles(TRANS_DIPOLES,NStates,BOMD,ISTATE)

    if ( DYN_PROPERTIES["MD_STEP"] >= 1 ):
        DYN_PROPERTIES["TRANS_DIPOLES_OLD"] = DYN_PROPERTIES["TRANS_DIPOLES_NEW"]
    DYN_PROPERTIES["TRANS_DIPOLES_NEW"] = TRANS_DIPOLES

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