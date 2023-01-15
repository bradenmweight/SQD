import numpy as np
import subprocess as sp
import os

def read_Energies(DIAG_ENERGIES,NStates):

    if ( NStates == 1 ):
        os.chdir("GS_NEW/")
        DIAG_ENERGIES = float( sp.check_output( "grep 'SCF Done' geometry.out | tail -n 1" ,shell=True).split()[4] )
        #np.savetxt(f"../DIAG_ENERGIES.dat", np.array([DIAG_ENERGIES]), header=f"Diagonal Energies (a.u.)", fmt="%2.8f" )

    else:
        os.chdir("TD_NEW_S1/")
        DIAG_ENERGIES[0] = float( sp.check_output( "grep 'SCF Done' geometry.out | tail -n 1" ,shell=True).split()[4] )
    
    if ( NStates >= 2 ):
        tmp = sp.check_output( "grep 'Excited State' geometry.out" ,shell=True).split(b"\n")
        for count,line in enumerate(tmp):
            if ( count+1 >= NStates ):
                break
            DIAG_ENERGIES[count+1] = float( line.split()[4] )/27.2114 + DIAG_ENERGIES[0]
        #np.savetxt(f"../DIAG_ENERGIES.dat", DIAG_ENERGIES, header=f"Diagonal Energies (a.u.)", fmt="%2.8f" )
    
    os.chdir("../")

    return DIAG_ENERGIES


def main(DYN_PROPERTIES):
    
    NStates = DYN_PROPERTIES["NStates"]

    DIAG_ENERGIES = np.zeros(( NStates )) # Diagonal gradients
    DIAG_ENERGIES = read_Energies(DIAG_ENERGIES,NStates)

    if ( DYN_PROPERTIES["MD_STEP"] >= 1 ):
        DYN_PROPERTIES["DIAG_ENERGIES_OLD"] = DYN_PROPERTIES["DIAG_ENERGIES_NEW"]
    DYN_PROPERTIES["DIAG_ENERGIES_NEW"] = DIAG_ENERGIES

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