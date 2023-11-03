import numpy as np
import subprocess as sp
import os

def read_Energies(DIAG_ENERGIES,NStates,BOMD,ISTATE):

    if ( BOMD == True ):
        if ( ISTATE == 0 and NStates == 1 ):
            # Read GS
            os.chdir(f"GS_NEW/")
            DIAG_ENERGIES[0] = float( sp.check_output( "grep 'SCF Done' geometry.out | tail -n 1" ,shell=True).split()[4] )
            os.chdir("../")
        elif ( ISTATE == 0 and NStates >= 2 ):
            # Read GS
            os.chdir(f"GS_NEW/")
            DIAG_ENERGIES[0] = float( sp.check_output( "grep 'SCF Done' geometry.out | tail -n 1" ,shell=True).split()[4] )
            os.chdir("../")

            # Read Excited States
            os.chdir(f"TD_NEW_S1/")
            tmp = sp.check_output( "grep 'Excited State' geometry.out" ,shell=True).split(b"\n")
            for count,line in enumerate(tmp):
                if ( count+1 >= NStates ):
                    break
                DIAG_ENERGIES[count+1] = float( line.split()[4] )/27.2114 + DIAG_ENERGIES[0]
            os.chdir("../")
        
        elif ( ISTATE != 0 and NStates >= 2 ):
            # Read Excited States
            os.chdir(f"TD_NEW_S{ISTATE}/")
            DIAG_ENERGIES[0] = float( sp.check_output( "grep 'SCF Done' geometry.out | tail -n 1" ,shell=True).split()[4] )
            tmp = sp.check_output( "grep 'Excited State' geometry.out" ,shell=True).split(b"\n")
            for count,line in enumerate(tmp):
                if ( count+1 >= NStates ):
                    break
                DIAG_ENERGIES[count+1] = float( line.split()[4] )/27.2114 + DIAG_ENERGIES[0]
            os.chdir("../")
    
    else:
        os.chdir("GS_NEW/")
        DIAG_ENERGIES[0] = float( sp.check_output( "grep 'SCF Done' geometry.out | tail -n 1" ,shell=True).split()[4] )
        os.chdir("../")
        os.chdir("TD_NEW_S1/")
        tmp = sp.check_output( "grep 'Excited State' geometry.out" ,shell=True).split(b"\n")
        for count,line in enumerate(tmp):
            if ( count+1 >= NStates ):
                break
            DIAG_ENERGIES[count+1] = float( line.split()[4] )/27.2114 + DIAG_ENERGIES[0]
        os.chdir("../")
    
    #print("Energy = ", DIAG_ENERGIES[:])
    

    return DIAG_ENERGIES


def check_state_labels(ENERGIES):
    NStates   = len(ENERGIES)
    ORDER_OLD = [ j for j in range(NStates) ]
    ORDER_NEW = [ j for j in range(NStates) ]
    for j in range( NStates ):
        for k in range( j+1, NStates ):
            if ( ENERGIES[j] > ENERGIES[k] ): # THIS WOULD BE WRONG, IF TRUE, SINCE j < k
                ENERGIES[j], ENERGIES[k] = ENERGIES[k], ENERGIES[j]
    return ENERGIES

def main(DYN_PROPERTIES):
    
    NStates = DYN_PROPERTIES["NStates"]
    BOMD    = DYN_PROPERTIES["BOMD"]
    ISTATE  = DYN_PROPERTIES["ISTATE"]

    DIAG_ENERGIES = np.zeros(( NStates )) # Diagonal gradients
    DIAG_ENERGIES = read_Energies(DIAG_ENERGIES,NStates,BOMD,ISTATE)

    DIAG_ENERGIES = check_state_labels( DIAG_ENERGIES )

    if ( DYN_PROPERTIES["MD_STEP"] >= 1 ):
        DYN_PROPERTIES["DIAG_ENERGIES_OLD"] = DYN_PROPERTIES["DIAG_ENERGIES_NEW"]
    DYN_PROPERTIES["DIAG_ENERGIES_NEW"] = DIAG_ENERGIES

    # #### FOR DEBUGGING ####
    # FILE01 = open(f"{DYN_PROPERTIES['SQD_RUNNING_DIR']}/AD_ENERGY.dat","a")
    # FILE01.write(f"{DYN_PROPERTIES['MD_STEP']}\t")
    # FILE01.write( "\t".join(map(str,DIAG_ENERGIES)) )
    # FILE01.write("\n")
    # FILE01.close()
    # #######################

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