import numpy as np
import subprocess as sp
import os, sys
import time
import multiprocessing as mp
from scipy.linalg import svd
import time

import get_cartesian_gradients
import get_diagonal_electronic_energies

def check_geometry(Atom_labels,Atom_coords_new):

    assert ( isinstance(Atom_labels, list) ), "Atoms labels needs to be a list" 
    assert ( isinstance(Atom_labels[0], str) ), "Atoms labels needs to be a list of strings"
    assert ( isinstance(Atom_coords_new, type(np.array([])) ) ) , "Atom coordinates need to be numpy array"

def clean_directory(NStates, MD_STEP, CPA_FLAG, BOMD_FLAG, ISTATE):

    if ( MD_STEP == 0 ): sp.call( "rm -rf GS* TD* DIMER OVERLAP" ,shell=True)

    if ( BOMD_FLAG ):

        if ( ISTATE == 0 and NStates == 1 ):
            # Remove old-old step
            sp.call( f"rm -rf GS_OLD" ,shell=True)

            # Move new to old
            sp.call( f"mv GS_NEW GS_OLD" ,shell=True)
            
            # Make fresh directory
            sp.call( f"mkdir GS_NEW" ,shell=True)

        elif ( ISTATE == 0 and NStates >= 2 ):
            # Remove old-old step
            sp.call( f"rm -rf GS_OLD" ,shell=True)
            states = [j for j in range(1,NStates)] * (not CPA_FLAG) + [1] * (CPA_FLAG)
            for state in states:
                sp.call( f"rm -rf TD_OLD_S{state}" ,shell=True)
            # Move new to old
            if ( MD_STEP >= 1 ):
                sp.call( "mv GS_NEW GS_OLD" ,shell=True)
                states = [j for j in range(1,NStates)] * (not CPA_FLAG) + [1] * (CPA_FLAG)
                for state in states:
                    sp.call( f"mv TD_NEW_S{state} TD_OLD_S{state}" ,shell=True)
            # Make fresh directory
            sp.call( "mkdir -p GS_NEW" ,shell=True)
            states = [j for j in range(1,NStates)] * (not CPA_FLAG) + [1] * (CPA_FLAG)
            for state in states:
                sp.call( f"mkdir -p TD_NEW_S{state}" ,shell=True)

        elif ( ISTATE != 0 and NStates >= 2 ):
            # Remove old-old step
            sp.call( f"rm -rf TD_OLD_S{ISTATE}" ,shell=True)
            states = [j for j in range(1,NStates)] * (not CPA_FLAG) + [1] * (CPA_FLAG)
            for state in states:
                sp.call( f"rm -rf TD_OLD_S{state}" ,shell=True)
            # Move new to old
            if ( MD_STEP >= 1 ): 
                states = [j for j in range(1,NStates)] * (not CPA_FLAG) + [1] * (CPA_FLAG)
                for state in states:
                    sp.call( f"mv TD_NEW_S{state} TD_OLD_S{state}" ,shell=True)
            # Make fresh directory
            states = [j for j in range(1,NStates)] * (not CPA_FLAG) + [1] * (CPA_FLAG)
            for state in states:
                sp.call( f"mkdir -p TD_NEW_S{state}" ,shell=True)

    else:

        # Remove old-old step
        sp.call( "rm -rf GS_OLD" ,shell=True)
        if ( NStates >= 2 ): 
            states = [j for j in range(1,NStates)] * (not CPA_FLAG) + [1] * (CPA_FLAG)
            for state in states:
                sp.call( f"rm -rf TD_OLD_S{state}" ,shell=True)
        if ( MD_STEP >= 2 and NStates >= 2 ): sp.call( "rm -rf DIMER" ,shell=True)

        # Move new to old
        if ( MD_STEP >= 1 ): 
            sp.call( "mv GS_NEW GS_OLD" ,shell=True)
            states = [j for j in range(1,NStates)] * (not CPA_FLAG) + [1] * (CPA_FLAG)
            for state in states:
                sp.call( f"mv TD_NEW_S{state} TD_OLD_S{state}" ,shell=True)

        # Make fresh directories
        sp.call( "mkdir -p GS_NEW" ,shell=True)
        if ( NStates >= 2 ): 
            states = [j for j in range(1,NStates)] * (not CPA_FLAG) + [1] * (CPA_FLAG)
            for state in states:
                sp.call( f"mkdir -p TD_NEW_S{state}" ,shell=True)
        if ( MD_STEP >= 1 and NStates >= 2 ): sp.call( "mkdir -p DIMER" ,shell=True)


def generate_inputs(DYN_PROPERTIES):

    NStates         = DYN_PROPERTIES["NStates"]
    Atom_labels     = DYN_PROPERTIES["Atom_labels"]
    Atom_coords_new = DYN_PROPERTIES["Atom_coords_new"]
    FUNCTIONAL      = DYN_PROPERTIES["FUNCTIONAL"]
    BASIS_SET       = DYN_PROPERTIES["BASIS_SET"]
    MEM             = DYN_PROPERTIES["MEMORY"]
    NCPUS_G16       = DYN_PROPERTIES["NCPUS_G16"]
    CHARGE          = DYN_PROPERTIES["CHARGE"]
    MULTIPLICITY    = DYN_PROPERTIES["MULTIPLICITY"]
    MD_STEP         = DYN_PROPERTIES["MD_STEP"]
    RUN_ELEC_STRUC  = DYN_PROPERTIES["RUN_ELEC_STRUC"]
    SBATCH_G16      = DYN_PROPERTIES["SBATCH_G16"]
    TDDFT_CONVERG   = DYN_PROPERTIES["TDDFT_CONVERG"]
    BOMD            = DYN_PROPERTIES["BOMD"]
    ISTATE          = DYN_PROPERTIES["ISTATE"]


    def write_header(file01,MEM,NCPUS_G16):
        file01.write(f"%chk=geometry.chk\n")
        file01.write(f"%nprocshared={NCPUS_G16}\n")
        file01.write(f"%mem={MEM}GB\n\n")

    def write_geom(file01, Atom_labels, Atom_coords_new, MD_STEP, CHARGE, MULTIPLICITY, method=[None,None]):
        file01.write(f"MD Step {MD_STEP}\n\n")
        file01.write(f"{CHARGE} {MULTIPLICITY}\n")
        for at in range( len(Atom_labels) ):
            file01.write( "%s  %2.8f  %2.8f  %2.8f \n" % (Atom_labels[at],Atom_coords_new[at,0]*0.529,Atom_coords_new[at,1]*0.529,Atom_coords_new[at,2]*0.529) )
        if ( method[0] == "DIMER" ):
            coords_dimer = method[1]
            for at in range( len(Atom_labels) ):
                file01.write( "%s  %2.8f  %2.8f  %2.8f \n" % (Atom_labels[at],coords_dimer[at,0]*0.529,coords_dimer[at,1]*0.529,coords_dimer[at,2]*0.529) )

    if ( BOMD == True and not (ISTATE == 0 and NStates >= 2) ):
        if ( ISTATE == 0 ):
            # Ground State for New Geometry
            os.chdir("GS_NEW/")
            file01 = open("geometry.com","w")
            if ( MD_STEP >= 1 ): 
                file01.write(f"%oldchk=../GS_OLD/geometry.chk\n")
            write_header(file01,MEM,NCPUS_G16)
            ### MAIN LINES ###  
            if ( FUNCTIONAL in ["DFTB", "DFTBA"] ):
                file01.write(f"# {FUNCTIONAL} SCF=XQC FORCE nosym ")
            else:
                file01.write(f"# {FUNCTIONAL}/{BASIS_SET} SCF=XQC FORCE nosym ")
            
            if ( MD_STEP >= 1 ):
                file01.write("guess=read\n\n")
            else:
                file01.write("\n\n")
            write_geom(file01,Atom_labels,Atom_coords_new,MD_STEP,CHARGE,MULTIPLICITY)
            if ( FUNCTIONAL in ["DFTB", "DFTBA"] ):
                file01.write(f"\n@GAUSS_EXEDIR:dftba.prm\n")
            # TODO ADD MIXED BASIS HERE
            file01.write("\n\n\n\n\n\n\n\n")
            os.chdir("../")


        elif ( ISTATE != 0 ):
            # Excited State for New Geometry (Use converged wavefunctions from GS at same geometry)
            os.chdir(f"TD_NEW_S{ISTATE}/")
            file01 = open("geometry.com","w")
            if ( MD_STEP >= 1 ):
                file01.write(f"%oldchk=../TD_OLD_S{ISTATE}/geometry.chk\n")
            write_header(file01,MEM,NCPUS_G16)
            # Add functional and basis set (unless DFTB or DFTBA)
            # Include one additional excited state even though we only perform NAMD on NStates-1 excied states to converge TD-DFT
            if ( FUNCTIONAL in ["DFTB", "DFTBA"] ):
                if ( MD_STEP == 0 ):
                    file01.write(f"# {FUNCTIONAL} SCF=XQC TD=(singlets,Conver={TDDFT_CONVERG},nstates={NStates},root={ISTATE}) FORCE nosym IOp(9/40=3)\n\n")
                if ( MD_STEP >= 1 ):
                    file01.write(f"# {FUNCTIONAL} SCF=XQC TD=(read,singlets,Conver={TDDFT_CONVERG},nstates={NStates},root={ISTATE}) FORCE nosym IOp(9/40=3) guess=read\n\n")
            else:
                if ( MD_STEP == 0 ):
                    file01.write(f"# {FUNCTIONAL}/{BASIS_SET} SCF=XQC TD=(singlets,Conver={TDDFT_CONVERG},nstates={NStates},root={ISTATE}) FORCE nosym IOp(9/40=3)\n\n")
                if ( MD_STEP >= 1 ):
                    file01.write(f"# {FUNCTIONAL}/{BASIS_SET} SCF=XQC TD=(read,singlets,Conver={TDDFT_CONVERG},nstates={NStates},root={ISTATE}) FORCE nosym IOp(9/40=3) guess=read\n\n")
            
            write_geom(file01,Atom_labels,Atom_coords_new,MD_STEP,CHARGE,MULTIPLICITY)
            if ( FUNCTIONAL in ["DFTB", "DFTBA"] ):
                file01.write(f"\n@GAUSS_EXEDIR:dftba.prm\n")
            # TODO ADD MIXED BASIS HERE
            file01.write("\n\n\n\n\n\n\n\n")
            os.chdir("../")



    else:

        # Ground State for New Geometry
        os.chdir("GS_NEW/")
        file01 = open("geometry.com","w")
        if ( MD_STEP >= 1 ): 
            file01.write(f"%oldchk=../GS_OLD/geometry.chk\n")
        write_header(file01,MEM,NCPUS_G16)
        ### MAIN LINES ###  
        if ( FUNCTIONAL in ["DFTB", "DFTBA"] ):
            file01.write(f"# {FUNCTIONAL} SCF=XQC FORCE nosym ")
        else:
            file01.write(f"# {FUNCTIONAL}/{BASIS_SET} SCF=XQC FORCE nosym ")
        
        if ( MD_STEP >= 1 ):
            file01.write("guess=read\n\n")
        else:
            file01.write("\n\n")
        write_geom(file01,Atom_labels,Atom_coords_new,MD_STEP,CHARGE,MULTIPLICITY)
        if ( FUNCTIONAL in ["DFTB", "DFTBA"] ):
            file01.write(f"\n@GAUSS_EXEDIR:dftba.prm\n")
        # TODO ADD MIXED BASIS HERE
        file01.write("\n\n\n\n\n\n\n\n")
        os.chdir("../")

        # Excited State for New Geometry (Use converged wavefunctions from GS at same geometry)
        if ( NStates >= 2 ):
            os.chdir("TD_NEW_S1/")
            file01 = open("geometry.com","w")
            if ( MD_STEP == 0 ):
                file01.write(f"%oldchk=../GS_NEW/geometry.chk\n")
            elif ( MD_STEP >= 1 ):
                file01.write(f"%oldchk=../TD_OLD_S1/geometry.chk\n")
            write_header(file01,MEM,NCPUS_G16)
            # Add functional and basis set (unless DFTB or DFTBA)
            # Include one additional excited state even though we only perform NAMD on NStates-1 excied states to converge TD-DFT
            if ( DYN_PROPERTIES["CPA"] == True ):
                if ( FUNCTIONAL in ["DFTB", "DFTBA"] ):
                    if ( MD_STEP == 0 ):
                        file01.write(f"# {FUNCTIONAL} SCF=XQC TD=(singlets,nstates={NStates},root=1) nosym IOp(9/40=3) guess=read\n\n")
                    elif ( MD_STEP >= 1 ):
                        file01.write(f"# {FUNCTIONAL} SCF=XQC TD=(read,singlets,nstates={NStates},root=1) nosym IOp(9/40=3) guess=read\n\n")
                else:
                    if ( MD_STEP == 0 ):
                        file01.write(f"# {FUNCTIONAL}/{BASIS_SET} SCF=XQC TD=(singlets,nstates={NStates},root=1) nosym IOp(9/40=3) guess=read\n\n")
                    elif ( MD_STEP >= 1 ):
                        file01.write(f"# {FUNCTIONAL}/{BASIS_SET} SCF=XQC TD=(read,singlets,nstates={NStates},root=1) nosym IOp(9/40=3) guess=read\n\n")
                
            else:
                if ( FUNCTIONAL in ["DFTB", "DFTBA"] ):
                    if ( MD_STEP == 0 ):
                        file01.write(f"# {FUNCTIONAL} SCF=XQC TD=(singlets,Conver={TDDFT_CONVERG},nstates={NStates},root=1) FORCE nosym IOp(9/40=3) guess=read\n\n")
                    if ( MD_STEP >= 1 ):
                        file01.write(f"# {FUNCTIONAL} SCF=XQC TD=(read,singlets,Conver={TDDFT_CONVERG},nstates={NStates},root=1) FORCE nosym IOp(9/40=3) guess=read\n\n")
                else:
                    if ( MD_STEP == 0 ):
                        file01.write(f"# {FUNCTIONAL}/{BASIS_SET} SCF=XQC TD=(singlets,Conver={TDDFT_CONVERG},nstates={NStates},root=1) FORCE nosym IOp(9/40=3) guess=read\n\n")
                    if ( MD_STEP >= 1 ):
                        file01.write(f"# {FUNCTIONAL}/{BASIS_SET} SCF=XQC TD=(read,singlets,Conver={TDDFT_CONVERG},nstates={NStates},root=1) FORCE nosym IOp(9/40=3) guess=read\n\n")
                
            write_geom(file01,Atom_labels,Atom_coords_new,MD_STEP,CHARGE,MULTIPLICITY)
            if ( FUNCTIONAL in ["DFTB", "DFTBA"] ):
                file01.write(f"\n@GAUSS_EXEDIR:dftba.prm\n")
            # TODO ADD MIXED BASIS HERE
            file01.write("\n\n\n\n\n\n\n\n")
            os.chdir("../")

        if ( NStates >= 3 and DYN_PROPERTIES["CPA"] == False ):
            for state in range( 2, NStates ):
                # Excited State for New Geometry (Use converged wavefunctions from TD root=1)
                os.chdir(f"TD_NEW_S{state}/")
                file01 = open("geometry.com","w")
                #file01.write(f"%oldchk=../TD_NEW_S1/geometry.chk\n")
                if ( MD_STEP == 0 ):
                    file01.write(f"%oldchk=../GS_NEW/geometry.chk\n")
                elif ( MD_STEP >= 1 ):
                    file01.write(f"%oldchk=../TD_OLD_S1/geometry.chk\n")
                write_header(file01,MEM,NCPUS_G16)
                # Add functional and basis set (unless DFTB or DFTBA)
                if ( FUNCTIONAL in ["DFTB", "DFTBA"] ):
                    if ( MD_STEP == 0 ):
                        file01.write(f"# {FUNCTIONAL} SCF=XQC TD=(singlets,Conver={TDDFT_CONVERG},nstates={NStates},root={state}) FORCE nosym IOp(9/40=3) guess=read\n\n")
                    elif ( MD_STEP >= 1 ):
                        file01.write(f"# {FUNCTIONAL} SCF=XQC TD=(read,singlets,Conver={TDDFT_CONVERG},nstates={NStates},root={state}) FORCE nosym IOp(9/40=3) guess=read\n\n")
                else:
                    if ( MD_STEP == 0 ):
                        file01.write(f"# {FUNCTIONAL}/{BASIS_SET} SCF=XQC TD=(singlets,Conver={TDDFT_CONVERG},nstates={NStates},root={state}) FORCE nosym IOp(9/40=3) guess=read\n\n")
                    elif ( MD_STEP >= 1 ):
                        file01.write(f"# {FUNCTIONAL}/{BASIS_SET} SCF=XQC TD=(read,singlets,Conver={TDDFT_CONVERG},nstates={NStates},root={state}) FORCE nosym IOp(9/40=3) guess=read\n\n")

                write_geom(file01,Atom_labels,Atom_coords_new,MD_STEP,CHARGE,MULTIPLICITY)
                if ( FUNCTIONAL in ["DFTB", "DFTBA"] ):
                    file01.write(f"\n@GAUSS_EXEDIR:dftba.prm\n")
                # TODO ADD MIXED BASIS HERE
                file01.write("\n\n\n\n\n\n\n\n")
                os.chdir("../")
        
        # Dimer method for atomic orbital overlaps
        if ( MD_STEP >= 1 and NStates >= 2 and BOMD == False ):
            Atom_coords_old = DYN_PROPERTIES["Atom_coords_old"]
            os.chdir("DIMER/")
            file01 = open("geometry.com","w")
            write_header(file01,MEM,NCPUS_G16)
            ### MAIN LINES ###
            # Add functional and basis set (unless DFTB or DFTBA)
            if ( FUNCTIONAL in ['AM1', 'PM3', 'PM6', 'PM7'] ): # By default, gaussian does not compute AO overlap for semi-empirical Hamiltonians
                file01.write(f"# {FUNCTIONAL}/{BASIS_SET} IOp(3/41=2000) nosymm iop(2/12=3,3/33=1) guess=only\n\n") ### MAIN LINE ###
            elif( FUNCTIONAL in ["DFTB", "DFTBA"] ):
                file01.write(f"# {FUNCTIONAL} IOp(3/41=2000) nosymm iop(2/12=3,3/33=1) guess=only\n\n")
            else:
                file01.write(f"# {FUNCTIONAL}/{BASIS_SET} nosymm iop(2/12=3,3/33=1) guess=only\n\n") ### MAIN LINE ###
            
            write_geom(file01,Atom_labels,Atom_coords_old,MD_STEP,CHARGE,MULTIPLICITY,method=["DIMER",Atom_coords_new])
            if ( DYN_PROPERTIES["FUNCTIONAL"] in ["DFTB", "DFTBA"] ):
                file01.write(f"\n@GAUSS_EXEDIR:dftba.prm\n")
            # TODO ADD MIXED BASIS HERE
            file01.write("\n\n\n\n\n\n\n\n")
            os.chdir("../")


def submit(RUN_ELEC_STRUC, SBATCH_G16, MD_STEP, BOMD, ISTATE, directory=None):
    if ( RUN_ELEC_STRUC == "SUBMIT_SBATCH" ):
        sp.call(f"cp {SBATCH_G16} .", shell=True)
        sp.call(f"sbatch {SBATCH_G16.split('/')[-1]}", shell=True)
        t0 = time.time()
        sleep_time = 0
        sleep_limit = 60 * 20 # 20 minutes, if electronic structure takes longer, we should not do NAMD
        sleep_check = 1 # Check gaussian output every {} seconds
        while ( True ):
            time.sleep(sleep_check)
            sleep_time += sleep_check # Add 1 seconds to sleep timer
            if ( sleep_limit > sleep_limit ):
                print(f"\tWARNING! Gaussian did not finish normally in the following directory:\n{os.getcwd()}", )
                exit()
            try:
                check1 = open("geometry.out","r").readlines()[-1]
                check2 = open("geometry.out","r").readlines()[-4]
            except (FileNotFoundError, IndexError) as errors:
                continue
            if ( check1.split()[:4] == "Normal termination of Gaussian".split() ):
                print("\tGaussian terminated normally in %2.2f s. (%s)" % (time.time() - t0, os.getcwd().split("/")[-1]) )
                sp.call(f"formchk geometry.chk > /dev/null 2>&1", shell=True)
                break
            elif ( check2.split()[:2] == "Error termination".split() ):
                print("\tGaussian crashed after %2.2f s. (%s)" % (time.time() - t0, os.getcwd().split("/")[-1]) )
                error = open("geometry.out","r").readlines()[-5] # Is this where all errors can be found ?
                print("Looking for possible error:\n", error)


    elif( RUN_ELEC_STRUC == "USE_CURRENT_NODE" ):
        t0 = time.time()
        sp.call(f"g16 < geometry.com > geometry.out", shell=True, stdout=True)
        #sp.call(f"g16 < geometry.com | tee geometry.out", shell=True, stdout=True)
        check = open("geometry.out","r").readlines()[-1]
        if ( check.split()[:4] == "Normal termination of Gaussian".split() ):
            print("\tGaussian terminated normally in %2.2f s. (%s)" % (time.time() - t0, os.getcwd().split("/")[-1]) )
            sp.call(f"formchk geometry.chk > /dev/null 2>&1", shell=True)
        else:
            print(f"\tWARNING! Gaussian did not finish normally in the following directory:\n{os.getcwd()}", )
    elif ( RUN_ELEC_STRUC == "TEST" ):
        # This is a test. Do nothing 
        print(f"Testing. I will not submit/run electronic structure calculations for step {MD_STEP}.")
    else:
        print(f"Error: 'RUN_ELEC_STRUC' was set to '{RUN_ELEC_STRUC}'. Not sure what to do.")


def run_ES_FORCE_parallel( inputs ):
    state, RUN_ELEC_STRUC, SBATCH_G16, MD_STEP, BOMD, ISTATE = inputs
    print(f"Starting forces for state {state}")
    os.chdir(f"TD_NEW_S{state}/")
    submit(RUN_ELEC_STRUC, SBATCH_G16, MD_STEP, BOMD, ISTATE)
    os.chdir("../")

def run_ES_FORCE_serial( NStates, RUN_ELEC_STRUC, SBATCH_G16, MD_STEP, BOMD, ISTATE, CPA_FLAG ):
    
    for state in range( 1, NStates ): # ALL STATES IN SERIAL
        if ( state >= 2 and CPA_FLAG == True ):
            continue
        if ( BOMD == True ):
            if ( state == ISTATE ):
                os.chdir(f"TD_NEW_S{state}/")
                submit(RUN_ELEC_STRUC, SBATCH_G16, MD_STEP, BOMD, ISTATE)
                os.chdir("../")
        else:
            os.chdir(f"TD_NEW_S{state}/")
            submit(RUN_ELEC_STRUC, SBATCH_G16, MD_STEP, BOMD, ISTATE)
            os.chdir("../")

def submit_jobs(DYN_PROPERTIES):
    NStates         = DYN_PROPERTIES["NStates"]
    RUN_ELEC_STRUC  = DYN_PROPERTIES["RUN_ELEC_STRUC"]
    SBATCH_G16      = DYN_PROPERTIES["SBATCH_G16"]
    MD_STEP         = DYN_PROPERTIES["MD_STEP"]
    BOMD            = DYN_PROPERTIES["BOMD"]
    ISTATE          = DYN_PROPERTIES["ISTATE"]

    print(f"Submitting electronic structure for step {MD_STEP}.")

    if ( BOMD == True and not (ISTATE == 0 and NStates >= 2) ):
        if ( ISTATE == 0 ):
            os.chdir("GS_NEW/")
            submit(RUN_ELEC_STRUC, SBATCH_G16, MD_STEP, BOMD, ISTATE)
            os.chdir("../")
        if ( ISTATE != 0 ):
            run_ES_FORCE_serial( NStates, RUN_ELEC_STRUC, SBATCH_G16, MD_STEP, BOMD, ISTATE, DYN_PROPERTIES["CPA"] )
    else:

        # GS must a serial job
        os.chdir("GS_NEW/")
        submit(RUN_ELEC_STRUC, SBATCH_G16, MD_STEP, BOMD, ISTATE)
        os.chdir("../")

        # Excited State
        if ( DYN_PROPERTIES["PARALLEL_FORCES"] == True and DYN_PROPERTIES["CPA"] == False ):
            state_List = [] 
            #for state in range( 2, NStates ): # Skip force for final excited state. We don't include in NAMD.
            for state in range( 1, NStates ): # Skip force for final excited state. We don't include in NAMD.
                state_List.append([ state, RUN_ELEC_STRUC, SBATCH_G16, MD_STEP, BOMD, ISTATE ])
            with mp.Pool(processes=DYN_PROPERTIES["NCPUS_NAMD"]) as pool:
                pool.map(run_ES_FORCE_parallel,state_List)
        else:
            run_ES_FORCE_serial( NStates, RUN_ELEC_STRUC, SBATCH_G16, MD_STEP, BOMD, ISTATE, DYN_PROPERTIES["CPA"] )

        # DIMER
        if ( MD_STEP >= 1 and NStates >= 2 and BOMD == False ):
            os.chdir("DIMER/")
            submit(RUN_ELEC_STRUC, SBATCH_G16, MD_STEP, BOMD, ISTATE)
            os.chdir("../")







def correct_phase(OVERLAP_ORTHO,MD_STEP,S_OLD=None):
    # Alexey Akimov, J. Phys. Chem. Lett. 2018, 9, 20, 6096–6102
    NStates = len(OVERLAP_ORTHO)
    f = np.zeros(( NStates ))
    for state in range( NStates ):
        f[state] = OVERLAP_ORTHO[state,state] / np.abs(OVERLAP_ORTHO[state,state]) # \pm 1
        #print("Phase Factor f(S%d) = %2.6f" % (state,f[state]))

    OVERLAP_corrected = np.zeros(( NStates, NStates ))
    for j in range( NStates ):
        for k in range( NStates ):
            OVERLAP_corrected[j,k] = OVERLAP_ORTHO[j,k] * f[k] # or f[j] ? # Here, f* = f since all wavefunctions are real-valued

    if ( MD_STEP >= 3 ):
        for j in range( NStates ):
            for k in range( j+1, NStates ):
                tmp1 = OVERLAP_corrected[j,k]
                tmp2 = np.abs(OVERLAP_corrected[j,k])
                tmp3 = S_OLD[j,k]
                tmp4 = np.abs(S_OLD[j,k])
                if ( np.isclose(tmp2,0.0) == False and np.isclose(tmp4,0.0) == False and int(tmp1/tmp2) != int(tmp3/tmp4) ):
                    if ( abs(OVERLAP_corrected[j,k] - S_OLD[j,k]) > 1e-3 ): # This is arbitrary threshold.
                        print( "OLD:\n", np.round(OVERLAP_corrected,8) )
                        OVERLAP_corrected[j,k] *= -1
                        OVERLAP_corrected[k,j] *= -1
                        print("I found a sign error upon phase correcting.")
                        print(f"Index Flip: {j} <--> {k}")
                        print( "NEW:\n", np.round(OVERLAP_corrected,8) )
                        print("[S.T @ S] after sign fix:")
                        print(np.round(OVERLAP_corrected.T @ OVERLAP_corrected,8))
                        OVERLAP_corrected = get_Lowdin_SVD(OVERLAP_corrected)
                        print("[S.T @ S] after sign fix and second orthogonalization:")
                        print(np.round(OVERLAP_corrected.T @ OVERLAP_corrected,8))


                

    #if ( MD_STEP >= 3 ):
    #    print("S(t)")
    #    print(S_OLD)
    #    print("S(t+dt)")
    #    print(OVERLAP_corrected)
    #    print("S(t)^-1 @ S(t+dt):")
    #    print( np.linalg.inv(S_OLD) @ OVERLAP_corrected )


    """
    if ( MD_STEP >= 3 ):
        for j in range( NStates ):
            for k in range( NStates ):
                if ( int(f[k]) == -1 and j != k  ):
                    sign_old = S_OLD[j,k] / np.abs(S_OLD[j,k])
                    sign_new = OVERLAP_corrected[j,k] / np.abs(OVERLAP_corrected[j,k])
                    if ( sign_old != sign_new and abs(OVERLAP_corrected[j,k]-S_OLD[j,k]) > 1e-5 ):
                        print(f"I found a sign issue in OVERLAP.") 
                        print(f"The sign changed from {sign_old} to {sign_new} with a large magnitude of {abs(OVERLAP_corrected[j,k]-S_OLD[j,k])}")
                        print(f"I changed the sign of index {j}-{k} to fix at MD step {MD_STEP}.")
                        OVERLAP_corrected[j,k] *= -1
            
            #print("Phase Corrected S_(%d-%d) = %2.8f" % (j,k,OVERLAP_corrected[j,k]))
            #print("Phase Non-corrected S_(%d-%d) = %2.8f " % (j,k,OVERLAP_ORTHO[j,k]))
    """

    return OVERLAP_corrected

def get_Lowdin_SVD(OVERLAP):
    """
    S = U @ diag(\lambda_i) @ V.T
    S_Ortho = U @ V.T
    """
    U, vals, VT = svd(OVERLAP)

    S_Ortho = U @ VT

    print("Orthogonalized [S.T @ S] :")
    print(np.round(S_Ortho.T @ S_Ortho,8))
    #print("Check orthogonalization. S.T @ S")
    #print("Saving to ortho_check.dat")
    #np.savetxt("ortho_check.dat", S_Ortho.T @ S_Ortho, fmt="%1.8f" )
    # It does work. ~BMW

    return S_Ortho

def calc_NACT(DYN_PROPERTIES):
    """
    Hammes-Schiffer and Tully, J. Chem. Phys., 101, 6, 1994
    NACT_{jk} \\approx (<j(t0)|k(t1)> - <j(t1)|k(t0)>) / (2*dt)
    """
    OVERLAP = DYN_PROPERTIES["OVERLAP_NEW"]
    dtI     = DYN_PROPERTIES["dtI"]

    #print("Original Overlap")
    #print(OVERLAP)

    OVERLAP_ORTHO = get_Lowdin_SVD(OVERLAP) * 1.0
    DYN_PROPERTIES["OVERLAP_NEW_uncorrected"] = OVERLAP_ORTHO * 1.0

    #print("Orthogonalized OVERLAP:")
    #print(OVERLAP_ORTHO)

    # Save uncorrected properties for debugging purposes. Remove later.
    NACT_uncorrected = (OVERLAP_ORTHO - OVERLAP_ORTHO.T) / 2 / dtI
    DYN_PROPERTIES["NACT_NEW_uncorrected"] = NACT_uncorrected * 1.0
    
    # Correct the phase of the OVERLAP matrix
    if ( DYN_PROPERTIES["MD_STEP"] >= 3 ):
        OVERLAP_CORR = correct_phase(OVERLAP_ORTHO,DYN_PROPERTIES["MD_STEP"],S_OLD=DYN_PROPERTIES["OVERLAP_OLD"]) * 1.0
    else:
        OVERLAP_CORR = correct_phase(OVERLAP_ORTHO,DYN_PROPERTIES["MD_STEP"]) * 1.0

    #print("Phase-corrected OVERLAP:")
    #print(OVERLAP_CORR)

    if ( DYN_PROPERTIES["CHECK_TRIVIAL_CROSSING"] == True ):
        OVERLAP_CORR, DYN_PROPERTIES = check_for_trivial_crossing(OVERLAP_CORR, DYN_PROPERTIES)

    NACT = (OVERLAP_CORR - OVERLAP_CORR.T) / 2 / dtI



    if ( DYN_PROPERTIES["MD_STEP"] >= 2 ):
        DYN_PROPERTIES["NACT_OLD"] = DYN_PROPERTIES["NACT_NEW"] * 1.0
        DYN_PROPERTIES["OVERLAP_OLD"] = DYN_PROPERTIES["OVERLAP_NEW"] * 1.0
    DYN_PROPERTIES["NACT_NEW"] = NACT * 1.0
    DYN_PROPERTIES["OVERLAP_NEW"] = OVERLAP_CORR * 1.0
    
    return DYN_PROPERTIES


def get_approx_NACR( DYN_PROPERTIES ):
    """
    Shu, ..., Truhlar, J. Chem. Theory Comput. 2022, 18, 3, 1320-1328
    """

    NAtoms = DYN_PROPERTIES["NAtoms"]
    NStates = DYN_PROPERTIES["NStates"]
    V     = DYN_PROPERTIES["Atom_velocs_new"]
    Ead   = DYN_PROPERTIES["DIAG_ENERGIES_NEW"]
    dEad  = DYN_PROPERTIES["DIAG_GRADIENTS"]
    NACT  = DYN_PROPERTIES["NACT_NEW"]

    alpha = np.zeros(( NStates, NStates ))
    g = np.zeros(( NStates, NStates, NAtoms, 3 ))
    G = np.zeros(( NStates, NStates, NAtoms, 3 ))
    
    if ( np.allclose(V, np.zeros((NAtoms,3))) ):
        print("\t Velocities are ZERO. Skipping approximate NACR calculation.")
        print("\t Setting NACR to zero. If this happens at time-steps later than 0, something is wrong.")
        if ( DYN_PROPERTIES["MD_STEP"] >= 2 ):
            DYN_PROPERTIES["NACR_APPROX_OLD"] = DYN_PROPERTIES["NACR_APPROX_NEW"]
        DYN_PROPERTIES["NACR_APPROX_NEW"] = np.zeros(( NStates, NStates, NAtoms, 3 ))
        return DYN_PROPERTIES

    for j in range( NStates ):
        for k in range( NStates ):
            g[j,k,:,:] = dEad[j,:,:] - dEad[k,:,:]
            alpha[j,k] = NACT[j,k] - np.einsum("ad,ad->", V[:,:], g[j,k,:,:])
            alpha[j,k] /= np.einsum("ad,ad->", V[:,:], V[:,:] )
            G[j,k,:,:] = g[j,k,:,:] + alpha[j,k] * V[:,:]

    if ( DYN_PROPERTIES["MD_STEP"] >= 2 ):
        DYN_PROPERTIES["NACR_APPROX_OLD"] = DYN_PROPERTIES["NACR_APPROX_NEW"]
    DYN_PROPERTIES["NACR_APPROX_NEW"] = G[:,:,:,:]
    
    return DYN_PROPERTIES


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
                for key, value in DYN_PROPERTIES.items() :
                    print (key)

    return OVERLAP_CORR, DYN_PROPERTIES


                
def main(DYN_PROPERTIES):

    # HOW TO HANDLE THIS PATH BETTER ?
    sys.path.append(f"{DYN_PROPERTIES['SQD_HOME_PATH']}/src/WFN_OVERLAP/PYTHON/")
    import G16_NAC

    NStates         = DYN_PROPERTIES["NStates"] # Total number of electronic states
    NAtoms          = DYN_PROPERTIES["NAtoms"] # Total number of electronic states
    Atom_labels     = DYN_PROPERTIES["Atom_labels"]
    Atom_coords_new = DYN_PROPERTIES["Atom_coords_new"]
    MD_STEP         = DYN_PROPERTIES["MD_STEP"]
    BOMD            = DYN_PROPERTIES["BOMD"]
    ISTATE          = DYN_PROPERTIES["ISTATE"]
    
    if ( not os.path.exists("G16") ):
        sp.call("mkdir G16", shell=True)
    os.chdir("G16")
    
    check_geometry(Atom_labels,Atom_coords_new)
    clean_directory(NStates,MD_STEP,DYN_PROPERTIES["CPA"],DYN_PROPERTIES["BOMD"],ISTATE)
    generate_inputs(DYN_PROPERTIES)
    submit_jobs(DYN_PROPERTIES)


    DYN_PROPERTIES = get_cartesian_gradients.main(DYN_PROPERTIES)
    DYN_PROPERTIES = get_diagonal_electronic_energies.main(DYN_PROPERTIES)
    if ( MD_STEP >= 1 ):
        if ( NStates >= 2 and BOMD == False ):
            T0 = time.time()
            DYN_PROPERTIES = G16_NAC.main(DYN_PROPERTIES) # Provides OVERLAP
            print( f"\tGET_OVERLAP OVERALL TIME (G16_TD.py):", round(time.time() - T0,2), "s" )
            T0 = time.time()
            DYN_PROPERTIES = calc_NACT(DYN_PROPERTIES) # Provides NACT
            print( f"\tGET_NACT OVERALL TIME (G16_TD.py):", round(time.time() - T0,2), "s" )
            T0 = time.time()
            if ( DYN_PROPERTIES["CPA"] == False ): # Only need NACR for off-diagonal forces, which only come from G.S. in CPA
                DYN_PROPERTIES = get_approx_NACR(DYN_PROPERTIES) # Provides NACR from NACT and OVERLAP
            print( f"\tGET_APPROX NACR OVERALL TIME (G16_TD.py):", round(time.time() - T0,2), "s" )
        else:
            DYN_PROPERTIES["OVERLAP_OLD"] = np.zeros(( NStates, NStates ))
            DYN_PROPERTIES["OVERLAP_NEW"] = np.zeros(( NStates, NStates ))
            DYN_PROPERTIES["NACR_APPROX_OLD"] = np.zeros(( NStates, NStates, NAtoms, 3 ))
            DYN_PROPERTIES["NACR_APPROX_NEW"] = np.zeros(( NStates, NStates, NAtoms, 3 ))
            DYN_PROPERTIES["NACT_NEW"] = np.zeros(( NStates, NStates ))

    os.chdir("../")

    return DYN_PROPERTIES

def read_XYZ(filename):
    XYZ_File = open(filename,"r").readlines()
    NAtoms = int(XYZ_File[0])
    Atom_labels = []
    Atom_coords_new = np.zeros(( NAtoms, 3 ))
    for count, line in enumerate(XYZ_File[2:]):
        t = line.split()
        Atom_labels.append( t[0] )
        Atom_coords_new[count,:] = np.array([ float(t[1]), float(t[2]), float(t[3]) ])

    return Atom_labels, Atom_coords_new

if ( __name__ == "__main__" ):

    # THIS CODE NEEDS TO BE RUN TWICE WHEN USING THE EXAMPLE

    sp.call("module load gaussian", shell=True)
    sp.call("module load intel", shell=True)

    Atom_labels, Atom_coords_new = read_XYZ("geometry_new.xyz")
    Atom_labels, Atom_velocs_new = read_XYZ("velocities_new.xyz")

    DYN_PROPERTIES = {"Atom_labels":Atom_labels, "Atom_coords_new":Atom_coords_new }
    Atom_coords_old = Atom_coords_new * 1.0
    Atom_coords_old[0,0] += 0.1
    DYN_PROPERTIES["Atom_coords_old"] = DYN_PROPERTIES["Atom_coords_new"]
    
    DYN_PROPERTIES["Atom_velocs_new"] = DYN_PROPERTIES["Atom_coords_new"]

    DYN_PROPERTIES["NStates"] = 4
    DYN_PROPERTIES["NAtoms"] = len(DYN_PROPERTIES["Atom_labels"])
    DYN_PROPERTIES["FUNCTIONAL"] = "BLYP"
    DYN_PROPERTIES["CHARGE"] = 0
    DYN_PROPERTIES["MULTIPLICITY"] = 1
    DYN_PROPERTIES["BASIS_SET"] = "sto-3g"
    DYN_PROPERTIES["MEMORY"] = 5
    DYN_PROPERTIES["NCPUS"] = 1
    DYN_PROPERTIES["MD_STEP"] = 1
    DYN_PROPERTIES["dtI"] = 0.1 # fs
    DYN_PROPERTIES["RUN_ELEC_STRUC"] = "SUBMIT_SBATCH" # "SUBMIT_SBATCH", "USE_CURRENT_NODE", "TEST"
    DYN_PROPERTIES["SBATCH_G16"] = "/scratch/bweight/software/many_molecule_many_mode_NAMD/src/ELECTRONIC_STRUCTURE_CONTROL/EXAMPLE/submit.gaussian" # For "SUBMIT_SBATCH" in previous only
    DYN_PROPERTIES = main(DYN_PROPERTIES)
    print( "NACR_APPROX_NEW (S1/S2):" )
    print( (DYN_PROPERTIES["NACR_APPROX_NEW"])[1,2] )
