import numpy as np
import subprocess as sp
import os, sys
import time
import multiprocessing as mp
from scipy.linalg import svd
import time

import get_cartesian_gradients
import get_diagonal_electronic_energies

sys.path.append("/scratch/bweight/software/many_molecule_many_mode_NAMD/src/WFN_OVERLAP/PYTHON/")
import G16_NAC

def check_geometry(Atom_labels,Atom_coords_new):

    assert ( isinstance(Atom_labels, list) ), "Atoms labels needs to be a list" 
    assert ( isinstance(Atom_labels[0], str) ), "Atoms labels needs to be a list of strings"
    assert ( isinstance(Atom_coords_new, type(np.array([])) ) ) , "Atoms coordinates need to be numpy array"

def clean_directory(NStates, MD_STEP):

    sp.call( "rm -rf GS_OLD" ,shell=True)
    if ( NStates >= 2 ): sp.call( "rm -rf TD_OLD_S1" ,shell=True)
    if ( MD_STEP >= 2 and NStates >= 2 ): sp.call( "rm -rf DIMER" ,shell=True)
    if ( NStates >= 3 ):
        for state in range( 2, NStates ):
            sp.call( f"rm -rf TD_OLD_S{state}" ,shell=True)

    
    if ( MD_STEP >= 1 ): 
        sp.call( "mv GS_NEW GS_OLD" ,shell=True)
        if ( NStates >= 2 ): sp.call( "mv TD_NEW_S1 TD_OLD_S1" ,shell=True)
        if ( NStates >= 3 ):
            for state in range( 2, NStates ):
                sp.call( f"mv TD_NEW_S{state} TD_OLD_S{state}" ,shell=True)

    sp.call( "mkdir -p GS_NEW" ,shell=True)
    if ( NStates >= 2 ): sp.call( "mkdir -p TD_NEW_S1" ,shell=True)
    if ( MD_STEP >= 1 and NStates >= 2 ): sp.call( "mkdir -p DIMER" ,shell=True)

    if ( NStates >= 3 ):
        for state in range( 2, NStates ):
            sp.call( f"mkdir -p TD_NEW_S{state}" ,shell=True)

    # ADD DIRECTORY TO KEEP TRACK OF PREVIOUS PHASE INFORMATION. ADD LATER. # TODO


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

    #assert(FUNCTIONAL.upper() not in ['AM1', 'PM3', 'PM6', 'PM7'] ), "CI Overlap code does not work with semi-empirical Hamiltonians."
    # Let's keep semi-empirical for BOMD where we don't need to overlaps of electronic wavefunctions.
    # Does the overlap work if we orthogonalize the overlap ??? Can test empirically first...

    def write_header(file01,MEM,NCPUS_G16):
        file01.write(f"%chk=geometry.chk\n")
        file01.write(f"%nprocshared={NCPUS_G16}\n")#file01.write(f"%nprocshared=1\n")
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
        file01.write("\n\n\n\n\n\n\n\n")

    # Ground State for New Geometry
    os.chdir("GS_NEW/")
    file01 = open("geometry.com","w")
    if ( MD_STEP >= 1 ): 
        file01.write(f"%oldchk=../GS_OLD/geometry.chk\n")
    write_header(file01,MEM,NCPUS_G16)
    file01.write(f"# {FUNCTIONAL}/{BASIS_SET} SCF=XQC FORCE nosym pop=full ") ### MAIN LINE ###
    if ( MD_STEP >= 1 ):
        file01.write("guess=read\n\n")
    else:
        file01.write("\n\n")
    write_geom(file01,Atom_labels,Atom_coords_new,MD_STEP,CHARGE,MULTIPLICITY)
    os.chdir("../")

    # Excited State for New Geometry (Use converged wavefunctions from GS at same geometry)
    if ( NStates >= 2 ):
        os.chdir("TD_NEW_S1/")
        file01 = open("geometry.com","w")
        file01.write(f"%oldchk=../GS_NEW/geometry.chk\n")
        write_header(file01,MEM,NCPUS_G16)
        # Include one additional excited state even though we only perform NAMD on NStates-1 excied states to converge TD-DFT
        file01.write(f"# {FUNCTIONAL}/{BASIS_SET} SCF=XQC TD=(singlets,nstates={NStates},root=1) FORCE nosym pop=full guess=read\n\n") ### MAIN LINE ###
        write_geom(file01,Atom_labels,Atom_coords_new,MD_STEP,CHARGE,MULTIPLICITY)
        os.chdir("../")

    if ( NStates >= 3 ):
        for state in range( 2, NStates ):
            # Excited State for New Geometry (Use converged wavefunctions from TD root=1)
            os.chdir(f"TD_NEW_S{state}/")
            file01 = open("geometry.com","w")
            #file01.write(f"%oldchk=../TD_NEW_S1/geometry.chk\n")
            file01.write(f"%oldchk=../GS_NEW/geometry.chk\n")
            write_header(file01,MEM,NCPUS_G16)
            # Include one additional excited state even though we only perform NAMD on NStates-1 excied states to converge TD-DFT
            #file01.write(f"# {FUNCTIONAL}/{BASIS_SET} SCF=XQC TD=(read,singlets,nstates={NStates},root={state}) FORCE nosym pop=full guess=read\n\n") ### MAIN LINE ###
            #file01.write(f"# {FUNCTIONAL}/{BASIS_SET} SCF=XQC TD=(read,root={state}) FORCE nosym pop=full guess=read\n\n") ### MAIN LINE ###
            file01.write(f"# {FUNCTIONAL}/{BASIS_SET} SCF=XQC TD=(singlets,nstates={NStates},root={state}) FORCE nosym pop=full guess=read\n\n") ### MAIN LINE ###
            write_geom(file01,Atom_labels,Atom_coords_new,MD_STEP,CHARGE,MULTIPLICITY)
            os.chdir("../")
    
    # Dimer method for atomic orbital overlaps
    if ( MD_STEP >= 1 and NStates >= 2 ):
        Atom_coords_old = DYN_PROPERTIES["Atom_coords_old"]
        os.chdir("DIMER/")
        file01 = open("geometry.com","w")
        write_header(file01,MEM,NCPUS_G16)
        if ( FUNCTIONAL.upper() in ['AM1', 'PM3', 'PM6', 'PM7'] ): # By default, gaussian does not compute AO overlap for semi-empirical Hamiltonians
            # THIS GIVES BAD OVERLAPS STILL. DO NOT USE. NEED TO FIX.
            file01.write(f"# {FUNCTIONAL}/{BASIS_SET} IOp(3/41=2000) nosymm iop(2/12=3,3/33=1) guess=only pop=full\n\n") ### MAIN LINE ###
        else:
            file01.write(f"# {FUNCTIONAL}/{BASIS_SET} nosymm iop(2/12=3,3/33=1) guess=only pop=full\n\n") ### MAIN LINE ###
        write_geom(file01,Atom_labels,Atom_coords_old,MD_STEP,CHARGE,MULTIPLICITY,method=["DIMER",Atom_coords_new])
        os.chdir("../")


def submit(RUN_ELEC_STRUC, SBATCH_G16, MD_STEP, directory=None):
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
    state, RUN_ELEC_STRUC, SBATCH_G16, MD_STEP = inputs
    print(f"Starting forces for state {state}")
    os.chdir(f"TD_NEW_S{state}/")
    submit(RUN_ELEC_STRUC, SBATCH_G16, MD_STEP)
    os.chdir("../")

def run_ES_FORCE_serial( RUN_ELEC_STRUC, SBATCH_G16, MD_STEP ):
    #for state in range( 2, NStates ): # ONLY 2,NStates in PARALLEL
    for state in range( 1, NStates ): # ALL STATES IN PARALLEL
        os.chdir(f"TD_NEW_S{state}/")
        submit(RUN_ELEC_STRUC, SBATCH_G16, MD_STEP)
        os.chdir("../")

def submit_jobs(DYN_PROPERTIES):
    NStates         = DYN_PROPERTIES["NStates"]
    RUN_ELEC_STRUC  = DYN_PROPERTIES["RUN_ELEC_STRUC"]
    SBATCH_G16      = DYN_PROPERTIES["SBATCH_G16"]
    MD_STEP         = DYN_PROPERTIES["MD_STEP"]

    print(f"Submitting electronic structure for step {MD_STEP}.")

    # GS must a serial job
    os.chdir("GS_NEW/")
    submit(RUN_ELEC_STRUC, SBATCH_G16, MD_STEP)
    os.chdir("../")

    if ( DYN_PROPERTIES["PARALLEL_FORCES"] == True ):
        state_List = [] 
        #for state in range( 2, NStates ): # Skip force for final excited state. We don't include in NAMD.
        for state in range( 1, NStates ): # Skip force for final excited state. We don't include in NAMD.
            state_List.append([ state, RUN_ELEC_STRUC, SBATCH_G16, MD_STEP ])
        with mp.Pool(processes=DYN_PROPERTIES["NCPUS_NAMD"]) as pool:
            pool.map(run_ES_FORCE_parallel,state_List)
    else:
        run_ES_FORCE_serial( RUN_ELEC_STRUC, SBATCH_G16, MD_STEP )




    """ # THIS IS FOR PARALLELIZING NSTATES >= 3
    if ( NStates >= 2 ):
        os.chdir("TD_NEW_S1/")
        submit(RUN_ELEC_STRUC, SBATCH_G16, MD_STEP)
        os.chdir("../")

    ### ADD PARALLELIZATION HERE FOR COMPUTING ALL THE EXCITED STATE FORCES IF NSTATES >= 3 ###
    if ( NStates >= 3 ):
        if ( DYN_PROPERTIES["PARALLEL_FORCES"] == True ):
            state_List = [] 
            for state in range( 2, NStates ): # Skip force for final excited state. We don't include in NAMD.
                state_List.append([ state, RUN_ELEC_STRUC, SBATCH_G16, MD_STEP ])
            with mp.Pool(processes=DYN_PROPERTIES["NCPUS_NAMD"]) as pool:
                pool.map(run_ES_FORCE_parallel,state_List)
        else:
            run_ES_FORCE_serial( RUN_ELEC_STRUC, SBATCH_G16, MD_STEP )
    """


    if ( MD_STEP >= 1 and NStates >= 2 ):
        os.chdir("DIMER/")
        submit(RUN_ELEC_STRUC, SBATCH_G16, MD_STEP)
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
                if ( int( S_OLD[j,k]/np.abs(S_OLD[j,k]) ) != int(OVERLAP_corrected[j,k]/np.abs(OVERLAP_corrected[j,k])) ):
                    if ( abs(OVERLAP_corrected[j,k] - S_OLD[j,k]) > 1e-3 ): # This is arbitrary threshold.
                        OVERLAP_corrected[j,k] *= -1
                        OVERLAP_corrected[k,j] *= -1
                        print("I found a sign error upon phase correcting.")
                        print(f"Index Flip: {j} <--> {k}")

                

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
    NACT = (OVERLAP_CORR - OVERLAP_CORR.T) / 2 / dtI

    #print("Phase-corrected OVERLAP:")
    #print(OVERLAP_CORR)





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


def main(DYN_PROPERTIES):

    NStates = DYN_PROPERTIES["NStates"] # Total number of electronic states
    NAtoms = DYN_PROPERTIES["NAtoms"] # Total number of electronic states
    Atom_labels = DYN_PROPERTIES["Atom_labels"]
    Atom_coords_new = DYN_PROPERTIES["Atom_coords_new"]
    MD_STEP = DYN_PROPERTIES["MD_STEP"]
    
    if ( not os.path.exists("G16") ):
        sp.call("mkdir G16", shell=True)
    os.chdir("G16")
    
    check_geometry(Atom_labels,Atom_coords_new)
    clean_directory(NStates,MD_STEP)
    generate_inputs(DYN_PROPERTIES)
    submit_jobs(DYN_PROPERTIES)


    DYN_PROPERTIES = get_cartesian_gradients.main(DYN_PROPERTIES)
    DYN_PROPERTIES = get_diagonal_electronic_energies.main(DYN_PROPERTIES)
    if ( MD_STEP >= 1 ):
        if ( NStates >= 2 ):
            T0 = time.time()
            DYN_PROPERTIES = G16_NAC.main(DYN_PROPERTIES) # Provides OVERLAP
            print( f"\tGET_OVERLAP OVERALL TIME (G16_TD.py):", round(time.time() - T0,2), "s" )
            T0 = time.time()
            DYN_PROPERTIES = calc_NACT(DYN_PROPERTIES) # Provides NACT
            print( f"\tGET_NACT OVERALL TIME (G16_TD.py):", round(time.time() - T0,2), "s" )
            T0 = time.time()
            DYN_PROPERTIES = get_approx_NACR(DYN_PROPERTIES) # Provides NACR from NACT and OVERLAP
            print( f"\tGET_APPROX NACR OVERALL TIME (G16_TD.py):", round(time.time() - T0,2), "s" )
        else:
            DYN_PROPERTIES["OVERLAP_OLD"] = 0
            DYN_PROPERTIES["OVERLAP_NEW"] = 0
            DYN_PROPERTIES["NACR_APPROX_OLD"] = np.zeros(( NStates, NStates, NAtoms, 3 ))
            DYN_PROPERTIES["NACR_APPROX_NEW"] = np.zeros(( NStates, NStates, NAtoms, 3 ))

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
