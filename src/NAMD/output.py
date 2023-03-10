import numpy as np
import subprocess as sp

import properties

# MOVE THIS FUNCTION TO NEW FILE OUTPUT.py
def save_data(DYN_PROPERTIES):

    NStates = DYN_PROPERTIES["NStates"]

    if ( DYN_PROPERTIES["MD_STEP"] == 0 ): 
        sp.call("rm -rf MD_OUTPUT ",shell=True)
        sp.call("mkdir MD_OUTPUT ",shell=True)

    with open("MD_OUTPUT/trajectory.xyz","a") as file01:
        file01.write(f"{DYN_PROPERTIES['NAtoms']}\n")
        file01.write(f"MD Step {DYN_PROPERTIES['MD_STEP']} Units = [Angstroms]\n")
        Atom_labels = DYN_PROPERTIES["Atom_labels"]
        Atom_coords = DYN_PROPERTIES["Atom_coords_new"] * 0.529 # 0.529 Ang./Bohr
        for count, atom in enumerate( Atom_labels ):
            #file01.write(f"{atom}\t{Atom_coords[count,0]*0.529}\t{Atom_coords[count,1]*0.529}\t{Atom_coords[count,2]*0.529}\n")
            file01.write(f"{atom}  " + " ".join(map("{:2.8f}".format,Atom_coords[count,:]))  + "\n")

    with open("MD_OUTPUT/velocity.xyz","a") as file01:
        file01.write(f"{DYN_PROPERTIES['NAtoms']}\n")
        file01.write(f"MD Step {DYN_PROPERTIES['MD_STEP']} Units = [Angstroms / fs]\n")
        Atom_labels = DYN_PROPERTIES["Atom_labels"]
        Atom_velocs = DYN_PROPERTIES["Atom_velocs_new"] * 0.529 * 41.341 # 0.529 Ang./Bohr, 41.341 a.u. / fs
        for count, atom in enumerate( Atom_labels ):
            file01.write(f"{atom}  " + " ".join(map("{:2.8f}".format,Atom_velocs[count,:]))  + "\n") # Ang / fs

    with open("MD_OUTPUT/forces.xyz","a") as file01:
        file01.write(f"{DYN_PROPERTIES['NAtoms']}\n")
        file01.write(f"MD Step {DYN_PROPERTIES['MD_STEP']} Units = [eV / Ang.]\n")
        Atom_labels = DYN_PROPERTIES["Atom_labels"]
        Atom_forces = DYN_PROPERTIES["FORCE_NEW"] / 0.529 * 27.2114 # 0.529 Ang./Bohr, 27.2114 eV / Hartree
        for count, atom in enumerate( Atom_labels ):
            file01.write(f"{atom}  " + " ".join(map("{:2.8f}".format,Atom_forces[count,:]))  + "\n") # eV / Ang.


    with open("MD_OUTPUT/PES.dat","a") as file01:
        if ( NStates >= 2 ):
            file01.write( f"{DYN_PROPERTIES['MD_STEP']}  " +  " ".join(map("{:2.8f}".format,DYN_PROPERTIES["DIAG_ENERGIES_NEW"]*27.2114 )) + "\n" )
        else:
            file01.write( f"{DYN_PROPERTIES['MD_STEP']}  {np.round(DYN_PROPERTIES['DIAG_ENERGIES_NEW']*27.2114,8)}\n" )

    with open("MD_OUTPUT/Energy.dat","a") as file01:
        
        DYN_PROPERTIES = properties.compute_KE(DYN_PROPERTIES)
        DYN_PROPERTIES = properties.compute_PE(DYN_PROPERTIES)

        KE = DYN_PROPERTIES["KE"] * 27.2114
        PE = DYN_PROPERTIES["PE"] * 27.2114
        TE = KE + PE

        file01.write(f"{DYN_PROPERTIES['MD_STEP']}  " + "%2.6f  %2.6f  %2.6f\n" % (KE,PE,TE))

    with open("MD_OUTPUT/Temperature.dat","a") as file01:
        
        DYN_PROPERTIES = properties.compute_KE(DYN_PROPERTIES)
        NAtoms = DYN_PROPERTIES['NAtoms']

        KE = DYN_PROPERTIES["KE"] * 27.2114 # a.u. --> eV

        T = (2/3) * KE / NAtoms * (300 / 0.025) # eV --> K

        file01.write(f"{DYN_PROPERTIES['MD_STEP']}  " + "%2.4f\n" % (T))


    if ( NStates >= 2 ):

        with open("MD_OUTPUT/mapping_re.dat","a") as file01:
            file01.write( f"{DYN_PROPERTIES['MD_STEP']}  " +  " ".join(map("{:2.8f}".format,DYN_PROPERTIES["MAPPING_VARS"].real )) + "\n" )

        with open("MD_OUTPUT/mapping_im.dat","a") as file01:
            file01.write( f"{DYN_PROPERTIES['MD_STEP']}  " +  " ".join(map("{:2.8f}".format,DYN_PROPERTIES["MAPPING_VARS"].imag )) + "\n" )

        with open("MD_OUTPUT/Population.dat","a") as file01:
            POP = np.real(properties.get_density_matrix( DYN_PROPERTIES )[np.diag_indices(NStates)])
            PSUM = np.sum(POP)
            file01.write( f"{DYN_PROPERTIES['MD_STEP']}  " +  " ".join(map("{:2.8f}".format,POP )) + "  %2.8f" % (PSUM) + "\n" )

        with open("MD_OUTPUT/Coherence_re.dat","a") as file01:
            if ( DYN_PROPERTIES['MD_STEP'] == 0 ): 
                file01.write(f"# Step " + " ".join([f'{j}-{k}' for j in range(NStates) for k in range(j+1,NStates)]) + "\n" )
            RHO = np.real(properties.get_density_matrix( DYN_PROPERTIES ))
            COH = np.array([RHO[j,k] for j in range(NStates) for k in range(j+1,NStates)]).real
            file01.write( f"{DYN_PROPERTIES['MD_STEP']}  " +  " ".join(map("{:2.8f}".format,COH )) + "\n" )

        with open("MD_OUTPUT/Coherence_im.dat","a") as file01:
            if ( DYN_PROPERTIES['MD_STEP'] == 0 ): 
                file01.write(f"# Step " + " ".join([f'{j}-{k}' for j in range(NStates) for k in range(j+1,NStates)]) + "\n" )
            RHO = np.real(properties.get_density_matrix( DYN_PROPERTIES ))
            COH = np.array([RHO[j,k] for j in range(NStates) for k in range(j+1,NStates)]).imag
            file01.write( f"{DYN_PROPERTIES['MD_STEP']}  " +  " ".join(map("{:2.8f}".format,COH )) + "\n" )

        with open("MD_OUTPUT/Overlap.dat","a") as file01:
            if ( DYN_PROPERTIES['MD_STEP'] == 0 ): 
                file01.write(f"# Step " + " ".join([f'{j}-{k}' for j in range(NStates) for k in range(j,NStates)]) + "\n" )
            if ( DYN_PROPERTIES['MD_STEP'] >= 1 ): 
                OVERLAP = DYN_PROPERTIES['OVERLAP_NEW'] * 1.0
                OVERLAP = np.array([OVERLAP[j,k] for j in range(NStates) for k in range(j,NStates)])
                file01.write( f"{DYN_PROPERTIES['MD_STEP']}  " +  " ".join(map("{:2.8f}".format,OVERLAP )) + "\n" )

        with open("MD_OUTPUT/NACT.dat","a") as file01: 
            if ( DYN_PROPERTIES['MD_STEP'] == 0 ): 
                file01.write(f"# Step " + " ".join([f'{j}-{k}' for j in range(NStates) for k in range(j,NStates)]) + " Units = [meV]\n" )
            if ( DYN_PROPERTIES['MD_STEP'] >= 1 ): 
                NACT = DYN_PROPERTIES['NACT_NEW'] * 27.2114 * 1000 # 1 / a.u.t. --> meV
                NACT = np.array([NACT[j,k] for j in range(NStates) for k in range(j,NStates)])
                file01.write( f"{DYN_PROPERTIES['MD_STEP']}  " +  " ".join(map("{:2.8f}".format,NACT )) + "\n" )

        with open("MD_OUTPUT/Overlap_uncorrected.dat","a") as file01:
            if ( DYN_PROPERTIES['MD_STEP'] == 0 ): 
                file01.write(f"# Step " + " ".join([f'{j}-{k}' for j in range(NStates) for k in range(j,NStates)]) + "\n" )
            if ( DYN_PROPERTIES['MD_STEP'] >= 1 ): 
                OVERLAP = DYN_PROPERTIES['OVERLAP_NEW_uncorrected'] * 1.0
                OVERLAP = np.array([OVERLAP[j,k] for j in range(NStates) for k in range(j,NStates)])
                file01.write( f"{DYN_PROPERTIES['MD_STEP']}  " +  " ".join(map("{:2.8f}".format,OVERLAP )) + "\n" )

        with open("MD_OUTPUT/NACT_uncorrected.dat","a") as file01:
            if ( DYN_PROPERTIES['MD_STEP'] == 0 ): 
                file01.write(f"# Step " + " ".join([f'{j}-{k}' for j in range(NStates) for k in range(j,NStates)]) + "\n" )
            if ( DYN_PROPERTIES['MD_STEP'] >= 1 ): 
                NACT = DYN_PROPERTIES['NACT_NEW_uncorrected'] * 27.2114 * 1000 # 1 / a.u.t. --> meV
                NACT = np.array([NACT[j,k] for j in range(NStates) for k in range(j,NStates)])
                file01.write( f"{DYN_PROPERTIES['MD_STEP']}  " +  " ".join(map("{:2.8f}".format,NACT )) + "\n" )

        #with open("MD_OUTPUT/NACR.xyz","a") as file01:
        #    file01.write(f"{DYN_PROPERTIES['NAtoms']}\n")
        #    file01.write(f"MD Step {DYN_PROPERTIES['MD_STEP']} Units = [1 / Angstroms]\n")
        #    Atom_labels = DYN_PROPERTIES["Atom_labels"]
        #    Atom_NACR = DYN_PROPERTIES["NACR_APPROX_NEW"] / 0.529 # 0.529 Ang./Bohr
        #    for count, atom in enumerate( Atom_labels ):
        #        file01.write(f"{atom}  " + " ".join(map("{:2.8f}".format,Atom_NACR[count,:]))  + "\n") # Ang / fs
