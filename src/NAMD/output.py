import numpy as np
import subprocess as sp

import properties

def save_data(DYN_PROPERTIES):

    NStates = DYN_PROPERTIES["NStates"]
    TIME    = round(DYN_PROPERTIES["MD_STEP"] * DYN_PROPERTIES["dtI"] / 41.341 ,4) # Print in fs


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
        file01.write(f"MD Step {TIME} Units = [eV / Ang.]\n")
        Atom_labels = DYN_PROPERTIES["Atom_labels"]
        Atom_forces = DYN_PROPERTIES["FORCE_NEW"] / 0.529 * 27.2114 # 0.529 Ang./Bohr, 27.2114 eV / Hartree
        for count, atom in enumerate( Atom_labels ):
            file01.write(f"{atom}  " + " ".join(map("{:2.8f}".format,Atom_forces[count,:]))  + "\n") # eV / Ang.


    with open("MD_OUTPUT/PES.dat","a") as file01:
        if ( NStates >= 2 ):
            file01.write( f"{TIME}  " +  " ".join(map("{:2.8f}".format,DYN_PROPERTIES["DIAG_ENERGIES_NEW"]*27.2114 )) + "\n" )
        else:
            file01.write( f"{TIME}  {np.round(DYN_PROPERTIES['DIAG_ENERGIES_NEW']*27.2114,8)}\n" )

    with open("MD_OUTPUT/Energy.dat","a") as file01:
        
        DYN_PROPERTIES = properties.compute_KE(DYN_PROPERTIES)
        DYN_PROPERTIES = properties.compute_PE(DYN_PROPERTIES)

        KE = DYN_PROPERTIES["KE"] * 27.2114
        PE = DYN_PROPERTIES["PE"] * 27.2114
        TE = KE + PE

        file01.write(f"{TIME}  " + "%2.6f  %2.6f  %2.6f\n" % (KE,PE,TE))

    with open("MD_OUTPUT/Temperature.dat","a") as file01:
        T = properties.compute_Temperature(DYN_PROPERTIES) # k
        file01.write(f"{TIME}  " + "%2.4f\n" % (T))


    if ( NStates >= 2 or ( DYN_PROPERTIES["BOMD"] == True and DYN_PROPERTIES["ISTATE"] != 0 ) ):

        with open("MD_OUTPUT/Population.dat","a") as file01:
            if ( DYN_PROPERTIES['NAMD_METHOD'] in ["GFSH"] ):
                AS = DYN_PROPERTIES['ACTIVE_STATE']
                POP     = np.zeros(( NStates ))
                POP[AS] = 1.0
                file01.write( f"{TIME}  " +  " ".join(map("{:2.2f}".format,POP )) + "  %2.2f" % (np.sum(POP)) + "\n" )
            else:
                POP = np.real(properties.get_density_matrix( DYN_PROPERTIES )[np.diag_indices(NStates)])
                file01.write( f"{TIME}  " +  " ".join(map("{:2.8f}".format,POP )) + "  %2.8f" % (np.sum(POP)) + "\n" )

        if ( DYN_PROPERTIES['NAMD_METHOD'] in ["GFSH"] ):
            with open("MD_OUTPUT/Active_State.dat","a") as file01:
                AS      = DYN_PROPERTIES['ACTIVE_STATE']
                file01.write( f"{TIME} {AS}\n" )


        if ( DYN_PROPERTIES["BOMD"] == False ):

            with open("MD_OUTPUT/mapping_re.dat","a") as file01:
                file01.write( f"{TIME}  " +  " ".join(map("{:2.8f}".format,DYN_PROPERTIES["MAPPING_VARS"].real )) + "\n" )

            with open("MD_OUTPUT/mapping_im.dat","a") as file01:
                file01.write( f"{TIME}  " +  " ".join(map("{:2.8f}".format,DYN_PROPERTIES["MAPPING_VARS"].imag )) + "\n" )

            with open("MD_OUTPUT/Coherence_re.dat","a") as file01:
                if ( DYN_PROPERTIES['MD_STEP'] == 0 ): 
                    file01.write(f"# Step " + " ".join([f'{j}-{k}' for j in range(NStates) for k in range(j+1,NStates)]) + "\n" )
                RHO = np.real(properties.get_density_matrix( DYN_PROPERTIES ))
                COH = np.array([RHO[j,k] for j in range(NStates) for k in range(j+1,NStates)])
                file01.write( f"{TIME}  " +  " ".join(map("{:2.8f}".format,COH )) + "\n" )

            with open("MD_OUTPUT/Coherence_im.dat","a") as file01:
                if ( DYN_PROPERTIES['MD_STEP'] == 0 ): 
                    file01.write(f"# Step " + " ".join([f'{j}-{k}' for j in range(NStates) for k in range(j+1,NStates)]) + "\n" )
                RHO = np.imag(properties.get_density_matrix( DYN_PROPERTIES ))
                COH = np.array([RHO[j,k] for j in range(NStates) for k in range(j+1,NStates)])
                file01.write( f"{TIME}  " +  " ".join(map("{:2.8f}".format,COH )) + "\n" )

            with open("MD_OUTPUT/Overlap.dat","a") as file01:
                if ( DYN_PROPERTIES['MD_STEP'] == 0 ): 
                    file01.write(f"# Step " + " ".join([f'{j}-{k}' for j in range(NStates) for k in range(j,NStates)]) + "\n" )
                if ( DYN_PROPERTIES['MD_STEP'] >= 1 ): 
                    OVERLAP = DYN_PROPERTIES['OVERLAP_NEW'] * 1.0
                    OVERLAP = np.array([OVERLAP[j,k] for j in range(NStates) for k in range(j,NStates)])
                    file01.write( f"{TIME}  " +  " ".join(map("{:2.8f}".format,OVERLAP )) + "\n" )

            with open("MD_OUTPUT/NACT.dat","a") as file01: 
                if ( DYN_PROPERTIES['MD_STEP'] == 0 ): 
                    file01.write(f"# Step " + " ".join([f'{j}-{k}' for j in range(NStates) for k in range(j,NStates)]) + " Units = [meV]\n" )
                if ( DYN_PROPERTIES['MD_STEP'] >= 1 ): 
                    NACT = DYN_PROPERTIES['NACT_NEW'] * 27.2114 * 1000 # 1 / a.u.t. --> meV
                    NACT = np.array([NACT[j,k] for j in range(NStates) for k in range(j,NStates)])
                    file01.write( f"{TIME}  " +  " ".join(map("{:2.8f}".format,NACT )) + "\n" )

            #with open("MD_OUTPUT/Overlap_uncorrected.dat","a") as file01:
            #    if ( DYN_PROPERTIES['MD_STEP'] == 0 ): 
            #        file01.write(f"# Step " + " ".join([f'{j}-{k}' for j in range(NStates) for k in range(j,NStates)]) + "\n" )
            #    if ( DYN_PROPERTIES['MD_STEP'] >= 1 ): 
            #        OVERLAP = DYN_PROPERTIES['OVERLAP_NEW_uncorrected'] * 1.0
            #        OVERLAP = np.array([OVERLAP[j,k] for j in range(NStates) for k in range(j,NStates)])
            #        file01.write( f"{TIME}  " +  " ".join(map("{:2.8f}".format,OVERLAP )) + "\n" )

            #with open("MD_OUTPUT/NACT_uncorrected.dat","a") as file01:
            #    if ( DYN_PROPERTIES['MD_STEP'] == 0 ): 
            #        file01.write(f"# Step " + " ".join([f'{j}-{k}' for j in range(NStates) for k in range(j,NStates)]) + "\n" )
            #    if ( DYN_PROPERTIES['MD_STEP'] >= 1 ): 
            #        NACT = DYN_PROPERTIES['NACT_NEW_uncorrected'] * 27.2114 * 1000 # 1 / a.u.t. --> meV
            #        NACT = np.array([NACT[j,k] for j in range(NStates) for k in range(j,NStates)])
            #        file01.write( f"{TIME}  " +  " ".join(map("{:2.8f}".format,NACT )) + "\n" )

            if ( DYN_PROPERTIES['PRINT_NACR'] ):
                with open("MD_OUTPUT/NACR.xyz","a") as file01:
                    file01.write(f"{DYN_PROPERTIES['NAtoms']}\n")
                    file01.write(f"MD Step {TIME} Units = [1 / Angstroms]\n")
                    Atom_labels = DYN_PROPERTIES["Atom_labels"]
                    Atom_NACR = DYN_PROPERTIES["NACR_APPROX_NEW"] / 0.529 # 0.529 Ang./Bohr
                    for count, atom in enumerate( Atom_labels ):
                        file01.write(f"{atom}  " + " ".join(map("{:2.8f}".format,Atom_NACR[count,:]))  + "\n") # Ang / fs

