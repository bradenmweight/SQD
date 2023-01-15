import numpy as np
import subprocess as sp
#from build.cioverlap import wf_overlap
from matplotlib import pyplot as plt
import os, re

def CI_overlap(DYN_PROPERTIES):

    # Save previous step
    #if (DYN_PROPERTIES["MD_STEP"] >= 1 ):
    #    DYN_PROPERTIES["AO_OVERLAP_OLD"] = DYN_PROPERTIES["AO_OVERLAP_NEW"] * 1.0
    #    DYN_PROPERTIES["MO_COEFF_OLD"]   = DYN_PROPERTIES["MO_COEFF_NEW"]   * 1.0
    #    DYN_PROPERTIES["CI_COEFF_OLD"]   = DYN_PROPERTIES["CI_COEFF_NEW"]   * 1.0

    # Read AO overlap
    os.chdir("DIMER/")
    DYN_PROPERTIES["AO_OVERLAP_NEW"] = read_ao_overlap('./geometry.rwf')
    os.chdir("../")
    
    # Read MO coefficients (WILL STORE FOR PRODUCION. REMOVE LATER)
    os.chdir("GS_OLD")
    DYN_PROPERTIES["MO_COEFF_OLD"] = read_mo_coef("./geometry.rwf") 
    os.chdir("../")

    # Read MO coefficients NEW
    os.chdir("GS_NEW")
    DYN_PROPERTIES["MO_COEFF_NEW"] = read_mo_coef("./geometry.rwf") 
    os.chdir("../")
    
    # Read CI coefficients OLD (WILL STORE FOR PRODUCION. REMOVE LATER)
    os.chdir("TD_OLD_S1")
    tmp, orb_dict = read_xy_coef("./geometry.rwf") # GS not included
    #DYN_PROPERTIES["CI_COEFF_NEW"] = np.append( np.zeros( tmp.shape[1:] ), tmp ) # GS is now included
    DYN_PROPERTIES["CI_COEFF_OLD"] = tmp # GS not included
    os.chdir("../")

    # Read CI coefficients NEW
    os.chdir("TD_NEW_S1")
    tmp, orb_dict = read_xy_coef("./geometry.rwf") # GS not included
    #DYN_PROPERTIES["CI_COEFF_NEW"] = np.append( np.zeros( tmp.shape[1:] ), tmp ) # GS is now included
    DYN_PROPERTIES["CI_COEFF_NEW"] = tmp # GS not included
    os.chdir("../")

    # Remove MO coefficients not involved in the TD-DFT correlation
    # Reshape the MO arrays to be (Nocc-NFC,NAOs,Nvir-NFV,NAOs)
    NAOs = orb_dict["NORBS"]
    KEEP = [orb_dict["N_CORE_FROZEN"],orb_dict["NORBS"]-orb_dict["N_VIRT_FROZEN"]]
    DYN_PROPERTIES["MO_COEFF_OLD"] = (DYN_PROPERTIES["MO_COEFF_OLD"])[KEEP[0]:KEEP[1]]
    DYN_PROPERTIES["MO_COEFF_NEW"] = (DYN_PROPERTIES["MO_COEFF_NEW"])[KEEP[0]:KEEP[1]]

    # Print all sizes
    print("AO:", DYN_PROPERTIES["AO_OVERLAP_NEW"].shape)
    print("MO (OLD):", DYN_PROPERTIES["MO_COEFF_OLD"].shape)
    print("MO (NEW):", DYN_PROPERTIES["MO_COEFF_NEW"].shape)
    print("CI (OLD):", DYN_PROPERTIES["CI_COEFF_OLD"].shape)
    print("CI (NEW):", DYN_PROPERTIES["CI_COEFF_NEW"].shape)

    # Calculate wavefunction overlap
    #wf_overlap(self, molecule, istep, dt)

    No = orb_dict["NOCC"] - orb_dict["N_CORE_FROZEN"]
    Nv = orb_dict["NVIR"] - orb_dict["N_VIRT_FROZEN"]

    normOLD = np.einsum("Jab,Jab->J", DYN_PROPERTIES["CI_COEFF_OLD"], DYN_PROPERTIES["CI_COEFF_OLD"])
    normNEW = np.einsum("Jab,Jab->J", DYN_PROPERTIES["CI_COEFF_NEW"], DYN_PROPERTIES["CI_COEFF_NEW"])
    #print( "CI Norms:", normOLD, normNEW )

    # Renormalize CI coefficients
    DYN_PROPERTIES["CI_COEFF_OLD"] = np.einsum( "Jab,J->Jab", DYN_PROPERTIES["CI_COEFF_OLD"], 1/np.sqrt(normOLD) )
    DYN_PROPERTIES["CI_COEFF_NEW"] = np.einsum( "Jab,J->Jab", DYN_PROPERTIES["CI_COEFF_NEW"], 1/np.sqrt(normNEW) )

    normOLD = np.einsum("Jab,Jab->J", DYN_PROPERTIES["CI_COEFF_OLD"], DYN_PROPERTIES["CI_COEFF_OLD"])
    normNEW = np.einsum("Jab,Jab->J", DYN_PROPERTIES["CI_COEFF_NEW"], DYN_PROPERTIES["CI_COEFF_NEW"])
    print( "CI Norms:", normOLD, normNEW )

    # Compute the overlap between MOs
    S_MO = np.einsum( "an,bm,nm->ab",   \
                    DYN_PROPERTIES["MO_COEFF_OLD"], \
                    DYN_PROPERTIES["MO_COEFF_NEW"], \
                    DYN_PROPERTIES["AO_OVERLAP_NEW"]
     )
    #print( "MO OVERLAP:", S_MO )
    
    S_JK = np.einsum( "Jaj,ab,jk,Kbk->JK",   \
                    DYN_PROPERTIES["CI_COEFF_OLD"][:,1:,:], \
                    S_MO[:No,:No], \
                    S_MO[No:,No:], \
                    DYN_PROPERTIES["CI_COEFF_NEW"][:,1:,:]
     )

    S_JK[0,0] = np.linalg.det(S_MO[:No,:No])

    print(S_JK)


    return DYN_PROPERTIES


def get_AO_Matrix_BRADEN(filename):
    FCHK_LINES = open(filename,"r").readlines()

    N_AO_BASIS_TOTAL = 0
    N_MO_BASIS_TOTAL = 0
    AO_overlap = []
    for count, line in enumerate(FCHK_LINES):
        t = line.split()
        if ( N_AO_BASIS_TOTAL == 0 and len(t) >= 2 and t[0] == "NBasis"):
            N_AO_BASIS_TOTAL = int(t[2])
            N_MO_BASIS_TOTAL = N_AO_BASIS_TOTAL
            AO_overlap = np.zeros(( N_AO_BASIS_TOTAL, N_AO_BASIS_TOTAL ), dtype=float)
        
        if ( len(t) == 3 and t == "*** Overlap *** ".split() ):
            counter = count
            while ( True ):
                counter += 1

                if ( FCHK_LINES[counter].split() == " *** Kinetic Energy ***".split() ):
                    break
            
                # READ AO ORBITAL OVERLAP
                try: # Will work if AO labels found on this line
                    AO1_Labels = [ int(j) for j in FCHK_LINES[counter].split() ]
                    N_COLUMNS = len( AO1_Labels )
                except ValueError: # Will return ValueError if AO data is found
                    s = FCHK_LINES[counter].split()
                    ao2 = int( s[0] )
                    for col in range( len(s)-1 ):
                        ao1 = AO1_Labels[col]
                        AO_overlap[ao1-1,ao2-1] = float( s[col+1].replace("D", "e") )
                        AO_overlap[ao2-1,ao1-1] = AO_overlap[ao1-1,ao2-1]

    OFF_BLOCK = AO_overlap[:N_AO_BASIS_TOTAL//2,N_AO_BASIS_TOTAL//2:]

    #np.savetxt("AO_OVERLAP_TOTAL.dat", AO_overlap)
    #np.savetxt("AO_OVERLAP_OFF_BLOCK.dat", OFF_BLOCK)


    #plt.imshow( np.abs(AO_overlap), origin='lower', cmap="hot_r" )
    #plt.colorbar(pad=0.01)
    #plt.xlabel("AO Basis Index",fontsize=15)
    #plt.ylabel("AO Basis Index",fontsize=15)
    #plt.tight_layout()
    #plt.savefig("AO_OVERLAP_DIMER.jpg",dpi=600)
    #plt.clf()

    #plt.imshow( np.abs(OFF_BLOCK), origin='lower', cmap="hot_r" )
    #plt.colorbar(pad=0.01)
    #plt.xlabel("AO Basis Index",fontsize=15)
    #plt.ylabel("AO Basis Index",fontsize=15)
    #plt.tight_layout()
    #plt.savefig("AO_OVERLAP_OFF_BLOCK.jpg",dpi=600)
    #plt.clf()

    # Print Fortran-friendly AO Overlap matrix as j k (AO)_jk (3-Column Format)
    #f = open("AO_OVERLAP_TOTAL_ijv.dat","w")
    #f.write("#  Overlap between R and R+dR : i_AO_R, j_AO_R+dR, S_ij\n")
    #for j in range( N_AO_BASIS_TOTAL ):
    #    for k in range( N_AO_BASIS_TOTAL ):
    #        f.write( f"{j+1} {k+1} {AO_overlap[j,k]}\n" )

    #f = open("AO_OVERLAP_OFF_BLOCK_ijv.dat","w")
    #f.write("#  Overlap between R and R+dR : i_AO_R, j_AO_R+dR, S_ij\n")
    #for j in range( N_AO_BASIS_TOTAL//2 ):
    #    for k in range( N_AO_BASIS_TOTAL//2, N_AO_BASIS_TOTAL ):
    #        f.write( f"{j+1} {k+1} {AO_overlap[j,k]}\n" )

    return OFF_BLOCK # Return t,t+dt overlap



def read_ao_overlap(RWF_FILE):
        """ 
        Read a rwf file to obtain ao_overlap data
        """
        os.system(f"rwfdump {RWF_FILE} ao_overlap.dat 514R")

        with open('ao_overlap.dat', "r") as f:
            log = f.read()

        tmp = re.findall('[-]?\d+\.\d+D[+-]\d\d', log)
        tmp = [float(x.replace('D', 'e')) for x in tmp]
 
        NBasis = int( sp.check_output("grep 'NBasis' geometry.out | head -n 1 | awk '{print $3}'",shell=True) )//2
        #print("NBasis", NBasis)
        tmp_ovr = np.zeros((NBasis * 2, NBasis * 2))
 
        cnt = 0
        for ibasis in range(NBasis * 2):
            for jbasis in range(ibasis + 1):
                tmp_ovr[ibasis, jbasis] = tmp[cnt]
                cnt += 1

        tmp_ovr += np.transpose(tmp_ovr) - np.diag(np.diag(tmp_ovr))

        # Slicing the components between t and t+dt
        return tmp_ovr[:NBasis, NBasis:]

def read_mo_coef(RWF_FILE):
        """ 
        Read a rwf file to obtain mo_coef data
        """
        os.system(f"rwfdump {RWF_FILE} mo_coef.dat 524R")

        with open('mo_coef.dat', "r") as f:
            log = f.read()

        tmp = re.findall('[-]?\d+\.\d+D[+-]\d\d', log)
        tmp = np.array([x.replace('D','e') for x in tmp], dtype=np.float64)

        NBasis = int(sp.check_output("grep 'NBasis' geometry.out | tail -n 1 | awk '{print $2}'", shell=True) )

        tmp_mo = tmp.reshape(NBasis, NBasis)

        #return tmp_mo[NFC:NBasis] # Why we are only returing after NFC ?
        return tmp_mo

def read_xy_coef(RWF_FILE):
        """ 
        Read a rwf file to obtain xy_coef data
        """
        os.system(f"rwfdump {RWF_FILE} xy_coef.dat 635R")

        with open(f'xy_coef.dat', "r") as f:
            log = f.read()

        NBasis, NAE, NFC, NFV = np.array(sp.check_output("grep 'NFC' geometry.out | tail -n 1 | awk '{print $2,$4,$8,$10}'", shell=True).split(), dtype=int)
        NOA,NOB,NVA,NVB = np.array(sp.check_output("grep 'NVA=' geometry.out | tail -n 1 | awk '{print $4,$6,$8,$10}'", shell=True).split(), dtype=int)
        NORB = NOA + NVA # Alpha occupied + unoccupied
        NStates = int(sp.check_output("grep 'initial guesses have been made' geometry.out | tail -n 1 | awk '{print $1}'", shell=True))
        NStates //= 4 # Gaussian saves 4x the number of requested roots

        NOCC,NVIR = NAE, NBasis-NAE
        #print(f"N Orbs (OCC,VIR,TOT) = ({NOCC},{NVIR},{NOCC+NVIR})")
        #print(f"N Frozen (OCC,VIR)   = ({NFC},{NFV})")
        #print(f"N Active (OCC,VIR)   = ({NOA},{NVA})")

        orb_dict = { "NORBS":NBasis, "NOCC":NOCC, "NVIR":NVIR, "N_CORE_FROZEN":NFC, "N_VIRT_FROZEN":NFV }

        tmp = re.findall('[-]?\d+\.\S+[+-]\d+', log)

        # Drop the first 12 dummy elements. Why are these here ??? ~ BMW
        tmp = tmp[12:]
 
        # Gaussian09 deals with 4 times as much roots as the input NStates value.
        # the nr. of excitation function => nocc \times nvirt
        # spin degrees of freedom => 2
        # X+Y, X-Y => 2
        roots = NStates * 4
        num_coef = 4 * (NOA * NVA) * roots
 
        tmp = tmp[:num_coef]
        tmp = np.array([x.replace('D','e') for x in tmp], dtype=np.float64)
        XplusY, XminusY = tmp.reshape(2, roots, 2, NOA * NVA)# Something like X+-Y, roots, spin, mo's ~ BMW
        x = 0.5 * (XplusY + XminusY)
 
        # ADD CHECK FOR ALPHA=BETA ~ BMW
        assert( np.allclose( x[:NStates, 0, :], x[:NStates, 1, :] ) ), "Alpha and beta X+Y and X-Y were not the same. Needs to be spin-restricted calculation."

        # Drop beta part and unrequested excited states
        x = x[:NStates, 0, :]
        x = x.reshape(NStates, NOA, NVA)

        # We need to add the ground state to the tensor
        X_GS = np.zeros(( NOA + 1, NVA )) # Define ground state as an extra configuration at (0,0). We add an extra occupied orbital
        X_GS[0,0] = 1

        x_TOT = np.zeros(( NStates+1, NOA+1,  NVA )) # NRoots + G.S.
        x_TOT[0,:,:] = X_GS
        x_TOT[1:,1:,:] = x

        #print(x_TOT[0])
        #print(np.round(x_TOT[1],5))

        return x_TOT, orb_dict

if (__name__ == "__main__"):
    #print( get_AO_Matrix_BRADEN("geometry.out") )
    #print( read_ao_overlap("geometry.rwf") ) # This gives same result as mine. This is good.
    #print( read_mo_coef("geometry.rwf") )
    #x = read_xy_coef("geometry.rwf") # THIS RETURNS [(X+Y) + (X-Y)]/2  in the shape of (NStates, NOCC, NVIRT)

    DYN_PROPERTIES = {"MD_STEP":1}
    CI_overlap(DYN_PROPERTIES)