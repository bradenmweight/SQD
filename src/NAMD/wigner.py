import numpy as np
from random import random
from scipy.special import laguerre
import subprocess as sp
import os

def get_Globals():
    global TEMP, NSamples, AU_to_K, AU_to_CMinv
    TEMP = 300 # K
    NSamples = 10 # Number of initial conditions

    AU_to_K     = 315777. # K / Hartree
    #K_to_eV     = 0.025/300 # eV / K
    AU_to_CMinv = 219474.6 # cm-1 / Hartree

def read_eigenvectors():

    """
    Run gaussian 16 job similar to this:
        # OPT FREQ=saveNM B3LYP/6-311++G**
    And generate formatted checkpoint file like this:
        formchk geometry.chk
    Need three files:
        1. geometry.com
        2. geometry.out
        3. geometry.fchk
    """

    COM_FILE  = open("geometry.com","r").readlines()
    OUT_FILE  = open("geometry.out","r").readlines()
    FCHK_FILE = open("geometry.fchk","r").readlines()

    NAtoms = int( FCHK_FILE[2].split()[-1] )
    write_text_file = True # For large systems, turn off.

    atom_labels = []
    geom      = np.zeros(( NAtoms, 3 ))
    for count, line in enumerate( COM_FILE ):
        if ( "0 1".split() == line.split() ):
            for at in range( NAtoms ):
                t = COM_FILE[ count+1+at ].split()
                atom_labels.append( t[0] )
                geom[at,:] = np.array( t[1:] ).astype(float)

    NM_EIGV = []
    for count, line in enumerate( FCHK_FILE ):
        t = line.split()
        if ( "Vib-Modes" in t ):
            length = int( t[-1] )
            NColumns = len(FCHK_FILE[ count+1 ].split())
            for n in range( int(np.ceil(length/NColumns)) ):
                s = FCHK_FILE[ count+1+n ].split()
                for j in range( len(s) ):
                    NM_EIGV.append( float(s[j]) )

    NM_EIGV = np.array( NM_EIGV )
    NModes = len(NM_EIGV) // NAtoms // 3
    print("Number of Atoms:", NAtoms)
    print("Number of Normal Modes (3*N - 6):", NModes, f"({3*NAtoms - 6})" )
    NM_EIGV = NM_EIGV.reshape(( NModes, NAtoms, 3 ))
    np.save( "NM_EIGV.dat.npy", NM_EIGV )
    if( write_text_file == True ):
        f = open("NM_EIGV.dat","w")
        #f.write("" + "\n")
        for mode in range(NModes):
            f.write(f"Mode #{mode+1}" + "\n")
            for at in range(NAtoms):
                f.write( " ".join(map("{:1.5f}".format, NM_EIGV[mode,at,:])) + "\n" )
        f.close()


    # Read Data
    NM_FREQ = []
    for line in OUT_FILE:
        t = line.split()
        if ( len(t) == 5 ):
            if ( t[0] == "Frequencies" and t[1] == "--" ):
                NM_FREQ.append( float(t[2]) ) # cm^-1
                NM_FREQ.append( float(t[3]) )
                NM_FREQ.append( float(t[4]) )

    NM_FREQ = np.array(NM_FREQ)
    np.savetxt("NM_FREQ.dat", np.array(NM_FREQ) )

    return atom_labels, geom, NM_FREQ, NM_EIGV

def set_masses(atom_labels):
    mass_amu_to_au = 1837/1.007 # au / amu
    masses_amu = { "H":1.007, 
                   "He":4.00260,
                   "Li":6.941,
                   "Be":9.01218,
                   "B":10.81,
                   "C":12.011,
                   "N":14.0067,
                   "O":15.9994,
                   "F":18.998403,
                   "Ne":20.179,
                   "Cl":35.453,
                   "Fe":55.845,
                   "Cu":63.546}

    masses = []
    for at in atom_labels:
        masses.append( masses_amu[at] )
    return np.array(masses) * mass_amu_to_au

def determine_state(mode_freq):
    print("I DONT KNOW WHAT IS HAPPENING HERE!!!!!!!!")
    print("Check the units for all things !!!!!!!!")
    thresh  = 0.9999 # Keep all n below population threshold

    #print( mode_freq/AU_to_CMinv*AU_to_K / TEMP )
    #exit()
    exponent   = mode_freq/AU_to_CMinv*AU_to_K / TEMP # TEMP is global
  
    if ( exponent > 800 ):
        exponent = 600
        print (f"The partition function is too close to zero to use, since its exponent was ~{int(exponent)}" )
        print (f"Either because of too low temperature (T = {TEMP} K) or very high frequency {mode_freq} cm^-1.")
        print ("It was has been reset to %e" % (np.exp(-exponent/2.) / ( 1. - np.exp(-exponent) ) ) )
  
    partition_function = np.exp(-exponent/2.) / ( 1. - np.exp(-exponent) )
    n     = -1
    sum_p = 0.
    prob  = []
    # calculate probabilities until sum is larger than threshold
    while ( True ):
        n += 1
        p = np.exp(-exponent*(n+1./2.))/partition_function
        prob.append(p)
        sum_p += prob[n]
        if (sum_p >= thresh):
            break
    n = -1
    probability = 0.
    # generate random number that is smaller than threshold
    while ( True ):
        random_state = random()
        if random_state < sum_p:
            break
    # determine state number by comparing with random number
    while ( True ):
        n += 1
        probability += prob[n]
        if ( probability >= random_state ):
            return n
            break
def wigner_prob(Q,P,mode_freq):
    """
    Q: Unitless Position
    P: Unitless Momentum
    n: Vibrational Mode Number (QHO)
    """

    n = determine_state(mode_freq)

    if ( n >= 500 ):
        print(f"Vibrational mode is n = {n}. This too large to handle with harmonic approximation.")
        print("We will simply discard this mode from consideration.")
        n = -1

    R2 = 2.0 * (P**2 + Q**2)
    prob = (-1.0)**n * laguerre(n)(R2) * np.exp(-R2/2.0)
    return (prob,n)

def shift_COM(GEOM,masses):

    NAtoms = len(GEOM)
    COM = np.zeros((3))
    for at in range(NAtoms):
        for d in range(3):
            COM[d] += masses[at] * GEOM[at,d]

    M_Total = np.sum( masses )

    COM = COM / M_Total
    GEOM -= COM

    return GEOM

def momentOfInertia(W_GEOM,W_VELOC,masses):
    NAtoms = len(W_GEOM)
    MOI     = np.zeros((3,3))
    one_mat = np.identity(3)

    for at in range(NAtoms):
        MOI += masses[at] * \
            ( np.dot(W_GEOM[at,:],W_GEOM[at,:])*one_mat - \
              np.outer(W_GEOM[at],W_GEOM[at]))
    return MOI

def myCrossProduct(a,b):
    return np.array([a[1]*b[2]-a[2]*b[1],\
                        a[2]*b[0]-a[0]*b[2],\
                        a[0]*b[1]-a[1]*b[0]])

def angularMomentum(W_GEOM,W_VELOC,masses):
    NAtoms = len(W_GEOM)
    ANG_MOM = np.array([0.0,0.0,0.0])
    for at in range(NAtoms):
        ANG_MOM += masses[at] * myCrossProduct(W_GEOM[at],W_VELOC[at])
    return ANG_MOM

def getAngularVelocity(W_GEOM,W_VELOC,masses):
    ANG_MOM = angularMomentum(W_GEOM,W_VELOC,masses)
    MOI     = momentOfInertia(W_GEOM,W_VELOC,masses)
    ANG_VELOC = np.dot(np.linalg.inv(MOI),ANG_MOM)
    return ANG_VELOC

def remove_rotations(W_GEOM,W_VELOC,masses):
    NAtoms = len(W_GEOM)
    ANG_VELOC = getAngularVelocity(W_GEOM,W_VELOC,masses)
    for at in range(NAtoms):
        W_VELOC[at,:] -= myCrossProduct(ANG_VELOC[:],W_GEOM[at,:])
    return W_VELOC

def get_one_initial_condition(atom_labels, GEOM, masses, NM_FREQ, NM_EIGV):
    NAtoms = len(atom_labels)
    NModes = len(NM_EIGV)

    W_GEOM  = np.zeros(( NAtoms, 3 ))
    W_VELOC = np.zeros(( NAtoms, 3 ))

    E_sample = 0.0
    for mode in range(NModes):
        while (True):
            # Sample unitless QHO coordinates from -10 to 10
            random_Q = 10 * (2*random() - 1)
            random_P = 10 * (2*random() - 1)
            prob,n = wigner_prob( random_Q, random_P, NM_FREQ[mode] )
            if ( prob > 1.0 or prob < 0 ):
                print(f"Wigner probability is P = {prob}. Something is wrong.")
                print(f"Recalculating the probability.")
                continue
                #exit()
            elif( prob > random() ):
                #print(prob,n)
                break
        
        # Convert from unitless Q,P to atomic units (AU) X:Bohr and V:Bohr/a.u.
        random_Q *= 0.529 / np.sqrt( NM_FREQ[mode] / AU_to_CMinv ) # X ~ sqrt[h/w]
        random_P *= np.sqrt( NM_FREQ[mode] / AU_to_CMinv ) # V ~ sqrt[ hw ] 

        # Compute QHO energy for this mode and add to total for this sample
        E_sample += 0.50000000 * (NM_FREQ[mode] / AU_to_CMinv)**2 * random_Q**2

        # Get new coordinates and velocities
        for count, atom in enumerate(atom_labels):
            for d in range(3):
                W_GEOM[count,d]  += GEOM[count,d] + 0.529 * random_Q * NM_EIGV[mode,count,d] * 1/np.sqrt(masses[count]) 
                W_VELOC[count,d] += 0.529 * 41.341 * random_P * NM_EIGV[mode,count,d] * 1/np.sqrt(masses[count]) # P ~ sqrt[ hwm / 2 ] -> v ~ sqrt[ hw / (2m) ] 

        W_GEOM  = shift_COM(W_GEOM,masses)
        W_VELOC = remove_rotations(W_GEOM,W_VELOC,masses)

        return W_GEOM, W_VELOC, E_sample

def write_XYZs(atom_labels,W_GEOM,W_VELOC,E_sample,sample):
    NAtoms = len(atom_labels)
    OUTPUT_DIR = "WIGNER_INITIAL_CONDITIONS"
    if ( sample == 0 and os.path.exists(OUTPUT_DIR) ):
        sp.call(f"rm -r {OUTPUT_DIR}", shell=True)
        sp.call(f"mkdir {OUTPUT_DIR}", shell=True)
    COORD_FILE = open(f"{OUTPUT_DIR}/INIT_GEOM_{sample}.xyz","w")
    COORD_FILE.write(f"{NAtoms}\n")
    COORD_FILE.write("E = %2.6f\n" % (E_sample) )
    VELOC_FILE = open(f"{OUTPUT_DIR}/INIT_VELOC_{sample}.xyz","w")
    VELOC_FILE.write(f"{NAtoms}\n")
    VELOC_FILE.write("E = %2.6f\n" % (E_sample) )
    for at,atom in enumerate( atom_labels ):
        COORD_FILE.write(f"{atom}  " + "  ".join(map("{:1.8f}".format, W_GEOM[at,:])) + "\n")
        VELOC_FILE.write(f"{atom}  " + "  ".join(map("{:1.8f}".format,W_VELOC[at,:])) + "\n")
    COORD_FILE.close()
    VELOC_FILE.close()

def get_initial_conditions(atom_labels, GEOM, masses, NM_FREQ, NM_EIGV):

    for sample in range( NSamples ):
        print(f"Sample {sample} of {NSamples}")
        W_GEOM, W_VELOC, E_sample = get_one_initial_condition(atom_labels, GEOM, masses, NM_FREQ, NM_EIGV)
        write_XYZs(atom_labels,W_GEOM,W_VELOC,E_sample,sample)

def check_orthogonality(atom_labels, GEOM, masses, NM_FREQ, NM_EIGV):
    print("Checking normal mode orthogonality.")
    NAtoms = len(atom_labels)
    NModes = len(NM_FREQ)
    AU_to_AMU = 1.00784/1837 # AU -> AMU
    #masses *= AU_to_AMU

    # Try to make these modes orthogonal
    
    # Weight each mode by 1/sqrt(m)
    #M = np.einsum("a,b->ab", 1/np.sqrt(masses), 1/np.sqrt(masses))
    #Ojk = np.einsum("xad,ab,ybd->xy", NM_EIGV, M, NM_EIGV)
    #print(Ojk[np.diag_indices(NModes)])
    #print(np.round(Ojk[:8,:8],4))
#
    ## Do nothing
    #Ojk = np.einsum("xad,yad->xy", NM_EIGV, NM_EIGV) # Default
    #print(Ojk[np.diag_indices(NModes)])
    #print(np.round(Ojk[:8,:8],4))
 #
    ## Try weighting each mode by 1/m
    #M = np.einsum("a,b->ab", 1/masses, masses)
    #Ojk = np.einsum("xad,ab,ybd->xy", NM_EIGV, M, NM_EIGV)
    #print(Ojk[np.diag_indices(NModes)])
    #print(np.round(Ojk[:8,:8],4))
    #        
    #        if ( mj == mk and Ojk < 0.98 ):
    #            print("%d %d %2.10f" % (mj,mk,Ojk) )
    #            print("%d %d %2.10f" % (mj,mk,Ojk) )
    #        if ( mj != mk and Ojk > 0.05 ):
    #            print("%d %d %2.10f" % (mj,mk,Ojk) )
    #            print("%d %d %2.10f" % (mj,mk,Ojk) )

    #exit()

def main(): 
    get_Globals()
    atom_labels, GEOM, NM_FREQ, NM_EIGV = read_eigenvectors()
    masses = set_masses(atom_labels)
    check_orthogonality(atom_labels, GEOM, masses, NM_FREQ, NM_EIGV)
    get_initial_conditions(atom_labels, GEOM, masses, NM_FREQ, NM_EIGV)

if ( __name__ == "__main__" ):
    main()