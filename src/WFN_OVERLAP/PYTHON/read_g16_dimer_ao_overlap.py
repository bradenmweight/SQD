import numpy as np
import sys
from matplotlib import pyplot as plt

# G16 INPUT: "# B3LYP/sto-3g nosymm iop(2/12=3,3/33=1) guess=only pop=full"

def get_AO_Matrix(filename):
    OUT_LINES = open(filename,"r").readlines()

    N_AO_BASIS_TOTAL = 0
    N_MO_BASIS_TOTAL = 0
    AO_overlap = []
    for count, line in enumerate(OUT_LINES):
        t = line.split()
        if ( N_AO_BASIS_TOTAL == 0 and len(t) >= 2 and t[0] == "NBasis"):
            try:
                N_AO_BASIS_TOTAL = int(t[2])
            except ValueError: # Not sure why this happens, but "= 1234" --> "=1234" sometimes.... ~BMW
                N_AO_BASIS_TOTAL = int(t[1].split("=")[1])
            N_MO_BASIS_TOTAL = N_AO_BASIS_TOTAL
            AO_overlap = np.zeros(( N_AO_BASIS_TOTAL, N_AO_BASIS_TOTAL ), dtype=float)
        
        if ( len(t) == 3 and t == "*** Overlap *** ".split() ):
            counter = count
            while ( True ):
                counter += 1

                if ( OUT_LINES[counter].split() == " *** Kinetic Energy ***".split() ):
                    break
            
                # READ AO ORBITAL OVERLAP
                try: # Will work if AO labels found on this line
                    AO1_Labels = [ int(j) for j in OUT_LINES[counter].split() ]
                    N_COLUMNS = len( AO1_Labels )
                except ValueError: # Will return ValueError if AO data is found
                    s = OUT_LINES[counter].split()
                    ao2 = int( s[0] )
                    for col in range( len(s)-1 ):
                        ao1 = AO1_Labels[col]
                        AO_overlap[ao1-1,ao2-1] = float( s[col+1].replace("D", "e") )
                        AO_overlap[ao2-1,ao1-1] = AO_overlap[ao1-1,ao2-1]

    OFF_BLOCK = AO_overlap[:N_AO_BASIS_TOTAL//2,N_AO_BASIS_TOTAL//2:]

    np.savetxt("AO_OVERLAP_TOTAL.dat", AO_overlap)
    np.savetxt("AO_OVERLAP_OFF_BLOCK.dat", OFF_BLOCK)


    plt.imshow( np.abs(AO_overlap), origin='lower', cmap="hot_r" )
    plt.colorbar(pad=0.01)
    plt.xlabel("AO Basis Index",fontsize=15)
    plt.ylabel("AO Basis Index",fontsize=15)
    plt.tight_layout()
    plt.savefig("AO_OVERLAP_DIMER.jpg",dpi=600)
    plt.clf()

    plt.imshow( np.abs(OFF_BLOCK), origin='lower', cmap="hot_r" )
    plt.colorbar(pad=0.01)
    plt.xlabel("AO Basis Index",fontsize=15)
    plt.ylabel("AO Basis Index",fontsize=15)
    plt.tight_layout()
    plt.savefig("AO_OVERLAP_OFF_BLOCK.jpg",dpi=600)
    plt.clf()

    # Print Fortran-friendly AO Overlap matrix as j k (AO)_jk (3-Column Format)
    f = open("AO_OVERLAP_TOTAL_ijv.dat","w")
    f.write("#  Overlap between R and R+dR : i_AO_R, j_AO_R+dR, S_ij\n")
    for j in range( N_AO_BASIS_TOTAL ):
        for k in range( N_AO_BASIS_TOTAL ):
            f.write( f"{j+1} {k+1} {AO_overlap[j,k]}\n" )

    f = open("AO_OVERLAP_OFF_BLOCK_ijv.dat","w")
    f.write("#  Overlap between R and R+dR : i_AO_R, j_AO_R+dR, S_ij\n")
    for j in range( N_AO_BASIS_TOTAL//2 ):
        for k in range( N_AO_BASIS_TOTAL//2, N_AO_BASIS_TOTAL ):
            f.write( f"{j+1} {k+1} {AO_overlap[j,k]}\n" )

    return AO_overlap

if ( __name__ == "__main__" ):
    filename = sys.argv[1]
    get_AO_Matrix(filename)