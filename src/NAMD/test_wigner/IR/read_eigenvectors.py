import numpy as np
import os
from matplotlib import pyplot as plt
import matplotlib as mpl
from scipy.interpolate import interp1d

COM_FILE = open("geometry.com","r").readlines()
FCHK_FILE = open("geometry.fchk","r").readlines()

NAtoms = 6
write_text_file = True # For large systems, turn off.

at_labels = []
geom      = np.zeros(( NAtoms, 3 ))
for count, line in enumerate( COM_FILE ):
    if ( "0 1".split() == line.split() ):
        for at in range( NAtoms ):
            t = COM_FILE[ count+1+at ].split()
            at_labels.append( t[0] )
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
print(len(NM_EIGV))
NModes = len(NM_EIGV) // NAtoms // 3
print(NModes)
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
f = open("geometry.out","r") # Open File
FREQ = []
for line in f:
    t = line.split()
    if ( len(t) == 5 ):
        if ( t[0] == "Frequencies" and t[1] == "--" ):
            FREQ.append( float(t[2]) ) # cm^-1
            FREQ.append( float(t[3]) )
            FREQ.append( float(t[4]) )

np.savetxt("NM_FREQ.dat", np.array(FREQ) )

