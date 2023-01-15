%oldchk=../GS_NEW/geometry.chk
%chk=geometry.chk
%nprocshared=12
%mem=5GB

# WB97XD/6-311G SCF=XQC TD=(singlets,nstates=3,root=1) FORCE nosym pop=full guess=read

MD Step 80

0 1
C  -0.33455903  -0.58263897  0.00750152 
H  0.18155580  -1.53363858  -0.02692242 
H  -1.42197947  -0.61318297  -0.02576156 
C  0.34215039  0.59203317  0.00555010 
H  -0.20788267  1.50280225  -0.03064309 
H  1.40146170  0.57570384  -0.02920670 








