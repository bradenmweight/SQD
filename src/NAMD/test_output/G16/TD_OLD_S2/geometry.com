%oldchk=../TD_NEW_S1/geometry.chk
%chk=geometry.chk
%nprocshared=12
%mem=5GB

# WB97XD/6-311G SCF=XQC TD=(read,singlets,nstates=3,root=2) FORCE nosym pop=full guess=read

MD Step 79

0 1
C  -0.33422811  -0.58212908  0.00715798 
H  0.18061101  -1.53634318  -0.02502026 
H  -1.42536774  -0.61649184  -0.02389277 
C  0.34236032  0.59217815  0.00525314 
H  -0.21027823  1.50341235  -0.02871510 
H  1.40174467  0.57329802  -0.02722661 








