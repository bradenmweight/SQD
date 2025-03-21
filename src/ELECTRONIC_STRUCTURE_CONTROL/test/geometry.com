%chk=geometry.chk
%mem=1GB
%nprocshared=1

#P B3LYP/6-31G
#P TD=(singlets,nstates=5)

Title

0 1
Li 0.0 0.0 0.0
F  0.0 0.0 3.5




