
%chk=geometry.chk
%rwf=geometry.rwf
%mem=1GB
%nprocshared=1

#P B3LYP/cc-pVTZ TD=(singlets,nstates=4)

TitleMe

0 1
H 0.0 0.0 0.0
F 0.0 0.0 1.5

                 


