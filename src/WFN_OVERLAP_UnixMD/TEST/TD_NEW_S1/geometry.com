%oldchk=../GS_NEW/geometry.chk
%chk=geometry.chk
%rwf=geometry.rwf
%nprocshared=1
%mem=5GB

# WB97XD/STO-3G SCF=XQC TD=(singlets,nstates=2,root=1) FORCE nosym pop=full guess=read

MD Step 1000

0 1
C  2.74013463  1.87525770  2.08311481 
H  -1.06507604  -1.03155638  2.48749307 
H  -1.09272681  -0.68990011  1.54417368 
C  1.53904851  3.05056607  2.05905601 
H  1.74189868  3.26331677  3.16784543 
H  1.55951981  3.73872426  1.25195983 








