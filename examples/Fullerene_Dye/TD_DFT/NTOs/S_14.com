%oldchk="../geometry.chk"
%chk=S_14.chk
%mem=5GB
%nprocshared=1

#p AM1/chkbasis
#p Geom=AllCheck Guess=(Read,Only) Density=(Check,Transition=14) Pop=(Minimal,NTO,SaveNTO)

TitleMe

0 1




