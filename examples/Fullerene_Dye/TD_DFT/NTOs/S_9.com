%oldchk="../geometry.chk"
%chk=S_9.chk
%mem=5GB
%nprocshared=1

#p AM1/chkbasis
#p Geom=AllCheck Guess=(Read,Only) Density=(Check,Transition=9) Pop=(Minimal,NTO,SaveNTO)

TitleMe

0 1



