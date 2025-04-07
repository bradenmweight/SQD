# p   "MD_OUTPUT/NACT.dat" u 1:($8/100)  title 'd12/100' @wl
# rep "MD_OUTPUT/NACT.dat" u 1:($12/100) title 'd23/100' @wl
# rep "MD_OUTPUT/NACT.dat" u 1:($15/100) title 'd34/100' @wl
# rep "MD_OUTPUT/PES.dat" u 1:($3-$2) title 'E2' @wp
# rep "MD_OUTPUT/PES.dat" u 1:($4-$2) title 'E3' @wp
# rep "MD_OUTPUT/PES.dat" u 1:($5-$2) title 'E4' @wp
# rep "MD_OUTPUT/PES.dat" u 1:($6-$2) title 'E5' @wp


p "MD_OUTPUT/Average_PES.dat" u 1:2 @wp
rep "MD_OUTPUT/PES.dat" u 1:2 @wl
rep "MD_OUTPUT/PES.dat" u 1:3 @wl
rep "MD_OUTPUT/PES.dat" u 1:4 @wl
rep "MD_OUTPUT/PES.dat" u 1:5 @wl
rep "MD_OUTPUT/PES.dat" u 1:6 @wl
while (1) {pause 10; @see}
