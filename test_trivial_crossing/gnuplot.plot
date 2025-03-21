
p "PES.dat" u 1:($3-$2)   title "E1" @wp
rep "PES.dat" u 1:($4-$2) title "E2" @wp
rep "Population.dat" u 1:2 title "P0" @wl
rep "Population.dat" u 1:3 title "P1" @wl
rep "Population.dat" u 1:4 title "P2" @wl
rep "Overlap.dat" u 1:3 title "S01" @wl
rep "Overlap.dat" u 1:4 title "S02" @wl
rep "Overlap.dat" u 1:6 title "S12" @wl

while (1) {@see; pause 10}