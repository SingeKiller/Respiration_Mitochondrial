set datafile separator "\t"
set term pngcairo size 1200,900

outdir = "plots"
system sprintf("mkdir -p %s", outdir)
set output sprintf("%s/resultats_graph13.png", outdir)
set multiplot layout 3,2 title "Bertram model Ca_c variation"

#A: Ca_c
set title "Ca_c(t)"
set xlabel "temps"; set ylabel "Ca_c"
set yrange [0:0.6]
plot "resultats/variation_Ca_c.txt" u 1:3 w l t "Ca_c"

#B: Ca_m
set title "Ca_m(t)"
set ylabel "Ca_m"
set yrange [0:6]
plot "resultats/variation_Ca_c.txt" u 1:5 w l t "Ca_m"


#C: NADH_m
set title "NADH(t)"
set xlabel "temps"; set ylabel "NADH_m"
set yrange [100:500]
plot "resultats/variation_Ca_c.txt" u 1:4 w l t "NADH_m"


#D: dpsi
set title "dpsi(t)"
set xlabel "temps"; set ylabel "dpsi"
set yrange [155:165]
plot "resultats/variation_Ca_c.txt" u 1:6 w l t "deltaPsi"


#E: J_o
set title "J_o(t)"
set xlabel "temps"; set ylabel "Jo"
set yrange [0.3:0.5]
plot "resultats/variation_Ca_c.txt" u 1:8 w l t "J_o"

#F: ATP
set title "ATP(t)"
set xlabel "temps"; set ylabel "ATP"
set yrange [3000:4000]
plot "resultats/variation_Ca_c.txt" u 1:7 w l t "ATP"

unset multiplot
