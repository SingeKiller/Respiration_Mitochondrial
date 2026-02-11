set datafile separator "\t"
set term pngcairo size 1200,900

outdir = "plots"
system sprintf("mkdir -p %s", outdir)
set output sprintf("%s/resultats_graph.png", outdir)
set multiplot layout 3,2 title "Bertram model FBP variation"

#A: FBP
set title "FBP(t)"
set xlabel "temps"; set ylabel "FBP"
set yrange  [0 : 20]
plot "resultats/variation_FBP.txt" u 1:2 w l t "FBP=5", \
	"resultats/variation_FBP_10.txt" u 1:2 w l t "FBP=10",\
	"resultats/variation_FBP_15.txt" u 1:2 w l t "FBP=15"

#B: Ca_m
set title "Ca_m(t)"
set ylabel "Ca_m"
set yrange [ 0 : 0.3]
plot "resultats/variation_FBP.txt" u 1:5 w l t "FBP=5", \
	"resultats/variation_FBP_10.txt" u 1:5 w l t "FBP=10",\
	"resultats/variation_FBP_15.txt" u 1:5 w l t "FBP=15"


#C: NADH_m
set title "NADH(t)"
set xlabel "temps"; set ylabel "NADH_m"
set yrange [0 : 2000]
plot "resultats/variation_FBP.txt" u 1:4 w l t "FBP=5", \
	"resultats/variation_FBP_10.txt" u 1:4 w l t "FBP=10",\
	"resultats/variation_FBP_15.txt" u 1:4 w l t "FBP=15"


#D: dpsi
set title "dpsi(t)"
set xlabel "temps"; set ylabel "dpsi"
set yrange [ 140 : 180]
plot "resultats/variation_FBP.txt" u 1:6 w l t "FBP=5", \
	"resultats/variation_FBP_10.txt" u 1:6 w l t "FBP=10",\
	"resultats/variation_FBP_15.txt" u 1:6 w l t "FBP=15"


#E: J_o
set title "J_o(t)"
set xlabel "temps"; set ylabel "Jo"
set yrange [0 : 0.8]
plot "resultats/variation_FBP.txt" u 1:8 w l t "FBP=5", \
	"resultats/variation_FBP_10.txt" u 1:8 w l t "FBP=10",\
	"resultats/variation_FBP_15.txt" u 1:8 w l t "FBP=15"


#F: ATP
set title "ATP(t)"
set xlabel "temps"; set ylabel "ATP"
set yrange [ 0 : 6000]
plot "resultats/variation_FBP.txt" u 1:7 w l t "FBP=5", \
	"resultats/variation_FBP_10.txt" u 1:7 w l t "FBP=10",\
	"resultats/variation_FBP_15.txt" u 1:7 w l t "FBP=15"

unset multiplot
