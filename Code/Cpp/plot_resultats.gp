reset
set datafile separator "\t"
set key left top
set grid

# entrée \ sortie
input = "resultats/resultat.txt"
outdir = "plots"

system(sprintf("mkdir %s 2> NUL", outdir))

set key autotitle columnhead

set terminal pngcairo size 1400,900 font ",12"
set output sprintf("%s/resultats_multiplot.png", outdir)

set multiplot layout 2,3 title "Simulation mitochondrie (resultat.txt)"

set xlabel "t"
set ylabel "NADH_m"
plot input every ::1 using 1:2 with lines lw 2

set xlabel "t"
set ylabel "ADP_m"
plot input every ::1 using 1:3 with lines lw 2

set xlabel "t"
set ylabel "deltaPsi"
plot input every ::1 using 1:4 with lines lw 2

set xlabel "t"
set ylabel "Ca_m"
plot input every ::1 using 1:5 with lines lw 2

set xlabel "t"
set ylabel "Ca_c"
plot input every ::1 using 1:6 with lines lw 2

set xlabel "t"
set ylabel "FBP"
plot input every ::1 using 1:7 with lines lw 2

unset multiplot
set output

set terminal pngcairo size 1200,700 font ",12"

set output sprintf("%s/NADH_m.png", outdir)
set title "NADH_m vs t"
set xlabel "t"
set ylabel "NADH_m"
plot input every ::1 using 1:2 with lines lw 2
set output

set output sprintf("%s/ADP_m.png", outdir)
set title "ADP_m vs t"
set xlabel "t"
set ylabel "ADP_m"
plot input every ::1 using 1:3 with lines lw 2
set output

set output sprintf("%s/deltaPsi.png", outdir)
set title "deltaPsi vs t"
set xlabel "t"
set ylabel "deltaPsi"
plot input every ::1 using 1:4 with lines lw 2
set output

set output sprintf("%s/Ca_m.png", outdir)
set title "Ca_m vs t"
set xlabel "t"
set ylabel "Ca_m"
plot input every ::1 using 1:5 with lines lw 2
set output

# Done
set terminal pop
