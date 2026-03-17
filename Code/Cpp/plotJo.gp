set datafile separator "\t"
set term pngcairo size 1200,900

outdir = "plots"
system sprintf("mkdir -p %s", outdir)
set output sprintf("%s/plotJo.png", outdir)

set multiplot layout 2,1 title "BE: Jo(t) et Jo normalise"

set title "Jo(t)"
set xlabel "temps (min)"
set ylabel "J_o"
set xrange [14:25]
set grid
plot "resultats/jo_BE_temps.txt" u ($1+14):2 w l lw 2 t "Jo"

set title "Jo normalise(t)"
set xlabel "temps (min)"
set ylabel "J_o_norm"
set yrange [0:1]
set xrange [14:25]
set grid
plot "resultats/jo_normalized.txt" u ($1+14):2 w l lw 2 t "Jo normalise"

unset multiplot
