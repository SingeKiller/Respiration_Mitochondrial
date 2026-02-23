set datafile separator "\t"
set term pngcairo size 1200,900

outdir = "plots"
system sprintf("mkdir -p %s", outdir)
set output sprintf("%s/solver_convergence.png", outdir)

set multiplot layout 2,2 title "Test solver: convergence"

# A: x(t) num vs exact 
set title "x(t) : RK4 vs Euler implicite vs exact"
set xlabel "t"
set ylabel "x"
plot "resultats/solver_test_dt3.txt" u 1:2 w l t "RK4", \
     "resultats/solver_test_dt3.txt" u 1:3 w l t "Euler implicite", \
     "resultats/solver_test_dt3.txt" u 1:4 w l t "valeur exacte"

# B: log pour pente (dt1)
set title "log(x/x0) : dt1"
set xlabel "t"
set ylabel "ln(x/x0)"
plot "resultats/solver_test_dt1.txt" u 1:5 w l t "ln RK4", \
     "resultats/solver_test_dt1.txt" u 1:6 w l t "ln Euler implicite", \
     "resultats/solver_test_dt1.txt" u 1:7 w l t "ln exact"

# C: erreur max vs dt (log-log)
set title "Erreur max vs dt (log-log)"
set xlabel "dt"
set ylabel "max erreur"
set logscale xy
set format x "10^{%L}"
set format y "10^{%L}"
set grid
set xrange [2.5e-4:1.0e-3]
set yrange [1e-10:1e-2]
plot "resultats/solver_errors.txt" u 1:2 w lp t "max erreur RK4", \
     "resultats/solver_errors.txt" u 1:4 w lp t "max erreur Euler implicite"

# D: erreur moyenne vs dt (log-log)
set title "Erreur Norme 2 vs dt "
set xlabel "dt"
set ylabel "norme 2_err"
set logscale xy
set format x "10^{%L}"
set format y "10^{%L}"
set grid
set xrange [2.5e-4:1.0e-3]
set yrange [1e-10:1e-2]
plot "resultats/solver_errors.txt" u 1:3 w lp t "rms RK4", \
     "resultats/solver_errors.txt" u 1:5 w lp t "rms Euler implicite"

unset multiplot
unset logscale
