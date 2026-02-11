set datafile separator "\t"
set term pngcairo size 1200,900

outdir = "plots"
system sprintf("mkdir -p %s", outdir)
set output sprintf("%s/solver_convergence.png", outdir)

set multiplot layout 2,2 title "Test solver: convergence"

# A: x(t) num vs exact 
set title "x(t) : num vs exact "
set xlabel "t"
set ylabel "x"
plot "resultats/solver_test_dt3.txt" u 1:2 w l t "solver", \
     "resultats/solver_test_dt3.txt" u 1:3 w l t "valeur exacte"

# B: log pour pente (dt1)
set title "log(x/x0) : dt1"
set xlabel "t"
set ylabel "ln(x/x0)"
plot "resultats/solver_test_dt1.txt" u 1:4 w l t "ln numerique", \
     "resultats/solver_test_dt1.txt" u 1:5 w l t "ln exacte"

# C: erreur max vs dt (log-log)
set title "Erreur max vs dt (log-log)"
set xlabel "dt"
set ylabel "max erreur"
set logscale xy
plot "resultats/solver_errors.txt" u 1:2 w lp t "max erreur"

# D: erreur moyenne vs dt (log-log)
set title "Erreur erreur moyenne vs dt (log-log)"
set xlabel "dt"
set ylabel "rms_err"
set logscale xy
plot "resultats/solver_errors.txt" u 1:3 w lp t "rms_err"

unset multiplot
unset logscale
