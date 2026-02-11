#include "formule.hpp"
#include "writer.hpp"
#include "solver.hpp"
#include <cmath>
#include <fstream>
#include<iostream>
#include <iomanip>
#include <map>
#include <vector>


std::vector<double> dxdt_exp(
    const std::vector<double>& x,
    entree /*I*/,
    const std::map<std::string, double>& params,
    double /*t*/
) {
    const double k = params.at("k");
    return { k * x[0] };
}

void teste_solver_convergence() {
    const double x0 = 1.0;
    const double k = 100.0;
    const double t_final = 0.001; // plus court pour eviter overflow avec k eleve

    std::map<std::string, double> test;
    test["k"] = k;

    entree I;
    I.Ca_c = 0.0;
    I.FBP = 0.0;

    const std::vector<double> dts = {0.001, 0.0005, 0.00025};
    std::ofstream err_out("resultats/solver_errors.txt");
    // ici on a un fichier d'erreur qui va nous permettre de suivre la convergence du solver en fonction du pas de temps
    err_out << "dt\t" << "max_err\t" << "rms_err" << std::endl;

    // boucle sur les différents pas de temps
    for (std::size_t i = 0; i < dts.size(); ++i) {
        const double dt_test = dts[i];
        std::vector<double> x(1, x0);
        double max_err = 0.0;
        double sum_sq = 0.0;
        std::size_t n = 0;

        // on génère un fichier de sortie pour chaque pas de temps
        std::string fichier = "resultats/solver_test_dt" + std::to_string(i + 1) + ".txt";
        std::ofstream out(fichier);
        out << "t\t" << "x_num\t" << "x_exact\t" << "ln_num\t" << "ln_exact\t" << "err" << std::endl;

        for (double t = 0.0; t <= t_final + 1e-12; t += dt_test) {
            const double x_exact = x0 * std::exp(k * t);
            const double ln_num = std::log(x[0] / x0);
            const double ln_exact = std::log(x_exact / x0);
            const double err = std::abs(x[0] - x_exact);

            out << std::setprecision(10)
                << t << "\t" << x[0] << "\t" << x_exact << "\t" << ln_num << "\t" << ln_exact
                << "\t" << err << std::endl;

            if (err > max_err) {
                max_err = err;
            }
            sum_sq += err * err;
            ++n;

            x = solver(x, I, test, t, dt_test, dxdt_exp);
        }

        const double rms_err = std::sqrt(sum_sq / n);
        
        err_out << std::setprecision(10) << dt_test << "\t" << max_err << "\t" << rms_err << std::endl;
        out.close();
    }
    err_out.close();
}