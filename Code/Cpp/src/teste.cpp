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
    err_out << "dt\t" << "max_err_rk4\t" << "rms_err_rk4\t" << "max_err_ie\t" << "rms_err_ie" << std::endl;

    // boucle sur les différents pas de temps
    for (std::size_t i = 0; i < dts.size(); ++i) {
        const double dt_test = dts[i];
        std::vector<double> x_rk4(1, x0);
        std::vector<double> x_ie(1, x0);
        double max_err_rk4 = 0.0;
        double sum_sq_rk4 = 0.0;
        double max_err_ie = 0.0;
        double sum_sq_ie = 0.0;
        std::size_t n = 0;

        // on génère un fichier de sortie pour chaque pas de temps
        std::string fichier = "resultats/solver_test_dt" + std::to_string(i + 1) + ".txt";
        std::ofstream out(fichier);
        out << "t\t" << "x_rk4\t" << "x_ie\t" << "x_exact\t" << "ln_rk4\t" << "ln_ie\t" << "ln_exact\t" << "err_rk4\t" << "err_ie" << std::endl;

        for (double t = 0.0; t <= t_final + 1e-12; t += dt_test) {
            const double x_exact = x0 * std::exp(k * t);
            const double ln_rk4 = std::log(x_rk4[0] / x0);
            const double ln_ie = std::log(x_ie[0] / x0);
            const double ln_exact = std::log(x_exact / x0);
            const double err_rk4 = std::abs(x_rk4[0] - x_exact);
            const double err_ie = std::abs(x_ie[0] - x_exact);

            out << std::setprecision(10)
                << t << "\t" << x_rk4[0] << "\t" << x_ie[0] << "\t" << x_exact << "\t"
                << ln_rk4 << "\t" << ln_ie << "\t" << ln_exact << "\t" << err_rk4 << "\t" << err_ie << std::endl;

            if (err_rk4 > max_err_rk4) {
                max_err_rk4 = err_rk4;
            }
            if (err_ie > max_err_ie) {
                max_err_ie = err_ie;
            }
            sum_sq_rk4 += err_rk4 * err_rk4;
            sum_sq_ie += err_ie * err_ie;
            ++n;

            x_rk4 = solver(x_rk4, I, test, t, dt_test, dxdt_exp);
            x_ie[0] = x_ie[0] / (1.0 - k * dt_test);
        }

        const double rms_err_rk4 = std::sqrt(sum_sq_rk4 / n);
        const double rms_err_ie = std::sqrt(sum_sq_ie / n);

        err_out << std::setprecision(10) << dt_test << "\t" << max_err_rk4 << "\t" << rms_err_rk4
                << "\t" << max_err_ie << "\t" << rms_err_ie << std::endl;
        out.close();
    }
    err_out.close();
}