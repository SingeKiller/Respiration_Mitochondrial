#include "extend.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <stdexcept>


double BE::ATP_c(const std::vector<double>& x, const ModelParams& params) {
    if (x.size() <= ADP_c) {
        throw std::invalid_argument("Etat insuffisant: ADP_c requis pour BE");
    }
    double ATP_c = params.Atot_c - x[ADP_c];
    if (ATP_c < 0) {
        ATP_c = 0;
    }
    return ATP_c;
}

double BE::J_ANT(const std::vector<double>& x, const ModelParams& params) {
    sync_state(x);
    if (x.size() <= ADP_c) {
        throw std::invalid_argument("Etat insuffisant: ADP_c requis pour J_ANT equation 35");
    }
    
    double ATP_m = this->ATP_M(x, params);
    double ADP_m = adp_value();
    if (ADP_m < 1e-16) {
        ADP_m = 1e-16;
    }
    const double RAT_m = ATP_m / ADP_m;

    double ADP_c_val = x[ADP_c]; // ADP_c en uM

    if (ADP_c_val < 1e-12) {
        ADP_c_val = 1e-12;
    }

    const double RAT_c = ATP_c(x, params) / ADP_c_val; // ATP_c calculé à partir de ADP_c

    double denom = (1.0 + params.kANT_c_den * RAT_c) *(RAT_m + params.kANT_m_den) * exp(-0.5 * params.FRT * dpsi_value());

    if (std::abs(denom) < 1e-20) {
        denom = (denom < 0.0) ? -1e-20 : 1e-20; // garde-fou numerique plus souple
    }

    return params.VmaxANT * (RAT_m - params.kANT_c_num * RAT_c) / denom;
}


// Flux externe d'ADP_c: somme cumulative des ajouts

double BE::J_ADP_ext_sum(const ModelParams& params, double t) {
    if (params.tau_plus <= 0.0) {
        return 0.0;
    }

    // Convention: t, t_plus_i et tau_plus sont en ms.
    const double tau_ms = params.tau_plus;
    double somme = 0.0;

    for (const ADP_ajout& ajout : params.adp_ajouts) {
        double ADP_c_plus_i = ajout.ADP_c_plus;
        double t_plus_i = ajout.t_plus ;

        if (t > t_plus_i) {
            const double dt_ms = t - t_plus_i;
            const double terme = ADP_c_plus_i * (1.0 - exp(-dt_ms / tau_ms));
            somme += terme;
        }
    }

    return somme;
}

// Derivee temporelle de J_ADP_ext_sum, utilisee dans d[ADP_c]/dt
double BE::J_ADP_ext_dot(const ModelParams& params, double t) {
    if (params.tau_plus <= 0.0) {
        return 0.0;
    }

    const double tau_ms = params.tau_plus;
    double somme = 0.0;

    for (const ADP_ajout& ajout : params.adp_ajouts) {
        const double ADP_c_plus_i = ajout.ADP_c_plus;
        const double t_plus_i = ajout.t_plus;

        if (t > t_plus_i) {
            const double dt_ms = t - t_plus_i;
            const double terme = (ADP_c_plus_i / tau_ms) * exp(-dt_ms / tau_ms);
            somme += terme;
        }
    }

    return somme;
}

std::vector<double> BE::dxdt(const std::vector<double>& x,
                              entree I, 
                              const ModelParams& params,
                              double t) {

    double jgpdh = this->J_GPDH(I.FBP, params);
    double jpdh = this->J_PDH(x,params,jgpdh);
    double jo = this->J_o(x,params);
    double jhres = this->J_Hres(x,params);
    double jhatp = this->J_Hatp(x,params);
    double jf1f0 = this->J_F1F0(x,params);
    double jhleak = this->J_Hleak(x,params);
    const double jant_m = this->J_ANT(x,params);
    const double jant_c = params.gamma * jant_m;
    double juni = this->J_uni(x,params,I.Ca_c);
    double jnaca = this->J_NaCa(x,params,I.Ca_c);

    // Flux source exterieur: derivee temporelle de J_ADP_ext_sum
    double j_adp_ext = this->J_ADP_ext_dot(params, t);
    
    std::vector<double> dx(taille_etendue);
    dx[NADH] = (jpdh - jo);
    dx[ADP] = (jant_m - jf1f0);
    dx[DPSI] = (jhres - jhatp - jant_m - jhleak - jnaca - 2*juni) / params.Cm;
    dx[CA_m] = params.fm * (juni - jnaca);
    // Schema Bertram etendu pour le compartiment cytosolique.
    dx[ADP_c] = - jant_c + j_adp_ext;
    
    return dx;
}

void BE::simulate_Jo_timecourse(std::vector<double> etat_initial,
                                const ModelParams& params,
                                entree I,
                                double dt_ms,
                                double tfinal_min,
                                const std::string& output_file,
                                const std::string& flux_output_file) {
    simulate_Jo_timecourse(
        etat_initial,
        params,
        I,
        dt_ms,
        tfinal_min,
        0.0,
        0.0,
        0.0,
        0.0,
        output_file,
        flux_output_file);
}

void BE::simulate_Jo_timecourse(std::vector<double> etat_initial,
                                const ModelParams& params,
                                entree I,
                                double dt_ms,
                                double tfinal_min,
                                double pulse1_time_min,
                                double pulse1_adp_c,
                                double pulse2_time_min,
                                double pulse2_adp_c,
                                const std::string& output_file,
                                const std::string& flux_output_file) {

    if (dt_ms <= 0.0) {
        throw std::invalid_argument("dt_ms doit etre > 0");
    }

    if (etat_initial.size() < taille_etendue) {
        etat_initial.resize(taille_etendue, 0.0);
    }

    // Rend la simulation robuste: si aucun ajout n'est pre-rempli,
    // on construit la sommation avec les pulses explicites passes au main.
    ModelParams params_effective = params;
    if (params_effective.adp_ajouts.empty()) {
        if (pulse1_adp_c > 0.0) {
            params_effective.adp_ajouts.push_back({pulse1_adp_c, pulse1_time_min * 60000.0});
        }
        if (pulse2_adp_c > 0.0) {
            params_effective.adp_ajouts.push_back({pulse2_adp_c, pulse2_time_min * 60000.0});
        }
    }

    std::filesystem::path out_path(output_file);
    if (!out_path.parent_path().empty()) {
        std::filesystem::create_directories(out_path.parent_path());
    }

    std::ofstream sortie(output_file.c_str());
    if (!sortie.is_open()) {
        throw std::runtime_error("Impossible d'ouvrir le fichier de sortie Jo");
    }

    const std::string effective_flux_output_file = 
            flux_output_file.empty() ? "resultats/fluxes_temps.txt" : flux_output_file;

    std::filesystem::path flux_out_path(effective_flux_output_file);
    if (!flux_out_path.parent_path().empty()) {
        std::filesystem::create_directories(flux_out_path.parent_path());
    }

    std::ofstream flux_sortie(effective_flux_output_file.c_str());
    if (!flux_sortie.is_open()) {
        throw std::runtime_error("Impossible d'ouvrir le fichier de sortie des flux");
    }
    // Keep enough decimals to preserve monotonic time stamps after ms->min conversion.
    sortie << std::fixed << std::setprecision(10);
    flux_sortie << std::fixed << std::setprecision(10);
    flux_sortie << "t_min\tJ_ANT\tJ_ADP_ext_sum\tJ_ADP_ext_dot\tADP_c" << std::endl;

    const std::function<std::vector<double>(const std::vector<double>&, entree, const ModelParams&, double)> dxdt_be =
        [this](const std::vector<double>& x, entree input, const ModelParams& p, double t){
            return this->dxdt(x, input, p, t);
        };
    sortie << "t_min\tJ_o" << std::endl;

    const double tfinal_ms = tfinal_min * 60000.0;
    // Internal sub-stepping improves robustness for stiff transients while
    // preserving the requested output cadence at dt_ms.
    const double max_internal_dt_ms = 1.0;
    const int n_substeps = std::max(1, static_cast<int>(std::ceil(dt_ms / max_internal_dt_ms)));
    const double dt_internal_ms = dt_ms / static_cast<double>(n_substeps);

    const double pulse1_time_ms = pulse1_time_min * 60000.0;
    const double pulse2_time_ms = pulse2_time_min * 60000.0;


    std::vector<double> etat = etat_initial;
    const double adp_c_min = 1e-12;
    const double adp_c_max = (params_effective.Atot_c > adp_c_min) ? 10.0 * params_effective.Atot_c : 1e12;
    etat[ADP_c] = std::clamp(etat[ADP_c], adp_c_min, adp_c_max);

    bool pulse1_done = false;
    bool pulse2_done = false;
    const bool use_discrete_pulses = params_effective.adp_ajouts.empty();

    for (double t = 0.0; t <= tfinal_ms; t += dt_ms){
        if (use_discrete_pulses) {
            if (!pulse1_done && t >= pulse1_time_ms) {
                etat[ADP_c] += pulse1_adp_c;
                pulse1_done = true;
            }
            if (!pulse2_done && t >= pulse2_time_ms) {
                etat[ADP_c] += pulse2_adp_c;
                pulse2_done = true;
            }
        }

        sortie << ms_to_min(t) << "\t" << this->J_o(etat, params_effective) << std::endl;
        flux_sortie << ms_to_min(t)
                   << "\t" << this->J_ANT(etat, params_effective)
                   << "\t" << this->J_ADP_ext_sum(params_effective, t)
                   << "\t" << this->J_ADP_ext_dot(params_effective, t)
                   << "\t" << etat[ADP_c]
                   << std::endl;
        double t_internal = t;
        for (int k = 0; k < n_substeps; ++k) {
            etat = solver(etat, I, params_effective, t_internal, dt_internal_ms, dxdt_be);
            t_internal += dt_internal_ms;
        }
        if (!std::isfinite(etat[ADP_c])) {
            etat[ADP_c] = adp_c_min;
        }
        etat[ADP_c] = std::clamp(etat[ADP_c], adp_c_min, adp_c_max);
    }
}

