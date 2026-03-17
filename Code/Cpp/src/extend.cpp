#include "extend.hpp"

#include <filesystem>
#include <fstream>
#include <functional>
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

    double ADP_c_val = x[ADP_c];
    if (ADP_c_val < 1e-16) {
        ADP_c_val = 1e-16;
    }
    const double RAT_c = ATP_c(x, params) / ADP_c_val;

    double denom = (1.0 + params.kANT_c_den * RAT_c) *
                   (RAT_m + params.kANT_m_den) *
                   exp(-0.5 * params.FRT * dpsi_value());
    if (std::abs(denom) < 1e-16) {
        denom = (denom < 0.0) ? -1e-16 : 1e-16;
    }

    return params.VmaxANT * (RAT_m - params.kANT_c_num * RAT_c) / denom;
}


// Calcul de d/dt(J_ADP,ext)
double BE::J_ADP_ext(const ModelParams& params, double t) {
    if (params.tau_plus <= 0.0) {
        return 0.0;
    }
    double somme = 0.0;

    for (const ADP_ajout& ajout : params.adp_ajouts) {
        double ADP_c_plus_i = ajout.ADP_c_plus;
        double t_plus_i = ajout.t_plus;
        
        if (t > t_plus_i) {
            double terme = (ADP_c_plus_i / params.tau_plus) * exp(-(t - t_plus_i) / params.tau_plus);
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
    double jant = this->J_ANT(x,params);
    double juni = this->J_uni(x,params,I.Ca_c);
    double jnaca = this->J_NaCa(x,params,I.Ca_c);

    // Nouveau flux exterieur: derivee de J_ADP,ext
    double j_adp_ext = this->J_ADP_ext(params, t);
    
    std::vector<double> dx(taille_etendue);
    dx[NADH] = (jpdh - jo);
    dx[ADP] = (jant - jf1f0);
    dx[DPSI] = (jhres - jhatp - jant - jhleak - jnaca - 2*juni) / params.Cm;
    dx[CA_m] = params.fm * (juni - jnaca);
    dx[ADP_c] = -jant + j_adp_ext;
    
    return dx;
}

void BE::simulate_Jo_timecourse(std::vector<double> etat_initial,
                                const ModelParams& params,
                                entree I,
                                double dt_ms,
                                double tfinal_min,
                                const std::string& output_file) {

    if (dt_ms <= 0.0) {
        throw std::invalid_argument("dt_ms doit etre > 0");
    }

    if (etat_initial.size() < taille_etendue) {
        etat_initial.resize(taille_etendue, 0.0);
    }

    std::filesystem::path out_path(output_file);
    if (!out_path.parent_path().empty()) {
        std::filesystem::create_directories(out_path.parent_path());
    }

    std::ofstream sortie(output_file.c_str());
    if (!sortie.is_open()) {
        throw std::runtime_error("Impossible d'ouvrir le fichier de sortie Jo");
    }

    const std::function<std::vector<double>(const std::vector<double>&, entree, const ModelParams&, double)> dxdt_be =
        [this](const std::vector<double>& x, entree input, const ModelParams& p, double t){
            return this->dxdt(x, input, p, t);
        };

    sortie << "t_min\tJ_o" << std::endl;

    const double tfinal_ms = tfinal_min * 60000.0;
    std::vector<double> etat = etat_initial;

    for (double t = 0.0; t <= tfinal_ms; t += dt_ms){
        sortie << ms_to_min(t) << "\t" << this->J_o(etat, params) << std::endl;
        etat = solver(etat, I, params, t, dt_ms, dxdt_be);
    }
}

void BE::simulate_o2_normalized_integrated_with_pulses(std::vector<double> etat_initial,
                                const ModelParams& params,
                                entree I,
                                double dt_ms,
                                double tfinal_min,
                                double pulse1_time_min,
                                double pulse1_adp_c,
                                double pulse2_time_min,
                                double pulse2_adp_c,
                                const std::string& output_file) {

    if (dt_ms <= 0.0) {
        throw std::invalid_argument("dt_ms doit etre > 0");
    }

    if (etat_initial.size() < taille_etendue) {
        etat_initial.resize(taille_etendue, 0.0);
    }

    std::filesystem::path out_path(output_file);
    if (!out_path.parent_path().empty()) {
        std::filesystem::create_directories(out_path.parent_path());
    }

    std::ofstream sortie(output_file.c_str());
    if (!sortie.is_open()) {
        throw std::runtime_error("Impossible d'ouvrir le fichier de sortie O2 normalise");
    }

    const std::function<std::vector<double>(const std::vector<double>&, entree, const ModelParams&, double)> dxdt_be =
        [this](const std::vector<double>& x, entree input, const ModelParams& p, double t){
            return this->dxdt(x, input, p, t);
        };

    ModelParams sim_params = params;
    sim_params.adp_ajouts.clear();

    const double tfinal_ms = tfinal_min * 60000.0;
    const double pulse1_time_ms = pulse1_time_min * 60000.0;
    const double pulse2_time_ms = pulse2_time_min * 60000.0;

    std::vector<double> etat = etat_initial;
    std::vector<double> t_vals_min;
    std::vector<double> jo_vals;
    bool pulse1_done = false;
    bool pulse2_done = false;

    for (double t = 0.0; t <= tfinal_ms; t += dt_ms) {
        if (!pulse1_done && t >= pulse1_time_ms) {
            etat[ADP_c] += pulse1_adp_c;
            pulse1_done = true;
        }
        if (!pulse2_done && t >= pulse2_time_ms) {
            etat[ADP_c] += pulse2_adp_c;
            pulse2_done = true;
        }

        t_vals_min.push_back(ms_to_min(t));
        jo_vals.push_back(this->J_o(etat, sim_params));
        etat = solver(etat, I, sim_params, t, dt_ms, dxdt_be);
    }

    if (jo_vals.empty()) {
        throw std::runtime_error("Aucune donnee Jo calculee pour integration");
    }

    std::vector<double> integ(jo_vals.size(), 0.0);
    for (std::size_t i = 1; i < jo_vals.size(); ++i) {
        const double dt_min = t_vals_min[i] - t_vals_min[i - 1];
        const double trap = 0.5 * (jo_vals[i] + jo_vals[i - 1]) * dt_min;
        integ[i] = integ[i - 1] + trap;
    }

    std::vector<double> o2_proxy(jo_vals.size(), 1.0);
    for (std::size_t i = 0; i < jo_vals.size(); ++i) {
        o2_proxy[i] = 1.0 - integ[i];
    }

    const double o2_start = o2_proxy.front();
    const double o2_end = o2_proxy.back();
    const double denom = o2_start - o2_end;

    sortie << "t_min\tO2_sim_norm" << std::endl;
    for (std::size_t i = 0; i < o2_proxy.size(); ++i) {
        double o2_norm = 1.0;
        if (std::abs(denom) > 1e-16) {
            o2_norm = (o2_proxy[i] - o2_end) / denom;
        }
        sortie << t_vals_min[i] << "\t" << o2_norm << std::endl;
    }
}