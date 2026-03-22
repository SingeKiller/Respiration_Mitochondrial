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

    double ADP_c_val = x[ADP_c]; // ADP_c en uM

    if (ADP_c_val < 1e-16) {
        ADP_c_val = 1e-16;
    }

    const double RAT_c = ATP_c(x, params) / ADP_c_val; // ATP_c calculé à partir de ADP_c

    double denom = (1.0 + params.kANT_c_den * RAT_c) *(RAT_m + params.kANT_m_den) * exp(-0.5 * params.FRT * dpsi_value());

    if (std::abs(denom) < 1e-16) {
        denom = (denom < 0.0) ? -1e-16 : 1e-16; // division zéros ? 
    }

    return params.VmaxANT * (RAT_m - params.kANT_c_num * RAT_c) / denom;
}


// Calcul de d/dt(J_ADP,ext)

double BE::J_ADP_ext(const ModelParams& params, double t) {
    if (params.tau_plus <= 0.0) {
        return 10e-16;
    }
    double somme = 0.0;
    
    double a = params.tau_plus *60000.0; // conversion en ms
    for (const ADP_ajout& ajout : params.adp_ajouts) {
        double ADP_c_plus_i = ajout.ADP_c_plus;
        double t_plus_i = ajout.t_plus * 60000.0;
        
        if (t > t_plus_i) {
            double terme = (ADP_c_plus_i)*(1 -  exp(-(t - t_plus_i) / a));
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
    dx[ADP_c] = - jant + j_adp_ext;
    
    return dx;
}

void BE::simulate_Jo_timecourse(std::vector<double> etat_initial,
                                const ModelParams& params,
                                entree I,
                                double dt_ms,
                                double tfinal_min,
                                const std::string& output_file) {
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
        output_file);
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

    const double pulse1_time_ms = pulse1_time_min * 60000.0;
    const double pulse2_time_ms = pulse2_time_min * 60000.0;
    

    std::vector<double> etat = etat_initial;
    bool pulse1_done = false;
    bool pulse2_done = false;

    
    for (double t = 0.0; t <= tfinal_ms; t += dt_ms){
         if (!pulse1_done && t >= pulse1_time_ms) {
            etat[ADP_c] += pulse1_adp_c;
            pulse1_done = true;
        }
        if (!pulse2_done && t >= pulse2_time_ms) {
            etat[ADP_c] += pulse2_adp_c;
            pulse2_done = true;
        }

        sortie << ms_to_min(t) << "\t" << this->J_o(etat, params) << std::endl;
        etat = solver(etat, I, params, t, dt_ms, dxdt_be);
    }
}

