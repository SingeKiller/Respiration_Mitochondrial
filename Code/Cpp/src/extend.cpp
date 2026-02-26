#include "extend.hpp"

// Calcul de  ATP avec  la relation de conservation
double BE::ATP_c(const std::vector<double>& x, const std::map<std::string,double>& params) {
    double ATP_c = params.at("Atot_c") - x[ADP_c];
    if (ATP_c < 0) {
        ATP_c = 0;
    }
    return ATP_c;
}

// Calcul de J_ADP_ext
double BE::J_ADP_ext(const std::vector<double>& x, 
                     const std::map<std::string,double>& params,
                     double t) {
    
    double ADP_c_0 = params.at("ADP_c_0");   
    double tau_plus = params.at("tau_plus"); 
    double somme = ADP_c_0;
    int n_ajouts = params.at("n_ajouts");
    
    for (int i = 0; i < n_ajouts; i++) {

        double ADP_c_plus_i = params.at("ADP_c_plus_" + std::to_string(i));
        double t_plus_i = params.at("t_plus_" + std::to_string(i));
        
        if (t > t_plus_i) {
            double terme = ADP_c_plus_i * (1.0 - exp(-(t - t_plus_i) / tau_plus));
            somme += terme;
        }
    }
    
    return somme;
}

// dxdt pour le modèle étendu avec ADP_c
std::vector<double> BE::dxdt(const std::vector<double>& x,
                              entree I, 
                              std::map<std::string,double> params, 
                              double t) {
    params["FBP"] = I.FBP;

    // Calcul de tous les flux hérités de Bertram

    double jgpdh = this->J_GPDH(params);
    double jpdh = this->J_PDH(x,params,jgpdh);
    double jo = this->J_o(x,params);
    double jhres = this->J_Hres(x,params);
    double jhatp = this->J_Hatp(x,params);
    double jf1f0 = this->J_F1F0(x,params);
    double jhleak = this->J_Hleak(x,params);
    double jant = this->J_ANT(x,params);
    double juni = this->J_uni(x,params,I.Ca_c);
    double jnaca = this->J_NaCa(x,params,I.Ca_c);
    // le nouveaux Flux
    double j_adp_ext = this->J_ADP_ext(x, params, t);
    
    std::vector<double> dx(taille_etendue);
    dx[NADH] = (jpdh - jo);
    dx[ADP] = (jant - jf1f0);
    dx[DPSI] = (jhres - jhatp - jant - jhleak - jnaca - 2*juni) / params.at("Cm");
    dx[CA_m] = params.at("fm") * (juni - jnaca);
    dx[ADP_c] = -jant + j_adp_ext;
    
    return dx;
}