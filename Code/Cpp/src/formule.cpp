#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include "formule.hpp"
#include "writer.hpp"
#include "solver.hpp"


std::vector<double> dxdt_(
    const std::vector<double>& x,
    entree I,
    const ModelParams& params,
    double t
){
    Bertram model;
    return model.dxdt(x, I, params, t);
}

void Bertram::sync_state(const std::vector<double>& x) {
    if (x.size() <= CA_m) {
        throw std::invalid_argument("Etat insuffisant: NADH, ADP, DPSI, CA_m requis");
    }
    nadh_ = x[NADH];
    adp_ = x[ADP];
    dpsi_ = x[DPSI];
    ca_m_ = x[CA_m];
}

double Bertram::NAD_m(const std::vector<double>& x, const ModelParams& params){
    sync_state(x);
    double NAD_m = params.NADtot - nadh_value();

    if (NAD_m < 1e-6){
        NAD_m = 1e-6;
    }
    return NAD_m;
}

double Bertram::ATP_M(const std::vector<double>& x, const ModelParams& params){
    sync_state(x);
    double ATP_m = params.Atot - adp_value();

    if (ATP_m < 0){
        ATP_m = 0;
    }
    return ATP_m;
}

// Calcul de tous les flux 

double Bertram::J_GPDH(double FBP, const ModelParams& params){
    if (FBP < 0){
        FBP = 0;
    }
    // JGPDH = kGPDH * FBP  (FBP en μM)
    double J_GPDH = params.kGPDH * sqrt(FBP);
    return J_GPDH;
}



double Bertram::J_PDH(const std::vector<double>& x, const ModelParams& params,double jgpdh){
    sync_state(x);
    double nad_m = this->NAD_m(x, params);
    double inibition = params.p[1]/(params.p[2] + nadh_value()/nad_m);
    double activation = ca_m_value()/(params.p[3] + ca_m_value());
    return inibition * activation * jgpdh;
}



double  Bertram::J_o(const std::vector<double>& x, const ModelParams& params){
    sync_state(x);
    double term1 = (params.p[4] * nadh_value()) /(params.p[5] + nadh_value());
    double term2 = 1 / (1 + exp(( dpsi_value() - params.p[6])/params.p[7]));
    return term1 * term2 ;
}



double Bertram::J_Hres(const std::vector<double>& x, const ModelParams& params){
    sync_state(x);
    double term1 = (params.p[8] * nadh_value()) /(params.p[9] + nadh_value());
    double term2 = 1 / (1 + exp(( dpsi_value() - params.p[10])/params.p[11]));
    return term1 * term2;
}


double Bertram::J_Hatp(const std::vector<double>& x, const ModelParams& params){
    sync_state(x);
    double ATP_m = this->ATP_M(x, params);
    double term1 = params.p[13] /(params.p[13] + ATP_m);
    double term2 = params.p[12] / (1 + exp((params.p[14] - dpsi_value() )/params.p[15]));
    return term1 * term2 ;
}



double Bertram::J_F1F0(const std::vector<double>& x, const ModelParams& params){
    sync_state(x);
    double ATP_m = this->ATP_M(x, params);
    double term1 = params.p[13] /(params.p[13] + ATP_m);
    double term2 = params.p[16] / (1 + exp((params.p[14]-dpsi_value() )/params.p[15]));
    return term1 * term2;
}



double Bertram::J_Hleak(const std::vector<double>& x, const ModelParams& params){
    sync_state(x);
    double Jhleak= params.p[17] * dpsi_value() + params.p[18];
    if (Jhleak < 0){
        Jhleak = 0;
    }
    return Jhleak;
}



double Bertram::J_ANT(const std::vector<double>& x, const ModelParams& params){
    sync_state(x);
    double ATP_m = this->ATP_M(x, params);
    double adp_m = adp_value();
    if (adp_m < 1e-16){
        adp_m = 1e-16;
    }
    double RAT_m = ATP_m / adp_m;
    double term1 = params.p[19];
    double term2 = RAT_m / (params.p[20] + RAT_m);
    double term3 = exp((.5 * dpsi_value() * params.FRT));
    return term1 * term2 * term3;
}



double Bertram::J_uni(const std::vector<double>& x, const ModelParams& params,double Ca_c){
    sync_state(x);
    if (Ca_c < 0){
        Ca_c = 0;
    }
    double J_uni = ((params.p[21] * dpsi_value()) - params.p[22])*Ca_c*Ca_c;
    if (J_uni < 0){
        J_uni = 0;
    }
    return J_uni;
}


double Bertram::J_NaCa(const std::vector<double>& x, const ModelParams& params,double Ca_c){
    sync_state(x);
    double ca_m = ca_m_value();
    if (ca_m < 0){
        ca_m = 0;
    }

    if (Ca_c < 1e-6){
        Ca_c = 1e-6;
    }

    double Jnaca = params.p[23]*(ca_m / Ca_c)* exp(params.p[24] * dpsi_value());
    return Jnaca;
}


std::vector<double> Bertram::dxdt(const std::vector<double>& x,entree I, const ModelParams& params, double t){
    (void)t;
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
    
    std::vector<double> dx(taille);
    dx[NADH] =  ( jpdh - jo);
    dx[ADP] =  (jant - jf1f0);
    dx[DPSI] = (jhres - jhatp - jant - jhleak - jnaca - 2*juni)/params.Cm;
    dx[CA_m] = params.fm * (juni - jnaca);
    return dx;
}

void Bertram::write(std::ofstream& sortie, double t, double FBP, double Ca_c, double NADH_m, double Ca_m, double DPSI, double ATP_m, double J_o){
    sortie << ms_to_min(t) << "\t" 
            << FBP << "\t" 
            << Ca_c << "\t" 
            << NADH_m << "\t" 
            << Ca_m << "\t" 
            << DPSI << "\t" 
            << ATP_m << "\t" 
            << J_o << std::endl;
}