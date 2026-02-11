#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <cmath>
#include "formule.hpp"
#include "writer.hpp"
#include "solver.hpp"


double Bertram::NAD_m(const std::vector<double>& x, const std::map<std::string, double>& params){
    double NAD_m = params.at("NADtot") - x[NADH];

    if (NAD_m < 1e-6){
        NAD_m = 1e-6;
    }
    return NAD_m;
}

double Bertram::ATP_M(const std::vector<double>& x, const std::map<std::string,double>& params){
    double ATP_m = params.at("Atot") - x[ADP];

    if (ATP_m < 0){
        ATP_m = 0;
    }
    return ATP_m;
}

// Calcul de tous les flux 

double Bertram::J_GPDH(std::map<std::string, double> params){
    if (params["FBP"] < 0){
        params["FBP"] = 0;
    }
    // JGPDH = kGPDH * FBP  (FBP en μM)
    double J_GPDH = params["kGPDH"] * sqrt(params["FBP"]);
    return J_GPDH;
}



double Bertram::J_PDH(const std::vector<double>& x, const std::map<std::string, double>& params,double jgpdh){
    double nad_m = this->NAD_m(x, params);
    double inibition = params.at("p1")/(params.at("p2") + x[NADH]/nad_m);
    double activation = x[CA_m]/(params.at("p3") + x[CA_m]);
    return inibition * activation * jgpdh;
}



double  Bertram::J_o(const std::vector<double>& x, const std::map<std::string, double>& params){
    double term1 = (params.at("p4") * x[NADH]) /(params.at("p5") + x[NADH]);
    double term2 = 1 / (1 + exp(( x[DPSI] - params.at("p6"))/params.at("p7")));
    return term1 * term2 ;
}



double Bertram::J_Hres(const std::vector<double>& x, const std::map<std::string, double>& params){
    double term1 = (params.at("p8") * x[NADH]) /(params.at("p9") + x[NADH]);
    double term2 = 1 / (1 + exp(( x[DPSI] - params.at("p10"))/params.at("p11")));
    return term1 * term2;
}


double Bertram::J_Hatp(const std::vector<double>& x, const std::map<std::string, double>& params){
    double ATP_m = this->ATP_M(x, params);
    double term1 = params.at("p13") /(params.at("p13") + ATP_m);
    double term2 = params.at("p12") / (1 + exp((params.at("p14") - x[DPSI] )/params.at("p15")));
    return term1 * term2 ;
}



double Bertram::J_F1F0(const std::vector<double>& x, const std::map<std::string, double>& params){
    double ATP_m = this->ATP_M(x, params);
    double term1 = params.at("p13") /(params.at("p13") + ATP_m);
    double term2 = params.at("p16") / (1 + exp((params.at("p14")-x[DPSI] )/params.at("p15")));
    return term1 * term2;
}



double Bertram::J_Hleak(const std::vector<double>& x, const std::map<std::string, double>& params){
    double Jhleak= params.at("p17") * x[DPSI] + params.at("p18");
    if (Jhleak < 0){
        Jhleak = 0;
    }
    return Jhleak;
}



double Bertram::J_ANT(const std::vector<double>& x, const std::map<std::string, double>& params){
    double ATP_m = this->ATP_M(x, params);
    double adp_m = x[ADP];
    if (adp_m < 1e-16){
        adp_m = 1e-16;
    }
    double RAT_m = ATP_m / adp_m;
    double term1 = params.at("p19");
    double term2 = RAT_m / (params.at("p20") + RAT_m);
    double term3 = exp((.5 * x[DPSI] * params.at("FRT")));
    return term1 * term2 * term3;
}



double Bertram::J_uni(const std::vector<double>& x, const std::map<std::string, double>& params,double Ca_c){
    if (Ca_c < 0){
        Ca_c = 0;
    }
    double J_uni = ((params.at("p21") * x[DPSI]) - params.at("p22"))*Ca_c*Ca_c;
    if (J_uni < 0){
        J_uni = 0;
    }
    return J_uni;
}


double Bertram::J_NaCa(const std::vector<double>& x, const std::map<std::string, double>& params,double Ca_c){
    double ca_m = x[CA_m];
    if (ca_m < 0){
        ca_m = 0;
    }

    if (Ca_c < 1e-6){
        Ca_c = 1e-6;
    }

    double Jnaca = params.at("p23")*(ca_m / Ca_c)* exp(params.at("p24") * x[DPSI]);
    return Jnaca;
}


std::vector<double> Bertram::dxdt(const std::vector<double>& x,entree I, std::map<std::string,double> params, double t){
    params["FBP"] = I.FBP;
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
    
    std::vector<double> dx(taille);
    dx[NADH] =  ( jpdh - jo);
    dx[ADP] =  (jant - jf1f0);
    dx[DPSI] = (jhres - jhatp - jant - jhleak - jnaca - 2*juni)/params.at("Cm");
    dx[CA_m] = params.at("fm") * (juni - jnaca);
    return dx;
}