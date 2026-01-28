#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <cmath>
#include "writer.hpp"


// Initialiation du Calcium cytosolique

double Ca_c(double t){
    double base = 0.0001;
    double pulse = 0.0006;
    if ((t>= 2 and t<= 2.5) or (t>= 3 and t<= 3.5) or (t>= 4 and t<= 4.5)){
        return pulse;
    }
    return base;
}

// Calcul de tous les flux

double J_GDH(std::map<std::string, double> params){
    if (params["FBP"] < 0){
        params["FBP"] = 0;
    }
    double J_GPDH = params["kGPDH"] * sqrt(params["FBP"]);
    return J_GPDH;
}



double J_PDH(etat x, std::map<std::string, double> params){
    double NAD_m = params["NAD_tot"] - x.NADH_m;
    if (NAD_m <= 0){
        NAD_m = 1e-16; //On eviter les valeurs négatif et nulles
    }

    double inibition = params["p1"]/(params["p1"]*x.NADH_m/NAD_m);
    double activation = x.Ca_m/(params["p3"] + x.Ca_m);
    return inibition * activation * params["JGPDH"];
}



double J_o(etat x, std::map<std::string, double> params){
    double term1 = (params["p4"] * x.NADH_m) /(params["p5"] + x.NADH_m);
    double term2 = 1 / (1 + exp(( x.deltaPsi - params["p6"])/params["p7"]));
    return term1 * term2 ;
}



double J_Hres(etat x, std::map<std::string, double> params){
    double term1 = (params["p8"] * x.NADH_m) /(params["p9"] + x.NADH_m);
    double term2 = 1 / (1 + exp(( x.deltaPsi - params["p10"])/params["p11"]));
    return term1 * term2;
}


double J_Hatp(etat x, std::map<std::string, double> params){
    double ATP_m = params["ADP_tot"] - x.ADP_m;

    if (ATP_m <= 0){
        ATP_m = 1e-16; //On eviter les valeurs négatif et nulles
    }

    double term1 = params["p13"] /(params["p13"] + ATP_m);
    double term2 = params["p12"] / (1 + exp((params["p14"] - x.deltaPsi )/params["p15"]));
    return term1 * term2 ;
}



double J_F1F0(etat x, std::map<std::string, double> params){
    double ATP_m = params["ADP_tot"] - x.ADP_m;

    if (ATP_m <= 0){
        ATP_m = 1e-16; //On eviter les valeurs négatif et nulles
    }

    double term1 = params["p13"] /(params["p13"] + ATP_m);
    double term2 = params["p16"] / (1 + exp((params["p14"] - x.deltaPsi )/params["p15"]));
    return term1 * term2 ;
}



double J_Hleak(etat x, std::map<std::string, double> params){
    return params["p17"] * x.deltaPsi + params["p18"];
}



double J_ANT(etat x, std::map<std::string, double> params){
    double ATP_m = params["ADP_tot"] - x.ADP_m;
    double RAT_m = ATP_m / x.ADP_m;
    
    if (RAT_m <= 0){
        RAT_m = 1e-16; //On eviter les valeurs négatif et nulles
    }
    double term1 = params["p19"];
    double term2 = RAT_m / (params["p20"] + RAT_m);
    double term3 = exp((.5 * x.deltaPsi * params["FRT"]));
    return term1 * term2 * term3;
}



double J_uni(etat x, std::map<std::string, double> params,double t){
    double ca_c = Ca_c(t);
    double J_uni = ((params["p21"] * x.deltaPsi) - params["p22"])*ca_c*ca_c;
    return J_uni;
}


double J_NaCa(etat x, std::map<std::string, double> params,double t){
    double ca_c = Ca_c(t);
    double Jnaca = params["p23"]*(x.Ca_m / ca_c)* exp(params["p24"] * x.deltaPsi);
    return Jnaca;
}
