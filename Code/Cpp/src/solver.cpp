#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <cmath>
#include "writer.hpp"
#include "formule.hpp"


// Initialiation du Calcium cytosolique en fonction du temps

double Ca_c(double t){
    double base = 0.0001;
    double pulse = 0.0006;
    if ((t>= 2 and t<= 2.5) or (t>= 3 and t<= 3.5) or (t>= 4 and t<= 4.5)){
        return pulse;
    }
    return base;
}

//calcule des EDO differentes ( NADH(t) , ADP(t), deltaPsi(t), Ca_m(t) )

double dNADH_dt(etat x, std::map<std::string, double> params,double t){
    return params["g"] * ( J_PDH(x, params) - J_o(x, params));
}

double dADP_m_dt(etat x, std::map<std::string, double> params, double t){
    return params["g"]*(J_ANT(x,params)-J_F1F0(x,params));
}

double ddeltapsi_dt(etat x, std::map<std::string, double> params, double t){
    return (J_Hres(x,params)-J_Hatp(x,params)-J_ANT(x,params) + J_Hleak(x,params) - J_NaCa(x,params,t)-2*J_uni(x,params,t))/params["Cm"];
}

double dCa_m_dt(etat x, std::map<std::string,double> params, double t){
    return params["fm"]*(J_uni(x,params,t)-J_NaCa(x,params,t));
}


// modification de l'état au cours du temps 

etat detat_dt(etat x, std::map<std::string,double> params, double t){
    etat dx_dt;
    dx_dt.NADH_m = dNADH_dt(x,params,t);
    dx_dt.ADP_m = dADP_m_dt(x,params,t);
    dx_dt.deltaPsi = ddeltapsi_dt(x,params,t);
    dx_dt.Ca_m = dCa_m_dt(x,params,t);
    return dx_dt;
}

etat solver(etat x, std::map<std::string,double> params, double t,double dt){

    // calcul de la pente
    etat s1 = detat_dt(x,params,t);
    etat x2 = {x.NADH_m + s1.NADH_m*dt/2,
               x.ADP_m + s1.ADP_m*dt/2,
               x.deltaPsi + s1.deltaPsi*dt/2,
               x.Ca_m + s1.Ca_m*dt/2};


    // calcul du milieu
    etat s2 = detat_dt(x2,params,t + dt/2);
    etat x3 = {x.NADH_m + s2.NADH_m*dt/2,
               x.ADP_m + s2.ADP_m*dt/2,
               x.deltaPsi + s2.deltaPsi*dt/2,
               x.Ca_m + s2.Ca_m*dt/2};

    // calcul du milieu 
    etat s3 = detat_dt(x3,params,t + dt/2);
    etat x4 = {x.NADH_m + s3.NADH_m*dt,
               x.ADP_m + s3.ADP_m*dt,
               x.deltaPsi + s3.deltaPsi*dt,
               x.Ca_m + s3.Ca_m*dt};

    // calcul de la fin

    etat s4 = detat_dt(x4,params,t + dt);

    // mise a jour de l'état final
    x.NADH_m = x.NADH_m + (dt/6)*(s1.NADH_m + 2*s2.NADH_m + 2*s3.NADH_m + s4.NADH_m);
    x.ADP_m = x.ADP_m + (dt/6)*(s1.ADP_m + 2*s2.ADP_m + 2*s3.ADP_m + s4.ADP_m);
    x.deltaPsi = x.deltaPsi + (dt/6)*(s1.deltaPsi + 2*s2.deltaPsi + 2*s3.deltaPsi + s4.deltaPsi);
    x.Ca_m = x.Ca_m + (dt/6)*(s1.Ca_m + 2*s2.Ca_m + 2*s3.Ca_m + s4.Ca_m);

    return x;
}
