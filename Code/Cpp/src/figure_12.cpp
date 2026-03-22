#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>

#include "solver.hpp"
#include "formule.hpp"
#include "extend.hpp"
#include "teste.hpp"
#include "writer.hpp"
#include "namespace.hpp"


entree FBP_var(double t, const ModelParams& params, double FBP_0){
    entree I;
    double t_min = ms_to_min(t);
    if (t_min>= 1 and t_min<=6){
        I.FBP = FBP_0;
    }
    else{
        I.FBP = 1;
    }

    I.Ca_c = params.const_Ca_c;
    return I;
}


entree Ca_c_var(double t, const ModelParams& params, double FBP_0){
    entree I;
    I.FBP = FBP_0;
    double t_min = ms_to_min(t);
    double init = params.const_Ca_c;
    double amplitude = 0.5; // μM
    int pulse = ((t_min>= 2.0) and (t_min < 2.5)) or
                ((t_min>=3.0) and (t_min<3.5)) or
                ((t_min>=4.0) and (t_min<4.5));
    
    I.Ca_c = pulse ? amplitude : init;  // c'est un IF else synthetique
    return I;
}

void variation_FBP( double dt,
                    double tfinal,
                    std::string nomfichier, 
                    ModelParams params,
                    std::vector<double> vect_etat,
                    double FBP_0){

                        std::ofstream sortie(nomfichier.c_str());
                        sortie << "t" << "\t"<< "FBP" << "\t" << "Ca_c" << "\t" << "NADH_m" << "\t" << "Ca_m" << "\t" << "deltaPsi" << "\t" << "ATP_m"<<"\t"<< "J_o" << std::endl;
                        Bertram X;
                        double t = 0.;
                        double tfinal_ms = tfinal * 60000.; // conversion en ms
                        while (t<=tfinal_ms){
                            entree I = FBP_var(t,params,FBP_0);
                            double ATP_m = X.ATP_M(vect_etat,params);
                            double jo = X.J_o(vect_etat,params);

                            sortie << ms_to_min(t) << "\t" 
                            << I.FBP << "\t" 
                            << I.Ca_c << "\t" 
                            << vect_etat[NADH] << "\t" 
                            << vect_etat[CA_m] << "\t" 
                            << vect_etat[DPSI] << "\t" 
                            << ATP_m << "\t" 
                            << jo << std::endl;
                            
                            vect_etat = solver(vect_etat, I, params, t, dt, dxdt_);
                            t += dt;
                        }
                        sortie.close();
                    }
                
void variation_Ca_c( double dt,
                    double tfinal,
                    std::string nomfichier, 
                    ModelParams params,
                    std::vector<double> vect_etat
                    ,double fbp_0
                    ,double test){

                        std::ofstream sortie(nomfichier.c_str());
                        sortie << "t" << "\t"<< "FBP" << "\t" << "Ca_c" << "\t" << "NADH_m" << "\t" << "Ca_m" << "\t" << "deltaPsi" << "\t" << "ATP_m"<<"\t"<< "J_o" << std::endl;
                        Bertram X;
                        double t = 0.;
                        if (test == 1.){
                            params.p[3] = 0;
                        }
                        if (test == 2.){
                            params.p[21] = 0.02;
                        }
                        double tfinal_ms = tfinal * 60000.; // conversion en ms
                        while (t<=tfinal_ms){
                            entree I = Ca_c_var(t,params,fbp_0);
                            double ATP_m = X.ATP_M(vect_etat,params);
                            double jo = X.J_o(vect_etat,params);

                            sortie << ms_to_min(t) << "\t" 
                            << I.FBP << "\t" 
                            << I.Ca_c << "\t" 
                            << vect_etat[NADH] << "\t" 
                            << vect_etat[CA_m] << "\t" 
                            << vect_etat[DPSI] << "\t" 
                            << ATP_m << "\t" 
                            << jo << std::endl;
                            
                            vect_etat = solver(vect_etat, I, params, t, dt, dxdt_);
                            t += dt;
                        }
                        sortie.close();
                    }