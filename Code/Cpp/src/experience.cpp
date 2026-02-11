#include "formule.hpp"
#include "writer.hpp"
#include "solver.hpp"
#include "experience.hpp"
#include <cmath>
#include <fstream>
#include<iostream>
#include <iomanip>
#include <map>
#include <vector>

experience::experience() {};
void experience::variation_FBP( double dt,
                    double tfinal,
                    std::string nomfichier, 
                    std::map<std::string,double> params,
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


void experience::variation_Ca_c( double dt,
                    double tfinal,
                    std::string nomfichier, 
                    std::map<std::string,double> params,
                    std::vector<double> vect_etat
                    ,double fbp_0
                    ,double test){
                        std::ofstream sortie(nomfichier.c_str());
                        sortie << "t" << "\t"<< "FBP" << "\t" << "Ca_c" << "\t" << "NADH_m" << "\t" << "Ca_m" << "\t" << "deltaPsi" << "\t" << "ATP_m"<<"\t"<< "J_o" << std::endl;
                        Bertram X;
                        double t = 0.;
                        if (test == 1.){
                            params["p3"]= 0;
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
experience::~experience() {};
