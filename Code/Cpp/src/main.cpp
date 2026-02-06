#include"solver.hpp"
#include"writer.hpp"
#include"formule.hpp"
#include <filesystem>


double ms_to_min(double t_ms){
    return t_ms / 60000.0;
}

entree FBP_var(double t,std::map<std::string,double> params, double FBP_0){
    entree I;
    double t_min = ms_to_min(t);
    if (t_min>= 1 and t_min<=6){
        I.FBP = FBP_0;
    }
    else{
        I.FBP = 1;
    }

    I.Ca_c = params.at("const_Ca_c"); 
    return I;
}


entree Ca_c_var(double t,std::map<std::string,double> params){
    entree I;
    I.FBP = params.at("const_FBP");
    double t_min = ms_to_min(t);
    double init = params.at("const_Ca_c");
    double amplitude = 0.5; // μM
    int pulse = ((t_min>= 2.0) and (t_min < 2.5)) or
                ((t_min>=3.0) and (t_min<3.5)) or
                ((t_min>=4.0) and (t_min<4.5));
    
    I.Ca_c = pulse ? amplitude : init;  // c'est un IF else synthetique
    return I;
}

std::vector<double> dxdt_(
    const std::vector<double>& x,
    entree I,
    const std::map<std::string,double>& params,
    double t
){
    Bertram model;
    return model.dxdt(x, I, params, t);
}


void variation_FBP( double dt,
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
                
void variation_Ca_c( double dt,
                    double tfinal,
                    std::string nomfichier, 
                    std::map<std::string,double> params,
                    std::vector<double> vect_etat){

                        std::ofstream sortie(nomfichier.c_str());
                        sortie << "t" << "\t"<< "FBP" << "\t" << "Ca_c" << "\t" << "NADH_m" << "\t" << "Ca_m" << "\t" << "deltaPsi" << "\t" << "ATP_m"<<"\t"<< "J_o" << std::endl;
                        Bertram X;
                        double t = 0.;
                        double tfinal_ms = tfinal * 60000.; // conversion en ms
                        while (t<=tfinal_ms){
                            entree I = Ca_c_var(t,params);
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

int main(int argc, char** argv){

    if (argc != taille + 2){
        std::cerr << "Usage: tfinal; NADH_m; ADP_m; Dpsi; Ca_m" << std::endl;
        return 1;
    }
    std::filesystem::create_directories("resultats");
    std::vector<double> vect_etat(taille);
    vect_etat[NADH] = std::stod(argv[2]);
    vect_etat[ADP] = std::stod(argv[3]);
    vect_etat[DPSI] = std::stod(argv[4]);
    vect_etat[CA_m] = std::stod(argv[5]);
    double tfinal = std::stod(argv[1]); // en min
    double dt = 1.0; // pas de temps en ms
    // paramètres du modèle

    std::map<std::string,double>params = parameters("parameters.txt");
    std::filesystem::create_directories("resultats");
    std::cout << params["p4"] << vect_etat[3] <<std::endl;  //test d'affichage
    //on génère un fichier d'ecriture de sortie


    std::filesystem::create_directories("resultats");
    variation_FBP(dt,
                    tfinal,
                    std::string("resultats/variation_FBP.txt"),
                    params,
                    vect_etat,
                    5.);
    variation_FBP(dt,
                    tfinal,
                    std::string("resultats/variation_FBP_10.txt"),
                    params,
                    vect_etat,
                    10.);
    variation_FBP(dt,
                    tfinal,
                    std::string("resultats/variation_FBP_15.txt"),
                    params,
                    vect_etat,
                    15.);

    variation_Ca_c(dt,
                    tfinal,
                    std::string("resultats/variation_Ca_c.txt"),
                    params,
                    vect_etat);
    return 0;  
}