#include "solver.hpp"
#include "formule.hpp"
#include "extend.hpp"
#include "teste.hpp"
#include "writer.hpp"
#include "namespace.hpp"
#include "figure_12.hpp"


#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>


int main(int argc, char** argv){

    double dt = 0.01;
    double tfinal = 25.0;

//Teste des arguments pour lancer les tests de convergence du solver

    if (argc == 2 and std::string(argv[1]) == "test"){
        std::filesystem::create_directories("resultats");
        std::filesystem::create_directories("plots");
        teste_solver_convergence();
        std::system("gnuplot plot_solver_convergence.gp");
        return 0;
    }

//Teste des arguments pour lancer le parametrage de l'individue 

    constexpr int expected_argc_legacy = 42;
    constexpr int expected_argc_extended = 44;
    if (argc != expected_argc_legacy and argc != expected_argc_extended) {
        std::cerr << "Usage:" << std::endl;
        std::cerr << "./main p1 p2 ... p24 FRT NADtot Atot Cm fm kGPDH ADP_c_0 Atot_c "
                     "const_FBP const_Ca_c dt_ms tfinal_min pulse1_time_min pulse1_adp_c "
                     "pulse2_time_min pulse2_adp_c [jo_norm_start_min o2_norm_start_min] output_filename" << std::endl;
        return 1;
    }

// Initialisation des paramètres à partir des arguments de la ligne de commande
    ModelParams params;

    int idx = 1;

    for (int p = 1; p <= 24; ++p) {
        params.p[p] = std::stod(argv[idx++]);
    }

    params.FRT = std::stod(argv[idx++]);
    params.NADtot = std::stod(argv[idx++]);
    params.Atot = std::stod(argv[idx++]);
    params.Cm = std::stod(argv[idx++]);
    params.fm = std::stod(argv[idx++]);
    params.kGPDH = std::stod(argv[idx++]);

    params.ADP_c_0 = std::stod(argv[idx++]);
    params.Atot_c = std::stod(argv[idx++]);
    params.const_FBP = std::stod(argv[idx++]);
    params.const_Ca_c = std::stod(argv[idx++]);
    params.dt_ms = std::stod(argv[idx++]);
    params.tfinal_min = std::stod(argv[idx++]);

    const double pulse1_time_min = std::stod(argv[idx++]);
    const double pulse1_adp_c = std::stod(argv[idx++]);
    const double pulse2_time_min = std::stod(argv[idx++]);
    const double pulse2_adp_c = std::stod(argv[idx++]);

    double jo_norm_start_min = 14.0;
    double o2_norm_start_min = 14.0;
    if (argc == expected_argc_extended) {
        jo_norm_start_min = std::stod(argv[idx++]);
        o2_norm_start_min = std::stod(argv[idx++]);
    }

    std::string output_filename = argv[idx++];

    params.ADP_c_plus_0 = pulse1_adp_c;
    params.ADP_c_plus_1 = pulse2_adp_c;
    params.adp_ajouts.clear();
    params.adp_ajouts.push_back({pulse1_adp_c, pulse1_time_min * 60000.0});
    params.adp_ajouts.push_back({pulse2_adp_c, pulse2_time_min * 60000.0});

    dt = params.dt_ms;
    tfinal = params.tfinal_min;

    std::filesystem::create_directories("resultats");
  
// Création d'une instancee BE
    BE X;

    std::vector<double> etat_be(taille_etendue, 0.0);
    etat_be[NADH] = 35.7411;
    etat_be[ADP] = 13351;
    etat_be[DPSI] = 152.271;
    etat_be[CA_m] = 0.0369;
    etat_be[ADP_c] = params.ADP_c_0;

// Ibe corresponds 
    entree Ibe;
    Ibe.Ca_c = params.const_Ca_c;
    Ibe.FBP = params.const_FBP;

    const std::string o2_norm_path = std::string("resultats/") + output_filename;
    const std::string jo_raw_path = "resultats/jo_temps.txt";
    const std::string jo_norm_path = "resultats/jo_normalized.txt";
    const std::string flux_raw_path = "resultats/fluxes_temps.txt";
    const std::string o2_from_jo_path = o2_norm_path;

    X.simulate_Jo_timecourse(
        etat_be,
        params,
        Ibe,
        dt,
        tfinal,
        pulse1_time_min, // temps du premier pulse en min
        pulse1_adp_c, // ADP_c du premier pulse en uM
        pulse2_time_min, // temps du second pulse en min
        pulse2_adp_c, // ADP_c du second pulse en uM
        jo_raw_path,
        flux_raw_path);
    
    write_normalized_jo(jo_raw_path, jo_norm_path, jo_norm_start_min);
    write_o2_from_jo(jo_raw_path, o2_from_jo_path, o2_norm_start_min);
    

    // variation_FBP(dt,
    //                 tfinal,
    //                 std::string("resultats/variation_FBP.txt"),
    //                 params,
    //                 std::vector<double>{35.7411, 13351.0, 152.271, 0.0369},
    //                 5.);
    // variation_FBP(dt,
    //                 tfinal,
    //                 std::string("resultats/variation_FBP_10.txt"),
    //                 params,
    //                 std::vector<double>{35.7411, 13351.0, 152.271, 0.0369},
    //                 10.);
    // variation_FBP(dt,
    //                 tfinal,
    //                 std::string("resultats/variation_FBP_15.txt"),
    //                 params,
    //                 std::vector<double>{35.7411, 13351.0, 152.271, 0.0369},
    //                 15.);

    // variation_Ca_c(dt,
    //                 tfinal,
    //                 std::string("resultats/variation_Ca_c.txt"),
    //                 params,
    //                 std::vector<double>{35.7411, 13351.0, 152.271, 0.0369},5.,
    //                 0.);

    // variation_Ca_c(dt,
    //                 tfinal,
    //                 std::string("resultats/variation_Ca_c_5.txt"),
    //                 params,
    //                 std::vector<double>{35.7411, 13351.0, 152.271, 0.0369},
    //                 5.,
    //                 1.);
    
    // variation_Ca_c(dt,
    //                 tfinal,
    //                 std::string("resultats/variation_Ca_c_10.txt"),
    //                 params,
    //                 std::vector<double>{35.7411, 13351.0, 152.271, 0.0369},
    //                 5.,
    //                 2.);

    return 0; 
}