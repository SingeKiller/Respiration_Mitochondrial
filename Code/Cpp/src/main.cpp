#include "solver.hpp"
#include "formule.hpp"
#include "extend.hpp"
#include "teste.hpp"
#include "writer.hpp"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>


namespace {

double get_or(const std::map<std::string, double>& raw, const std::string& key, double fallback) {
    const auto it = raw.find(key);
    return (it != raw.end()) ? it->second : fallback;
}

void extra_params(ModelParams& params, const std::map<std::string, double>& raw) {
    params.gamma = get_or(raw, "gamma", params.gamma);
    params.Atot_c = get_or(raw, "Atot_c", params.Atot_c);
    params.ADP_c_0 = get_or(raw, "ADP_c_0", params.ADP_c_0);
    params.tau_plus = get_or(raw, "tau_plus", params.tau_plus);
    params.const_FBP = get_or(raw, "const_FBP", params.const_FBP);
    params.const_Ca_c = get_or(raw, "const_Ca_c", params.const_Ca_c);
    params.dt_ms = get_or(raw, "dt_ms", params.dt_ms);
    params.tfinal_min = get_or(raw, "tfinal_min", 25.0);
    params.VmaxANT = get_or(raw, "VmaxANT", params.p[19]);
    params.kANT_c_num = get_or(raw, "kANT_c_num", params.kANT_c_num);
    params.kANT_c_den = get_or(raw, "kANT_c_den", params.kANT_c_den);
    params.kANT_m_den = get_or(raw, "kANT_m_den", params.kANT_m_den);

    const int n_ajouts = static_cast<int>(get_or(raw, "n_ajouts", 0.0));
    params.adp_ajouts.clear();
    for (int i = 0; i < n_ajouts; ++i) {
        const std::string adp_key = "ADP_c_plus_" + std::to_string(i);
        const std::string t_key = "t_plus_" + std::to_string(i);
        const auto it_adp = raw.find(adp_key);
        const auto it_t = raw.find(t_key);
        if (it_adp != raw.end() and it_t != raw.end()) {
            params.adp_ajouts.push_back({it_adp->second, it_t->second});
        }
    }
}

void write_normalized_jo(const std::string& in_file, const std::string& out_file) {
    std::ifstream in(in_file.c_str());
    if (!in.is_open()) {
        throw std::runtime_error("Impossible d'ouvrir le fichier Jo: " + in_file);
    }

    std::string header;
    std::getline(in, header);

    std::vector<double> t_vals;
    std::vector<double> jo_vals;
    double t = 0.0;
    double jo = 0.0;
    while (in >> t >> jo) {
        t_vals.push_back(t);
        jo_vals.push_back(jo);
    }

    if (jo_vals.empty()) {
        throw std::runtime_error("Aucune donnee Jo trouvee dans: " + in_file);
    }

    double min_jo = jo_vals[0];
    double max_jo = jo_vals[0];
    for (double v : jo_vals) {
        if (v < min_jo) {
            min_jo = v;
        }
        if (v > max_jo) {
            max_jo = v;
        }
    }

    const double denom = max_jo - min_jo;
    std::ofstream out(out_file.c_str());
    if (!out.is_open()) {
        throw std::runtime_error("Impossible d'ecrire le fichier Jo normalise: " + out_file);
    }

    out << "t_min\tJ_o_norm" << std::endl;
    for (std::size_t i = 0; i < jo_vals.size(); ++i) {
        double jo_norm = 0.0;
        if (denom > 1e-16) {
            jo_norm = (jo_vals[i] - min_jo) / denom;
        }
        out << t_vals[i] << "\t" << jo_norm << std::endl;
    }
}

void write_o2_from_jo(const std::string& in_file, const std::string& out_file) {
    std::ifstream in(in_file.c_str());
    if (!in.is_open()) {
        throw std::runtime_error("Impossible d'ouvrir le fichier Jo pour integration O2: " + in_file);
    }

    std::string header;
    std::getline(in, header);

    std::vector<double> t_vals;
    std::vector<double> jo_vals;
    double t = 0.0;
    double jo = 0.0;
    while (in >> t >> jo) {
        t_vals.push_back(t);
        jo_vals.push_back(jo);
    }

    if (jo_vals.empty()) {
        throw std::runtime_error("Aucune donnee Jo trouvee pour integration O2: " + in_file);
    }

    std::vector<double> integ(jo_vals.size(), 0.0);
    for (std::size_t i = 1; i < jo_vals.size(); ++i) {
        const double dt_min = t_vals[i] - t_vals[i - 1];
        const double trap = 0.5 * (jo_vals[i] + jo_vals[i - 1]) * dt_min;
        integ[i] = integ[i - 1] + trap;
    }

   
    std::vector<double> o2_proxy(jo_vals.size(), 1.0);
    for (std::size_t i = 0; i < jo_vals.size(); ++i) {
        o2_proxy[i] = 1.0 - integ[i];
    }

    const double o2_start = o2_proxy.front();
    const double o2_end = o2_proxy.back();
    const double denom = o2_start - o2_end;

    std::ofstream out(out_file.c_str());
    if (!out.is_open()) {
        throw std::runtime_error("Impossible d'ecrire le fichier O2 simule: " + out_file);
    }

    out << "t_min\tO2_sim_norm" << std::endl;
    for (std::size_t i = 0; i < o2_proxy.size(); ++i) {
        double o2_norm = 1.0;
        if (std::abs(denom) > 1e-16) {
            o2_norm = (o2_proxy[i] - o2_end) / denom;
        }
        out << t_vals[i] << "\t" << o2_norm << std::endl;
    }
}

}


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



int main(int argc, char** argv){

    double dt = 1.0; // ms
    double tfinal = 7.0; // min

    if (argc == 2 and std::string(argv[1]) == "test"){
        std::filesystem::create_directories("resultats");
        std::filesystem::create_directories("plots");
        teste_solver_convergence();
        std::system("gnuplot plot_solver_convergence.gp");
        return 0;
    }

    constexpr int expected_argc = 32;
    if (argc != expected_argc) {
        std::cerr << "Usage:" << std::endl;
        std::cerr << "  ./main p1 p2 ... p24 FRT NADtot Atot Cm fm kGPDH output_filename" << std::endl;
        return 1;
    }

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
    const std::string output_filename = argv[idx++];

    const std::map<std::string, double> raw = parameters("parameters.txt");
    extra_params(params, raw);
    dt = params.dt_ms;
    tfinal = params.tfinal_min;

    std::filesystem::create_directories("resultats");

    
    BE X;

    std::vector<double> etat_be(taille_etendue, 0.0);
    etat_be[NADH] = 35.7411;
    etat_be[ADP] = 13351;
    etat_be[DPSI] = 152.271;
    etat_be[CA_m] = 0.0369;
    etat_be[ADP_c] = params.ADP_c_0;

    entree Ibe;
    Ibe.Ca_c = params.const_Ca_c;
    Ibe.FBP = params.const_FBP;

    

    const std::string o2_norm_path = std::string("resultats/") + output_filename;
    X.simulate_o2_normalized_integrated_with_pulses(
        etat_be,
        params,
        Ibe,
        dt,
        tfinal,
        2.5,
        330.0,
        10.0,
        1000.0,
        o2_norm_path);
    
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