#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>

constexpr std::size_t taille = 4; // taille modèle Bertram original
constexpr std::size_t taille_etendue = 5; // taille modèle étendu 
constexpr std::size_t NADH = 0;
constexpr std::size_t ADP = 1;
constexpr std::size_t DPSI = 2;
constexpr std::size_t CA_m = 3;
constexpr std::size_t ADP_c = 4;

// parametre de ;odification du milieu

struct entree{
    double Ca_c;
    double FBP;
};

struct ADP_ajout {
    double ADP_c_plus;
    double t_plus;
};

struct ModelParams {
    double p[25] = {0.0};
    double FRT = 0.0;
    double gamma = 0.001;
    double NADtot = 0.0;
    double Atot = 0.0;
    double Cm = 1.8;
    double fm = 0.01;
    double kGPDH = 0.0;

    // Parametres ANT (equation 35)
    double VmaxANT = 0.35;
    double kANT_c_num = 0.8;
    double kANT_c_den = 0.11;
    double kANT_m_den = 7.2;
    double Atot_c = 15000.0;
    double ADP_c_0 = 0.0;
    double tau_plus = 1000.0;
    
    std::vector<ADP_ajout> adp_ajouts;
    double const_FBP = 5.0;
    double const_Ca_c = 0.1;
    double dt_ms = 1.0;
    double tfinal_min = 10.0;
    double pre_equil_min = 2.0;
    double k = 0.0;
};

inline double ms_to_min(double t_ms){
    return t_ms / 60000.0;
}

std::vector<double> dxdt_(
    const std::vector<double>& x,
    entree I,
    const ModelParams& params,
    double t);


// Calcul de tous les flux
class Bertram {
    public:
        
        Bertram(){}
        double NAD_m(const std::vector<double>& x, const ModelParams& params);
        double ATP_M(const std::vector<double>& x, const ModelParams& params);

        // flux independant du temps
        double J_GPDH(double FBP, const ModelParams& params);
        double J_PDH(const std::vector<double>& x, const ModelParams& params, double jgpdh);
        double J_o(const std::vector<double>& x, const ModelParams& params);
        double J_Hres(const std::vector<double>& x, const ModelParams& params);
        double J_Hatp(const std::vector<double>& x, const ModelParams& params);
        double J_F1F0(const std::vector<double>& x, const ModelParams& params);
        double J_Hleak(const std::vector<double>& x, const ModelParams& params);
        virtual double J_ANT(const std::vector<double>& x, const ModelParams& params);

        // flux en fonction du temps
        double J_uni(const std::vector<double>& x, const ModelParams& params, double Ca_c);
        double J_NaCa(const std::vector<double>& x, const ModelParams& params, double Ca_c);

        virtual std::vector<double> dxdt(const std::vector<double>& x,entree I, const ModelParams& params, double t);
        
        void write(std::ofstream& sortie, double t, double FBP, double Ca_c, double NADH_m, double Ca_m, double DPSI, double ATP_m, double J_o);

    protected:
        double nadh_value() const { return nadh_; }
        double adp_value() const { return adp_; }
        double dpsi_value() const { return dpsi_; }
        double ca_m_value() const { return ca_m_; }
        void sync_state(const std::vector<double>& x);

    private:
        double nadh_ = 35.7411;
        double adp_ = 13351;
        double dpsi_ = 152.271;
        double ca_m_ = 0.0369789;
    };

