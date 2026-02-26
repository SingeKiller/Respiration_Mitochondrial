#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <cmath>
#include "writer.hpp"

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

double ms_to_min(double t_ms){
    return t_ms / 60000.0;
}

std::vector<double> dxdt_(
    const std::vector<double>& x,
    entree I,
    const std::map<std::string,double>& params,
    double t);


// Calcul de tous les flux
class Bertram {
    public:
        
        Bertram(){}
        double NAD_m(const std::vector<double>& x, const std::map<std::string, double>& params);
        double ATP_M(const std::vector<double>& x, const std::map<std::string, double>& params);

        // flux independant du temps
        double J_GPDH(std::map<std::string, double> params);
        double J_PDH(const std::vector<double>& x, const std::map<std::string, double>& params, double jgpdh);
        double J_o(const std::vector<double>& x, const std::map<std::string, double>& params);
        double J_Hres(const std::vector<double>& x, const std::map<std::string, double>& params);
        double J_Hatp(const std::vector<double>& x, const std::map<std::string, double>& params);
        double J_F1F0(const std::vector<double>& x, const std::map<std::string, double>& params);
        double J_Hleak(const std::vector<double>& x, const std::map<std::string, double>& params);
        double J_ANT(const std::vector<double>& x, const std::map<std::string, double>& params);

        // flux en fonction du temps
        double J_uni(const std::vector<double>& x, const std::map<std::string, double>& params, double Ca_c);
        double J_NaCa(const std::vector<double>& x, const std::map<std::string, double>& params, double Ca_c);

        virtual std::vector<double> dxdt(const std::vector<double>& x,entree I, std::map<std::string,double> params, double t);
        
        void write(std::ofstream& sortie, double t, double FBP, double Ca_c, double NADH_m, double Ca_m, double DPSI, double ATP_m, double J_o);
    };

