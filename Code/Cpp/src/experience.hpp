#pragma once
#include "formule.hpp"
#include "writer.hpp"
#include "solver.hpp"
#include <cmath>
#include <fstream>
#include<iostream>
#include <iomanip>
#include <map>
#include <vector>

class experience {
    public:
        experience() {}
        void variation_FBP( double dt,
                            double tfinal,
                            std::string nomfichier, 
                            std::map<std::string,double> params,
                            std::vector<double> vect_etat,
                            double FBP_0);

        void variation_Ca_c( double dt,
                            double tfinal,
                            std::string nomfichier, 
                            std::map<std::string,double> params,
                            std::vector<double> vect_etat,
                            double fbp_0,
                            double test);
    private:
        double ms_to_min(double t_ms);
        entree FBP_var(double t,std::map<std::string,double> params, double FBP_0);
        entree Ca_c_var(double t,std::map<std::string,double> params, double FBP_0);
        std::vector<double> dxdt_(
            const std::vector<double>& x,
            entree I,
            const std::map<std::string,double>& params,
            double t
        );
};