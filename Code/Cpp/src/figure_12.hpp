#pragma once
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


entree FBP_var(double t, const ModelParams& params, double FBP_0);
entree Ca_c_var(double t, const ModelParams& params, double FBP_0);
void variation_FBP( double dt,
                    double tfinal,
                    std::string nomfichier, 
                    ModelParams params,
                    std::vector<double> vect_etat,
                    double FBP_0);
void variation_Ca_c( double dt,
                    double tfinal,
                    std::string nomfichier, 
                    ModelParams params,
                    std::vector<double> vect_etat
                    ,double fbp_0
                    ,double test);