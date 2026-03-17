#pragma once
#include "formule.hpp"
#include "writer.hpp"
#include "solver.hpp"
#include <cmath>
#include <fstream>
#include<iostream>
#include <iomanip>
#include <vector>


std::vector<double> dxdt_exp(const std::vector<double>& x,entree /*I*/,const ModelParams& params,double /*t*/);
void teste_solver_convergence();