#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <cmath>
#include "writer.hpp"


// Initialiation du Calcium cytosolique en fonction du temps

double Ca_c(double t);

// Equation differentiels ordinaires
double dNADH_dt(etat x, std::map<std::string, double> params,double t);
double dADP_m_dt(etat x, std::map<std::string, double> params, double t);
double ddeltapsi_dt(etat x, std::map<std::string, double> params, double t);
double dCa_m_dt(etat x, std::map<std::string,double> params, double t);

// Evolution de l'etat au cours du temps a chaque pas de temps t

etat detat_dt(etat x, std::map<std::string,double> params, double t);
etat solver(etat x, std::map<std::string,double> params, double t,double dt);