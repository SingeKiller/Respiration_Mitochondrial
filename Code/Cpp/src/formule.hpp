#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <cmath>
#include "writer.hpp"

// Calcul de tous les flux

double J_GDH(std::map<std::string, double> params);
double J_PDH(etat x, std::map<std::string, double> params);
double J_o(etat x, std::map<std::string, double> params);
double J_Hres(etat x, std::map<std::string, double> params);
double J_Hatp(etat x, std::map<std::string, double> params);
double J_F1F0(etat x, std::map<std::string, double> params);
double J_Hleak(etat x, std::map<std::string, double> params);
double J_ANT(etat x, std::map<std::string, double> params);

// flux en fonction du temps
double J_uni(etat x, std::map<std::string, double> params, double t);
double J_NaCa(etat x, std::map<std::string, double> params, double t);