#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <cmath>


struct etat{
    double NADH_m;
    double ADP_m;
    double deltaPsi;
    double Ca_m;
};


std::map<std::string, double> parameters(const std::string& nomfichier);