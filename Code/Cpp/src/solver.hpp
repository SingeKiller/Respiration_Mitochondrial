#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <cmath>
#include <functional>
#include "formule.hpp"


// Initialiation du Calcium cytosolique en fonction du temps

double Ca_c(double t);


std::vector<double> add(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double> mult(const std::vector<double>& a, double k);

// Evolution de l'etat au cours du temps a chaque pas de temps t

std::vector<double> solver(
	std::vector<double> x,
    entree I,
	const ModelParams& params,
	double t,
	double dt,
	const std::function<std::vector<double>(const std::vector<double>&, entree, const ModelParams&, double)>& dxdt
);