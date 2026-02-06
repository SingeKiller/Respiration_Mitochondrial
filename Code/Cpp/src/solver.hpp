#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <cmath>
#include "formule.hpp"
#include "writer.hpp"


// Initialiation du Calcium cytosolique en fonction du temps

double Ca_c(double t);


std::vector<double> add(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double> mult(const std::vector<double>& a, double k);

// Evolution de l'etat au cours du temps a chaque pas de temps t

std::vector<double> solver(
	std::vector<double> x,
    entree I,
	const std::map<std::string,double>& params,
	double t,
	double dt,
	std::vector<double> (*dxdt)(const std::vector<double>&, entree, const std::map<std::string,double>&, double)
);