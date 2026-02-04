#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <cmath>
#include "solver.hpp"
#include "writer.hpp"
#include "formule.hpp"


std::vector<double> add(const std::vector<double>& a, const std::vector<double>& b){
    std::vector<double> result(a.size());
    for (std::size_t i = 0; i < a.size(); ++i){
        result[i] = a[i] + b[i];
    }
    return result;
}

std::vector<double> mult(const std::vector<double>& a, double k){
    std::vector<double> result(a.size());
    for (std::size_t i = 0; i < a.size(); ++i){
        result[i] = a[i] * k;
    }
    return result;
}


//;odification pour entree une fomntion en parametre , revoir la structure 
std::vector<double> solver(
                        std::vector<double> x,
                        entree I,
                        const std::map<std::string,double>& params,
                        double t,
                        double dt,
                        std::vector<double> (*dxdt)(const std::vector<double>&, entree, const std::map<std::string,double>&, double)){

    // calcul de la pente
    std::vector<double> s1 = dxdt(x,I,params,t);
    std::vector<double> s2 = dxdt(add(x, mult(s1, dt/2)),I,params,t + dt/2);
    std::vector<double> s3 = dxdt(add(x, mult(s2, dt/2)),I,params,t + dt/2);
    std::vector<double> s4 = dxdt(add(x, mult(s3, dt)),I,params,t + dt);

    // mise a jour de l'état final
    std::vector<double> xn = add(x, mult(add(s1,add(mult(s2,2), add(mult(s3,2), s4))), dt/6));

    return xn;
}
