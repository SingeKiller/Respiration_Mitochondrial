#pragma once

#include "formule.hpp"
#include "solver.hpp"
#include "writer.hpp"

// Structure pour les ajouts d'ADP
struct ADP_ajout {
    double ADP_c_plus;  // [ADP]+c,i : concentration ajoutée
    double t_plus;      // t+i : temps de l'ajout
};

class BE : public Bertram {
    public : 
        BE() : Bertram() {}
        
        // Reutilisation de dxdt pour le modèle étendu
        std::vector<double> dxdt(const std::vector<double>& x,
                                entree I, 
                                std::map<std::string,double> params, 
                                double t) override;
        
        // Fonction pour calculer J_ADP,ext
        double J_ADP_ext(const std::vector<double>& x,
                        const std::map<std::string,double>& params,
                        double t);
        
        // Calcul de [ATP]c
        double ATP_c(const std::vector<double>& x,const std::map<std::string,double>& params);
};
