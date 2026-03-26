#pragma once

#include "formule.hpp"
#include "solver.hpp"

#include <string>

class BE : public Bertram {
    public : 
        BE() : Bertram() {}
        
        // Reutilisation de dxdt pour le modèle étendu
        std::vector<double> dxdt(const std::vector<double>& x,
                                entree I, 
                                const ModelParams& params,
                                double t) override;

        // Flux ANT selon l'equation 35
        double J_ANT(const std::vector<double>& x, const ModelParams& params) override;
        
        // Somme cumulative des ajouts externes d'ADP
        double J_ADP_ext_sum(const ModelParams& params, double t);

        // Derivee temporelle de la somme, utilisee comme flux source dans d[ADP_c]/dt
        double J_ADP_ext_dot(const ModelParams& params, double t);
                                        
        // Calcul de [ATP]c
        double ATP_c(const std::vector<double>& x,const ModelParams& params);

        // Lance une simulation et ecrit Jo(t) dans un fichier
        void simulate_Jo_timecourse(std::vector<double> etat_initial,
                const ModelParams& params,
                entree I,
                double dt_ms,
                double tfinal_min,
            const std::string& output_file,
            const std::string& flux_output_file = "");

        // Variante avec 2 pulses ADP_c explicites
        void simulate_Jo_timecourse(std::vector<double> etat_initial,
                        const ModelParams& params,
                        entree I,
                        double dt_ms,
                        double tfinal_min,
                        double pulse1_time_min,
                        double pulse1_adp_c,
                        double pulse2_time_min,
                        double pulse2_adp_c,
                        const std::string& output_file,
                        const std::string& flux_output_file = "");
};
