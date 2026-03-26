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

namespace {

double get_or(const std::map<std::string, double>& raw, const std::string& key, double fallback) {
    const auto it = raw.find(key);
    return (it != raw.end()) ? it->second : fallback;
}

void extra_params(ModelParams& params, const std::map<std::string, double>& raw) {
    params.gamma = get_or(raw, "gamma", params.gamma);
    params.Atot_c = get_or(raw, "Atot_c", params.Atot_c);
    params.ADP_c_0 = get_or(raw, "ADP_c_0", params.ADP_c_0);
    params.ADP_c_plus_0 = get_or(raw, "ADP_c_plus_0", params.ADP_c_plus_0);
    params.ADP_c_plus_1 = get_or(raw, "ADP_c_plus_1", params.ADP_c_plus_1);
    params.tau_plus = get_or(raw, "tau_plus", params.tau_plus);
    params.const_FBP = get_or(raw, "const_FBP", params.const_FBP);
    params.const_Ca_c = get_or(raw, "const_Ca_c", params.const_Ca_c);
    params.dt_ms = get_or(raw, "dt_ms", params.dt_ms);
    params.tfinal_min = get_or(raw, "tfinal_min", 25.0);
    params.VmaxANT = get_or(raw, "VmaxANT", params.p[19]);
    params.kANT_c_num = get_or(raw, "kANT_c_num", params.kANT_c_num);
    params.kANT_c_den = get_or(raw, "kANT_c_den", params.kANT_c_den);
    params.kANT_m_den = get_or(raw, "kANT_m_den", params.kANT_m_den);

    const int n_ajouts = static_cast<int>(get_or(raw, "n_ajouts", 0.0));
    params.adp_ajouts.clear();
    for (int i = 0; i < n_ajouts; ++i) {
        const std::string adp_key = "ADP_c_plus_" + std::to_string(i);
        const std::string t_key = "t_plus_" + std::to_string(i);
        const auto it_adp = raw.find(adp_key);
        const auto it_t = raw.find(t_key);
        if (it_adp != raw.end() and it_t != raw.end()) {
            params.adp_ajouts.push_back({it_adp->second, it_t->second});
        }
    }
}

void write_normalized_jo(const std::string& in_file, const std::string& out_file) {
    std::ifstream in(in_file.c_str());
    if (!in.is_open()) {
        throw std::runtime_error("Impossible d'ouvrir le fichier Jo: " + in_file);
    }

    std::string header;
    std::getline(in, header);

    std::vector<double> t_vals;
    std::vector<double> jo_vals;
    double t = 0.0;
    double jo = 0.0;
    while (in >> t >> jo) {
        t_vals.push_back(t);
        jo_vals.push_back(jo);
    }

    if (jo_vals.empty()) {
        throw std::runtime_error("Aucune donnee Jo trouvee dans: " + in_file);
    }

    // Normalize using the analysis window from 14 to 20 min.
    const double norm_start_min = 14.0;
    const double norm_end_min = 20.0;
    bool found_window_point = false;
    double min_jo = 0.0;
    double max_jo = 0.0;
    for (std::size_t i = 0; i < jo_vals.size(); ++i) {
        if (t_vals[i] < norm_start_min || t_vals[i] > norm_end_min) {
            continue;
        }
        if (!found_window_point) {
            min_jo = jo_vals[i];
            max_jo = jo_vals[i];
            found_window_point = true;
        } else {
            if (jo_vals[i] < min_jo) {
                min_jo = jo_vals[i];
            }
            if (jo_vals[i] > max_jo) {
                max_jo = jo_vals[i];
            }
        }
    }

    // Fallback if no sample exists in the [14, 20] min window.
    if (!found_window_point) {
        min_jo = jo_vals[0];
        max_jo = jo_vals[0];
        for (double v : jo_vals) {
            if (v < min_jo) {
                min_jo = v;
            }
            if (v > max_jo) {
                max_jo = v;
            }
        }
    }

    const double denom = max_jo - min_jo;
    std::ofstream out(out_file.c_str());
    if (!out.is_open()) {
        throw std::runtime_error("Impossible d'ecrire le fichier Jo normalise: " + out_file);
    }

    out << "t_min\tJ_o_norm" << std::endl;
    for (std::size_t i = 0; i < jo_vals.size(); ++i) {
        double jo_norm = 0.0;
        if (denom > 1e-16) {
            jo_norm = (jo_vals[i] - min_jo) / denom;
        }
        out << t_vals[i] << "\t" << jo_norm << std::endl;
    }
}

void write_o2_from_jo(const std::string& in_file, const std::string& out_file) {
    std::ifstream in(in_file.c_str());
    if (!in.is_open()) {
        throw std::runtime_error("Impossible d'ouvrir le fichier Jo pour integration O2: " + in_file);
    }

    std::string header;
    std::getline(in, header);

    std::vector<double> t_vals;
    std::vector<double> jo_vals;
    double t = 0.0;
    double jo = 0.0;
    while (in >> t >> jo) {
        t_vals.push_back(t);
        jo_vals.push_back(jo);
    }

    if (jo_vals.empty()) {
        throw std::runtime_error("Aucune donnee Jo trouvee pour integration O2: " + in_file);
    }

    std::vector<double> integ(jo_vals.size(), 0.0);
    
    for (std::size_t i = 1; i < jo_vals.size(); ++i) {
        const double dt_min = t_vals[i] - t_vals[i - 1];
        const double trap = 0.5 * (jo_vals[i] + jo_vals[i - 1]) * dt_min;
        integ[i] = integ[i - 1] + trap;
    }

   
    std::vector<double> o2_proxy(jo_vals.size(), 1.0);

    for (std::size_t i = 0; i < jo_vals.size(); ++i) {
        o2_proxy[i] = 1.0 - integ[i];
    }

    std::ofstream out(out_file.c_str());
    if (!out.is_open()) {
        throw std::runtime_error("Impossible d'ecrire le fichier O2 simule: " + out_file);
    }

    const double output_start_min = 14.0;
    bool has_output_window = false;
    std::size_t window_start_idx = 0;
    for (std::size_t i = 0; i < t_vals.size(); ++i) {
        if (t_vals[i] >= output_start_min) {
            has_output_window = true;
            window_start_idx = i;
            break;
        }
    }

    const double o2_start = has_output_window ? o2_proxy[window_start_idx] : o2_proxy.front();
    const double o2_end = o2_proxy.back();
    const double denom = o2_start - o2_end;

    out << "t_min\tO2_sim_norm" << std::endl;
    for (std::size_t i = 0; i < o2_proxy.size(); ++i) {
        if (has_output_window && t_vals[i] < output_start_min) {
            continue;
        }
        double o2_norm = 1.0;
        if (std::abs(denom) > 1e-16) {
            o2_norm = (o2_proxy[i] - o2_end) / denom;
        }
        out << t_vals[i] << "\t" << o2_norm << std::endl;
    }
}

}