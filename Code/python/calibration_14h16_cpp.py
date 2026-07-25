#!/usr/bin/env python3

from pathlib import Path
import subprocess
import time
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

CPP_ARGS = [
    "400", "1", "0.01", "0.6", "100", "177", "5", "7", "100", "177",
    "5", "120", "10000", "190", "8.5", "35", "0.002", "-0.03", "0.35",
    "2", "0.01", "1.1", "0.001", "0.016", "0.037", "10000", "15000",
    "1.8", "0.01", "0.0005", "0", "0.1", "15000", "0.1", "1000", "21",
    "16.15", "330", "19.15", "1000", "o2_test_normalized.txt",
]

NOMS_PARAMETRES_CALIBRES = [f"P{i}" for i in range(1, 25)] + ["FRT", "NADtot", "Atot","Cm", "fm", "kGPDH"]
INDICES_PARAMETRES_SUPPLEMENTAIRES = {
    "FRT": 24,
    "NADtot": 25,
    "Atot": 26,
    "Cm": 27,
    "fm": 28,    
    "kGPDH": 29,
}

JEUX_CALIBRATION = {
    #"14h16": {
    #    "fichier": "data14h16.npy",
    #    "fenetres": [[14.6, 15.5], [16.25, 16.8], [17.6, 18.5], [19.7, 20.5]],
    #    "fenetre_data": [14.0, 21.0],
    #    "debut_jo_norm": 14.0,
    #    "debut_o2_norm": 14.0,
    #    "ajouts_adp_mM": [0.33, 1.0],
    #    "pulse_times_min": [16.1981, 18.9969],
    #},
    "15h12": {
        "fichier": "data15h12.npy",
        "fenetres": [[10.8, 12.85], [13.4, 15.0], [16., 17.0], [17.4, 17.8]],
        "fenetre_data": [10.0, 17.9],
        "debut_jo_norm": 10.8,
        "debut_o2_norm": 10.8,
        "ajouts_adp_mM": [0.66, 1.0],
        "pulse_times_min": [13.1431, 17.2056],
    },
    "15h37": {
        "fichier": "data15h37.npy",
        "fenetres": [[10.6, 11.35], [12.25, 12.65], [13.05, 13.6], [14.25, 15.7]],
        "fenetre_data": [10.0, 15.7],
        "debut_jo_norm": 10.6,
        "debut_o2_norm": 10.6,
        "ajouts_adp_mM": [0.17, 1.0],
        "pulse_times_min": [11.9725, 14.0256],
    },
    "16h00": {
        "fichier": "data16h00.npy",
        "fenetres": [[10.75, 11.8], [12.45, 13.8], [14.2, 15.05], [15.4, 16.3]],
        "fenetre_data": [10.0, 16.3],
        "debut_jo_norm": 10.75,
        "debut_o2_norm": 10.75,
        "ajouts_adp_mM": [0.66, 1.0],
        "pulse_times_min": [12.1288, 15.2069],
    },
    "17h03": {
        "fichier": "data17h03.npy",
        "fenetres": [[10.3, 11.5], [12.15, 12.7], [13.7, 14.85], [15.35, 17.0]],
        "fenetre_data": [10.0, 18.0],
        "debut_jo_norm": 10.3,
        "debut_o2_norm": 10.3,
        "ajouts_adp_mM": [0.17, 1.0],
        "pulse_times_min": [11.8156, 15.0250],
    },
}
OUTPUT_NAME = "o2_test_normalized.txt"

class ConfigurationCalibration:
    def __init__(
        self,
        population_size=20,
        generations=5,
        elite_size=25,
        init_low=0.70,
        init_high=1.30,
        mutation_low=0.97,
        mutation_high=1.03,
        poids_ratios=1.0,
        poids_points_l2=0.2,
        random_seed=42,
        checkpoint_every=5,
        cpp_timeout_sec=600,
    ):
        self.population_size = population_size
        self.generations = generations
        self.elite_size = elite_size
        self.init_low = init_low
        self.init_high = init_high
        self.mutation_low = mutation_low
        self.mutation_high = mutation_high
        self.poids_ratios = poids_ratios
        self.poids_points_l2 = poids_points_l2
        self.random_seed = random_seed
        self.checkpoint_every = checkpoint_every
        self.cpp_timeout_sec = cpp_timeout_sec


def seed_generation(base_seed, generation):
    return int(base_seed) + int(generation)

def resoudre_chemins(fichier_experimental):
    racine = Path(__file__).resolve().parents[1] 
    dossier_cpp = racine / "Cpp"
    dossier_py = racine / "python"
    dossier_sortie = dossier_py / "plots"
    dossier_sortie.mkdir(parents=True, exist_ok=True)
    executable = dossier_cpp / "main.exe"
    
    if not executable.exists():
        executable = dossier_cpp / "main"
        
    if not executable.exists():
        raise FileNotFoundError("Executable C++ introuvable (main.exe ou main).")
    
    chemin_experimental = dossier_py / "normalized_respiration_data" / fichier_experimental
    
    if not chemin_experimental.exists():
        raise FileNotFoundError(f"Fichier experimental introuvable: {chemin_experimental}")

    chemin_modele = dossier_cpp / "resultats" / OUTPUT_NAME
    return {"racine": racine,"dossier_cpp": dossier_cpp,"dossier_py": dossier_py,"dossier_sortie": dossier_sortie,"executable": executable,"experimental": chemin_experimental,"sortie_modele": chemin_modele,
    }

def charger_deux_colonnes(chemin):
    try:
        tableau = np.loadtxt(chemin, skiprows=1)
    except ValueError:
        tableau = np.loadtxt(chemin, skiprows=0)

    if tableau.ndim == 1:
        tableau = tableau.reshape(1, -1)
    if tableau.shape[1] < 2:
        raise ValueError(f"Le fichier {chemin} doit avoir au moins 2 colonnes")

    return tableau[:, 0], tableau[:, 1]


def extraire_fenetre_temps_depuis_npy(dossier_py, fichier_experimental):
    chemin = dossier_py / "normalized_respiration_data" / fichier_experimental
    if not chemin.exists():
        raise FileNotFoundError(f"Fichier experimental introuvable: {chemin}")

    courbe = np.load(chemin)
    if courbe.ndim != 2 or courbe.shape[0] < 1 or courbe.shape[1] < 1:
        raise ValueError(f"Format invalide pour {chemin}: tableau attendu [2, N]")

    temps = courbe[0, :]
    t_min = float(np.min(temps))
    t_max = float(np.max(temps))
    if t_min >= t_max:
        raise ValueError(f"Fenetre temps invalide dans {chemin}: min={t_min}, max={t_max}")

    return [t_min, t_max]


def extraire_pentes_segments(temps, signal, fenetres):
    pentes = []
    for t_min, t_max in fenetres:
        masque = (temps >= t_min) & (temps <= t_max)
        x_seg = temps[masque]
        y_seg = signal[masque]
        if x_seg.size < 3:
            raise ValueError(f"Segment [{t_min}, {t_max}] insuffisant ({x_seg.size} points)")
        regression = stats.linregress(x_seg, y_seg)
        pentes.append(float(regression.slope))
    return np.array(pentes, dtype=float)


def valider_couverture_fenetres(temps, fenetres):
    t_min = float(np.min(temps))
    t_max = float(np.max(temps))
    for debut, fin in fenetres:
        if debut < t_min or fin > t_max:
            raise ValueError(
                f"Couverture temporelle sim insuffisante: sim=[{t_min:.4g}, {t_max:.4g}] "
                f"ne couvre pas segment [{debut}, {fin}]"
            )


def calculer_ratios_pentes(pentes):

    epsilon = 1e-12
    p1, p2, p3, p4 = pentes
    ratios = np.array([
        p1 / (p2 + epsilon),
        p1 / (p3 + epsilon),
        p1 / (p4 + epsilon),
        p2 / (p3 + epsilon),
        p2 / (p4 + epsilon),
        p3 / (p4 + epsilon),
    ], dtype=float)
    return ratios

def calculer_cout_quadratique(pentes_exp, pentes_sim, ratios_exp, ratios_sim, poids_ratios):
    ecart_pentes = pentes_exp - pentes_sim
    ecart_ratios = ratios_exp - ratios_sim
    cout_pentes = float(np.dot(ecart_pentes, ecart_pentes))
    cout_ratios = float(np.dot(ecart_ratios, ecart_ratios))
    return cout_pentes + poids_ratios * cout_ratios


def calculer_norme_l2_points(temps_exp, o2_exp, temps_sim, o2_sim):
    debut_commun = max(float(np.min(temps_exp)), float(np.min(temps_sim)))
    fin_commun = min(float(np.max(temps_exp)), float(np.max(temps_sim)))
    if debut_commun >= fin_commun:
        return float("inf")

    masque_exp = (temps_exp >= debut_commun) & (temps_exp <= fin_commun)
    temps_exp_commun = temps_exp[masque_exp]
    o2_exp_commun = o2_exp[masque_exp]

    if temps_exp_commun.size < 2:
        return float("inf")

    o2_sim_interp = np.interp(temps_exp_commun, temps_sim, o2_sim)
    diff = o2_exp_commun - o2_sim_interp
    return float(np.linalg.norm(diff, ord=2))


def construire_vecteur_base_calibrable():
    base = [float(CPP_ARGS[i]) for i in range(24)]
    
    for nom in ["FRT", "NADtot", "Atot", "Cm", "fm", "kGPDH"]:
        base.append(CPP_ARGS[INDICES_PARAMETRES_SUPPLEMENTAIRES[nom]])
        
    return np.array(base, dtype=float)


def ajouts_adp_um_depuis_mm(ajouts_adp_mm):
    return ajouts_adp_mm[0] * 1000.0, ajouts_adp_mm[1] * 1000.0


def preparer_arguments_cpp(
    individu,
    debut_jo_norm,
    debut_o2_norm,
    ajouts_adp_mm,
    tfinal_min=None,
    pulse_times_min=None,
):
    if individu.size != len(NOMS_PARAMETRES_CALIBRES):
        raise ValueError(
            f"Taille individu invalide: {individu.size}, attendu {len(NOMS_PARAMETRES_CALIBRES)}"
        )

    args = CPP_ARGS.copy()

    for i in range(24):
        args[i] = f"{float(individu[i]):.12g}"

    indices_individu = {nom: idx for idx, nom in enumerate(NOMS_PARAMETRES_CALIBRES)}
    for nom, idx_cpp in INDICES_PARAMETRES_SUPPLEMENTAIRES.items():
        args[idx_cpp] = f"{float(individu[indices_individu[nom]]):.12g}"

    ajout1_um, ajout2_um = ajouts_adp_um_depuis_mm(ajouts_adp_mm)
    args[37] = f"{ajout1_um:.12g}"
    args[39] = f"{ajout2_um:.12g}"

    if pulse_times_min is not None:
        if len(pulse_times_min) != 2:
            raise ValueError("pulse_times_min doit contenir exactement 2 temps (min).")
        args[36] = f"{float(pulse_times_min[0]):.12g}"
        args[38] = f"{float(pulse_times_min[1]):.12g}"

    if tfinal_min is not None:
        # C++ argument order: ... dt_ms (idx 34), tfinal_min (idx 35), pulses...
        args[35] = f"{float(tfinal_min):.12g}"

    args_sans_sortie = args[:40]
    fichier_sortie = args[40]
    return args_sans_sortie + [f"{debut_jo_norm:.12g}", f"{debut_o2_norm:.12g}", fichier_sortie]


def executer_modele_cpp(executable, dossier_cpp, arguments, timeout_sec):
    subprocess.run(
        [str(executable), *arguments],
        cwd=str(dossier_cpp),
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        timeout=float(timeout_sec),
    )


def population_initiale(aleatoire, cfg):
    base = construire_vecteur_base_calibrable()
    echelles = echantillonner_facteurs_normaux_tronques(
        aleatoire,
        cfg.init_low,
        cfg.init_high,
        size=(cfg.population_size, base.size),
    )
    return base * echelles


def echantillonner_facteurs_normaux_tronques(aleatoire, low, high, size):
    if not (low < 1.0 < high):
        raise ValueError(f"Bornes invalides pour une loi normale centree en 1: low={low}, high={high}")

    sigma = min(1.0 - low, high - 1.0) / 3.0
    a = (low - 1.0) / sigma
    b = (high - 1.0) / sigma
    return stats.truncnorm.rvs(a, b, loc=1.0, scale=sigma, size=size, random_state=aleatoire)


def faire_evoluer_population(aleatoire, parents, cfg):
    enfants = []
    nb_enfants = max(0, cfg.population_size - len(parents))
    while len(enfants) < nb_enfants:
        i, j = aleatoire.integers(0, len(parents), size=2)
        enfant = (parents[i] + parents[j]) / 2.0
        mutation = echantillonner_facteurs_normaux_tronques(
            aleatoire,
            cfg.mutation_low,
            cfg.mutation_high,
            size=parents.shape[1],
        )
        enfant = enfant * mutation
        enfants.append(enfant)
    if len(enfants) == 0:
        return np.empty((0, parents.shape[1]), dtype=float)
    return np.array(enfants, dtype=float)

 
def borner_intervalle_mutation(low, high):
    low_borne = max(1e-6, float(low))
    high_borne = min(3.0, float(high))
    if not (low_borne < 1.0 < high_borne):
        low_borne = min(low_borne, 0.999)
        high_borne = max(high_borne, 1.001)
    return low_borne, high_borne


def lancer_calibration(
    cfg,
    nom_jeu,
    fichier_experimental,
    fenetres,
    debut_jo_norm,
    debut_o2_norm,
    ajouts_adp_mm,
    tfinal_min,
    pulse_times_min,
):
    chemins = resoudre_chemins(fichier_experimental)

    courbe_exp = np.load(chemins["experimental"])
    temps_exp = courbe_exp[0, :]
    o2_exp = courbe_exp[1, :]
    pentes_exp = extraire_pentes_segments(temps_exp, o2_exp, fenetres)
    ratios_exp = calculer_ratios_pentes(pentes_exp)

    seed_init = seed_generation(cfg.random_seed, 0)
    aleatoire_init = np.random.default_rng(seed_init)
    population = population_initiale(aleatoire_init, cfg)
    historique = []
    stagnation_count = 0
    mutation_low_actuelle = float(cfg.mutation_low)
    mutation_high_actuelle = float(cfg.mutation_high)

    meilleur_global = None

    for gen in range(cfg.generations):
        scores = np.full(cfg.population_size, np.inf, dtype=float)
        pentes_sim_stock = [None] * cfg.population_size
        ratios_sim_stock = [None] * cfg.population_size
        print(
            f"Generation {gen:02d} | evaluation en cours ({cfg.population_size} individus)",
            flush=True,
        )

        for idx, individu in enumerate(population):
            print(
                f"Generation {gen:02d} | individu {idx + 1:02d}/{cfg.population_size:02d}",
                flush=True,
            )
            try:
                executer_modele_cpp(
                    chemins["executable"],
                    chemins["dossier_cpp"],
                    preparer_arguments_cpp(
                        individu,
                        debut_jo_norm,
                        debut_o2_norm,
                        ajouts_adp_mm,
                        tfinal_min=tfinal_min,
                        pulse_times_min=pulse_times_min,
                    ),
                    cfg.cpp_timeout_sec,
                )
                temps_sim, o2_sim = charger_deux_colonnes(chemins["sortie_modele"])
                valider_couverture_fenetres(temps_sim, fenetres)
                pentes_sim = extraire_pentes_segments(temps_sim, o2_sim, fenetres)
                ratios_sim = calculer_ratios_pentes(pentes_sim)
                cout_quadratique = calculer_cout_quadratique(
                    pentes_exp,
                    pentes_sim,
                    ratios_exp,
                    ratios_sim,
                    cfg.poids_ratios,
                )
                norme_l2_points = calculer_norme_l2_points(temps_exp, o2_exp, temps_sim, o2_sim)
                score = cout_quadratique + cfg.poids_points_l2 * norme_l2_points
                scores[idx] = score
                pentes_sim_stock[idx] = pentes_sim
                ratios_sim_stock[idx] = ratios_sim
            except subprocess.TimeoutExpired:
                print(
                    f"Generation {gen:02d} | individu {idx + 1:02d} timeout > {cfg.cpp_timeout_sec}s, penalise",
                    flush=True,
                )
                scores[idx] = 1e9
                pentes_sim_stock[idx] = np.full_like(pentes_exp, np.nan)
                ratios_sim_stock[idx] = np.full_like(ratios_exp, np.nan)
            except Exception as exc:
                print(
                    f"Generation {gen:02d} | individu {idx + 1:02d} erreur: {type(exc).__name__}: {exc}",
                    flush=True,
                )
                scores[idx] = 1e9
                pentes_sim_stock[idx] = np.full_like(pentes_exp, np.nan)
                ratios_sim_stock[idx] = np.full_like(ratios_exp, np.nan)

        ordre = np.argsort(scores)
        meilleur_idx = int(ordre[0])
        meilleur_score = float(scores[meilleur_idx])
        meilleur_individu = population[meilleur_idx].copy()
        meilleures_pentes_sim = pentes_sim_stock[meilleur_idx]
        meilleurs_ratios_sim = ratios_sim_stock[meilleur_idx]

        if meilleur_global is None or meilleur_score < meilleur_global["fitness"]:
            meilleur_global = {
                "generation": gen,
                "fitness": meilleur_score,
                "individual": meilleur_individu,
                "sim_slopes": np.array(meilleures_pentes_sim, dtype=float),
                "exp_slopes": pentes_exp.copy(),
                "sim_ratios": np.array(meilleurs_ratios_sim, dtype=float),
                "exp_ratios": ratios_exp.copy(),
            }
            stagnation_count = 0
        else:
            stagnation_count += 1

        historique.append(
            {
                "generation": gen,
                "best_fitness": meilleur_score,
                "seed_generation": seed_generation(cfg.random_seed, gen),
            }
        )

        print(
            f"Generation {gen:02d} | best={meilleur_score:.6e} | seed={seed_generation(cfg.random_seed, gen)}",
            flush=True,
        )

        if (gen + 1) % cfg.checkpoint_every == 0:
            ecrire_checkpoint_intermediaire(
                chemins["dossier_sortie"],
                nom_jeu,
                historique,
                meilleur_global,
                gen,
            )

        if cfg.elite_size <= 1:
            nb_parents = int(round(cfg.population_size * float(cfg.elite_size)))
        else:
            nb_parents = int(round(cfg.population_size * float(cfg.elite_size) / 100.0))
        nb_parents = max(2, min(cfg.population_size - 1, nb_parents))
        parents = population[ordre[:nb_parents]]
        seed_enfants = seed_generation(cfg.random_seed, gen + 1)
        aleatoire_enfants = np.random.default_rng(seed_enfants)
        if stagnation_count > 3:
            mutation_low_actuelle -= 0.02
            mutation_high_actuelle += 0.02
            mutation_low_actuelle, mutation_high_actuelle = borner_intervalle_mutation(
                mutation_low_actuelle,
                mutation_high_actuelle,
            )
            stagnation_count = 0
            print(
                f"Generation {gen:02d} | stagnation detectee, mutation elargie a [{mutation_low_actuelle:.4f}, {mutation_high_actuelle:.4f}]",
                flush=True,
            )

        mutation_low_originale = cfg.mutation_low
        mutation_high_originale = cfg.mutation_high
        cfg.mutation_low = mutation_low_actuelle
        cfg.mutation_high = mutation_high_actuelle
        enfants = faire_evoluer_population(aleatoire_enfants, parents, cfg)
        cfg.mutation_low = mutation_low_originale
        cfg.mutation_high = mutation_high_originale
        population = np.vstack([parents, enfants])

    assert meilleur_global is not None
    return {
        "paths": chemins,
        "best": meilleur_global,
        "config": cfg,
        "nom_jeu": nom_jeu,
        "fichier_experimental": fichier_experimental,
        "fenetres": fenetres,
        "debut_jo_norm": debut_jo_norm,
        "debut_o2_norm": debut_o2_norm,
        "ajouts_adp_mm": ajouts_adp_mm,
        "tfinal_min": tfinal_min,
        "pulse_times_min": pulse_times_min,
    }, historique


def ecrire_sorties(resultat, historique):
    dossier_sortie = resultat["paths"]["dossier_sortie"]
    best = resultat["best"]
    cfg = resultat["config"]

    nom_jeu = resultat["nom_jeu"]
    chemin_historique = dossier_sortie / f"calibration_{nom_jeu}_generation.csv"
    chemin_parametres = dossier_sortie / f"calibration_{nom_jeu}_parametre.csv"

    en_tete = "generation,best_fitness\n"
    lignes = [
        f"{h['generation']},{h['best_fitness']}\n"
        for h in historique
    ]
    with chemin_historique.open("w", encoding="utf-8") as f:
        f.write(en_tete)
        f.writelines(lignes)

    en_tete_parametres = NOMS_PARAMETRES_CALIBRES + ["best_fitness", "best_generation"]
    valeurs_parametres = [str(float(best["individual"][i])) for i in range(len(NOMS_PARAMETRES_CALIBRES))]
    valeurs_parametres.extend([str(float(best["fitness"])), str(int(best["generation"]))])
    with chemin_parametres.open("w", encoding="utf-8") as f:
        f.write(",".join(en_tete_parametres) + "\n")
        f.write(",".join(valeurs_parametres) + "\n")

    return chemin_historique, chemin_parametres


def ecrire_checkpoint_intermediaire(dossier_sortie, nom_jeu, historique, best, generation):
    chemin_historique = dossier_sortie / f"calibration_{nom_jeu}_checkpoint_generation.csv"
    chemin_parametres = dossier_sortie / f"calibration_{nom_jeu}_checkpoint_parametre.csv"
    chemin_log = dossier_sortie / f"calibration_{nom_jeu}_checkpoint_log.csv"

    en_tete = "generation,best_fitness,seed_generation\n"
    lignes = [
        f"{h['generation']},{h['best_fitness']},{h.get('seed_generation', '')}\n"
        for h in historique
    ]
    with chemin_historique.open("w", encoding="utf-8") as f:
        f.write(en_tete)
        f.writelines(lignes)

    en_tete_parametres = NOMS_PARAMETRES_CALIBRES + ["best_fitness", "best_generation"]
    valeurs_parametres = [str(float(best["individual"][i])) for i in range(len(NOMS_PARAMETRES_CALIBRES))]
    valeurs_parametres.extend([str(float(best["fitness"])), str(int(best["generation"]))])
    with chemin_parametres.open("w", encoding="utf-8") as f:
        f.write(",".join(en_tete_parametres) + "\n")
        f.write(",".join(valeurs_parametres) + "\n")

    log_existe = chemin_log.exists()
    with chemin_log.open("a", encoding="utf-8") as f:
        if not log_existe:
            f.write("checkpoint_generation,best_generation,best_fitness\n")
        f.write(f"{generation},{int(best['generation'])},{float(best['fitness'])}\n")


def lire_csv_meilleurs_parametres(chemin_parametres):
    data = np.loadtxt(chemin_parametres, delimiter=",", skiprows=1)
    if data.ndim == 1:
        ligne = data
    else:
        ligne = data[0, :]

    expected = len(NOMS_PARAMETRES_CALIBRES)
    legacy = 24
    if ligne.size >= expected:
        return np.array(ligne[:expected], dtype=float)
    if ligne.size >= legacy:
        params = np.array(ligne[:legacy], dtype=float)
        # Backward compatibility: old CSV files may only contain P1..P24.
        fallback_plus = np.array([
            float(CPP_ARGS[INDICES_PARAMETRES_SUPPLEMENTAIRES["FRT"]]),
            float(CPP_ARGS[INDICES_PARAMETRES_SUPPLEMENTAIRES["NADtot"]]),
            float(CPP_ARGS[INDICES_PARAMETRES_SUPPLEMENTAIRES["Atot"]]),
            float(CPP_ARGS[INDICES_PARAMETRES_SUPPLEMENTAIRES["Cm"]]),
            float(CPP_ARGS[INDICES_PARAMETRES_SUPPLEMENTAIRES["fm"]]),
            float(CPP_ARGS[INDICES_PARAMETRES_SUPPLEMENTAIRES["kGPDH"]]),
        ], dtype=float)
        return np.concatenate([params, fallback_plus])
    raise ValueError("CSV des parametres invalide: colonnes calibrables manquantes.")


def definir_intervalle_trace(
    temps_exp,
    temps_sim,
    fenetre_trace,
):
    debut_commun = max(float(np.min(temps_exp)), float(np.min(temps_sim)))
    fin_commun = min(float(np.max(temps_exp)), float(np.max(temps_sim)))
    if debut_commun >= fin_commun:
        raise ValueError("Pas de recouvrement temporel entre l'experimental et la simulation.")

    if fenetre_trace is None:
        return debut_commun, fin_commun

    debut_demande, fin_demande = float(fenetre_trace[0]), float(fenetre_trace[1])
    debut_effectif = max(debut_demande, debut_commun)
    fin_effective = min(fin_demande, fin_commun)
    if debut_effectif >= fin_effective:
        raise ValueError("Intervalle de trace invalide apres intersection avec les donnees disponibles.")
    return debut_effectif, fin_effective


def relancer_cpp_depuis_csv_et_tracer(
    resultat,
    chemin_parametres,
    fenetre_data,
):
    chemins = resultat["paths"]
    individu = lire_csv_meilleurs_parametres(chemin_parametres)

    executer_modele_cpp(
        chemins["executable"],
        chemins["dossier_cpp"],
        preparer_arguments_cpp(
            individu,
            resultat["debut_jo_norm"],
            resultat["debut_o2_norm"],
            resultat["ajouts_adp_mm"],
            tfinal_min=float(resultat["tfinal_min"]),
            pulse_times_min=list(resultat["pulse_times_min"]),
        ),
        resultat["config"].cpp_timeout_sec,
    )

    temps_sim, o2_sim = charger_deux_colonnes(chemins["sortie_modele"])
    courbe_exp = np.load(chemins["experimental"])
    temps_exp = courbe_exp[0, :]
    o2_exp = courbe_exp[1, :]

    t_min, t_max = definir_intervalle_trace(temps_exp, temps_sim, fenetre_data)
    masque_sim = (temps_sim >= t_min) & (temps_sim <= t_max)
    masque_exp = (temps_exp >= t_min) & (temps_exp <= t_max)

    temps_sim_plot = temps_sim[masque_sim]
    temps_exp_plot = temps_exp[masque_exp]
    o2_sim_plot = o2_sim[masque_sim]
    o2_exp_plot = o2_exp[masque_exp]

    figure, axe = plt.subplots(figsize=(10, 5))
    axe.plot(temps_exp_plot, o2_exp_plot, color="tab:blue", lw=2.0, label="O2 normalisé experimental")
    axe.plot(temps_sim_plot, o2_sim_plot, color="tab:red", lw=1.8, label="O2 normalisé modele C++ (depuis CSV)")
    axe.set_xlim(t_min, t_max)
    axe.set_xlabel("Temps (min)")
    axe.set_ylabel("O2 normalisé")
    titre_fenetre = f"{t_min:.2f}-{t_max:.2f}"
    axe.set_title(f"Comparaison O2 normalisé modele C++ vs donnees experimentales ({titre_fenetre} min)")
    axe.grid(alpha=0.3)
    axe.legend()
    figure.tight_layout()

    nom_jeu = resultat["nom_jeu"]
    borne_min_nom = str(round(t_min, 2)).replace(".", "p")
    borne_max_nom = str(round(t_max, 2)).replace(".", "p")
    chemin_figure = chemins["dossier_sortie"] / f"calibration_{nom_jeu}_o2_comparison_{borne_min_nom}_{borne_max_nom}.png"
    figure.savefig(chemin_figure, dpi=150)
    plt.close(figure)
    return chemin_figure


def calibrer_un_jeu(
    cfg,
    nom_jeu,
    fichier_experimental,
    fenetres,
    fenetre_data,
    debut_jo_norm,
    debut_o2_norm,
    ajouts_adp_mm,
    tfinal_min,
    pulse_times_min,
):
    resultat, historique = lancer_calibration(
        cfg,
        nom_jeu,
        fichier_experimental,
        fenetres,
        debut_jo_norm,
        debut_o2_norm,
        ajouts_adp_mm,
        tfinal_min,
        pulse_times_min,
    )
    chemin_historique, chemin_parametres = ecrire_sorties(resultat, historique)
    chemin_comparaison = relancer_cpp_depuis_csv_et_tracer(
        resultat,
        chemin_parametres,
        fenetre_data,
    )

    print(f"[{nom_jeu}] Calibration terminee.", flush=True)
    print(f"[{nom_jeu}] Fenetre data source: {fenetre_data}", flush=True)
    print(f"[{nom_jeu}] tfinal C++: {tfinal_min} min", flush=True)
    print(f"[{nom_jeu}] Pulses ADP (min): {pulse_times_min}", flush=True)
    print(f"[{nom_jeu}] Debut jo_norm={debut_jo_norm} min, debut o2_norm={debut_o2_norm} min", flush=True)
    print(f"[{nom_jeu}] Ajouts ADP_c (mM) utilises: {ajouts_adp_mm}", flush=True)


def main():
    random_seed_run = int(time.time_ns() % (2**32))
    cfg = ConfigurationCalibration(
        population_size=100,
        generations=25,
        elite_size=25,
        poids_ratios=1.0,
        poids_points_l2=0.2,
        random_seed=random_seed_run,
    )

    print(f"Seed de base du run: {cfg.random_seed}", flush=True)

    for nom_jeu, metadonnees in JEUX_CALIBRATION.items():
        chemins_jeu = resoudre_chemins(metadonnees["fichier"])
        fenetre_npy = extraire_fenetre_temps_depuis_npy(chemins_jeu["dossier_py"], metadonnees["fichier"])

        fenetre_data = metadonnees.get("fenetre_data", fenetre_npy)
        tfinal_jeu = float(metadonnees.get("tfinal_min", fenetre_npy[1]))
        pulses_jeu = list(metadonnees.get("pulse_times_min", []))

        print(f"[{nom_jeu}] Fenetre calculee depuis NPY: {fenetre_npy}", flush=True)
        print(f"[{nom_jeu}] tfinal par defaut (max NPY): {fenetre_npy[1]:.2f} min", flush=True)
        print(f"[{nom_jeu}] Pulses ADP utilises (min): {pulses_jeu}", flush=True)
        print(f"[{nom_jeu}] Temps d'ajout ADP utilises (min): {pulses_jeu}", flush=True)
        print(f"[{nom_jeu}] Nombre de generation {cfg.generations}", flush=True)

        calibrer_un_jeu(
            cfg,
            nom_jeu,
            metadonnees["fichier"],
            metadonnees["fenetres"],
            fenetre_data,
            float(metadonnees.get("debut_jo_norm", 14.0)),
            float(metadonnees.get("debut_o2_norm", 14.0)),
            list(metadonnees.get("ajouts_adp_mM", [0.33, 1.0])),
            tfinal_jeu,
            pulses_jeu,
        )


if __name__ == "__main__":
    main()
