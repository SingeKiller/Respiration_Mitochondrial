#!/usr/bin/env python3

from pathlib import Path
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


CPP_READY = True

CPP_ARGS = [
    "400", "1", "0.01", "0.6", "100", "177", "5", "7", "100", "177",
    "5", "120", "10000", "190", "8.5", "35", "0.002", "-0.03", "0.35",
    "2", "0.01", "1.1", "0.001", "0.016", "0.037", "10000", "15000",
    "1.8", "0.01", "0.0005", "0", "0.1", "15000", "0.1", "1000", "21",
    "16.15", "330", "19.15", "1000", "o2_test_normalized.txt"
]



SECTIONS = {
    "data14h16.npy": [[14.6, 15.5], [16.25, 16.8], [17.6, 18.5], [19.7, 20.5]],
    "data15h12.npy": [[10.8, 12.8], [13.3, 14.35], [16.0, 17.0], [17.25, 17.8]],
    "data15h37.npy": [[10.6, 11.3], [12.3, 12.6], [13.1, 13.55], [14.3, 15.7]],
    "data16h00.npy": [[10.75, 11.75], [12.5, 13.75], [14.25, 15.0], [15.4, 16.3]],
    "data17h03.npy": [[10.3, 11.45], [12.2, 12.85], [13.75, 14.8], [15.4, 17.0]],
}


def load_two_cols(path: Path): 
    arr = np.loadtxt(path, skiprows=1)
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    return arr[:, 0], arr[:, 1]


def window_0_25(t: np.ndarray, y: np.ndarray):
    mask = (t >= 0.0) & (t <= 25.0)
    return t[mask], y[mask]


def run_cpp_if_enabled(root: Path):
    if not CPP_READY:
        return

    cpp_dir = root / "Cpp"
    exe = cpp_dir / "main.exe"
    if not exe.exists():
        exe = cpp_dir / "main"

    subprocess.run([str(exe), *CPP_ARGS], check=True, cwd=str(cpp_dir))


def regression_segments(courbe: np.ndarray, decoupage: list[list[float]]):
    droite_a_b = []
    segments = []

    for t_min, t_max in decoupage:
        mask = (courbe[0, :] > t_min) & (courbe[0, :] < t_max)
        x = courbe[0, :][mask]
        y = courbe[1, :][mask]
        reg = stats.linregress(x, y)
        droite = reg.slope * x + reg.intercept
        droite_a_b.append([float(reg.slope), float(reg.intercept)])
        segments.append((x, y, droite))

    return droite_a_b, segments


def plot_regression_overview(py_dir: Path, out_dir: Path):
    data_dir = py_dir / "normalized_respiration_data"
    coeffs_by_file = {}

    for file_name, windows in SECTIONS.items():
        p = data_dir / file_name
        if not p.exists():
            raise FileNotFoundError(f"Fichier introuvable: {p}")

        courbe = np.load(p)
        coeffs, segments = regression_segments(courbe, windows)
        coeffs_by_file[file_name] = coeffs

        fig, axes = plt.subplots(1, len(windows) + 1, figsize=(18, 4.5))

        axes[0].plot(courbe[0, :], courbe[1, :], color="blue", linestyle=":")
        for i, (x, y, droite) in enumerate(segments, start=1):
            axes[0].plot(x, droite, color="red")
            axes[i].plot(x, y, color="blue", linestyle=":")
            axes[i].plot(x, droite, color="red")
            axes[i].set_title(f"Segment {i}")
            axes[i].grid(alpha=0.25)

        axes[0].set_title(f"Courbe globale - {file_name}")
        axes[0].grid(alpha=0.25)
        fig.tight_layout()

        out_png = out_dir / f"regression_{Path(file_name).stem}.png"
        fig.savefig(out_png, dpi=150)
        plt.close(fig)

    return coeffs_by_file


def plot_cpp_outputs(root: Path, out_dir: Path):
    res = root / "Cpp" / "resultats"

    jo_raw_candidates = [res / "jo_temps.txt"]
    jo_raw_path = None
    for cand in jo_raw_candidates:
        if cand.exists():
            jo_raw_path = cand
            break
    if jo_raw_path is None:
        raise FileNotFoundError(
            f"Fichier introuvable: aucun des candidats {[str(p) for p in jo_raw_candidates]}"
        )

    jo_norm_path = res / "jo_normalized.txt"
    o2_candidates = [
        res / "o2_avec_pulse_normalized.txt",  # ancien nom
        res / "o2_test_normalized.txt",
        res / "o2_tmp_dt10_precise_normalized.txt",
        res / "o2_tmp_dt10_substep_normalized.txt",
        res / "o2_tmp_dt1ms_normalized.txt",
    ]
    o2_path = None
    for cand in o2_candidates:
        if cand.exists():
            o2_path = cand
            break

    missing = [p for p in [jo_raw_path, jo_norm_path] if not p.exists()]
    if missing:
        raise FileNotFoundError(f"Fichier introuvable: {missing[0]}")
    if o2_path is None:
        raise FileNotFoundError(
            f"Fichier introuvable: aucun des candidats {[str(p) for p in o2_candidates]}"
        )

    t1, y1 = load_two_cols(jo_raw_path)
    t2, y2 = load_two_cols(jo_norm_path)
    t3, y3 = load_two_cols(o2_path)

    t1, y1 = window_0_25(t1, y1)
    t2, y2 = window_0_25(t2, y2)
    t3, y3 = window_0_25(t3, y3)

    fig, axes = plt.subplots(1, 1, figsize=(10, 8))
    
    plt.plot(t3, y3, color="tab:green", lw=1.5)
    plt.ylabel("O2_norm")
    plt.xlabel("Temps (min)")
    plt.xlim(14.0, 21.0)
    plt.grid(alpha=0.25)

    fig.tight_layout()
    out_png = out_dir / "Sortie_Cpp.png"
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    return out_png


def main():
    plt.ion()
    plt.close('all')
    root = Path(__file__).resolve().parents[1]
    out_dir = root / "python" / "plots"
    py_dir = root / "python"
    out_dir.mkdir(parents=True, exist_ok=True)

    run_cpp_if_enabled(root)
    coeffs_by_file = plot_regression_overview(py_dir, out_dir)
    out_png = plot_cpp_outputs(root, out_dir)

    print(f"Plot enregistre: {out_png}")
    for file_name, coeffs in coeffs_by_file.items():
        print(f"Regressions {file_name}: {coeffs}")


if __name__ == "__main__": 
    main()

"""
===================================================
algo génétiques
=================================================

"""


def initialisation():
    # on initialise une population en focntion des parametre prédéfini dans bertram
    genes_base = np.array(CPP_ARGS[:24], dtype=float)
    population = []
    for _ in range(100):
        mutation_initiale = np.random.uniform(0.5, 1.5, size=24)
        individu = genes_base * mutation_initiale
        population.append(individu)

    population = np.array(population)
    return population


def creer_nouvelle_generation(parents, taille_cible=100):
    nouveaux_individus = []
    
    while len(nouveaux_individus) < taille_cible:
        #on prends deux parent on fait la moyenne et on ajoute une modification génétique
        indices = np.random.randint(0, len(parents), size=2)
        p1, p2 = parents[indices[0]], parents[indices[1]]
        enfant = (p1 + p2) / 2
        
        coefficients = np.random.uniform(0.95, 1.05, size=24)
        enfant = enfant * coefficients
        
        nouveaux_individus.append(enfant)
        
    return np.array(nouveaux_individus)

def preparer_commande(individu_genes):
    #on ajoute les paremtres en str pour avoir totu les parammetres
    #à mettre dans un boucle 
    genes_str = individu_genes.astype(str).tolist()
    params_fixes = CPP_ARGS[24:]
    commande_finale = genes_str + params_fixes
    return commande_finale

def analyser_courbe_relatif(chemin_fichier, decoupage):
    """
    Lit le fichier O2, calcule les angles de chaque segment et renvoie le tableau des angles RELATIFS (différences).
    """
    t, y = load_two_cols(chemin_fichier)
    courbe = np.vstack((t, y))
    angles_absolus = []
    #calcule des angles par segment
    for t_min, t_max in decoupage:
        masque = (courbe[0, :] > t_min) & (courbe[0, :] < t_max)
        reg = stats.linregress(courbe[0, :][masque], courbe[1, :][masque])
        angles_absolus.append(np.arctan(reg.slope))
    # [b-a, c-b]
    angles_relatifs = np.diff(angles_absolus)
    return np.array(angles_relatifs)

def calculer_score_fitness(angles_reels, angles_simules):
    ecarts = np.abs(angles_reels - angles_simules)
    score = np.sum(ecarts)
    return float(score)

def filtrer_meilleurs(population, scores):
    indices_tries = np.argsort(scores)
    indices_gagnants = indices_tries[:50]
    meilleurs_individus = population[indices_gagnants]
    
    return meilleurs_individus


def boucle_optimisation(nb_generations=100):
    # --- INITIALISATION ---
    population = initialisation() # Ton tableau (100, 24)
    cpp_dir = root / "Cpp"
    exe = cpp_dir / "main.exe"
    fichier_o2 = cpp_dir / "resultats" / "o2_test_normalized.txt"

    for gen in range(nb_generations):
        scores = []
        for individu in population:
            args_cpp = preparer_commande(individu)
            
            subprocess.run([str(exe), *args_cpp], cwd=str(cpp_dir), stdout=subprocess.DEVNULL)
            # au cas ou  il y a un pb on elimine d'office les parametre avec des pbs
            try:
                angles_simules = analyser_courbe_relatif(fichier_o2, decoupage_cible)
                score = calculer_score_fitness(angles_cibles, angles_simules)
            except:
                score = 999
            
            scores.append(score)

        scores = np.array(scores)
        
        #on séléctionne et on les repoduit 
        meilleurs = filtrer_meilleurs(population, scores)
        population = creer_nouvelle_generation(meilleurs, taille_cible=100)
        
        print(f"Génération {gen} : Meilleur score = {np.min(scores):.4f}")
    #résultat finale
    return population[np.argmin(scores)]