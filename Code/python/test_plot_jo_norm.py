#!/usr/bin/env python3

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def load_two_cols(path: Path):
    arr = np.loadtxt(path, skiprows=1)
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    return arr[:, 0], arr[:, 1]


def window_0_25(t: np.ndarray, y: np.ndarray):
    mask = (t >= 0.0) & (t <= 25.0)
    return t[mask], y[mask]


def main():
    root = Path(__file__).resolve().parents[1]
    res = root / "Cpp" / "resultats"
    out_dir = root / "python" / "plots"
    out_dir.mkdir(parents=True, exist_ok=True)

    jo_raw_candidates = [res / "jo_BE_temps.txt", res / "jo_temps.txt"]
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
    o2_path = res / "o2_pulses_2p5_10_normalized.txt"

    for p in [jo_raw_path, jo_norm_path, o2_path]:
        if not p.exists():
            raise FileNotFoundError(f"Fichier introuvable: {p}")

    t1, y1 = load_two_cols(jo_raw_path)
    t2, y2 = load_two_cols(jo_norm_path)
    t3, y3 = load_two_cols(o2_path)

    t1, y1 = window_0_25(t1, y1)
    t2, y2 = window_0_25(t2, y2)
    t3, y3 = window_0_25(t3, y3)

    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)

    axes[0].plot(t1, y1, color="black", lw=1.3)
    axes[0].set_ylabel("Jo")
    axes[0].set_title("Sorties C++ sur 0-25 min")
    axes[0].grid(alpha=0.25)

    axes[1].plot(t2, y2, color="tab:blue", lw=1.3)
    axes[1].set_ylabel("Jo_norm")
    axes[1].grid(alpha=0.25)

    axes[2].plot(t3, y3, color="tab:green", lw=1.5)
    axes[2].set_ylabel("O2_norm")
    axes[2].set_xlabel("Temps (min)")
    axes[2].set_xlim(14.0, 25.0)
    axes[2].grid(alpha=0.25)

    fig.tight_layout()
    out_png = out_dir / "cpp_outputs_0_25.png"
    fig.savefig(out_png, dpi=150)
    print(f"Plot enregistre: {out_png}")


if __name__ == "__main__":
    main()
