#!/usr/bin/env python3

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def load_two_cols(path: Path):
    # For very large outputs (dt small, long tfinal), sample lines to keep plotting responsive.
    stride = 1
    if path.stat().st_size > 50_000_000:
        stride = 20

    rows = []
    with path.open("r", encoding="utf-8") as f:
        next(f, None)  # header
        for i, line in enumerate(f):
            if i % stride != 0:
                continue
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            rows.append((float(parts[0]), float(parts[1])))

    if not rows:
        raise ValueError(f"Aucune donnee lisible dans {path}")

    arr = np.asarray(rows, dtype=float)
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    return arr[:, 0], arr[:, 1]


def load_flux(path: Path):
    stride = 1
    if path.stat().st_size > 50_000_000:
        stride = 20

    rows = []
    with path.open("r", encoding="utf-8") as f:
        next(f, None)  # header
        for i, line in enumerate(f):
            if i % stride != 0:
                continue
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            rows.append([float(v) for v in parts])

    if not rows:
        raise ValueError(f"Aucune donnee lisible dans {path}")

    return np.asarray(rows, dtype=float)


def window_10_25(t: np.ndarray, y: np.ndarray):
    mask = (t >= 10.0) & (t <= 25.0)
    return t[mask], y[mask]


def collapse_duplicate_times(t: np.ndarray, y: np.ndarray):
    """Ensure strictly increasing time axis by averaging duplicate samples."""
    if t.size <= 1:
        return t, y

    order = np.argsort(t, kind="stable")
    t_sorted = t[order]
    y_sorted = y[order]

    unique_t, inverse = np.unique(t_sorted, return_inverse=True)
    if unique_t.size == t_sorted.size:
        return t_sorted, y_sorted

    sums = np.bincount(inverse, weights=y_sorted)
    counts = np.bincount(inverse)
    y_collapsed = sums / counts
    return unique_t, y_collapsed


def collapse_duplicate_times_multi(t: np.ndarray, *ys: np.ndarray):
    """Collapse duplicate times once and apply the same aggregation to many series."""
    if t.size <= 1:
        return (t, *ys)

    order = np.argsort(t, kind="stable")
    t_sorted = t[order]
    unique_t, inverse = np.unique(t_sorted, return_inverse=True)

    if unique_t.size == t_sorted.size:
        return (t_sorted, *(y[order] for y in ys))

    counts = np.bincount(inverse)
    collapsed = []
    for y in ys:
        y_sorted = y[order]
        sums = np.bincount(inverse, weights=y_sorted)
        collapsed.append(sums / counts)

    return (unique_t, *collapsed)


def main():
    root = Path(__file__).resolve().parents[1]
    res = root / "Cpp" / "resultats"
    out_dir = root / "python" / "plots"
    out_dir.mkdir(parents=True, exist_ok=True)
    pulse_times_min = [10.0, 18.0]

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
    o2_candidates = sorted(
        res.glob("o2*_normalized.txt"),
        key=lambda p: p.stat().st_mtime,
        reverse=True,
    )
    if not o2_candidates:
        raise FileNotFoundError("Aucun fichier O2 normalise trouve dans resultats")
    o2_path = res / "o2_test_normalized.txt"
    flux_path = res / "fluxes_temps.txt"

    for p in [jo_raw_path, jo_norm_path, o2_path, flux_path]:
        if not p.exists():
            raise FileNotFoundError(f"Fichier introuvable: {p}")

    t1, y1 = load_two_cols(jo_raw_path)
    t2, y2 = load_two_cols(jo_norm_path)
    t3, y3 = load_two_cols(o2_path)
    flux = load_flux(flux_path)
    if flux.ndim == 1:
        flux = flux.reshape(1, -1)
    t4 = flux[:, 0]
    j_ant = flux[:, 1]
    if flux.shape[1] >= 5:
        j_adp_ext_sum = flux[:, 2]
        j_adp_ext_dot = flux[:, 3]
        adp_c = flux[:, 4]
    else:
        # Compatibilite anciens fichiers: seule une courbe J_ADP_ext exportee.
        j_adp_ext_sum = flux[:, 2]
        j_adp_ext_dot = flux[:, 2]
        adp_c = np.full_like(t4, np.nan)

    t1, y1 = window_10_25(t1, y1)
    t2, y2 = window_10_25(t2, y2)
    t3, y3 = window_10_25(t3, y3)
    t1, y1 = collapse_duplicate_times(t1, y1)
    t2, y2 = collapse_duplicate_times(t2, y2)
    t3, y3 = collapse_duplicate_times(t3, y3)
    if t3.size >= 3:
        o2_slope = -np.gradient(y3, t3)
    else:
        o2_slope = np.zeros_like(y3)
    t4_all = t4
    t4, j_ant = window_10_25(t4_all, j_ant)
    _, j_adp_ext_sum = window_10_25(t4_all, j_adp_ext_sum)
    _, j_adp_ext_dot = window_10_25(t4_all, j_adp_ext_dot)
    _, adp_c = window_10_25(t4_all, adp_c)
    t4, j_ant, j_adp_ext_sum, j_adp_ext_dot, adp_c = collapse_duplicate_times_multi(
        t4, j_ant, j_adp_ext_sum, j_adp_ext_dot, adp_c
    )

    fig, axes = plt.subplots(7, 1, figsize=(11, 18), sharex=True)

    axes[0].plot(t1, y1, color="black", lw=1.3)
    axes[0].set_ylabel("Jo")
    axes[0].set_xlim(14.0, 20.0)
    axes[0].set_title("Sorties C++ sur 14 à 20 min")
    axes[0].grid(alpha=0.25)

    axes[1].plot(t2, y2, color="tab:blue", lw=1.3)
    axes[1].set_ylabel("Jo_norm")
    axes[1].set_xlim(14.0, 20.0)
    axes[1].grid(alpha=0.25)

    axes[2].plot(t3, y3, color="tab:green", lw=1.5)
    axes[2].set_ylabel("O2_norm")
    axes[2].set_xlim(14.0, 20.0)
    axes[2].grid(alpha=0.25)

    axes[3].plot(t4, j_ant, color="tab:red", lw=1.3)
    axes[3].set_ylabel("J_ANT")
    axes[3].grid(alpha=0.25)

    axes[4].plot(t4, j_adp_ext_sum, color="tab:orange", lw=1.3)
    axes[4].set_ylabel("J_ADP_ext_sum")
    axes[4].set_xlim(14.0, 20.0)
    axes[4].grid(alpha=0.25)

    axes[5].plot(t4, j_adp_ext_dot, color="saddlebrown", lw=1.3)
    axes[5].set_ylabel("J_ADP_ext_dot")
    axes[5].set_xlim(14.0, 20.0)
    axes[5].grid(alpha=0.25)

    axes[6].plot(t4, adp_c, color="tab:purple", lw=1.3)
    axes[6].set_ylabel("ADP_c")
    axes[6].set_xlabel("Temps (min)")
    axes[6].set_xlim(14.0, 20.0)
    axes[6].grid(alpha=0.25)

    for ax in axes:
        for pulse_t in pulse_times_min:
            ax.axvline(pulse_t, color="gray", ls=":", lw=0.9, alpha=0.6)

    fig.tight_layout()
    out_png = out_dir / "cpp_outputs_10_25.png"
    fig.savefig(out_png, dpi=150)
    print(f"Plot enregistre: {out_png}")


if __name__ == "__main__":
    main()
