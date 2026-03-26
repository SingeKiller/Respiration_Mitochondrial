#!/usr/bin/env python3

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def load_o2_model(path: Path):
    rows = []
    with path.open("r", encoding="utf-8") as f:
        next(f, None)  # skip header
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            rows.append((float(parts[0]), float(parts[1])))

    if not rows:
        raise ValueError(f"Aucune donnee O2 lisible dans {path}")

    arr = np.asarray(rows, dtype=float)
    return arr[:, 0], arr[:, 1]


def load_npy_trace(path: Path):
    data = np.load(path)
    if data.ndim != 2:
        raise ValueError(f"Format inattendu pour {path}: shape={data.shape}")

    if data.shape[0] == 2:
        t = data[0]
        y = data[1]
    elif data.shape[1] == 2:
        t = data[:, 0]
        y = data[:, 1]
    else:
        raise ValueError(f"Impossible d'interpreter {path}: shape={data.shape}")

    return np.asarray(t, dtype=float), np.asarray(y, dtype=float)


def window_trace(t: np.ndarray, y: np.ndarray, t_min: float, t_max: float):
    mask = (t >= t_min) & (t <= t_max)
    return t[mask], y[mask]


def rebase_time_and_normalize(t: np.ndarray, y: np.ndarray, t_min: float, t_max: float):
    if t.size < 2:
        return t, y

    # Ensure monotonic time for robust remapping.
    order = np.argsort(t, kind="stable")
    t = t[order]
    y = y[order]

    t0 = float(t[0])
    t1 = float(t[-1])
    if abs(t1 - t0) < 1e-12:
        t_rebased = np.full_like(t, t_min, dtype=float)
    else:
        t_rebased = t_min + (t - t0) * (t_max - t_min) / (t1 - t0)

    y_start = float(y[0])
    y_end = float(y[-1])
    denom = y_start - y_end
    if abs(denom) < 1e-12:
        y_norm = np.ones_like(y, dtype=float)
    else:
        y_norm = (y - y_end) / denom

    return t_rebased, y_norm


def main():
    py_dir = Path(__file__).resolve().parent
    root = py_dir.parent
    res_dir = root / "Cpp" / "resultats"
    out_dir = py_dir / "plots"
    out_dir.mkdir(parents=True, exist_ok=True)

    model_candidates = [
        res_dir / "o2_test_normalized.txt",
        res_dir / "o2_tmp_dt10_substep_normalized.txt",
        res_dir / "o2_tmp_dt10_precise_normalized.txt",
        res_dir / "o2_tmp_dt1ms_normalized.txt",
    ]
    model_path = None
    for candidate in model_candidates:
        if candidate.exists():
            model_path = candidate
            break
    if model_path is None:
        matches = sorted(res_dir.glob("o2*_normalized.txt"), key=lambda p: p.stat().st_mtime, reverse=True)
        if not matches:
            raise FileNotFoundError("Aucun fichier O2 normalise trouve dans Cpp/resultats")
        model_path = matches[0]

    data_paths = [
        py_dir / "data14h16.npy",
        py_dir / "data15h12.npy",
        py_dir / "data15h37.npy",
        py_dir / "data16h00.npy",
        py_dir / "data17h03.npy",
    ]
    for data_path in data_paths:
        if not data_path.exists():
            raise FileNotFoundError(f"Fichier introuvable: {data_path}")

    t_start = 14.0
    t_end = 21.0

    t_model, y_model = load_o2_model(model_path)
    t_model, y_model = window_trace(t_model, y_model, t_start, t_end)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(t_model, y_model, lw=2.0, color="tab:green", label=f"Modele O2_norm ({model_path.name})")

    data_colors = ["black", "tab:blue", "tab:orange", "tab:red", "tab:purple"]
    for data_path, color in zip(data_paths, data_colors):
        t_data, y_data = load_npy_trace(data_path)
        t_data, y_data = window_trace(t_data, y_data, t_start, t_end)
        t_data, y_data = rebase_time_and_normalize(t_data, y_data, t_start, t_end)
        ax.plot(t_data, y_data, lw=1.6, alpha=0.9, color=color, label=f"Donnee rebassee {data_path.name}")

    ax.set_title("O2 modele + 5 donnees rebassees (14-21, normalisees)")
    ax.set_xlabel("Temps (min)")
    ax.set_ylabel("Concentration O2 normalisee")
    ax.set_xlim(t_start, t_end)
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()

    out_png = out_dir / "o2_norm_vs_5_data_overlay.png"
    fig.savefig(out_png, dpi=150)
    print(f"Graphe enregistre: {out_png}")

    plt.show()


if __name__ == "__main__":
    main()
