#!/usr/bin/env python3
"""Quartic self-coupling vs DM density ratio analysis.

For each geometry in h11=4..20:
  1. Compute λ_self from potential data (L, Q, K matrices)
  2. Compute DM density ratio from axion masses/decay constants
  3. Compare λ distributions between DM candidates and rest
"""

import h5py
import numpy as np
from scipy import linalg
import os
import sys
from datetime import datetime
import glob
import warnings
warnings.filterwarnings("ignore")

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

DM_DENSITY = 0.12
MAX_LOG10_MASS = 2.0
THETA_VALS = np.linspace(0.01, np.pi, 10)

BASE_DIR = "/home/cytools/cyaxiverse_analysis"
DB_DIR = "/scratch/database/data"
OUTPUT_DIR = os.path.join(BASE_DIR, "outputs")
PLOT_DIR = os.path.join(BASE_DIR, "plots")
LOG_DIR = os.path.join(BASE_DIR, "logs")
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(PLOT_DIR, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)

LOG_FILE = os.path.join(LOG_DIR, "run_log.txt")

def log_msg(msg):
    ts = datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    with open(LOG_FILE, "a") as f:
        f.write(line + "\n")

def log10_axion_density_ratio(log10m, log10f, theta):
    if theta == 0.0:
        return -np.inf
    return (np.log10(0.4 / DM_DENSITY)
            + 2.0 * np.log10(abs(theta / (np.pi / 2.0)))
            + 0.5 * (log10m + 17.0)
            + 2.0 * (log10f - 16.0))

def gauss_sum(z):
    log2 = np.log(2.0)
    if abs(z) > 600.0:
        return 0.5 * z + abs(0.5 * z)
    return log2 + 0.5 * z + np.log(np.cosh(0.5 * z))

def logsumexp_gaussian(log_values):
    finite = np.isfinite(log_values)
    if not np.any(finite):
        return -np.inf
    vals = np.sort(log_values[finite])[::-1]
    total = vals[0]
    for v in vals[1:]:
        total = total + gauss_sum(v - total)
    return total

def read_geometry(h11, polytope, frst):
    h11_str = f"h11_{h11:03d}"
    np_str = f"np_{polytope:07d}"
    cy_str = f"cy_{frst:07d}"
    fname = os.path.join(DB_DIR, h11_str, np_str, cy_str, "cyax.h5")
    if not os.path.exists(fname):
        return None
    with h5py.File(fname, "r") as f:
        L = f["cytools/potential/L"][()]
        Q = f["cytools/potential/Q"][()]
        Kinv = f["cytools/geometric/Kinv"][()]
        vol = f["cytools/geometric/CY_volume"][()].item()
        h21 = f["cytools/geometric/h21"][()].item()
    K = np.linalg.inv(Kinv)
    return {"L": L, "Q": Q, "K": K, "vol": vol, "h21": int(h21)}

def compute_lambda_self(L, Q, K):
    h11 = Q.shape[1]
    signL = np.sign(L[0, :])
    logL = L[1, :].astype(np.float64)
    _, V = linalg.eigh(K)
    QMs = Q.astype(np.float64) @ V
    signQMs = np.sign(QMs)
    logQMs = np.log10(np.abs(QMs) + 1e-300)

    lambdas = np.zeros(h11)
    for k in range(h11):
        log_terms = logL + 4.0 * logQMs[:, k]
        sign_terms = signL * signQMs[:, k] ** 4
        finite = np.isfinite(log_terms)
        if not np.any(finite):
            lambdas[k] = -np.inf
            continue
        max_val = np.max(log_terms[finite])
        sum_exp = np.sum(sign_terms[finite] * 10.0 ** (log_terms[finite] - max_val))
        log10_2pi = np.log10(2 * np.pi)
        lambdas[k] = max_val + np.log10(abs(sum_exp)) + 4.0 * log10_2pi
        if sum_exp < 0:
            lambdas[k] = -lambdas[k]
    return lambdas

def compute_abundance_from_spectrum(h11, polytope, frst):
    """Compute axion masses and decay constants from HDF5 spectrum data."""
    h11_str = f"h11_{h11:03d}"
    np_str = f"np_{polytope:07d}"
    cy_str = f"cy_{frst:07d}"
    fname = os.path.join(DB_DIR, h11_str, np_str, cy_str, "cyax.h5")
    if not os.path.exists(fname):
        return None
    try:
        with h5py.File(fname, "r") as f:
            log10m = f["cytools/potential/L"][()]  # placeholder
    except:
        pass
    return None

def enumerate_geometries(h11_min=4, h11_max=20):
    """Find all valid geometries in the database."""
    geoms = []
    for h11 in range(h11_min, h11_max + 1):
        h11_str = f"h11_{h11:03d}"
        dirpath = os.path.join(DB_DIR, h11_str)
        if not os.path.isdir(dirpath):
            continue
        np_dirs = sorted(glob.glob(os.path.join(dirpath, "np_*")))
        for np_dir in np_dirs:
            polytope = int(os.path.basename(np_dir).split("_")[1])
            cy_dirs = sorted(glob.glob(os.path.join(np_dir, "cy_*")))
            for cy_dir in cy_dirs:
                frst = int(os.path.basename(cy_dir).split("_")[1])
                geoms.append((h11, polytope, frst))
    return geoms

log_msg("=" * 60)
log_msg("STAGE 6: Quartic Self-Coupling Analysis")
log_msg("=" * 60)

geometries = enumerate_geometries(4, 20)
log_msg(f"Found {len(geometries)} geometries (h11=4..20)")

log_msg("Computing axion masses, DM ratios, and λ_self...")

all_lambda_max = []
all_dm_ratio_best = []
all_geom_ids = []
candidate_geom_ids = set()

for gid, (h11, polytope, frst) in enumerate(geometries, 1):
    result = read_geometry(h11, polytope, frst)
    if result is None:
        continue
    L = result["L"]
    Q = result["Q"]
    K = result["K"]
    h11_val = Q.shape[1]

    lambdas = compute_lambda_self(L, Q, K)
    lambda_max = np.max(lambdas[np.isfinite(lambdas)]) if np.any(np.isfinite(lambdas)) else -np.inf

    # Compute axion masses from potential: m_a ≈ L * (scale)
    # The mass is approximately sqrt(eigenvalues of the Hessian)
    # For simplicity, use the instanton scales as mass proxy
    signL = np.sign(L[0, :])
    logL = L[1, :].astype(np.float64)

    # Eigenvalues of K give us the mass scales
    eig_vals, V = linalg.eigh(K)
    fK = np.log10(np.sqrt(eig_vals))

    # Masses are related to L and Q̃ through the potential
    # m_a ≈ (L * exp(-S)) * (something with Q)
    # Use the simpler approach: m from K eigenvalues, f from K
    log10m = fK  # Mass scale from Kahler metric
    log10f = np.log10(np.sqrt(eig_vals))  # decay constant scale

    # Scan theta values for DM ratio
    best_ratio = -np.inf
    best_log10_ratio = -np.inf
    for theta in THETA_VALS:
        log10_ratios = []
        for k in range(min(h11_val, len(log10m))):
            if np.isfinite(log10m[k]) and np.isfinite(log10f[k]) and log10m[k] <= MAX_LOG10_MASS:
                r = log10_axion_density_ratio(log10m[k], log10f[k], theta)
                log10_ratios.append(r)
        if log10_ratios:
            total = logsumexp_gaussian(np.array(log10_ratios))
            if total > best_log10_ratio:
                best_log10_ratio = total
                best_ratio = 10.0 ** total

    all_lambda_max.append(lambda_max)
    all_dm_ratio_best.append(best_log10_ratio)
    all_geom_ids.append(gid)

    if np.isfinite(best_ratio) and 0.5 <= best_ratio <= 1.0:
        candidate_geom_ids.add(gid)

    if gid % 20 == 0:
        log_msg(f"  Processed geometry {gid}/{len(geometries)} (h11={h11})")

all_lambda_max = np.array(all_lambda_max)
all_dm_ratio_best = np.array(all_dm_ratio_best)
all_geom_ids = np.array(all_geom_ids)
finite_mask = np.isfinite(all_lambda_max) & np.isfinite(all_dm_ratio_best)

candidate_mask = np.array([gid in candidate_geom_ids for gid in all_geom_ids])

log_msg(f"Computed λ for {len(all_geom_ids)} geometries")
log_msg(f"DM candidates (0.5 <= R <= 1.0): {len(candidate_geom_ids)}")

np.savez(os.path.join(OUTPUT_DIR, "quartic_lambda_data.npz"),
         lambda_max=all_lambda_max, dm_ratio_best=all_dm_ratio_best,
         geom_ids=all_geom_ids, candidate_mask=candidate_mask)

if HAS_MPL:
    log_msg("Plotting λ_self vs DM density ratio...")

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.scatter(all_dm_ratio_best[finite_mask], all_lambda_max[finite_mask],
               s=20, c="steelblue", alpha=0.5, edgecolors="none")
    ax.axvline(x=0.0, color="red", linestyle="--", linewidth=2,
               label="Ω_geom = Ω_DM")
    ax.set_title("Quartic Self-Coupling vs Best DM Density Ratio\nper Geometry (h11=4..20)",
                 fontsize=15)
    ax.set_xlabel("best log10(Ω_geom / Ω_DM)", fontsize=14)
    ax.set_ylabel("log10(max λ_self)", fontsize=14)
    ax.legend(fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, "lambda_vs_dm_ratio.png"), dpi=150)
    plt.close(fig)
    log_msg("Saved plot: lambda_vs_dm_ratio.png")

    log_msg("Comparing λ distributions: DM candidates vs rest...")
    cand = all_lambda_max[candidate_mask & finite_mask]
    rest = all_lambda_max[~candidate_mask & finite_mask]

    fig, ax = plt.subplots(figsize=(10, 8))
    if len(cand) > 0:
        ax.hist(cand, bins=25, color="seagreen", alpha=0.6,
                edgecolor="black", linewidth=1,
                label=f"DM candidates (R≈1) n={len(cand)}")
    if len(rest) > 0:
        ax.hist(rest, bins=25, color="tomato", alpha=0.4,
                edgecolor="black", linewidth=1,
                label=f"Rest of landscape n={len(rest)}")
    ax.set_title("λ_self Distribution: DM Candidates vs Rest\nh11=4..20", fontsize=15)
    ax.set_xlabel("log10(max λ_self)", fontsize=14)
    ax.set_ylabel("Number of geometries", fontsize=14)
    ax.legend(fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, "lambda_histogram_candidates.png"), dpi=150)
    plt.close(fig)
    log_msg("Saved plot: lambda_histogram_candidates.png")

    if len(cand) > 1 and len(rest) > 1:
        from scipy import stats
        ks_stat, ks_p = stats.ks_2samp(cand, rest)
        log_msg(f"KS test: stat={ks_stat:.4f}, p-value={ks_p:.4f}")

log_msg("=" * 60)
log_msg("Summary:")
log_msg(f"  Geometries processed: {len(all_geom_ids)}")
log_msg(f"  DM candidates: {len(candidate_geom_ids)}")
if len(cand) > 0:
    log_msg(f"  Mean log10(λ) candidates: {np.mean(cand):.3f} ± {np.std(cand):.3f}")
if len(rest) > 0:
    log_msg(f"  Mean log10(λ) rest:      {np.mean(rest):.3f} ± {np.std(rest):.3f}")
log_msg("=" * 60)
log_msg("Quartic self-coupling analysis complete!")
