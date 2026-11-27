"""
Theorem D16 — Attacking the Two Residual Zones of Lee's Conjecture on m = pq:
  Zone 1 (Very dense signed graphs, rho in [m^{-1/3}, 1]):
    D16.1  Spectral Tail Dichotomy & Hanson-Wright Elimination of Dense Expanders:
           For any A_UW in {-1, 0, +1}^{n/2 x n/2} and non-wrapping frequency 1 <= |t| <= T_nw,
           either there exists r <= O(ln^2 m) with E_tail(r) >= 16 m^2 ln m / t^2 and
           E_tail(r) / sigma_{r+1}^2 >= C_0 ln m (giving |S_{1/2}(t)| < 1/m via Hanson-Wright),
           or A_UW admits a real rank-r* approximation ||A_UW - L||_F^2 < 16 m^2 ln m / t^2.
           Benchmarked at rho = 0.5 and rho = 1.0 on random +-1, random {0,1}, Sylvester-Hadamard,
           and Paley signed graphs.
  Zone 2 (Arbitrary modular quadratics a_ij, b_i in Z/mZ at t = 1):
    D16.2  Dyadic Coefficient Exclusion & Independent-Set Reduction:
           At every dyadic scale E_k = {{i,j} : 2^{k-1} <= |t a_ij mod+- m| < 2^k},
           |S_{1/2}(t)| <= exp(-(1/4) nu_ind(G_k) (2^{k-1}/m)^2) <= exp(-(|E_k|/(8 Delta_k^2))(2^{k-1}/m)^2).
           In particular, for large coefficients (|t a_ij mod+- m| >= m/4), |S_{1/2}(t)| > 1/m forces
           nu_ind(G_t(m/4)) <= 64 ln m and alpha(G_t(m/4)) >= n - 128 Delta_D ln m, killing all
           sparse/medium-degree large-coefficient modular quadratics at n = O(Delta ln m) (independent of m!).
"""

import json
import math
import os
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_PATH = os.path.join(BASE_DIR, "results", "theorem_d16_zone1_zone2_attack.json")
FIGURE_PATH = os.path.join(BASE_DIR, "figures", "theorem_d16_zone1_zone2_attack.png")


def sylvester_hadamard(k):
    """Build 2^k x 2^k Sylvester-Hadamard matrix in {-1, +1} without scipy."""
    H = np.array([[1]], dtype=np.int8)
    for _ in range(k):
        H = np.block([[H, H], [H, -H]])
    return H


def paley_signed_matrix(p):
    """Build p x p Paley signed adjacency matrix (Legendre symbol (i-j | p)) for prime p = 1 mod 4."""
    qr = set((x * x) % p for x in range(1, p))
    A = np.zeros((p, p), dtype=np.int8)
    for i in range(p):
        for j in range(p):
            if i != j:
                A[i, j] = 1 if ((i - j) % p) in qr else -1
    return A


def spectral_tail_stats(A, m):
    """Compute singular values, best tail ratio E_tail(r)/sigma_{r+1}^2 subject to E_tail(r) >= m^2."""
    s = np.linalg.svd(A.astype(float), compute_uv=False)
    s2 = s**2
    frob2 = float(np.sum(s2))
    # Tail energy E_tail(r) = sum_{j >= r} s2[j] (0-indexed)
    tail = np.cumsum(s2[::-1])[::-1]
    best_ratio = 0.0
    best_r = 0
    for r in range(len(s2)):
        if tail[r] >= m**2 and s2[r] > 1e-12:
            ratio = float(tail[r] / s2[r])
            if ratio > best_ratio:
                best_ratio = ratio
                best_r = r
    return {
        "frob2_over_m2": round(frob2 / (m**2), 3),
        "r_stab_raw": round(float(frob2**2 / np.sum(s2**2)), 2),
        "best_tail_r": int(best_r),
        "best_tail_ratio": round(best_ratio, 2),
    }


def measure_decoupled_phi(A, m, n_samples=2500, seed=42):
    """Measure max_{1 <= t <= m/2} E_y[prod_w |cos(pi t (A^T y)_w / m)|] for y in {-1, 0, +1}^u."""
    rng = np.random.default_rng(seed)
    u, w = A.shape
    y = rng.integers(0, 2, size=(n_samples, u), dtype=np.int8) - rng.integers(0, 2, size=(n_samples, u), dtype=np.int8)
    Z = y.astype(np.int32) @ A.astype(np.int32)
    max_phi = 0.0
    worst_t = 1
    for t in range(1, m // 2 + 1):
        c = np.abs(np.cos((np.pi * t / m) * Z))
        # Use log-sum for numerical stability when products are tiny
        with np.errstate(divide="ignore"):
            log_p = np.sum(np.log(np.maximum(c, 1e-300)), axis=1)
        phi = float(np.mean(np.exp(np.maximum(log_p, -700.0))))
        if phi > max_phi:
            max_phi = phi
            worst_t = t
    return max_phi, worst_t


def benchmark_zone1_dense_signed():
    """Zone 1 benchmarks: rho = 0.5 and rho = 1.0 dense signed graphs vs 1/m^2 threshold."""
    m = 35
    rng = np.random.default_rng(16)
    configs = [
        ("Rand +-1 (rho=1, 140x140)", rng.choice([-1, 1], size=(140, 140)).astype(np.int8)),
        ("Rand 0/1 (rho=0.5, 140x140)", rng.integers(0, 2, size=(140, 140), dtype=np.int8)),
        ("Hadamard H_128 (rho=1)", sylvester_hadamard(7)),
        ("Hadamard H_256 (rho=1)", sylvester_hadamard(8)),
        ("Paley P_197 (rho=1)", paley_signed_matrix(197)),
        ("Rank-1 J_140 (control)", np.ones((140, 140), dtype=np.int8)),
    ]
    results = []
    for label, A in configs:
        st = spectral_tail_stats(A, m)
        max_phi, worst_t = measure_decoupled_phi(A, m, n_samples=2500)
        killed = bool(max_phi < 1.0 / (m**2))
        row = {
            "label": label,
            "shape": list(A.shape),
            "m": m,
            **st,
            "max_phi": max_phi,
            "worst_t": int(worst_t),
            "threshold_1_over_m2": 1.0 / (m**2),
            "killed_below_1_over_m2": killed,
        }
        results.append(row)
        print(f"[Zone 1] {label:28s} | E/m^2={st['frob2_over_m2']:5.1f} | tail_r={st['best_tail_r']} (ratio={st['best_tail_ratio']:5.1f}) | "
              f"max_phi={max_phi:.3e} (t={worst_t:2d}) | <1/m^2? {killed}")
    return results


def greedy_induced_matching_size(n, edges):
    """Greedy lower bound on induced matching number nu_ind(G)."""
    adj = [set() for _ in range(n)]
    for u, v in edges:
        adj[u].add(v)
        adj[v].add(u)
    alive = np.ones(n, dtype=bool)
    matching = 0
    for u, v in edges:
        if alive[u] and alive[v]:
            matching += 1
            to_kill = {u, v} | adj[u] | adj[v]
            for w in to_kill:
                alive[w] = False
    return matching


def benchmark_zone2_dyadic_modular():
    """
    Zone 2 benchmarks: arbitrary large-coefficient modular quadratics a_ij ~ Unif(Z/mZ)
    on sparse/medium graphs (Delta in {4, 8, 16}) at n << m (e.g. n = 160 for m = 899!).
    Verify:
      1. Greedy induced matching in G_t(m/4) exceeds 16 ln m -> certified |S_{1/2}(t)| << 1/m
      2. Empirical Pr_{1/2}[P == 0 mod m] ~ 1/m and max_t |S_{1/2}(t)| << 1/m
    """
    rng = np.random.default_rng(160)
    results = []
    for m, n, deg in [(35, 160, 4), (143, 280, 6), (323, 380, 6), (899, 460, 6)]:
        # Build random approx deg-regular graph
        edges = set()
        for _ in range(deg):
            perm = rng.permutation(n)
            for i in range(0, n - 1, 2):
                u, v = int(perm[i]), int(perm[i + 1])
                if u > v:
                    u, v = v, u
                edges.add((u, v))
        edges = sorted(edges)
        eu = np.array([e[0] for e in edges], dtype=np.int32)
        ev = np.array([e[1] for e in edges], dtype=np.int32)
        # Arbitrary modular coefficients a_ij in {1..m-1} and arbitrary b_i in {0..m-1}
        a_coeff = rng.integers(1, m, size=len(edges), dtype=np.int64)
        b_coeff = rng.integers(0, m, size=n, dtype=np.int64)

        # Check induced matching in G_t(m/4) at t = 1 and worst over 20 random unit frequencies
        test_ts = [1] + [int(t) for t in rng.integers(2, m // 2 + 1, size=19) if math.gcd(int(t), m) == 1]
        min_nu_ind = 10**9
        for t in test_ts:
            # centered residue |t * a_ij mod+- m| >= m/4
            res = (t * a_coeff) % m
            res = np.minimum(res, m - res)
            large_idx = np.where(res >= m / 4.0)[0]
            large_edges = [(int(eu[i]), int(ev[i])) for i in large_idx]
            nu = greedy_induced_matching_size(n, large_edges)
            if nu < min_nu_ind:
                min_nu_ind = nu

        # Certified Band-3 upper bound on |S_{1/2}(t)| from nu_ind on G_t(m/4):
        # sin(pi/4)^2 = 1/2, pre-Weyl factor (1/4) * sin^2 >= 1/8 => exp(-nu / 8)
        cert_bound = math.exp(-min_nu_ind / 8.0)

        # Empirical Pr_{1/2}[P == 0 mod m] and max_t |S_{1/2}(t)|
        n_samp = 12000
        X = rng.integers(0, 2, size=(n_samp, n), dtype=np.int8)
        P_vals = (X @ b_coeff) % m
        P_vals = (P_vals + ((X[:, eu] * X[:, ev]) @ a_coeff)) % m
        pr_zero = float(np.mean(P_vals == 0))
        phases = (2.0 * np.pi / m) * P_vals
        max_S_emp = max(float(np.abs(np.mean(np.exp(1j * t * phases)))) for t in test_ts)

        row = {
            "m": m, "n": n, "n_over_m": round(n / m, 3), "target_deg": deg, "n_edges": len(edges),
            "min_nu_ind_large": int(min_nu_ind),
            "certified_S_ub": cert_bound,
            "target_1_over_m": 1.0 / m,
            "certified_killed": bool(cert_bound < 1.0 / m),
            "pr_zero_emp": round(pr_zero, 5),
            "max_S_emp": round(max_S_emp, 5),
        }
        results.append(row)
        print(f"[Zone 2] m={m:3d}, n={n:3d} (n/m={n/m:.2f}), deg={deg} | min nu_ind(G_t(m/4))={min_nu_ind:2d} | "
              f"cert |S|<={cert_bound:.2e} vs 1/m={1/m:.2e} ({row['certified_killed']}) | Pr0={pr_zero:.5f}")
    return results


def make_figure(z1, z2):
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.0), dpi=160)
    plt.rcParams["font.size"] = 9.2

    # Panel A: Zone 1 decoupled max_phi vs 1/m^2
    ax = axes[0, 0]
    labels = [r["label"].replace(" (rho=1, 140x140)", "\n(140x140)").replace(" (rho=0.5, 140x140)", "\n(140x140)").replace(" (rho=1)", "").replace(" (control)", "\n(ctrl)") for r in z1]
    vals = [max(r["max_phi"], 1e-12) for r in z1]
    colors = ["#2ca02c" if r["killed_below_1_over_m2"] else "#d62728" for r in z1]
    x = np.arange(len(z1))
    ax.bar(x, vals, color=colors, alpha=0.88)
    ax.axhline(1.0 / (35**2), color="#111111", ls="--", lw=1.6, label=r"Elimination threshold $1/m^2$ ($m=35$)")
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7.5)
    ax.set_ylabel(r"$\max_{1 \leq t \leq m/2} \Phi(A_{UW}, t/m)$ (log)")
    ax.set_title("A. Zone 1: Decoupled Fourier Decay on Dense Signed Matrices")
    ax.grid(alpha=0.25, axis="y")
    ax.legend(fontsize=8)

    # Panel B: Zone 1 Spectral Tail Ratio vs Rank-1 Control
    ax = axes[0, 1]
    ratios = [max(r["best_tail_ratio"], 0.1) for r in z1]
    ax.bar(x, ratios, color=["#1f77b4" if r["best_tail_ratio"] >= 10 else "#d62728" for r in z1], alpha=0.88)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7.5)
    ax.set_ylabel(r"Best spectral tail ratio $E_{\mathrm{tail}}(r) / \sigma_{r+1}^2$")
    ax.set_title("B. Zone 1: Hanson-Wright Spectral Tail Ratio (D16.1)")
    ax.grid(alpha=0.25, axis="y")

    # Panel C: Zone 2 Certified |S| <= exp(-nu_ind/8) vs 1/m even at n << m!
    ax = axes[1, 0]
    x2 = np.arange(len(z2))
    l2 = [f"m={r['m']}\nn={r['n']} (n/m={r['n_over_m']:.2f})" for r in z2]
    w = 0.28
    ax.bar(x2 - w, [r["certified_S_ub"] for r in z2], width=w, color="#2ca02c", label=r"Certified $\exp(-\nu_{\mathrm{ind}}/8)$")
    ax.bar(x2, [r["target_1_over_m"] for r in z2], width=w, color="#1f77b4", label=r"Threshold $1/m$")
    ax.bar(x2 + w, [r["max_S_emp"] for r in z2], width=w, color="#ff7f0e", label=r"Empirical $\max_t |S_{1/2}(t)|$")
    ax.set_yscale("log")
    ax.set_xticks(x2)
    ax.set_xticklabels(l2, fontsize=8)
    ax.set_ylabel("Fourier modulus (log)")
    ax.set_title(r"C. Zone 2: Large-Coeff Modular Quadratics Killed at $n \ll m$ (D16.2)")
    ax.grid(alpha=0.25, axis="y")
    ax.legend(fontsize=7.8)

    # Panel D: Zone 2 Empirical Pr[P == 0 mod m] vs 1/m
    ax = axes[1, 1]
    ax.bar(x2 - 0.17, [r["pr_zero_emp"] for r in z2], width=0.34, color="#2ca02c", label=r"Empirical $\Pr_{1/2}[P \equiv 0\ (\mathrm{mod}\ m)]$")
    ax.bar(x2 + 0.17, [r["target_1_over_m"] for r in z2], width=0.34, color="#1f77b4", label=r"Uniform baseline $1/m$")
    ax.set_yscale("log")
    ax.set_xticks(x2)
    ax.set_xticklabels(l2, fontsize=8)
    ax.set_ylabel("Zero probability (log)")
    ax.set_title(r"D. Zone 2: Uniform Equidistribution at $n = O(\Delta \ln m) \ll m$")
    ax.grid(alpha=0.25, axis="y")
    ax.legend(fontsize=8)

    fig.suptitle("Theorem D16: Spectral Tail Dichotomy (Zone 1) & Dyadic Coefficient Exclusion (Zone 2)", fontsize=12, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(FIGURE_PATH)
    plt.close(fig)


def main():
    t0 = time.time()
    z1 = benchmark_zone1_dense_signed()
    z2 = benchmark_zone2_dyadic_modular()
    make_figure(z1, z2)
    out = {
        "theorem": "D16",
        "elapsed_sec": round(time.time() - t0, 2),
        "zone1_dense_signed_spectral_tail": z1,
        "zone2_dyadic_modular_exclusion": z2,
    }
    with open(RESULTS_PATH, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print(f"Saved {RESULTS_PATH} and {FIGURE_PATH} in {out['elapsed_sec']}s.")


if __name__ == "__main__":
    main()
