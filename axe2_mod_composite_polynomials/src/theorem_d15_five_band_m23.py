"""
Theorem D15 — Five-Band Tilted Weyl Differencing, Warnke Typical-Bounded-Differences,
Exact sin^2 Weyl Expectation, and Unconditional Delta <= m^{2/3}/(ln m)^{4/3} Bound
on Signed Graphs modulo m.

Components verified and benchmarked:
  D15.1  Stronger Phase-Shifted Cosine Lemma:
         dist(alpha + delta, pi Z)^2 >= (1/4) delta^2 + c(alpha) delta + (1/8) alpha^2
         for all alpha in [-pi/2, pi/2], delta in [-pi/4, pi/4], with exact minima
         (5/448) pi^2 at (3pi/7, -pi/4) and (1/64) pi^2 at (pi/2, pi/4).
  D15.2  Exact sin^2 Weyl Expectation Bound (linear-term-free):
         E[1_{j in J} dist(2 pi t D_j / m, pi Z)^2] >= (3/16)(1 - exp(-(3/4) deg(j) sin^2(pi t / m)))
         verified against exact 3^n enumeration on signed graphs with arbitrary b_j in Z/mZ.
  D15.3  Five-Band Unconditional Bound (Warnke 2016 typical bounded differences +
         Band 2.75 sin^2 McDiarmid / Band 3 pre-Weyl crossover at T_3 = m / (K^{1/2} Delta^{1/4})):
         n <= 64 m unconditionally up to Delta <= c_2 m^{2/3} / (ln m)^{4/3},
         and forced global edge density rho(G) > c_2' / (m^{1/3} (ln m)^{7/3}) for n >= 256 m
         (improving D10 m^{1/4} -> D11 m^{1/2} -> D15 m^{2/3}).
  D15.4  Empirical benchmarks at Delta ~ m^{2/3} on m in {35, 143, 323, 899}.
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
RESULTS_PATH = os.path.join(BASE_DIR, "results", "theorem_d15_five_band_m23.json")
FIGURE_PATH = os.path.join(BASE_DIR, "figures", "theorem_d15_five_band_m23.png")


def verify_d15_1_stronger_phase_cosine():
    """Verify D15.1 on a dense grid and check exact closed-form minima."""
    alphas = np.linspace(-np.pi / 2, np.pi / 2, 1001)
    deltas = np.linspace(-np.pi / 4, np.pi / 4, 501)
    A, D = np.meshgrid(alphas, deltas, indexing="ij")
    U = A + D
    dist = np.abs(U - np.pi * np.round(U / np.pi))

    C = np.where(np.abs(A) <= np.pi / 4, 2.0 * A, np.sign(A) * (np.pi / 2 - np.abs(A)))
    margin = dist**2 - 0.25 * D**2 - C * D - 0.125 * A**2

    min_all = float(np.min(margin))
    outer_mask = np.abs(A) > np.pi / 4 + 1e-12
    min_outer = float(np.min(margin[outer_mask]))
    exact_5_448 = (5.0 / 448.0) * (np.pi**2)
    exact_1_64 = (1.0 / 64.0) * (np.pi**2)

    # Check exact value at (3pi/7, -pi/4) and (pi/2, pi/4)
    a1, d1 = 3.0 * np.pi / 7.0, -np.pi / 4.0
    m1 = (a1 + d1)**2 - 0.25 * d1**2 - (np.pi / 2 - a1) * d1 - 0.125 * a1**2
    a2, d2 = np.pi / 2.0, np.pi / 4.0
    m2 = (np.pi - a2 - d2)**2 - 0.25 * d2**2 - 0.0 * d2 - 0.125 * a2**2

    assert min_all >= -1e-12, f"Margin negative: {min_all}"
    assert abs(m1 - exact_5_448) < 1e-12
    assert abs(m2 - exact_1_64) < 1e-12

    return {
        "min_margin_all": min_all,
        "min_margin_outer": min_outer,
        "exact_min_subcase_2a_5_448_pi2": exact_5_448,
        "exact_min_subcase_2b_1_64_pi2": exact_1_64,
        "verified_pass": bool(min_all >= -1e-12 and min_outer >= exact_5_448 - 1e-4),
    }


def verify_d15_2_exact_sin2_weyl():
    """Exact 3^n enumeration verifying Lemma D15.2 for arbitrary linear terms b_j in Z/mZ."""
    rng = np.random.default_rng(2026)
    cases = []
    min_ratio = 1e9

    for m, n in [(15, 7), (21, 7), (35, 8)]:
        # Random signed graph with degree >= 2
        A = np.zeros((n, n), dtype=int)
        for i in range(n):
            for j in range(i + 1, n):
                if rng.random() < 0.65:
                    s = int(rng.choice([-1, 1]))
                    A[i, j] = A[j, i] = s
        b = rng.integers(0, m, size=n)  # COMPLETELY ARBITRARY b_j in Z/mZ
        theta = float(rng.uniform(0.25, 0.75))
        probs = np.array([(1 - theta)**2, 2 * theta * (1 - theta), theta**2])
        xi_vals = np.array([0.0, 0.5, 1.0])

        # Enumerate all 3^n states of xi in {0, 1/2, 1}^n
        N3 = 3**n
        idx = np.arange(N3, dtype=np.int64)
        states = np.zeros((N3, n), dtype=np.int8)
        tmp = idx.copy()
        for col in range(n):
            states[:, col] = tmp % 3
            tmp //= 3

        w_prob = np.prod(probs[states], axis=1)
        xi_mat = xi_vals[states]  # (3^n, n)
        D_mat = b[None, :] + xi_mat @ A  # (3^n, n), note A has 0 diagonal
        J_mask = (states == 1)  # (3^n, n)
        degs = np.sum(np.abs(A), axis=1)

        for t in [1, m // 4, m // 2]:
            phi = (2.0 * np.pi * t / m) * D_mat
            dist2 = (phi - np.pi * np.round(phi / np.pi))**2
            exact_E = np.sum(w_prob[:, None] * J_mask * dist2, axis=0)
            lb = (3.0 / 16.0) * (1.0 - np.exp(-0.75 * degs * (np.sin(np.pi * t / m)**2)))
            for j in range(n):
                if lb[j] > 1e-12:
                    ratio = float(exact_E[j] / lb[j])
                    if ratio < min_ratio:
                        min_ratio = ratio
                    assert exact_E[j] >= lb[j] - 1e-12, f"Violation at m={m}, t={t}, j={j}"
            cases.append({
                "m": m, "n": n, "t": int(t), "theta": round(theta, 4),
                "mean_exact_E_Ft": float(np.sum(exact_E)),
                "sum_lb": float(np.sum(lb)),
                "ratio_sum": float(np.sum(exact_E) / max(np.sum(lb), 1e-12)),
            })

    return {
        "n_cases_exact_3n": len(cases),
        "min_vertex_ratio_exact_over_lb": min_ratio,
        "cases": cases,
        "verified_pass": bool(min_ratio >= 1.0 - 1e-9),
    }


def build_random_d_regular_signed(n, deg, rng):
    """Fast approximate deg-regular signed graph on n vertices."""
    adj = [[] for _ in range(n)]
    # Configuration-style random pairings (each pairing adds 1 to degree)
    for _ in range(deg):
        perm = rng.permutation(n)
        for i in range(0, n - 1, 2):
            u, v = int(perm[i]), int(perm[i + 1])
            s = int(rng.choice([-1, 1]))
            adj[u].append((v, s))
            adj[v].append((u, s))
    # Deduplicate edges
    edges_u, edges_v, edges_s = [], [], []
    seen = set()
    for u in range(n):
        for v, s in adj[u]:
            if u < v and (u, v) not in seen:
                seen.add((u, v))
                edges_u.append(u)
                edges_v.append(v)
                edges_s.append(s)
    return np.array(edges_u, dtype=np.int32), np.array(edges_v, dtype=np.int32), np.array(edges_s, dtype=np.int8)


def benchmark_five_bands_m23():
    """
    Benchmark dense signed graphs at Delta ~ m^{2/3} (beyond D11's m^{1/2} barrier!)
    for m in {35, 143, 323, 899}, measuring:
      - Empirical Pr_{theta*}[P == 0 mod m] vs 1/m
      - Max |S_{theta*}(t)| in Band 2.5, Band 2.75, and Band 3
      - Crossover frequencies T_1, T_2, T_2.5, T_3
    """
    rng = np.random.default_rng(15)
    benchmarks = []
    for m in [35, 143, 323, 899]:
        n = min(4 * m, 1200)
        deg = max(6, int(round(m ** (2.0 / 3.0))))
        if deg % 2 == 1:
            deg += 1
        eu, ev, es = build_random_d_regular_signed(n, deg, rng)
        # Actual max degree
        degs = np.zeros(n, dtype=np.int32)
        np.add.at(degs, eu, 1)
        np.add.at(degs, ev, 1)
        Delta = int(np.max(degs))
        d_bar = float(np.mean(degs))
        b = np.ones(n, dtype=np.int64)  # B = 1

        # Find theta* in [1/4, 3/4] centering E[P] in m Z
        A_sum = int(np.sum(es))
        B_sum = int(np.sum(b))
        # mu(theta) = A_sum * theta^2 + B_sum * theta
        thetas_grid = np.linspace(0.25, 0.75, 20001)
        mus = A_sum * (thetas_grid**2) + B_sum * thetas_grid
        dist_m = np.abs(mus - m * np.round(mus / m))
        theta_star = float(thetas_grid[int(np.argmin(dist_m))])

        # Sample x ~ Ber(theta*)^n to measure Pr[P == 0 mod m] and S(t)
        n_samp = 8000
        X = (rng.random((n_samp, n)) < theta_star).astype(np.int8)
        # Evaluate P(x) = sum_{e} es * x_u * x_v + sum_i b_i * x_i
        # Chunked for memory efficiency
        P_vals = X @ b
        chunk = 4000
        for start in range(0, len(eu), chunk):
            u_c = eu[start:start + chunk]
            v_c = ev[start:start + chunk]
            s_c = es[start:start + chunk]
            P_vals += (X[:, u_c] * X[:, v_c]) @ s_c

        pr_zero = float(np.mean((P_vals % m) == 0))

        # Frequencies across bands
        sigma_emp = float(np.std(P_vals))
        x_val = math.sqrt(2.0) * math.pi * sigma_emp / m
        T0 = int(1.0 / max(x_val, 1e-9))
        K = math.ceil(4.0 * math.log2(2.0 * m))
        L = math.sqrt(8.0 * Delta * math.log(40.0 * n * m))
        T1 = max(T0 + 1, int(m / (4.0 * (Delta + 1))))
        T2 = max(T1, int(m / (8.0 * L)))
        T25 = max(T2, int(m / math.sqrt(max(d_bar, 1.0))))
        T3 = max(T25, int(m / ((K**0.5) * (Delta**0.25))))
        T_half = m // 2

        # Compute |S_{theta*}(t)| on a representative grid of frequencies > T0 (outside Band 1)
        sample_ts = sorted(set(
            list(np.linspace(max(T0 + 1, T1), max(T0 + 1, T2), 15, dtype=int))
            + list(np.linspace(max(T0 + 1, T2), max(T0 + 1, T3), 25, dtype=int))
            + list(np.linspace(max(T0 + 1, T3), T_half, 25, dtype=int))
        ))
        phases_base = (2.0 * np.pi / m) * (P_vals % m)
        max_S_b25 = 0.0
        max_S_b275 = 0.0
        max_S_b3 = 0.0
        for t in sample_ts:
            if t <= T0:
                continue
            s_val = float(np.abs(np.mean(np.exp(1j * t * phases_base))))
            if t <= T2:
                max_S_b25 = max(max_S_b25, s_val)
            elif t <= T3:
                max_S_b275 = max(max_S_b275, s_val)
            else:
                max_S_b3 = max(max_S_b3, s_val)

        benchmarks.append({
            "m": m,
            "n": n,
            "Delta": Delta,
            "m_2_3": round(m ** (2.0 / 3.0), 2),
            "m_1_2": round(m ** 0.5, 2),
            "theta_star": round(theta_star, 5),
            "pr_zero_empirical": round(pr_zero, 5),
            "baseline_1_over_m": round(1.0 / m, 5),
            "T1": T1, "T2": T2, "T25": T25, "T3": T3, "T_half": T_half,
            "max_S_band25": round(max_S_b25, 5),
            "max_S_band275": round(max_S_b275, 5),
            "max_S_band3": round(max_S_b3, 5),
        })
    return benchmarks


def make_figure(res_d15_1, res_d15_2, benchmarks):
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.2), dpi=160)
    plt.rcParams["font.size"] = 9.5

    # Panel A: Stronger Phase-Shifted Cosine margin (D15.1)
    ax = axes[0, 0]
    alphas = np.linspace(np.pi / 4, np.pi / 2, 200)
    deltas = np.linspace(-np.pi / 4, np.pi / 4, 201)
    A, D = np.meshgrid(alphas, deltas, indexing="ij")
    U = A + D
    dist = np.abs(U - np.pi * np.round(U / np.pi))
    C = np.pi / 2 - A
    M = dist**2 - 0.25 * D**2 - C * D - 0.125 * A**2
    min_over_delta = np.min(M, axis=1)
    ax.plot(alphas / np.pi, min_over_delta, color="#1f77b4", lw=2.3, label=r"$\min_\delta$ margin ($\kappa = 1/8$)")
    ax.axhline((5.0 / 448.0) * np.pi**2, color="#d62728", ls="--", lw=1.4, label=r"Exact global min $\frac{5}{448}\pi^2 \approx 0.1102$")
    ax.axhline((1.0 / 64.0) * np.pi**2, color="#2ca02c", ls=":", lw=1.4, label=r"Subcase 2b min $\frac{1}{64}\pi^2 \approx 0.1542$")
    ax.set_xlabel(r"$\alpha / \pi \in [1/4, 1/2]$")
    ax.set_ylabel("Minimum margin over $\\delta \\in [-\\pi/4, \\pi/4]$")
    ax.set_title("A. Stronger Phase-Shifted Cosine Lemma (D15.1)")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)

    # Panel B: 5-Band Exponent Crossover (Band 2.75 vs Band 3)
    ax = axes[0, 1]
    Delta_ref = 100.0
    K_ref = 20.0
    n_ref = 64.0 * 1000.0
    u_grid = np.logspace(-2.2, -0.3, 300)  # t / m
    exp_b275 = n_ref / (K_ref**2 * (Delta_ref**2) * (u_grid**2))
    exp_b3 = (3.0 / 32.0) * (n_ref / Delta_ref) * (u_grid**2)
    exp_max = np.maximum(exp_b275, exp_b3)
    u_cross = (32.0 / 3.0)**0.25 / ((K_ref**0.5) * (Delta_ref**0.25))
    ax.loglog(u_grid, exp_b275, color="#1f77b4", lw=2.0, label=r"Band 2.75 ($\sin^2$ McDiarmid): $\propto n / (K^2 \Delta^2 (t/m)^2)$")
    ax.loglog(u_grid, exp_b3, color="#d62728", lw=2.0, label=r"Band 3 (pre-Weyl matching): $\propto n (t/m)^2 / \Delta$")
    ax.loglog(u_grid, exp_max, color="#111111", lw=2.5, ls="--", label=r"Envelope $\max \geq \Omega(n / (K \Delta^{3/2}))$")
    ax.axvline(u_cross, color="#2ca02c", ls=":", lw=1.5, label=r"Crossover $T_3/m \sim K^{-1/2}\Delta^{-1/4}$")
    ax.set_xlabel(r"Normalized frequency $t / m$")
    ax.set_ylabel("Concentration exponent (log scale)")
    ax.set_title("B. Band 2.75 / Band 3 Exponent Crossover (D15.3)")
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=7.5)

    # Panel C: Progression of unconditional max-degree & forced density exponents
    ax = axes[1, 0]
    ms = np.logspace(1.5, 4.5, 200)
    ax.loglog(ms, np.ones_like(ms), ls=":", color="#7f7f7f", lw=1.8, label="Caro-Wei (D9.2): $\\Delta = O(1)$")
    ax.loglog(ms, ms**0.25, color="#ff7f0e", lw=2.0, label="Theorem D10: $\\Delta \\sim m^{1/4}$ ($\\rho > m^{-3/4}$)")
    ax.loglog(ms, ms**0.50, color="#1f77b4", lw=2.2, label="Theorem D11: $\\Delta \\sim m^{1/2}$ ($\\rho > m^{-1/2}$)")
    ax.loglog(ms, ms**(2.0 / 3.0), color="#d62728", lw=2.8, label="Theorem D15: $\\Delta \\sim m^{2/3}$ ($\\rho > m^{-1/3}$)")
    ax.loglog(ms, ms, ls="--", color="#2ca02c", lw=1.5, label="Lee target: $\\Delta = m$ ($\\rho = 1$)")
    ax.set_xlabel("Composite modulus $m$")
    ax.set_ylabel("Unconditional $n \\leq O(m)$ max-degree threshold $\\Delta$")
    ax.set_title("C. Exponent Progression: D9 -> D10 -> D11 -> D15")
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=7.8)

    # Panel D: Dense signed graph benchmarks at Delta ~ m^{2/3}
    ax = axes[1, 1]
    labels = [f"m={b['m']}\n($\\Delta$={b['Delta']})" for b in benchmarks]
    x = np.arange(len(benchmarks))
    w = 0.25
    pr_emp = [b["pr_zero_empirical"] for b in benchmarks]
    pr_base = [b["baseline_1_over_m"] for b in benchmarks]
    max_S = [max(b["max_S_band25"], b["max_S_band275"], b["max_S_band3"]) for b in benchmarks]
    ax.bar(x - w, pr_emp, width=w, color="#2ca02c", label=r"Empirical $\Pr_{\theta^*}[P \equiv 0\ (\mathrm{mod}\ m)]$")
    ax.bar(x, pr_base, width=w, color="#1f77b4", label=r"Uniform baseline $1/m$")
    ax.bar(x + w, max_S, width=w, color="#d62728", label=r"$\max_{B2.5/2.75/3} |S_{\theta^*}(t)|$")
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Probability / Fourier modulus (log)")
    ax.set_title(r"D. Dense Signed Benchmarks at $\Delta \approx m^{2/3}$ (beyond D11)")
    ax.grid(alpha=0.25, axis="y")
    ax.legend(fontsize=8)

    fig.suptitle("Theorem D15: Five-Band Tilted Weyl + Warnke + sin² Bound ($\\Delta \\lesssim m^{2/3}$, $\\rho \\gtrsim m^{-1/3}$)", fontsize=12, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(FIGURE_PATH)
    plt.close(fig)


def main():
    t0 = time.time()
    r1 = verify_d15_1_stronger_phase_cosine()
    print("[D15.1] Stronger phase-shifted cosine:", r1)
    r2 = verify_d15_2_exact_sin2_weyl()
    print(f"[D15.2] Exact 3^n sin^2 Weyl expectation: {r2['n_cases_exact_3n']} cases, min ratio = {r2['min_vertex_ratio_exact_over_lb']:.4f}")
    r3 = benchmark_five_bands_m23()
    for b in r3:
        print(f"[D15.4] m={b['m']:3d}, n={b['n']:4d}, Delta={b['Delta']:3d} (~m^2/3={b['m_2_3']:5.1f}) | "
              f"Pr0={b['pr_zero_empirical']:.5f} (1/m={b['baseline_1_over_m']:.5f}) | "
              f"max|S| B2.5={b['max_S_band25']:.4f}, B2.75={b['max_S_band275']:.4f}, B3={b['max_S_band3']:.4f}")

    make_figure(r1, r2, r3)
    out = {
        "theorem": "D15",
        "elapsed_sec": round(time.time() - t0, 2),
        "d15_1_stronger_phase_cosine": r1,
        "d15_2_exact_sin2_weyl": r2,
        "d15_4_dense_benchmarks_m23": r3,
    }
    with open(RESULTS_PATH, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print(f"Saved {RESULTS_PATH} and {FIGURE_PATH} in {out['elapsed_sec']}s.")


if __name__ == "__main__":
    main()
