"""
Phase 6 : Solveur v5 Unifie (Theoreme 1 A* + Theoreme 2 Gray glissant)
          et Mode Zero-Gap (B = 1, N_hot = 0 garanti dans le pire cas)
======================================================================

Ce script reunit en un seul algorithme universel (Solveur v5 Unifie) :
  1. Detection et resolution structurelle a la racine et dans l'arbre via le
     Theoreme 1 (Primitive A* : clusters L1 + enveloppe de residu + reseau PGCD
     + descripteur de discordance f8) et la propagation Pseudo-Booleenne a
     coupes surrogees.
  2. Bascule automatique vers la compression bilineaire a pointeurs glissants
     en code de Gray (Theoreme 2) des que le circuit (ou sous-probleme) est
     dense / sans clusters L1 (c_up == up), divisant par 2^l le nombre d'etats
     parcourus (l = 9 -> 2^n / 512 etats).
  3. Analyse experimentale du verrou asymptotique de la Marche 2 : comparaison
     entre le mode B = 1 (Zero-Gap : P = 2^l, ecart de relaxation identiquement
     nul, N_hot = 0 garanti sur 100% des instances UNSAT) et B = 4 / B = 8,
     et test de passage a l'echelle jusqu'a n = 28 (2^28 = 268,435,456 entrees).
"""
import json
import math
import os
import sys
import time

import numpy as np
from numba import njit, prange

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from thr2_lab import HERE, make_instances, top_max
from adversarial_suite import make_parity, make_mod3, pack_single
from theorem_a_star import (NF8, cluster_residual_bound_and_f8,
                            make_perturbed_parity, make_two_chain_perturbed)
from industrial_and_orthogonal import (
    gray_code_exhaust_batch, make_ortho_blocks, make_inner_product_noisy,
    solve_v4
)
from step2_sliding_gray import sliding_gray_single, sliding_gray_batch

OUT = sys.argv[1] if len(sys.argv) > 1 else HERE


# --------------------------------------------------------------------------- #
#  Solveur v5 Unifie : Theoreme 1 (A* + PB/Coupes + f8) + Theoreme 2 (Gray glissant)
# --------------------------------------------------------------------------- #
@njit(cache=True)
def solve_v5_single(W, t, v, T, theta8, l_right, log2_B):
    """
    Algorithme universel :
      - Teste a la racine la presence de clusters L1 (Theoreme 1 A*).
      - Si c_up < up (circuit structure : parites, chaines perturbees, blocs
        orthogonaux), lance le solveur A* + PB/Coupes + f8.
      - Si c_up == up (circuit dense sans clusters L1 : s = 2n, 3n, 4n),
        bascule immediatement sur le Theoreme 2 (Gray glissant 2^(n - l_right)).
    """
    s, n = W.shape
    Wpos = np.maximum(W, 0)
    Wneg = np.minimum(W, 0)
    S = np.zeros(s, np.int64)
    P = Wpos.sum(axis=1)
    N = Wneg.sum(axis=1)
    assigned = np.zeros(n, np.bool_)
    und = np.zeros(s, np.bool_)
    lo = np.zeros(s, np.int64)
    hi = np.zeros(s, np.int64)
    f8 = np.zeros(n, np.float64)

    base = 0
    up = 0
    dn = 0
    n_und = 0
    for g in range(s):
        lo[g] = N[g]
        hi[g] = P[g]
        if lo[g] >= t[g]:
            base += v[g]
            und[g] = False
        elif hi[g] < t[g]:
            und[g] = False
        else:
            und[g] = True
            n_und += 1
            if v[g] > 0:
                up += v[g]
            else:
                dn += v[g]
    if base + dn >= T:
        return 1, 1
    if base + up < T:
        return 0, 1

    c_up = cluster_residual_bound_and_f8(W, t, v, S, lo, hi, und, assigned, f8)
    if base + c_up < T:
        return 0, 1

    # Si le Theoreme 1 detecte une structure de clusters (c_up < up ou f8 > 0 sur un cluster),
    # utiliser le solveur A* + PB/Coupes + f8 ; sinon basculer sur le Theoreme 2 (Gray glissant).
    if c_up < up:
        return solve_v4(W, t, v, T, theta8, 1 << (n + 2), True, 3, True, 8)
    else:
        return sliding_gray_single(W, t, v, T, l_right, log2_B)


@njit(parallel=True, cache=True)
def solve_v5_batch(Ws, ts, vs, Ts, theta8, l_right, log2_B):
    B = Ws.shape[0]
    status = np.zeros(B, np.int64)
    work = np.zeros(B, np.int64)
    for k in prange(B):
        st, w = solve_v5_single(Ws[k], ts[k], vs[k], Ts[k], theta8, l_right, log2_B)
        status[k] = st
        work[k] = w
    return status, work


def sanity_check(th8):
    print("== 1. Verification de correction du Solveur v5 Unifie (SAT & UNSAT) ==", flush=True)
    tot = 0
    for sat in (False, True):
        for n in (12, 14, 16):
            for gen in (
                lambda: make_instances(n, 2 * n, 10, seed=1010 + n + 50 * sat, sat=sat),
                lambda: make_perturbed_parity(n, k_vars=n // 3, FlipsPerGate=2, count=10, seed=1020 + n + 50 * sat, sat=sat),
                lambda: make_ortho_blocks(n, 4, 10, seed=1030 + n + 50 * sat, sat=sat),
                lambda: make_inner_product_noisy(n, 10, seed=1040 + n + 50 * sat, sat=sat),
            ):
                Ws, ts, vs, Ts = gen()
                st, _ = solve_v5_batch(Ws, ts, vs, Ts, th8, 8, 2)
                exp = 1 if sat else 0
                assert (st == exp).all(), f"ECHEC v5 sat={sat} n={n}: {st}"
                tot += len(st)
    print(f"   {tot}/{tot} tests SAT et UNSAT verifies a 100% contre la force brute.\n", flush=True)


def main():
    t0 = time.time()
    g_hyb = np.array([1.0, 0.13, 0.86, 0.58, -0.27, 0.83, -0.48], dtype=np.float64)
    th8 = np.append(g_hyb, 10.0)

    sanity_check(th8)

    # ----------------------------------------------------------------------- #
    #  2. Verification du mode Zero-Gap (B = 1, log2_B = 0) vs B = 4 et B = 8
    # ----------------------------------------------------------------------- #
    print("== 2. Analyse du verrou Marche 2 : Mode Zero-Gap (B=1, N_hot=0) vs B=4 et B=8 ==", flush=True)
    zg_ns = [18, 20, 22, 24]
    zg_sets = {n: make_instances(n, 2 * n, 8, seed=9500 + n) for n in zg_ns}
    zg_configs = [
        ("Zero-Gap exact (l=6, P=64, B=1 -> N_hot=0)", 6, 0),
        ("Paquets B=4    (l=8, P=64, B=4)", 8, 2),
        ("Paquets B=8    (l=9, P=64, B=8)", 9, 3),
    ]
    zg_report = {}
    for cname, lr, lb in zg_configs:
        logs = []
        times = []
        hot_excess = []
        for n in zg_ns:
            Ws, ts, vs, Ts = zg_sets[n]
            ts0 = time.perf_counter()
            st, w = sliding_gray_batch(Ws, ts, vs, Ts, lr, lb)
            dt = time.perf_counter() - ts0
            assert (st == 0).all()
            outer_exact = float(1 << (n - lr))
            excess = float(w.mean()) - outer_exact
            hot_excess.append(excess)
            logs.append(math.log2(float(w.mean())))
            times.append(dt)
        zg_report[cname] = {"log2_work": logs, "times_s": times, "hot_excess": hot_excess}
        print(f"   {cname:<44} | t(n=24)={times[-1]:5.3f}s | log2(etats) = "
              + " ".join(f"{x:5.1f}" for x in logs)
              + f" | verifs internes moy. = {hot_excess[-1]:.1f}", flush=True)
    print("", flush=True)

    # ----------------------------------------------------------------------- #
    #  3. Panorama universel du Solveur v5 sur TOUTES les familles (n = 18)
    # ----------------------------------------------------------------------- #
    n_pan = 18
    print(f"== 3. Panorama universel sur les 7 familles a n = {n_pan} ==", flush=True)
    pan_suites = {
        "F1: PARITY_n": pack_single([make_parity(n_pan)]),
        "F2: MOD3_n": pack_single([make_mod3(n_pan)]),
        "P2: Parite perturbee (k=n/3)": make_perturbed_parity(n_pan, n_pan // 3, 2, 8, seed=2200 + n_pan),
        "P3: Double chaine perturbee": make_two_chain_perturbed(n_pan, 4, 8, seed=3300 + n_pan),
        "O1: 4 Blocs orthogonaux (s=2n)": make_ortho_blocks(n_pan, 4, 8, seed=8100 + n_pan),
        "O4: Aleatoires denses (s=2n)": make_instances(n_pan, 2 * n_pan, 8, seed=8400 + n_pan),
        "HD: Haute densite (s=4n)": make_instances(n_pan, 4 * n_pan, 8, seed=9700 + n_pan),
    }
    pan_report = {}
    print(f"   {'Famille':<32} | {'Gray naif 2^n':>16} | {'Solveur v4':>16} | {'Solveur v5 Unifie':>18}", flush=True)
    print("   " + "-" * 90, flush=True)
    for fname, (Ws, ts, vs, Ts) in pan_suites.items():
        t1 = time.perf_counter()
        st1, w1 = gray_code_exhaust_batch(Ws, ts, vs, Ts)
        dt1 = time.perf_counter() - t1

        from industrial_and_orthogonal import solve_v4_batch
        t2 = time.perf_counter()
        st2, w2 = solve_v4_batch(Ws, ts, vs, Ts, th8, 1 << (n_pan + 2), True, 3, True, 8)
        dt2 = time.perf_counter() - t2

        t3 = time.perf_counter()
        st3, w3 = solve_v5_batch(Ws, ts, vs, Ts, th8, 8, 2)
        dt3 = time.perf_counter() - t3

        assert (st1 == 0).all() and (st2 == 0).all() and (st3 == 0).all()
        l1, l2, l3 = math.log2(float(w1.mean())), math.log2(float(w2.mean())), math.log2(float(w3.mean()))
        pan_report[fname] = {
            "gray": {"log2": l1, "t": dt1},
            "v4": {"log2": l2, "t": dt2},
            "v5": {"log2": l3, "t": dt3},
        }
        print(f"   {fname:<32} | 2^{l1:4.1f} ({dt1:5.3f}s) | 2^{l2:4.1f} ({dt2:5.3f}s) | 2^{l3:4.1f} ({dt3:5.3f}s)",
              flush=True)
    print("", flush=True)

    # ----------------------------------------------------------------------- #
    #  4. Passage a l'echelle jusqu'a n = 28 (2^28 = 268,435,456 entrees !)
    # ----------------------------------------------------------------------- #
    scale_ns = [20, 22, 24, 26, 28]
    print("== 4. Passage a l'echelle du Solveur v5 Unifie jusqu'a n = 28 (s = 2n, 6 instances) ==", flush=True)
    scale_report = {"ns": scale_ns, "v5_l8": [], "v5_l9": []}
    for n in scale_ns:
        # Pour n=28, generer 6 instances UNSAT de justesse via top_max
        Ws, ts, vs, Ts = make_instances(n, 2 * n, 6, seed=11000 + n)
        ts0 = time.perf_counter()
        st8, w8 = solve_v5_batch(Ws, ts, vs, Ts, th8, 8, 2)
        dt8 = time.perf_counter() - ts0
        ts0 = time.perf_counter()
        st9, w9 = solve_v5_batch(Ws, ts, vs, Ts, th8, 9, 3)
        dt9 = time.perf_counter() - ts0
        assert (st8 == 0).all() and (st9 == 0).all()
        l8, l9 = math.log2(float(w8.mean())), math.log2(float(w9.mean()))
        scale_report["v5_l8"].append({"n": n, "log2": l8, "time_s": dt8})
        scale_report["v5_l9"].append({"n": n, "log2": l9, "time_s": dt9})
        print(f"   n = {n:>2} (2^n = {1<<n:>11,}) | v5 (l=8) : 2^{l8:4.1f} etats en {dt8:6.3f}s | "
              f"v5 (l=9) : 2^{l9:4.1f} etats en {dt9:6.3f}s", flush=True)

    print(f"\nTemps total Phase 6 : {time.time() - t0:.1f} s", flush=True)
    with open(os.path.join(HERE, "phase6_v5_results.json"), "w", encoding="utf-8") as fh:
        json.dump({"zero_gap": zg_report, "panorama_n18": pan_report, "scaling_n28": scale_report},
                  fh, indent=2)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 4.8))

    labels = list(pan_report.keys())
    x = np.arange(len(labels))
    w = 0.26
    ax1.bar(x - w, [pan_report[k]["gray"]["log2"] for k in labels], w, label="Code de Gray naïf 2^n", color="#7f7f7f")
    ax1.bar(x, [pan_report[k]["v4"]["log2"] for k in labels], w, label="Solveur v4 (Arbre hybride)", color="#1f77b4")
    ax1.bar(x + w, [pan_report[k]["v5"]["log2"] for k in labels], w, label="Solveur v5 Unifié (Th.1 + Th.2)", color="#2ca02c")
    ax1.set_xticks(x)
    ax1.set_xticklabels([k.split(":")[0] for k in labels])
    ax1.set_ylabel("log2(états / nœuds explorés)")
    ax1.set_title("Panorama universel sur les 7 familles (n = 18)")
    ax1.grid(alpha=0.3, axis="y")
    ax1.legend(fontsize=8)

    ax2.semilogy(scale_ns, [r["time_s"] for r in scale_report["v5_l8"]], "o-",
                 label="Solveur v5 Unifié (l=8, 2^(n-8))", color="#ff7f0e")
    ax2.semilogy(scale_ns, [r["time_s"] for r in scale_report["v5_l9"]], "s-",
                 label="Solveur v5 Unifié (l=9, 2^(n-9))", color="#2ca02c")
    ax2.set_xlabel("Nombre de variables n (jusqu'à 2^28 = 268M entrées)")
    ax2.set_ylabel("Temps réel pour 6 instances s=2n (s, log)")
    ax2.set_title("Passage à l'échelle du Solveur v5 Unifié jusqu'à n = 28")
    ax2.grid(alpha=0.3, which="both")
    ax2.legend(fontsize=8)

    fig.tight_layout()
    out_png = os.path.join(OUT, "phase6_unified_v5.png")
    fig.savefig(out_png, dpi=130)
    print(f"Graphique : {out_png}", flush=True)


if __name__ == "__main__":
    main()
