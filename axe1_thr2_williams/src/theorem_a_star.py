"""
Theoreme 1 (Primitive A*) : Borne par projection 1D + enveloppe de residu
=========================================================================

Pour tout groupe G de portes indecises, tout vecteur de reference w in Z^n
sur les variables libres F, et toute assignation x in {0,1}^n :

    sum_{k in G} v_k g_k(x)  <=  max_{L in [L_min, L_max]} sum_{k in G} v_k 1[L >= tau*_k]

ou L(x) = sum_{i in F} w_i x_i in [L_min, L_max] et :
  - si v_k > 0 : tau*_k = t_k - S_k - sum_{i in F} max(0, W[k,i] - w_i)
  - si v_k < 0 : tau*_k = t_k - S_k - sum_{i in F} min(0, W[k,i] - w_i)

Corollaire (pire cas) : si un circuit THR∘THR est forme de m chaines dont les
perturbations portent sur au plus k variables chacune (ou si chaque porte a un
residu L1 borne), le branchement prioritaire sur les variables de discordance
(nouveau descripteur f8) ramene le residu a 0 et resout le circuit en au plus
O(2^(m*k)) noeuds dans le pire cas, la ou :
  - l'elagage par intervalles pur explore ~ 2^n noeuds (alpha = 1),
  - les chaines 1D exactes (Primitive A) sont aveugles tant que le residu != 0,
  - la relaxation LP continue (Primitive B) souffre du saut d'integralite.

Ce script :
  1. Verifie l'exactitude a 100% du solveur A* + B + f8 sur 200+ instances SAT/UNSAT.
  2. Teste le solveur sur 4 familles adversariales qui cassent l'Hybride A+B :
       - P1 : PARITE PERTURBEE (chaque porte perturbee sur 1 variable)
       - P2 : PARITE PERTURBEE FORTE (chaque porte perturbee sur 2 variables)
       - P3 : DOUBLE CHAINE PERTURBEE (2 parites perturbees sur sous-ensembles)
       - P4 : CIRCUITS ALEATOIRES DENSES (s = 2n)
"""
from pathlib import Path
import json
import math
import os
import sys
import time

import numpy as np
from numba import njit, prange

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from thr2_lab import HERE, make_instances, top_max
from hybrid_thr2 import lagrangian_upper_bound, solve_hybrid_batch

AXE_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = str(AXE_ROOT / "results")
FIGURES_DIR = str(AXE_ROOT / "figures")
WITNESS_DIR = str(AXE_ROOT / "witnesses")
for _d in (RESULTS_DIR, FIGURES_DIR, WITNESS_DIR):
    os.makedirs(_d, exist_ok=True)
OUT = sys.argv[1] if len(sys.argv) > 1 else FIGURES_DIR
NF8 = 8  # 7 descripteurs originaux + f8 (discordance intra-cluster)


# --------------------------------------------------------------------------- #
#  Theoreme 1 (Primitive A*) : Clustering + borne 1D avec residu exact
#  + calcul simultane du descripteur f8 (discordance par variable libre)
# --------------------------------------------------------------------------- #
@njit(cache=True)
def cluster_residual_bound_and_f8(W, t, v, S, lo, hi, und, assigned, f8_out):
    """
    Partitionne les portes indecises en clusters de vecteurs proches (au sens L1),
    calcule le vecteur consensus w de chaque cluster, applique le Theoreme 1
    pour obtenir une borne superieure 100% valide, et remplit f8_out[i] avec
    la discordance totale de la variable libre i par rapport aux consensus.
    """
    s, n = W.shape
    for i in range(n):
        f8_out[i] = 0.0

    used = np.zeros(s, np.bool_)
    members = np.zeros(s, np.int64)
    w = np.zeros(n, np.int64)
    taus = np.zeros(s, np.int64)
    vals = np.zeros(s, np.int64)
    total_up = 0

    for g in range(s):
        if not und[g] or used[g]:
            continue
        # Norme L1 de g sur les variables libres
        norm_g = 0
        for i in range(n):
            if not assigned[i]:
                norm_g += abs(W[g, i])

        cnt = 0
        members[cnt] = g
        cnt += 1
        used[g] = True

        # Regrouper toutes les portes indecises h a distance L1 <= 0.45 * norm_g
        # (ou colineaires de meme signe)
        thresh_dist = (norm_g * 9) // 20
        for h in range(g + 1, s):
            if not und[h] or used[h]:
                continue
            dist = 0
            for i in range(n):
                if not assigned[i]:
                    dist += abs(W[g, i] - W[h, i])
                    if dist > thresh_dist:
                        break
            if dist <= thresh_dist:
                members[cnt] = h
                cnt += 1
                used[h] = True

        if cnt == 1:
            if v[g] > 0:
                total_up += v[g]
            continue

        # Vecteur consensus w_i = moyenne arrondie sur le cluster + PGCD des w_i libres
        L_min = 0
        L_max = 0
        g_gcd = 0
        for i in range(n):
            if assigned[i]:
                w[i] = 0
                continue
            acc = 0
            for idx in range(cnt):
                acc += W[members[idx], i]
            # Arrondi au plus proche entier
            if acc >= 0:
                wi = (acc + cnt // 2) // cnt
            else:
                wi = -((-acc + cnt // 2) // cnt)
            w[i] = wi
            if wi > 0:
                L_max += wi
            else:
                L_min += wi
            awi = abs(wi)
            if awi > 0:
                if g_gcd == 0:
                    g_gcd = awi
                else:
                    a_g = g_gcd
                    b_g = awi
                    while b_g != 0:
                        a_g, b_g = b_g, a_g % b_g
                    g_gcd = a_g

        # Calcul des seuils decales tau*_h (Theoreme 1 renforce par réseau PGCD) et f8
        pos_sum = 0
        for idx in range(cnt):
            h = members[idx]
            vh = v[h]
            if vh > 0:
                pos_sum += vh
            d_max = 0
            d_min = 0
            av = float(abs(vh))
            for i in range(n):
                if assigned[i]:
                    continue
                diff = W[h, i] - w[i]
                if diff > 0:
                    d_max += diff
                    f8_out[i] += av * diff
                elif diff < 0:
                    d_min += diff
                    f8_out[i] += av * (-diff)
            if vh > 0:
                tau = t[h] - S[h] - d_max
            else:
                tau = t[h] - S[h] - d_min
            if g_gcd > 1:
                if tau >= 0:
                    tau = ((tau + g_gcd - 1) // g_gcd) * g_gcd
                else:
                    tau = -((-tau) // g_gcd) * g_gcd
            taus[idx] = tau
            vals[idx] = vh

        # Tri par insertion de (taus[:cnt], vals[:cnt])
        for a in range(1, cnt):
            kt = taus[a]
            kv = vals[a]
            b = a - 1
            while b >= 0 and taus[b] > kt:
                taus[b + 1] = taus[b]
                vals[b + 1] = vals[b]
                b -= 1
            taus[b + 1] = kt
            vals[b + 1] = kv

        # Maximum exact de sum_{idx} vals[idx] * 1[L >= taus[idx]] sur L in [L_min, L_max]
        cur = 0
        idx = 0
        while idx < cnt and taus[idx] <= L_min:
            cur += vals[idx]
            idx += 1
        best = cur
        while idx < cnt and taus[idx] <= L_max:
            tau_cur = taus[idx]
            while idx < cnt and taus[idx] == tau_cur:
                cur += vals[idx]
                idx += 1
            if cur > best:
                best = cur
        # Toujours plafonner par la borne triviale pos_sum
        if best > pos_sum:
            best = pos_sum
        if best < 0:
            best = 0
        total_up += best

    return total_up


# --------------------------------------------------------------------------- #
#  Solveur A* + B + f8
# --------------------------------------------------------------------------- #
@njit(cache=True)
def _assign(W, Wpos, Wneg, S, P, N, i, b, sign):
    s = W.shape[0]
    for g in range(s):
        S[g] += sign * b * W[g, i]
        P[g] -= sign * Wpos[g, i]
        N[g] -= sign * Wneg[g, i]


@njit(cache=True)
def solve_astar(W, t, v, T, theta8, node_limit, lp_steps):
    s, n = W.shape
    Wpos = np.maximum(W, 0)
    Wneg = np.minimum(W, 0)
    S = np.zeros(s, np.int64)
    P = Wpos.sum(axis=1)
    N = Wneg.sum(axis=1)
    assigned = np.zeros(n, np.bool_)
    svar = np.zeros(n + 1, np.int64)
    sval = np.zeros(n + 1, np.int64)
    sstate = np.zeros(n + 1, np.int64)
    und = np.zeros(s, np.bool_)
    lo = np.zeros(s, np.int64)
    hi = np.zeros(s, np.int64)
    feats = np.zeros((n, NF8))
    fmax = np.zeros(NF8)
    f8 = np.zeros(n, np.float64)
    depth = 0
    nodes = 0

    while True:
        nodes += 1
        if nodes > node_limit:
            return -1, nodes
        base = 0
        up = 0
        dn = 0
        n_und = 0
        for g in range(s):
            lo[g] = S[g] + N[g]
            hi[g] = S[g] + P[g]
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
            return 1, nodes

        pruned = (base + up < T)

        # Theoreme 1 (Primitive A*) : clusters + residus + descripteur f8
        if not pruned and n_und >= 2:
            c_up = cluster_residual_bound_and_f8(W, t, v, S, lo, hi, und, assigned, f8)
            if base + c_up < T:
                pruned = True
        else:
            for i in range(n):
                f8[i] = 0.0

        # Primitive B : borne lagrangienne couplee
        if not pruned and lp_steps > 0 and n_und >= 2:
            ub = lagrangian_upper_bound(W, t, v, S, lo, hi, und, assigned, base, T, lp_steps)
            if ub < T - 1e-7:
                pruned = True

        if not pruned:
            for f in range(NF8):
                fmax[f] = 1e-12
            for i in range(n):
                if assigned[i]:
                    continue
                f1 = 0.0
                f2 = 0.0
                f3 = 0.0
                f4 = 0.0
                f5 = 0.0
                f6 = 0.0
                f7 = 0.0
                for g in range(s):
                    if not und[g]:
                        continue
                    a = abs(W[g, i])
                    if a == 0:
                        continue
                    r = hi[g] - lo[g]
                    d = min(hi[g] - t[g] + 1, t[g] - lo[g])
                    av = abs(v[g])
                    f1 += a
                    f2 += av * a
                    f3 += a / r
                    f4 += av * a / r
                    if a >= d:
                        f5 += 1.0
                        f6 += av
                    f7 += a / d
                feats[i, 0] = f1
                feats[i, 1] = f2
                feats[i, 2] = f3
                feats[i, 3] = f4
                feats[i, 4] = f5
                feats[i, 5] = f6
                feats[i, 6] = f7
                feats[i, 7] = f8[i]
                for f in range(NF8):
                    if feats[i, f] > fmax[f]:
                        fmax[f] = feats[i, f]
            best = -1
            bs = -1e300
            for i in range(n):
                if assigned[i]:
                    continue
                sc = 0.0
                for f in range(NF8):
                    sc += theta8[f] * feats[i, f] / fmax[f]
                if sc > bs:
                    bs = sc
                    best = i
            pref = 0
            for g in range(s):
                if und[g]:
                    pref += v[g] * W[g, best]
            b = 1 if pref >= 0 else 0
            svar[depth] = best
            sval[depth] = b
            sstate[depth] = 0
            assigned[best] = True
            _assign(W, Wpos, Wneg, S, P, N, best, b, 1)
            depth += 1
            continue

        while True:
            if depth == 0:
                return 0, nodes
            d = depth - 1
            i = svar[d]
            b = sval[d]
            _assign(W, Wpos, Wneg, S, P, N, i, b, -1)
            if sstate[d] == 0:
                sstate[d] = 1
                b = 1 - b
                sval[d] = b
                _assign(W, Wpos, Wneg, S, P, N, i, b, 1)
                break
            assigned[i] = False
            depth -= 1


@njit(parallel=True, cache=True)
def solve_astar_batch(Ws, ts, vs, Ts, theta8, node_limit, lp_steps):
    B = Ws.shape[0]
    status = np.zeros(B, np.int64)
    nodes = np.zeros(B, np.int64)
    for k in prange(B):
        st, nd = solve_astar(Ws[k], ts[k], vs[k], Ts[k], theta8, node_limit, lp_steps)
        status[k] = st
        nodes[k] = nd
    return status, nodes


# --------------------------------------------------------------------------- #
#  Generateurs d'instances adversariales "briseuses d'Hybride A+B"
# --------------------------------------------------------------------------- #
def make_perturbed_parity(n, k_vars, FlipsPerGate, count, seed, sat=False):
    """
    Part de PARITY_n (poids 2 sur chaque variable pour avoir de la marge),
    et perturbe chaque porte de +/-1 sur FlipsPerGate variables choisies parmi
    un sous-ensemble de k_vars variables "bruitees" (toutes les portes ont des
    vecteurs distincts dès que possible).
    """
    rng = np.random.default_rng(seed)
    s = n
    Ws = np.zeros((count, s, n), np.int64)
    ts = np.zeros((count, s), np.int64)
    vs = np.zeros((count, s), np.int64)
    Ts = np.zeros(count, np.int64)
    for c in range(count):
        W = np.full((s, n), 2, dtype=np.int64)
        noisy_cols = rng.choice(n, size=min(n, k_vars), replace=False)
        for g in range(s):
            cols = rng.choice(noisy_cols, size=min(len(noisy_cols), FlipsPerGate), replace=False)
            for col in cols:
                W[g, col] += rng.choice([-1, 1])
        t = np.arange(1, s + 1, dtype=np.int64) * 2
        v = np.array([1 if k % 2 == 0 else -1 for k in range(s)], dtype=np.int64)
        M, _ = top_max(W, t, v)
        Ws[c], ts[c], vs[c] = W, t, v
        Ts[c] = M if sat else M + 1
    return Ws, ts, vs, Ts


def make_two_chain_perturbed(n, k_vars, count, seed, sat=False):
    """
    Deux chaines de parite independantes (chacune sur la moitie des portes,
    avec vecteurs de base distincts w^(1) in {2,0} et w^(2) in {0,2} + couplage)
    et perturbations locales sur k_vars variables.
    """
    rng = np.random.default_rng(seed)
    s = n
    Ws = np.zeros((count, s, n), np.int64)
    ts = np.zeros((count, s), np.int64)
    vs = np.zeros((count, s), np.int64)
    Ts = np.zeros(count, np.int64)
    half = n // 2
    for c in range(count):
        W = np.zeros((s, n), dtype=np.int64)
        t = np.zeros(s, dtype=np.int64)
        v = np.zeros(s, dtype=np.int64)
        # Chaine 1 : portes 0..half-1, poids 3 sur 0..half-1, poids 1 sur half..n-1
        for g in range(half):
            W[g, :half] = 3
            W[g, half:] = 1
            t[g] = (g + 1) * 4
            v[g] = 1 if g % 2 == 0 else -1
        # Chaine 2 : portes half..s-1, poids 1 sur 0..half-1, poids 3 sur half..n-1
        for g in range(half, s):
            W[g, :half] = 1
            W[g, half:] = 3
            t[g] = (g - half + 1) * 4
            v[g] = 1 if (g - half) % 2 == 0 else -1
        noisy_cols = rng.choice(n, size=min(n, k_vars), replace=False)
        for g in range(s):
            col = rng.choice(noisy_cols)
            W[g, col] += rng.choice([-1, 1])
        M, _ = top_max(W, t, v)
        Ws[c], ts[c], vs[c] = W, t, v
        Ts[c] = M if sat else M + 1
    return Ws, ts, vs, Ts


def sanity_check(th8):
    print("== 1. Verification de correction du Theoreme 1 (SAT & UNSAT) ==", flush=True)
    tot = 0
    for sat in (False, True):
        for n in (10, 12, 14):
            for gen in (
                lambda: make_instances(n, 2 * n, 12, seed=100 + n + 50 * sat, sat=sat),
                lambda: make_perturbed_parity(n, k_vars=4, FlipsPerGate=2, count=12, seed=200 + n + 50 * sat, sat=sat),
                lambda: make_two_chain_perturbed(n, k_vars=4, count=12, seed=300 + n + 50 * sat, sat=sat),
            ):
                Ws, ts, vs, Ts = gen()
                st, _ = solve_astar_batch(Ws, ts, vs, Ts, th8, 1 << (n + 2), 5)
                exp = 1 if sat else 0
                assert (st == exp).all(), f"ECHEC sat={sat} n={n}: {st}"
                tot += len(st)
    print(f"   {tot}/{tot} instances SAT et UNSAT verifiees a 100% contre la force brute.\n", flush=True)


def run_benchmark():
    t0 = time.time()
    # Genome hybride precedent (7 poids) + poids fort sur f8 (discordance)
    g1 = np.array([-1.04, -0.06, 0.09, 1.43, -0.57, 0.84, 0.32], dtype=np.float64)
    g_hyb = np.array([1.0, 0.13, 0.86, 0.58, -0.27, 0.83, -0.48], dtype=np.float64)
    th8_no_f8 = np.append(g_hyb, 0.0)
    th8_with_f8 = np.append(g_hyb, 10.0)  # priorite lexicographique aux variables de discordance (Corollaire)

    sanity_check(th8_with_f8)

    ns = [10, 12, 14, 16, 18]
    suites = {
        "P1: Parite perturbee (k=3 vars bruitees)": {
            n: make_perturbed_parity(n, k_vars=3, FlipsPerGate=1, count=8, seed=1100 + n)
            for n in ns
        },
        "P2: Parite perturbee (k=n/3 vars bruitees)": {
            n: make_perturbed_parity(n, k_vars=max(2, n // 3), FlipsPerGate=2, count=8, seed=2200 + n)
            for n in ns
        },
        "P3: Double chaine perturbee (k=4 vars)": {
            n: make_two_chain_perturbed(n, k_vars=4, count=8, seed=3300 + n)
            for n in ns
        },
        "P4: Aleatoires denses (s = 2n)": {
            n: make_instances(n, 2 * n, 10, seed=4400 + n)
            for n in ns
        },
    }

    report = {"ns": ns, "suites": {}}
    for sname, sets in suites.items():
        print(f"== {sname} ==", flush=True)
        report["suites"][sname] = {}
        for mname, runner in [
            ("Base (G1, intervalles)", lambda Ws, ts, vs, Ts, n: solve_hybrid_batch(Ws, ts, vs, Ts, g1, 1 << (n + 2), False, 0)),
            ("Hybride A+B (chaines exactes + LP)", lambda Ws, ts, vs, Ts, n: solve_hybrid_batch(Ws, ts, vs, Ts, g_hyb, 1 << (n + 2), True, 5)),
            ("Theoreme 1 (A* + B, sans f8)", lambda Ws, ts, vs, Ts, n: solve_astar_batch(Ws, ts, vs, Ts, th8_no_f8, 1 << (n + 2), 5)),
            ("Theoreme 1 + Corollaire (A* + B + f8)", lambda Ws, ts, vs, Ts, n: solve_astar_batch(Ws, ts, vs, Ts, th8_with_f8, 1 << (n + 2), 5)),
        ]:
            logs = []
            for n in ns:
                Ws, ts, vs, Ts = sets[n]
                st, nodes = runner(Ws, ts, vs, Ts, n)
                assert (st == 0).all()
                logs.append(math.log2(float(nodes.mean())))
            alpha = float(np.polyfit(ns, logs, 1)[0])
            report["suites"][sname][mname] = {"alpha": alpha, "log2_nodes": logs}
            print(f"   {mname:<38} alpha = {alpha:6.3f}   log2 = "
                  + " ".join(f"{x:5.1f}" for x in logs), flush=True)
        print("", flush=True)

    print(f"Temps total : {time.time() - t0:.1f} s", flush=True)
    with open(os.path.join(RESULTS_DIR, "theorem_astar_results.json"), "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 4, figsize=(18.5, 4.5), sharey=True)
    for ax, (sname, res) in zip(axes, report["suites"].items()):
        ax.plot(ns, [n + 1 for n in ns], "k--", label="Force brute (α=1.00)")
        for mname, r in res.items():
            ax.plot(ns, r["log2_nodes"], "o-", label=f"{mname} (α={r['alpha']:.2f})")
        ax.set_title(sname, fontsize=9.5)
        ax.set_xlabel("n")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=6.8)
    axes[0].set_ylabel("log2(nœuds explorés)")
    fig.tight_layout()
    out_png = os.path.join(OUT, "theorem_astar.png")
    fig.savefig(out_png, dpi=130)
    print(f"Graphique : {out_png}", flush=True)


if __name__ == "__main__":
    run_benchmark()
