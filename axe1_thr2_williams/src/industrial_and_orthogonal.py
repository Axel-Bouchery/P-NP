"""
Phase 4 : Confrontation aux solveurs de reference (Code de Gray compile &
Solveur Pseudo-Booleen a coupes surrogees) et attaque du regime orthogonal s >> n
=================================================================================

Nouvelles briques exactes ajoutees a notre solveur (Solveur Complet v4) :
  1. Propagation unitaire inverse (Sommet -> Portes -> Variables) :
     Si base + up - |v_k| < T, la porte indecise k est FORCEE (a 1 si v_k > 0,
     a 0 si v_k < 0). Toute porte forcee propage des assignations unitaires
     sur les variables libres x_i telles que hi_k - |W_{k,i}| < t_k (ou
     lo_k + |W_{k,i}| >= t_k).
  2. Coupe surrogee (Cutting Plane / Surrogate Constraint) :
     Toutes les portes forcees k in F_1 U F_0 sont combinees en une unique
     inegalite lineaire ponderee sum_{i in F} W^surr_i x_i >= t^surr.
     Si le membre gauche maximal est < t^surr -> conflit immediat (elagage).
     Sinon, cette coupe peut elle-meme forcer des variables x_i qu'aucune
     porte individuelle ne pouvait forcer !
  3. Evaluation exacte des sous-cubes feuilles en code de Gray (Primitive C) :
     Des que le nombre de variables libres tombe a u <= u_leaf, le sous-cube
     {0,1}^u restant est evalue exactement par parcours incremental en code de
     Gray sur les seules portes indecises.
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
from hybrid_thr2 import lagrangian_upper_bound
from theorem_a_star import (NF8, cluster_residual_bound_and_f8,
                            make_perturbed_parity, solve_astar_batch)

AXE_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = str(AXE_ROOT / "results")
FIGURES_DIR = str(AXE_ROOT / "figures")
WITNESS_DIR = str(AXE_ROOT / "witnesses")
for _d in (RESULTS_DIR, FIGURES_DIR, WITNESS_DIR):
    os.makedirs(_d, exist_ok=True)
OUT = sys.argv[1] if len(sys.argv) > 1 else FIGURES_DIR


# --------------------------------------------------------------------------- #
#  Reference 1 : Enumeration exhaustive en code de Gray compilee (style libFES)
# --------------------------------------------------------------------------- #
@njit(cache=True)
def gray_code_exhaust_single(W, t, v, T):
    s, n = W.shape
    total = 1 << n
    acc = np.zeros(s, np.int64)
    x = 0
    for k in range(total):
        if k > 0:
            i = 0
            while not (k >> i) & 1:
                i += 1
            x ^= 1 << i
            sg = 1 if (x >> i) & 1 else -1
            for g in range(s):
                acc[g] += sg * W[g, i]
        top = 0
        for g in range(s):
            if acc[g] >= t[g]:
                top += v[g]
        if top >= T:
            return 1, k + 1
    return 0, total


@njit(parallel=True, cache=True)
def gray_code_exhaust_batch(Ws, ts, vs, Ts):
    B = Ws.shape[0]
    status = np.zeros(B, np.int64)
    steps = np.zeros(B, np.int64)
    for k in prange(B):
        st, nd = gray_code_exhaust_single(Ws[k], ts[k], vs[k], Ts[k])
        status[k] = st
        steps[k] = nd
    return status, steps


# --------------------------------------------------------------------------- #
#  Propagation unitaire inverse + Coupe surrogee (Cutting Plane)
# --------------------------------------------------------------------------- #
@njit(cache=True)
def _assign(W, Wpos, Wneg, S, P, N, i, b, sign):
    s = W.shape[0]
    for g in range(s):
        S[g] += sign * b * W[g, i]
        P[g] -= sign * Wpos[g, i]
        N[g] -= sign * Wneg[g, i]


@njit(cache=True)
def propagate_and_cut(W, Wpos, Wneg, t, v, T, S, P, N, lo, hi, und, assigned,
                      trail_var, trail_val, t0_idx, ws):
    """
    Execute en boucle jusqu'au point fixe :
      - calcul de [lo, hi], base, up, dn
      - detection de conflit (base + up < T) ou succes (base + dn >= T)
      - identification des portes forcees (base + up - |v_g| < T)
      - propagation unitaire depuis chaque porte forcee
      - construction de la coupe surrogee (combinaison des portes forcees)
        et propagation unitaire depuis cette coupe.
    Renvoie (code, new_t_idx) ou code :
      -1 = conflit (elaguer), +1 = SAT trouve, 0 = continuer.
    """
    s, n = W.shape
    t_idx = t0_idx

    while True:
        base = 0
        up = 0
        dn = 0
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
                if v[g] > 0:
                    up += v[g]
                else:
                    dn += v[g]
        if base + dn >= T:
            return 1, t_idx
        if base + up < T:
            return -1, t_idx

        slack = base + up - T
        new_assign = False
        for i in range(n):
            ws[i] = 0
        t_surr = 0
        n_forced_gates = 0

        for g in range(s):
            if not und[g]:
                continue
            vg = v[g]
            avg = abs(vg)
            if avg <= slack:
                continue
            n_forced_gates += 1
            if vg > 0:
                req = t[g] - S[g]
                rem_max = P[g]
                t_surr += avg * req
                for i in range(n):
                    if assigned[i]:
                        continue
                    wgi = W[g, i]
                    ws[i] += avg * wgi
                    if wgi > 0 and rem_max - wgi < req:
                        assigned[i] = True
                        trail_var[t_idx] = i
                        trail_val[t_idx] = 1
                        t_idx += 1
                        _assign(W, Wpos, Wneg, S, P, N, i, 1, 1)
                        new_assign = True
                        break
                    elif wgi < 0 and rem_max + wgi < req:
                        assigned[i] = True
                        trail_var[t_idx] = i
                        trail_val[t_idx] = 0
                        t_idx += 1
                        _assign(W, Wpos, Wneg, S, P, N, i, 0, 1)
                        new_assign = True
                        break
            else:
                req = S[g] - t[g] + 1
                rem_max = -N[g]
                t_surr += avg * req
                for i in range(n):
                    if assigned[i]:
                        continue
                    wgi = -W[g, i]
                    ws[i] += avg * wgi
                    if wgi > 0 and rem_max - wgi < req:
                        assigned[i] = True
                        trail_var[t_idx] = i
                        trail_val[t_idx] = 1
                        t_idx += 1
                        _assign(W, Wpos, Wneg, S, P, N, i, 1, 1)
                        new_assign = True
                        break
                    elif wgi < 0 and rem_max + wgi < req:
                        assigned[i] = True
                        trail_var[t_idx] = i
                        trail_val[t_idx] = 0
                        t_idx += 1
                        _assign(W, Wpos, Wneg, S, P, N, i, 0, 1)
                        new_assign = True
                        break
            if new_assign:
                break

        if new_assign:
            continue

        if n_forced_gates >= 2:
            surr_max = 0
            for i in range(n):
                if not assigned[i] and ws[i] > 0:
                    surr_max += ws[i]
            if surr_max < t_surr:
                return -1, t_idx
            for i in range(n):
                if assigned[i]:
                    continue
                wsi = ws[i]
                if wsi > 0 and surr_max - wsi < t_surr:
                    assigned[i] = True
                    trail_var[t_idx] = i
                    trail_val[t_idx] = 1
                    t_idx += 1
                    _assign(W, Wpos, Wneg, S, P, N, i, 1, 1)
                    new_assign = True
                    break
                elif wsi < 0 and surr_max + wsi < t_surr:
                    assigned[i] = True
                    trail_var[t_idx] = i
                    trail_val[t_idx] = 0
                    t_idx += 1
                    _assign(W, Wpos, Wneg, S, P, N, i, 0, 1)
                    new_assign = True
                    break
            if new_assign:
                continue

        break

    return 0, t_idx


# --------------------------------------------------------------------------- #
#  Primitive C : Evaluation exacte de sous-cube feuille en code de Gray
# --------------------------------------------------------------------------- #
@njit(cache=True)
def eval_leaf_gray(W, t, v, T, S, und, base, free_vars, u, u_gates, n_und):
    acc = np.empty(n_und, np.int64)
    thr = np.empty(n_und, np.int64)
    val = np.empty(n_und, np.int64)
    for idx in range(n_und):
        g = u_gates[idx]
        acc[idx] = S[g]
        thr[idx] = t[g]
        val[idx] = v[g]

    total = 1 << u
    x = 0
    for k in range(total):
        if k > 0:
            bit = 0
            while not (k >> bit) & 1:
                bit += 1
            x ^= 1 << bit
            var_i = free_vars[bit]
            sg = 1 if (x >> bit) & 1 else -1
            for idx in range(n_und):
                acc[idx] += sg * W[u_gates[idx], var_i]
        top = base
        for idx in range(n_und):
            if acc[idx] >= thr[idx]:
                top += val[idx]
        if top >= T:
            return 1
    return 0


# --------------------------------------------------------------------------- #
#  Solveur v4 (PB-Propagation + Coupe surrogee + A* + LP B + Feuilles Gray C)
# --------------------------------------------------------------------------- #
@njit(cache=True)
def solve_v4(W, t, v, T, theta8, node_limit, use_astar, lp_steps, use_pb_cut, u_leaf):
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
    strail_start = np.zeros(n + 1, np.int64)
    strail_end = np.zeros(n + 1, np.int64)
    trail_var = np.zeros(n + 1, np.int64)
    trail_val = np.zeros(n + 1, np.int64)
    t_idx = 0

    und = np.zeros(s, np.bool_)
    lo = np.zeros(s, np.int64)
    hi = np.zeros(s, np.int64)
    feats = np.zeros((n, NF8))
    fmax = np.zeros(NF8)
    f8 = np.zeros(n, np.float64)
    ws = np.zeros(n, np.int64)
    free_vars = np.zeros(n, np.int64)
    u_gates = np.zeros(s, np.int64)
    depth = 0
    nodes = 0

    has_clusters = True
    while True:
        nodes += 1
        if nodes > node_limit:
            return -1, nodes

        t0_idx = t_idx
        if use_pb_cut:
            code, t_idx = propagate_and_cut(
                W, Wpos, Wneg, t, v, T, S, P, N, lo, hi, und, assigned,
                trail_var, trail_val, t0_idx, ws
            )
            if code == 1:
                return 1, nodes
            pruned = (code == -1)
        else:
            pruned = False

        base = 0
        up = 0
        dn = 0
        n_und = 0
        if not pruned:
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
                    u_gates[n_und] = g
                    n_und += 1
                    if v[g] > 0:
                        up += v[g]
                    else:
                        dn += v[g]
            if base + dn >= T:
                return 1, nodes
            if base + up < T:
                pruned = True

        # Compter les variables libres restantes
        u = 0
        if not pruned:
            for i in range(n):
                if not assigned[i]:
                    free_vars[u] = i
                    u += 1
            if u == 0:
                pruned = True

        # Primitive A* (Theoreme 1) avec declenchement adaptatif :
        # toujours teste a la racine (depth == 0) ; si aucun cluster n'ameliore `up`
        # a la racine (circuit aleatoire dense), on economise le cout O(s^2 n) dans l'arbre.
        if not pruned and use_astar and n_und >= 2 and (depth == 0 or has_clusters):
            c_up = cluster_residual_bound_and_f8(W, t, v, S, lo, hi, und, assigned, f8)
            if depth == 0:
                has_clusters = (c_up < up)
            if base + c_up < T:
                pruned = True
        else:
            for i in range(n):
                f8[i] = 0.0

        # Primitive B (LP lagrangien couple) : actif uniquement au-dessus des feuilles
        # (u > u_leaf + 2), la ou elaguer un noeud economise un gros sous-arbre.
        if not pruned and lp_steps > 0 and n_und >= 2 and u > u_leaf + 2:
            ub = lagrangian_upper_bound(W, t, v, S, lo, hi, und, assigned, base, T, lp_steps)
            if ub < T - 1e-7:
                pruned = True

        # Primitive C : sous-cube feuille en code de Gray lorsque u <= u_leaf
        if not pruned and u <= u_leaf:
            if eval_leaf_gray(W, t, v, T, S, und, base, free_vars, u, u_gates, n_und) == 1:
                return 1, nodes
            pruned = True

        if not pruned:
            strail_start[depth] = t0_idx
            strail_end[depth] = t_idx

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

        # Annuler les propagations faites a ce noeud avant de remonter
        while t_idx > t0_idx:
            t_idx -= 1
            vi = trail_var[t_idx]
            vb = trail_val[t_idx]
            _assign(W, Wpos, Wneg, S, P, N, vi, vb, -1)
            assigned[vi] = False

        # Retour arriere
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
            # Annuler aussi la trail du parent d
            ts_start = strail_start[d]
            while t_idx > ts_start:
                t_idx -= 1
                vi = trail_var[t_idx]
                vb = trail_val[t_idx]
                _assign(W, Wpos, Wneg, S, P, N, vi, vb, -1)
                assigned[vi] = False
            depth -= 1


@njit(parallel=True, cache=True)
def solve_v4_batch(Ws, ts, vs, Ts, theta8, node_limit, use_astar, lp_steps, use_pb_cut, u_leaf):
    B = Ws.shape[0]
    status = np.zeros(B, np.int64)
    nodes = np.zeros(B, np.int64)
    for k in prange(B):
        st, nd = solve_v4(Ws[k], ts[k], vs[k], Ts[k], theta8, node_limit,
                          use_astar, lp_steps, use_pb_cut, u_leaf)
        status[k] = st
        nodes[k] = nd
    return status, nodes


# --------------------------------------------------------------------------- #
#  Generateurs : Blocs orthogonaux et Produit scalaire modulo 2 (Inner Product)
# --------------------------------------------------------------------------- #
def make_ortho_blocks(n, n_blocks, count, seed, sat=False):
    """
    Divise les n variables en n_blocks sous-espaces orthogonaux disjoints.
    Chaque bloc porte une chaine de parite locale + un bruit faible global.
    """
    rng = np.random.default_rng(seed)
    s = 2 * n
    Ws = np.zeros((count, s, n), np.int64)
    ts = np.zeros((count, s), np.int64)
    vs = np.zeros((count, s), np.int64)
    Ts = np.zeros(count, np.int64)
    for c in range(count):
        W = np.zeros((s, n), dtype=np.int64)
        t = np.zeros(s, dtype=np.int64)
        v = np.zeros(s, dtype=np.int64)
        for g in range(s):
            b = g % n_blocks
            # Variables du bloc b
            cols = np.arange(b, n, n_blocks)
            W[g, cols] = 4
            # Perturbation locale sur 1 variable hors bloc
            j = rng.integers(n)
            W[g, j] += rng.choice([-1, 1])
            k_in_b = (g // n_blocks) + 1
            t[g] = min(len(cols), k_in_b) * 4
            v[g] = 1 if (g // n_blocks) % 2 == 0 else -1
        M, _ = top_max(W, t, v)
        Ws[c], ts[c], vs[c] = W, t, v
        Ts[c] = M if sat else M + 1
    return Ws, ts, vs, Ts


def make_inner_product_noisy(n, count, seed, sat=False):
    """
    Produit scalaire bruite : m = n/2 paires disjointes (x_{2j}, x_{2j+1})
    avec portes AND/XOR locales + portes de couplage aleatoires (s = 2n).
    C'est l'archetype de fonction a petits supports orthogonaux.
    """
    rng = np.random.default_rng(seed)
    s = 2 * n
    Ws = np.zeros((count, s, n), np.int64)
    ts = np.zeros((count, s), np.int64)
    vs = np.zeros((count, s), np.int64)
    Ts = np.zeros(count, np.int64)
    pairs = n // 2
    for c in range(count):
        W = np.zeros((s, n), dtype=np.int64)
        t = np.zeros(s, dtype=np.int64)
        v = np.zeros(s, dtype=np.int64)
        # Portes 0..n-1 : paires locales (XOR / AND sur chaque paire + couplage faible)
        for p in range(pairs):
            i1, i2 = 2 * p, 2 * p + 1
            # Porte >= 1 sur la paire
            W[2 * p, i1] = 4
            W[2 * p, i2] = 4
            t[2 * p] = 4
            v[2 * p] = 2
            # Porte >= 2 sur la paire (avec signe -4 -> XOR de la paire)
            W[2 * p + 1, i1] = 4
            W[2 * p + 1, i2] = 4
            t[2 * p + 1] = 8
            v[2 * p + 1] = -4
        # Portes n..2n-1 : contraintes globales aleatoires couplant toutes les paires
        W[n:, :] = rng.integers(-3, 4, (n, n))
        t[n:] = np.floor(W[n:, :].sum(axis=1) / 2).astype(np.int64)
        v[n:] = rng.choice([-2, -1, 1, 2], size=n)
        M, _ = top_max(W, t, v)
        Ws[c], ts[c], vs[c] = W, t, v
        Ts[c] = M if sat else M + 1
    return Ws, ts, vs, Ts


def sanity_check(th8):
    print("== 1. Verification de correction du Solveur v4 (SAT & UNSAT) ==", flush=True)
    tot = 0
    for sat in (False, True):
        for n in (10, 12, 14):
            for gen in (
                lambda: make_instances(n, 2 * n, 12, seed=510 + n + 50 * sat, sat=sat),
                lambda: make_ortho_blocks(n, 4, 12, seed=610 + n + 50 * sat, sat=sat),
                lambda: make_inner_product_noisy(n, 12, seed=710 + n + 50 * sat, sat=sat),
            ):
                Ws, ts, vs, Ts = gen()
                for cfg in (
                    (False, 0, True, 0),   # Solveur Pseudo-Booleen pur (PB + coupes surrogees)
                    (True, 5, True, 5),    # Solveur Complet v4 (A* + LP B + PB/Coupes + Feuilles Gray C)
                ):
                    st, _ = solve_v4_batch(Ws, ts, vs, Ts, th8, 1 << (n + 2), *cfg)
                    exp = 1 if sat else 0
                    assert (st == exp).all(), f"ECHEC sat={sat} n={n} cfg={cfg}: {st}"
                    tot += len(st)
    print(f"   {tot}/{tot} tests SAT et UNSAT verifies a 100% contre la force brute.\n", flush=True)


def main():
    t0 = time.time()
    g_hyb = np.array([1.0, 0.13, 0.86, 0.58, -0.27, 0.83, -0.48], dtype=np.float64)
    th8 = np.append(g_hyb, 10.0)
    th_pb = np.array([0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float64)  # activite f2 standard

    sanity_check(th8)

    ns = [12, 14, 16, 18, 20]
    suites = {
        "O1: 4 Blocs orthogonaux perturbes (s=2n)": {
            n: make_ortho_blocks(n, 4, 8, seed=8100 + n) for n in ns
        },
        "O2: Produit scalaire + couplage (s=2n)": {
            n: make_inner_product_noisy(n, 8, seed=8200 + n) for n in ns
        },
        "O3: Haute densite aleatoire (s=3n)": {
            n: make_instances(n, 3 * n, 8, seed=8300 + n) for n in ns
        },
        "O4:Aleatoires s=2n (comparaison temps reel)": {
            n: make_instances(n, 2 * n, 10, seed=8400 + n) for n in ns
        },
    }

    solvers = [
        ("Ref 1 : Code de Gray compile (type libFES)",
         lambda Ws, ts, vs, Ts, n: gray_code_exhaust_batch(Ws, ts, vs, Ts)),
        ("Ref 2 : Solveur Pseudo-Booleen (PB + Coupes)",
         lambda Ws, ts, vs, Ts, n: solve_v4_batch(Ws, ts, vs, Ts, th_pb, 1 << (n + 2), False, 0, True, 0)),
        ("Solveur v3 : Theoreme 1 (A* + LP B + f8)",
         lambda Ws, ts, vs, Ts, n: solve_v4_batch(Ws, ts, vs, Ts, th8, 1 << (n + 2), True, 5, False, 0)),
        ("Solveur v4 : Complet (A* + B + PB/Coupes + Gray C)",
         lambda Ws, ts, vs, Ts, n: solve_v4_batch(Ws, ts, vs, Ts, th8, 1 << (n + 2), True, 3, True, 8)),
    ]

    report = {"ns": ns, "suites": {}}
    for sname, sets in suites.items():
        print(f"== {sname} ==", flush=True)
        report["suites"][sname] = {}
        for mname, fn in solvers:
            logs = []
            t_tot = 0.0
            t_last = 0.0
            for n in ns:
                Ws, ts, vs, Ts = sets[n]
                ts0 = time.perf_counter()
                st, nodes = fn(Ws, ts, vs, Ts, n)
                dt = time.perf_counter() - ts0
                assert (st == 0).all()
                logs.append(math.log2(float(nodes.mean())))
                t_tot += dt
                t_last = dt
            alpha = float(np.polyfit(ns, logs, 1)[0])
            report["suites"][sname][mname] = {
                "alpha": alpha, "log2_nodes": logs, "time_total_s": t_tot, "time_n20_s": t_last
            }
            print(f"   {mname:<48} alpha = {alpha:6.3f} | t(n=20)={t_last:5.3f}s | log2 = "
                  + " ".join(f"{x:5.1f}" for x in logs), flush=True)
        print("", flush=True)

    # Test de passage a l'echelle en temps reel jusqu'a n = 26 sur s = 2n
    big_ns = [18, 20, 22, 24, 26]
    print("== Duel en temps reel sur n = 18..26 (s = 2n, 8 instances par taille) ==", flush=True)
    big_sets = {n: make_instances(n, 2 * n, 8, seed=9500 + n) for n in big_ns}
    report["wallclock_duel"] = {"ns": big_ns, "methods": {}}
    for mname, fn in solvers:
        times = []
        logs = []
        for n in big_ns:
            Ws, ts, vs, Ts = big_sets[n]
            ts0 = time.perf_counter()
            st, nodes = fn(Ws, ts, vs, Ts, n)
            dt = time.perf_counter() - ts0
            assert (st == 0).all()
            times.append(dt)
            logs.append(math.log2(float(nodes.mean())))
        alpha = float(np.polyfit(big_ns, logs, 1)[0])
        report["wallclock_duel"]["methods"][mname] = {
            "alpha": alpha, "times_s": times, "log2_nodes": logs
        }
        print(f"   {mname:<48} alpha = {alpha:6.3f} | temps(s) = "
              + " ".join(f"{t:6.3f}" for t in times), flush=True)

    print(f"\nTemps total Phase 4 : {time.time() - t0:.1f} s", flush=True)
    with open(os.path.join(RESULTS_DIR, "phase4_results.json"), "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 4, figsize=(19, 4.6))
    suite_keys = [
        "O1: 4 Blocs orthogonaux perturbes (s=2n)",
        "O2: Produit scalaire + couplage (s=2n)",
        "O3: Haute densite aleatoire (s=3n)",
    ]
    for ax, sname in zip(axes[:3], suite_keys):
        res = report["suites"][sname]
        for mname, r in res.items():
            ls = "k--" if "Ref 1" in mname else "o-"
            ax.plot(ns, r["log2_nodes"], ls, label=f"{mname.split(':')[0]} (α={r['alpha']:.2f})")
        ax.set_title(sname, fontsize=9.2)
        ax.set_xlabel("n")
        ax.set_ylabel("log2(nœuds / étapes)")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7.0)

    ax = axes[3]
    for mname, r in report["wallclock_duel"]["methods"].items():
        ls = "k--" if "Ref 1" in mname else "o-"
        ax.semilogy(big_ns, [max(1e-4, t) for t in r["times_s"]], ls,
                    label=f"{mname.split(':')[0]} ({r['times_s'][-1]:.2f}s @ n=26)")
    ax.set_title("Temps réel mesuré (n = 18..26, s = 2n)", fontsize=9.2)
    ax.set_xlabel("n")
    ax.set_ylabel("Temps réel (secondes, échelle log)")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=7.0)

    fig.tight_layout()
    out_png = os.path.join(OUT, "phase4_duel.png")
    fig.savefig(out_png, dpi=130)
    print(f"Graphique : {out_png}", flush=True)


if __name__ == "__main__":
    main()
