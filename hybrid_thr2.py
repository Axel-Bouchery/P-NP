"""
Partie 2 : Solveur hybride pour THR∘THR-SAT + Evolution
=======================================================

Ajoute deux primitives mathematiquement exactes au Branch-and-Bound :
  - Primitive A (Groupement en chaines 1D) :
      Les portes indecises partageant le meme vecteur de poids sur les variables
      libres forment une fonction en escalier d'une seule somme L(x). On calcule
      son MAXIMUM EXACT sur [L_min, L_max] -> elimine instantanement PARITY,
      MOD3, et WEIGHTED-PARITY (1 seul noeud !).
  - Primitive B (Borne lagrangienne / enveloppe lineaire couplee) :
      Pour chaque porte indecise k :
        v_k > 0 : g_k(x) <= (1 - lam_k) + lam_k * (S_k(x) - lo_k) / (t_k - lo_k)
        v_k < 0 : v_k g_k(x) <= lam_k * v_k * (S_k(x) - (t_k - 1)) / (hi_k - t_k + 1)
      En sommant sur k, les variables libres x_i se couplent entre toutes les portes :
        sum v_k g_k(x) <= C(lam) + sum_{i libre} max(0, c_i(lam)).
      On optimise lam in [0, 1]^s par quelques pas de sous-gradient projetes.
      Toute valeur de lam in [0,1]^s donne une borne superieure 100% valide.

Le script :
  1. Verifie l'exactitude (100% d'accord avec la force brute sur SAT et UNSAT).
  2. Evolue un nouveau genome hybride (7 poids de branchement + frequence LP).
  3. Compare Solveur de base vs Solveur hybride sur :
       - instances aleatoires s = n et s = 2n (jusqu'a n = 26)
       - banc d'essai adversarial (PARITY, MOD3, WEIGHTED-PARITY, LARGE-WEIGHTS).
"""
import json
import math
import os
import sys
import time

import numpy as np
from numba import njit, prange

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from thr2_lab import HERE, NF, make_instances, top_max
from adversarial_suite import (make_parity, make_mod3, make_weighted_parity,
                               make_large_weights, pack_single)

OUT = sys.argv[1] if len(sys.argv) > 1 else HERE


# --------------------------------------------------------------------------- #
#  Primitive A : borne exacte par regroupement des portes de meme vecteur libre
# --------------------------------------------------------------------------- #
@njit(cache=True)
def chain_upper_bound(W, t, v, S, lo, hi, und, assigned):
    """
    Regroupe les portes indecises ayant exactement les memes poids sur les
    variables libres (signature de hachage + verification exacte) et calcule
    la somme des maxima exacts de chaque chaine 1D.
    """
    s, n = W.shape
    used = np.zeros(s, np.bool_)
    total_up = 0
    taus = np.zeros(s, np.int64)
    vals = np.zeros(s, np.int64)

    for g in range(s):
        if not und[g] or used[g]:
            continue
        # Chercher toutes les portes indecises ayant le meme vecteur libre que g
        cnt = 0
        taus[cnt] = t[g] - S[g]
        vals[cnt] = v[g]
        cnt += 1
        used[g] = True

        for h in range(g + 1, s):
            if not und[h] or used[h]:
                continue
            same = True
            for i in range(n):
                if not assigned[i] and W[g, i] != W[h, i]:
                    same = False
                    break
            if same:
                taus[cnt] = t[h] - S[h]
                vals[cnt] = v[h]
                cnt += 1
                used[h] = True

        if cnt == 1:
            if vals[0] > 0:
                total_up += vals[0]
        else:
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
            # Balayage sur L in [L_min, L_max] :
            # avant le premier seuil, somme = 0 ; apres chaque seuil distinct, on ajoute
            # tous les vals partageant ce seuil.
            cur = 0
            best = 0  # L < min(taus) est toujours atteignable car toutes les portes sont indecises (lo < t)
            idx = 0
            while idx < cnt:
                tau_cur = taus[idx]
                while idx < cnt and taus[idx] == tau_cur:
                    cur += vals[idx]
                    idx += 1
                if cur > best:
                    best = cur
            total_up += best
    return total_up


# --------------------------------------------------------------------------- #
#  Primitive B : borne lagrangienne couplee (enveloppe lineaire des seuils)
# --------------------------------------------------------------------------- #
@njit(cache=True)
def lagrangian_upper_bound(W, t, v, S, lo, hi, und, assigned, base, T, n_steps):
    """
    Calcule une borne superieure 100% valide sur sum_{k in U} v_k g_k(x) + base
    en couplant les portes via les variables partagees x_i.
    Optimise lam_k in [0, 1] par sous-gradient projete.
    """
    s, n = W.shape
    lam = np.full(s, 0.5, np.float64)
    c = np.zeros(n, np.float64)
    best_ub = 1e300

    # Precalcul des pentes alpha_k et constantes beta_k pour chaque porte indecise
    # v_k g_k(x) <= (1 - lam_k) * u0_k + lam_k * (beta_k + sum_i alpha_k * W[k,i] * x_i)
    alpha = np.zeros(s, np.float64)
    beta = np.zeros(s, np.float64)
    u0 = np.zeros(s, np.float64)
    for k in range(s):
        if not und[k]:
            continue
        if v[k] > 0:
            denom = float(t[k] - lo[k])  # >= 1
            u0[k] = float(v[k])
            alpha[k] = float(v[k]) / denom
            beta[k] = alpha[k] * float(S[k] - lo[k])
        else:
            denom = float(hi[k] - t[k] + 1)  # >= 1
            u0[k] = 0.0
            alpha[k] = float(v[k]) / denom
            beta[k] = alpha[k] * float(S[k] - t[k] + 1)

    step = 0.6
    for it in range(n_steps):
        const_term = float(base)
        for i in range(n):
            c[i] = 0.0
        for k in range(s):
            if not und[k]:
                continue
            lk = lam[k]
            const_term += (1.0 - lk) * u0[k] + lk * beta[k]
            ak = lk * alpha[k]
            for i in range(n):
                if not assigned[i]:
                    c[i] += ak * W[k, i]
        ub = const_term
        for i in range(n):
            if not assigned[i] and c[i] > 0.0:
                ub += c[i]
        if ub < best_ub:
            best_ub = ub
            if best_ub < T - 1e-7:
                return best_ub
        # Pas de sous-gradient sur lam_k
        if it + 1 < n_steps:
            for k in range(s):
                if not und[k]:
                    continue
                # derivee de ub par rapport a lam_k a x*_i = 1[c[i] > 0]
                dot = 0.0
                for i in range(n):
                    if not assigned[i] and c[i] > 0.0:
                        dot += W[k, i]
                grad = -u0[k] + beta[k] + alpha[k] * dot
                # normalisation douce par |v[k]|
                av = abs(float(v[k]))
                if av > 0.0:
                    lam[k] -= step * (grad / av)
                if lam[k] < 0.0:
                    lam[k] = 0.0
                elif lam[k] > 1.0:
                    lam[k] = 1.0
            step *= 0.65
    return best_ub


# --------------------------------------------------------------------------- #
#  Solveur hybride exact (numba)
# --------------------------------------------------------------------------- #
@njit(cache=True)
def _assign(W, Wpos, Wneg, S, P, N, i, b, sign):
    s = W.shape[0]
    for g in range(s):
        S[g] += sign * b * W[g, i]
        P[g] -= sign * Wpos[g, i]
        N[g] -= sign * Wneg[g, i]


@njit(cache=True)
def solve_hybrid(W, t, v, T, theta, node_limit, use_chain, lp_steps):
    """
    Solveur exact avec :
      - elagage par intervalles (toujours actif)
      - Primitive A : groupement en chaines 1D (si use_chain)
      - Primitive B : borne lagrangienne couplee (si lp_steps > 0)
    """
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
    feats = np.zeros((n, NF))
    fmax = np.zeros(NF)
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

        # Primitive A : chaines 1D exactes (seulement si pas deja elague)
        if not pruned and use_chain and n_und >= 2:
            c_up = chain_upper_bound(W, t, v, S, lo, hi, und, assigned)
            if base + c_up < T:
                pruned = True

        # Primitive B : borne lagrangienne couplee
        if not pruned and lp_steps > 0 and n_und >= 2:
            ub = lagrangian_upper_bound(W, t, v, S, lo, hi, und, assigned, base, T, lp_steps)
            if ub < T - 1e-7:
                pruned = True

        if not pruned:
            # Choix de variable par genome
            for f in range(NF):
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
                for f in range(NF):
                    if feats[i, f] > fmax[f]:
                        fmax[f] = feats[i, f]
            best = -1
            bs = -1e300
            for i in range(n):
                if assigned[i]:
                    continue
                sc = 0.0
                for f in range(NF):
                    sc += theta[f] * feats[i, f] / fmax[f]
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
            depth -= 1


@njit(parallel=True, cache=True)
def solve_hybrid_batch(Ws, ts, vs, Ts, theta, node_limit, use_chain, lp_steps):
    B = Ws.shape[0]
    status = np.zeros(B, np.int64)
    nodes = np.zeros(B, np.int64)
    for k in prange(B):
        st, nd = solve_hybrid(Ws[k], ts[k], vs[k], Ts[k], theta, node_limit, use_chain, lp_steps)
        status[k] = st
        nodes[k] = nd
    return status, nodes


def sanity_check():
    print("== 1. Verification de correction du solveur hybride (SAT & UNSAT) ==", flush=True)
    th = np.array([-1.04, -0.06, 0.09, 1.43, -0.57, 0.84, 0.32], dtype=np.float64)
    total = 0
    for sat in (False, True):
        for n in (8, 10, 12, 14):
            Ws, ts, vs, Ts = make_instances(n, 2 * n, 16, seed=1234 + n + 99 * sat, sat=sat)
            st, _ = solve_hybrid_batch(Ws, ts, vs, Ts, th, 1 << (n + 2), True, 6)
            expect = 1 if sat else 0
            assert (st == expect).all(), f"ECHEC sur sat={sat}, n={n}: {st}"
            total += len(st)
    # Verifier aussi sur PARITY, MOD3, WEIGHTED-PARITY en versions SAT (T=M) et UNSAT (T=M+1)
    for n in (10, 12):
        for fn in (make_parity, make_mod3):
            W, t, v, T = fn(n)
            for delta, exp in ((0, 0), (-1, 1)):
                st, _ = solve_hybrid(W, t, v, T + delta, th, 1 << (n + 2), True, 6)
                assert st == exp, f"ECHEC sur {fn.__name__} n={n} delta={delta}: {st}"
                total += 1
    print(f"   {total}/{total} tests passes avec 100% d'exactitude.\n", flush=True)


def measure_hybrid(theta, sets, ns, use_chain, lp_steps):
    logs = []
    times = []
    for n in ns:
        Ws, ts, vs, Ts = sets[n]
        t0 = time.perf_counter()
        st, nodes = solve_hybrid_batch(Ws, ts, vs, Ts, theta, 1 << (n + 2), use_chain, lp_steps)
        dt = time.perf_counter() - t0
        assert (st == 0).all(), f"Erreur statut n={n}: {st}"
        logs.append(math.log2(float(nodes.mean())))
        times.append(dt)
    alpha = float(np.polyfit(ns, logs, 1)[0])
    return alpha, logs, sum(times)


def evolve_hybrid(train, ns, gens, pop, rng, lp_steps=5):
    elite = [np.zeros(NF)] + [np.eye(NF)[k] for k in range(NF)]
    scored = [(measure_hybrid(th, train, ns, True, lp_steps)[0], th) for th in elite]
    scored.sort(key=lambda z: z[0])
    sigma = 0.6
    for gen in range(gens):
        parents = [th for _, th in scored[:4]]
        children = []
        for j in range(pop):
            if j % 4 == 3:
                a, b = rng.choice(len(parents), 2, replace=False)
                w = rng.random(NF)
                child = w * parents[a] + (1 - w) * parents[b]
            else:
                child = parents[j % len(parents)] + sigma * rng.normal(size=NF)
            children.append(child)
        scored += [(measure_hybrid(th, train, ns, True, lp_steps)[0], th) for th in children]
        scored.sort(key=lambda z: z[0])
        scored = scored[:8]
        sigma *= 0.92
        print(f"   generation {gen + 1:>2} : meilleur alpha hybride = {scored[0][0]:.4f}", flush=True)
    return scored[0][1], scored[0][0]


def main():
    t0 = time.time()
    rng = np.random.default_rng(2026)
    sanity_check()

    # Charger Genome 1 de base
    res_path = os.path.join(HERE, "thr2_results.json")
    if os.path.exists(res_path):
        with open(res_path, "r", encoding="utf-8") as fh:
            g1 = np.array(json.load(fh)["genome1"], dtype=np.float64)
    else:
        g1 = np.array([-1.04, -0.06, 0.09, 1.43, -0.57, 0.84, 0.32], dtype=np.float64)

    print("== 2. Evolution du genome hybride (sur s = 2n, UNSAT de justesse) ==", flush=True)
    train_ns = [10, 12, 14, 16]
    train = {n: make_instances(n, 2 * n, 10, seed=2000 + n) for n in train_ns}
    g_hyb, a_tr = evolve_hybrid(train, train_ns, gens=12, pop=12, rng=rng, lp_steps=5)
    print(f"   genome hybride : {np.round(g_hyb, 2).tolist()}  (alpha train = {a_tr:.4f})\n", flush=True)

    configs = [
        ("Base (G1, intervalles seuls)", g1, False, 0),
        ("Hybride A (G1 + chaines 1D)", g1, True, 0),
        ("Hybride B (G1 + LP couple)", g1, False, 5),
        ("Hybride A+B (G_hyb + chaines + LP)", g_hyb, True, 5),
    ]

    results = {}

    # Evaluation sur aleatoires s=n et s=2n (n = 12..24)
    val_ns = [12, 14, 16, 18, 20, 22, 24]
    for label, mult in [("Aleatoires s = n", 1), ("Aleatoires s = 2n", 2)]:
        print(f"== 3. Comparaison sur {label} (n = 12..24) ==", flush=True)
        sets = {n: make_instances(n, mult * n, 12, seed=9000 + 17 * n + mult) for n in val_ns}
        results[label] = {"ns": val_ns, "m": {}}
        for cname, th, uc, lps in configs:
            a, logs, dt = measure_hybrid(th, sets, val_ns, uc, lps)
            results[label]["m"][cname] = {"alpha": a, "log2_nodes": logs, "time_s": dt}
            print(f"   {cname:<36} alpha = {a:.3f}  ({dt:5.2f}s)  log2 = "
                  + " ".join(f"{x:5.1f}" for x in logs), flush=True)
        print("", flush=True)

    # Evaluation sur les 4 familles adversariales structurees (n = 10..18)
    adv_ns = [10, 12, 14, 16, 18]
    adv_sets = {
        "PARITY_n": {n: pack_single([make_parity(n)]) for n in adv_ns},
        "MOD3_n": {n: pack_single([make_mod3(n)]) for n in adv_ns},
        "WEIGHTED-PARITY": {
            n: pack_single([make_weighted_parity(n, seed=300 + 17 * n + k) for k in range(6)])
            for n in adv_ns
        },
        "LARGE-WEIGHTS (s=2n)": {
            n: make_large_weights(n, 2 * n, count=6, seed=700 + 19 * n, R=1000)
            for n in adv_ns
        },
    }
    print("== 4. Comparaison sur les 4 familles adversariales (n = 10..18) ==", flush=True)
    results["adversarial"] = {"ns": adv_ns, "families": {}}
    for fname, sets in adv_sets.items():
        print(f"  -- {fname} --", flush=True)
        results["adversarial"]["families"][fname] = {}
        for cname, th, uc, lps in configs:
            a, logs, dt = measure_hybrid(th, sets, adv_ns, uc, lps)
            results["adversarial"]["families"][fname][cname] = {
                "alpha": a, "log2_nodes": logs, "time_s": dt
            }
            print(f"   {cname:<36} alpha = {a:.3f}  log2 = "
                  + " ".join(f"{x:5.1f}" for x in logs), flush=True)

    print(f"\nTemps total Partie 2 : {time.time() - t0:.1f} s", flush=True)
    with open(os.path.join(HERE, "hybrid_results.json"), "w", encoding="utf-8") as fh:
        json.dump({"genome_hybrid": g_hyb.tolist(), "results": results}, fh, indent=2)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(17, 4.8))
    for ax, label in zip(axes[:2], ["Aleatoires s = n", "Aleatoires s = 2n"]):
        ns = results[label]["ns"]
        ax.plot(ns, [n + 1 for n in ns], "k--", label="Force brute (α=1.00)")
        for cname, r in results[label]["m"].items():
            ax.plot(ns, r["log2_nodes"], "o-", label=f"{cname} (α={r['alpha']:.2f})")
        ax.set_title(label)
        ax.set_xlabel("n")
        ax.set_ylabel("log2(nœuds)")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7.2)

    ax = axes[2]
    ax.plot(adv_ns, [n + 1 for n in adv_ns], "k--", label="Force brute (α=1.00)")
    for fname in ("PARITY_n", "LARGE-WEIGHTS (s=2n)"):
        r_base = results["adversarial"]["families"][fname]["Base (G1, intervalles seuls)"]
        r_hyb = results["adversarial"]["families"][fname]["Hybride A+B (G_hyb + chaines + LP)"]
        ax.plot(adv_ns, r_base["log2_nodes"], "o--", label=f"Base / {fname} (α={r_base['alpha']:.2f})")
        ax.plot(adv_ns, r_hyb["log2_nodes"], "s-", label=f"Hybride / {fname} (α={r_hyb['alpha']:.2f})")
    ax.set_title("Impact sur familles adversariales")
    ax.set_xlabel("n")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7.2)

    fig.tight_layout()
    out_png = os.path.join(OUT, "hybrid_comparison.png")
    fig.savefig(out_png, dpi=130)
    print(f"Graphique : {out_png}", flush=True)


if __name__ == "__main__":
    main()
