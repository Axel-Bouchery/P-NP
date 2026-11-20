"""
Laboratoire : recherche evolutionnaire de strategies SAT pour circuits THR∘THR
=============================================================================

Circuit THR∘THR (reseau de neurones binaire a 2 couches) :
  porte g du bas  : g(x) = [ sum_i W[g,i] x_i >= t[g] ]
  porte du sommet : C(x) = [ sum_g v[g] g(x) >= T ]
Question : existe-t-il x dans {0,1}^n tel que C(x) = 1 ?

Probleme ouvert (frontiere de la methode de Williams) : un algorithme en
2^n / n^omega(1) pour THR∘THR-SAT prouverait NEXP ⊄ THR∘THR.

Solveur : exploration en profondeur avec elagage par intervalles.
  - pour chaque porte du bas, l'intervalle [lo, hi] de sa somme ponderee
    sur les completions possibles => porte decidee a 0, a 1, ou indecise ;
  - pour la porte du sommet, idem => on coupe la branche si le circuit
    ne peut plus valoir 1, on s'arrete si il vaut 1 quelle que soit la suite.
Le "genome" = 7 poids theta qui definissent la regle de choix de variable.
L'evolution cherche le genome minimisant l'exposant alpha (noeuds ~ 2^(alpha n)).
Force brute : alpha = 1.
"""
import json
import math
import os
import sys
import time

import numpy as np
from numba import njit, prange

HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault("MPLCONFIGDIR", os.path.join(HERE, ".mpl"))

NF = 7
FEATURE_NAMES = [
    "f1 somme |W| (portes indecises)",
    "f2 somme |v|·|W|",
    "f3 somme |W|/largeur intervalle",
    "f4 somme |v|·|W|/largeur",
    "f5 nb portes decidables seule",
    "f6 somme |v| portes decidables seule",
    "f7 somme |W|/distance au seuil",
]


# --------------------------------------------------------------------------- #
#  Solveur (numba)
# --------------------------------------------------------------------------- #
@njit(cache=True)
def _assign(W, Wpos, Wneg, S, P, N, i, b, sign):
    s = W.shape[0]
    for g in range(s):
        S[g] += sign * b * W[g, i]
        P[g] -= sign * Wpos[g, i]
        N[g] -= sign * Wneg[g, i]


@njit(cache=True)
def solve(W, t, v, T, theta, node_limit):
    """Renvoie (statut, noeuds) ; statut 1 = SAT, 0 = UNSAT, -1 = limite atteinte."""
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
            return 1, nodes
        if base + up >= T:
            # ---- choix de la variable (le "genome") ----
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
            # valeur essayee en premier : celle qui pousse le sommet vers 1
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
        # ---- retour arriere ----
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
def solve_batch(Ws, ts, vs, Ts, theta, node_limit):
    B = Ws.shape[0]
    status = np.zeros(B, np.int64)
    nodes = np.zeros(B, np.int64)
    for k in prange(B):
        st, nd = solve(Ws[k], ts[k], vs[k], Ts[k], theta, node_limit)
        status[k] = st
        nodes[k] = nd
    return status, nodes


# --------------------------------------------------------------------------- #
#  Generation d'instances "insatisfiables de justesse" (et "SAT rares")
# --------------------------------------------------------------------------- #
@njit(parallel=True, cache=True)
def top_max(W, t, v):
    s, n = W.shape
    nblocks = 256
    total = 1 << n
    per = (total + nblocks - 1) // nblocks
    mx = np.full(nblocks, -(1 << 40), np.int64)
    cnt = np.zeros(nblocks, np.int64)
    for blk in prange(nblocks):
        k0 = blk * per
        k1 = min(total, (blk + 1) * per)
        if k0 >= k1:
            continue
        x = k0 ^ (k0 >> 1)                     # parcours en code de Gray
        acc = np.zeros(s, np.int64)
        for g in range(s):
            for i in range(n):
                if (x >> i) & 1:
                    acc[g] += W[g, i]
        for k in range(k0, k1):
            if k > k0:
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
            if top > mx[blk]:
                mx[blk] = top
                cnt[blk] = 1
            elif top == mx[blk]:
                cnt[blk] += 1
    M = mx.max()
    c = 0
    for blk in range(nblocks):
        if mx[blk] == M:
            c += cnt[blk]
    return M, c


def make_instances(n, s, count, seed, sat=False, R=3, density=0.6):
    rng = np.random.default_rng(seed)
    Ws = np.zeros((count, s, n), np.int64)
    ts = np.zeros((count, s), np.int64)
    vs = np.zeros((count, s), np.int64)
    Ts = np.zeros(count, np.int64)
    for k in range(count):
        W = rng.integers(-R, R + 1, (s, n)) * (rng.random((s, n)) < density)
        t = np.floor(W.sum(axis=1) / 2).astype(np.int64) + rng.integers(-1, 2, s)
        v = rng.integers(1, R + 1, s) * rng.choice([-1, 1], s)
        M, _ = top_max(W.astype(np.int64), t, v.astype(np.int64))
        Ws[k], ts[k], vs[k] = W, t, v
        Ts[k] = M if sat else M + 1          # SAT : solutions rares ; UNSAT : rate d'un point
    return Ws, ts, vs, Ts


# --------------------------------------------------------------------------- #
#  Mesure de l'exposant
# --------------------------------------------------------------------------- #
def measure(theta, sets, ns):
    """log2(noeuds moyens) pour chaque n, et pente alpha (moindres carres)."""
    theta = np.asarray(theta, np.float64)
    logs = []
    for n in ns:
        Ws, ts, vs, Ts = sets[n]
        _, nodes = solve_batch(Ws, ts, vs, Ts, theta, 1 << (n + 1))
        logs.append(math.log2(nodes.mean()))
    alpha = np.polyfit(ns, logs, 1)[0]
    return alpha, logs


def evolve(train, ns, gens, pop, rng, log):
    elite = [np.zeros(NF)] + [np.eye(NF)[k] for k in range(NF)]   # graines simples
    scored = [(measure(th, train, ns)[0], th) for th in elite]
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
        scored += [(measure(th, train, ns)[0], th) for th in children]
        scored.sort(key=lambda z: z[0])
        scored = scored[:8]
        sigma *= 0.93
        log(f"   generation {gen + 1:>2} : meilleur alpha = {scored[0][0]:.4f}")
    return scored[0][1], scored[0][0]


def main():
    out_dir = sys.argv[1] if len(sys.argv) > 1 else HERE
    rng = np.random.default_rng(7)
    lines = []

    def log(msg):
        print(msg, flush=True)
        lines.append(msg)

    t0 = time.time()
    train_ns = [10, 12, 14, 16]
    log("== Generation des instances d'entrainement (s = n portes, UNSAT de justesse) ==")
    train = {n: make_instances(n, n, 12, seed=1000 + n) for n in train_ns}

    log("== Evolution ==")
    best, alpha_train = evolve(train, train_ns, gens=20, pop=16, rng=rng, log=log)
    log(f"   genome trouve : {np.round(best, 3).tolist()}")
    log(f"   (temps : {time.time() - t0:.0f} s)\n")

    candidates = {
        "Ordre fixe (elagage seul)": np.zeros(NF),
        "Heuristique f2 (a la main)": np.eye(NF)[1],
        "Heuristique f7 (a la main)": np.eye(NF)[6],
        "Evolue": best,
    }

    results = {}
    for label, fam_s, sat in [("UNSAT de justesse, s = n", 1, False),
                              ("UNSAT de justesse, s = 2n", 2, False),
                              ("SAT rares, s = n", 1, True)]:
        val_ns = [10, 12, 14, 16, 18, 20]
        log(f"== Validation sur instances NOUVELLES : {label} ==")
        sets = {n: make_instances(n, fam_s * n, 12, seed=5000 + 31 * n + 7 * fam_s + sat, sat=sat)
                for n in val_ns}
        fam = {}
        for name, th in candidates.items():
            a, logs = measure(th, sets, val_ns)
            fam[name] = {"alpha": a, "log2_nodes": logs}
            log(f"   {name:<28} alpha = {a:.3f}   log2(noeuds) = "
                + " ".join(f"{x:5.1f}" for x in logs))
        log(f"   {'Force brute':<28} alpha = 1.000   log2(noeuds) = "
            + " ".join(f"{n + 1:5.1f}" for n in val_ns) + "\n")
        results[label] = {"ns": val_ns, "methods": fam}

    log(f"Temps total : {time.time() - t0:.0f} s")
    with open(os.path.join(HERE, "thr2_results.json"), "w", encoding="utf-8") as fh:
        json.dump({"genome": best.tolist(), "alpha_train": alpha_train,
                   "features": FEATURE_NAMES, "results": results}, fh, indent=2)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6), sharey=True)
    colors = {"Ordre fixe (elagage seul)": "#ff7f0e", "Heuristique f2 (a la main)": "#9467bd",
              "Heuristique f7 (a la main)": "#8c564b", "Evolue": "#2ca02c"}
    for ax, (label, res) in zip(axes, results.items()):
        ns = res["ns"]
        ax.plot(ns, [n + 1 for n in ns], "k--", label="Force brute (α = 1)")
        for name, r in res["methods"].items():
            ax.plot(ns, r["log2_nodes"], "o-", color=colors[name],
                    label=f"{name} (α = {r['alpha']:.2f})")
        ax.set_title(label)
        ax.set_xlabel("n")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7.5)
    axes[0].set_ylabel("log2(nœuds explorés)")
    fig.tight_layout()
    out = os.path.join(out_dir, "thr2_evolution.png")
    fig.savefig(out, dpi=130)
    log(f"Graphique : {out}")


if __name__ == "__main__":
    main()
