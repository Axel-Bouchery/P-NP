"""
Partie 1 : Banc d'essai adversarial structure pour THR∘THR-SAT
===============================================================

Teste le solveur de base (elagage par intervalles + Genome 1) sur 4 familles
algebriques / pathologiques encodees exactement en THR∘THR (UNSAT de justesse) :
  F1. PARITY_n        : somme(x_i) mod 2, max = 1, seuil T = 2
  F2. MOD3_n          : somme(x_i) mod 3, max = 2, seuil T = 3
  F3. WEIGHTED-PARITY : somme(a_i x_i) avec poids a_i in {1,2,3} et seuils
                        alternes (non symetrique, mais annulations exactes)
  F4. LARGE-WEIGHTS   : circuits aleatoires s = 2n a grands poids (R = 1000)

Objectif : mesurer exactement ou l'elagage par intervalles pur s'effondre.
"""
import json
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from thr2_lab import HERE, NF, solve_batch, top_max

OUT = sys.argv[1] if len(sys.argv) > 1 else HERE


def make_parity(n):
    """PARITY_n : g_k(x) = [sum x_i >= k], v_k = (-1)^(k-1), max = 1, T = 2."""
    W = np.ones((n, n), dtype=np.int64)
    t = np.arange(1, n + 1, dtype=np.int64)
    v = np.array([1 if k % 2 == 0 else -1 for k in range(n)], dtype=np.int64)
    T = 2  # max exact = 1 -> UNSAT de justesse
    return W, t, v, T


def make_mod3(n):
    """MOD3_n : g_k(x) = [sum x_i >= k], v_k = (+1, +1, -2, ...), max = 2, T = 3."""
    W = np.ones((n, n), dtype=np.int64)
    t = np.arange(1, n + 1, dtype=np.int64)
    pattern = [1, 1, -2]
    v = np.array([pattern[k % 3] for k in range(n)], dtype=np.int64)
    T = 3  # max exact = 2 -> UNSAT de justesse
    return W, t, v, T


def make_weighted_parity(n, seed):
    """Parite ponderee : L(x) = sum a_i x_i (a_i in {1,2,3}), s = n seuils alternes."""
    rng = np.random.default_rng(seed)
    a = rng.integers(1, 4, size=n, dtype=np.int64)
    s = n
    W = np.tile(a, (s, 1))
    max_L = int(a.sum())
    # Seuils regulierement espaces dans [1, max_L]
    t = np.linspace(1, max_L, s, dtype=np.int64)
    t = np.maximum.accumulate(t)
    v = np.array([1 if k % 2 == 0 else -1 for k in range(s)], dtype=np.int64)
    M, _ = top_max(W, t, v)
    return W, t, v, int(M) + 1


def make_large_weights(n, s, count, seed, R=1000, density=0.6):
    """Circuits aleatoires a grands poids R = 1000, UNSAT de justesse."""
    rng = np.random.default_rng(seed)
    Ws = np.zeros((count, s, n), np.int64)
    ts = np.zeros((count, s), np.int64)
    vs = np.zeros((count, s), np.int64)
    Ts = np.zeros(count, np.int64)
    for k in range(count):
        W = rng.integers(-R, R + 1, (s, n)) * (rng.random((s, n)) < density)
        t = np.floor(W.sum(axis=1) / 2).astype(np.int64) + rng.integers(-R // 10, R // 10 + 1, s)
        v = rng.integers(1, R + 1, s) * rng.choice([-1, 1], s)
        M, _ = top_max(W.astype(np.int64), t, v.astype(np.int64))
        Ws[k], ts[k], vs[k], Ts[k] = W, t, v, M + 1
    return Ws, ts, vs, Ts


def pack_single(items):
    Ws = np.stack([x[0] for x in items])
    ts = np.stack([x[1] for x in items])
    vs = np.stack([x[2] for x in items])
    Ts = np.array([x[3] for x in items], dtype=np.int64)
    return Ws, ts, vs, Ts


def main():
    t0 = time.time()
    # Charger Genome 1 si disponible
    res_path = os.path.join(HERE, "thr2_results.json")
    if os.path.exists(res_path):
        with open(res_path, "r", encoding="utf-8") as fh:
            g1 = np.array(json.load(fh)["genome1"], dtype=np.float64)
    else:
        g1 = np.array([-1.04, -0.06, 0.09, 1.43, -0.57, 0.84, 0.32], dtype=np.float64)

    methods = {
        "Ordre fixe": np.zeros(NF),
        "Heuristique f7": np.eye(NF)[6],
        "Genome 1": g1,
    }

    ns = [10, 12, 14, 16, 18]
    families = {}
    print("== Generation du banc d'essai adversarial structure ==", flush=True)
    families["F1: PARITY_n"] = {n: pack_single([make_parity(n)]) for n in ns}
    families["F2: MOD3_n"] = {n: pack_single([make_mod3(n)]) for n in ns}
    families["F3: WEIGHTED-PARITY"] = {
        n: pack_single([make_weighted_parity(n, seed=300 + 17 * n + k) for k in range(6)])
        for n in ns
    }
    families["F4: LARGE-WEIGHTS (s=2n, R=1000)"] = {
        n: make_large_weights(n, 2 * n, count=6, seed=700 + 19 * n, R=1000)
        for n in ns
    }

    report = {"ns": ns, "families": {}}
    for fname, sets in families.items():
        print(f"\n== {fname} ==", flush=True)
        report["families"][fname] = {}
        for mname, th in methods.items():
            logs = []
            for n in ns:
                Ws, ts, vs, Ts = sets[n]
                st, nodes = solve_batch(Ws, ts, vs, Ts, th, 1 << (n + 2))
                assert (st == 0).all(), f"Erreur statut sur {fname} n={n}: {st}"
                logs.append(math.log2(float(nodes.mean())))
            alpha = float(np.polyfit(ns, logs, 1)[0])
            report["families"][fname][mname] = {"alpha": alpha, "log2_nodes": logs}
            print(f"   {mname:<18} alpha = {alpha:.3f}   log2(noeuds) = "
                  + " ".join(f"{x:5.1f}" for x in logs), flush=True)

    print(f"\nTemps total Partie 1 : {time.time() - t0:.1f} s", flush=True)
    with open(os.path.join(HERE, "adversarial_results.json"), "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 4, figsize=(18, 4.4), sharey=True)
    for ax, (fname, res) in zip(axes, report["families"].items()):
        ax.plot(ns, [n + 1 for n in ns], "k--", label="Force brute (α=1.00)")
        for mname, r in res.items():
            ax.plot(ns, r["log2_nodes"], "o-", label=f"{mname} (α={r['alpha']:.2f})")
        ax.set_title(fname, fontsize=9.5)
        ax.set_xlabel("n")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7.5)
    axes[0].set_ylabel("log2(nœuds explorés)")
    fig.tight_layout()
    out_png = os.path.join(OUT, "adversarial_suite.png")
    fig.savefig(out_png, dpi=130)
    print(f"Graphique : {out_png}", flush=True)


if __name__ == "__main__":
    main()
