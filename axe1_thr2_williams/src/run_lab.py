"""
Experience complete :
  A. evolution d'une regle de branchement sur instances aleatoires
  B. validation sur instances nouvelles, n jusqu'a 26
  C. course aux armements : un "adversaire" fabrique les pires instances
     pour le solveur, le solveur re-evolue, l'adversaire recommence.
     -> estime le PIRE CAS (c'est lui qui compte pour la theorie), pas la moyenne.
"""
from pathlib import Path
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from thr2_lab import (HERE, NF, FEATURE_NAMES, evolve, make_instances, measure,
                      solve, top_max)

AXE_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = str(AXE_ROOT / "results")
FIGURES_DIR = str(AXE_ROOT / "figures")
WITNESS_DIR = str(AXE_ROOT / "witnesses")
for _d in (RESULTS_DIR, FIGURES_DIR, WITNESS_DIR):
    os.makedirs(_d, exist_ok=True)
OUT = sys.argv[1] if len(sys.argv) > 1 else FIGURES_DIR
LINES = []


def log(msg=""):
    print(msg, flush=True)
    LINES.append(msg)


def random_circuit(n, s, rng, R=3, density=0.6):
    W = (rng.integers(-R, R + 1, (s, n)) * (rng.random((s, n)) < density)).astype(np.int64)
    t = (np.floor(W.sum(axis=1) / 2).astype(np.int64) + rng.integers(-1, 2, s)).astype(np.int64)
    v = (rng.integers(1, R + 1, s) * rng.choice([-1, 1], s)).astype(np.int64)
    return W, t, v


def harden(n, s, theta, steps, rng, R=3):
    """Recherche locale d'une instance UNSAT-de-justesse maximisant les noeuds du solveur."""
    W, t, v = random_circuit(n, s, rng)
    T = top_max(W, t, v)[0] + 1
    best = solve(W, t, v, T, theta, 1 << (n + 1))[1]
    for _ in range(steps):
        W2, t2, v2 = W.copy(), t.copy(), v.copy()
        r = rng.random()
        g = rng.integers(s)
        if r < 0.6:
            W2[g, rng.integers(n)] = rng.integers(-R, R + 1)
        elif r < 0.8:
            t2[g] += rng.choice([-1, 1])
        else:
            v2[g] = rng.integers(1, R + 1) * rng.choice([-1, 1])
        T2 = top_max(W2, t2, v2)[0] + 1
        nd = solve(W2, t2, v2, T2, theta, 1 << (n + 1))[1]
        if nd >= best:
            W, t, v, T, best = W2, t2, v2, T2, nd
    return W, t, v, T


def hard_sets(theta, ns, count, steps, rng):
    sets = {}
    for n in ns:
        inst = [harden(n, n, theta, steps, rng) for _ in range(count)]
        sets[n] = tuple(np.stack([x[k] for x in inst]) for k in range(3)) + (
            np.array([x[3] for x in inst], np.int64),)
    return sets


def report(name, theta, sets, ns):
    a, logs = measure(theta, sets, ns)
    log(f"   {name:<34} alpha = {a:.3f}   log2(noeuds) = " + " ".join(f"{x:5.1f}" for x in logs))
    return {"alpha": float(a), "log2_nodes": [float(x) for x in logs]}


def main():
    rng = np.random.default_rng(2026)
    t0 = time.time()
    res = {}

    # ---------------- A. evolution sur instances aleatoires ----------------
    train_ns = [10, 12, 14, 16, 18]
    log("== A. Evolution sur instances aleatoires (UNSAT de justesse, s = n) ==")
    train = {n: make_instances(n, n, 12, seed=1000 + n) for n in train_ns}
    th1, a1 = evolve(train, train_ns, gens=15, pop=16, rng=rng, log=log)
    log(f"   genome 1 : {np.round(th1, 2).tolist()}   ({time.time() - t0:.0f} s)\n")

    base = {"Ordre fixe (elagage seul)": np.zeros(NF),
            "Heuristique f2 (a la main)": np.eye(NF)[1],
            "Heuristique f7 (a la main)": np.eye(NF)[6],
            "Evolue (genome 1)": th1}

    # ---------------- B. validation ----------------
    val_ns = [12, 14, 16, 18, 20, 22, 24, 26]
    for label, mult in [("aleatoires s = n", 1), ("aleatoires s = 2n", 2)]:
        log(f"== B. Validation, instances NOUVELLES {label} ==")
        sets = {n: make_instances(n, mult * n, 12, seed=9000 + 17 * n + mult) for n in val_ns}
        res[label] = {"ns": val_ns, "m": {k: report(k, th, sets, val_ns) for k, th in base.items()}}
        log(f"   ({time.time() - t0:.0f} s)\n")

    # ---------------- C. course aux armements ----------------
    adv_ns = [10, 12, 14, 16]
    log("== C1. L'adversaire fabrique les pires instances contre le genome 1 ==")
    H1 = hard_sets(th1, adv_ns, count=6, steps=250, rng=rng)
    r_c1 = {k: report(k, th, H1, adv_ns) for k, th in base.items()}
    log(f"   ({time.time() - t0:.0f} s)\n")

    log("== C2. Le solveur re-evolue sur ces instances difficiles ==")
    th2, _ = evolve(H1, adv_ns, gens=12, pop=16, rng=rng, log=log)
    log(f"   genome 2 : {np.round(th2, 2).tolist()}")
    r_c2 = report("Genome 2 sur instances H1", th2, H1, adv_ns)
    log("")

    log("== C3. L'adversaire contre-attaque contre le genome 2 ==")
    H2 = hard_sets(th2, adv_ns, count=6, steps=250, rng=rng)
    r_c3 = {"Genome 1 sur H2": report("Genome 1 sur H2", th1, H2, adv_ns),
            "Genome 2 sur H2": report("Genome 2 sur H2", th2, H2, adv_ns),
            "Ordre fixe sur H2": report("Ordre fixe sur H2", np.zeros(NF), H2, adv_ns)}
    rand_adv = {n: make_instances(n, n, 12, seed=4000 + n) for n in adv_ns}
    r_rand = report("Genome 1 sur aleatoires (memes n)", th1, rand_adv, adv_ns)
    log(f"\nTemps total : {time.time() - t0:.0f} s")

    res["arms_race"] = {"ns": adv_ns, "C1": r_c1, "C2": r_c2, "C3": r_c3, "random": r_rand}
    with open(os.path.join(RESULTS_DIR, "thr2_results.json"), "w", encoding="utf-8") as fh:
        json.dump({"genome1": th1.tolist(), "genome2": th2.tolist(), "features": FEATURE_NAMES,
                   "results": res, "log": LINES}, fh, indent=2)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(17, 4.8))
    for ax, label in zip(axes[:2], ["aleatoires s = n", "aleatoires s = 2n"]):
        ns = res[label]["ns"]
        ax.plot(ns, [n + 1 for n in ns], "k--", label="Force brute (α = 1)")
        for k, r in res[label]["m"].items():
            ax.plot(ns, r["log2_nodes"], "o-", label=f"{k} (α = {r['alpha']:.2f})")
        ax.set_title(f"Instances {label}")
        ax.set_xlabel("n")
        ax.set_ylabel("log2(nœuds)")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7.5)
    ax = axes[2]
    ax.plot(adv_ns, [n + 1 for n in adv_ns], "k--", label="Force brute (α = 1)")
    curves = [("Génome 1 / aléatoires", r_rand), ("Génome 1 / pires (H1)", r_c1["Evolue (genome 1)"]),
              ("Génome 2 / H1", r_c2), ("Génome 2 / pires (H2)", r_c3["Genome 2 sur H2"])]
    for k, r in curves:
        ax.plot(adv_ns, r["log2_nodes"], "o-", label=f"{k} (α = {r['alpha']:.2f})")
    ax.set_title("Course aux armements : pire cas")
    ax.set_xlabel("n")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7.5)
    fig.tight_layout()
    out = os.path.join(OUT, "thr2_evolution.png")
    fig.savefig(out, dpi=130)
    log(f"Graphique : {out}")


if __name__ == "__main__":
    main()
