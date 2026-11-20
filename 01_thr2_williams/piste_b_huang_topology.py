"""
Piste B : Topologie de l'hypercube et spectre de Huang (2019) sur THR∘THR
=========================================================================

Explore deux invariants topologiques de l'hypercube Q_n = {0,1}^n :
  1. L'influence totale Inf(f) = E_x[s(f, x)] et la borne exacte de frontiere :
         Inf(C) <= s * Inf(MAJ_n) ~ 0.798 * s * sqrt(n)
     qui prouve en 1 ligne s*_{THR∘THR}(PARITY_n) >= Omega(sqrt(n)).
  2. La matrice d'adjacence signee de Huang (2019) A_n in {-1, 0, +1}^(2^n x 2^n)
     (valeurs propres +/- sqrt(n)) :
       - lambda_H(f) = max(lambda_max(A_n[S_1]), lambda_max(A_n[S_0]))
       - s_max(f)    = sensibilite maximale max_x |{i : f(x) != f(x ^ e_i)}|
       - lambda_bd(f)= rayon spectral du graphe des aretes frontieres E_d(f)
                       (aretes ou f(x) != f(x ^ e_i)).

Experiences :
  1. Sur l'atlas exact des 222 classes NPN a n = 4 variables : comparaison de
     Inf(f), s_max(f), lambda_H(f) et lambda_bd(f) selon la taille exacte s*.
  2. Sur n = 6, 8, 10 : comparaison entre circuits THR∘THR de taille s = 1..2n,
     PARITY_n, Produit Scalaire IP_n et fonctions booleennes aleatoires.
"""
import json
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from thr2_lab import HERE, make_instances

OUT = sys.argv[1] if len(sys.argv) > 1 else HERE


def huang_matrix(n):
    """Matrice d'adjacence signee de Huang A_n (taille 2^n x 2^n, valeurs propres +/- sqrt(n))."""
    A = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=np.float64)
    for k in range(2, n + 1):
        m = 1 << (k - 1)
        I = np.eye(m, dtype=np.float64)
        A = np.block([[A, I], [I, -A]])
    return A


def hypercube_adj(n):
    """Matrice d'adjacence standard de l'hypercube Q_n."""
    N = 1 << n
    Q = np.zeros((N, N), dtype=np.float64)
    for x in range(N):
        for i in range(n):
            Q[x, x ^ (1 << i)] = 1.0
    return Q


def eval_circuit_tt(W, t, v, T):
    """Table de verite (vecteur 0/1 de longueur 2^n) d'un circuit THR∘THR."""
    s, n = W.shape
    N = 1 << n
    idx = np.arange(N, dtype=np.int64)[:, None]
    bits = ((idx >> np.arange(n, dtype=np.int64)[None, :]) & 1).astype(np.int64)
    sums = bits @ W.T  # (N, s)
    top = ((sums >= t[None, :]).astype(np.int64) @ v)
    return (top >= T).astype(np.int8)


def topological_invariants(tt, n, A_huang, Q_adj):
    """
    Calcule les 4 invariants topologiques d'une fonction booleenne tt in {0,1}^(2^n) :
      - Inf(f)      : sensibilite moyenne (influence totale)
      - s_max(f)    : sensibilite maximale
      - lambda_H(f) : max valeur propre de A_huang restreinte a S_1 ou S_0
      - lambda_bd(f): rayon spectral du graphe des aretes frontieres (ou f bascule)
    """
    N = 1 << n
    # Sensibilite en chaque sommet x
    sens = np.zeros(N, dtype=np.int32)
    for i in range(n):
        flip = tt ^ tt[np.arange(N) ^ (1 << i)]
        sens += flip
    inf_val = float(sens.mean())
    s_max = int(sens.max())

    # Valeur propre maximale de Huang sur S_1 et S_0
    idx1 = np.nonzero(tt == 1)[0]
    idx0 = np.nonzero(tt == 0)[0]
    l1 = float(np.linalg.eigvalsh(A_huang[np.ix_(idx1, idx1)]).max()) if len(idx1) > 0 else 0.0
    l0 = float(np.linalg.eigvalsh(A_huang[np.ix_(idx0, idx0)]).max()) if len(idx0) > 0 else 0.0
    lam_h = max(l1, l0)

    # Graphe des aretes frontieres : B[x, y] = Q[x, y] * 1[tt[x] != tt[y]]
    diff_mask = (tt[:, None] != tt[None, :]).astype(np.float64)
    B = Q_adj * diff_mask
    lam_bd = float(np.linalg.eigvalsh(B).max())

    return inf_val, s_max, lam_h, lam_bd


def main():
    t0 = time.time()

    # ----------------------------------------------------------------------- #
    #  1. Analyse sur les 222 classes NPN de l'atlas exact n = 4
    # ----------------------------------------------------------------------- #
    print("== 1. Invariants topologiques (Huang 2019 + Frontiere) sur les 222 classes NPN (n = 4) ==", flush=True)
    atlas = np.load(os.path.join(HERE, "atlas_n4.npz"))
    s_star = atlas["s_star"]
    reps = atlas["reps"]
    A4 = huang_matrix(4)
    Q4 = hypercube_adj(4)

    by_s = {k: {"inf": [], "smax": [], "lam_h": [], "lam_bd": []} for k in (1, 2, 3, 4)}
    for r in reps:
        tt = np.array([(int(r) >> x) & 1 for x in range(16)], dtype=np.int8)
        inf_v, smax_v, lh_v, lbd_v = topological_invariants(tt, 4, A4, Q4)
        k = int(s_star[r])
        by_s[k]["inf"].append(inf_v)
        by_s[k]["smax"].append(smax_v)
        by_s[k]["lam_h"].append(lh_v)
        by_s[k]["lam_bd"].append(lbd_v)

    print(f"   {'s*':>3} | {'Classes':>7} | {'Inf moy (max)':>14} | {'s_max moy (max)':>15} | "
          f"{'Huang λ_H moy (min)':>19} | {'Frontiere λ_bd moy (max)':>24}", flush=True)
    print("   " + "-" * 94, flush=True)
    n4_summary = {}
    for k in (1, 2, 3, 4):
        d = by_s[k]
        row = {
            "n_classes": len(d["inf"]),
            "inf_mean": float(np.mean(d["inf"])), "inf_max": float(np.max(d["inf"])),
            "smax_mean": float(np.mean(d["smax"])), "smax_max": int(np.max(d["smax"])),
            "lam_h_mean": float(np.mean(d["lam_h"])), "lam_h_min": float(np.min(d["lam_h"])),
            "lam_bd_mean": float(np.mean(d["lam_bd"])), "lam_bd_max": float(np.max(d["lam_bd"])),
        }
        n4_summary[k] = row
        print(f"   {k:>3} | {row['n_classes']:>7} | {row['inf_mean']:6.3f} ({row['inf_max']:5.2f}) | "
              f"{row['smax_mean']:7.2f} ({row['smax_max']:5d}) | "
              f"{row['lam_h_mean']:10.3f} ({row['lam_h_min']:5.2f}) | "
              f"{row['lam_bd_mean']:15.3f} ({row['lam_bd_max']:5.2f})", flush=True)

    # ----------------------------------------------------------------------- #
    #  2. Passage a n = 8 (N = 256 sommets) : THR∘THR(s) vs PARITY_8, IP_8, Aleatoire
    # ----------------------------------------------------------------------- #
    n = 8
    N = 1 << n
    print(f"\n== 2. Topologie de Q_{n} (N={N} sommets) : circuits THR∘THR(s) vs fonctions dures ==", flush=True)
    A8 = huang_matrix(n)
    Q8 = hypercube_adj(n)
    rng = np.random.default_rng(2026)

    # Borne exacte Inf(MAJ_n) pour n = 8
    inf_maj8 = n * math.comb(n - 1, (n - 1) // 2) / (1 << (n - 1))
    print(f"   (Borne theorique 1 porte : Inf(MAJ_8) = {inf_maj8:.3f}, sqrt(n) = {math.sqrt(n):.3f})\n", flush=True)
    print(f"   {'Famille (n = 8)':<30} | {'Inf(f)':>7} | {'s_max':>6} | {'Huang λ_H':>9} | {'Frontiere λ_bd':>14}",
          flush=True)
    print("   " + "-" * 78, flush=True)

    n8_rows = {}
    # Circuits THR∘THR equilibres (seuil median pour eviter les fonctions constantes)
    for s in (1, 2, 4, 8, 16, 32):
        vals = []
        for k_inst in range(8):
            W = (rng.integers(-3, 4, (s, n)) * (rng.random((s, n)) < 0.7)).astype(np.int64)
            t = np.floor(W.sum(axis=1) / 2).astype(np.int64)
            v = (rng.integers(1, 4, s) * rng.choice([-1, 1], s)).astype(np.int64)
            # Evaluer les sommes du sommet pour prendre le seuil median (fonction equilibree)
            idx = np.arange(N, dtype=np.int64)[:, None]
            bits = ((idx >> np.arange(n, dtype=np.int64)[None, :]) & 1).astype(np.int64)
            top = ((bits @ W.T >= t[None, :]).astype(np.int64) @ v)
            T = int(np.median(top))
            tt = (top >= T).astype(np.int8)
            if tt.sum() == 0 or tt.sum() == N:
                tt = (top > T - 1).astype(np.int8)
            vals.append(topological_invariants(tt, n, A8, Q8))
        row = {
            "inf": float(np.mean([x[0] for x in vals])),
            "smax": float(np.mean([x[1] for x in vals])),
            "lam_h": float(np.mean([x[2] for x in vals])),
            "lam_bd": float(np.mean([x[3] for x in vals])),
        }
        name = f"THR∘THR aleatoire (s={s:>2})"
        n8_rows[name] = row
        print(f"   {name:<30} | {row['inf']:>7.3f} | {row['smax']:>6.2f} | {row['lam_h']:>9.3f} | {row['lam_bd']:>14.3f}",
              flush=True)

    # Fonctions de reference dures sur n = 8
    tt_parity = np.array([bin(x).count("1") & 1 for x in range(N)], dtype=np.int8)
    # Produit scalaire IP_8(x_L, x_R) = sum_{i=0..3} x_i * x_{i+4} mod 2
    tt_ip = np.array([bin((x & 0xF) & (x >> 4)).count("1") & 1 for x in range(N)], dtype=np.int8)
    # Fonction booleenne uniforme aleatoire (moyenne sur 8 tirages)
    rand_vals = [topological_invariants(rng.integers(0, 2, size=N, dtype=np.int8), n, A8, Q8) for _ in range(8)]

    refs = [
        ("Produit scalaire IP_8", topological_invariants(tt_ip, n, A8, Q8)),
        ("Fonction aleatoire uniforme", tuple(float(np.mean([x[i] for x in rand_vals])) for i in range(4))),
        ("PARITY_8", topological_invariants(tt_parity, n, A8, Q8)),
    ]
    for name, (inf_v, smax_v, lh_v, lbd_v) in refs:
        n8_rows[name] = {"inf": float(inf_v), "smax": float(smax_v), "lam_h": float(lh_v), "lam_bd": float(lbd_v)}
        print(f"   {name:<30} | {inf_v:>7.3f} | {smax_v:>6.2f} | {lh_v:>9.3f} | {lbd_v:>14.3f}", flush=True)

    with open(os.path.join(HERE, "piste_b_huang_results.json"), "w", encoding="utf-8") as fh:
        json.dump({"n4_npn_atlas": n4_summary, "n8_topology": n8_rows, "inf_maj8": inf_maj8}, fh, indent=2)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 4.8))

    # Panel 1 : Atlas n=4
    ks = [1, 2, 3, 4]
    ax1.plot(ks, [n4_summary[k]["inf_mean"] for k in ks], "o-", label="Influence moyenne Inf(f)")
    ax1.plot(ks, [n4_summary[k]["lam_bd_mean"] for k in ks], "s-", label="Rayon spectral frontière λ_bd(f)")
    ax1.plot(ks, [n4_summary[k]["lam_h_mean"] for k in ks], "^--", label="Valeur propre Huang λ_H(f)")
    ax1.set_xticks(ks)
    ax1.set_xlabel("Taille minimale exacte s* (n = 4, 222 classes NPN)")
    ax1.set_ylabel("Invariant topologique sur Q_4")
    ax1.set_title("Topologie de Q_4 vs taille minimale THR∘THR")
    ax1.grid(alpha=0.3)
    ax1.legend(fontsize=8)

    # Panel 2 : n = 8 (THR∘THR(s) vs fonctions dures)
    s_vals = [1, 2, 4, 8, 16, 32]
    ax2.plot(s_vals, [n8_rows[f"THR∘THR aleatoire (s={s:>2})"]["inf"] for s in s_vals], "o-",
             color="#1f77b4", label="THR∘THR(s) : Influence Inf(f)")
    ax2.plot(s_vals, [n8_rows[f"THR∘THR aleatoire (s={s:>2})"]["lam_bd"] for s in s_vals], "s-",
             color="#2ca02c", label="THR∘THR(s) : Rayon frontière λ_bd(f)")
    ax2.axhline(n8_rows["Produit scalaire IP_8"]["lam_bd"], color="#ff7f0e", ls="--",
                label=f"IP_8 : λ_bd = {n8_rows['Produit scalaire IP_8']['lam_bd']:.2f}")
    ax2.axhline(n8_rows["Fonction aleatoire uniforme"]["lam_bd"], color="#d62728", ls=":",
                label=f"Aléatoire : λ_bd = {n8_rows['Fonction aleatoire uniforme']['lam_bd']:.2f}")
    ax2.axhline(n8_rows["PARITY_8"]["lam_bd"], color="black", ls="-.",
                label=f"PARITY_8 : λ_bd = {n8_rows['PARITY_8']['lam_bd']:.2f}")
    ax2.set_xscale("log", base=2)
    ax2.set_xlabel("Nombre de portes cachées s (n = 8, échelle log2)")
    ax2.set_ylabel("Invariant topologique sur Q_8")
    ax2.set_title("Séparation topologique sur Q_8 : THR∘THR(s) vs fonctions dures")
    ax2.grid(alpha=0.3, which="both")
    ax2.legend(fontsize=7.5)

    fig.tight_layout()
    out_png = os.path.join(OUT, "piste_b_huang_topology.png")
    fig.savefig(out_png, dpi=130)
    print(f"\nTemps total Piste B : {time.time() - t0:.2f} s", flush=True)
    print(f"Graphique : {out_png}", flush=True)


if __name__ == "__main__":
    main()
