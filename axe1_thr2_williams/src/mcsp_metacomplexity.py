"""
Piste 3 : Meta-complexite (MCSP) et Barriere des Preuves Naturelles en miniature
================================================================================

Exploite l'atlas exact `atlas_n4.npz` des 65,536 fonctions booleennes a n = 4
variables pour :
  1. Calculer 5 invariants mathematiques / spectraux pour chacune des 65,536
     tables de verite :
       - I1 : Norme spectrale de Walsh-Hadamard L1(f^) = sum_S |f^(S)|
              (complexite spectrale de Bruck-Smolensky)
       - I2 : Energie de Chow C(f) = sum_{|S| <= 1} |f^(S)|^2
              (par le theoreme de Chow, C(f) caracterise les portes s* = 1)
       - I3 : Influence totale / sensibilite moyenne Inf(f)
              (perimetre de la frontiere sur l'hypercube {0,1}^4)
       - I4 : Degre algebrique deg_{GF(2)}(f) et nombre de monomes ANF
       - I5 : Entropie de sous-cubes H_sub(f) (proxy de Kolmogorov)
  2. Evaluer le probleme MCSP (Minimum Circuit Size Problem) : quelle est la
     capacite d'un predicat "naturel" (au sens de Razborov-Rudich : constructif
     et large) a separer les fonctions faciles (s* <= 2) des fonctions dures ?
  3. Identifier les "Fonctions Cameleons" (temoins concrets de la barriere des
     Preuves Naturelles) : des fonctions qui possedent un petit circuit THR∘THR
     (s* = 2) mais dont le spectre de Fourier, l'influence et le degre GF(2)
     imitent ceux des fonctions les plus dures (s* maximal) !
"""
from pathlib import Path
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
AXE_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = str(AXE_ROOT / "results")
FIGURES_DIR = str(AXE_ROOT / "figures")
WITNESS_DIR = str(AXE_ROOT / "witnesses")
for _d in (RESULTS_DIR, FIGURES_DIR, WITNESS_DIR):
    os.makedirs(_d, exist_ok=True)
OUT = sys.argv[1] if len(sys.argv) > 1 else FIGURES_DIR


def compute_all_invariants():
    """Calcule les 5 invariants de meta-complexite pour les 65,536 fonctions de n = 4."""
    N = 65536
    # Matrice des tables de verite : bits[f, x] in {0, 1}, pm[f, x] = (-1)^bits[f, x]
    f_idx = np.arange(N, dtype=np.int32)[:, None]
    x_idx = np.arange(16, dtype=np.int32)[None, :]
    bits = ((f_idx >> x_idx) & 1).astype(np.int32)
    pm = 1.0 - 2.0 * bits.astype(np.float64)

    # 1) Transformee de Walsh-Hadamard rapide sur {0,1}^4
    H = pm.copy()
    for i in range(4):
        H = H.reshape(N, -1, 2, 1 << i)
        a = H[:, :, 0, :].copy()
        b = H[:, :, 1, :].copy()
        H[:, :, 0, :] = a + b
        H[:, :, 1, :] = a - b
        H = H.reshape(N, 16)
    H /= 16.0  # f^(S) in [-1, 1], sum_S f^(S)^2 = 1 (Parseval)

    # I1 : Norme spectrale L1 = sum_S |f^(S)|
    spectral_l1 = np.abs(H).sum(axis=1)

    # I2 : Energie de Chow = sum_{|S| <= 1} f^(S)^2  (proportion d'energie affine)
    deg1_mask = np.array([bin(s).count("1") <= 1 for s in range(16)], dtype=bool)
    chow_energy = (H[:, deg1_mask] ** 2).sum(axis=1)

    # I3 : Influence totale Inf(f) = sum_{i=0..3} Pr_x[f(x) != f(x ^ 2^i)]
    # Par identite de Fourier : Inf(f) = sum_S |S| * f^(S)^2
    weights_S = np.array([bin(s).count("1") for s in range(16)], dtype=np.float64)
    influence = ((H ** 2) * weights_S[None, :]).sum(axis=1)

    # I4 : Transformee de Mobius sur GF(2) -> Polynome ANF (Reed-Muller)
    anf = bits.copy()
    for i in range(4):
        anf = anf.reshape(N, -1, 2, 1 << i)
        anf[:, :, 1, :] ^= anf[:, :, 0, :]
        anf = anf.reshape(N, 16)
    anf_monomials = anf.sum(axis=1)
    anf_deg = np.max(anf * weights_S[None, :].astype(np.int32), axis=1)

    # I5 : Diversite de sous-cubes (nombre de motifs 4-bits distincts sur les 6 paires de variables)
    # Proxy de complexite de Kolmogorov / compressibilite de la table de verite
    sub_diversity = np.zeros(N, dtype=np.float64)
    # 4 sous-cubes obtenus en fixant chaque variable x_i in {0, 1} (8 bits chacun)
    for i in range(4):
        m0 = [x for x in range(16) if not ((x >> i) & 1)]
        m1 = [x for x in range(16) if (x >> i) & 1]
        w8 = (1 << np.arange(8, dtype=np.int32))
        s0 = bits[:, m0] @ w8
        s1 = bits[:, m1] @ w8
        # Distance de Hamming entre les deux moities + asymetrie
        sub_diversity += np.abs(bits[:, m0] - bits[:, m1]).sum(axis=1) / 4.0

    return spectral_l1, chow_energy, influence, anf_deg, anf_monomials, sub_diversity


def main():
    t0 = time.time()
    atlas_path = os.path.join(WITNESS_DIR, "atlas_n4.npz")
    if not os.path.exists(atlas_path):
        import exact_thr2_atlas
        exact_thr2_atlas.main()
    data = np.load(atlas_path)
    s_star = data["s_star"]
    reps = data["reps"]

    print("== 1. Calcul des 5 invariants de meta-complexite sur les 65,536 tables de verite ==", flush=True)
    l1, chow, inf, deg, mon, subdiv = compute_all_invariants()

    max_s = int(s_star.max())
    print("\n== 2. Profil moyen des invariants selon la taille minimale exacte s*_{THR∘THR}(f) ==", flush=True)
    print(f"   {'s*':>3} | {'Nb fonct.':>9} | {'Walsh L1':>9} | {'Chow E_1':>9} | "
          f"{'Influence':>9} | {'Deg GF(2)':>9} | {'Monomes ANF':>11}", flush=True)
    print("   " + "-" * 74, flush=True)

    stats_by_s = {}
    for k in range(1, max_s + 1):
        mask = (s_star == k)
        cnt = int(mask.sum())
        if cnt == 0:
            continue
        row = {
            "count": cnt,
            "walsh_l1_mean": float(l1[mask].mean()),
            "walsh_l1_max": float(l1[mask].max()),
            "chow_mean": float(chow[mask].mean()),
            "chow_min": float(chow[mask].min()),
            "influence_mean": float(inf[mask].mean()),
            "deg_gf2_mean": float(deg[mask].mean()),
            "anf_monomials_mean": float(mon[mask].mean()),
        }
        stats_by_s[k] = row
        print(f"   {k:>3} | {cnt:>9,} | {row['walsh_l1_mean']:>9.3f} | {row['chow_mean']:>9.3f} | "
              f"{row['influence_mean']:>9.3f} | {row['deg_gf2_mean']:>9.3f} | {row['anf_monomials_mean']:>11.3f}",
              flush=True)

    # ----------------------------------------------------------------------- #
    #  3. Test des "Preuves Naturelles" (Razborov-Rudich) et MCSP :
    #     Quelle fraction des fonctions s* >= 3 un critere spectral/combinatoire
    #     peut-il certifier SANS aucun faux positif dans s* <= 2 ?
    # ----------------------------------------------------------------------- #
    print("\n== 3. Analyse de la barriere des Preuves Naturelles (Razborov-Rudich) ==", flush=True)
    mask_le2 = (s_star <= 2)
    mask_ge3 = (s_star >= 3)
    max_l1_s1 = float(l1[s_star == 1].max())
    max_l1_s2 = float(l1[mask_le2].max())
    min_chow_s1 = float(chow[s_star == 1].min())
    min_chow_s2 = float(chow[mask_le2].min())
    max_inf_s2 = float(inf[mask_le2].max())

    print(f"   Borne spectrale exacte pour s* = 1 : max Walsh L1 = {max_l1_s1:.3f}, min Chow E_1 = {min_chow_s1:.3f}",
          flush=True)
    print(f"   Borne spectrale exacte pour s* <= 2 : max Walsh L1 = {max_l1_s2:.3f}, min Chow E_1 = {min_chow_s2:.3f}, "
          f"max Influence = {max_inf_s2:.3f}", flush=True)

    cert_l1_over_s1 = float((l1[s_star >= 2] > max_l1_s1 + 1e-9).mean()) * 100
    cert_l1_over_s2 = float((l1[mask_ge3] > max_l1_s2 + 1e-9).mean()) * 100
    print(f"   -> Le critere naturel 'Walsh L1 > {max_l1_s1:.2f}' certifie s* >= 2 sur {cert_l1_over_s1:.1f}% des fonctions s* >= 2.",
          flush=True)
    print(f"   -> En revanche, 'Walsh L1 > {max_l1_s2:.2f}' ne certifie s* >= 3 que sur {cert_l1_over_s2:.1f}% des fonctions s* >= 3 !",
          flush=True)

    # ----------------------------------------------------------------------- #
    #  4. Identification des "Fonctions Cameleons" (s* = 2 mais profil s* = max)
    # ----------------------------------------------------------------------- #
    print("\n== 4. Decouverte des 'Fonctions Cameleons' (s* = 2 imitant les fonctions les plus dures) ==", flush=True)
    # Chercher des representants NPN avec s* = 2 ayant Chow E_1 = 0 (aucune correlation lineaire,
    # comme PARITY_4 !), degre GF(2) = 4, et norme spectrale L1 ou Influence maximale !
    chameleons = []
    for r in reps:
        if s_star[r] == 2 and chow[r] < 1e-9 and deg[r] == 4:
            chameleons.append((float(l1[r]), float(inf[r]), int(r)))
    chameleons.sort(reverse=True)
    for wl1, winf, r in chameleons[:5]:
        print(f"   Cameleon NPN 0x{r:04x} : taille réelle s* = 2, MAIS Chow E_1 = 0.000 (orthogonal a toute fonction lineaire), "
              f"deg_GF2 = 4, Walsh L1 = {wl1:.3f}, Influence = {winf:.2f}", flush=True)

    # ----------------------------------------------------------------------- #
    #  5. Evaluation d'un predicteur MCSP (classification s* par les invariants)
    # ----------------------------------------------------------------------- #
    # Regle lineaire multi-invariants par moindres carres pour predire s*(f)
    X = np.column_stack([np.ones(65536), l1, chow, inf, deg, mon, subdiv])
    y = s_star.astype(np.float64)
    beta, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    y_pred = X @ beta
    corr = float(np.corrcoef(y, y_pred)[0, 1])
    acc_exact = float((np.rint(y_pred) == y).mean()) * 100
    print(f"\n== 5. Predicteur MCSP spectral/combinatoire ==", flush=True)
    print(f"   Correlation r(s*, s*_pred) = {corr:.3f} | Precision exacte d'arrondi = {acc_exact:.1f}%", flush=True)

    with open(os.path.join(RESULTS_DIR, "mcsp_results.json"), "w", encoding="utf-8") as fh:
        json.dump({
            "stats_by_s": stats_by_s,
            "max_l1_s1": max_l1_s1,
            "max_l1_s2": max_l1_s2,
            "cert_l1_over_s1_pct": cert_l1_over_s1,
            "cert_l1_over_s2_pct": cert_l1_over_s2,
            "chameleons": [
                {"tt_hex": f"0x{r:04x}", "walsh_l1": wl1, "influence": winf}
                for wl1, winf, r in chameleons[:8]
            ],
            "mcsp_corr": corr,
            "mcsp_acc_pct": acc_exact,
        }, fh, indent=2)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(17, 4.8))

    # Panel 1 : Distribution exacte de s*(f) sur les 65,536 fonctions et 222 classes NPN
    ks = sorted(stats_by_s.keys())
    counts = [stats_by_s[k]["count"] for k in ks]
    ax = axes[0]
    bars = ax.bar([str(k) for k in ks], counts, color=["#2ca02c", "#1f77b4", "#ff7f0e", "#d62728"][:len(ks)])
    for bar, c in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width() / 2, c + 800, f"{c:,}\n({100*c/65536:.1f}%)",
                ha="center", va="bottom", fontsize=8.5)
    ax.set_ylim(0, max(counts) * 1.18)
    ax.set_xlabel("Taille minimale exacte s* (nombre de portes THR cachées)")
    ax.set_ylabel("Nombre de fonctions booléennes (sur 65 536)")
    ax.set_title("Atlas exact THR∘THR pour n = 4 (65 536 fonctions)")
    ax.grid(alpha=0.3, axis="y")

    # Panel 2 : Nuage (Energie de Chow E_1 vs Norme spectrale Walsh L1) sur les 222 classes NPN
    ax = axes[1]
    colors = {1: "#2ca02c", 2: "#1f77b4", 3: "#ff7f0e", 4: "#d62728", 5: "#9467bd"}
    for k in ks:
        r_k = [r for r in reps if s_star[r] == k]
        ax.scatter(chow[r_k], l1[r_k], c=colors.get(k, "black"), label=f"s* = {k} ({len(r_k)} classes NPN)",
                   alpha=0.8, s=36, edgecolors="k", linewidths=0.3)
    if chameleons:
        r_cham = [r for _, _, r in chameleons[:3]]
        ax.scatter(chow[r_cham], l1[r_cham], facecolors="none", edgecolors="red", s=130, linewidths=1.8,
                   label="Caméléons (s*=2 mais E₁=0)")
    ax.set_xlabel("Énergie affine de Chow E₁(f) = ∑_{|S|≤1} f̂(S)²")
    ax.set_ylabel("Norme spectrale de Walsh-Hadamard L₁(f̂)")
    ax.set_title("Espace de méta-complexité MCSP et Caméléons")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7.5)

    # Panel 3 : Evolution des invariants moyens avec s*
    ax = axes[2]
    ax.plot(ks, [stats_by_s[k]["walsh_l1_mean"] for k in ks], "o-", label="Norme spectrale Walsh L₁")
    ax.plot(ks, [stats_by_s[k]["influence_mean"] for k in ks], "s-", label="Influence totale Inf(f)")
    ax.plot(ks, [stats_by_s[k]["chow_mean"] * 3.0 for k in ks], "^--", label="3 × Énergie de Chow E₁")
    ax.set_xticks(ks)
    ax.set_xlabel("Taille minimale exacte s*")
    ax.set_ylabel("Valeur moyenne de l'invariant")
    ax.set_title("Invariants spectraux vs complexité minimale s*")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    fig.tight_layout()
    out_png = os.path.join(OUT, "atlas_and_mcsp.png")
    fig.savefig(out_png, dpi=130)
    print(f"\nTemps total Piste 3 (MCSP) : {time.time() - t0:.1f} s", flush=True)
    print(f"Graphique : {out_png}", flush=True)


if __name__ == "__main__":
    main()
