"""
Piste A : Spectre des valeurs singulieres (SVD) de la matrice M = base + U·V
          et faisabilite de CAPP (approximation en rang faible) pour THR∘THR
============================================================================

Question centrale :
  Dans le Theoreme 3, la matrice M (de taille N x N avec N = 2^(n/2)) d'un
  circuit THR∘THR a un rang exact K <= R * n * s, qui grandit lineairement
  avec le nombre de portes s.
  MAIS pour CAPP (estimer Pr[M >= T] a +/- 0.05 pres, ce qui suffit pour la
  methode de Williams), on n'a pas besoin du rang exact : il suffit d'une
  approximation de rang r (M_r) en norme L2.

  Ce script mesure :
    1. Comment evolue le rang effectif K_95% et K_99% (capture de 95%/99% de
       l'energie de M centree) quand s passe de n a 2n, 4n, 8n, 16n (n = 16, N = 256).
       Le rang effectif sature-t-il a O(n) independamment de s ?
    2. Quelle erreur CAPP |Pr[M_r >= T] - Pr[M >= T]| obtient-on pour r = 4, 8, 16, 32 ?
    3. Quelles familles adversariales (PARITY, MOD3, Produit scalaire, Grands poids)
       resistent ou cedent a la compression SVD ?
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from thr2_lab import HERE, make_instances
from adversarial_suite import make_parity, make_mod3, make_large_weights
from industrial_and_orthogonal import make_inner_product_noisy

OUT = sys.argv[1] if len(sys.argv) > 1 else HERE


def build_exact_matrix_M(W, t, v):
    """Construit la matrice exacte M[x_L, x_R] = sum_k v_k 1[A_k(x_L) + B_k(x_R) >= t_k]."""
    s, n = W.shape
    l = n // 2
    n_left = n - l
    N_L = 1 << n_left
    N_R = 1 << l

    idx_L = np.arange(N_L, dtype=np.int64)[:, None]
    bits_L = ((idx_L >> np.arange(n_left, dtype=np.int64)[None, :]) & 1).astype(np.int64)
    A = bits_L @ W[:, :n_left].T  # (N_L, s)

    idx_R = np.arange(N_R, dtype=np.int64)[:, None]
    bits_R = ((idx_R >> np.arange(l, dtype=np.int64)[None, :]) & 1).astype(np.int64)
    B = bits_R @ W[:, n_left:].T  # (N_R, s)

    M = np.zeros((N_L, N_R), dtype=np.float64)
    for k in range(s):
        # (N_L, 1) + (1, N_R) >= t[k]
        M += v[k] * ((A[:, k:k+1] + B[:, k:k+1].T) >= t[k])
    return M


def analyze_matrix_svd_capp(M):
    """Calcule la SVD de M centree, les rangs effectifs et les erreurs CAPP."""
    mean_M = M.mean()
    Mc = M - mean_M
    U, S, Vt = np.linalg.svd(Mc, full_matrices=False)
    energy = S ** 2
    tot_energy = float(energy.sum())
    if tot_energy < 1e-12:
        return {"rank_exact": 0, "k90": 0, "k95": 0, "k99": 0, "capp_err": {}, "sv": S[:32].tolist()}

    cum = np.cumsum(energy) / tot_energy
    rank_exact = int((S > 1e-7 * S[0]).sum())
    k90 = int(np.searchsorted(cum, 0.90)) + 1
    k95 = int(np.searchsorted(cum, 0.95)) + 1
    k99 = int(np.searchsorted(cum, 0.99)) + 1

    # Seuil CAPP median (Pr[M >= T] ~ 0.5) et seuil queue (Pr[M >= T] ~ 0.1)
    T_med = float(np.median(M))
    T_q90 = float(np.quantile(M, 0.90))
    p_med = float((M >= T_med).mean())
    p_q90 = float((M >= T_q90).mean())

    capp_err = {}
    for r in (2, 4, 8, 16, 32):
        r_eff = min(r, len(S))
        Mr = mean_M + (U[:, :r_eff] * S[:r_eff]) @ Vt[:r_eff, :]
        err_med = abs(float((Mr >= T_med).mean()) - p_med)
        err_q90 = abs(float((Mr >= T_q90).mean()) - p_q90)
        capp_err[r] = {"err_med": err_med, "err_q90": err_q90}

    return {
        "rank_exact": rank_exact,
        "k90": k90,
        "k95": k95,
        "k99": k99,
        "capp_err": capp_err,
        "sv_norm": (S[:40] / S[0]).tolist(),
    }


def main():
    t0 = time.time()
    n = 16  # N = 2^8 = 256 -> matrices 256 x 256
    N = 1 << (n // 2)
    print(f"== 1. Evolution du rang effectif (SVD) quand le nombre de portes s explose (n={n}, N={N}) ==", flush=True)
    print(f"   {'Portes s':>8} | {'Rang exact':>10} | {'K_90%':>6} | {'K_95%':>6} | {'K_99%':>6} | "
          f"{'Err CAPP (r=8)':>14} | {'Err CAPP (r=16)':>15}", flush=True)
    print("   " + "-" * 82, flush=True)

    s_sweep = {}
    for mult in (1, 2, 4, 8, 16):
        s = mult * n
        Ws, ts, vs, _ = make_instances(n, s, 3, seed=5500 + s)
        res_list = [analyze_matrix_svd_capp(build_exact_matrix_M(Ws[i], ts[i], vs[i])) for i in range(3)]
        row = {
            "s": s,
            "rank_exact": float(np.mean([r["rank_exact"] for r in res_list])),
            "k90": float(np.mean([r["k90"] for r in res_list])),
            "k95": float(np.mean([r["k95"] for r in res_list])),
            "k99": float(np.mean([r["k99"] for r in res_list])),
            "err_r8": float(np.mean([r["capp_err"][8]["err_med"] for r in res_list])),
            "err_r16": float(np.mean([r["capp_err"][16]["err_med"] for r in res_list])),
            "err_q90_r16": float(np.mean([r["capp_err"][16]["err_q90"] for r in res_list])),
            "sv_norm": res_list[0]["sv_norm"],
        }
        s_sweep[s] = row
        print(f"   {s:>8} | {row['rank_exact']:>10.1f} | {row['k90']:>6.1f} | {row['k95']:>6.1f} | {row['k99']:>6.1f} | "
              f"{row['err_r8']:>14.4f} | {row['err_r16']:>15.4f}", flush=True)

    print(f"\n== 2. Comparaison sur familles adversariales (n={n}, N={N}) ==", flush=True)
    print(f"   {'Famille':<30} | {'Rang exact':>10} | {'K_90%':>6} | {'K_95%':>6} | {'K_99%':>6} | "
          f"{'Err CAPP (r=16)':>15}", flush=True)
    print("   " + "-" * 86, flush=True)

    W_par, t_par, v_par, _ = make_parity(n)
    W_m3, t_m3, v_m3, _ = make_mod3(n)
    W_ip, t_ip, v_ip, _ = make_inner_product_noisy(n, 1, seed=6600 + n)
    W_lw, t_lw, v_lw, _ = make_large_weights(n, 2 * n, 1, seed=6700 + n, R=1000)
    W_rnd, t_rnd, v_rnd, _ = make_instances(n, 2 * n, 1, seed=6800 + n)

    fams = {
        "Aleatoire (s=2n, R=3)": (W_rnd[0], t_rnd[0], v_rnd[0]),
        "Grands poids (s=2n, R=1000)": (W_lw[0], t_lw[0], v_lw[0]),
        "Produit scalaire + couplage": (W_ip[0], t_ip[0], v_ip[0]),
        "PARITY_n": (W_par, t_par, v_par),
        "MOD3_n": (W_m3, t_m3, v_m3),
    }
    fam_report = {}
    for fname, (W, t, v) in fams.items():
        r = analyze_matrix_svd_capp(build_exact_matrix_M(W, t, v))
        fam_report[fname] = {
            "rank_exact": r["rank_exact"],
            "k90": r["k90"],
            "k95": r["k95"],
            "k99": r["k99"],
            "err_r16": r["capp_err"][16]["err_med"],
            "sv_norm": r["sv_norm"],
        }
        print(f"   {fname:<30} | {r['rank_exact']:>10d} | {r['k90']:>6d} | {r['k95']:>6d} | {r['k99']:>6d} | "
              f"{r['capp_err'][16]['err_med']:>15.4f}", flush=True)

    # 3. Dependance en n : n = 10, 12, 14, 16, 18 (s = 4n)
    print("\n== 3. Passage a l'echelle du rang effectif K_95% en fonction de n (s = 4n) ==", flush=True)
    n_sweep = []
    for n_cur in (10, 12, 14, 16, 18):
        N_cur = 1 << (n_cur // 2)
        Ws, ts, vs, _ = make_instances(n_cur, 4 * n_cur, 2, seed=7100 + n_cur)
        r_list = [analyze_matrix_svd_capp(build_exact_matrix_M(Ws[i], ts[i], vs[i])) for i in range(2)]
        k95 = float(np.mean([r["k95"] for r in r_list]))
        k99 = float(np.mean([r["k99"] for r in r_list]))
        re = float(np.mean([r["rank_exact"] for r in r_list]))
        err16 = float(np.mean([r["capp_err"][16]["err_med"] for r in r_list]))
        n_sweep.append({"n": n_cur, "N": N_cur, "rank_exact": re, "k95": k95, "k99": k99, "err_r16": err16})
        print(f"   n = {n_cur:>2} (N = {N_cur:>3}, s = {4*n_cur:>2}) | Rang exact = {re:5.1f} | "
              f"K_95% = {k95:5.1f} | K_99% = {k99:5.1f} | Err CAPP (r=16) = {err16:.4f}", flush=True)

    with open(os.path.join(HERE, "piste_a_svd_results.json"), "w", encoding="utf-8") as fh:
        json.dump({"s_sweep_n16": s_sweep, "families_n16": fam_report, "n_sweep_s4n": n_sweep}, fh, indent=2)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(17.5, 4.7))

    # Panel 1 : Rang exact vs K_95% vs K_99% en fonction de s
    ax = axes[0]
    sv = sorted(s_sweep.keys())
    ax.plot(sv, [s_sweep[s]["rank_exact"] for s in sv], "k--o", label="Rang exact (sur N=256)")
    ax.plot(sv, [s_sweep[s]["k99"] for s in sv], "s-", color="#ff7f0e", label="Rang effectif K_99%")
    ax.plot(sv, [s_sweep[s]["k95"] for s in sv], "^-", color="#2ca02c", label="Rang effectif K_95%")
    ax.plot(sv, [s_sweep[s]["k90"] for s in sv], "d-", color="#1f77b4", label="Rang effectif K_90%")
    ax.set_xlabel("Nombre de portes cachées s (n = 16, N = 256)")
    ax.set_ylabel("Rang de la matrice M")
    ax.set_title("Saturation du rang effectif quand s explose")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    # Panel 2 : Decroissance des valeurs singulieres par famille
    ax = axes[1]
    for fname, r in fam_report.items():
        ax.semilogy(range(1, len(r["sv_norm"]) + 1), r["sv_norm"], "o-", markersize=3.5, label=fname)
    ax.set_xlabel("Indice de la valeur singulière j")
    ax.set_ylabel("σ_j / σ_1 (échelle log)")
    ax.set_title("Décroissance spectrale de M selon la famille (n = 16)")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=7.5)

    # Panel 3 : Evolution de K_95% vs N = 2^(n/2) en fonction de n (s = 4n)
    ax = axes[2]
    ns_s = [r["n"] for r in n_sweep]
    ax.semilogy(ns_s, [r["N"] for r in n_sweep], "k--o", label="Dimension N = 2^(n/2)")
    ax.semilogy(ns_s, [r["rank_exact"] for r in n_sweep], "s-", color="#d62728", label="Rang exact")
    ax.semilogy(ns_s, [r["k95"] for r in n_sweep], "^-", color="#2ca02c", label="Rang effectif K_95% (CAPP)")
    ax.set_xlabel("Nombre de variables n (avec s = 4n)")
    ax.set_ylabel("Rang / Dimension (log)")
    ax.set_title("Effondrement spectral en fonction de n (s = 4n)")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8)

    fig.tight_layout()
    out_png = os.path.join(OUT, "piste_a_svd_capp.png")
    fig.savefig(out_png, dpi=130)
    print(f"\nTemps total Piste A : {time.time() - t0:.1f} s", flush=True)
    print(f"Graphique : {out_png}", flush=True)


if __name__ == "__main__":
    main()
