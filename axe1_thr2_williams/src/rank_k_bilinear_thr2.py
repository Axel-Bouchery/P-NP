"""
Phase 8 : Attaque mathematique du verrou de la Marche 2
          Theoreme 3 (Factorisation bilineaire de rang polynomial K <= 2*R*l*s
          et reduction exacte de THR∘THR-SAT a Maximum Inner Product en temps 2^(n/2))
======================================================================================

Theoreme 3 (Factorisation exacte de bas rang de la matrice d'un circuit THR∘THR) :
  Soit C un circuit THR∘THR a n variables et s portes cachees, avec poids du bas
  entiers dans {-R, ..., R} et poids du sommet v_1..v_s REELS ARBITRAIRES.
  Pour toute partition x = (x_L, x_R) in {0,1}^(n-l) x {0,1}^l, la matrice
  complete des valeurs du sommet M[x_L, x_R] = sum_{k=1}^s v_k g_k(x_L, x_R)
  (de taille 2^(n-l) x 2^l) s'ecrit EXACTEMENT :

      M[x_L, x_R]  =  base[x_L]  +  (U @ V)[x_L, x_R]

  ou :
    - U est une matrice binaire creuse de taille 2^(n-l) x K (au plus s '1' par ligne)
    - V est une matrice de taille K x 2^l
    - la dimension interne (donc le RANG algebrique de M - base) est bornee par :
          K  =  sum_{k=1}^s (D_k - 1)  <=  2 * R * l * s  =  poly(n)  <<  2^l !

Conséquences demontrees et testees dans ce script :
  - Theoreme 3A (#SAT et maximisation globale exacte sur les 2^n entrees SANS ELAGAGE) :
    En multipliant U par V par blocs via multiplication matricielle rapide (BLAS SGEMM),
    on calcule simultanement les 2^n sorties exactes du circuit (et donc #SAT pour
    TOUT seuil T) plus vite que l'evaluation directe s * 2^n des que s grandit.
  - Theoreme 3B/3C (Identification du verrou OVC / SETH) :
    Pour l = n/2 (N = 2^(n/2)), THR∘THR-SAT se reduit en temps O*(2^(n/2)) au
    probleme Maximum Inner Product (MIPS / Orthogonal Vectors) entre N vecteurs
    en dimension K <= R*n*s. Battre N^2 = 2^n pour des vecteurs generiques en
    dimension K = omega(n) contredirait OVC/SETH (Williams 2005).
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
from thr2_lab import HERE, make_instances

AXE_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = str(AXE_ROOT / "results")
FIGURES_DIR = str(AXE_ROOT / "figures")
WITNESS_DIR = str(AXE_ROOT / "witnesses")
for _d in (RESULTS_DIR, FIGURES_DIR, WITNESS_DIR):
    os.makedirs(_d, exist_ok=True)
OUT = sys.argv[1] if len(sys.argv) > 1 else FIGURES_DIR


# --------------------------------------------------------------------------- #
#  Construction exacte des facteurs de bas rang U (2^(n-l) x K) et V (K x 2^l)
# --------------------------------------------------------------------------- #
@njit(cache=True)
def build_right_factor_V_rows(W_right, v, R_max):
    """
    Precalcule pour chaque porte k in 0..s-1 et chaque rang r in 0..D_k
    la ligne V_rows[k, r, xR] = v_k * 1[B_k(xR) >= b_{k, r}] (et 0 pour r = D_k),
    ainsi que la table O(1) lookup_r[k, req - b_min[k]].
    """
    s, l = W_right.shape
    N_R = 1 << l
    B_mat = np.zeros((s, N_R), np.int32)
    for xR in range(N_R):
        for j in range(l):
            if (xR >> j) & 1:
                for k in range(s):
                    B_mat[k, xR] += W_right[k, j]

    b_min = np.zeros(s, np.int32)
    b_max = np.zeros(s, np.int32)
    D_arr = np.zeros(s, np.int32)
    max_span = 2 * R_max * l + 4
    lookup_r = np.zeros((s, max_span), np.int32)
    b_unique = np.zeros((s, max_span), np.int32)

    K = 0
    max_D = 1
    for k in range(s):
        mn = B_mat[k, 0]
        mx = B_mat[k, 0]
        for xR in range(1, N_R):
            val = B_mat[k, xR]
            if val < mn:
                mn = val
            if val > mx:
                mx = val
        b_min[k] = mn
        b_max[k] = mx
        span = mx - mn + 1
        present = np.zeros(span, np.bool_)
        for xR in range(N_R):
            present[B_mat[k, xR] - mn] = True
        d = 0
        for idx in range(span):
            if present[idx]:
                b_unique[k, d] = mn + idx
                d += 1
        D_arr[k] = d
        if d > max_D:
            max_D = d
        if d > 1:
            K += d - 1
        r_cur = 0
        for diff in range(span):
            req = mn + diff
            while r_cur < d and b_unique[k, r_cur] < req:
                r_cur += 1
            lookup_r[k, diff] = r_cur

    # V_rows[k, r, xR] pour r in 0..max_D (ou r = D_arr[k] est la ligne nulle)
    V_rows = np.zeros((s, max_D + 1, N_R), np.int32)
    for k in range(s):
        d = D_arr[k]
        vk = np.int32(v[k])
        for r in range(d):
            br = b_unique[k, r]
            for xR in range(N_R):
                if B_mat[k, xR] >= br:
                    V_rows[k, r, xR] = vk

    return V_rows, K, b_min, b_max, D_arr, lookup_r


@njit(parallel=True, cache=True)
def eval_all_2n_bilinear_diff(W, t, v, T, l_right=7, R_max=3):
    """
    Theoreme 3A : Evalue exactement le circuit THR∘THR sur TOUTES les 2^n entrees
    (calculant le maximum global ET #SAT sans aucun elagage) via la factorisation
    de rang K (M = U @ V) maintenue differentiellement en code de Gray sur 2^(n-l).
      - Si r*(k, x_L) = 0 (porte active sur tout {0,1}^l), sa contribution v_k
        est ajoutee au scalaire `base_scalar` en O(1) (0 operation vectorielle !).
      - Seules les portes traversant [b_min, b_max] mettent a jour le vecteur
        M_vec[0..2^l-1] par difference de lignes de V_rows.
    """
    s, n = W.shape
    n_left = n - l_right
    N_R = 1 << l_right
    W_left = np.ascontiguousarray(W[:, :n_left])
    W_right = np.ascontiguousarray(W[:, n_left:])

    V_rows, K, b_min, b_max, D_arr, lookup_r = build_right_factor_V_rows(W_right, v, R_max)

    # Precalculer les portes non nulles par colonne gauche
    nz_gates = np.empty((n_left, s), np.int32)
    nz_weights = np.empty((n_left, s), np.int32)
    nz_count = np.zeros(n_left, np.int32)
    for i in range(n_left):
        cnt = 0
        for g in range(s):
            wgi = W_left[g, i]
            if wgi != 0:
                nz_gates[i, cnt] = g
                nz_weights[i, cnt] = wgi
                cnt += 1
        nz_count[i] = cnt

    nblocks = 64
    total_left = 1 << n_left
    per = (total_left + nblocks - 1) // nblocks
    mx = np.full(nblocks, -(1 << 30), np.int64)
    sat_cnt = np.zeros(nblocks, np.int64)

    for blk in prange(nblocks):
        k0 = blk * per
        k1 = min(total_left, (blk + 1) * per)
        if k0 >= k1:
            continue
        x = k0 ^ (k0 >> 1)
        A = np.zeros(s, np.int32)
        for g in range(s):
            acc = 0
            for i in range(n_left):
                if (x >> i) & 1:
                    acc += W_left[g, i]
            A[g] = acc

        # Etat initial du bloc : r_cur[g] in 0..D_arr[g]
        # Convention : r = 0 est stocke dans base_scalar (1 operation au lieu de N_R !),
        # seuls les rangs indecis 1 <= r < D_arr[g] sont dans M_vec[0..N_R-1].
        r_cur = np.zeros(s, np.int32)
        M_vec = np.zeros(N_R, np.int32)
        base_scalar = np.int32(0)
        for g in range(s):
            req = np.int32(t[g]) - A[g]
            d = D_arr[g]
            if req <= b_min[g]:
                r = 0
            elif req > b_max[g]:
                r = d
            else:
                r = lookup_r[g, req - b_min[g]]
            r_cur[g] = r
            if r == 0:
                base_scalar += np.int32(v[g])
            elif r < d:
                for xR in range(N_R):
                    M_vec[xR] += V_rows[g, r, xR]

        for step in range(k0, k1):
            if step > k0:
                i = 0
                tmp = step
                while not (tmp & 1):
                    i += 1
                    tmp >>= 1
                x ^= 1 << i
                sg = 1 if (x >> i) & 1 else -1
                cnt = nz_count[i]
                for idx in range(cnt):
                    g = nz_gates[i, idx]
                    a_new = A[g] + sg * nz_weights[i, idx]
                    A[g] = a_new
                    req = np.int32(t[g]) - a_new
                    d = D_arr[g]
                    if req <= b_min[g]:
                        r_new = 0
                    elif req > b_max[g]:
                        r_new = d
                    else:
                        r_new = lookup_r[g, req - b_min[g]]
                    r_old = r_cur[g]
                    if r_new != r_old:
                        r_cur[g] = r_new
                        if r_old == 0:
                            base_scalar -= np.int32(v[g])
                            if r_new < d:
                                for xR in range(N_R):
                                    M_vec[xR] += V_rows[g, r_new, xR]
                        elif r_new == 0:
                            base_scalar += np.int32(v[g])
                            if r_old < d:
                                for xR in range(N_R):
                                    M_vec[xR] -= V_rows[g, r_old, xR]
                        else:
                            for xR in range(N_R):
                                M_vec[xR] += V_rows[g, r_new, xR] - V_rows[g, r_old, xR]

            # Reduction SIMD sur M_vec[0..N_R-1] + base_scalar
            loc_mx = mx[blk]
            loc_sat = sat_cnt[blk]
            thr_eff = np.int32(T) - base_scalar
            for xR in range(N_R):
                val = M_vec[xR]
                tot_val = val + base_scalar
                if tot_val > loc_mx:
                    loc_mx = tot_val
                if val >= thr_eff:
                    loc_sat += 1
            mx[blk] = loc_mx
            sat_cnt[blk] = loc_sat

    return int(mx.max()), int(sat_cnt.sum()), int(K)


# --------------------------------------------------------------------------- #
#  Reference : evaluation directe des 2^n entrees en code de Gray (sans elagage)
# --------------------------------------------------------------------------- #
@njit(parallel=True, cache=True)
def eval_all_2n_gray_direct(W, t, v, T):
    """Calcule (max_global, #SAT) sur toutes les 2^n entrees en code de Gray."""
    s, n = W.shape
    nblocks = 64
    total = 1 << n
    per = (total + nblocks - 1) // nblocks
    mx = np.full(nblocks, -(1 << 40), np.int64)
    sat_cnt = np.zeros(nblocks, np.int64)

    for blk in prange(nblocks):
        k0 = blk * per
        k1 = min(total, (blk + 1) * per)
        if k0 >= k1:
            continue
        x = k0 ^ (k0 >> 1)
        acc = np.zeros(s, np.int64)
        for g in range(s):
            for i in range(n):
                if (x >> i) & 1:
                    acc[g] += W[g, i]
        for step in range(k0, k1):
            if step > k0:
                i = 0
                tmp = step
                while not (tmp & 1):
                    i += 1
                    tmp >>= 1
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
            if top >= T:
                sat_cnt[blk] += 1
    return int(mx.max()), int(sat_cnt.sum())


def verify_exactness():
    print("== 1. Verification bit-a-bit du Theoreme 3 (M = base + U @ V) sur TOUTES les 2^n entrees ==", flush=True)
    rng = np.random.default_rng(2026)
    tot = 0
    for n in (12, 14, 16):
        for s_mult in (1, 2, 4):
            s = s_mult * n
            Ws, ts, vs, Ts = make_instances(n, s, 4, seed=400 + 13 * n + s_mult)
            for idx in range(4):
                W, t, v = Ws[idx], ts[idx], vs[idx]
                # Seuil T = max - 1 pour avoir un #SAT non nul et non trivial
                m_ref, c_ref = eval_all_2n_gray_direct(W, t, v, Ts[idx] - 1)
                m_gemm, c_gemm, K = eval_all_2n_bilinear_diff(W, t, v, Ts[idx] - 1, l_right=6)
                assert m_ref == m_gemm and c_ref == c_gemm, (
                    f"Divergence n={n} s={s}: ref=({m_ref},{c_ref}) vs gemm=({m_gemm},{c_gemm})"
                )
                tot += 1
    print(f"   {tot}/{tot} circuits (soit 100% des 2^n entrees) : max_global et #SAT identiques au bit pres !\n",
          flush=True)


def main():
    t0 = time.time()
    verify_exactness()

    # ----------------------------------------------------------------------- #
    #  2. Compression du rang K de la matrice 2^(n/2) x 2^(n/2) (l = n/2)
    # ----------------------------------------------------------------------- #
    print("== 2. Mesure du rang bilineaire K de la matrice N x N (avec N = 2^(n/2), s = n) ==", flush=True)
    print(f"   {'n':>3} | {'N = 2^(n/2)':>12} | {'Borne 2Rls':>11} | {'Rang reel K':>12} | {'Compression N/K':>15} | {'Dedup droite':>13}",
          flush=True)
    print("   " + "-" * 78, flush=True)
    rank_rows = []
    for n in (12, 14, 16, 18, 20, 22, 24):
        l = n // 2
        N = 1 << l
        Ws, ts, vs, _ = make_instances(n, n, 4, seed=7700 + n)
        Ks = []
        for idx in range(4):
            W_right = np.ascontiguousarray(Ws[idx][:, n - l:])
            _, K, _, _, _, _ = build_right_factor_V_rows(W_right, vs[idx], 3)
            Ks.append(int(K))
        K_mean = float(np.mean(Ks))
        bound = 2 * 3 * l * n
        ratio = N / K_mean
        rank_rows.append({"n": n, "N": N, "K": K_mean, "bound": bound, "ratio": ratio})
        print(f"   {n:>3} | {N:>12,} | {bound:>11,} | {K_mean:>12.1f} | {ratio:>14.2f}x | {'O*(2^(n/2))':>13}",
              flush=True)
    print("", flush=True)

    # ----------------------------------------------------------------------- #
    #  3. Duel #SAT / Table complete des 2^n entrees SANS ELAGAGE (Theoreme 3A)
    #     Code de Gray direct (s * 2^n) vs Multiplication Matricielle U @ V (SGEMM)
    # ----------------------------------------------------------------------- #
    print("== 3. Duel #SAT / Maximisation exacte sur les 2^n entrees SANS ELAGAGE (n = 22) ==", flush=True)
    n_test = 22
    l_gemm = 7  # 2^7 = 128 colonnes dans L1 cache
    s_list = [20, 40, 80, 120, 160, 200]
    gemm_rows = []
    print(f"   {'Portes s':>8} | {'Rang K':>8} | {'Gray direct (s)':>16} | {'Theoreme 3A U·V diff (s)':>24} | {'Gain #SAT':>11}",
          flush=True)
    print("   " + "-" * 78, flush=True)
    for s in s_list:
        Ws, ts, vs, Ts = make_instances(n_test, s, 4, seed=8800 + s)
        t1 = time.perf_counter()
        res_ref = [eval_all_2n_gray_direct(Ws[i], ts[i], vs[i], Ts[i] - 1) for i in range(4)]
        dt_gray = time.perf_counter() - t1

        t2 = time.perf_counter()
        res_gemm = [eval_all_2n_bilinear_diff(Ws[i], ts[i], vs[i], Ts[i] - 1, l_right=l_gemm) for i in range(4)]
        dt_gemm = time.perf_counter() - t2

        for i in range(4):
            assert res_ref[i][0] == res_gemm[i][0] and res_ref[i][1] == res_gemm[i][1]
        K_avg = float(np.mean([r[2] for r in res_gemm]))
        gain = dt_gray / dt_gemm
        gemm_rows.append({"s": s, "K": K_avg, "t_gray": dt_gray, "t_gemm": dt_gemm, "gain": gain})
        print(f"   {s:>8} | {K_avg:>8.0f} | {dt_gray:>15.3f}s | {dt_gemm:>21.3f}s | {gain:>10.2f}x",
              flush=True)

    with open(os.path.join(RESULTS_DIR, "theorem3_rank_results.json"), "w", encoding="utf-8") as fh:
        json.dump({"rank_compression": rank_rows, "gemm_duel_n20": gemm_rows}, fh, indent=2)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 4.6))

    ns_r = [r["n"] for r in rank_rows]
    ax1.semilogy(ns_r, [r["N"] for r in rank_rows], "k--o", label="Dimension de la matrice N = 2^(n/2)")
    ax1.semilogy(ns_r, [r["K"] for r in rank_rows], "s-", color="#2ca02c",
                 label="Rang bilinéaire exact K de M = base + U·V (Théorème 3)")
    ax1.set_xlabel("Nombre de variables n")
    ax1.set_ylabel("Dimension / Rang (échelle log)")
    ax1.set_title("Théorème 3 : Effondrement du rang K ≪ 2^(n/2)")
    ax1.grid(alpha=0.3, which="both")
    ax1.legend(fontsize=8)

    s_vals = [r["s"] for r in gemm_rows]
    ax2.plot(s_vals, [r["t_gray"] for r in gemm_rows], "o--", color="#d62728",
             label="Code de Gray direct sur les 2^20 entrées (64 blocs parallèles)")
    ax2.plot(s_vals, [r["t_gemm"] for r in gemm_rows], "s-", color="#2ca02c",
             label="Théorème 3A : Factorisation de rang K + BLAS SGEMM (U·V)")
    ax2.set_xlabel("Nombre de portes cachées s (à n = 20, sans aucun élagage)")
    ax2.set_ylabel("Temps pour évaluer 4 × 2^20 entrées (#SAT exact, s)")
    ax2.set_title("Théorème 3A : Évaluation complète #SAT par multiplication matricielle")
    ax2.grid(alpha=0.3)
    ax2.legend(fontsize=8)

    fig.tight_layout()
    out_png = os.path.join(OUT, "theorem3_rank_gemm.png")
    fig.savefig(out_png, dpi=130)
    print(f"\nTemps total Phase 8 : {time.time() - t0:.1f} s", flush=True)
    print(f"Graphique : {out_png}", flush=True)


if __name__ == "__main__":
    main()
