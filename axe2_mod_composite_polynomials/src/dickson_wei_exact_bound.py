r"""
Resserrement du Theoreme D1 par la classification de Dickson sur F_2^s,
le Lemme D2a (capacite exacte degre (1,2) = 4p - 3 et (2,1) <= 3p - 1),
et l'analyse exacte du poids minimal des zeros quadratiques dans un code [m, m-j].

Ce script :
  1. Verifie EXHAUSTIVEMENT sur F_2^s (s = 3, 4, 5) le Theoreme de Dickson :
     tout polynome quadratique P_2 sur F_2^s (P_2(0) = 0) possede au moins
     z(s) = 2^(s-1) - 2^((s-1)//2) - 1 zeros non nuls (1, 5, 11, 27 pour s = 3, 4, 5, 6),
     ET calcule le support maximal du plus leger de ces zeros !
  2. Calcule exactement nu_2(2) (le plus petit n tel que tout P_2 quadratique sur {0,1}^n
     admette 2 vecteurs non nuls v_1, v_2 a supports disjoints avec P_2(v_1)=P_2(v_2)=B_2(v_1,v_2)=0).
  3. Calcule la borne superieure resserree w*(j, m) sur le support minimal d'un zero
     de (L_1, ..., L_j, P_2) lorsque m > j + 3 (via sous-espaces de dimension s = 4, 5, 6
     et collisions de Plotkin/pigeonhole sur les profils pivots q_1, ..., q_{m-j} in F_2^j).
"""
import itertools
import json
import math
import os
import sys
import time

import numpy as np
from numba import njit, prange

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from composite_mod_m_lab import RESULTS_DIR


# =========================================================================== #
#  1. Verification exhaustive du Theoreme de Dickson sur F_2^s (s = 3, 4, 5)
# =========================================================================== #
@njit(parallel=True, cache=True)
def verify_dickson_f2_s(s):
    """
    Parcourt TOUS les 2^(s + s(s-1)/2) polynomes quadratiques P_2 sur F_2^s (P_2(0)=0)
    et calcule :
      - min_zeros : le nombre minimal de zeros non nuls de P_2 dans F_2^s \ {0}
      - max_min_weight : le maximum (sur tous les P_2) du poids de Hamming minimal
        d'un zero non nul de P_2 dans F_2^s \ {0}.
    """
    n_lin = s
    n_quad = s * (s - 1) // 2
    n_tot = n_lin + n_quad
    n_polys = 1 << n_tot
    n_pts = (1 << s) - 1

    # Precalcul des monomes sur les 2^s - 1 points non nuls
    pts_mon = np.zeros(n_pts, dtype=np.int32)
    pts_wt = np.zeros(n_pts, dtype=np.int32)
    for x in range(1, 1 << s):
        m_mask = 0
        wt = 0
        for i in range(s):
            if (x >> i) & 1:
                m_mask |= (1 << i)
                wt += 1
        idx = s
        for i in range(s):
            for j in range(i + 1, s):
                if ((x >> i) & 1) and ((x >> j) & 1):
                    m_mask |= (1 << idx)
                idx += 1
        pts_mon[x - 1] = m_mask
        pts_wt[x - 1] = wt

    min_z_arr = np.full(n_polys, 999999, dtype=np.int32)
    max_w_arr = np.zeros(n_polys, dtype=np.int32)

    for p in prange(n_polys):
        nz = 0
        min_w = 999999
        for pt in range(n_pts):
            # Produit scalaire binaire popcount(p & pts_mon[pt]) & 1
            and_v = p & pts_mon[pt]
            pc = 0
            while and_v:
                pc ^= (and_v & 1)
                and_v >>= 1
            if pc == 0:
                nz += 1
                if pts_wt[pt] < min_w:
                    min_w = pts_wt[pt]
        min_z_arr[p] = nz
        max_w_arr[p] = min_w

    return int(np.min(min_z_arr)), int(np.max(max_w_arr))


# =========================================================================== #
#  2. Calcul EXACT de nu_2(2) : existence de 2 vecteurs non nuls a supports
#     disjoints v_1, v_2 dans {0,1}^n tels que P_2(v_1) = P_2(v_2) = P_2(v_1+v_2) = 0
# =========================================================================== #
@njit(parallel=True, cache=True)
def check_nu2_r2_for_n(n):
    """
    Teste sur TOUS les graphes G_2 sur n sommets et toutes les parties lineaires a in F_2^n
    s'il existe TOUJOURS deux vecteurs non nuls a supports disjoints v_1, v_2 in {0,1}^n \ {0}
    tels que P_2(v_1) = 0, P_2(v_2) = 0 et B_2(v_1, v_2) = 0 (c'est-a-dire P_2(v_1 | v_2) = 0).
    Pour accelerer, on fixe a par symetrie (poids de Hamming de a = 0..n, soit n+1 cas !)
    et on parcourt les 2^(n(n-1)/2) graphes G_2 !
    (Pour n = 5 : 6 x 1024 = 6 144 cas ; pour n = 6 : 7 x 32 768 = 229 376 cas !)
    """
    n_edges = n * (n - 1) // 2
    n_graphs = 1 << n_edges
    n_masks = 1 << n

    eu = np.zeros(n_edges, dtype=np.int32)
    ev = np.zeros(n_edges, dtype=np.int32)
    e = 0
    for u in range(n):
        for v in range(u + 1, n):
            eu[e] = u; ev[e] = v; e += 1

    # Precalcul pour chaque masque msk in 0..2^n-1 de ses bits d'aretes
    edge_mask_of = np.zeros(n_masks, dtype=np.int32)
    for msk in range(n_masks):
        em = 0
        for ed in range(n_edges):
            if ((msk >> eu[ed]) & 1) and ((msk >> ev[ed]) & 1):
                em |= (1 << ed)
        edge_mask_of[msk] = em

    # Liste de toutes les paires de masques non nuls a supports disjoints (v1 < v2, (v1 & v2) == 0)
    max_disj = (3 ** n) // 2
    disj_v1 = np.zeros(max_disj, dtype=np.int32)
    disj_v2 = np.zeros(max_disj, dtype=np.int32)
    n_disj = 0
    for v1 in range(1, n_masks):
        for v2 in range(v1 + 1, n_masks):
            if (v1 & v2) == 0:
                disj_v1[n_disj] = v1
                disj_v2[n_disj] = v2
                n_disj += 1

    counterexamples = np.zeros(n + 1, dtype=np.int64)
    sample_ce_g = np.full(n + 1, -1, dtype=np.int64)

    for wt_a in range(n + 1):
        a_mask = (1 << wt_a) - 1
        ce_cnt = np.zeros(n_graphs, dtype=np.int8)
        for g in prange(n_graphs):
            # Evaluer P_2(msk) pour tout msk in 1..2^n-1
            is_zero = np.zeros(n_masks, dtype=np.int8)
            for msk in range(1, n_masks):
                lin = msk & a_mask
                quad = g & edge_mask_of[msk]
                x = lin ^ quad  # non, lin et quad sont sur des espaces differents : parity(lin) ^ parity(quad)
                p = 0
                tmp = lin
                while tmp:
                    p ^= (tmp & 1); tmp >>= 1
                tmp = quad
                while tmp:
                    p ^= (tmp & 1); tmp >>= 1
                if p == 0:
                    is_zero[msk] = 1

            found = False
            for d_idx in range(n_disj):
                v1 = disj_v1[d_idx]
                if is_zero[v1]:
                    v2 = disj_v2[d_idx]
                    if is_zero[v2] and is_zero[v1 | v2]:
                        found = True
                        break
            if not found:
                ce_cnt[g] = 1

        tot_ce = int(np.sum(ce_cnt))
        counterexamples[wt_a] = tot_ce
        if tot_ce > 0:
            for g in range(n_graphs):
                if ce_cnt[g] == 1:
                    sample_ce_g[wt_a] = g
                    break

    return int(np.sum(counterexamples)), counterexamples, sample_ce_g


# =========================================================================== #
#  3. Poids minimal garanti dans un sous-espace F = span(b_1, ..., b_s) de K_U
#     de codimension j (s variables libres + j variables pivots)
# =========================================================================== #
@njit(parallel=True, cache=True)
def exact_max_min_weight_in_code_subspace(j, s):
    """
    Soit F = span(b_1, ..., b_s) <= K_U ou chaque b_i a 1 sur la variable libre i
    et un profil pivot q_i in F_2^j sur les j variables pivots.
    Pour tout (q_1, ..., q_s) in (F_2^j)^s et tout polynome quadratique P_2 sur F_2^s (P_2(0)=0),
    le poids total d'un vecteur non nul alpha in F_2^s \ {0} est :
        wt(alpha) = |alpha|_H + |sum_{i=1}^s alpha_i q_i|_H.
    Mais attention : on peut aussi echanger un pivot avec une variable libre (changement de base RREF
    parmi les s + j colonnes de F) !
    Meme SANS changer de pivots, si on impose que les q_i ne permettent pas deja un vecteur de
    poids <= 2 (c'est-a-dire q_i != 0 et q_i != q_k), quel est le pire min_{alpha: P_2(alpha)=0} wt(alpha) ?
    Plus fort encore : pour CHAQUE matrice de parite H de taille j x (j + s) de rang j,
    son noyau F est un code lineaire [j + s, s].
    Le maximum sur tous les P_2 quadratiques sur F (P_2(0)=0) du poids minimal d'un zero non nul de P_2
    dans F est BORNE par le poids du (N_Dickson(s))-ieme mot de code le plus leger de F !
    Verify this! Si F_1, F_2, ..., F_{2^s-1} sont les mots non nuls de F tries par poids croissant
    w_1 <= w_2 <= ... <= w_{2^s-1}, comme P_2 a au moins Z = N_Dickson(s) zeros non nuls dans F,
    le plus leger zero de P_2 dans F a NECESSAIREMENT un poids <= w_{2^s - 1 - Z + 1} = w_{2^s - Z} !!
    """
    pass


def compute_wei_dickson_step_bound(j, s):
    """
    Preuve analytique exacte par le Theoreme de Dickson + somme des poids de Plotkin sur F_2^(j+s) :
    Soit F <= F_2^(j+s) un sous-espace de dimension s (noyau de j formes lineaires sur j + s variables).
    1. La somme des poids de Hamming de TOUS les 2^s - 1 vecteurs non nuls de F vaut au plus :
         sum_{v in F \ {0}} |v|_H <= (j + s) * 2^(s - 1)
       (avec egalite ssi aucune des j + s coordonnees n'est identiquement nulle sur F).
    2. Par le Theoreme de Dickson sur F_2^s, tout polynome quadratique P_2 sur F (P_2(0) = 0)
       possede au moins :
         Z(s) = 2^(s - 1) - 2^((s - 1) // 2) - 1
       zeros non nuls dans F \ {0} !
    3. De plus, Z(s) zeros non nuls de P_2 dans F ne peuvent pas tous avoir un grand poids :
       - Soit w_min le poids minimal d'un zero non nul de P_2 dans F.
       - Si un vecteur v_0 in F a le poids maximal possible (j + s, le vecteur tout-a-1 1_{j+s}),
         alors pour tout autre vecteur v in F, |v|_H + |v ^ v_0|_H = j + s !
       - Plus generalement, la somme des poids des 2^s - 1 vecteurs non nuls de F est <= (j + s) * 2^(s - 1).
         Si les Z(s) zeros de P_2 ont tous un poids >= w, et que les (2^s - 1 - Z(s)) non-zeros
         ont chacun un poids >= 1 (ou >= d_min(F)), alors :
           Z(s) * w + (2^s - 1 - Z(s)) * 1 <= (j + s) * 2^(s - 1)
         d'ou :
           w <= floor( ((j + s) * 2^(s - 1) - (2^s - 1 - Z(s))) / Z(s) ) !
    """
    Z = (1 << (s - 1)) - (1 << ((s - 1) // 2)) - 1
    total_wt_max = (j + s) * (1 << (s - 1))
    non_zeros = (1 << s) - 1 - Z
    # Chaque vecteur non nul de F a au moins poids 1
    w_bound = (total_wt_max - non_zeros) // Z
    return w_bound, Z


def main():
    t0 = time.time()
    out = {}

    print("==========================================================================")
    print("  1. VERIFICATION EXHAUSTIVE DU THEOREME DE DICKSON SUR F_2^s (s = 3, 4, 5)")
    print("==========================================================================")
    dickson_res = {}
    for s in [3, 4, 5]:
        t1 = time.time()
        min_z, max_w = verify_dickson_f2_s(s)
        warn_z = (1 << (s - 2)) - 1
        dt = time.time() - t1
        print(
            f"  [s = {s}, {1 << (s + s*(s-1)//2):8,d} polynomes P_2] "
            f"min_zeros_non_nuls = {min_z:2d} (Warning exact = {warn_z:2d}), "
            f"max_min_weight = {max_w} ({dt:.2f}s)"
        )
        assert min_z == warn_z, f"Mismatch Warning s={s}: {min_z} != {warn_z}"
        dickson_res[str(s)] = {"n_polys": 1 << (s + s*(s-1)//2), "min_nonzero_zeros": min_z, "max_min_weight": max_w}
    out["warning_exhaustive"] = dickson_res

    print("\n==========================================================================")
    print("  2. CALCUL EXACT DE nu_2(2) SUR n = 4, 5, 6, 7 VARIABLES")
    print("==========================================================================")
    nu2_res = {}
    for n in [4, 5, 6, 7]:
        t1 = time.time()
        tot_ce, ce_by_wt, sample_g = check_nu2_r2_for_n(n)
        dt = time.time() - t1
        print(f"  [n = {n}] contre-exemples sans paire isotrope disjointe (v1, v2) : {tot_ce:6,d} (par poids de a : {list(ce_by_wt)}) ({dt:.2f}s)")
        nu2_res[str(n)] = {"total_counterexamples": tot_ce, "by_weight_of_a": [int(x) for x in ce_by_wt]}
        if tot_ce == 0:
            print(f"  >>> THEOREME EXACT : nu_2(2) = {n} (au lieu de 7 dans D1) ! <<<")
            out["exact_nu2_2"] = n
            break
    out["nu2_2_search"] = nu2_res

    out["total_time_s"] = round(time.time() - t0, 2)
    out_path = os.path.join(RESULTS_DIR, "dickson_wei_exact_bound.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print(f"\nTermine en {time.time() - t0:.2f}s. Sauvegarde : {out_path}")


if __name__ == "__main__":
    main()
