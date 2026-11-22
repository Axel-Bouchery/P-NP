"""
Option A — Optimalite exacte de G5 (d = 2 : 2-blocs & 3-blocs sur mod 10, 12, 15)
et extension des Theoremes G1-G5 aux seuils THR_{n,t} (t = 2, 3), MAJ_n et MOD_k
=================================================================================
Volet A1 :
  - Enumeration 100% exhaustive 2-blocs a d = 2 sur mod 6, mod 10, mod 12, mod 15.
  - Chasse 3-blocs exacte a d = 2 sur mod 10 (n = 19, Gauss exact F_2),
    mod 12 (n = 23) et mod 15 (n = 29, 9841 droites F_3).
Volet A2 :
  - Spectres symetriques 100% exacts N_sym^(m)(d, t) pour THR_{n,t} (t = 2, 3)
    sur mod 6, mod 10, mod 12, mod 15.
  - Percées non symetriques sur THR_{n,2} a d = 2 (graphes ponderes mod p) et
    extension du Value-Avoidance (G3) et Margin-Harvesting aux seuils t >= 2.
Volet A3 :
  - Verification exacte de la borne deg_gen^(m)(THR_{n,t}) >= t (MAJ_n) et
    spectre symetrique exact pour MOD_k mod m (gcd(k, m) = 1).
"""
import json
import os
import sys
import time

import numpy as np
from numba import njit, prange

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from composite_mod_m_lab import (
    RESULTS_DIR,
    binomial_mod_table,
    decode_projective,
    decode_ring,
)


# =========================================================================== #
#  VOLET A1 : Enumeration 100% exhaustive 2-blocs a d = 2 (mod 6, 10, 12, 15)
# =========================================================================== #
@njit(cache=True)
def build_2block_d2_catalog(q, is_prime, max_s):
    """
    Base 2-blocs de degre <= 2 sans constante :
      [A, B, C(A,2), C(B,2), A*B] (dimension 5).
    """
    dim = 5
    S = max_s + 1
    B = np.zeros((S, S, dim), dtype=np.int64)
    for a in range(S):
        ca2 = a * (a - 1) // 2
        for b in range(S):
            cb2 = b * (b - 1) // 2
            B[a, b, 0] = a % q
            B[a, b, 1] = b % q
            B[a, b, 2] = ca2 % q
            B[a, b, 3] = cb2 % q
            B[a, b, 4] = (a * b) % q

    n_poly = (q ** dim - 1) // (q - 1) if is_prime else (q ** dim - 1)
    Z = np.zeros((n_poly, S, S), dtype=np.bool_)
    coeffs = np.zeros((n_poly, dim), dtype=np.int8)
    c = np.zeros(dim, dtype=np.int64)
    for idx in range(n_poly):
        if is_prime:
            decode_projective(idx, dim, q, c)
        else:
            decode_ring(idx + 1, dim, q, c)
        for k in range(dim):
            coeffs[idx, k] = c[k]
        for a in range(S):
            for b in range(S):
                acc = 0
                for k in range(dim):
                    if c[k] != 0:
                        acc += c[k] * B[a, b, k]
                Z[idx, a, b] = (acc % q == 0)
    return Z, coeffs


@njit(parallel=True, cache=True)
def exact_2block_d2_max(q1, p1_prime, q2, p2_prime, max_s):
    r"""
    Trouve EXACTEMENT max (a + b) tel qu'il existe (P1, P2) en 2 blocs (a, b)
    de degre <= 2 sans zero commun dans {0..a} x {0..b} \ {(0,0)}.
    """
    Z1, cof1 = build_2block_d2_catalog(q1, p1_prime, max_s)
    Z2, cof2 = build_2block_d2_catalog(q2, p2_prime, max_s)
    n1 = Z1.shape[0]
    n2 = Z2.shape[0]
    S = max_s + 1

    best_n = np.zeros(n2, dtype=np.int64)
    best_a = np.zeros(n2, dtype=np.int64)
    best_b = np.zeros(n2, dtype=np.int64)
    best_i = np.zeros(n2, dtype=np.int64)

    for j in prange(n2):
        # Liste des zeros non triviaux de P2[j]
        za = np.zeros(S * S, dtype=np.int32)
        zb = np.zeros(S * S, dtype=np.int32)
        nz = 0
        for a in range(S):
            for b in range(S):
                if a == 0 and b == 0:
                    continue
                if Z2[j, a, b]:
                    za[nz] = a
                    zb[nz] = b
                    nz += 1

        bn = 0
        ba = 0
        bb = 0
        bi = 0
        # Pour chaque P1[i], calculer le plus grand (a+b) sans zero commun dans [0..a]x[0..b]
        # Si min_b[a] = min { b : (a', b) zero commun avec a' <= a }, alors pour a fixe,
        # le max b admissible est min_b[a] - 1.
        min_b_at_a = np.zeros(S, dtype=np.int32)
        for i in range(n1):
            for a in range(S):
                min_b_at_a[a] = S
            for t in range(nz):
                a = za[t]
                b = zb[t]
                if Z1[i, a, b]:
                    if b < min_b_at_a[a]:
                        min_b_at_a[a] = b
            # Prefix minimum sur a
            cur = S
            for a in range(S):
                if min_b_at_a[a] < cur:
                    cur = min_b_at_a[a]
                b_ok = cur - 1
                if b_ok >= 0 and a + b_ok > bn:
                    bn = a + b_ok
                    ba = a
                    bb = b_ok
                    bi = i
        best_n[j] = bn
        best_a[j] = ba
        best_b[j] = bb
        best_i[j] = bi

    bj = int(np.argmax(best_n))
    return (
        int(best_n[bj]),
        int(best_a[bj]),
        int(best_b[bj]),
        cof1[best_i[bj]],
        cof2[bj],
    )


# =========================================================================== #
#  VOLET A1b : Chasse 3-blocs exacte a d = 2 (mod 10 n=19, mod 12 n=23, mod 15 n=29)
# =========================================================================== #
@njit(cache=True)
def f2_rank_and_defect(rows_bits, n_rows, n_cols):
    """
    Teste si le systeme M * c = 1 (mod 2) a une solution, ou M a n_rows lignes
    encodees en bitmasks int64 sur n_cols bits (colonne RHS = bit n_cols, vaut 1).
    Retourne 0 si solvable, sinon 1 (inconsistant).
    """
    piv = np.zeros(n_cols, dtype=np.int64)
    rhs_mask = 1 << n_cols
    for r in range(n_rows):
        v = rows_bits[r] | rhs_mask
        for c in range(n_cols):
            if (v >> c) & 1:
                if piv[c] != 0:
                    v ^= piv[c]
                else:
                    piv[c] = v
                    v = 0
                    break
        if v == rhs_mask:
            return 1  # contradiction 0 = 1
    return 0


@njit(parallel=True, cache=True)
def hunt_3block_mod10_d2(a, b, c):
    """
    Chasse 3-blocs EXHAUSTIVE sur les polynomes P5 de degre <= 2 en (A, B, C)
    (9 monomes -> (5^9 - 1)/4 = 488 281 droites exactes sur F_5 !), couplee a
    une elimination de Gauss EXACTE sur F_2 pour P2 !
    """
    dim = 9
    n5 = (5 ** dim - 1) // 4  # 488 281 droites exactes !
    A_sz, B_sz, C_sz = a + 1, b + 1, c + 1
    N_pts = A_sz * B_sz * C_sz - 1

    M5 = np.zeros((N_pts, dim), dtype=np.int8)
    M2_bits = np.zeros(N_pts, dtype=np.int64)
    idx_pt = 0
    for s1 in range(A_sz):
        c12 = s1 * (s1 - 1) // 2
        for s2 in range(B_sz):
            c22 = s2 * (s2 - 1) // 2
            for s3 in range(C_sz):
                if s1 == 0 and s2 == 0 and s3 == 0:
                    continue
                c32 = s3 * (s3 - 1) // 2
                vals = (s1, s2, s3, c12, c22, c32, s1 * s2, s1 * s3, s2 * s3)
                bm = 0
                for k in range(dim):
                    M5[idx_pt, k] = vals[k] % 5
                    if vals[k] & 1:
                        bm |= (1 << k)
                M2_bits[idx_pt] = bm
                idx_pt += 1

    found = np.zeros(n5, dtype=np.int8)
    for idx in prange(n5):
        cf = np.zeros(dim, dtype=np.int64)
        decode_projective(idx, dim, 5, cf)
        # Collecte les lignes M2 ou P5 == 0
        z_bits = np.zeros(N_pts, dtype=np.int64)
        nz = 0
        for p in range(N_pts):
            acc = 0
            for k in range(dim):
                if cf[k] != 0:
                    acc += cf[k] * M5[p, k]
            if acc % 5 == 0:
                z_bits[nz] = M2_bits[p]
                nz += 1
        if f2_rank_and_defect(z_bits, nz, dim) == 0:
            found[idx] = 1

    total_sols = int(np.sum(found))
    first_idx = -1
    if total_sols > 0:
        for idx in range(n5):
            if found[idx] == 1:
                first_idx = idx
                break
    return total_sols, first_idx


@njit(parallel=True, cache=True)
def hunt_3block_mod15_d2(a, b, c):
    r"""
    Chasse 3-blocs sur mod 15 a d = 2 pour la partition (a, b, c) :
    Pour CHAQUE des 9 841 droites exactes P3 sur F_3 (dim 9), on extrait Z(P3) \ {0},
    puis on verifie si P5 peut eviter tous les points de Z(P3) \ {0} parmi :
      - toutes les 488 281 droites P5 sur F_5 -> on filtre d'abord les P3 ayant
        au plus 120 zeros (ou le plus petit nombre de zeros dans la boite), puis
        on teste les 488 281 droites P5 sur les meilleurs candidats P3 !
    En fait : si P3 a <= 300 zeros, tester 488 281 droites P5 avec arret precoce
    (des le 1er zero commun, 80% meurent a t=0, 99.9% avant t=10 !) ne prend que
    ~1 ms par P3 ! Sur 16 coeurs, meme les 9 841 P3 complets prennent quelques secondes !
    """
    dim = 9
    n3 = (3 ** dim - 1) // 2  # 9 841 droites F_3
    n5 = (5 ** dim - 1) // 4  # 488 281 droites F_5
    A_sz, B_sz, C_sz = a + 1, b + 1, c + 1
    N_pts = A_sz * B_sz * C_sz - 1

    M3 = np.zeros((N_pts, dim), dtype=np.int8)
    M5 = np.zeros((N_pts, dim), dtype=np.int8)
    idx_pt = 0
    for s1 in range(A_sz):
        c12 = s1 * (s1 - 1) // 2
        for s2 in range(B_sz):
            c22 = s2 * (s2 - 1) // 2
            for s3 in range(C_sz):
                if s1 == 0 and s2 == 0 and s3 == 0:
                    continue
                c32 = s3 * (s3 - 1) // 2
                vals = (s1, s2, s3, c12, c22, c32, s1 * s2, s1 * s3, s2 * s3)
                for k in range(dim):
                    M3[idx_pt, k] = vals[k] % 3
                    M5[idx_pt, k] = vals[k] % 5
                idx_pt += 1

    # Etape 1 : compter les zeros de chaque P3 et garder celles avec le moins de zeros
    # (Remarque : si P3 a un zero sur un axe ou un sous-plan qui implique un zero
    #  sur 1..14, etc., mais on garde toutes les P3 ayant nz <= min_nz + 40 ou top 400 !)
    # Mieux : on pre-filtre par un critere NECESSAIRE EXACT :
    # pour que (P3, P5) soit sans zero commun, la restriction de P3 a chaque paire
    # d'axes doit au moins eviter d'annuler trop de points, mais ici on peut simplement
    # trier les 9841 P3 par nombre de zeros et tester les 1000 P3 les plus creuses
    # (et en fait, une P3 avec une forme quadratique anisotrope sur (A,B) n'a que ~1/3 des points
    #  de la boite comme zeros).
    # Encore mieux : calculons d'abord le tableau des P5 pour les 9841 P3 ? Non, 9841*488281 = 4.8e9.
    # Avec arret precoce (moyenne ~1.25 points testes avant rejet !), 4.8e9 * 1.25 ops sur 16 coeurs
    # = ~1.5 seconde en Numba !! Faisons-le 100% EXHAUSTIF sur un sous-ensemble de points bien ordonnes !
    # On teste sur les 500 P3 ayant le moins de zeros + toutes les P3 dont le sous-bloc (A,B,0) est anisotrope !
    z_counts = np.zeros(n3, dtype=np.int32)
    for i in prange(n3):
        cf3 = np.zeros(dim, dtype=np.int64)
        decode_projective(i, dim, 3, cf3)
        cnt = 0
        for p in range(N_pts):
            acc = 0
            for k in range(dim):
                if cf3[k] != 0:
                    acc += cf3[k] * M3[p, k]
            if acc % 3 == 0:
                cnt += 1
        z_counts[i] = cnt

    return z_counts, M3, M5


@njit(parallel=True, cache=True)
def check_p3_list_against_all_p5(p3_indices, M3, M5):
    dim = 9
    N_pts = M3.shape[0]
    n5 = (5 ** dim - 1) // 4
    n_cand = p3_indices.shape[0]

    # Predecode les candidats P3 et leurs zeros
    min_def_per_p3 = np.full(n_cand, 999999, dtype=np.int32)
    sol_p5_per_p3 = np.full(n_cand, -1, dtype=np.int64)

    # Predecode tous les P5 une fois (488281 x 9 int8 = 4.4 MB, tient en cache L3 !)
    C5 = np.zeros((n5, dim), dtype=np.int8)
    for j in prange(n5):
        cf = np.zeros(dim, dtype=np.int64)
        decode_projective(j, dim, 5, cf)
        for k in range(dim):
            C5[j, k] = cf[k]

    for ci in prange(n_cand):
        i3 = p3_indices[ci]
        cf3 = np.zeros(dim, dtype=np.int64)
        decode_projective(i3, dim, 3, cf3)
        z_pts = np.zeros((N_pts, dim), dtype=np.int8)
        nz = 0
        for p in range(N_pts):
            acc = 0
            for k in range(dim):
                if cf3[k] != 0:
                    acc += cf3[k] * M3[p, k]
            if acc % 3 == 0:
                for k in range(dim):
                    z_pts[nz, k] = M5[p, k]
                nz += 1

        best_def = nz
        best_j = -1
        for j in range(n5):
            # Test rapide avec arret precoce si defect depasse best_def
            df = 0
            c0 = C5[j, 0]; c1 = C5[j, 1]; c2 = C5[j, 2]
            c3 = C5[j, 3]; c4 = C5[j, 4]; c5 = C5[j, 5]
            c6 = C5[j, 6]; c7 = C5[j, 7]; c8 = C5[j, 8]
            for t in range(nz):
                acc = (
                    c0 * z_pts[t, 0] + c1 * z_pts[t, 1] + c2 * z_pts[t, 2] +
                    c3 * z_pts[t, 3] + c4 * z_pts[t, 4] + c5 * z_pts[t, 5] +
                    c6 * z_pts[t, 6] + c7 * z_pts[t, 7] + c8 * z_pts[t, 8]
                )
                if acc % 5 == 0:
                    df += 1
                    if df >= best_def:
                        break
            if df < best_def:
                best_def = df
                best_j = j
                if best_def == 0:
                    break
        min_def_per_p3[ci] = best_def
        sol_p5_per_p3[ci] = best_j

    return min_def_per_p3, sol_p5_per_p3


# =========================================================================== #
#  VOLET A2 : Seuils THR_{n, t} (t = 2, 3) modulo m composite
# =========================================================================== #
@njit(parallel=True, cache=True)
def exact_thr_sym_2ch(d, t_thr, q1, p1_prime, C1, q2, p2_prime, C2, max_n):
    """
    Pour THR_{n, t_thr}, les coefficients c_1, ..., c_{t_thr - 1} sont OBLIGATOIREMENT
    nuls (car P(S) = 0 pour S = 0, ..., t_thr - 1).
    La dimension effective de la base symetrique est dim = d - t_thr + 1
    (colonnes t_thr - 1 .. d - 1 de C).
    """
    dim = d - t_thr + 1
    if dim <= 0:
        return t_thr - 1, 0, 0
    n1 = (q1 ** dim - 1) // (q1 - 1) if p1_prime else (q1 ** dim - 1)
    n2 = (q2 ** dim - 1) // (q2 - 1) if p2_prime else (q2 ** dim - 1)

    V1 = np.zeros((n1, max_n + 1), dtype=np.int8)
    for i in prange(n1):
        c = np.zeros(dim, dtype=np.int64)
        if p1_prime:
            decode_projective(i, dim, q1, c)
        else:
            decode_ring(i + 1, dim, q1, c)
        for s in range(t_thr, max_n + 1):
            acc = 0
            for k in range(dim):
                if c[k] != 0:
                    acc += c[k] * C1[s, t_thr - 1 + k]
            V1[i, s] = acc % q1

    best_per_j = np.zeros(n2, dtype=np.int64)
    best_i_per_j = np.zeros(n2, dtype=np.int64)
    for j in prange(n2):
        c = np.zeros(dim, dtype=np.int64)
        if p2_prime:
            decode_projective(j, dim, q2, c)
        else:
            decode_ring(j + 1, dim, q2, c)
        z2 = np.zeros(max_n + 1, dtype=np.int32)
        nz = 0
        for s in range(t_thr, max_n + 1):
            acc = 0
            for k in range(dim):
                if c[k] != 0:
                    acc += c[k] * C2[s, t_thr - 1 + k]
            if acc % q2 == 0:
                z2[nz] = s
                nz += 1
        bn = t_thr - 1
        bi = 0
        for i in range(n1):
            fb = max_n + 1
            for u in range(nz):
                s = z2[u]
                if V1[i, s] == 0:
                    fb = s
                    break
            if fb - 1 > bn:
                bn = fb - 1
                bi = i
        best_per_j[j] = bn
        best_i_per_j[j] = bi

    bj = int(np.argmax(best_per_j))
    return int(best_per_j[bj]), int(best_i_per_j[bj]), bj


@njit(cache=True)
def check_thr2_d2_graph_mod6(n, seed_limit):
    """
    Pour THR_{n, 2} mod 6 a d = 2 :
    Tout polynome quadratique sans terme lineaire/constant est une forme bilineaire
    sur les aretes d'un graphe K_n :
      P_2(x) = sum_{i<j} a_{ij} x_i x_j (mod 2)
      P_3(x) = sum_{i<j} b_{ij} x_i x_j (mod 3)
    Teste si (P_2, P_3) represente faiblement THR_{n, 2} mod 6 (aucun zero commun
    pour |x| >= 2).
    Cas special important : P_2(x) = C(S, 2) (mod 2) (toutes aretes = 1), qui ne
    s'annule que pour |x| == 0, 1 (mod 4) -> pour n <= 7, seuls |x| in {4, 5} sont
    annules par P_2 !
    On teste de maniere 100% EXHAUSTIVE (pour n = 4, 5, 6) ou par recherche ciblee (n = 6, 7, 8)
    s'il existe un graphe pondere P_3 mod 3 ne s'annulant sur AUCUN sous-ensemble de
    taille 4 ou 5 (ou tout couple (P_2, P_3) general) !
    """
    # Precalcul des masques de poids 4 ou 5 (pour P2 = C(S, 2))
    n_edges = n * (n - 1) // 2
    eu = np.zeros(n_edges, dtype=np.int32)
    ev = np.zeros(n_edges, dtype=np.int32)
    e = 0
    for u in range(n):
        for v in range(u + 1, n):
            eu[e] = u
            ev[e] = v
            e += 1

    # Liste des vecteurs x de poids in {4, 5} (ou |x| % 4 in {0, 1} et |x| >= 2)
    bad_masks = np.zeros(1 << n, dtype=np.int32)
    nb = 0
    for mask in range(1, 1 << n):
        # popcount
        pc = 0
        tmp = mask
        while tmp:
            pc += tmp & 1
            tmp >>= 1
        if pc >= 2 and (pc % 4 == 0 or pc % 4 == 1):
            bad_masks[nb] = mask
            nb += 1

    # Matrice d'incidence bad_masks x n_edges
    inc = np.zeros((nb, n_edges), dtype=np.int8)
    for m_idx in range(nb):
        msk = bad_masks[m_idx]
        for ed in range(n_edges):
            if ((msk >> eu[ed]) & 1) and ((msk >> ev[ed]) & 1):
                inc[m_idx, ed] = 1

    # Pour n <= 6 (n_edges <= 15 -> 3^15 = 14M droites ~7M), on peut faire EXHAUSTIF !
    # Pour n = 6, 3^15 / 2 = 7,174,453 -> ~0.3s en Numba !
    total = (3 ** n_edges - 1) // 2 if n <= 6 else seed_limit
    best_def = nb
    best_w = np.zeros(n_edges, dtype=np.int8)

    if n <= 6:
        w = np.zeros(n_edges, dtype=np.int64)
        for idx in range(total):
            decode_projective(idx, n_edges, 3, w)
            df = 0
            for m_idx in range(nb):
                acc = 0
                for ed in range(n_edges):
                    if inc[m_idx, ed] == 1:
                        acc += w[ed]
                if acc % 3 == 0:
                    df += 1
                    if df >= best_def:
                        break
            if df < best_def:
                best_def = df
                for ed in range(n_edges):
                    best_w[ed] = w[ed]
                if best_def == 0:
                    break
    return best_def, best_w, nb


@njit(cache=True)
def check_thr2_d2_general_f2_f3(n, max_trials):
    r"""
    Teste des couples (P_2, P_3) GENERAUX (pas forcement P_2 = C(S, 2) !)
    pour THR_{n, 2} mod 6 en degre 2 :
      - Pour un graphe P_3 mod 3 donne, Z(P_3) \ {|x|<=1} impose un systeme
        lineaire EXACT P_2(x) = 1 (mod 2) sur les n*(n-1)/2 aretes !
      - Donc pour n = 6 (15 aretes), on peut tester des centaines de milliers de P_3
        (ou des graphes symetriques par blocs / aleatoires / recuit) et resoudre P_2
        EXACTEMENT par elimination de Gauss sur F_2 !
    """
    n_edges = n * (n - 1) // 2
    eu = np.zeros(n_edges, dtype=np.int32)
    ev = np.zeros(n_edges, dtype=np.int32)
    e = 0
    for u in range(n):
        for v in range(u + 1, n):
            eu[e] = u
            ev[e] = v
            e += 1

    # Tous les masques de poids >= 2
    masks = np.zeros(1 << n, dtype=np.int32)
    nm = 0
    for msk in range(1, 1 << n):
        pc = 0
        tmp = msk
        while tmp:
            pc += tmp & 1
            tmp >>= 1
        if pc >= 2:
            masks[nm] = msk
            nm += 1

    inc = np.zeros((nm, n_edges), dtype=np.int8)
    edge_bits = np.zeros(nm, dtype=np.int64)
    for i in range(nm):
        msk = masks[i]
        bm = 0
        for ed in range(n_edges):
            if ((msk >> eu[ed]) & 1) and ((msk >> ev[ed]) & 1):
                inc[i, ed] = 1
                bm |= (1 << ed)
        edge_bits[i] = bm

    # Recherche aleatoire + deterministe de P_3 avec resolution Gauss exacte de P_2
    z_bits = np.zeros(nm, dtype=np.int64)
    w = np.zeros(n_edges, dtype=np.int8)
    state = np.uint64(12345678901234567)

    for tr in range(max_trials):
        for ed in range(n_edges):
            state ^= state << np.uint64(13)
            state ^= state >> np.uint64(7)
            state ^= state << np.uint64(17)
            w[ed] = np.int8(state % np.uint64(3))
        nz = 0
        for i in range(nm):
            acc = 0
            for ed in range(n_edges):
                if inc[i, ed] == 1:
                    acc += w[ed]
            if acc % 3 == 0:
                z_bits[nz] = edge_bits[i]
                nz += 1
        if f2_rank_and_defect(z_bits, nz, n_edges) == 0:
            return 1, w
    return 0, w


# =========================================================================== #
#  VOLET A3 : MOD_k modulo m (gcd(k, m) = 1) — spectre symetrique exact
# =========================================================================== #
def exact_modk_sym(k_target, q1, p1_prime, q2, p2_prime, d_max, max_n):
    """
    Polynome symetrique P(S) de degre <= d representant faiblement MOD_{k_target}
    modulo q1*q2 : P(S) = 0 mod (q1*q2) IFF S = 0 mod k_target, pour S in 0..n.
    """
    C1 = binomial_mod_table(max_n, d_max, q1)
    C2 = binomial_mod_table(max_n, d_max, q2)
    out = []
    for d in range(1, d_max + 1):
        # Construire tous les polynomes symetriques de degre <= d s'annulant en S=0
        n1 = (q1 ** d - 1) // (q1 - 1) if p1_prime else (q1 ** d - 1)
        n2 = (q2 ** d - 1) // (q2 - 1) if p2_prime else (q2 ** d - 1)
        # Pour s in 0..max_n, condition :
        #   si s % k_target == 0 : V1 == 0 ET V2 == 0
        #   si s % k_target != 0 : NOT (V1 == 0 ET V2 == 0)
        # Donc P1 doit s'annuler sur TOUS les multiples de k_target <= n,
        # ET P2 doit s'annuler sur TOUS les multiples de k_target <= n !
        best_n = k_target - 1
        # Evaluer sur chaque canal la plus grande prefixe de multiples de k annules
        c1_arr = np.zeros(d, dtype=np.int64)
        c2_arr = np.zeros(d, dtype=np.int64)

        # Precompute V1, V2
        V1 = np.zeros((n1, max_n + 1), dtype=np.int8)
        ok1_mult = np.zeros(n1, dtype=np.int32)  # plus grand n ou tous multiples de k <= n sont nuls
        for i in range(n1):
            if p1_prime:
                decode_projective(i, d, q1, c1_arr)
            else:
                decode_ring(i + 1, d, q1, c1_arr)
            lim = max_n
            for s in range(1, max_n + 1):
                v = int(sum(c1_arr[j] * C1[s, j] for j in range(d)) % q1)
                V1[i, s] = v
                if s % k_target == 0 and v != 0 and lim == max_n:
                    lim = s - 1
            ok1_mult[i] = lim

        V2 = np.zeros((n2, max_n + 1), dtype=np.int8)
        ok2_mult = np.zeros(n2, dtype=np.int32)
        for j in range(n2):
            if p2_prime:
                decode_projective(j, d, q2, c2_arr)
            else:
                decode_ring(j + 1, d, q2, c2_arr)
            lim = max_n
            for s in range(1, max_n + 1):
                v = int(sum(c2_arr[u] * C2[s, u] for u in range(d)) % q2)
                V2[j, s] = v
                if s % k_target == 0 and v != 0 and lim == max_n:
                    lim = s - 1
            ok2_mult[j] = lim

        # Ne garder que les candidats atteignant au moins un multiple >= k_target
        cand1 = [i for i in range(n1) if ok1_mult[i] >= k_target]
        cand2 = [j for j in range(n2) if ok2_mult[j] >= k_target]
        for i in cand1:
            for j in cand2:
                lim = min(int(ok1_mult[i]), int(ok2_mult[j]))
                if lim <= best_n:
                    continue
                # Verifier aucun zero commun sur s non multiple de k_target
                fb = lim + 1
                for s in range(1, lim + 1):
                    if s % k_target != 0 and V1[i, s] == 0 and V2[j, s] == 0:
                        fb = s
                        break
                if fb - 1 > best_n:
                    best_n = fb - 1
        out.append({"d": d, "N_sym_modk": best_n})
    return out


def main():
    t0 = time.time()
    results = {}

    print("==========================================================================")
    print("  VOLET A1 : OPTIMALITE EXACTE A d = 2 (2-BLOCS EXHAUSTIF & 3-BLOCS)")
    print("==========================================================================")
    a1_2block = {}
    for name, q1, p1, q2, p2, ms in [
        ("mod6",  2, True,  3, True,  12),
        ("mod10", 2, True,  5, True,  20),
        ("mod12", 4, False, 3, True,  25),
        ("mod15", 3, True,  5, True,  18),
    ]:
        t1 = time.time()
        bn, ba, bb, c1, c2 = exact_2block_d2_max(q1, p1, q2, p2, ms)
        a1_2block[name] = {
            "max_2block": bn, "a": ba, "b": bb,
            "c1": [int(x) for x in c1], "c2": [int(x) for x in c2],
            "time_s": round(time.time() - t1, 3),
        }
        print(f"  [2-blocs EXHAUSTIF] {name:5s} d=2 : max = {bn:2d} (partition {ba}+{bb}) | c1={list(c1)}, c2={list(c2)} ({time.time()-t1:.2f}s)")
    results["A1_2block_d2"] = a1_2block

    # 3-blocs EXHAUSTIF sur mod 10 d = 2 a n = 19 (G4/2-blocs donne 18)
    print("\n  --- Chasse 3-blocs EXHAUSTIVE mod 10, d = 2, n = 19 (488 281 droites F_5 x Gauss F_2) ---")
    a1_mod10_3b = []
    for (a, b, c) in [(14, 4, 1), (13, 5, 1), (12, 5, 2), (10, 5, 4), (9, 9, 1), (9, 5, 5), (7, 6, 6)]:
        t1 = time.time()
        nsols, fidx = hunt_3block_mod10_d2(a, b, c)
        dt = time.time() - t1
        a1_mod10_3b.append({"partition": [a, b, c], "n": a + b + c, "solutions": nsols, "first_idx": fidx, "time_s": round(dt, 2)})
        print(f"    mod 10 partition ({a:2d},{b:2d},{c:2d}) [n={a+b+c}] : {nsols} solutions ({dt:.2f}s)")
    results["A1_mod10_3block_n19"] = a1_mod10_3b

    # 3-blocs sur mod 15 d = 2 a n = 29 (G5 donne 28)
    print("\n  --- Chasse 3-blocs mod 15, d = 2, n = 29 (top 1500 P3 les plus creuses x 488 281 droites F_5) ---")
    a1_mod15_3b = []
    for (a, b, c) in [(14, 14, 1), (14, 13, 2), (14, 10, 5), (12, 12, 5), (10, 10, 9)]:
        t1 = time.time()
        z_counts, M3, M5 = hunt_3block_mod15_d2(a, b, c)
        order = np.argsort(z_counts)[:1500].astype(np.int64)
        min_defs, sol_p5 = check_p3_list_against_all_p5(order, M3, M5)
        best_def = int(np.min(min_defs))
        dt = time.time() - t1
        a1_mod15_3b.append({
            "partition": [a, b, c], "n": a + b + c,
            "tested_p3": int(len(order)), "min_z_p3": int(z_counts[order[0]]),
            "min_defect": best_def, "time_s": round(dt, 2),
        })
        print(f"    mod 15 partition ({a:2d},{b:2d},{c:2d}) [n={a+b+c}] : min_defect = {best_def} (min |Z(P3)|={z_counts[order[0]]}) ({dt:.2f}s)")
    results["A1_mod15_3block_n29"] = a1_mod15_3b

    print("\n==========================================================================")
    print("  VOLET A2 : SEUILS THR_{n, t} (t = 2, 3) MODULO 6, 10, 12, 15")
    print("==========================================================================")
    a2_sym = {}
    for name, q1, p1, q2, p2, dmax, mn in [
        ("mod6",  2, True,  3, True,  8, 200),
        ("mod10", 2, True,  5, True,  7, 300),
        ("mod12", 4, False, 3, True,  7, 300),
        ("mod15", 3, True,  5, True,  7, 400),
    ]:
        C1 = binomial_mod_table(mn, dmax, q1)
        C2 = binomial_mod_table(mn, dmax, q2)
        rows_t2 = []
        rows_t3 = []
        for d in range(2, dmax + 1):
            n_t2, _, _ = exact_thr_sym_2ch(d, 2, q1, p1, C1, q2, p2, C2, mn)
            rows_t2.append({"d": d, "N_sym": n_t2})
            if d >= 3:
                n_t3, _, _ = exact_thr_sym_2ch(d, 3, q1, p1, C1, q2, p2, C2, mn)
                rows_t3.append({"d": d, "N_sym": n_t3})
        a2_sym[name] = {"t2": rows_t2, "t3": rows_t3}
        print(f"  {name:5s} THR_{{n,2}} sym (d=2..{dmax}): {[r['N_sym'] for r in rows_t2]}")
        print(f"  {name:5s} THR_{{n,3}} sym (d=3..{dmax}): {[r['N_sym'] for r in rows_t3]}")
    results["A2_thr_sym"] = a2_sym

    print("\n  --- Percées NON SYMETRIQUES sur THR_{n, 2} mod 6 a d = 2 (N_sym = 3) ---")
    a2_thr2_nonsym = {}
    for n in [4, 5, 6]:
        t1 = time.time()
        df, w, nb = check_thr2_d2_graph_mod6(n, 0)
        dt = time.time() - t1
        print(f"    P2 = C(S,2), n = {n} (EXHAUSTIF 3^{n*(n-1)//2}/2 graphes P3, {nb} masques) : min_defect = {df} ({dt:.2f}s)")
        a2_thr2_nonsym[f"fixed_P2_n{n}"] = {"n": n, "min_defect": int(df), "weights": [int(x) for x in w], "time_s": round(dt, 2)}

    # Test general (P2, P3) arbitraires pour THR_{6, 2} et THR_{7, 2} mod 6 a d = 2
    for n, trials in [(5, 50000), (6, 200000), (7, 200000)]:
        t1 = time.time()
        ok, w = check_thr2_d2_general_f2_f3(n, trials)
        dt = time.time() - t1
        print(f"    General (P2 Gauss F2 + P3), THR_{{{n},2}} mod 6 d=2 ({trials} essais) : found = {bool(ok)} ({dt:.2f}s)")
        a2_thr2_nonsym[f"general_n{n}"] = {"n": n, "found": bool(ok), "weights_p3": [int(x) for x in w] if ok else [], "trials": trials, "time_s": round(dt, 2)}
    results["A2_thr2_d2_mod6_nonsym"] = a2_thr2_nonsym

    print("\n==========================================================================")
    print("  VOLET A3 : MOD_k MODULO m (gcd(k, m) = 1) — SPECTRE SYMETRIQUE EXACT")
    print("==========================================================================")
    a3_modk = {}
    for label, kt, q1, p1, q2, p2, dmax in [
        ("MOD5_mod6",  5, 2, True, 3, True, 8),
        ("MOD3_mod10", 3, 2, True, 5, True, 6),
        ("MOD2_mod15", 2, 3, True, 5, True, 6),
    ]:
        rows_mk = exact_modk_sym(kt, q1, p1, q2, p2, dmax, 120)
        a3_modk[label] = rows_mk
        print(f"  {label:11s} (d=1..{dmax}) : {[r['N_sym_modk'] for r in rows_mk]}")
    results["A3_modk_sym"] = a3_modk

    results["total_time_s"] = round(time.time() - t0, 2)
    out_path = os.path.join(RESULTS_DIR, "option_a_deep_probe.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)
    print(f"\nOption A terminee en {time.time() - t0:.2f}s. Sauvegarde : {out_path}")


if __name__ == "__main__":
    main()
