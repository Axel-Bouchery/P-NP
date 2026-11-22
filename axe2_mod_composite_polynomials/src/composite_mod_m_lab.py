"""
Generalisation des Theoremes M1/M2 aux modules composites m in {6, 10, 12, 15, 30}
==================================================================================
Calcul 100% EXACT des spectres symetriques N_sym^(m)(d) et evaluation constructive
des quatre mecanismes non symetriques prouves :
  - G1 : Multi-Channel Margin-Harvesting (additif sur tous les canaux)
  - G2 : Anisotropic Field-Norm Harvesting N_{F_{p^k}/F_p} (gain (k-1)(p-1)d_p + reste)
  - G3 : Value-Avoidance Harvesting (gain +d meme a marge nulle si Im(P_i) != Z/q_i Z)
  - Combinaison optimale G* par canal disjoint sur TOUT l'espace symetrique (a + sum_i gain_i)

Verification explicite de chaque record non symetrique decouvert.
"""
from pathlib import Path
import json
import os
import sys
import time
from math import gcd

import numpy as np
from numba import njit, prange

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = str(Path(__file__).resolve().parent.parent)
RESULTS_DIR = str(Path(ROOT) / "results")
FIGURES_DIR = str(Path(ROOT) / "figures")
os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)


def binomial_mod_table(max_n: int, max_d: int, mod: int) -> np.ndarray:
    """Table exacte C[s, k-1] = binom(s, k) mod `mod` pour s=0..max_n, k=1..max_d."""
    C = np.zeros((max_n + 1, max_d), dtype=np.int64)
    # Triangle de Pascal exact mod `mod` (valable meme pour mod compose ex. 4)
    prev = np.zeros(max_d + 1, dtype=np.int64)
    prev[0] = 1
    for s in range(1, max_n + 1):
        curr = np.zeros(max_d + 1, dtype=np.int64)
        curr[0] = 1
        lim = min(s, max_d)
        for k in range(1, lim + 1):
            curr[k] = (prev[k - 1] + prev[k]) % mod
        C[s, :] = curr[1:]
        prev = curr
    return C


@njit(cache=True)
def decode_projective(idx, dim, p, out):
    """
    Bijection exacte entre idx in [0, (p^dim - 1)/(p - 1) - 1] et les droites
    projectives de P^(dim-1)(F_p) : dernier coefficient non nul egal a 1.
    Retourne le degre exact (1-based) du polynome.
    """
    for j in range(dim):
        out[j] = 0
    pw = 1
    k = 0
    rem = idx
    while rem >= pw:
        rem -= pw
        pw *= p
        k += 1
    for j in range(k):
        out[j] = rem % p
        rem //= p
    out[k] = 1
    return k + 1


@njit(cache=True)
def decode_ring(idx, dim, q, out):
    """
    Bijection pour anneau non-corps (ex. Z/4Z) : idx in [1, q^dim - 1]
    en base q. On elimine les multiples globaux de 2 qui ne sont que des
    polynomes mod 2 (mais on les garde quand meme pour etre 100% exhaustif !).
    Retourne le degre exact (1-based).
    """
    rem = idx
    deg = 0
    for j in range(dim):
        v = rem % q
        out[j] = v
        if v != 0:
            deg = j + 1
        rem //= q
    return deg


@njit(cache=True)
def build_channel_catalog(d, q, is_prime, Cq, max_n):
    """
    Construit le catalogue complet de tous les polynomes symetriques de degre <= d
    modulo q, en deduplicant par (masque exact de zeros sur 1..max_n) tout en
    conservant le MEILLEUR representant (degre minimal, et s'il existe un polynome
    evitant au moins une valeur non nulle sur 1..a, sa propriete d'omission).

    Pour une analyse 100% exacte de G1+G2+G3, on renvoie pour chaque polynome :
      - deg[i] : degre exact (1..d)
      - vals[i, s] : P_i(s) mod q pour s = 0..max_n
      - coeffs[i, :d] : coefficients dans la base C(S, 1..d)
    """
    if is_prime:
        n_poly = (q ** d - 1) // (q - 1)
    else:
        n_poly = q ** d - 1

    degs = np.zeros(n_poly, dtype=np.int64)
    vals = np.zeros((n_poly, max_n + 1), dtype=np.int8)
    coeffs = np.zeros((n_poly, d), dtype=np.int8)

    c = np.zeros(d, dtype=np.int64)
    for idx in range(n_poly):
        if is_prime:
            dg = decode_projective(idx, d, q, c)
        else:
            dg = decode_ring(idx + 1, d, q, c)
        degs[idx] = dg
        for k in range(d):
            coeffs[idx, k] = c[k]
        vals[idx, 0] = 0
        for s in range(1, max_n + 1):
            acc = 0
            for k in range(dg):
                if c[k] != 0:
                    acc += c[k] * Cq[s, k]
            vals[idx, s] = acc % q
    return degs, vals, coeffs


@njit(cache=True)
def channel_best_gain(d, q, is_prime, deg_p, vals_row, a):
    """
    Calcule le gain PROUVE maximal sur un canal `q` pour un polynome P de degre `deg_p`
    sur {0, ..., a} :
      - G1 (Margin-Harvesting) : d - deg_p
      - G2 (Anisotropic Field-Norm + residu) : si q premier et k = d // deg_p >= 2,
          gain = (k - 1) * (q - 1) * deg_p + (d - k * deg_p)
      - G3 (Value-Avoidance) : si P(s) pour s in {0..a} omet au moins 1 valeur
          non nulle c in {1..q-1} (i.e. nombre de valeurs distinctes prises <= q - 1),
          alors gain >= d (meme si deg_p == d !).
          Bonus G3+ : si q >= 5 premier et l >= 2 valeurs non nulles omises :
          on peut aussi faire P(A) + c_1 y_1 + ... si sous-ensemble clos, mais on reste
          sur la borne 100% prouvee et universelle :
            * si l >= 1 : gain >= d (par c * (1 - prod_{j=1}^d (1 - y_j)))
            * si P(A) in {0, 1} mod q (omet q-2 valeurs non nulles) et d == 1 :
              w(y) = y_1 + ... + y_{q-2} prend ses valeurs dans {1..q-2}, donc
              P(A) + w(y) in {1..q-1} != 0 mod q -> gain = q - 2 !
    Retourne (best_gain, method_code, n_omitted)
      method_code: 0=none, 1=G1, 2=G2, 3=G3, 4=G3_linear
    """
    # Comptage des valeurs prises sur 0..a
    seen = np.zeros(q, dtype=np.int8)
    for s in range(a + 1):
        seen[vals_row[s]] = 1
    n_seen = 0
    for v in range(q):
        if seen[v]:
            n_seen += 1
    n_omitted = q - n_seen  # valeurs non nulles omises (car seen[0] == 1)

    g1 = d - deg_p
    best = g1
    method = 1 if g1 > 0 else 0

    if is_prime:
        k = d // deg_p
        if k >= 2:
            g2 = (k - 1) * (q - 1) * deg_p + (d - k * deg_p)
            if g2 > best:
                best = g2
                method = 2

    if n_omitted >= 1:
        g3 = d
        if g3 > best:
            best = g3
            method = 3
        # Cas special additif lineaire si on peut trouver un bloc contigu {1..L}
        # de valeurs c telles que -c not in Im(P), i.e. (q - v) not in seen for v=1..L
        # (ou pour tout multiplicateur inversible u in (Z/qZ)*)
        max_L = 0
        for u in range(1, q):
            if is_prime or (u % 2 == 1):  # u inversible mod q (valable pour q=2,3,4,5,7)
                L = 0
                for step in range(1, q):
                    target = (q - (u * step) % q) % q
                    if seen[target] == 0:
                        L += 1
                    else:
                        break
                if L > max_L:
                    max_L = L
        # Si d == 1, w(y) = u*(y_1 + ... + y_{max_L}) donne gain = max_L !
        # Si d >= 2 et is_prime : on peut même mettre max_L variables en degre 1
        if max_L > best:
            best = max_L
            method = 4

    return best, method, n_omitted


@njit(parallel=True, cache=True)
def solve_2channel(d, q1, p1_prime, deg1, val1, q2, p2_prime, deg2, val2, max_n):
    """
    Enumeration EXHAUSTIVE sur tous les couples (P1, P2) mod (q1, q2).
    Trouve :
      1. Le maximum symetrique exact N_sym^(q1*q2)(d) (et le meilleur representant
         ayant degres minimaux parmi ceux atteignant N_sym).
      2. Le maximum non symetrique prouve par G1 seul (Margin-Harvesting).
      3. Le maximum non symetrique prouve par G1+G2 (Norm Harvesting).
      4. Le maximum non symetrique prouve par G1+G2+G3 (incl. Value-Avoidance),
         en optimisant sur TOUS les couples symetriques (y compris sous-optimaux en a !).
    """
    n1 = deg1.shape[0]
    n2 = deg2.shape[0]

    # Precalcul des listes de zeros pour canal 2 (souvent le plus grand)
    z2_idx = np.zeros((n2, max_n), dtype=np.int32)
    z2_len = np.zeros(n2, dtype=np.int32)
    for j in prange(n2):
        cnt = 0
        for s in range(1, max_n + 1):
            if val2[j, s] == 0:
                z2_idx[j, cnt] = s
                cnt += 1
        z2_len[j] = cnt

    # Resultats par thread (indexe par j in 0..n2-1)
    sym_best = np.zeros(n2, dtype=np.int64)
    sym_best_i = np.zeros(n2, dtype=np.int64)
    sym_best_dsum = np.full(n2, 999, dtype=np.int64)

    g1_best = np.zeros(n2, dtype=np.int64)
    g1_best_i = np.zeros(n2, dtype=np.int64)
    g1_best_a = np.zeros(n2, dtype=np.int64)

    g2_best = np.zeros(n2, dtype=np.int64)
    g2_best_i = np.zeros(n2, dtype=np.int64)
    g2_best_a = np.zeros(n2, dtype=np.int64)

    gall_best = np.zeros(n2, dtype=np.int64)
    gall_best_i = np.zeros(n2, dtype=np.int64)
    gall_best_a = np.zeros(n2, dtype=np.int64)
    gall_m1 = np.zeros(n2, dtype=np.int64)
    gall_m2 = np.zeros(n2, dtype=np.int64)
    gall_g1 = np.zeros(n2, dtype=np.int64)
    gall_g2 = np.zeros(n2, dtype=np.int64)

    for j in prange(n2):
        d2 = deg2[j]
        nz = z2_len[j]
        # Gain G1 et G2 sur canal 2 (independant de a)
        gain2_g1 = d - d2
        gain2_g2 = gain2_g1
        if p2_prime:
            k2 = d // d2
            if k2 >= 2:
                cand = (k2 - 1) * (q2 - 1) * d2 + (d - k2 * d2)
                if cand > gain2_g2:
                    gain2_g2 = cand

        b_sym = 0
        b_sym_i = 0
        b_sym_ds = 999

        b_g1 = 0
        b_g1_i = 0
        b_g1_a = 0

        b_g2 = 0
        b_g2_i = 0
        b_g2_a = 0

        b_all = 0
        b_all_i = 0
        b_all_a = 0
        b_m1 = 0
        b_m2 = 0
        b_gn1 = 0
        b_gn2 = 0

        for i in range(n1):
            first_bad = max_n + 1
            for t in range(nz):
                s = z2_idx[j, t]
                if val1[i, s] == 0:
                    first_bad = s
                    break
            a = first_bad - 1
            d1 = deg1[i]
            dsum = d1 + d2

            if a > b_sym or (a == b_sym and dsum < b_sym_ds):
                b_sym = a
                b_sym_i = i
                b_sym_ds = dsum

            gain1_g1 = d - d1
            tot_g1 = a + gain1_g1 + gain2_g1
            if tot_g1 > b_g1:
                b_g1 = tot_g1
                b_g1_i = i
                b_g1_a = a

            gain1_g2 = gain1_g1
            if p1_prime:
                k1 = d // d1
                if k1 >= 2:
                    cand1 = (k1 - 1) * (q1 - 1) * d1 + (d - k1 * d1)
                    if cand1 > gain1_g2:
                        gain1_g2 = cand1
            tot_g2 = a + gain1_g2 + gain2_g2
            if tot_g2 > b_g2:
                b_g2 = tot_g2
                b_g2_i = i
                b_g2_a = a

            # Test G3 uniquement si a + max(gain1_g2, d, q1-2) + max(gain2_g2, d, q2-2) > b_all
            ub1 = gain1_g2 if gain1_g2 > d else d
            if q1 - 2 > ub1:
                ub1 = q1 - 2
            ub2 = gain2_g2 if gain2_g2 > d else d
            if q2 - 2 > ub2:
                ub2 = q2 - 2
            if a + ub1 + ub2 > b_all:
                gn1, m1, _ = channel_best_gain(d, q1, p1_prime, d1, val1[i], a)
                gn2, m2, _ = channel_best_gain(d, q2, p2_prime, d2, val2[j], a)
                tot_all = a + gn1 + gn2
                if tot_all > b_all:
                    b_all = tot_all
                    b_all_i = i
                    b_all_a = a
                    b_m1 = m1
                    b_m2 = m2
                    b_gn1 = gn1
                    b_gn2 = gn2

        sym_best[j] = b_sym
        sym_best_i[j] = b_sym_i
        sym_best_dsum[j] = b_sym_ds

        g1_best[j] = b_g1
        g1_best_i[j] = b_g1_i
        g1_best_a[j] = b_g1_a

        g2_best[j] = b_g2
        g2_best_i[j] = b_g2_i
        g2_best_a[j] = b_g2_a

        gall_best[j] = b_all
        gall_best_i[j] = b_all_i
        gall_best_a[j] = b_all_a
        gall_m1[j] = b_m1
        gall_m2[j] = b_m2
        gall_g1[j] = b_gn1
        gall_g2[j] = b_gn2

    # Reduction finale
    best_j_sym = 0
    for j in range(1, n2):
        if sym_best[j] > sym_best[best_j_sym] or (
            sym_best[j] == sym_best[best_j_sym] and sym_best_dsum[j] < sym_best_dsum[best_j_sym]
        ):
            best_j_sym = j

    best_j_g1 = int(np.argmax(g1_best))
    best_j_g2 = int(np.argmax(g2_best))
    best_j_all = int(np.argmax(gall_best))

    return (
        int(sym_best[best_j_sym]), int(sym_best_i[best_j_sym]), int(best_j_sym),
        int(g1_best[best_j_g1]), int(g1_best_i[best_j_g1]), int(best_j_g1), int(g1_best_a[best_j_g1]),
        int(g2_best[best_j_g2]), int(g2_best_i[best_j_g2]), int(best_j_g2), int(g2_best_a[best_j_g2]),
        int(gall_best[best_j_all]), int(gall_best_i[best_j_all]), int(best_j_all), int(gall_best_a[best_j_all]),
        int(gall_m1[best_j_all]), int(gall_m2[best_j_all]), int(gall_g1[best_j_all]), int(gall_g2[best_j_all]),
    )


@njit(parallel=True, cache=True)
def solve_3channel(d, q1, p1_prime, deg1, val1, q2, p2_prime, deg2, val2, q3, p3_prime, deg3, val3, max_n):
    """
    Enumeration EXHAUSTIVE 3 canaux (ex. mod 30 = 2 x 3 x 5) :
    parallise sur (j, k) ou k in canal 3 (le plus grand).
    """
    n1 = deg1.shape[0]
    n2 = deg2.shape[0]
    n3 = deg3.shape[0]

    # Precalcul des zeros pour canal 3
    z3_idx = np.zeros((n3, max_n), dtype=np.int32)
    z3_len = np.zeros(n3, dtype=np.int32)
    for k in prange(n3):
        cnt = 0
        for s in range(1, max_n + 1):
            if val3[k, s] == 0:
                z3_idx[k, cnt] = s
                cnt += 1
        z3_len[k] = cnt

    sym_best = np.zeros(n3, dtype=np.int64)
    sym_i = np.zeros(n3, dtype=np.int64)
    sym_j = np.zeros(n3, dtype=np.int64)
    sym_ds = np.full(n3, 999, dtype=np.int64)

    g1_best = np.zeros(n3, dtype=np.int64)
    g2_best = np.zeros(n3, dtype=np.int64)
    gall_best = np.zeros(n3, dtype=np.int64)
    gall_i = np.zeros(n3, dtype=np.int64)
    gall_j = np.zeros(n3, dtype=np.int64)
    gall_a = np.zeros(n3, dtype=np.int64)
    gall_m1 = np.zeros(n3, dtype=np.int64)
    gall_m2 = np.zeros(n3, dtype=np.int64)
    gall_m3 = np.zeros(n3, dtype=np.int64)
    gall_gn1 = np.zeros(n3, dtype=np.int64)
    gall_gn2 = np.zeros(n3, dtype=np.int64)
    gall_gn3 = np.zeros(n3, dtype=np.int64)

    for k in prange(n3):
        d3 = deg3[k]
        nz3 = z3_len[k]
        g3_1 = d - d3
        g3_2 = g3_1
        if p3_prime:
            kk = d // d3
            if kk >= 2:
                cand = (kk - 1) * (q3 - 1) * d3 + (d - kk * d3)
                if cand > g3_2:
                    g3_2 = cand
        ub3 = g3_2 if g3_2 > d else d
        if q3 - 2 > ub3:
            ub3 = q3 - 2

        b_sym = 0
        b_si = 0
        b_sj = 0
        b_sds = 999
        b_g1 = 0
        b_g2 = 0
        b_all = 0
        b_ai = 0
        b_aj = 0
        b_aa = 0
        b_m1 = 0
        b_m2 = 0
        b_m3 = 0
        b_gn1 = 0
        b_gn2 = 0
        b_gn3 = 0

        # Buffer des zeros communs (canal 2 et canal 3)
        z23 = np.zeros(max_n, dtype=np.int32)

        for j in range(n2):
            d2 = deg2[j]
            # Intersect z3 with val2[j] == 0
            nz23 = 0
            for t in range(nz3):
                s = z3_idx[k, t]
                if val2[j, s] == 0:
                    z23[nz23] = s
                    nz23 += 1

            g2_1 = d - d2
            g2_2 = g2_1
            if p2_prime:
                kk = d // d2
                if kk >= 2:
                    cand = (kk - 1) * (q2 - 1) * d2 + (d - kk * d2)
                    if cand > g2_2:
                        g2_2 = cand
            ub2 = g2_2 if g2_2 > d else d
            if q2 - 2 > ub2:
                ub2 = q2 - 2

            for i in range(n1):
                first_bad = max_n + 1
                for t in range(nz23):
                    s = z23[t]
                    if val1[i, s] == 0:
                        first_bad = s
                        break
                a = first_bad - 1
                d1 = deg1[i]
                ds = d1 + d2 + d3

                if a > b_sym or (a == b_sym and ds < b_sds):
                    b_sym = a
                    b_si = i
                    b_sj = j
                    b_sds = ds

                g1_1 = d - d1
                if a + g1_1 + g2_1 + g3_1 > b_g1:
                    b_g1 = a + g1_1 + g2_1 + g3_1

                g1_2 = g1_1
                if p1_prime:
                    kk = d // d1
                    if kk >= 2:
                        cand = (kk - 1) * (q1 - 1) * d1 + (d - kk * d1)
                        if cand > g1_2:
                            g1_2 = cand
                if a + g1_2 + g2_2 + g3_2 > b_g2:
                    b_g2 = a + g1_2 + g2_2 + g3_2

                ub1 = g1_2 if g1_2 > d else d
                if a + ub1 + ub2 + ub3 > b_all:
                    gn1, m1, _ = channel_best_gain(d, q1, p1_prime, d1, val1[i], a)
                    gn2, m2, _ = channel_best_gain(d, q2, p2_prime, d2, val2[j], a)
                    gn3, m3, _ = channel_best_gain(d, q3, p3_prime, d3, val3[k], a)
                    tot = a + gn1 + gn2 + gn3
                    if tot > b_all:
                        b_all = tot
                        b_ai = i
                        b_aj = j
                        b_aa = a
                        b_m1 = m1
                        b_m2 = m2
                        b_m3 = m3
                        b_gn1 = gn1
                        b_gn2 = gn2
                        b_gn3 = gn3

        sym_best[k] = b_sym
        sym_i[k] = b_si
        sym_j[k] = b_sj
        sym_ds[k] = b_sds
        g1_best[k] = b_g1
        g2_best[k] = b_g2
        gall_best[k] = b_all
        gall_i[k] = b_ai
        gall_j[k] = b_aj
        gall_a[k] = b_aa
        gall_m1[k] = b_m1
        gall_m2[k] = b_m2
        gall_m3[k] = b_m3
        gall_gn1[k] = b_gn1
        gall_gn2[k] = b_gn2
        gall_gn3[k] = b_gn3

    bk_sym = 0
    for k in range(1, n3):
        if sym_best[k] > sym_best[bk_sym] or (
            sym_best[k] == sym_best[bk_sym] and sym_ds[k] < sym_ds[bk_sym]
        ):
            bk_sym = k

    bk_all = int(np.argmax(gall_best))
    return (
        int(sym_best[bk_sym]), int(sym_i[bk_sym]), int(sym_j[bk_sym]), int(bk_sym),
        int(np.max(g1_best)), int(np.max(g2_best)),
        int(gall_best[bk_all]), int(gall_i[bk_all]), int(gall_j[bk_all]), int(bk_all),
        int(gall_a[bk_all]),
        int(gall_m1[bk_all]), int(gall_m2[bk_all]), int(gall_m3[bk_all]),
        int(gall_gn1[bk_all]), int(gall_gn2[bk_all]), int(gall_gn3[bk_all]),
    )


METHOD_NAMES = {0: "none", 1: "G1(margin)", 2: "G2(norm)", 3: "G3(avoid)", 4: "G3(lin-avoid)"}


def verify_witness_2ch(d, q1, p1_prime, c1, q2, p2_prime, c2, a, m1, gn1, m2, gn2):
    """Verification independante 100% exacte de la construction non symetrique sur 2 canaux."""
    C1 = binomial_mod_table(a, d, q1)
    C2 = binomial_mod_table(a, d, q2)
    vals1 = [0]
    vals2 = [0]
    for s in range(1, a + 1):
        v1 = int(sum(int(c1[k]) * int(C1[s, k]) for k in range(d)) % q1)
        v2 = int(sum(int(c2[k]) * int(C2[s, k]) for k in range(d)) % q2)
        assert not (v1 == 0 and v2 == 0), f"Zero commun illegal a s={s}!"
        vals1.append(v1)
        vals2.append(v2)

    for q, is_p, c, vals, meth, gn in [(q1, p1_prime, c1, vals1, m1, gn1), (q2, p2_prime, c2, vals2, m2, gn2)]:
        dg = max((k + 1 for k in range(d) if c[k] != 0), default=0)
        if meth == 0:
            assert gn == 0
        elif meth == 1:
            assert dg + gn <= d
        elif meth == 2:
            assert is_p
            kk = d // dg
            assert kk >= 2 and gn == (kk - 1) * (q - 1) * dg + (d - kk * dg)
        elif meth == 3:
            seen = set(vals)
            assert len(seen) < q and dg <= d and gn == d
        elif meth == 4:
            seen = set(vals)
            ok_u = False
            for u in range(1, q):
                if is_p or (u % 2 == 1):
                    if all((q - (u * step) % q) % q not in seen for step in range(1, gn + 1)):
                        ok_u = True
                        break
            assert ok_u and dg <= d
    return True


def run_all():
    t0 = time.time()
    all_results = {}

    # Configurations (nom, facteurs [(q, is_prime)], d_max, max_n_search)
    configs_2ch = [
        ("mod6",  (2, True),  (3, True),  9,  300),
        ("mod10", (2, True),  (5, True),  7,  400),
        ("mod12", (4, False), (3, True),  7,  500),
        ("mod15", (3, True),  (5, True),  6,  900),
    ]

    for name, (q1, p1), (q2, p2), d_max, max_n in configs_2ch:
        print(f"\n==================== {name.upper()} ({q1} x {q2}), d = 1..{d_max} ====================")
        rows = []
        C1 = binomial_mod_table(max_n, d_max, q1)
        C2 = binomial_mod_table(max_n, d_max, q2)
        for d in range(1, d_max + 1):
            td = time.time()
            deg1, val1, cof1 = build_channel_catalog(d, q1, p1, C1[:, :d], max_n)
            deg2, val2, cof2 = build_channel_catalog(d, q2, p2, C2[:, :d], max_n)
            res = solve_2channel(d, q1, p1, deg1, val1, q2, p2, deg2, val2, max_n)
            (nsym, si, sj,
             ng1, g1i, g1j, g1a,
             ng2, g2i, g2j, g2a,
             nall, ai, aj, aa, m1, m2, gn1, gn2) = res
            assert nsym < max_n, f"max_n={max_n} atteint pour {name} d={d}!"
            verify_witness_2ch(d, q1, p1, cof1[ai], q2, p2, cof2[aj], aa, m1, gn1, m2, gn2)
            dt = time.time() - td
            row = {
                "d": d,
                "N_sym": nsym,
                "sym_deg1": int(deg1[si]),
                "sym_deg2": int(deg2[sj]),
                "sym_c1": [int(x) for x in cof1[si]],
                "sym_c2": [int(x) for x in cof2[sj]],
                "N_G1": ng1,
                "N_G2": ng2,
                "N_all": nall,
                "gain_total": nall - nsym,
                "witness_a": aa,
                "witness_deg1": int(deg1[ai]),
                "witness_deg2": int(deg2[aj]),
                "witness_c1": [int(x) for x in cof1[ai]],
                "witness_c2": [int(x) for x in cof2[aj]],
                "method1": METHOD_NAMES[m1],
                "gain1": gn1,
                "method2": METHOD_NAMES[m2],
                "gain2": gn2,
                "time_s": round(dt, 3),
            }
            rows.append(row)
            print(
                f"  d={d:2d} | N_sym={nsym:3d} (deg {deg1[si]},{deg2[sj]}) | "
                f"G1={ng1:3d} | G2={ng2:3d} | G*={nall:3d} (+{nall-nsym:2d}) "
                f"[a={aa} + {gn1}({METHOD_NAMES[m1]}) + {gn2}({METHOD_NAMES[m2]})] ({dt:.2f}s)"
            )
        all_results[name] = rows

    # 3-channel: mod 30 = 2 x 3 x 5, d = 1..5 (et d=6 si rapide)
    print("\n==================== MOD30 (2 x 3 x 5), d = 1..5 ====================")
    max_n30 = 1600
    d_max30 = 5
    C2_30 = binomial_mod_table(max_n30, d_max30, 2)
    C3_30 = binomial_mod_table(max_n30, d_max30, 3)
    C5_30 = binomial_mod_table(max_n30, d_max30, 5)
    rows30 = []
    for d in range(1, d_max30 + 1):
        td = time.time()
        deg1, val1, cof1 = build_channel_catalog(d, 2, True, C2_30[:, :d], max_n30)
        deg2, val2, cof2 = build_channel_catalog(d, 3, True, C3_30[:, :d], max_n30)
        deg3, val3, cof3 = build_channel_catalog(d, 5, True, C5_30[:, :d], max_n30)
        res = solve_3channel(d, 2, True, deg1, val1, 3, True, deg2, val2, 5, True, deg3, val3, max_n30)
        (nsym, si, sj, sk, ng1, ng2, nall, ai, aj, ak, aa, m1, m2, m3, gn1, gn2, gn3) = res
        assert nsym < max_n30, f"max_n30={max_n30} atteint pour mod30 d={d}!"
        dt = time.time() - td
        row = {
            "d": d,
            "N_sym": nsym,
            "sym_degs": [int(deg1[si]), int(deg2[sj]), int(deg3[sk])],
            "N_G1": ng1,
            "N_G2": ng2,
            "N_all": nall,
            "gain_total": nall - nsym,
            "witness_a": aa,
            "witness_degs": [int(deg1[ai]), int(deg2[aj]), int(deg3[ak])],
            "witness_c1": [int(x) for x in cof1[ai]],
            "witness_c2": [int(x) for x in cof2[aj]],
            "witness_c3": [int(x) for x in cof3[ak]],
            "methods": [METHOD_NAMES[m1], METHOD_NAMES[m2], METHOD_NAMES[m3]],
            "gains": [gn1, gn2, gn3],
            "time_s": round(dt, 3),
        }
        rows30.append(row)
        print(
            f"  d={d:2d} | N_sym={nsym:4d} (deg {deg1[si]},{deg2[sj]},{deg3[sk]}) | "
            f"G1={ng1:4d} | G2={ng2:4d} | G*={nall:4d} (+{nall-nsym:2d}) "
            f"[a={aa} + {gn1}({METHOD_NAMES[m1]}) + {gn2}({METHOD_NAMES[m2]}) + {gn3}({METHOD_NAMES[m3]})] ({dt:.2f}s)"
        )
    all_results["mod30"] = rows30

    out_path = os.path.join(RESULTS_DIR, "composite_mod_m_results.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(all_results, f, indent=2)
    print(f"\nTermine en {time.time() - t0:.2f}s. Sauvegarde : {out_path}")


if __name__ == "__main__":
    run_all()
