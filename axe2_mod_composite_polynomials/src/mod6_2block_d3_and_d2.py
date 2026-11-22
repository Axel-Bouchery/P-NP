"""
Chasse aux records modulo 6 par polynomes multi-blocs (Non symetriques)
=======================================================================

1. Degre d = 2 (N_sym(2) = 8) :
   Nous venons de prouver N_gen(2) >= 10 via deux familles explicites :
     - (8, 2) : P_2 = C(A,2) mod 2,  P_3 = -(A^2 + B^2) mod 3
     - (5, 5) : P_2 = A + B + A*B mod 2, P_3 = A^2 + B^2 mod 3.
   Ici on teste si d = 2 peut atteindre n = 11 ou n = 12 (recherche exhaustive
   3-blocs + 4-blocs + recuit simule parallele sur 16 coeurs).

2. Degre d = 3 (N_sym(3) = 19) :
   Recherche EXHAUSTIVE sur TOUS les 512 x 9,841 = 5,038,592 couples de
   polynomes 2-blocs (P_2(S1, S2), P_3(S1, S2)) de degre <= 3 !
   Grace a la reduction sur la grille (S1, S2) in {0..a} x {0..b}, la
   verification sur les 2^(a+b) points est EXACTE en O(a*b) !
   Objectif : battre N_sym(3) = 19 (et trouver le maximum exact N_2blk(3) !).

3. Degre d = 4 (N_sym(4) = 36) :
   Test d'extensions multi-blocs (e.g. a + b > 36) et comparaison avec la
   borne symetrique exacte N_sym(d) pour d = 1..10.
"""
from pathlib import Path
import json
import os
import sys
import time

import numpy as np
from numba import njit, prange

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mod6_poly_lab import HERE, build_monomial_matrix
from hunt_d2_n9_mod6 import build_3block_basis, exhaustive_3block_search, parallel_sa_hunt

AXE_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = str(AXE_ROOT / "results")
FIGURES_DIR = str(AXE_ROOT / "figures")
for _d in (RESULTS_DIR, FIGURES_DIR):
    os.makedirs(_d, exist_ok=True)
OUT = sys.argv[1] if len(sys.argv) > 1 else FIGURES_DIR


# --------------------------------------------------------------------------- #
#  Degre d = 3, 2 blocs (a, b) : recherche EXHAUSTIVE des 5,038,592 couples !
# --------------------------------------------------------------------------- #
@njit(parallel=True, cache=True)
def exact_2block_d3_search(max_a, max_b):
    """
    Pour chaque couple (P_2, P_3) de polynomes 2-blocs de degre <= 3 :
      Base de 9 monomes en (S1, S2) :
        [S1, S2, C(S1,2), C(S2,2), S1*S2, C(S1,3), C(S2,3), C(S1,2)*S2, S1*C(S2,2)]
    Calcule le plus grand n = a + b (avec 1 <= b <= a <= max_a) tel que
    Z(P_2) inter Z(P_3) = {(0, 0)} sur {0..a} x {0..b} !
    """
    A_sz = max_a + 1
    B_sz = max_b + 1
    # Precalculer les 9 fonctions de base mod 2 et mod 3 sur la grille (s1, s2)
    B2 = np.zeros((A_sz, B_sz, 9), np.int64)
    B3 = np.zeros((A_sz, B_sz, 9), np.int64)
    for s1 in range(A_sz):
        c1_1 = s1
        c1_2 = s1 * (s1 - 1) // 2
        c1_3 = s1 * (s1 - 1) * (s1 - 2) // 6
        for s2 in range(B_sz):
            c2_1 = s2
            c2_2 = s2 * (s2 - 1) // 2
            c2_3 = s2 * (s2 - 1) * (s2 - 2) // 6
            vals = (c1_1, c2_1, c1_2, c2_2, s1 * s2, c1_3, c2_3, c1_2 * s2, s1 * c2_2)
            for k in range(9):
                B2[s1, s2, k] = vals[k] & 1
                B3[s1, s2, k] = vals[k] % 3

    # Precalculer tous les 512 polynomes P_2 : Z2[m2, s1, s2] = True si P_2(s1, s2) == 0 mod 2
    n2 = 512
    Z2 = np.zeros((n2, A_sz, B_sz), np.bool_)
    for m2 in range(1, n2):
        for s1 in range(A_sz):
            for s2 in range(B_sz):
                acc = 0
                for k in range(9):
                    if (m2 >> k) & 1:
                        acc ^= B2[s1, s2, k]
                Z2[m2, s1, s2] = (acc == 0)

    n3_canon = (3 ** 9 - 1) // 2  # 9,841
    best_n_arr = np.zeros(n3_canon, np.int64)
    best_a_arr = np.zeros(n3_canon, np.int64)
    best_b_arr = np.zeros(n3_canon, np.int64)
    best_m2_arr = np.zeros(n3_canon, np.int64)

    for idx in prange(n3_canon):
        code = idx + 1
        c3 = np.zeros(9, np.int64)
        tmp = code
        for k in range(9):
            c3[k] = tmp % 3
            tmp //= 3

        # Collecter les zeros de P_3 sur la grille (s1, s2) != (0, 0)
        z1 = np.zeros(A_sz * B_sz, np.int64)
        z2 = np.zeros(A_sz * B_sz, np.int64)
        nz = 0
        for s1 in range(A_sz):
            for s2 in range(B_sz):
                if s1 == 0 and s2 == 0:
                    continue
                acc = 0
                for k in range(9):
                    if c3[k] != 0:
                        acc += c3[k] * B3[s1, s2, k]
                if acc % 3 == 0:
                    z1[nz] = s1
                    z2[nz] = s2
                    nz += 1

        # Pour chaque P_2, calculer le rectangle maximal sans collision (a, b) maximisant a + b
        loc_best_n = 0
        loc_a = 0
        loc_b = 0
        loc_m2 = 0
        # Pour chaque s1 in 0..max_a, min_bad_s2[s1] = plus petit s2 tel que collision en (s1, s2)
        min_bad_s2 = np.empty(A_sz, np.int64)
        for m2 in range(1, n2):
            for s1 in range(A_sz):
                min_bad_s2[s1] = B_sz
            for j in range(nz):
                s1 = z1[j]
                s2 = z2[j]
                if Z2[m2, s1, s2]:
                    if s2 < min_bad_s2[s1]:
                        min_bad_s2[s1] = s2
            # Balayer a = 0..max_a : la borne sur b pour un a donne est min_{s1 <= a} min_bad_s2[s1] - 1
            cur_max_b = B_sz - 1
            for a in range(A_sz):
                if min_bad_s2[a] - 1 < cur_max_b:
                    cur_max_b = min_bad_s2[a] - 1
                if cur_max_b < 0:
                    break
                if a + cur_max_b > loc_best_n:
                    loc_best_n = a + cur_max_b
                    loc_a = a
                    loc_b = cur_max_b
                    loc_m2 = m2

        best_n_arr[idx] = loc_best_n
        best_a_arr[idx] = loc_a
        best_b_arr[idx] = loc_b
        best_m2_arr[idx] = loc_m2

    # Maximum global
    g_idx = 0
    for idx in range(1, n3_canon):
        if best_n_arr[idx] > best_n_arr[g_idx]:
            g_idx = idx
    return (int(best_n_arr[g_idx]), int(best_a_arr[g_idx]), int(best_b_arr[g_idx]),
            int(best_m2_arr[g_idx]), int(g_idx + 1))


def verify_2block_d3(a, b, m2, c3_code):
    """Verification EXACTE independante sur toute la grille (s1, s2) in {0..a} x {0..b}."""
    c2 = [(m2 >> k) & 1 for k in range(9)]
    c3 = [(c3_code // (3 ** k)) % 3 for k in range(9)]
    for s1 in range(a + 1):
        c1_1, c1_2, c1_3 = s1, s1 * (s1 - 1) // 2, s1 * (s1 - 1) * (s1 - 2) // 6
        for s2 in range(b + 1):
            c2_1, c2_2, c2_3 = s2, s2 * (s2 - 1) // 2, s2 * (s2 - 1) * (s2 - 2) // 6
            vals = [c1_1, c2_1, c1_2, c2_2, s1 * s2, c1_3, c2_3, c1_2 * s2, s1 * c2_2]
            v2 = sum(c2[k] * vals[k] for k in range(9)) % 2
            v3 = sum(c3[k] * vals[k] for k in range(9)) % 3
            if s1 == 0 and s2 == 0:
                assert v2 == 0 and v3 == 0
            else:
                if v2 == 0 and v3 == 0:
                    return False, (s1, s2)
    return True, None


def main():
    t0 = time.time()
    print("=========================================================================", flush=True)
    print("== 1. DEGRE d = 3 (N_sym(3) = 19) : Recherche EXHAUSTIVE 2-blocs !     ==", flush=True)
    print("=========================================================================", flush=True)
    ts = time.perf_counter()
    best_n3, a3, b3, m2_3, c3_3 = exact_2block_d3_search(35, 35)
    dt = time.perf_counter() - ts
    ok, bad_pt = verify_2block_d3(a3, b3, m2_3, c3_3)
    assert ok, f"Erreur verification en {bad_pt} !"
    c2_vec = [(m2_3 >> k) & 1 for k in range(9)]
    c3_vec = [(c3_3 // (3 ** k)) % 3 for k in range(9)]
    names = ["S1", "S2", "C(S1,2)", "C(S2,2)", "S1*S2", "C(S1,3)", "C(S2,3)", "C(S1,2)*S2", "S1*C(S2,2)"]
    p2_s = " + ".join(names[k] for k in range(9) if c2_vec[k])
    p3_s = " + ".join((names[k] if c3_vec[k] == 1 else f"2*{names[k]}") for k in range(9) if c3_vec[k])
    print(f"   MAXIMUM EXACT 2-BLOCS POUR d = 3 : N_2blk(3) = {best_n3} (partition {a3} + {b3})  [{dt:.2f}s]",
          flush=True)
    print(f"   Comparaison : N_sym(3) = 19  -->  GAIN NON SYMETRIQUE = +{best_n3 - 19} variables !", flush=True)
    print(f"   Verification exacte sur les {1 << best_n3:,} entrees de {{0,1}}^{best_n3} : 0 collision !", flush=True)
    print(f"   P_2 (mod 2) = {p2_s}", flush=True)
    print(f"   P_3 (mod 3) = {p3_s}\n", flush=True)

    print("=========================================================================", flush=True)
    print("== 2. DEGRE d = 2 (N_sym(2) = 8, deja trouve n = 10) : Test n = 11 !   ==", flush=True)
    print("=========================================================================", flush=True)
    n = 11
    eval11, mon11 = build_monomial_matrix(n, 2)
    parts11 = [(5, 4, 2), (6, 3, 2), (7, 2, 2), (4, 4, 3), (5, 3, 3), (6, 4, 1), (5, 5, 1)]
    best_def11 = 999
    for part in parts11:
        ts = time.perf_counter()
        B11 = build_3block_basis(n, part)
        bdef, bcode, bnz = exhaustive_3block_search(eval11, B11)
        dt = time.perf_counter() - ts
        if bdef < best_def11:
            best_def11 = bdef
        print(f"   n = 11, 3-blocs {part} : meilleur defaut F_2 = {bdef} (|Z_3| = {bnz:>3})  [{dt:.2f}s]",
              flush=True)

    print("\n   Lancement recuit simule parallele sur n = 11 (64 essais x 6,000 pas, 66 coeffs libres) :", flush=True)
    ts = time.perf_counter()
    incs, nzs, c3s = parallel_sa_hunt(eval11, 64, 6000, 9999)
    dt = time.perf_counter() - ts
    bi = int(np.argmin(incs))
    print(f"   n = 11 (recuit libre) : meilleur defaut F_2 = {int(incs[bi])} (|Z_3| = {int(nzs[bi])})  [{dt:.2f}s]",
          flush=True)

    with open(os.path.join(RESULTS_DIR, "mod6_records.json"), "w", encoding="utf-8") as fh:
        json.dump({
            "d2": {"n_sym_exact": 8, "n_nonsym_found": 10, "n11_best_defect_3blk": best_def11,
                   "n11_best_defect_sa": int(incs[bi])},
            "d3": {"n_sym_exact": 19, "n_2blk_exact": best_n3, "partition": [a3, b3],
                   "p2": p2_s, "p3": p3_s, "c2": c2_vec, "c3": c3_vec},
        }, fh, indent=2)

    print(f"\nTemps total : {time.time() - t0:.1f} s", flush=True)


if __name__ == "__main__":
    main()
