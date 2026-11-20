"""
Exploration EXACTE et SANS BIAIS des polynomes modulo 6 (d = 1..10 symetriques,
d = 2 non symetrique a n = 11, et d = 3 non symetrique 2-blocs)
===============================================================================

Correction cruciale de l'enumeration canonique sur F_3^D \ {0} modulo {+1, -1} :
  Chaque droite projective de P^(D-1)(F_3) a un unique representant dont le
  dernier coefficient non nul vaut 1 (nombre exact = sum_{k=0}^{D-1} 3^k = (3^D - 1)/2).
  Verification : pour D = 3, les 13 droites incluent bien [0, 2, 1] -> N_sym(3) = 27 !
"""
import json
import os
import sys
import time

import numpy as np
from numba import njit, prange

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mod6_poly_lab import (HERE, binomial_mod_table, build_monomial_matrix,
                           f2_defect_and_solve, format_poly)
from hunt_d2_n9_mod6 import build_3block_basis

OUT = sys.argv[1] if len(sys.argv) > 1 else HERE


@njit(cache=True)
def decode_canonical_f3(idx, dim, out):
    """
    Bijection exacte entre idx in [0, (3^dim - 1)/2 - 1] et les (3^dim - 1)/2
    droites projectives de F_3^dim (dernier coefficient non nul egal a 1).
    """
    for j in range(dim):
        out[j] = 0
    p3 = 1
    k = 0
    rem = idx
    while rem >= p3:
        rem -= p3
        p3 *= 3
        k += 1
    for j in range(k):
        out[j] = rem % 3
        rem //= 3
    out[k] = 1


# --------------------------------------------------------------------------- #
#  1. Spectre symetrique 100% EXACT N_sym(d) pour d = 1..10
# --------------------------------------------------------------------------- #
@njit(parallel=True, cache=True)
def exact_n_sym_fixed(d, C2, C3):
    max_n = C2.shape[0] - 1
    n2 = 1 << d
    val2_zero = np.zeros((n2, max_n + 1), np.bool_)
    for mask in range(1, n2):
        for s in range(1, max_n + 1):
            acc = 0
            for k in range(d):
                if (mask >> k) & 1:
                    acc ^= C2[s, k]
            val2_zero[mask, s] = (acc == 0)

    n3_canon = (3 ** d - 1) // 2
    best_per_t = np.zeros(n3_canon, np.int64)
    best_m2_per_t = np.zeros(n3_canon, np.int64)

    for idx in prange(n3_canon):
        c3 = np.zeros(d, np.int64)
        decode_canonical_f3(idx, d, c3)
        z3 = np.zeros(max_n, np.int64)
        nz = 0
        for s in range(1, max_n + 1):
            acc = 0
            for k in range(d):
                if c3[k] != 0:
                    acc += c3[k] * C3[s, k]
            if acc % 3 == 0:
                z3[nz] = s
                nz += 1
        best_n = 0
        best_m2 = 0
        for mask in range(1, n2):
            first_bad = max_n + 1
            for j in range(nz):
                s = z3[j]
                if val2_zero[mask, s]:
                    first_bad = s
                    break
            n_ok = first_bad - 1
            if n_ok > best_n:
                best_n = n_ok
                best_m2 = mask
        best_per_t[idx] = best_n
        best_m2_per_t[idx] = best_m2

    best_idx = 0
    for idx in range(1, n3_canon):
        if best_per_t[idx] > best_per_t[best_idx]:
            best_idx = idx
    c3_best = np.zeros(d, np.int64)
    decode_canonical_f3(best_idx, d, c3_best)
    return int(best_per_t[best_idx]), int(best_m2_per_t[best_idx]), c3_best


# --------------------------------------------------------------------------- #
#  2. Recherche EXHAUSTIVE 2-blocs pour d = 3 (TOUTES les 9,841 droites F_3 x 511 F_2)
# --------------------------------------------------------------------------- #
@njit(parallel=True, cache=True)
def exact_2block_d3_fixed(max_a, max_b):
    A_sz = max_a + 1
    B_sz = max_b + 1
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

    n3_canon = (3 ** 9 - 1) // 2  # 9,841 droites exactes
    best_n_arr = np.zeros(n3_canon, np.int64)
    best_a_arr = np.zeros(n3_canon, np.int64)
    best_b_arr = np.zeros(n3_canon, np.int64)
    best_m2_arr = np.zeros(n3_canon, np.int64)

    for idx in prange(n3_canon):
        c3 = np.zeros(9, np.int64)
        decode_canonical_f3(idx, 9, c3)

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

        loc_best_n = 0
        loc_a = 0
        loc_b = 0
        loc_m2 = 0
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

    g_idx = 0
    for idx in range(1, n3_canon):
        if best_n_arr[idx] > best_n_arr[g_idx]:
            g_idx = idx
    c3_best = np.zeros(9, np.int64)
    decode_canonical_f3(g_idx, 9, c3_best)
    return (int(best_n_arr[g_idx]), int(best_a_arr[g_idx]), int(best_b_arr[g_idx]),
            int(best_m2_arr[g_idx]), c3_best)


# --------------------------------------------------------------------------- #
#  3. Recherche EXHAUSTIVE 3-blocs corrigee pour d = 2, n = 11
# --------------------------------------------------------------------------- #
@njit(parallel=True, cache=True)
def exhaustive_3block_fixed(eval_mat, basis3):
    N, M = eval_mat.shape
    n_canon = (3 ** 9 - 1) // 2
    defects = np.full(n_canon, 999999, np.int64)
    z3_sizes = np.zeros(n_canon, np.int64)
    for idx in prange(n_canon):
        c = np.zeros(9, np.int64)
        decode_canonical_f3(idx, 9, c)
        nz = 0
        for r in range(N):
            acc = 0
            for k in range(9):
                if c[k] != 0:
                    acc += c[k] * basis3[r, k]
            if acc % 3 == 0:
                nz += 1
        z3_sizes[idx] = nz
        if nz == 0:
            defects[idx] = 0
            continue
        sub = np.empty((nz, M), np.uint8)
        pos = 0
        for r in range(N):
            acc = 0
            for k in range(9):
                if c[k] != 0:
                    acc += c[k] * basis3[r, k]
            if acc % 3 == 0:
                for col in range(M):
                    sub[pos, col] = eval_mat[r, col]
                pos += 1
        inc, _ = f2_defect_and_solve(sub)
        defects[idx] = inc
    best_idx = 0
    for idx in range(1, n_canon):
        if defects[idx] < defects[best_idx] or (defects[idx] == defects[best_idx] and z3_sizes[idx] < z3_sizes[best_idx]):
            best_idx = idx
    c_best = np.zeros(9, np.int64)
    decode_canonical_f3(best_idx, 9, c_best)
    return int(defects[best_idx]), c_best, int(z3_sizes[best_idx])


def main():
    t0 = time.time()
    print("=========================================================================", flush=True)
    print("== 1. Spectre symetrique 100% EXACT N_sym(d) pour d = 1..10 (CORRIGE)  ==", flush=True)
    print("=========================================================================", flush=True)
    max_d = 10
    max_n = 450
    C2 = binomial_mod_table(max_n, max_d, 2)
    C3 = binomial_mod_table(max_n, max_d, 3)

    print(f"   {'deg d':>5} | {'Periode max p2*p3-1':>20} | {'EXACT N_sym(d)':>14} | {'N_sym / d^2':>11} | {'temps':>6}",
          flush=True)
    print("   " + "-" * 68, flush=True)
    sym_table = {}
    for d in range(1, max_d + 1):
        ts = time.perf_counter()
        n_opt, m2, c3_vec = exact_n_sym_fixed(d, C2[:, :d], C3[:, :d])
        dt = time.perf_counter() - ts
        p2 = 1
        while p2 <= d:
            p2 *= 2
        p3 = 1
        while p3 <= d:
            p3 *= 3
        c2_vec = [(m2 >> k) & 1 for k in range(d)]
        sym_table[d] = {
            "period_bound": p2 * p3 - 1, "p2": p2, "p3": p3,
            "n_sym_exact": n_opt, "ratio_d2": n_opt / (d * d),
            "c2": c2_vec, "c3": [int(x) for x in c3_vec],
        }
        print(f"   {d:>5} | {p2*p3-1:>10} ({p2:>2}x{p3:>2}-1) | {n_opt:>14} | {n_opt/(d*d):>11.2f} | {dt:>5.2f}s",
              flush=True)
        if d <= 7:
            print(f"           c2 = {c2_vec}   |   c3 = {[int(x) for x in c3_vec]}", flush=True)

    print("\n=========================================================================", flush=True)
    print("== 2. DEGRE d = 3 : Recherche EXHAUSTIVE 2-blocs (CORRIGEE, 5M paires) ==", flush=True)
    print("=========================================================================", flush=True)
    ts = time.perf_counter()
    best_n3, a3, b3, m2_3, c3_3 = exact_2block_d3_fixed(42, 42)
    dt = time.perf_counter() - ts
    c2_vec = [(m2_3 >> k) & 1 for k in range(9)]
    c3_vec = [int(x) for x in c3_3]
    names = ["S1", "S2", "C(S1,2)", "C(S2,2)", "S1*S2", "C(S1,3)", "C(S2,3)", "C(S1,2)*S2", "S1*C(S2,2)"]
    p2_s = " + ".join(names[k] for k in range(9) if c2_vec[k])
    p3_s = " + ".join((names[k] if c3_vec[k] == 1 else f"2*{names[k]}") for k in range(9) if c3_vec[k])
    # Verification independante sur toute la grille (a3 + 1) x (b3 + 1)
    for s1 in range(a3 + 1):
        c1_1, c1_2, c1_3 = s1, s1 * (s1 - 1) // 2, s1 * (s1 - 1) * (s1 - 2) // 6
        for s2 in range(b3 + 1):
            c2_1, c2_2, c2_3 = s2, s2 * (s2 - 1) // 2, s2 * (s2 - 1) * (s2 - 2) // 6
            vals = [c1_1, c2_1, c1_2, c2_2, s1 * s2, c1_3, c2_3, c1_2 * s2, s1 * c2_2]
            v2 = sum(c2_vec[k] * vals[k] for k in range(9)) % 2
            v3 = sum(c3_vec[k] * vals[k] for k in range(9)) % 3
            if s1 == 0 and s2 == 0:
                assert v2 == 0 and v3 == 0
            else:
                assert not (v2 == 0 and v3 == 0), f"Collision en ({s1}, {s2})!"
    print(f"   MAXIMUM EXACT 2-BLOCS POUR d = 3 : N_2blk(3) = {best_n3} (partition {a3} + {b3})  [{dt:.2f}s]",
          flush=True)
    print(f"   Comparaison avec N_sym(3) = {sym_table[3]['n_sym_exact']}  -->  GAIN = +{best_n3 - sym_table[3]['n_sym_exact']} variables !",
          flush=True)
    print(f"   Verification exacte sur 2^{best_n3} = {1<<best_n3:,} entrees : 0 collision !", flush=True)
    print(f"   P_2 (mod 2) = {p2_s}", flush=True)
    print(f"   P_3 (mod 3) = {p3_s}", flush=True)

    print("\n=========================================================================", flush=True)
    print("== 3. DEGRE d = 2, n = 11 : Recherche EXHAUSTIVE 3-blocs (CORRIGEE)    ==", flush=True)
    print("=========================================================================", flush=True)
    eval11, mon11 = build_monomial_matrix(11, 2)
    d2_n11 = {}
    for part in [(5, 4, 2), (6, 3, 2), (7, 2, 2), (8, 2, 1), (5, 5, 1), (4, 4, 3)]:
        ts = time.perf_counter()
        B11 = build_3block_basis(11, part)
        bdef, bc, bnz = exhaustive_3block_fixed(eval11, B11)
        dt = time.perf_counter() - ts
        d2_n11[str(part)] = {"defect": bdef, "z3": bnz, "c3": [int(x) for x in bc]}
        print(f"   n = 11, 3-blocs {part} : meilleur defaut F_2 = {bdef} (|Z_3| = {bnz:>3})  [{dt:.2f}s]", flush=True)
        if bdef == 0:
            val3 = (B11 @ bc) % 3
            z3_idx = np.nonzero(val3 == 0)[0]
            inc, sol2 = f2_defect_and_solve(eval11[z3_idx])
            assert inc == 0
            v2_all = (eval11 @ sol2) % 2
            assert int(((v2_all == 0) & (val3 == 0)).sum()) == 0
            print(f"   *** RECORD ! SOLUTION EXACTE POUR d = 2, n = 11 TROUVEE ({part}) ! ***", flush=True)
            print(f"   P_3 (3-blocs {part}) = {[int(x) for x in bc]}", flush=True)
            print(f"   P_2 (mod 2) = {format_poly(mon11, sol2, 2)}", flush=True)

    with open(os.path.join(HERE, "mod6_exact_complete.json"), "w", encoding="utf-8") as fh:
        json.dump({
            "symmetric_exact": sym_table,
            "d3_2block_exact": {
                "n_2blk": best_n3, "partition": [a3, b3], "p2": p2_s, "p3": p3_s,
                "c2": c2_vec, "c3": c3_vec,
            },
            "d2_n11_3block": d2_n11,
        }, fh, indent=2)

    print(f"\nTemps total : {time.time() - t0:.1f} s", flush=True)


if __name__ == "__main__":
    main()
