"""
Chasse intensive sur d = 2, n = 9 (et n = 10) modulo 6 :
Un polynome NON SYMETRIQUE de degre 2 peut-il battre la borne symetrique exacte N_sym(2) = 8 ?
==============================================================================================

1. Recherche EXHAUSTIVE sur tous les polynomes P_3 in F_3[x] a 2 blocs (a, b)
   et a 3 blocs (3, 3, 3) (9,841 polynomes P_3 canoniques), chacun couple a
   l'elimination de Gauss EXACTE sur F_2 (qui teste simultanement les 2^45
   polynomes P_2 possibles, symetriques ou non !).
2. Recherche massive par recuit simule parallele sur 16 coeurs (128 essais x
   15,000 pas) sur l'espace totalement libre des 45 coefficients de P_3 in F_3.
3. Inspection algebrique du meilleur candidat (si defaut = 0 -> PREUVE par
   construction explicite que N_gen(2) >= 9 > N_sym(2) = 8 ; si defaut = 1 ->
   analyse de la dependance lineaire F_2 qui bloque !).
"""
import json
import os
import sys
import time

import numpy as np
from numba import njit, prange

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mod6_poly_lab import (HERE, build_monomial_matrix, f2_defect_and_solve,
                           format_poly, verify_polynomial_pair)

OUT = sys.argv[1] if len(sys.argv) > 1 else HERE


# --------------------------------------------------------------------------- #
#  1. Recherche EXHAUSTIVE 3-blocs (b0, b1, b2) pour P_3 + Gauss exact sur F_2
# --------------------------------------------------------------------------- #
def build_3block_basis(n, sizes):
    """
    Construit les 9 polynomes de base de degre <= 2 symetriques sur 3 blocs
    de tailles (b0, b1, b2) (somme = n) evalues sur x in 1..2^n - 1 :
      S0, S1, S2, C(S0,2), C(S1,2), C(S2,2), S0*S1, S0*S2, S1*S2  (mod 3).
    """
    b0, b1, b2 = sizes
    N = (1 << n) - 1
    xs = np.arange(1, 1 << n, dtype=np.int64)
    S = np.zeros((3, N), dtype=np.int64)
    off = 0
    for grp, sz in enumerate((b0, b1, b2)):
        for bit in range(off, off + sz):
            S[grp] += (xs >> bit) & 1
        off += sz
    B = np.zeros((N, 9), dtype=np.int64)
    B[:, 0] = S[0] % 3
    B[:, 1] = S[1] % 3
    B[:, 2] = S[2] % 3
    B[:, 3] = (S[0] * (S[0] - 1) // 2) % 3
    B[:, 4] = (S[1] * (S[1] - 1) // 2) % 3
    B[:, 5] = (S[2] * (S[2] - 1) // 2) % 3
    B[:, 6] = (S[0] * S[1]) % 3
    B[:, 7] = (S[0] * S[2]) % 3
    B[:, 8] = (S[1] * S[2]) % 3
    return B


@njit(parallel=True, cache=True)
def exhaustive_3block_search(eval_mat, basis3):
    """
    Teste EXHAUSTIVEMENT les (3^9 - 1)/2 = 9,841 polynomes 3-blocs P_3 in F_3
    et resout pour chacun le systeme F_2 EXACT sur tous les 2^M polynomes P_2.
    """
    N, M = eval_mat.shape
    n_canon = (3 ** 9 - 1) // 2  # 9,841
    defects = np.full(n_canon, 999999, np.int64)
    z3_sizes = np.zeros(n_canon, np.int64)

    for idx in prange(n_canon):
        code = idx + 1
        c = np.zeros(9, np.int64)
        tmp = code
        for k in range(9):
            c[k] = tmp % 3
            tmp //= 3

        # Evaluer P_3(x) mod 3 sur x = 1..2^n - 1
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

        # Construire la sous-matrice Z_3 pour Gauss sur F_2
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
    best_def = defects[0]
    for idx in range(1, n_canon):
        if defects[idx] < best_def or (defects[idx] == best_def and z3_sizes[idx] < z3_sizes[best_idx]):
            best_def = defects[idx]
            best_idx = idx
    return int(best_def), int(best_idx + 1), int(z3_sizes[best_idx])


# --------------------------------------------------------------------------- #
#  2. Recuit simule massif parallele (128 essais sur 16 coeurs) avec 1-flip + 2-flip
# --------------------------------------------------------------------------- #
@njit(cache=True)
def sa_worker(eval_mat, n_steps, seed):
    np.random.seed(seed)
    N, M = eval_mat.shape
    c3 = np.random.randint(0, 3, size=M)

    val3 = np.zeros(N, np.int64)
    for col in range(M):
        cc = c3[col]
        if cc != 0:
            for r in range(N):
                if eval_mat[r, col]:
                    val3[r] += cc
    for r in range(N):
        val3[r] %= 3

    def eval_def(v3):
        nz = 0
        for r in range(N):
            if v3[r] == 0:
                nz += 1
        if nz == 0:
            return 0, 0
        sub = np.empty((nz, M), np.uint8)
        pos = 0
        for r in range(N):
            if v3[r] == 0:
                for col in range(M):
                    sub[pos, col] = eval_mat[r, col]
                pos += 1
        inc, _ = f2_defect_and_solve(sub)
        return inc, nz

    cur_inc, cur_nz = eval_def(val3)
    best_inc, best_nz = cur_inc, cur_nz
    best_c3 = c3.copy()
    if best_inc == 0:
        return 0, best_nz, best_c3

    temp = 1.2
    for step in range(n_steps):
        # 25% de 2-flips pour sortir des minima locaux a defaut = 1
        two_flip = (np.random.randint(0, 4) == 0)
        col1 = np.random.randint(0, M)
        d1 = 1 if np.random.randint(0, 2) == 0 else 2
        col2 = np.random.randint(0, M) if two_flip else -1
        if col2 == col1:
            col2 = -1
        d2 = 1 if np.random.randint(0, 2) == 0 else 2

        for r in range(N):
            if eval_mat[r, col1]:
                val3[r] = (val3[r] + d1) % 3
        if col2 >= 0:
            for r in range(N):
                if eval_mat[r, col2]:
                    val3[r] = (val3[r] + d2) % 3

        inc, nz = eval_def(val3)
        diff = (inc - cur_inc) + 0.001 * (nz - cur_nz)
        if diff <= 0.0 or np.random.random() < np.exp(-diff / max(temp, 0.03)):
            c3[col1] = (c3[col1] + d1) % 3
            if col2 >= 0:
                c3[col2] = (c3[col2] + d2) % 3
            cur_inc = inc
            cur_nz = nz
            if inc < best_inc or (inc == best_inc and nz < best_nz):
                best_inc = inc
                best_nz = nz
                best_c3 = c3.copy()
                if best_inc == 0:
                    return 0, best_nz, best_c3
        else:
            id1 = (3 - d1) % 3
            for r in range(N):
                if eval_mat[r, col1]:
                    val3[r] = (val3[r] + id1) % 3
            if col2 >= 0:
                id2 = (3 - d2) % 3
                for r in range(N):
                    if eval_mat[r, col2]:
                        val3[r] = (val3[r] + id2) % 3
        temp *= 0.9996
    return best_inc, best_nz, best_c3


@njit(parallel=True, cache=True)
def parallel_sa_hunt(eval_mat, n_trials, n_steps, base_seed):
    M = eval_mat.shape[1]
    incs = np.zeros(n_trials, np.int64)
    nzs = np.zeros(n_trials, np.int64)
    c3s = np.zeros((n_trials, M), np.int64)
    for t in prange(n_trials):
        inc, nz, c3 = sa_worker(eval_mat, n_steps, base_seed + 97 * t)
        incs[t] = inc
        nzs[t] = nz
        c3s[t] = c3
    return incs, nzs, c3s


def main():
    t0 = time.time()
    n = 9
    eval_mat, monoms = build_monomial_matrix(n, 2)
    M = len(monoms)
    print(f"== CHASSE INTENSIVE : d = 2, n = {n} (N_sym(2) = 8, M = {M} monomes, 2^9-1 = 511 points) ==\n",
          flush=True)

    print("1. Recherche EXHAUSTIVE 3-blocs sur P_3 (9,841 polynomes x Gauss F_2 exact) :", flush=True)
    partitions = [(3, 3, 3), (4, 3, 2), (5, 2, 2), (5, 3, 1), (6, 2, 1), (7, 1, 1), (4, 4, 1)]
    found_3blk = None
    for part in partitions:
        ts = time.perf_counter()
        B = build_3block_basis(n, part)
        bdef, bcode, bnz = exhaustive_3block_search(eval_mat, B)
        dt = time.perf_counter() - ts
        c = [(bcode // (3 ** k)) % 3 for k in range(9)]
        print(f"   Partition {part} : meilleur defaut F_2 = {bdef} (|Z_3| = {bnz:>3}, c = {c})  [{dt:.2f}s]",
              flush=True)
        if bdef == 0 and found_3blk is None:
            found_3blk = (part, c, B)

    if found_3blk is not None:
        part, c, B = found_3blk
        # Reconstruire la solution exacte P_2 et verifier sur les 512 points !
        val3 = (B @ np.array(c, dtype=np.int64)) % 3
        z3_idx = np.nonzero(val3 == 0)[0]
        inc, sol2 = f2_defect_and_solve(eval_mat[z3_idx])
        assert inc == 0
        # Convertir c (sur base 3-blocs) en coefficients sur monoms
        # Verification directe sur les 512 points :
        v2_all = (eval_mat @ sol2) % 2
        collisions = int(((v2_all == 0) & (val3 == 0)).sum())
        print(f"\n   *** DECOUVERTE ! SOLUTION EXACTE POUR d = 2, n = 9 TROUVEE (partition {part}) ! ***", flush=True)
        print(f"   Verification exhaustive sur les 511 points non nuls : {collisions} collisions !", flush=True)
        print(f"   P_3 (3-blocs {part}) : coeffs = {c}", flush=True)
        print(f"   P_2 (mod 2) = {format_poly(monoms, sol2, 2)}\n", flush=True)
        assert collisions == 0

        # Tester immediatement n = 10 avec 3-blocs !!
        print("== BONUS : Test immediat de d = 2, n = 10 par recherche exhaustive 3-blocs ! ==", flush=True)
        eval10, mon10 = build_monomial_matrix(10, 2)
        for part10 in [(4, 3, 3), (4, 4, 2), (5, 3, 2), (6, 2, 2), (5, 4, 1)]:
            ts = time.perf_counter()
            B10 = build_3block_basis(10, part10)
            bdef10, bcode10, bnz10 = exhaustive_3block_search(eval10, B10)
            dt = time.perf_counter() - ts
            c10 = [(bcode10 // (3 ** k)) % 3 for k in range(9)]
            print(f"   n = 10, partition {part10} : meilleur defaut F_2 = {bdef10} (|Z_3| = {bnz10:>3})  [{dt:.2f}s]",
                  flush=True)
    else:
        print("\n2. Recuit simule massif parallele sur 16 coeurs (128 essais x 12,000 pas, 45 coeffs libres) :",
              flush=True)
        ts = time.perf_counter()
        incs, nzs, c3s = parallel_sa_hunt(eval_mat, 128, 12000, 7777)
        dt = time.perf_counter() - ts
        best_i = int(np.argmin(incs))
        print(f"   Termine en {dt:.2f}s | meilleur defaut F_2 = {int(incs[best_i])} "
              f"(|Z_3| = {int(nzs[best_i])}, nb d'essais a defaut <= 1 : {int((incs <= 1).sum())}/128)",
              flush=True)
        if incs[best_i] == 0:
            c3 = c3s[best_i]
            val3 = (eval_mat @ c3) % 3
            z3_idx = np.nonzero(val3 == 0)[0]
            inc, sol2 = f2_defect_and_solve(eval_mat[z3_idx])
            assert inc == 0 and verify_polynomial_pair(n, 2, monoms, sol2, c3)
            print("\n   *** DECOUVERTE ! POLYNOME NON SYMETRIQUE EXACT TROUVE POUR d = 2, n = 9 ! ***", flush=True)
            print(f"   P_2 (mod 2) = {format_poly(monoms, sol2, 2)}", flush=True)
            print(f"   P_3 (mod 3) = {format_poly(monoms, c3, 3)}", flush=True)

    print(f"\nTemps total : {time.time() - t0:.1f} s", flush=True)


if __name__ == "__main__":
    main()
