"""
Chasse aux constructions NON SYMETRIQUES depassant N_sym(d) (et sondage > d^2)
==============================================================================

Exploite la reduction sur grille de poids (S_1, ..., S_k) + elimination de Gauss
EXACTE sur F_2 (qui teste simultanement les 2^M polynomes P_2 pour chaque P_3) :

  1. Degre d = 3, 3 blocs (S_1, S_2, S_3) : objectif battre N_sym(3) = 27
     (M = 19 monomes de degre <= 3 en 3 blocs).
  2. Degre d = 4, 2 blocs (S_1, S_2) : objectif battre N_sym(4) = 36
     (M = 14 monomes de degre <= 4 en 2 blocs).
  3. Degre d = 5, 2 blocs (S_1, S_2) : objectif battre N_sym(5) = 53
     (M = 20 monomes de degre <= 5 en 2 blocs).
  4. Degre d = 6, 2 blocs (S_1, S_2) : objectif battre N_nonsym(6) = 63
     (M = 27 monomes de degre <= 6 en 2 blocs).

Chaque candidat trouve avec defaut F_2 = 0 est verifie 100% EXACTEMENT sur toute
la grille des poids (ce qui prouve 0 collision sur les 2^n entrees booleennes !).
"""
from pathlib import Path
import json
import os
import sys
import time

import numpy as np
from numba import njit, prange

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mod6_poly_lab import HERE, f2_defect_and_solve

AXE_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = str(AXE_ROOT / "results")
FIGURES_DIR = str(AXE_ROOT / "figures")
for _d in (RESULTS_DIR, FIGURES_DIR):
    os.makedirs(_d, exist_ok=True)
OUT = sys.argv[1] if len(sys.argv) > 1 else FIGURES_DIR


def binom_int(n, k):
    if k < 0 or k > n:
        return 0
    num = 1
    den = 1
    for i in range(1, k + 1):
        num *= (n - i + 1)
        den *= i
    return num // den


def build_2block_grid_matrices(a, b, d):
    """
    Construit les matrices d'evaluation mod 2 et mod 3 des M = (d+1)(d+2)/2 - 1
    monomes binomiaux C(S1, k1) * C(S2, k2) (1 <= k1 + k2 <= d) sur tous les
    (a + 1)*(b + 1) - 1 couples non nuls (s1, s2) in {0..a} x {0..b}.
    """
    exps = []
    for deg in range(1, d + 1):
        for k1 in range(deg, -1, -1):
            k2 = deg - k1
            exps.append((k1, k2))
    M = len(exps)
    N = (a + 1) * (b + 1) - 1
    mat2 = np.zeros((N, M), dtype=np.uint8)
    mat3 = np.zeros((N, M), dtype=np.int64)
    row = 0
    for s1 in range(a + 1):
        for s2 in range(b + 1):
            if s1 == 0 and s2 == 0:
                continue
            for col, (k1, k2) in enumerate(exps):
                val = binom_int(s1, k1) * binom_int(s2, k2)
                mat2[row, col] = val & 1
                mat3[row, col] = val % 3
            row += 1
    return mat2, mat3, exps


def build_3block_grid_matrices(a, b, c, d):
    """Monomes binomiaux C(S1, k1)*C(S2, k2)*C(S3, k3) sur {0..a}x{0..b}x{0..c}\{(0,0,0)}."""
    exps = []
    for deg in range(1, d + 1):
        for k1 in range(deg, -1, -1):
            for k2 in range(deg - k1, -1, -1):
                k3 = deg - k1 - k2
                exps.append((k1, k2, k3))
    M = len(exps)
    N = (a + 1) * (b + 1) * (c + 1) - 1
    mat2 = np.zeros((N, M), dtype=np.uint8)
    mat3 = np.zeros((N, M), dtype=np.int64)
    row = 0
    for s1 in range(a + 1):
        for s2 in range(b + 1):
            for s3 in range(c + 1):
                if s1 == 0 and s2 == 0 and s3 == 0:
                    continue
                for col, (k1, k2, k3) in enumerate(exps):
                    val = binom_int(s1, k1) * binom_int(s2, k2) * binom_int(s3, k3)
                    mat2[row, col] = val & 1
                    mat3[row, col] = val % 3
                row += 1
    return mat2, mat3, exps


@njit(cache=True)
def sa_grid_worker(mat2, mat3, n_steps, seed):
    np.random.seed(seed)
    N, M = mat2.shape
    c3 = np.random.randint(0, 3, size=M)
    val3 = np.zeros(N, np.int64)
    for col in range(M):
        cc = c3[col]
        if cc != 0:
            for r in range(N):
                val3[r] += cc * mat3[r, col]
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
                    sub[pos, col] = mat2[r, col]
                pos += 1
        inc, _ = f2_defect_and_solve(sub)
        return inc, nz

    cur_inc, cur_nz = eval_def(val3)
    best_inc, best_nz = cur_inc, cur_nz
    best_c3 = c3.copy()
    if best_inc == 0:
        return 0, best_nz, best_c3

    temp = 1.5
    for step in range(n_steps):
        two_flip = (np.random.randint(0, 3) == 0)
        col1 = np.random.randint(0, M)
        d1 = 1 if np.random.randint(0, 2) == 0 else 2
        col2 = np.random.randint(0, M) if two_flip else -1
        if col2 == col1:
            col2 = -1
        d2 = 1 if np.random.randint(0, 2) == 0 else 2

        for r in range(N):
            val3[r] = (val3[r] + d1 * mat3[r, col1]) % 3
        if col2 >= 0:
            for r in range(N):
                val3[r] = (val3[r] + d2 * mat3[r, col2]) % 3

        inc, nz = eval_def(val3)
        diff = (inc - cur_inc) + 0.0005 * (nz - cur_nz)
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
                val3[r] = (val3[r] + id1 * mat3[r, col1]) % 3
            if col2 >= 0:
                id2 = (3 - d2) % 3
                for r in range(N):
                    val3[r] = (val3[r] + id2 * mat3[r, col2]) % 3
        temp *= 0.9994
    return best_inc, best_nz, best_c3


@njit(parallel=True, cache=True)
def parallel_grid_hunt(mat2, mat3, n_trials, n_steps, seed0):
    M = mat2.shape[1]
    incs = np.zeros(n_trials, np.int64)
    nzs = np.zeros(n_trials, np.int64)
    c3s = np.zeros((n_trials, M), np.int64)
    for t in prange(n_trials):
        inc, nz, c3 = sa_grid_worker(mat2, mat3, n_steps, seed0 + 131 * t)
        incs[t] = inc
        nzs[t] = nz
        c3s[t] = c3
    return incs, nzs, c3s


def verify_and_format(mat2, mat3, exps, c3):
    val3 = (mat3 @ c3) % 3
    z3_idx = np.nonzero(val3 == 0)[0]
    inc, sol2 = f2_defect_and_solve(mat2[z3_idx])
    assert inc == 0
    val2 = (mat2 @ sol2) % 2
    assert int(((val2 == 0) & (val3 == 0)).sum()) == 0
    return sol2, int(len(z3_idx))


def hunt_case(label, mat2, mat3, exps, n_trials=128, n_steps=6000, seed0=100):
    ts = time.perf_counter()
    incs, nzs, c3s = parallel_grid_hunt(mat2, mat3, n_trials, n_steps, seed0)
    dt = time.perf_counter() - ts
    bi = int(np.argmin(incs))
    best_inc = int(incs[bi])
    if best_inc == 0:
        sol2, nz = verify_and_format(mat2, mat3, exps, c3s[bi])
        print(f"   [RECORD TROUVE !] {label:<34} -> DEFAUT = 0 (VERIFIE 100%, |Z_3|={nz:>4}) [{dt:.2f}s]",
              flush=True)
        return True, {"c2": sol2.tolist(), "c3": c3s[bi].tolist(), "exps": exps, "z3": nz}
    else:
        print(f"   [bloque]          {label:<34} -> meilleur defaut F_2 = {best_inc:>2} (|Z_3|={int(nzs[bi]):>4}) [{dt:.2f}s]",
              flush=True)
        return False, {"best_defect": best_inc}


def main():
    t0 = time.time()
    records = {}

    print("== 1. DEGRE d = 3, 3 BLOCS (objectif : battre N_sym(3) = 27) ==", flush=True)
    for part in [(18, 8, 2), (20, 6, 2), (16, 10, 2), (14, 10, 4), (10, 10, 8), (24, 2, 2), (26, 1, 1)]:
        n = sum(part)
        m2, m3, exps = build_3block_grid_matrices(*part, 3)
        ok, info = hunt_case(f"d=3, n={n}, 3-blocs {part}", m2, m3, exps, n_trials=96, n_steps=5000, seed0=3000 + n)
        if ok:
            records[f"d3_n{n}_{part}"] = info

    print("\n== 2. DEGRE d = 4, 2 BLOCS (objectif : battre N_sym(4) = 36) ==", flush=True)
    for n_target in (37, 38, 40, 44):
        hit = False
        for b in (1, 2, 4, 8, 12, n_target // 2):
            a = n_target - b
            if a < b:
                continue
            m2, m3, exps = build_2block_grid_matrices(a, b, 4)
            ok, info = hunt_case(f"d=4, n={n_target}, 2-blocs ({a},{b})", m2, m3, exps,
                                 n_trials=128, n_steps=6000, seed0=4000 + 37 * n_target + b)
            if ok:
                records[f"d4_n{n_target}_{a}_{b}"] = info
                hit = True
                break
        if not hit and n_target >= 40:
            break

    print("\n== 3. DEGRE d = 5, 2 BLOCS (objectif : battre N_sym(5) = 53) ==", flush=True)
    for n_target in (54, 56, 60):
        hit = False
        for b in (1, 2, 4, 8, 16, n_target // 2):
            a = n_target - b
            if a < b:
                continue
            m2, m3, exps = build_2block_grid_matrices(a, b, 5)
            ok, info = hunt_case(f"d=5, n={n_target}, 2-blocs ({a},{b})", m2, m3, exps,
                                 n_trials=96, n_steps=5000, seed0=5000 + 41 * n_target + b)
            if ok:
                records[f"d5_n{n_target}_{a}_{b}"] = info
                hit = True
                break
        if not hit and n_target >= 56:
            break

    with open(os.path.join(RESULTS_DIR, "super_quad_mod6_results.json"), "w", encoding="utf-8") as fh:
        json.dump(records, fh, indent=2)
    print(f"\nTemps total : {time.time() - t0:.1f} s | Nouveaux records trouves : {list(records.keys())}",
          flush=True)


if __name__ == "__main__":
    main()
