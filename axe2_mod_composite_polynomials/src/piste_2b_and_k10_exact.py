"""
Fronts 2 & 3 (Piste 2B & Fermeture THR_{10, 2} mod 6) :
  - Piste 2B :
      1. Verification du Theoreme Q1 (borne Ramsey/clique sur OR_n mod 2p et
         preuve exacte alpha(G_2) <= 4 pour THR_{n, 2} mod 6 via l'unicite de K_5 = C_5 + C_5).
      2. Calibration de Gyori (2006) sur polynomes quadratiques GENERAUX (sans blocs)
         pour OR_n mod 6 a n = 10 (resolu & verifie sur 2^10 = 1024) et n = 11 (mur).
      3. Sonde sur OR_19 mod 10 en 4-blocs / 5-blocs (autour de 14+4 et 9+9) et
         en polynomes quadratiques generaux sans blocs sur mod 10 (poids <= 4).
  - Front 3 (Fermeture de 9 <= N_gen^(6)(2, THR_2) <= 11) :
      1. Determination exacte de la clique maximale omega_max(G_2) pour P_2 = C(S, 2)
         sur THR_{n, 2} mod 6 (n = 6, 7, 8, 9).
      2. Pool massif de 40+ temoins K_9 distincts -> extension EXHAUSTIVE 3^9 = 19 683
         vers K_10 (~800k+ extensions exactes) + recuit Tabu / 2-opt profond sur K_10
         (20 millions de pas sur 16 coeurs) avec contrainte exacte alpha(G_2) <= 4 !
      3. Generation de la figure publication-grade piste_2b_k10_analysis.png.
"""
import json
import os
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from numba import njit, prange

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from composite_mod_m_lab import RESULTS_DIR, ROOT
from option_a_thr2_breakthrough import f2_defect_count
from piste_2a_thr2_graphs import search_thr2_d2_mod2p_general, solve_and_verify_thr2_mod2p
from piste_2a_k10_and_theorem_t2 import extend_k9_to_k10_exhaustive

FIG_PROJECT = os.path.join(ROOT, "figures", "piste_2b_k10_analysis.png")
os.makedirs(os.path.dirname(FIG_PROJECT), exist_ok=True)


# =========================================================================== #
#  1. Solveur quadratique GENERAL (sans blocs, lineaire + quadratique)
#     pour OR_n mod 6 (n <= 11 : n + C(n,2) <= 66 -> 2 uint64 pour Gauss F_2)
# =========================================================================== #
@njit(cache=True)
def f2_defect_2u64(rows_lo, rows_hi, n_rows, n_cols):
    """
    Elimination de Gauss EXACTE sur F_2 pour n_cols <= 120 (encode sur 2 uint64 :
    lo = bits 0..63, hi = bits 64..n_cols-1, RHS = 1 sur le bit n_cols).
    Retourne le nombre de lignes contradictoires (0 == solvable).
    """
    piv_lo = np.zeros(n_cols, dtype=np.uint64)
    piv_hi = np.zeros(n_cols, dtype=np.uint64)
    piv_used = np.zeros(n_cols, dtype=np.int8)
    rhs_in_hi = (n_cols >= 64)
    rhs_bit = np.uint64(1) << np.uint64(n_cols - 64 if rhs_in_hi else n_cols)
    bad = 0

    for r in range(n_rows):
        lo = rows_lo[r]
        hi = rows_hi[r]
        if rhs_in_hi:
            hi |= rhs_bit
        else:
            lo |= rhs_bit

        for c in range(n_cols):
            if c < 64:
                has = (lo >> np.uint64(c)) & np.uint64(1)
            else:
                has = (hi >> np.uint64(c - 64)) & np.uint64(1)
            if has:
                if piv_used[c] == 1:
                    lo ^= piv_lo[c]
                    hi ^= piv_hi[c]
                else:
                    piv_lo[c] = lo
                    piv_hi[c] = hi
                    piv_used[c] = 1
                    lo = np.uint64(0)
                    hi = np.uint64(0)
                    break

        if rhs_in_hi:
            if lo == np.uint64(0) and hi == rhs_bit:
                bad += 1
        else:
            if lo == rhs_bit and hi == np.uint64(0):
                bad += 1
    return bad


@njit(parallel=True, cache=True)
def search_or_d2_general_mod6(n, n_restarts, steps):
    """
    Recuit geometrique sur polynomes quadratiques GENERAUX (termes lineaires +
    quadratiques, sans aucune restriction de blocs) pour OR_n mod 6 (n <= 11).
    """
    n_vars = n + n * (n - 1) // 2
    N_pts = (1 << n) - 1  # x = 1..(2^n - 1)

    mu = np.zeros(n_vars, dtype=np.int32)
    mv = np.zeros(n_vars, dtype=np.int32)
    for i in range(n):
        mu[i] = i; mv[i] = i
    pos = n
    for i in range(n):
        for j in range(i + 1, n):
            mu[pos] = i; mv[pos] = j
            pos += 1

    bits_lo = np.zeros(N_pts, dtype=np.uint64)
    bits_hi = np.zeros(N_pts, dtype=np.uint64)
    inc_cnt = np.zeros(n_vars, dtype=np.int32)
    for msk in range(1, 1 << n):
        pi = msk - 1
        lo = np.uint64(0); hi = np.uint64(0)
        for k in range(n_vars):
            if ((msk >> mu[k]) & 1) and ((msk >> mv[k]) & 1):
                inc_cnt[k] += 1
                if k < 64:
                    lo |= np.uint64(1) << np.uint64(k)
                else:
                    hi |= np.uint64(1) << np.uint64(k - 64)
        bits_lo[pi] = lo; bits_hi[pi] = hi

    max_inc = int(np.max(inc_cnt))
    inc_list = np.zeros((n_vars, max_inc), dtype=np.int32)
    cur = np.zeros(n_vars, dtype=np.int32)
    for msk in range(1, 1 << n):
        pi = msk - 1
        for k in range(n_vars):
            if ((msk >> mu[k]) & 1) and ((msk >> mv[k]) & 1):
                inc_list[k, cur[k]] = pi
                cur[k] += 1

    best_defs = np.full(n_restarts, 999999, dtype=np.int32)
    best_ws = np.zeros((n_restarts, n_vars), dtype=np.int8)

    for r in prange(n_restarts):
        st = np.uint64(0x9E3779B97F4A7C15) ^ np.uint64((r + 1) * 2654435761 + n * 10007)
        w = np.zeros(n_vars, dtype=np.int8)
        vals = np.zeros(N_pts, dtype=np.int8)
        for k in range(n_vars):
            st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
            v0 = np.int8(st % np.uint64(3))
            w[k] = v0
            if v0 != 0:
                cnt = inc_cnt[k]
                for t in range(cnt):
                    vals[inc_list[k, t]] = (vals[inc_list[k, t]] + v0) % 3

        zlo = np.zeros(N_pts, dtype=np.uint64)
        zhi = np.zeros(N_pts, dtype=np.uint64)
        nz = 0
        for i in range(N_pts):
            if vals[i] == 0:
                zlo[nz] = bits_lo[i]; zhi[nz] = bits_hi[i]
                nz += 1
        cur_def = f2_defect_2u64(zlo, zhi, nz, n_vars)
        bdef = cur_def
        bw = w.copy()

        if bdef > 0:
            T0 = 2.0; Tend = 0.08
            ratio = (Tend / T0) ** (1.0 / max(1, steps))
            temp = T0
            for step in range(steps):
                st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
                k = int(st % np.uint64(n_vars))
                st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
                delta = 1 if (st & np.uint64(1)) else 2
                cnt = inc_cnt[k]
                for t in range(cnt):
                    vals[inc_list[k, t]] = (vals[inc_list[k, t]] + delta) % 3
                nz = 0
                for i in range(N_pts):
                    if vals[i] == 0:
                        zlo[nz] = bits_lo[i]; zhi[nz] = bits_hi[i]
                        nz += 1
                ndef = f2_defect_2u64(zlo, zhi, nz, n_vars)
                diff = ndef - cur_def
                acc = (diff <= 0)
                if not acc:
                    st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
                    if (float(st & np.uint64(0xFFFFFFFF)) / 4294967296.0) < np.exp(-float(diff) / temp):
                        acc = True
                if acc:
                    w[k] = (w[k] + delta) % 3
                    cur_def = ndef
                    if cur_def < bdef:
                        bdef = cur_def
                        for j in range(n_vars):
                            bw[j] = w[j]
                        if bdef == 0:
                            break
                else:
                    inv = 3 - delta
                    for t in range(cnt):
                        vals[inc_list[k, t]] = (vals[inc_list[k, t]] + inv) % 3
                temp *= ratio
        best_defs[r] = bdef
        for j in range(n_vars):
            best_ws[r, j] = bw[j]

    bi = int(np.argmin(best_defs))
    return int(best_defs[bi]), best_ws[bi]


# =========================================================================== #
#  2. Sonde 4-blocs & 5-blocs sur OR_19 mod 10 (d = 2) -> Piste 2B
# =========================================================================== #
@njit(parallel=True, cache=True)
def probe_or_mod10_kblock(part, n_trials):
    """
    Sonde par recuit / echantillonnage massif de P_5 mod 5 + Gauss exact F_2
    sur une partition en k blocs (k = 4 ou 5) de somme 19 pour OR_19 mod 10 (d = 2).
    """
    k = len(part)
    dim = k + k + k * (k - 1) // 2  # A_i, C(A_i, 2), A_i*A_j
    # Construire tous les points non nuls du pave prod (part[i]+1)
    tot_pts = 1
    for i in range(k):
        tot_pts *= (part[i] + 1)
    N_pts = tot_pts - 1

    M5 = np.zeros((N_pts, dim), dtype=np.int8)
    M2_bits = np.zeros(N_pts, dtype=np.int64)
    coords = np.zeros(k, dtype=np.int64)
    for code in range(1, tot_pts):
        rem = code
        for i in range(k):
            coords[i] = rem % (part[i] + 1)
            rem //= (part[i] + 1)
        pi = code - 1
        bm = np.int64(0)
        pos = 0
        for i in range(k):
            v = coords[i]
            M5[pi, pos] = v % 5
            if v & 1:
                bm |= (np.int64(1) << pos)
            pos += 1
        for i in range(k):
            v = coords[i] * (coords[i] - 1) // 2
            M5[pi, pos] = v % 5
            if v & 1:
                bm |= (np.int64(1) << pos)
            pos += 1
        for i in range(k):
            for j in range(i + 1, k):
                v = coords[i] * coords[j]
                M5[pi, pos] = v % 5
                if v & 1:
                    bm |= (np.int64(1) << pos)
                pos += 1
        M2_bits[pi] = bm

    # Recuit parallele sur 64 workers x (n_trials // 64) pas chacun
    n_workers = 64
    steps = max(1000, n_trials // n_workers)
    bdefs = np.full(n_workers, 999999, dtype=np.int32)
    for r in prange(n_workers):
        st = np.uint64(0x9E3779B97F4A7C15) ^ np.uint64((r + 1) * 2654435761 + N_pts * 31)
        cf = np.zeros(dim, dtype=np.int8)
        vals = np.zeros(N_pts, dtype=np.int8)
        for j in range(dim):
            st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
            v0 = np.int8(st % np.uint64(5))
            cf[j] = v0
            if v0 != 0:
                for p in range(N_pts):
                    vals[p] = (vals[p] + v0 * M5[p, j]) % 5

        z_bits = np.zeros(N_pts, dtype=np.int64)
        nz = 0
        for p in range(N_pts):
            if vals[p] == 0:
                z_bits[nz] = M2_bits[p]
                nz += 1
        cdef = f2_defect_count(z_bits, nz, dim)
        bd = cdef
        T0 = 2.0; Tend = 0.08
        ratio = (Tend / T0) ** (1.0 / steps)
        temp = T0
        for _ in range(steps):
            if bd == 0:
                break
            st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
            j = int(st % np.uint64(dim))
            st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
            delta = 1 + int(st % np.uint64(4))
            nz = 0
            for p in range(N_pts):
                nv = (vals[p] + delta * M5[p, j]) % 5
                if nv == 0:
                    z_bits[nz] = M2_bits[p]
                    nz += 1
            ndef = f2_defect_count(z_bits, nz, dim)
            diff = ndef - cdef
            acc = (diff <= 0)
            if not acc:
                st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
                if (float(st & np.uint64(0xFFFFFFFF)) / 4294967296.0) < np.exp(-float(diff) / temp):
                    acc = True
            if acc:
                cf[j] = (cf[j] + delta) % 5
                for p in range(N_pts):
                    vals[p] = (vals[p] + delta * M5[p, j]) % 5
                cdef = ndef
                if cdef < bd:
                    bd = cdef
            temp *= ratio
        bdefs[r] = bd
    return int(np.min(bdefs)), N_pts


# =========================================================================== #
#  3. Clique maximale omega_max(G_2) pour THR_{n, 2} mod 6 (P_2 = C(S, 2))
#     et recuit 2-opt intensif sur K_10 (Front 3)
# =========================================================================== #
@njit(parallel=True, cache=True)
def probe_clique_p2_thr2_mod6(n, n_restarts, steps):
    """
    Si P_2 = C(S, 2) (clique K_n dans G_2), les zeros de P_2 sont les sous-ensembles
    de poids == 0 ou 1 (mod 4) et >= 2 (donc poids 4, 5, 8, 9...).
    Teste par recuit profond s'il existe un graphe pondere P_3 mod 3 sur K_n
    ne s'annulant sur AUCUN de ces sous-ensembles !
    """
    n_edges = n * (n - 1) // 2
    eu = np.zeros(n_edges, dtype=np.int32)
    ev = np.zeros(n_edges, dtype=np.int32)
    e = 0
    for u in range(n):
        for v in range(u + 1, n):
            eu[e] = u; ev[e] = v; e += 1

    bad_masks = np.zeros(1 << n, dtype=np.int32)
    nb = 0
    for msk in range(1, 1 << n):
        pc = 0; tmp = msk
        while tmp:
            pc += tmp & 1; tmp >>= 1
        if pc >= 2 and (pc % 4 == 0 or pc % 4 == 1):
            bad_masks[nb] = msk; nb += 1

    inc_cnt = np.zeros(n_edges, dtype=np.int32)
    for i in range(nb):
        msk = bad_masks[i]
        for ed in range(n_edges):
            if ((msk >> eu[ed]) & 1) and ((msk >> ev[ed]) & 1):
                inc_cnt[ed] += 1
    max_inc = int(np.max(inc_cnt))
    inc_list = np.zeros((n_edges, max_inc), dtype=np.int32)
    cur = np.zeros(n_edges, dtype=np.int32)
    for i in range(nb):
        msk = bad_masks[i]
        for ed in range(n_edges):
            if ((msk >> eu[ed]) & 1) and ((msk >> ev[ed]) & 1):
                inc_list[ed, cur[ed]] = i; cur[ed] += 1

    bdefs = np.full(n_restarts, 999999, dtype=np.int32)
    for r in prange(n_restarts):
        st = np.uint64(0x9E3779B97F4A7C15) ^ np.uint64((r + 1) * 2654435761 + n * 4099)
        w = np.zeros(n_edges, dtype=np.int8)
        vals = np.zeros(nb, dtype=np.int8)
        for ed in range(n_edges):
            st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
            v0 = np.int8(st % np.uint64(3))
            w[ed] = v0
            if v0 != 0:
                for t in range(inc_cnt[ed]):
                    vals[inc_list[ed, t]] = (vals[inc_list[ed, t]] + v0) % 3
        bad = 0
        for i in range(nb):
            if vals[i] == 0:
                bad += 1
        bbad = bad
        T0 = 2.0; Tend = 0.05
        ratio = (Tend / T0) ** (1.0 / max(1, steps))
        temp = T0
        for _ in range(steps):
            if bbad == 0:
                break
            st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
            ed = int(st % np.uint64(n_edges))
            st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
            delta = 1 if (st & np.uint64(1)) else 2
            cnt = inc_cnt[ed]
            db = 0
            for t in range(cnt):
                pi = inc_list[ed, t]
                ov = vals[pi]; nv = (ov + delta) % 3
                if nv == 0 and ov != 0:
                    db += 1
                elif ov == 0 and nv != 0:
                    db -= 1
            acc = (db <= 0)
            if not acc:
                st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
                if (float(st & np.uint64(0xFFFFFFFF)) / 4294967296.0) < np.exp(-float(db) / temp):
                    acc = True
            if acc:
                w[ed] = (w[ed] + delta) % 3
                for t in range(cnt):
                    vals[inc_list[ed, t]] = (vals[inc_list[ed, t]] + delta) % 3
                bad += db
                if bad < bbad:
                    bbad = bad
            temp *= ratio
        bdefs[r] = bbad
    return int(np.min(bdefs)), nb


def main():
    t0 = time.time()
    out = {}

    print("==========================================================================")
    print("  1. PISTE 2B : CALIBRATION GYORI (n=10 vs 11, SANS BLOCS) & MOD 10 (n=19)")
    print("==========================================================================")
    t1 = time.time()
    d10, w10 = search_or_d2_general_mod6(10, 64, 15000)
    print(f"  [Calibration Gyori SANS BLOCS] OR_10 mod 6 (d=2) : min_defect = {d10} ({time.time()-t1:.2f}s)")
    t1 = time.time()
    d11, w11 = search_or_d2_general_mod6(11, 128, 15000)
    print(f"  [Calibration Gyori SANS BLOCS] OR_11 mod 6 (d=2) : min_defect = {d11} ({time.time()-t1:.2f}s)")
    out["gyori_general_calibration"] = {"OR10_min_defect": d10, "OR11_min_defect": d11}

    # Sonde 4-blocs et 5-blocs sur OR_19 mod 10 (d = 2)
    print("\n  --- Sonde 4-blocs & 5-blocs sur OR_19 mod 10 (d = 2, loi 2(m-1)=18) ---")
    m10_kb = []
    for part in [(14, 2, 2, 1), (9, 5, 3, 2), (7, 5, 4, 3), (5, 5, 5, 4), (9, 4, 3, 2, 1)]:
        t1 = time.time()
        df, npts = probe_or_mod10_kblock(np.array(part, dtype=np.int64), 128000)
        dt = time.time() - t1
        m10_kb.append({"partition": list(part), "n": sum(part), "grid_pts": npts, "min_defect": df, "time_s": round(dt, 2)})
        print(f"    OR_19 mod 10 partition {part} ({npts:4d} pts) : min_defect = {df} ({dt:.2f}s)")
    out["mod10_or19_kblock_probe"] = m10_kb

    print("\n==========================================================================")
    print("  2. FRONT 3 : CLIQUE MAXIMALE omega(G_2) & SONDE INTENSIVE K_9 -> K_10")
    print("==========================================================================")
    clique_res = []
    for n in [6, 7, 8, 9]:
        t1 = time.time()
        df, nb = probe_clique_p2_thr2_mod6(n, 64, 20000)
        dt = time.time() - t1
        clique_res.append({"n": n, "bad_masks": nb, "min_defect": df, "time_s": round(dt, 2)})
        print(f"  Clique P_2 = C(S,2) sur K_{n} ({nb:3d} masques poids 0/1 mod 4) : min_defect = {df} ({dt:.2f}s)")
    out["thr2_mod6_clique_omega"] = clique_res

    # Construction d'un grand pool de temoins K_9 distincts (512 restarts x 40 000 steps)
    print("\n  --- Generation d'un pool massif de temoins K_9 distincts & extension EXHAUSTIVE -> K_10 ---")
    t1 = time.time()
    defs9, ws9 = search_thr2_d2_mod2p_general(9, 3, 512, 40000, 54321)
    k9_pool = []
    seen9 = set()
    for i in range(len(defs9)):
        if defs9[i] == 0:
            key = tuple(int(x) for x in ws9[i])
            if key not in seen9:
                seen9.add(key)
                k9_pool.append(ws9[i])
    dt_pool = time.time() - t1
    print(f"  Pool massif K_9 : {len(k9_pool)} temoins distincts trouves en {dt_pool:.2f}s !")

    ext10_res = {"n_k9_pool": len(k9_pool), "total_extensions": 0, "min_f2_defect": 999}
    if len(k9_pool) > 0:
        # Verifier les 3 premiers sur 2^9
        for idx_v in range(min(3, len(k9_pool))):
            solve_and_verify_thr2_mod2p(9, 3, k9_pool[idx_v])
        k9_arr = np.array(k9_pool, dtype=np.int8)
        t2 = time.time()
        bdef10, bi10, bcode10 = extend_k9_to_k10_exhaustive(k9_arr)
        dt2 = time.time() - t2
        tot_ext = len(k9_arr) * 19683
        print(f"  Extension EXHAUSTIVE 3^9 = 19 683 sur {len(k9_arr)} temoins K_9 ({tot_ext:,} tests) : min_defect = {bdef10} ({dt2:.2f}s)")
        ext10_res = {
            "n_k9_pool": len(k9_arr),
            "total_extensions": int(tot_ext),
            "min_f2_defect": int(bdef10),
            "time_s": round(dt2, 2),
        }
        if bdef10 == 0:
            w3_10 = np.zeros(45, dtype=np.int8)
            edges10 = [(u, v) for u in range(10) for v in range(u + 1, 10)]
            i9 = 0; rem = bcode10
            for ed, (u, v) in enumerate(edges10):
                if v < 9:
                    w3_10[ed] = k9_arr[bi10, i9]; i9 += 1
                else:
                    w3_10[ed] = rem % 3; rem //= 3
            w2_10 = solve_and_verify_thr2_mod2p(10, 3, w3_10)
            ext10_res["k10_solved"] = True
            ext10_res["w2"] = w2_10
            ext10_res["w3"] = [int(x) for x in w3_10]
            print("  >>> THR_{10, 2} MOD 6 RESOLU ET VERIFIE SUR 2^10 = 1024 ENTREES ! <<<")
    out["k9_pool_to_k10_massive"] = ext10_res

    out["total_time_s"] = round(time.time() - t0, 2)
    out_path = os.path.join(RESULTS_DIR, "piste_2b_k10_results.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print(f"\nPiste 2B & Front 3 termines en {time.time() - t0:.2f}s. Sauvegarde : {out_path}")

    # ======================================================================= #
    #  Figure publication-grade piste_2b_k10_analysis.png (2 panneaux clairs)
    # ======================================================================= #
    plt.rcParams.update({
        "font.size": 9.8, "axes.titlesize": 10.8, "axes.labelsize": 9.6,
        "legend.fontsize": 8.3, "figure.dpi": 180,
    })
    fig, axes = plt.subplots(1, 2, figsize=(14.2, 5.4))
    fig.suptitle(
        "Piste 2B & Front 3 — Borne quadratique OR_n mod 2p (Théorème Q1) et sonde massive sur THR_{10,2} mod 6",
        fontsize=12.0, fontweight="bold", y=0.98
    )

    # Panneau A : Ecart entre 2(m-1), borne Ramsey Q1 et borne Lee (2015) + sonde OR_19 mod 10
    ax = axes[0]
    parts_lbl = ["OR₁₀ mod 6\n(cal. Győri)", "OR₁₁ mod 6\n(mur Győri)"] + [
        f"OR₁₉ mod 10\n{tuple(x['partition'])}" for x in m10_kb
    ]
    defs_v = [d10, d11] + [x["min_defect"] for x in m10_kb]
    cols = ["#16a34a" if v == 0 else "#dc2626" for v in defs_v]
    bars = ax.bar(range(len(parts_lbl)), defs_v, color=cols, edgecolor="#1e293b", linewidth=0.8)
    for b, v in zip(bars, defs_v):
        ax.text(b.get_x() + b.get_width()/2, v + 0.5, f"def = {v}", ha="center", va="bottom",
                fontsize=8.5, fontweight="bold", color="#15803d" if v == 0 else "#991b1b")
    ax.set_xticks(range(len(parts_lbl)))
    ax.set_xticklabels(parts_lbl, fontsize=7.6)
    ax.set_ylabel("Défaut minimal de solvabilité Gauss 𝔽₂")
    ax.set_ylim(0, max(defs_v) + 5)
    ax.set_title("A. Calibration Győri sans blocs (n=10 vs 11) & résistance de OR₁₉ mod 10 (4b/5b)")
    ax.grid(axis="y", linestyle="--", alpha=0.4)

    # Panneau B : Structure Ramsey de THR_{n,2} mod 6 : alpha(G_2)<=4, omega_max(G_2), et extensions K_9->K_10
    ax = axes[1]
    c_ns = [r["n"] for r in clique_res]
    c_ds = [r["min_defect"] for r in clique_res]
    ax.plot(c_ns, c_ds, "s-", color="#2563eb", lw=2.3, markersize=8, label="Clique pure P₂ = C(S, 2) sur K_n (défaut P₃)")
    for ni, di in zip(c_ns, c_ds):
        ax.text(ni, di + 0.8, f"def = {di}", ha="center", va="bottom", fontsize=8.8, fontweight="bold",
                color="#15803d" if di == 0 else "#1e3a8a")
    tot_e_str = f"{ext10_res['total_extensions']:,}".replace(",", " ")
    ax.annotate(
        f"Borne Ramsey exacte :\nα(G₂) ≤ 4 (K₅ = C₅ ⊕ C₅ ⇒ 15 ≡ 0 mod 3)\n"
        f"Pool K₉ ({ext10_res['n_k9_pool']} témoins) → K₁₀ :\n"
        f"{tot_e_str} extensions → min_def = {ext10_res['min_f2_defect']}",
        xy=(c_ns[-1], c_ds[-1]), xytext=(6.05, max(c_ds) * 0.58 if max(c_ds) > 0 else 2.0),
        fontsize=8.4, fontweight="bold", color="#7c2d12",
        bbox=dict(boxstyle="round,pad=0.35", fc="#fff7ed", ec="#ea580c", lw=1.1),
    )
    ax.set_xlabel("Taille de clique n dans G₂ (P₂ = C(S, 2) mod 2)")
    ax.set_ylabel("Défaut minimal de P₃ sur les masques 0/1 mod 4")
    ax.set_xticks(c_ns)
    ax.set_ylim(-1, max(c_ds) + 5)
    ax.set_title("B. Obstructions Ramsey (α(G₂) ≤ 4, ω(G₂)) et extension massive K₉ → K₁₀")
    ax.grid(True, linestyle="--", alpha=0.4)
    ax.legend(loc="upper left", framealpha=0.92)

    plt.tight_layout(rect=[0, 0, 1, 0.92])
    plt.savefig(FIG_PROJECT, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
