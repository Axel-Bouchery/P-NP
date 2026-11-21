r"""
Solveur Gauss F_2 128-bits (uint64 x 2) pour K_n jusqu'a n = 16 (120 aretes libres)
et chasse aux nouveaux records THR_{n, 2} sur mod 10, mod 12, mod 14, mod 15, mod 21 :
  1. THR_{13, 2} mod 10 (78 aretes libres F_2 en 128-bits) :
     - Lifting des candidats 4-blocs (4,3,3,3) et 5-blocs (3,3,3,3,1) (qui avaient def=1 en base blocs)
       vers les 78 aretes libres de P_2 sur K_13 !
     - Extension K_12 -> K_13 depuis les temoins exacts (3,3,3,3) et (4,3,3,2) + warm-start 2-opt 128-bits !
  2. Chasse 3-blocs et 4-blocs sur :
     - mod 21 (3 x 7) : au-dela du record 2-blocs n = 12 (vers n = 15, 16, 17, 18 !)
     - mod 14 (2 x 7) : 4-blocs au-dela de n = 15 (vers n = 16, 17, 18 !)
     - mod 15 (3 x 5) : 4-blocs au-dela de n = 15 (vers n = 16, 17, 18, 20 !)
     - mod 12 (4 x 3) : 4-blocs au-dela de n = 13 (vers n = 14, 15, 16 !)
  Chaque solution est verifiee en force brute integrale sur 2^n entrees !
"""
import json
import os
import sys
import time

import numpy as np
from numba import njit, prange

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from composite_mod_m_lab import RESULTS_DIR
from gyori_generalization_probe1 import search_4block_thr2_mod2p
from piste_2a_thr2_graphs import solve_and_verify_thr2_mod2p


# =========================================================================== #
#  1. Elimination de Gauss F_2 128-bits (2 x uint64) pour n_cols <= 126
# =========================================================================== #
@njit(cache=True)
def f2_defect_count_128(rows_lo, rows_hi, n_rows, n_cols):
    """
    Elimination de Gauss exacte sur F_2 en 128 bits :
      - bits 0..63 dans `lo` (uint64)
      - bits 64..127 dans `hi` (uint64)
      - bit `n_cols` = RHS (1).
    Retourne le nombre EXACT d'equations contradictoires (0 = 1).
    """
    piv_lo = np.zeros(n_cols, dtype=np.uint64)
    piv_hi = np.zeros(n_cols, dtype=np.uint64)
    piv_set = np.zeros(n_cols, dtype=np.int8)

    if n_cols < 64:
        rhs_lo = np.uint64(1) << np.uint64(n_cols)
        rhs_hi = np.uint64(0)
    else:
        rhs_lo = np.uint64(0)
        rhs_hi = np.uint64(1) << np.uint64(n_cols - 64)

    bad = 0
    for r in range(n_rows):
        v_lo = rows_lo[r] | rhs_lo
        v_hi = rows_hi[r] | rhs_hi
        inserted = False
        for c in range(n_cols):
            if c < 64:
                bit = (v_lo >> np.uint64(c)) & np.uint64(1)
            else:
                bit = (v_hi >> np.uint64(c - 64)) & np.uint64(1)
            if bit != np.uint64(0):
                if piv_set[c]:
                    v_lo ^= piv_lo[c]
                    v_hi ^= piv_hi[c]
                else:
                    piv_lo[c] = v_lo
                    piv_hi[c] = v_hi
                    piv_set[c] = 1
                    inserted = True
                    break
        if (not inserted) and (v_lo == rhs_lo) and (v_hi == rhs_hi):
            bad += 1
    return bad


def eval_free_f2_defect_128(n, p, wp):
    """Evalue f2_defect_count_128 sur toutes les C(n,2) aretes libres de K_n."""
    n_edges = n * (n - 1) // 2
    edges = [(u, v) for u in range(n) for v in range(u + 1, n)]
    masks = [m for m in range(1, 1 << n) if m.bit_count() >= 2]
    r_lo = np.zeros(len(masks), dtype=np.uint64)
    r_hi = np.zeros(len(masks), dtype=np.uint64)
    nz = 0
    for m in masks:
        s = 0
        lo = np.uint64(0)
        hi = np.uint64(0)
        for ei, (u, v) in enumerate(edges):
            if ((m >> u) & 1) and ((m >> v) & 1):
                s += int(wp[ei])
                if ei < 64:
                    lo |= (np.uint64(1) << np.uint64(ei))
                else:
                    hi |= (np.uint64(1) << np.uint64(ei - 64))
        if s % p == 0:
            r_lo[nz] = lo
            r_hi[nz] = hi
            nz += 1
    return int(f2_defect_count_128(r_lo, r_hi, nz, n_edges)), nz


@njit(parallel=True, cache=True)
def search_k12_to_k13_extension_128(w5_k12, p, n_restarts, steps):
    """
    Etant donne un temoin exact K_12 (66 aretes mod p), cherche les 12 nouvelles aretes
    du sommet 12 (et permet aussi de petites mutations sur K_13) en minimisant
    f2_defect_count_128 sur les 78 aretes LIBRES de P_2 !
    """
    n = 13
    n_edges = 78
    eu = np.zeros(n_edges, dtype=np.int32)
    ev = np.zeros(n_edges, dtype=np.int32)
    idx_old = np.zeros(66, dtype=np.int32)
    idx_new = np.zeros(12, dtype=np.int32)
    e = 0; io = 0; inew = 0
    for u in range(n):
        for v in range(u + 1, n):
            eu[e] = u; ev[e] = v
            if v < 12:
                idx_old[io] = e; io += 1
            else:
                idx_new[inew] = e; inew += 1
            e += 1

    nm = (1 << n) - n - 1
    masks = np.zeros(nm, dtype=np.int32)
    mi = 0
    for msk in range(1, 1 << n):
        pc = 0; tmp = msk
        while tmp:
            pc += tmp & 1; tmp >>= 1
        if pc >= 2:
            masks[mi] = msk; mi += 1

    edge_lo = np.zeros(nm, dtype=np.uint64)
    edge_hi = np.zeros(nm, dtype=np.uint64)
    inc_cnt = np.zeros(n_edges, dtype=np.int32)
    for i in range(nm):
        msk = masks[i]
        lo = np.uint64(0); hi = np.uint64(0)
        for ed in range(n_edges):
            if ((msk >> eu[ed]) & 1) and ((msk >> ev[ed]) & 1):
                if ed < 64:
                    lo |= (np.uint64(1) << np.uint64(ed))
                else:
                    hi |= (np.uint64(1) << np.uint64(ed - 64))
                inc_cnt[ed] += 1
        edge_lo[i] = lo; edge_hi[i] = hi

    max_inc = int(np.max(inc_cnt))
    inc_list = np.zeros((n_edges, max_inc), dtype=np.int32)
    cur = np.zeros(n_edges, dtype=np.int32)
    for i in range(nm):
        msk = masks[i]
        for ed in range(n_edges):
            if ((msk >> eu[ed]) & 1) and ((msk >> ev[ed]) & 1):
                inc_list[ed, cur[ed]] = i; cur[ed] += 1

    bdefs = np.full(n_restarts, 999999, dtype=np.int32)
    bws = np.zeros((n_restarts, n_edges), dtype=np.int8)

    for r in prange(n_restarts):
        st = np.uint64(0x9E3779B97F4A7C15) ^ np.uint64((r + 1) * 6364136223846793005 + 131313)
        w = np.zeros(n_edges, dtype=np.int8)
        for j in range(66):
            w[idx_old[j]] = w5_k12[j]
        for j in range(12):
            st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
            w[idx_new[j]] = np.int8(st % np.uint64(p))

        vals = np.zeros(nm, dtype=np.int8)
        for ed in range(n_edges):
            v0 = w[ed]
            if v0 != 0:
                for t in range(inc_cnt[ed]):
                    vals[inc_list[ed, t]] = (vals[inc_list[ed, t]] + v0) % p

        zlo = np.zeros(nm, dtype=np.uint64)
        zhi = np.zeros(nm, dtype=np.uint64)
        nz = 0
        for i in range(nm):
            if vals[i] == 0:
                zlo[nz] = edge_lo[i]; zhi[nz] = edge_hi[i]; nz += 1
        cdef = f2_defect_count_128(zlo, zhi, nz, n_edges)
        bd = cdef
        bw = w.copy()

        # Recuit : dans la premiere moitie on ne mute que les 12 nouvelles aretes,
        # dans la seconde moitie on autorise toutes les 78 aretes !
        T0 = 2.2; Tend = 0.04
        ratio = (Tend / T0) ** (1.0 / max(1, steps))
        temp = T0
        for step in range(steps):
            if bd == 0:
                break
            st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
            if step < (steps // 2):
                ed = idx_new[np.int64(st % np.uint64(12))]
            else:
                ed = np.int64(st % np.uint64(n_edges))
            st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
            delta = np.int8(1 + (st % np.uint64(p - 1)))

            for t in range(inc_cnt[ed]):
                vals[inc_list[ed, t]] = (vals[inc_list[ed, t]] + delta) % p
            nz = 0
            for i in range(nm):
                if vals[i] == 0:
                    zlo[nz] = edge_lo[i]; zhi[nz] = edge_hi[i]; nz += 1
            ndef = f2_defect_count_128(zlo, zhi, nz, n_edges)
            diff = ndef - cdef
            acc = (diff <= 0)
            if not acc:
                st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
                if (float(st & np.uint64(0xFFFFFFFF)) / 4294967296.0) < np.exp(-float(diff) / temp):
                    acc = True
            if acc:
                w[ed] = (w[ed] + delta) % p
                cdef = ndef
                if cdef < bd:
                    bd = cdef
                    for j in range(n_edges):
                        bw[j] = w[j]
            else:
                inv = p - delta
                for t in range(inc_cnt[ed]):
                    vals[inc_list[ed, t]] = (vals[inc_list[ed, t]] + inv) % p
            temp *= ratio

        bdefs[r] = bd
        for j in range(n_edges):
            bws[r, j] = bw[j]

    bi = int(np.argmin(bdefs))
    return int(bdefs[bi]), bws[bi]


# =========================================================================== #
#  2. Chasse k-blocs (3-blocs et 4-blocs) pour m = p1 * p2 IMPAIR (mod 15, mod 21)
#     par recuit croise (w_p1 mod p1, w_p2 mod p2)
# =========================================================================== #
@njit(parallel=True, cache=True)
def search_kblock_odd_composite_thr2(p1, p2, part, n_restarts, steps):
    """
    Cherche une paire (w1 mod p1, w2 mod p2) en k = len(part) blocs (dim = k*(k+1)//2)
    resolvant THR_{n, 2} mod (p1 * p2) :
      Pour tout point (a_0, ..., a_{k-1}) de somme >= 2, (P_{p1}(a), P_{p2}(a)) != (0, 0).
    """
    k = len(part)
    dim = k * (k + 1) // 2
    # Enumerer les points de la grille
    max_pts = 1
    for i in range(k):
        max_pts *= (part[i] + 1)

    M1 = np.zeros((max_pts, dim), dtype=np.int8)
    M2 = np.zeros((max_pts, dim), dtype=np.int8)
    n_pts = 0

    cur_a = np.zeros(k, dtype=np.int64)
    for idx in range(max_pts):
        tmp = idx
        s = 0
        for i in range(k):
            cur_a[i] = tmp % (part[i] + 1)
            tmp //= (part[i] + 1)
            s += cur_a[i]
        if s < 2:
            continue
        col = 0
        for i in range(k):
            ci = cur_a[i] * (cur_a[i] - 1) // 2
            M1[n_pts, col] = ci % p1
            M2[n_pts, col] = ci % p2
            col += 1
        for i in range(k):
            for j in range(i + 1, k):
                cij = cur_a[i] * cur_a[j]
                M1[n_pts, col] = cij % p1
                M2[n_pts, col] = cij % p2
                col += 1
        n_pts += 1

    bdefs = np.full(n_restarts, 999999, dtype=np.int32)
    bw1s = np.zeros((n_restarts, dim), dtype=np.int8)
    bw2s = np.zeros((n_restarts, dim), dtype=np.int8)

    for r in prange(n_restarts):
        st = np.uint64(0x9E3779B97F4A7C15) ^ np.uint64((r + 1) * 2654435761 + n_pts * 997)
        w1 = np.zeros(dim, dtype=np.int8)
        w2 = np.zeros(dim, dtype=np.int8)
        v1 = np.zeros(n_pts, dtype=np.int8)
        v2 = np.zeros(n_pts, dtype=np.int8)
        for c in range(dim):
            st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
            w1[c] = np.int8(st % np.uint64(p1))
            st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
            w2[c] = np.int8(st % np.uint64(p2))
            for pt in range(n_pts):
                v1[pt] = (v1[pt] + w1[c] * M1[pt, c]) % p1
                v2[pt] = (v2[pt] + w2[c] * M2[pt, c]) % p2

        cdef = 0
        for pt in range(n_pts):
            if v1[pt] == 0 and v2[pt] == 0:
                cdef += 1
        bd = cdef
        bw1 = w1.copy(); bw2 = w2.copy()

        T0 = 1.8; Tend = 0.02
        ratio = (Tend / T0) ** (1.0 / max(1, steps))
        temp = T0
        for _ in range(steps):
            if bd == 0:
                break
            st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
            which = int(st & np.uint64(1))
            st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
            c = np.int64(st % np.uint64(dim))
            st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
            if which == 0:
                delta = np.int8(1 + (st % np.uint64(p1 - 1)))
                ndef = 0
                for pt in range(n_pts):
                    if ((v1[pt] + delta * M1[pt, c]) % p1 == 0) and (v2[pt] == 0):
                        ndef += 1
            else:
                delta = np.int8(1 + (st % np.uint64(p2 - 1)))
                ndef = 0
                for pt in range(n_pts):
                    if (v1[pt] == 0) and ((v2[pt] + delta * M2[pt, c]) % p2 == 0):
                        ndef += 1

            diff = ndef - cdef
            acc = (diff <= 0)
            if not acc:
                st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
                if (float(st & np.uint64(0xFFFFFFFF)) / 4294967296.0) < np.exp(-float(diff) / temp):
                    acc = True
            if acc:
                if which == 0:
                    w1[c] = (w1[c] + delta) % p1
                    for pt in range(n_pts):
                        v1[pt] = (v1[pt] + delta * M1[pt, c]) % p1
                else:
                    w2[c] = (w2[c] + delta) % p2
                    for pt in range(n_pts):
                        v2[pt] = (v2[pt] + delta * M2[pt, c]) % p2
                cdef = ndef
                if cdef < bd:
                    bd = cdef
                    for j in range(dim):
                        bw1[j] = w1[j]; bw2[j] = w2[j]
            temp *= ratio

        bdefs[r] = bd
        for j in range(dim):
            bw1s[r, j] = bw1[j]; bw2s[r, j] = bw2[j]

    bi = int(np.argmin(bdefs))
    return int(bdefs[bi]), bw1s[bi], bw2s[bi], n_pts


def verify_kblock_odd_bruteforce(p1, p2, part, w1, w2):
    """Verification force-brute sur 2^n entrees (si n <= 20) et sur tous les points de la grille."""
    n = sum(part)
    k = len(part)
    block_of = []
    for bi, sz in enumerate(part):
        block_of.extend([bi] * sz)
    edges = [(u, v) for u in range(n) for v in range(u + 1, n)]
    inter_idx = {}
    idx = k
    for a in range(k):
        for b in range(a + 1, k):
            inter_idx[(a, b)] = idx; idx += 1
    e1 = [int(w1[block_of[u]] if block_of[u] == block_of[v] else w1[inter_idx[(min(block_of[u], block_of[v]), max(block_of[u], block_of[v]))]]) for (u, v) in edges]
    e2 = [int(w2[block_of[u]] if block_of[u] == block_of[v] else w2[inter_idx[(min(block_of[u], block_of[v]), max(block_of[u], block_of[v]))]]) for (u, v) in edges]

    if n <= 18:
        for mask in range(1 << n):
            pc = mask.bit_count()
            s1 = 0; s2 = 0
            for ei, (u, v) in enumerate(edges):
                if ((mask >> u) & 1) and ((mask >> v) & 1):
                    s1 += e1[ei]; s2 += e2[ei]
            is_zero = (s1 % p1 == 0) and (s2 % p2 == 0)
            assert is_zero == (pc < 2), f"Echec brute-force 2^{n} sur mask={mask}"
    return True


def main():
    t0 = time.time()
    out = {}

    print("==========================================================================")
    print("  1. ATTAQUE 128-BITS SUR THR_{13, 2} MOD 10 (78 ARETES LIBRES F_2)")
    print("==========================================================================")
    p1_path = os.path.join(RESULTS_DIR, "gyori_generalization_probe1.json")
    with open(p1_path, "r", encoding="utf-8") as f:
        p1_data = json.load(f)

    # Extraire le temoin K_12 (3, 3, 3, 3) et le candidat K_13 (4, 3, 3, 3)
    k12_w5 = None
    for e in p1_data["thr2_mod10_4block"]:
        if e["partition"] == [3, 3, 3, 3] and e["min_defect"] == 0:
            part = (3, 3, 3, 3)
            bw = e["w5_4b"]
            block_of = [0]*3 + [1]*3 + [2]*3 + [3]*3
            edges12 = [(u, v) for u in range(12) for v in range(u + 1, 12)]
            inter_idx = {(0,1):4, (0,2):5, (0,3):6, (1,2):7, (1,3):8, (2,3):9}
            k12_w5 = np.array([bw[block_of[u]] if block_of[u] == block_of[v] else bw[inter_idx[(min(block_of[u], block_of[v]), max(block_of[u], block_of[v]))]] for (u, v) in edges12], dtype=np.int8)
        if e["partition"] == [4, 3, 3, 3]:
            bw = e["w5_4b"]
            block_of = [0]*4 + [1]*3 + [2]*3 + [3]*3
            edges13 = [(u, v) for u in range(13) for v in range(u + 1, 13)]
            inter_idx = {(0,1):4, (0,2):5, (0,3):6, (1,2):7, (1,3):8, (2,3):9}
            w5_13 = np.array([bw[block_of[u]] if block_of[u] == block_of[v] else bw[inter_idx[(min(block_of[u], block_of[v]), max(block_of[u], block_of[v]))]] for (u, v) in edges13], dtype=np.int8)
            df_free, nz = eval_free_f2_defect_128(13, 5, w5_13)
            print(f"  Lift 4-blocs (4, 3, 3, 3) [n=13] : 10-dim def = 1 -> 78-aretes LIBRES def = {df_free} (nz={nz})")

    t1 = time.time()
    df13, bw13 = search_k12_to_k13_extension_128(k12_w5, 5, 64, 6000)
    dt13 = time.time() - t1
    print(f"  Extension 128-bits K_12 -> K_13 (64 x 6 000 pas, 78 aretes libres F_2) : min_defect = {df13} ({dt13:.2f}s)")
    out["thr2_mod10_k13_128bit"] = {"min_defect": df13, "time_s": round(dt13, 2)}
    if df13 == 0:
        w2_13 = solve_and_verify_thr2_mod2p(13, 5, bw13)
        out["thr2_mod10_k13_128bit"]["verified_2_to_13"] = True
        print("    >>> PERCEE HISTORIQUE : THR_{13, 2} MOD 10 RESOLU ET VERIFIE SUR 2^13 = 8 192 ENTREES ! <<<")

    print("\n==========================================================================")
    print("  2. CHASSE AUX NOUVEAUX RECORDS SUR MOD 21 (3 x 7), MOD 14 (2 x 7), MOD 15 (3 x 5)")
    print("==========================================================================")
    # (a) mod 21 (3 x 7) : 3-blocs et 4-blocs au-dela de n = 12 (2-blocs = 6+6 = 12)
    mod21_res = []
    for part in [(6, 6, 3), (6, 6, 5), (6, 6, 6), (6, 6, 6, 2), (6, 6, 6, 3), (6, 6, 6, 6)]:
        t1 = time.time()
        df, w1, w2, npts = search_kblock_odd_composite_thr2(3, 7, np.array(part, dtype=np.int64), 128, 35000)
        dt = time.time() - t1
        n = sum(part)
        print(f"  [THR_2 mod 21, {len(part)}-blocs {part}, n={n:2d}] min_defect = {df} ({dt:.2f}s)")
        entry = {"partition": list(part), "n": n, "min_defect": df, "w3": [int(x) for x in w1], "w7": [int(x) for x in w2], "time_s": round(dt, 2)}
        if df == 0:
            verify_kblock_odd_bruteforce(3, 7, part, w1, w2)
            entry["verified_bruteforce"] = True
            print(f"    >>> NOUVEAU RECORD MONDIAL : THR_{{{n}, 2}} MOD 21 RESOLU ({n/6:.2f}x le symetrique 6) ! <<<")
        mod21_res.append(entry)
    out["thr2_mod21_hunt"] = mod21_res

    # (b) mod 15 (3 x 5) : 4-blocs au-dela de n = 15 (3-blocs = 5+5+5 = 15)
    mod15_res = []
    for part in [(5, 5, 5, 1), (5, 5, 4, 2), (5, 5, 5, 2), (5, 5, 5, 3), (5, 5, 5, 5)]:
        t1 = time.time()
        df, w1, w2, npts = search_kblock_odd_composite_thr2(3, 5, np.array(part, dtype=np.int64), 128, 35000)
        dt = time.time() - t1
        n = sum(part)
        print(f"  [THR_2 mod 15, 4-blocs {part}, n={n:2d}] min_defect = {df} ({dt:.2f}s)")
        entry = {"partition": list(part), "n": n, "min_defect": df, "w3": [int(x) for x in w1], "w5": [int(x) for x in w2], "time_s": round(dt, 2)}
        if df == 0:
            verify_kblock_odd_bruteforce(3, 5, part, w1, w2)
            entry["verified_bruteforce"] = True
            print(f"    >>> NOUVEAU RECORD MONDIAL : THR_{{{n}, 2}} MOD 15 RESOLU ({n/5:.2f}x le symetrique 5) ! <<<")
        mod15_res.append(entry)
    out["thr2_mod15_hunt"] = mod15_res

    # (c) mod 14 (2 x 7) : 4-blocs au-dela de n = 15 (3-blocs = 2+7+6 = 15)
    mod14_res = []
    for part in [(5, 4, 4, 3), (6, 4, 3, 3), (4, 4, 4, 4), (5, 4, 4, 4), (6, 5, 4, 3), (5, 5, 4, 4)]:
        t1 = time.time()
        df, bw, npts = search_4block_thr2_mod2p(7, np.array(part, dtype=np.int64), 128, 30000)
        dt = time.time() - t1
        n = sum(part)
        print(f"  [THR_2 mod 14, 4-blocs {part}, n={n:2d}] min_defect = {df} ({dt:.2f}s)")
        entry = {"partition": list(part), "n": n, "min_defect": df, "w7": [int(x) for x in bw], "time_s": round(dt, 2)}
        if df == 0:
            block_of = []
            for bi, sz in enumerate(part):
                block_of.extend([bi] * sz)
            edges = [(u, v) for u in range(n) for v in range(u + 1, n)]
            inter_idx = {(0,1):4, (0,2):5, (0,3):6, (1,2):7, (1,3):8, (2,3):9}
            w7_full = np.array([bw[block_of[u]] if block_of[u] == block_of[v] else bw[inter_idx[(min(block_of[u], block_of[v]), max(block_of[u], block_of[v]))]] for (u, v) in edges], dtype=np.int8)
            if n <= 17:
                w2_ok = solve_and_verify_thr2_mod2p(n, 7, w7_full)
                entry["verified_2_to_n"] = True
            print(f"    >>> NOUVEAU RECORD MONDIAL : THR_{{{n}, 2}} MOD 14 RESOLU ({n/7:.2f}x le symetrique 7) ! <<<")
        mod14_res.append(entry)
    out["thr2_mod14_hunt"] = mod14_res

    out["total_time_s"] = round(time.time() - t0, 2)
    out_path = os.path.join(RESULTS_DIR, "thr2_mod10_128bit_and_multi_mod.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print(f"\nTermine en {time.time() - t0:.2f}s. Sauvegarde : {out_path}")


if __name__ == "__main__":
    main()
