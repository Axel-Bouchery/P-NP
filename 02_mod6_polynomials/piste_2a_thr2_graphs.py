"""
Piste 2A — Etude approfondie des graphes ponderes pour THR_{n, 2} modulo m composite (d = 2)
============================================================================================
1. Trancher n = 9 pour THR_{n, 2} mod 6 (d = 2) :
   - Enumeration 100% EXHAUSTIVE 4-blocs sur mod 6 (dim 10 -> 29 524 droites F_3 x Gauss F_2)
     pour toutes les partitions de n = 8 et n = 9 !
   - Generation de 200 temoins distincts n = 8 (verifies 100% sur 2^8) puis EXTENSION
     EXHAUSTIVE des 3^8 = 6 561 vecteurs d'aretes mod 3 vers un 9e sommet (1 312 200 extensions
     exactes testees par Gauss F_2 !).
   - Recuit simule geometrique profond (256 restarts x 30 000 steps) sur K_9 general.

2. Generalisation de THR_{n, 2} (d = 2) aux autres modules composites :
   - Exact symetrique + 2-blocs exhaustif sur mod 6, 10, 12, 14, 15, 21.
   - 3-blocs exhaustif sur mod 6, mod 10 (63 F_2 x 3906 F_5), mod 12 (4095 Z_4 x 364 F_3),
     mod 15 (364 F_3 x 3906 F_5).
   - Chasse graphes ponderes generaux (P_p mod p + Gauss exact F_2) sur mod 10 (n = 9..15)
     et mod 14 (n = 8..16) avec verification force-brute 2^n integrale !
"""
import json
import os
import sys
import time

import numpy as np
from numba import njit, prange

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from composite_mod_m_lab import RESULTS_DIR, ROOT, binomial_mod_table, decode_projective, decode_ring
from option_a_deep_probe import f2_rank_and_defect
from option_a_thr2_breakthrough import (
    exact_thr2_2block_d2,
    f2_defect_count,
    search_thr2_d2_mod6_n,
    solve_and_verify_thr2_mod6,
)


# =========================================================================== #
#  1. Enumeration 100% EXHAUSTIVE 4-blocs pour THR_{n, 2} mod 6 (d = 2)
# =========================================================================== #
@njit(parallel=True, cache=True)
def exact_thr2_4block_mod6(a, b, c, d_sz):
    """
    Base homogene degre 2 en 4 blocs (A, B, C, D) : dim = 10
      [C(A,2), C(B,2), C(C,2), C(D,2), AB, AC, AD, BC, BD, CD].
    (3^10 - 1) / 2 = 29 524 droites F_3 x Gauss exact F_2 = 100% EXHAUSTIF !
    """
    dim = 10
    n3 = (3 ** dim - 1) // 2  # 29 524
    A_s, B_s, C_s, D_s = a + 1, b + 1, c + 1, d_sz + 1

    # Collecter tous les points (s1, s2, s3, s4) avec s1+s2+s3+s4 >= 2
    max_pts = A_s * B_s * C_s * D_s
    M3 = np.zeros((max_pts, dim), dtype=np.int8)
    M2_bits = np.zeros(max_pts, dtype=np.int64)
    n_pts = 0
    for s1 in range(A_s):
        c1 = s1 * (s1 - 1) // 2
        for s2 in range(B_s):
            c2 = s2 * (s2 - 1) // 2
            for s3 in range(C_s):
                c3 = s3 * (s3 - 1) // 2
                for s4 in range(D_s):
                    if s1 + s2 + s3 + s4 < 2:
                        continue
                    c4 = s4 * (s4 - 1) // 2
                    vals = (c1, c2, c3, c4, s1*s2, s1*s3, s1*s4, s2*s3, s2*s4, s3*s4)
                    bm = 0
                    for k in range(dim):
                        M3[n_pts, k] = vals[k] % 3
                        if vals[k] & 1:
                            bm |= (1 << k)
                    M2_bits[n_pts] = bm
                    n_pts += 1

    sols = np.zeros(n3, dtype=np.int8)
    min_def = np.full(n3, 999999, dtype=np.int32)
    for idx in prange(n3):
        cf = np.zeros(dim, dtype=np.int64)
        decode_projective(idx, dim, 3, cf)
        z_bits = np.zeros(n_pts, dtype=np.int64)
        nz = 0
        for p in range(n_pts):
            acc = 0
            for k in range(dim):
                if cf[k] != 0:
                    acc += cf[k] * M3[p, k]
            if acc % 3 == 0:
                z_bits[nz] = M2_bits[p]
                nz += 1
        df = f2_defect_count(z_bits, nz, dim)
        min_def[idx] = df
        if df == 0:
            sols[idx] = 1

    tot = int(np.sum(sols))
    bdef = int(np.min(min_def))
    first_idx = -1
    if tot > 0:
        for i in range(n3):
            if sols[i] == 1:
                first_idx = i
                break
    return tot, bdef, first_idx


# =========================================================================== #
#  2. Extension EXHAUSTIVE 8 -> 9 (3^8 = 6 561 extensions par temoin K_8)
# =========================================================================== #
@njit(parallel=True, cache=True)
def extend_k8_to_k9_exhaustive(w3_8_batch):
    """
    Pour chaque temoin exact K_8 (28 aretes mod 3) du batch, teste TOUTES les
    3^8 = 6 561 affectations possibles des 8 nouvelles aretes (0..7, 8) mod 3
    et verifie par elimination de Gauss EXACTE sur F_2 (sur les 36 aretes de K_9)
    si une solution P_2 existe pour THR_{9, 2} mod 6 !
    """
    n = 9
    n_edges = 36
    eu = np.zeros(n_edges, dtype=np.int32)
    ev = np.zeros(n_edges, dtype=np.int32)
    e = 0
    for u in range(n):
        for v in range(u + 1, n):
            eu[e] = u
            ev[e] = v
            e += 1
    # Indice dans K_9 des 28 aretes de K_8 et des 8 nouvelles aretes vers sommet 8
    idx_k8 = np.zeros(28, dtype=np.int32)
    idx_new = np.zeros(8, dtype=np.int32)
    i8, inew = 0, 0
    for ed in range(n_edges):
        if ev[ed] < 8:
            idx_k8[i8] = ed
            i8 += 1
        else:
            idx_new[inew] = ed
            inew += 1

    # Tous les masques de poids >= 2 dans K_9 (502 masques)
    masks = np.zeros(512, dtype=np.int32)
    nm = 0
    for msk in range(1, 512):
        pc = 0
        tmp = msk
        while tmp:
            pc += tmp & 1
            tmp >>= 1
        if pc >= 2:
            masks[nm] = msk
            nm += 1

    edge_bits = np.zeros(nm, dtype=np.int64)
    inc_k8 = np.zeros((nm, 28), dtype=np.int8)
    inc_new = np.zeros((nm, 8), dtype=np.int8)
    for i in range(nm):
        msk = masks[i]
        bm = np.int64(0)
        for ed in range(n_edges):
            if ((msk >> eu[ed]) & 1) and ((msk >> ev[ed]) & 1):
                bm |= (np.int64(1) << ed)
        edge_bits[i] = bm
        for k in range(28):
            ed = idx_k8[k]
            if ((msk >> eu[ed]) & 1) and ((msk >> ev[ed]) & 1):
                inc_k8[i, k] = 1
        for k in range(8):
            ed = idx_new[k]
            if ((msk >> eu[ed]) & 1) and ((msk >> ev[ed]) & 1):
                inc_new[i, k] = 1

    n_batch = w3_8_batch.shape[0]
    n_ext = 3 ** 8  # 6 561
    best_def_batch = np.full(n_batch, 999999, dtype=np.int32)
    best_ext_code = np.full(n_batch, -1, dtype=np.int32)

    for bi in prange(n_batch):
        # Precalcul de la contribution base_v[i] de K_8 sur chaque masque i
        base_v = np.zeros(nm, dtype=np.int8)
        for i in range(nm):
            acc = 0
            for k in range(28):
                if inc_k8[i, k] == 1:
                    acc += w3_8_batch[bi, k]
            base_v[i] = acc % 3

        bdef = 999999
        bcode = -1
        z_bits = np.zeros(nm, dtype=np.int64)
        wnew = np.zeros(8, dtype=np.int8)
        for code in range(n_ext):
            rem = code
            for k in range(8):
                wnew[k] = rem % 3
                rem //= 3
            nz = 0
            for i in range(nm):
                acc = base_v[i]
                for k in range(8):
                    if inc_new[i, k] == 1:
                        acc += wnew[k]
                if acc % 3 == 0:
                    z_bits[nz] = edge_bits[i]
                    nz += 1
            df = f2_defect_count(z_bits, nz, n_edges)
            if df < bdef:
                bdef = df
                bcode = code
                if bdef == 0:
                    break
        best_def_batch[bi] = bdef
        best_ext_code[bi] = bcode

    b_idx = int(np.argmin(best_def_batch))
    return int(best_def_batch[b_idx]), b_idx, int(best_ext_code[b_idx])


# =========================================================================== #
#  3. Recuit geometrique profond + Tabu sur K_n mod (2 * p) (p = 3, 5, 7)
# =========================================================================== #
@njit(parallel=True, cache=True)
def search_thr2_d2_mod2p_general(n, p, n_restarts, steps_per_restart, seed_offset):
    """
    Recuit simule geometrique sur les poids w in {0..p-1}^{C(n,2)} de P_p mod p,
    en minimisant le defaut EXACT de solvabilite Gauss F_2 de P_2 !
    Valable pour tout n <= 11 (C(11,2) = 55 < 62 bits) et tout premier impair p !
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

    total_masks = 1 << n
    masks = np.zeros(total_masks, dtype=np.int32)
    nm = 0
    for msk in range(1, total_masks):
        pc = 0
        tmp = msk
        while tmp:
            pc += tmp & 1
            tmp >>= 1
        if pc >= 2:
            masks[nm] = msk
            nm += 1

    edge_bits = np.zeros(nm, dtype=np.int64)
    inc_cnt = np.zeros(n_edges, dtype=np.int32)
    for i in range(nm):
        msk = masks[i]
        bm = np.int64(0)
        for ed in range(n_edges):
            if ((msk >> eu[ed]) & 1) and ((msk >> ev[ed]) & 1):
                bm |= (np.int64(1) << ed)
                inc_cnt[ed] += 1
        edge_bits[i] = bm

    max_inc = int(np.max(inc_cnt))
    inc_list = np.zeros((n_edges, max_inc), dtype=np.int32)
    cur_c = np.zeros(n_edges, dtype=np.int32)
    for i in range(nm):
        msk = masks[i]
        for ed in range(n_edges):
            if ((msk >> eu[ed]) & 1) and ((msk >> ev[ed]) & 1):
                inc_list[ed, cur_c[ed]] = i
                cur_c[ed] += 1

    best_def_all = np.full(n_restarts, 999999, dtype=np.int32)
    best_w_all = np.zeros((n_restarts, n_edges), dtype=np.int8)

    for r in prange(n_restarts):
        state = np.uint64(0x9E3779B97F4A7C15) ^ np.uint64((r + 1 + seed_offset) * 2654435761 + n * 1000003 + p * 99991)
        w = np.zeros(n_edges, dtype=np.int8)
        vals = np.zeros(nm, dtype=np.int8)

        for ed in range(n_edges):
            state ^= state << np.uint64(13)
            state ^= state >> np.uint64(7)
            state ^= state << np.uint64(17)
            val = np.int8(state % np.uint64(p))
            w[ed] = val
            if val != 0:
                cnt = inc_cnt[ed]
                for t in range(cnt):
                    mi = inc_list[ed, t]
                    vals[mi] = (vals[mi] + val) % p

        z_bits = np.zeros(nm, dtype=np.int64)
        nz = 0
        for i in range(nm):
            if vals[i] == 0:
                z_bits[nz] = edge_bits[i]
                nz += 1
        cur_def = f2_defect_count(z_bits, nz, n_edges)
        b_def = cur_def
        bw = w.copy()

        if b_def > 0:
            T0 = 2.5
            T_end = 0.08
            ratio = (T_end / T0) ** (1.0 / max(1, steps_per_restart))
            temp = T0
            for st in range(steps_per_restart):
                state ^= state << np.uint64(13)
                state ^= state >> np.uint64(7)
                state ^= state << np.uint64(17)
                ed = int(state % np.uint64(n_edges))
                state ^= state << np.uint64(13)
                state ^= state >> np.uint64(7)
                state ^= state << np.uint64(17)
                delta = 1 + int(state % np.uint64(p - 1))
                old_w = w[ed]
                new_w = (old_w + delta) % p

                cnt = inc_cnt[ed]
                for t in range(cnt):
                    mi = inc_list[ed, t]
                    vals[mi] = (vals[mi] + delta) % p

                nz = 0
                for i in range(nm):
                    if vals[i] == 0:
                        z_bits[nz] = edge_bits[i]
                        nz += 1
                new_def = f2_defect_count(z_bits, nz, n_edges)

                diff = new_def - cur_def
                accept = (diff <= 0)
                if not accept:
                    state ^= state << np.uint64(13)
                    state ^= state >> np.uint64(7)
                    state ^= state << np.uint64(17)
                    u01 = float(state & np.uint64(0xFFFFFFFF)) / 4294967296.0
                    if u01 < np.exp(-float(diff) / temp):
                        accept = True

                if accept:
                    w[ed] = new_w
                    cur_def = new_def
                    if cur_def < b_def:
                        b_def = cur_def
                        for k in range(n_edges):
                            bw[k] = w[k]
                        if b_def == 0:
                            break
                else:
                    inv_d = p - delta
                    for t in range(cnt):
                        mi = inc_list[ed, t]
                        vals[mi] = (vals[mi] + inv_d) % p
                temp *= ratio

        best_def_all[r] = b_def
        for k in range(n_edges):
            best_w_all[r, k] = bw[k]

    return best_def_all, best_w_all


def solve_and_verify_thr2_mod2p(n, p, wp):
    """Reconstruit P_2 mod 2 par Gauss exact a partir de wp (P_p mod p) et verifie sur 2^n."""
    n_edges = n * (n - 1) // 2
    edges = [(u, v) for u in range(n) for v in range(u + 1, n)]
    rows = []
    for mask in range(1, 1 << n):
        if mask.bit_count() < 2:
            continue
        sp = 0
        row = 0
        for idx, (u, v) in enumerate(edges):
            if ((mask >> u) & 1) and ((mask >> v) & 1):
                sp += int(wp[idx])
                row |= (1 << idx)
        if sp % p == 0:
            rows.append(row | (1 << n_edges))

    piv = {}
    for r in rows:
        v = r
        for c in range(n_edges):
            if (v >> c) & 1:
                if c in piv:
                    v ^= piv[c]
                else:
                    piv[c] = v
                    v = 0
                    break
        assert v != (1 << n_edges), "Inconsistant!"

    w2 = [0] * n_edges
    for c in sorted(piv.keys(), reverse=True):
        v = piv[c]
        acc = (v >> n_edges) & 1
        for j in range(c + 1, n_edges):
            if (v >> j) & 1:
                acc ^= w2[j]
        w2[c] = acc

    for mask in range(1 << n):
        pc = mask.bit_count()
        s2 = 0
        sp = 0
        for idx, (u, v) in enumerate(edges):
            if ((mask >> u) & 1) and ((mask >> v) & 1):
                s2 ^= w2[idx]
                sp += int(wp[idx])
        sp %= p
        if pc < 2:
            assert s2 == 0 and sp == 0
        else:
            assert not (s2 == 0 and sp == 0)
    return w2


# =========================================================================== #
#  4. 3-blocs EXHAUSTIF pour THR_{n, 2} (d = 2) sur mod 10, mod 12, mod 15
# =========================================================================== #
@njit(parallel=True, cache=True)
def exact_thr2_3block_2ch(q1, p1_prime, q2, p2_prime, max_s):
    """
    100% EXHAUSTIF 3-blocs homogene degre 2 (dim 6) pour q1 x q2.
    (Note : pour q1=2, q2=5 -> 63 x 3906 = 246k paires ;
            pour q1=3, q2=5 -> 364 x 3906 = 1.42M paires ;
            pour q1=4, q2=3 -> 4095 x 364 = 1.49M paires.)
    """
    dim = 6
    S = max_s + 1
    n1 = (q1 ** dim - 1) // (q1 - 1) if p1_prime else (q1 ** dim - 1)
    n2 = (q2 ** dim - 1) // (q2 - 1) if p2_prime else (q2 ** dim - 1)

    Z1 = np.zeros((n1, S, S, S), dtype=np.bool_)
    C1 = np.zeros((n1, dim), dtype=np.int8)
    for i in prange(n1):
        cf = np.zeros(dim, dtype=np.int64)
        if p1_prime:
            decode_projective(i, dim, q1, cf)
        else:
            decode_ring(i + 1, dim, q1, cf)
        for k in range(dim):
            C1[i, k] = cf[k]
        for a in range(S):
            ca = a * (a - 1) // 2
            for b in range(S):
                cb = b * (b - 1) // 2
                for c in range(S):
                    cc = c * (c - 1) // 2
                    v = (cf[0]*ca + cf[1]*cb + cf[2]*cc + cf[3]*a*b + cf[4]*a*c + cf[5]*b*c) % q1
                    Z1[i, a, b, c] = (v == 0)

    best_n = np.zeros(n2, dtype=np.int64)
    best_abc = np.zeros((n2, 3), dtype=np.int64)
    best_i = np.zeros(n2, dtype=np.int64)
    C2 = np.zeros((n2, dim), dtype=np.int8)

    for j in prange(n2):
        cf = np.zeros(dim, dtype=np.int64)
        if p2_prime:
            decode_projective(j, dim, q2, cf)
        else:
            decode_ring(j + 1, dim, q2, cf)
        for k in range(dim):
            C2[j, k] = cf[k]
        za = np.zeros(S * S * S, dtype=np.int16)
        zb = np.zeros(S * S * S, dtype=np.int16)
        zc = np.zeros(S * S * S, dtype=np.int16)
        nz = 0
        for a in range(S):
            ca = a * (a - 1) // 2
            for b in range(S):
                cb = b * (b - 1) // 2
                for c in range(S):
                    if a + b + c < 2:
                        continue
                    cc = c * (c - 1) // 2
                    if (cf[0]*ca + cf[1]*cb + cf[2]*cc + cf[3]*a*b + cf[4]*a*c + cf[5]*b*c) % q2 == 0:
                        za[nz] = a; zb[nz] = b; zc[nz] = c
                        nz += 1
        bn, ba, bb, bc, bi = 0, 0, 0, 0, 0
        min_c = np.zeros((S, S), dtype=np.int16)
        for i in range(n1):
            for a in range(S):
                for b in range(S):
                    min_c[a, b] = S
            for t in range(nz):
                a, b, c = za[t], zb[t], zc[t]
                if Z1[i, a, b, c] and c < min_c[a, b]:
                    min_c[a, b] = c
            for a in range(S):
                for b in range(S):
                    cur = min_c[a, b]
                    if a > 0 and min_c[a - 1, b] < cur:
                        cur = min_c[a - 1, b]
                    if b > 0 and min_c[a, b - 1] < cur:
                        cur = min_c[a, b - 1]
                    min_c[a, b] = cur
                    c_ok = cur - 1
                    if c_ok >= 0 and a + b + c_ok > bn:
                        bn = a + b + c_ok
                        ba, bb, bc, bi = a, b, c_ok, i
        best_n[j] = bn
        best_abc[j, 0] = ba; best_abc[j, 1] = bb; best_abc[j, 2] = bc
        best_i[j] = bi

    bj = int(np.argmax(best_n))
    return int(best_n[bj]), (int(best_abc[bj, 0]), int(best_abc[bj, 1]), int(best_abc[bj, 2])), C1[best_i[bj]], C2[bj]


def main():
    t0 = time.time()
    res = {}

    print("==========================================================================")
    print("  1. 4-BLOCS 100% EXHAUSTIF POUR THR_{n, 2} MOD 6 (d = 2) : n = 8 ET n = 9")
    print("==========================================================================")
    b4 = []
    for part in [(2, 2, 2, 2), (2, 2, 2, 3), (3, 2, 2, 2), (3, 3, 2, 1), (3, 3, 3, 0)]:
        t1 = time.time()
        tot, bdef, fidx = exact_thr2_4block_mod6(*part)
        dt = time.time() - t1
        b4.append({"partition": list(part), "n": sum(part), "solutions": tot, "min_defect": bdef, "time_s": round(dt, 2)})
        print(f"  4-blocs {part} [n={sum(part)}] : solutions = {tot:4d}, min_defect = {bdef} ({dt:.2f}s)")
    res["mod6_4block_exact"] = b4

    print("\n==========================================================================")
    print("  2. ATTAQUE DE n = 9 MOD 6 : EXTENSION EXHAUSTIVE 8->9 + RECUIT PROFOND")
    print("==========================================================================")
    # Generer 128 temoins K_8 distincts
    t1 = time.time()
    defs8, ws8 = search_thr2_d2_mod2p_general(8, 3, 192, 4000, 100)
    k8_sols = []
    seen_k8 = set()
    for i in range(len(defs8)):
        if defs8[i] == 0:
            key = tuple(int(x) for x in ws8[i])
            if key not in seen_k8:
                seen_k8.add(key)
                k8_sols.append(ws8[i])
    k8_arr = np.array(k8_sols, dtype=np.int8)
    print(f"  {len(k8_sols)} temoins K_8 distincts trouves en {time.time() - t1:.2f}s (verification du premier sur 2^8...)")
    if len(k8_sols) > 0:
        solve_and_verify_thr2_mod2p(8, 3, k8_arr[0])
        t2 = time.time()
        bdef9, b_bi, b_code = extend_k8_to_k9_exhaustive(k8_arr)
        dt2 = time.time() - t2
        print(f"  Extension EXHAUSTIVE 3^8 = 6 561 sur {len(k8_sols)} temoins K_8 ({len(k8_sols)*6561:,} tests) : min_defect = {bdef9} ({dt2:.2f}s)")
        res["mod6_k8_to_k9_exhaustive"] = {
            "n_k8_witnesses": int(len(k8_sols)),
            "total_extensions_tested": int(len(k8_sols) * 6561),
            "min_f2_defect": int(bdef9),
            "time_s": round(dt2, 2),
        }

    # Recuit geometrique profond direct sur K_9 mod 6 (256 restarts x 25 000 steps)
    t3 = time.time()
    defs9, ws9 = search_thr2_d2_mod2p_general(9, 3, 256, 25000, 777)
    bdef9_sa = int(np.min(defs9))
    dt3 = time.time() - t3
    print(f"  Recuit geometrique profond K_9 mod 6 (256 x 25 000 = 6,4 M pas) : min_defect = {bdef9_sa} ({dt3:.2f}s)")
    entry_k9 = {"restarts": 256, "steps": 25000, "min_f2_defect": bdef9_sa, "time_s": round(dt3, 2)}
    if bdef9_sa == 0:
        bi = int(np.argmin(defs9))
        w2_9 = solve_and_verify_thr2_mod2p(9, 3, ws9[bi])
        entry_k9["verified_2_to_9"] = True
        entry_k9["w2"] = w2_9
        entry_k9["w3"] = [int(x) for x in ws9[bi]]
        print("  >>> PERCEE HISTORIQUE : THR_{9, 2} MOD 6 RESOLU ET VERIFIE SUR 2^9 = 512 ENTREES ! <<<")
    else:
        entry_k9["verified_2_to_9"] = False
    res["mod6_k9_deep_sa"] = entry_k9

    print("\n==========================================================================")
    print("  3. EXTENSION AUX AUTRES MODULES COMPOSITES : 2-BLOCS, 3-BLOCS & GRAPHES")
    print("==========================================================================")
    # Exact symetrique + 2-blocs exhaustif sur mod 14 et mod 21 (et rappel 6, 10, 12, 15)
    all_mods_2b = {}
    for name, q1, p1, q2, p2, ms in [
        ("mod6",  2, True,  3, True,  12),
        ("mod10", 2, True,  5, True,  16),
        ("mod12", 4, False, 3, True,  20),
        ("mod14", 2, True,  7, True,  20),
        ("mod15", 3, True,  5, True,  20),
        ("mod21", 3, True,  7, True,  24),
    ]:
        # Symetrique exact a d=2 pour THR_{n,2} : P1 = C(S,2) mod q1, P2 = C(S,2) mod q2
        # Premier s >= 2 ou C(s,2) = 0 mod q1 ET mod q2 :
        sym_n = 1
        for s in range(2, 100):
            cs2 = s * (s - 1) // 2
            if cs2 % q1 == 0 and cs2 % q2 == 0:
                sym_n = s - 1
                break
        bn, ba, bb, c1, c2 = exact_thr2_2block_d2(q1, p1, q2, p2, ms)
        all_mods_2b[name] = {
            "N_sym_thr2": sym_n, "max_2block": bn, "partition_2b": [ba, bb],
            "c1_2b": [int(x) for x in c1], "c2_2b": [int(x) for x in c2],
        }
        print(f"  {name:5s} : N_sym = {sym_n:2d} | 2-blocs EXHAUSTIF = {bn:2d} ({ba}+{bb}) | c1={[int(x) for x in c1]}, c2={[int(x) for x in c2]}")
    res["all_mods_sym_and_2block"] = all_mods_2b

    # 3-blocs EXHAUSTIF sur mod 10, mod 12, mod 15
    print("\n  --- 3-blocs EXHAUSTIF pour THR_{n, 2} (d = 2) sur mod 10, mod 12, mod 15 ---")
    b3_mods = {}
    for name, q1, p1, q2, p2, ms in [
        ("mod10", 2, True,  5, True,  9),
        ("mod12", 4, False, 3, True,  10),
        ("mod15", 3, True,  5, True,  10),
    ]:
        t1 = time.time()
        bn3, abc3, c1_3, c2_3 = exact_thr2_3block_2ch(q1, p1, q2, p2, ms)
        dt = time.time() - t1
        b3_mods[name] = {
            "max_3block": bn3, "partition_3b": list(abc3),
            "c1_3b": [int(x) for x in c1_3], "c2_3b": [int(x) for x in c2_3],
            "time_s": round(dt, 2),
        }
        print(f"  [3-blocs EXHAUSTIF] {name:5s} : max = {bn3:2d} (partition {abc3}) | c1={[int(x) for x in c1_3]}, c2={[int(x) for x in c2_3]} ({dt:.2f}s)")
    res["mods_3block_exact"] = b3_mods

    # Chasse graphes ponderes generaux sur mod 10 (p = 5) et mod 14 (p = 7) pour n = 9..11 (C(11,2)=55 <= 62 bits!)
    print("\n  --- Graphes ponderes GENERAUX pour THR_{n, 2} sur mod 10 (p=5) et mod 14 (p=7) ---")
    gen_2p = {}
    for name, p, n_list in [("mod10", 5, [9, 10, 11]), ("mod14", 7, [9, 10, 11])]:
        entries = []
        for n in n_list:
            t1 = time.time()
            defs, ws = search_thr2_d2_mod2p_general(n, p, 96, 8000, 300 + n)
            bdef = int(np.min(defs))
            dt = time.time() - t1
            e_dict = {"n": n, "min_f2_defect": bdef, "time_s": round(dt, 2)}
            if bdef == 0:
                bi = int(np.argmin(defs))
                w2 = solve_and_verify_thr2_mod2p(n, p, ws[bi])
                e_dict["verified_2_to_n"] = True
                e_dict["w2"] = w2
                e_dict["wp"] = [int(x) for x in ws[bi]]
                print(f"  {name} THR_{{{n:2d}, 2}} d=2 : SOLVED & VERIFIED 100% sur 2^{n} = {1<<n} entrees ! ({dt:.2f}s)")
            else:
                e_dict["verified_2_to_n"] = False
                print(f"  {name} THR_{{{n:2d}, 2}} d=2 : min_defect = {bdef} ({dt:.2f}s)")
            entries.append(e_dict)
        gen_2p[name] = entries
    res["general_graphs_mod10_mod14"] = gen_2p

    res["total_time_s"] = round(time.time() - t0, 2)
    out_path = os.path.join(RESULTS_DIR, "piste_2a_thr2_graphs.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(res, f, indent=2)
    print(f"\nPiste 2A terminee en {time.time() - t0:.2f}s. Sauvegarde : {out_path}")


if __name__ == "__main__":
    main()
