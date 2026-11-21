r"""
Etape 1 & 2 — Generalisation de Gyori (2006) et Theoreme D1 (Contraction isotrope)
==================================================================================
1. Sous-graphe V_0 = {v : b_v = 0 (mod p)} dans OR_n mod 2p (d = 2) :
   Pour tout v in V_0, P_p(e_v) = 0 impose P_2(e_v) = a_v = 1 (mod 2).
   Donc pour tout S subseteq V_0 (|S| >= 2) avec P_p(S) = 0 (mod p), on doit avoir :
       sum_{u < v in S} a_{uv} = 1 + |S| (mod 2).
   Pour p = 3 (mod 6) : nous enumerons 100% EXHAUSTIVEMENT (mod S_k x {+-1}) tous
   les sous-graphes V_0 valides pour k = 1, 2, 3, ... jusqu'a extinction EXACTE !
   Cela donne une borne superieure EXACTE |V_0| <= k_max sur le nombre de variables
   sans terme lineaire mod 3 dans TOUTE solution de OR_n mod 6 !

2. Meme classification canonique EXACTE pour le sous-graphe V_{12} = {v : b_v = 1 (mod 3)}
   (quitte a utiliser l'astuce projectivisee ou l'enumeration des sommets b_v in {1, 2})
   afin de classifier OR_n mod 6 de maniere 100% exacte !

3. Sonde sur THR_{n, 2} mod 10 (d = 2) :
   Actuellement 11 <= N_gen^(10)(2, THR_2) <= 19 (3-blocs = 11).
   Testons 4-blocs, 5-blocs et recuit 2-opt general sur K_12, K_13, K_14, K_15 mod 10
   pour voir si THR_{n, 2} mod 10 depasse n = 11 !
"""
import itertools
import json
import os
import sys
import time

import numpy as np
from numba import njit, prange

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from composite_mod_m_lab import RESULTS_DIR
from option_a_thr2_breakthrough import f2_defect_count
from piste_2a_thr2_graphs import solve_and_verify_thr2_mod2p
from thr2_mod6_exact_closure import canonicalize_exact_s_k_pm1


# =========================================================================== #
#  1. Classification EXACTE du sous-graphe V_0 (b_v = 0, a_v = 1) de OR_n mod 6
# =========================================================================== #
@njit(parallel=True, cache=True)
def extend_v0_or_mod6_exact(prev_list, k):
    """
    Dans V_0 (ou b_v = 0 et a_v = 1 pour tout v) :
    Pour tout sous-ensemble S subseteq {0..k-1} de taille |S| >= 2 tel que
    P_3(S) = sum_{u<v in S} b_{uv} == 0 (mod 3), l'equation sur F_2 est :
        sum_{u<v in S} a_{uv} == (1 + |S|) & 1  (mod 2).
    Attention : pour |S| pair (ex. |S| = 2, 4, 6...), le membre de droite vaut 1 !
                pour |S| impair (ex. |S| = 3, 5, 7...), le membre de droite vaut 0 !
    En particulier, pour |S| = 2 ({u, v}), P_3({u, v}) = b_{uv}, donc b_{uv} = 0
    impose a_{uv} = 1 ; et si un triangle {u, v, w} a b_{uv}=b_{vw}=b_{wu}=0,
    alors a_{uv}+a_{vw}+a_{wu} = 1+1+1 = 1 != 0 (mod 2) -> K_3 de 0-aretes INTERDIT !
    """
    n_prev = prev_list.shape[0]
    n_ext = 3 ** (k - 1)
    n_edges = k * (k - 1) // 2
    n_prev_edges = (k - 1) * (k - 2) // 2
    rhs_mask = np.int64(1) << n_edges

    eu = np.zeros(n_edges, dtype=np.int32)
    ev = np.zeros(n_edges, dtype=np.int32)
    idx_prev = np.zeros(n_prev_edges, dtype=np.int32)
    idx_new = np.zeros(k - 1, dtype=np.int32)
    e = 0; ip = 0; inew = 0
    for u in range(k):
        for v in range(u + 1, k):
            eu[e] = u; ev[e] = v
            if v < k - 1:
                idx_prev[ip] = e; ip += 1
            else:
                idx_new[inew] = e; inew += 1
            e += 1

    masks = np.zeros(1 << k, dtype=np.int32)
    rhs_of = np.zeros(1 << k, dtype=np.int64)
    nm = 0
    for msk in range(1, 1 << k):
        pc = 0; tmp = msk
        while tmp:
            pc += tmp & 1; tmp >>= 1
        if pc >= 2:
            masks[nm] = msk
            # RHS = (1 + |S|) % 2 : vaut 1 si |S| pair, 0 si |S| impair
            rhs_of[nm] = rhs_mask if ((pc & 1) == 0) else np.int64(0)
            nm += 1

    edge_bits = np.zeros(nm, dtype=np.int64)
    mask_ed_cnt = np.zeros(nm, dtype=np.int32)
    mask_ed = np.zeros((nm, n_edges), dtype=np.int32)
    for i in range(nm):
        msk = masks[i]
        bm = np.int64(0)
        cnt = 0
        for ed in range(n_edges):
            if ((msk >> eu[ed]) & 1) and ((msk >> ev[ed]) & 1):
                bm |= (np.int64(1) << ed)
                mask_ed[i, cnt] = ed
                cnt += 1
        edge_bits[i] = bm
        mask_ed_cnt[i] = cnt

    valid_counts = np.zeros(n_prev, dtype=np.int32)
    for pi in prange(n_prev):
        w = np.zeros(n_edges, dtype=np.int8)
        for j in range(n_prev_edges):
            w[idx_prev[j]] = prev_list[pi, j]
        piv = np.zeros(n_edges, dtype=np.int64)
        cnt_ok = 0
        for ec in range(n_ext):
            tmp = ec
            for j in range(k - 1):
                w[idx_new[j]] = tmp % 3
                tmp //= 3
            for c in range(n_edges):
                piv[c] = 0
            ok = True
            for i in range(nm):
                s = 0
                mc = mask_ed_cnt[i]
                for t in range(mc):
                    s += w[mask_ed[i, t]]
                if s % 3 == 0:
                    v = edge_bits[i] | rhs_of[i]
                    for c in range(n_edges):
                        if (v >> c) & 1:
                            if piv[c] != 0:
                                v ^= piv[c]
                            else:
                                piv[c] = v
                                v = 0
                                break
                    if v == rhs_mask:
                        ok = False
                        break
            if ok:
                cnt_ok += 1
        valid_counts[pi] = cnt_ok

    total_ok = int(np.sum(valid_counts))
    offsets = np.zeros(n_prev, dtype=np.int64)
    acc = 0
    for pi in range(n_prev):
        offsets[pi] = acc
        acc += valid_counts[pi]

    out_w = np.zeros((total_ok, n_edges), dtype=np.int8)
    if total_ok == 0:
        return out_w, valid_counts

    for pi in prange(n_prev):
        if valid_counts[pi] == 0:
            continue
        w = np.zeros(n_edges, dtype=np.int8)
        for j in range(n_prev_edges):
            w[idx_prev[j]] = prev_list[pi, j]
        piv = np.zeros(n_edges, dtype=np.int64)
        pos = offsets[pi]
        for ec in range(n_ext):
            tmp = ec
            for j in range(k - 1):
                w[idx_new[j]] = tmp % 3
                tmp //= 3
            for c in range(n_edges):
                piv[c] = 0
            ok = True
            for i in range(nm):
                s = 0
                mc = mask_ed_cnt[i]
                for t in range(mc):
                    s += w[mask_ed[i, t]]
                if s % 3 == 0:
                    v = edge_bits[i] | rhs_of[i]
                    for c in range(n_edges):
                        if (v >> c) & 1:
                            if piv[c] != 0:
                                v ^= piv[c]
                            else:
                                piv[c] = v
                                v = 0
                                break
                    if v == rhs_mask:
                        ok = False
                        break
            if ok:
                for j in range(n_edges):
                    out_w[pos, j] = w[j]
                pos += 1
    return out_w, valid_counts


# =========================================================================== #
#  2. Classification EXACTE de V_1 = {v : b_v = 1 (mod 3)} dans OR_n mod 6
#     (a_v in {0, 1} LIBRE sur F_2, a_{uv} in {0, 1} LIBRE sur F_2 !)
# =========================================================================== #
@njit(parallel=True, cache=True)
def extend_v1_or_mod6_exact(prev_list, k):
    """
    Dans V_1 (ou b_v = 1 mod 3 pour tout v in {0..k-1}) :
    - Pour tout S != empty, P_3(S) = (|S| + sum_{u<v in S} b_{uv}) mod 3.
    - Les variables F_2 sont : k variables lineaires a_v (colonnes 0..k-1)
      + C(k,2) variables quadratiques a_{uv} (colonnes k..k+C(k,2)-1).
    - Pour tout S != empty tel que P_3(S) == 0 (mod 3), on impose :
        sum_{v in S} a_v + sum_{u<v in S} a_{uv} == 1 (mod 2).
    Note : les permutations S_k agissent sur V_1 en preservant b_v = 1 (sans negation +-1,
    car -1 envoie b_v = 1 sur b_v = 2 !).
    """
    n_prev = prev_list.shape[0]
    n_ext = 3 ** (k - 1)
    n_edges = k * (k - 1) // 2
    n_cols = k + n_edges
    n_prev_edges = (k - 1) * (k - 2) // 2
    rhs_mask = np.int64(1) << n_cols

    eu = np.zeros(n_edges, dtype=np.int32)
    ev = np.zeros(n_edges, dtype=np.int32)
    idx_prev = np.zeros(n_prev_edges, dtype=np.int32)
    idx_new = np.zeros(k - 1, dtype=np.int32)
    e = 0; ip = 0; inew = 0
    for u in range(k):
        for v in range(u + 1, k):
            eu[e] = u; ev[e] = v
            if v < k - 1:
                idx_prev[ip] = e; ip += 1
            else:
                idx_new[inew] = e; inew += 1
            e += 1

    nm = (1 << k) - 1
    masks = np.zeros(nm, dtype=np.int32)
    pop_mod3 = np.zeros(nm, dtype=np.int8)
    row_bits = np.zeros(nm, dtype=np.int64)
    mask_ed_cnt = np.zeros(nm, dtype=np.int32)
    mask_ed = np.zeros((nm, n_edges), dtype=np.int32)

    for i in range(nm):
        msk = i + 1
        masks[i] = msk
        pc = 0; tmp = msk
        while tmp:
            pc += tmp & 1; tmp >>= 1
        pop_mod3[i] = pc % 3
        # Bits lineaires (0..k-1) + bits quadratiques (k..k+n_edges-1)
        bm = np.int64(msk)
        cnt = 0
        for ed in range(n_edges):
            if ((msk >> eu[ed]) & 1) and ((msk >> ev[ed]) & 1):
                bm |= (np.int64(1) << (k + ed))
                mask_ed[i, cnt] = ed
                cnt += 1
        row_bits[i] = bm
        mask_ed_cnt[i] = cnt

    valid_counts = np.zeros(n_prev, dtype=np.int32)
    for pi in prange(n_prev):
        w = np.zeros(n_edges, dtype=np.int8)
        for j in range(n_prev_edges):
            w[idx_prev[j]] = prev_list[pi, j]
        piv = np.zeros(n_cols, dtype=np.int64)
        cnt_ok = 0
        for ec in range(n_ext):
            tmp = ec
            for j in range(k - 1):
                w[idx_new[j]] = tmp % 3
                tmp //= 3
            for c in range(n_cols):
                piv[c] = 0
            ok = True
            for i in range(nm):
                s = int(pop_mod3[i])
                mc = mask_ed_cnt[i]
                for t in range(mc):
                    s += w[mask_ed[i, t]]
                if s % 3 == 0:
                    v = row_bits[i] | rhs_mask
                    for c in range(n_cols):
                        if (v >> c) & 1:
                            if piv[c] != 0:
                                v ^= piv[c]
                            else:
                                piv[c] = v
                                v = 0
                                break
                    if v == rhs_mask:
                        ok = False
                        break
            if ok:
                cnt_ok += 1
        valid_counts[pi] = cnt_ok

    total_ok = int(np.sum(valid_counts))
    offsets = np.zeros(n_prev, dtype=np.int64)
    acc = 0
    for pi in range(n_prev):
        offsets[pi] = acc
        acc += valid_counts[pi]

    out_w = np.zeros((total_ok, n_edges), dtype=np.int8)
    if total_ok == 0:
        return out_w, valid_counts

    for pi in prange(n_prev):
        if valid_counts[pi] == 0:
            continue
        w = np.zeros(n_edges, dtype=np.int8)
        for j in range(n_prev_edges):
            w[idx_prev[j]] = prev_list[pi, j]
        piv = np.zeros(n_cols, dtype=np.int64)
        pos = offsets[pi]
        for ec in range(n_ext):
            tmp = ec
            for j in range(k - 1):
                w[idx_new[j]] = tmp % 3
                tmp //= 3
            for c in range(n_cols):
                piv[c] = 0
            ok = True
            for i in range(nm):
                s = int(pop_mod3[i])
                mc = mask_ed_cnt[i]
                for t in range(mc):
                    s += w[mask_ed[i, t]]
                if s % 3 == 0:
                    v = row_bits[i] | rhs_mask
                    for c in range(n_cols):
                        if (v >> c) & 1:
                            if piv[c] != 0:
                                v ^= piv[c]
                            else:
                                piv[c] = v
                                v = 0
                                break
                    if v == rhs_mask:
                        ok = False
                        break
            if ok:
                for j in range(n_edges):
                    out_w[pos, j] = w[j]
                pos += 1
    return out_w, valid_counts


def canonicalize_exact_s_k_only(arr, k):
    """Deduplique EXACTEMENT `arr` modulo S_k (sans negation +-1, car b_v = 1 est fixe)."""
    from thr2_mod6_exact_closure import canonical_key_exact_single
    edges = [(u, v) for u in range(k) for v in range(u + 1, k)]
    seen = {}
    keep_idx = []
    for g in range(arr.shape[0]):
        w = arr[g]
        A = [[0] * k for _ in range(k)]
        for ei, (u, v) in enumerate(edges):
            val = int(w[ei])
            A[u][v] = A[v][u] = val
        ckey = canonical_key_exact_single(A, k)
        if ckey not in seen:
            seen[ckey] = g
            keep_idx.append(g)
    if not keep_idx:
        return np.zeros((0, arr.shape[1]), dtype=np.int8), []
    return arr[keep_idx].copy(), list(seen.keys())


# =========================================================================== #
#  3. Recherche 4-blocs et graphes ponderes sur THR_{n, 2} mod 10 (n = 12, 13, 14)
# =========================================================================== #
@njit(parallel=True, cache=True)
def search_4block_thr2_mod2p(p, part, n_restarts, steps):
    """
    4-blocs pour THR_{n, 2} mod 2p en degre 2 : dim = 4 + 6 = 10.
    """
    dim = 10
    S0, S1, S2, S3 = part[0] + 1, part[1] + 1, part[2] + 1, part[3] + 1
    max_pts = S0 * S1 * S2 * S3
    Mp = np.zeros((max_pts, dim), dtype=np.int8)
    M2_bits = np.zeros(max_pts, dtype=np.int64)
    n_pts = 0
    for a in range(S0):
        ca = a * (a - 1) // 2
        for b in range(S1):
            cb = b * (b - 1) // 2
            for c in range(S2):
                cc = c * (c - 1) // 2
                for d in range(S3):
                    if a + b + c + d < 2:
                        continue
                    cd = d * (d - 1) // 2
                    vals = (ca, cb, cc, cd, a * b, a * c, a * d, b * c, b * d, c * d)
                    bm = np.int64(0)
                    for k in range(dim):
                        Mp[n_pts, k] = vals[k] % p
                        if vals[k] & 1:
                            bm |= (np.int64(1) << k)
                    M2_bits[n_pts] = bm
                    n_pts += 1

    bdefs = np.full(n_restarts, 999999, dtype=np.int32)
    bws = np.zeros((n_restarts, dim), dtype=np.int8)
    for r in prange(n_restarts):
        st = np.uint64(0x9E3779B97F4A7C15) ^ np.uint64((r + 1) * 2654435761 + n_pts * 131)
        w = np.zeros(dim, dtype=np.int8)
        vals = np.zeros(n_pts, dtype=np.int8)
        for k in range(dim):
            st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
            v0 = np.int8(st % np.uint64(p))
            w[k] = v0
            if v0 != 0:
                for pt in range(n_pts):
                    vals[pt] = (vals[pt] + v0 * Mp[pt, k]) % p
        zb = np.zeros(n_pts, dtype=np.int64)
        nz = 0
        for pt in range(n_pts):
            if vals[pt] == 0:
                zb[nz] = M2_bits[pt]; nz += 1
        cdef = f2_defect_count(zb, nz, dim)
        bd = cdef
        bw = w.copy()
        T0 = 2.4; Tend = 0.04
        ratio = (Tend / T0) ** (1.0 / max(1, steps))
        temp = T0
        for _ in range(steps):
            if bd == 0:
                break
            st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
            k = np.int64(st % np.uint64(dim))
            st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
            delta = np.int8(1 + (st % np.uint64(p - 1)))
            nz = 0
            for pt in range(n_pts):
                nv = (vals[pt] + delta * Mp[pt, k]) % p
                if nv == 0:
                    zb[nz] = M2_bits[pt]; nz += 1
            ndef = f2_defect_count(zb, nz, dim)
            diff = ndef - cdef
            acc = (diff <= 0)
            if not acc:
                st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
                if (float(st & np.uint64(0xFFFFFFFF)) / 4294967296.0) < np.exp(-float(diff) / temp):
                    acc = True
            if acc:
                w[k] = (w[k] + delta) % p
                for pt in range(n_pts):
                    vals[pt] = (vals[pt] + delta * Mp[pt, k]) % p
                cdef = ndef
                if cdef < bd:
                    bd = cdef
                    for j in range(dim):
                        bw[j] = w[j]
            temp *= ratio
        bdefs[r] = bd
        for j in range(dim):
            bws[r, j] = bw[j]

    bi = int(np.argmin(bdefs))
    return int(bdefs[bi]), bws[bi], n_pts


def main():
    t0 = time.time()
    out = {}

    print("==========================================================================")
    print("  1. BORNE EXACTE SUR |V_0| = #{v : b_v = 0 mod 3} DANS OR_n MOD 6 (d = 2)")
    print("==========================================================================")
    cur0 = np.array([[0], [1]], dtype=np.int8)
    v0_levels = {"2": {"canon_count": 2}}
    print("  [V_0, k = 2] : 2 classes canoniques (mod S_2 x +-1)")
    for k in range(3, 12):
        t1 = time.time()
        ext_w, vcnt = extend_v0_or_mod6_exact(cur0, k)
        uniq_w, _ = canonicalize_exact_s_k_pm1(ext_w, k)
        dt = time.time() - t1
        print(
            f"  [V_0, k = {k:2d}] {cur0.shape[0]:5d} parents x 3^{k-1} -> {ext_w.shape[0]:6d} valides -> "
            f"{uniq_w.shape[0]:5d} canoniques (dead={int(np.sum(vcnt==0))}/{cur0.shape[0]}) ({dt:.2f}s)"
        )
        v0_levels[str(k)] = {
            "parents": int(cur0.shape[0]),
            "valid_extensions": int(ext_w.shape[0]),
            "canon_count": int(uniq_w.shape[0]),
            "time_s": round(dt, 2),
        }
        cur0 = uniq_w
        if uniq_w.shape[0] == 0:
            print(f"  >>> THEOREME EXACT : |V_0| <= {k - 1} DANS TOUTE SOLUTION DE OR_n MOD 6 ! <<<")
            out["exact_max_V0_mod6"] = k - 1
            break
    out["v0_levels_mod6"] = v0_levels

    print("\n==========================================================================")
    print("  2. CLASSIFICATION EXACTE DE V_1 = #{v : b_v = 1 mod 3} DANS OR_n MOD 6")
    print("==========================================================================")
    cur1 = np.array([[0], [1], [2]], dtype=np.int8)
    v1_levels = {"2": {"canon_count": 3}}
    print("  [V_1, k = 2] : 3 classes canoniques (mod S_2)")
    for k in range(3, 8):
        t1 = time.time()
        ext_w, vcnt = extend_v1_or_mod6_exact(cur1, k)
        uniq_w, _ = canonicalize_exact_s_k_only(ext_w, k)
        dt = time.time() - t1
        print(
            f"  [V_1, k = {k:2d}] {cur1.shape[0]:5d} parents x 3^{k-1} -> {ext_w.shape[0]:6d} valides -> "
            f"{uniq_w.shape[0]:5d} canoniques (dead={int(np.sum(vcnt==0))}/{cur1.shape[0]}) ({dt:.2f}s)"
        )
        v1_levels[str(k)] = {
            "parents": int(cur1.shape[0]),
            "valid_extensions": int(ext_w.shape[0]),
            "canon_count": int(uniq_w.shape[0]),
            "time_s": round(dt, 2),
        }
        cur1 = uniq_w
        if uniq_w.shape[0] == 0 or uniq_w.shape[0] > 15000:
            break
    out["v1_levels_mod6"] = v1_levels

    print("\n==========================================================================")
    print("  3. ATTAQUE SUR THR_{n, 2} MOD 10 EN 4-BLOCS ET GRAPHES (n = 12, 13, 14)")
    print("==========================================================================")
    mod10_4b = []
    for part in [
        (4, 4, 2, 2), (4, 3, 3, 2), (3, 3, 3, 3),  # n = 12
        (4, 4, 3, 2), (4, 3, 3, 3),                 # n = 13
        (4, 4, 4, 2), (4, 4, 3, 3),                 # n = 14
        (4, 4, 4, 3),                               # n = 15
    ]:
        t1 = time.time()
        df, bw, npts = search_4block_thr2_mod2p(5, np.array(part, dtype=np.int64), 128, 25000)
        dt = time.time() - t1
        print(f"  [THR_2 mod 10, 4-blocs {part}, n={sum(part):2d}] min_defect = {df} ({dt:.2f}s)")
        entry = {"partition": list(part), "n": sum(part), "min_defect": df, "w5_4b": [int(x) for x in bw], "time_s": round(dt, 2)}
        if df == 0:
            # Reconstruire w5 sur K_n et verifier sur 2^n !
            n = sum(part)
            block_of = []
            for bi, sz in enumerate(part):
                block_of.extend([bi] * sz)
            edges = [(u, v) for u in range(n) for v in range(u + 1, n)]
            inter_idx = {(0,1):4, (0,2):5, (0,3):6, (1,2):7, (1,3):8, (2,3):9}
            w5_full = np.zeros(len(edges), dtype=np.int8)
            for ei, (u, v) in enumerate(edges):
                bu, bv = block_of[u], block_of[v]
                w5_full[ei] = bw[bu] if bu == bv else bw[inter_idx[(min(bu, bv), max(bu, bv))]]
            w2_ok = solve_and_verify_thr2_mod2p(n, 5, w5_full)
            entry["verified_2_to_n"] = True
            entry["w2"] = w2_ok
            print(f"    >>> PERCEE THR_{{{n}, 2}} MOD 10 RESOLU EN 4-BLOCS ET VERIFIE SUR 2^{n} = {1<<n} ENTREES ! <<<")
        mod10_4b.append(entry)
    out["thr2_mod10_4block"] = mod10_4b

    out["total_time_s"] = round(time.time() - t0, 2)
    out_path = os.path.join(RESULTS_DIR, "gyori_generalization_probe1.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print(f"\nTermine en {time.time() - t0:.2f}s. Sauvegarde : {out_path}")


if __name__ == "__main__":
    main()
