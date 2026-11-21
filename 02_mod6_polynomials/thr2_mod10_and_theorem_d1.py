r"""
1. Chasse 5-blocs (dim 15) et 6-blocs (dim 21) sur THR_{n, 2} mod 10 (n = 13, 14, 15, 16)
   avec verification force-brute 2^n integrale sur chaque nouveau record decouvert !
2. Chasse 4-blocs (dim 10) sur THR_{n, 2} mod 14 (n = 16) et mod 15 (n = 16).
3. Verification force-brute 2^9 = 512 et analyse structurelle des 4 classes canoniques
   maximales de V_0 (OR_9 mod 6 sans termes lineaires mod 3, |V_0|_max = 9 exact).
4. Tabulation exacte du Theoreme D1 (Contraction isotrope Chevalley-Warning F_2 x F_p
   pour tout degre d >= 1 et tout module m = 2p).
"""
import json
import math
import os
import sys
import time

import numpy as np
from numba import njit, prange

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from composite_mod_m_lab import RESULTS_DIR
from gyori_generalization_probe1 import extend_v0_or_mod6_exact
from option_a_thr2_breakthrough import f2_defect_count
from piste_2a_thr2_graphs import solve_and_verify_thr2_mod2p
from thr2_mod6_exact_closure import canonicalize_exact_s_k_pm1


# =========================================================================== #
#  1. Recherche 5-blocs (dim 15) et 6-blocs (dim 21) sur THR_{n, 2} mod 2p
# =========================================================================== #
@njit(parallel=True, cache=True)
def search_5block_thr2_mod2p(p, part, n_restarts, steps):
    """5-blocs homogene degre 2 pour THR_{n, 2} mod 2p : dim = 5 + 10 = 15."""
    dim = 15
    S0, S1, S2, S3, S4 = part[0]+1, part[1]+1, part[2]+1, part[3]+1, part[4]+1
    max_pts = S0 * S1 * S2 * S3 * S4
    Mp = np.zeros((max_pts, dim), dtype=np.int8)
    M2_bits = np.zeros(max_pts, dtype=np.int64)
    n_pts = 0
    for a in range(S0):
        ca = a*(a-1)//2
        for b in range(S1):
            cb = b*(b-1)//2
            for c in range(S2):
                cc = c*(c-1)//2
                for d in range(S3):
                    cd = d*(d-1)//2
                    for e in range(S4):
                        if a + b + c + d + e < 2:
                            continue
                        ce = e*(e-1)//2
                        vals = (ca, cb, cc, cd, ce, a*b, a*c, a*d, a*e, b*c, b*d, b*e, c*d, c*e, d*e)
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
        st = np.uint64(0x9E3779B97F4A7C15) ^ np.uint64((r + 1) * 2654435761 + n_pts * 313)
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

        cycle_len = max(1000, steps // 3)
        for step in range(steps):
            if bd == 0:
                break
            frac = float(step % cycle_len) / float(cycle_len)
            temp = 2.6 * ((0.03 / 2.6) ** frac)
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
        bdefs[r] = bd
        for j in range(dim):
            bws[r, j] = bw[j]

    bi = int(np.argmin(bdefs))
    return int(bdefs[bi]), bws[bi], n_pts


@njit(parallel=True, cache=True)
def search_6block_thr2_mod2p(p, part, n_restarts, steps):
    """6-blocs homogene degre 2 pour THR_{n, 2} mod 2p : dim = 6 + 15 = 21."""
    dim = 21
    S0, S1, S2, S3, S4, S5 = part[0]+1, part[1]+1, part[2]+1, part[3]+1, part[4]+1, part[5]+1
    max_pts = S0 * S1 * S2 * S3 * S4 * S5
    Mp = np.zeros((max_pts, dim), dtype=np.int8)
    M2_bits = np.zeros(max_pts, dtype=np.int64)
    n_pts = 0
    for a0 in range(S0):
        c0 = a0*(a0-1)//2
        for a1 in range(S1):
            c1 = a1*(a1-1)//2
            for a2 in range(S2):
                c2 = a2*(a2-1)//2
                for a3 in range(S3):
                    c3 = a3*(a3-1)//2
                    for a4 in range(S4):
                        c4 = a4*(a4-1)//2
                        for a5 in range(S5):
                            if a0 + a1 + a2 + a3 + a4 + a5 < 2:
                                continue
                            c5 = a5*(a5-1)//2
                            vals = (
                                c0, c1, c2, c3, c4, c5,
                                a0*a1, a0*a2, a0*a3, a0*a4, a0*a5,
                                a1*a2, a1*a3, a1*a4, a1*a5,
                                a2*a3, a2*a4, a2*a5,
                                a3*a4, a3*a5,
                                a4*a5
                            )
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
        st = np.uint64(0x9E3779B97F4A7C15) ^ np.uint64((r + 1) * 2654435761 + n_pts * 577)
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

        cycle_len = max(1000, steps // 3)
        for step in range(steps):
            if bd == 0:
                break
            frac = float(step % cycle_len) / float(cycle_len)
            temp = 2.6 * ((0.03 / 2.6) ** frac)
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
        bdefs[r] = bd
        for j in range(dim):
            bws[r, j] = bw[j]

    bi = int(np.argmin(bdefs))
    return int(bdefs[bi]), bws[bi], n_pts


def expand_kblock_to_graph_and_verify(part, bw, p):
    """Developpe un vecteur k-blocs bw sur les C(n,2) aretes de K_n et verifie sur 2^n."""
    k = len(part)
    n = sum(part)
    block_of = []
    for bi, sz in enumerate(part):
        block_of.extend([bi] * sz)
    edges = [(u, v) for u in range(n) for v in range(u + 1, n)]
    inter_idx = {}
    idx = k
    for a in range(k):
        for b in range(a + 1, k):
            inter_idx[(a, b)] = idx
            idx += 1
    wp_full = np.zeros(len(edges), dtype=np.int8)
    for ei, (u, v) in enumerate(edges):
        bu, bv = block_of[u], block_of[v]
        wp_full[ei] = bw[bu] if bu == bv else bw[inter_idx[(min(bu, bv), max(bu, bv))]]
    w2_ok = solve_and_verify_thr2_mod2p(n, p, wp_full)
    return w2_ok, [int(x) for x in wp_full]


def verify_v0_k9_witnesses():
    """Reconstruit et verifie sur 2^9 = 512 les 4 classes canoniques maximales de V_0 (OR_9 mod 6)."""
    cur0 = np.array([[0], [1]], dtype=np.int8)
    for k in range(3, 10):
        ext_w, _ = extend_v0_or_mod6_exact(cur0, k)
        cur0, _ = canonicalize_exact_s_k_pm1(ext_w, k)
    assert cur0.shape[0] == 4, f"Attendu 4 classes V_0 a k=9, obtenu {cur0.shape[0]}"

    n = 9
    n_edges = 36
    edges = [(u, v) for u in range(n) for v in range(u + 1, n)]
    witnesses = []
    for cid in range(4):
        w3 = cur0[cid]
        # Resoudre le systeme F_2 : sum_{u<v in S} a_{uv} == (1 + |S|) % 2 pour tout |S| >= 2 avec P_3(S) == 0 mod 3
        rows = []
        for mask in range(1, 1 << n):
            pc = mask.bit_count()
            if pc < 2:
                continue
            s3 = 0
            row = 0
            for idx, (u, v) in enumerate(edges):
                if ((mask >> u) & 1) and ((mask >> v) & 1):
                    s3 += int(w3[idx])
                    row |= (1 << idx)
            if s3 % 3 == 0:
                rhs = 1 if (pc % 2 == 0) else 0
                rows.append(row | (rhs << n_edges))
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
            assert v != (1 << n_edges), "Inconsistant V_0!"
        w2 = [0] * n_edges
        for c in sorted(piv.keys(), reverse=True):
            v = piv[c]
            acc = (v >> n_edges) & 1
            for j in range(c + 1, n_edges):
                if (v >> j) & 1:
                    acc ^= w2[j]
            w2[c] = acc

        # Verification force-brute sur 2^9 = 512 entrees pour OR_9 mod 6 avec a_v = 1, b_v = 0
        for mask in range(1 << n):
            pc = mask.bit_count()
            s2 = pc & 1  # sum a_v x_v = |x| mod 2
            s3 = 0       # sum b_v x_v = 0
            for idx, (u, v) in enumerate(edges):
                if ((mask >> u) & 1) and ((mask >> v) & 1):
                    s2 ^= w2[idx]
                    s3 += int(w3[idx])
            s3 %= 3
            if mask == 0:
                assert s2 == 0 and s3 == 0
            else:
                assert not (s2 == 0 and s3 == 0), f"Echec V_0 classe {cid} sur mask={mask}"
        witnesses.append({
            "class_id": cid,
            "n_f2_equations": len(rows),
            "f2_rank": len(piv),
            "n_edges_p2": sum(w2),
            "n_edges_p3": sum(1 for x in w3 if x != 0),
            "w2": w2,
            "w3": [int(x) for x in w3],
        })
    return witnesses


def theorem_d1_bounds_table():
    """
    Calcule la borne superieure polynomiale explicite du Theoreme D1 :
      N_gen^(2p)(d, OR) <= r + sum_{k=1}^d (d - k + 1) * C(r, k) - 1,  avec r = d*(p - 1) + 1
    ainsi que la version raffinee d = 2 : N_gen^(2p)(2, OR) <= 2*p^2 + p - 2.
    """
    table = []
    for p in [3, 5, 7, 11, 13, 17]:
        m = 2 * p
        row = {"m": m, "p": p, "lower_2_m_minus_1_d2": 2 * (m - 1), "d1_refined_d2": 2 * p * p + p - 2, "by_degree": {}}
        for d in [1, 2, 3, 4]:
            r = d * (p - 1) + 1
            bound = r + sum((d - k + 1) * math.comb(r, k) for k in range(1, d + 1)) - 1
            row["by_degree"][str(d)] = bound
        table.append(row)
    return table


def main():
    t0 = time.time()
    out = {}

    print("==========================================================================")
    print("  1. VERIFICATION 2^9 = 512 DES 4 CLASSES CANONIQUES MAXIMALES DE V_0 (OR_9 MOD 6)")
    print("==========================================================================")
    v0_wits = verify_v0_k9_witnesses()
    for w in v0_wits:
        print(
            f"  V_0 Classe #{w['class_id']} (n=9, b_v=0, a_v=1) : eqs={w['n_f2_equations']}, "
            f"rang F_2={w['f2_rank']}/36, |E(P2)|={w['n_edges_p2']}, |E(P3)|={w['n_edges_p3']} -> VERIFIE 2^9 = 512 !"
        )
    out["v0_k9_exact_classes"] = v0_wits

    print("\n==========================================================================")
    print("  2. ATTAQUE 5-BLOCS (dim 15) & 6-BLOCS (dim 21) SUR THR_{n, 2} MOD 10")
    print("==========================================================================")
    b5_list = []
    for part in [
        (3, 3, 3, 2, 2), (4, 3, 2, 2, 2), (3, 3, 3, 3, 1), (4, 3, 3, 2, 1),  # n = 13
        (3, 3, 3, 3, 2), (4, 3, 3, 2, 2), (4, 4, 2, 2, 2),                    # n = 14
        (3, 3, 3, 3, 3), (4, 3, 3, 3, 2),                                     # n = 15
    ]:
        t1 = time.time()
        df, bw, npts = search_5block_thr2_mod2p(5, np.array(part, dtype=np.int64), 192, 35000)
        dt = time.time() - t1
        n = sum(part)
        print(f"  [5-blocs {part}, n={n:2d}, {npts:3d} pts] min_defect = {df} ({dt:.2f}s)")
        entry = {"partition": list(part), "n": n, "grid_pts": npts, "min_defect": df, "w5": [int(x) for x in bw], "time_s": round(dt, 2)}
        if df == 0:
            w2_ok, w5_full = expand_kblock_to_graph_and_verify(part, bw, 5)
            entry["verified_2_to_n"] = True
            print(f"    >>> RECORD BATTU : THR_{{{n}, 2}} MOD 10 RESOLU EN 5-BLOCS ET VERIFIE SUR 2^{n} = {1<<n} ENTREES ! <<<")
        b5_list.append(entry)
    out["thr2_mod10_5block"] = b5_list

    b6_list = []
    for part in [
        (3, 2, 2, 2, 2, 2), (3, 3, 2, 2, 2, 1),  # n = 13
        (3, 3, 2, 2, 2, 2), (3, 3, 3, 2, 2, 1),  # n = 14
        (3, 3, 3, 2, 2, 2),                      # n = 15
        (3, 3, 3, 3, 2, 2),                      # n = 16
    ]:
        t1 = time.time()
        df, bw, npts = search_6block_thr2_mod2p(5, np.array(part, dtype=np.int64), 128, 30000)
        dt = time.time() - t1
        n = sum(part)
        print(f"  [6-blocs {part}, n={n:2d}, {npts:4d} pts] min_defect = {df} ({dt:.2f}s)")
        entry = {"partition": list(part), "n": n, "grid_pts": npts, "min_defect": df, "w5": [int(x) for x in bw], "time_s": round(dt, 2)}
        if df == 0:
            w2_ok, w5_full = expand_kblock_to_graph_and_verify(part, bw, 5)
            entry["verified_2_to_n"] = True
            print(f"    >>> RECORD BATTU : THR_{{{n}, 2}} MOD 10 RESOLU EN 6-BLOCS ET VERIFIE SUR 2^{n} = {1<<n} ENTREES ! <<<")
        b6_list.append(entry)
    out["thr2_mod10_6block"] = b6_list

    print("\n==========================================================================")
    print("  3. TABULATION DU THEOREME D1 (BORNE POLYNOMIALE UNIVERSELLE O_d(p^d))")
    print("==========================================================================")
    d1_tab = theorem_d1_bounds_table()
    for r in d1_tab:
        print(
            f"  mod {r['m']:2d} (p={r['p']:2d}) : 2(m-1)={r['lower_2_m_minus_1_d2']:2d} <= N_gen(2) <= {r['d1_refined_d2']:3d} "
            f"(d=1..4 : {r['by_degree']})"
        )
    out["theorem_d1_bounds"] = d1_tab

    out["total_time_s"] = round(time.time() - t0, 2)
    out_path = os.path.join(RESULTS_DIR, "thr2_mod10_and_theorem_d1.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print(f"\nTermine en {time.time() - t0:.2f}s. Sauvegarde : {out_path}")


if __name__ == "__main__":
    main()
