"""
Attaque chirurgicale pour fermer THR_{n, 2} mod 6 en degre 2 (9 <= N <= 11)
===========================================================================
1. Enumeration 100% EXHAUSTIVE des familles symetriques sur K_10 :
   - Graphes circulants Z_10 (dim 5 -> 121 droites F_3 x Gauss F_2)
   - Graphes Z_5-invariants sur 5 + 5 sommets (dim 9 -> 9 841 droites F_3 x Gauss F_2,
     incluant Petersen et tous les bicirculants de K_10 !)
   - Graphes Z_3-invariants sur 3 + 3 + 3 + 1 sommets (recuit + echantillonnage)
   - 5-blocs (2, 2, 2, 2, 1) [n=9] et (2, 2, 2, 2, 2) [n=10] (dim 15)
2. Recuit 2-opt (1-flip + 2-flip + re-heating) massif sur K_10 general (45 aretes)
   pour tester si n = 10 existe sans symetrie.
3. Comptage canonique EXACT des graphes w_3 valides (f2_defect == 0) pour
   k = 2, 3, 4, 5, 6 afin d'evaluer la faisabilite d'un certificat d'impossibilite
   par branchement canonique (et tester directement la faisabilite de n = 11 via
   la classification de Gyori).
"""
import itertools
import json
import os
import sys
import time

import numpy as np
from numba import njit, prange

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from composite_mod_m_lab import RESULTS_DIR, decode_projective
from option_a_thr2_breakthrough import f2_defect_count
from piste_2a_thr2_graphs import solve_and_verify_thr2_mod2p


# =========================================================================== #
#  1. Enumeration 100% EXHAUSTIVE sur orbites d'aretes de K_10
# =========================================================================== #
@njit(parallel=True, cache=True)
def exact_orbit_enum_k10(orbit_of_edge, n_orbits):
    """
    Etant donne une partition des 45 aretes de K_10 en `n_orbits` orbites
    (ex. dim 5 pour Z_10, dim 9 pour Z_5 sur 5+5), enumere 100% EXHAUSTIVEMENT
    toutes les (3^n_orbits - 1)/2 droites P_3 sur F_3 et verifie la solvabilite
    EXACTE de P_2 sur F_2 (sur TOUTES les 45 aretes libres de K_10, donc encore
    plus general que P_2 invariant !).
    """
    n = 10
    n_edges = 45
    eu = np.zeros(n_edges, dtype=np.int32)
    ev = np.zeros(n_edges, dtype=np.int32)
    e = 0
    for u in range(n):
        for v in range(u + 1, n):
            eu[e] = u; ev[e] = v; e += 1

    # Tous les 1013 masques de poids >= 2
    masks = np.zeros(1024, dtype=np.int32)
    nm = 0
    for msk in range(1, 1024):
        pc = 0; tmp = msk
        while tmp:
            pc += tmp & 1; tmp >>= 1
        if pc >= 2:
            masks[nm] = msk; nm += 1

    edge_bits = np.zeros(nm, dtype=np.int64)
    # Matrice d'incidence masque x orbite (mod 3)
    M_orb = np.zeros((nm, n_orbits), dtype=np.int8)
    for i in range(nm):
        msk = masks[i]
        bm = np.int64(0)
        for ed in range(n_edges):
            if ((msk >> eu[ed]) & 1) and ((msk >> ev[ed]) & 1):
                bm |= (np.int64(1) << ed)
                M_orb[i, orbit_of_edge[ed]] = (M_orb[i, orbit_of_edge[ed]] + 1) % 3
        edge_bits[i] = bm

    n3 = (3 ** n_orbits - 1) // 2
    defs = np.full(n3, 999999, dtype=np.int32)
    for idx in prange(n3):
        cf = np.zeros(n_orbits, dtype=np.int64)
        decode_projective(idx, n_orbits, 3, cf)
        z_bits = np.zeros(nm, dtype=np.int64)
        nz = 0
        for i in range(nm):
            acc = 0
            for o in range(n_orbits):
                if cf[o] != 0:
                    acc += cf[o] * M_orb[i, o]
            if acc % 3 == 0:
                z_bits[nz] = edge_bits[i]
                nz += 1
        defs[idx] = f2_defect_count(z_bits, nz, n_edges)

    bi = int(np.argmin(defs))
    cf_b = np.zeros(n_orbits, dtype=np.int64)
    decode_projective(bi, n_orbits, 3, cf_b)
    return int(defs[bi]), cf_b, int(np.sum(defs == 0))


# =========================================================================== #
#  2. Recherche sur 5-blocs (2,2,2,2,1) [n=9] et (2,2,2,2,2) [n=10] (dim 15)
# =========================================================================== #
@njit(parallel=True, cache=True)
def search_5block_thr2_mod6(part, n_restarts, steps):
    """
    Base homogene degre 2 en 5 blocs : dim = 5 + 10 = 15.
    (A_i in 0..part[i], somme >= 2).
    """
    dim = 15
    S0, S1, S2, S3, S4 = part[0]+1, part[1]+1, part[2]+1, part[3]+1, part[4]+1
    max_pts = S0 * S1 * S2 * S3 * S4
    M3 = np.zeros((max_pts, dim), dtype=np.int8)
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
                            M3[n_pts, k] = vals[k] % 3
                            if vals[k] & 1:
                                bm |= (np.int64(1) << k)
                        M2_bits[n_pts] = bm
                        n_pts += 1

    bdefs = np.full(n_restarts, 999999, dtype=np.int32)
    bws = np.zeros((n_restarts, dim), dtype=np.int8)
    for r in prange(n_restarts):
        st = np.uint64(0x9E3779B97F4A7C15) ^ np.uint64((r + 1) * 2654435761 + n_pts * 101)
        w = np.zeros(dim, dtype=np.int8)
        vals = np.zeros(n_pts, dtype=np.int8)
        for k in range(dim):
            st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
            v0 = np.int8(st % np.uint64(3))
            w[k] = v0
            if v0 != 0:
                for p in range(n_pts):
                    vals[p] = (vals[p] + v0 * M3[p, k]) % 3
        zb = np.zeros(n_pts, dtype=np.int64)
        nz = 0
        for p in range(n_pts):
            if vals[p] == 0:
                zb[nz] = M2_bits[p]; nz += 1
        cdef = f2_defect_count(zb, nz, dim)
        bd = cdef
        bw = w.copy()
        T0 = 2.5; Tend = 0.05
        ratio = (Tend / T0) ** (1.0 / max(1, steps))
        temp = T0
        for _ in range(steps):
            if bd == 0:
                break
            st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
            k = int(st % np.uint64(dim))
            st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
            delta = 1 if (st & np.uint64(1)) else 2
            nz = 0
            for p in range(n_pts):
                nv = (vals[p] + delta * M3[p, k]) % 3
                if nv == 0:
                    zb[nz] = M2_bits[p]; nz += 1
            ndef = f2_defect_count(zb, nz, dim)
            diff = ndef - cdef
            acc = (diff <= 0)
            if not acc:
                st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
                if (float(st & np.uint64(0xFFFFFFFF)) / 4294967296.0) < np.exp(-float(diff) / temp):
                    acc = True
            if acc:
                w[k] = (w[k] + delta) % 3
                for p in range(n_pts):
                    vals[p] = (vals[p] + delta * M3[p, k]) % 3
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


# =========================================================================== #
#  3. Recuit 2-opt + re-heating massif sur K_10 general (45 aretes)
# =========================================================================== #
@njit(parallel=True, cache=True)
def search_k10_2opt_deep(n_restarts, steps):
    """
    Recuit simule avec mouvements 1-arete (65%) et 2-aretes simultanes (35%)
    + reheating cyclique sur K_10 (45 aretes, 1013 masques de poids >= 2).
    """
    n = 10
    n_edges = 45
    eu = np.zeros(n_edges, dtype=np.int32)
    ev = np.zeros(n_edges, dtype=np.int32)
    e = 0
    for u in range(n):
        for v in range(u + 1, n):
            eu[e] = u; ev[e] = v; e += 1

    masks = np.zeros(1024, dtype=np.int32)
    nm = 0
    for msk in range(1, 1024):
        pc = 0; tmp = msk
        while tmp:
            pc += tmp & 1; tmp >>= 1
        if pc >= 2:
            masks[nm] = msk; nm += 1

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
    cur = np.zeros(n_edges, dtype=np.int32)
    for i in range(nm):
        msk = masks[i]
        for ed in range(n_edges):
            if ((msk >> eu[ed]) & 1) and ((msk >> ev[ed]) & 1):
                inc_list[ed, cur[ed]] = i; cur[ed] += 1

    bdefs = np.full(n_restarts, 999999, dtype=np.int32)
    bws = np.zeros((n_restarts, n_edges), dtype=np.int8)

    for r in prange(n_restarts):
        st = np.uint64(0x9E3779B97F4A7C15) ^ np.uint64((r + 1) * 6364136223846793005 + 1010101)
        w = np.zeros(n_edges, dtype=np.int8)
        vals = np.zeros(nm, dtype=np.int8)
        for ed in range(n_edges):
            st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
            v0 = np.int8(st % np.uint64(3))
            w[ed] = v0
            if v0 != 0:
                for t in range(inc_cnt[ed]):
                    vals[inc_list[ed, t]] = (vals[inc_list[ed, t]] + v0) % 3

        zb = np.zeros(nm, dtype=np.int64)
        nz = 0
        for i in range(nm):
            if vals[i] == 0:
                zb[nz] = edge_bits[i]; nz += 1
        cdef = f2_defect_count(zb, nz, n_edges)
        bd = cdef
        bw = w.copy()

        cycle_len = max(1000, steps // 3)
        for step in range(steps):
            if bd == 0:
                break
            frac = float(step % cycle_len) / float(cycle_len)
            temp = 2.8 * ((0.06 / 2.8) ** frac)

            st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
            do_two = (st % np.uint64(100)) < np.uint64(35)

            st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
            e1 = np.int64(st % np.uint64(n_edges))
            st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
            d1 = 1 if (st & np.uint64(1)) else 2

            for t in range(inc_cnt[e1]):
                vals[inc_list[e1, t]] = (vals[inc_list[e1, t]] + d1) % 3

            e2 = np.int64(0)
            d2 = 0
            if do_two:
                st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
                e2 = np.int64(st % np.uint64(n_edges))
                if e2 == e1:
                    e2 = (e1 + 1) % n_edges
                st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
                d2 = 1 if (st & np.uint64(1)) else 2
                for t in range(inc_cnt[e2]):
                    vals[inc_list[e2, t]] = (vals[inc_list[e2, t]] + d2) % 3

            nz = 0
            for i in range(nm):
                if vals[i] == 0:
                    zb[nz] = edge_bits[i]; nz += 1
            ndef = f2_defect_count(zb, nz, n_edges)
            diff = ndef - cdef
            acc = (diff <= 0)
            if not acc:
                st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
                if (float(st & np.uint64(0xFFFFFFFF)) / 4294967296.0) < np.exp(-float(diff) / temp):
                    acc = True

            if acc:
                w[e1] = (w[e1] + d1) % 3
                if do_two:
                    w[e2] = (w[e2] + d2) % 3
                cdef = ndef
                if cdef < bd:
                    bd = cdef
                    for j in range(n_edges):
                        bw[j] = w[j]
            else:
                inv1 = 3 - d1
                for t in range(inc_cnt[e1]):
                    vals[inc_list[e1, t]] = (vals[inc_list[e1, t]] + inv1) % 3
                if do_two:
                    inv2 = 3 - d2
                    for t in range(inc_cnt[e2]):
                        vals[inc_list[e2, t]] = (vals[inc_list[e2, t]] + inv2) % 3

        bdefs[r] = bd
        for j in range(n_edges):
            bws[r, j] = bw[j]

    bi = int(np.argmin(bdefs))
    return int(bdefs[bi]), bws[bi]


# =========================================================================== #
#  4. Comptage canonique EXACT des graphes w_3 valides pour k = 2, 3, 4, 5, 6
# =========================================================================== #
def canonical_w3_key(n, w):
    """Forme canonique exacte d'un graphe pondere w in {0,1,2}^{C(n,2)} sous S_n x {+1,-1}."""
    edges = [(u, v) for u in range(n) for v in range(u + 1, n)]
    mat = [[0] * n for _ in range(n)]
    for idx, (u, v) in enumerate(edges):
        mat[u][v] = mat[v][u] = int(w[idx])
    best = None
    for perm in itertools.permutations(range(n)):
        seq1 = tuple(mat[perm[u]][perm[v]] for (u, v) in edges)
        seq2 = tuple((3 - x) % 3 for x in seq1)
        cand = seq1 if seq1 < seq2 else seq2
        if best is None or cand < best:
            best = cand
    return best


def count_canonical_valid_w3(max_k=5):
    """
    Enumere EXHAUSTIVEMENT tous les graphes w_3 valides (admettant un P_2 mod 2)
    a isomorphisme S_k x {+1,-1} pres, sommet par sommet de k = 2 jusqu'a max_k.
    """
    counts = {}
    # k = 2 : w in {(0,), (1,)} (car (2,) ~ (1,) par negation)
    cur_list = [(0,), (1,)]
    counts[2] = len(cur_list)

    for k in range(3, max_k + 1):
        t1 = time.time()
        n_edges = k * (k - 1) // 2
        edges = [(u, v) for u in range(k) for v in range(u + 1, k)]
        idx_prev = [i for i, (u, v) in enumerate(edges) if v < k - 1]
        idx_new = [i for i, (u, v) in enumerate(edges) if v == k - 1]

        # Precalcul des masques de poids >= 2 sur k sommets
        masks = [m for m in range(1, 1 << k) if m.bit_count() >= 2]
        edge_bits = np.zeros(len(masks), dtype=np.int64)
        inc_e = []
        for mi, m in enumerate(masks):
            bm = 0
            es = []
            for ei, (u, v) in enumerate(edges):
                if ((m >> u) & 1) and ((m >> v) & 1):
                    bm |= (1 << ei)
                    es.append(ei)
            edge_bits[mi] = bm
            inc_e.append(es)

        seen_canon = set()
        next_list = []
        for w_prev in cur_list:
            for ext in itertools.product((0, 1, 2), repeat=k - 1):
                w = [0] * n_edges
                for pos, ei in enumerate(idx_prev):
                    w[ei] = w_prev[pos]
                for pos, ei in enumerate(idx_new):
                    w[ei] = ext[pos]
                # Test exact f2_defect == 0
                zb = np.zeros(len(masks), dtype=np.int64)
                nz = 0
                for mi, es in enumerate(inc_e):
                    s = 0
                    for ei in es:
                        s += w[ei]
                    if s % 3 == 0:
                        zb[nz] = edge_bits[mi]
                        nz += 1
                if f2_defect_count(zb, nz, n_edges) == 0:
                    key = canonical_w3_key(k, w)
                    if key not in seen_canon:
                        seen_canon.add(key)
                        next_list.append(key)
        cur_list = next_list
        counts[k] = len(cur_list)
        print(f"  [Canonique EXACT] k = {k} : {len(cur_list):6d} graphes w_3 valides (mod S_k x +-1) ({time.time()-t1:.2f}s)")
    return counts


def main():
    t0 = time.time()
    out = {}

    print("==========================================================================")
    print("  1. FAMILLES SYMETRIQUES EXHAUSTIVES SUR K_10 (THR_{10,2} mod 6, d = 2)")
    print("==========================================================================")
    edges10 = [(u, v) for u in range(10) for v in range(u + 1, 10)]

    # (a) Circulants Z_10 : orbite = min((v-u)%10, (u-v)%10) - 1 in 0..4
    orb_z10 = np.array([min((v - u) % 10, (u - v) % 10) - 1 for (u, v) in edges10], dtype=np.int32)
    df, cf, nsol = exact_orbit_enum_k10(orb_z10, 5)
    print(f"  [Z_10 circulant EXHAUSTIF, 121 droites F_3] : min_defect = {df}, sols = {nsol}")
    out["k10_Z10_circulant"] = {"min_defect": df, "solutions": nsol}

    # (b) Z_5-invariants sur 5 + 5 sommets (A = 0..4, B = 5..9) :
    #     orbites : dist 1 dans A (0), dist 2 dans A (1), dist 1 dans B (2), dist 2 dans B (3),
    #     inter (v-5 - u) % 5 in {0,1,2,3,4} -> 5 orbites (4..8) => total 9 orbites !
    orb_z5 = []
    for (u, v) in edges10:
        if u < 5 and v < 5:
            d = min((v - u) % 5, (u - v) % 5)
            orb_z5.append(0 if d == 1 else 1)
        elif u >= 5 and v >= 5:
            d = min((v - u) % 5, (u - v) % 5)
            orb_z5.append(2 if d == 1 else 3)
        else:
            orb_z5.append(4 + ((v - 5 - u) % 5))
    orb_z5 = np.array(orb_z5, dtype=np.int32)
    t1 = time.time()
    df, cf, nsol = exact_orbit_enum_k10(orb_z5, 9)
    print(f"  [Z_5 bicirculant EXHAUSTIF (incl. Petersen), 9 841 droites F_3] : min_defect = {df}, sols = {nsol} ({time.time()-t1:.2f}s)")
    out["k10_Z5_bicirculant"] = {"min_defect": df, "solutions": nsol, "best_cf": [int(x) for x in cf]}
    if df == 0:
        w3_10 = np.array([cf[o] for o in orb_z5], dtype=np.int8)
        w2_10 = solve_and_verify_thr2_mod2p(10, 3, w3_10)
        out["k10_SOLVED"] = {"source": "Z5_bicirculant", "w2": w2_10, "w3": [int(x) for x in w3_10]}
        print("  >>> INCROYABLE : THR_{10, 2} MOD 6 RESOLU PAR BICIRCULANT Z_5 ET VERIFIE SUR 2^10 = 1024 ! <<<")

    print("\n==========================================================================")
    print("  2. 5-BLOCS SUR n = 9 ET n = 10 : (2,2,2,2,1), (2,2,2,2,2), (3,2,2,2,1)")
    print("==========================================================================")
    b5_res = []
    for part in [(2, 2, 2, 2, 1), (2, 2, 2, 2, 2), (3, 2, 2, 2, 1), (3, 3, 2, 1, 1)]:
        t1 = time.time()
        df, bw, npts = search_5block_thr2_mod6(np.array(part, dtype=np.int64), 128, 30000)
        dt = time.time() - t1
        b5_res.append({"partition": list(part), "n": sum(part), "grid_pts": npts, "min_defect": df, "w3_5b": [int(x) for x in bw], "time_s": round(dt, 2)})
        print(f"  5-blocs {part} [n={sum(part)}, {npts:3d} pts] : min_defect = {df} ({dt:.2f}s)")
    out["mod6_5block"] = b5_res

    print("\n==========================================================================")
    print("  3. RECUIT 2-OPT + RE-HEATING MASSIF SUR K_10 GENERAL (512 x 40k = 20,5 M pas)")
    print("==========================================================================")
    t1 = time.time()
    df10, bw10 = search_k10_2opt_deep(512, 40000)
    dt10 = time.time() - t1
    print(f"  K_10 2-opt massif (512 x 40 000 = 20,48 M pas) : min_defect = {df10} ({dt10:.2f}s)")
    out["k10_2opt_massive"] = {"restarts": 512, "steps": 40000, "min_defect": df10, "time_s": round(dt10, 2)}
    if df10 == 0:
        w2_10 = solve_and_verify_thr2_mod2p(10, 3, bw10)
        out["k10_SOLVED"] = {"source": "2opt_massive", "w2": w2_10, "w3": [int(x) for x in bw10]}
        print("  >>> INCROYABLE : THR_{10, 2} MOD 6 RESOLU PAR 2-OPT ET VERIFIE SUR 2^10 = 1024 ! <<<")

    print("\n==========================================================================")
    print("  4. COMPTAGE CANONIQUE EXACT DES GRAPHES w_3 VALIDES (k = 2..5)")
    print("==========================================================================")
    out["canonical_w3_counts"] = count_canonical_valid_w3(max_k=5)

    out["total_time_s"] = round(time.time() - t0, 2)
    out_path = os.path.join(RESULTS_DIR, "thr2_mod6_n10_attack.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print(f"\nTermine en {time.time() - t0:.2f}s. Sauvegarde : {out_path}")


if __name__ == "__main__":
    main()
