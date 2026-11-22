r"""
Probe 2 sur THR_{10, 2} mod 6 :
1. Evaluation 45-aretes libres des meilleurs 5-blocs w_3 trouves sur K_10
   + optimisation 5-blocs w_3 avec P_2 LIBRE (45 aretes)
   + warm-start 2-opt general K_10 depuis ces candidats (et depuis K_9 -> K_10).
2. Preuve EXACTE de omega(G_2) <= 6 (backtracking exhaustif sur K_7 : aucun w_3
   ne satisfait P_3(T) != 0 mod 3 pour tous |T| in {4, 5}, donc aucune 7-clique dans G_2).
3. Enumeration canonique EXACTE des graphes w_3 valides a k = 6 (en Numba rapide)
   et estimation du branchement k = 6 -> 7 pour evaluer un certificat d'impossibilite exact.
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


# =========================================================================== #
#  1. Lifting 5-blocs w_3 -> P_2 LIBRE sur toutes les C(n,2) aretes
# =========================================================================== #
def build_block_orbit_map(part):
    n = sum(part)
    block_of = []
    for bi, sz in enumerate(part):
        block_of.extend([bi] * sz)
    edges = [(u, v) for u in range(n) for v in range(u + 1, n)]
    # 5 intra (0..4) + 10 inter (5..14)
    inter_idx = {}
    idx = 5
    for a in range(5):
        for b in range(a + 1, 5):
            inter_idx[(a, b)] = idx
            idx += 1
    orb = np.zeros(len(edges), dtype=np.int32)
    for ei, (u, v) in enumerate(edges):
        bu, bv = block_of[u], block_of[v]
        if bu == bv:
            orb[ei] = bu
        else:
            orb[ei] = inter_idx[(min(bu, bv), max(bu, bv))]
    return n, edges, orb


def eval_free_f2_defect(n, w3):
    n_edges = n * (n - 1) // 2
    edges = [(u, v) for u in range(n) for v in range(u + 1, n)]
    masks = [m for m in range(1, 1 << n) if m.bit_count() >= 2]
    zb = np.zeros(len(masks), dtype=np.int64)
    nz = 0
    for m in masks:
        s = 0
        bm = 0
        for ei, (u, v) in enumerate(edges):
            if ((m >> u) & 1) and ((m >> v) & 1):
                s += int(w3[ei])
                bm |= (1 << ei)
        if s % 3 == 0:
            zb[nz] = bm
            nz += 1
    return int(f2_defect_count(zb, nz, n_edges)), nz


@njit(parallel=True, cache=True)
def search_5block_w3_free_p2(orb, n, n_restarts, steps):
    """
    Optimise w_3 sur 15 parametres (5-blocs), mais evalue f2_defect avec
    P_2 COMPLETEMENT LIBRE sur les C(n,2) aretes !
    """
    n_edges = n * (n - 1) // 2
    dim = 15
    eu = np.zeros(n_edges, dtype=np.int32)
    ev = np.zeros(n_edges, dtype=np.int32)
    e = 0
    for u in range(n):
        for v in range(u + 1, n):
            eu[e] = u; ev[e] = v; e += 1

    masks = np.zeros(1 << n, dtype=np.int32)
    nm = 0
    for msk in range(1, 1 << n):
        pc = 0; tmp = msk
        while tmp:
            pc += tmp & 1; tmp >>= 1
        if pc >= 2:
            masks[nm] = msk; nm += 1

    edge_bits = np.zeros(nm, dtype=np.int64)
    M_orb = np.zeros((nm, dim), dtype=np.int8)
    for i in range(nm):
        msk = masks[i]
        bm = np.int64(0)
        for ed in range(n_edges):
            if ((msk >> eu[ed]) & 1) and ((msk >> ev[ed]) & 1):
                bm |= (np.int64(1) << ed)
                M_orb[i, orb[ed]] = (M_orb[i, orb[ed]] + 1) % 3
        edge_bits[i] = bm

    bdefs = np.full(n_restarts, 999999, dtype=np.int32)
    bws = np.zeros((n_restarts, dim), dtype=np.int8)

    for r in prange(n_restarts):
        st = np.uint64(0x9E3779B97F4A7C15) ^ np.uint64((r + 1) * 2654435761 + n * 7919)
        w = np.zeros(dim, dtype=np.int8)
        vals = np.zeros(nm, dtype=np.int8)
        for k in range(dim):
            st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
            v0 = np.int8(st % np.uint64(3))
            w[k] = v0
            if v0 != 0:
                for p in range(nm):
                    vals[p] = (vals[p] + v0 * M_orb[p, k]) % 3
        zb = np.zeros(nm, dtype=np.int64)
        nz = 0
        for p in range(nm):
            if vals[p] == 0:
                zb[nz] = edge_bits[p]; nz += 1
        cdef = f2_defect_count(zb, nz, n_edges)
        bd = cdef
        bw = w.copy()
        T0 = 2.5; Tend = 0.04
        ratio = (Tend / T0) ** (1.0 / max(1, steps))
        temp = T0
        for _ in range(steps):
            if bd == 0:
                break
            st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
            k = np.int64(st % np.uint64(dim))
            st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
            delta = 1 if (st & np.uint64(1)) else 2
            nz = 0
            for p in range(nm):
                nv = (vals[p] + delta * M_orb[p, k]) % 3
                if nv == 0:
                    zb[nz] = edge_bits[p]; nz += 1
            ndef = f2_defect_count(zb, nz, n_edges)
            diff = ndef - cdef
            acc = (diff <= 0)
            if not acc:
                st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
                if (float(st & np.uint64(0xFFFFFFFF)) / 4294967296.0) < np.exp(-float(diff) / temp):
                    acc = True
            if acc:
                w[k] = (w[k] + delta) % 3
                for p in range(nm):
                    vals[p] = (vals[p] + delta * M_orb[p, k]) % 3
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
    return int(bdefs[bi]), bws[bi]


@njit(parallel=True, cache=True)
def warmstart_k10_2opt(seeds, steps):
    """
    Warm-start 2-opt SA sur K_10 (45 aretes libres) depuis un ensemble de graines w_3.
    """
    n = 10
    n_edges = 45
    n_seeds = seeds.shape[0]
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

    bdefs = np.full(n_seeds, 999999, dtype=np.int32)
    bws = np.zeros((n_seeds, n_edges), dtype=np.int8)

    for r in prange(n_seeds):
        st = np.uint64(0x9E3779B97F4A7C15) ^ np.uint64((r + 1) * 6364136223846793005 + 424242)
        w = seeds[r].copy()
        vals = np.zeros(nm, dtype=np.int8)
        for ed in range(n_edges):
            v0 = w[ed]
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
            temp = 1.6 * ((0.03 / 1.6) ** frac)

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
#  2. PREUVE EXACTE DE omega(G_2) <= 6 (Backtracking exhaustif sur K_7)
# =========================================================================== #
@njit(cache=True)
def exact_prove_omega_g2_le_6():
    """
    Teste EXHAUSTIVEMENT s'il existe w_3 in {0,1,2}^{21} sur K_7 tel que
    P_3(T) != 0 (mod 3) pour TOUS les 56 sous-ensembles T de K_7 de taille 4 ou 5.
    Ordre des aretes : sommet par sommet (k = 2, 3, 4, 5, 6, 7).
    Pour chaque nouveau sommet v = 3..6, des que toutes les aretes vers {0..v-1}
    sont fixees, on elague immediatement sur tous les sous-ensembles T de taille 4 ou 5
    contenant v et contenus dans {0..v}.
    De plus, on exploite la symetrie :
    - Par notre lemme algebrique C_2(U) = 0 mod 3, on n'en a meme pas besoin pour la preuve :
      on fait un backtracking 100% SANS HYPOTHESE, uniquement avec brisure de symetrie
      exacte sur les 6 aretes de K_4 (les 27+ classes d'isomorphisme de K_4 valide !).
    En fait, les 6 aretes de K_4 ont seulement 3^6 = 729 configs, dont 486 ont somme != 0 mod 3 !
    Avec negation (+-1), seulement 243 -> elagage instantane !
    """
    # Aretes ordonnees par max vertex v :
    # v=1: (0,1) [0]
    # v=2: (0,2),(1,2) [1,2]
    # v=3: (0,3),(1,3),(2,3) [3,4,5] -> check T={0,1,2,3}
    # v=4: (0,4)..(3,4) [6..9] -> check 4 subsets of size 4 + 1 subset of size 5
    # v=5: (0,5)..(4,5) [10..14] -> check 10 subsets of size 4 + 5 subsets of size 5
    # v=6: (0,6)..(5,6) [15..20] -> check 20 subsets of size 4 + 15 subsets of size 5
    edge_id = np.full((7, 7), -1, dtype=np.int32)
    idx = 0
    for v in range(1, 7):
        for u in range(v):
            edge_id[u, v] = idx
            edge_id[v, u] = idx
            idx += 1

    # Precalcul des contraintes qui s'activent exactement a l'arete `e`
    # Une contrainte = masque de taille 4 ou 5, activee quand sa derniere arete (max_e) est assignee.
    c_masks = np.zeros(56, dtype=np.int32)
    c_trigger = np.zeros(56, dtype=np.int32)
    c_edges = np.zeros((56, 10), dtype=np.int32)
    c_ned = np.zeros(56, dtype=np.int32)
    nc = 0
    for m in range(1, 128):
        pc = 0; tmp = m
        while tmp:
            pc += tmp & 1; tmp >>= 1
        if pc == 4 or pc == 5:
            c_masks[nc] = m
            me = -1
            cnt = 0
            for u in range(7):
                if (m >> u) & 1:
                    for v in range(u + 1, 7):
                        if (m >> v) & 1:
                            eid = edge_id[u, v]
                            c_edges[nc, cnt] = eid
                            cnt += 1
                            if eid > me:
                                me = eid
            c_ned[nc] = cnt
            c_trigger[nc] = me
            nc += 1

    # Grouper les contraintes par arete declencheuse
    trig_cnt = np.zeros(21, dtype=np.int32)
    trig_list = np.zeros((21, 56), dtype=np.int32)
    for c in range(nc):
        t = c_trigger[c]
        trig_list[t, trig_cnt[t]] = c
        trig_cnt[t] += 1

    # Backtracking iteratif avec brisure de symetrie exacte :
    # On peut imposer w[0] in {0, 1} par negation (+-1) (si w[0]==0, la negation garde 0,
    # sinon w[0]==2 s'envoie sur 1). Donc 100% EXHAUSTIF !
    w = np.zeros(21, dtype=np.int8)
    choice = np.zeros(21, dtype=np.int8)
    depth = 0
    nodes = np.int64(0)
    sols = np.int64(0)
    sol_w = np.zeros(21, dtype=np.int8)

    while depth >= 0:
        max_v = 2 if depth == 0 else 3
        if choice[depth] >= max_v:
            choice[depth] = 0
            depth -= 1
            continue
        w[depth] = choice[depth]
        choice[depth] += 1
        nodes += 1

        # Verifier toutes les contraintes declenchees a `depth`
        ok = True
        tc = trig_cnt[depth]
        for k in range(tc):
            cid = trig_list[depth, k]
            s = 0
            ned = c_ned[cid]
            for j in range(ned):
                s += w[c_edges[cid, j]]
            if s % 3 == 0:
                ok = False
                break
        if not ok:
            continue

        if depth == 20:
            sols += 1
            for j in range(21):
                sol_w[j] = w[j]
            # Continuer pour avoir le compte exact !
        else:
            depth += 1
            choice[depth] = 0

    return int(sols), int(nodes), sol_w


# =========================================================================== #
#  3. Comptage canonique EXACT des graphes w_3 valides a k = 6 (et sonde k = 7)
# =========================================================================== #
@njit(parallel=True, cache=True)
def extend_valid_w3_batch(prev_list, k):
    """
    Etend une liste de graphes w_3 valides sur k-1 sommets a k sommets
    (3^(k-1) extensions chacun) en filtrant EXACTEMENT par f2_defect == 0.
    Retourne les indices (p_idx, ext_code) valides.
    """
    n_prev = prev_list.shape[0]
    n_ext = 3 ** (k - 1)
    n_edges = k * (k - 1) // 2
    n_prev_edges = (k - 1) * (k - 2) // 2

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
    nm = 0
    for msk in range(1, 1 << k):
        pc = 0; tmp = msk
        while tmp:
            pc += tmp & 1; tmp >>= 1
        if pc >= 2:
            masks[nm] = msk; nm += 1

    edge_bits = np.zeros(nm, dtype=np.int64)
    # Matrice d'incidence masque x arete
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
    # On fait un premier passage pour compter, puis on alloue exact
    for pi in prange(n_prev):
        w = np.zeros(n_edges, dtype=np.int8)
        for j in range(n_prev_edges):
            w[idx_prev[j]] = prev_list[pi, j]
        zb = np.zeros(nm, dtype=np.int64)
        cnt_ok = 0
        for ec in range(n_ext):
            tmp = ec
            for j in range(k - 1):
                w[idx_new[j]] = tmp % 3
                tmp //= 3
            nz = 0
            for i in range(nm):
                s = 0
                mc = mask_ed_cnt[i]
                for t in range(mc):
                    s += w[mask_ed[i, t]]
                if s % 3 == 0:
                    zb[nz] = edge_bits[i]
                    nz += 1
            if f2_defect_count(zb, nz, n_edges) == 0:
                cnt_ok += 1
        valid_counts[pi] = cnt_ok

    total_ok = int(np.sum(valid_counts))
    offsets = np.zeros(n_prev, dtype=np.int64)
    acc = 0
    for pi in range(n_prev):
        offsets[pi] = acc
        acc += valid_counts[pi]

    out_w = np.zeros((total_ok, n_edges), dtype=np.int8)
    for pi in prange(n_prev):
        if valid_counts[pi] == 0:
            continue
        w = np.zeros(n_edges, dtype=np.int8)
        for j in range(n_prev_edges):
            w[idx_prev[j]] = prev_list[pi, j]
        zb = np.zeros(nm, dtype=np.int64)
        pos = offsets[pi]
        for ec in range(n_ext):
            tmp = ec
            for j in range(k - 1):
                w[idx_new[j]] = tmp % 3
                tmp //= 3
            nz = 0
            for i in range(nm):
                s = 0
                mc = mask_ed_cnt[i]
                for t in range(mc):
                    s += w[mask_ed[i, t]]
                if s % 3 == 0:
                    zb[nz] = edge_bits[i]
                    nz += 1
            if f2_defect_count(zb, nz, n_edges) == 0:
                for j in range(n_edges):
                    out_w[pos, j] = w[j]
                pos += 1
    return out_w, valid_counts


@njit(cache=True)
def canonicalize_batch_k(arr, k, perms):
    """
    Calcule la forme canonique exacte (entier int64 en base 3, <= 3^21) pour chaque
    graphe de `arr` sous S_k x {+1, -1}.
    """
    n_graphs = arr.shape[0]
    n_edges = k * (k - 1) // 2
    n_perms = perms.shape[0]
    # Table (u, v) -> edge_idx
    eidx = np.zeros((k, k), dtype=np.int32)
    e = 0
    eu = np.zeros(n_edges, dtype=np.int32)
    ev = np.zeros(n_edges, dtype=np.int32)
    for u in range(k):
        for v in range(u + 1, k):
            eidx[u, v] = e; eidx[v, u] = e
            eu[e] = u; ev[e] = v
            e += 1
    # Precalculer pour chaque permutation le mapping edge_out -> edge_in
    perm_emap = np.zeros((n_perms, n_edges), dtype=np.int32)
    for p in range(n_perms):
        for ed in range(n_edges):
            perm_emap[p, ed] = eidx[perms[p, eu[ed]], perms[p, ev[ed]]]

    keys = np.zeros(n_graphs, dtype=np.int64)
    for g in range(n_graphs):
        best = np.int64(9223372036854775807)
        for p in range(n_perms):
            v1 = np.int64(0)
            v2 = np.int64(0)
            for ed in range(n_edges):
                x = np.int64(arr[g, perm_emap[p, ed]])
                v1 = v1 * 3 + x
                v2 = v2 * 3 + ((3 - x) % 3)
            c = v1 if v1 < v2 else v2
            if c < best:
                best = c
        keys[g] = best
    return keys


def main():
    t0 = time.time()
    out = {}

    print("==========================================================================")
    print("  1. PREUVE EXACTE : omega(G_2) <= 6 (backtracking exhaustif sur K_7)")
    print("==========================================================================")
    t1 = time.time()
    sols_k7, nodes_k7, sol_w7 = exact_prove_omega_g2_le_6()
    dt7 = time.time() - t1
    print(f"  [K_7 clique dans G_2] solutions EXACTES = {sols_k7}, noeuds explores = {nodes_k7:,} ({dt7:.3f}s)")
    out["exact_omega_g2_le_6"] = {
        "k7_clique_solutions": sols_k7,
        "nodes_explored": nodes_k7,
        "time_s": round(dt7, 3),
        "proven_omega_le_6": (sols_k7 == 0),
    }

    print("\n==========================================================================")
    print("  2. LIFTING 5-BLOCS w_3 -> P_2 LIBRE (45 ARETES) + WARM-START K_10")
    print("==========================================================================")
    atk_path = os.path.join(RESULTS_DIR, "thr2_mod6_n10_attack.json")
    with open(atk_path, "r", encoding="utf-8") as f:
        atk = json.load(f)

    seed_list = []
    lifted_evals = []
    for entry in atk["mod6_5block"]:
        part = tuple(entry["partition"])
        w5 = entry["w3_5b"]
        n, edges, orb = build_block_orbit_map(part)
        w3_full = np.array([w5[o] for o in orb], dtype=np.int8)
        df_free, nz = eval_free_f2_defect(n, w3_full)
        print(f"  Lift 5-blocs {part} [n={n}] : 15-dim def = {entry['min_defect']} -> P_2 LIBRE ({len(edges)} aretes) def = {df_free} (nz={nz})")
        lifted_evals.append({"partition": list(part), "n": n, "def_15d": entry["min_defect"], "def_free_p2": df_free, "nz": nz})
        if n == 10:
            seed_list.append(w3_full)
        if df_free == 0:
            w2_ok = solve_and_verify_thr2_mod2p(n, 3, w3_full)
            print(f"    >>> RESOLU SUR n = {n} PAR LIFTING 5-BLOCS ! VERIFIE SUR 2^{n} ! <<<")
            out[f"solved_lift_n{n}"] = {"partition": list(part), "w2": w2_ok, "w3": [int(x) for x in w3_full]}
    out["lifted_5block_evals"] = lifted_evals

    print("\n  -- Optimisation 5-blocs w_3 directement sous P_2 LIBRE (45 aretes) --")
    free_5b_res = []
    for part in [(2, 2, 2, 2, 2), (3, 2, 2, 2, 1), (4, 2, 2, 1, 1), (3, 3, 2, 1, 1)]:
        t1 = time.time()
        n, edges, orb = build_block_orbit_map(part)
        df, bw5 = search_5block_w3_free_p2(orb, n, 64, 10000)
        dt = time.time() - t1
        print(f"  5-blocs w_3 + P_2 LIBRE {part} : min_defect = {df} ({dt:.2f}s)")
        w3_full = np.array([bw5[o] for o in orb], dtype=np.int8)
        seed_list.append(w3_full)
        free_5b_res.append({"partition": list(part), "min_defect": df, "w3_5b": [int(x) for x in bw5], "time_s": round(dt, 2)})
        if df == 0:
            w2_ok = solve_and_verify_thr2_mod2p(10, 3, w3_full)
            out["k10_SOLVED"] = {"source": f"5block_free_p2_{part}", "w2": w2_ok, "w3": [int(x) for x in w3_full]}
            print("    >>> THR_{10, 2} MOD 6 RESOLU ET VERIFIE SUR 2^10 = 1024 ! <<<")
    out["k10_5block_w3_free_p2"] = free_5b_res

    # Warm-start 2-opt sur K_10 general depuis nos meilleures graines (x32 copies perturbees chacune)
    if "k10_SOLVED" not in out:
        rng = np.random.default_rng(20261120)
        warm_seeds = []
        for s in seed_list:
            warm_seeds.append(s.copy())
            for _ in range(31):
                c = s.copy()
                # perturber 1 a 4 aretes
                n_mut = int(rng.integers(1, 5))
                idxs = rng.choice(45, size=n_mut, replace=False)
                for ii in idxs:
                    c[ii] = (c[ii] + int(rng.integers(1, 3))) % 3
                warm_seeds.append(c)
        warm_arr = np.array(warm_seeds, dtype=np.int8)
        t1 = time.time()
        df_ws, bw_ws = warmstart_k10_2opt(warm_arr, 25000)
        dt_ws = time.time() - t1
        print(f"  Warm-start 2-opt K_10 ({warm_arr.shape[0]} graines x 25k pas) : min_defect = {df_ws} ({dt_ws:.2f}s)")
        out["k10_warmstart_2opt"] = {"n_seeds": int(warm_arr.shape[0]), "steps": 25000, "min_defect": df_ws, "time_s": round(dt_ws, 2)}
        if df_ws == 0:
            w2_ok = solve_and_verify_thr2_mod2p(10, 3, bw_ws)
            out["k10_SOLVED"] = {"source": "warmstart_2opt", "w2": w2_ok, "w3": [int(x) for x in bw_ws]}
            print("    >>> THR_{10, 2} MOD 6 RESOLU PAR WARM-START 2-OPT ET VERIFIE SUR 2^10 ! <<<")

    print("\n==========================================================================")
    print("  3. ENUMERATION CANONIQUE EXACTE k = 2..6 + SONDE BRANCHEMENT k = 6 -> 7")
    print("==========================================================================")
    cur = np.array([[0], [1]], dtype=np.int8)
    canon_stats = {"2": {"canon_count": 2}}
    for k in [3, 4, 5, 6]:
        t1 = time.time()
        ext_w, vcnt = extend_valid_w3_batch(cur, k)
        perms = np.array(list(itertools.permutations(range(k))), dtype=np.int32)
        keys = canonicalize_batch_k(ext_w, k, perms)
        uniq_keys, first_idx = np.unique(keys, return_index=True)
        cur = ext_w[first_idx].copy()
        dt = time.time() - t1
        avg_ext = float(np.mean(vcnt))
        max_ext = int(np.max(vcnt))
        zero_ext = int(np.sum(vcnt == 0))
        print(
            f"  [EXACT] k = {k} : {ext_w.shape[0]:7d} ext. valides -> {cur.shape[0]:6d} canoniques "
            f"(avg_ext={avg_ext:.2f}/{3**(k-1)}, max={max_ext}, dead_ends={zero_ext}/{len(vcnt)}) ({dt:.2f}s)"
        )
        canon_stats[str(k)] = {
            "valid_extensions": int(ext_w.shape[0]),
            "canon_count": int(cur.shape[0]),
            "avg_valid_ext": round(avg_ext, 3),
            "max_valid_ext": max_ext,
            "dead_ends": zero_ext,
            "parent_count": int(len(vcnt)),
            "time_s": round(dt, 2),
        }

    # Sonde statistique sur 200 graphes canoniques de k = 6 vers k = 7 (3^6 = 729 extensions chacun)
    rng = np.random.default_rng(777)
    sample_sz = min(256, cur.shape[0])
    sample_idx = rng.choice(cur.shape[0], size=sample_sz, replace=False)
    t1 = time.time()
    ext7, vcnt7 = extend_valid_w3_batch(cur[sample_idx], 7)
    dt7s = time.time() - t1
    print(
        f"  [SONDE k=6->7, {sample_sz} parents x 729] : {ext7.shape[0]:6d} ext. valides "
        f"(avg={np.mean(vcnt7):.2f}/729, max={np.max(vcnt7)}, dead_ends={int(np.sum(vcnt7==0))}/{sample_sz}) ({dt7s:.2f}s)"
    )
    canon_stats["probe_6_to_7"] = {
        "sample_parents": sample_sz,
        "valid_extensions": int(ext7.shape[0]),
        "avg_valid_ext": round(float(np.mean(vcnt7)), 3),
        "max_valid_ext": int(np.max(vcnt7)),
        "dead_ends": int(np.sum(vcnt7 == 0)),
        "est_total_ext_k7": int(round(float(np.mean(vcnt7)) * cur.shape[0])),
        "time_s": round(dt7s, 2),
    }
    out["canonical_exact_k2_to_6"] = canon_stats

    out["total_time_s"] = round(time.time() - t0, 2)
    out_path = os.path.join(RESULTS_DIR, "thr2_mod6_n10_probe2.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print(f"\nTermine en {time.time() - t0:.2f}s. Sauvegarde : {out_path}")


if __name__ == "__main__":
    main()
