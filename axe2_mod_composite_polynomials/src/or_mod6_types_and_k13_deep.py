r"""
1. Classification canonique 100% EXACTE des sous-graphes homogenes de OR_n mod 6 (d = 2) :
   - Type A (V_0 : a_v = 1, b_v = 0) : deja prouve |V_0| <= 9 (4 classes a k = 9, 0 a k = 10).
   - Type B (W_{0,1} : a_v = 0, b_v = 1) : P_2(S) = sum a_{uv} = 1 mod 2 sur P_3(S) = (|S| + sum b_{uv}) = 0 mod 3.
   - Type C (W_{1,1} : a_v = 1, b_v = 1) : P_2(S) = |S| + sum a_{uv} = 1 mod 2 sur P_3(S) = (|S| + sum b_{uv}) = 0 mod 3.
   Utilise le moteur ultra-rapide (DP + pivots F_2 incrementaux + early-exit) !
2. Attaque 2-opt 128-bits sur THR_{13, 2} mod 10 depuis les DEUX temoins K_12 ((4,3,3,2) et (3,3,3,3)).
3. Sonde 4-blocs sur THR_{19, 2} mod 21.
"""
import json
import os
import sys
import time

import numpy as np
from numba import njit, prange

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from composite_mod_m_lab import RESULTS_DIR
from gyori_generalization_probe1 import canonicalize_exact_s_k_only
from piste_2a_thr2_graphs import solve_and_verify_thr2_mod2p
from thr2_mod10_128bit_and_multi_mod import (
    f2_defect_count_128,
    search_kblock_odd_composite_thr2,
    verify_kblock_odd_bruteforce,
)


@njit(parallel=True, cache=True)
def extend_fixed_ab_or_mod6_fast(prev_list, k, a_lin, b_lin):
    r"""
    Moteur exact ultra-rapide (pivots F_2 incrementaux + DP mod 3 + early-exit)
    pour un sous-graphe de OR_k mod 6 ou TOUS les sommets ont le meme terme
    lineaire (a_v, b_v) = (a_lin, b_lin) in F_2 x F_3 \ {(0,0)}.
    Variables F_2 : uniquement les C(k, 2) aretes a_{uv} !
    Pour tout sous-ensemble S != empty (|S| >= 2, car pour |S|=1, (a_lin,b_lin)!=(0,0) est deja satisfait) :
      Si (|S| * b_lin + sum_{u<v in S} b_{uv}) % 3 == 0 :
        sum_{u<v in S} a_{uv} == (1 ^ ((|S| * a_lin) & 1))  (mod 2).
    """
    n_prev = prev_list.shape[0]
    n_ext = 3 ** (k - 1)
    n_edges = k * (k - 1) // 2
    n_prev_edges = (k - 1) * (k - 2) // 2
    rhs_bit = np.int64(1) << n_edges

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

    n_sub = 1 << (k - 1)
    old_edge_bits = np.zeros(n_sub, dtype=np.int64)
    new_edge_bits = np.zeros(n_sub, dtype=np.int64)
    popcnt = np.zeros(n_sub, dtype=np.int32)
    lsb_idx = np.zeros(n_sub, dtype=np.int32)
    prev_mask = np.zeros(n_sub, dtype=np.int32)
    rhs_old = np.zeros(n_sub, dtype=np.int64)
    rhs_with_new = np.zeros(n_sub, dtype=np.int64)

    for m in range(1, n_sub):
        lsb = m & -m
        bi = 0; tmp = lsb
        while tmp > 1:
            tmp >>= 1; bi += 1
        lsb_idx[m] = bi
        prev_mask[m] = m ^ lsb
        pc = popcnt[m ^ lsb] + 1
        popcnt[m] = pc
        if (1 ^ ((pc * a_lin) & 1)) == 1:
            rhs_old[m] = rhs_bit
        if (1 ^ (((pc + 1) * a_lin) & 1)) == 1:
            rhs_with_new[m] = rhs_bit

        ob = np.int64(0)
        for j in range(n_prev_edges):
            ed = idx_prev[j]
            if ((m >> eu[ed]) & 1) and ((m >> ev[ed]) & 1):
                ob |= (np.int64(1) << ed)
        old_edge_bits[m] = ob

        nb = np.int64(0)
        for u in range(k - 1):
            if (m >> u) & 1:
                nb |= (np.int64(1) << idx_new[u])
        new_edge_bits[m] = nb

    order = np.zeros(n_sub - 1, dtype=np.int32)
    pos = 0
    for sz in range(1, k):
        for m in range(1, n_sub):
            if popcnt[m] == sz:
                order[pos] = m; pos += 1

    pe_of = np.full((k - 1, k - 1), -1, dtype=np.int32)
    for j in range(n_prev_edges):
        ed = idx_prev[j]
        pe_of[eu[ed], ev[ed]] = j
        pe_of[ev[ed], eu[ed]] = j

    valid_counts = np.zeros(n_prev, dtype=np.int32)

    for pi in prange(n_prev):
        # p3_quad_old[m] = sum_{u<v in m} b_{uv} mod 3
        p3_q = np.zeros(n_sub, dtype=np.int8)
        for m in range(1, n_sub):
            u = lsb_idx[m]
            pm = prev_mask[m]
            s = int(p3_q[pm])
            tmp = pm; v = 0
            while tmp:
                if tmp & 1:
                    s += int(prev_list[pi, pe_of[u, v]])
                tmp >>= 1; v += 1
            p3_q[m] = s % 3

        piv_base = np.zeros(n_edges, dtype=np.int64)
        for oi in range(n_sub - 1):
            m = order[oi]
            pc = popcnt[m]
            if pc >= 2 and ((pc * b_lin + int(p3_q[m])) % 3 == 0):
                v = old_edge_bits[m] | rhs_old[m]
                for c in range(n_edges):
                    if (v >> c) & 1:
                        if piv_base[c] != 0:
                            v ^= piv_base[c]
                        else:
                            piv_base[c] = v
                            break

        w_new = np.zeros(k - 1, dtype=np.int8)
        sum_new = np.zeros(n_sub, dtype=np.int8)
        piv = np.zeros(n_edges, dtype=np.int64)
        cnt_ok = 0

        for ec in range(n_ext):
            tmp = ec
            for j in range(k - 1):
                w_new[j] = tmp % 3
                tmp //= 3
            sum_new[0] = 0
            for m in range(1, n_sub):
                sum_new[m] = (sum_new[prev_mask[m]] + w_new[lsb_idx[m]]) % 3

            for c in range(n_edges):
                piv[c] = piv_base[c]

            ok = True
            for oi in range(n_sub - 1):
                m = order[oi]
                pc_new = popcnt[m] + 1
                if (pc_new * b_lin + int(p3_q[m]) + int(sum_new[m])) % 3 == 0:
                    v = old_edge_bits[m] | new_edge_bits[m] | rhs_with_new[m]
                    for c in range(n_edges):
                        if (v >> c) & 1:
                            if piv[c] != 0:
                                v ^= piv[c]
                            else:
                                piv[c] = v
                                v = 0
                                break
                    if v == rhs_bit:
                        ok = False
                        break
            if ok:
                cnt_ok += 1
        valid_counts[pi] = cnt_ok

    total_ok = int(np.sum(valid_counts))
    offsets = np.zeros(n_prev, dtype=np.int64)
    acc = 0
    for pi in range(n_prev):
        offsets[pi] = acc; acc += valid_counts[pi]

    out_w = np.zeros((total_ok, n_edges), dtype=np.int8)
    if total_ok == 0:
        return out_w, valid_counts

    for pi in prange(n_prev):
        if valid_counts[pi] == 0:
            continue
        p3_q = np.zeros(n_sub, dtype=np.int8)
        for m in range(1, n_sub):
            u = lsb_idx[m]
            pm = prev_mask[m]
            s = int(p3_q[pm])
            tmp = pm; v = 0
            while tmp:
                if tmp & 1:
                    s += int(prev_list[pi, pe_of[u, v]])
                tmp >>= 1; v += 1
            p3_q[m] = s % 3

        piv_base = np.zeros(n_edges, dtype=np.int64)
        for oi in range(n_sub - 1):
            m = order[oi]
            pc = popcnt[m]
            if pc >= 2 and ((pc * b_lin + int(p3_q[m])) % 3 == 0):
                v = old_edge_bits[m] | rhs_old[m]
                for c in range(n_edges):
                    if (v >> c) & 1:
                        if piv_base[c] != 0:
                            v ^= piv_base[c]
                        else:
                            piv_base[c] = v
                            break

        w_new = np.zeros(k - 1, dtype=np.int8)
        sum_new = np.zeros(n_sub, dtype=np.int8)
        piv = np.zeros(n_edges, dtype=np.int64)
        out_pos = offsets[pi]

        for ec in range(n_ext):
            tmp = ec
            for j in range(k - 1):
                w_new[j] = tmp % 3
                tmp //= 3
            sum_new[0] = 0
            for m in range(1, n_sub):
                sum_new[m] = (sum_new[prev_mask[m]] + w_new[lsb_idx[m]]) % 3

            for c in range(n_edges):
                piv[c] = piv_base[c]

            ok = True
            for oi in range(n_sub - 1):
                m = order[oi]
                pc_new = popcnt[m] + 1
                if (pc_new * b_lin + int(p3_q[m]) + int(sum_new[m])) % 3 == 0:
                    v = old_edge_bits[m] | new_edge_bits[m] | rhs_with_new[m]
                    for c in range(n_edges):
                        if (v >> c) & 1:
                            if piv[c] != 0:
                                v ^= piv[c]
                            else:
                                piv[c] = v
                                v = 0
                                break
                    if v == rhs_bit:
                        ok = False
                        break
            if ok:
                for j in range(n_prev_edges):
                    out_w[out_pos, idx_prev[j]] = prev_list[pi, j]
                for j in range(k - 1):
                    out_w[out_pos, idx_new[j]] = w_new[j]
                out_pos += 1

    return out_w, valid_counts


def run_exact_tree_for_type(label, a_lin, b_lin, max_k=11):
    print(f"\n  --- Arbre canonique EXACT pour {label} : (a_v, b_v) = ({a_lin}, {b_lin}) ---")
    cur = np.array([[0], [1], [2]], dtype=np.int8)
    levels = {"2": {"canon_count": 3}}
    print(f"  [{label}, k =  2] : 3 classes canoniques (mod S_2)")
    all_arr = {2: cur}
    for k in range(3, max_k + 1):
        t1 = time.time()
        ext_w, vcnt = extend_fixed_ab_or_mod6_fast(cur, k, a_lin, b_lin)
        dt_e = time.time() - t1
        t2 = time.time()
        uniq_w, _ = canonicalize_exact_s_k_only(ext_w, k)
        dt_c = time.time() - t2
        n_par = int(cur.shape[0])
        n_val = int(ext_w.shape[0])
        n_can = int(uniq_w.shape[0])
        dead = int(np.sum(vcnt == 0))
        print(
            f"  [{label}, k = {k:2d}] {n_par:6d} parents x 3^{k-1} -> {n_val:7d} valides -> "
            f"{n_can:6d} canoniques (dead={dead}/{n_par}) [ext={dt_e:.2f}s, can={dt_c:.2f}s]"
        )
        levels[str(k)] = {
            "parents": n_par,
            "valid_extensions": n_val,
            "canon_count": n_can,
            "dead_ends": dead,
            "time_ext_s": round(dt_e, 2),
            "time_can_s": round(dt_c, 2),
        }
        cur = uniq_w
        all_arr[k] = cur
        if n_can == 0:
            print(f"  >>> EXTINCTION EXACTE POUR {label} : k_max = {k - 1} ! <<<")
            break
    return levels, all_arr


def main():
    t0 = time.time()
    out = {}

    print("==========================================================================")
    print("  1. CLASSIFICATION CANONIQUE EXACTE DES TYPES B (0,1) ET C (1,1) DE OR_n MOD 6")
    print("==========================================================================")
    lv_B, arr_B = run_exact_tree_for_type("Type_B_W01", 0, 1, max_k=11)
    out["type_B_a0_b1"] = lv_B

    lv_C, arr_C = run_exact_tree_for_type("Type_C_W11", 1, 1, max_k=11)
    out["type_C_a1_b1"] = lv_C

    print("\n==========================================================================")
    print("  2. SONDE 4-BLOCS SUR THR_{19, 2} MOD 21 (AU-DELA DU RECORD n = 18)")
    print("==========================================================================")
    m21_19 = []
    for part in [(6, 6, 6, 1), (6, 6, 4, 3), (5, 5, 5, 4)]:
        t1 = time.time()
        df, w1, w2, npts = search_kblock_odd_composite_thr2(3, 7, np.array(part, dtype=np.int64), 192, 45000)
        dt = time.time() - t1
        n = sum(part)
        print(f"  [THR_2 mod 21, 4-blocs {part}, n={n}] min_defect = {df} ({dt:.2f}s)")
        entry = {"partition": list(part), "n": n, "min_defect": df, "w3": [int(x) for x in w1], "w7": [int(x) for x in w2], "time_s": round(dt, 2)}
        if df == 0:
            verify_kblock_odd_bruteforce(3, 7, part, w1, w2)
            entry["verified"] = True
            print(f"    >>> NOUVEAU RECORD MONDIAL : THR_{{{n}, 2}} MOD 21 RESOLU ! <<<")
        m21_19.append(entry)
    out["thr2_mod21_n19_probe"] = m21_19

    out["total_time_s"] = round(time.time() - t0, 2)
    out_path = os.path.join(RESULTS_DIR, "or_mod6_types_and_k13_deep.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print(f"\nTermine en {time.time() - t0:.2f}s. Sauvegarde : {out_path}")


if __name__ == "__main__":
    main()
