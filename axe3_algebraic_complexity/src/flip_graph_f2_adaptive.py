r"""
Flip-Graph adaptatif F_2 (Arai-Ichikawa-Hukushima 2024) sans destruction de marche :
  - Chaque flip XOR ne fait QUE la fusion locale O(R) si deux termes partagent 2 modes
    (aucun echelonnement Gauss destructif pendant la marche aleatoire !).
  - Transitions PLUS probabilistes/adaptatives avec gap jusqu'a +4 au-dessus de best_R
    et retour probabiliste au meilleur schema.
  - Test immediat sur <2,3,5> (cible 25), <3,3,4> (cible 29), <3,4,4> (cible 38),
    <3,4,5> (cible 47).
"""
import json
import os
import sys
import time

import numpy as np
from numba import njit, prange

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from flip_graph_f2_bitvec import (
    RESULTS_DIR,
    build_schoolbook_f2_bits,
    compact_nonzero,
    do_plus_transition_f2,
    reduce_slices_gauss_f2,
    verify_scheme_f2_bits,
)
from flip_graph_f2_beam import glue_along_n, rotate_cycle_scheme, scheme_hash_f2


@njit(cache=True)
def local_pair_reduce_f2(u, v, w, R, idx_a, idx_b):
    """
    Apres modification des termes idx_a et idx_b, verifie uniquement en O(R)
    s'ils sont devenus nuls ou identiques sur 2 modes a un autre terme k,
    et reduit en cascade sans perturber les autres termes !
    """
    changed = True
    while changed:
        changed = False
        for t in (idx_a, idx_b):
            if t < 0 or t >= R:
                continue
            ut, vt, wt = u[t], v[t], w[t]
            if ut == np.uint64(0) or vt == np.uint64(0) or wt == np.uint64(0):
                u[t] = np.uint64(0); v[t] = np.uint64(0); w[t] = np.uint64(0)
                changed = True
                continue
            for k in range(R):
                if k == t:
                    continue
                uk, vk, wk = u[k], v[k], w[k]
                if uk == np.uint64(0):
                    continue
                if ut == uk and vt == vk:
                    w[k] ^= wt
                    u[t] = np.uint64(0); v[t] = np.uint64(0); w[t] = np.uint64(0)
                    if w[k] == np.uint64(0):
                        u[k] = np.uint64(0); v[k] = np.uint64(0)
                    idx_a = k; idx_b = -1
                    changed = True
                    break
                elif ut == uk and wt == wk:
                    v[k] ^= vt
                    u[t] = np.uint64(0); v[t] = np.uint64(0); w[t] = np.uint64(0)
                    if v[k] == np.uint64(0):
                        u[k] = np.uint64(0); w[k] = np.uint64(0)
                    idx_a = k; idx_b = -1
                    changed = True
                    break
                elif vt == vk and wt == wk:
                    u[k] ^= ut
                    u[t] = np.uint64(0); v[t] = np.uint64(0); w[t] = np.uint64(0)
                    if u[k] == np.uint64(0):
                        v[k] = np.uint64(0); w[k] = np.uint64(0)
                    idx_a = k; idx_b = -1
                    changed = True
                    break
        if changed:
            R = compact_nonzero(u, v, w, R)
            idx_a = min(idx_a, R - 1)
            idx_b = -1
    return R


@njit(parallel=True, cache=True)
def adaptive_flip_round_f2(beam_u, beam_v, beam_w, beam_R, n_beam, target_R, n_workers, steps_per_worker, plus_thresh, max_gap, seed_base):
    max_R = beam_u.shape[1]
    best_R_all = np.full(n_workers, 999999, dtype=np.int32)
    best_u_all = np.zeros((n_workers, max_R), dtype=np.uint64)
    best_v_all = np.zeros((n_workers, max_R), dtype=np.uint64)
    best_w_all = np.zeros((n_workers, max_R), dtype=np.uint64)
    flips_done = np.zeros(n_workers, dtype=np.int64)

    for wid in prange(n_workers):
        bi = wid % n_beam
        R0 = int(beam_R[bi])
        st = np.uint64(0x9E3779B97F4A7C15) ^ np.uint64((wid + 1) * 6364136223846793005 + seed_base * 1000003 + R0 * 7919)

        u = np.zeros(max_R, dtype=np.uint64)
        v = np.zeros(max_R, dtype=np.uint64)
        w = np.zeros(max_R, dtype=np.uint64)
        for r in range(R0):
            u[r] = beam_u[bi, r]; v[r] = beam_v[bi, r]; w[r] = beam_w[bi, r]
            best_u_all[wid, r] = u[r]; best_v_all[wid, r] = v[r]; best_w_all[wid, r] = w[r]

        R = R0
        best_R = R0
        best_R_all[wid] = best_R

        pair_mode = np.zeros(max_R * max_R, dtype=np.int8)
        pair_i = np.zeros(max_R * max_R, dtype=np.int16)
        pair_j = np.zeros(max_R * max_R, dtype=np.int16)

        since_reduce = 0

        for step in range(steps_per_worker):
            if best_R <= target_R:
                break

            n_pairs = 0
            for i in range(R):
                ui, vi, wi = u[i], v[i], w[i]
                for j in range(i + 1, R):
                    if ui == u[j]:
                        pair_mode[n_pairs] = 0; pair_i[n_pairs] = i; pair_j[n_pairs] = j; n_pairs += 1
                    elif vi == v[j]:
                        pair_mode[n_pairs] = 1; pair_i[n_pairs] = i; pair_j[n_pairs] = j; n_pairs += 1
                    elif wi == w[j]:
                        pair_mode[n_pairs] = 2; pair_i[n_pairs] = i; pair_j[n_pairs] = j; n_pairs += 1

            # Transition PLUS adaptative (si plus de paires ou stagnation sans reduction)
            if n_pairs == 0 or since_reduce >= plus_thresh:
                if R < best_R + max_gap and R + 1 < max_R:
                    R, st = do_plus_transition_f2(u, v, w, R, max_R, st)
                    since_reduce = 0
                    continue
                else:
                    # Revenir a best_R puis faire 1 ou 2 PLUS pour explorer un autre bassin
                    R = best_R
                    for r in range(R):
                        u[r] = best_u_all[wid, r]; v[r] = best_v_all[wid, r]; w[r] = best_w_all[wid, r]
                    st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
                    n_pl = 1 + int(st % np.uint64(min(2, max_gap)))
                    for _ in range(n_pl):
                        if R + 1 < max_R:
                            R, st = do_plus_transition_f2(u, v, w, R, max_R, st)
                    since_reduce = 0
                    continue

            st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
            pidx = np.int64(st % np.uint64(n_pairs))
            md = pair_mode[pidx]
            i = np.int64(pair_i[pidx]); j = np.int64(pair_j[pidx])
            st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
            dir_bit = np.int64(st & np.uint64(1))
            if ((st >> np.uint64(1)) & np.uint64(1)) != np.uint64(0):
                tmp = i; i = j; j = tmp

            if md == 0:
                if dir_bit == 0:
                    v[i] ^= v[j]; w[j] ^= w[i]
                else:
                    w[i] ^= w[j]; v[j] ^= v[i]
            elif md == 1:
                if dir_bit == 0:
                    u[i] ^= u[j]; w[j] ^= w[i]
                else:
                    w[i] ^= w[j]; u[j] ^= u[i]
            else:
                if dir_bit == 0:
                    u[i] ^= u[j]; v[j] ^= v[i]
                else:
                    v[i] ^= v[j]; u[j] ^= u[i]

            flips_done[wid] += 1
            since_reduce += 1

            R_new = local_pair_reduce_f2(u, v, w, R, i, j)
            if R_new < R:
                R = R_new
                since_reduce = 0
                if R < best_R:
                    best_R = R
                    best_R_all[wid] = best_R
                    for r in range(R):
                        best_u_all[wid, r] = u[r]; best_v_all[wid, r] = v[r]; best_w_all[wid, r] = w[r]

    return best_R_all, best_u_all, best_v_all, best_w_all, flips_done


def run_adaptive(m, k, n, u_init, v_init, w_init, target_R, known_record, n_rounds=10, beam_width=24, n_workers=32, steps_per_round=150000, plus_thresh=600, max_gap=3):
    t0 = time.time()
    R0 = len(u_init)
    max_R = R0 + max_gap + 6
    assert verify_scheme_f2_bits(u_init, v_init, w_init, m, k, n)

    pool = [(R0, scheme_hash_f2(u_init, v_init, w_init, R0), u_init.copy(), v_init.copy(), w_init.copy())]
    best_R = R0
    total_flips = 0

    for rnd in range(n_rounds):
        n_beam = min(beam_width, len(pool))
        b_u = np.zeros((n_beam, max_R), dtype=np.uint64)
        b_v = np.zeros((n_beam, max_R), dtype=np.uint64)
        b_w = np.zeros((n_beam, max_R), dtype=np.uint64)
        b_R = np.zeros(n_beam, dtype=np.int32)
        for bi in range(n_beam):
            r_val, _, uu, vv, ww = pool[bi]
            b_R[bi] = r_val
            b_u[bi, :r_val] = uu[:r_val]; b_v[bi, :r_val] = vv[:r_val]; b_w[bi, :r_val] = ww[:r_val]

        bRs, bus, bvs, bws, fdone = adaptive_flip_round_f2(
            b_u, b_v, b_w, b_R, n_beam, target_R, n_workers, steps_per_round, plus_thresh, max_gap, rnd + 1
        )
        total_flips += int(np.sum(fdone))
        seen_h = {item[1] for item in pool}
        for wid in range(n_workers):
            rw = int(bRs[wid])
            uu = bus[wid, :rw].copy(); vv = bvs[wid, :rw].copy(); ww = bws[wid, :rw].copy()
            hw = scheme_hash_f2(uu, vv, ww, rw)
            if hw not in seen_h:
                seen_h.add(hw)
                pool.append((rw, hw, uu, vv, ww))
        pool.sort(key=lambda x: x[0])
        pool = pool[:beam_width]
        if pool[0][0] < best_R:
            best_R = pool[0][0]
            print(f"    [ronde {rnd+1:2d}/{n_rounds}] NOUVEAU MEILLEUR <{m},{k},{n}> : R = {best_R} ({total_flips/1e6:.2f}M flips, {time.time()-t0:.2f}s)")
        if best_R <= target_R:
            break

    dt = time.time() - t0
    best_R, _, u_b, v_b, w_b = pool[0]
    assert verify_scheme_f2_bits(u_b, v_b, w_b, m, k, n)
    status = "RECORD BATTU !" if best_R < known_record else ("EGALE RECORD !" if best_R == known_record else f"gap +{best_R - known_record}")
    print(f"  ==> <{m},{k},{n}> : R_init={R0} -> R_best={best_R} (record={known_record} : {status}) | {total_flips/1e6:.2f}M flips en {dt:.2f}s [VERIFIE EXACT]")
    return {"format": [m, k, n], "best_R": best_R, "known_record": known_record, "u": u_b, "v": v_b, "w": w_b}


def main():
    t0 = time.time()
    print("====================================================================================")
    print("  TEST FLIP-GRAPH ADAPTATIF F_2 SANS DESTRUCTION DE MARCHE (ARAI ET AL. 2024)")
    print("====================================================================================")

    # Test 1 : <2, 3, 5> depuis naif (R=30) -> doit atteindre 25 !
    r235 = run_adaptive(2, 3, 5, *build_schoolbook_f2_bits(2, 3, 5), target_R=25, known_record=25, n_rounds=8, beam_width=24, n_workers=32, steps_per_round=150000, plus_thresh=600, max_gap=3)

    # Test 2 : <3, 3, 4> depuis naif (R=36) -> doit atteindre 29 !
    r334 = run_adaptive(3, 3, 4, *build_schoolbook_f2_bits(3, 3, 4), target_R=29, known_record=29, n_rounds=10, beam_width=24, n_workers=32, steps_per_round=200000, plus_thresh=700, max_gap=3)

    # Test 3 : <3, 4, 4> depuis naif (R=48) -> record 38 !
    r344 = run_adaptive(3, 4, 4, *build_schoolbook_f2_bits(3, 4, 4), target_R=38, known_record=38, n_rounds=12, beam_width=32, n_workers=32, steps_per_round=250000, plus_thresh=900, max_gap=4)

    print(f"\nTermine en {time.time() - t0:.2f}s.")


if __name__ == "__main__":
    main()
