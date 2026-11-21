r"""
Flip-Graph F_2 INCREMENAL + Detection de tranches Gauss NON-DESTRUCTIVE
=======================================================================
Corrige le dilemme critique observe :
  - `reduce_slices_gauss_f2` en place detruit la marche aleatoire.
  - `local_pair_reduce_f2` seul est aveugle aux dependances lineaires d'ordre >= 3
    dans une tranche (ex. v_1 ^ v_2 ^ v_3 == 0 a u fixe).
Solution : tous les 64 flips, on copie (u, v, w) dans des buffers scratch, on lance
`reduce_slices_gauss_f2` sur la COPIE, et on ne remplace (u, v, w) QUE si R_new < R !
Ainsi : 0 perturbation de la marche + 100% capture des dependances d'ordre >= 3 !

Test immediat sur <2, 3, 5> (cible 25), <3, 5, 5> (cible <= 58), <4, 4, 6> (cible < 73).
"""
import json
import os
import sys
import time

import numpy as np
from numba import njit, prange

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from flip_graph_f2_adaptive import local_pair_reduce_f2
from flip_graph_f2_beam import scheme_hash_f2
from flip_graph_f2_bitvec import (
    RESULTS_DIR,
    WITNESS_DIR,
    build_schoolbook_f2_bits,
    do_plus_transition_f2,
    reduce_slices_gauss_f2,
    verify_scheme_f2_bits,
)
from flip_graph_f2_fast_inc import rebuild_pairs_full, update_pairs_two


@njit(parallel=True, cache=True)
def nondest_flip_round_f2(beam_u, beam_v, beam_w, beam_R, n_beam, target_R, n_workers, steps_per_worker, plus_thresh, max_gap, seed_base):
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
        su = np.zeros(max_R, dtype=np.uint64)
        sv = np.zeros(max_R, dtype=np.uint64)
        sw = np.zeros(max_R, dtype=np.uint64)

        for r in range(R0):
            u[r] = beam_u[bi, r]; v[r] = beam_v[bi, r]; w[r] = beam_w[bi, r]
            best_u_all[wid, r] = u[r]; best_v_all[wid, r] = v[r]; best_w_all[wid, r] = w[r]

        R = R0
        best_R = R0
        best_R_all[wid] = best_R

        p_mode = np.zeros(max_R * max_R, dtype=np.int8)
        p_i = np.zeros(max_R * max_R, dtype=np.int16)
        p_j = np.zeros(max_R * max_R, dtype=np.int16)
        n_pairs = rebuild_pairs_full(u, v, w, R, p_mode, p_i, p_j)

        since_reduce = 0

        for step in range(steps_per_worker):
            if best_R <= target_R:
                break

            # Verification Gauss de tranche NON-DESTRUCTIVE tous les 64 pas
            if (step & 63) == 0:
                for r in range(R):
                    su[r] = u[r]; sv[r] = v[r]; sw[r] = w[r]
                R_g = reduce_slices_gauss_f2(su, sv, sw, R)
                if R_g < R:
                    R = R_g
                    for r in range(R):
                        u[r] = su[r]; v[r] = sv[r]; w[r] = sw[r]
                    since_reduce = 0
                    n_pairs = rebuild_pairs_full(u, v, w, R, p_mode, p_i, p_j)
                    if R < best_R:
                        best_R = R
                        best_R_all[wid] = best_R
                        for r in range(R):
                            best_u_all[wid, r] = u[r]; best_v_all[wid, r] = v[r]; best_w_all[wid, r] = w[r]
                        if best_R <= target_R:
                            break

            if n_pairs == 0 or since_reduce >= plus_thresh:
                if R < best_R + max_gap and R + 1 < max_R:
                    R, st = do_plus_transition_f2(u, v, w, R, max_R, st)
                    n_pairs = rebuild_pairs_full(u, v, w, R, p_mode, p_i, p_j)
                    since_reduce = 0
                    continue
                else:
                    R = best_R
                    for r in range(R):
                        u[r] = best_u_all[wid, r]; v[r] = best_v_all[wid, r]; w[r] = best_w_all[wid, r]
                    st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
                    n_pl = 1 + int(st % np.uint64(min(3, max_gap)))
                    for _ in range(n_pl):
                        if R + 1 < max_R:
                            R, st = do_plus_transition_f2(u, v, w, R, max_R, st)
                    n_pairs = rebuild_pairs_full(u, v, w, R, p_mode, p_i, p_j)
                    since_reduce = 0
                    continue

            st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
            pidx = np.int64(st % np.uint64(n_pairs))
            md = p_mode[pidx]
            i = np.int64(p_i[pidx]); j = np.int64(p_j[pidx])
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
                n_pairs = rebuild_pairs_full(u, v, w, R, p_mode, p_i, p_j)
                if R < best_R:
                    best_R = R
                    best_R_all[wid] = best_R
                    for r in range(R):
                        best_u_all[wid, r] = u[r]; best_v_all[wid, r] = v[r]; best_w_all[wid, r] = w[r]
            else:
                n_pairs = update_pairs_two(u, v, w, R, i, j, p_mode, p_i, p_j, n_pairs)

    return best_R_all, best_u_all, best_v_all, best_w_all, flips_done


def run_nondest(m, k, n, init_schemes, target_R, known_record, n_rounds=16, beam_width=40, n_workers=48, steps_per_round=800000, plus_thresh=6000, max_gap=4):
    t0 = time.time()
    max_R0 = max(len(s[0]) for s in init_schemes)
    max_R = max_R0 + max_gap + 6
    pool = []
    seen_h = set()
    for (ui, vi, wi) in init_schemes:
        assert verify_scheme_f2_bits(ui, vi, wi, m, k, n)
        r0 = len(ui)
        h = scheme_hash_f2(ui, vi, wi, r0)
        if h not in seen_h:
            seen_h.add(h)
            pool.append((r0, h, ui.copy(), vi.copy(), wi.copy()))
    pool.sort(key=lambda x: x[0])
    best_R = pool[0][0]
    R_start = best_R
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

        bRs, bus, bvs, bws, fdone = nondest_flip_round_f2(
            b_u, b_v, b_w, b_R, n_beam, target_R, n_workers, steps_per_round, plus_thresh, max_gap, rnd + 1
        )
        total_flips += int(np.sum(fdone))
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
    speed = total_flips / max(1e-6, dt) / 1e6
    status = "RECORD BATTU !" if best_R < known_record else ("EGALE RECORD !" if best_R == known_record else f"gap +{best_R - known_record}")
    print(f"  ==> <{m},{k},{n}> : R_init={R_start} -> R_best={best_R} (record={known_record} : {status}) | {total_flips/1e6:.2f}M flips en {dt:.2f}s ({speed:.2f} M flips/s) [VERIFIE EXACT]")
    return {"format": [m, k, n], "best_R": best_R, "known_record": known_record, "u": u_b, "v": v_b, "w": w_b}


def main():
    t0 = time.time()
    print("====================================================================================")
    print("  FLIP-GRAPH F_2 NON-DESTRUCTIF + GAUSS SUR COPIE : <2,3,5> & <4,4,6> & <4,5,5>")
    print("====================================================================================")

    # 1. Test sur <2, 3, 5> (doit enfin atteindre 25 !)
    r235 = run_nondest(2, 3, 5, [build_schoolbook_f2_bits(2, 3, 5)], target_R=25, known_record=25, n_rounds=8, beam_width=32, n_workers=48, steps_per_round=400000, plus_thresh=4000, max_gap=4)

    # 2. <4, 4, 6> depuis notre temoin R = 73 (cible < 73 !)
    d446 = np.load(os.path.join(WITNESS_DIR, "f2_best_4x4x6_R73.npz"))
    r446 = run_nondest(4, 4, 6, [(d446["u"], d446["v"], d446["w"])], target_R=66, known_record=73, n_rounds=16, beam_width=48, n_workers=48, steps_per_round=800000, plus_thresh=8000, max_gap=4)
    np.savez(os.path.join(WITNESS_DIR, f"f2_best_4x4x6_R{r446['best_R']}.npz"), u=r446["u"], v=r446["v"], w=r446["w"], m=4, k=4, n=6, R=r446["best_R"])

    # 3. <4, 5, 5> depuis notre temoin R = 75 (cible <= 73 !)
    d455 = np.load(os.path.join(WITNESS_DIR, "f2_best_4x5x5_R75.npz"))
    r455 = run_nondest(4, 5, 5, [(d455["u"], d455["v"], d455["w"])], target_R=72, known_record=73, n_rounds=16, beam_width=48, n_workers=48, steps_per_round=800000, plus_thresh=8000, max_gap=4)
    np.savez(os.path.join(WITNESS_DIR, f"f2_best_4x5x5_R{r455['best_R']}.npz"), u=r455["u"], v=r455["v"], w=r455["w"], m=4, k=5, n=5, R=r455["best_R"])

    print(f"\nTermine en {time.time() - t0:.2f}s.")


if __name__ == "__main__":
    main()
