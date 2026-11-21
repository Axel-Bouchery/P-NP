r"""
Chasse aux records mondiaux sur F_2 : Flip-Graph bit-vectoriel a Faisceau (Beam Search)
+ Collage exact de blocs (Strassen, Laderman, sous-blocs optimaux).

Principe :
  1. Collage exact de tenseurs T_{<m, k, n1>} (+) T_{<m, k, n2>} -> T_{<m, k, n1+n2>}
     (et sur les axes m et k par permutation cyclique des modes (u, v, w)).
     Exemple : <3,4,5> peut demarrer depuis <3,4,2> + <3,4,3> ou <3,3,5> + <3,1,5>
     a un rang deja tres inferieur au schema naif (60).
  2. Beam Search multi-rondes (16 coeurs Numba) :
     - Pool elitiste de `beam_width` schemas distincts (dedoublonnes par hash 64-bit
       ordre-invariant du multi-ensemble des termes) ;
     - A chaque ronde, les 16 workers partent de membres du beam, executent des
       marches de flips XOR + reduction de tranches Gauss + lookahead triangle XOR
       + transitions PLUS adaptatives (Arai et al. 2024) ;
     - Tout nouveau schema distinct au meilleur rang (ou R_min + 1) rejoint le beam.
  3. Verification exacte 100% independante sur T_{<m,k,n>} mod 2 pour chaque record.
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
    WITNESS_DIR,
    build_schoolbook_f2_bits,
    do_plus_transition_f2,
    reduce_slices_gauss_f2,
    try_greedy_lookahead_reduction,
    verify_scheme_f2_bits,
)


def transpose_scheme_axes(u, v, w, m, k, n, target_order):
    """
    Reindexe un schema (u, v, w) de <m, k, n> selon une permutation de (0:m, 1:k, 2:n) :
    - Si on veut coller le long du 3e axe, on peut permuter cycliquement ou transposer
      la multiplication matricielle : T_{<m,k,n>} <-> T_{<k,n,m>} <-> T_{<n,m,k>}
      (par rotation des modes u->v->w->u) et T_{<m,k,n>} <-> T_{<n,k,m>} (transposition A*B -> B^T*A^T).
    Ici on fournit la rotation cyclique simple et la transposition simple.
    """
    pass


def glue_along_n(u1, v1, w1, m, k, n1, u2, v2, w2, n2):
    """
    Colle exactement un schema de <m, k, n1> et un schema de <m, k, n2>
    en un schema exact de <m, k, n1 + n2> en placant B = [B1 | B2] et C = [C1 | C2].
    - u : indices i*k + t in 0..m*k-1 (identiques !)
    - v1 : indices t*n1 + j -> t*(n1+n2) + j
    - v2 : indices t*n2 + j -> t*(n1+n2) + (n1 + j)
    - w1 : indices i*n1 + j -> i*(n1+n2) + j
    - w2 : indices i*n2 + j -> i*(n1+n2) + (n1 + j)
    """
    n_tot = n1 + n2
    assert m * k <= 64 and k * n_tot <= 64 and m * n_tot <= 64

    def remap_cols(arr, rows, cols_old, cols_new, col_offset):
        out = np.zeros(len(arr), dtype=np.uint64)
        for idx, val in enumerate(arr):
            iv = int(val)
            nv = 0
            for r in range(rows):
                for c in range(cols_old):
                    if (iv >> (r * cols_old + c)) & 1:
                        nv |= 1 << (r * cols_new + col_offset + c)
            out[idx] = np.uint64(nv)
        return out

    v1_r = remap_cols(v1, k, n1, n_tot, 0)
    v2_r = remap_cols(v2, k, n2, n_tot, n1)
    w1_r = remap_cols(w1, m, n1, n_tot, 0)
    w2_r = remap_cols(w2, m, n2, n_tot, n1)

    u_out = np.concatenate([u1, u2]).astype(np.uint64)
    v_out = np.concatenate([v1_r, v2_r]).astype(np.uint64)
    w_out = np.concatenate([w1_r, w2_r]).astype(np.uint64)
    return u_out, v_out, w_out


def rotate_cycle_scheme(u, v, w, m, k, n):
    """
    Rotation cyclique exacte de trace Tr(A B C) = Tr(B C A) :
    Transforme un schema de <m, k, n> en un schema exact de <k, n, m> :
      - A (m x k) devient le 3e mode (C^T : m x k -> transpose k x m ?).
    Verification directe :
      T_{<m,k,n>}((i,t), (t',j), (i',j')) = delta_{t,t'} delta_{j,j'} delta_{i,i'}.
      En renommant (u, v, w) -> (v, w^T, u^T) :
      v est (k x n), w^T est (n x m), u^T est (m x k)? Non : pour <k, n, m>,
      U' doit etre (k x n), V' (n x m), W' (k x m).
      Avec w (m x n) -> w^T (n x m) et u (m x k) -> u^T (k x m), on obtient bien <k, n, m> !
    """
    def transpose_bits(arr, rows, cols):
        out = np.zeros(len(arr), dtype=np.uint64)
        for idx, val in enumerate(arr):
            iv = int(val)
            nv = 0
            for r in range(rows):
                for c in range(cols):
                    if (iv >> (r * cols + c)) & 1:
                        nv |= 1 << (c * rows + r)
            out[idx] = np.uint64(nv)
        return out

    u_new = v.copy()
    v_new = transpose_bits(w, m, n)
    w_new = transpose_bits(u, m, k)
    return u_new, v_new, w_new, k, n, m


def scheme_hash_f2(u, v, w, R):
    """Hash 64-bit ordre-invariant d'un schema exact (u, v, w) de rang R."""
    h_acc = 0
    for r in range(R):
        x = (int(u[r]) * 0x9E3779B97F4A7C15) ^ (int(v[r]) * 0xBF58476D1CE4E5B9) ^ (int(w[r]) * 0x94D049BB133111EB)
        x = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9) & 0xFFFFFFFFFFFFFFFF
        x = ((x ^ (x >> 27)) * 0x94D049BB133111EB) & 0xFFFFFFFFFFFFFFFF
        h_acc ^= x ^ (x >> 31)
    return int(h_acc & 0xFFFFFFFFFFFFFFFF)


@njit(parallel=True, cache=True)
def beam_round_f2_parallel(beam_u, beam_v, beam_w, beam_R, n_beam, target_R, n_workers, steps_per_worker, plus_int, max_gap, seed_base):
    """
    Execute une ronde parallele sur `n_workers` : chaque worker demarre d'un
    membre du faisceau `beam` (modulo n_beam), effectue `steps_per_worker` flips
    avec reduction de tranches, lookahead et PLUS adaptatifs, et retourne son
    meilleur schema trouve.
    """
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
            u[r] = beam_u[bi, r]
            v[r] = beam_v[bi, r]
            w[r] = beam_w[bi, r]
            best_u_all[wid, r] = u[r]
            best_v_all[wid, r] = v[r]
            best_w_all[wid, r] = w[r]

        R = R0
        best_R = R0
        best_R_all[wid] = best_R

        # Si on part de plus haut que la premiere ronde, effectuer 1 PLUS d'entree sur la moitie des workers
        if (wid & 1) == 1 and R + 1 < max_R:
            R, st = do_plus_transition_f2(u, v, w, R, max_R, st)

        pair_mode = np.zeros(max_R * max_R, dtype=np.int8)
        pair_i = np.zeros(max_R * max_R, dtype=np.int16)
        pair_j = np.zeros(max_R * max_R, dtype=np.int16)
        since_improve = 0
        since_plus = 0

        for step in range(steps_per_worker):
            if best_R <= target_R:
                break

            if (step & 15) == 0:
                R_new, hit = try_greedy_lookahead_reduction(u, v, w, R)
                if hit:
                    R = reduce_slices_gauss_f2(u, v, w, R_new)
                    while True:
                        R2, hit2 = try_greedy_lookahead_reduction(u, v, w, R)
                        if not hit2:
                            break
                        R = reduce_slices_gauss_f2(u, v, w, R2)
                    if R < best_R:
                        best_R = R
                        best_R_all[wid] = best_R
                        for r in range(R):
                            best_u_all[wid, r] = u[r]
                            best_v_all[wid, r] = v[r]
                            best_w_all[wid, r] = w[r]
                        since_improve = 0
                        since_plus = 0
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

            if n_pairs == 0 or since_plus >= plus_int:
                if R < best_R + max_gap and R + 1 < max_R:
                    R, st = do_plus_transition_f2(u, v, w, R, max_R, st)
                    since_plus = 0
                    continue
                else:
                    R = best_R
                    for r in range(R):
                        u[r] = best_u_all[wid, r]; v[r] = best_v_all[wid, r]; w[r] = best_w_all[wid, r]
                    if R + 1 < max_R:
                        R, st = do_plus_transition_f2(u, v, w, R, max_R, st)
                    since_plus = 0
                    since_improve = 0
                    continue

            st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
            pidx = np.int64(st % np.uint64(n_pairs))
            md = pair_mode[pidx]
            i = np.int64(pair_i[pidx]); j = np.int64(pair_j[pidx])
            st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
            dir_bit = np.int64(st & np.uint64(1))
            coin = np.int64((st >> np.uint64(1)) & np.uint64(1))
            if coin != 0:
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
            since_plus += 1
            since_improve += 1

            if (
                u[i] == np.uint64(0) or v[i] == np.uint64(0) or w[i] == np.uint64(0)
                or u[j] == np.uint64(0) or v[j] == np.uint64(0) or w[j] == np.uint64(0)
                or (step & 7) == 0
            ):
                R_new = reduce_slices_gauss_f2(u, v, w, R)
                if R_new < R:
                    R = R_new
                    since_plus = 0
                    if R < best_R:
                        best_R = R
                        best_R_all[wid] = best_R
                        for r in range(R):
                            best_u_all[wid, r] = u[r]; best_v_all[wid, r] = v[r]; best_w_all[wid, r] = w[r]
                        since_improve = 0

            if since_improve > plus_int * 6 and R > best_R:
                R = best_R
                for r in range(R):
                    u[r] = best_u_all[wid, r]; v[r] = best_v_all[wid, r]; w[r] = best_w_all[wid, r]
                since_improve = 0
                since_plus = plus_int

    return best_R_all, best_u_all, best_v_all, best_w_all, flips_done


def run_beam_hunt(m, k, n, u_init, v_init, w_init, target_R, known_record, n_rounds=12, beam_width=16, n_workers=32, steps_per_round=120000, plus_int=350, max_gap=2):
    """Execute un Beam Search complet sur <m, k, n> et verifie 100% exactement le meilleur schema."""
    t0 = time.time()
    R0 = len(u_init)
    max_R = R0 + max_gap + 6
    assert verify_scheme_f2_bits(u_init, v_init, w_init, m, k, n), f"Init invalide pour <{m},{k},{n}>!"

    pool = [(R0, scheme_hash_f2(u_init, v_init, w_init, R0), u_init.copy(), v_init.copy(), w_init.copy())]
    global_best_R = R0
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
            b_u[bi, :r_val] = uu[:r_val]
            b_v[bi, :r_val] = vv[:r_val]
            b_w[bi, :r_val] = ww[:r_val]

        bRs, bus, bvs, bws, fdone = beam_round_f2_parallel(
            b_u, b_v, b_w, b_R, n_beam, target_R, n_workers, steps_per_round, plus_int, max_gap, rnd + 1
        )
        total_flips += int(np.sum(fdone))

        # Fusionner les resultats dans le pool
        seen_h = {item[1] for item in pool}
        for wid in range(n_workers):
            rw = int(bRs[wid])
            uu = bus[wid, :rw].copy()
            vv = bvs[wid, :rw].copy()
            ww = bws[wid, :rw].copy()
            hw = scheme_hash_f2(uu, vv, ww, rw)
            if hw not in seen_h:
                seen_h.add(hw)
                pool.append((rw, hw, uu, vv, ww))

        pool.sort(key=lambda x: x[0])
        pool = pool[:beam_width]
        if pool[0][0] < global_best_R:
            global_best_R = pool[0][0]
            print(f"    [ronde {rnd+1:2d}/{n_rounds}] NOUVEAU MEILLEUR <{m},{k},{n}> : R = {global_best_R} ({total_flips/1e6:.2f}M flips, {time.time()-t0:.2f}s)")
        if global_best_R <= target_R:
            break

    dt = time.time() - t0
    best_R, _, u_b, v_b, w_b = pool[0]
    ok = verify_scheme_f2_bits(u_b, v_b, w_b, m, k, n)
    assert ok, f"Echec verification finale <{m},{k},{n}> R={best_R}!"
    status = "RECORD BATTU !" if best_R < known_record else ("EGALE RECORD !" if best_R == known_record else f"gap +{best_R - known_record}")
    print(
        f"  ==> <{m},{k},{n}> : R_init={R0} -> R_best={best_R} (record={known_record} : {status}) | "
        f"{total_flips/1e6:.2f}M flips en {dt:.2f}s [VERIFIE EXACT]"
    )
    return {
        "format": [m, k, n],
        "R_init": R0,
        "best_R": best_R,
        "known_record": known_record,
        "beat_record": (best_R < known_record),
        "matched_record": (best_R == known_record),
        "total_flips": total_flips,
        "time_s": round(dt, 2),
        "u": [int(x) for x in u_b],
        "v": [int(x) for x in v_b],
        "w": [int(x) for x in w_b],
    }


def main():
    t0 = time.time()
    print("====================================================================================")
    print("  CHASSE BEAM SEARCH F_2 + COLLAGE DE BLOCS SUR LES FORMATS A RECORDS")
    print("====================================================================================")
    results = {}

    # 1. Petits blocs de base (R exacts) pour assembler les grands formats !
    u222, v222, w222 = build_schoolbook_f2_bits(2, 2, 2)
    r222 = run_beam_hunt(2, 2, 2, u222, v222, w222, target_R=7, known_record=7, n_rounds=2, beam_width=4, n_workers=16, steps_per_round=10000, plus_int=150, max_gap=1)
    u222, v222, w222 = np.array(r222["u"], np.uint64), np.array(r222["v"], np.uint64), np.array(r222["w"], np.uint64)

    u223, v223, w223 = build_schoolbook_f2_bits(2, 2, 3)
    r223 = run_beam_hunt(2, 2, 3, u223, v223, w223, target_R=11, known_record=11, n_rounds=2, beam_width=4, n_workers=16, steps_per_round=15000, plus_int=150, max_gap=1)
    u223, v223, w223 = np.array(r223["u"], np.uint64), np.array(r223["v"], np.uint64), np.array(r223["w"], np.uint64)

    # <2, 3, 5> : collage <2, 3, 2> (11) + <2, 3, 3> (15) -> init = 26, cible <= 25 !
    u233, v233, w233 = build_schoolbook_f2_bits(2, 3, 3)
    r233 = run_beam_hunt(2, 3, 3, u233, v233, w233, target_R=15, known_record=15, n_rounds=3, beam_width=8, n_workers=16, steps_per_round=40000, plus_int=250, max_gap=2)
    u233, v233, w233 = np.array(r233["u"], np.uint64), np.array(r233["v"], np.uint64), np.array(r233["w"], np.uint64)

    # Rotation cyclique <2,2,3> -> <2,3,2> (non : cycle <2,2,3> -> <2,3,2> est 2 rotations : <2,2,3>-><2,3,2>-><3,2,2> -> 2x)
    u232, v232, w232, m1, k1, n1 = rotate_cycle_scheme(u223, v223, w223, 2, 2, 3)  # -> <2, 3, 2> !
    assert (m1, k1, n1) == (2, 3, 2)
    u235_init, v235_init, w235_init = glue_along_n(u232, v232, w232, 2, 3, 2, u233, v233, w233, 3)
    r235 = run_beam_hunt(2, 3, 5, u235_init, v235_init, w235_init, target_R=25, known_record=25, n_rounds=6, beam_width=16, n_workers=32, steps_per_round=80000, plus_int=300, max_gap=2)
    results["2_3_5"] = r235

    # <3, 3, 3> -> R = 23
    u333, v333, w333 = build_schoolbook_f2_bits(3, 3, 3)
    r333 = run_beam_hunt(3, 3, 3, u333, v333, w333, target_R=23, known_record=23, n_rounds=4, beam_width=12, n_workers=16, steps_per_round=80000, plus_int=300, max_gap=2)
    u333, v333, w333 = np.array(r333["u"], np.uint64), np.array(r333["v"], np.uint64), np.array(r333["w"], np.uint64)

    # <3, 3, 4> : collage <3, 3, 2> (15) + <3, 3, 2> (15) = 30 -> cible 29 (voire 28 !)
    u332, v332, w332, m2, k2, n2 = rotate_cycle_scheme(u233, v233, w233, 2, 3, 3)  # -> <3, 3, 2>
    assert (m2, k2, n2) == (3, 3, 2)
    u334_init, v334_init, w334_init = glue_along_n(u332, v332, w332, 3, 3, 2, u332, v332, w332, 2)
    r334 = run_beam_hunt(3, 3, 4, u334_init, v334_init, w334_init, target_R=28, known_record=29, n_rounds=8, beam_width=16, n_workers=32, steps_per_round=100000, plus_int=350, max_gap=2)
    results["3_3_4"] = r334

    # <3, 4, 4> : rotation <3,3,4> -> <3,4,3> + <3,4,1> (=12) ; ou collage direct <3,4,2> (20) + <3,4,2> (20) = 40 -> record 38 !
    # <2,3,4> (R=20) -> rotate -> <3,4,2> (R=20) !
    u234 = np.array(run_beam_hunt(2, 3, 4, *build_schoolbook_f2_bits(2, 3, 4), target_R=20, known_record=20, n_rounds=4, beam_width=12, n_workers=32, steps_per_round=60000, plus_int=300, max_gap=2)["u"], np.uint64)
    # Re-recuperer u234, v234, w234 proprement
    r234 = run_beam_hunt(2, 3, 4, *build_schoolbook_f2_bits(2, 3, 4), target_R=20, known_record=20, n_rounds=4, beam_width=12, n_workers=32, steps_per_round=60000, plus_int=300, max_gap=2)
    u234, v234, w234 = np.array(r234["u"], np.uint64), np.array(r234["v"], np.uint64), np.array(r234["w"], np.uint64)
    u342, v342, w342, m3, k3, n3 = rotate_cycle_scheme(u234, v234, w234, 2, 3, 4)  # -> <3, 4, 2>
    assert (m3, k3, n3) == (3, 4, 2)
    u344_init, v344_init, w344_init = glue_along_n(u342, v342, w342, 3, 4, 2, u342, v342, w342, 2)
    r344 = run_beam_hunt(3, 4, 4, u344_init, v344_init, w344_init, target_R=37, known_record=38, n_rounds=10, beam_width=24, n_workers=32, steps_per_round=120000, plus_int=400, max_gap=2)
    results["3_4_4"] = r344

    # <3, 4, 5> : collage <3,4,4> (r344) + <3,4,1> (12) -> init <= 50, cible < 47 !
    u344_b, v344_b, w344_b = np.array(r344["u"], np.uint64), np.array(r344["v"], np.uint64), np.array(r344["w"], np.uint64)
    u341, v341, w341 = build_schoolbook_f2_bits(3, 4, 1)
    u345_init, v345_init, w345_init = glue_along_n(u344_b, v344_b, w344_b, 3, 4, 4, u341, v341, w341, 1)
    r345 = run_beam_hunt(3, 4, 5, u345_init, v345_init, w345_init, target_R=46, known_record=47, n_rounds=14, beam_width=24, n_workers=32, steps_per_round=150000, plus_int=450, max_gap=2)
    results["3_4_5"] = r345

    out = {"results": results, "total_time_s": round(time.time() - t0, 2)}
    out_path = os.path.join(RESULTS_DIR, "flip_graph_f2_beam_records.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print(f"\nTermine en {time.time() - t0:.2f}s. Sauvegarde : {out_path}")


if __name__ == "__main__":
    main()
