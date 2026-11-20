"""
Piste 3 (Phase 6) : Graphe de Flips GUIDE PAR PARCIMONIE + WARM-START EXACT (F_3 & F_5)
========================================================================================

Ameliorations majeures tirees du diagnostic de Phase 5 :
  1. Selection guidee par parcimonie (K=8 candidats par pas, biais vers delta_nnz <= 0) :
     empeche la diffusion vers le bruit dense dans F_p^N (nnz reste ~110 au lieu de 240,
     multipliant par 5-10x le nombre de paires colineaires et de reductions !).
  2. Warm-start EXACT depuis Laderman (<3,3,3> R=23, nnz=123, valide sur Z donc F_3 et F_5)
     et Strassen (<2,2,2> R=7) embarque par blocs :
       * <2, 3, 4> : depart exact R=22 -> calibration 20 -> chasse 19 (F_3, F_5)
       * <3, 3, 3> : depart exact R=23 (Laderman) -> chasse 22 (F_3, F_5)
       * <2, 4, 4> : depart exact R=28 (4x Strassen) -> calibration 26 -> chasse 25 (F_3, F_5)
       * <3, 3, 4> : depart exact R=32 (Laderman + 9) -> calibration 29 -> chasse 28 (F_3, F_5)
  3. Pool d'elites etendu (64 schemas distincts) + reduction exacte de tranches + PLUS (+1/+2).
"""
import concurrent.futures
import json
import os
import sys
import time

import numpy as np
from numba import njit

from flip_graph_fp import (
    verify_scheme_fp,
    make_inv_table,
    update_term_keys,
    reduce_fast_and_slices,
    do_flip,
    do_plus,
    scheme_signature,
)

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = sys.argv[1] if len(sys.argv) > 1 else HERE


def build_strassen_z():
    """Schema exact de Strassen (1969) pour <2, 2, 2> (R=7) sur Z."""
    # Indices: a_00=0, a_01=1, a_10=2, a_11=3 ; meme ordre pour B et C
    U = np.array([
        [1, 0, 0, 1],   # m1 = (a00 + a11)(b00 + b11)
        [0, 0, 1, 1],   # m2 = (a10 + a11) b00
        [1, 0, 0, 0],   # m3 = a00 (b01 - b11)
        [0, 0, 0, 1],   # m4 = a11 (b10 - b00)
        [1, 1, 0, 0],   # m5 = (a00 + a01) b11
        [-1, 0, 1, 0],  # m6 = (a10 - a00)(b00 + b01)
        [0, 1, 0, -1],  # m7 = (a01 - a11)(b10 + b11)
    ], dtype=np.int16)
    V = np.array([
        [1, 0, 0, 1],
        [1, 0, 0, 0],
        [0, 1, 0, -1],
        [-1, 0, 1, 0],
        [0, 0, 0, 1],
        [1, 1, 0, 0],
        [0, 0, 1, 1],
    ], dtype=np.int16)
    W = np.array([
        [1, 0, 0, 1],   # m1 -> c00, c11
        [0, 0, 1, -1],  # m2 -> c10, -c11
        [0, 1, 0, 1],   # m3 -> c01, c11
        [1, 0, 1, 0],   # m4 -> c00, c10
        [-1, 1, 0, 0],  # m5 -> -c00, c01
        [0, 0, 0, 1],   # m6 -> c11
        [1, 0, 0, 0],   # m7 -> c00
    ], dtype=np.int16)
    return U, V, W


def build_laderman_z():
    """Schema exact de Laderman (1976) pour <3, 3, 3> (R=23) sur Z."""
    def a(i, j): return (i - 1) * 3 + (j - 1)
    def b(i, j): return (i - 1) * 3 + (j - 1)
    def c(i, j): return (i - 1) * 3 + (j - 1)

    U = np.zeros((23, 9), dtype=np.int16)
    V = np.zeros((23, 9), dtype=np.int16)
    W = np.zeros((23, 9), dtype=np.int16)

    def set_row(M, r, pairs):
        for idx, val in pairs:
            M[r - 1, idx] = val

    # m1..m23 (1-indexed)
    set_row(U, 1, [(a(1,1),1),(a(1,2),1),(a(1,3),1),(a(2,1),-1),(a(2,2),-1),(a(3,2),-1),(a(3,3),-1)])
    set_row(V, 1, [(b(2,2),1)])

    set_row(U, 2, [(a(1,1),1),(a(2,1),-1)])
    set_row(V, 2, [(b(1,2),-1),(b(2,2),1)])

    set_row(U, 3, [(a(2,2),1)])
    set_row(V, 3, [(b(1,1),-1),(b(1,2),1),(b(2,1),1),(b(2,2),-1),(b(2,3),-1),(b(3,1),-1),(b(3,3),1)])

    set_row(U, 4, [(a(1,1),-1),(a(2,1),1),(a(2,2),1)])
    set_row(V, 4, [(b(1,1),1),(b(1,2),-1),(b(2,2),1)])

    set_row(U, 5, [(a(2,1),1),(a(2,2),1)])
    set_row(V, 5, [(b(1,1),-1),(b(1,2),1)])

    set_row(U, 6, [(a(1,1),1)])
    set_row(V, 6, [(b(1,1),1)])

    set_row(U, 7, [(a(1,1),-1),(a(3,1),1),(a(3,2),1)])
    set_row(V, 7, [(b(1,1),1),(b(1,3),-1),(b(2,3),1)])

    set_row(U, 8, [(a(1,1),-1),(a(3,1),1)])
    set_row(V, 8, [(b(1,3),1),(b(2,3),-1)])

    set_row(U, 9, [(a(3,1),1),(a(3,2),1)])
    set_row(V, 9, [(b(1,1),-1),(b(1,3),1)])

    set_row(U, 10, [(a(1,1),1),(a(1,2),1),(a(1,3),1),(a(2,2),-1),(a(2,3),-1),(a(3,1),-1),(a(3,2),-1)])
    set_row(V, 10, [(b(2,3),1)])

    set_row(U, 11, [(a(3,2),1)])
    set_row(V, 11, [(b(1,1),-1),(b(1,3),1),(b(2,1),1),(b(2,2),-1),(b(2,3),-1),(b(3,1),-1),(b(3,2),1)])

    set_row(U, 12, [(a(1,3),-1),(a(3,2),1),(a(3,3),1)])
    set_row(V, 12, [(b(2,2),1),(b(3,1),1),(b(3,2),-1)])

    set_row(U, 13, [(a(1,3),1),(a(3,3),-1)])
    set_row(V, 13, [(b(2,2),1),(b(3,2),-1)])

    set_row(U, 14, [(a(1,3),1)])
    set_row(V, 14, [(b(3,1),1)])

    set_row(U, 15, [(a(3,2),1),(a(3,3),1)])
    set_row(V, 15, [(b(3,1),-1),(b(3,2),1)])

    set_row(U, 16, [(a(1,3),-1),(a(2,2),1),(a(2,3),1)])
    set_row(V, 16, [(b(2,3),1),(b(3,1),1),(b(3,3),-1)])

    set_row(U, 17, [(a(1,3),1),(a(2,3),-1)])
    set_row(V, 17, [(b(2,3),1),(b(3,3),-1)])

    set_row(U, 18, [(a(2,2),1),(a(2,3),1)])
    set_row(V, 18, [(b(3,1),-1),(b(3,3),1)])

    set_row(U, 19, [(a(1,2),1)]); set_row(V, 19, [(b(2,1),1)])
    set_row(U, 20, [(a(2,3),1)]); set_row(V, 20, [(b(3,2),1)])
    set_row(U, 21, [(a(2,1),1)]); set_row(V, 21, [(b(1,3),1)])
    set_row(U, 22, [(a(3,1),1)]); set_row(V, 22, [(b(1,2),1)])
    set_row(U, 23, [(a(3,3),1)]); set_row(V, 23, [(b(3,3),1)])

    outputs = {
        c(1,1): [6, 14, 19],
        c(1,2): [1, 4, 5, 6, 12, 14, 15],
        c(1,3): [6, 7, 9, 10, 14, 16, 18],
        c(2,1): [2, 3, 4, 6, 14, 16, 17],
        c(2,2): [2, 4, 5, 6, 20],
        c(2,3): [14, 16, 17, 18, 21],
        c(3,1): [6, 7, 8, 11, 12, 13, 14],
        c(3,2): [12, 13, 14, 15, 22],
        c(3,3): [6, 7, 8, 9, 23],
    }
    for c_idx, ms in outputs.items():
        for m_id in ms:
            W[m_id - 1, c_idx] = 1
    return U, V, W


def embed_block_schemes(m, k, n, p):
    """
    Construit un schema EXACT parcimonieux de depart par decoupage en blocs Strassen / Laderman :
      - <3, 3, 3> : Laderman R = 23
      - <2, 3, 4> : 2x Strassen (<2,2,2> sur cols 0..1 et 2..3) + naif (ligne k=2) -> 7 + 7 + 8 = 22
      - <2, 4, 4> : 4x Strassen (<2,2,2> sur 2x2 blocs) -> 28
      - <3, 3, 4> : Laderman (<3,3,3> sur cols 0..2) + naif (<3,3,1> sur col 3) -> 23 + 9 = 32
    """
    N1, N2, N3 = m * k, k * n, m * n
    Us, Vs, Ws = [], [], []

    def add_sub(U_sub, V_sub, W_sub, rows_i, mids_t, cols_j):
        R_s = U_sub.shape[0]
        ms, ks, ns = len(rows_i), len(mids_t), len(cols_j)
        U_big = np.zeros((R_s, N1), dtype=np.int16)
        V_big = np.zeros((R_s, N2), dtype=np.int16)
        W_big = np.zeros((R_s, N3), dtype=np.int16)
        for r in range(R_s):
            for ii, gi in enumerate(rows_i):
                for tt, gt in enumerate(mids_t):
                    U_big[r, gi * k + gt] = U_sub[r, ii * ks + tt]
            for tt, gt in enumerate(mids_t):
                for jj, gj in enumerate(cols_j):
                    V_big[r, gt * n + gj] = V_sub[r, tt * ns + jj]
            for ii, gi in enumerate(rows_i):
                for jj, gj in enumerate(cols_j):
                    W_big[r, gi * n + gj] = W_sub[r, ii * ns + jj]
        Us.append(U_big); Vs.append(V_big); Ws.append(W_big)

    def add_naive(rows_i, mids_t, cols_j):
        R_s = len(rows_i) * len(mids_t) * len(cols_j)
        U_big = np.zeros((R_s, N1), dtype=np.int16)
        V_big = np.zeros((R_s, N2), dtype=np.int16)
        W_big = np.zeros((R_s, N3), dtype=np.int16)
        r = 0
        for gi in rows_i:
            for gt in mids_t:
                for gj in cols_j:
                    U_big[r, gi * k + gt] = 1
                    V_big[r, gt * n + gj] = 1
                    W_big[r, gi * n + gj] = 1
                    r += 1
        Us.append(U_big); Vs.append(V_big); Ws.append(W_big)

    if (m, k, n) == (3, 3, 3):
        U, V, W = build_laderman_z()
        return (U % p).astype(np.int16), (V % p).astype(np.int16), (W % p).astype(np.int16)
    elif (m, k, n) == (2, 3, 4):
        U2, V2, W2 = build_strassen_z()
        add_sub(U2, V2, W2, [0, 1], [0, 1], [0, 1])
        add_sub(U2, V2, W2, [0, 1], [0, 1], [2, 3])
        add_naive([0, 1], [2], [0, 1, 2, 3])
    elif (m, k, n) == (2, 4, 4):
        U2, V2, W2 = build_strassen_z()
        for t_blk in ([0, 1], [2, 3]):
            for j_blk in ([0, 1], [2, 3]):
                add_sub(U2, V2, W2, [0, 1], t_blk, j_blk)
    elif (m, k, n) == (3, 3, 4):
        U3, V3, W3 = build_laderman_z()
        add_sub(U3, V3, W3, [0, 1, 2], [0, 1, 2], [0, 1, 2])
        add_naive([0, 1, 2], [0, 1, 2], [3])
    else:
        add_naive(list(range(m)), list(range(k)), list(range(n)))

    U = (np.vstack(Us) % p).astype(np.int16)
    V = (np.vstack(Vs) % p).astype(np.int16)
    W = (np.vstack(Ws) % p).astype(np.int16)
    return U, V, W


@njit(cache=True)
def eval_flip_delta_nnz(U, V, W, lw, r, s, mode, orient, alpha, N1, N2, N3, p, inv):
    """Calcule en O(N) la variation de poids de Hamming delta_nnz d'un flip candidat (et s'il annule un vecteur !)."""
    delta = 0
    if mode == 0:
        if orient == 0:
            z1 = 0; z2 = 0
            for i in range(N2):
                old = 1 if V[r, i] != 0 else 0
                nv = (V[r, i] + alpha * V[s, i]) % p
                nw = 1 if nv != 0 else 0
                delta += (nw - old); z1 += nw
            for i in range(N3):
                old = 1 if W[s, i] != 0 else 0
                nv = (W[s, i] - alpha * W[r, i]) % p
                nw = 1 if nv != 0 else 0
                delta += (nw - old); z2 += nw
            if z1 == 0 or z2 == 0:
                return -1000
        else:
            z1 = 0; z2 = 0
            for i in range(N3):
                old = 1 if W[r, i] != 0 else 0
                nv = (W[r, i] + alpha * W[s, i]) % p
                nw = 1 if nv != 0 else 0
                delta += (nw - old); z1 += nw
            for i in range(N2):
                old = 1 if V[s, i] != 0 else 0
                nv = (V[s, i] - alpha * V[r, i]) % p
                nw = 1 if nv != 0 else 0
                delta += (nw - old); z2 += nw
            if z1 == 0 or z2 == 0:
                return -1000
    elif mode == 1:
        if orient == 0:
            z1 = 0; z2 = 0
            for i in range(N1):
                old = 1 if U[r, i] != 0 else 0
                nv = (U[r, i] + alpha * U[s, i]) % p
                nw = 1 if nv != 0 else 0
                delta += (nw - old); z1 += nw
            for i in range(N3):
                old = 1 if W[s, i] != 0 else 0
                nv = (W[s, i] - alpha * W[r, i]) % p
                nw = 1 if nv != 0 else 0
                delta += (nw - old); z2 += nw
            if z1 == 0 or z2 == 0:
                return -1000
        else:
            z1 = 0; z2 = 0
            for i in range(N3):
                old = 1 if W[r, i] != 0 else 0
                nv = (W[r, i] + alpha * W[s, i]) % p
                nw = 1 if nv != 0 else 0
                delta += (nw - old); z1 += nw
            for i in range(N1):
                old = 1 if U[s, i] != 0 else 0
                nv = (U[s, i] - alpha * U[r, i]) % p
                nw = 1 if nv != 0 else 0
                delta += (nw - old); z2 += nw
            if z1 == 0 or z2 == 0:
                return -1000
    else:
        c = (lw[s] * inv[lw[r]]) % p
        ac = (alpha * c) % p
        if orient == 0:
            z1 = 0; z2 = 0
            for i in range(N1):
                old = 1 if U[r, i] != 0 else 0
                nv = (U[r, i] + ac * U[s, i]) % p
                nw = 1 if nv != 0 else 0
                delta += (nw - old); z1 += nw
            for i in range(N2):
                old = 1 if V[s, i] != 0 else 0
                nv = (V[s, i] - alpha * V[r, i]) % p
                nw = 1 if nv != 0 else 0
                delta += (nw - old); z2 += nw
            if z1 == 0 or z2 == 0:
                return -1000
        else:
            z1 = 0; z2 = 0
            for i in range(N2):
                old = 1 if V[r, i] != 0 else 0
                nv = (V[r, i] + ac * V[s, i]) % p
                nw = 1 if nv != 0 else 0
                delta += (nw - old); z1 += nw
            for i in range(N1):
                old = 1 if U[s, i] != 0 else 0
                nv = (U[s, i] - alpha * U[r, i]) % p
                nw = 1 if nv != 0 else 0
                delta += (nw - old); z2 += nw
            if z1 == 0 or z2 == 0:
                return -1000
    return delta


@njit(cache=True)
def guided_flip_worker(
    U_init, V_init, W_init, p, max_steps, burst_steps,
    plus_depth, seed, target_rank, pool_cap=64, k_cand=8
):
    np.random.seed(seed)
    R0 = U_init.shape[0]
    N1 = U_init.shape[1]
    N2 = V_init.shape[1]
    N3 = W_init.shape[1]
    R_max = R0 + plus_depth + 4
    inv = make_inv_table(p)

    U = np.zeros((R_max, N1), dtype=np.int16)
    V = np.zeros((R_max, N2), dtype=np.int16)
    W = np.zeros((R_max, N3), dtype=np.int16)
    ku = np.zeros(R_max, dtype=np.int64)
    kv = np.zeros(R_max, dtype=np.int64)
    kw = np.zeros(R_max, dtype=np.int64)
    lw = np.zeros(R_max, dtype=np.int16)

    for r in range(R0):
        for i in range(N1): U[r, i] = U_init[r, i] % p
        for i in range(N2): V[r, i] = V_init[r, i] % p
        for i in range(N3): W[r, i] = W_init[r, i] % p
        update_term_keys(U, V, W, ku, kv, kw, lw, r, N1, N2, N3, p, inv)
    R = reduce_fast_and_slices(U, V, W, ku, kv, kw, lw, R0, N1, N2, N3, p, inv, True)

    R_best = R
    pool_U = np.zeros((pool_cap, R_max, N1), dtype=np.int16)
    pool_V = np.zeros((pool_cap, R_max, N2), dtype=np.int16)
    pool_W = np.zeros((pool_cap, R_max, N3), dtype=np.int16)
    pool_sig = np.zeros(pool_cap, dtype=np.uint64)
    pool_size = 1
    for r in range(R):
        for i in range(N1): pool_U[0, r, i] = U[r, i]
        for i in range(N2): pool_V[0, r, i] = V[r, i]
        for i in range(N3): pool_W[0, r, i] = W[r, i]
    pool_sig[0] = scheme_signature(ku, kv, kw, R)

    max_pairs = R_max * R_max
    pr = np.zeros(max_pairs, dtype=np.int16)
    ps = np.zeros(max_pairs, dtype=np.int16)
    pm = np.zeros(max_pairs, dtype=np.int8)

    steps_since_best = 0
    total_flips = 0
    total_plus = 0
    total_reds = 0
    restart_cnt = 0

    for step in range(max_steps):
        if R_best <= target_rank:
            break

        n_pairs = 0
        for r in range(R):
            k1 = ku[r]; k2 = kv[r]; k3 = kw[r]
            for s in range(r + 1, R):
                if k1 == ku[s]:
                    pr[n_pairs] = r; ps[n_pairs] = s; pm[n_pairs] = 0; n_pairs += 1
                if k2 == kv[s]:
                    pr[n_pairs] = r; ps[n_pairs] = s; pm[n_pairs] = 1; n_pairs += 1
                if k3 == kw[s]:
                    pr[n_pairs] = r; ps[n_pairs] = s; pm[n_pairs] = 2; n_pairs += 1

        if n_pairs == 0 or steps_since_best >= burst_steps:
            idx = np.random.randint(0, pool_size)
            R = R_best
            for r in range(R):
                for i in range(N1): U[r, i] = pool_U[idx, r, i]
                for i in range(N2): V[r, i] = pool_V[idx, r, i]
                for i in range(N3): W[r, i] = pool_W[idx, r, i]
                update_term_keys(U, V, W, ku, kv, kw, lw, r, N1, N2, N3, p, inv)
            restart_cnt += 1
            do_jump = (n_pairs == 0) or (restart_cnt % 3 != 0)
            if do_jump:
                n_jumps = 1 if (plus_depth <= 1 or restart_cnt % 4 != 0) else plus_depth
                for _ in range(n_jumps):
                    if R + 1 < R_max and R >= 2:
                        r1 = np.random.randint(0, R)
                        r2 = np.random.randint(0, R - 1)
                        if r2 >= r1:
                            r2 += 1
                        pmode = np.random.randint(0, 3)
                        alpha = np.random.randint(1, p)
                        R = do_plus(U, V, W, ku, kv, kw, lw, R, r1, r2, pmode, alpha, N1, N2, N3, p, inv)
                        R = reduce_fast_and_slices(U, V, W, ku, kv, kw, lw, R, N1, N2, N3, p, inv, True)
                        total_plus += 1
            steps_since_best = 0
            continue

        # Selection guidee par parcimonie (75% meilleur parmi k_cand, 25% exploration libre)
        if np.random.randint(0, 4) != 0:
            best_score = 100000
            best_r = pr[0]; best_s = ps[0]; best_m = pm[0]; best_o = 0; best_a = 1
            for _ in range(k_cand):
                pick = np.random.randint(0, n_pairs)
                cr = pr[pick]; cs = ps[pick]; cm = pm[pick]
                co = np.random.randint(0, 2)
                ca = np.random.randint(1, p)
                sc = eval_flip_delta_nnz(U, V, W, lw, cr, cs, cm, co, ca, N1, N2, N3, p, inv)
                if sc < best_score:
                    best_score = sc
                    best_r = cr; best_s = cs; best_m = cm; best_o = co; best_a = ca
                    if sc <= -1000:
                        break
            r = best_r; s = best_s; mode = best_m; orient = best_o; alpha = best_a
        else:
            pick = np.random.randint(0, n_pairs)
            r = pr[pick]; s = ps[pick]; mode = pm[pick]
            orient = np.random.randint(0, 2)
            alpha = np.random.randint(1, p)

        do_flip(U, V, W, ku, kv, kw, lw, r, s, mode, orient, alpha, N1, N2, N3, p, inv)
        total_flips += 1

        check_sl = (total_flips & 3) == 0
        R_old = R
        R = reduce_fast_and_slices(U, V, W, ku, kv, kw, lw, R, N1, N2, N3, p, inv, check_sl)
        if R < R_old:
            total_reds += (R_old - R)
            if R < R_best:
                R_best = R
                pool_size = 1
                for rr in range(R):
                    for i in range(N1): pool_U[0, rr, i] = U[rr, i]
                    for i in range(N2): pool_V[0, rr, i] = V[rr, i]
                    for i in range(N3): pool_W[0, rr, i] = W[rr, i]
                pool_sig[0] = scheme_signature(ku, kv, kw, R)
                steps_since_best = 0
            elif R == R_best:
                sig = scheme_signature(ku, kv, kw, R)
                seen = False
                for q in range(pool_size):
                    if pool_sig[q] == sig:
                        seen = True; break
                if not seen:
                    slot = pool_size if pool_size < pool_cap else np.random.randint(1, pool_cap)
                    if pool_size < pool_cap:
                        pool_size += 1
                    for rr in range(R):
                        for i in range(N1): pool_U[slot, rr, i] = U[rr, i]
                        for i in range(N2): pool_V[slot, rr, i] = V[rr, i]
                        for i in range(N3): pool_W[slot, rr, i] = W[rr, i]
                    pool_sig[slot] = sig
                steps_since_best = 0
            else:
                steps_since_best += 1
        else:
            steps_since_best += 1

    # Trouver le schema le plus parcimonieux dans le pool au meilleur rang R_best
    best_slot = 0
    best_nnz = 10**9
    for q in range(pool_size):
        nnz = 0
        for rr in range(R_best):
            for i in range(N1):
                if pool_U[q, rr, i] != 0: nnz += 1
            for i in range(N2):
                if pool_V[q, rr, i] != 0: nnz += 1
            for i in range(N3):
                if pool_W[q, rr, i] != 0: nnz += 1
        if nnz < best_nnz:
            best_nnz = nnz
            best_slot = q

    U_out = pool_U[best_slot, :R_best, :].copy()
    V_out = pool_V[best_slot, :R_best, :].copy()
    W_out = pool_W[best_slot, :R_best, :].copy()
    return R_best, U_out, V_out, W_out, total_flips, total_plus, total_reds, pool_size, best_nnz


def _worker_v2(args):
    m, k, n, p, max_steps, burst, pdepth, seed, target_rank = args
    U0, V0, W0 = embed_block_schemes(m, k, n, p)
    R_b, Ub, Vb, Wb, nf, np_, nr, psize, nnz = guided_flip_worker(
        U0, V0, W0, p, max_steps, burst, pdepth, seed, target_rank
    )
    ok = verify_scheme_fp(Ub, Vb, Wb, m, k, n, p)
    return {
        "seed": seed, "R_best": int(R_b), "nnz": int(nnz), "verified": bool(ok),
        "flips": int(nf), "plus": int(np_), "reds": int(nr), "pool_size": int(psize),
        "factors": (Ub, Vb, Wb),
    }


def run_v2_campaign(pool, label, m, k, n, p, n_runs, max_steps, burst, pdepth, target_rank, world_rec, seed_base):
    ts = time.perf_counter()
    tasks = [
        (m, k, n, p, max_steps, burst, pdepth, seed_base + w, target_rank)
        for w in range(n_runs)
    ]
    outs = list(pool.map(_worker_v2, tasks))
    dt = time.perf_counter() - ts
    all_ok = all(o["verified"] for o in outs)
    assert all_ok, f"Verification failed on {label}!"

    outs.sort(key=lambda z: (z["R_best"], z["nnz"]))
    best = outs[0]
    tot_flips = sum(o["flips"] for o in outs)
    tot_plus = sum(o["plus"] for o in outs)
    ranks = [o["R_best"] for o in outs]
    rate_m = (tot_flips / dt) / 1e6

    saved = None
    if best["R_best"] <= world_rec:
        fname = f"flipv2_{label}_R{best['R_best']}.npz"
        saved = os.path.join(HERE, fname)
        Ub, Vb, Wb = best["factors"]
        np.savez(saved, U=Ub, V=Vb, W=Wb, shape=np.array([m, k, n]), p=np.array([p]))

    tag = "*** RECORD BATTU! ***" if best["R_best"] < world_rec else ("== RECORD EGALE ==" if best["R_best"] == world_rec else "")
    print(
        f"   {label:<10} | <{m},{k},{n}> F_{p} | best R={best['R_best']:>2} (rec {world_rec:>2}) | nnz={best['nnz']:>3} | "
        f"pool={best['pool_size']:>2} | rangs={ranks[:6]} | {tot_flips/1e6:>5.1f}M flips ({rate_m:>4.1f}M/s) | {dt:>5.1f}s {tag}",
        flush=True,
    )
    return {
        "label": label, "shape": [m, k, n], "p": p, "world_rec": world_rec, "target_rank": target_rank,
        "best_R": best["R_best"], "best_nnz": best["nnz"], "all_ranks": ranks,
        "total_flips": tot_flips, "total_plus": tot_plus, "time_s": dt,
        "verified_100pct": all_ok, "saved_npz": saved,
    }


def main():
    t0 = time.time()
    print("== Verification exacte des schemas de depart (Strassen & Laderman) sur F_3 et F_5 ==\n", flush=True)
    for (m, k, n) in [(2, 3, 4), (3, 3, 3), (2, 4, 4), (3, 3, 4)]:
        for p in (3, 5):
            U0, V0, W0 = embed_block_schemes(m, k, n, p)
            ok = verify_scheme_fp(U0, V0, W0, m, k, n, p)
            nnz0 = int(np.count_nonzero(U0) + np.count_nonzero(V0) + np.count_nonzero(W0))
            assert ok, f"Bad warm-start for {(m,k,n)} F_{p}"
            if p == 3:
                print(f"   Warm-start <{m},{k},{n}> : R0={U0.shape[0]:>2}, nnz={nnz0:>3} | verifie exact F_3 & F_5 = True", flush=True)

    # Warmup JIT rapide
    U0, V0, W0 = embed_block_schemes(2, 3, 4, 3)
    guided_flip_worker(U0, V0, W0, 3, 2000, 300, 1, 1, 19)

    n_workers = min(10, os.cpu_count() or 4)
    print(f"\n== Phase 6 : Flip-Graph GUIDE PAR PARCIMONIE + Warm-Start ({n_workers} coeurs Numba JIT) ==\n", flush=True)

    campaigns = [
        # (label, m, k, n, p, n_runs, max_steps, burst, pdepth, target_rank, world_rec, seed_base)
        ("2x3x4_F3", 2, 3, 4, 3, 10, 2_000_000, 1200, 2, 19, 20, 10000),
        ("2x3x4_F5", 2, 3, 4, 5, 10, 2_000_000, 1200, 2, 19, 20, 11000),
        ("3x3x3_F3", 3, 3, 3, 3, 10, 4_000_000, 1600, 2, 22, 23, 20000),
        ("3x3x3_F5", 3, 3, 3, 5, 10, 4_000_000, 1600, 2, 22, 23, 21000),
        ("2x4x4_F3", 2, 4, 4, 3, 10, 2_500_000, 1500, 2, 25, 26, 30000),
        ("2x4x4_F5", 2, 4, 4, 5, 10, 2_500_000, 1500, 2, 25, 26, 31000),
        ("3x3x4_F3", 3, 3, 4, 3, 10, 3_000_000, 1800, 2, 28, 29, 40000),
        ("3x3x4_F5", 3, 3, 4, 5, 10, 3_000_000, 1800, 2, 28, 29, 41000),
    ]

    report = {"n_workers": n_workers, "campaigns": {}}
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as pool:
        for cfg in campaigns:
            res = run_v2_campaign(pool, *cfg)
            report["campaigns"][cfg[0]] = res

    out_json = os.path.join(HERE, "flip_graph_v2_results.json")
    with open(out_json, "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14.0, 5.0))

    fmts = ["2x3x4", "3x3x3", "2x4x4", "3x3x4"]
    starts = [22, 23, 28, 32]
    world = [20, 23, 26, 29]
    bf3 = [report["campaigns"][f"{f}_F3"]["best_R"] for f in fmts]
    bf5 = [report["campaigns"][f"{f}_F5"]["best_R"] for f in fmts]
    x = np.arange(len(fmts))
    w = 0.2

    ax1.bar(x - 1.5 * w, starts, width=w, color="#cccccc", label="Warm-start exact (blocs)")
    ax1.bar(x - 0.5 * w, world, width=w, color="#2ca02c", label="Record mondial exact")
    ax1.bar(x + 0.5 * w, bf3, width=w, color="#1f77b4", label="Flip-Graph guidé F_3")
    ax1.bar(x + 1.5 * w, bf5, width=w, color="#d62728", label="Flip-Graph guidé F_5")
    ax1.set_xticks(x)
    ax1.set_xticklabels([f"<{f.replace('x',',')}>" for f in fmts])
    ax1.set_ylabel("Rang exact (multiplications scalaires)")
    ax1.set_title("Flip-Graph guidé par parcimonie + Warm-Start (F_3 & F_5)")
    ax1.grid(axis="y", alpha=0.3)
    ax1.legend(fontsize=8)

    # Panel 2 : Comparaison Phase 5 (non-guide) vs Phase 6 (guide par parcimonie) sur F_3 et F_5
    with open(os.path.join(HERE, "flip_graph_fp_results.json"), "r", encoding="utf-8") as fh:
        p5 = json.load(fh)["campaigns"]
    p5_f3 = [p5[f"{f}_F3"]["best_R"] - r for f, r in zip(fmts, world)]
    p5_f5 = [p5[f"{f}_F5"]["best_R"] - r for f, r in zip(fmts, world)]
    p6_f3 = [b - r for b, r in zip(bf3, world)]
    p6_f5 = [b - r for b, r in zip(bf5, world)]

    ax2.plot(x, p5_f3, "o:", color="#1f77b4", alpha=0.45, markersize=7, label="Phase 5 non-guidé (F_3)")
    ax2.plot(x, p5_f5, "s:", color="#d62728", alpha=0.45, markersize=7, label="Phase 5 non-guidé (F_5)")
    ax2.plot(x, p6_f3, "o-", color="#1f77b4", linewidth=2.2, markersize=9, label="Phase 6 guidé parcimonie (F_3)")
    ax2.plot(x, p6_f5, "s-", color="#d62728", linewidth=2.2, markersize=9, label="Phase 6 guidé parcimonie (F_5)")
    ax2.axhline(0, color="#2ca02c", linestyle="--", linewidth=1.8, label="Record mondial égalé (0)")
    ax2.axhline(-1, color="#9467bd", linestyle="--", linewidth=1.4, label="Nouveau record (-1)")
    ax2.set_xticks(x)
    ax2.set_xticklabels([f"<{f.replace('x',',')}>" for f in fmts])
    ax2.set_ylabel("Meilleur rang trouvé − Record mondial")
    ax2.set_title("Impact du guidage par parcimonie (F_3 & F_5)")
    ax2.grid(alpha=0.3)
    ax2.legend(fontsize=7.8)

    fig.tight_layout()
    out_png = os.path.join(OUT, "flip_graph_v2_guided.png")
    fig.savefig(out_png, dpi=130)
    print(f"\nTemps total : {time.time() - t0:.1f} s", flush=True)
    print(f"Graphique : {out_png}", flush=True)


if __name__ == "__main__":
    main()
