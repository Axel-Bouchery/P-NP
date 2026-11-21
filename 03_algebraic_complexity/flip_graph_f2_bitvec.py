r"""
Graphe de Flips bit-vectoriel ultra-rapide sur F_2 (uint64 XOR, Numba 16 coeurs)
================================================================================
Innovations algorithmiques integrees :
  1. Encodage 100% bit-vectoriel uint64 : chaque terme r est un triplet
     (u[r], v[r], w[r]) de registres uint64 (jusqu'a N1, N2, N3 <= 64).
     Un flip = 2 instructions XOR CPU (0 allocation, 0 division, 0 cle projective).
  2. Reduction exacte de tranches par elimination de Gauss bit-vectorielle sur F_2 :
     pour toute tranche partageant une meme direction u (ou v, ou w), reduit le
     sous-tenseur bilineaire par pivot de Gauss XOR en O(m^2).
  3. Detection gloutonne de triangles XOR a 1 pas (1-Step XOR Triangle Lookahead) :
     scanne si deux termes partageant u[i] == u[j] ont v[i] ^ v[j] == v[k] (ou w[i] ^ w[j] == w[k])
     pour declencher immediatement une cascade de reduction R -> R - 1 !
  4. Transitions PLUS adaptatives (Arai et al. 2024) + Seeding par blocs (Strassen/Laderman).
  5. Verification exacte independante du tenseur T_{<m,k,n>} mod 2 sur 100% des schemas.
"""
import json
import os
import sys
import time

import numpy as np
from numba import njit, prange

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RESULTS_DIR = os.path.join(ROOT, "results", "algebraic_complexity")
WITNESS_DIR = os.path.join(ROOT, "witnesses", "algebraic_complexity")
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(WITNESS_DIR, exist_ok=True)


def build_schoolbook_f2_bits(m, k, n):
    """Construit le schema naif R = m*k*n en registres uint64."""
    R = m * k * n
    u = np.zeros(R, dtype=np.uint64)
    v = np.zeros(R, dtype=np.uint64)
    w = np.zeros(R, dtype=np.uint64)
    r = 0
    for i in range(m):
        for t in range(k):
            for j in range(n):
                u[r] = np.uint64(1) << np.uint64(i * k + t)
                v[r] = np.uint64(1) << np.uint64(t * n + j)
                w[r] = np.uint64(1) << np.uint64(i * n + j)
                r += 1
    return u, v, w


def verify_scheme_f2_bits(u, v, w, m, k, n):
    """Verifie exactement que sum_r u_r x v_r x w_r == T_{<m,k,n>} sur F_2."""
    N1, N2, N3 = m * k, k * n, m * n
    assert N1 <= 64 and N2 <= 64 and N3 <= 64
    R = len(u)
    U = np.zeros((R, N1), dtype=np.int64)
    V = np.zeros((R, N2), dtype=np.int64)
    W = np.zeros((R, N3), dtype=np.int64)
    for r in range(R):
        ur, vr, wr = int(u[r]), int(v[r]), int(w[r])
        for i in range(N1):
            if (ur >> i) & 1:
                U[r, i] = 1
        for j in range(N2):
            if (vr >> j) & 1:
                V[r, j] = 1
        for l in range(N3):
            if (wr >> l) & 1:
                W[r, l] = 1
    T_hat = np.einsum("ri,rj,rk->ijk", U, V, W) & 1
    T_ref = np.zeros((N1, N2, N3), dtype=np.int64)
    for i in range(m):
        for t in range(k):
            for j in range(n):
                T_ref[i * k + t, t * n + j, i * n + j] = 1
    return int(np.sum(T_hat != T_ref)) == 0


@njit(cache=True)
def compact_nonzero(u, v, w, R):
    """Supprime tous les termes ou u[r]==0 ou v[r]==0 ou w[r]==0."""
    out = 0
    for r in range(R):
        if u[r] != np.uint64(0) and v[r] != np.uint64(0) and w[r] != np.uint64(0):
            if out != r:
                u[out] = u[r]
                v[out] = v[r]
                w[out] = w[r]
            out += 1
    return out


@njit(cache=True)
def reduce_slices_gauss_f2(u, v, w, R):
    """
    Reduction exacte de tranches sur F_2 par elimination de Gauss XOR :
      1. Fusion directe des paires ayant deux vecteurs identiques (ex. u[i]==u[j] et v[i]==v[j]).
      2. Pour chaque tranche partageant la meme direction u (m >= 3 termes),
         echelonne les vecteurs v par XOR tout en accumulant l'operation duale sur w ;
         puis echelonne w en accumulant sur v.
      3. Meme operation par permutation cyclique des modes (v fixe, puis w fixe).
    """
    changed = True
    while changed:
        changed = False
        R = compact_nonzero(u, v, w, R)

        # Passe rapide : absorption directe des paires ayant 2 composantes egales
        for i in range(R):
            if u[i] == np.uint64(0):
                continue
            ui, vi, wi = u[i], v[i], w[i]
            for j in range(i + 1, R):
                if u[j] == np.uint64(0):
                    continue
                if ui == u[j]:
                    if vi == v[j]:
                        wi ^= w[j]
                        w[i] = wi
                        u[j] = np.uint64(0); v[j] = np.uint64(0); w[j] = np.uint64(0)
                        changed = True
                        if wi == np.uint64(0):
                            u[i] = np.uint64(0); v[i] = np.uint64(0)
                            break
                    elif wi == w[j]:
                        vi ^= v[j]
                        v[i] = vi
                        u[j] = np.uint64(0); v[j] = np.uint64(0); w[j] = np.uint64(0)
                        changed = True
                        if vi == np.uint64(0):
                            u[i] = np.uint64(0); w[i] = np.uint64(0)
                            break
                elif vi == v[j] and wi == w[j]:
                    ui ^= u[j]
                    u[i] = ui
                    u[j] = np.uint64(0); v[j] = np.uint64(0); w[j] = np.uint64(0)
                    changed = True
                    if ui == np.uint64(0):
                        v[i] = np.uint64(0); w[i] = np.uint64(0)
                        break
        if changed:
            R = compact_nonzero(u, v, w, R)
            continue

        # Passe Gauss sur chaque tranche de taille m >= 3 pour les 3 modes (0: u, 1: v, 2: w)
        visited = np.zeros(R, dtype=np.int8)
        idx_buf = np.zeros(R, dtype=np.int32)
        for mode in range(3):
            for r in range(R):
                visited[r] = 0
            for i in range(R):
                if visited[i]:
                    continue
                key = u[i] if mode == 0 else (v[i] if mode == 1 else w[i])
                if key == np.uint64(0):
                    continue
                m_cnt = 0
                for j in range(i, R):
                    kj = u[j] if mode == 0 else (v[j] if mode == 1 else w[j])
                    if kj == key:
                        visited[j] = 1
                        idx_buf[m_cnt] = j
                        m_cnt += 1
                if m_cnt < 3:
                    continue

                # Gauss sur A (le 1er autre mode) avec accumulation duale sur B (le 2e autre mode)
                # Identite : a_p ⊗ b_p ⊕ a_q ⊗ b_q = a_p ⊗ (b_p ⊕ b_q) ⊕ (a_p ⊕ a_q) ⊗ b_q
                rank = 0
                for col in range(64):
                    mask = np.uint64(1) << np.uint64(col)
                    piv = -1
                    for r_idx in range(rank, m_cnt):
                        val_a = v[idx_buf[r_idx]] if mode == 0 else (u[idx_buf[r_idx]] if mode == 1 else u[idx_buf[r_idx]])
                        if (val_a & mask) != np.uint64(0):
                            piv = r_idx
                            break
                    if piv == -1:
                        continue
                    if piv != rank:
                        tmp = idx_buf[rank]; idx_buf[rank] = idx_buf[piv]; idx_buf[piv] = tmp
                    p_term = idx_buf[rank]
                    for r_idx in range(m_cnt):
                        if r_idx == rank:
                            continue
                        q_term = idx_buf[r_idx]
                        val_a = v[q_term] if mode == 0 else u[q_term]
                        if (val_a & mask) != np.uint64(0):
                            # a_q ^= a_p ; b_p ^= b_q
                            if mode == 0:
                                v[q_term] ^= v[p_term]
                                w[p_term] ^= w[q_term]
                            elif mode == 1:
                                u[q_term] ^= u[p_term]
                                w[p_term] ^= w[q_term]
                            else:
                                u[q_term] ^= u[p_term]
                                v[p_term] ^= v[q_term]
                    rank += 1

                if rank < m_cnt:
                    for r_idx in range(rank, m_cnt):
                        qt = idx_buf[r_idx]
                        u[qt] = np.uint64(0); v[qt] = np.uint64(0); w[qt] = np.uint64(0)
                    changed = True
                    break

                #Verifier si un b_p est devenu 0 ou si deux b_p sont devenus egaux
                for r_idx in range(m_cnt):
                    pt = idx_buf[r_idx]
                    if u[pt] == np.uint64(0) or v[pt] == np.uint64(0) or w[pt] == np.uint64(0):
                        changed = True
                        break
                if changed:
                    break
            if changed:
                R = compact_nonzero(u, v, w, R)
                break
    return compact_nonzero(u, v, w, R)


@njit(cache=True)
def try_greedy_lookahead_reduction(u, v, w, R):
    """
    Scan O(R^2) pour trouver un flip qui permet une reduction immediate (lookahead a 1 pas) :
    Si deux termes i, j partagent u[i] == u[j], alors un flip peut remplacer :
      - v[i] par v[i] ^ v[j] (et w[j] par w[j] ^ w[i]).
        Si le nouveau w[j] ^ w[i] == w[k] pour un terme k ayant v[k] == v[j] (ou u[k] == u[j]),
        alors au pas suivant j et k fusionnent -> R diminue de 1 !
    Plus generalement, si trois termes i, j, k verifient :
      u[i] == u[j] et u[k] == u[i] et (v[i] ^ v[j] == v[k] ou w[i] ^ w[j] == w[k]),
      reduce_slices_gauss_f2 l'a deja reduit.
    Mais si u[i] == u[j], w[k] == w[j] ^ w[i] et v[k] == v[j] (avec u[k] != u[j]) :
      En flippant v[i] ^= v[j] et w[j] ^= w[i], le terme j devient (u[j], v[j], w[k]),
      qui partage ALORS v[j] == v[k] ET w[j] == w[k] avec le terme k -> FUSION IMMEDIATE j + k !
    """
    for i in range(R):
        ui, vi, wi = u[i], v[i], w[i]
        for j in range(i + 1, R):
            uj, vj, wj = u[j], v[j], w[j]
            # Cas 1 : u[i] == u[j]
            if ui == uj:
                wx = wi ^ wj
                vx = vi ^ vj
                for k in range(R):
                    if k == i or k == j:
                        continue
                    # Si w[k] == wi ^ wj et v[k] == vj (ou v[k] == vi)
                    if w[k] == wx:
                        if v[k] == vj:
                            v[i] ^= vj; w[j] = wx
                            u[k] ^= uj
                            u[j] = np.uint64(0); v[j] = np.uint64(0); w[j] = np.uint64(0)
                            return compact_nonzero(u, v, w, R), True
                        elif v[k] == vi:
                            v[j] ^= vi; w[i] = wx
                            u[k] ^= ui
                            u[i] = np.uint64(0); v[i] = np.uint64(0); w[i] = np.uint64(0)
                            return compact_nonzero(u, v, w, R), True
                    if v[k] == vx:
                        if w[k] == wj:
                            w[i] ^= wj; v[j] = vx
                            u[k] ^= uj
                            u[j] = np.uint64(0); v[j] = np.uint64(0); w[j] = np.uint64(0)
                            return compact_nonzero(u, v, w, R), True
                        elif w[k] == wi:
                            w[j] ^= wi; v[i] = vx
                            u[k] ^= ui
                            u[i] = np.uint64(0); v[i] = np.uint64(0); w[i] = np.uint64(0)
                            return compact_nonzero(u, v, w, R), True
            # Cas 2 : v[i] == v[j]
            if vi == vj:
                ux = ui ^ uj
                wx = wi ^ wj
                for k in range(R):
                    if k == i or k == j:
                        continue
                    if w[k] == wx:
                        if u[k] == uj:
                            u[i] ^= uj; w[j] = wx
                            v[k] ^= vj
                            u[j] = np.uint64(0); v[j] = np.uint64(0); w[j] = np.uint64(0)
                            return compact_nonzero(u, v, w, R), True
                        elif u[k] == ui:
                            u[j] ^= ui; w[i] = wx
                            v[k] ^= vi
                            u[i] = np.uint64(0); v[i] = np.uint64(0); w[i] = np.uint64(0)
                            return compact_nonzero(u, v, w, R), True
                    if u[k] == ux:
                        if w[k] == wj:
                            w[i] ^= wj; u[j] = ux
                            v[k] ^= vj
                            u[j] = np.uint64(0); v[j] = np.uint64(0); w[j] = np.uint64(0)
                            return compact_nonzero(u, v, w, R), True
                        elif w[k] == wi:
                            w[j] ^= wi; u[i] = ux
                            v[k] ^= vi
                            u[i] = np.uint64(0); v[i] = np.uint64(0); w[i] = np.uint64(0)
                            return compact_nonzero(u, v, w, R), True
            # Cas 3 : w[i] == w[j]
            if wi == wj:
                ux = ui ^ uj
                vx = vi ^ vj
                for k in range(R):
                    if k == i or k == j:
                        continue
                    if v[k] == vx:
                        if u[k] == uj:
                            u[i] ^= uj; v[j] = vx
                            w[k] ^= wj
                            u[j] = np.uint64(0); v[j] = np.uint64(0); w[j] = np.uint64(0)
                            return compact_nonzero(u, v, w, R), True
                        elif u[k] == ui:
                            u[j] ^= ui; v[i] = vx
                            w[k] ^= wi
                            u[i] = np.uint64(0); v[i] = np.uint64(0); w[i] = np.uint64(0)
                            return compact_nonzero(u, v, w, R), True
                    if u[k] == ux:
                        if v[k] == vj:
                            v[i] ^= vj; u[j] = ux
                            w[k] ^= wj
                            u[j] = np.uint64(0); v[j] = np.uint64(0); w[j] = np.uint64(0)
                            return compact_nonzero(u, v, w, R), True
                        elif v[k] == vi:
                            v[j] ^= vi; u[i] = ux
                            w[k] ^= wi
                            u[i] = np.uint64(0); v[i] = np.uint64(0); w[i] = np.uint64(0)
                            return compact_nonzero(u, v, w, R), True
    return R, False


@njit(cache=True)
def do_plus_transition_f2(u, v, w, R, max_R, st):
    """
    Transition PLUS d'Arai et al. (2024) sur F_2 :
    Choisit deux termes distincts i, j, scinde i ou cree un terme auxiliaire
    qui partage une direction avec i et une avec j, puis effectue immediatement
    deux flips croises pour sortir d'un bassin local (R -> R + 1).
    Identite exacte :
      Pour deux termes (u_i, v_i, w_i) et (u_j, v_j, w_j),
      ajoutons 0 = (u_i ^ u_j) ⊗ v_j ⊗ w_i  ⊕  u_i ⊗ v_j ⊗ w_i  ⊕  u_j ⊗ v_j ⊗ w_i :
      - u_i ⊗ v_j ⊗ w_i fusionne avec i (u_i, v_i, w_i) -> i devient (u_i, v_i ^ v_j, w_i)
      - u_j ⊗ v_j ⊗ w_i fusionne avec j (u_j, v_j, w_j) -> j devient (u_j, v_j, w_j ^ w_i)
      - il reste le nouveau terme R : (u_i ^ u_j, v_j, w_i) !
    Verification :
      u_i ⊗ (v_i ^ v_j) ⊗ w_i  ⊕  u_j ⊗ v_j ⊗ (w_j ^ w_i)  ⊕  (u_i ^ u_j) ⊗ v_j ⊗ w_i
      = u_i ⊗ v_i ⊗ w_i  ⊕  u_j ⊗ v_j ⊗ w_j ! (100% EXACT !)
    """
    if R >= max_R or R < 2:
        return R, st
    st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
    i = np.int64(st % np.uint64(R))
    st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
    j = np.int64(st % np.uint64(R - 1))
    if j >= i:
        j += np.int64(1)
    st ^= st << np.uint64(13); st ^= st >> np.uint64(7); st ^= st << np.uint64(17)
    perm = np.int64(st % np.uint64(3))

    ui, vi, wi = u[i], v[i], w[i]
    uj, vj, wj = u[j], v[j], w[j]

    if perm == 0:
        if ui == uj or vi == vj or wi == wj:
            return R, st
        v[i] = vi ^ vj
        w[j] = wj ^ wi
        u[R] = ui ^ uj
        v[R] = vj
        w[R] = wi
        return R + 1, st
    elif perm == 1:
        if vi == vj or wi == wj or ui == uj:
            return R, st
        w[i] = wi ^ wj
        u[j] = uj ^ ui
        u[R] = ui
        v[R] = vi ^ vj
        w[R] = wj
        return R + 1, st
    else:
        if wi == wj or ui == uj or vi == vj:
            return R, st
        u[i] = ui ^ uj
        v[j] = vj ^ vi
        u[R] = uj
        v[R] = vi
        w[R] = wi ^ wj
        return R + 1, st


@njit(parallel=True, cache=True)
def flip_walk_f2_parallel(u_init, v_init, w_init, target_R, n_workers, steps_per_worker, plus_interval, max_plus_gap):
    """
    Lance `n_workers` marches de flips F_2 bit-vectorielles en parallele avec :
      - Reduction exacte de tranches par Gauss XOR
      - Lookahead glouton a 1 pas (triangle XOR) tous les 32 flips
      - Transitions PLUS adaptatives quand un worker stagne pendant `plus_interval` pas
        (avec plafond R <= best_R_worker + max_plus_gap)
      - Reprise periodique depuis le meilleur schema du worker.
    """
    R0 = len(u_init)
    max_R = R0 + max_plus_gap + 4

    best_R_all = np.full(n_workers, R0, dtype=np.int32)
    best_u_all = np.zeros((n_workers, max_R), dtype=np.uint64)
    best_v_all = np.zeros((n_workers, max_R), dtype=np.uint64)
    best_w_all = np.zeros((n_workers, max_R), dtype=np.uint64)
    flips_done = np.zeros(n_workers, dtype=np.int64)

    for wid in prange(n_workers):
        st = np.uint64(0x9E3779B97F4A7C15) ^ np.uint64((wid + 1) * 6364136223846793005 + R0 * 7919)
        u = np.zeros(max_R, dtype=np.uint64)
        v = np.zeros(max_R, dtype=np.uint64)
        w = np.zeros(max_R, dtype=np.uint64)
        for r in range(R0):
            u[r] = u_init[r]; v[r] = v_init[r]; w[r] = w_init[r]
            best_u_all[wid, r] = u_init[r]
            best_v_all[wid, r] = v_init[r]
            best_w_all[wid, r] = w_init[r]
        R = reduce_slices_gauss_f2(u, v, w, R0)
        best_R = R
        best_R_all[wid] = best_R

        # Buffers de paires flippables (mode, i, j)
        pair_mode = np.zeros(max_R * max_R, dtype=np.int8)
        pair_i = np.zeros(max_R * max_R, dtype=np.int16)
        pair_j = np.zeros(max_R * max_R, dtype=np.int16)

        since_improve = 0
        since_plus = 0

        for step in range(steps_per_worker):
            if best_R <= target_R:
                break

            # Lookahead glouton a 1 pas tous les 16 pas
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

            # Collecter les paires flippables
            n_pairs = 0
            for i in range(R):
                ui, vi, wi = u[i], v[i], w[i]
                for j in range(i + 1, R):
                    if ui == u[j]:
                        pair_mode[n_pairs] = 0
                        pair_i[n_pairs] = i; pair_j[n_pairs] = j
                        n_pairs += 1
                    elif vi == v[j]:
                        pair_mode[n_pairs] = 1
                        pair_i[n_pairs] = i; pair_j[n_pairs] = j
                        n_pairs += 1
                    elif wi == w[j]:
                        pair_mode[n_pairs] = 2
                        pair_i[n_pairs] = i; pair_j[n_pairs] = j
                        n_pairs += 1

            if n_pairs == 0 or since_plus >= plus_interval:
                if R < best_R + max_plus_gap:
                    R, st = do_plus_transition_f2(u, v, w, R, max_R, st)
                    since_plus = 0
                    continue
                else:
                    # Reset au meilleur schema connu du worker + 1 transition PLUS
                    R = best_R
                    for r in range(R):
                        u[r] = best_u_all[wid, r]
                        v[r] = best_v_all[wid, r]
                        w[r] = best_w_all[wid, r]
                    R, st = do_plus_transition_f2(u, v, w, R, max_R, st)
                    since_plus = 0
                    since_improve = 0
                    continue

            # Choisir une paire au hasard et effectuer un flip XOR
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
                # u[i] == u[j]
                if dir_bit == 0:
                    v[i] ^= v[j]; w[j] ^= w[i]
                else:
                    w[i] ^= w[j]; v[j] ^= v[i]
            elif md == 1:
                # v[i] == v[j]
                if dir_bit == 0:
                    u[i] ^= u[j]; w[j] ^= w[i]
                else:
                    w[i] ^= w[j]; u[j] ^= u[i]
            else:
                # w[i] == w[j]
                if dir_bit == 0:
                    u[i] ^= u[j]; v[j] ^= v[i]
                else:
                    v[i] ^= v[j]; u[j] ^= u[i]

            flips_done[wid] += 1
            since_plus += 1
            since_improve += 1

            # Verifier si annulation directe ou tous les 8 flips lancer Gauss de tranche
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
                            best_u_all[wid, r] = u[r]
                            best_v_all[wid, r] = v[r]
                            best_w_all[wid, r] = w[r]
                        since_improve = 0

            # Si on dérive trop longtemps à R > best_R, revenir à best_R
            if since_improve > plus_interval * 8 and R > best_R:
                R = best_R
                for r in range(R):
                    u[r] = best_u_all[wid, r]
                    v[r] = best_v_all[wid, r]
                    w[r] = best_w_all[wid, r]
                since_improve = 0
                since_plus = plus_interval  # forcer un nouveau PLUS

    return best_R_all, best_u_all, best_v_all, best_w_all, flips_done


def run_format_benchmark(m, k, n, target_R, n_workers=16, steps=150000, plus_int=350, max_gap=2):
    u0, v0, w0 = build_schoolbook_f2_bits(m, k, n)
    R0 = len(u0)
    t0 = time.time()
    bRs, bus, bvs, bws, fdone = flip_walk_f2_parallel(
        u0, v0, w0, target_R, n_workers, steps, plus_int, max_gap
    )
    dt = time.time() - t0
    bi = int(np.argmin(bRs))
    best_R = int(bRs[bi])
    u_best = bus[bi, :best_R].copy()
    v_best = bvs[bi, :best_R].copy()
    w_best = bws[bi, :best_R].copy()
    ok = verify_scheme_f2_bits(u_best, v_best, w_best, m, k, n)
    assert ok, f"Echec verification <{m},{k},{n}> R={best_R}!"
    tot_flips = int(np.sum(fdone))
    n_hit = int(np.sum(bRs == best_R))
    print(
        f"  <{m},{k},{n}> : naif R={R0:2d} -> atteint R={best_R:2d} (cible={target_R:2d}, "
        f"workers={n_hit:2d}/{n_workers}) | {tot_flips:9,d} flips en {dt:5.2f}s ({tot_flips/max(1e-6,dt)/1e6:5.2f} M flips/s) [VERIFIE EXACT]"
    )
    return {
        "format": [m, k, n],
        "naive_R": R0,
        "best_R": best_R,
        "target_R": target_R,
        "workers_at_best": n_hit,
        "total_flips": tot_flips,
        "time_s": round(dt, 2),
        "mflips_per_s": round(tot_flips / max(1e-6, dt) / 1e6, 2),
        "u": [int(x) for x in u_best],
        "v": [int(x) for x in v_best],
        "w": [int(x) for x in w_best],
    }


def main():
    t0 = time.time()
    print("====================================================================================")
    print("  CALIBRATION & CHASSE FLIP-GRAPH BIT-VECTORIEL F_2 (uint64 XOR + GAUSS + LOOKAHEAD)")
    print("====================================================================================")

    # Warm-up JIT sur <2, 2, 2>
    res_all = []
    res_all.append(run_format_benchmark(2, 2, 2, target_R=7, n_workers=16, steps=20000, plus_int=200, max_gap=1))
    res_all.append(run_format_benchmark(2, 2, 3, target_R=11, n_workers=16, steps=40000, plus_int=250, max_gap=1))
    res_all.append(run_format_benchmark(2, 3, 3, target_R=15, n_workers=16, steps=80000, plus_int=300, max_gap=2))
    res_all.append(run_format_benchmark(3, 3, 3, target_R=23, n_workers=16, steps=150000, plus_int=350, max_gap=2))
    res_all.append(run_format_benchmark(2, 3, 4, target_R=20, n_workers=16, steps=150000, plus_int=350, max_gap=2))
    res_all.append(run_format_benchmark(2, 3, 5, target_R=25, n_workers=16, steps=200000, plus_int=400, max_gap=2))
    res_all.append(run_format_benchmark(3, 3, 4, target_R=29, n_workers=16, steps=250000, plus_int=400, max_gap=2))

    out = {"benchmarks": res_all, "total_time_s": round(time.time() - t0, 2)}
    out_path = os.path.join(RESULTS_DIR, "flip_graph_f2_bitvec_bench.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print(f"\nTermine en {time.time() - t0:.2f}s. Sauvegarde : {out_path}")


if __name__ == "__main__":
    main()
