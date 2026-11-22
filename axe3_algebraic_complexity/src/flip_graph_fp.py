"""
Piste 3 (Phase 5) : Graphe de Flips combinatoire EXACT sur F_3 et F_5 (Numba JIT)
=================================================================================

Moteur Kauers-Moosbauer (2023) + transitions PLUS (Arai et al. 2024 /
Moosbauer-Poole 2025) + REDUCTION EXACTE DE TRANCHES PAR DEPENDANCE LINEAIRE
(elimine le facteur combinatoire sur F_3 et F_5 : des que m >= 2 termes partagent
une meme direction u, v ou w, toute dependance lineaire parmi les deux autres
modes est reduite instantanement en O(m^2 * N) !).

Verification exacte independante sur 100% des schemas retournes.
"""
from pathlib import Path
import concurrent.futures
import json
import os
import sys
import time

import numpy as np
from numba import njit

HERE = os.path.dirname(os.path.abspath(__file__))
AXE_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = str(AXE_ROOT / "results")
FIGURES_DIR = str(AXE_ROOT / "figures")
WITNESS_DIR = str(AXE_ROOT / "witnesses")
for _d in (RESULTS_DIR, FIGURES_DIR, WITNESS_DIR):
    os.makedirs(_d, exist_ok=True)
OUT = sys.argv[1] if len(sys.argv) > 1 else FIGURES_DIR


def build_schoolbook_fp(m, k, n, p):
    N1, N2, N3 = m * k, k * n, m * n
    R = m * k * n
    U = np.zeros((R, N1), dtype=np.int16)
    V = np.zeros((R, N2), dtype=np.int16)
    W = np.zeros((R, N3), dtype=np.int16)
    r = 0
    for i in range(m):
        for t in range(k):
            for j in range(n):
                U[r, i * k + t] = 1
                V[r, t * n + j] = 1
                W[r, i * n + j] = 1
                r += 1
    return U, V, W


def verify_scheme_fp(U, V, W, m, k, n, p):
    N1, N2, N3 = m * k, k * n, m * n
    T_hat = np.einsum("ri,rj,rk->ijk", U.astype(np.int64), V.astype(np.int64), W.astype(np.int64)) % p
    T_ref = np.zeros((N1, N2, N3), dtype=np.int64)
    for i in range(m):
        for t in range(k):
            for j in range(n):
                T_ref[i * k + t, t * n + j, i * n + j] = 1
    return int(np.sum(T_hat != T_ref)) == 0


@njit(cache=True)
def make_inv_table(p):
    inv = np.zeros(p, dtype=np.int16)
    for a in range(1, p):
        for b in range(1, p):
            if (a * b) % p == 1:
                inv[a] = b
                break
    return inv


@njit(cache=True)
def canon_and_key_row(vec, N, p, inv):
    lead = 0
    for i in range(N):
        v = vec[i] % p
        vec[i] = v
        if lead == 0 and v != 0:
            lead = v
    if lead == 0:
        return np.int64(0), np.int16(0)
    il = inv[lead]
    key = np.int64(0)
    for i in range(N):
        val = (vec[i] * il) % p
        vec[i] = val
        key = key * p + val
    return key, np.int16(lead)


@njit(cache=True)
def proj_key_readonly(vec, N, p, inv):
    lead = 0
    for i in range(N):
        v = vec[i] % p
        vec[i] = v
        if lead == 0 and v != 0:
            lead = v
    if lead == 0:
        return np.int64(0), np.int16(0)
    il = inv[lead]
    key = np.int64(0)
    for i in range(N):
        val = (vec[i] * il) % p
        key = key * p + val
    return key, np.int16(lead)


@njit(cache=True)
def update_term_keys(U, V, W, ku, kv, kw, lw, r, N1, N2, N3, p, inv):
    k1, l1 = canon_and_key_row(U[r], N1, p, inv)
    if l1 == 0:
        ku[r] = 0; kv[r] = 0; kw[r] = 0; lw[r] = 0
        return True
    k2, l2 = canon_and_key_row(V[r], N2, p, inv)
    if l2 == 0:
        ku[r] = 0; kv[r] = 0; kw[r] = 0; lw[r] = 0
        return True
    scale = (l1 * l2) % p
    if scale != 1:
        for i in range(N3):
            W[r, i] = (W[r, i] * scale) % p
    k3, l3 = proj_key_readonly(W[r], N3, p, inv)
    if l3 == 0:
        ku[r] = 0; kv[r] = 0; kw[r] = 0; lw[r] = 0
        return True
    ku[r] = k1
    kv[r] = k2
    kw[r] = k3
    lw[r] = l3
    return False


@njit(cache=True)
def delete_term(U, V, W, ku, kv, kw, lw, r, R, N1, N2, N3):
    last = R - 1
    if r != last:
        for i in range(N1): U[r, i] = U[last, i]
        for i in range(N2): V[r, i] = V[last, i]
        for i in range(N3): W[r, i] = W[last, i]
        ku[r] = ku[last]; kv[r] = kv[last]; kw[r] = kw[last]; lw[r] = lw[last]
    return last


@njit(cache=True)
def find_linear_dep(M, m, N, p, inv, coeffs_out):
    """
    Teste si les m lignes de M (m x N) sur F_p sont lineairement dependantes.
    Si oui, renvoie l'indice d'une ligne dependante `dep_row` et remplit `coeffs_out[:m]`
    tel que M[dep_row] == sum_{j != dep_row} coeffs_out[j] * M[j] (mod p).
    Sinon renvoie -1.
    """
    # Elimination de Gauss sur A = [M | I_m]
    A = np.zeros((m, N + m), dtype=np.int16)
    for i in range(m):
        for c in range(N):
            A[i, c] = M[i, c] % p
        A[i, N + i] = 1

    piv_row = 0
    for col in range(N):
        if piv_row >= m:
            break
        sel = -1
        for r in range(piv_row, m):
            if A[r, col] != 0:
                sel = r
                break
        if sel == -1:
            continue
        if sel != piv_row:
            for c in range(N + m):
                tmp = A[piv_row, c]; A[piv_row, c] = A[sel, c]; A[sel, c] = tmp
        iv = inv[A[piv_row, col]]
        for c in range(col, N + m):
            A[piv_row, c] = (A[piv_row, c] * iv) % p
        for r in range(m):
            if r != piv_row and A[r, col] != 0:
                f = A[r, col]
                for c in range(col, N + m):
                    A[r, c] = (A[r, c] - f * A[piv_row, c]) % p
        piv_row += 1

    if piv_row == m:
        return -1

    # La ligne piv_row a partie gauche nulle -> sum_j A[piv_row, N + j] * M[j] == 0
    dep = -1
    for j in range(m):
        if A[piv_row, N + j] != 0:
            dep = j
            break
    if dep == -1:
        return -1
    # M[dep] = sum_{j != dep} (-A[piv_row, N+j] * inv[A[piv_row, N+dep]]) * M[j]
    inv_d = inv[A[piv_row, N + dep]]
    for j in range(m):
        if j == dep:
            coeffs_out[j] = 0
        else:
            coeffs_out[j] = (-A[piv_row, N + j] * inv_d) % p
    return dep


@njit(cache=True)
def reduce_one_slice(U, V, W, ku, kv, kw, lw, R, mode, target_key, N1, N2, N3, p, inv):
    """
    Collecte tous les termes partageant `target_key` dans `mode` (0=U, 1=V, 2=W).
    Si m >= 2 et qu'un des deux autres modes est lineairement dependant, reduit R -> R-1 !
    Renvoie (nouveau R, reduced_bool).
    """
    idxs = np.empty(R, dtype=np.int16)
    m = 0
    if mode == 0:
        for r in range(R):
            if ku[r] == target_key:
                idxs[m] = r; m += 1
    elif mode == 1:
        for r in range(R):
            if kv[r] == target_key:
                idxs[m] = r; m += 1
    else:
        for r in range(R):
            if kw[r] == target_key:
                idxs[m] = r; m += 1

    if m < 2:
        return R, False

    coeffs = np.zeros(m, dtype=np.int16)

    if mode == 0:
        # Tous ont u_r identique -> tester dependance sur V, puis sur W
        Mv = np.empty((m, N2), dtype=np.int16)
        for j in range(m):
            r = idxs[j]
            for c in range(N2): Mv[j, c] = V[r, c]
        dep = find_linear_dep(Mv, m, N2, p, inv, coeffs)
        if dep != -1:
            r_dep = idxs[dep]
            for j in range(m):
                if j != dep and coeffs[j] != 0:
                    rj = idxs[j]; cj = coeffs[j]
                    for c in range(N3): W[rj, c] = (W[rj, c] + cj * W[r_dep, c]) % p
                    update_term_keys(U, V, W, ku, kv, kw, lw, rj, N1, N2, N3, p, inv)
            R = delete_term(U, V, W, ku, kv, kw, lw, r_dep, R, N1, N2, N3)
            return R, True

        Mw = np.empty((m, N3), dtype=np.int16)
        for j in range(m):
            r = idxs[j]
            for c in range(N3): Mw[j, c] = W[r, c]
        dep = find_linear_dep(Mw, m, N3, p, inv, coeffs)
        if dep != -1:
            r_dep = idxs[dep]
            for j in range(m):
                if j != dep and coeffs[j] != 0:
                    rj = idxs[j]; cj = coeffs[j]
                    for c in range(N2): V[rj, c] = (V[rj, c] + cj * V[r_dep, c]) % p
                    update_term_keys(U, V, W, ku, kv, kw, lw, rj, N1, N2, N3, p, inv)
            R = delete_term(U, V, W, ku, kv, kw, lw, r_dep, R, N1, N2, N3)
            return R, True

    elif mode == 1:
        # Tous ont v_r identique -> tester dependance sur U, puis sur W
        Mu = np.empty((m, N1), dtype=np.int16)
        for j in range(m):
            r = idxs[j]
            for c in range(N1): Mu[j, c] = U[r, c]
        dep = find_linear_dep(Mu, m, N1, p, inv, coeffs)
        if dep != -1:
            r_dep = idxs[dep]
            for j in range(m):
                if j != dep and coeffs[j] != 0:
                    rj = idxs[j]; cj = coeffs[j]
                    for c in range(N3): W[rj, c] = (W[rj, c] + cj * W[r_dep, c]) % p
                    update_term_keys(U, V, W, ku, kv, kw, lw, rj, N1, N2, N3, p, inv)
            R = delete_term(U, V, W, ku, kv, kw, lw, r_dep, R, N1, N2, N3)
            return R, True

        Mw = np.empty((m, N3), dtype=np.int16)
        for j in range(m):
            r = idxs[j]
            for c in range(N3): Mw[j, c] = W[r, c]
        dep = find_linear_dep(Mw, m, N3, p, inv, coeffs)
        if dep != -1:
            r_dep = idxs[dep]
            for j in range(m):
                if j != dep and coeffs[j] != 0:
                    rj = idxs[j]; cj = coeffs[j]
                    for c in range(N1): U[rj, c] = (U[rj, c] + cj * U[r_dep, c]) % p
                    update_term_keys(U, V, W, ku, kv, kw, lw, rj, N1, N2, N3, p, inv)
            R = delete_term(U, V, W, ku, kv, kw, lw, r_dep, R, N1, N2, N3)
            return R, True

    else:
        # Tous ont w_r colineaire : w_j = lw_j * w_canon
        # Somme_j (lw_j * u_j) (x) v_j (x) w_canon
        # Dependance sur (lw_j * u_j) -> absorbee dans v_j ; ou sur (lw_j * v_j) -> absorbee dans u_j
        Mu = np.empty((m, N1), dtype=np.int16)
        for j in range(m):
            r = idxs[j]; sc = lw[r]
            for c in range(N1): Mu[j, c] = (U[r, c] * sc) % p
        dep = find_linear_dep(Mu, m, N1, p, inv, coeffs)
        if dep != -1:
            r_dep = idxs[dep]
            for j in range(m):
                if j != dep and coeffs[j] != 0:
                    rj = idxs[j]; cj = coeffs[j]
                    # (lw_dep u_dep) = sum cj (lw_j u_j) => somme = sum_j (lw_j u_j) (x) (v_j + cj v_dep) (x) w_canon
                    for c in range(N2): V[rj, c] = (V[rj, c] + cj * V[r_dep, c]) % p
                    update_term_keys(U, V, W, ku, kv, kw, lw, rj, N1, N2, N3, p, inv)
            R = delete_term(U, V, W, ku, kv, kw, lw, r_dep, R, N1, N2, N3)
            return R, True

        Mv = np.empty((m, N2), dtype=np.int16)
        for j in range(m):
            r = idxs[j]; sc = lw[r]
            for c in range(N2): Mv[j, c] = (V[r, c] * sc) % p
        dep = find_linear_dep(Mv, m, N2, p, inv, coeffs)
        if dep != -1:
            r_dep = idxs[dep]
            for j in range(m):
                if j != dep and coeffs[j] != 0:
                    rj = idxs[j]; cj = coeffs[j]
                    for c in range(N1): U[rj, c] = (U[rj, c] + cj * U[r_dep, c]) % p
                    update_term_keys(U, V, W, ku, kv, kw, lw, rj, N1, N2, N3, p, inv)
            R = delete_term(U, V, W, ku, kv, kw, lw, r_dep, R, N1, N2, N3)
            return R, True

    return R, False


@njit(cache=True)
def reduce_fast_and_slices(U, V, W, ku, kv, kw, lw, R, N1, N2, N3, p, inv, check_slices):
    """
    1. Supprime tout terme nul et fusionne toute double colinearite (O(R^2) ultra rapide).
    2. Si check_slices=True, verifie aussi toute dependance lineaire multi-termes dans
       chaque tranche colineaire (capture TOUTES les reductions d'ordre m >= 3 sur F_3/F_5 !).
    """
    changed = True
    while changed:
        changed = False
        r = 0
        while r < R:
            if lw[r] == 0:
                R = delete_term(U, V, W, ku, kv, kw, lw, r, R, N1, N2, N3)
                changed = True
            else:
                r += 1

        # Doubles colinearites rapides
        found = False
        for r in range(R):
            k1 = ku[r]; k2 = kv[r]; k3 = kw[r]
            for s in range(r + 1, R):
                if k1 == ku[s] and k2 == kv[s]:
                    for i in range(N3): W[r, i] = (W[r, i] + W[s, i]) % p
                    R = delete_term(U, V, W, ku, kv, kw, lw, s, R, N1, N2, N3)
                    update_term_keys(U, V, W, ku, kv, kw, lw, r, N1, N2, N3, p, inv)
                    changed = True; found = True; break
                elif k1 == ku[s] and k3 == kw[s]:
                    c = (lw[s] * inv[lw[r]]) % p
                    for i in range(N2): V[r, i] = (V[r, i] + c * V[s, i]) % p
                    R = delete_term(U, V, W, ku, kv, kw, lw, s, R, N1, N2, N3)
                    update_term_keys(U, V, W, ku, kv, kw, lw, r, N1, N2, N3, p, inv)
                    changed = True; found = True; break
                elif k2 == kv[s] and k3 == kw[s]:
                    c = (lw[s] * inv[lw[r]]) % p
                    for i in range(N1): U[r, i] = (U[r, i] + c * U[s, i]) % p
                    R = delete_term(U, V, W, ku, kv, kw, lw, s, R, N1, N2, N3)
                    update_term_keys(U, V, W, ku, kv, kw, lw, r, N1, N2, N3, p, inv)
                    changed = True; found = True; break
            if found:
                break

        if not changed and check_slices:
            # Tester dependances lineaires dans les tranches de taille >= 3
            for mode in range(3):
                for r in range(R):
                    tk = ku[r] if mode == 0 else (kv[r] if mode == 1 else kw[r])
                    # Eviter de tester 2 fois la meme cle : verifier qu'elle n'apparait pas avant r
                    seen_before = False
                    cnt = 1
                    for j in range(r):
                        kj = ku[j] if mode == 0 else (kv[j] if mode == 1 else kw[j])
                        if kj == tk:
                            seen_before = True; break
                    if seen_before:
                        continue
                    for j in range(r + 1, R):
                        kj = ku[j] if mode == 0 else (kv[j] if mode == 1 else kw[j])
                        if kj == tk:
                            cnt += 1
                    if cnt >= 3:
                        R, red = reduce_one_slice(U, V, W, ku, kv, kw, lw, R, mode, tk, N1, N2, N3, p, inv)
                        if red:
                            changed = True
                            break
                if changed:
                    break
    return R


@njit(cache=True)
def do_flip(U, V, W, ku, kv, kw, lw, r, s, mode, orient, alpha, N1, N2, N3, p, inv):
    if mode == 0:
        if orient == 0:
            for i in range(N2): V[r, i] = (V[r, i] + alpha * V[s, i]) % p
            for i in range(N3): W[s, i] = (W[s, i] - alpha * W[r, i]) % p
        else:
            for i in range(N3): W[r, i] = (W[r, i] + alpha * W[s, i]) % p
            for i in range(N2): V[s, i] = (V[s, i] - alpha * V[r, i]) % p
    elif mode == 1:
        if orient == 0:
            for i in range(N1): U[r, i] = (U[r, i] + alpha * U[s, i]) % p
            for i in range(N3): W[s, i] = (W[s, i] - alpha * W[r, i]) % p
        else:
            for i in range(N3): W[r, i] = (W[r, i] + alpha * W[s, i]) % p
            for i in range(N1): U[s, i] = (U[s, i] - alpha * U[r, i]) % p
    else:
        c = (lw[s] * inv[lw[r]]) % p
        ac = (alpha * c) % p
        if orient == 0:
            for i in range(N1): U[r, i] = (U[r, i] + ac * U[s, i]) % p
            for i in range(N2): V[s, i] = (V[s, i] - alpha * V[r, i]) % p
        else:
            for i in range(N2): V[r, i] = (V[r, i] + ac * V[s, i]) % p
            for i in range(N1): U[s, i] = (U[s, i] - alpha * U[r, i]) % p
    update_term_keys(U, V, W, ku, kv, kw, lw, r, N1, N2, N3, p, inv)
    update_term_keys(U, V, W, ku, kv, kw, lw, s, N1, N2, N3, p, inv)


@njit(cache=True)
def do_plus(U, V, W, ku, kv, kw, lw, R, r, s, perm_mode, alpha, N1, N2, N3, p, inv):
    rn = R
    if perm_mode == 0:
        for i in range(N1): U[rn, i] = U[r, i]
        for i in range(N2): V[rn, i] = (alpha * (V[r, i] - V[s, i])) % p
        for i in range(N3): W[rn, i] = W[s, i]
        for i in range(N3): W[r, i] = (W[r, i] - alpha * W[s, i]) % p
        for i in range(N1): U[s, i] = (alpha * U[r, i] + U[s, i]) % p
    elif perm_mode == 1:
        for i in range(N2): V[rn, i] = V[r, i]
        for i in range(N3): W[rn, i] = (alpha * (W[r, i] - W[s, i])) % p
        for i in range(N1): U[rn, i] = U[s, i]
        for i in range(N1): U[r, i] = (U[r, i] - alpha * U[s, i]) % p
        for i in range(N2): V[s, i] = (alpha * V[r, i] + V[s, i]) % p
    else:
        for i in range(N3): W[rn, i] = W[r, i]
        for i in range(N1): U[rn, i] = (alpha * (U[r, i] - U[s, i])) % p
        for i in range(N2): V[rn, i] = V[s, i]
        for i in range(N2): V[r, i] = (V[r, i] - alpha * V[s, i]) % p
        for i in range(N3): W[s, i] = (alpha * W[r, i] + W[s, i]) % p

    update_term_keys(U, V, W, ku, kv, kw, lw, r, N1, N2, N3, p, inv)
    update_term_keys(U, V, W, ku, kv, kw, lw, s, N1, N2, N3, p, inv)
    update_term_keys(U, V, W, ku, kv, kw, lw, rn, N1, N2, N3, p, inv)
    return R + 1


@njit(cache=True)
def scheme_signature(ku, kv, kw, R):
    h = np.uint64(0)
    for r in range(R):
        x = np.uint64(ku[r]) * np.uint64(1000003) ^ np.uint64(kv[r]) * np.uint64(1000033) ^ np.uint64(kw[r])
        x = (x ^ (x >> np.uint64(30))) * np.uint64(0xbf58476d1ce4e5b9)
        h += x
    return h


@njit(cache=True)
def flip_graph_worker(
    U_init, V_init, W_init, p, max_steps, burst_steps,
    plus_depth, seed, target_rank, pool_cap=32
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
            # Alterner : marche pure a R_best (si flips disponibles) vs saut PLUS (+1 ou +2)
            restart_cnt += 1
            do_jump = (n_pairs == 0) or (restart_cnt % 3 != 0)
            if do_jump:
                n_jumps = 1 if (plus_depth <= 1 or restart_cnt % 5 != 0) else plus_depth
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

        pick = np.random.randint(0, n_pairs)
        r = pr[pick]; s = ps[pick]; mode = pm[pick]
        orient = np.random.randint(0, 2)
        alpha = np.random.randint(1, p)
        do_flip(U, V, W, ku, kv, kw, lw, r, s, mode, orient, alpha, N1, N2, N3, p, inv)
        total_flips += 1

        # Reduction rapide a chaque pas + reduction complete de tranches tous les 4 pas !
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
                        seen = True
                        break
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

    U_out = pool_U[0, :R_best, :].copy()
    V_out = pool_V[0, :R_best, :].copy()
    W_out = pool_W[0, :R_best, :].copy()
    return R_best, U_out, V_out, W_out, total_flips, total_plus, total_reds, pool_size


def _run_one_worker(args):
    m, k, n, p, max_steps, burst, pdepth, seed, target_rank = args
    U0, V0, W0 = build_schoolbook_fp(m, k, n, p)
    R_b, Ub, Vb, Wb, nf, np_, nr, psize = flip_graph_worker(
        U0, V0, W0, p, max_steps, burst, pdepth, seed, target_rank
    )
    ok = verify_scheme_fp(Ub, Vb, Wb, m, k, n, p)
    return {
        "seed": seed, "R_best": int(R_b), "verified": bool(ok),
        "flips": int(nf), "plus": int(np_), "reds": int(nr), "pool_size": int(psize),
        "factors": (Ub, Vb, Wb),
    }


def run_flip_campaign(pool, label, m, k, n, p, n_workers_runs, max_steps, burst, pdepth, target_rank, seed_base):
    ts = time.perf_counter()
    tasks = [
        (m, k, n, p, max_steps, burst, pdepth, seed_base + w, target_rank)
        for w in range(n_workers_runs)
    ]
    outs = list(pool.map(_run_one_worker, tasks))
    dt = time.perf_counter() - ts
    all_ok = all(o["verified"] for o in outs)
    assert all_ok, f"Verification failed on {label}!"

    outs.sort(key=lambda z: z["R_best"])
    best = outs[0]
    tot_flips = sum(o["flips"] for o in outs)
    tot_plus = sum(o["plus"] for o in outs)
    ranks = [o["R_best"] for o in outs]
    rate_m = (tot_flips / dt) / 1e6

    saved = None
    if best["R_best"] <= target_rank + 1:
        fname = f"flip_{label}_R{best['R_best']}.npz"
        saved = os.path.join(WITNESS_DIR, fname)
        Ub, Vb, Wb = best["factors"]
        np.savez(saved, U=Ub, V=Vb, W=Wb, shape=np.array([m, k, n]), p=np.array([p]))

    print(
        f"   {label:<12} | <{m},{k},{n}> F_{p} | best R={best['R_best']:>2} (cible {target_rank:>2}) | "
        f"rangs={ranks[:6]} | {tot_flips/1e6:>5.1f}M flips ({rate_m:>4.1f}M/s) | +{tot_plus:>5} PLUS | {dt:>5.1f}s | exact={all_ok}",
        flush=True,
    )
    return {
        "label": label, "shape": [m, k, n], "p": p, "target_rank": target_rank,
        "best_R": best["R_best"], "all_ranks": ranks,
        "total_flips": tot_flips, "total_plus": tot_plus,
        "flips_per_sec": tot_flips / dt, "time_s": dt,
        "verified_100pct": all_ok, "saved_npz": saved,
    }


def main():
    t0 = time.time()
    print("== Compilation Numba JIT + Sanity check exact <2,2,2> (Strassen) sur F_3 et F_5 ==\n", flush=True)
    for p in (3, 5):
        U0, V0, W0 = build_schoolbook_fp(2, 2, 2, p)
        Rb, Ub, Vb, Wb, nf, np_, nr, _ = flip_graph_worker(U0, V0, W0, p, 50000, 400, 1, 42 + p, 7)
        ok = verify_scheme_fp(Ub, Vb, Wb, 2, 2, 2, p)
        print(f"   Sanity <2,2,2> F_{p} : 8 -> R={Rb} en {nf} flips | verifie exact = {ok}", flush=True)
        assert ok and Rb == 7

    n_workers = min(10, os.cpu_count() or 4)
    print(f"\n== Chasse combinatoire EXACT Flip-Graph + Reduction de Tranches + PLUS ({n_workers} coeurs Numba JIT) ==\n", flush=True)

    campaigns = [
        ("2x3x4_F3", 2, 3, 4, 3, 10,  1_200_000, 1000, 2, 19, 1000),
        ("2x3x4_F5", 2, 3, 4, 5, 10,  1_200_000, 1000, 2, 19, 1100),
        ("3x3x3_F3", 3, 3, 3, 3, 10,  2_000_000, 1500, 2, 22, 2000),
        ("3x3x3_F5", 3, 3, 3, 5, 10,  2_000_000, 1500, 2, 22, 2100),
        ("2x4x4_F3", 2, 4, 4, 3, 10,  1_200_000, 1200, 2, 25, 3000),
        ("2x4x4_F5", 2, 4, 4, 5, 10,  1_200_000, 1200, 2, 25, 3100),
        ("3x3x4_F3", 3, 3, 4, 3, 10,  1_500_000, 1500, 2, 28, 4000),
        ("3x3x4_F5", 3, 3, 4, 5, 10,  1_500_000, 1500, 2, 28, 4100),
        ("2x4x5_F3", 2, 4, 5, 3, 10,  1_200_000, 1600, 2, 31, 5000),
        ("2x4x5_F5", 2, 4, 5, 5, 10,  1_200_000, 1600, 2, 31, 5100),
        ("3x4x4_F3", 3, 4, 4, 3, 10,  1_200_000, 1800, 2, 37, 6000),
        ("3x4x4_F5", 3, 4, 4, 5, 10,  1_200_000, 1800, 2, 37, 6100),
    ]

    report = {"n_workers": n_workers, "campaigns": {}}
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as pool:
        for cfg in campaigns:
            res = run_flip_campaign(pool, *cfg)
            report["campaigns"][cfg[0]] = res

    out_json = os.path.join(RESULTS_DIR, "flip_graph_fp_results.json")
    with open(out_json, "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14.0, 5.0))

    fmts = ["2x3x4", "3x3x3", "2x4x4", "3x3x4", "2x4x5", "3x4x4"]
    schoolbook = [24, 27, 32, 36, 40, 48]
    world_rec = [20, 23, 26, 29, 32, 38]
    best_f3 = [report["campaigns"][f"{f}_F3"]["best_R"] for f in fmts]
    best_f5 = [report["campaigns"][f"{f}_F5"]["best_R"] for f in fmts]
    x = np.arange(len(fmts))
    w = 0.22

    ax1.bar(x - 1.5 * w, schoolbook, width=w, color="#cccccc", label="Naïf (m·k·n)")
    ax1.bar(x - 0.5 * w, world_rec, width=w, color="#2ca02c", label="Record mondial exact")
    ax1.bar(x + 0.5 * w, best_f3, width=w, color="#1f77b4", label="Notre Flip-Graph F_3")
    ax1.bar(x + 1.5 * w, best_f5, width=w, color="#d62728", label="Notre Flip-Graph F_5")
    ax1.set_xticks(x)
    ax1.set_xticklabels([f"<{f.replace('x',',')}>" for f in fmts])
    ax1.set_ylabel("Rang exact (multiplications scalaires)")
    ax1.set_title("Graphe de Flips exact + Tranches sur F_3 & F_5 vs Records")
    ax1.grid(axis="y", alpha=0.3)
    ax1.legend(fontsize=8)

    gap_f3 = [b - r for b, r in zip(best_f3, world_rec)]
    gap_f5 = [b - r for b, r in zip(best_f5, world_rec)]
    ax2.plot(x - 0.08, gap_f3, "o-", color="#1f77b4", markersize=9, linewidth=2, label="Écart au record sur F_3")
    ax2.plot(x + 0.08, gap_f5, "s--", color="#d62728", markersize=9, linewidth=2, label="Écart au record sur F_5")
    ax2.axhline(0, color="#2ca02c", linestyle=":", linewidth=2, label="Record mondial égalé (écart = 0)")
    ax2.axhline(-1, color="#9467bd", linestyle="--", linewidth=1.5, label="NOUVEAU RECORD MONDIAL (écart = -1)")
    ax2.set_xticks(x)
    ax2.set_xticklabels([f"<{f.replace('x',',')}>" for f in fmts])
    ax2.set_ylabel("Meilleur rang trouvé − Record mondial")
    ax2.set_title("Distance combinatoire exacte au record mondial (zéro bord)")
    ax2.grid(alpha=0.3)
    ax2.legend(fontsize=8)

    fig.tight_layout()
    out_png = os.path.join(OUT, "flip_graph_fp_hunt.png")
    fig.savefig(out_png, dpi=130)
    print(f"\nTemps total : {time.time() - t0:.1f} s", flush=True)
    print(f"Graphique : {out_png}", flush=True)


if __name__ == "__main__":
    main()
