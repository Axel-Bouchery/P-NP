"""
Piste 3 : Rang tensoriel exact et rang de bord de la multiplication de matrices 3x3
===================================================================================

Solveur Gauss-Newton amorti / Levenberg-Marquardt EXACT (forme fermee hermitienne)
sur le tenseur de multiplication matricielle T_<n,n,n> (n=2 calibration Strassen,
n=3 cible Laderman/Smirnov) sur R (float64) et C (complex128).

Deux modes structures :
  1. Mode libre (sans contrainte de symetrie) : 3 * n^2 * R parametres.
  2. Mode Z_3-cyclique (invariance Tr(ABC) = Tr(BCA)) : R = F + 3*C
     (division par 3 du nombre de parametres : 198 params pour R=22, 171 pour R=19 !).

Et deux regimes physiques :
  - Regime BORD (mu = 0, norme libre) : chasse au rang de bord R_border(<3,3,3>) <= 20, 19, 18.
  - Regime EXACT (regularisation L_{2/3} recuite + equilibrage AM-GM -> blocage des
    degenerescences de bord, puis polissage mu = 0) : chasse au rang exact R(<3,3,3>) = 22
    sur R et C (sachant que R_{F_2} = R_Z = 23 vient d'etre prouve optimal !).
"""
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = sys.argv[1] if len(sys.argv) > 1 else HERE


def build_matmul_tensor(n):
    """
    Tenseur T_<n,n,n> de forme (n^2, n^2, n^2) pour Tr(A B C) :
      T[(i,k), (k,j), (j,i)] = 1  (forme cycliquement symetrique !).
    Verification : sum_{i,j,k} A_{i,k} B_{k,j} C_{j,i} = Tr(A B C),
    ce qui est exactement equivalent a C^T = A B et strictement invariant par
    permutation cyclique (A, B, C) -> (B, C, A) !
    """
    N = n * n
    T = np.zeros((N, N, N), dtype=np.float64)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                a_idx = i * n + k
                b_idx = k * n + j
                c_idx = j * n + i
                T[a_idx, b_idx, c_idx] = 1.0
    return T


def build_cyclic_projection(N, F, C):
    """
    Construit la matrice de projection lineaire P (taille 3*N*R x N*(F + 3*C))
    telle que theta_full = P @ theta_cyc impose la Z_3-symetrie cyclique :
      U = [A_fix, U_c, V_c, W_c]
      V = [A_fix, V_c, W_c, U_c]
      W = [A_fix, W_c, U_c, V_c]
    avec R = F + 3*C.
    Convention de vectorisation de theta_full : [vec_row(U), vec_row(V), vec_row(W)]
    ou chaque facteur (N, R) est aplati en ordre C (i*R + r).
    """
    R = F + 3 * C
    n_full = 3 * N * R
    n_cyc = N * (F + 3 * C)
    P = np.zeros((n_full, n_cyc), dtype=np.float64)

    # Indices dans theta_cyc :
    #   A_fix : i*F + f                         (0 .. N*F - 1)
    #   U_c   : N*F + i*C + c                   (N*F .. N*F + N*C - 1)
    #   V_c   : N*F + N*C + i*C + c
    #   W_c   : N*F + 2*N*C + i*C + c
    off_U = 0
    off_V = N * R
    off_W = 2 * N * R

    for i in range(N):
        for f in range(F):
            col = i * F + f
            P[off_U + i * R + f, col] = 1.0
            P[off_V + i * R + f, col] = 1.0
            P[off_W + i * R + f, col] = 1.0
        for c in range(C):
            col_u = N * F + i * C + c
            col_v = N * F + N * C + i * C + c
            col_w = N * F + 2 * N * C + i * C + c
            # U : [..., U_c, V_c, W_c]
            P[off_U + i * R + (F + c), col_u] = 1.0
            P[off_U + i * R + (F + C + c), col_v] = 1.0
            P[off_U + i * R + (F + 2 * C + c), col_w] = 1.0
            # V : [..., V_c, W_c, U_c]
            P[off_V + i * R + (F + c), col_v] = 1.0
            P[off_V + i * R + (F + C + c), col_w] = 1.0
            P[off_V + i * R + (F + 2 * C + c), col_u] = 1.0
            # W : [..., W_c, U_c, V_c]
            P[off_W + i * R + (F + c), col_w] = 1.0
            P[off_W + i * R + (F + C + c), col_u] = 1.0
            P[off_W + i * R + (F + 2 * C + c), col_v] = 1.0
    return P


def cp_reconstruct(U, V, W):
    return np.einsum("ir,jr,kr->ijk", U, V, W)


def balance_factors(U, V, W):
    """Equilibrage de jauge AM-GM exact : ||u_r|| = ||v_r|| = ||w_r||."""
    nu = np.linalg.norm(U, axis=0)
    nv = np.linalg.norm(V, axis=0)
    nw = np.linalg.norm(W, axis=0)
    prod = nu * nv * nw
    mask = prod > 1e-30
    s = np.ones_like(prod)
    s[mask] = np.cbrt(prod[mask])
    Un = U.copy(); Vn = V.copy(); Wn = W.copy()
    Un[:, mask] *= (s[mask] / nu[mask])[None, :]
    Vn[:, mask] *= (s[mask] / nv[mask])[None, :]
    Wn[:, mask] *= (s[mask] / nw[mask])[None, :]
    return Un, Vn, Wn


def normal_equations_cp(U, V, W, T_unfold):
    """
    Calcule en forme fermee EXACTE (R ou C hermitien) :
      g = J^H (T_hat - T)  de taille 3*N*R
      H = J^H J            de taille (3*N*R) x (3*N*R)
      cost = 0.5 * ||T_hat - T||_F^2
    ou T_unfold = (T1, T2, T3) sont les 3 depliages (N, N^2) de T et ||T||_F^2.
    """
    N, R = U.shape
    T1, T2, T3, T_norm2 = T_unfold
    dtype = U.dtype

    GU = U.conj().T @ U  # (R, R)
    GV = V.conj().T @ V
    GW = W.conj().T @ W

    GVW = GV * GW
    GUW = GU * GW
    GUV = GU * GV

    # MTTKRP via produits matriciels (ultra rapide)
    # T1 est (N, N*N) indexe par (i, j*N + k) -> on veut sum_{j,k} T_{i,j,k} conj(V_{j,r}) conj(W_{k,r})
    # Produit de Khatri-Rao conj(V) (*) conj(W) de forme (N*N, R)
    VW_c = (V.conj()[:, None, :] * W.conj()[None, :, :]).reshape(N * N, R)
    UW_c = (U.conj()[:, None, :] * W.conj()[None, :, :]).reshape(N * N, R)
    UV_c = (U.conj()[:, None, :] * V.conj()[None, :, :]).reshape(N * N, R)

    M1 = T1 @ VW_c  # (N, R)
    M2 = T2 @ UW_c
    M3 = T3 @ UV_c

    # Note : sum_s U_{i,s} GVW_{r,s} = (U @ GVW.T)_{i,r} = (U @ GVW.conj())_{i,r}
    gU = U @ GVW.T - M1
    gV = V @ GUW.T - M2
    gW = W @ GUV.T - M3
    g = np.concatenate([gU.ravel(), gV.ravel(), gW.ravel()])

    # Cout exact : 0.5 * (||T_hat||^2 - 2 Re<T_hat, T> + ||T||^2)
    # ||T_hat||^2 = sum_{r,s} GU_{r,s} GV_{r,s} GW_{r,s} = real(sum(GU * GVW))
    norm_hat2 = float(np.real(np.sum(GU * GVW)))
    inner = float(np.real(np.sum(U.conj() * M1)))
    # Note: M1 = T1 @ (V^* (*) W^*), donc sum U^* M1 = sum T_{ijk} U^* V^* W^* = conj(<T_hat, T>), meme partie reelle!
    diff2 = max(0.0, norm_hat2 - 2.0 * inner + T_norm2)
    cost = 0.5 * diff2

    # Construction vectorisee exacte de H = J^H J de taille (3*N*R, 3*N*R)
    NR = N * R
    H = np.zeros((3 * NR, 3 * NR), dtype=dtype)

    # Blocs diagonaux : I_N kron GVW, I_N kron GUW, I_N kron GUV
    for i in range(N):
        sl = slice(i * R, (i + 1) * R)
        H[sl, sl] = GVW
        H[NR + sl.start:NR + sl.stop, NR + sl.start:NR + sl.stop] = GUW
        H[2 * NR + sl.start:2 * NR + sl.stop, 2 * NR + sl.start:2 * NR + sl.stop] = GUV

    # Blocs hors-diagonaux :
    # H_UV[(i,r), (j,s)] = conj(V_{j,r}) * U_{i,s} * GW_{r,s}
    H_UV = np.einsum("jr,is,rs->irjs", V.conj(), U, GW).reshape(NR, NR)
    H_UW = np.einsum("kr,is,rs->irks", W.conj(), U, GV).reshape(NR, NR)
    H_VW = np.einsum("kr,js,rs->jrks", W.conj(), V, GU).reshape(NR, NR)

    H[:NR, NR:2 * NR] = H_UV
    H[NR:2 * NR, :NR] = H_UV.conj().T
    H[:NR, 2 * NR:] = H_UW
    H[2 * NR:, :NR] = H_UW.conj().T
    H[NR:2 * NR, 2 * NR:] = H_VW
    H[2 * NR:, NR:2 * NR] = H_VW.conj().T

    return cost, g, H


def prepare_unfoldings(T):
    N = T.shape[0]
    T1 = T.reshape(N, N * N)
    T2 = np.transpose(T, (1, 0, 2)).reshape(N, N * N)
    T3 = np.transpose(T, (2, 0, 1)).reshape(N, N * N)
    return (T1, T2, T3, float(np.sum(T * T)))


def verify_normal_equations():
    """Verification exacte du gradient et du Hessien Gauss-Newton sur C (complex128)."""
    rng = np.random.default_rng(123)
    T = build_matmul_tensor(2)
    Tu = prepare_unfoldings(T)
    U = rng.normal(size=(4, 5)) + 1j * rng.normal(size=(4, 5))
    V = rng.normal(size=(4, 5)) + 1j * rng.normal(size=(4, 5))
    W = rng.normal(size=(4, 5)) + 1j * rng.normal(size=(4, 5))
    cost, g, H = normal_equations_cp(U, V, W, Tu)
    # Cout direct
    T_hat = cp_reconstruct(U, V, W)
    c_direct = 0.5 * float(np.sum(np.abs(T_hat - T) ** 2))
    # Jacobienne explicite par differences finies complexes
    N, R = U.shape
    theta = np.concatenate([U.ravel(), V.ravel(), W.ravel()])
    J = np.zeros((N ** 3, 3 * N * R), dtype=np.complex128)
    eps = 1e-7
    for p in range(3 * N * R):
        th2 = theta.copy(); th2[p] += eps
        U2 = th2[:N * R].reshape(N, R)
        V2 = th2[N * R:2 * N * R].reshape(N, R)
        W2 = th2[2 * N * R:].reshape(N, R)
        J[:, p] = (cp_reconstruct(U2, V2, W2) - T_hat).ravel() / eps
    g_ref = J.conj().T @ (T_hat - T).ravel()
    H_ref = J.conj().T @ J
    return (
        abs(cost - c_direct),
        float(np.linalg.norm(g - g_ref) / np.linalg.norm(g_ref)),
        float(np.linalg.norm(H - H_ref) / np.linalg.norm(H_ref)),
    )


def solve_cp_lm(
    T, R, is_complex=False, seed=0, max_iter=450,
    cyclic_FC=None, mu_schedule=None, max_norm_stop=1e5,
    init_scale=0.55
):
    """
    Solveur Levenberg-Marquardt exact sur T (float64 ou complex128).
    - cyclic_FC = (F, C) avec F + 3*C = R pour le sous-espace Z_3-cyclique, ou None (libre).
    - mu_schedule : liste de (n_iters, mu) pour recuit L_{2/3} (mode EXACT) ou None (mode BORD mu=0).
    """
    N = T.shape[0]
    Tu = prepare_unfoldings(T)
    T_norm = math.sqrt(Tu[3])
    rng = np.random.default_rng(seed)

    P_cyc = None
    if cyclic_FC is not None:
        F, C = cyclic_FC
        assert F + 3 * C == R
        P_cyc = build_cyclic_projection(N, F, C)
        n_par = P_cyc.shape[1]
        if is_complex:
            theta = (rng.normal(0.0, init_scale, n_par) + 1j * rng.normal(0.0, init_scale, n_par)) / math.sqrt(2)
        else:
            theta = rng.normal(0.0, init_scale, n_par)
    else:
        n_par = 3 * N * R
        if is_complex:
            theta = (rng.normal(0.0, init_scale, n_par) + 1j * rng.normal(0.0, init_scale, n_par)) / math.sqrt(2)
        else:
            theta = rng.normal(0.0, init_scale, n_par)

    NR = N * R

    def to_UVW(th):
        full = (P_cyc @ th) if (P_cyc is not None) else th
        return full[:NR].reshape(N, R), full[NR:2 * NR].reshape(N, R), full[2 * NR:].reshape(N, R)

    def eval_data_cost(th):
        U, V, W = to_UVW(th)
        GU = U.conj().T @ U; GV = V.conj().T @ V; GW = W.conj().T @ W
        GVW = GV * GW
        VW_c = (V.conj()[:, None, :] * W.conj()[None, :, :]).reshape(N * N, R)
        M1 = Tu[0] @ VW_c
        norm_hat2 = float(np.real(np.sum(GU * GVW)))
        inner = float(np.real(np.sum(U.conj() * M1)))
        return 0.5 * max(0.0, norm_hat2 - 2.0 * inner + Tu[3])

    if mu_schedule is None:
        stages = [(max_iter, 0.0)]
    else:
        stages = mu_schedule

    total_it = 0
    best_err = 1e9
    best_th = theta.copy()
    history = []

    for n_st, mu in stages:
        lam = 1e-2 if mu == 0.0 else max(1e-3, mu)
        for st_it in range(n_st):
            # Equilibrage AM-GM periodique (en mode libre uniquement, conserve exactement T_hat)
            if P_cyc is None and (st_it % 15 == 0) and st_it > 0:
                U, V, W = to_UVW(theta)
                U, V, W = balance_factors(U, V, W)
                theta = np.concatenate([U.ravel(), V.ravel(), W.ravel()])

            U, V, W = to_UVW(theta)
            c_data, g_full, H_full = normal_equations_cp(U, V, W, Tu)
            err = math.sqrt(2.0 * c_data) / T_norm
            norm_f = float(np.linalg.norm(np.concatenate([U.ravel(), V.ravel(), W.ravel()])))

            if err < best_err:
                best_err = err
                best_th = theta.copy()

            if total_it % 40 == 0 or err < 1e-13:
                history.append((total_it, err, norm_f, mu))
            if err < 1e-13 and mu == 0.0:
                break
            if norm_f > max_norm_stop:
                break

            if P_cyc is not None:
                g = P_cyc.T @ g_full
                H = P_cyc.T @ H_full @ P_cyc
                # Poids de ridge = 3 pour les termes cycliques (equivaut a ||U||^2+||V||^2+||W||^2)
                pen_w = 3.0 * mu
            else:
                g = g_full
                H = H_full
                pen_w = mu

            c_tot = c_data + 0.5 * pen_w * float(np.real(np.vdot(theta, theta)))
            if pen_w > 0.0:
                g = g + pen_w * theta
                H = H + pen_w * np.eye(n_par, dtype=H.dtype)

            diag_H = np.maximum(np.real(np.diag(H)), 1e-6)
            accepted = False
            for _ in range(12):
                Hd = H + lam * np.diag(diag_H)
                try:
                    step = -np.linalg.solve(Hd, g)
                except np.linalg.LinAlgError:
                    lam *= 10.0
                    continue
                th_new = theta + step
                c_new = eval_data_cost(th_new) + 0.5 * pen_w * float(np.real(np.vdot(th_new, th_new)))
                if math.isfinite(c_new) and c_new < c_tot:
                    theta = th_new
                    lam = max(lam * 0.35, 1e-14)
                    accepted = True
                    break
                else:
                    lam = min(lam * 4.5, 1e12)
            total_it += 1
            if not accepted:
                break

    # Evaluation directe exacte finale (sans erreur d'arrondi de la formule ||T_hat||^2)
    U, V, W = to_UVW(best_th)
    if P_cyc is None:
        U, V, W = balance_factors(U, V, W)
    T_hat = cp_reconstruct(U, V, W)
    diff = T_hat - T
    rel_err = float(np.linalg.norm(diff) / T_norm)
    max_err = float(np.max(np.abs(diff)))
    frob_norm = float(math.sqrt(
        np.sum(np.abs(U) ** 2) + np.sum(np.abs(V) ** 2) + np.sum(np.abs(W) ** 2)
    ))
    max_entry = float(max(np.max(np.abs(U)), np.max(np.abs(V)), np.max(np.abs(W))))

    return {
        "seed": seed,
        "rel_err": rel_err,
        "max_err": max_err,
        "frob_norm": frob_norm,
        "max_entry": max_entry,
        "iters": total_it,
        "history": history,
    }, (U, V, W)


def run_benchmark(label, T, R, is_complex, n_trials, max_iter=450,
                  cyclic_FC=None, mu_schedule=None, seed_base=0, save_if_below=1e-9):
    ts = time.perf_counter()
    trials = []
    best_factors = None
    best_err = 1e9
    for tr in range(n_trials):
        res, facs = solve_cp_lm(
            T, R, is_complex=is_complex, seed=seed_base + tr,
            max_iter=max_iter, cyclic_FC=cyclic_FC, mu_schedule=mu_schedule
        )
        trials.append(res)
        if res["rel_err"] < best_err:
            best_err = res["rel_err"]
            best_factors = facs
    trials.sort(key=lambda z: z["rel_err"])
    best = trials[0]
    med = float(np.median([z["rel_err"] for z in trials]))
    dt = time.perf_counter() - ts
    n_low = sum(1 for z in trials if z["rel_err"] < 1e-4)

    saved_path = None
    if best["rel_err"] < save_if_below:
        fname = f"best_{label.replace(' ', '_').replace('/', '_')}.npz"
        saved_path = os.path.join(HERE, fname)
        np.savez(saved_path, U=best_factors[0], V=best_factors[1], W=best_factors[2])

    print(
        f"   {label:<28} | R={R:>2} | {n_trials:>2} ess. | <1e-4:{n_low:>2} | "
        f"best_rel={best['rel_err']:>9.2e} | best_max={best['max_err']:>9.2e} | "
        f"||F||={best['frob_norm']:>7.1f} | max|x|={best['max_entry']:>6.1f} | {dt:>5.1f}s",
        flush=True,
    )
    return {
        "label": label, "R": R, "is_complex": is_complex, "cyclic_FC": cyclic_FC,
        "best": best, "median_rel_err": med, "n_below_1e-4": n_low,
        "trials": trials, "time_s": dt, "saved_npz": saved_path,
    }


def main():
    t0 = time.time()
    ec, eg, eh = verify_normal_equations()
    print(f"== Sanity check Gauss-Newton complexe : err_cost={ec:.2e}, err_g={eg:.2e}, err_H={eh:.2e} ==\n", flush=True)
    assert eg < 1e-6 and eh < 1e-6

    T2 = build_matmul_tensor(2)
    T3 = build_matmul_tensor(3)

    report = {"sanity": {"ec": ec, "eg": eg, "eh": eh}, "runs": {}}

    print("== PHASE A : Calibration exacte Strassen <2,2,2> et Laderman/Smirnov <3,3,3> ==\n", flush=True)
    exact_sched = [(80, 1e-2), (80, 1e-3), (80, 1e-4), (160, 0.0)]

    cal_configs = [
        ("Cal_2x2_R_libre",      T2, 7,  False, 6,  200, None,   None,        100),
        ("Cal_2x2_R_cyc_1+2x3",  T2, 7,  False, 6,  200, (1, 2), None,        200),
        ("Cal_2x2_C_bord_R6",    T2, 6,  True,  6,  250, None,   None,        300),
        ("Cal_3x3_R_exact_R23",  T3, 23, False, 6,  400, None,   exact_sched, 400),
        ("Cal_3x3_R_bord_R20",   T3, 20, False, 8,  450, None,   None,        500),
    ]
    for cfg in cal_configs:
        res = run_benchmark(*cfg)
        report["runs"][cfg[0]] = res

    print("\n== PHASE B : Frontiere OUVERTE 1 — Rang de bord R_border(<3,3,3>) <= 19, 18 ==\n", flush=True)
    bord_configs = [
        ("Bord_3x3_R_libre_R19",     T3, 19, False, 14, 550, None,   None, 1000),
        ("Bord_3x3_C_libre_R19",     T3, 19, True,  10, 500, None,   None, 1100),
        ("Bord_3x3_R_cyc_1+6x3_R19", T3, 19, False, 20, 500, (1, 6), None, 1200),
        ("Bord_3x3_C_cyc_1+6x3_R19", T3, 19, True,  16, 500, (1, 6), None, 1300),
        ("Bord_3x3_R_cyc_4+5x3_R19", T3, 19, False, 16, 500, (4, 5), None, 1400),
        ("Bord_3x3_R_libre_R18",     T3, 18, False,  8, 500, None,   None, 1500),
    ]
    for cfg in bord_configs:
        res = run_benchmark(*cfg)
        report["runs"][cfg[0]] = res

    print("\n== PHASE C : Frontiere OUVERTE 2 — Rang EXACT borne R(<3,3,3>) = 22 (R & C) ==\n", flush=True)
    # Recuit L_{2/3} pour forcer une solution bornee (non-degeneree) a R=22
    exact_sched_strong = [(100, 3e-2), (100, 5e-3), (100, 5e-4), (200, 0.0)]
    exact_configs = [
        ("Exact_3x3_R_libre_R22",     T3, 22, False, 14, 500, None,   exact_sched_strong, 2000),
        ("Exact_3x3_C_libre_R22",     T3, 22, True,  10, 500, None,   exact_sched_strong, 2100),
        ("Exact_3x3_R_cyc_1+7x3_R22", T3, 22, False, 20, 500, (1, 7), exact_sched_strong, 2200),
        ("Exact_3x3_C_cyc_1+7x3_R22", T3, 22, True,  16, 500, (1, 7), exact_sched_strong, 2300),
        ("Exact_3x3_R_cyc_4+6x3_R22", T3, 22, False, 16, 500, (4, 6), exact_sched_strong, 2400),
    ]
    for cfg in exact_configs:
        res = run_benchmark(*cfg)
        report["runs"][cfg[0]] = res

    out_json = os.path.join(HERE, "matmul3_tensor_results.json")
    with open(out_json, "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2)

    # Graphique de synthese
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14.0, 5.0))

    # Panel 1 : Erreur relative vs rang R sur <3,3,3> (bord vs exact borne)
    bord_pts = [
        (18, report["runs"]["Bord_3x3_R_libre_R18"]["best"]["rel_err"]),
        (19, min(report["runs"]["Bord_3x3_R_libre_R19"]["best"]["rel_err"],
                 report["runs"]["Bord_3x3_R_cyc_1+6x3_R19"]["best"]["rel_err"],
                 report["runs"]["Bord_3x3_R_cyc_4+5x3_R19"]["best"]["rel_err"])),
        (20, report["runs"]["Cal_3x3_R_bord_R20"]["best"]["rel_err"]),
        (23, report["runs"]["Cal_3x3_R_exact_R23"]["best"]["rel_err"]),
    ]
    bord_c_19 = min(report["runs"]["Bord_3x3_C_libre_R19"]["best"]["rel_err"],
                    report["runs"]["Bord_3x3_C_cyc_1+6x3_R19"]["best"]["rel_err"])
    exact_22_r = min(report["runs"]["Exact_3x3_R_libre_R22"]["best"]["rel_err"],
                     report["runs"]["Exact_3x3_R_cyc_1+7x3_R22"]["best"]["rel_err"],
                     report["runs"]["Exact_3x3_R_cyc_4+6x3_R22"]["best"]["rel_err"])
    exact_22_c = min(report["runs"]["Exact_3x3_C_libre_R22"]["best"]["rel_err"],
                     report["runs"]["Exact_3x3_C_cyc_1+7x3_R22"]["best"]["rel_err"])

    ax1.semilogy([p[0] for p in bord_pts], [max(p[1], 1e-15) for p in bord_pts],
                 "o-", color="#1f77b4", linewidth=2, markersize=8, label="Rang de bord libre sur R")
    ax1.semilogy([19], [max(bord_c_19, 1e-15)], "s", color="#d62728", markersize=9,
                 label="Rang de bord R=19 sur C")
    ax1.semilogy([22], [max(exact_22_r, 1e-15)], "^", color="#2ca02c", markersize=9,
                 label="Rang exact recuit R=22 sur R")
    ax1.semilogy([22], [max(exact_22_c, 1e-15)], "D", color="#9467bd", markersize=8,
                 label="Rang exact recuit R=22 sur C")
    ax1.axvline(20, color="#ff7f0e", linestyle=":", alpha=0.8, label="Smirnov R_bord <= 20")
    ax1.axvline(23, color="#2ca02c", linestyle="--", alpha=0.7, label="Laderman R = 23 (optimal sur F_2)")
    ax1.set_xticks([18, 19, 20, 21, 22, 23])
    ax1.set_xlabel("Rang tensoriel R pour <3, 3, 3>")
    ax1.set_ylabel("Meilleure erreur relative ||T_hat - T|| / ||T|| (log)")
    ax1.set_title("Multiplication de matrices 3x3 : Rang exact vs Rang de bord")
    ax1.grid(alpha=0.3, which="both")
    ax1.legend(fontsize=8)

    # Panel 2 : Nuage Norme ||U,V,W||_F vs Erreur relative (distinguer EXACT vs BORD !)
    plot_groups = [
        ("Cal_3x3_R_exact_R23", "#2ca02c", "*", "Exact R=23 (Laderman)"),
        ("Cal_3x3_R_bord_R20",  "#ff7f0e", "o", "Bord R=20 (Smirnov)"),
        ("Bord_3x3_R_libre_R19","#1f77b4", "^", "Bord R=19 (R libre)"),
        ("Bord_3x3_C_libre_R19","#d62728", "s", "Bord R=19 (C libre)"),
        ("Exact_3x3_R_libre_R22","#17becf","v", "Exact recuit R=22 (R)"),
        ("Exact_3x3_C_libre_R22","#9467bd","D", "Exact recuit R=22 (C)"),
    ]
    for key, col, mk, lbl in plot_groups:
        trs = report["runs"][key]["trials"]
        ax2.loglog(
            [t["frob_norm"] for t in trs],
            [max(t["rel_err"], 1e-15) for t in trs],
            mk, color=col, label=lbl, markersize=6.5, alpha=0.85,
        )
    ax2.set_xlabel("Norme de Frobenius ||U, V, W||_F (log)")
    ax2.set_ylabel("Erreur relative ||T_hat - T|| / ||T|| (log)")
    ax2.set_title("Signature Exact (norme bornée) vs Bord (norme divergente)")
    ax2.grid(alpha=0.3, which="both")
    ax2.legend(fontsize=7.8)

    fig.tight_layout()
    out_png = os.path.join(OUT, "matmul3_tensor_hunt.png")
    fig.savefig(out_png, dpi=130)
    print(f"\nTemps total : {time.time() - t0:.1f} s", flush=True)
    print(f"Graphique : {out_png}", flush=True)


if __name__ == "__main__":
    main()
