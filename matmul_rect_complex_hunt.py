"""
Piste 3 (Phase 3) : Chasse sur formats rectangulaires ouverts sur C (complex128)
================================================================================

Formats cibles <m, p, n> (records mondiaux exacts sur C, Sedoglavic / Smirnov / AlphaEvolve) :
  - <2, 3, 4> : record = 20  -> cible ouverte R = 19 sur C (cal. 20)
  - <2, 3, 5> : record = 25  -> cible ouverte R = 24 sur C (cal. 25)
  - <2, 4, 4> : record = 26  -> cible ouverte R = 25 sur C (cal. 26)
  - <3, 3, 4> : record = 29  -> cible ouverte R = 28 sur C (cal. 29)

Outillage :
  1. Solveur Gauss-Newton amorti (Levenberg-Marquardt) HERMITIEN EXACT rectangulaire
     (N1 = m*p, N2 = p*n, N3 = m*n) en complex128, parallelise sur 8 coeurs.
  2. Recuit L_{2/3} + equilibrage AM-GM + detecteur de pente alpha = d log(err)/d log(mu)
     (separe instantanement rang exact borne alpha ~ 1.0 vs rang de bord alpha ~ 0.30..0.42).
  3. Double initialisation complexe (gaussienne circulaire + cyclotomique {0, +-1, +-i, w})
     + compression warm-start depuis nos solutions exactes R_rec -> R_rec - 1.
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import concurrent.futures
import json
import math
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = sys.argv[1] if len(sys.argv) > 1 else HERE


def build_rect_matmul_tensor(m, p, n):
    """
    Tenseur T_<m,p,n> de forme (m*p, p*n, m*n) tel que
      (A B)_{i,j} = sum_k A_{i,k} B_{k,j}.
    """
    N1, N2, N3 = m * p, p * n, m * n
    T = np.zeros((N1, N2, N3), dtype=np.float64)
    for i in range(m):
        for k in range(p):
            for j in range(n):
                T[i * p + k, k * n + j, i * n + j] = 1.0
    return T


def prepare_rect_unfoldings(T):
    N1, N2, N3 = T.shape
    T1 = T.reshape(N1, N2 * N3)
    T2 = np.transpose(T, (1, 0, 2)).reshape(N2, N1 * N3)
    T3 = np.transpose(T, (2, 0, 1)).reshape(N3, N1 * N2)
    return (T1, T2, T3, float(np.sum(T * T)))


def cp_reconstruct(U, V, W):
    return np.einsum("ir,jr,kr->ijk", U, V, W)


def balance_factors(U, V, W):
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


def normal_equations_rect_cp(U, V, W, Tu):
    """
    Equations normales Gauss-Newton hermitiennes EXACTES pour un tenseur rectangulaire
    N1 x N2 x N3 en float64 ou complex128.
    """
    N1, R = U.shape
    N2 = V.shape[0]
    N3 = W.shape[0]
    T1, T2, T3, T_norm2 = Tu
    dtype = U.dtype

    GU = U.conj().T @ U
    GV = V.conj().T @ V
    GW = W.conj().T @ W

    GVW = GV * GW
    GUW = GU * GW
    GUV = GU * GV

    VW_c = (V.conj()[:, None, :] * W.conj()[None, :, :]).reshape(N2 * N3, R)
    UW_c = (U.conj()[:, None, :] * W.conj()[None, :, :]).reshape(N1 * N3, R)
    UV_c = (U.conj()[:, None, :] * V.conj()[None, :, :]).reshape(N1 * N2, R)

    M1 = T1 @ VW_c
    M2 = T2 @ UW_c
    M3 = T3 @ UV_c

    gU = U @ GVW.T - M1
    gV = V @ GUW.T - M2
    gW = W @ GUV.T - M3
    g = np.concatenate([gU.ravel(), gV.ravel(), gW.ravel()])

    norm_hat2 = float(np.real(np.sum(GU * GVW)))
    inner = float(np.real(np.sum(U.conj() * M1)))
    cost = 0.5 * max(0.0, norm_hat2 - 2.0 * inner + T_norm2)

    S1 = N1 * R
    S2 = N2 * R
    S3 = N3 * R
    P = S1 + S2 + S3
    H = np.zeros((P, P), dtype=dtype)

    for i in range(N1):
        sl = slice(i * R, (i + 1) * R)
        H[sl, sl] = GVW
    for j in range(N2):
        sl = slice(S1 + j * R, S1 + (j + 1) * R)
        H[sl, sl] = GUW
    for k in range(N3):
        sl = slice(S1 + S2 + k * R, S1 + S2 + (k + 1) * R)
        H[sl, sl] = GUV

    H_UV = np.einsum("jr,is,rs->irjs", V.conj(), U, GW).reshape(S1, S2)
    H_UW = np.einsum("kr,is,rs->irks", W.conj(), U, GV).reshape(S1, S3)
    H_VW = np.einsum("kr,js,rs->jrks", W.conj(), V, GU).reshape(S2, S3)

    H[:S1, S1:S1 + S2] = H_UV
    H[S1:S1 + S2, :S1] = H_UV.conj().T
    H[:S1, S1 + S2:] = H_UW
    H[S1 + S2:, :S1] = H_UW.conj().T
    H[S1:S1 + S2, S1 + S2:] = H_VW
    H[S1 + S2:, S1:S1 + S2] = H_VW.conj().T

    return cost, g, H


def verify_rect_complex():
    rng = np.random.default_rng(7)
    T = build_rect_matmul_tensor(2, 3, 2)
    Tu = prepare_rect_unfoldings(T)
    N1, N2, N3 = T.shape
    R = 5
    U = rng.normal(size=(N1, R)) + 1j * rng.normal(size=(N1, R))
    V = rng.normal(size=(N2, R)) + 1j * rng.normal(size=(N2, R))
    W = rng.normal(size=(N3, R)) + 1j * rng.normal(size=(N3, R))
    cost, g, H = normal_equations_rect_cp(U, V, W, Tu)
    T_hat = cp_reconstruct(U, V, W)
    c_dir = 0.5 * float(np.sum(np.abs(T_hat - T) ** 2))
    theta = np.concatenate([U.ravel(), V.ravel(), W.ravel()])
    P = theta.size
    J = np.zeros((N1 * N2 * N3, P), dtype=np.complex128)
    eps = 1e-7
    S1, S2 = N1 * R, N2 * R
    for q in range(P):
        th2 = theta.copy(); th2[q] += eps
        U2 = th2[:S1].reshape(N1, R)
        V2 = th2[S1:S1 + S2].reshape(N2, R)
        W2 = th2[S1 + S2:].reshape(N3, R)
        J[:, q] = (cp_reconstruct(U2, V2, W2) - T_hat).ravel() / eps
    g_ref = J.conj().T @ (T_hat - T).ravel()
    H_ref = J.conj().T @ J
    return (
        abs(cost - c_dir),
        float(np.linalg.norm(g - g_ref) / np.linalg.norm(g_ref)),
        float(np.linalg.norm(H - H_ref) / np.linalg.norm(H_ref)),
    )


def sample_init(N1, N2, N3, R, seed, mode="gauss", warm_factors=None):
    rng = np.random.default_rng(seed)
    if warm_factors is not None:
        # Compression warm-start : partir d'une solution exacte R+1, supprimer le terme de plus petite norme
        # ou fusionner deux termes proches, puis perturber dans C !
        U0, V0, W0 = warm_factors
        U0, V0, W0 = balance_factors(U0, V0, W0)
        R0 = U0.shape[1]
        if R0 > R:
            # Supprimer (R0 - R) colonnes : soit les plus petites (tr%2==0), soit aleatoires
            norms = np.linalg.norm(U0, axis=0)
            if seed % 2 == 0:
                keep = np.argsort(norms)[R0 - R:]
            else:
                keep = rng.choice(R0, size=R, replace=False)
            U = U0[:, keep].astype(np.complex128)
            V = V0[:, keep].astype(np.complex128)
            W = W0[:, keep].astype(np.complex128)
        else:
            U, V, W = U0.astype(np.complex128), V0.astype(np.complex128), W0.astype(np.complex128)
        noise = 0.12 if (seed % 3 != 0) else 0.25
        U = U + noise * (rng.normal(size=(N1, R)) + 1j * rng.normal(size=(N1, R))) / math.sqrt(2)
        V = V + noise * (rng.normal(size=(N2, R)) + 1j * rng.normal(size=(N2, R))) / math.sqrt(2)
        W = W + noise * (rng.normal(size=(N3, R)) + 1j * rng.normal(size=(N3, R))) / math.sqrt(2)
        return U, V, W

    if mode == "cyclotomic":
        # Alphabet complexe {0, +-1, +-i, +-w, +-w^2} + bruit doux
        w = np.exp(2j * math.pi / 3.0)
        alph = np.array([0.0, 0.0, 0.0, 1.0, -1.0, 1j, -1j, w, -w, w**2, -w**2], dtype=np.complex128)
        scale = 0.45
        U = scale * rng.choice(alph, size=(N1, R)) + 0.15 * (rng.normal(size=(N1, R)) + 1j * rng.normal(size=(N1, R)))
        V = scale * rng.choice(alph, size=(N2, R)) + 0.15 * (rng.normal(size=(N2, R)) + 1j * rng.normal(size=(N2, R)))
        W = scale * rng.choice(alph, size=(N3, R)) + 0.15 * (rng.normal(size=(N3, R)) + 1j * rng.normal(size=(N3, R)))
        return U, V, W

    # Gaussienne complexe circulaire standard
    sc = 0.52
    U = sc * (rng.normal(size=(N1, R)) + 1j * rng.normal(size=(N1, R))) / math.sqrt(2)
    V = sc * (rng.normal(size=(N2, R)) + 1j * rng.normal(size=(N2, R))) / math.sqrt(2)
    W = sc * (rng.normal(size=(N3, R)) + 1j * rng.normal(size=(N3, R))) / math.sqrt(2)
    return U, V, W


def solve_rect_complex_lm(m, p, n, R, seed, mu_schedule, init_mode="gauss",
                          warm_factors=None, max_norm_stop=5e4):
    T = build_rect_matmul_tensor(m, p, n)
    Tu = prepare_rect_unfoldings(T)
    T_norm = math.sqrt(Tu[3])
    N1, N2, N3 = T.shape
    S1, S2, S3 = N1 * R, N2 * R, N3 * R
    P = S1 + S2 + S3

    U, V, W = sample_init(N1, N2, N3, R, seed, mode=init_mode, warm_factors=warm_factors)
    theta = np.concatenate([U.ravel(), V.ravel(), W.ravel()])

    def to_UVW(th):
        return (
            th[:S1].reshape(N1, R),
            th[S1:S1 + S2].reshape(N2, R),
            th[S1 + S2:].reshape(N3, R),
        )

    def eval_data_cost(th):
        Uc, Vc, Wc = to_UVW(th)
        GU = Uc.conj().T @ Uc; GV = Vc.conj().T @ Vc; GW = Wc.conj().T @ Wc
        GVW = GV * GW
        VW_c = (Vc.conj()[:, None, :] * Wc.conj()[None, :, :]).reshape(N2 * N3, R)
        M1 = Tu[0] @ VW_c
        norm_hat2 = float(np.real(np.sum(GU * GVW)))
        inner = float(np.real(np.sum(Uc.conj() * M1)))
        return 0.5 * max(0.0, norm_hat2 - 2.0 * inner + Tu[3])

    total_it = 0
    best_err = 1e9
    best_th = theta.copy()
    history = []
    mu_end_stats = {}

    for n_st, mu in mu_schedule:
        lam = 1e-2 if mu == 0.0 else max(1e-3, mu)
        for st_it in range(n_st):
            if st_it % 15 == 0 and st_it > 0:
                Uc, Vc, Wc = balance_factors(*to_UVW(theta))
                theta = np.concatenate([Uc.ravel(), Vc.ravel(), Wc.ravel()])
            Uc, Vc, Wc = to_UVW(theta)
            c_data, g, H = normal_equations_rect_cp(Uc, Vc, Wc, Tu)
            err = math.sqrt(2.0 * c_data) / T_norm
            nf = float(np.linalg.norm(theta))
            if err < best_err:
                best_err = err
                best_th = theta.copy()
            if total_it % 40 == 0 or err < 1e-13:
                history.append((total_it, err, nf, mu))
            if err < 1e-13 and mu == 0.0:
                break
            if nf > max_norm_stop:
                break

            c_tot = c_data + 0.5 * mu * float(np.real(np.vdot(theta, theta)))
            if mu > 0.0:
                g = g + mu * theta
                H = H + mu * np.eye(P, dtype=np.complex128)
            diag_H = np.maximum(np.real(np.diag(H)), 1e-6)
            accepted = False
            for _ in range(12):
                try:
                    step = -np.linalg.solve(H + lam * np.diag(diag_H), g)
                except np.linalg.LinAlgError:
                    lam *= 10.0
                    continue
                th_new = theta + step
                c_new = eval_data_cost(th_new) + 0.5 * mu * float(np.real(np.vdot(th_new, th_new)))
                if math.isfinite(c_new) and c_new < c_tot:
                    theta = th_new
                    lam = max(lam * 0.35, 1e-14)
                    accepted = True
                    break
                lam = min(lam * 4.5, 1e12)
            total_it += 1
            if not accepted:
                break
        if mu > 0.0:
            Uc, Vc, Wc = to_UVW(theta)
            e_stage = math.sqrt(2.0 * eval_data_cost(theta)) / T_norm
            mu_end_stats[mu] = (e_stage, float(np.linalg.norm(theta)))

    # Calcul de la pente alpha sur les deux plus petits mu > 0
    mus = sorted(mu_end_stats.keys())
    alpha = None
    if len(mus) >= 2:
        m1, m2 = mus[0], mus[1]
        e1, e2 = mu_end_stats[m1][0], mu_end_stats[m2][0]
        if e1 > 1e-12 and e2 > 1e-12:
            alpha = math.log(e2 / e1) / math.log(m2 / m1)

    Ub, Vb, Wb = balance_factors(*to_UVW(best_th))
    diff = cp_reconstruct(Ub, Vb, Wb) - T
    rel_err = float(np.linalg.norm(diff) / T_norm)
    max_err = float(np.max(np.abs(diff)))
    frob = float(math.sqrt(
        np.sum(np.abs(Ub) ** 2) + np.sum(np.abs(Vb) ** 2) + np.sum(np.abs(Wb) ** 2)
    ))
    max_x = float(max(np.max(np.abs(Ub)), np.max(np.abs(Vb)), np.max(np.abs(Wb))))

    return {
        "seed": seed,
        "init_mode": init_mode if warm_factors is None else "warm_compress",
        "rel_err": rel_err,
        "max_err": max_err,
        "frob_norm": frob,
        "max_entry": max_x,
        "alpha_slope": alpha,
        "mu_end_stats": {str(k): [round(v[0], 6), round(v[1], 2)] for k, v in mu_end_stats.items()},
        "iters": total_it,
        "history": history,
    }, (Ub, Vb, Wb)


def _worker(args):
    m, p, n, R, seed, sched, mode, warm_facs = args
    res, facs = solve_rect_complex_lm(m, p, n, R, seed, sched, init_mode=mode, warm_factors=warm_facs)
    keep = facs if (res["rel_err"] < 2e-2 or warm_facs is None and res["rel_err"] < 1e-6) else None
    return res, keep


def run_group(pool, label, m, p, n, R, n_trials, sched, seed_base, warm_facs=None, save_always_best=False):
    ts = time.perf_counter()
    tasks = []
    for tr in range(n_trials):
        if warm_facs is not None and tr < n_trials // 3:
            tasks.append((m, p, n, R, seed_base + tr, sched, "warm", warm_facs))
        else:
            mode = "gauss" if (tr % 2 == 0) else "cyclotomic"
            tasks.append((m, p, n, R, seed_base + tr, sched, mode, None))
    outs = list(pool.map(_worker, tasks))
    trials = [o[0] for o in outs]
    best_i = int(np.argmin([t["rel_err"] for t in trials]))
    best = trials[best_i]
    best_facs = outs[best_i][1]

    trials_sorted = sorted(trials, key=lambda z: z["rel_err"])
    med = float(np.median([z["rel_err"] for z in trials]))
    n_1e2 = sum(1 for z in trials if z["rel_err"] < 1e-2)
    n_1e6 = sum(1 for z in trials if z["rel_err"] < 1e-6)
    wall = 1.0 / math.sqrt(m * p * n)
    n_below_wall = sum(1 for z in trials if z["rel_err"] < 0.85 * wall)
    dt = time.perf_counter() - ts

    saved = None
    if best_facs is not None and (save_always_best or best["rel_err"] < 1e-2):
        fname = f"rect_{label}.npz"
        saved = os.path.join(HERE, fname)
        np.savez(saved, U=best_facs[0], V=best_facs[1], W=best_facs[2])

    a_s = f"{best['alpha_slope']:.2f}" if best["alpha_slope"] is not None else "  --"
    print(
        f"   {label:<22} | <{m},{p},{n}> R={R:>2} | {n_trials:>2} ess. | <mur:{n_below_wall:>2} | <1e-2:{n_1e2:>2} | <1e-6:{n_1e6:>2} | "
        f"best={best['rel_err']:>9.2e} | ||F||={best['frob_norm']:>6.1f} | max|x|={best['max_entry']:>5.1f} | a={a_s:>4} | {dt:>5.1f}s",
        flush=True,
    )
    return {
        "label": label, "shape": [m, p, n], "R": R, "discrete_wall_1": wall,
        "best": best, "median_rel_err": med,
        "n_below_wall": n_below_wall, "n_below_1e2": n_1e2, "n_below_1e6": n_1e6,
        "trials": trials_sorted, "time_s": dt, "saved_npz": saved,
    }, best_facs


def main():
    t0 = time.time()
    ec, eg, eh = verify_rect_complex()
    print(f"== Sanity check rectangulaire complexe <2,3,2> : ec={ec:.2e}, eg={eg:.2e}, eh={eh:.2e} ==\n", flush=True)
    assert eg < 1e-6 and eh < 1e-6

    n_workers = min(8, os.cpu_count() or 4)
    print(f"== Chasse RECTANGULAIRE sur C ({n_workers} workers, recuit L_{{2/3}} + cyclotomique + warm-start) ==\n", flush=True)

    sched_cal = [(60, 2e-2), (60, 3e-3), (60, 3e-4), (140, 0.0)]
    sched_hunt = [(70, 3e-2), (70, 5e-3), (70, 5e-4), (70, 5e-5), (200, 0.0)]

    report = {"sanity": {"ec": ec, "eg": eg, "eh": eh}, "runs": {}}

    targets = [
        # (name, m, p, n, R_record, R_target, n_cal, n_hunt)
        ("2x3x4", 2, 3, 4, 20, 19, 6, 30),
        ("2x3x5", 2, 3, 5, 25, 24, 6, 24),
        ("2x4x4", 2, 4, 4, 26, 25, 6, 24),
        ("3x3x4", 3, 3, 4, 29, 28, 6, 24),
    ]

    with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as pool:
        for idx, (name, m, p, n, R_rec, R_tgt, n_cal, n_hunt) in enumerate(targets):
            print(f"-- Format <{m}, {p}, {n}> (record mondial = {R_rec}, mur 1/sqrt({m*p*n}) = {1/math.sqrt(m*p*n):.4f}) --", flush=True)
            cal_res, cal_facs = run_group(
                pool, f"Cal_{name}_R{R_rec}", m, p, n, R_rec, n_cal, sched_cal,
                seed_base=10000 + 1000 * idx, warm_facs=None, save_always_best=True
            )
            report["runs"][cal_res["label"]] = cal_res
            warm = cal_facs if (cal_facs is not None and cal_res["best"]["rel_err"] < 1e-6) else None

            hunt_res, _ = run_group(
                pool, f"Hunt_{name}_R{R_tgt}", m, p, n, R_tgt, n_hunt, sched_hunt,
                seed_base=20000 + 1000 * idx, warm_facs=warm, save_always_best=False
            )
            report["runs"][hunt_res["label"]] = hunt_res

    out_json = os.path.join(HERE, "matmul_rect_complex_results.json")
    with open(out_json, "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2)

    # Graphique de synthese
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14.0, 5.0))

    names = [t[0] for t in targets]
    x = np.arange(len(names))
    cal_errs = [max(report["runs"][f"Cal_{t[0]}_R{t[4]}"]["best"]["rel_err"], 1e-14) for t in targets]
    hunt_errs = [max(report["runs"][f"Hunt_{t[0]}_R{t[5]}"]["best"]["rel_err"], 1e-14) for t in targets]
    walls = [1.0 / math.sqrt(t[1] * t[2] * t[3]) for t in targets]

    ax1.semilogy(x - 0.12, cal_errs, "o", color="#2ca02c", markersize=10, label="Calibration au record R_rec (exact C)")
    ax1.semilogy(x + 0.12, hunt_errs, "s", color="#d62728", markersize=10, label="Chasse ouverte R_rec - 1 sur C")
    ax1.semilogy(x, walls, "^--", color="#555555", alpha=0.7, label="Mur discret 1/sqrt(m·p·n)")
    ax1.set_xticks(x)
    ax1.set_xticklabels([f"<{t[1]},{t[2]},{t[3]}>\n({t[4]} → {t[5]})" for t in targets])
    ax1.set_ylabel("Meilleure erreur relative ||T_hat - T|| / ||T|| (log)")
    ax1.set_title("Chasse rectangulaire sur C : Calibration vs Record - 1")
    ax1.grid(alpha=0.3, which="both")
    ax1.legend(fontsize=8.2)

    cols = ["#1f77b4", "#ff7f0e", "#9467bd", "#d62728"]
    mks = ["o", "^", "s", "D"]
    for t, col, mk in zip(targets, cols, mks):
        key = f"Hunt_{t[0]}_R{t[5]}"
        trs = report["runs"][key]["trials"]
        ax2.loglog(
            [z["frob_norm"] for z in trs],
            [max(z["rel_err"], 1e-14) for z in trs],
            mk, color=col, alpha=0.8, markersize=6,
            label=f"<{t[1]},{t[2]},{t[3]}> R={t[5]} (C)"
        )
    ax2.set_xlabel("Norme de Frobenius ||U, V, W||_F (log)")
    ax2.set_ylabel("Erreur relative ||T_hat - T|| / ||T|| (log)")
    ax2.set_title("Signature Exact (norme bornée) vs Bord (norme divergente) sur C")
    ax2.grid(alpha=0.3, which="both")
    ax2.legend(fontsize=8.2)

    fig.tight_layout()
    out_png = os.path.join(OUT, "matmul_rect_complex_hunt.png")
    fig.savefig(out_png, dpi=130)
    print(f"\nTemps total : {time.time() - t0:.1f} s", flush=True)
    print(f"Graphique : {out_png}", flush=True)


if __name__ == "__main__":
    main()
