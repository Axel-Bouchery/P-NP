"""
Piste 3 (Phase 4) : Optimisation continue sous BARRIERE L_inf STRICTE sur C
===========================================================================

Principe mathematique (methode de barriere circulaire L_inf + equilibrage L_inf) :
  Pour fermer completement les vallees de rang de bord (ou ||F|| -> infini)
  sans introduire AUCUN biais sur les vraies solutions exactes bornees :
    L(z) = 0.5 * ||T_hat(z) - T||_F^2 + 0.5 * gamma * sum_q max(0, |z_q| - M)^2.

  Proprietes exactes :
    1. Si une solution exacte satisfait |z_q| <= M partout (tous les schemas
       connus comme Strassen, Laderman, Smirnov, AlphaEvolve satisfont max|z| <= 1.5
       apres equilibrage L_inf), la barriere, son gradient et son Hessien sont
       IDENTIQUEMENT NULS -> convergence quadratique vers 1e-14 sans jamais couper gamma !
    2. Le majorant hermitien H_bar = gamma * diag(1_{|z_q| > M}) est un majorant
       MM exact du vrai Hessien reel 2x2 de 0.5*max(0,|z|-M)^2 (courbures radiale 1,
       tangentielle 1 - M/|z| in (0,1)), conservant la resolution hermitienne P x P
       4x plus rapide sur C !
    3. Equilibrage de jauge L_inf exact : pour chaque terme r, on pose
       s_r = (||u_r||_inf * ||v_r||_inf * ||w_r||_inf)^{1/3}, minimisant exactement
       max_q |z_q| sur la jauge de mise a l'echelle sans modifier T_hat.

Test tournesol (calibration critique) :
  - <3, 3, 4> a R=29 : sans barriere, la Phase 3 etait aspiree par la vallee de bord
    (err=7.3e-4, max|x|=4.0). Avec barriere L_inf, le solveur doit trouver la VRAIE
    solution EXACTE bornee de Smirnov (err < 1e-12, max|x| <= M) !
  - Puis chasse stricte sous barriere L_inf aux records exacts ouverts sur C :
      * <3, 3, 3> a R=22 (M=1.8 et M=2.5, alphabet AlphaEvolve {0, +-1/2, +-i/2, +-1, +-i})
      * <3, 3, 4> a R=28
      * <2, 4, 4> a R=25
      * <2, 3, 5> a R=24
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

from matmul_rect_complex_hunt import (
    build_rect_matmul_tensor,
    prepare_rect_unfoldings,
    normal_equations_rect_cp,
    cp_reconstruct,
)

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = sys.argv[1] if len(sys.argv) > 1 else HERE


def balance_linf(U, V, W):
    """
    Equilibrage de jauge L_inf exact :
      ||u_r||_inf = ||v_r||_inf = ||w_r||_inf = (a_r * b_r * c_r)^(1/3).
    Minimise exactement max(|u_r|, |v_r|, |w_r|) sans modifier T_hat !
    """
    au = np.max(np.abs(U), axis=0)
    av = np.max(np.abs(V), axis=0)
    aw = np.max(np.abs(W), axis=0)
    prod = au * av * aw
    mask = prod > 1e-30
    s = np.ones_like(prod)
    s[mask] = np.cbrt(prod[mask])
    Un = U.copy(); Vn = V.copy(); Wn = W.copy()
    Un[:, mask] *= (s[mask] / au[mask])[None, :]
    Vn[:, mask] *= (s[mask] / av[mask])[None, :]
    Wn[:, mask] *= (s[mask] / aw[mask])[None, :]
    return Un, Vn, Wn


def barrier_cost_grad_diag(theta, M, gamma):
    """
    Barriere circulaire L_inf exacte sur C :
      cost = 0.5 * gamma * sum max(0, |z_q| - M)^2
      grad = gamma * max(0, |z_q| - M) / |z_q| * z_q   (derivee Wirtinger 2*dPhi/dz*)
      diag_H = gamma * 1_{|z_q| > M}                   (majorant hermitien exact)
    """
    if gamma <= 0.0:
        return 0.0, np.zeros_like(theta), np.zeros(theta.shape, dtype=np.float64)
    abs_z = np.abs(theta)
    excess = np.maximum(0.0, abs_z - M)
    active = abs_z > M
    cost = 0.5 * gamma * float(np.sum(excess * excess))
    grad = np.zeros_like(theta)
    if np.any(active):
        grad[active] = gamma * (excess[active] / abs_z[active]) * theta[active]
    diag_h = np.where(active, gamma, 0.0)
    return cost, grad, diag_h


def sample_barrier_init(N1, N2, N3, R, seed, mode="alphaevolve"):
    rng = np.random.default_rng(seed)
    if mode == "alphaevolve":
        # Alphabet complexe d'AlphaEvolve : {0, +-1/2, +-i/2, +-1, +-i} + bruit doux
        alph = np.array([
            0.0, 0.0, 0.0, 0.0,
            0.5, -0.5, 0.5j, -0.5j,
            1.0, -1.0, 1.0j, -1.0j,
        ], dtype=np.complex128)
        U = rng.choice(alph, size=(N1, R)) + 0.12 * (rng.normal(size=(N1, R)) + 1j * rng.normal(size=(N1, R)))
        V = rng.choice(alph, size=(N2, R)) + 0.12 * (rng.normal(size=(N2, R)) + 1j * rng.normal(size=(N2, R)))
        W = rng.choice(alph, size=(N3, R)) + 0.12 * (rng.normal(size=(N3, R)) + 1j * rng.normal(size=(N3, R)))
        return balance_linf(U, V, W)
    elif mode == "cyclotomic":
        w = np.exp(2j * math.pi / 3.0)
        alph = np.array([0.0, 0.0, 0.0, 1.0, -1.0, 1j, -1j, w, -w, w**2, -w**2], dtype=np.complex128)
        U = 0.5 * rng.choice(alph, size=(N1, R)) + 0.12 * (rng.normal(size=(N1, R)) + 1j * rng.normal(size=(N1, R)))
        V = 0.5 * rng.choice(alph, size=(N2, R)) + 0.12 * (rng.normal(size=(N2, R)) + 1j * rng.normal(size=(N2, R)))
        W = 0.5 * rng.choice(alph, size=(N3, R)) + 0.12 * (rng.normal(size=(N3, R)) + 1j * rng.normal(size=(N3, R)))
        return balance_linf(U, V, W)
    else:
        sc = 0.55
        U = sc * (rng.normal(size=(N1, R)) + 1j * rng.normal(size=(N1, R))) / math.sqrt(2)
        V = sc * (rng.normal(size=(N2, R)) + 1j * rng.normal(size=(N2, R))) / math.sqrt(2)
        W = sc * (rng.normal(size=(N3, R)) + 1j * rng.normal(size=(N3, R))) / math.sqrt(2)
        return balance_linf(U, V, W)


def solve_barrier_lm(m, p, n, R, seed, M_bound=2.0, init_mode="alphaevolve",
                     stages=None, is_real=False):
    """
    Levenberg-Marquardt exact sous barriere L_inf stricte permanente |z_q| <= M_bound.
    """
    T = build_rect_matmul_tensor(m, p, n)
    Tu = prepare_rect_unfoldings(T)
    T_norm = math.sqrt(Tu[3])
    N1, N2, N3 = T.shape
    S1, S2, S3 = N1 * R, N2 * R, N3 * R
    P = S1 + S2 + S3

    U, V, W = sample_barrier_init(N1, N2, N3, R, seed, mode=init_mode)
    if is_real:
        U, V, W = U.real.astype(np.float64), V.real.astype(np.float64), W.real.astype(np.float64)
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

    if stages is None:
        # (n_iters, mu_ridge, gamma_barrier)
        # Note : gamma_barrier reste ACTIF jusqu'a la fin !
        stages = [
            (70,  2e-2, 10.0),
            (70,  2e-3, 25.0),
            (80,  0.0,  40.0),
            (220, 0.0,  80.0),
        ]

    idx_diag = np.arange(P)
    total_it = 0
    best_err = 1e9
    best_th = theta.copy()
    history = []

    for n_st, mu, gamma in stages:
        lam = 1e-2 if mu == 0.0 else max(1e-3, mu)
        for st_it in range(n_st):
            # Equilibrage L_inf periodique : reduit max|z| sans changer T_hat !
            if st_it % 12 == 0 and st_it > 0:
                Uc, Vc, Wc = balance_linf(*to_UVW(theta))
                theta = np.concatenate([Uc.ravel(), Vc.ravel(), Wc.ravel()])

            Uc, Vc, Wc = to_UVW(theta)
            c_data, g_data, H = normal_equations_rect_cp(Uc, Vc, Wc, Tu)
            c_bar, g_bar, d_bar = barrier_cost_grad_diag(theta, M_bound, gamma)
            c_rid = 0.5 * mu * float(np.real(np.vdot(theta, theta)))
            c_tot = c_data + c_bar + c_rid

            err = math.sqrt(2.0 * c_data) / T_norm
            max_z = float(np.max(np.abs(theta)))
            nf = float(np.linalg.norm(theta))

            # Ne retenir comme best que si la contrainte barriere est quasi respectee (ou meilleur global)
            if max_z <= M_bound * 1.15 and err < best_err:
                best_err = err
                best_th = theta.copy()

            if total_it % 40 == 0 or (err < 1e-13 and max_z <= M_bound * 1.05):
                history.append((total_it, err, nf, max_z, gamma))
            if err < 1e-13 and mu == 0.0 and max_z <= M_bound * 1.02:
                break

            g = g_data + g_bar + (mu * theta if mu > 0.0 else 0.0)
            H[idx_diag, idx_diag] += d_bar + mu
            diag_H = np.maximum(np.real(H[idx_diag, idx_diag]), 1e-6)

            accepted = False
            for _ in range(12):
                Hd = H.copy()
                Hd[idx_diag, idx_diag] += lam * diag_H
                try:
                    step = -np.linalg.solve(Hd, g)
                except np.linalg.LinAlgError:
                    lam *= 10.0
                    continue
                th_new = theta + step
                cn_data = eval_data_cost(th_new)
                cn_bar, _, _ = barrier_cost_grad_diag(th_new, M_bound, gamma)
                cn_rid = 0.5 * mu * float(np.real(np.vdot(th_new, th_new)))
                c_new = cn_data + cn_bar + cn_rid
                if math.isfinite(c_new) and c_new < c_tot:
                    theta = th_new
                    lam = max(lam * 0.35, 1e-14)
                    accepted = True
                    break
                lam = min(lam * 4.5, 1e12)
            total_it += 1
            if not accepted:
                break

    if best_err > 1e8:
        best_th = theta.copy()

    Ub, Vb, Wb = balance_linf(*to_UVW(best_th))
    diff = cp_reconstruct(Ub, Vb, Wb) - T
    rel_err = float(np.linalg.norm(diff) / T_norm)
    max_err = float(np.max(np.abs(diff)))
    frob = float(math.sqrt(
        np.sum(np.abs(Ub) ** 2) + np.sum(np.abs(Vb) ** 2) + np.sum(np.abs(Wb) ** 2)
    ))
    max_x = float(max(np.max(np.abs(Ub)), np.max(np.abs(Vb)), np.max(np.abs(Wb))))
    n_over_M = int(np.sum(np.abs(np.concatenate([Ub.ravel(), Vb.ravel(), Wb.ravel()])) > M_bound * 1.01))

    return {
        "seed": seed,
        "M_bound": M_bound,
        "init_mode": init_mode,
        "rel_err": rel_err,
        "max_err": max_err,
        "frob_norm": frob,
        "max_entry": max_x,
        "n_over_M": n_over_M,
        "iters": total_it,
        "history": history,
    }, (Ub, Vb, Wb)


def _worker(args):
    m, p, n, R, seed, M_bound, mode, stages, is_real = args
    res, facs = solve_barrier_lm(m, p, n, R, seed, M_bound=M_bound, init_mode=mode, stages=stages, is_real=is_real)
    keep = facs if res["rel_err"] < 5e-2 else None
    return res, keep


def run_barrier_group(pool, label, m, p, n, R, n_trials, M_list, seed_base, stages=None, is_real=False):
    ts = time.perf_counter()
    modes = ["alphaevolve", "cyclotomic", "gauss"]
    tasks = [
        (m, p, n, R, seed_base + tr, M_list[tr % len(M_list)], modes[tr % len(modes)], stages, is_real)
        for tr in range(n_trials)
    ]
    outs = list(pool.map(_worker, tasks))
    trials = [o[0] for o in outs]
    best_i = int(np.argmin([t["rel_err"] for t in trials]))
    best = trials[best_i]
    best_facs = outs[best_i][1]

    trials_sorted = sorted(trials, key=lambda z: z["rel_err"])
    med = float(np.median([z["rel_err"] for z in trials]))
    wall = 1.0 / math.sqrt(m * p * n)
    n_under_wall = sum(1 for z in trials if z["rel_err"] < 0.85 * wall)
    n_1e3 = sum(1 for z in trials if z["rel_err"] < 1e-3)
    n_1e8 = sum(1 for z in trials if z["rel_err"] < 1e-8)
    dt = time.perf_counter() - ts

    saved = None
    if best_facs is not None and best["rel_err"] < 1e-2:
        fname = f"barrier_{label}.npz"
        saved = os.path.join(HERE, fname)
        np.savez(saved, U=best_facs[0], V=best_facs[1], W=best_facs[2])

    print(
        f"   {label:<24} | <{m},{p},{n}> R={R:>2} | {n_trials:>2} ess. | <mur:{n_under_wall:>2} | <1e-3:{n_1e3:>2} | <1e-8:{n_1e8:>2} | "
        f"best={best['rel_err']:>9.2e} | med={med:>8.2e} | ||F||={best['frob_norm']:>5.1f} | max|x|={best['max_entry']:>4.2f} (M={best['M_bound']}) | {dt:>5.1f}s",
        flush=True,
    )
    return {
        "label": label, "shape": [m, p, n], "R": R, "M_list": M_list,
        "discrete_wall": wall, "best": best, "median_rel_err": med,
        "n_under_wall": n_under_wall, "n_under_1e3": n_1e3, "n_under_1e8": n_1e8,
        "trials": trials_sorted, "time_s": dt, "saved_npz": saved,
    }


def main():
    t0 = time.time()
    n_workers = min(8, os.cpu_count() or 4)
    print(f"== Phase 4 : Optimisation sous BARRIERE L_inf STRICTE sur C ({n_workers} workers) ==\n", flush=True)

    stages_cal = [(60, 2e-2, 10.0), (60, 2e-3, 25.0), (80, 0.0, 50.0), (180, 0.0, 100.0)]
    stages_hunt = [(70, 2e-2, 10.0), (70, 2e-3, 25.0), (90, 0.0, 50.0), (240, 0.0, 100.0)]

    report = {"n_workers": n_workers, "runs": {}}

    with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as pool:
        print("-- 1. TEST TOURNESOL (Calibration sous barriere L_inf) : fermer les vallees de bord --", flush=True)
        # Sur <3,3,4> R=29, Phase 3 (sans barriere) s'etait fait pieger en bord (7.3e-4, max|x|=4.0) !
        # La barriere L_inf doit debloquer la vraie solution EXACTE de Smirnov (<1e-12, max|x|<=2) !
        for cfg in [
            ("CalBar_3x3x3_R23", 3, 3, 3, 23, 6,  [1.6, 2.0], 30000, stages_cal),
            ("CalBar_3x3x4_R29", 3, 3, 4, 29, 12, [1.5, 2.0], 31000, stages_cal),
        ]:
            report["runs"][cfg[0]] = run_barrier_group(pool, *cfg)

        print("\n-- 2. CHASSE AUX RECORDS EXACTS SOUS BARRIERE L_inf STRICTE (vallees de bord fermees !) --", flush=True)
        for cfg in [
            ("Bar_3x3x3_R22_C", 3, 3, 3, 22, 32, [1.8, 2.5], 40000, stages_hunt),
            ("Bar_2x3x5_R24_C", 2, 3, 5, 24, 24, [1.8, 2.5], 41000, stages_hunt),
            ("Bar_2x4x4_R25_C", 2, 4, 4, 25, 24, [1.8, 2.5], 42000, stages_hunt),
            ("Bar_3x3x4_R28_C", 3, 3, 4, 28, 24, [1.8, 2.5], 43000, stages_hunt),
        ]:
            report["runs"][cfg[0]] = run_barrier_group(pool, *cfg)

    out_json = os.path.join(HERE, "matmul_barrier_exact_results.json")
    with open(out_json, "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2)

    # Graphique comparatif : Effet de la barriere L_inf (Cal <3,3,4> R=29 + chasse exacte)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14.0, 5.0))

    groups_p1 = [
        ("CalBar_3x3x3_R23", "<3,3,3>\nR=23 (cal)", "#2ca02c"),
        ("CalBar_3x3x4_R29", "<3,3,4>\nR=29 (cal)", "#17becf"),
        ("Bar_3x3x3_R22_C",  "<3,3,3>\nR=22 (C)",   "#d62728"),
        ("Bar_2x3x5_R24_C",  "<2,3,5>\nR=24 (C)",   "#ff7f0e"),
        ("Bar_2x4x4_R25_C",  "<2,4,4>\nR=25 (C)",   "#9467bd"),
        ("Bar_3x3x4_R28_C",  "<3,3,4>\nR=28 (C)",   "#1f77b4"),
    ]
    for idx, (key, lbl, col) in enumerate(groups_p1):
        trs = report["runs"][key]["trials"]
        errs = [max(z["rel_err"], 1e-14) for z in trs]
        xs = np.full(len(errs), idx) + np.linspace(-0.2, 0.2, len(errs))
        ax1.semilogy(xs, errs, "o", color=col, alpha=0.8, markersize=5.5)
        w = report["runs"][key]["discrete_wall"]
        ax1.hlines(w, idx - 0.3, idx + 0.3, colors="#555555", linestyles="--", linewidth=1.2)
    ax1.set_xticks(range(len(groups_p1)))
    ax1.set_xticklabels([g[1] for g in groups_p1], fontsize=8.5)
    ax1.set_ylabel("Erreur relative ||T_hat - T|| / ||T|| (log)")
    ax1.set_title("Barrière L_inf stricte : Calibration exacte vs Chasse R_rec - 1")
    ax1.grid(alpha=0.3, which="both")

    for key, lbl, col in groups_p1:
        trs = report["runs"][key]["trials"]
        ax2.semilogy(
            [z["max_entry"] for z in trs],
            [max(z["rel_err"], 1e-14) for z in trs],
            "o", color=col, label=lbl.replace("\n", " "), alpha=0.8, markersize=6,
        )
    ax2.axvline(1.8, color="#888888", linestyle=":", label="Barrière M=1.8")
    ax2.axvline(2.5, color="#444444", linestyle="--", label="Barrière M=2.5")
    ax2.set_xlabel("Entrée maximale max|z_q| (après équilibrage L_inf)")
    ax2.set_ylabel("Erreur relative (log)")
    ax2.set_title("Verrouillage des vallées de bord : max|z_q| vs Erreur")
    ax2.grid(alpha=0.3, which="both")
    ax2.legend(fontsize=7.5)

    fig.tight_layout()
    out_png = os.path.join(OUT, "matmul_barrier_exact_hunt.png")
    fig.savefig(out_png, dpi=130)
    print(f"\nTemps total : {time.time() - t0:.1f} s", flush=True)
    print(f"Graphique : {out_png}", flush=True)


if __name__ == "__main__":
    main()
