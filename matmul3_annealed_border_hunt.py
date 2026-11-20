"""
Piste 3 (Phase 2) : Chasse parallele par recuit L_{2/3} au rang de bord R=20, 19
                    et detection de pente exacte (alpha ~ 1 vs 3/8) a R=22 sur <3,3,3>
======================================================================================

Decouverte de Phase 1 :
  - Sans recuit L_{2/3} (mu=0), Gauss-Newton reste piege dans les bassins discrets
    sqrt(k/27) (0.19245, 0.27217, 0.33333, 0.38490).
  - Avec recuit L_{2/3}, il franchit ces murs (6/6 exacts a R=23, et chute a 5e-4 a R=22).
  - De plus, la pente alpha = d log(err) / d log(mu) separe mathematiquement :
      * Solution EXACTE bornee : alpha = 1.00  (err ~ mu, ||F|| -> const)
      * Degenerescence de BORD d'ordre e : alpha = 3e / (6e + 2) (= 0.375 pour e=1).

On deploie ici le recuit L_{2/3} en parallele (8 workers, 1 thread BLAS/worker) sur :
  1. Calibration R_bord : R = 21 (Schonhage) et R = 20 (Smirnov, R & C).
  2. Frontiere OUVERTE : R_bord = 19 sur R et C (libre et Z_3-cyclique).
  3. Detecteur de pente alpha sur R = 22 (R & C) pour reperer tout bassin exact borne.
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

from matmul3_tensor_lab import build_matmul_tensor, solve_cp_lm

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = sys.argv[1] if len(sys.argv) > 1 else HERE


def _worker_trial(args):
    R, is_complex, seed, cyclic_FC, sched, init_scale = args
    T3 = build_matmul_tensor(3)
    res, facs = solve_cp_lm(
        T3, R, is_complex=is_complex, seed=seed,
        cyclic_FC=cyclic_FC, mu_schedule=sched,
        max_norm_stop=5e4, init_scale=init_scale
    )
    # Estimer la pente alpha = d log(err) / d log(mu) sur les deux derniers paliers mu > 0
    hist = res["history"]
    mu_errs = {}
    for _, err, norm_f, mu in hist:
        if mu > 0:
            mu_errs[mu] = (err, norm_f)
    mus_sorted = sorted(mu_errs.keys())
    alpha = None
    if len(mus_sorted) >= 2:
        m1, m2 = mus_sorted[0], mus_sorted[1]
        e1, e2 = mu_errs[m1][0], mu_errs[m2][0]
        if e1 > 1e-12 and e2 > 1e-12:
            alpha = math.log(e2 / e1) / math.log(m2 / m1)
    res["alpha_slope"] = alpha
    # Ne renvoyer les facteurs lourds que si erreur tres basse
    keep_facs = facs if res["rel_err"] < 1e-2 else None
    return res, keep_facs


def run_parallel_group(pool, label, R, is_complex, n_trials, cyclic_FC, sched, seed_base, init_scale=0.55):
    ts = time.perf_counter()
    tasks = [
        (R, is_complex, seed_base + tr, cyclic_FC, sched, init_scale)
        for tr in range(n_trials)
    ]
    outs = list(pool.map(_worker_trial, tasks))
    trials = [o[0] for o in outs]
    best_idx = int(np.argmin([t["rel_err"] for t in trials]))
    best = trials[best_idx]
    best_facs = outs[best_idx][1]

    trials_sorted = sorted(trials, key=lambda z: z["rel_err"])
    med = float(np.median([z["rel_err"] for z in trials]))
    n_under_015 = sum(1 for z in trials if z["rel_err"] < 0.15)
    n_under_1e2 = sum(1 for z in trials if z["rel_err"] < 1e-2)
    dt = time.perf_counter() - ts

    saved = None
    if best_facs is not None and best["rel_err"] < 1e-2:
        fname = f"phase2_{label}.npz"
        saved = os.path.join(HERE, fname)
        np.savez(saved, U=best_facs[0], V=best_facs[1], W=best_facs[2])

    a_str = f"{best['alpha_slope']:.2f}" if best["alpha_slope"] is not None else "  --"
    print(
        f"   {label:<26} | R={R:>2} | {n_trials:>2} ess. | <0.15:{n_under_015:>2} | <1e-2:{n_under_1e2:>2} | "
        f"best={best['rel_err']:>9.2e} | med={med:>8.2e} | ||F||={best['frob_norm']:>6.1f} | "
        f"alpha={a_str:>4} | {dt:>5.1f}s",
        flush=True,
    )
    return {
        "label": label, "R": R, "is_complex": is_complex, "cyclic_FC": cyclic_FC,
        "best": best, "median_rel_err": med,
        "n_under_015": n_under_015, "n_under_1e2": n_under_1e2,
        "trials": trials_sorted, "time_s": dt, "saved_npz": saved,
    }


def main():
    t0 = time.time()
    n_workers = min(8, os.cpu_count() or 4)
    print(f"== Chasse PARALLELE ({n_workers} workers) par recuit L_{{2/3}} sur <3,3,3> ==\n", flush=True)

    # Recuit L_{2/3} progressif : casse les bassins sqrt(k/27) puis relache a mu=0 pour polir le bord
    sched_border = [(70, 4e-2), (70, 8e-3), (70, 1.5e-3), (70, 2e-4), (240, 0.0)]
    sched_exact = [(80, 3e-2), (80, 5e-3), (80, 5e-4), (80, 5e-5), (180, 0.0)]

    report = {"n_workers": n_workers, "runs": {}}

    with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as pool:
        print("-- 1. Calibration recuit L_{2/3} sur R_bord = 21 (Schonhage) et R_bord = 20 (Smirnov) --", flush=True)
        for cfg in [
            ("Anneal_R_libre_R21",  21, False, 16, None,   sched_border, 3000),
            ("Anneal_R_libre_R20",  20, False, 24, None,   sched_border, 3100),
            ("Anneal_C_libre_R20",  20, True,  12, None,   sched_border, 3200),
        ]:
            report["runs"][cfg[0]] = run_parallel_group(pool, *cfg)

        print("\n-- 2. FRONTIERE OUVERTE : Rang de bord R_border(<3,3,3>) <= 19 (record = 20) --", flush=True)
        for cfg in [
            ("Anneal_R_libre_R19",     19, False, 24, None,   sched_border, 4000),
            ("Anneal_C_libre_R19",     19, True,  16, None,   sched_border, 4100),
            ("Anneal_R_cyc1+6_R19",    19, False, 24, (1, 6), sched_border, 4200),
            ("Anneal_C_cyc1+6_R19",    19, True,  16, (1, 6), sched_border, 4300),
            ("Anneal_R_cyc4+5_R19",    19, False, 24, (4, 5), sched_border, 4400),
        ]:
            report["runs"][cfg[0]] = run_parallel_group(pool, *cfg)

        print("\n-- 3. Detecteur de pente alpha (exact ~1.0 vs bord ~0.375) a R = 22 --", flush=True)
        for cfg in [
            ("Slope_R_libre_R22",      22, False, 24, None,   sched_exact, 5000),
            ("Slope_C_libre_R22",      22, True,  16, None,   sched_exact, 5100),
        ]:
            report["runs"][cfg[0]] = run_parallel_group(pool, *cfg)

    out_json = os.path.join(HERE, "matmul3_phase2_results.json")
    with open(out_json, "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2)

    # Graphique de synthese Phase 2
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14.0, 5.0))

    # Panel 1 : Meilleure erreur vs R avec recuit L_{2/3} (R=19, 20, 21, 22)
    r_vals = [19, 20, 21, 22]
    best_r = [
        min(report["runs"]["Anneal_R_libre_R19"]["best"]["rel_err"],
            report["runs"]["Anneal_R_cyc1+6_R19"]["best"]["rel_err"],
            report["runs"]["Anneal_R_cyc4+5_R19"]["best"]["rel_err"]),
        report["runs"]["Anneal_R_libre_R20"]["best"]["rel_err"],
        report["runs"]["Anneal_R_libre_R21"]["best"]["rel_err"],
        report["runs"]["Slope_R_libre_R22"]["best"]["rel_err"],
    ]
    best_c = [
        min(report["runs"]["Anneal_C_libre_R19"]["best"]["rel_err"],
            report["runs"]["Anneal_C_cyc1+6_R19"]["best"]["rel_err"]),
        report["runs"]["Anneal_C_libre_R20"]["best"]["rel_err"],
        None,
        report["runs"]["Slope_C_libre_R22"]["best"]["rel_err"],
    ]
    ax1.semilogy(r_vals, [max(e, 1e-14) for e in best_r], "o-", color="#1f77b4", linewidth=2,
                 markersize=8, label="Recuit L_{2/3} sur R")
    c_xs = [r for r, e in zip(r_vals, best_c) if e is not None]
    c_ys = [max(e, 1e-14) for e in best_c if e is not None]
    ax1.semilogy(c_xs, c_ys, "s--", color="#d62728", linewidth=2, markersize=8, label="Recuit L_{2/3} sur C")
    ax1.axhline(1.0 / math.sqrt(27.0), color="#555555", linestyle=":", label="Mur discret 1/sqrt(27) ≈ 0.192")
    ax1.axhline(math.sqrt(2.0 / 27.0), color="#888888", linestyle="--", alpha=0.6, label="Mur discret sqrt(2/27) ≈ 0.272")
    ax1.set_xticks([19, 20, 21, 22])
    ax1.set_xlabel("Rang tensoriel R sur <3, 3, 3>")
    ax1.set_ylabel("Meilleure erreur relative ||T_hat - T|| / ||T|| (log)")
    ax1.set_title("Recuit L_{2/3} : franchissement des murs discrets (R=19..22)")
    ax1.grid(alpha=0.3, which="both")
    ax1.legend(fontsize=8)

    # Panel 2 : Nuage (||F|| vs err) + pente alpha sur tous les essais R=20, 21, 22, 19
    groups = [
        ("Slope_R_libre_R22",  "#2ca02c", "v", "R=22 (R)"),
        ("Slope_C_libre_R22",  "#9467bd", "D", "R=22 (C)"),
        ("Anneal_R_libre_R21", "#17becf", "^", "R=21 (R)"),
        ("Anneal_R_libre_R20", "#ff7f0e", "o", "R=20 (R)"),
        ("Anneal_C_libre_R20", "#e377c2", "P", "R=20 (C)"),
        ("Anneal_R_libre_R19", "#1f77b4", "x", "R=19 (R)"),
        ("Anneal_C_libre_R19", "#d62728", "s", "R=19 (C)"),
    ]
    for key, col, mk, lbl in groups:
        trs = report["runs"][key]["trials"]
        ax2.loglog(
            [t["frob_norm"] for t in trs],
            [max(t["rel_err"], 1e-14) for t in trs],
            mk, color=col, label=lbl, markersize=6, alpha=0.8
        )
    ax2.set_xlabel("Norme ||U, V, W||_F (log)")
    ax2.set_ylabel("Erreur relative (log)")
    ax2.set_title("Signature de bord (err ~ ||F||^{-3e}) vs murs discrets")
    ax2.grid(alpha=0.3, which="both")
    ax2.legend(fontsize=7.8)

    fig.tight_layout()
    out_png = os.path.join(OUT, "matmul3_phase2_annealed.png")
    fig.savefig(out_png, dpi=130)
    print(f"\nTemps total : {time.time() - t0:.1f} s", flush=True)
    print(f"Graphique : {out_png}", flush=True)


if __name__ == "__main__":
    main()
