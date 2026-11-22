"""
Piste 1 (VP =? VNP / GCT) - Sonde profonde sur la frontiere k=6 (perm_3)
========================================================================

Trois experiences calibrees en Parseval-LM (racines de l'unite) :
  1. Taux de succes calibre de l'initialisation diag_r sur k=7 (R) vs
     chasse intensive sur k=6 (24 essais sur R + 16 essais sur C).
  2. Homotopie de bord 7x7 -> 6x6 : part de Grenet 7x7 EXACT (err = 1e-14)
     avec penalite croissante mu * ||(u, v, w)||^2 sur la 7e ligne/colonne
     (sur R et sur C) pour pousser la variete exacte V_7(perm_3) vers le bord det_6 !
  3. Extraction symbolique des 6 monomes de perm_3 et des monomes parasites
     sur le meilleur polynome obtenu a k=6.
"""
from pathlib import Path
import json
import math
import os
import sys
import time

import numpy as np
import torch

from vp_vnp_parseval_lm import (
    perm3_batch,
    sample_sphere,
    build_parseval_batch,
    eval_det_and_cofactor_c,
    degree_profile,
    solve_parseval_lm,
)

HERE = os.path.dirname(os.path.abspath(__file__))
AXE_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = str(AXE_ROOT / "results")
FIGURES_DIR = str(AXE_ROOT / "figures")
WITNESS_DIR = str(AXE_ROOT / "witnesses")
for _d in (RESULTS_DIR, FIGURES_DIR, WITNESS_DIR):
    os.makedirs(_d, exist_ok=True)
OUT = sys.argv[1] if len(sys.argv) > 1 else FIGURES_DIR


def build_grenet_diag_gauge():
    """
    Construit la matrice exacte de Grenet 7x7 en jauge A_0 = diag(1,1,1,1,1,1,0)
    (en permuttant la ligne 0 vers la ligne 6 et en changeant le signe de 6 lignes).
    """
    # Grenet originale :
    # row 0: [x11, x12, x13, 0, 0, 0, 0]
    # row 1: [-1, 0, 0, x22, x23, 0, 0]
    # row 2: [0, -1, 0, x21, 0, x23, 0]
    # row 3: [0, 0, -1, 0, x21, x22, 0]
    # row 4: [0, 0, 0, -1, 0, 0, x33]
    # row 5: [0, 0, 0, 0, -1, 0, x32]
    # row 6: [0, 0, 0, 0, 0, -1, x31]
    # Deplacons row 0 a la fin (6 transpositions -> signe +1) et multiplions rows 0..5 par -1 (signe (-1)^6 = +1) !
    # Alors A_0 = diag(1, 1, 1, 1, 1, 1, 0) !
    A = np.zeros((10, 7, 7), dtype=np.float64)
    for i in range(6):
        A[0, i, i] = 1.0
    # Indices m = 1..9 correspondent a X.reshape(9) :
    # 1:x11, 2:x12, 3:x13, 4:x21, 5:x22, 6:x23, 7:x31, 8:x32, 9:x33
    # Nouvelles lignes 0..5 = - (anciennes lignes 1..6), nouvelle ligne 6 = ancienne ligne 0
    A[5, 0, 3] = -1.0; A[6, 0, 4] = -1.0
    A[4, 1, 3] = -1.0; A[6, 1, 5] = -1.0
    A[4, 2, 4] = -1.0; A[5, 2, 5] = -1.0
    A[9, 3, 6] = -1.0
    A[8, 4, 6] = -1.0
    A[7, 5, 6] = -1.0
    A[1, 6, 0] = 1.0; A[2, 6, 1] = 1.0; A[3, 6, 2] = 1.0
    return A


def homotopy_7_to_6(is_complex=False, seed=0, pivot_idx=0):
    """
    Homotopie de bord V_7(perm_3) -> det_6 :
      On fixe le bloc 6x6 = indices != pivot_idx, et on penalise mu * ||bord(pivot_idx)||^2
      (tous les elements hors-diagonale de la ligne et colonne pivot_idx, plus A_lin[pivot, pivot]).
      Quand le bord -> 0 et A_0[pivot, pivot] = 1, det_7(M) = det_6(L_6) !
    """
    rng = np.random.default_rng(seed)
    A_init = build_grenet_diag_gauge()
    # Permuter pivot_idx avec 0 pour la simplicite : ici pivot_idx in {0..5} a deja A_0[p,p] = 1 !
    # On applique une petite rotation SL_7 aleatoire sur le bloc 6x6 et un bruit complexe si C
    if is_complex:
        A_c = A_init.astype(np.complex128) + 1e-4 * (
            rng.normal(size=(10, 7, 7)) + 1j * rng.normal(size=(10, 7, 7))
        )
        # Garder A_0[pivot, pivot] = 1
        A_c[0, pivot_idx, pivot_idx] = 1.0
    else:
        A_c = (A_init + 1e-4 * rng.normal(size=(10, 7, 7))).astype(np.complex128)
        A_c[0, pivot_idx, pivot_idx] = 1.0

    # Verification initiale
    X_tr, y_tr = sample_sphere(200, is_complex, seed=55000 + seed, radius=2.2)
    feats_tr, targ_tr, _ = build_parseval_batch(X_tr, y_tr, 7, t_radius=1.0)
    X_va, y_va = sample_sphere(320, is_complex, seed=88888, radius=2.2)
    targ_norm = float(torch.sqrt(torch.mean(torch.abs(targ_tr) ** 2)).item())
    N_eq = feats_tr.shape[0]

    # Masque des parametres du "bord" de la 7e dimension (ligne/col pivot sauf A_0[p,p]=1)
    border_mask = np.zeros((10, 7, 7), dtype=np.float64)
    border_mask[:, pivot_idx, :] = 1.0
    border_mask[:, :, pivot_idx] = 1.0
    border_mask[0, pivot_idx, pivot_idx] = 0.0  # on fixe A_0[p,p] = 1.0

    # Parametres libres = tous les (10,7,7) sauf A_0[p,p] (fixe a 1)
    free_mask = np.ones((10, 7, 7), dtype=bool)
    free_mask[0, pivot_idx, pivot_idx] = False
    free_idx = np.where(free_mask.ravel())[0]
    bmask_free = torch.tensor(border_mask.ravel()[free_idx], dtype=torch.float64)

    n_free = len(free_idx)
    if is_complex:
        theta = torch.tensor(
            np.concatenate([A_c.real.ravel()[free_idx], A_c.imag.ravel()[free_idx]]),
            dtype=torch.float64,
        )
        b_weights = torch.cat([bmask_free, bmask_free])
    else:
        theta = torch.tensor(A_c.real.ravel()[free_idx], dtype=torch.float64)
        b_weights = bmask_free

    free_idx_t = torch.tensor(free_idx, dtype=torch.long)
    keep6 = [i for i in range(7) if i != pivot_idx]

    def unpack(th):
        A_flat = torch.zeros(10 * 49, dtype=torch.complex128)
        A_flat[pivot_idx * 7 + pivot_idx] = 1.0 + 0.0j
        if is_complex:
            A_flat[free_idx_t] = torch.complex(th[:n_free], th[n_free:])
        else:
            A_flat[free_idx_t] = th.to(torch.complex128)
        return A_flat.reshape(10, 7, 7)

    # Continuation sur mu (poids de la penalite sur le bord de la 7e dimension)
    mus = [0.0, 1e-4, 1e-3, 1e-2, 5e-2, 2e-1, 1.0, 5.0, 25.0, 150.0, 1000.0]
    trace = []

    for mu in mus:
        lam = 1e-3
        sqrt_mu = math.sqrt(mu)
        w_pen = sqrt_mu * b_weights

        def eval_total_cost(th):
            A = unpack(th)
            L = torch.einsum("bm,mij->bij", feats_tr, A)
            r = torch.linalg.det(L) - targ_tr
            c_data = 0.5 * float(torch.sum(torch.abs(r) ** 2).item())
            c_pen = 0.5 * float(torch.sum((w_pen * th) ** 2).item())
            return c_data + c_pen

        cost = eval_total_cost(theta)
        for _ in range(35):
            A = unpack(theta)
            d, G = eval_det_and_cofactor_c(A, feats_tr)
            r = d - targ_tr
            Jall = (feats_tr[:, :, None, None] * G[:, None, :, :]).reshape(N_eq, 10 * 49)[:, free_idx_t]
            if is_complex:
                J_top = torch.cat([Jall.real, -Jall.imag], dim=1)
                J_bot = torch.cat([Jall.imag,  Jall.real], dim=1)
            else:
                J_top = Jall.real
                J_bot = Jall.imag
            J_data = torch.cat([J_top, J_bot], dim=0)
            r_data = torch.cat([r.real, r.imag], dim=0)

            JTJ = J_data.T @ J_data + torch.diag(w_pen ** 2)
            g = J_data.T @ r_data + (w_pen ** 2) * theta
            diag_H = torch.clamp(torch.diagonal(JTJ), min=1e-6)

            accepted = False
            for _ in range(10):
                H = JTJ + lam * torch.diag(diag_H)
                try:
                    step = -torch.linalg.solve(H, g)
                except RuntimeError:
                    lam *= 10.0
                    continue
                th_new = theta + step
                c_new = eval_total_cost(th_new)
                if math.isfinite(c_new) and c_new < cost:
                    theta = th_new
                    cost = c_new
                    lam = max(lam * 0.35, 1e-13)
                    accepted = True
                    break
                else:
                    lam = min(lam * 4.0, 1e12)
            if not accepted:
                break

        A_now = unpack(theta)
        err7, _ = degree_profile(A_now, X_va, y_va, 7)
        A6 = A_now[:, keep6, :][:, :, keep6]
        err6, _ = degree_profile(A6, X_va, y_va, 6)
        b_norm = float(torch.linalg.norm(b_weights * theta).item())
        n6 = float(torch.linalg.norm(A6).item())
        trace.append({
            "mu": mu, "err7": err7, "err6": err6,
            "border_norm": b_norm, "norm_L6": n6,
        })
    return trace


def extract_degree3_monomials(A_c, k):
    """
    Extrait par regression lineaire exacte (165 monomes de degre 3 en 9 variables)
    le polynome homogene P_3(X) calcule par det_k(L(X)).
    """
    # Enumerer les 165 monomes de degre 3 en 9 variables (x_0 .. x_8)
    monoms = []
    for i in range(9):
        for j in range(i, 9):
            for l in range(j, 9):
                monoms.append((i, j, l))
    assert len(monoms) == 165

    X_va, y_va = sample_sphere(350, True, seed=77777, radius=2.2)
    B = X_va.shape[0]
    M = k + 1
    feats, _, _ = build_parseval_batch(X_va, y_va, k, t_radius=1.0)
    with torch.no_grad():
        L = torch.einsum("bm,mij->bij", feats, A_c)
        vals = torch.linalg.det(L).reshape(B, M)
        p3_vals = (torch.fft.fft(vals, dim=1)[:, 3] / M).numpy()

    X_flat = X_va.reshape(B, 9).numpy()
    Phi = np.empty((B, 165), dtype=np.complex128)
    for idx, (i, j, l) in enumerate(monoms):
        Phi[:, idx] = X_flat[:, i] * X_flat[:, j] * X_flat[:, l]
    coeffs, _, _, _ = np.linalg.lstsq(Phi, p3_vals, rcond=None)
    fit_err = float(np.linalg.norm(Phi @ coeffs - p3_vals) / np.linalg.norm(p3_vals))

    # Monomes exacts de perm_3 : (0,4,8), (0,5,7), (1,3,8), (1,5,6), (2,3,7), (2,4,6)
    perm_set = {(0, 4, 8), (0, 5, 7), (1, 3, 8), (1, 5, 6), (2, 3, 7), (2, 4, 6)}
    perm_coeffs = {}
    parasite_l2 = 0.0
    top_parasites = []
    for idx, m in enumerate(monoms):
        c = coeffs[idx]
        if m in perm_set:
            perm_coeffs[str(m)] = [float(c.real), float(c.imag)]
        else:
            parasite_l2 += abs(c) ** 2
            top_parasites.append((abs(c), m, complex(c)))
    top_parasites.sort(reverse=True)
    return {
        "lstsq_fit_err": fit_err,
        "perm3_coeffs": perm_coeffs,
        "parasite_l2_norm": math.sqrt(parasite_l2),
        "top5_parasites": [(str(m), float(a)) for a, m, _ in top_parasites[:5]],
    }


def main():
    t0 = time.time()
    torch.set_num_threads(min(8, os.cpu_count() or 4))

    print("== 1. Calibration diag_r sur k=7 (R) + Chasse intensive k=6 (R & C) ==\n", flush=True)
    results = {}

    for field, k, n_trials, max_it in [
        ("R", 7, 10, 140),
        ("R", 6, 24, 160),
        ("C", 6, 16, 160),
    ]:
        is_c = (field == "C")
        ts = time.perf_counter()
        trials = []
        n_exact = 0
        for tr in range(n_trials):
            res = solve_parseval_lm(
                k, is_complex=is_c, seed=5000 + 100 * k + (50 if is_c else 0) + tr,
                max_iter=max_it, B_mats=210, init_mode="diag_r", fix_diag=True
            )
            trials.append(res)
            if res["va_err"] < 1e-3:
                n_exact += 1
        trials.sort(key=lambda z: z["va_err"])
        best = trials[0]
        med = float(np.median([z["va_err"] for z in trials]))
        dt = time.perf_counter() - ts
        key = f"{field}_k{k}"
        results[key] = {
            "n_trials": n_trials, "n_below_1e-3": n_exact,
            "best_va_err": best["va_err"], "median_va_err": med,
            "best_norm": best["norm"], "best_seed": best["seed"],
            "all_va_errs": [round(z["va_err"], 6) for z in trials],
            "all_norms": [round(z["norm"], 2) for z in trials],
            "time_s": dt,
        }
        print(f"   {key:<5} | {n_trials:>2} essais diag_r | <1e-3 : {n_exact:>2}/{n_trials} | "
              f"best_va = {best['va_err']:.3e} (||A||={best['norm']:.1f}) | med = {med:.3e} | {dt:.1f}s",
              flush=True)

    print("\n== 2. Homotopie de bord exacte V_7(perm_3) -> det_6 (Grenet -> bord) ==\n", flush=True)
    homotopies = {}
    for field, is_c, p_idx in [("R", False, 0), ("R", False, 3), ("C", True, 0)]:
        ts = time.perf_counter()
        tr = homotopy_7_to_6(is_complex=is_c, seed=42 + p_idx, pivot_idx=p_idx)
        key = f"{field}_piv{p_idx}"
        homotopies[key] = tr
        print(f"   Homotopie {key} ({time.perf_counter()-ts:.1f}s) :", flush=True)
        for row in tr:
            print(f"     mu={row['mu']:>8.1e} | err7={row['err7']:.2e} | err6={row['err6']:.3e} | "
                  f"||bord||={row['border_norm']:.3e} | ||L6||={row['norm_L6']:.2f}", flush=True)
    results["homotopies"] = homotopies

    # 3. Extraction symbolique sur la meilleure solution R_k6 (relance rapide de son seed)
    print("\n== 3. Inspection des monomes de degre 3 (meilleur R_k6) ==\n", flush=True)
    best_seed = results["R_k6"]["best_seed"]
    # Rejouer avec retour de A_best
    rng = np.random.default_rng(best_seed)
    # On evalue directement via une petite fonction d'extraction sur Grenet et sur R_k6 best_seed
    A_grenet = torch.tensor(build_grenet_diag_gauge(), dtype=torch.complex128)
    mon_g = extract_degree3_monomials(A_grenet, 7)
    print(f"   Sanity check monomes Grenet 7x7 : perm3 = {list(mon_g['perm3_coeffs'].values())[:3]}..., "
          f"parasites L2 = {mon_g['parasite_l2_norm']:.2e}", flush=True)

    # Rejouer best_seed pour R_k6 en recuperant theta
    X_tr, y_tr = sample_sphere(210, False, seed=20000 + best_seed, radius=2.2)
    feats_tr, targ_tr, _ = build_parseval_batch(X_tr, y_tr, 6, t_radius=1.0)
    d_re = np.zeros(6); d_re[:5] = 1.0; d_re += 0.02 * rng.normal(size=6)
    a_re = rng.normal(0.0, 0.25 / math.sqrt(6), (9, 6, 6))
    theta = torch.tensor(np.concatenate([d_re, a_re.ravel()]), dtype=torch.float64)
    idx_k = torch.arange(6)
    N_eq = feats_tr.shape[0]
    lam = 1e-2
    def unpack6(th):
        A = torch.zeros((10, 6, 6), dtype=torch.complex128)
        A[0, idx_k, idx_k] = th[:6].to(torch.complex128)
        A[1:] = th[6:].reshape(9, 6, 6).to(torch.complex128)
        return A
    def c6(th):
        L = torch.einsum("bm,mij->bij", feats_tr, unpack6(th))
        return 0.5 * float(torch.sum(torch.abs(torch.linalg.det(L) - targ_tr) ** 2).item())
    cost = c6(theta)
    for _ in range(160):
        A = unpack6(theta)
        d, G = eval_det_and_cofactor_c(A, feats_tr)
        r = d - targ_tr
        J0 = G[:, idx_k, idx_k]
        J1 = (feats_tr[:, 1:, None, None] * G[:, None, :, :]).reshape(N_eq, 9 * 36)
        J = torch.cat([torch.cat([J0.real, J1.real], dim=1), torch.cat([J0.imag, J1.imag], dim=1)], dim=0)
        rv = torch.cat([r.real, r.imag], dim=0)
        JTJ = J.T @ J; g = J.T @ rv
        dH = torch.clamp(torch.diagonal(JTJ), min=1e-6)
        acc = False
        for _ in range(12):
            try:
                st = -torch.linalg.solve(JTJ + lam * torch.diag(dH), g)
            except RuntimeError:
                lam *= 10.0; continue
            cn = c6(theta + st)
            if math.isfinite(cn) and cn < cost:
                theta = theta + st; cost = cn; lam = max(lam * 0.35, 1e-14); acc = True; break
            lam = min(lam * 4.0, 1e12)
        if not acc:
            break
    mon_k6 = extract_degree3_monomials(unpack6(theta), 6)
    results["monomials_best_R_k6"] = mon_k6
    print(f"   Meilleur R_k6 (seed {best_seed}) : coeffs sur les 6 monomes de perm_3 :", flush=True)
    for m, val in mon_k6["perm3_coeffs"].items():
        print(f"      {m} : {val[0]:+.4f}", flush=True)
    print(f"   Norme L2 des 159 monomes parasites : {mon_k6['parasite_l2_norm']:.4f} "
          f"(top 3 : {mon_k6['top5_parasites'][:3]})", flush=True)

    out_json = os.path.join(RESULTS_DIR, "vp_vnp_k6_deep_results.json")
    with open(out_json, "w", encoding="utf-8") as fh:
        json.dump(results, fh, indent=2)

    # Graphique homotopie + distribution des essais
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 4.8))

    for key, col, mk, lbl in [
        ("R_piv0", "#1f77b4", "o-", "Homotopie R (pivot couche 2)"),
        ("R_piv3", "#ff7f0e", "^--", "Homotopie R (pivot couche 3)"),
        ("C_piv0", "#d62728", "s:", "Homotopie C (pivot couche 2)"),
    ]:
        tr = homotopies[key]
        bn = [max(r["border_norm"], 1e-6) for r in tr]
        e7 = [max(r["err7"], 1e-15) for r in tr]
        e6 = [max(r["err6"], 1e-15) for r in tr]
        ax1.loglog(bn, e7, mk, color=col, linewidth=2, markersize=6, label=f"{lbl} : err7")
        ax1.loglog(bn, e6, "--", color=col, alpha=0.5, linewidth=1.4, label=f"{lbl} : err6")
    ax1.invert_xaxis()
    ax1.set_xlabel("Norme ||bord 7e dim|| (de droite à gauche : 7x7 → 6x6)")
    ax1.set_ylabel("Erreur relative de validation Parseval (log)")
    ax1.set_title("Homotopie de bord V_7(perm_3) → det_6 : mur géométrique")
    ax1.grid(alpha=0.3, which="both")
    ax1.legend(fontsize=7.5)

    labels = ["R k=7\n(10 diag_r)", "R k=6\n(24 diag_r)", "C k=6\n(16 diag_r)"]
    keys = ["R_k7", "R_k6", "C_k6"]
    cols = ["#2ca02c", "#1f77b4", "#d62728"]
    for idx, (kname, col) in enumerate(zip(keys, cols)):
        errs = [max(e, 1e-14) for e in results[kname]["all_va_errs"]]
        xs = np.full(len(errs), idx) + np.linspace(-0.18, 0.18, len(errs))
        ax2.semilogy(xs, errs, "o", color=col, alpha=0.8, markersize=6)
    ax2.axhline(math.sqrt(2.0 / 6.0), color="#555555", linestyle="--", alpha=0.7,
                label="Seuil 4/6 monômes = sqrt(2/6) ≈ 0.577")
    ax2.axhline(math.sqrt(1.0 / 6.0), color="#9467bd", linestyle=":", alpha=0.8,
                label="Seuil 5/6 monômes = sqrt(1/6) ≈ 0.408")
    ax2.set_xticks([0, 1, 2])
    ax2.set_xticklabels(labels)
    ax2.set_ylabel("Erreur relative de validation Parseval (log)")
    ax2.set_title("Calibration statistique : taux de succès k=7 vs mur k=6")
    ax2.grid(alpha=0.3, which="both")
    ax2.legend(fontsize=8)

    fig.tight_layout()
    out_png = os.path.join(OUT, "vp_vnp_k6_homotopy.png")
    fig.savefig(out_png, dpi=130)
    print(f"\nTemps total : {time.time() - t0:.1f} s", flush=True)
    print(f"Graphique : {out_png}", flush=True)


if __name__ == "__main__":
    main()
