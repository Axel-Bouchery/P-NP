"""
Piste 1 (VP =? VNP / GCT) : Solveur Levenberg-Marquardt de Parseval (racines de l'unite)
========================================================================================

Principe exact :
  Pour L(t, X) = A_0 + t * sum_{m=1}^9 x_m A_m, on a :
      det_k(L(t, X)) = sum_{d=0}^k t^d P_d(X).
  En evaluant sur les (k+1) racines de l'unite t_j = exp(2i*pi*j/(k+1)), le theoreme
  de Parseval garantit :
      (1/(k+1)) sum_j |det_k(L(t_j, X)) - t_j^3 perm_3(X)|^2
          = |P_3(X) - perm_3(X)|^2 + sum_{d != 3} |P_d(X)|^2.

  Avantages decisifs :
    1. Decouplage orthogonal EXACT de tous les degres d = 0..k (zero compensation !).
    2. Echantillonnage uniforme sur la sphere ||X||_F = r (zero outlier ||X||^7 !).
    3. Jacobienne exacte par cofacteurs vectorises + Levenberg-Marquardt (float64 / complex128).
    4. Validation sur points independants + decomposition spectrale par degre d = 0..k.
"""
from pathlib import Path
import json
import math
import os
import sys
import time

import numpy as np
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
AXE_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = str(AXE_ROOT / "results")
FIGURES_DIR = str(AXE_ROOT / "figures")
WITNESS_DIR = str(AXE_ROOT / "witnesses")
for _d in (RESULTS_DIR, FIGURES_DIR, WITNESS_DIR):
    os.makedirs(_d, exist_ok=True)
OUT = sys.argv[1] if len(sys.argv) > 1 else FIGURES_DIR


def perm3_batch(X):
    return (
        X[:, 0, 0] * (X[:, 1, 1] * X[:, 2, 2] + X[:, 1, 2] * X[:, 2, 1])
        + X[:, 0, 1] * (X[:, 1, 0] * X[:, 2, 2] + X[:, 1, 2] * X[:, 2, 0])
        + X[:, 0, 2] * (X[:, 1, 0] * X[:, 2, 1] + X[:, 1, 1] * X[:, 2, 0])
    )


def sample_sphere(B, is_complex, seed, radius=2.2):
    """Echantillonne B matrices 3x3 sur la sphere ||X||_F = radius."""
    rng = np.random.default_rng(seed)
    if is_complex:
        X = rng.normal(size=(B, 3, 3)) + 1j * rng.normal(size=(B, 3, 3))
    else:
        X = rng.normal(size=(B, 3, 3))
    norms = np.linalg.norm(X.reshape(B, 9), axis=1, keepdims=True).reshape(B, 1, 1)
    X = X * (radius / norms)
    Xt = torch.tensor(X, dtype=torch.complex128)  # toujours complex128 pour racines de l'unite
    y = perm3_batch(Xt)
    return Xt, y


def build_parseval_batch(X, y, k, t_radius=1.0):
    """
    Pour chaque X^(b) (B matrices) et chaque racine t_j = t_radius * exp(2i*pi*j/(k+1)),
    construit les features [B*(k+1), 10] et la cible t_j^3 * perm_3(X^(b)).
    """
    B = X.shape[0]
    M = k + 1
    angles = 2.0 * math.pi * np.arange(M) / M
    ts = torch.tensor(t_radius * np.exp(1j * angles), dtype=torch.complex128)  # [M]
    # feats_x : [B, 9] -> [B, M, 9] avec facteur t_j
    X_flat = X.reshape(B, 9)
    t_grid = ts[None, :, None]  # [1, M, 1]
    X_scaled = X_flat[:, None, :] * t_grid  # [B, M, 9]
    ones = torch.ones((B, M, 1), dtype=torch.complex128)
    feats = torch.cat([ones, X_scaled], dim=2).reshape(B * M, 10)
    targets = (y[:, None] * (ts[None, :] ** 3)).reshape(B * M)
    return feats, targets, ts


def eval_det_and_cofactor_c(A_c, feats_c):
    L = torch.einsum("bm,mij->bij", feats_c, A_c)
    d = torch.linalg.det(L)
    try:
        Linv = torch.linalg.inv(L)
        G = d.unsqueeze(-1).unsqueeze(-1) * Linv.transpose(-1, -2)
        if not (torch.isfinite(G.real).all() and torch.isfinite(G.imag).all()):
            raise RuntimeError("Singular")
    except RuntimeError:
        U, S, Vh = torch.linalg.svd(L)
        k = S.shape[-1]
        cof_s = torch.ones_like(S)
        for r in range(k):
            mask = [l for l in range(k) if l != r]
            cof_s[:, r] = torch.prod(S[:, mask], dim=-1)
        det_uv = torch.linalg.det(U) * torch.linalg.det(Vh)
        G = det_uv.unsqueeze(-1).unsqueeze(-1) * torch.einsum(
            "bir,br,brj->bij", U.conj(), cof_s.to(torch.complex128), Vh.conj()
        )
    return d, G


def unpack_theta(theta, k, is_complex, fix_diag=True):
    if fix_diag:
        if is_complex:
            n0 = k
            n1 = 9 * k * k
            d0 = torch.complex(theta[:n0], theta[n0:2 * n0])
            a1 = torch.complex(theta[2 * n0:2 * n0 + n1], theta[2 * n0 + n1:]).reshape(9, k, k)
        else:
            d0 = theta[:k].to(torch.complex128)
            a1 = theta[k:].reshape(9, k, k).to(torch.complex128)
        A = torch.zeros((10, k, k), dtype=torch.complex128)
        idx = torch.arange(k)
        A[0, idx, idx] = d0
        A[1:] = a1
        return A
    else:
        n_all = 10 * k * k
        if is_complex:
            return torch.complex(theta[:n_all], theta[n_all:]).reshape(10, k, k)
        return theta.reshape(10, k, k).to(torch.complex128)


def degree_profile(A_c, X_va, y_va, k):
    """Decompose det_k(A_0 + t L_1(X)) = sum_{d=0}^k t^d P_d(X) par FFT exacte."""
    B = X_va.shape[0]
    M = k + 1
    feats, _, _ = build_parseval_batch(X_va, y_va, k, t_radius=1.0)
    with torch.no_grad():
        L = torch.einsum("bm,mij->bij", feats, A_c)
        vals = torch.linalg.det(L).reshape(B, M)
        # FFT le long de dim=1 -> coeff de t^d est fft[:, d] / M (attention signe exp(-2ipi d j / M))
        coeffs = torch.fft.fft(vals, dim=1) / M  # [B, M] : coeffs[:, d] = P_d(X)
        y_norm = float(torch.sqrt(torch.mean(torch.abs(y_va) ** 2)).item())
        prof = {}
        for d in range(M):
            if d == 3:
                err_d = float(torch.sqrt(torch.mean(torch.abs(coeffs[:, 3] - y_va) ** 2)).item()) / y_norm
            else:
                err_d = float(torch.sqrt(torch.mean(torch.abs(coeffs[:, d]) ** 2)).item()) / y_norm
            prof[f"d{d}"] = err_d
        total_err = math.sqrt(sum(v * v for v in prof.values()))
        return total_err, prof


def solve_parseval_lm(k, is_complex=False, seed=0, max_iter=250, B_mats=240,
                      init_mode="diag_r", fix_diag=True, max_norm_stop=1e7):
    rng = np.random.default_rng(seed)
    X_tr, y_tr = sample_sphere(B_mats, is_complex, seed=20000 + seed, radius=2.2)
    feats_tr, targ_tr, _ = build_parseval_batch(X_tr, y_tr, k, t_radius=1.0)
    X_va, y_va = sample_sphere(400, is_complex, seed=88888, radius=2.2)
    targ_norm = float(torch.sqrt(torch.mean(torch.abs(targ_tr) ** 2)).item())

    N_eq = feats_tr.shape[0]
    idx_k = torch.arange(k)

    # Initialisation calibree
    scale_lin = 0.25 / math.sqrt(k)
    if fix_diag:
        d_re = np.zeros(k, dtype=np.float64)
        if init_mode == "diag_r":
            d_re[:k - 1] = 1.0
            d_re += 0.02 * rng.normal(size=k)
        elif init_mode == "diag_full":
            d_re[:] = 1.0 + 0.1 * rng.normal(size=k)
        else:
            d_re[:k - 2] = 1.0
            d_re += 0.05 * rng.normal(size=k)
        if is_complex:
            d_im = 0.02 * rng.normal(size=k)
            a_re = rng.normal(0.0, scale_lin / math.sqrt(2), (9, k, k))
            a_im = rng.normal(0.0, scale_lin / math.sqrt(2), (9, k, k))
            theta = torch.tensor(np.concatenate([d_re, d_im, a_re.ravel(), a_im.ravel()]), dtype=torch.float64)
        else:
            a_re = rng.normal(0.0, scale_lin, (9, k, k))
            theta = torch.tensor(np.concatenate([d_re, a_re.ravel()]), dtype=torch.float64)
    else:
        A0 = np.zeros((k, k), dtype=np.float64)
        for i in range(k - 1):
            A0[i + 1, i] = -1.0
        A0 += 0.02 * rng.normal(size=(k, k))
        A1 = rng.normal(0.0, scale_lin, (9, k, k))
        Ar = np.concatenate([A0[None], A1], axis=0)
        if is_complex:
            Ai = rng.normal(0.0, scale_lin * 0.5, (10, k, k))
            theta = torch.tensor(np.concatenate([Ar.ravel(), Ai.ravel()]), dtype=torch.float64)
        else:
            theta = torch.tensor(Ar.ravel(), dtype=torch.float64)

    def residuals_and_jac(th):
        A_c = unpack_theta(th, k, is_complex, fix_diag=fix_diag)
        d, G = eval_det_and_cofactor_c(A_c, feats_tr)
        r = d - targ_tr  # [N_eq] complex
        if fix_diag:
            J0 = G[:, idx_k, idx_k]  # [N_eq, k]
            J1 = (feats_tr[:, 1:, None, None] * G[:, None, :, :]).reshape(N_eq, 9 * k * k)
            if is_complex:
                J_top = torch.cat([J0.real, -J0.imag, J1.real, -J1.imag], dim=1)
                J_bot = torch.cat([J0.imag,  J0.real, J1.imag,  J1.real], dim=1)
            else:
                J_top = torch.cat([J0.real, J1.real], dim=1)
                J_bot = torch.cat([J0.imag, J1.imag], dim=1)
        else:
            Jall = (feats_tr[:, :, None, None] * G[:, None, :, :]).reshape(N_eq, 10 * k * k)
            if is_complex:
                J_top = torch.cat([Jall.real, -Jall.imag], dim=1)
                J_bot = torch.cat([Jall.imag,  Jall.real], dim=1)
            else:
                J_top = Jall.real
                J_bot = Jall.imag
        J_real = torch.cat([J_top, J_bot], dim=0)
        r_real = torch.cat([r.real, r.imag], dim=0)
        return r_real, J_real

    def eval_cost(th):
        A_c = unpack_theta(th, k, is_complex, fix_diag=fix_diag)
        L = torch.einsum("bm,mij->bij", feats_tr, A_c)
        r = torch.linalg.det(L) - targ_tr
        return 0.5 * float(torch.sum(torch.abs(r) ** 2).item())

    lam = 1e-2
    cost = eval_cost(theta)
    history = []

    for it in range(max_iter):
        r_vec, J = residuals_and_jac(theta)
        cost = 0.5 * float(torch.dot(r_vec, r_vec).item())
        tr_rel = math.sqrt(2.0 * cost / N_eq) / targ_norm
        anorm = float(torch.linalg.norm(theta).item())
        if it % 25 == 0 or tr_rel < 1e-12:
            history.append((it, tr_rel, anorm))
        if tr_rel < 1e-12 or anorm > max_norm_stop:
            break

        JTJ = J.T @ J
        g = J.T @ r_vec
        diag_JTJ = torch.clamp(torch.diagonal(JTJ), min=1e-6)

        accepted = False
        for _ in range(12):
            H = JTJ + lam * torch.diag(diag_JTJ)
            try:
                step = -torch.linalg.solve(H, g)
            except RuntimeError:
                lam *= 10.0
                continue
            th_new = theta + step
            c_new = eval_cost(th_new)
            if math.isfinite(c_new) and c_new < cost:
                theta = th_new
                cost = c_new
                lam = max(lam * 0.35, 1e-14)
                accepted = True
                break
            else:
                lam = min(lam * 4.0, 1e12)
        if not accepted:
            break

    A_best = unpack_theta(theta, k, is_complex, fix_diag=fix_diag)
    va_rel, prof = degree_profile(A_best, X_va, y_va, k)
    tr_rel = math.sqrt(2.0 * cost / N_eq) / targ_norm
    A_norm = float(torch.linalg.norm(A_best).item())
    A_max = float(torch.max(torch.abs(A_best)).item())
    return {
        "seed": seed,
        "init_mode": init_mode,
        "tr_err": tr_rel,
        "va_err": va_rel,
        "deg_profile": prof,
        "norm": A_norm,
        "max_entry": A_max,
        "iters": it + 1,
        "history": history,
    }


def main():
    t0 = time.time()
    torch.set_num_threads(min(8, os.cpu_count() or 4))

    print("== Chasse GCT par Levenberg-Marquardt de PARSEVAL (racines de l'unite) ==\n", flush=True)
    print(f"   {'Corps':<4} | {'k':>2} | {'Essais':>6} | {'Best err TRAIN':>14} | {'Best err VALID':>14} | "
          f"{'Med err VALID':>14} | {'||A||_F':>9} | {'max|A|':>9} | {'t(s)':>6}", flush=True)
    print("   " + "-" * 102, flush=True)

    configs = [
        ("R", 7, 8,  280, 220),   # Calibration : doit trouver err ~ 1e-13 sur R !
        ("R", 4, 4,  180, 200),
        ("R", 5, 6,  220, 220),
        ("R", 6, 12, 320, 240),   # Frontiere ouverte sur R
        ("C", 5, 6,  220, 220),
        ("C", 6, 12, 350, 240),   # Frontiere ouverte GCT sur C !
    ]

    report = {}
    modes = ["diag_r", "diag_full", "diag_r2"]
    for field, k, n_trials, max_it, B_m in configs:
        is_c = (field == "C")
        ts = time.perf_counter()
        trials = []
        for tr in range(n_trials):
            res = solve_parseval_lm(
                k, is_complex=is_c, seed=100 * k + (50 if is_c else 0) + tr,
                max_iter=max_it, B_mats=B_m, init_mode=modes[tr % len(modes)], fix_diag=True
            )
            trials.append(res)
        trials.sort(key=lambda z: z["va_err"])
        best = trials[0]
        med = float(np.median([z["va_err"] for z in trials]))
        dt = time.perf_counter() - ts
        key = f"{field}_k{k}"
        report[key] = {
            "field": field, "k": k, "best": best,
            "median_va_err": med, "trials": trials, "time_s": dt,
        }
        print(f"   {field:<4} | {k:>2} | {n_trials:>6} | {best['tr_err']:>14.3e} | {best['va_err']:>14.3e} | "
              f"{med:>14.3e} | {best['norm']:>9.2f} | {best['max_entry']:>9.2f} | {dt:>5.1f}s", flush=True)
        prof_str = " ".join(f"{dk}={dv:.1e}" for dk, dv in best["deg_profile"].items())
        print(f"       -> Profil degres (best {key}) : {prof_str}", flush=True)

    out_json = os.path.join(RESULTS_DIR, "vp_vnp_perm3_results.json")
    with open(out_json, "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 4.8))

    ks_r = [4, 5, 6, 7]
    errs_r = [max(report[f"R_k{k}"]["best"]["va_err"], 1e-15) for k in ks_r]
    ax1.semilogy(ks_r, errs_r, "o-", color="#1f77b4", linewidth=2, markersize=8, label="Sur R (float64)")
    ks_c = [5, 6]
    errs_c = [max(report[f"C_k{k}"]["best"]["va_err"], 1e-15) for k in ks_c]
    ax1.semilogy(ks_c, errs_c, "s--", color="#d62728", linewidth=2, markersize=8, label="Sur C (complex128, GCT)")
    ax1.axvline(7, color="#2ca02c", linestyle=":", alpha=0.7, label="dc(perm_3) = 7 (Alper-Bogart-Velasco)")
    ax1.set_xticks([4, 5, 6, 7])
    ax1.set_xlabel("Taille k de la matrice affine L(X)")
    ax1.set_ylabel("Erreur relative de validation Parseval (log)")
    ax1.set_title("VP vs VNP : Complexité déterminantale & de bord de perm_3 (Parseval-LM)")
    ax1.grid(alpha=0.3, which="both")
    ax1.legend(fontsize=8.5)

    for key, col, mk, lbl in [
        ("R_k5", "#ff7f0e", "^", "k=5 sur R"),
        ("R_k6", "#1f77b4", "o", "k=6 sur R"),
        ("C_k6", "#d62728", "s", "k=6 sur C"),
        ("R_k7", "#2ca02c", "*", "k=7 sur R (exact)"),
    ]:
        trs = report[key]["trials"]
        ax2.loglog(
            [t["norm"] for t in trs],
            [max(t["va_err"], 1e-15) for t in trs],
            mk, color=col, label=lbl, markersize=7, alpha=0.85
        )
    ax2.set_xlabel("Norme ||A||_F des paramètres (échelle log)")
    ax2.set_ylabel("Erreur relative de validation Parseval (échelle log)")
    ax2.set_title("Signature de bord (GCT) : ||A||_F vs Erreur")
    ax2.grid(alpha=0.3, which="both")
    ax2.legend(fontsize=8.5)

    fig.tight_layout()
    out_png = os.path.join(OUT, "vp_vnp_perm3_border.png")
    fig.savefig(out_png, dpi=130)
    print(f"\nTemps total : {time.time() - t0:.1f} s", flush=True)
    print(f"Graphique : {out_png}", flush=True)


if __name__ == "__main__":
    main()
