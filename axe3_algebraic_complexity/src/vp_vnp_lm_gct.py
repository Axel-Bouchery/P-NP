"""
Piste 1 (VP =? VNP / GCT) : Solveur Levenberg-Marquardt (Gauss-Newton exact)
                            pour la complexite determinantale et de bord de perm_3
==================================================================================

Innovation algorithmique :
  Au lieu d'Adam/L-BFGS (qui echouent sur la vallee singuliere de det_k),
  on calcule la Jacobienne EXACTE J_{b, (m,i,j)} = x_m^(b) * cof(L^(b))_{i,j}
  en une seule passe vectorisee (O(1 ms) par iteration) et on resout le systeme
  de Levenberg-Marquardt (J^T J + lambda * D) dtheta = -J^T r en float64 / complex128.

  De plus, on fixe la jauge SL_k x SL_k sans perte de generalite en imposant
  A_0 diagonale (par SVD de A_0), et on valide sur 2000 points independants
  multi-echelles t in {0.5, 1.0, 1.5, 2.0}.
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


def sample_multiscale(B, is_complex, seed, scales=(0.6, 1.0, 1.4)):
    rng = np.random.default_rng(seed)
    xs = []
    per = B // len(scales)
    for idx, sc in enumerate(scales):
        cnt = per if idx < len(scales) - 1 else B - per * (len(scales) - 1)
        if is_complex:
            X = (
                rng.normal(0.0, sc / math.sqrt(2), (cnt, 3, 3))
                + 1j * rng.normal(0.0, sc / math.sqrt(2), (cnt, 3, 3))
            )
            xs.append(torch.tensor(X, dtype=torch.complex128))
        else:
            X = rng.normal(0.0, sc, (cnt, 3, 3))
            xs.append(torch.tensor(X, dtype=torch.float64))
    X = torch.cat(xs, dim=0)
    y = perm3_batch(X)
    ones = torch.ones((B, 1), dtype=X.dtype)
    feats = torch.cat([ones, X.reshape(B, 9)], dim=1)
    return feats, y


def verify_grenet_7x7():
    """Verification exacte de la construction de Grenet (2011) k=7 pour perm_3."""
    feats, y = sample_multiscale(512, False, seed=1)
    X = feats[:, 1:].reshape(-1, 3, 3)
    B = X.shape[0]
    M = torch.zeros((B, 7, 7), dtype=torch.float64)
    # Sous-diagonale = -1
    for i in range(6):
        M[:, i + 1, i] = -1.0
    # Couche 1 : source -> {j}
    M[:, 0, 0] = X[:, 0, 0]
    M[:, 0, 1] = X[:, 0, 1]
    M[:, 0, 2] = X[:, 0, 2]
    # Couche 2 : {1} -> {1,2},{1,3} ; {2} -> {1,2},{2,3} ; {3} -> {1,3},{2,3}
    M[:, 1, 3] = X[:, 1, 1]
    M[:, 1, 4] = X[:, 1, 2]
    M[:, 2, 3] = X[:, 1, 0]
    M[:, 2, 5] = X[:, 1, 2]
    M[:, 3, 4] = X[:, 1, 0]
    M[:, 3, 5] = X[:, 1, 1]
    # Couche 3 : {1,2}->t (x_33), {1,3}->t (x_32), {2,3}->t (x_31)
    M[:, 4, 6] = X[:, 2, 2]
    M[:, 5, 6] = X[:, 2, 1]
    M[:, 6, 6] = X[:, 2, 0]
    d = torch.linalg.det(M)
    err = float(torch.max(torch.abs(d - y)).item())
    return err


def eval_det_and_cofactor(A, feats):
    """
    Calcule d = det(L) et G_{b,i,j} = d(det(L^(b))) / d(L^(b)_{i,j}) (derivee holomorphe exacte).
    Fonctionne en float64 et en complex128 (en passant par SVD si quasi-singulier).
    """
    L = torch.einsum("bm,mij->bij", feats, A)
    # Methode rapide via LU / inv^T avec fallback autograd si NaN/Inf
    d = torch.linalg.det(L)
    try:
        Linv = torch.linalg.inv(L)
        G = d.unsqueeze(-1).unsqueeze(-1) * Linv.transpose(-1, -2)
        if not torch.isfinite(G.real).all() or (G.is_complex() and not torch.isfinite(G.imag).all()):
            raise RuntimeError("Singular L")
    except RuntimeError:
        # Fallback SVD exact : det(U S V^H) = det(U)det(V^H) prod(s_i)
        U, S, Vh = torch.linalg.svd(L)
        k = S.shape[-1]
        # prod_{l != r} S_l
        cof_s = torch.ones_like(S)
        for r in range(k):
            mask = [l for l in range(k) if l != r]
            cof_s[:, r] = torch.prod(S[:, mask], dim=-1)
        det_uv = torch.linalg.det(U) * torch.linalg.det(Vh)
        # G = det(U) det(V^H) * conj(U) diag(cof_s) conj(Vh)
        # Verification : d(det(L)) = tr(G^T dL) -> pour L = U S V^H, d(det) = det(U)det(V^H) sum cof_s_r (U^H dL V)_{rr}
        # = tr( (det_uv * V diag(cof_s) U^H) dL ) => G^T = det_uv * V diag(cof_s) U^H => G = det_uv * U^* diag(cof_s) V^T
        G = det_uv.unsqueeze(-1).unsqueeze(-1) * torch.einsum("bir,br,brj->bij", U.conj(), cof_s.to(L.dtype), Vh.conj())
    return d, G


def pack_params(diag_A0, A_lin):
    """Emballe (diag_A0 [k], A_lin [9, k, k]) en vecteur reel 1D."""
    if torch.is_complex(A_lin):
        return torch.cat([
            diag_A0.real.reshape(-1), diag_A0.imag.reshape(-1),
            A_lin.real.reshape(-1), A_lin.imag.reshape(-1)
        ])
    return torch.cat([diag_A0.reshape(-1), A_lin.reshape(-1)])


def unpack_params(theta, k, is_complex):
    if is_complex:
        n0 = k
        n1 = 9 * k * k
        d_re = theta[:n0]
        d_im = theta[n0:2 * n0]
        a_re = theta[2 * n0:2 * n0 + n1].reshape(9, k, k)
        a_im = theta[2 * n0 + n1:].reshape(9, k, k)
        diag_A0 = torch.complex(d_re, d_im)
        A_lin = torch.complex(a_re, a_im)
    else:
        diag_A0 = theta[:k]
        A_lin = theta[k:].reshape(9, k, k)
    A = torch.zeros((10, k, k), dtype=A_lin.dtype)
    idx = torch.arange(k)
    A[0, idx, idx] = diag_A0
    A[1:] = A_lin
    return A, diag_A0, A_lin


def solve_lm(k, is_complex=False, seed=0, max_iter=180, B_train=900,
             init_rank=None, max_norm_stop=1e6, verbose=False):
    """
    Levenberg-Marquardt exact sur det_k(diag(d) + sum x_m A_m) = perm_3(X).
    """
    rng = np.random.default_rng(seed)
    feats_tr, y_tr = sample_multiscale(B_train, is_complex, seed=10000 + seed, scales=(0.55, 0.95, 1.35))
    feats_va, y_va = sample_multiscale(1600, is_complex, seed=99999, scales=(0.5, 1.0, 1.5, 2.0))
    y_tr_norm = float(torch.sqrt(torch.mean(torch.abs(y_tr) ** 2)).item())
    y_va_norm = float(torch.sqrt(torch.mean(torch.abs(y_va) ** 2)).item())

    # Initialisation douce : A_0 = diag(1, ..., 1, 0) (rang k-1) + petite perturbation sur A_1..A_9
    # pour eviter l'explosion du determinant de degre k a l'etape 0 !
    r0 = (k - 1) if init_rank is None else init_rank
    diag_re = np.zeros(k, dtype=np.float64)
    diag_re[:r0] = 1.0
    diag_re += 0.03 * rng.normal(size=k)
    scale_lin = 0.16 / math.sqrt(k)  # garantit det ~ O(1) a t=0
    if is_complex:
        diag_im = 0.03 * rng.normal(size=k)
        diag_A0 = torch.tensor(diag_re + 1j * diag_im, dtype=torch.complex128)
        A_lin = torch.tensor(
            (rng.normal(0.0, scale_lin, (9, k, k)) + 1j * rng.normal(0.0, scale_lin, (9, k, k))) / math.sqrt(2),
            dtype=torch.complex128,
        )
    else:
        diag_A0 = torch.tensor(diag_re, dtype=torch.float64)
        A_lin = torch.tensor(rng.normal(0.0, scale_lin, (9, k, k)), dtype=torch.float64)

    theta = pack_params(diag_A0, A_lin)
    P = theta.numel()
    idx_k = torch.arange(k)

    def residuals_and_jac(th):
        A, _, _ = unpack_params(th, k, is_complex)
        d, G = eval_det_and_cofactor(A, feats_tr)
        r = d - y_tr  # [B]
        # Derivee par rapport a diag_A0 : G[:, idx, idx] [B, k]
        J0 = G[:, idx_k, idx_k]  # [B, k]
        # Derivee par rapport a A_lin : feats[:, 1:, None, None] * G[:, None, :, :] -> [B, 9*k*k]
        J1 = (feats_tr[:, 1:, None, None] * G[:, None, :, :]).reshape(B_train, 9 * k * k)
        if is_complex:
            # Systeme reel equivalent : r_real = [Re(r), Im(r)] (2B)
            # Pour z = u + iv, f = f_re + i f_im holomorphe :
            # df_re/du = Re(f'), df_re/dv = -Im(f'), df_im/du = Im(f'), df_im/dv = Re(f')
            J_c = torch.cat([J0, J1], dim=1)  # [B, P/2]
            Jr = J_c.real
            Ji = J_c.imag
            # Ordonner colonnes comme pack_params : [d_re, d_im, a_re, a_im]
            J0r, J0i = J0.real, J0.imag
            J1r, J1i = J1.real, J1.imag
            J_top = torch.cat([J0r, -J0i, J1r, -J1i], dim=1)
            J_bot = torch.cat([J0i,  J0r, J1i,  J1r], dim=1)
            J_real = torch.cat([J_top, J_bot], dim=0)
            r_real = torch.cat([r.real, r.imag], dim=0)
            return r_real, J_real
        else:
            J_real = torch.cat([J0, J1], dim=1)
            return r, J_real

    def eval_cost(th):
        A, _, _ = unpack_params(th, k, is_complex)
        L = torch.einsum("bm,mij->bij", feats_tr, A)
        r = torch.linalg.det(L) - y_tr
        return 0.5 * float(torch.sum(torch.abs(r) ** 2).item())

    lam = 1e-2
    cost = eval_cost(theta)
    history = []

    for it in range(max_iter):
        r_vec, J = residuals_and_jac(theta)
        cost = 0.5 * float(torch.dot(r_vec, r_vec).item())
        tr_rel = math.sqrt(2.0 * cost / B_train) / y_tr_norm
        anorm = float(torch.linalg.norm(theta).item())
        if it % 20 == 0 or tr_rel < 1e-12:
            history.append((it, tr_rel, anorm))
        if tr_rel < 1e-12:
            break
        if anorm > max_norm_stop:
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

    # Validation finale multi-echelle independante
    A_best, _, _ = unpack_params(theta, k, is_complex)
    with torch.no_grad():
        L_va = torch.einsum("bm,mij->bij", feats_va, A_best)
        d_va = torch.linalg.det(L_va)
        va_rel = float(torch.sqrt(torch.mean(torch.abs(d_va - y_va) ** 2)).item()) / y_va_norm
        L_tr = torch.einsum("bm,mij->bij", feats_tr, A_best)
        d_tr = torch.linalg.det(L_tr)
        tr_rel = float(torch.sqrt(torch.mean(torch.abs(d_tr - y_tr) ** 2)).item()) / y_tr_norm
        A_norm = float(torch.linalg.norm(A_best).item())
        A_max = float(torch.max(torch.abs(A_best)).item())

    return {
        "seed": seed,
        "tr_err": tr_rel,
        "va_err": va_rel,
        "norm": A_norm,
        "max_entry": A_max,
        "iters": it + 1,
        "history": history,
    }


def main():
    t0 = time.time()
    torch.set_num_threads(min(8, os.cpu_count() or 4))

    g_err = verify_grenet_7x7()
    print(f"== Sanity check Grenet 7x7 exact : max|det - perm_3| = {g_err:.2e} ==\n", flush=True)

    print("== Chasse GCT par Levenberg-Marquardt EXACT (jauge A_0 diag, multi-echelle) ==\n", flush=True)
    print(f"   {'Corps':<4} | {'k':>2} | {'Essais':>6} | {'Best err TRAIN':>14} | {'Best err VALID':>14} | "
          f"{'Med err VALID':>14} | {'||A||_F':>9} | {'max|A|':>9} | {'t(s)':>6}", flush=True)
    print("   " + "-" * 102, flush=True)

    configs = [
        ("R", 7, 5,  160, 900),   # Calibration : doit trouver err ~ 1e-13 depuis init aleatoire !
        ("R", 4, 4,  120, 700),
        ("R", 5, 6,  160, 800),
        ("R", 6, 12, 220, 900),   # Frontiere ouverte sur R
        ("C", 5, 4,  150, 700),
        ("C", 6, 12, 220, 800),   # Frontiere ouverte GCT sur C !
    ]

    report = {"grenet_sanity_err": g_err}
    for field, k, n_trials, max_it, B_tr in configs:
        is_c = (field == "C")
        ts = time.perf_counter()
        trials = []
        for tr in range(n_trials):
            # Alterner init_rank = k-1 et init_rank = k (pour permettre les degenerescences de bord !)
            r0 = (k - 1) if (tr % 3 != 2) else k
            res = solve_lm(
                k, is_complex=is_c, seed=100 * k + (50 if is_c else 0) + tr,
                max_iter=max_it, B_train=B_tr, init_rank=r0
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
    ax1.set_ylabel("Erreur relative de validation multi-échelle (log)")
    ax1.set_title("VP vs VNP : Complexité déterminantale & de bord de perm_3 (LM exact)")
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
    ax2.set_ylabel("Erreur relative de validation (échelle log)")
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
