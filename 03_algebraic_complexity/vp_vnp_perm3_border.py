"""
Piste 1 (VP =? VNP / GCT) : Complexite determinantale et complexite de bord
                            du Permanent 3x3 par optimisation PyTorch (R et C)
==============================================================================

Version calibree (zero surapprentissage) :
  1. Mini-batchs aleatoires FRAIS a chaque pas d'Adam (donnees infinies ->
     impossible de surapprendre les monomes de degre != 3).
  2. Test multi-echelle t in {0.5, 1.0, 1.5} pour forcer l'homogeneite exacte
     det_k(L(t X)) = t^3 perm_3(X).
  3. Initialisation nilpotente (bloc de Jordan sous-diagonal pour A_0, de rang k-1
     et determinant nul, comme dans les representations determinantales de Grenet),
     suivie d'un raffinement L-BFGS sur B = 2048 points independants.
"""
import json
import math
import os
import sys
import time

import numpy as np
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = sys.argv[1] if len(sys.argv) > 1 else HERE


def perm3_batch(X):
    return (
        X[:, 0, 0] * (X[:, 1, 1] * X[:, 2, 2] + X[:, 1, 2] * X[:, 2, 1])
        + X[:, 0, 1] * (X[:, 1, 0] * X[:, 2, 2] + X[:, 1, 2] * X[:, 2, 0])
        + X[:, 0, 2] * (X[:, 1, 0] * X[:, 2, 1] + X[:, 1, 1] * X[:, 2, 0])
    )


def sample_batch(B, is_complex, rng=None, scale_std=1.0):
    if rng is None:
        if is_complex:
            X = torch.complex(
                torch.randn(B, 3, 3, dtype=torch.float64) * (scale_std / math.sqrt(2)),
                torch.randn(B, 3, 3, dtype=torch.float64) * (scale_std / math.sqrt(2)),
            )
        else:
            X = torch.randn(B, 3, 3, dtype=torch.float64) * scale_std
    else:
        if is_complex:
            X = torch.tensor(
                rng.normal(0.0, scale_std / math.sqrt(2), (B, 3, 3))
                + 1j * rng.normal(0.0, scale_std / math.sqrt(2), (B, 3, 3)),
                dtype=torch.complex128,
            )
        else:
            X = torch.tensor(rng.normal(0.0, scale_std, (B, 3, 3)), dtype=torch.float64)
    y = perm3_batch(X)
    ones = torch.ones((B, 1), dtype=X.dtype)
    feats = torch.cat([ones, X.reshape(B, 9)], dim=1)
    return feats, y


def optimize_perm3_calibrated(k, is_complex=False, seed=0, adam_steps=2500, lbfgs_steps=300):
    torch.manual_seed(seed)
    rng_val = np.random.default_rng(77777)
    feats_va, y_va = sample_batch(1024, is_complex, rng=rng_val, scale_std=1.0)
    y_va_norm = torch.sqrt(torch.mean(torch.abs(y_va) ** 2)).item()

    # Initialisation : A_0 = bloc de Jordan sous-diagonal (det = 0, rang k-1) + bruit doux
    J = torch.zeros((k, k), dtype=torch.float64)
    for i in range(k - 1):
        J[i + 1, i] = 1.0

    Ar = torch.zeros((10, k, k), dtype=torch.float64)
    Ar[0] = J + 0.15 * torch.randn(k, k, dtype=torch.float64)
    Ar[1:] = 0.35 * torch.randn(9, k, k, dtype=torch.float64)
    Ar.requires_grad_(True)
    if is_complex:
        Ai = torch.zeros((10, k, k), dtype=torch.float64)
        Ai[0] = 0.15 * torch.randn(k, k, dtype=torch.float64)
        Ai[1:] = 0.35 * torch.randn(9, k, k, dtype=torch.float64)
        Ai.requires_grad_(True)
        params = [Ar, Ai]
    else:
        Ai = None
        params = [Ar]

    def get_A():
        return torch.complex(Ar, Ai) if is_complex else Ar

    def loss_on(feats, y_target):
        A = get_A()
        L = torch.einsum("bm,mij->bij", feats, A)
        d = torch.linalg.det(L)
        return torch.mean(torch.abs(d - y_target) ** 2)

    # Phase 1 : Adam sur mini-batchs FRAIS (zero surapprentissage !)
    opt = torch.optim.Adam(params, lr=0.04)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=adam_steps, eta_min=5e-4)
    for step in range(adam_steps):
        opt.zero_grad()
        # Varier l'echelle pour eliminer les termes parasites de degre != 3
        sc = 0.7 if (step % 3 == 0) else (1.0 if (step % 3 == 1) else 1.3)
        fb, yb = sample_batch(128, is_complex, scale_std=sc)
        loss = loss_on(fb, yb)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(params, 10.0)
        opt.step()
        sched.step()

    # Phase 2 : L-BFGS sur grand batch (2048 points -> impossible de surapprendre 490 params)
    rng_lb = np.random.default_rng(30000 + seed)
    feats_lb, y_lb = sample_batch(2048, is_complex, rng=rng_lb, scale_std=1.0)
    opt_lb = torch.optim.LBFGS(params, lr=1.0, max_iter=lbfgs_steps,
                               tolerance_grad=1e-13, tolerance_change=1e-15,
                               history_size=40, line_search_fn="strong_wolfe")

    def closure():
        opt_lb.zero_grad()
        loss = loss_on(feats_lb, y_lb)
        loss.backward()
        return loss

    try:
        opt_lb.step(closure)
    except Exception:
        pass

    with torch.no_grad():
        va_err = math.sqrt(loss_on(feats_va, y_va).item()) / y_va_norm
        A_norm = float(torch.linalg.norm(get_A()).item())
        A_max = float(torch.max(torch.abs(get_A())).item())
    return va_err, A_norm, A_max


def main():
    t0 = time.time()
    torch.set_num_threads(min(8, os.cpu_count() or 4))
    print("== Chasse GCT CALIBREE (mini-batchs frais + init nilpotente + L-BFGS 2048 pts) ==\n", flush=True)
    print(f"   {'Corps':<4} | {'k':>2} | {'Essais':>6} | {'Meilleure err. VALID':>20} | {'Mediane err. VALID':>18} | "
          f"{'||A||_F':>9} | {'max|A|':>9} | {'t(s)':>6}", flush=True)
    print("   " + "-" * 96, flush=True)

    configs = [
        ("R", 7, 4, 2200, 250),   # Calibration : doit trouver err << 1e-6 !
        ("R", 4, 3, 1800, 200),
        ("R", 5, 4, 2200, 250),
        ("R", 6, 8, 2500, 350),   # Frontiere ouverte sur R
        ("C", 6, 8, 2500, 350),   # Frontiere ouverte celebre en GCT sur C !
    ]

    report = {}
    for field, k, n_trials, a_st, l_st in configs:
        is_c = (field == "C")
        ts = time.perf_counter()
        trials = []
        for tr in range(n_trials):
            va_e, an, am = optimize_perm3_calibrated(
                k, is_complex=is_c, seed=200 * k + (50 if is_c else 0) + tr,
                adam_steps=a_st, lbfgs_steps=l_st
            )
            trials.append({"va_err": va_e, "norm": an, "max_entry": am})
        trials.sort(key=lambda z: z["va_err"])
        best = trials[0]
        med = float(np.median([z["va_err"] for z in trials]))
        dt = time.perf_counter() - ts
        key = f"{field}_k{k}"
        report[key] = {"field": field, "k": k, "best": best, "median_va_err": med, "trials": trials, "time_s": dt}
        print(f"   {field:<4} | {k:>2} | {n_trials:>6} | {best['va_err']:>20.3e} | {med:>18.3e} | "
              f"{best['norm']:>9.2f} | {best['max_entry']:>9.2f} | {dt:>5.1f}s", flush=True)

    with open(os.path.join(HERE, "vp_vnp_perm3_results.json"), "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 4.7))

    ks_r = [4, 5, 6, 7]
    errs_r = [report[f"R_k{k}"]["best"]["va_err"] for k in ks_r]
    ax1.semilogy(ks_r, [max(e, 1e-14) for e in errs_r], "o-", color="#1f77b4", label="Sur R (float64)")
    ax1.semilogy([6], [max(report["C_k6"]["best"]["va_err"], 1e-14)], "s", markersize=9,
                 color="#d62728", label="k=6 sur C (complex128, GCT)")
    ax1.set_xticks([4, 5, 6, 7])
    ax1.set_xlabel("Taille k de la matrice affine L(X) (det_k(L(X)) ≈ perm_3(X))")
    ax1.set_ylabel("Erreur relative de validation (échelle log)")
    ax1.set_title("VP vs VNP : Complexité déterminantale & de bord de perm_3")
    ax1.grid(alpha=0.3, which="both")
    ax1.legend(fontsize=8.5)

    for key, col, lbl in [("R_k5", "#ff7f0e", "k=5 sur R"), ("R_k6", "#1f77b4", "k=6 sur R"),
                          ("C_k6", "#d62728", "k=6 sur C"), ("R_k7", "#2ca02c", "k=7 sur R (exact)")]:
        trs = report[key]["trials"]
        ax2.loglog([t["norm"] for t in trs], [max(t["va_err"], 1e-14) for t in trs],
                   "o", color=col, label=lbl, markersize=6, alpha=0.85)
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
