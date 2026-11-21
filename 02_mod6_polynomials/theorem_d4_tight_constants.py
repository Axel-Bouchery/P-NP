r"""
Axe II — Theoreme D4 (Resserrement des Constantes du Theoreme D3 par Couplage B_2 / F_p) :
1. En degre d = 2 modulo 2p :
   - Borne generale : 4p - 2 <= N_gen^{(2p)}(2, OR) <= 8p - 7 (ratio < 2.00x pour TOUT p >= 3 !),
     et <= 16 pour p = 3 (avec |w_1|_H <= 3 sur m >= 9 prouve exhaustivement).
   - Pour THR_{n, t} : N_gen^{(2p)}(2, THR_t) <= 8p - 8 + t (<= 34 pour THR_{n,2} mod 10).
   - Pour toute famille avec au moins un c_i = 0 (dont Type A) : N_gen <= 6p - 4 (ratio < 1.50x !).
   - Pour THR_{n, 2} avec G_p != K_n (au moins une non-arete mod p) : N_gen <= 6p - 3 (<= 27 sur mod 10).
2. A tout degre d >= 2 modulo 2p :
   - Partition par la (d-1)-ieme derivee discrete lineaire Delta_{w_1,...,w_{d-1}} P_2(x) in F_2
     sautant l'etape terminale j = d :
     C_d^{D4} = (3d + 5) * 2^{d-2} - (d + 1)  (vs C_d^{D3} = (d + 1)(2^d - 1), gain -> -25%).
3. Fermeture exhaustive exacte du degre mixte (1, 2) sur mod 6 : N_gen^{(6)}((1, 2), OR) = 7 = 3p - 2.
4. Amelioration directe sur mod 70 = 2 x 5 x 7 (R = 34 au lieu de 37 : 4471 -> 3796).
"""

import json
import math
import shutil
import time
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np


def d0_step(d: int, j: int) -> int:
    return sum((d - r + 1) * math.comb(j, r - 1) for r in range(1, min(d, j + 1) + 1))


def compute_d4_bound_2p(p: int, d: int = 2) -> int:
    """
    Calcule la borne exacte du Theoreme D4 en degre d >= 2 sur m = 2p :
    - Etapes j = 0 .. d - 2 : on extrait d - 1 blocs w_1, ..., w_{d-1} par Chevalley-Warning sur F_p
      (+ Lemme D3a d'inversion de Mobius des que 2 <= D_0(j) <= p - 1).
    - Etape finale j = d - 1 : sur les variables restantes U = U_0 U U_1 partitionnees par la
      forme lineaire L(x) = Delta_{w_1,...,w_{d-1}} P_2(x) in F_2, les k = |U_0| + floor(|U_1|/2)
      blocs dans ker(L) satisfont k <= D_0(d - 1) * (p - 1), d'ou |U| <= 2 * D_0(d - 1) * (p - 1) + 1.
    """
    cum_w = 0
    max_n_needed = 0
    for j in range(d - 1):
        d0 = d0_step(d, j)
        if 2 <= d0 <= p - 1:
            m_req = d0 * (p - 1) + 2
            w_used = d0 * (p - 1)
        else:
            m_req = d0 * (p - 1) + 1
            w_used = d0 * (p - 1) + 1
        if cum_w + m_req > max_n_needed:
            max_n_needed = cum_w + m_req
        cum_w += w_used

    d0_last = d0_step(d, d - 1)
    # |U| <= 2 * d0_last * (p - 1) + 1 => contradiction des que |U| >= 2 * d0_last * (p - 1) + 2
    u_contradiction = 2 * d0_last * (p - 1) + 2
    if cum_w + u_contradiction > max_n_needed:
        max_n_needed = cum_w + u_contradiction
    return max_n_needed - 1


def compute_d3_bound_2p(p: int, d: int = 2) -> int:
    R = d + 1
    cum_w = 0
    max_n_needed = 0
    for j in range(R - 1):
        d0 = d0_step(d, j)
        if 2 <= d0 <= p - 1:
            m_req = d0 * (p - 1) + 2
            w_used = d0 * (p - 1)
        else:
            m_req = d0 * (p - 1) + 1
            w_used = d0 * (p - 1) + 1
        if cum_w + m_req > max_n_needed:
            max_n_needed = cum_w + m_req
        cum_w += w_used
    m_last = d0_step(d, R - 1) * (p - 1) + 1
    if cum_w + m_last > max_n_needed:
        max_n_needed = cum_w + m_last
    return max_n_needed - 1


def compute_pq_from_R(q: int, R: int) -> int:
    cum_w = 0
    max_n = 0
    for j in range(R - 1):
        d0 = j + 2
        if 2 <= d0 <= q - 1:
            m_req = d0 * (q - 1) + 2
            w_used = d0 * (q - 1)
        else:
            m_req = d0 * (q - 1) + 1
            w_used = d0 * (q - 1) + 1
        if cum_w + m_req > max_n:
            max_n = cum_w + m_req
        cum_w += w_used
    m_last = (R + 1) * (q - 1) + 1
    if cum_w + m_last > max_n:
        max_n = cum_w + m_last
    return max_n - 1


def make_figure(summary: dict, fig_path: Path, artifact_fig_path: Path):
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.titleweight": "bold",
        "axes.labelsize": 10,
        "figure.facecolor": "white",
        "axes.facecolor": "#f8fafc",
        "axes.edgecolor": "#cbd5e1",
        "grid.color": "#e2e8f0",
        "grid.linestyle": "--",
        "grid.alpha": 0.8,
    })

    fig, axes = plt.subplots(2, 2, figsize=(14.8, 10.6), dpi=180)
    fig.suptitle(
        "Axe II — Théorème D4 : Resserrement des Constantes par Couplage Bilinéaire $(B_2, B_p)$\n"
        "$N_{\\mathrm{gen}}^{(2p)}(2, \\mathrm{OR}) \\leq 8p - 7$ (ratio $< 2{,}00\\times$), Type A $\\leq 6p - 4$ ($< 1{,}50\\times$) et gain $-25\\%$ à tout degré $d$",
        fontsize=13,
        fontweight="bold",
        y=0.98,
        color="#0f172a",
    )

    # Panel (a): Ratio Upper / (4p - 2) showing D3 (2.25x) -> D4 (<2.00x) -> Type A (<1.50x) -> Lower (1.00x)
    ax = axes[0, 0]
    p_vals = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 53, 101]
    r_d3 = [compute_d3_bound_2p(p, 2) / (4 * p - 2) for p in p_vals]
    r_d4 = [(16 if p == 3 else (8 * p - 7)) / (4 * p - 2) for p in p_vals]
    r_typeA = [(6 * p - 4) / (4 * p - 2) for p in p_vals]

    ax.plot(p_vals, r_d3, color="#94a3b8", lw=2.0, ls="--", marker="s", markersize=4.5,
            label="Théorème D3 : $9p - 9$ (asymptote $2{,}25\\times$)")
    ax.plot(p_vals, r_d4, color="#dc2626", lw=2.6, marker="D", markersize=5.5, zorder=5,
            label="Théorème D4 (général) : $8p - 7$ (STRICTEMENT $< 2{,}00\\times$)")
    ax.plot(p_vals, r_typeA, color="#2563eb", lw=2.3, marker="o", markersize=5, zorder=5,
            label="Théorème D4 ($\\exists\\,c_i=0$ / Type A) : $6p - 4$ (STRICTEMENT $< 1{,}50\\times$)")
    ax.axhline(2.00, color="#991b1b", lw=1.4, ls=":", label="Barrière symbolique $2{,}00\\times$")
    ax.axhline(1.50, color="#1d4ed8", lw=1.3, ls=":", label="Asymptote Type A : $1{,}50\\times$")
    ax.axhline(1.00, color="#059669", lw=1.8, ls="--", label="Borne inf. constructive $4p - 2$ ($1{,}00\\times$)")
    ax.fill_between(p_vals, [1.0] * len(p_vals), r_d4, color="#ef4444", alpha=0.09)
    ax.fill_between(p_vals, [1.0] * len(p_vals), r_typeA, color="#3b82f6", alpha=0.11)

    for pm in [3, 5, 11, 23, 101]:
        idx = p_vals.index(pm)
        ax.annotate(f"{r_d4[idx]:.2f}x", (pm, r_d4[idx]),
                    textcoords="offset points", xytext=(0, 7),
                    ha="center", fontsize=8.0, fontweight="bold", color="#991b1b")

    ax.set_title("(a) Franchissement du facteur $2{,}00\\times$ sur $N_{\\mathrm{gen}}^{(2p)}(2, \\mathrm{OR}) / (4p - 2)$")
    ax.set_xlabel("Nombre premier impair $p$ (échelle log)")
    ax.set_ylabel("Ratio Borne Supérieure / $(4p - 2)$")
    ax.set_xscale("log")
    ax.set_ylim(0.88, 2.42)
    ax.grid(True, which="both", alpha=0.45)
    ax.legend(loc="lower right", fontsize=7.6, framealpha=0.95)

    # Panel (b): Concrete bounds on mod 6, 10, 14, 22, 26 (D1 -> D2 -> D3 -> D4 -> Type A -> LB)
    ax = axes[0, 1]
    p_bar = [3, 5, 7, 11]
    labels_b = [f"mod {2*p}\n(p={p})" for p in p_bar]
    x_b = np.arange(len(p_bar))
    w = 0.16

    v_lb = [4 * p - 2 for p in p_bar]
    v_ta = [6 * p - 4 for p in p_bar]
    v_d4 = [16 if p == 3 else (8 * p - 7) for p in p_bar]
    v_d3 = [compute_d3_bound_2p(p, 2) for p in p_bar]
    v_d1 = [2 * p * p + p - 2 for p in p_bar]

    b_lb = ax.bar(x_b - 2 * w, v_lb, w, color="#059669", label="Borne inf. $4p - 2$")
    b_ta = ax.bar(x_b - 1 * w, v_ta, w, color="#2563eb", label="D4 ($\\exists c_i=0$) : $6p - 4$")
    b_d4 = ax.bar(x_b + 0 * w, v_d4, w, color="#dc2626", label="D4 (général) : $8p - 7$")
    b_d3 = ax.bar(x_b + 1 * w, v_d3, w, color="#f59e0b", label="D3 : $9p - 9$")
    b_d1 = ax.bar(x_b + 2 * w, v_d1, w, color="#94a3b8", label="D1 : $2p^2 + p - 2$")

    for bars, col in [(b_lb, "#065f46"), (b_ta, "#1e3a8a"), (b_d4, "#991b1b"), (b_d1, "#475569")]:
        for bar in bars:
            h = int(bar.get_height())
            ax.annotate(f"{h}", (bar.get_x() + bar.get_width() / 2, h),
                        textcoords="offset points", xytext=(0, 3),
                        ha="center", fontsize=7.2, fontweight="bold", color=col)

    ax.set_title("(b) Progression des bornes concrètes D1 $\\to$ D3 $\\to$ D4 sur `mod 6, 10, 14, 22`")
    ax.set_xticks(x_b)
    ax.set_xticklabels(labels_b)
    ax.set_ylabel("Nombre maximal de variables $n$")
    ax.set_ylim(0, 285)
    ax.grid(True, axis="y")
    ax.legend(loc="upper left", fontsize=8.0, framealpha=0.95)

    # Panel (c): Higher-degree linear constant C_d^{D4} vs C_d^{D3} (d = 2..7)
    ax = axes[1, 0]
    degs = [2, 3, 4, 5, 6, 7]
    c_d3 = [(d + 1) * ((1 << d) - 1) for d in degs]
    c_d4 = [(3 * d + 5) * (1 << (d - 2)) - (d + 1) for d in degs]
    gain_pct = [100.0 * (a - b) / a for a, b in zip(c_d3, c_d4)]

    ax.plot(degs, c_d3, color="#94a3b8", lw=2.2, ls="--", marker="s", markersize=5.5,
            label="Théorème D3 : $C_d^{\\mathrm{D3}} = (d+1)(2^d - 1)$")
    ax.plot(degs, c_d4, color="#7c3aed", lw=2.6, marker="o", markersize=6,
            label="Théorème D4 : $C_d^{\\mathrm{D4}} = (3d+5)2^{d-2} - (d+1)$")
    ax.fill_between(degs, c_d4, c_d3, color="#8b5cf6", alpha=0.16,
                    label="Réduction relative $\\Delta C_d \\to -25\\%$")

    for d_k, y4, g in zip(degs, c_d4, gain_pct):
        ax.annotate(f"{y4}\n(-{g:.1f}%)", (d_k, y4),
                    textcoords="offset points", xytext=(0, -24),
                    ha="center", fontsize=7.6, fontweight="bold", color="#5b21b6")

    ax.set_title("(c) Constante linéaire $C_d$ en $p$ pour $N_{\\mathrm{gen}}^{(2p)}(d, \\mathrm{OR}) \\leq C_d(p-1) + 1$")
    ax.set_xlabel("Degré polynomial $d$")
    ax.set_ylabel("Constante linéaire $C_d$ (échelle log)")
    ax.set_yscale("log")
    ax.set_ylim(5, 1600)
    ax.set_xticks(degs)
    ax.grid(True, which="both", alpha=0.45)
    ax.legend(loc="upper left", fontsize=8.2, framealpha=0.95)

    # Panel (d): Exact closure of mixed degree (1, 2) on mod 6 (N_gen = 7 exact)
    ax = axes[1, 1]
    k_vals = [2, 3, 4, 5, 6, 7, 8]
    cls_counts = [6, 28, 119, 272, 109, 13, 0]
    colors_k = ["#3b82f6"] * 5 + ["#059669", "#dc2626"]
    bars = ax.bar(k_vals, cls_counts, color=colors_k, alpha=0.88, width=0.62)

    for k, h in zip(k_vals, cls_counts):
        lbl = f"{h}"
        if k == 7:
            lbl = "13\n(MAX EXACT)"
        elif k == 8:
            lbl = "0 / 4 992\n(EXTINCTION)"
        ax.annotate(lbl, (k, max(h, 2)),
                    textcoords="offset points", xytext=(0, 5),
                    ha="center", fontsize=8.0, fontweight="bold",
                    color=("#065f46" if k == 7 else ("#991b1b" if k == 8 else "#1e3a8a")))

    ax.set_title("(d) Fermeture exacte du degré mixte $(1, 2)$ sur `mod 6` : $N_{\\mathrm{gen}}^{(6)}((1, 2), \\mathrm{OR}) = 7$")
    ax.set_xlabel("Nombre de variables $k$ (polynômes $P_3$ sans zéro non nul de poids pair)")
    ax.set_ylabel("Classes canoniques (mod $S_k \\times \\{\\pm 1\\}$)")
    ax.set_xticks(k_vals)
    ax.set_ylim(0, 325)
    ax.grid(True, axis="y")

    plt.tight_layout(rect=[0, 0, 1, 0.94])
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(fig_path, dpi=180)
    plt.close(fig)
    shutil.copy2(fig_path, artifact_fig_path)
    print(f"Figure sauvegardee : {fig_path} et copiee vers {artifact_fig_path}")


def main():
    t0 = time.time()
    print("=" * 90)
    print("AXE II — THEOREME D4 : RESSERREMENT DES CONSTANTES (8p - 7, 6p - 4, C_d^{D4} -25%)")
    print("=" * 90)

    primes = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 53, 101]
    print(f"{'p':>3s} | {'2p':>4s} | {'4p-2':>5s} | {'D3 (9p-9)':>9s} | {'D4 (8p-7)':>9s} | {'D4 TypeA':>8s} | {'D4 THR2':>8s} | {'Ratio D4':>9s} | {'Ratio TA':>9s}")
    print("-" * 90)

    table_d2 = []
    for p in primes:
        lb = 4 * p - 2
        d3 = compute_d3_bound_2p(p, 2)
        d4 = 16 if p == 3 else compute_d4_bound_2p(p, 2)
        d4_ta = 6 * p - 4
        d4_thr2 = 8 * p - 6 if p >= 5 else 17
        d4_thr2_non_complete = 6 * p - 3
        table_d2.append({
            "p": p,
            "modulus_2p": 2 * p,
            "lower_4p_minus_2": lb,
            "d3_upper": d3,
            "d4_upper_or": d4,
            "d4_upper_typeA": d4_ta,
            "d4_upper_thr2_general": d4_thr2,
            "d4_upper_thr2_non_complete_Gp": d4_thr2_non_complete,
            "ratio_d4_over_lb": round(d4 / lb, 4),
            "ratio_typeA_over_lb": round(d4_ta / lb, 4),
        })
        print(f"{p:3d} | {2*p:4d} | {lb:5d} | {d3:9d} | {d4:9d} | {d4_ta:8d} | {d4_thr2:8d} | "
              f"{d4/lb:8.3f}x | {d4_ta/lb:8.3f}x")

    print("\nConstantes lineaires C_d en degre d = 2..7 :")
    higher_deg = []
    for d in range(2, 8):
        c3 = (d + 1) * ((1 << d) - 1)
        c4 = (3 * d + 5) * (1 << (d - 2)) - (d + 1)
        assert c4 == sum(d0_step(d, j) for j in range(d - 1)) + 2 * d0_step(d, d - 1)
        gain = round(100.0 * (c3 - c4) / c3, 2)
        higher_deg.append({"d": d, "C_d_D3": c3, "C_d_D4": c4, "saving_abs": c3 - c4, "saving_pct": gain})
        print(f"   d = {d} : C_d^{'{D3}'} = {c3:4d}  -->  C_d^{'{D4}'} = {c4:4d}  (gain -{c3 - c4:3d}, -{gain:.1f}%)")

    mod70_d3 = compute_pq_from_R(7, 37)
    mod70_d4 = compute_pq_from_R(7, 34)
    print(f"\nImpact sur mod 70 (2x5x7) : R=37 -> {mod70_d3}  ==>  R=34 (via D4 mod 10 <= 33) -> {mod70_d4} (-{100*(mod70_d3-mod70_d4)/mod70_d3:.1f}%)")

    summary = {
        "theorem_d4": {
            "d2_general_2p": "4p - 2 <= N_gen^{(2p)}(2, OR) <= 8p - 7 (ratio < 2.00x for ALL odd primes p >= 3; <= 16 for p = 3)",
            "d2_thr_t_2p": "N_gen^{(2p)}(2, THR_t) <= 8p - 8 + t (and <= 6p - 3 for THR_{n,2} whenever G_p != K_n)",
            "d2_typeA_ci_zero": "N_gen^{(2p)}(2, OR with exists c_i = 0) <= 6p - 4 (ratio < 1.50x)",
            "higher_d_2p": "N_gen^{(2p)}(d, OR) <= [(3d + 5)*2^{d-2} - (d + 1)]*(p - 1) + 1 (-25% asymptotic vs D3)",
            "mixed_1_2_exact_mod6": {
                "N_gen_exact": 7,
                "formula_3p_minus_2": 7,
                "canonical_counts_k2_to_k8": [6, 28, 119, 272, 109, 13, 0],
                "extinction_tests_k8": 4992,
                "general_bounds_mixed_1_2": "3p - 2 <= N_gen^{(2p)}((1, 2), OR) <= 4p - 3 (with N_sym = 2p - 1)",
            },
            "mod70_improvement": {"d3": mod70_d3, "d4": mod70_d4},
        },
        "table_d2": table_d2,
        "higher_degrees": higher_deg,
        "elapsed_s": round(time.time() - t0, 3),
    }

    out_json = Path("results/mod6_polynomials/theorem_d4_tight_constants.json")
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    fig_path = Path("figures/theorem_d4_tight_constants.png")
    artifact_fig_path = Path(r"C:\Users\Axoul\.gemini\antigravity\brain\52cfbf17-f603-488a-9fbc-24b50a9ed2c4\theorem_d4_tight_constants.png")
    make_figure(summary, fig_path, artifact_fig_path)
    print(f"Synthese JSON : {out_json} ({summary['elapsed_s']}s)")


if __name__ == "__main__":
    main()
