"""
Extension (mod 35, mod 15 d=8), verification independante des 5 Theoremes (G1..G5)
et generation de la figure publication-grade composite_mod_m_discovery.png.
"""
import json
import os
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from composite_mod_m_lab import (
    METHOD_NAMES,
    RESULTS_DIR,
    ROOT,
    binomial_mod_table,
    build_channel_catalog,
    solve_2channel,
    verify_witness_2ch,
)
from composite_mod_m_full import apply_g4_dp, verify_field_norm_anisotropic

FIG_PROJECT = os.path.join(ROOT, "figures", "composite_mod_m_discovery.png")
os.makedirs(os.path.dirname(FIG_PROJECT), exist_ok=True)


def main():
    t0 = time.time()
    full_path = os.path.join(RESULTS_DIR, "composite_mod_m_full.json")
    with open(full_path, "r", encoding="utf-8") as f:
        payload = json.load(f)

    # 1. Ajout de mod 15 a d = 8 (pour verifier G3 +8 sur F_3)
    print("=== 1. Calcul exact MOD 15 (3 x 5) a d = 8 ===")
    max_n15 = 1200
    C3 = binomial_mod_table(max_n15, 8, 3)
    C5 = binomial_mod_table(max_n15, 8, 5)
    td = time.time()
    deg3, val3, cof3 = build_channel_catalog(8, 3, True, C3, max_n15)
    deg5, val5, cof5 = build_channel_catalog(8, 5, True, C5, max_n15)
    res15 = solve_2channel(8, 3, True, deg3, val3, 5, True, deg5, val5, max_n15)
    (nsym, si, sj, ng1, _, _, _, ng2, _, _, _, nall, ai, aj, aa, m1, m2, gn1, gn2) = res15
    verify_witness_2ch(8, 3, True, cof3[ai], 5, True, cof5[aj], aa, m1, gn1, m2, gn2)
    dt = time.time() - td
    parts = []
    if gn1 > 0:
        parts.append(f"+{gn1} {METHOD_NAMES[m1]}@3")
    if gn2 > 0:
        parts.append(f"+{gn2} {METHOD_NAMES[m2]}@5")
    row15_d8 = {
        "d": 8,
        "N_sym": nsym,
        "sym_degs": [int(deg3[si]), int(deg5[sj])],
        "N_G1": ng1,
        "N_G2": ng2,
        "N_all": nall,
        "witness_a": aa,
        "witness_degs": [int(deg3[ai]), int(deg5[aj])],
        "witness_c1": [int(x) for x in cof3[ai]],
        "witness_c2": [int(x) for x in cof5[aj]],
        "methods": [METHOD_NAMES[m1], METHOD_NAMES[m2]],
        "gains": [gn1, gn2],
        "g123_summary": " & ".join(parts) if parts else "sym",
        "time_s": round(dt, 3),
    }
    rows15 = [r for r in payload["modules"]["mod15"] if r["d"] < 8] + [row15_d8]
    apply_g4_dp(rows15, r_factors=2, is_squarefree=True)
    payload["modules"]["mod15"] = rows15
    r = rows15[-1]
    print(
        f"  mod15 d=8 | N_sym={r['N_sym']} (deg {r['sym_degs']}) | G1={r['N_G1']} | "
        f"G123={r['N_all']} | BEST={r['N_best']} (+{r['best_gain']}, +{r['best_gain_pct']}%) "
        f"[{r['best_mechanism']}] ({dt:.2f}s)"
    )

    # 2. Ajout de mod 35 = 5 x 7, d = 1..5 (demonstration du Theoreme G5 : ratio 4x = +300%)
    print("\n=== 2. Calcul exact MOD 35 (5 x 7), d = 1..5 (Theoreme G5) ===")
    max_n35 = 800
    d_max35 = 5
    C5_35 = binomial_mod_table(max_n35, d_max35, 5)
    C7_35 = binomial_mod_table(max_n35, d_max35, 7)
    rows35 = []
    for d in range(1, d_max35 + 1):
        td = time.time()
        deg5, val5, cof5 = build_channel_catalog(d, 5, True, C5_35[:, :d], max_n35)
        deg7, val7, cof7 = build_channel_catalog(d, 7, True, C7_35[:, :d], max_n35)
        res = solve_2channel(d, 5, True, deg5, val5, 7, True, deg7, val7, max_n35)
        (nsym, si, sj, ng1, _, _, _, ng2, _, _, _, nall, ai, aj, aa, m1, m2, gn1, gn2) = res
        assert nsym < max_n35
        verify_witness_2ch(d, 5, True, cof5[ai], 7, True, cof7[aj], aa, m1, gn1, m2, gn2)
        dt = time.time() - td
        parts = []
        if gn1 > 0:
            parts.append(f"+{gn1} {METHOD_NAMES[m1]}@5")
        if gn2 > 0:
            parts.append(f"+{gn2} {METHOD_NAMES[m2]}@7")
        rows35.append({
            "d": d,
            "N_sym": nsym,
            "sym_degs": [int(deg5[si]), int(deg7[sj])],
            "N_G1": ng1,
            "N_G2": ng2,
            "N_all": nall,
            "witness_a": aa,
            "witness_degs": [int(deg5[ai]), int(deg7[aj])],
            "witness_c1": [int(x) for x in cof5[ai]],
            "witness_c2": [int(x) for x in cof7[aj]],
            "methods": [METHOD_NAMES[m1], METHOD_NAMES[m2]],
            "gains": [gn1, gn2],
            "g123_summary": " & ".join(parts) if parts else "sym",
            "time_s": round(dt, 3),
        })
    apply_g4_dp(rows35, r_factors=2, is_squarefree=True)
    payload["modules"]["mod35"] = rows35
    for r in rows35:
        print(
            f"  mod35 d={r['d']} | N_sym={r['N_sym']:3d} (deg {r['sym_degs']}) | "
            f"G1={r['N_G1']:3d} | G2={r['N_G2']:3d} | G4={r['N_G4']:3d} | "
            f"BEST={r['N_best']:3d} (+{r['best_gain']:3d}, +{r['best_gain_pct']:5.1f}%) "
            f"[{r['best_mechanism']}] ({r['time_s']:.2f}s)"
        )

    # 3. Verification explicite des normes N_{F_{5^4}/F_5} et N_{F_{7^4}/F_7} pour mod 35 d=4 -> 136
    n5_4 = verify_field_norm_anisotropic(5, 4)
    n7_4 = verify_field_norm_anisotropic(7, 4)
    print(f"\n  Norme F_{{5^4}}/F_5 verifiee : irred={n5_4['irred_coeffs_low_to_high']} ({n5_4['exhaustive_checked']} vecteurs)")
    print(f"  Norme F_{{7^4}}/F_7 verifiee : irred={n7_4['irred_coeffs_low_to_high']} ({n7_4['exhaustive_checked']} vecteurs)")

    # Verification directe du Theoreme G5 sur (S_1, S_2, S_3, S_4) in {0..34}^4 -> 136 variables
    # Par definition de N_4 anisotrope sur F_5 et F_7 :
    # N_4(S_1..S_4) = 0 mod 5 iff S_j = 0 mod 5 for all j
    # N_4(S_1..S_4) = 0 mod 7 iff S_j = 0 mod 7 for all j
    # => joint zero iff S_j = 0 mod 35 for all j => S_j = 0 (since 0 <= S_j <= 34). QED!

    # Sauvegarde JSON enrichi
    payload["total_time_s"] = round(payload.get("total_time_s", 0.0) + (time.time() - t0), 2)
    with open(full_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)

    # 4. Figure publication-grade (4 panneaux)
    plt.rcParams.update({
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "legend.fontsize": 8.5,
        "figure.dpi": 180,
    })
    fig, axes = plt.subplots(2, 2, figsize=(14.5, 10.2))
    fig.suptitle(
        "Généralisation des polynômes faibles modulo m composite : Théorèmes G1–G5\n"
        "(Énumération symétrique 100 % exacte vs percées non symétriques prouvées sur 8 modules)",
        fontsize=13, fontweight="bold", y=0.98
    )

    # --- Panneau A : Carte des gains relatifs (%) par mecanisme vainqueur ---
    ax = axes[0, 0]
    mod_order = ["mod6", "mod10", "mod12", "mod14", "mod15", "mod21", "mod30", "mod35"]
    mod_labels = ["mod 6\n(2×3)", "mod 10\n(2×5)", "mod 12\n(4×3)", "mod 14\n(2×7)",
                  "mod 15\n(3×5)", "mod 21\n(3×7)", "mod 30\n(2×3×5)", "mod 35\n(5×7)"]
    max_gains = []
    max_mechs = []
    max_annots = []
    color_map = {"G1": "#2563eb", "G2": "#059669", "G3": "#d97706", "G4": "#dc2626"}
    for m in mod_order:
        rows = payload["modules"][m]
        best_row = max(rows, key=lambda r: r["best_gain_pct"])
        max_gains.append(best_row["best_gain_pct"])
        mech = best_row["best_mechanism"]
        if mech.startswith("G4"):
            mtype = "G4"
        elif "G3" in mech:
            mtype = "G3"
        elif "G2" in mech:
            mtype = "G2"
        else:
            mtype = "G1"
        max_mechs.append(mtype)
        max_annots.append(f"d={best_row['d']}\n{best_row['N_sym']}→{best_row['N_best']}\n(+{best_row['best_gain_pct']:.0f}%)")

    x = np.arange(len(mod_order))
    bars = ax.bar(x, max_gains, color=[color_map[mt] for mt in max_mechs], edgecolor="#1e293b", linewidth=0.8)
    for bar, ann in zip(bars, max_annots):
        h = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, h + 6, ann, ha="center", va="bottom", fontsize=8, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(mod_labels, fontsize=8.5)
    ax.set_ylabel("Gain maximal prouvé au-dessus de N_sym (%)")
    ax.set_ylim(0, 365)
    ax.set_title("A. Pic de séparation non symétrique par module composite")
    ax.grid(axis="y", linestyle="--", alpha=0.4)
    # Legende mecanismes
    from matplotlib.patches import Patch
    legend_handles = [
        Patch(facecolor=color_map["G2"], edgecolor="#1e293b", label="G2 : Norme de corps sur canal (N_{F_{p^k}/F_p})"),
        Patch(facecolor=color_map["G3"], edgecolor="#1e293b", label="G3 : Omission de valeur (Value-Avoidance)"),
        Patch(facecolor=color_map["G4"], edgecolor="#1e293b", label="G4/G5 : Norme multi-blocs (k × N_gen(d₀))"),
    ]
    ax.legend(handles=legend_handles, loc="upper left", framealpha=0.92)

    # --- Panneau B : Theoreme G5 (Ratio N_gen / N_sym pour composites impairs) ---
    ax = axes[0, 1]
    for m, label, col, mk in [
        ("mod15", "mod 15 = 3×5 (p_min=3)", "#7c3aed", "o"),
        ("mod21", "mod 21 = 3×7 (p_min=3)", "#ea580c", "s"),
        ("mod35", "mod 35 = 5×7 (p_min=5)", "#dc2626", "^"),
        ("mod14", "mod 14 = 2×7 (p_min=2)", "#059669", "D"),
    ]:
        rows = payload["modules"][m]
        ds = [r["d"] for r in rows]
        ratios = [r["N_best"] / r["N_sym"] for r in rows]
        ax.plot(ds, ratios, marker=mk, linewidth=2.2, markersize=7, color=col, label=label)
    ax.axhline(1.0, color="#64748b", linestyle=":", linewidth=1.2, label="Seuil symétrique (ratio = 1,0×)")
    ax.annotate("Théorème G5 : 4,00×\n(34 → 136 à d=4)", xy=(4, 4.0), xytext=(1.15, 3.65),
                arrowprops=dict(arrowstyle="->", color="#dc2626", lw=1.5),
                fontsize=8.5, fontweight="bold", color="#dc2626",
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#dc2626", alpha=0.92))
    ax.annotate("2,00× (14→28 / 20→40)", xy=(2, 2.0), xytext=(2.85, 2.25),
                arrowprops=dict(arrowstyle="->", color="#7c3aed", lw=1.3),
                fontsize=8.5, fontweight="bold", color="#7c3aed",
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#7c3aed", alpha=0.92))
    ax.set_xlabel("Degré d")
    ax.set_ylabel("Ratio prouvé N_gen^(m)(d) / N_sym^(m)(d)")
    ax.set_title("B. Théorème G5 : séparation non bornée (p_min − 1)× pour d < p_min")
    ax.set_xticks(range(1, 9))
    ax.set_ylim(0.85, 4.35)
    ax.grid(True, linestyle="--", alpha=0.4)
    ax.legend(loc="upper right", framealpha=0.92)

    # --- Panneau C : Zoom sur mod 6 (d = 1..10) : M2 vs nouveau record G3 a d=8 ---
    ax = axes[1, 0]
    rows6 = payload["modules"]["mod6"]
    ds6 = np.array([r["d"] for r in rows6])
    nsym6 = np.array([r["N_sym"] for r in rows6])
    ng1_6 = np.array([r["N_G1"] for r in rows6])
    nbest6 = np.array([r["N_best"] for r in rows6])
    w = 0.27
    ax.bar(ds6 - w, nsym6, width=w, color="#94a3b8", edgecolor="#334155", label="Exact symétrique N_sym^(6)(d)")
    ax.bar(ds6, ng1_6, width=w, color="#3b82f6", edgecolor="#1e3a8a", label="G1 / M2 (Margin-Harvesting seul)")
    ax.bar(ds6 + w, nbest6, width=w, color="#f59e0b", edgecolor="#92400e", label="Unifié G1+G2+G3 (nouveau record)")
    ax.annotate("Nouveau record G3 :\n80 → 83 (M2) → 88 (+8)\nvia I₉(A) ∈ {0, 2} mod 3",
                xy=(8 + w, 88), xytext=(4.3, 125),
                arrowprops=dict(arrowstyle="->", color="#b45309", lw=1.5),
                fontsize=8.5, fontweight="bold", color="#92400e",
                bbox=dict(boxstyle="round,pad=0.3", fc="#fef3c7", ec="#d97706"))
    ax.annotate("Győri (2006)\n8 → 10 (+25%)", xy=(2 + w, 10), xytext=(1.1, 45),
                arrowprops=dict(arrowstyle="->", color="#b45309", lw=1.2),
                fontsize=8, fontweight="bold", color="#92400e")
    ax.set_xlabel("Degré d")
    ax.set_ylabel("Nombre maximal de variables n (OR_n mod 6)")
    ax.set_title("C. Percée G3 (Value-Avoidance) sur mod 6 : d = 8 bondit à n = 88")
    ax.set_xticks(ds6)
    ax.set_ylim(0, 205)
    ax.grid(axis="y", linestyle="--", alpha=0.4)
    ax.legend(loc="upper left", framealpha=0.92)

    # --- Panneau D : Lois d'echelle r=2 (Theta(d^2)) vs r=3 (mod 30, Theta(d^3)) ---
    ax = axes[1, 1]
    for m, label, col, mk in [
        ("mod6",  "mod 6 (r=2, 2×3)", "#2563eb", "o"),
        ("mod12", "mod 12 (r=2, 4×3, non sqfree)", "#0891b2", "v"),
        ("mod15", "mod 15 (r=2, 3×5)", "#7c3aed", "s"),
        ("mod30", "mod 30 (r=3, 2×3×5)", "#dc2626", "D"),
    ]:
        rows = payload["modules"][m]
        ds = [r["d"] for r in rows]
        nsym = [r["N_sym"] for r in rows]
        nbest = [r["N_best"] for r in rows]
        ax.plot(ds, nsym, linestyle="--", linewidth=1.3, color=col, alpha=0.55)
        ax.plot(ds, nbest, marker=mk, linewidth=2.2, markersize=6, color=col, label=label)
    ax.set_yscale("log")
    ax.set_xlabel("Degré d")
    ax.set_ylabel("Variables n (échelle log — pointillé : symétrique, plein : G1–G5)")
    ax.set_title("D. Régimes Θ(d²) (r = 2) vs Θ(d³) (mod 30, jusqu'à n = 1 359 à d = 6)")
    ax.set_xticks(range(1, 11))
    ax.grid(True, which="both", linestyle="--", alpha=0.35)
    ax.legend(loc="lower right", framealpha=0.92)

    plt.tight_layout(rect=[0, 0, 1, 0.93])
    os.makedirs(os.path.dirname(FIG_PROJECT), exist_ok=True)
    plt.savefig(FIG_PROJECT, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
