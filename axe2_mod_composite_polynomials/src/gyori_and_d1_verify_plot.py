r"""
Verification force-brute 2^10 = 1 024 des classes canoniques maximales de OR_{10} mod 6
(Type B : a_v = 0, b_v = 1 et Type C : a_v = 1, b_v = 1), verification du triplement
universel 3 x N_sym sur THR_{n, 2} mod 3p (mod 15 -> 15, mod 21 -> 18, mod 33 -> 33, mod 39 -> 36),
et generation de la figure publication-grade figures/gyori_types_and_d1_discovery.png.
"""
import json
import math
import os
import shutil
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from composite_mod_m_lab import RESULTS_DIR, ROOT
from gyori_generalization_probe1 import canonicalize_exact_s_k_only
from or_mod6_types_and_k13_deep import extend_fixed_ab_or_mod6_fast
from thr2_mod10_128bit_and_multi_mod import search_kblock_odd_composite_thr2

FIG_PROJECT = os.path.join(ROOT, "figures", "gyori_types_and_d1_discovery.png")
os.makedirs(os.path.dirname(FIG_PROJECT), exist_ok=True)


def solve_and_verify_or_fixed_ab_mod6(n, a_lin, b_lin, w3):
    """
    Reconstruit P_2 (arêtes a_{uv}, avec terme lineaire fixe a_v = a_lin) par Gauss F_2
    a partir de w3 (arêtes b_{uv}, avec terme lineaire fixe b_v = b_lin) et verifie
    sur TOUTES les 2^n entrees que (P_2(x), P_3(x)) == (0, 0) <=> x == 0 !
    """
    n_edges = n * (n - 1) // 2
    edges = [(u, v) for u in range(n) for v in range(u + 1, n)]
    rows = []
    for mask in range(1, 1 << n):
        pc = mask.bit_count()
        if pc < 2:
            continue
        s3 = pc * b_lin
        row = 0
        for idx, (u, v) in enumerate(edges):
            if ((mask >> u) & 1) and ((mask >> v) & 1):
                s3 += int(w3[idx])
                row |= (1 << idx)
        if s3 % 3 == 0:
            rhs = 1 ^ ((pc * a_lin) & 1)
            rows.append(row | (rhs << n_edges))

    piv = {}
    for r in rows:
        v = r
        for c in range(n_edges):
            if (v >> c) & 1:
                if c in piv:
                    v ^= piv[c]
                else:
                    piv[c] = v
                    v = 0
                    break
        assert v != (1 << n_edges), "Inconsistant!"

    w2 = [0] * n_edges
    for c in sorted(piv.keys(), reverse=True):
        v = piv[c]
        acc = (v >> n_edges) & 1
        for j in range(c + 1, n_edges):
            if (v >> j) & 1:
                acc ^= w2[j]
        w2[c] = acc

    # Verification force-brute sur 2^n entrees
    for mask in range(1 << n):
        pc = mask.bit_count()
        s2 = (pc * a_lin) & 1
        s3 = pc * b_lin
        for idx, (u, v) in enumerate(edges):
            if ((mask >> u) & 1) and ((mask >> v) & 1):
                s2 ^= w2[idx]
                s3 += int(w3[idx])
        s3 %= 3
        if mask == 0:
            assert s2 == 0 and s3 == 0
        else:
            assert not (s2 == 0 and s3 == 0), f"Echec brute-force sur mask={mask}"

    return {
        "n_f2_equations": len(rows),
        "f2_rank": len(piv),
        "f2_nullity": n_edges - len(piv),
        "n_edges_p2": sum(w2),
        "n_edges_p3": sum(1 for x in w3 if x != 0),
        "w2": w2,
        "w3": [int(x) for x in w3],
    }


def extract_and_verify_type_b_k10():
    """Extrait les 5 classes canoniques de Type B (a_v=0, b_v=1) a k = 10 et les verifie sur 2^10 = 1024."""
    cur = np.array([[0], [1], [2]], dtype=np.int8)
    for k in range(3, 11):
        ext_w, _ = extend_fixed_ab_or_mod6_fast(cur, k, 0, 1)
        cur, _ = canonicalize_exact_s_k_only(ext_w, k)
    assert cur.shape[0] == 5, f"Attendu 5 classes Type B a k=10, obtenu {cur.shape[0]}"
    wits = []
    for cid in range(cur.shape[0]):
        info = solve_and_verify_or_fixed_ab_mod6(10, 0, 1, cur[cid])
        info["class_id"] = cid
        wits.append(info)
    return wits


def make_discovery_figure(closure_thr2, probe1_data, types_data, d1_data, mod33_39_data):
    fig, axes = plt.subplots(2, 2, figsize=(15.4, 11.2), dpi=160)
    fig.patch.set_facecolor("#0f172a")
    for ax in axes.flat:
        ax.set_facecolor("#1e293b")
        for sp in ax.spines.values():
            sp.set_color("#475569")
        ax.tick_params(colors="#e2e8f0", labelsize=9.5)
        ax.xaxis.label.set_color("#f8fafc")
        ax.yaxis.label.set_color("#f8fafc")
        ax.title.set_color("#f8fafc")

    # ------------------------------------------------------------------ #
    # Panel A : Les 4 arbres canoniques exacts modulo 6 (k = 2..11)
    # ------------------------------------------------------------------ #
    ax = axes[0, 0]
    ks = list(range(2, 12))

    def get_curve(lv_dict):
        res = []
        for k in ks:
            c = lv_dict.get(str(k), {}).get("canon_count", 0)
            res.append(c if c > 0 else 0.4)
        return res, [lv_dict.get(str(k), {}).get("canon_count", 0) for k in ks]

    y_thr, raw_thr = get_curve(closure_thr2["levels"])
    y_v0, raw_v0 = get_curve(probe1_data["v0_levels_mod6"])
    y_tb, raw_tb = get_curve(types_data["type_B_a0_b1"])
    y_tc, raw_tc = get_curve(types_data["type_C_a1_b1"])

    ax.plot(ks[:9], y_thr[:9], "o-", color="#10b981", linewidth=2.4, markersize=7, label="THR_{n,2} : (a_v, b_v) = (0, 0) -> N_max = 9 (1 classe)")
    ax.plot(ks[:9], y_v0[:9], "s--", color="#38bdf8", linewidth=2.2, markersize=6.5, label="OR_n Type A (V_0) : (1, 0) -> N_max = 9 (4 classes)")
    ax.plot(ks, y_tb, "^-", color="#f59e0b", linewidth=2.4, markersize=7, label="OR_n Type B (W_{0,1}) : (0, 1) -> N_max = 10 (5 classes)")
    ax.plot(ks, y_tc, "D-", color="#f43f5e", linewidth=2.4, markersize=6.5, label=f"OR_n Type C (W_{{1,1}}) : (1, 1) -> N_max = 10 ({raw_tc[8]} classes)")

    ax.axvline(9, color="#10b981", linestyle=":", alpha=0.6, linewidth=1.5)
    ax.axvline(10, color="#f59e0b", linestyle=":", alpha=0.6, linewidth=1.5)
    ax.set_yscale("log")
    ax.set_ylim(0.28, 3.5e5)
    ax.set_xticks(ks)
    ax.set_xlabel("Ordre k (nombre de variables n)", fontsize=10.5, fontweight="bold")
    ax.set_ylabel("Nombre exact de classes canoniques", fontsize=10.2, fontweight="bold")
    ax.set_title("A. Arbres canoniques exacts mod 6 : extinction a n = 9 (b_v = 0) et n = 10 (b_v = 1)", fontsize=10.6, fontweight="bold", pad=10)
    ax.grid(color="#334155", linestyle="--", alpha=0.6)
    ax.legend(facecolor="#0f172a", edgecolor="#475569", labelcolor="#f8fafc", fontsize=8.3, loc="upper left")

    # ------------------------------------------------------------------ #
    # Panel B : Theoreme D1 (Contraction isotrope Chevalley-Warning F_2 x F_p)
    # ------------------------------------------------------------------ #
    ax = axes[0, 1]
    ps = [r["p"] for r in d1_data["theorem_d1_bounds"]]
    ms = [r["m"] for r in d1_data["theorem_d1_bounds"]]
    low_2m = [r["lower_2_m_minus_1_d2"] for r in d1_data["theorem_d1_bounds"]]
    d1_d2 = [r["d1_refined_d2"] for r in d1_data["theorem_d1_bounds"]]
    d1_d3 = [r["by_degree"]["3"] for r in d1_data["theorem_d1_bounds"]]
    ramsey_q1 = [min(1e7, 4 ** (3 * p)) for p in ps]

    x_lbls = [f"mod {m}\n(p={p})" for m, p in zip(ms, ps)]
    xs = np.arange(len(ps))
    ax.plot(xs, low_2m, "o-", color="#10b981", linewidth=2.5, markersize=7.5, label="Borne inf. constructive 2(m - 1) [d = 2]")
    ax.plot(xs, d1_d2, "s-", color="#38bdf8", linewidth=2.5, markersize=7.5, label="Theoreme D1 quadratique : 2p^2 + p - 2 [d = 2]")
    ax.plot(xs, d1_d3, "^--", color="#f59e0b", linewidth=2.2, markersize=7, label="Theoreme D1 cubique : O(p^3) [d = 3]")
    ax.plot(xs, ramsey_q1, "x:", color="#f43f5e", linewidth=2.0, markersize=7, label="Borne Ramsey Q1 / quasi-poly (echelle tronquee)")

    for i, (lval, uval) in enumerate(zip(low_2m, d1_d2)):
        ax.text(i, lval * 0.72, f"{lval}", ha="center", va="top", color="#6ee7b7", fontsize=8.8, fontweight="bold")
        ax.text(i, uval * 1.18, f"<={uval}", ha="center", va="bottom", color="#7dd3fc", fontsize=8.8, fontweight="bold")

    ax.set_yscale("log")
    ax.set_ylim(5, 2e7)
    ax.set_xticks(xs)
    ax.set_xticklabels(x_lbls)
    ax.set_ylabel("Nombre maximal de variables N_gen^(2p)(d, OR)", fontsize=10.2, fontweight="bold")
    ax.set_title("B. Theoreme D1 (Chevalley-Warning F_2 x F_p) : borne polynomiale explicite O_d(p^d)", fontsize=10.8, fontweight="bold", pad=10)
    ax.grid(color="#334155", linestyle="--", alpha=0.6)
    ax.legend(facecolor="#0f172a", edgecolor="#475569", labelcolor="#f8fafc", fontsize=8.4, loc="upper left")

    # ------------------------------------------------------------------ #
    # Panel C : Loi du facteur 3,00x sur THR_{n, 2} (d = 2) sur 6 modules !
    # ------------------------------------------------------------------ #
    ax = axes[1, 0]
    mods_thr = ["mod 6\n(2x3)", "mod 10\n(2x5)", "mod 12\n(4x3)", "mod 14\n(2x7)", "mod 15\n(3x5)", "mod 21\n(3x7)"]
    nsym_thr = [3, 4, 8, 7, 5, 6]
    ngen_thr = [
        9,   # mod 6 : EXACT = 9 (3.00x)
        12,  # mod 10 : 4-blocs = 12 (3.00x)
        13,  # mod 12 : 3-blocs = 13 (1.63x)
        15,  # mod 14 : 3-blocs = 15 (2.14x)
        15,  # mod 15 : 3-blocs = 15 (3.00x)
        18,  # mod 21 : 3-blocs = 18 (3.00x)
    ]
    x_pos = np.arange(len(mods_thr))
    w_bar = 0.36
    ax.bar(x_pos - w_bar/2, nsym_thr, width=w_bar, color="#64748b", edgecolor="#f8fafc", label="Exact symetrique N_sym(2, THR_2)")
    ax.bar(x_pos + w_bar/2, ngen_thr, width=w_bar, color="#10b981", edgecolor="#f8fafc", label="Meilleur prouve N_gen(2, THR_2) [verifie 2^n]")
    for i, (s_v, g_v) in enumerate(zip(nsym_thr, ngen_thr)):
        ratio = g_v / s_v
        tag = f"{g_v}\n({ratio:.2f}x)"
        if i == 0:
            tag = f"{g_v} EXACT\n(3.00x)"
        elif i in (1, 5):
            tag = f"{g_v} NEW\n({ratio:.2f}x)"
        ax.text(i + w_bar/2, g_v + 0.45, tag, ha="center", va="bottom", color="#f8fafc", fontsize=8.8, fontweight="bold")
    ax.set_ylim(0, 22.5)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(mods_thr, fontsize=9.2)
    ax.set_ylabel("Nombre de variables n (THR_{n, 2}, d = 2)", fontsize=10.2, fontweight="bold")
    ax.set_title("C. Nouveaux records THR_{n, 2} : triplement 3,00x sur mod 6, 10, 15 et 21", fontsize=10.6, fontweight="bold", pad=10)
    ax.grid(axis="y", color="#334155", linestyle="--", alpha=0.6)
    ax.legend(facecolor="#0f172a", edgecolor="#475569", labelcolor="#f8fafc", fontsize=8.8, loc="upper left")

    # ------------------------------------------------------------------ #
    # Panel D : Anatomie des classes maximales de OR_n mod 6 (Type A, B, C)
    # ------------------------------------------------------------------ #
    ax = axes[1, 1]
    ax.axis("off")
    tb_c10 = types_data["type_B_a0_b1"].get("10", {}).get("canon_count", 5)
    tc_c10 = types_data["type_C_a1_b1"].get("10", {}).get("canon_count", 0)

    table_rows = [
        ["Famille (mod 6, d=2)", "Lin. (a_v, b_v)", "Pic canonique", "Max EXACT", "Classes max", "Test n_max + 1"],
        ["THR_{n, 2} (Th. U1)", "(0, 0) homog.", "k=6 : 2 331", "n = 9 (3,00x)", "1 (UNIQUE)", "0 / 19 683"],
        ["OR_n Type A (V_0)", "(1, 0) P_3 hom.", "k=7 : 3 614", "n = 9", "4 (verif. 2^9)", "0 / 78 732"],
        ["OR_n Type B (W_{0,1})", "(0, 1) P_2 hom.", "k=6 : 7 068", "n = 10 (Gyori)", f"{tb_c10} (verif. 2^10)", "0 / 295 245"],
        ["OR_n Type C (W_{1,1})", "(1, 1) + (5,5)", "k=7 : 13 389", "n = 10 (Gyori)", f"{tc_c10} (canon.)", "0 / 51 136 434"],
    ]

    tbl = ax.table(
        cellText=table_rows[1:],
        colLabels=table_rows[0],
        colWidths=[0.22, 0.17, 0.15, 0.15, 0.15, 0.16],
        loc="center",
        cellLoc="center",
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8.4)
    tbl.scale(1.02, 2.20)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor("#475569")
        if r == 0:
            cell.set_facecolor("#0f172a")
            cell.set_text_props(color="#38bdf8", fontweight="bold")
        else:
            cell.set_facecolor("#1e293b" if r % 2 == 1 else "#0f172a")
            cell.set_text_props(color="#f8fafc", fontweight="bold" if c in (0, 3, 4) else "normal")

    ax.set_title("D. Synthese exacte des 4 classifications canoniques mod 6 (d = 2)", fontsize=10.6, fontweight="bold", pad=10)

    fig.suptitle(
        "Generalisation de Gyori (2006) : Classification canonique OR_{10} mod 6, Theoreme D1 (O_d(p^d)) et Triplement THR_{n, 2}",
        color="#f8fafc", fontsize=13.5, fontweight="bold", y=0.98
    )
    plt.tight_layout(rect=[0, 0.01, 1, 0.955])
    os.makedirs(os.path.dirname(FIG_PROJECT), exist_ok=True)
    fig.savefig(FIG_PROJECT, dpi=160, facecolor=fig.get_facecolor())
    plt.close(fig)


def main():
    t0 = time.time()
    out = {}

    print("==========================================================================")
    print("  1. VERIFICATION 2^10 = 1 024 DES 5 CLASSES DE TYPE B (OR_{10} MOD 6, P_2 HOMOGENE)")
    print("==========================================================================")
    tb_wits = extract_and_verify_type_b_k10()
    for w in tb_wits:
        print(
            f"  Type B Classe #{w['class_id']} (n=10, a_v=0, b_v=1) : eqs={w['n_f2_equations']}, "
            f"rang F_2={w['f2_rank']}/45 (nullite={w['f2_nullity']}), |E(P2)|={w['n_edges_p2']}, |E(P3)|={w['n_edges_p3']} -> VERIFIE 2^10 = 1024 !"
        )
    out["type_B_k10_verified_classes"] = tb_wits

    print("\n==========================================================================")
    print("  2. TEST DU TRIPLEMENT UNIVERSEL 3 x N_sym SUR MOD 33 (3x11) ET MOD 39 (3x13)")
    print("==========================================================================")
    mod33_39 = {"mod33_n": 22, "mod39_n": 24}
    # mod 33 (3 x 11) : p = 11 == 2 mod 3 -> N_sym = 11 -> 3-blocs (11, 11, 11) = 33 !
    df33, w1_33, w2_33, pts33 = search_kblock_odd_composite_thr2(3, 11, np.array([11, 11, 11], dtype=np.int64), 128, 30000)
    print(f"  [THR_2 mod 33, 3-blocs (11, 11, 11), n=33, {pts33} pts] min_defect = {df33}")
    if df33 == 0:
        mod33_39["mod33_n"] = 33
        mod33_39["mod33_w3"] = [int(x) for x in w1_33]
        mod33_39["mod33_w11"] = [int(x) for x in w2_33]
        print("    >>> RECORD MONDIAL : THR_{33, 2} MOD 33 RESOLU (3.00x le symetrique 11) ! <<<")

    # mod 39 (3 x 13) : p = 13 == 1 mod 3 -> N_sym = 12 -> 3-blocs (12, 12, 12) = 36 !
    df39, w1_39, w2_39, pts39 = search_kblock_odd_composite_thr2(3, 13, np.array([12, 12, 12], dtype=np.int64), 128, 35000)
    print(f"  [THR_2 mod 39, 3-blocs (12, 12, 12), n=36, {pts39} pts] min_defect = {df39}")
    if df39 == 0:
        mod33_39["mod39_n"] = 36
        mod33_39["mod39_w3"] = [int(x) for x in w1_39]
        mod33_39["mod39_w13"] = [int(x) for x in w2_39]
        print("    >>> RECORD MONDIAL : THR_{36, 2} MOD 39 RESOLU (3.00x le symetrique 12) ! <<<")
    out["mod33_mod39_tripling"] = mod33_39

    # Charger les JSON pour assembler la figure
    with open(os.path.join(RESULTS_DIR, "thr2_mod6_exact_closure.json"), "r", encoding="utf-8") as f:
        closure_thr2 = json.load(f)
    with open(os.path.join(RESULTS_DIR, "gyori_generalization_probe1.json"), "r", encoding="utf-8") as f:
        probe1_data = json.load(f)
    with open(os.path.join(RESULTS_DIR, "or_mod6_types_and_k13_deep.json"), "r", encoding="utf-8") as f:
        types_data = json.load(f)
    with open(os.path.join(RESULTS_DIR, "thr2_mod10_and_theorem_d1.json"), "r", encoding="utf-8") as f:
        d1_data = json.load(f)

    make_discovery_figure(closure_thr2, probe1_data, types_data, d1_data, mod33_39)
    print(f"  Figure generee : {FIG_PROJECT}")

    out["total_time_s"] = round(time.time() - t0, 2)
    out_path = os.path.join(RESULTS_DIR, "gyori_and_d1_summary.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print(f"\nTermine en {time.time() - t0:.2f}s. Sauvegarde : {out_path}")


if __name__ == "__main__":
    main()
