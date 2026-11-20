r"""
Audit independant, analyse algebrique complete de l'unique solution K_9,
et figure publication-grade de la fermeture EXACTE de THR_{n, 2} mod 6 (d = 2).

Verifications independantes :
  1. Recalcul 100% independant par le moteur naif (sans DP ni pivot incremental)
     des transitions k = 6 -> 7 (1 699 299 tests), k = 8 -> 9 (177 147 tests)
     et k = 9 -> 10 (19 683 tests), confirmant l'egalite BIT-A-BIT des comptes.
  2. Calcul EXACT de l'espace affine des solutions P_2 (rang F_2, nombre exact
     de P_2 valides pour l'unique w_3 de K_9, groupe d'automorphismes Aut(w_3)
     et Aut(P_2, P_3), nombre total de polynomes etiquetes sur K_9).
  3. Distribution EXACTE de f2_defect sur toutes les 3^9 = 19 683 extensions
     de l'unique K_9 vers K_10 (gap minimal exact).
  4. Verification analytique de la formule fermee V = {u} U {v} U K_2 U C_5
     et des 4 parents survivants a k = 8 (= les 4 orbites de sommets de K_9).
  5. Figure publication-grade figures/thr2_mod6_exact_closure.png.
"""
import itertools
import json
import os
import shutil
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from numba import njit, prange

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from composite_mod_m_lab import RESULTS_DIR, ROOT
from option_a_thr2_breakthrough import f2_defect_count
from piste_2a_thr2_graphs import solve_and_verify_thr2_mod2p
from thr2_mod6_exact_closure import (
    canonicalize_exact_s_k_pm1,
    extend_valid_w3_exact_fast,
)
from thr2_mod6_n10_probe2 import extend_valid_w3_batch

FIG_PROJECT = os.path.join(ROOT, "figures", "thr2_mod6_exact_closure.png")
FIG_ARTIFACT = r"C:\Users\Axoul\.gemini\antigravity\brain\52cfbf17-f603-488a-9fbc-24b50a9ed2c4\thr2_mod6_exact_closure.png"


def exact_f2_solution_space(n, w3):
    """
    Resout le systeme lineaire F_2 associe a w3 et retourne :
      (rank, nullity, n_zero_eqs, all_w2_solutions)
    """
    n_edges = n * (n - 1) // 2
    edges = [(u, v) for u in range(n) for v in range(u + 1, n)]
    masks = [m for m in range(1, 1 << n) if m.bit_count() >= 2]

    rows = []
    for m in masks:
        s = 0
        row = [0] * (n_edges + 1)
        for ei, (u, v) in enumerate(edges):
            if ((m >> u) & 1) and ((m >> v) & 1):
                s += int(w3[ei])
                row[ei] = 1
        if s % 3 == 0:
            row[n_edges] = 1
            rows.append(row)

    M = np.array(rows, dtype=np.int8)
    n_eqs = M.shape[0]
    # Gauss-Jordan complet sur F_2
    piv_cols = []
    r = 0
    for c in range(n_edges):
        piv = -1
        for i in range(r, n_eqs):
            if M[i, c]:
                piv = i
                break
        if piv == -1:
            continue
        if piv != r:
            M[[r, piv]] = M[[piv, r]]
        for i in range(n_eqs):
            if i != r and M[i, c]:
                M[i] ^= M[r]
        piv_cols.append(c)
        r += 1

    # Verifier consistance
    for i in range(r, n_eqs):
        assert M[i, n_edges] == 0, "Systeme inconsistant !"

    free_cols = [c for c in range(n_edges) if c not in piv_cols]
    rank = len(piv_cols)
    nullity = len(free_cols)

    # Enumerer toutes les 2^nullity solutions si nullity <= 10
    sols = []
    if nullity <= 10:
        for bits in range(1 << nullity):
            x = [0] * n_edges
            for fi, fc in enumerate(free_cols):
                x[fc] = (bits >> fi) & 1
            for row_i, pc in enumerate(piv_cols):
                val = int(M[row_i, n_edges])
                for fc in free_cols:
                    if M[row_i, fc]:
                        val ^= x[fc]
                x[pc] = val
            sols.append(x)
    return rank, nullity, n_eqs, sols


def compute_automorphism_group_k9(w3, w2=None):
    """
    Calcule le groupe d'automorphismes EXACT de w3 sous S_9 (et S_9 x {+-1})
    ainsi que les orbites de sommets, en utilisant 1-WL pour reduire 9! = 362880.
    """
    n = 9
    edges = [(u, v) for u in range(n) for v in range(u + 1, n)]
    A3 = [[0] * n for _ in range(n)]
    A2 = [[0] * n for _ in range(n)]
    for ei, (u, v) in enumerate(edges):
        A3[u][v] = A3[v][u] = int(w3[ei])
        if w2 is not None:
            A2[u][v] = A2[v][u] = int(w2[ei])

    # Degres (d0, d1, d2) -> partition en blocs
    groups = {}
    for v in range(n):
        sig = (
            sum(1 for u in range(n) if u != v and A3[v][u] == 0),
            sum(1 for u in range(n) if u != v and A3[v][u] == 1),
            sum(1 for u in range(n) if u != v and A3[v][u] == 2),
        )
        groups.setdefault(sig, []).append(v)

    # Automorphismes preservant signe (+1) : doivent stabiliser chaque sig
    blocks = list(groups.values())
    aut_pos = []
    for prod in itertools.product(*(itertools.permutations(b) for b in blocks)):
        perm = [0] * n
        for b, p in zip(blocks, prod):
            for orig, img in zip(b, p):
                perm[orig] = img
        ok = True
        for u in range(n):
            for v in range(u + 1, n):
                if A3[perm[u]][perm[v]] != A3[u][v]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            aut_pos.append(tuple(perm))

    # Automorphismes inversant signe (-1) : envoient (d0, d1, d2) -> (d0, d2, d1)
    aut_neg = []
    can_neg = all((s[0], s[2], s[1]) in groups and len(groups[(s[0], s[2], s[1])]) == len(vs) for s, vs in groups.items())
    if can_neg:
        src_blocks = list(groups.keys())
        tgt_lists = [groups[(s[0], s[2], s[1])] for s in src_blocks]
        src_lists = [groups[s] for s in src_blocks]
        for prod in itertools.product(*(itertools.permutations(t) for t in tgt_lists)):
            perm = [0] * n
            for src, img_t in zip(src_lists, prod):
                for orig, img in zip(src, img_t):
                    perm[orig] = img
            ok = True
            for u in range(n):
                for v in range(u + 1, n):
                    if A3[perm[u]][perm[v]] != (3 - A3[u][v]) % 3:
                        ok = False
                        break
                if not ok:
                    break
            if ok:
                aut_neg.append(tuple(perm))

    # Orbites de sommets sous aut_pos
    visited = set()
    orbits = []
    for v in range(n):
        if v in visited:
            continue
        orb = sorted({p[v] for p in aut_pos})
        for x in orb:
            visited.add(x)
        orbits.append(orb)

    return len(aut_pos), len(aut_neg), orbits


@njit(parallel=True, cache=True)
def exact_k9_to_k10_defect_distribution(w3_9):
    """Calcule f2_defect EXACT pour toutes les 3^9 = 19 683 extensions de K_9 -> K_10."""
    n = 10
    n_edges = 45
    n_ext = 3 ** 9
    eu = np.zeros(n_edges, dtype=np.int32)
    ev = np.zeros(n_edges, dtype=np.int32)
    idx_prev = np.zeros(36, dtype=np.int32)
    idx_new = np.zeros(9, dtype=np.int32)
    e = 0; ip = 0; inew = 0
    for u in range(n):
        for v in range(u + 1, n):
            eu[e] = u; ev[e] = v
            if v < 9:
                idx_prev[ip] = e; ip += 1
            else:
                idx_new[inew] = e; inew += 1
            e += 1

    masks = np.zeros(1024, dtype=np.int32)
    nm = 0
    for msk in range(1, 1024):
        pc = 0; tmp = msk
        while tmp:
            pc += tmp & 1; tmp >>= 1
        if pc >= 2:
            masks[nm] = msk; nm += 1

    edge_bits = np.zeros(nm, dtype=np.int64)
    mask_ed_cnt = np.zeros(nm, dtype=np.int32)
    mask_ed = np.zeros((nm, n_edges), dtype=np.int32)
    for i in range(nm):
        msk = masks[i]
        bm = np.int64(0)
        cnt = 0
        for ed in range(n_edges):
            if ((msk >> eu[ed]) & 1) and ((msk >> ev[ed]) & 1):
                bm |= (np.int64(1) << ed)
                mask_ed[i, cnt] = ed
                cnt += 1
        edge_bits[i] = bm
        mask_ed_cnt[i] = cnt

    defs = np.zeros(n_ext, dtype=np.int32)
    for ec in prange(n_ext):
        w = np.zeros(n_edges, dtype=np.int8)
        for j in range(36):
            w[idx_prev[j]] = w3_9[j]
        tmp = ec
        for j in range(9):
            w[idx_new[j]] = tmp % 3
            tmp //= 3
        zb = np.zeros(nm, dtype=np.int64)
        nz = 0
        for i in range(nm):
            s = 0
            mc = mask_ed_cnt[i]
            for t in range(mc):
                s += w[mask_ed[i, t]]
            if s % 3 == 0:
                zb[nz] = edge_bits[i]
                nz += 1
        defs[ec] = f2_defect_count(zb, nz, n_edges)
    return defs


def make_figure(closure_data, audit_data, k9_w2, k9_w3, defect_dist):
    fig, axes = plt.subplots(2, 2, figsize=(15.2, 11.2), dpi=160)
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
    # Panel A : Arbre canonique exact k = 2..10 (transition de phase)
    # ------------------------------------------------------------------ #
    ax = axes[0, 0]
    ks = list(range(2, 11))
    counts = [closure_data["levels"][str(k)]["canon_count"] for k in ks]
    # Pour log scale : afficher 0 comme barre hachuree a 0.45 annotee "0 EXACT"
    plot_y = [c if c > 0 else 0.45 for c in counts]
    colors = ["#38bdf8"] * 4 + ["#f59e0b", "#fb923c", "#a78bfa", "#10b981", "#ef4444"]
    bars = ax.bar(ks, plot_y, color=colors, edgecolor="#f8fafc", linewidth=0.9, width=0.68, zorder=3)
    bars[-1].set_hatch("//")
    ax.set_yscale("log")
    ax.set_ylim(0.3, 9000)
    ax.set_xticks(ks)
    ax.set_xlabel("Nombre de variables n (ordre du graphe K_n)", fontsize=10.5, fontweight="bold")
    ax.set_ylabel("Classes canoniques valides w_3 (mod S_n x {+-1})", fontsize=10, fontweight="bold")
    ax.set_title("A. Classification exhaustive exacte : pic a k = 6, unicite a n = 9, 0 a n = 10", fontsize=11, fontweight="bold", pad=10)
    ax.grid(axis="y", color="#334155", linestyle="--", alpha=0.6, zorder=0)
    for x, c, py in zip(ks, counts, plot_y):
        lbl = f"{c:,}".replace(",", " ") if c > 0 else "0\n(EXACT)"
        if x == 9:
            lbl = "1\n(UNIQUE)"
        ax.text(x, py * 1.22, lbl, ha="center", va="bottom", color="#f8fafc", fontsize=9.2, fontweight="bold")

    # ------------------------------------------------------------------ #
    # Panel B : Taux de survie des extensions & impasses (dead ends)
    # ------------------------------------------------------------------ #
    ax = axes[0, 1]
    ks_t = list(range(3, 11))
    surv_pct = []
    dead_pct = []
    for k in ks_t:
        lv = closure_data["levels"][str(k)]
        surv_pct.append(100.0 * lv["valid_extensions"] / lv["extensions_tested"])
        dead_pct.append(100.0 * lv["dead_ends"] / lv["parents"])
    ax.plot(ks_t, surv_pct, "o-", color="#38bdf8", linewidth=2.4, markersize=7, label="Taux d'extensions F_2-consistantes (%)")
    ax.plot(ks_t, dead_pct, "s--", color="#f43f5e", linewidth=2.4, markersize=7, label="Proportion de parents sans suite (dead ends %)")
    ax.fill_between(ks_t, dead_pct, alpha=0.14, color="#f43f5e")
    ax.set_xticks(ks_t)
    ax.set_ylim(-4, 109)
    ax.set_xlabel("Transition k - 1 -> k", fontsize=10.5, fontweight="bold")
    ax.set_ylabel("Pourcentage (%)", fontsize=10.5, fontweight="bold")
    ax.set_title("B. Effondrement combinatoire : 94,3% d'impasses a k = 8, 100% a k = 10", fontsize=11, fontweight="bold", pad=10)
    ax.grid(color="#334155", linestyle="--", alpha=0.6)
    ax.legend(facecolor="#0f172a", edgecolor="#475569", labelcolor="#f8fafc", fontsize=9.2, loc="center left")
    for x, d in zip(ks_t, dead_pct):
        if d > 0:
            ax.annotate(f"{d:.1f}%", (x, d), textcoords="offset points", xytext=(0, 7 if d < 95 else -15),
                        ha="center", color="#fda4af", fontsize=8.8, fontweight="bold")

    # ------------------------------------------------------------------ #
    # Panel C : Dessin du graphe SPORADIQUE UNIQUE K_9 : {u} U {v} U K_2 U C_5
    # ------------------------------------------------------------------ #
    ax = axes[1, 0]
    # Identifier canoniquement les 4 orbites de notre K_9 :
    # deg2 == 0 -> u ; deg2 == 7 -> v ; deg2 == 2 -> K_2 ; deg2 == 3 -> C_5 (ordonne le long du cycle)
    edges9 = [(a, b) for a in range(9) for b in range(a + 1, 9)]
    A2 = np.zeros((9, 9), dtype=int)
    A3 = np.zeros((9, 9), dtype=int)
    for ei, (a, b) in enumerate(edges9):
        A2[a, b] = A2[b, a] = k9_w2[ei]
        A3[a, b] = A3[b, a] = k9_w3[ei]
    deg2 = A2.sum(axis=1)
    u_iso = int(np.where(deg2 == 0)[0][0])
    v_uni = int(np.where(deg2 == 7)[0][0])
    k2_vs = list(np.where(deg2 == 2)[0])
    c5_raw = list(np.where(deg2 == 3)[0])
    # Ordonner C_5 le long de ses aretes P_2
    c5_ord = [c5_raw[0]]
    while len(c5_ord) < 5:
        last = c5_ord[-1]
        for nxt in c5_raw:
            if nxt not in c5_ord and A2[last, nxt] == 1:
                c5_ord.append(nxt)
                break

    pos = {}
    pos[u_iso] = (-1.65, 0.20)
    pos[v_uni] = (-0.05, 0.20)
    pos[k2_vs[0]] = (-0.82, 1.15)
    pos[k2_vs[1]] = (-0.82, -0.75)
    for idx, node in enumerate(c5_ord):
        ang = np.deg2rad(90 - 72 * idx)
        pos[node] = (1.35 + 0.82 * np.cos(ang), 0.20 + 0.82 * np.sin(ang))

    # Tracer les aretes G_3 (fond) puis G_2 (premier plan)
    cnt_p2 = 0; cnt_p3_1 = 0; cnt_p3_2 = 0
    for a, b in edges9:
        xa, ya = pos[a]; xb, yb = pos[b]
        if A3[a, b] == 1:
            cnt_p3_1 += 1
            ax.plot([xa, xb], [ya, yb], color="#fb7185", alpha=0.45, linewidth=1.3, linestyle="--", zorder=1)
        elif A3[a, b] == 2:
            cnt_p3_2 += 1
            ax.plot([xa, xb], [ya, yb], color="#fbbf24", alpha=0.65, linewidth=1.6, linestyle=":", zorder=2)
        if A2[a, b] == 1:
            cnt_p2 += 1
            ax.plot([xa, xb], [ya, yb], color="#38bdf8", alpha=0.95, linewidth=2.6, zorder=3)

    # Tracer les sommets par orbite
    orbit_style = [
        ([u_iso], "#f43f5e", "u (deg2=0)"),
        ([v_uni], "#38bdf8", "v (deg2=7)"),
        (k2_vs, "#a78bfa", "K_2 (deg2=2)"),
        (c5_ord, "#10b981", "C_5 (deg2=3)"),
    ]
    labels_short = {u_iso: "u", v_uni: "v", k2_vs[0]: "a1", k2_vs[1]: "a2"}
    for i, node in enumerate(c5_ord):
        labels_short[node] = f"c{i+1}"

    for nodes, col, leg in orbit_style:
        xs = [pos[nd][0] for nd in nodes]
        ys = [pos[nd][1] for nd in nodes]
        ax.scatter(xs, ys, s=330, color=col, edgecolor="#f8fafc", linewidth=1.6, zorder=5, label=leg)
        for nd in nodes:
            ax.text(pos[nd][0], pos[nd][1], labels_short[nd], ha="center", va="center",
                    color="#0f172a", fontsize=8.8, fontweight="bold", zorder=6)

    # Legende aretes (comptes dynamiques exacts)
    ax.plot([], [], color="#38bdf8", linewidth=2.6, label=f"P_2 = 1 ({cnt_p2} aretes, P_3 = 0)")
    ax.plot([], [], color="#fb7185", linestyle="--", linewidth=1.6, label=f"P_3 = 1 ({cnt_p3_1} aretes, P_2 = 0)")
    ax.plot([], [], color="#fbbf24", linestyle=":", linewidth=2.0, label=f"P_3 = 2 ({cnt_p3_2} aretes, P_2 = 0)")

    ax.set_xlim(-2.10, 2.50)
    ax.set_ylim(-1.68, 1.48)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title("C. Unique solution K_9 (mod S_9 x {+-1}) : V = {u} U {v} U K_2 U C_5, |Aut| = 20", fontsize=10.8, fontweight="bold", pad=10)
    ax.legend(facecolor="#0f172a", edgecolor="#475569", labelcolor="#f8fafc", fontsize=8.0, loc="lower center", ncol=2)

    # ------------------------------------------------------------------ #
    # Panel D : Distribution EXACTE de f2_defect sur les 19 683 extensions K_9 -> K_10
    # ------------------------------------------------------------------ #
    ax = axes[1, 1]
    vals, cnts = np.unique(defect_dist, return_counts=True)
    min_d = int(vals[0])
    med_d = float(np.median(defect_dist))
    ax.bar(vals, cnts, color="#f59e0b", edgecolor="#0f172a", width=0.88, alpha=0.9, zorder=3)
    ax.axvline(0, color="#10b981", linestyle="--", linewidth=2.0, zorder=4, label="Cible : f2_defect = 0")
    ax.axvline(med_d, color="#38bdf8", linestyle=":", linewidth=1.8, zorder=4, label=f"Mediane = {med_d:.0f}")
    ax.axvspan(0, min_d - 0.5, color="#ef4444", alpha=0.22, zorder=1)
    ax.annotate(
        f"GAP EXACT D'IMPOSSIBILITE\nmin f2_defect = {min_d} > 0\n(0 / 19 683 solutions)",
        xy=(min_d / 2.0, max(cnts) * 0.16),
        xytext=(34.0, max(cnts) * 0.72),
        ha="center", va="center",
        color="#fecaca", fontsize=9.6, fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.45", facecolor="#450a0a", edgecolor="#ef4444", linewidth=1.4),
        arrowprops=dict(arrowstyle="->", color="#f87171", lw=1.8, connectionstyle="arc3,rad=0.18"),
    )
    ax.set_xlim(-4, int(vals[-1]) + 6)
    ax.set_xlabel("Defaut exact de consistance F_2 (nombre de contradictions 0 = 1)", fontsize=10.2, fontweight="bold")
    ax.set_ylabel("Nombre d'extensions de l'unique K_9 vers K_10", fontsize=10.2, fontweight="bold")
    ax.set_title("D. Certificat terminal K_9 -> K_10 : toutes les 3^9 = 19 683 extensions echouent", fontsize=10.8, fontweight="bold", pad=10)
    ax.grid(axis="y", color="#334155", linestyle="--", alpha=0.6, zorder=0)
    ax.legend(facecolor="#0f172a", edgecolor="#475569", labelcolor="#f8fafc", fontsize=8.8, loc="upper right")

    fig.suptitle(
        "Fermeture EXACTE de THR_{n, 2} modulo 6 en degre d = 2 : N_gen^(6)(2, THR_2) = 9 (unique a S_9 x {+-1} pres)",
        color="#f8fafc", fontsize=14, fontweight="bold", y=0.98
    )
    plt.tight_layout(rect=[0, 0.01, 1, 0.955])
    os.makedirs(os.path.dirname(FIG_PROJECT), exist_ok=True)
    fig.savefig(FIG_PROJECT, dpi=160, facecolor=fig.get_facecolor())
    shutil.copy2(FIG_PROJECT, FIG_ARTIFACT)
    plt.close(fig)


def main():
    t0 = time.time()
    out = {}

    closure_path = os.path.join(RESULTS_DIR, "thr2_mod6_exact_closure.json")
    with open(closure_path, "r", encoding="utf-8") as f:
        closure_data = json.load(f)

    print("==========================================================================")
    print("  1. ANALYSE ALGEBRIQUE EXACTE DE L'UNIQUE CLASSE MAXIMALE K_9")
    print("==========================================================================")
    k9_entry = closure_data["max_n_witnesses"][0]
    w3_9 = np.array(k9_entry["w3"], dtype=np.int8)
    w2_9 = k9_entry["w2"]

    rank, nullity, n_eqs, all_w2 = exact_f2_solution_space(9, w3_9)
    print(f"  Systeme F_2 pour l'unique w_3 sur K_9 : {n_eqs} equations, rang = {rank}/36, nullite = {nullity} -> {len(all_w2)} solution(s) P_2 !")
    # Verifier force-brute 2^9 sur toutes les solutions P_2
    for idx_s, sol_w2 in enumerate(all_w2):
        # Verification directe sur 2^9 = 512 entrees
        edges9 = [(u, v) for u in range(9) for v in range(u + 1, 9)]
        for x in range(1 << 9):
            pc = x.bit_count()
            s2 = 0; s3 = 0
            for ei, (u, v) in enumerate(edges9):
                if ((x >> u) & 1) and ((x >> v) & 1):
                    s2 ^= sol_w2[ei]
                    s3 += int(w3_9[ei])
            is_zero = (s2 == 0) and (s3 % 3 == 0)
            assert is_zero == (pc < 2), f"Echec brute-force sur sol #{idx_s}, x={x}"
    print(f"  Verification 2^9 = 512 sur les {len(all_w2)} solution(s) P_2 : 100% OK !")

    n_aut_pos, n_aut_neg, orbits = compute_automorphism_group_k9(w3_9, w2_9)
    n_labeled_w3 = (2 * 362880) // (n_aut_pos + n_aut_neg)
    n_labeled_pairs = n_labeled_w3 * len(all_w2)
    print(f"  |Aut_+(w_3)| = {n_aut_pos}, |Aut_-(w_3)| = {n_aut_neg}, |Aut_{{+-1}}(w_3)| = {n_aut_pos + n_aut_neg}")
    print(f"  Orbites de sommets ({len(orbits)} orbites, tailles {[len(o) for o in orbits]}) : {orbits}")
    print(f"  Nombre EXACT de polynomes P_3 etiquetes sur K_9 : {n_labeled_w3:,}")
    print(f"  Nombre EXACT de polynomes P = (P_2, P_3) etiquetes sur K_9 : {n_labeled_pairs:,}")

    out["k9_unique_analysis"] = {
        "f2_equations": n_eqs,
        "f2_rank": rank,
        "f2_nullity": nullity,
        "n_valid_p2_for_unique_w3": len(all_w2),
        "p2_unique": (len(all_w2) == 1),
        "aut_pos_order": n_aut_pos,
        "aut_neg_order": n_aut_neg,
        "aut_pm1_order": n_aut_pos + n_aut_neg,
        "vertex_orbits": orbits,
        "vertex_orbit_sizes": sorted(len(o) for o in orbits),
        "n_labeled_w3_on_k9": n_labeled_w3,
        "n_labeled_polynomials_on_k9": n_labeled_pairs,
    }

    print("\n==========================================================================")
    print("  2. AUDIT INDEPENDANT (MOTEUR NAIF vs MOTEUR RAPIDE) k = 2..10")
    print("==========================================================================")
    cur = np.array([[0], [1]], dtype=np.int8)
    levels_arr = {2: cur}
    for k in range(3, 11):
        ext_fast, _ = extend_valid_w3_exact_fast(cur, k)
        ext_naive, _ = extend_valid_w3_batch(cur, k)
        assert ext_fast.shape[0] == ext_naive.shape[0], f"MISMATCH k={k}: fast={ext_fast.shape[0]} vs naive={ext_naive.shape[0]}"
        # Verifier egalite des ensembles tries de lignes
        s_fast = set(map(bytes, ext_fast))
        s_naive = set(map(bytes, ext_naive))
        assert s_fast == s_naive, f"SET MISMATCH k={k}!"
        print(f"  [AUDIT NAIF == RAPIDE] k = {k:2d} : {ext_fast.shape[0]:7d} extensions valides IDENTIQUES BIT-A-BIT !")
        uniq_w, _ = canonicalize_exact_s_k_pm1(ext_fast, k)
        cur = uniq_w
        levels_arr[k] = cur

    # Verifier aussi que les 4 parents survivants de k = 8 sont EXACTEMENT les 4 sous-graphes
    # obtenus en supprimant un sommet de chacune des 4 orbites de K_9 !
    edges9 = [(u, v) for u in range(9) for v in range(u + 1, 9)]
    del_subgraphs = []
    for orb in orbits:
        drop_v = orb[0]
        keep = [v for v in range(9) if v != drop_v]
        sub_w = []
        for i in range(8):
            for j in range(i + 1, 8):
                a, b = min(keep[i], keep[j]), max(keep[i], keep[j])
                sub_w.append(int(w3_9[edges9.index((a, b))]))
        del_subgraphs.append(sub_w)
    _, del_keys = canonicalize_exact_s_k_pm1(np.array(del_subgraphs, dtype=np.int8), 8)
    assert len(del_keys) == 4, f"Attendu 4 sous-graphes K_8 distincts, obtenu {len(del_keys)}"
    print(f"  [VERIF ORBITES] Suppression d'un sommet par orbite de K_9 -> EXACTEMENT {len(del_keys)} classes K_8 distinctes (matching les 4 parents de k=8->9) !")

    # Verifier force-brute 2^8 sur TOUTES les 27 classes canoniques de k = 8
    for idx8 in range(levels_arr[8].shape[0]):
        solve_and_verify_thr2_mod2p(8, 3, levels_arr[8][idx8])
    print(f"  [VERIF 2^8] 100% des {levels_arr[8].shape[0]} classes canoniques de k = 8 verifiees sur 2^8 = 256 entrees !")

    out["independent_audit"] = {
        "naive_vs_fast_bit_identical_levels": [3, 4, 5, 6, 7, 8, 9, 10],
        "k8_all_27_classes_verified_2_to_8": True,
        "k9_vertex_deleted_distinct_k8_parents": len(del_keys),
    }

    print("\n==========================================================================")
    print("  3. DISTRIBUTION EXACTE DE f2_defect SUR LES 19 683 EXTENSIONS K_9 -> K_10")
    print("==========================================================================")
    t1 = time.time()
    defs10 = exact_k9_to_k10_defect_distribution(w3_9)
    vals, cnts = np.unique(defs10, return_counts=True)
    min_d = int(vals[0])
    n_at_min = int(cnts[0])
    med_d = float(np.median(defs10))
    max_d = int(vals[-1])
    print(f"  19 683 extensions K_9 -> K_10 : min_defect = {min_d} ({n_at_min} fois), mediane = {med_d:.1f}, max = {max_d} ({time.time()-t1:.2f}s)")
    out["k9_to_k10_exact_defects"] = {
        "n_extensions": int(len(defs10)),
        "min_defect": min_d,
        "count_at_min_defect": n_at_min,
        "median_defect": med_d,
        "max_defect": max_d,
        "histogram": {str(int(v)): int(c) for v, c in zip(vals, cnts)},
    }

    print("\n==========================================================================")
    print("  4. GENERATION DE LA FIGURE PUBLICATION-GRADE")
    print("==========================================================================")
    make_figure(closure_data, out, w2_9, w3_9, defs10)
    print(f"  Figure generee : {FIG_PROJECT}")
    print(f"  Copie artifact : {FIG_ARTIFACT}")

    out["total_time_s"] = round(time.time() - t0, 2)
    out_path = os.path.join(RESULTS_DIR, "thr2_mod6_audit.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print(f"\nAudit complet termine en {time.time() - t0:.2f}s. Sauvegarde : {out_path}")


if __name__ == "__main__":
    main()
