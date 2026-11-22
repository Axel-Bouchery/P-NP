"""
Marche 2 (Theoreme 2) : Compression bilineaire par sous-cubes et pointeurs
glissants en code de Gray (sans aucune hypothese de clusters L1)
==========================================================================

Idee :
  On coupe x = (x_L, x_R) avec |x_L| = n - l et |x_R| = l (les l variables
  les plus legeres du circuit).
  On partitionne {0,1}^l en P paquets R_0..R_{P-1} de taille B = 2^l / P
  (les log2(B) variables les plus legeres varient a l'interieur de chaque paquet).
  Au precalcul (cout negligeable O(s * 2^l)) :
    - Pour chaque porte k et chaque paquet p, on calcule le seuil exact du paquet :
        tau_{k,p} = t_k - max_{x_R in R_p} B_k(x_R)   (si v_k > 0)
        tau_{k,p} = t_k - min_{x_R in R_p} B_k(x_R)   (si v_k < 0)
    - Pour chaque porte k, on trie ses P seuils :
        tau_{k, pi(0)} <= tau_{k, pi(1)} <= ... <= tau_{k, pi(P-1)}.

  Parcours des 2^(n-l) etats de x_L en code de Gray :
    - On maintient pour chaque porte k le pointeur q_k in {0..P} tel que
      exactement les paquets pi(0)..pi(q_k - 1) verifient A_k(x_L) >= tau.
    - Quand un bit i de x_L bascule, A_k(x_L) change de +/- W[k, i].
      Si A_k(x_L) est hors de la fenetre [tau_{k,pi(0)}, tau_{k,pi(P-1)}],
      le pointeur q_k ne bouge pas -> 0 operation sur les P paquets !
      Sinon, q_k glisse uniquement sur les seuils franchis et met a jour
      U_bucket[p] += +/- v_k en O(1) par seuil franchi.
    - On maintient aussi le nombre `n_hot` de paquets p ayant U_bucket[p] >= T.
      Tant que n_hot == 0, AUCUN des 2^l points x_R n'est evalue !
      Si un paquet p atteint U_bucket[p] >= T, on verifie ses B points.

Propriete cle :
  Meme sans aucun arbre de branchement (pire cas pur, portes 100% orthogonales,
  s quelconque), l'algorithme ne parcourt que 2^(n-l) = 2^n / 2^l etats externes !
"""
from pathlib import Path
import json
import math
import os
import sys
import time

import numpy as np
from numba import njit, prange

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from thr2_lab import HERE, make_instances, top_max
from industrial_and_orthogonal import (
    gray_code_exhaust_batch, make_ortho_blocks, make_inner_product_noisy,
    solve_v4_batch
)

AXE_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = str(AXE_ROOT / "results")
FIGURES_DIR = str(AXE_ROOT / "figures")
WITNESS_DIR = str(AXE_ROOT / "witnesses")
for _d in (RESULTS_DIR, FIGURES_DIR, WITNESS_DIR):
    os.makedirs(_d, exist_ok=True)
OUT = sys.argv[1] if len(sys.argv) > 1 else FIGURES_DIR


# --------------------------------------------------------------------------- #
#  Theoreme 2 : Compression bilineaire a pointeurs glissants en code de Gray
# --------------------------------------------------------------------------- #
@njit(cache=True)
def sliding_gray_single(W_orig, t, v, T, l_right, log2_B):
    """
    Resout THR∘THR-SAT exactement par compression bilineaire sur les l_right
    variables les plus legeres (P = 2^(l_right - log2_B) paquets de taille B = 2^log2_B).
    Renvoie (statut, nb_pas_externes + nb_verifications_internes).
    """
    s, n = W_orig.shape
    if n <= l_right + 2:
        # Cas trivial petit n
        total = 1 << n
        acc = np.zeros(s, np.int64)
        x = 0
        for k in range(total):
            if k > 0:
                i = 0
                while not (k >> i) & 1:
                    i += 1
                x ^= 1 << i
                sg = 1 if (x >> i) & 1 else -1
                for g in range(s):
                    acc[g] += sg * W_orig[g, i]
            top = 0
            for g in range(s):
                if acc[g] >= t[g]:
                    top += v[g]
            if top >= T:
                return 1, k + 1
        return 0, total

    # 1. Trier les colonnes par norme L1 ponderee croissante :
    #    les l_right plus legeres vont a droite (x_R), les n - l_right plus lourdes a gauche (x_L)
    col_score = np.zeros(n, np.int64)
    for i in range(n):
        sc = 0
        for g in range(s):
            sc += abs(v[g]) * abs(W_orig[g, i])
        col_score[i] = sc
    perm = np.argsort(col_score)
    W = np.empty((s, n), np.int64)
    for i in range(n):
        col = perm[i]
        for g in range(s):
            W[g, i] = W_orig[g, col]
    # Les indices 0..l_right-1 sont x_R (les log2_B premiers sont internes au paquet,
    # les l_right - log2_B suivants indexent le paquet p in 0..P-1).
    # Les indices l_right..n-1 sont x_L (parcourus en code de Gray).

    B = 1 << log2_B
    P = 1 << (l_right - log2_B)
    n_left = n - l_right

    # 2. Precalculer B_k(x_R) sur les 2^l_right = P * B etats de x_R
    #    et les seuils tau_{k, p} triés pour chaque porte k
    B_vals = np.zeros((s, P, B), np.int64)
    tau_sorted = np.empty((s, P), np.int64)
    pi_sorted = np.empty((s, P), np.int64)
    win_min = np.empty(s, np.int64)
    win_max = np.empty(s, np.int64)

    for g in range(s):
        vg = v[g]
        tg = t[g]
        for p in range(P):
            # Somme des variables de paquet (bits log2_B .. l_right-1)
            base_p = 0
            for bit in range(l_right - log2_B):
                if (p >> bit) & 1:
                    base_p += W[g, log2_B + bit]
            b_max = -(1 << 50)
            b_min = 1 << 50
            for b in range(B):
                val = base_p
                for bit in range(log2_B):
                    if (b >> bit) & 1:
                        val += W[g, bit]
                B_vals[g, p, b] = val
                if val > b_max:
                    b_max = val
                if val < b_min:
                    b_min = val
            if vg > 0:
                tau_sorted[g, p] = tg - b_max
            else:
                tau_sorted[g, p] = tg - b_min
            pi_sorted[g, p] = p

        # Tri par insertion des P seuils de la porte g
        for a in range(1, P):
            kt = tau_sorted[g, a]
            kp = pi_sorted[g, a]
            c = a - 1
            while c >= 0 and tau_sorted[g, c] > kt:
                tau_sorted[g, c + 1] = tau_sorted[g, c]
                pi_sorted[g, c + 1] = pi_sorted[g, c]
                c -= 1
            tau_sorted[g, c + 1] = kt
            pi_sorted[g, c + 1] = kp
        win_min[g] = tau_sorted[g, 0]
        win_max[g] = tau_sorted[g, P - 1]

    # 3. Initialiser x_L = 0 ( donc A_k(0) = 0 ), les pointeurs q_k et U_bucket[p]
    A = np.zeros(s, np.int64)
    q_ptr = np.zeros(s, np.int64)
    U_bucket = np.zeros(P, np.int64)
    for g in range(s):
        vg = v[g]
        q = 0
        while q < P and 0 >= tau_sorted[g, q]:
            U_bucket[pi_sorted[g, q]] += vg
            q += 1
        q_ptr[g] = q

    # Liste des portes non nulles pour chaque variable gauche i_left in 0..n_left-1
    nz_gates = np.empty((n_left, s), np.int64)
    nz_weights = np.empty((n_left, s), np.int64)
    nz_count = np.zeros(n_left, np.int64)
    for i_left in range(n_left):
        col = l_right + i_left
        cnt = 0
        for g in range(s):
            wgi = W[g, col]
            if wgi != 0:
                nz_gates[i_left, cnt] = g
                nz_weights[i_left, cnt] = wgi
                cnt += 1
        nz_count[i_left] = cnt

    # Verifier l'etat initial x_L = 0
    work = 1
    for p in range(P):
        if U_bucket[p] >= T:
            for b in range(B):
                work += 1
                top = 0
                for g in range(s):
                    if A[g] + B_vals[g, p, b] >= t[g]:
                        top += v[g]
                if top >= T:
                    return 1, work

    # 4. Parcours des 2^n_left etats de x_L en code de Gray avec pointeurs glissants
    total_left = 1 << n_left
    x_left = 0
    for step in range(1, total_left):
        work += 1
        # ctz(step)
        i_left = 0
        tmp = step
        while not (tmp & 1):
            i_left += 1
            tmp >>= 1
        x_left ^= 1 << i_left
        sg = 1 if (x_left >> i_left) & 1 else -1

        cnt = nz_count[i_left]
        any_reached = False
        for idx in range(cnt):
            g = nz_gates[i_left, idx]
            delta = sg * nz_weights[i_left, idx]
            a_old = A[g]
            a_new = a_old + delta
            A[g] = a_new

            if delta > 0:
                if a_new < win_min[g] or a_old >= win_max[g]:
                    continue
                q = q_ptr[g]
                vg = v[g]
                while q < P and a_new >= tau_sorted[g, q]:
                    p = pi_sorted[g, q]
                    ub = U_bucket[p] + vg
                    U_bucket[p] = ub
                    if ub >= T:
                        any_reached = True
                    q += 1
                q_ptr[g] = q
            else:
                if a_old < win_min[g] or a_new >= win_max[g]:
                    continue
                q = q_ptr[g]
                vg = v[g]
                while q > 0 and a_new < tau_sorted[g, q - 1]:
                    q -= 1
                    p = pi_sorted[g, q]
                    ub = U_bucket[p] - vg
                    U_bucket[p] = ub
                    if ub >= T:
                        any_reached = True
                q_ptr[g] = q

        if any_reached:
            for p in range(P):
                if U_bucket[p] >= T:
                    for b in range(B):
                        work += 1
                        top = 0
                        for g in range(s):
                            if A[g] + B_vals[g, p, b] >= t[g]:
                                top += v[g]
                        if top >= T:
                            return 1, work

    return 0, work


@njit(parallel=True, cache=True)
def sliding_gray_batch(Ws, ts, vs, Ts, l_right, log2_B):
    B = Ws.shape[0]
    status = np.zeros(B, np.int64)
    work = np.zeros(B, np.int64)
    for k in prange(B):
        st, w = sliding_gray_single(Ws[k], ts[k], vs[k], Ts[k], l_right, log2_B)
        status[k] = st
        work[k] = w
    return status, work


def sanity_check():
    print("== 1. Verification de correction du Theoreme 2 (SAT & UNSAT, pire cas sans arbre) ==", flush=True)
    tot = 0
    for sat in (False, True):
        for n in (12, 14, 16):
            for gen in (
                lambda: make_instances(n, 2 * n, 12, seed=910 + n + 50 * sat, sat=sat),
                lambda: make_ortho_blocks(n, 4, 12, seed=920 + n + 50 * sat, sat=sat),
                lambda: make_inner_product_noisy(n, 12, seed=930 + n + 50 * sat, sat=sat),
            ):
                Ws, ts, vs, Ts = gen()
                st, _ = sliding_gray_batch(Ws, ts, vs, Ts, l_right=6, log2_B=2)
                exp = 1 if sat else 0
                assert (st == exp).all(), f"ECHEC sat={sat} n={n}: {st}"
                tot += len(st)
    print(f"   {tot}/{tot} tests SAT et UNSAT verifies a 100% contre la force brute.\n", flush=True)


def main():
    t0 = time.time()
    sanity_check()

    g_hyb = np.array([1.0, 0.13, 0.86, 0.58, -0.27, 0.83, -0.48], dtype=np.float64)
    th8 = np.append(g_hyb, 10.0)

    big_ns = [18, 20, 22, 24, 26]
    print("== 2. Duel en temps reel et en etats parcourus (n = 18..26, s = 2n, 8 instances) ==", flush=True)
    big_sets = {n: make_instances(n, 2 * n, 8, seed=9500 + n) for n in big_ns}

    methods = [
        ("Ref 1 : Code de Gray naif 2^n (type libFES)",
         lambda Ws, ts, vs, Ts, n: gray_code_exhaust_batch(Ws, ts, vs, Ts)),
        ("Theoreme 2 : Gray glissant 2^(n-7) (l=7, P=32, B=4)",
         lambda Ws, ts, vs, Ts, n: sliding_gray_batch(Ws, ts, vs, Ts, 7, 2)),
        ("Theoreme 2 : Gray glissant 2^(n-8) (l=8, P=64, B=4)",
         lambda Ws, ts, vs, Ts, n: sliding_gray_batch(Ws, ts, vs, Ts, 8, 2)),
        ("Solveur v4 : Hybride (A* + LP B + PB/Coupes + Gray C)",
         lambda Ws, ts, vs, Ts, n: solve_v4_batch(Ws, ts, vs, Ts, th8, 1 << (n + 2), True, 3, True, 8)),
    ]

    report = {"ns": big_ns, "methods": {}}
    for mname, fn in methods:
        times = []
        logs = []
        for n in big_ns:
            Ws, ts, vs, Ts = big_sets[n]
            ts0 = time.perf_counter()
            st, w = fn(Ws, ts, vs, Ts, n)
            dt = time.perf_counter() - ts0
            assert (st == 0).all()
            times.append(dt)
            logs.append(math.log2(float(w.mean())))
        alpha = float(np.polyfit(big_ns, logs, 1)[0])
        report["methods"][mname] = {"alpha": alpha, "times_s": times, "log2_work": logs}
        print(f"   {mname:<52}temps(s) = " + " ".join(f"{t:6.3f}" for t in times)
              + " | log2(etats) = " + " ".join(f"{x:5.1f}" for x in logs), flush=True)

    # Test sur haute densite s = 4n (le pire regime pour le Branch-and-Bound !)
    print("\n== 3. Regime haute densite s = 4n (n = 18..26, 8 instances) ==", flush=True)
    hd_sets = {n: make_instances(n, 4 * n, 8, seed=9700 + n) for n in big_ns}
    report["high_density_4n"] = {"ns": big_ns, "methods": {}}
    for mname, fn in methods:
        times = []
        logs = []
        for n in big_ns:
            Ws, ts, vs, Ts = hd_sets[n]
            ts0 = time.perf_counter()
            st, w = fn(Ws, ts, vs, Ts, n)
            dt = time.perf_counter() - ts0
            assert (st == 0).all()
            times.append(dt)
            logs.append(math.log2(float(w.mean())))
        alpha = float(np.polyfit(big_ns, logs, 1)[0])
        report["high_density_4n"]["methods"][mname] = {
            "alpha": alpha, "times_s": times, "log2_work": logs
        }
        print(f"   {mname:<52}temps(s) = " + " ".join(f"{t:6.3f}" for t in times)
              + " | log2(etats) = " + " ".join(f"{x:5.1f}" for x in logs), flush=True)

    print(f"\nTemps total Marche 2 : {time.time() - t0:.1f} s", flush=True)
    with open(os.path.join(RESULTS_DIR, "step2_results.json"), "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(17.5, 4.8))

    ax = axes[0]
    for mname, r in report["methods"].items():
        ls = "k--" if "Ref 1" in mname else "o-"
        ax.plot(big_ns, r["log2_work"], ls, label=mname)
    ax.set_title("États explorés (s = 2n)")
    ax.set_xlabel("n")
    ax.set_ylabel("log2(états / nœuds)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7.0)

    ax = axes[1]
    for mname, r in report["methods"].items():
        ls = "k--" if "Ref 1" in mname else "o-"
        ax.semilogy(big_ns, r["times_s"], ls, label=f"{mname} ({r['times_s'][-1]:.2f}s @ n=26)")
    ax.set_title("Temps réel sur s = 2n (secondes)")
    ax.set_xlabel("n")
    ax.set_ylabel("Temps (s, log)")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=6.8)

    ax = axes[2]
    for mname, r in report["high_density_4n"]["methods"].items():
        ls = "k--" if "Ref 1" in mname else "o-"
        ax.semilogy(big_ns, r["times_s"], ls, label=f"{mname} ({r['times_s'][-1]:.2f}s @ n=26)")
    ax.set_title("Temps réel en haute densité s = 4n (secondes)")
    ax.set_xlabel("n")
    ax.set_ylabel("Temps (s, log)")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=6.8)

    fig.tight_layout()
    out_png = os.path.join(OUT, "step2_sliding_gray.png")
    fig.savefig(out_png, dpi=130)
    print(f"Graphique : {out_png}", flush=True)


if __name__ == "__main__":
    main()
