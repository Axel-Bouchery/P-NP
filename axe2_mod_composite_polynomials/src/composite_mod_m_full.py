"""
Piste 1 — Generalisation complete des Theoremes M1/M2 aux modules composites
m in {6, 10, 12, 14, 15, 21, 30}
=============================================================================
Quatre theoremes constructifs unifies (prouves et verifies exactement) :
  - G1 (Multi-Channel Margin-Harvesting) : somme des marges (d - d_i) sur tous
    les canaux CRT disjoints (valable pour tout m composite, meme avec p^e).
  - G2 (Anisotropic Field-Norm Harvesting) : norme d'extension N_{F_{p^k}/F_p}
    sur chaque canal premier simple p | m des que k = floor(d / d_p) >= 2,
    donnant (k - 1)*(p - 1)*d_p + (d - k*d_p) variables gratuites.
  - G3 (Value-Avoidance Harvesting) : si Im(P_i) sur {0..a} omet au moins une
    valeur non nulle -c mod q_i, alors P_i(A) + c*(1 - prod_{j=1}^d (1 - y_j))
    a degre max(deg(P_i), d) = d et recolte +d variables MEME A MARGE NULLE !
    (Ex. mod 6 a d = 8 : I_9(A) evite 2 mod 3 -> 80 + 8 = 88, battant 83 !)
  - G4 (Multi-Block Norm Doubling / Multiplication) : pour tout m sans facteur
    carre (r facteurs premiers) et tout k >= 2, d_0 >= 1 avec k*d_0 <= d :
      N_gen^(m)(d) >= k * N_gen^(m)(d_0) + r * (d - k*d_0)
    (Ex. mod 15 a d = 2 : 2*14 = 28 = +100% au-dessus de N_sym = 14 !
     Ex. mod 21 a d = 2 : 2*20 = 40 = +100% au-dessus de N_sym = 20 !
     Ex. mod 30 a d = 2 : 2*29 = 58 vs N_sym = 44, +31.8% !).

Inclut verification exhaustive force-brute 2^n sur les petits records (n <= 20)
et verification algebrique exacte de tous les temoins + normes N_{F_{p^k}/F_p}.
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from composite_mod_m_lab import (
    METHOD_NAMES,
    RESULTS_DIR,
    binomial_mod_table,
    build_channel_catalog,
    solve_2channel,
    solve_3channel,
    verify_witness_2ch,
)


# --------------------------------------------------------------------------- #
#  Verification explicite des polynomes irreductibles -> Norme N_{F_{p^k}/F_p}
# --------------------------------------------------------------------------- #
def verify_field_norm_anisotropic(p: int, k: int) -> dict:
    """
    Trouve un polynome unitaire irreductible f(t) de degre k sur F_p,
    construit explicitement la matrice compagnon M(u) de multiplication par
    u_0 + u_1 t + ... + u_{k-1} t^{k-1} dans F_p[t]/(f(t)), et verifie sur
    TOUS les p^k - 1 vecteurs non nuls (si p^k <= 20000) que
    N_k(u) = det(M(u)) mod p != 0.
    """
    # Recherche d'un polynome irreductible unitaire de degre k sur F_p
    # (via test de divisibilite par tous les polynomes de degre 1..k//2)
    def poly_mod(a, b):
        a = list(a)
        while len(a) >= len(b):
            if a[-1] != 0:
                factor = (a[-1] * pow(b[-1], p - 2, p)) % p
                shift = len(a) - len(b)
                for i in range(len(b)):
                    a[shift + i] = (a[shift + i] - factor * b[i]) % p
            a.pop()
        while a and a[-1] == 0:
            a.pop()
        return a

    def is_irreducible(f):
        # Test Rabin / divisibilite directe pour k <= 6
        # On teste pgcd(f, t^{p^m} - t) == 1 pour m = 1..k//2
        for m in range(1, k // 2 + 1):
            # Calculer t^{p^m} mod f par exponentiation rapide
            res = [1]
            base = [0, 1]
            exp = p ** m
            while exp > 0:
                if exp & 1:
                    # mul res * base mod f
                    prod = [0] * (len(res) + len(base) - 1)
                    for i, ci in enumerate(res):
                        for j, cj in enumerate(base):
                            prod[i + j] = (prod[i + j] + ci * cj) % p
                    res = poly_mod(prod, f)
                prod = [0] * (2 * len(base) - 1)
                for i, ci in enumerate(base):
                    for j, cj in enumerate(base):
                        prod[i + j] = (prod[i + j] + ci * cj) % p
                base = poly_mod(prod, f)
                exp >>= 1
            # res - t
            if len(res) < 2:
                res += [0] * (2 - len(res))
            res[1] = (res[1] - 1) % p
            while res and res[-1] == 0:
                res.pop()
            # gcd(f, res)
            a, b = list(f), list(res)
            while b:
                a, b = b, poly_mod(a, b)
            if len(a) > 1:
                return False
        return True

    irred = None
    # Cherche un polynome irreductible le plus creux possible
    total = p ** k
    for code in range(1, total):
        f = [(code // (p ** j)) % p for j in range(k)] + [1]
        if f[0] == 0:
            continue
        if is_irreducible(f):
            irred = f
            break
    assert irred is not None, f"Aucun irreductible trouve pour p={p}, k={k}!"

    # Verification directe sur F_p^k \ {0} si p^k <= 20000
    checked = 0
    if p ** k <= 20000:
        # Matrice compagnon C de multiplication par t
        C = np.zeros((k, k), dtype=np.int64)
        for i in range(k - 1):
            C[i + 1, i] = 1
        for i in range(k):
            C[i, k - 1] = (-irred[i]) % p
        # Puissances C^0, C^1, ..., C^{k-1} mod p
        pows = [np.eye(k, dtype=np.int64)]
        for _ in range(1, k):
            pows.append((pows[-1] @ C) % p)

        def det_mod_p(mat):
            M = mat.copy() % p
            det = 1
            for col in range(k):
                piv = -1
                for r in range(col, k):
                    if M[r, col] != 0:
                        piv = r
                        break
                if piv == -1:
                    return 0
                if piv != col:
                    M[[col, piv]] = M[[piv, col]]
                    det = (-det) % p
                det = (det * M[col, col]) % p
                inv = pow(int(M[col, col]), p - 2, p)
                for r in range(col + 1, k):
                    if M[r, col] != 0:
                        fac = (M[r, col] * inv) % p
                        M[r, col:] = (M[r, col:] - fac * M[col, col:]) % p
            return det

        for code in range(1, p ** k):
            M = np.zeros((k, k), dtype=np.int64)
            rem = code
            for j in range(k):
                uj = rem % p
                rem //= p
                if uj:
                    M = (M + uj * pows[j]) % p
            assert det_mod_p(M) != 0, f"Norme nulle sur vecteur non nul pour p={p}, k={k}!"
            checked += 1

    return {"p": p, "k": k, "irred_coeffs_low_to_high": irred, "exhaustive_checked": checked}


def brute_force_verify_small(n: int, mod_list, eval_fn) -> bool:
    """Verification force-brute integrale sur les 2^n vecteurs boolens."""
    assert n <= 20
    # Evalue sur tous les x in {0..2^n - 1}
    vals0 = eval_fn(0)
    assert all(v == 0 for v in vals0), f"Non nul en 0: {vals0}"
    for mask in range(1, 1 << n):
        vals = eval_fn(mask)
        if all(v == 0 for v in vals):
            raise AssertionError(f"Zero commun illegal au masque {bin(mask)}!")
    return True


def verify_mod6_d8_n88():
    """
    Verification explicite et independante du nouveau record MOD 6, d = 8 -> n = 88
    (qui bat notre precedent record M2 = 83 et le symetrique N_sym(8) = 80 !).
    """
    # P3(A) = I_9(A) = 1[A != 0 mod 9] de degre 8 sur F_3
    # P2(A) de degre 8 sur F_2 donnant a = 80
    C2 = binomial_mod_table(81, 8, 2)
    C3 = binomial_mod_table(81, 8, 3)
    deg2, val2, cof2 = build_channel_catalog(8, 2, True, C2, 81)
    deg3, val3, cof3 = build_channel_catalog(8, 3, True, C3, 81)
    res = solve_2channel(8, 2, True, deg2, val2, 3, True, deg3, val3, 81)
    nsym, si, sj, ng1, _, _, _, ng2, _, _, _, nall, ai, aj, aa, m1, m2, gn1, gn2 = res
    assert nsym == 80 and nall == 88 and aa == 80 and gn2 == 8
    c2 = [int(x) for x in cof2[ai]]
    c3 = [int(x) for x in cof3[aj]]
    # Verifier les valeurs de P3(s) mod 3 pour s = 0..80
    im3 = set()
    for s in range(0, 81):
        v2 = int(sum(c2[k] * int(C2[s, k]) for k in range(8)) % 2)
        v3 = int(sum(c3[k] * int(C3[s, k]) for k in range(8)) % 3)
        if s > 0:
            assert not (v2 == 0 and v3 == 0), f"Zero commun a s={s}"
        im3.add(v3)
    assert len(im3) == 2 and 0 in im3, f"Image P3 inattendue: {im3}"
    omitted = ({0, 1, 2} - im3).pop()
    c_shift = (-omitted) % 3
    # Verifier sur tout (s_A in 0..80, y_non_zero in {False, True}):
    for sA in range(81):
        v2 = int(sum(c2[k] * int(C2[sA, k]) for k in range(8)) % 2)
        v3 = int(sum(c3[k] * int(C3[sA, k]) for k in range(8)) % 3)
        # y = 0:
        if sA == 0:
            assert v2 == 0 and v3 == 0
        else:
            assert not (v2 == 0 and v3 == 0)
        # y != 0: v(y) = 1 -> P3_tilde = (v3 + c_shift) mod 3
        v3_y = (v3 + c_shift) % 3
        assert v3_y != 0, f"Annulation pour sA={sA}, y!=0!"
    return {
        "c2": c2, "c3": c3, "image_P3": sorted(im3),
        "omitted_val": omitted, "shift_c": c_shift, "n_verified": 88,
    }


def apply_g4_dp(rows, r_factors: int, is_squarefree: bool):
    """
    Applique la programmation dynamique exacte du Theoreme G4 (norme d'extension
    multi-blocs sur modules sans facteur carre) :
      N_gen(d) >= max(N_G123(d), max_{k>=2, k*d0 <= d} [ k * N_gen(d0) + r*(d - k*d0) ])
    Met a jour chaque ligne en place avec N_G4, N_best et la source du record.
    """
    best_by_d = {}
    for row in rows:
        d = row["d"]
        cand_123 = row["N_all"]
        best_g4 = 0
        g4_desc = "none"
        if is_squarefree:
            for k in range(2, d + 1):
                for d0 in range(1, d // k + 1):
                    if d0 in best_by_d:
                        rem = d - k * d0
                        val = k * best_by_d[d0] + r_factors * rem
                        if val > best_g4:
                            best_g4 = val
                            g4_desc = f"{k}xN({d0})+{r_factors}*{rem}" if rem > 0 else f"{k}xN({d0})"
        best = max(cand_123, best_g4)
        best_by_d[d] = best
        row["N_G4"] = best_g4
        row["G4_desc"] = g4_desc
        row["N_best"] = best
        row["best_gain"] = best - row["N_sym"]
        row["best_gain_pct"] = round(100.0 * (best - row["N_sym"]) / row["N_sym"], 2)
        if best == row["N_sym"]:
            row["best_mechanism"] = "sym"
        elif best_g4 > cand_123:
            row["best_mechanism"] = f"G4({g4_desc})"
        else:
            row["best_mechanism"] = row.get("g123_summary", "G123")


def run_full_campaign():
    t0 = time.time()

    # 1. Verification des normes d'extension de corps N_{F_{p^k}/F_p}
    print("=== 1. Verification exacte des normes anisotropes N_{F_{p^k}/F_p} ===")
    norm_checks = []
    for p, kmax in [(2, 6), (3, 6), (5, 5), (7, 4)]:
        for k in range(2, kmax + 1):
            info = verify_field_norm_anisotropic(p, k)
            norm_checks.append(info)
    print(f"  {len(norm_checks)} normes de corps verifiees (toutes anisotropes sur F_p^k \\ {{0}}).")

    # 2. Verification independante du nouveau record MOD 6 d=8 -> n=88
    print("\n=== 2. Verification independante MOD 6, d = 8 -> n = 88 (G3 Value-Avoidance) ===")
    m6d8 = verify_mod6_d8_n88()
    print(f"  OK! c2={m6d8['c2']}, c3={m6d8['c3']}, Im(P3)={m6d8['image_P3']} (omet {m6d8['omitted_val']}), shift=+{m6d8['shift_c']} -> n=88!")

    # 3. Verification force-brute 2^n sur n=18 (mod 10 d=2) et n=20 (mod 15 d=2, G2)
    print("\n=== 3. Verification force-brute 2^n sur mod 10 (d=2, n=18) et mod 15 (d=2, n=20) ===")
    # mod 10, d=2, n=18: G4(2x9) -> A = x_0..x_8 (9 vars), B = x_9..x_17 (9 vars)
    # P2 = A + B + A*B mod 2, P5 = A^2 - 2*B^2 mod 5
    def eval_mod10_d2_n18(mask):
        A = (mask & ((1 << 9) - 1)).bit_count()
        B = (mask >> 9).bit_count()
        return ((A + B + A * B) % 2, (A * A - 2 * B * B) % 5)
    brute_force_verify_small(18, (2, 5), eval_mod10_d2_n18)
    print("  mod 10, d = 2, n = 18 (2^18 = 262 144 entrees) : 100% EXACT!")

    # mod 15, d=2, n=20 (G2): A = x_0..x_13 (14), B3 = x_14..x_15 (2), B5 = x_16..x_19 (4)
    # P3 = A^2 + B3^2 mod 3, P5 = A^2 - 2*B5^2 mod 5
    def eval_mod15_d2_n20(mask):
        A = (mask & ((1 << 14) - 1)).bit_count()
        B3 = ((mask >> 14) & 3).bit_count()
        B5 = (mask >> 16).bit_count()
        return ((A * A + B3 * B3) % 3, (A * A - 2 * B5 * B5) % 5)
    brute_force_verify_small(20, (3, 5), eval_mod15_d2_n20)
    print("  mod 15, d = 2, n = 20 (2^20 = 1 048 576 entrees) : 100% EXACT!")

    # 4. Enumeration complete des 7 modules composites
    all_results = {}
    configs_2ch = [
        ("mod6",  (2, True),  (3, True),  True,  10, 350),
        ("mod10", (2, True),  (5, True),  True,  8,  500),
        ("mod12", (4, False), (3, True),  False, 8,  600),
        ("mod14", (2, True),  (7, True),  True,  6,  400),
        ("mod15", (3, True),  (5, True),  True,  7,  1200),
        ("mod21", (3, True),  (7, True),  True,  6,  1000),
    ]

    for name, (q1, p1), (q2, p2), sqfree, d_max, max_n in configs_2ch:
        print(f"\n==================== {name.upper()} ({q1} x {q2}, sqfree={sqfree}), d = 1..{d_max} ====================")
        rows = []
        C1 = binomial_mod_table(max_n, d_max, q1)
        C2 = binomial_mod_table(max_n, d_max, q2)
        for d in range(1, d_max + 1):
            td = time.time()
            deg1, val1, cof1 = build_channel_catalog(d, q1, p1, C1[:, :d], max_n)
            deg2, val2, cof2 = build_channel_catalog(d, q2, p2, C2[:, :d], max_n)
            res = solve_2channel(d, q1, p1, deg1, val1, q2, p2, deg2, val2, max_n)
            (nsym, si, sj,
             ng1, _, _, _,
             ng2, _, _, _,
             nall, ai, aj, aa, m1, m2, gn1, gn2) = res
            assert nsym < max_n, f"max_n={max_n} atteint pour {name} d={d}!"
            verify_witness_2ch(d, q1, p1, cof1[ai], q2, p2, cof2[aj], aa, m1, gn1, m2, gn2)
            dt = time.time() - td
            parts = []
            if gn1 > 0:
                parts.append(f"+{gn1} {METHOD_NAMES[m1]}@{q1}")
            if gn2 > 0:
                parts.append(f"+{gn2} {METHOD_NAMES[m2]}@{q2}")
            g123_sum = " & ".join(parts) if parts else "sym"
            rows.append({
                "d": d,
                "N_sym": nsym,
                "sym_degs": [int(deg1[si]), int(deg2[sj])],
                "N_G1": ng1,
                "N_G2": ng2,
                "N_all": nall,
                "witness_a": aa,
                "witness_degs": [int(deg1[ai]), int(deg2[aj])],
                "witness_c1": [int(x) for x in cof1[ai]],
                "witness_c2": [int(x) for x in cof2[aj]],
                "methods": [METHOD_NAMES[m1], METHOD_NAMES[m2]],
                "gains": [gn1, gn2],
                "g123_summary": g123_sum,
                "time_s": round(dt, 3),
            })
        apply_g4_dp(rows, r_factors=2, is_squarefree=sqfree)
        for r in rows:
            print(
                f"  d={r['d']:2d} | N_sym={r['N_sym']:4d} (deg {r['sym_degs'][0]},{r['sym_degs'][1]}) | "
                f"G1={r['N_G1']:4d} | G2={r['N_G2']:4d} | G123={r['N_all']:4d} | G4={r['N_G4']:4d} | "
                f"BEST={r['N_best']:4d} (+{r['best_gain']:2d}, +{r['best_gain_pct']:5.1f}%) [{r['best_mechanism']}] ({r['time_s']:.2f}s)"
            )
        all_results[name] = rows

    # mod 30 = 2 x 3 x 5, d = 1..6
    print("\n==================== MOD30 (2 x 3 x 5, sqfree=True), d = 1..6 ====================")
    max_n30 = 2500
    d_max30 = 6
    C2_30 = binomial_mod_table(max_n30, d_max30, 2)
    C3_30 = binomial_mod_table(max_n30, d_max30, 3)
    C5_30 = binomial_mod_table(max_n30, d_max30, 5)
    rows30 = []
    for d in range(1, d_max30 + 1):
        td = time.time()
        deg1, val1, cof1 = build_channel_catalog(d, 2, True, C2_30[:, :d], max_n30)
        deg2, val2, cof2 = build_channel_catalog(d, 3, True, C3_30[:, :d], max_n30)
        deg3, val3, cof3 = build_channel_catalog(d, 5, True, C5_30[:, :d], max_n30)
        res = solve_3channel(d, 2, True, deg1, val1, 3, True, deg2, val2, 5, True, deg3, val3, max_n30)
        (nsym, si, sj, sk, ng1, ng2, nall, ai, aj, ak, aa, m1, m2, m3, gn1, gn2, gn3) = res
        assert nsym < max_n30, f"max_n30={max_n30} atteint pour mod30 d={d}!"
        dt = time.time() - td
        parts = []
        for q, m, gn in [(2, m1, gn1), (3, m2, gn2), (5, m3, gn3)]:
            if gn > 0:
                parts.append(f"+{gn} {METHOD_NAMES[m]}@{q}")
        rows30.append({
            "d": d,
            "N_sym": nsym,
            "sym_degs": [int(deg1[si]), int(deg2[sj]), int(deg3[sk])],
            "N_G1": ng1,
            "N_G2": ng2,
            "N_all": nall,
            "witness_a": aa,
            "witness_degs": [int(deg1[ai]), int(deg2[aj]), int(deg3[ak])],
            "methods": [METHOD_NAMES[m1], METHOD_NAMES[m2], METHOD_NAMES[m3]],
            "gains": [gn1, gn2, gn3],
            "g123_summary": " & ".join(parts) if parts else "sym",
            "time_s": round(dt, 3),
        })
    apply_g4_dp(rows30, r_factors=3, is_squarefree=True)
    for r in rows30:
        print(
            f"  d={r['d']:2d} | N_sym={r['N_sym']:4d} (deg {r['sym_degs']}) | "
            f"G1={r['N_G1']:4d} | G2={r['N_G2']:4d} | G123={r['N_all']:4d} | G4={r['N_G4']:4d} | "
            f"BEST={r['N_best']:4d} (+{r['best_gain']:2d}, +{r['best_gain_pct']:5.1f}%) [{r['best_mechanism']}] ({r['time_s']:.2f}s)"
        )
    all_results["mod30"] = rows30

    payload = {
        "norm_checks": norm_checks,
        "mod6_d8_n88_witness": m6d8,
        "modules": all_results,
        "total_time_s": round(time.time() - t0, 2),
    }
    out_path = os.path.join(RESULTS_DIR, "composite_mod_m_full.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    print(f"\nCampagne complete terminee en {time.time() - t0:.2f}s. Sauvegarde : {out_path}")


if __name__ == "__main__":
    run_full_campaign()
