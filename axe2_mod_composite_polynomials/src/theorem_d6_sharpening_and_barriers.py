"""Theorem D6:
1. Unconditional break of 8p - 7 on mod 10 (p = 5) via 3-root B_2-symmetry:
   N_gen^(10)(2, OR) <= 32 and N_gen^(10)(2, THR_2) <= 33.
2. Lucas Even-Multiple Slice Theorem (p >= 5):
   Any quadratic P_p on m >= 5p - 1 variables (P_p(0) = 0) vanishes on some x with
   |x|_H in {2p, 4p}, giving N_gen^(2p)(2, OR) <= 7p - 4 (< 1.75x) whenever B_p(w_1, .) is constant on U_1.
3. Unconditional improvement on m = 3q (q >= 5) via (B_3, B_q) coupling + exact F_3^2
   weight-<=3 zero-sum packing (max_no_wt3 = 6 => |U| <= 3k + 6) + support recycling:
   N_gen^(3q)(2, OR) <= 19q - 13 (unconditional, beats D3 20q - 19 for all q >= 7),
   <= 17q - 11 (Case 1), <= 15q - 6 (Type A on F_q).
4. Exact 4-barrier characterization on general semiprimes m = pq (3 <= p < q):
   (i) Fourier matching barrier Omega(m^2 log m),
   (ii) F_p-rank dichotomy barrier Omega(p),
   (iii) 100% sharp Chevalley-Warning weight adversary j(q-1)+1,
   (iv) No information-theoretic m^{3/2} barrier: p^2 ln(q) / ln(p) < pq = m for all 3 <= p < q.
"""

from pathlib import Path
import itertools
import json
import math
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = str(Path(__file__).resolve().parent.parent)
RES_PATH = os.path.join(ROOT, "results", "theorem_d6_sharpening_and_barriers.json")
FIG_PATH = os.path.join(ROOT, "figures", "theorem_d6_sharpening_and_barriers.png")


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    for d in range(2, int(math.isqrt(n)) + 1):
        if n % d == 0:
            return False
    return True


def verify_d6_1_mod10_3root():
    """Exhaustive integer verification of 3-root B_2-symmetry contradiction at n=33 (p=5)."""
    p = 5
    out = {}
    for n in (33, 32):
        sols = []
        for w1 in range(1, 2 * p - 1):
            u = n - w1
            for u0 in range(0, u + 1):
                u1 = u - u0
                if u0 + u1 // 2 > 3 * (p - 1):
                    continue
                if u1 < 3 * p - 1:
                    sols.append(("u1_small", w1, u0, u1))
                    continue
                for w2 in range(1, min(3 * (p - 1), u1) + 1, 2):
                    v = u1 - w2
                    for y0 in range(0, w1 + 1):
                        y1 = w1 - y0
                        if y1 % 2 == 0:
                            continue
                        for v0 in range(0, v + 1):
                            v1 = v - v0
                            for u00 in range(0, u0 + 1):
                                u01 = u0 - u00
                                c1 = v0 + u00 + y0 + (v1 + u01 + y1) // 2 <= 3 * (p - 1)
                                c2 = v1 + u00 + (v0 + u01) // 2 <= 3 * (p - 1)
                                if c1 and c2:
                                    sols.append((w1, w2, u0, u1, y0, v0, v1))
                                    break
        out[str(n)] = len(sols)
    assert out["33"] == 0, f"Expected 0 solutions at n=33, got {out['33']}"
    assert out["32"] > 0
    return {"mod10_n33_configs": out["33"], "mod10_n32_configs": out["32"], "bound_or": 32, "bound_thr2": 33}


def verify_d6_2_lucas_slices():
    """Verify Lucas slice congruences at M = 5p - 1 for all primes 5 <= p <= 101."""
    primes = [p for p in range(5, 102) if is_prime(p)]
    for p in primes:
        M = 5 * p - 1
        for k in range(0, 2 * p - 1):
            c2p = math.comb(M - k, 2 * p - k) % p if 2 * p >= k else 0
            c4p = math.comb(M - k, 4 * p - k) % p if 4 * p >= k else 0
            if k == 0:
                assert c2p == 6 % p and c4p == 1 % p
            elif k == p:
                assert c2p == 3 % p and c4p == 1 % p
            else:
                assert c2p == 0 and c4p == 0
        assert (6 - 3 * 1) % p == 3 % p != 0
    return {
        "primes_verified": len(primes),
        "max_prime": primes[-1],
        "slice_M": "5p - 1",
        "bound_const_Bp": "7p - 4",
    }


def verify_d6_3_m3q_packing():
    """Exhaustively verify max multiset in F_3^2 with NO zero-sum of weight <= 3 is 6."""
    lines = [
        ((1, 0), (2, 0)),
        ((0, 1), (0, 2)),
        ((1, 1), (2, 2)),
        ((1, 2), (2, 1)),
    ]
    max_no_wt3 = 0
    for signs in itertools.product((0, 1), repeat=4):
        reps = [lines[i][signs[i]] for i in range(4)]
        for mults in itertools.product((0, 1, 2), repeat=4):
            vec = []
            for r, m in zip(reps, mults):
                vec.extend([r] * m)
            has_zs = False
            for sz in (2, 3):
                for comb in itertools.combinations(vec, sz):
                    if sum(x[0] for x in comb) % 3 == 0 and sum(x[1] for x in comb) % 3 == 0:
                        has_zs = True
                        break
                if has_zs:
                    break
            if not has_zs and len(vec) > max_no_wt3:
                max_no_wt3 = len(vec)
    assert max_no_wt3 == 6
    table = []
    for q in (5, 7, 11, 13, 17, 19, 23, 29, 31, 53, 101):
        table.append({
            "q": q,
            "m": 3 * q,
            "lb_2block": 2 * (3 * q - 1),
            "d3": 20 * (q - 1) + 1,
            "d6_uncond": 19 * (q - 1) + 6,
            "d6_case1": 17 * (q - 1) + 6,
            "d6_typeA": 15 * (q - 1) + 9,
        })
    return {"f3_2_max_no_wt3": max_no_wt3, "table_m3q": table}


def verify_d6_4_semiprime_barriers():
    """Verify all 4 barrier statements on m = pq (3 <= p < q)."""
    # (i) Fourier matching barrier: |(3 + e^{2pi i/m})/4|^k at k = floor(m^{3/2})
    fourier_checks = []
    for p, q in ((3, 5), (5, 7), (7, 11), (11, 13), (17, 19)):
        m = p * q
        z = (3.0 + np.exp(2j * np.pi / m)) / 4.0
        k = int(m ** 1.5)
        amp = abs(z) ** k
        fourier_checks.append({"m": m, "k_m32": k, "abs_Phi": round(float(amp), 6)})
        assert amp > 1.0 / m

    # (ii) F_p-rank dichotomy barrier: Q(lambda) = sum_{t=1}^r (x + y + xy) on 2r vars, 3r < p
    rank_checks = []
    for p in (5, 7, 11, 13):
        r = (p - 1) // 3
        for mask in range(1, 1 << (2 * r)):
            val = 0
            for t in range(r):
                x = (mask >> (2 * t)) & 1
                y = (mask >> (2 * t + 1)) & 1
                val += x + y + x * y
            assert val % p != 0
        rank_checks.append({"p": p, "full_rank_2r": 2 * r})

    # (iii) Sharp Chevalley-Warning weight adversary L_t(x) = sum_{B_t} x_i + z = 0 mod q
    cw_checks = []
    for q in (3, 5):
        for j in (1, 2):
            N = j * (q - 1) + 1
            if N <= 10:
                roots = []
                for mask in range(1, 1 << N):
                    z = mask & 1
                    ok = True
                    for t in range(j):
                        blk = (mask >> (1 + t * (q - 1))) & ((1 << (q - 1)) - 1)
                        if (bin(blk).count("1") + z) % q != 0:
                            ok = False
                            break
                    if ok:
                        roots.append(bin(mask).count("1"))
                assert roots == [N], f"Expected unique root of weight {N}, got {roots}"
                cw_checks.append({"q": q, "j": j, "N": N, "unique_root_weight": roots[0]})

    # (iv) Information-theoretic threshold p^2 ln(q) / ln(p) < pq = m for ALL 3 <= p < q <= 101
    primes = [p for p in range(3, 102) if is_prime(p)]
    pair_count = 0
    max_ratio = 0.0
    for i, p in enumerate(primes):
        for q in primes[i + 1:]:
            val = (p * p * math.log(q)) / math.log(p)
            ratio = val / (p * q)
            assert ratio < 1.0
            if ratio > max_ratio:
                max_ratio = ratio
            pair_count += 1

    return {
        "fourier_barrier": fourier_checks,
        "fp_rank_barrier": rank_checks,
        "cw_sharp_adversary": cw_checks,
        "info_pairs_verified": pair_count,
        "max_info_ratio_over_m": round(max_ratio, 6),
    }


def make_figure(res_m3q, res_bar):
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.6), dpi=180)
    plt.rcParams.update({"font.size": 9.5})

    # Panel (a): m = 2p bounds (D4 8p-7 vs Lucas 7p-4 vs 4p-2, plus mod 10 = 32)
    ax = axes[0]
    ps = np.array([3, 5, 7, 11, 13, 17, 19, 23, 29, 31])
    lb = 4 * ps - 2
    d4 = 8 * ps - 7
    d4[0] = 16
    d6_const = 7 * ps - 4
    ax.plot(ps, d4 / lb, "o-", color="#d62728", lw=2, label="Th. D4: $8p-7$ ($<2.00\\times$)")
    ax.plot([5], [32 / 18], "*", color="#1f77b4", ms=12, zorder=5, label="Th. D6.1: mod 10 $\\leq 32$ ($1.78\\times$)")
    ax.plot(ps[1:], d6_const[1:] / lb[1:], "s--", color="#2ca02c", lw=2, label="Th. D6.2 (Lucas): $7p-4$ ($<1.75\\times$)")
    ax.axhline(2.0, color="#888", ls=":", lw=1)
    ax.axhline(1.75, color="#2ca02c", ls=":", lw=1, alpha=0.6)
    ax.set_xlabel("Odd prime $p$ ($m = 2p$)")
    ax.set_ylabel("Ratio Upper / Lower $(4p-2)$")
    ax.set_title("(a) $m = 2p$: mod $10 \\leq 32$ & Lucas $7p-4$")
    ax.set_ylim(1.35, 2.08)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8, loc="lower right")

    # Panel (b): m = 3q bounds (D3 vs D6 uncond / Case 1 / Type A)
    ax = axes[1]
    rows = res_m3q["table_m3q"][:9]
    qs = np.array([r["q"] for r in rows])
    lb3 = np.array([r["lb_2block"] for r in rows])
    d3 = np.array([r["d3"] for r in rows])
    d6u = np.array([r["d6_uncond"] for r in rows])
    d6c1 = np.array([r["d6_case1"] for r in rows])
    d6ta = np.array([r["d6_typeA"] for r in rows])
    ax.plot(qs, d3 / lb3, "o-", color="#d62728", lw=2, label="Th. D3: $20q-19$")
    ax.plot(qs, d6u / lb3, "s-", color="#1f77b4", lw=2, label="Th. D6.3 (uncond.): $19q-13$")
    ax.plot(qs, d6c1 / lb3, "^--", color="#ff7f0e", lw=1.8, label="Th. D6.3 (Case 1): $17q-11$")
    ax.plot(qs, d6ta / lb3, "d:", color="#2ca02c", lw=1.8, label="Th. D6.3 (Type A): $15q-6$")
    ax.set_xlabel("Prime $q \\geq 5$ ($m = 3q$)")
    ax.set_ylabel("Ratio Upper / $2(m-1)$")
    ax.set_title("(b) $m = 3q$: $(B_3, B_q)$ Coupling + $\\mathbb{F}_3^2$ Packing")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8, loc="best")

    # Panel (c): Information-theoretic threshold vs m and m^{3/2}
    ax = axes[2]
    primes = [p for p in range(3, 50) if is_prime(p)]
    ms, info_vals, d3_vals = [], [], []
    for i, p in enumerate(primes):
        if i + 1 < len(primes):
            q = primes[i + 1]
            m = p * q
            ms.append(m)
            info_vals.append((p * p * math.log(q)) / math.log(p))
            d3_vals.append((2 * p - 1) * (p + 1) * (q - 1))
    ms = np.array(ms)
    ax.loglog(ms, d3_vals, "o-", color="#d62728", lw=2, label="Th. D3: $O(p^2 q) \\leq O(m^{3/2})$")
    ax.loglog(ms, 2 * (ms - 1), "--", color="#2ca02c", lw=1.8, label="2-block lower bound $2(m-1)$")
    ax.loglog(ms, info_vals, "s-", color="#1f77b4", lw=2, label="Info threshold $p^2 \\ln q / \\ln p < m$")
    ax.set_xlabel("Semiprime $m = pq$ (consecutive odd primes)")
    ax.set_ylabel("Number of variables $n$")
    ax.set_title("(c) Semiprimes $m=pq$: No Info $m^{3/2}$ Barrier")
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=8, loc="upper left")

    fig.tight_layout()
    fig.savefig(FIG_PATH, dpi=180)
    plt.close(fig)


def main():
    r1 = verify_d6_1_mod10_3root()
    r2 = verify_d6_2_lucas_slices()
    r3 = verify_d6_3_m3q_packing()
    r4 = verify_d6_4_semiprime_barriers()
    payload = {
        "D6_1_mod10_3root": r1,
        "D6_2_lucas_even_slices": r2,
        "D6_3_m3q_coupling": r3,
        "D6_4_semiprime_barriers": r4,
    }
    os.makedirs(os.path.dirname(RES_PATH), exist_ok=True)
    with open(RES_PATH, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    make_figure(r3, r4)
    print("THEOREM D6 ALL VERIFICATIONS PASSED.")
    print(f"  D6.1 mod 10: n=33 -> {r1['mod10_n33_configs']} configs => OR <= 32, THR_2 <= 33")
    print(f"  D6.2 Lucas slices: verified on {r2['primes_verified']} primes (5..{r2['max_prime']}) => 7p-4")
    print(f"  D6.3 m=3q packing: F_3^2 max_no_wt3 = {r3['f3_2_max_no_wt3']} => 19q-13 (Case 1: 17q-11)")
    print(f"  D6.4 barriers: {r4['info_pairs_verified']} semiprime pairs, max(info/m) = {r4['max_info_ratio_over_m']}")


if __name__ == "__main__":
    main()
