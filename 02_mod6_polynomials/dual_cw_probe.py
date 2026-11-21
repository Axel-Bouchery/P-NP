r"""
Sonde exacte pour le Theoreme D3 (borne LINEAIRE en p pour N_gen^{(2p)}(d, OR)) :
1. Sur F_3 (p = 3), test exact / recherche sur m = 5, 6 :
   Quel est le poids minimal maximal mu_1(3, m) d'un zero non nul w_1 in {0,1}^m \ {0}
   d'un polynome quadratique P_3 (P_3(0) = 0) ?
2. Pour m = 5 sur F_3 : existe-t-il un P_3 sans zero de poids 1, 2, 3 ? (donc min_wt >= 4)
   Ou sans zero de poids 1, 2, 3, 4 ? (donc min_wt == 5)
3. Verification de la table comparative D1 vs D2 vs D3 (9p - 7) et modules pq (15, 21, 35).
"""

import time
import random
from itertools import combinations


def probe_p3_min_zero_weight_m5():
    """
    Sur m = 5 variables sur F_3, il y a 5 termes lineaires (c_i in {1, 2} s'il n'y a aucun zero de poids 1)
    et 10 termes quadratiques d_{ij} in {0, 1, 2}.
    A symetrie S_5 pres, le vecteur c est determine par k = #{i : c_i = 2} in {0, 1, 2}
    (car on peut multiplier P_3 par -1 mod 3, donc k <= 2).
    Pour chaque k in {0, 1, 2}, il n'y a que 3^{10} = 59 049 polynomes quadratiques !
    C'est 100% EXHAUSTIF en moins de 1 seconde !
    """
    m = 5
    N = 1 << m
    pairs = [(i, j) for i in range(m) for j in range(i + 1, m)]
    # Pre-calculer les vecteurs de poids 1, 2, 3, 4, 5
    vecs_by_wt = {w: [x for x in range(1, N) if x.bit_count() == w] for w in range(1, m + 1)}
    wt2 = vecs_by_wt[2]
    wt3 = vecs_by_wt[3]
    wt4 = vecs_by_wt[4]
    wt5 = vecs_by_wt[5]

    # Pre-calculer les paires dans chaque x
    pair_idx = {p: idx for idx, p in enumerate(pairs)}
    x_pairs = {}
    for x in range(1, N):
        bits = [i for i in range(m) if (x >> i) & 1]
        x_pairs[x] = [pair_idx[(bits[a], bits[b])] for a in range(len(bits)) for b in range(a + 1, len(bits))]

    # Pour chaque k in {0, 1, 2} (k variables avec c_i = 2, 5-k avec c_i = 1)
    # Une paire (i, j) n'est PAS un zero de poids 2 ssi c_i + c_j + d_{ij} != 0 (mod 3),
    # donc d_{ij} ne peut prendre que 2 valeurs (au lieu de 3) !
    # Cela reduit l'espace de 3^{10} = 59 049 a seulement 2^{10} = 1 024 polynomes par k !!
    results = {}
    for k in range(3):
        c = [2 if i < k else 1 for i in range(m)]
        allowed_d = []
        for (i, j) in pairs:
            forb = (-c[i] - c[j]) % 3
            allowed_d.append([v for v in range(3) if v != forb])

        c_sum = {x: sum(c[i] for i in range(m) if (x >> i) & 1) % 3 for x in range(1, N)}

        max_min_wt = 2
        count_ge3 = 0
        count_ge4 = 0
        count_ge5 = 0
        example_ge4 = None
        example_ge5 = None

        for mask in range(1 << 10):
            d = [allowed_d[idx][(mask >> idx) & 1] for idx in range(10)]
            count_ge3 += 1

            # Tester poids 3
            has_wt3 = False
            for x in wt3:
                val = c_sum[x]
                for pidx in x_pairs[x]:
                    val += d[pidx]
                if val % 3 == 0:
                    has_wt3 = True
                    break
            if has_wt3:
                if max_min_wt < 3:
                    max_min_wt = 3
                continue

            count_ge4 += 1
            if example_ge4 is None:
                example_ge4 = (c[:], d[:])

            # Tester poids 4
            has_wt4 = False
            for x in wt4:
                val = c_sum[x]
                for pidx in x_pairs[x]:
                    val += d[pidx]
                if val % 3 == 0:
                    has_wt4 = True
                    break
            if has_wt4:
                if max_min_wt < 4:
                    max_min_wt = 4
                continue

            # Verifier poids 5
            x5 = wt5[0]
            val5 = (c_sum[x5] + sum(d)) % 3
            assert val5 == 0, f"Contre-exemple Chevalley-Warning?! c={c}, d={d}"
            count_ge5 += 1
            max_min_wt = 5
            if example_ge5 is None:
                example_ge5 = (c[:], d[:])

        print(f"[p=3, m=5, k={k}] ge3={count_ge3}, ge4={count_ge4}, ge5={count_ge5}, max_min_wt={max_min_wt}")
        results[k] = (count_ge4, count_ge5, max_min_wt, example_ge4, example_ge5)

    return results


def probe_p3_min_zero_weight_m6():
    """
    Sur m = 6 variables sur F_3 : existe-t-il un polynome quadratique P_3 sans zero de poids <= 3 ?
    Pour un tel P_3, chaque sous-ensemble de 5 variables (les 6 sous-graphes a 5 sommets)
    doit etre dans la liste des polynomes ge4 de m = 5 !
    """
    m = 5
    N = 1 << m
    pairs5 = [(i, j) for i in range(m) for j in range(i + 1, m)]
    pair5_idx = {p: idx for idx, p in enumerate(pairs5)}
    wt3 = [x for x in range(1, N) if x.bit_count() == 3]
    x_pairs = {x: [pair5_idx[(bits[a], bits[b])] for bits in [[i for i in range(m) if (x >> i) & 1]]
                   for a in range(len(bits)) for b in range(a + 1, len(bits))] for x in wt3}

    # Collecter tous les polynomes (c, d) sur m = 5 sans zero de poids 1, 2, 3
    bad5 = set()
    for c_code in range(1 << 5):
        c = [1 + ((c_code >> i) & 1) for i in range(5)]
        allowed_d = [[v for v in range(3) if v != (-c[i] - c[j]) % 3] for (i, j) in pairs5]
        c_sum = {x: sum(c[i] for i in range(5) if (x >> i) & 1) % 3 for x in wt3}
        for mask in range(1 << 10):
            d = tuple(allowed_d[idx][(mask >> idx) & 1] for idx in range(10))
            if not any((c_sum[x] + sum(d[pidx] for pidx in x_pairs[x])) % 3 == 0 for x in wt3):
                bad5.add((tuple(c), d))

    print(f"[p=3, m=5] Total exact (c, d) sans zero de poids <= 3 : {len(bad5)}")

    # Etendre a m = 6 : essayer les 2 choix de c_5 in {1, 2} et les 2^5 = 32 choix autorises pour d_{i, 5}
    pairs6 = [(i, j) for i in range(6) for j in range(i + 1, 6)]
    pair6_idx = {p: idx for idx, p in enumerate(pairs6)}
    bad6_count = 0
    bad6_no_wt4 = 0
    for (c5, d5) in bad5:
        adj = [[0] * 6 for _ in range(6)]
        for idx, (i, j) in enumerate(pairs5):
            adj[i][j] = adj[j][i] = d5[idx]
        for c_new in (1, 2):
            c6 = c5 + (c_new,)
            allowed_ext = [[v for v in range(3) if v != (-c5[i] - c_new) % 3] for i in range(5)]
            for emask in range(1 << 5):
                for i in range(5):
                    v = allowed_ext[i][(emask >> i) & 1]
                    adj[i][5] = adj[5][i] = v
                # Verifier tous les triples contenant le sommet 5
                ok_wt3 = True
                for i in range(5):
                    for j in range(i + 1, 5):
                        if (c6[i] + c6[j] + c6[5] + adj[i][j] + adj[i][5] + adj[j][5]) % 3 == 0:
                            ok_wt3 = False
                            break
                    if not ok_wt3:
                        break
                if ok_wt3:
                    bad6_count += 1
                    # Verifier s'il a un zero de poids 4
                    has_wt4 = False
                    for comb in combinations(range(6), 4):
                        s = sum(c6[u] for u in comb)
                        for a in range(4):
                            for b in range(a + 1, 4):
                                s += adj[comb[a]][comb[b]]
                        if s % 3 == 0:
                            has_wt4 = True
                            break
                    if not has_wt4:
                        bad6_no_wt4 += 1

    print(f"[p=3, m=6] Total exact sans zero de poids <= 3 : {bad6_count}, sans zero de poids <= 4 : {bad6_no_wt4}")


def main():
    t0 = time.time()
    print("=" * 78)
    print("AXE II — SONDE EXACTE SUR LE POIDS DES ZEROS DE P_3 ET BORNE LINEAIRE D3")
    print("=" * 78)
    probe_p3_min_zero_weight_m5()
    probe_p3_min_zero_weight_m6()
    print(f"\nTermine en {time.time() - t0:.2f}s")


if __name__ == "__main__":
    main()
