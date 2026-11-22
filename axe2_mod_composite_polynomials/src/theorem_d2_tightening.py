"""
Axe II — Resserrer le Theoreme D1 (Theoreme D2) :
1. Verification exhaustive de l'indice de Witt (sous-espace totalement singulier V_d)
   et du poids minimax exact alpha(s) pour les polynomes quadratiques sur F_2^s (s = 3, 4, 5, 6).
2. Calcul exact de la pire suite de poids w*(j) à l'etape j pour un sous-espace F de K_U
   de dimension s = 3, 4, 5, 6, et optimisation dynamique du choix de s(j) à chaque etape j.
3. Combinaison avec :
   - nu_2(2) = 7 exact (et nu_2^{THR}(2) = 5 exact pour a = 0 / |a|_H < n)
   - Lemme D2a : N_gen^{(2p)}((1, 2), OR) = 4p - 3 exact
   - Lemme D2b : N_gen^{(2p)}((2, 1), OR) <= 3p - 1 exact
"""

import json
import math
import time
from itertools import combinations
from pathlib import Path


def analyze_quadratic_forms_f2(s: int):
    r"""
    Pour chaque polynome quadratique Q sur F_2^s avec Q(0) = 0 :
    - Calcule Z(Q) = {y in F_2^s \ {0} : Q(y) = 0}
    - Trouve la dimension maximale d d'un sous-espace vectoriel V_d <= Z(Q) U {0}
    - Calcule le poids moyen / minimax des vecteurs dans V_d \ {0} ou dans un sous-ensemble
      de Z(Q) de somme nulle.
    """
    n_quad = s + s * (s - 1) // 2
    pairs = [(i, j) for i in range(s) for j in range(i + 1, s)]

    # Pre-calculer les masques de chaque monome sur x in 1..(2^s - 1)
    N = 1 << s
    mono_masks = []
    for i in range(s):
        mask = 0
        for x in range(1, N):
            if (x >> i) & 1:
                mask |= (1 << x)
        mono_masks.append(mask)
    for (i, j) in pairs:
        mask = 0
        for x in range(1, N):
            if ((x >> i) & 1) and ((x >> j) & 1):
                mask |= (1 << x)
        mono_masks.append(mask)

    # Pre-calculer tous les sous-espaces vectoriels de dimension d = 1, 2, 3 de F_2^s
    # Un sous-espace de dim 2 est {0, a, b, a^b} (masque (1<<a)|(1<<b)|(1<<(a^b)))
    subspaces_d2 = []
    subspaces_d2_free_wt_sum = []
    seen_d2 = set()
    for a in range(1, N):
        for b in range(a + 1, N):
            c = a ^ b
            key = tuple(sorted((a, b, c)))
            if key not in seen_d2:
                seen_d2.add(key)
                mask = (1 << a) | (1 << b) | (1 << c)
                wt_sum = a.bit_count() + b.bit_count() + c.bit_count()
                subspaces_d2.append(mask)
                subspaces_d2_free_wt_sum.append(wt_sum)

    subspaces_d3 = []
    subspaces_d3_free_wt_sum = []
    if s >= 5:
        seen_d3 = set()
        for a in range(1, N):
            for b in range(a + 1, N):
                ab = a ^ b
                for c in range(b + 1, N):
                    if c == ab:
                        continue
                    pts = (a, b, c, ab, a ^ c, b ^ c, ab ^ c)
                    key = tuple(sorted(pts))
                    if key not in seen_d3:
                        seen_d3.add(key)
                        mask = 0
                        wt_sum = 0
                        for p in pts:
                            mask |= (1 << p)
                            wt_sum += p.bit_count()
                        subspaces_d3.append(mask)
                        subspaces_d3_free_wt_sum.append(wt_sum)

    full_mask = (1 << N) - 2  # bits 1..N-1

    min_zeros = N
    count_has_d2 = 0
    count_has_d3 = 0
    max_min_d2_free_wt_sum = 0
    max_min_d3_free_wt_sum = 0
    no_d2_zero_counts = set()

    # Pour chaque Q, calculer val_mask (1 là où Q(x) = 1), donc zero_mask = full_mask ^ val_mask
    # Utiliser un parcours en code de Gray pour aller très vite !
    val_mask = 0
    total_polys = 1 << n_quad

    def check_zero_mask(zmask):
        nonlocal min_zeros, count_has_d2, count_has_d3
        nonlocal max_min_d2_free_wt_sum, max_min_d3_free_wt_sum
        nz = zmask.bit_count()
        if nz < min_zeros:
            min_zeros = nz

        best_d2_wt = 999
        for idx, smask in enumerate(subspaces_d2):
            if (zmask & smask) == smask:
                wt = subspaces_d2_free_wt_sum[idx]
                if wt < best_d2_wt:
                    best_d2_wt = wt
        if best_d2_wt < 999:
            count_has_d2 += 1
            if best_d2_wt > max_min_d2_free_wt_sum:
                max_min_d2_free_wt_sum = best_d2_wt
        else:
            no_d2_zero_counts.add(nz)

        if s >= 5:
            best_d3_wt = 999
            for idx, smask in enumerate(subspaces_d3):
                if (zmask & smask) == smask:
                    wt = subspaces_d3_free_wt_sum[idx]
                    if wt < best_d3_wt:
                        best_d3_wt = wt
            if best_d3_wt < 999:
                count_has_d3 += 1
                if best_d3_wt > max_min_d3_free_wt_sum:
                    max_min_d3_free_wt_sum = best_d3_wt

    check_zero_mask(full_mask ^ val_mask)
    for step in range(1, total_polys):
        bit = (step & -step).bit_length() - 1
        val_mask ^= mono_masks[bit]
        check_zero_mask(full_mask ^ val_mask)

    return {
        "s": s,
        "total_polys": total_polys,
        "min_nonzero_zeros": min_zeros,
        "num_subspaces_d2": len(subspaces_d2),
        "count_has_d2": count_has_d2,
        "all_have_d2": (count_has_d2 == total_polys),
        "max_min_d2_free_wt_sum": max_min_d2_free_wt_sum,
        "no_d2_zero_counts": sorted(list(no_d2_zero_counts)),
        "num_subspaces_d3": len(subspaces_d3),
        "count_has_d3": count_has_d3,
        "all_have_d3": (count_has_d3 == total_polys) if s >= 5 else False,
        "max_min_d3_free_wt_sum": max_min_d3_free_wt_sum,
    }


def main():
    t0 = time.time()
    print("=" * 78)
    print("AXE II — THEOREME D2 : ANALYSE EXACTE DES SOUS-ESPACES SINGULIERS SUR F_2^s")
    print("=" * 78)

    results = {}
    for s in [3, 4, 5]:
        t1 = time.time()
        res = analyze_quadratic_forms_f2(s)
        dt = time.time() - t1
        results[str(s)] = res
        print(f"[s = {s}] ({dt:.2f}s) polys={res['total_polys']}, min_zeros={res['min_nonzero_zeros']}")
        print(f"   all_have_d2 = {res['all_have_d2']} ({res['count_has_d2']}/{res['total_polys']}), "
              f"max_min_d2_free_wt_sum = {res['max_min_d2_free_wt_sum']}, "
              f"no_d2_zero_counts = {res['no_d2_zero_counts']}")
        if s >= 5:
            print(f"   all_have_d3 = {res['all_have_d3']} ({res['count_has_d3']}/{res['total_polys']}), "
                  f"max_min_d3_free_wt_sum = {res['max_min_d3_free_wt_sum']}")

    out_path = Path(__file__).resolve().parent.parent / "results" / "theorem_d2_subspaces.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)
    print(f"\nSauvegarde : {out_path} ({time.time() - t0:.2f}s)")


if __name__ == "__main__":
    main()
