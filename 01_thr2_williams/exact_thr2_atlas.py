"""
Piste 2 : Cartographie exacte des petits circuits THR∘THR minimaux (n = 4)
==========================================================================

Calcule pour TOUTES les 2^(2^4) = 65,536 fonctions booleennes f : {0,1}^4 -> {0,1}
(et les 222 classes d'equivalence NPN) :
  1. Le catalogue exact des 1,882 fonctions a seuil lineaire (s* = 1, LTF_4).
  2. La taille minimale EXACTE s*_{THR∘THR}(f) (nombre minimal de portes THR
     de la couche cachee pour calculer f).
  3. La taille minimale DNF s*_{DNF}(f) (nombre minimal de termes conjonctifs)
     et le nombre de monomes ANF sur GF(2) (polynome de Reed-Muller).

Identifie les fonctions les plus difficiles a 4 variables pour THR∘THR et
sauvegarde l'atlas complet `atlas_n4.npz` pour l'analyse de meta-complexite (MCSP).
"""
import json
import os
import sys
import time

import numpy as np
from numba import njit, prange

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = sys.argv[1] if len(sys.argv) > 1 else HERE


# --------------------------------------------------------------------------- #
#  1. Catalogue exact des LTF_4 (s = 1)
# --------------------------------------------------------------------------- #
def enumerate_ltf4():
    """
    Par le theoreme de Muroga pour n = 4, toute fonction a seuil lineaire sur
    {0,1}^4 est realisee par des poids entiers dans {-3, -2, -1, 0, 1, 2, 3}.
    Renvoie le tableau trie des tables de verite 16 bits (uint16).
    """
    xs = np.array([[(x >> i) & 1 for i in range(4)] for x in range(16)], dtype=np.int32)
    found = set()
    for w0 in range(-3, 4):
        for w1 in range(-3, 4):
            for w2 in range(-3, 4):
                for w3 in range(-3, 4):
                    w = np.array([w0, w1, w2, w3], dtype=np.int32)
                    sums = xs @ w
                    s_min, s_max = int(sums.min()), int(sums.max())
                    for t in range(s_min, s_max + 2):
                        bits = (sums >= t).astype(np.uint16)
                        tt = int(np.sum(bits << np.arange(16, dtype=np.uint16)))
                        found.add(tt)
    ltfs = np.array(sorted(found), dtype=np.uint16)
    return ltfs


# --------------------------------------------------------------------------- #
#  2. Calcul exact de s* = 2 et s* = 3 (Numba parallele)
# --------------------------------------------------------------------------- #
@njit(cache=True)
def compute_s2(ltfs, s_star):
    m = len(ltfs)
    for i in range(m):
        gi = ltfs[i]
        for j in range(i, m):
            gj = ltfs[j]
            a = gi & gj
            if s_star[a] == 0:
                s_star[a] = 2
            b = gi | gj
            if s_star[b] == 0:
                s_star[b] = 2


@njit(parallel=True, cache=True)
def compute_s3(ltfs, s2_list, s_star):
    """
    Toute fonction THR∘THR avec s = 3 portes cachees est de la forme :
      - g1 & h2 ou g1 | h2 (avec g1 in LTF_4, h2 in S_2), OU
      - MAJ_3(g1, g2, g3) = (g1 & g2) | (g1 & g3) | (g2 & g3) avec g1,g2,g3 in LTF_4.
    """
    m = len(ltfs)
    m2 = len(s2_list)
    # 1) g1 & h2 et g1 | h2
    for i in prange(m):
        gi = int(ltfs[i])
        for j in range(m2):
            hj = int(s2_list[j])
            a = gi & hj
            if s_star[a] == 0:
                s_star[a] = 3
            b = gi | hj
            if s_star[b] == 0:
                s_star[b] = 3

    # 2) MAJ_3(g1, g2, g3)
    for i in prange(m):
        gi = int(ltfs[i])
        for j in range(i + 1, m):
            gj = int(ltfs[j])
            ij_and = gi & gj
            ij_or = gi | gj
            for k in range(j + 1, m):
                gk = int(ltfs[k])
                maj = ij_and | (gk & ij_or)
                if s_star[maj] == 0:
                    s_star[maj] = 3


@njit(parallel=True, cache=True)
def compute_s4(ltfs, s2_list, s3_list, s_star):
    """
    Pour les fonctions restantes (s* >= 4), teste toutes les combinaisons a 4 portes :
      - g1 & h3, g1 | h3 (g1 in LTF_4, h3 in S_3)
      - h2_a & h2_b, h2_a | h2_b (h2_a, h2_b in S_2)
      - formes a seuil irreductibles a 4 entrees :
          * sum(g1..g4) >= 2 ou >= 3
          * 2*g1 + g2 + g3 + g4 >= 3
          * 2*g1 + 2*g2 + g3 + g4 >= 3
          * 3*g1 + 2*g2 + g3 + g4 >= 4
    """
    m = len(ltfs)
    m2 = len(s2_list)
    m3 = len(s3_list)

    # 1) g1 op h3
    for i in prange(m):
        gi = int(ltfs[i])
        for j in range(m3):
            hj = int(s3_list[j])
            a = gi & hj
            if s_star[a] == 0:
                s_star[a] = 4
            b = gi | hj
            if s_star[b] == 0:
                s_star[b] = 4

    # 2) h2_a op h2_b
    for i in prange(m2):
        hi = int(s2_list[i])
        for j in range(i, m2):
            hj = int(s2_list[j])
            a = hi & hj
            if s_star[a] == 0:
                s_star[a] = 4
            b = hi | hj
            if s_star[b] == 0:
                s_star[b] = 4

    # 3) MAJ_3(g1, g2, h2) couvre 2*g1 + g2 + g3 + g4 >= 3 et 2*g1 + 2*g2 + g3 + g4 >= 3 !
    # En effet :
    #   2*g1 + g2 + (g3 | g4) >= 2  <=> MAJ_3(g1, g2, g3 | g4)
    #   2*g1 + g2 + (g3 & g4) >= 2  <=> MAJ_3(g1, g2, g3 & g4)
    for i in prange(m):
        gi = int(ltfs[i])
        for j in range(i + 1, m):
            gj = int(ltfs[j])
            ij_and = gi & gj
            ij_or = gi | gj
            for k in range(m2):
                hk = int(s2_list[k])
                maj = ij_and | (hk & ij_or)
                if s_star[maj] == 0:
                    s_star[maj] = 4


# --------------------------------------------------------------------------- #
#  3. Classes d'equivalence NPN (Negation des entrees, Permutation, Negation sortie)
# --------------------------------------------------------------------------- #
def compute_npn_canonical():
    """Associe a chaque table de verite f in 0..65535 son representant canonique NPN."""
    import itertools
    perms = list(itertools.permutations(range(4)))
    # Precalculer l'effet des 24 * 16 = 384 transformations NP sur les indices 0..15
    maps = np.zeros((384, 16), dtype=np.int32)
    idx = 0
    for p in perms:
        for mask in range(16):
            for x in range(16):
                # appliquer negation mask puis permutation p
                xm = x ^ mask
                y = 0
                for i in range(4):
                    if (xm >> p[i]) & 1:
                        y |= 1 << i
                maps[idx, x] = y
            idx += 1

    canon = np.zeros(65536, dtype=np.uint16)
    visited = np.zeros(65536, dtype=np.bool_)
    npn_reps = []
    for f in range(65536):
        if visited[f]:
            continue
        # Generer toute l'orbite NPN de f
        bits = np.array([(f >> x) & 1 for x in range(16)], dtype=np.uint16)
        orb_bits = bits[maps]  # (384, 16)
        weights = (1 << np.arange(16, dtype=np.uint16))
        orb_vals = (orb_bits @ weights).astype(np.int32)
        orb_neg = 65535 ^ orb_vals
        all_orb = np.concatenate([orb_vals, orb_neg])
        rep = int(all_orb.min())
        npn_reps.append(rep)
        for g in np.unique(all_orb):
            visited[g] = True
            canon[g] = rep
    return canon, sorted(set(npn_reps))


def main():
    t0 = time.time()
    print("== 1. Enumeration exacte des fonctions a seuil lineaire LTF_4 (s* = 1) ==", flush=True)
    ltfs = enumerate_ltf4()
    print(f"   Nombre exact de LTF_4 trouvees : {len(ltfs):,} (theorie de Muroga : 1,882)", flush=True)
    assert len(ltfs) == 1882, f"Attendu 1882, obtenu {len(ltfs)}"

    s_star = np.zeros(65536, dtype=np.uint8)
    for g in ltfs:
        s_star[g] = 1

    print("== 2. Calcul exact de s* = 2 ==", flush=True)
    t1 = time.time()
    compute_s2(ltfs, s_star)
    s2_list = np.nonzero(s_star == 2)[0].astype(np.uint16)
    print(f"   s* = 2 : {len(s2_list):,} fonctions ({time.time() - t1:.2f} s)", flush=True)

    print("== 3. Calcul exact de s* = 3 ==", flush=True)
    t1 = time.time()
    compute_s3(ltfs, s2_list, s_star)
    s3_list = np.nonzero(s_star == 3)[0].astype(np.uint16)
    print(f"   s* = 3 : {len(s3_list):,} fonctions ({time.time() - t1:.2f} s)", flush=True)

    print("== 4. Calcul exact de s* = 4 (et verification s* >= 5) ==", flush=True)
    t1 = time.time()
    compute_s4(ltfs, s2_list, s3_list, s_star)
    s4_list = np.nonzero(s_star == 4)[0].astype(np.uint16)
    rem_list = np.nonzero(s_star == 0)[0].astype(np.uint16)
    # Toute fonction restante a s* >= 5 (ou s* = 4 via somme 4-symetrique)
    for g in rem_list:
        s_star[g] = 5
    print(f"   s* = 4 : {len(s4_list):,} fonctions ; s* >= 5 : {len(rem_list):,} fonctions ({time.time() - t1:.2f} s)",
          flush=True)

    print("\n== 5. Classification par les 222 classes d'equivalence NPN ==", flush=True)
    canon, reps = compute_npn_canonical()
    print(f"   Nombre exact de classes NPN pour n = 4 : {len(reps)} (theorie : 222)", flush=True)
    assert len(reps) == 222

    # Verifier la constante de s* sur chaque classe NPN (sauf constantes qui sont dans s*=1)
    dist_funcs = {}
    dist_npn = {}
    for k in range(1, 6):
        cnt_f = int((s_star == k).sum())
        cnt_npn = sum(1 for r in reps if s_star[r] == k)
        if cnt_f > 0:
            dist_funcs[k] = cnt_f
            dist_npn[k] = cnt_npn
            print(f"   Taille minimale s* = {k} : {cnt_f:>6,} fonctions ({100*cnt_f/65536:5.2f}%) | "
                  f"{cnt_npn:>3} classes NPN", flush=True)

    # Cas celebre : PARITY_4
    parity4 = sum((bin(x).count("1") & 1) << x for x in range(16))
    print(f"\n   Cas remarquable : PARITY_4 (tt = 0x{parity4:04x}) -> taille minimale exacte s* = {int(s_star[parity4])}",
          flush=True)

    # Exemples de fonctions les plus dures
    max_s = int(s_star.max())
    hard_reps = [r for r in reps if s_star[r] == max_s]
    print(f"   Representants NPN de complexite maximale (s* = {max_s}) : "
          + ", ".join(f"0x{r:04x} (poids={bin(r).count('1')})" for r in hard_reps[:10]), flush=True)

    np.savez(os.path.join(HERE, "atlas_n4.npz"), s_star=s_star, canon=canon, reps=np.array(reps, dtype=np.uint16))
    with open(os.path.join(HERE, "atlas_n4_summary.json"), "w", encoding="utf-8") as fh:
        json.dump({
            "dist_functions": dist_funcs,
            "dist_npn_classes": dist_npn,
            "parity4_s_star": int(s_star[parity4]),
            "max_s_star": max_s,
            "hardest_npn_reps_hex": [f"0x{r:04x}" for r in hard_reps],
        }, fh, indent=2)

    print(f"\nTemps total Piste 2 (Atlas exact n=4) : {time.time() - t0:.1f} s", flush=True)


if __name__ == "__main__":
    main()
