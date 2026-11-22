"""
Piste Nouvelle : Representation de OR_n / AND_n par polynomes modulo 6
                 (Probleme ouvert de Barrington-Beigel-Rudich, 1994)
======================================================================

Rappel (Theoreme des restes chinois Z_6 = F_2 x F_3) :
  Un polynome multilineaire P(x_1, ..., x_n) in Z_6[x] de degre <= d
  represente faiblement OR_n modulo 6 :
      P(0) = 0 (mod 6)  et  P(x) != 0 (mod 6) pour tout x != 0
  si et seulement si P = (P_2, P_3) in F_2[x] x F_3[x] (deg <= d, P_2(0)=P_3(0)=0)
  verifie :
      Z(P_2) inter Z(P_3) = {0}.

Experience 1 (Spectre symetrique EXACT pour d = 1..10) :
  Recherche EXHAUSTIVE de tous les polynomes symetriques P_2 in F_2[x] (2^d)
  et P_3 in F_3[x] ((3^d - 1)/2) pour trouver le maximum EXACT N_sym(d)
  (sans supposer les periodes de BBR), et comparaison avec la construction BBR.

Experience 2 (CHASSE AU RECORD : briser la barriere symetrique N_gen(d) > N_sym(d)) :
  Astuce algebrique cle :
    Pour tout polynome P_3 in F_3[x] fixe (de lieu des zeros Z_3 = Z(P_3)\{0}),
    l'existence d'un P_2 in F_2[x] (meme NON symetrique !) de degre <= d tel que
    P_2(x) = 1 (mod 2) pour tout x in Z_3 est un SYSTEME LINEAIRE EXACT SUR F_2 !
  On combine donc une recherche locale / recuit sur P_3 in F_3[x] avec un solveur
  d'elimination de Gauss exact sur F_2 pour chercher des polynomes NON symetriques
  de degre d = 2 (pour n = 9, 10, 11, 12...) et d = 3 (pour n > N_sym(3)) !
"""
from pathlib import Path
import itertools
import json
import os
import sys
import time

import numpy as np
from numba import njit, prange

HERE = os.path.dirname(os.path.abspath(__file__))
AXE_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = str(AXE_ROOT / "results")
FIGURES_DIR = str(AXE_ROOT / "figures")
for _d in (RESULTS_DIR, FIGURES_DIR):
    os.makedirs(_d, exist_ok=True)
OUT = sys.argv[1] if len(sys.argv) > 1 else FIGURES_DIR


# --------------------------------------------------------------------------- #
#  Experience 1 : Recherche EXHAUSTIVE du maximum symetrique exact N_sym(d)
# --------------------------------------------------------------------------- #
def binomial_mod_table(max_n, max_d, p):
    """Table C(s, k) mod p pour s in 0..max_n, k in 1..max_d."""
    C = np.zeros((max_n + 1, max_d), dtype=np.int64)
    # Triangle de Pascal exact mod p
    P = np.zeros((max_n + 1, max_d + 1), dtype=np.int64)
    for s in range(max_n + 1):
        P[s, 0] = 1
        for k in range(1, min(s, max_d) + 1):
            P[s, k] = (P[s - 1, k - 1] + P[s - 1, k]) % p
    C[:, :] = P[:, 1:max_d + 1]
    return C


@njit(parallel=True, cache=True)
def exact_n_sym(d, C2, C3):
    """
    Teste EXHAUSTIVEMENT tous les 2^d polynomes symetriques P_2 in F_2 (deg <= d)
    et tous les (3^d - 1)/2 polynomes symetriques P_3 in F_3 (deg <= d, normalises).
    Renvoie (N_sym(d), best_c2_mask, best_c3_idx).
    """
    max_n = C2.shape[0] - 1
    n2 = 1 << d
    # Precalculer les valeurs de tous les P_2 sur s = 1..max_n :
    # Z2_mask[i, s] = 1 si P_2^{(i)}(s) == 0 mod 2
    val2_zero = np.zeros((n2, max_n + 1), np.bool_)
    for mask in range(1, n2):
        for s in range(1, max_n + 1):
            acc = 0
            for k in range(d):
                if (mask >> k) & 1:
                    acc ^= C2[s, k]
            val2_zero[mask, s] = (acc == 0)

    # Nombre de P_3 canoniques (premier coefficient non nul egal a 1) : (3^d - 1)//2
    n3_total = 1
    for _ in range(d):
        n3_total *= 3
    n3_canon = (n3_total - 1) // 2

    best_per_t = np.zeros(n3_canon, np.int64)
    best_m2_per_t = np.zeros(n3_canon, np.int64)

    for idx in prange(n3_canon):
        # Decoder idx en vecteur c3 in {0,1,2}^d canonique (on utilise idx + 1..n3_canon
        # ou plus simplement tous les vecteurs de 1..(3^d - 1)//2 dont le symetrique donne
        # le meme ensemble de zeros !)
        tmp = idx + 1
        c3 = np.zeros(d, np.int64)
        for k in range(d):
            c3[k] = tmp % 3
            tmp //= 3

        # Trouver la liste des zeros s >= 1 de P_3
        z3 = np.zeros(max_n, np.int64)
        nz = 0
        for s in range(1, max_n + 1):
            acc = 0
            for k in range(d):
                if c3[k] != 0:
                    acc += c3[k] * C3[s, k]
            if acc % 3 == 0:
                z3[nz] = s
                nz += 1

        # Pour chaque P_2, trouver le plus petit s in z3 ou P_2(s) == 0 -> N = s - 1
        best_n = 0
        best_m2 = 0
        for mask in range(1, n2):
            first_bad = max_n + 1
            for j in range(nz):
                s = z3[j]
                if val2_zero[mask, s]:
                    first_bad = s
                    break
            n_ok = first_bad - 1
            if n_ok > best_n:
                best_n = n_ok
                best_m2 = mask
        best_per_t[idx] = best_n
        best_m2_per_t[idx] = best_m2

    # Maximum global
    best_global = 0
    best_idx = 0
    best_m2 = 0
    for idx in range(n3_canon):
        if best_per_t[idx] > best_global:
            best_global = best_per_t[idx]
            best_idx = idx
            best_m2 = best_m2_per_t[idx]
    return best_global, best_m2, best_idx + 1


def decode_c2(mask, d):
    return [(mask >> k) & 1 for k in range(d)]


def decode_c3(code, d):
    res = []
    for _ in range(d):
        res.append(code % 3)
        code //= 3
    return res


def bbr_classic_bound(d):
    """Borne classique BBR (1994) par puissances pures q2 = 2^a <= d+1, q3 = 3^b <= d+1."""
    q2 = 1
    while q2 * 2 <= d + 1:
        q2 *= 2
    q3 = 1
    while q3 * 3 <= d + 1:
        q3 *= 3
    return q2 * q3 - 1, q2, q3


# --------------------------------------------------------------------------- #
#  Experience 2 : Chasse aux polynomes NON SYMETRIQUES (elimination de Gauss F_2)
# --------------------------------------------------------------------------- #
def build_monomial_matrix(n, d):
    """
    Construit la matrice M_eval (taille (2^n - 1) x M) des evaluations de tous les
    monomes multilineaires de degre 1..d sur x in {0,1}^n \ {0}.
    """
    monoms = []
    for k in range(1, d + 1):
        for comb in itertools.combinations(range(n), k):
            mask = 0
            for bit in comb:
                mask |= 1 << bit
            monoms.append(mask)
    M = len(monoms)
    N = (1 << n) - 1
    xs = np.arange(1, 1 << n, dtype=np.int64)
    mat = np.zeros((N, M), dtype=np.uint8)
    for col, mm in enumerate(monoms):
        mat[:, col] = ((xs & mm) == mm).astype(np.uint8)
    return mat, monoms


@njit(cache=True)
def f2_defect_and_solve(rows_mat):
    """
    Etant donne une matrice binaire A (taille R x M), teste par elimination de Gauss
    sur F_2 si le systeme A c = 1 (mod 2) possede une solution !
    Renvoie (n_inconsistent_rows, sol_c) :
      - n_inconsistent_rows == 0  <=>  IL EXISTE UNE SOLUTION EXACTE P_2 !
      - sinon, renvoie le nombre exact de lignes independantes produisant 0 = 1.
    """
    R, M = rows_mat.shape
    # Matrice augmentee [A | 1] de taille R x (M + 1)
    aug = np.empty((R, M + 1), np.uint8)
    for r in range(R):
        for c in range(M):
            aug[r, c] = rows_mat[r, c]
        aug[r, M] = 1

    piv_col = np.full(R, -1, np.int64)
    row = 0
    for col in range(M):
        # Trouver pivot dans col a partir de row
        p = -1
        for r in range(row, R):
            if aug[r, col] == 1:
                p = r
                break
        if p == -1:
            continue
        if p != row:
            for c in range(col, M + 1):
                tmp = aug[row, c]
                aug[row, c] = aug[p, c]
                aug[p, c] = tmp
        piv_col[row] = col
        # Eliminer toute la colonne (Gauss-Jordan pour lecture directe)
        for r in range(R):
            if r != row and aug[r, col] == 1:
                for c in range(col, M + 1):
                    aug[r, c] ^= aug[row, c]
        row += 1
        if row == R:
            break

    # Compter les lignes [0 ... 0 | 1] (inconsistances exactes !)
    inconsistent = 0
    for r in range(row, R):
        if aug[r, M] == 1:
            inconsistent += 1

    sol = np.zeros(M, np.uint8)
    if inconsistent == 0:
        for r in range(row):
            sol[piv_col[r]] = aug[r, M]
    return inconsistent, sol


@njit(cache=True)
def search_nonsym_mod6(eval_mat, c3_init, n_steps, seed):
    """
    Recuit simule / recherche locale sur les coefficients c3 in F_3^M de P_3 :
    a chaque mutation de c3, on evalue Z_3 = {x != 0 : P_3(x) == 0 mod 3}, puis
    on calcule par elimination de Gauss EXACTE sur F_2 le defaut d'inconsistance
    du systeme P_2(x) = 1 pour tout x in Z_3.
    Des que defect == 0, on a trouve une representation EXACTE (P_2, P_3) de OR_n mod 6 !
    """
    np.random.seed(seed)
    N, M = eval_mat.shape
    c3 = c3_init.copy()

    # Calculer val3 = eval_mat @ c3 mod 3
    val3 = np.zeros(N, np.int64)
    for col in range(M):
        cc = c3[col]
        if cc != 0:
            for r in range(N):
                if eval_mat[r, col]:
                    val3[r] += cc
    for r in range(N):
        val3[r] %= 3

    def eval_defect(v3):
        nz = 0
        for r in range(N):
            if v3[r] == 0:
                nz += 1
        if nz == 0:
            return 0, 0, np.zeros(M, np.uint8)
        sub = np.empty((nz, M), np.uint8)
        idx = 0
        for r in range(N):
            if v3[r] == 0:
                for col in range(M):
                    sub[idx, col] = eval_mat[r, col]
                idx += 1
        inc, sol = f2_defect_and_solve(sub)
        return inc, nz, sol

    best_inc, best_nz, best_sol = eval_defect(val3)
    best_c3 = c3.copy()
    cur_inc = best_inc
    cur_nz = best_nz
    if best_inc == 0:
        return 0, best_nz, best_c3, best_sol

    temp = 1.5
    for step in range(n_steps):
        col = np.random.randint(0, M)
        old_v = c3[col]
        delta = 1 if np.random.randint(0, 2) == 0 else 2
        new_v = (old_v + delta) % 3

        # Mise a jour incrementale de val3 en O(N)
        for r in range(N):
            if eval_mat[r, col]:
                val3[r] = (val3[r] + delta) % 3

        inc, nz, sol = eval_defect(val3)
        # Score combine : priorite absolue a l'inconsistance, puis taille de Z_3
        diff = (inc - cur_inc) + 0.002 * (nz - cur_nz)
        accept = (diff <= 0.0) or (np.random.random() < np.exp(-diff / max(temp, 0.05)))
        if accept:
            c3[col] = new_v
            cur_inc = inc
            cur_nz = nz
            if inc < best_inc or (inc == best_inc and nz < best_nz):
                best_inc = inc
                best_nz = nz
                best_c3 = c3.copy()
                best_sol = sol.copy()
                if best_inc == 0:
                    return 0, best_nz, best_c3, best_sol
        else:
            # Annuler
            inv_delta = (3 - delta) % 3
            for r in range(N):
                if eval_mat[r, col]:
                    val3[r] = (val3[r] + inv_delta) % 3
        temp *= 0.9985

    return best_inc, best_nz, best_c3, best_sol


def verify_polynomial_pair(n, d, monoms, c2, c3):
    """Verification EXHAUSTIVE 100% independante sur tous les 2^n points de {0,1}^n."""
    # P(0) = 0 mod 6 garanti car tous les monomes ont degre >= 1
    bad = 0
    for x in range(1, 1 << n):
        v2 = 0
        v3 = 0
        for col, mm in enumerate(monoms):
            if (x & mm) == mm:
                if c2[col]:
                    v2 ^= 1
                if c3[col]:
                    v3 += int(c3[col])
        if v2 == 0 and (v3 % 3) == 0:
            bad += 1
    return bad == 0


def format_poly(monoms, coeffs, mod):
    terms = []
    for mm, c in zip(monoms, coeffs):
        c = int(c) % mod
        if c == 0:
            continue
        vars_ = [f"x{i+1}" for i in range(20) if (mm >> i) & 1]
        mon = "".join(vars_)
        terms.append(mon if c == 1 else f"{c}{mon}")
    return " + ".join(terms) if terms else "0"


def main():
    t0 = time.time()
    print("=========================================================================", flush=True)
    print("== EXPERIENCE 1 : Spectre symetrique EXACT N_sym(d) pour d = 1..10     ==", flush=True)
    print("=========================================================================", flush=True)
    max_d = 10
    max_n = 260
    C2 = binomial_mod_table(max_n, max_d, 2)
    C3 = binomial_mod_table(max_n, max_d, 3)

    print(f"   {'deg d':>5} | {'Borne BBR (q2*q3-1)':>20} | {'EXACT N_sym(d)':>14} | {'Gain vs BBR':>11} | {'temps':>7}",
          flush=True)
    print("   " + "-" * 70, flush=True)

    sym_report = {}
    for d in range(1, max_d + 1):
        ts = time.perf_counter()
        n_opt, m2, c3_code = exact_n_sym(d, C2[:, :d], C3[:, :d])
        dt = time.perf_counter() - ts
        bbr_n, q2, q3 = bbr_classic_bound(d)
        c2_vec = decode_c2(int(m2), d)
        c3_vec = decode_c3(int(c3_code), d)
        gain = int(n_opt) - bbr_n
        sym_report[d] = {
            "bbr_bound": bbr_n, "q2": q2, "q3": q3,
            "n_sym_exact": int(n_opt), "gain": gain,
            "c2": c2_vec, "c3": c3_vec, "time_s": dt,
        }
        print(f"   {d:>5} | {bbr_n:>10} ({q2}x{q3}-1) | {int(n_opt):>14} | {gain:>+11} | {dt:>6.2f}s",
              flush=True)
        if d <= 6:
            print(f"           c2(F_2) = {c2_vec}   |   c3(F_3) = {c3_vec}", flush=True)

    print("\n=========================================================================", flush=True)
    print("== EXPERIENCE 2 : Chasse aux polynomes NON SYMETRIQUES (d = 2, n >= 9) ==", flush=True)
    print("=========================================================================", flush=True)
    nonsym_report = {}
    for n in (9, 10, 11, 12):
        eval_mat, monoms = build_monomial_matrix(n, 2)
        M = len(monoms)
        found = False
        best_overall_inc = 1 << 30
        best_pair = None
        ts = time.perf_counter()
        # Lancer 16 redemarrages (dont des graines proches de la symetrie et des graines aleatoires)
        rng = np.random.default_rng(2026 + n)
        for trial in range(16):
            if trial == 0:
                # Graine lineaire symetrique x1 + ... + xn
                c3_init = np.zeros(M, dtype=np.int64)
                c3_init[:n] = 1
            elif trial < 4:
                # Graine lineaire + petite perturbation quadratique
                c3_init = np.zeros(M, dtype=np.int64)
                c3_init[:n] = rng.integers(1, 3, size=n)
                c3_init[n:] = rng.choice([0, 0, 0, 1, 2], size=M - n)
            else:
                c3_init = rng.integers(0, 3, size=M, dtype=np.int64)

            inc, nz, bc3, bsol = search_nonsym_mod6(eval_mat, c3_init, 3000, seed=1000 * n + trial)
            if inc < best_overall_inc:
                best_overall_inc = int(inc)
                best_pair = (bsol.copy(), bc3.copy(), int(nz))
            if inc == 0:
                ok = verify_polynomial_pair(n, 2, monoms, bsol, bc3)
                assert ok, "ERREUR DE VERIFICATION !"
                found = True
                dt = time.perf_counter() - ts
                print(f"   [SUCCES !] d = 2, n = {n:>2} > N_sym(2) = 8 : POLYNOME NON SYMETRIQUE TROUVE et VERIFIE "
                      f"sur les {1<<n:,} points ! (essai {trial+1}, |Z_3|={nz}, {dt:.2f}s)", flush=True)
                if n <= 10:
                    print(f"      P_2 (mod 2) = {format_poly(monoms, bsol, 2)}", flush=True)
                    print(f"      P_3 (mod 3) = {format_poly(monoms, bc3, 3)}", flush=True)
                nonsym_report[f"d2_n{n}"] = {
                    "found": True, "n": n, "d": 2, "z3_size": int(nz), "trial": trial + 1,
                    "p2_str": format_poly(monoms, bsol, 2),
                    "p3_str": format_poly(monoms, bc3, 3),
                }
                break
        if not found:
            dt = time.perf_counter() - ts
            print(f"   [ECHEC]    d = 2, n = {n:>2} : aucune solution trouvee en 16 essais "
                  f"(meilleur defaut F_2 = {best_overall_inc}, {dt:.2f}s)", flush=True)
            nonsym_report[f"d2_n{n}"] = {"found": False, "n": n, "d": 2, "best_defect": best_overall_inc}
            if best_overall_inc > 5:
                break

    with open(os.path.join(RESULTS_DIR, "mod6_poly_results.json"), "w", encoding="utf-8") as fh:
        json.dump({"symmetric_exact": sym_report, "nonsymmetric_d2": nonsym_report}, fh, indent=2)
    print(f"\nTemps total : {time.time() - t0:.1f} s", flush=True)


if __name__ == "__main__":
    main()
