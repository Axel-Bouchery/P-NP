"""
Methode polynomiale de Williams, version complete (avec regroupement de variables)
==================================================================================

Probleme NP-complet attaque : MQ sur GF(2)
  Etant donne m polynomes quadratiques p_1..p_m en n variables booleennes,
  trouver x tel que p_i(x) = 0 pour tout i.
  (C'est le probleme a la base de la cryptographie "multivariee".)

Algorithme (Lokshtanov-Paturi-Tamaki-Williams-Yu 2017, meme idee que Williams 2011) :
  1. On coupe x = (y, z) avec |z| = l variables "regroupees".
  2. Razborov-Smolensky : on remplace le ET des m equations par
         R(x) = prod_{j=1..r} (1 + sum_{i in S_j} p_i(x))      (S_j aleatoires)
     R(x) = 1 si x est solution, et R(x) = 0 sinon avec proba 1 - 2^-r.
     R est un polynome de degre <= 2r : c'est le "polynome probabiliste".
  3. F(y) = sum_z R(y, z) mod 2  = parite du nombre d'extensions solutions de y.
     Sommer sur TOUS les z tue tout monome qui ne contient pas z_1...z_l,
     donc deg F <= 2r - l : F est un polynome de PETIT degre en n - l variables.
  4. Un polynome de degre D est determine par ses valeurs sur les points de
     poids de Hamming <= D (beaucoup moins que 2^(n-l)) : on ne calcule F que la,
     on interpole (transformee de Mobius), puis on evalue F sur les 2^(n-l)
     points d'un coup (transformee zeta).  ->  cout ~ 2^(n-l) au lieu de 2^n.
  5. Les y avec F(y) = 1 sont des candidats : on verifie (2^l essais chacun).
  6. Isolation (Valiant-Vazirani) : si toutes les parites sont paires, on ajoute
     des equations lineaires aleatoires et on recommence.

Garantie : si l'algorithme renvoie x, x est une vraie solution (verifiee).
           S'il repond "pas de solution", il peut se tromper avec faible proba.
"""
from pathlib import Path
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
AXE_ROOT = Path(__file__).resolve().parent.parent
FIGURES_DIR = str(AXE_ROOT / "figures")
os.makedirs(FIGURES_DIR, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", os.path.join(HERE, ".mpl"))

EVALS = 0  # compteur : nombre d'evaluations de polynomes quadratiques en un point


# --------------------------------------------------------------------------- #
#  Polynomes quadratiques sur GF(2) : (Q triangulaire stricte, L lineaire, c)
# --------------------------------------------------------------------------- #
def random_poly(n, rng):
    Q = np.triu(rng.integers(0, 2, (n, n)), 1).astype(np.float32)
    L = rng.integers(0, 2, n).astype(np.float32)
    return Q, L, int(rng.integers(0, 2))


def eval_poly(p, X):
    """p evalue sur les lignes de X (N x n, 0/1 en float32). Exact tant que n < 4000."""
    global EVALS
    Q, L, c = p
    EVALS += X.shape[0]
    v = ((X @ Q) * X).sum(axis=1) + X @ L + c
    return v.astype(np.int64) & 1


def add_polys(polys):
    n = polys[0][0].shape[0]
    Q = np.zeros((n, n), np.float32)
    L = np.zeros(n, np.float32)
    c = 0
    for q, l, cc in polys:
        Q += q
        L += l
        c += cc
    return np.mod(Q, 2), np.mod(L, 2), c & 1


def planted_instance(n, m, rng):
    """m equations aleatoires avec une solution cachee x*."""
    secret = rng.integers(0, 2, n).astype(np.float32)[None, :]
    polys = []
    for _ in range(m):
        Q, L, _ = random_poly(n, rng)
        val = int(((secret @ Q) * secret).sum() + (secret @ L).item()) & 1
        polys.append((Q, L, val))  # p(x*) = val + val = 0
    return polys


def random_instance(n, m, rng):
    return [random_poly(n, rng) for _ in range(m)]


def to_bits(idx, nbits):
    return ((idx[:, None] >> np.arange(nbits, dtype=np.int64)) & 1).astype(np.float32)


def surviving_rows(polys, X):
    """Indices des lignes de X ou tous les polynomes s'annulent (filtrage progressif)."""
    alive = np.arange(X.shape[0])
    for p in polys:
        if alive.size == 0:
            break
        alive = alive[eval_poly(p, X[alive]) == 0]
    return alive


def is_solution(polys, x):
    X = np.asarray(x, np.float32)[None, :]
    return all(eval_poly(p, X)[0] == 0 for p in polys)


# --------------------------------------------------------------------------- #
#  Force brute (avec filtrage progressif : ~2 evaluations par point en moyenne)
# --------------------------------------------------------------------------- #
def brute_force(polys, n, find_all=False, chunk=1 << 18):
    sols = []
    for start in range(0, 1 << n, chunk):
        idx = np.arange(start, min(start + chunk, 1 << n), dtype=np.int64)
        alive = surviving_rows(polys, to_bits(idx, n))
        if alive.size:
            if not find_all:
                return to_bits(idx[alive[:1]], n)[0]
            sols.extend(idx[alive].tolist())
    return sols if find_all else None


# --------------------------------------------------------------------------- #
#  Transformee de Mobius / zeta sur GF(2) (identiques mod 2) : O(k * 2^k)
# --------------------------------------------------------------------------- #
def mobius_gf2(a, k):
    for i in range(k):
        v = a.reshape(-1, 2, 1 << i)
        v[:, 1, :] ^= v[:, 0, :]
    return a


def popcount_table(k):
    idx = np.arange(1 << k, dtype=np.int64)
    pc = np.zeros(1 << k, np.uint8)
    for b in range(k):
        pc += ((idx >> b) & 1).astype(np.uint8)
    return pc


# --------------------------------------------------------------------------- #
#  Modele de cout (en evaluations de polynomes) et choix des parametres
# --------------------------------------------------------------------------- #
def low_weight_count(k, D):
    return sum(math.comb(k, i) for i in range(min(k, D) + 1))


def cost_poly(n, l, c):
    """Cout estime d'une passe de la methode polynomiale (en evaluations)."""
    mp, r = n - l, l + c
    D = min(mp, l + 2 * c)
    interpolation = 2.0 * low_weight_count(mp, D) * (1 << l)   # F sur poids <= D
    transforms = 2.0 * mp * 2.0 ** mp / (n * n)                  # 2 transformees bit a bit
    verification = 2.0 * 2.0 ** (n - c)                          # faux positifs ~ 2^-c
    return interpolation + transforms + verification


def cost_brute(n):
    return 2.0 * 2.0 ** n


def best_params(n):
    best = None
    for l in range(1, max(2, n // 3)):
        for c in range(1, l + 12):
            if l + 2 * c > n - l:
                break
            cst = cost_poly(n, l, c)
            if best is None or cst < best[0]:
                best = (cst, l, c)
    return best[1], best[2]


# --------------------------------------------------------------------------- #
#  Methode polynomiale : une passe
# --------------------------------------------------------------------------- #
def poly_method_pass(system, original, n, l, c, rng, chunk_rows=1 << 18):
    mp, r = n - l, l + c
    D = min(mp, l + 2 * c)
    nz = 1 << l

    # (2) polynome probabiliste de Razborov-Smolensky : r combinaisons aleatoires
    combos = []
    for _ in range(r):
        S = [p for p in system if rng.integers(0, 2)]
        combos.append(add_polys(S) if S else (np.zeros((n, n), np.float32),
                                                np.zeros(n, np.float32), 0))

    # (3)-(4) F(y) = sum_z R(y,z) sur les seuls points y de poids <= D
    pc = popcount_table(mp)
    low = np.nonzero(pc <= D)[0].astype(np.int64)
    F = np.zeros(1 << mp, np.uint8)
    zs = np.arange(nz, dtype=np.int64) << mp
    B = max(1, chunk_rows // nz)
    for s in range(0, low.size, B):
        ys = low[s:s + B]
        X = to_bits((ys[:, None] | zs[None, :]).ravel(), n)
        alive = surviving_rows(combos, X)              # R(y,z) = 1  <=>  combos nulles
        F[ys] = (np.bincount(alive >> l, minlength=ys.size) & 1).astype(np.uint8)

    # interpolation (Mobius), troncature au degre D, puis evaluation partout (zeta)
    mobius_gf2(F, mp)
    F[pc > D] = 0
    mobius_gf2(F, mp)

    # (5) verification des candidats sur le systeme ORIGINAL
    cand = np.nonzero(F)[0].astype(np.int64)
    for s in range(0, cand.size, B):
        ys = cand[s:s + B]
        idx = (ys[:, None] | zs[None, :]).ravel()
        alive = surviving_rows(original, to_bits(idx, n))
        if alive.size:
            return to_bits(idx[alive[:1]], n)[0], cand.size
    return None, cand.size


def poly_method(polys, n, rng, l=None, c=None, trials=2):
    if l is None:
        l, c = best_params(n)
    for k in range(0, n + 1):                           # (6) isolation
        iso = [(np.zeros((n, n), np.float32),
                rng.integers(0, 2, n).astype(np.float32),
                int(rng.integers(0, 2))) for _ in range(k)]
        for _ in range(trials):
            sol, _ = poly_method_pass(polys + iso, polys, n, l, c, rng)
            if sol is not None:
                return sol, (l, c, k)
    return None, (l, c, n)


# --------------------------------------------------------------------------- #
#  Experiences
# --------------------------------------------------------------------------- #
def correctness(rng):
    print("== 1. Verification de correction (comparaison avec la force brute) ==")
    ok = total = 0
    for n in range(10, 19, 2):
        for t in range(12):
            if t % 2 == 0:
                polys, kind = planted_instance(n, n, rng), "SAT plante"
            else:
                polys, kind = random_instance(n, n + 12, rng), "aleatoire"
            truth = brute_force(polys, n) is not None
            sol, _ = poly_method(polys, n, rng)
            got = sol is not None
            if got:
                assert is_solution(polys, sol), "solution invalide !"
            total += 1
            ok += truth == got
    print(f"   {ok}/{total} instances : meme reponse que la force brute "
          f"(toute solution renvoyee est verifiee)\n")


def benchmark(rng):
    global EVALS
    print("== 2. Mesures reelles (instances plantees, m = n equations) ==")
    print(f"{'n':>3} | {'l':>2} {'c':>2} | {'brute: evals':>13} {'temps':>7} | "
          f"{'polynomiale: evals':>18} {'temps':>7} | ratio evals")
    print("-" * 82)
    rows = []
    for n in range(14, 25, 2):
        polys = planted_instance(n, n, rng)
        EVALS = 0
        t0 = time.perf_counter()
        brute_force(polys, n, find_all=True)            # pire cas : tout parcourir
        tb, eb = time.perf_counter() - t0, EVALS

        l, c = best_params(n)
        EVALS = 0
        t0 = time.perf_counter()
        sol, (_, _, k) = poly_method(polys, n, rng, l, c)
        tp, ep = time.perf_counter() - t0, EVALS
        assert sol is not None and is_solution(polys, sol)
        rows.append((n, eb, ep, tb, tp))
        print(f"{n:>3} | {l:>2} {c:>2} | {eb:>13,} {tb:>6.2f}s | {ep:>18,} {tp:>6.2f}s | "
              f"{ep / eb:>6.2f}")
    print()
    return rows


def model():
    print("== 3. Modele de cout asymptotique (evaluations de polynomes) ==")
    ns = list(range(20, 201, 5))
    lb, lp, params = [], [], []
    cross = None
    for n in ns:
        l, c = best_params(n)
        b, p = math.log2(cost_brute(n)), math.log2(cost_poly(n, l, c))
        lb.append(b)
        lp.append(p)
        params.append((l, c))
        if cross is None and p < b:
            cross = n
    slope = (lp[-1] - lp[-9]) / (ns[-1] - ns[-9])
    for n in (40, 60, 100, 150, 200):
        i = ns.index(n)
        print(f"   n = {n:>3} : force brute 2^{lb[i]:6.1f}   polynomiale 2^{lp[i]:6.1f}"
              f"   (l = {params[i][0]}, c = {params[i][1]})")
    print(f"   Croisement (la methode polynomiale devient meilleure) : n ~ {cross}")
    print(f"   Exposant mesure a grand n : cout ~ 2^({slope:.3f} n)  "
          f"(theorie LPTWY : 2^(0.9 n) pour des quadratiques)\n")
    return ns, lb, lp, cross, slope


def plot(rows, ns, lb, lp, cross, slope, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 4.8))
    n_small = [r[0] for r in rows]
    a1.semilogy(n_small, [r[1] for r in rows], "o-", color="#d62728", label="Force brute")
    a1.semilogy(n_small, [r[2] for r in rows], "o-", color="#2ca02c", label="Méthode polynomiale")
    a1.set_title("Mesures réelles (n ≤ 24) : moins d'évaluations")
    a1.set_xlabel("Nombre de variables n")
    a1.set_ylabel("Évaluations de polynômes (log)")
    a1.grid(alpha=0.3, which="both")
    a1.legend()

    a2.plot(ns, lb, color="#d62728", label="Force brute : 2^n")
    a2.plot(ns, lp, color="#2ca02c", label=f"Méthode polynomiale : ~2^({slope:.2f}·n)")
    a2.set_title("Modèle asymptotique : l'écart se creuse exponentiellement")
    a2.set_xlabel("Nombre de variables n")
    a2.set_ylabel("log2(coût)")
    a2.grid(alpha=0.3)
    a2.legend()
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    print(f"Graphique : {out}")


def main():
    rng = np.random.default_rng(42)
    correctness(rng)
    rows = benchmark(rng)
    ns, lb, lp, cross, slope = model()
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(FIGURES_DIR, "polymethod.png")
    plot(rows, ns, lb, lp, cross, slope, out)


if __name__ == "__main__":
    main()
