"""
Demo du coeur technique de l'approche de Ryan Williams (2011) :
"un algorithme de SAT un peu meilleur que la force brute => une borne inferieure".

Dans la preuve NEXP ⊄ ACC0, l'etape algorithmique est :
  1. tout circuit ACC0 se reecrit en circuit SYM∘AND (une porte symetrique
     au sommet, des ET de litteraux en dessous) de taille s quasi-polynomiale
     (Yao, Beigel-Tarui) ;
  2. on evalue ce circuit sur TOUTES les 2^n entrees en temps ~ 2^n * poly(n),
     c'est-a-dire en temps quasi INDEPENDANT de la taille s du circuit.

Ce script illustre le point 2 :
  - methode naive  : pour chaque entree, on evalue chaque porte ET -> ~ s * 2^n
  - methode rapide : on developpe chaque ET en polynome multilineaire, puis une
    transformee "zeta" (somme sur les sous-ensembles) donne le nombre de ET
    satisfaits pour toutes les entrees en ~ n * 2^n (+ s * 2^k pour developper)

La porte symetrique du sommet est un seuil : le circuit accepte x si au moins
`theta` portes ET sont satisfaites. On resout SAT et meme #SAT (comptage).
"""
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault("MPLCONFIGDIR", os.path.join(HERE, ".mpl"))


def random_sym_and_circuit(n, s, k, rng):
    """s portes ET, chacune sur k litteraux aleatoires (positifs ou negatifs)."""
    terms = []
    for _ in range(s):
        vars_ = rng.choice(n, size=k, replace=False)
        signs = rng.integers(0, 2, size=k)
        pos = 0
        neg = 0
        for v, sg in zip(vars_, signs):
            if sg:
                pos |= 1 << int(v)
            else:
                neg |= 1 << int(v)
        terms.append((pos, neg))
    return terms


def counts_naive(n, terms):
    """Nombre de portes ET satisfaites pour chaque entree x : cout ~ s * 2^n."""
    X = np.arange(1 << n, dtype=np.int64)
    cnt = np.zeros(1 << n, dtype=np.int64)
    for pos, neg in terms:
        cnt += ((X & pos) == pos) & ((X & neg) == 0)
    return cnt


def counts_fast(n, terms):
    """Meme resultat via polynome multilineaire + transformee zeta : ~ s*2^k + n*2^n."""
    coef = np.zeros(1 << n, dtype=np.int64)
    for pos, neg in terms:
        # prod_{i in pos} x_i * prod_{j in neg} (1 - x_j)
        #   = sum_{S subset of neg} (-1)^{|S|} * x^{pos ∪ S}
        sub = neg
        while True:
            sign = -1 if bin(sub).count("1") & 1 else 1
            coef[pos | sub] += sign
            if sub == 0:
                break
            sub = (sub - 1) & neg
    # Transformee zeta : f(x) = sum_{T subset of x} coef[T]
    a = coef
    for i in range(n):
        a = a.reshape(-1, 2, 1 << i)
        a[:, 1, :] += a[:, 0, :]
        a = a.reshape(-1)
    return a


def main():
    n = 18
    k = 6
    sizes = [100, 500, 1000, 2000, 5000, 10000, 20000]
    rng = np.random.default_rng(2026)

    print(f"n = {n} variables  ->  2^n = {1 << n:,} entrees ; ET de largeur k = {k}\n")
    print(f"{'taille s':>9} | {'naif (s)':>9} | {'rapide (s)':>10} | {'gain':>6} | "
          f"{'max ET sat.':>11} | {'#SAT (seuil=max)':>16} | ok")
    print("-" * 82)

    rows = []
    for s in sizes:
        terms = random_sym_and_circuit(n, s, k, rng)

        t0 = time.perf_counter()
        c1 = counts_naive(n, terms)
        t_naive = time.perf_counter() - t0

        t0 = time.perf_counter()
        c2 = counts_fast(n, terms)
        t_fast = time.perf_counter() - t0

        same = bool(np.array_equal(c1, c2))
        theta = int(c2.max())                 # porte seuil : "au moins theta ET vrais"
        nb_sat = int((c2 >= theta).sum())     # #SAT du circuit seuil∘ET
        rows.append((s, t_naive, t_fast))
        print(f"{s:>9,} | {t_naive:>9.3f} | {t_fast:>10.3f} | {t_naive / t_fast:>5.1f}x | "
              f"{theta:>11} | {nb_sat:>16,} | {'OK' if same else 'ERREUR'}")
        if not same:
            sys.exit("Les deux methodes divergent !")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    s_vals = [r[0] for r in rows]
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.plot(s_vals, [r[1] for r in rows], "o-", color="#d62728",
            label="Naïf : évaluer chaque porte sur chaque entrée (~ s·2ⁿ)")
    ax.plot(s_vals, [r[2] for r in rows], "o-", color="#2ca02c",
            label="Polynôme + transformée zêta (~ n·2ⁿ + s·2ᵏ)")
    ax.set_xlabel("Taille du circuit s (nombre de portes ET)")
    ax.set_ylabel("Temps (secondes)")
    ax.set_title(f"Évaluer un circuit SYM∘AND sur les 2^{n} = {1 << n:,} entrées")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "williams_demo.png")
    fig.savefig(out, dpi=130)
    print(f"\nGraphique : {out}")


if __name__ == "__main__":
    main()
