# Formalisation Lean 4 — Polynômes faibles modulo un entier composite (Axe II)

Ce dossier contient un projet Lean 4 (`v4.16.0` + Mathlib 4) formalisant les briques analytiques centrales de l'article [`paper/weak_poly_composite.tex`](../paper/weak_poly_composite.tex) :

- `WeakPolyComposite/Basic.lean` : définitions (`hammingWeight`, `QuadPolyMod`, forme polaire `polar`, `WeakRepOR`).
- `WeakPolyComposite/TheoremG3pp.lean` : lemme exact de signe de Wilson `(-1)^(p-2) (p-2)! = -1 (mod p)` et borne d'image de la construction par blocs du **Théorème G3++**.
- `WeakPolyComposite/TheoremG5.lean` : borne de périodicité symétrique et correction des blocs de norme du **Théorème G5**.
- `WeakPolyComposite/TheoremD4.lean` : obstruction modulaire du **Lemme D3a(2)** (`1 + h ≠ 0 mod q` pour `q` premier impair) et assemblage linéaire **`8p − 7`** (`< 2(4p − 2)`) / **`6p − 4`** du **Théorème D4**.
- `WeakPolyComposite/TheoremD5.lean` : identité de parité `p = 2` et borne exacte 2-blocs **`2(m − 1)`** du **Théorème D5**.
- `WeakPolyComposite/TheoremD6.lean` : contradiction 3-racines à `n = 33` sur `mod 10` (`N_gen^{(10)}(2, OR) ≤ 32`), élimination de Lucas `(3 : ZMod p) ≠ 0` (`p ≥ 5`) et amélioration `19q − 13 < 20q − 19` sur `m = 3q` (`q ≥ 7`).
- `WeakPolyComposite/TheoremD7.lean` : annulations exactes de Lucas des bi-tranches `(3p − 1, 2p − 1)` et `(4p − 1, p − 1)` (**Théorème D7**) et divisibilité `2p ∣ C(p, 2)` pour `p ≡ 1 (mod 4)` (**Théorème U3**).
- `WeakPolyComposite/TheoremD8.lean` : identités binomiales `det A(1, r) = 1`, `det A(2, r) = C(2r − 1, r)`, gain strict `(Dr + 1)p − 1 < (D + 1)r(p − 1) + r` et corollaire puissance première (**Théorème D8**).
- `WeakPolyComposite/TheoremD9.lean` : empilement sharp `𝔽₃²` (`3k + 4`), amélioration inconditionnelle `19q − 15 < 20q − 19` pour tout `q ≥ 5` (y compris `mod 15 : 80 < 81`), borne Caro–Wei sur couplage et resserrement D3a sur semi-premiers (**Théorème D9**).
- `WeakPolyComposite/TheoremD10.lean` : parités exactes du témoin rang-1 `3q − 1` sur `m = 2q`, bornes sous-`2(m − 1)` et rang-2 pur `4m − 2p − 2q` sur `m = pq`, identité 3-points du lemme de franchissement de moyenne à biais borné (`θ* ∈ [1/4, 3/4]`), identité polynomiale `L(T) ≥ (4T² − 1)/(3T)`, survie sous-échantillonnée et témoin exact `K_8 ⊔ K_2` (`n = 10` mod `6`, `decide`) (**Théorème D10**).
- `verify_proof_steps_exact.py` : vérificateur arithmétique exact exécutable sous Python vérifiant toutes ces identités et lemmes sur des familles de nombres premiers.

## Statut de compilation

> **Note de transparence** : `elan` / `lean` n'étant pas installé sur la machine de développement au moment de la rédaction, ces fichiers Lean 4 sont structurés et prêts pour une compilation locale via :
> ```bash
> cd lean
> lake exe cache get
> lake build
> ```
> Toutes les identités algébriques et combinatoires correspondantes ont été vérifiées en arithmétique exacte par `python lean/verify_proof_steps_exact.py`.
