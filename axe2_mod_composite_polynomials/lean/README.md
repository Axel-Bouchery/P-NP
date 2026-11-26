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
- `WeakPolyComposite/TheoremD11.lean` : identités quadratiques exactes du lemme universel de cosinus déphasé (`(3/256) π² > 0` sur `∣α∣ ∈ (π/4, π/2]`, `α² + (3/4) δ² ≥ 0` sur `∣α∣ ≤ π/4`) et contraction de densité `72/121 < 2/3` du cœur quasi-régulier (**Théorème D11**).
- `WeakPolyComposite/TheoremD12.lean` : certificats exacts (`decide`) par tranches de Lucas sur `mod 6` à `d = 3, n = 28` (`P_3*` et tranches simples `w = 8, 16`), dual `d = 3, n = 40` (`I_4(S)`, 4 termes mod `3`), et clôture `d = 6, n = 63` (`I_9(S)`) (**Théorème D12**).
- `WeakPolyComposite/TheoremD13.lean` : constante pré-Weyl Bande 3 (`3/16 = 4 × 3/64`), annulation exacte de tout poids de cluster `λ` par repliement de paires à poids constant, et borne linéaire `< 8pq` sur graphes cocktail-party (`ρ = 1`) par rigidité de Valiant (**Théorème D13**).
- `WeakPolyComposite/TheoremD14.lean` : seuil de contraction `OR_3` (`17`), élimination inconditionnelle de `n ≥ 30` sur `mod 10` depuis l'escalier héréditaire de Lucas et amélioration du ratio `< 1,56×` (`28 / 18`) (**Théorème D14**).
- `WeakPolyComposite/TheoremU4.lean` : surjectivité exacte (`decide`) des décalages affines des formes normiques anisotropes sur `𝔽₃²` et `𝔽₅²` (lemme d'obstruction 3-blocs anisotrope), amélioration stricte `2(m − 1) < 2m + q − 3` (`18, 22, 28`), dénombrement exact des orbites `ℤ₁₉` (`19 × 27 594 + 1 = 2¹⁹ − 1`) et obstructions ponctuelles `C₁₉` / `P₁₉` (**Théorème U4**).
- `verify_proof_steps_exact.py` : vérificateur arithmétique exact exécutable sous Python vérifiant toutes ces identités et lemmes sur des familles de nombres premiers.

## Statut de compilation

> **Compilation** : `elan` / `lean` n'étant pas installé dans l'environnement d'exécution initial, ces fichiers Lean 4 sont prêts pour une compilation locale via :
> ```bash
> cd lean
> lake exe cache get
> lake build
> ```
> Toutes les identités algébriques et combinatoires correspondantes ont été vérifiées en arithmétique exacte par `python lean/verify_proof_steps_exact.py`.
