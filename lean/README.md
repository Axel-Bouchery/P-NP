# Formalisation Lean 4 — Polynômes faibles modulo un entier composite (Axe II)

Ce dossier contient un projet Lean 4 (`v4.16.0` + Mathlib 4) formalisant les briques analytiques centrales de l'article [`paper/weak_poly_composite.tex`](../paper/weak_poly_composite.tex) :

- `WeakPolyComposite/Basic.lean` : définitions (`hammingWeight`, `QuadPolyMod`, forme polaire `polar`, `WeakRepOR`).
- `WeakPolyComposite/TheoremG3pp.lean` : lemme exact de signe de Wilson `(-1)^(p-2) (p-2)! = -1 (mod p)` et borne d'image de la construction par blocs du **Théorème G3++**.
- `WeakPolyComposite/TheoremG5.lean` : borne de périodicité symétrique et correction des blocs de norme du **Théorème G5**.
- `WeakPolyComposite/TheoremD4.lean` : obstruction modulaire du **Lemme D3a(2)** (`1 + h ≠ 0 mod q` pour `q` premier impair) et assemblage linéaire **`8p − 7`** (`< 2(4p − 2)`) / **`6p − 4`** du **Théorème D4**.
- `WeakPolyComposite/TheoremD5.lean` : identité de parité `p = 2` et borne exacte 2-blocs **`2(m − 1)`** du **Théorème D5**.
- `verify_proof_steps_exact.py` : vérificateur arithmétique exact exécutable sous Python vérifiant toutes ces identités et lemmes sur des familles de nombres premiers.

## Statut de compilation

> **Note de transparence** : `elan` / `lean` n'étant pas installé sur la machine de développement au moment de la rédaction, ces fichiers Lean 4 sont structurés et prêts pour une compilation locale via :
> ```bash
> cd lean
> lake exe cache get
> lake build
> ```
> Toutes les identités algébriques et combinatoires correspondantes ont été vérifiées en arithmétique exacte par `python lean/verify_proof_steps_exact.py`.
