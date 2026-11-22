# Complexité booléenne, polynômes modulo un entier composite et complexité algébrique

Ce dépôt est organisé en **trois axes de recherche entièrement autonomes** rassemblant au total **78 scripts** Python/Numba/PyTorch (77 scripts de recherche + 1 vérificateur exact), **61 fichiers de résultats JSON**, **32 témoins algébriques (`.npz`)**, **39 figures**, un manuscrit LaTeX autonome (`paper/`) et un projet Lean 4 (`lean/`). Les questions P ≟ NP et VP ≟ VNP servent uniquement de cadre motivant : aucun résultat présenté ici ne les résout.

## Index des trois axes

| Axe | Dossier | Contenu | Documentation |
| :--- | :--- | :--- | :--- |
| **Axe I** | **[axe1_thr2_williams/](axe1_thr2_williams/README.md)** | Paradigme de Williams, solveurs exacts `THR ∘ THR-SAT` (Théorèmes 1–3, v5 jusqu'à `n = 28`), Atlas exhaustif `n = 4`, MCSP, CAPP déterministe CUR & topologie de Huang (`17` scripts, `13` JSON, `1` `.npz`, `13` figures) | [README.md](axe1_thr2_williams/README.md) · [RAPPORT.md](axe1_thr2_williams/RAPPORT.md) |
| **Axe II** | **[axe2_mod_composite_polynomials/](axe2_mod_composite_polynomials/README.md)** | Représentation faible modulo un entier composite `m` (`OR_n`, `THR_{n, t}`, `MAJ_n`, `MOD_k`) : 26 théorèmes (`M1–M2`, `G1–G5`, `G3++`, `T1–T2`, `Q1`, `P1`, `U1–U4`, `D1–D10`), manuscrit LaTeX et formalisation Lean 4 (`41` scripts, `33` JSON, `18` figures) | [README.md](axe2_mod_composite_polynomials/README.md) · [RAPPORT.md](axe2_mod_composite_polynomials/RAPPORT.md) · [paper/](axe2_mod_composite_polynomials/paper/weak_poly_composite.tex) · [lean/](axe2_mod_composite_polynomials/lean/README.md) |
| **Axe III** | **[axe3_algebraic_complexity/](axe3_algebraic_complexity/README.md)** | GCT sur `perm₃` (Parseval-LM, piège `SL₇`) et tenseurs de multiplication matricielle (recuit `L_{2/3}`, loi `α`, barrière `L_∞`, Graphe de Flips `𝔽₃/𝔽₅` et `uint64` sur `𝔽₂` sur 16 formats) (`20` scripts, `15` JSON, `31` `.npz`, `8` figures) | [README.md](axe3_algebraic_complexity/README.md) · [RAPPORT.md](axe3_algebraic_complexity/RAPPORT.md) |

---

## Extraire chaque axe en dépôt Git indépendant (`git filter-repo`)

Chaque dossier d'axe étant désormais autonome (`src/`, `results/`, `figures/`, `witnesses/`, `RAPPORT.md`, `README.md`, et `paper/` + `lean/` pour l'Axe II), il peut être extrait en un dépôt Git indépendant tout en préservant l'historique des commits via [`git-filter-repo`](https://github.com/newren/git-filter-repo) :

```bash
# 1. Extraire l'Axe I (THR o THR / Williams)
git clone <URL_DU_DEPOT> repo-axe1-thr2-williams
cd repo-axe1-thr2-williams
git filter-repo --subdirectory-filter axe1_thr2_williams

# 2. Extraire l'Axe II (Polynomes faibles modulo m composite)
git clone <URL_DU_DEPOT> repo-axe2-mod-composite
cd repo-axe2-mod-composite
git filter-repo --subdirectory-filter axe2_mod_composite_polynomials

# 3. Extraire l'Axe III (Complexite algebrique : perm_3 & MatMul)
git clone <URL_DU_DEPOT> repo-axe3-algebraic-complexity
cd repo-axe3-algebraic-complexity
git filter-repo --subdirectory-filter axe3_algebraic_complexity
```
