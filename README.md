# Complexité booléenne, polynômes modulo un entier composite et complexité algébrique

Ce dépôt est organisé en **trois axes de recherche entièrement autonomes** rassemblant au total **80 scripts** Python/Numba/PyTorch (79 scripts de recherche + 1 vérificateur exact), **63 fichiers de résultats JSON**, **32 témoins algébriques (`.npz`)**, **41 figures**, un manuscrit LaTeX autonome (`paper/`) et un projet Lean 4 (`lean/`). Les questions P ≟ NP et VP ≟ VNP servent uniquement de cadre motivant : aucun résultat présenté ici ne les résout.

## Index des trois axes

| Axe | Dossier | Contenu | Documentation |
| :--- | :--- | :--- | :--- |
| **Axe I** | **[axe1_thr2_williams/](axe1_thr2_williams/README.md)** | Paradigme de Williams, solveurs exacts `THR ∘ THR-SAT` (Théorèmes 1–3, v5 jusqu'à `n = 28`), Atlas exhaustif `n = 4`, MCSP, CAPP déterministe CUR & topologie de Huang (`17` scripts, `13` JSON, `1` `.npz`, `13` figures) | [README.md](axe1_thr2_williams/README.md) · [RAPPORT.md](axe1_thr2_williams/RAPPORT.md) |
| **Axe II** | **[axe2_mod_composite_polynomials/](axe2_mod_composite_polynomials/README.md)** | Représentation faible modulo un entier composite `m` (`OR_n`, `THR_{n, t}`, `MAJ_n`, `MOD_k`) : 28 théorèmes (`M1–M2`, `G1–G5`, `G3++`, `T1–T2`, `Q1`, `P1`, `U1–U4`, `D1–D12`), manuscrit LaTeX et formalisation Lean 4 (`43` scripts, `35` JSON, `20` figures) | [README.md](axe2_mod_composite_polynomials/README.md) · [RAPPORT.md](axe2_mod_composite_polynomials/RAPPORT.md) · [paper/](axe2_mod_composite_polynomials/paper/weak_poly_composite.tex) · [lean/](axe2_mod_composite_polynomials/lean/README.md) |
| **Axe III** | **[axe3_algebraic_complexity/](axe3_algebraic_complexity/README.md)** | GCT sur `perm₃` (Parseval-LM, piège `SL₇`) et tenseurs de multiplication matricielle (recuit `L_{2/3}`, loi `α`, barrière `L_∞`, Graphe de Flips `𝔽₃/𝔽₅` et `uint64` sur `𝔽₂` sur 16 formats) (`20` scripts, `15` JSON, `31` `.npz`, `8` figures) | [README.md](axe3_algebraic_complexity/README.md) · [RAPPORT.md](axe3_algebraic_complexity/RAPPORT.md) |

---

## Arborescence du projet

```text
.
├── README.md                              # Index général et commandes d'extraction Git
├── .gitignore
│
├── axe1_thr2_williams/                    # Axe I — Paradigme de Williams & THR ∘ THR
│   ├── README.md                          # Carte et commandes d'exécution de l'Axe I
│   ├── RAPPORT.md                         # Monographie autonome (Expériences E1–E10)
│   ├── src/                               # 17 scripts Python/Numba
│   ├── results/                           # 13 fichiers de résultats JSON
│   ├── witnesses/                         # 1 témoin exact (.npz : atlas_n4.npz)
│   └── figures/                           # 13 figures (.png)
│
├── axe2_mod_composite_polynomials/        # Axe II — Polynômes faibles modulo m composite
│   ├── README.md                          # Carte et commandes d'exécution de l'Axe II
│   ├── RAPPORT.md                         # Monographie autonome (Expériences E11–E15q)
│   ├── src/                               # 42 scripts Python/Numba
│   ├── results/                           # 35 fichiers de résultats JSON
│   ├── figures/                           # 20 figures (.png)
│   ├── paper/                             # Manuscrit LaTeX autonome (weak_poly_composite.tex, references.bib)
│   └── lean/                              # Formalisation Lean 4 (14 .lean + verify_proof_steps_exact.py)
│
└── axe3_algebraic_complexity/             # Axe III — GCT sur perm₃ & Multiplication matricielle
    ├── README.md                          # Carte et commandes d'exécution de l'Axe III
    ├── RAPPORT.md                         # Monographie autonome (Expériences E16–E23)
    ├── src/                               # 20 scripts Python/Numba/PyTorch
    ├── results/                           # 15 fichiers de résultats JSON
    ├── witnesses/                         # 31 témoins algébriques (.npz)
    └── figures/                           # 8 figures (.png)
```

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
