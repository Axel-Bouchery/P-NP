# Axe I — Paradigme de Williams, THR ∘ THR-SAT, Atlas n = 4, MCSP, CAPP & Topologie de Huang

Ce dossier autonome rassemble **17 scripts** Python/Numba (`src/`), **13 fichiers de résultats JSON** (`results/`), **1 témoin exact `.npz`** (`witnesses/atlas_n4.npz`), **13 figures** (`figures/`) et le rapport complet **[RAPPORT.md](RAPPORT.md)** dédiés aux circuits à seuil de profondeur 2 (`THR ∘ THR`) dans le paradigme algorithmique de Williams.

## Carte du dossier

| Script (`src/`) | Rôle | Données / Figures |
| :--- | :--- | :--- |
| [williams_demo.py](src/williams_demo.py) | Transformation zêta SYM ∘ AND (× 331,4 à n = 18) | [williams_demo.png](figures/williams_demo.png) |
| [polymethod_mq.py](src/polymethod_mq.py) | Méthode polynomiale LPTWY sur MQ / 𝔽₂ (60/60 exact) | [polymethod.png](figures/polymethod.png) |
| [thr2_lab.py](src/thr2_lab.py), [run_lab.py](src/run_lab.py), [test_thr2.py](src/test_thr2.py), [hybrid_thr2.py](src/hybrid_thr2.py), [adversarial_suite.py](src/adversarial_suite.py) | Solveur B&B évolutif (avec tests unitaires), suite adversariale et Hybride A+B | [thr2_results.json](results/thr2_results.json), [adversarial_results.json](results/adversarial_results.json), [hybrid_results.json](results/hybrid_results.json), [thr2_evolution.png](figures/thr2_evolution.png), [adversarial_suite.png](figures/adversarial_suite.png), [hybrid_comparison.png](figures/hybrid_comparison.png) |
| [theorem_a_star.py](src/theorem_a_star.py), [industrial_and_orthogonal.py](src/industrial_and_orthogonal.py), [step2_sliding_gray.py](src/step2_sliding_gray.py), [unified_v5_and_poly_top.py](src/unified_v5_and_poly_top.py) | Théorèmes 1 & 2 : réseau PGCD, duel industriel v4, Gray glissant + Solveur unifié v5 (n = 28 : 6 instances en 0,502 s) | [theorem_astar_results.json](results/theorem_astar_results.json), [phase4_results.json](results/phase4_results.json), [step2_results.json](results/step2_results.json), [phase6_v5_results.json](results/phase6_v5_results.json), [theorem_astar.png](figures/theorem_astar.png), [phase4_duel.png](figures/phase4_duel.png), [step2_sliding_gray.png](figures/step2_sliding_gray.png), [phase6_unified_v5.png](figures/phase6_unified_v5.png) |
| [exact_thr2_atlas.py](src/exact_thr2_atlas.py), [mcsp_metacomplexity.py](src/mcsp_metacomplexity.py) | Atlas exhaustif n = 4 (65 536 fonctions / 222 classes NPN) et échec des preuves naturelles L₁ | [atlas_n4.npz](witnesses/atlas_n4.npz), [atlas_n4_summary.json](results/atlas_n4_summary.json), [mcsp_results.json](results/mcsp_results.json), [atlas_and_mcsp.png](figures/atlas_and_mcsp.png) |
| [rank_k_bilinear_thr2.py](src/rank_k_bilinear_thr2.py), [piste_a_svd_capp.py](src/piste_a_svd_capp.py), [det_capp_nystrom.py](src/det_capp_nystrom.py), [piste_b_huang_topology.py](src/piste_b_huang_topology.py) | Théorème 3 (factorisation rang K / MIPS), CAPP déterministe Nyström (n = 26) et topologie de Huang | [theorem3_rank_results.json](results/theorem3_rank_results.json), [piste_a_svd_results.json](results/piste_a_svd_results.json), [det_capp_results.json](results/det_capp_results.json), [piste_b_huang_results.json](results/piste_b_huang_results.json), [theorem3_rank_gemm.png](figures/theorem3_rank_gemm.png), [piste_a_svd_capp.png](figures/piste_a_svd_capp.png), [piste_b_huang_topology.png](figures/piste_b_huang_topology.png) |

## Reproduire les calculs

Dépendances : Python 3.14, `numpy`, `numba`, `matplotlib`. Exemples d'exécution depuis la racine de cet axe (`axe1_thr2_williams/`) :

```powershell
python -u src/test_thr2.py
python -u src/williams_demo.py
python -u src/mcsp_metacomplexity.py
python -u src/unified_v5_and_poly_top.py
```
