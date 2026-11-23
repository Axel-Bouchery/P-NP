# Axe III — Complexité algébrique : GCT sur perm₃ & Multiplication matricielle

Ce dossier autonome rassemble **20 scripts** Python/Numba/PyTorch (`src/`), **15 fichiers de résultats JSON** (`results/`), **31 témoins algébriques `.npz`** (`witnesses/`), **8 figures** (`figures/`) et le rapport complet **[RAPPORT.md](RAPPORT.md)** dédiés à la complexité déterminantale de `perm₃` et au rang exact / de bord des tenseurs de multiplication matricielle `⟨m, p, n⟩`.

## Arborescence du dossier

```text
axe3_algebraic_complexity/
├── README.md          # Carte et commandes d'exécution de l'Axe III
├── RAPPORT.md         # Monographie autonome (Expériences E16–E23)
├── src/               # 20 scripts Python/Numba/PyTorch (Parseval-LM, Recuit L2/3, Graphe de Flips)
├── results/           # 15 fichiers de résultats JSON
├── witnesses/         # 31 témoins algébriques (.npz)
└── figures/           # 8 figures (.png)
```

## Carte du dossier

| Script (`src/`) | Rôle | Données & Témoins `.npz` / Figures |
| :--- | :--- | :--- |
| [vp_vnp_perm3_border.py](src/vp_vnp_perm3_border.py), [vp_vnp_lm_gct.py](src/vp_vnp_lm_gct.py), [vp_vnp_parseval_lm.py](src/vp_vnp_parseval_lm.py) | Diagnostic Adam/L-BFGS et LM gaussien puis solveur Levenberg–Marquardt de Parseval (`perm₃`, `k = 4..7`, redécouvre `k = 7` à `8,01 × 10⁻¹³`) | [vp_vnp_perm3_results.json](results/vp_vnp_perm3_results.json), [vp_vnp_perm3_border.png](figures/vp_vnp_perm3_border.png) |
| [vp_vnp_k6_deep_probe.py](src/vp_vnp_k6_deep_probe.py) | Sonde `diag_r` (`k = 7` vs `k = 6`), homotopie `7 → 6` et extraction des 165 monômes | [vp_vnp_k6_deep_results.json](results/vp_vnp_k6_deep_results.json), [vp_vnp_k6_homotopy.png](figures/vp_vnp_k6_homotopy.png) |
| [matmul3_tensor_lab.py](src/matmul3_tensor_lab.py) | Phase 1 : Gauss–Newton hermitien exact libre & ℤ₃-cyclique (`⟨2,2,2⟩`, `⟨3,3,3⟩`) | [matmul3_tensor_results.json](results/matmul3_tensor_results.json), [best_Cal_3x3_R_exact_R23.npz](witnesses/best_Cal_3x3_R_exact_R23.npz), [matmul3_tensor_hunt.png](figures/matmul3_tensor_hunt.png) |
| [matmul3_annealed_border_hunt.py](src/matmul3_annealed_border_hunt.py) | Phase 2 : Recuit L<sub>2/3</sub> parallèle + loi d'échelle empirique `α` (`R = 19..22`, bord 21 de Schönhage retrouvé ; la meilleure borne connue reste 20 par Smirnov, 2013) | [matmul3_phase2_results.json](results/matmul3_phase2_results.json), [phase2_Anneal_R_libre_R21.npz](witnesses/phase2_Anneal_R_libre_R21.npz), [matmul3_phase2_annealed.png](figures/matmul3_phase2_annealed.png) |
| [matmul_rect_complex_hunt.py](src/matmul_rect_complex_hunt.py) | Phase 3 : Formats rectangulaires sur ℂ (`⟨2,3,4⟩`, `⟨2,3,5⟩`, `⟨2,4,4⟩`, `⟨3,3,4⟩`) | [matmul_rect_complex_results.json](results/matmul_rect_complex_results.json), [matmul_rect_complex_hunt.png](figures/matmul_rect_complex_hunt.png) |
| [matmul_barrier_exact_hunt.py](src/matmul_barrier_exact_hunt.py) | Phase 4 : Barrière L<sub>∞</sub> circulaire stricte permanente + équilibrage L<sub>∞</sub> | [matmul_barrier_exact_results.json](results/matmul_barrier_exact_results.json), [matmul_barrier_exact_hunt.png](figures/matmul_barrier_exact_hunt.png) |
| [flip_graph_fp.py](src/flip_graph_fp.py) | Phase 5 : Graphe de Flips combinatoire exact + réduction de tranches + PLUS (𝔽₃, 𝔽₅, Numba 10 cœurs) | [flip_graph_fp_results.json](results/flip_graph_fp_results.json), [flip_3x3x3_F3_R23.npz](witnesses/flip_3x3x3_F3_R23.npz), [flip_graph_fp_hunt.png](figures/flip_graph_fp_hunt.png) |
| [flip_graph_v2_guided.py](src/flip_graph_v2_guided.py) | Phase 6 : Flip-Graph guidé par parcimonie + warm-start exact Strassen / Laderman | [flip_graph_v2_results.json](results/flip_graph_v2_results.json), [flip_graph_v2_guided.png](figures/flip_graph_v2_guided.png) |
| [flip_graph_f2_bitvec.py](src/flip_graph_f2_bitvec.py), [flip_graph_f2_adaptive.py](src/flip_graph_f2_adaptive.py), [flip_graph_f2_fast_inc.py](src/flip_graph_f2_fast_inc.py), [flip_graph_f2_nondestructive.py](src/flip_graph_f2_nondestructive.py), [flip_graph_f2_beam.py](src/flip_graph_f2_beam.py), [flip_graph_f2_world_hunt.py](src/flip_graph_f2_world_hunt.py), [flip_graph_f2_world_hunt2.py](src/flip_graph_f2_world_hunt2.py), [flip_graph_f2_446_hunt.py](src/flip_graph_f2_446_hunt.py), [flip_graph_f2_sa_hunt.py](src/flip_graph_f2_sa_hunt.py), [flip_graph_f2_447_448.py](src/flip_graph_f2_447_448.py) | Phase 7 : Graphe de Flips bit-vectoriel `uint64` sur 𝔽₂ (16,7 M flips/s, réduction de Gauss non destructive, Laderman `R = 23` en 0,28 s, reproduction de l'état de l'art connu sur 16 formats jusqu'à `⟨4, 4, 8⟩ = 94` sans nouveau record) | [flip_graph_f2_bitvec_bench.json](results/flip_graph_f2_bitvec_bench.json), [flip_graph_f2_beam_records.json](results/flip_graph_f2_beam_records.json), [flip_graph_f2_world_hunt.json](results/flip_graph_f2_world_hunt.json), [flip_graph_f2_world_hunt2.json](results/flip_graph_f2_world_hunt2.json), [flip_graph_f2_446_hunt.json](results/flip_graph_f2_446_hunt.json), [flip_graph_f2_sa_hunt.json](results/flip_graph_f2_sa_hunt.json), [flip_graph_f2_447_448.json](results/flip_graph_f2_447_448.json), `witnesses/f2_best_*.npz` |

## Reproduire les calculs

Dépendances : Python 3.14, `numpy`, `numba`, `torch`, `matplotlib`. Exemples d'exécution depuis la racine de cet axe (`axe3_algebraic_complexity/`) :

```powershell
python -u src/vp_vnp_parseval_lm.py
python -u src/flip_graph_fp.py
python -u src/flip_graph_f2_bitvec.py
```
