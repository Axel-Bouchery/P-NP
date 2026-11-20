# Exploration théorique et expérimentale des frontières de la complexité (P vs NP, ACC⁰, THR ∘ THR, VP vs VNP)

Ce dépôt rassemble l'ensemble des **46 scripts** Python/Numba/PyTorch, des **37 fichiers de résultats JSON**, des **22 témoins algébriques (`.npz`)** et des **28 figures** produits au cours de la campagne de recherche sur trois axes complémentaires.

## Monographie de recherche unifiée

- **[RAPPORT.md](RAPPORT.md)** — Rapport complet de recherche (problématique, définitions formelles, protocoles, démonstrations des Théorèmes 1–3, M1–M2, G1–G5, G3++, T1–T2, Q1, P1 et U1, tableaux, 28 figures, bilan des cinq murs structurels et bibliographie unifiée).

---

## Carte du dépôt

### Axe I : Paradigme de Williams, THR ∘ THR-SAT, Atlas n = 4, MCSP, CAPP & Topologie de Huang
| Script | Rôle | Données / Figures |
| :--- | :--- | :--- |
| [williams_demo.py](01_thr2_williams/williams_demo.py) | Transformation zêta SYM ∘ AND (× 331,4 à n = 18) | [williams_demo.png](figures/williams_demo.png) |
| [polymethod_mq.py](01_thr2_williams/polymethod_mq.py) | Méthode polynomiale LPTWY sur MQ / 𝔽₂ (60/60 exact) | [polymethod.png](figures/polymethod.png) |
| [thr2_lab.py](01_thr2_williams/thr2_lab.py), [hybrid_thr2.py](01_thr2_williams/hybrid_thr2.py), [adversarial_suite.py](01_thr2_williams/adversarial_suite.py) | Solveur B&B évolutif, suite adversariale et Hybride A+B | [thr2_results.json](results/thr2_williams/thr2_results.json), [hybrid_results.json](results/thr2_williams/hybrid_results.json) |
| [theorem_a_star.py](01_thr2_williams/theorem_a_star.py), [step2_sliding_gray.py](01_thr2_williams/step2_sliding_gray.py), [unified_v5_and_poly_top.py](01_thr2_williams/unified_v5_and_poly_top.py) | Théorèmes 1 & 2 : réseau PGCD + Gray glissant + Solveur unifié v5 (n = 28 en 0,502 s) | [theorem_astar_results.json](results/thr2_williams/theorem_astar_results.json), [step2_results.json](results/thr2_williams/step2_results.json), [phase6_v5_results.json](results/thr2_williams/phase6_v5_results.json) |
| [exact_thr2_atlas.py](01_thr2_williams/exact_thr2_atlas.py), [mcsp_metacomplexity.py](01_thr2_williams/mcsp_metacomplexity.py) | Atlas exhaustif n = 4 (65 536 fonctions / 222 classes NPN) et échec des preuves naturelles L₁ | `witnesses/thr2_williams/atlas_n4.npz`, [mcsp_results.json](results/thr2_williams/mcsp_results.json) |
| [rank_k_bilinear_thr2.py](01_thr2_williams/rank_k_bilinear_thr2.py), [piste_a_svd_capp.py](01_thr2_williams/piste_a_svd_capp.py), [det_capp_nystrom.py](01_thr2_williams/det_capp_nystrom.py), [piste_b_huang_topology.py](01_thr2_williams/piste_b_huang_topology.py) | Théorème 3 (factorisation rang K / MIPS), CAPP déterministe Nyström (n = 26) et topologie de Huang | [theorem3_rank_results.json](results/thr2_williams/theorem3_rank_results.json), [det_capp_results.json](results/thr2_williams/det_capp_results.json), [piste_b_huang_results.json](results/thr2_williams/piste_b_huang_results.json) |

### Axe II : Polynômes faibles modulo m composite (ORₙ, THR<sub>n, t</sub>, MAJₙ, MOD<sub>k</sub>)
| Script | Rôle | Données / Figures |
| :--- | :--- | :--- |
| [mod6_exact_complete.py](02_mod6_polynomials/mod6_exact_complete.py) | Table symétrique exacte N<sub>sym</sub>(1..10) par canonisation ℙ<sup>d−1</sup>(𝔽₃) | [mod6_exact_complete.json](results/mod6_polynomials/mod6_exact_complete.json) |
| [verify_mod6_theorems.py](02_mod6_polynomials/verify_mod6_theorems.py) | Preuve & vérification exhaustive des Théorèmes M1 (N<sub>gen</sub>(2) = 10) et M2 (*Margin-Harvesting* d = 6..9) | [mod6_theorems_summary.json](results/mod6_polynomials/mod6_theorems_summary.json), [mod6_discovery.png](figures/mod6_discovery.png) |
| [hunt_super_quadratic_mod6.py](02_mod6_polynomials/hunt_super_quadratic_mod6.py) | Chasse multi-blocs à d = 3, 4, 5 (mur `defect = 1` à N<sub>sym</sub>(d) + 1) | [super_quad_mod6_results.json](results/mod6_polynomials/super_quad_mod6_results.json) |
| [composite_mod_m_lab.py](02_mod6_polynomials/composite_mod_m_lab.py), [composite_mod_m_full.py](02_mod6_polynomials/composite_mod_m_full.py), [composite_mod_m_verify_and_plot.py](02_mod6_polynomials/composite_mod_m_verify_and_plot.py) | Théorèmes G1–G5 sur 8 modules composites (`mod 6, 10, 12, 14, 15, 21, 30, 35`), record `mod 6 d = 8 → 88`, séparation non bornée `(p_min − 1)×` | [composite_mod_m_full.json](results/mod6_polynomials/composite_mod_m_full.json), [composite_mod_m_discovery.png](figures/composite_mod_m_discovery.png) |
| [option_a_deep_probe.py](02_mod6_polynomials/option_a_deep_probe.py), [option_a_thr2_breakthrough.py](02_mod6_polynomials/option_a_thr2_breakthrough.py) | Optimalité 2b/3b N<sub>gen</sub><sup>(m)</sup>(2) = 2(m − 1), Théorème T1, percée THR<sub>n, 2</sub>, MAJₙ & MOD<sub>k</sub> | [option_a_deep_probe.json](results/mod6_polynomials/option_a_deep_probe.json), [option_a_thr2_breakthrough.json](results/mod6_polynomials/option_a_thr2_breakthrough.json), [option_a_discovery.png](figures/option_a_discovery.png) |
| [piste_2a_thr2_graphs.py](02_mod6_polynomials/piste_2a_thr2_graphs.py), [piste_2a_k10_and_theorem_t2.py](02_mod6_polynomials/piste_2a_k10_and_theorem_t2.py) | Piste 2A — Percée THR<sub>9, 2</sub> mod 6 (K₉ vérifié sur 2⁹ = 512, 3× le symétrique), Théorème T2, records 3-blocs `mod 10, 12, 14, 15` | [piste_2a_thr2_graphs.json](results/mod6_polynomials/piste_2a_thr2_graphs.json), [piste_2a_phase2.json](results/mod6_polynomials/piste_2a_phase2.json), [piste_2a_thr2_breakthrough.png](figures/piste_2a_thr2_breakthrough.png) |
| [piste_2c_2d_mva_and_prob.py](02_mod6_polynomials/piste_2c_2d_mva_and_prob.py) | Pistes 2C & 2D — Théorème G3++ (A<sub>p</sub>(d) = (p − 2)·d exact sur 𝔽₅/𝔽₇) et Théorème P1 (effondrement probabiliste symétrique ε = 100 % vs m<sup>−d</sup>) | [piste_2c_2d_results.json](results/mod6_polynomials/piste_2c_2d_results.json), [piste_2c_2d_mva_prob.png](figures/piste_2c_2d_mva_prob.png) |
| [piste_2b_and_k10_exact.py](02_mod6_polynomials/piste_2b_and_k10_exact.py) | Piste 2B & Front 3 — Théorème Q1 (Ramsey `α(G₂) ≤ 4`), résistance 4b/5b OR₁₉ mod 10, clique `ω(G₂)` et 236k extensions K₉ → K₁₀ | [piste_2b_k10_results.json](results/mod6_polynomials/piste_2b_k10_results.json), [piste_2b_k10_analysis.png](figures/piste_2b_k10_analysis.png) |
| [thr2_mod6_exact_closure.py](02_mod6_polynomials/thr2_mod6_exact_closure.py), [thr2_mod6_audit_and_plot.py](02_mod6_polynomials/thr2_mod6_audit_and_plot.py), [thr2_mod6_n10_probe2.py](02_mod6_polynomials/thr2_mod6_n10_probe2.py), [thr2_mod6_n10_attack.py](02_mod6_polynomials/thr2_mod6_n10_attack.py) | **Théorème U1** — Fermeture exhaustive exacte **N<sub>gen</sub><sup>(6)</sup>(2, THR₂) = 9** (unique à S₉ × {±1} près, `V = {u} ⊔ {v} ⊔ K₂ ⊔ C₅`, `∣Aut∣ = 20`, `ω(G₂) ≤ 6`, 0 sol. à n = 10 auditée bit-à-bit sur 4,78 M tests) | [thr2_mod6_exact_closure.json](results/mod6_polynomials/thr2_mod6_exact_closure.json), [thr2_mod6_audit.json](results/mod6_polynomials/thr2_mod6_audit.json), [thr2_mod6_exact_closure.png](figures/thr2_mod6_exact_closure.png) |

### Axe III : Complexité algébrique (GCT sur perm₃ & Multiplication matricielle)
| Script | Rôle | Données & Témoins `.npz` |
| :--- | :--- | :--- |
| [vp_vnp_parseval_lm.py](03_algebraic_complexity/vp_vnp_parseval_lm.py) | Solveur Levenberg–Marquardt de Parseval (perm₃, k = 4..7, redécouvre k = 7 à 8,01 × 10⁻¹³) | [vp_vnp_perm3_results.json](results/algebraic_complexity/vp_vnp_perm3_results.json) |
| [vp_vnp_k6_deep_probe.py](03_algebraic_complexity/vp_vnp_k6_deep_probe.py) | Sonde `diag_r` (k = 7 vs k = 6), homotopie 7 → 6 et extraction des 165 monômes | [vp_vnp_k6_deep_results.json](results/algebraic_complexity/vp_vnp_k6_deep_results.json) |
| [matmul3_tensor_lab.py](03_algebraic_complexity/matmul3_tensor_lab.py) | Phase 1 : Gauss–Newton hermitien exact libre & ℤ₃-cyclique (⟨2,2,2⟩, ⟨3,3,3⟩) | [matmul3_tensor_results.json](results/algebraic_complexity/matmul3_tensor_results.json), `witnesses/algebraic_complexity/best_Cal_3x3_R_exact_R23.npz` |
| [matmul3_annealed_border_hunt.py](03_algebraic_complexity/matmul3_annealed_border_hunt.py) | Phase 2 : Recuit L<sub>2/3</sub> parallèle + loi de pente α (R = 19..22, bord 21 retrouvé) | [matmul3_phase2_results.json](results/algebraic_complexity/matmul3_phase2_results.json), `witnesses/algebraic_complexity/phase2_Anneal_R_libre_R21.npz` |
| [matmul_rect_complex_hunt.py](03_algebraic_complexity/matmul_rect_complex_hunt.py) | Phase 3 : Formats rectangulaires sur ℂ (⟨2,3,4⟩, ⟨2,3,5⟩, ⟨2,4,4⟩, ⟨3,3,4⟩) | [matmul_rect_complex_results.json](results/algebraic_complexity/matmul_rect_complex_results.json) |
| [matmul_barrier_exact_hunt.py](03_algebraic_complexity/matmul_barrier_exact_hunt.py) | Phase 4 : Barrière L<sub>∞</sub> circulaire stricte permanente + équilibrage L<sub>∞</sub> | [matmul_barrier_exact_results.json](results/algebraic_complexity/matmul_barrier_exact_results.json) |
| [flip_graph_fp.py](03_algebraic_complexity/flip_graph_fp.py) | Phase 5 : Graphe de Flips combinatoire exact + réduction de tranches + PLUS (𝔽₃, 𝔽₅, Numba 10 cœurs) | [flip_graph_fp_results.json](results/algebraic_complexity/flip_graph_fp_results.json), `witnesses/algebraic_complexity/flip_3x3x3_F3_R23.npz` |
| [flip_graph_v2_guided.py](03_algebraic_complexity/flip_graph_v2_guided.py) | Phase 6 : Flip-Graph guidé par parcimonie + warm-start exact Strassen / Laderman | [flip_graph_v2_results.json](results/algebraic_complexity/flip_graph_v2_results.json) |

## Reproduire les calculs

Environnement : Python 3.14, `numpy`, `numba`, `torch`, `matplotlib` (CPU multi-cœurs). Exemples d'exécution depuis la racine du dépôt :

```powershell
python -u 02_mod6_polynomials/composite_mod_m_full.py
python -u 02_mod6_polynomials/option_a_thr2_breakthrough.py
python -u 03_algebraic_complexity/vp_vnp_parseval_lm.py
python -u 03_algebraic_complexity/flip_graph_fp.py
```
