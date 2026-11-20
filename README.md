# Exploration computationnelle & théorique autour de $P \stackrel{?}{=} NP$ et $VP \stackrel{?}{=} VNP$

Ce répertoire regroupe l'ensemble des scripts, solveurs exacts, résultats JSON et témoins algébriques (`.npz`) produits au cours des trois arcs de recherche.

## Rapports scientifiques associés

1. **Arc 1 — Paradigme de Williams, $\mathrm{THR} \circ \mathrm{THR}\text{-SAT}$, Atlas $n=4$, MCSP, CAPP & Topologie de Huang**  
   [rapport_recherche_thr2_sat.md](file:///C:/Users/Axoul/.gemini/antigravity/brain/52cfbf17-f603-488a-9fbc-24b50a9ed2c4/rapport_recherche_thr2_sat.md)
2. **Arc 2 — Représentation faible de $\mathrm{OR}_n$ par polynômes modulo 6 (Théorèmes M1 & M2, percées non symétriques)**  
   [rapport_polynomes_mod6.md](file:///C:/Users/Axoul/.gemini/antigravity/brain/52cfbf17-f603-488a-9fbc-24b50a9ed2c4/rapport_polynomes_mod6.md)
3. **Arc 3 — Complexité algébrique : GCT ($\mathrm{perm}_3$) & Multiplication de matrices (Phases 1 à 6 : continu, recuit $L_{2/3}$, barrière $L_\infty$, Graphe de Flips $\mathbb{F}_3/\mathbb{F}_5$)**  
   [rapport_complexite_algebrique.md](file:///C:/Users/Axoul/.gemini/antigravity/brain/52cfbf17-f603-488a-9fbc-24b50a9ed2c4/rapport_complexite_algebrique.md)

---

## Carte des fichiers

### Arc 1 : $\mathrm{THR} \circ \mathrm{THR}\text{-SAT}$, Atlas $n=4$, MCSP, CAPP & Huang
| Script | Rôle | Données / Sorties |
| :--- | :--- | :--- |
| [williams_demo.py](file:///C:/Users/Axoul/.gemini/antigravity/scratch/pnp-williams/01_thr2_williams/williams_demo.py) | Transformation zêta $\mathrm{SYM} \circ \mathrm{AND}$ ($\times 331{,}4$ à $n=18$) | `williams_demo.png` |
| [polymethod_mq.py](file:///C:/Users/Axoul/.gemini/antigravity/scratch/pnp-williams/01_thr2_williams/polymethod_mq.py) | Méthode polynomiale LPTWY sur $\mathrm{MQ}/\mathbb{F}_2$ (60/60 exact) | `polymethod.png` |
| [thr2_lab.py](file:///C:/Users/Axoul/.gemini/antigravity/scratch/pnp-williams/01_thr2_williams/thr2_lab.py), [hybrid_thr2.py](file:///C:/Users/Axoul/.gemini/antigravity/scratch/pnp-williams/01_thr2_williams/hybrid_thr2.py), [adversarial_suite.py](file:///C:/Users/Axoul/.gemini/antigravity/scratch/pnp-williams/01_thr2_williams/adversarial_suite.py) | Solveur B&B évolutif, suite adversariale et Hybride A+B | `thr2_results.json`, `hybrid_results.json` |
| [theorem_a_star.py](file:///C:/Users/Axoul/.gemini/antigravity/scratch/pnp-williams/01_thr2_williams/theorem_a_star.py), [step2_sliding_gray.py](file:///C:/Users/Axoul/.gemini/antigravity/scratch/pnp-williams/01_thr2_williams/step2_sliding_gray.py), [unified_v5_and_poly_top.py](file:///C:/Users/Axoul/.gemini/antigravity/scratch/pnp-williams/01_thr2_williams/unified_v5_and_poly_top.py) | Théorèmes 1 & 2 : réseau $\mathrm{PGCD}$ + Gray glissant + Solveur unifié v5 ($n=28$ en $0{,}5\text{ s}$) | `theorem_astar_results.json`, `step2_results.json`, `phase6_v5_results.json` |
| [exact_thr2_atlas.py](file:///C:/Users/Axoul/.gemini/antigravity/scratch/pnp-williams/01_thr2_williams/exact_thr2_atlas.py), [mcsp_metacomplexity.py](file:///C:/Users/Axoul/.gemini/antigravity/scratch/pnp-williams/01_thr2_williams/mcsp_metacomplexity.py) | Atlas exhaustif $n=4$ ($65\,536$ fonctions / 222 classes NPN) et échec des preuves naturelles $L_1$ | `atlas_n4.npz`, `mcsp_results.json` |
| [rank_k_bilinear_thr2.py](file:///C:/Users/Axoul/.gemini/antigravity/scratch/pnp-williams/01_thr2_williams/rank_k_bilinear_thr2.py), [piste_a_svd_capp.py](file:///C:/Users/Axoul/.gemini/antigravity/scratch/pnp-williams/01_thr2_williams/piste_a_svd_capp.py), [det_capp_nystrom.py](file:///C:/Users/Axoul/.gemini/antigravity/scratch/pnp-williams/01_thr2_williams/det_capp_nystrom.py), [piste_b_huang_topology.py](file:///C:/Users/Axoul/.gemini/antigravity/scratch/pnp-williams/01_thr2_williams/piste_b_huang_topology.py) | Théorème 3 (factorisation rang $K$ / MIPS), CAPP déterministe Nyström ($n=26$), et topologie de Huang | `theorem3_rank_results.json`, `det_capp_results.json`, `piste_b_huang_results.json` |

### Arc 2 : Polynômes modulo 6 ($\mathrm{OR}_n$)
| Script | Rôle | Données / Sorties |
| :--- | :--- | :--- |
| [mod6_exact_complete.py](file:///C:/Users/Axoul/.gemini/antigravity/scratch/pnp-williams/02_mod6_polynomials/mod6_exact_complete.py) | Table symétrique exacte $N_{\text{sym}}(1..10)$ par canonisation $\mathbb{P}^{d-1}(\mathbb{F}_3)$ | [mod6_exact_complete.json](file:///C:/Users/Axoul/.gemini/antigravity/scratch/pnp-williams/results/mod6_polynomials/mod6_exact_complete.json) |
| [verify_mod6_theorems.py](file:///C:/Users/Axoul/.gemini/antigravity/scratch/pnp-williams/02_mod6_polynomials/verify_mod6_theorems.py) | Preuve & vérification exhaustive des Théorèmes M1 ($N_{\text{gen}}(2)=10$) et M2 (Margin-Harvesting $d=6,7,8,9$) | [mod6_theorems_summary.json](file:///C:/Users/Axoul/.gemini/antigravity/scratch/pnp-williams/results/mod6_polynomials/mod6_theorems_summary.json) |
| [hunt_super_quadratic_mod6.py](file:///C:/Users/Axoul/.gemini/antigravity/scratch/pnp-williams/02_mod6_polynomials/hunt_super_quadratic_mod6.py) | Chasse multi-blocs à $d = 3, 4, 5$ (mur `defect = 1` à $N_{\text{sym}}(d)+1$) | [super_quad_mod6_results.json](file:///C:/Users/Axoul/.gemini/antigravity/scratch/pnp-williams/results/mod6_polynomials/super_quad_mod6_results.json) |

### Arc 3 : Complexité algébrique ($VP \stackrel{?}{=} VNP$ & Multiplication matricielle)
| Script | Rôle | Données & Témoins `.npz` |
| :--- | :--- | :--- |
| [vp_vnp_parseval_lm.py](file:///C:/Users/Axoul/.gemini/antigravity/scratch/pnp-williams/03_algebraic_complexity/vp_vnp_parseval_lm.py) | Solveur Levenberg–Marquardt de Parseval ($\mathrm{perm}_3$, $k=4..7$, redécouvre $k=7$ à $8 \times 10^{-13}$) | [vp_vnp_perm3_results.json](file:///C:/Users/Axoul/.gemini/antigravity/scratch/pnp-williams/results/algebraic_complexity/vp_vnp_perm3_results.json) |
| [vp_vnp_k6_deep_probe.py](file:///C:/Users/Axoul/.gemini/antigravity/scratch/pnp-williams/03_algebraic_complexity/vp_vnp_k6_deep_probe.py) | Sonde `diag_r` ($k=7$ vs $k=6$), homotopie $7 \to 6$ et extraction des 165 monômes | [vp_vnp_k6_deep_results.json](file:///C:/Users/Axoul/.gemini/antigravity/scratch/pnp-williams/results/algebraic_complexity/vp_vnp_k6_deep_results.json) |
| [matmul3_tensor_lab.py](file:///C:/Users/Axoul/.gemini/antigravity/scratch/pnp-williams/03_algebraic_complexity/matmul3_tensor_lab.py) | Phase 1 : Gauss–Newton hermitien exact libre & $\mathbb{Z}_3$-cyclique ($\langle 2,2,2 \rangle$, $\langle 3,3,3 \rangle$) | [matmul3_tensor_results.json](file:///C:/Users/Axoul/.gemini/antigravity/scratch/pnp-williams/results/algebraic_complexity/matmul3_tensor_results.json), `best_Cal_3x3_R_exact_R23.npz` |
| [matmul3_annealed_border_hunt.py](file:///C:/Users/Axoul/.gemini/antigravity/scratch/pnp-williams/03_algebraic_complexity/matmul3_annealed_border_hunt.py) | Phase 2 : Recuit $L_{2/3}$ parallèle + détecteur de pente $\alpha$ ($R=19..22$, bord 21 retrouvé) | [matmul3_phase2_results.json](file:///C:/Users/Axoul/.gemini/antigravity/scratch/pnp-williams/results/algebraic_complexity/matmul3_phase2_results.json), `phase2_Anneal_R_libre_R21.npz` |
| [matmul_rect_complex_hunt.py](file:///C:/Users/Axoul/.gemini/antigravity/scratch/pnp-williams/03_algebraic_complexity/matmul_rect_complex_hunt.py) | Phase 3 : Formats rectangulaires sur $\mathbb{C}$ ($\langle 2,3,4 \rangle, \langle 2,3,5 \rangle, \langle 2,4,4 \rangle, \langle 3,3,4 \rangle$) | [matmul_rect_complex_results.json](file:///C:/Users/Axoul/.gemini/antigravity/scratch/pnp-williams/results/algebraic_complexity/matmul_rect_complex_results.json), `rect_Hunt_*_R*.npz` |
| [matmul_barrier_exact_hunt.py](file:///C:/Users/Axoul/.gemini/antigravity/scratch/pnp-williams/03_algebraic_complexity/matmul_barrier_exact_hunt.py) | Phase 4 : Barrière $L_\infty$ circulaire stricte permanente + équilibrage $L_\infty$ | [matmul_barrier_exact_results.json](file:///C:/Users/Axoul/.gemini/antigravity/scratch/pnp-williams/results/algebraic_complexity/matmul_barrier_exact_results.json), `barrier_*.npz` |
| [flip_graph_fp.py](file:///C:/Users/Axoul/.gemini/antigravity/scratch/pnp-williams/03_algebraic_complexity/flip_graph_fp.py) | Phase 5 : Graphe de Flips combinatoire exact + réduction de tranches + PLUS ($\mathbb{F}_3, \mathbb{F}_5$, Numba 10 cœurs) | [flip_graph_fp_results.json](file:///C:/Users/Axoul/.gemini/antigravity/scratch/pnp-williams/results/algebraic_complexity/flip_graph_fp_results.json), `flip_3x3x3_F3_R23.npz` |
| [flip_graph_v2_guided.py](file:///C:/Users/Axoul/.gemini/antigravity/scratch/pnp-williams/03_algebraic_complexity/flip_graph_v2_guided.py) | Phase 6 : Flip-Graph guidé par parcimonie + warm-start exact Strassen / Laderman | [flip_graph_v2_results.json](file:///C:/Users/Axoul/.gemini/antigravity/scratch/pnp-williams/results/algebraic_complexity/flip_graph_v2_results.json), `flipv2_3x3x3_F*.npz` |

## Reproduire les calculs

Environnement : Python 3.14, `numpy`, `numba`, `torch`, `matplotlib` (CPU multi-cœurs). Chaque script accepte optionnellement le répertoire de sortie des graphiques en premier argument :

```powershell
python -u vp_vnp_parseval_lm.py
python -u matmul3_annealed_border_hunt.py
python -u flip_graph_fp.py
```

