# RAPPORT DE RECHERCHE — Complexité algébrique : GCT sur le permanent 3 × 3 et tenseurs de multiplication matricielle

**Nature du document** : Monographie de recherche théorique et expérimentale (Axe III)  
**Date** : Novembre 2026  
**Dépôt de reproductibilité** : [README.md](README.md)  
**Environnement matériel et logiciel** : Station de travail Windows 16 cœurs CPU, Python 3.14, NumPy 2.3, Numba 0.66 (JIT), PyTorch 2.14 (`float64` / `complex128`), Matplotlib 3.11.

---

## Résumé (*Abstract*)

Ce rapport présente les résultats d'une campagne de recherche numérique et combinatoire consacrée à la complexité algébrique exacte et de bord du **permanent $3 \times 3$** ($\mathrm{perm}_3$, *Geometric Complexity Theory*) et des **tenseurs de multiplication matricielle** $\langle m, p, n \rangle$ sur $\mathbb{R}$, $\mathbb{C}$, $\mathbb{F}_2$, $\mathbb{F}_3$ et $\mathbb{F}_5$.

Le fil directeur des huit expériences relie directement l'optimisation continue à la combinatoire exacte. Sur $\mathrm{perm}_3$, un solveur Levenberg–Marquardt de Parseval à jauge SVD diagonale redécouvre de zéro une représentation déterminantale exacte $7 \times 7$ à **$8{,}01 \times 10^{-13}$** en $33$ itérations, tout en isolant un blocage du solveur à $k = 6$ ($0/64$ essais, sans prouver l'inexistence d'une représentation de bord) et un piège de jauge $\mathrm{SL}_7$ sur l'homotopie $7 \to 6$. Transposée aux tenseurs de multiplication matricielle, l'approche continue (Gauss–Newton hermitien exact à recuit $L_{2/3}$) met en évidence une loi d'échelle empirique $\mathrm{err}(\mu) \propto \mu^\alpha$ avec $\alpha = p / (2p + 2)$ qui discrimine numériquement rang exact borné ($\alpha = 1{,}00$), bord d'ordre 1 ($\alpha = 3/8$) et bord d'ordre 2 ($\alpha = 3/10$). Elle redécouvre Strassen, Laderman, Schönhage ($R = 21$) et trois séparations strictes $\underline{R} < R_{\mathrm{rec}}$ ($\langle 2,3,5 \rangle$, $\langle 2,4,4 \rangle$, $\langle 3,3,4 \rangle$), mais l'essai d'une barrière $L_{\infty}$ stricte démontre ensuite expérimentalement que les vallées de bord capturent $100\,\%$ des trajectoires dès que $\underline{R} < R$.

Pour s'affranchir par construction de ces vallées de bord, la fin de l'étude bascule en combinatoire exacte dans le Graphe de Flips sur corps finis. L'intégration d'une **réduction exacte de tranches par élimination de Gauss** sur $\mathbb{F}_3$ et $\mathbb{F}_5$ ($2{,}6\text{ M flips/s}$, Strassen en $385$ flips, Laderman $R = 23$ sur $\mathbb{F}_3$ en $7{,}8\text{ s}$) quantifie l'explosion projective $|\mathbb{P}^{N-1}(\mathbb{F}_p)|$, tandis que son transfert dans une architecture bit-vectorielle `uint64` sur $\mathbb{F}_2$ (**$16{,}7\text{ M flips/s}$**, réduction de Gauss *non destructive*) redécouvre Laderman en **$0{,}28\text{ s}$** et reproduit l'état de l'art sur **16 formats** jusqu'à $\langle 4, 4, 8 \rangle = 94$.

---

## Table des matières

1. [Problème d'origine, contexte théorique et définitions formelles](#1-problème-dorigine-contexte-théorique-et-définitions-formelles)
   - [1.1 Complexité déterminantale, GCT et rang tensoriel](#11-complexité-déterminantale-gct-et-rang-tensoriel)
   - [1.2 Définitions formelles : complexité déterminantale, rang tensoriel et Graphe de Flips](#12-définitions-formelles--complexité-déterminantale-rang-tensoriel-et-graphe-de-flips)
2. [Méthodologie, architecture et protocoles de validation](#2-méthodologie-architecture-et-protocoles-de-validation)
3. [Bilan synthétique : résultats établis et obstructions structurelles](#3-bilan-synthétique--résultats-établis-et-obstructions-structurelles)
   - [3.1 Synthèse des résultats établis](#31-synthèse-des-résultats-établis)
   - [3.2 Limites et obstructions structurelles](#32-limites-et-obstructions-structurelles)
4. [Présentation détaillée des expériences, buts, résultats et graphiques (Expériences E16 à E23)](#4-présentation-détaillée-des-expériences-buts-résultats-et-graphiques-expériences-e16-à-e23)
5. [Conclusion générale et perspectives de recherche](#5-conclusion-générale-et-perspectives-de-recherche)
6. [Mentions, reproductibilité et bibliographie](#6-mentions-reproductibilité-et-bibliographie)

---

## 1. Problème d'origine, contexte théorique et définitions formelles

### 1.1 Complexité déterminantale, GCT et rang tensoriel

En complexité algébrique (Valiant, 1979), l'expressivité relative du permanent ($\mathrm{perm}_n$) vis-à-vis du déterminant et le coût bilinéaire de la multiplication matricielle se traduisent par des problèmes géométriques et combinatoires concrets. La géométrie algébrique et la théorie des tenseurs (*Geometric Complexity Theory* et rang tensoriel) étudient notamment la frontière entre complexité **exacte** et complexité **de bord** (fermeture d'orbite) pour le permanent $\mathrm{perm}_3$ et pour les tenseurs de multiplication matricielle rapide.

### 1.2 Définitions formelles : complexité déterminantale, rang tensoriel et Graphe de Flips

Trois objets structurent l'ensemble de cet axe. D'abord, la **complexité déterminantale exacte** $\mathrm{dc}(\mathrm{perm}_n)$ est le plus petit entier $k$ tel que $\mathrm{perm}_n(X) = \det_k(A_0 + L_1(X))$ pour une matrice affine de taille $k \times k$, tandis que la **complexité déterminantale de bord** $\overline{\mathrm{dc}}(\mathrm{perm}_n)$ est le plus petit $k$ tel que $\mathrm{perm}_n$ appartienne à la fermeture d'orbite $\overline{\mathrm{GL}_{k^2} \cdot \det_k}$ (limite de déterminants $k \times k$ lorsque $\varepsilon \to 0$). Pour $n = 3$, il est établi que $\mathrm{dc}(\mathrm{perm}_3) = 7$ (Grenet, 2011 ; Alper–Bogart–Velasco, 2017) et $5 \le \overline{\mathrm{dc}}(\mathrm{perm}_3) \le 7$, la question $\overline{\mathrm{dc}}(\mathrm{perm}_3) \stackrel{?}{=} 6$ restant ouverte.

Ensuite, la multiplication d'une matrice $m \times p$ par une matrice $p \times n$ correspond à un **tenseur trilinéaire** $T_{\langle m, p, n \rangle}$ d'ordre $3$. Son **rang exact** $R(\langle m, p, n \rangle)$ sur un corps $\mathbb{K}$ est le nombre minimal $R$ de produits bilinéaires de rang $1$ dont la somme égale exactement $T_{\langle m, p, n \rangle}$, tandis que son **rang de bord** $\underline{R}(\langle m, p, n \rangle)$ est le plus petit $R$ tel que $T_{\langle m, p, n \rangle}$ soit limite d'une famille $T(\varepsilon)$ de rang $\le R$. Pour $\langle 3, 3, 3 \rangle$, l'algorithme de Laderman (1976) fixe $R \le 23$ (prouvé optimal sur $\mathbb{F}_2$ et $\mathbb{Z}$ en 2026), laissant ouvert $19 \le R \le 23$ sur $\mathbb{Q}/\mathbb{R}/\mathbb{C}$ (Bläser, 2003), alors que $17 \le \underline{R}(\langle 3, 3, 3 \rangle) \le 20$ (Conner–Harper–Landsberg, 2023 ; Smirnov, 2013).

Enfin, le **Graphe de Flips** (Kauers et Moosbauer, 2023) est un graphe combinatoire défini sur les schémas bilinéaires exacts sur un corps fini $\mathbb{F}_p$, dont les transitions locales (*flips*, réductions de rang et transitions PLUS d'Arai et al., 2024) conservent exactement le tenseur tout en recherchant une baisse de rang sans jamais quitter la variété exacte.

---

## 2. Méthodologie, architecture et protocoles de validation

```mermaid
flowchart TD
    Root["Complexite algebrique : GCT perm_3 et MatMul"]
    Root --> A3a["GCT perm_3 : Parseval-LM (k=7 exact, k=6 mur, piege SL_7)"]
    Root --> A3b["MatMul continu : Recuit L2/3, loi alpha, barriere Linf stricte"]
    Root --> A3c["MatMul discret : Graphe de Flips F2/F3/F5 (uint64 16.7M/s) + Gauss non destructif"]
    A3b -.-> M4["Mur 4 : Domination volumetrique des vallees de bord"]
    A3c -.-> M5["Mur 5 : Explosion projective |P^(N-1)(F_p)|"]
```

Les **8 expériences (E16 à E23, 8 figures)** appliquent deux règles de validation strictes. En régime continu, aucun solveur n'est déployé sur une frontière ouverte avant d'avoir redécouvert, depuis une initialisation aléatoire pure, les références exactes à précision machine ($\det_7$ pour $\mathrm{perm}_3$ à $8{,}01 \times 10^{-13}$, Strassen $R = 7$ à $9{,}1 \times 10^{-13}$, Laderman $R = 23$ à $5{,}61 \times 10^{-13}$). En régime discret, chaque schéma tensoriel produit dans le Graphe de Flips ($\sim 400$ millions de flips sur $\mathbb{F}_3$ et $\mathbb{F}_5$, et plus de $500$ millions de flips sur $\mathbb{F}_2$) est vérifié bit-à-bit contre le tenseur exact $T_{\langle m, p, n \rangle}$.

---

## 3. Bilan synthétique : résultats établis et obstructions structurelles

### 3.1 Synthèse des résultats établis

| Résultat | Nature | Énoncé clé / Vérification |
| :--- | :--- | :--- |
| **Parseval-LM à jauge SVD** et diagnostic du piège $\mathrm{SL}_7$ sur $\mathrm{perm}_3$ | **Solveur numérique et diagnostic** | Redécouverte exacte $k = 7$ à **$8{,}01 \times 10^{-13}$** (33 iters) ; mur numérique $k = 6$ ($0/64$ essais, limite du solveur et non preuve d'inexistence) |
| **Loi d'échelle empirique $\alpha = p / (2p + 2)$** sous recuit $L_{2/3}$ + barrière $L_{\infty}$ | **Analyse numérique** | Comportement numérique du solveur continu ($\alpha = 1{,}00 / 0{,}375 / 0{,}300$) ; verrouillage de paroi observé sur Smirnov $R = 29$ |
| **Réduction exacte de tranches mod $p$ et Graphe de Flips `uint64` $\mathbb{F}_2$** | **Algorithme combinatoire exact** | Strassen sur $\mathbb{F}_3$ en **$385$ flips** ; Laderman $R = 23$ sur $\mathbb{F}_3$ en **$7{,}8\text{ s}$** ($2{,}6\text{ M flips/s}$) et sur $\mathbb{F}_2$ en **$0{,}28\text{ s}$** (**$16{,}7\text{ M flips/s}$**, état de l'art reproduit sur **16 formats** jusqu'à $\langle 4, 4, 8 \rangle = 94$) |

### 3.2 Limites et obstructions structurelles

Aucune représentation de bord de $\mathrm{perm}_3$ à $k = 6$ ($\overline{\mathrm{dc}}(\mathrm{perm}_3) \stackrel{?}{=} 6$) ni amélioration des bornes de rang exact en multiplication matricielle (en particulier $R(\langle 3, 3, 3 \rangle) < 23$ sur $\mathbb{R}/\mathbb{C}$) n'a été obtenue. Deux murs structurels complémentaires expliquent cette frontière :

| # | Mur structurel | Contexte | Manifestation quantitative mesurée |
| :---: | :--- | :--- | :--- |
| **4** | **Domination des vallées de bord** | GCT $\mathrm{perm}_3$ & MatMul continu | Dès que $\underline{R} < R$, les vallées de bord occupent presque tout le bassin d'attraction continu : sans barrière elles divergent ($\alpha \approx 0{,}31$ à $0{,}40$) ; sous barrière $L_{\infty}$ stricte, $100\,\%$ des trajectoires se plaquent contre la paroi $|z_q| = M$, masquant même le rang exact connu de Smirnov ($R = 29$). |
| **5** | **Explosion projective en corps fini** | Graphe de Flips $\mathbb{F}_3 / \mathbb{F}_5$ | L'élimination des vallées de bord via $\mathbb{F}_p$ ($p \ge 3$) multiplie le cardinal projectif $|\mathbb{P}^{N-1}(\mathbb{F}_p)|$ de $511$ ($\mathbb{F}_2$, $N = 9$) à $9\,841$ ($\mathbb{F}_3$, $N = 9$), $488\,281$ ($\mathbb{F}_5$, $N = 9$) et $6{,}1 \times 10^7$ ($\mathbb{F}_5$, $N = 12$), raréfiant les collisions créatrices de réductions. |

---

## 4. Présentation détaillée des expériences, buts, résultats et graphiques (Expériences E16 à E23)

#### Expérience E16 — Solveur Levenberg–Marquardt de Parseval sur perm₃ (k = 4..7)

L'exploration commence sur la question ouverte $\overline{\mathrm{dc}}(\mathrm{perm}_3) \stackrel{?}{=} 6$ : avant de sonder $k = 6$ sur $\mathbb{R}$ et sur $\mathbb{C}$, il faut disposer d'un solveur numérique capable de redécouvrir de zéro, sans biais, la représentation déterminantale exacte $7 \times 7$ de Grenet (2011). Les premiers essais par optimisation générique montrent vite pourquoi la formulation naïve échoue : Adam/L-BFGS ([vp_vnp_perm3_border.py](src/vp_vnp_perm3_border.py), $\mathrm{err}_{\mathrm{val}} = 0{,}937$) comme Levenberg–Marquardt gaussien ([vp_vnp_lm_gct.py](src/vp_vnp_lm_gct.py), $\mathrm{err}_{\mathrm{val}} = 25{,}3$) sont submergés par l'explosion des termes homogènes parasites de degré jusqu'à $7$.

Pour lever ce verrou, le solveur original [vp_vnp_parseval_lm.py](src/vp_vnp_parseval_lm.py) évalue $\det_k(A_0 + t_j L_1(X))$ aux $(k + 1)$ racines de l'unité $t_j = e^{2i\pi j/(k+1)}$, ce qui découple orthogonalement les composantes homogènes $P_0, \dots, P_k$ par identité de Parseval. En combinant ce découplage avec un échantillonnage sphérique ($\|X\|_F = 2{,}2$), une fixation de jauge $\mathrm{SL}_k \times \mathrm{SL}_k$ diagonalisant $A_0$ par SVD et le calcul analytique de la Jacobienne par cofacteurs vectorisés $G^{(b)} = \det(L^{(b)})(L^{(b)})^{-T}$, le solveur réussit la calibration $k = 7$ en **$33$ itérations** à **$8{,}01 \times 10^{-13}$** sur les $400$ points de validation. En revanche, dès que $k$ descend à $6$ (ou $5, 4$), l'erreur plafonne nettement au-dessus de $0{,}56$ aussi bien sur $\mathbb{R}$ que sur $\mathbb{C}$ :

| Corps | $k$ | Essais | Best err. TRAIN | Best err. VALID (400 pts) | Médiane VALID | $\|A\|_F$ | Degrés parasites ($d \ne 3$) |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| $\mathbb{R}$ | **7** (cal.) | $8$ | $5{,}55 \times 10^{-13}$ | **$8{,}01 \times 10^{-13}$** | $0{,}704$ | $5{,}17$ | Tous $\le 7{,}4 \times 10^{-13}$ (**$33$ iters**) |
| $\mathbb{R}$ | $4$ | $4$ | $0{,}444$ | $0{,}743$ | $0{,}786$ | $5{,}80$ | $d_2 = 0{,}14$, $d_4 = 0{,}19$ |
| $\mathbb{R}$ | $5$ | $6$ | $0{,}407$ | $0{,}781$ | $0{,}870$ | $5{,}90$ | Tous $< 3 \times 10^{-9}$ |
| $\mathbb{C}$ | $5$ | $6$ | $0{,}391$ | $0{,}731$ | $0{,}842$ | $5{,}49$ | Tous $< 2 \times 10^{-12}$ |
| $\mathbb{R}$ | **6** | $12$ | $0{,}256$ | **$0{,}566$** | $0{,}790$ | $15{,}09$ | $\le 3{,}2 \times 10^{-3}$ |
| $\mathbb{C}$ | **6** | $12$ | $0{,}273$ | **$0{,}640$** | $0{,}779$ | $11{,}28$ | $\le 8{,}4 \times 10^{-4}$ |

![Calibration k=7 et frontière k=6 sur perm₃ par Parseval-LM](figures/vp_vnp_perm3_border.png)

---

#### Expérience E17 — Sonde profonde k = 6, extraction des 165 monômes et piège de jauge SL₇ (7 → 6)

Le contraste observé dans l'Expérience E16 entre $k = 7$ et $k = 6$ appelle deux vérifications complémentaires ([vp_vnp_k6_deep_probe.py](src/vp_vnp_k6_deep_probe.py)) : confirmer statistiquement ce mur en initialisation `diag_r`, et tester si une homotopie pénalisée partant de la solution exacte $7 \times 7$ de Grenet peut déformer continûment $\det_7$ vers une limite de bord $\det_6$.

Sur le plan statistique, le fossé se confirme : **$4/10$ essais** convergent à $k = 7$ ($< 0{,}015$, jusqu'à $4{,}1 \times 10^{-4}$ en $140$ itérations), contre **$0/24$ sur $\mathbb{R}$** et **$0/16$ sur $\mathbb{C}$** à $k = 6$ (soit **$0/64$** au cumul sous le plancher $0{,}566 \approx \sqrt{2/6}$). L'extraction exacte des $165$ monômes cubiques du meilleur candidat `R_k6` montre qu'il ne capture que $(0{,}90;\; 0{,}51;\; 0{,}86;\; 0{,}85;\; 0{,}50;\; 0{,}84)$ sur les $6$ termes de $\mathrm{perm}_3$, tout en laissant une norme $L_2 = 1{,}128$ sur les $159$ monômes parasites (contre $7{,}33 \times 10^{-15}$ pour Grenet $7 \times 7$). Plus révélateur encore, sur l'homotopie pénalisée $V_7(\mathrm{perm}_3) \to \det_6$, la norme de bordure $\|\mathrm{bord}\|$ s'effondre bien de $1{,}73$ à $4{,}7 \times 10^{-3}$ pendant que $\|L_6\|_F$ grimpe à $67{,}4$, mais l'erreur `err6` reste rigoureusement figée à **$0{,}621$ ($\mathbb{R}$)** et **$0{,}607$ ($\mathbb{C}$)** : la trajectoire ne s'approche pas d'une limite de bord, elle se contente d'effectuer un recalibrage diagonal dans $\mathrm{SL}_7$ qui conserve simultanément $\det_7$ et $\det_6$.

![Homotopie 7 vers 6 et calibration statistique k=7 vs k=6](figures/vp_vnp_k6_homotopy.png)

---

#### Expérience E18 — Multiplication matricielle ⟨3, 3, 3⟩ : Gauss–Newton hermitien, recuit L₂/₃ et loi de pente α (Phases 1–2)

Le piège de jauge mis au jour sur $\mathrm{perm}_3$ montre qu'en optimisation continue, la divergence des normes de facteurs peut aussi bien signaler une véritable limite de bord qu'une dégénérescence de jauge stérile. Pour étudier ce phénomène sur un terrain où rang exact et rang de bord sont déjà finement cartographiés, nous transposons la méthode aux tenseurs de multiplication matricielle $\langle 2, 2, 2 \rangle$ et $\langle 3, 3, 3 \rangle$ ([matmul3_tensor_lab.py](src/matmul3_tensor_lab.py), [matmul3_annealed_border_hunt.py](src/matmul3_annealed_border_hunt.py), $8$ cœurs). L'holomorphie de $\widehat{T}(U, V, W)$ permet de résoudre les équations normales exactes en forme hermitienne $P \times P$ avec équilibrage AM-GM — ce qui transforme la pénalité Ridge en quasi-norme $L_{2/3}$ : $(3\mu/2)\sum_r \sigma_r^{2/3}$ — couplé à une projection linéaire $\mathbb{Z}_3$-cyclique.

Minimiser $C \cdot F^{-2p} + \mu \cdot F^2$ lorsque l'erreur décroît comme $\mathrm{err} \sim \|F\|^{-p}$ impose analytiquement une loi d'échelle qui distingue immédiatement la nature de la trajectoire :

$$
\mathrm{err}(\mu) \;\propto\; \mu^{\alpha}, \qquad \alpha \;=\; \frac{p}{2(p+1)} \;=\; \begin{cases} 1{,}000 & \text{(rang exact borné, } p = \infty\text{)}, \\ 0{,}375 & \text{(bord d'ordre 1, } p = 3\text{)}, \\ 0{,}300 & \text{(bord d'ordre 2, } p = 3/2\text{)}. \end{cases}
$$

Les mesures confirment cette signature à deux chiffres significatifs. En calibration exacte, **Strassen $\langle 2,2,2 \rangle$ ($R = 7$)** converge sur $6/6$ essais à **$9{,}1 \times 10^{-13}$** en $0{,}04\text{ s}$ ($4/6$ en mode $\mathbb{Z}_3$-cyclique $1 + 2 \times 3$, avec mur $R = 6$ confirmé à $1/\sqrt{8} = 0{,}3536$), et **Laderman $\langle 3,3,3 \rangle$ ($R = 23$)** converge sur **$6/6$ essais exacts bornés** à **$5{,}61 \times 10^{-13}$** ($\|F\| = 10{,}5$, $\max|x| = 1{,}25$, $\alpha \approx 1{,}00$). Dès que l'on descend à **$R = 22$**, le recuit $L_{2/3}$ franchit le plateau discret $1/\sqrt{27} = 0{,}1925$ sur **$29\,\%$ des essais $\mathbb{R}$** ($5{,}14 \times 10^{-4}$, $\alpha = 0{,}386$), **$56\,\%$ des essais $\mathbb{C}$** ($5{,}37 \times 10^{-4}$, $\alpha = 0{,}398$) et en symétrie cyclique $\mathbb{Z}_3$ ($4 + 6 \times 3$, $6{,}58 \times 10^{-4}$), mais **tous** portent la signature d'un **bord d'ordre 1** ($\alpha \approx 0{,}38$, $\mathrm{err} \propto \|F\|^{-3}$), sans aucun point exact borné. À **$R = 21$ (Schönhage)**, $1$ essai sur $16$ converge et atteint **$5{,}17 \times 10^{-3}$** en polissage ($\|F\| = 93{,}5$, $\alpha = 0{,}311$, $p = -1{,}58 \approx -3/2$, signature d'un **bord d'ordre 2**, [phase2_Anneal_R_libre_R21.npz](witnesses/phase2_Anneal_R_libre_R21.npz)). Plus bas encore, $0/36$ essais à $R = 20$ et $0/188$ à $R = 19$ restent bloqués sur les plateaux discrets $\sqrt{2/27}$ et $\sqrt{3/27}$.

![Phase 1 : Rang exact vs Rang de bord sur 3x3](figures/matmul3_tensor_hunt.png)

![Phase 2 : Recuit L2/3 et spectres de pente sur 3x3](figures/matmul3_phase2_annealed.png)

---

#### Expérience E19 — Chasse sur formats rectangulaires ouverts sur ℂ (Phase 3)

Armé du détecteur de pente $\alpha$ calibré dans l'Expérience E18, nous déployons le solveur hermitien à recuit $L_{2/3}$ sur quatre formats rectangulaires ouverts sur $\mathbb{C}$ ($\langle 2,3,4 \rangle$, $\langle 2,3,5 \rangle$, $\langle 2,4,4 \rangle$, $\langle 3,3,4 \rangle$) afin de repérer automatiquement des séparations exact/bord ou d'éventuelles solutions exactes bornées ([matmul_rect_complex_hunt.py](src/matmul_rect_complex_hunt.py), $126$ essais combinant initialisation gaussienne, initialisation cyclotomique $\{0, \pm 1, \pm i, \omega\}$ et *warm-start* compressé).

La calibration exacte est d'abord retrouvée avec $\alpha = 0{,}98$ et $\max|x| \approx 1{,}0$ à $1{,}1$ sur $\langle 2,3,4 \rangle$ ($R = 20$, $8{,}95 \times 10^{-10}$), $\langle 2,3,5 \rangle$ ($R = 25$, $1{,}45 \times 10^{-10}$) et $\langle 2,4,4 \rangle$ ($R = 26$, $5{,}31 \times 10^{-13}$). Juste en dessous de ces rangs records, le détecteur isole automatiquement **trois séparations strictes $\underline{R} < R_{\mathrm{rec}}$** confirmées par polissage : $\underline{R}(\langle 2,3,5 \rangle) \le 24$ ($3{,}26 \times 10^{-4}$, $\alpha = 0{,}376$, $p = -3{,}01$, cohérent avec Bini $10 + 14$), $\underline{R}(\langle 2,4,4 \rangle) \le 25$ ($6{,}26 \times 10^{-3}$, $\alpha = 0{,}314$, $p = -1{,}64$, bord d'ordre 2 révélé par l'initialisation cyclotomique) et $\underline{R}(\langle 3,3,4 \rangle) \le 28$ ($3{,}72 \times 10^{-4}$, $\alpha = 0{,}378$, $p = -3{,}04$, $5$ témoins indépendants cohérents avec $2 \times 14$). Seul le mur $\langle 2,3,4 \rangle$ à $R = 19$ résiste intégralement ($0/30$ essais sous $1/\sqrt{24} = 0{,}2041$).

![Phase 3 : Chasse rectangulaire sur C](figures/matmul_rect_complex_hunt.png)

---

#### Expérience E20 — Barrière L∞ circulaire stricte permanente et test tournesol ⟨3,3,4⟩ R = 29 (Phase 4)

Les Expériences E18 et E19 laissent ouverte une question décisive : si toutes les trajectoires convergent vers des solutions de bord dès que $\underline{R} < R$, est-ce parce qu'aucune solution exacte n'existe, ou simplement parce que les vallées de bord aspirent le solveur continu ? Pour trancher, nous soumettons le solveur à un **test tournesol** : imposer une paroi $L_{\infty}$ stricte empêchant la divergence des facteurs tout en testant un format où l'on **sait** qu'il existe à la fois une solution exacte et une vallée de bord — en l'occurrence $\langle 3,3,4 \rangle$ à $R = 29$ (record exact connu de Smirnov, alors que $\underline{R} \le 28$).

Dans [matmul_barrier_exact_hunt.py](src/matmul_barrier_exact_hunt.py), la pénalité $\Phi_M(z) = (\gamma/2)\sum_q \max(0, |z_q| - M)^2$ ($\gamma = 100$, majorant MM hermitien exact $\gamma \cdot \mathrm{diag}(\mathbf{1}_{|z_q| > M})$) est couplée à un équilibrage de jauge $L_{\infty}$ exact ($\|u_r\|_{\infty} = \|v_r\|_{\infty} = \|w_r\|_{\infty}$). Sur la calibration $\langle 3,3,3 \rangle$ ($R = 23$), les trajectoires restent bien intérieures ($6/6$ à $\max|z_q| = 0{,}80$ à $0{,}92 < M$, convergeant à $2{,}04 \times 10^{-11}$). En revanche, sur le test tournesol $\langle 3,3,4 \rangle$ ($R = 29$) comme sur les cibles $R_{\mathrm{rec}} - 1$ ($104$ essais), **$100\,\%$ des trajectoires** viennent s'écraser verticalement contre la paroi $|z_q| = M$ ($\mathrm{err} = 1{,}77 \times 10^{-3}$ à $5{,}79 \times 10^{-3}$, $23/23$ collés à la paroi sur $\langle 3,3,3 \rangle$ $R = 22$). Le résultat est sans appel : dès que $\underline{R} < R$, le volume de la vallée de bord se transforme en un minimum local profond plaqué au mur qui masque même les solutions exactes connues (**Mur n° 4**).

![Phase 4 : Barrière L_inf stricte et verrouillage sur la paroi](figures/matmul_barrier_exact_hunt.png)

---

#### Expérience E21 — Graphe de Flips combinatoire exact et réduction de tranches sur 𝔽₃ et 𝔽₅ (Phase 5)

Le verdict de l'Expérience E20 impose un changement de paradigme : pour chasser des schémas exacts sans être piégé par les vallées de bord, il faut quitter $\mathbb{R}$ et $\mathbb{C}$ et travailler en combinatoire exacte sur des corps finis $\mathbb{F}_p$, où aucune limite $\varepsilon \to 0$ n'existe. Nous déployons donc dans [flip_graph_fp.py](src/flip_graph_fp.py) un moteur de Graphe de Flips compilé en Numba JIT sur **$10$ cœurs** sur $\mathbb{F}_3$ et $\mathbb{F}_5$ (contournant l'obstruction $\mathbb{F}_2$, avec $1/2 \in \mathbb{F}_3, \mathbb{F}_5$ et $i = 2 \in \mathbb{F}_5$), avec encodage projectif sans perte sur `int64` ($5^{20} < 2^{47} < 2^{63}$) et $\sim 400$ millions de flips vérifiés exacts à $100\,\%$.

L'apport méthodologique central consiste à intégrer dans la marche aléatoire une **réduction exacte de tranches par élimination de Gauss mod $p$** sur chaque tranche partageant une direction $u$, $v$ ou $w$ ($m \ge 3$ termes). Cet opérateur ramène la redécouverte de Strassen $\langle 2,2,2 \rangle$ sur $\mathbb{F}_3$ de $15\,707$ à **$385$ flips**, et débloque $\mathbb{F}_5$ en **$8\,804$ flips**. Parti du schéma naïf $R = 27$ sur $\langle 3,3,3 \rangle$, le solveur redécouvre **Laderman $R = 23$ exact sur $\mathbb{F}_3$ en $7{,}8\text{ s}$** ($3/10$ workers, cadence **$2{,}6\text{ M flips/s}$**, [flip_3x3x3_F3_R23.npz](witnesses/flip_3x3x3_F3_R23.npz)). Toutefois, dès que la taille du format ou du corps augmente, la taille projective $|\mathbb{P}^{N-1}(\mathbb{F}_p)|$ explose — de $511$ ($\mathbb{F}_2, N = 9$) à $9\,841$ ($\mathbb{F}_3, N = 9$), $488\,281$ ($\mathbb{F}_5, N = 9$) et $6{,}1 \times 10^7$ ($\mathbb{F}_5, N = 12$) — raréfiant les collisions créatrices de réductions (**Mur n° 5**).

![Phase 5 : Graphe de Flips exact non guidé sur F_3 et F_5](figures/flip_graph_fp_hunt.png)

---

#### Expérience E22 — Guidage par parcimonie et warm-start exact dans le Graphe de Flips (Phase 6)

Face à la raréfaction des collisions mesurée sur $\mathbb{F}_3$ et $\mathbb{F}_5$ dans l'Expérience E21, une idée naturelle consiste à biaiser la marche aléatoire vers les schémas les plus creux ($\Delta\mathrm{nnz} \le 0$ à $75\,\%$) et à initialiser la recherche par *warm-start* depuis des schémas exacts déjà connus (Strassen, Laderman). L'expérience [flip_graph_v2_guided.py](src/flip_graph_v2_guided.py) met toutefois en évidence un résultat négatif instructif : loin de franchir le mur projectif, ce biais glouton de parcimonie piège la marche aléatoire à l'intérieur des blocs initiaux sans mélanger les sous-schémas, bloquant précisément les transitions PLUS nécessaires à la découverte de nouvelles factorisations.

![Phase 6 : Impact du guidage par parcimonie et warm-start exact](figures/flip_graph_v2_guided.png)

---

#### Expérience E23 — Graphe de Flips bit-vectoriel `uint64` sur 𝔽₂ à 16,7 M flips/s et réduction de Gauss non destructive sur 16 formats (Phase 7)

Les deux leçons des Expériences E21 et E22 — l'avantage décisif du petit cardinal projectif de $\mathbb{F}_2$ et la nécessité de préserver l'ergodicité de la marche — convergent dans cette dernière phase ([flip_graph_f2_bitvec.py](src/flip_graph_f2_bitvec.py), [flip_graph_f2_beam.py](src/flip_graph_f2_beam.py), [flip_graph_f2_adaptive.py](src/flip_graph_f2_adaptive.py), [flip_graph_f2_fast_inc.py](src/flip_graph_f2_fast_inc.py), [flip_graph_f2_nondestructive.py](src/flip_graph_f2_nondestructive.py), [flip_graph_f2_world_hunt.py](src/flip_graph_f2_world_hunt.py), [flip_graph_f2_world_hunt2.py](src/flip_graph_f2_world_hunt2.py), [flip_graph_f2_446_hunt.py](src/flip_graph_f2_446_hunt.py), [flip_graph_f2_sa_hunt.py](src/flip_graph_f2_sa_hunt.py), [flip_graph_f2_447_448.py](src/flip_graph_f2_447_448.py) ; données dans [flip_graph_f2_bitvec_bench.json](results/flip_graph_f2_bitvec_bench.json), [flip_graph_f2_beam_records.json](results/flip_graph_f2_beam_records.json), [flip_graph_f2_world_hunt.json](results/flip_graph_f2_world_hunt.json), [flip_graph_f2_world_hunt2.json](results/flip_graph_f2_world_hunt2.json), [flip_graph_f2_446_hunt.json](results/flip_graph_f2_446_hunt.json), [flip_graph_f2_sa_hunt.json](results/flip_graph_f2_sa_hunt.json) et [flip_graph_f2_447_448.json](results/flip_graph_f2_447_448.json)).

Sur le plan algorithmique, le transfert sur $\mathbb{F}_2$ en représentation bit-vectorielle `uint64` (opérations `XOR` natives) et le remplacement du balayage quadratique $O(R^2)$ des paires compatibles par une mise à jour incrémentale en $O(R)$ portent la cadence à **$10{,}5$–$16{,}7$ millions de flips/s** sur $12$ cœurs. L'expérience révèle parallèlement un piège dynamique subtil : appliquer l'élimination de Gauss sur les tranches *en place* tous les $8$ flips re-canonise périodiquement les vecteurs et enferme la marche aléatoire dans un cycle court (bloquant $\langle 3, 4, 4 \rangle$ à $R = 39$ et $\langle 4, 4, 4 \rangle$ à $R = 48$). En exécutant la réduction de Gauss tous les $64$ flips sur une **copie tampon non destructive** `(su, sv, sw)` qui n'est adoptée que si elle abaisse strictement le rang ($R_g < R$), la marche conserve toute son ergodicité tout en capturant instantanément les réductions multi-termes.

Ainsi équipé, le moteur reproduit l'état de l'art exact sur **16 formats sur $\mathbb{F}_2$**, chaque témoin étant certifié bit-à-bit contre $T_{\langle m, p, n \rangle}$ : sur les petits et moyens formats, **$\langle 2, 2, 2 \rangle = 7$** ($0{,}00\text{ s}$), **$\langle 2, 2, 3 \rangle = 11$** ($0{,}00\text{ s}$), **$\langle 2, 3, 3 \rangle = 15$** ($0{,}00\text{ s}$), **$\langle 2, 3, 4 \rangle = 20$** ($0{,}05\text{ s}$), **$\langle 2, 3, 5 \rangle = 25$** ($0{,}08\text{ s}$), **$\langle 2, 4, 4 \rangle = 26$** ($0{,}08\text{ s}$), **$\langle 3, 3, 3 \rangle = 23$** (**$0{,}28\text{ s}$**, $28\times$ plus rapide que sur $\mathbb{F}_3$) et **$\langle 3, 3, 4 \rangle = 29$** ($0{,}43\text{ s}$) ; sur les formats rectangulaires et carrés d'ordre $4$ à $8$ (archivés dans `witnesses/`), **$\langle 3, 4, 4 \rangle = 38$**, **$\langle 3, 4, 5 \rangle = 47$** ([f2_best_3x4x5_R47.npz](witnesses/f2_best_3x4x5_R47.npz)), **$\langle 3, 4, 6 \rangle = 56$** ([f2_best_3x4x6_R56.npz](witnesses/f2_best_3x4x6_R56.npz)), **$\langle 3, 4, 7 \rangle = 66$** ([f2_best_3x4x7_R66.npz](witnesses/f2_best_3x4x7_R66.npz)), **$\langle 4, 4, 4 \rangle = 47$** ([f2_best_4x4x4_R47.npz](witnesses/f2_best_4x4x4_R47.npz)), **$\langle 4, 4, 6 \rangle = 73$** ([f2_best_4x4x6_R73.npz](witnesses/f2_best_4x4x6_R73.npz)), **$\langle 4, 4, 7 \rangle = 85$** ([f2_best_4x4x7_R85.npz](witnesses/f2_best_4x4x7_R85.npz)) et **$\langle 4, 4, 8 \rangle = 94$** ([f2_best_4x4x8_R94.npz](witnesses/f2_best_4x4x8_R94.npz)). Sur les paliers profonds $\langle 3, 5, 5 \rangle$ ([f2_best_3x5x5_R60.npz](witnesses/f2_best_3x5x5_R60.npz), $R = 60$ vs $58$), $\langle 4, 4, 5 \rangle$ ([f2_best_4x4x5_R61.npz](witnesses/f2_best_4x4x5_R61.npz), $R = 61$ vs $60$) et $\langle 4, 5, 5 \rangle$ ([f2_best_4x5x5_R75.npz](witnesses/f2_best_4x5x5_R75.npz), $R = 75$ vs $73$), plus de $500$ millions de flips confirment la largeur très étroite des entonnoirs de transition PLUS au-delà de $R > 55$.

---

## 5. Conclusion générale et perspectives de recherche

Cette campagne met en place une chaîne cohérente d'outils reproductibles — le solveur **Parseval-LM**, la loi d'échelle empirique **$\alpha = p / (2p + 2)$** sous recuit $L_{2/3}$, le diagnostic du verrouillage sous barrière $L_{\infty}$ et la **réduction exacte de tranches non destructive** dans le Graphe de Flips sur $\mathbb{F}_2/\mathbb{F}_3/\mathbb{F}_5$ (jusqu'à $16{,}7\text{ M flips/s}$) — tout en caractérisant quantitativement les deux verrous complémentaires qui encadrent la recherche : domination volumétrique des vallées de bord en continu et explosion projective en discret.

Deux prolongements directs découlent de ce bilan. Le premier consiste à concevoir un **solveur en séries de Laurent explicites en $\varepsilon$** ($L(\varepsilon, X) = \sum_{j=0}^e \varepsilon^j L^{(j)}$) afin de supprimer l'explosion de norme flottante (Mur n° 4) en résolvant directement les équations de Parseval-LM ou de Gauss–Newton sur des coefficients bornés aux ordres $e = 2, 3$ pour attaquer $\overline{\mathrm{dc}}(\mathrm{perm}_3) \stackrel{?}{=} 6$ et $\underline{R}(\langle 3, 3, 3 \rangle) \stackrel{?}{=} 19$. Le second consiste à coupler le **Graphe de Flips `uint64` sur $\mathbb{F}_2$** au relèvement de Hensel mod $2^k$ et aux invariants de symétrie de groupe (actions diédrales ou cycliques) afin d'exploiter la cadence de $16{,}7\text{ M flips/s}$ dans des sous-variétés symétriques avant relèvement vers $\mathbb{Z}$ et $\mathbb{Q}$.

---

## 6. Mentions, reproductibilité et bibliographie

### 6.1 Mentions scientifiques et reproductibilité
Cette étude fournit des solveurs calibrés, des témoins exacts vérifiés bit-à-bit sur 16 formats (reproduisant les meilleures bornes connues sur $\mathbb{F}_2$ jusqu'à $\langle 4, 4, 8 \rangle = 94$) et une caractérisation quantitative des obstructions numériques et combinatoires. L'ensemble des **20 scripts** (`src/`), des **15 fichiers JSON** (`results/`), des **31 témoins algébriques `.npz`** (`witnesses/`) et des **8 figures** (`figures/`) est archivé dans ce dossier et indexé dans [README.md](README.md). Tous les calculs des Phases 5–7 ont été exécutés en précision exacte (`int64` / `uint64` / arithmétique modulaire entière), et les solveurs continus des Phases 1–4 en double précision (`float64` / `complex128`).

### 6.2 Bibliographie

- Alper, J., Bogart, T., & Velasco, M. (2017). A lower bound for the determinantal complexity of a hypersurface. *Foundations of Computational Mathematics*, 17(3), 829–836.
- Arai, Y., Ichikawa, Y., & Hukushima, K. (2024). Adaptive Flip Graph Algorithm for Matrix Multiplication. *Proceedings of ISSAC 2024*.
- Bini, D., Capovani, M., Romani, F., & Lotti, G. (1979). $O(n^{2.7799})$ complexity for $n \times n$ approximate matrix multiplication. *Information Processing Letters*, 8(5), 234–235.
- Bläser, M. (2003). On the complexity of the multiplication of matrices of small formats. *Journal of Complexity*, 19(1), 43–60.
- Conner, A., Harper, A., & Landsberg, J. M. (2023). New lower bounds for matrix multiplication and $\det_3$. *Forum of Mathematics, Pi*, 11, e17.
- Grenet, B. (2011). An upper bound for the permanent versus determinant problem. *Manuscript*.
- Kauers, M., & Moosbauer, J. (2023). Flip graphs for matrix multiplication. *Proceedings of ISSAC 2023*.
- Laderman, J. D. (1976). A noncommutative algorithm for multiplying $3 \times 3$ matrices using 23 multiplications. *Bull. Amer. Math. Soc.*, 82(1), 126–128.
- Moosbauer, J., & Poole, M. (2025). Flip graphs with PLUS transitions for fast matrix multiplication over finite fields and rationals. *Preprint*.
- Novikov, A., et al. / Google DeepMind (2025). AlphaEvolve: A coding agent for scientific and algorithmic discovery.
- Schönhage, A. (1981). Partial and total matrix multiplication. *SIAM J. Comput.*, 10(3), 434–455.
- Smirnov, A. V. (2013). The bilinear complexity and practical algorithms for matrix multiplication. *Computational Mathematics and Mathematical Physics*, 53(12), 1781–1795.
- Valiant, L. G. (1979). Completeness classes in algebra. *Proceedings of STOC 1979*, 249–261.
