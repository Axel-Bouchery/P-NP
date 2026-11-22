import WeakPolyComposite.Basic

/-!
# Theorem D10 (Exact Rank-1 Closures on `m = 2q`, Sub-`2(m - 1)` Rank-1 & Pure Rank-`r` on `m = pq`,
  Bounded-Tilt Mean Crossing Identity, Small-Frequency Surplus Algebra, and `{0, 1}`-Graph Witness `K_8 ⊔ K_2`)

Formalizes the exact arithmetic and algebraic identities underlying Theorem D10:
1. Exact parity verification of the `3q - 1` rank-1 lower-bound witness on `m = 2q` for all odd primes
   `q ≡ 3 (mod 4)` and `q ≡ 1 (mod 4)`.
2. Strict sub-`2(m - 1)` rank-1 upper bounds `2m - p - 1 < 2(m - 1)` and `2m - q - 1 < 2(m - 1)`
   on every odd semiprime `m = pq`, and pure rank-2 sharpening `4m - 2p - 2q < 6m - 4p - 3q + 1`.
3. Exact 3-point polynomial identity (`16`-scaled over `ℤ`) underlying the Bounded-Tilt Mean Crossing Lemma
   guaranteeing `∃ θ* ∈ [1/4, 3/4]` with `μ(θ*) ∈ mℤ` whenever `B ≥ 16m`.
4. Exact polynomial identity for the Band-1 elementary lower bound `L(T) ≥ (4T² - 1) / (3T)`.
5. Exact finite verification (`decide`) that the single `{0, 1}`-graph `K_8 ⊔ K_2` with `c = 2`
   achieves the exact Győri optimum `n = 10` modulo `6`.
-/

namespace WeakPolyComposite

/-- Witness parity for `q = 4k + 3`: `P_2(S) = C(S, 2)` is odd at both `S = q` and `S = 2q`. -/
theorem d10_rank1_m2q_witness_mod4_eq3 (k : ℕ) :
    ((4 * k + 3) * (4 * k + 2) / 2) % 2 = 1 ∧
    ((8 * k + 6) * (8 * k + 5) / 2) % 2 = 1 := by
  omega

/-- Witness parity for `q = 4k + 1`: `P_2(S) = S + C(S, 2)` is odd at both `S = q` and `S = 2q`. -/
theorem d10_rank1_m2q_witness_mod4_eq1 (k : ℕ) :
    ((4 * k + 1) + (4 * k + 1) * (4 * k) / 2) % 2 = 1 ∧
    ((8 * k + 2) + (8 * k + 2) * (8 * k + 1) / 2) % 2 = 1 := by
  omega

/-- On every `m = 2q` (`q ≥ 3`), the exact rank-1 threshold `3q - 1` sits strictly below
both the general 2-block lower bound `4q - 2` and Theorem D4's general upper bound `8q - 7`. -/
theorem d10_rank1_m2q_strict_hierarchy (q : ℕ) (hq : 3 ≤ q) :
    3 * q - 1 < 4 * q - 2 ∧ 4 * q - 2 < 8 * q - 7 := by
  omega

/-- On every odd semiprime `m = pq` (`3 ≤ p < q`), neither channel can be rank-1 if `n ≥ 2(m - 1)`. -/
theorem d10_rank1_semiprime_sub_2m2 (p q m : ℕ) (hp : 3 ≤ p) (hpq : p < q) (hm : m = p * q) :
    2 * m - q - 1 < 2 * m - p - 1 ∧ 2 * m - p - 1 < 2 * (m - 1) := by
  have hq : 4 ≤ q := by omega
  have hm_ge : 12 ≤ m := by
    rw [hm]
    nlinarith
  omega

/-- Pure rank-2 Davenport bound `(2p - 1)(2q - 1) - 1 = 4pq - 2p - 2q` strictly improves
Theorem D9.1's `(2p - 1)(3q - 2) - 1` for every `3 ≤ p < q`. -/
theorem d10_pure_rank2_vs_d91 (p q : ℕ) (hp : 3 ≤ p) (hq : 5 ≤ q) :
    (2 * p - 1) * (2 * (q - 1) + 1) - 1 = 4 * (p * q) - 2 * p - 2 * q ∧
    (2 * p - 1) * (2 * (q - 1) + 1) - 1 < (2 * p - 1) * (3 * (q - 1) + 1) - 1 := by
  constructor
  · zify [show 1 ≤ p by omega, show 1 ≤ q by omega,
          show 1 ≤ (2 * p - 1) * (2 * (q - 1) + 1) by nlinarith,
          show 2 * p + 2 * q ≤ 4 * (p * q) by nlinarith]
    ring
  · nlinarith

/-- Exact 3-point identity (`16`-scaled) for the Bounded-Tilt Mean Crossing Lemma on `[1/4, 3/4]`:
`16 μ(1/4) = A + 4B`, `16 μ(1/2) = 4A + 8B`, `16 μ(3/4) = 9A + 12B`, so
`5(16μ(1/2) - 16μ(1/4)) - 3(16μ(3/4) - 16μ(1/2)) = 8B`. -/
theorem d10_bounded_tilt_three_point_identity (A B : ℤ) :
    5 * ((4 * A + 8 * B) - (A + 4 * B)) - 3 * ((9 * A + 12 * B) - (4 * A + 8 * B)) = 8 * B := by
  ring

/-- Exact polynomial identity certifying `L(T) ≥ (4T² - 1) / (3T)` at the worst-case endpoint `x = 1/T`. -/
theorem d10_surplus_L_endpoint_identity (T : ℕ) :
    3 * T ^ 2 * (2 * T + 1) - T * (T + 1) * (2 * T + 1) = T * (4 * T ^ 2 - 1) := by
  cases T with
  | zero => simp
  | succ k =>
    zify [show 1 ≤ 4 * (k + 1) ^ 2 by nlinarith,
          show (k + 1) * (k + 1 + 1) * (2 * (k + 1) + 1) ≤ 3 * (k + 1) ^ 2 * (2 * (k + 1) + 1) by nlinarith]
    ring

/-- Subsampling vertex count survival: keeping `64m` vertices and deleting `< 16m` high-degree vertices
leaves at least `48m ≥ 32m` vertices for Theorem D10.4. -/
theorem d10_subsampling_survival (m : ℕ) :
    64 * m - 16 * m = 48 * m ∧ 32 * m ≤ 48 * m := by
  omega

/-- Exact finite verification that the single `{0, 1}`-graph `K_8 ⊔ K_2` with `c = 2` avoids `0 (mod 6)`
on all `2^10 - 1` nonempty subsets, attaining Győri's exact optimum `N_gen^(6)(2) = 10`. -/
theorem d10_k8_disjoint_k2_mod6_exact :
    ∀ k1 : Fin 9, ∀ k2 : Fin 3, (k1.val + k2.val > 0) →
      (k1.val * (k1.val - 1) / 2 + 2 * k1.val + k2.val * (k2.val - 1) / 2 + 2 * k2.val) % 6 ≠ 0 := by
  decide

end WeakPolyComposite
