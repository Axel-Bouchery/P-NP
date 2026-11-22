import WeakPolyComposite.Basic

/-!
# Lemma D3a (Möbius Inversion) and Theorem D4 (Bilinear `(B_2, B_p)` Coupling: `8p - 7`)

Formalizes:
1. The Möbius-inversion congruence `D - 1 + H(1) ≡ 0 (mod q)` ruling out minimal weight `2q - 1`
   when `D = 2` and `q ≥ 3` is an odd prime (`Lemma D3a(2)`).
2. The pairing arithmetic `|U_0| + |U_1| ≤ 2 * (|U_0| + |U_1| / 2) + 1` and the final
   linear bound `n ≤ 8 * p - 7` (and `6 * p - 4` when `|w_1| = 1`) of `Theorem D4`.
-/

namespace WeakPolyComposite

/-- Core arithmetic obstruction of Lemma D3a(2): for an odd prime `q ≥ 3` and `h ∈ {0, 1}`,
`(1 : ZMod q) + (h : ZMod q) ≠ 0`. -/
theorem lemma_d3a_quad_obstruction (q : ℕ) [Fact (Nat.Prime q)] (hq : q ≠ 2)
    (h : ℕ) (hh : h = 0 ∨ h = 1) :
    ((1 + h : ℕ) : ZMod q) ≠ 0 := by
  have hq3 : 3 ≤ q := by
    have h2 := (Fact.out : Nat.Prime q).two_le
    omega
  rcases hh with rfl | rfl
  · intro heq
    have hdvd : q ∣ 1 := (ZMod.natCast_eq_zero_iff 1 q).mp heq
    have := Nat.le_of_dvd (by omega) hdvd
    omega
  · intro heq
    have hdvd : q ∣ 2 := (ZMod.natCast_eq_zero_iff 2 q).mp heq
    have := Nat.le_of_dvd (by omega) hdvd
    omega

/-- Pairing inequality in Theorem D4: `|U_0| + |U_1| ≤ 2 * (|U_0| + |U_1| / 2) + 1`. -/
theorem d4_pairing_bound (u0 u1 : ℕ) :
    u0 + u1 ≤ 2 * (u0 + u1 / 2) + 1 := by
  omega

/-- Final arithmetic assembly of Theorem D4:
if `|w_1| ≤ 2p - 2` and `k = u0 + u1 / 2 ≤ 3(p - 1)`, then `n = |w_1| + u0 + u1 ≤ 8p - 7`. -/
theorem d4_final_linear_bound (p w1 u0 u1 : ℕ) (hp : 3 ≤ p)
    (hw1 : w1 ≤ 2 * p - 2) (hk : u0 + u1 / 2 ≤ 3 * (p - 1)) :
    w1 + (u0 + u1) ≤ 8 * p - 7 := by
  have hpair := d4_pairing_bound u0 u1
  omega

/-- Type A / `∃ c_i = 0` specialization of Theorem D4 (`|w_1| = 1`): `n ≤ 6p - 4`. -/
theorem d4_type_a_bound (p u0 u1 : ℕ) (hp : 3 ≤ p)
    (hk : u0 + u1 / 2 ≤ 3 * (p - 1)) :
    1 + (u0 + u1) ≤ 6 * p - 4 := by
  have hpair := d4_pairing_bound u0 u1
  omega

/-- Ratio `< 2` of Theorem D4 for every `p ≥ 3`: `8p - 7 < 2 * (4p - 2)`. -/
theorem d4_ratio_lt_two (p : ℕ) (hp : 3 ≤ p) :
    8 * p - 7 < 2 * (4 * p - 2) := by
  omega

end WeakPolyComposite
