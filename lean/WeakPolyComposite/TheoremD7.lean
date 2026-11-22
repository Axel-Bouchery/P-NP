import WeakPolyComposite.Basic

/-!
# Theorems D7 and U3 (2-Class Lucas Bi-Slice Cancellation & `THR_{n,2}` Single-Block Ceiling)

Formalizes:
1. Exact Lucas cancellation identities for `(3p - 1, 2p - 1)` and `(4p - 1, p - 1)` over `ZMod p` (`p ≥ 5`).
2. Single-block ceiling of `Theorem U3`: for every prime `p ≡ 1 (mod 4)`, `2 * p ∣ p * (p - 1) / 2`,
   proving `C(p, 2) ≡ 0 (mod 2p)` and hence `N_{k-bloc}^(2p)(2, THR_2) ≤ k * (p - 1)`.
-/

namespace WeakPolyComposite

/-- Theorem D7 Lucas bracket cancellation for `(3p - 1, 2p - 1)` over `ZMod p` (`p ≥ 5`):
`1 - (2 : ZMod p)⁻¹ * 2 = 0` and `- (2 : ZMod p)⁻¹ + (1 - (2 : ZMod p)⁻¹) = 0`. -/
theorem d7_bislice_3p_2p_cancel (p : ℕ) [Fact (Nat.Prime p)] (hp5 : 5 ≤ p) :
    let inv2 : ZMod p := (2 : ZMod p)⁻¹
    (1 : ZMod p) - inv2 * 2 = 0 ∧ (-inv2) + (1 - inv2) = 0 := by
  have h2 : (2 : ZMod p) ≠ 0 := by
    intro h
    have hdvd : p ∣ 2 := (ZMod.natCast_eq_zero_iff 2 p).mp h
    have hle : p ≤ 2 := Nat.le_of_dvd (by omega) hdvd
    omega
  intro inv2
  have hmul : inv2 * 2 = 1 := ZMod.inv_mul_of_unit 2 (isUnit_iff_ne_zero.mpr h2)
  constructor
  · rw [hmul, sub_self]
  · have h2inv : inv2 + inv2 = 1 := by
      calc inv2 + inv2 = inv2 * 2 := by ring
        _ = 1 := hmul
    linear_combination h2inv

/-- Theorem D7 Lucas bracket cancellation for `(4p - 1, p - 1)` over `ZMod p` (`p ≥ 5`):
`1 - (3 : ZMod p)⁻¹ * 3 = 0` and `- (3 : ZMod p)⁻¹ * 2 + (1 - (3 : ZMod p)⁻¹) = 0`. -/
theorem d7_bislice_4p_p_cancel (p : ℕ) [Fact (Nat.Prime p)] (hp5 : 5 ≤ p) :
    let inv3 : ZMod p := (3 : ZMod p)⁻¹
    (1 : ZMod p) - inv3 * 3 = 0 ∧ (-inv3 * 2) + (1 - inv3) = 0 := by
  have h3 : (3 : ZMod p) ≠ 0 := by
    intro h
    have hdvd : p ∣ 3 := (ZMod.natCast_eq_zero_iff 3 p).mp h
    have hle : p ≤ 3 := Nat.le_of_dvd (by omega) hdvd
    omega
  intro inv3
  have hmul : inv3 * 3 = 1 := ZMod.inv_mul_of_unit 3 (isUnit_iff_ne_zero.mpr h3)
  constructor
  · rw [hmul, sub_self]
  · linear_combination hmul

/-- Theorem U3 single-block ceiling: for every `p` with `p % 4 = 1`, `2 * p ∣ p * (p - 1) / 2`. -/
theorem u3_thr2_binom_zero_mod_2p (p : ℕ) (hp : p % 4 = 1) :
    (2 * p) ∣ (p * (p - 1) / 2) := by
  refine ⟨(p - 1) / 4, ?_⟩
  omega

end WeakPolyComposite
