import WeakPolyComposite.Basic

/-!
# Theorem D8 (Lucas–Vandermonde Weight-Constrained Chevalley–Warning Theorem)

Formalizes the exact arithmetic identities underlying Theorem D8:
1. At `D = 1`: `det A(1, r) = C(r, 0) = 1` (valid for every prime `p` and every `r ≥ 1`).
2. At `D = 2`: `det A(2, r) = C(2r, r) * 1 - C(2r - 1, r) * 1 = C(2r - 1, r)` via Pascal's
   identity and symmetry `C(2r - 1, r - 1) = C(2r - 1, r)`.
3. Strict variable saving over rigid `r`-tupling Chevalley–Warning:
   `(D * r + 1) * p - 1 < (D + 1) * r * (p - 1) + r` whenever `r ≥ 2` and `p ≥ D * r`
   (and for all `p ≥ 2` at `D = 1`), with exact saving `(r - 1) * p - D * r + 1`.
4. Exact prime-power check `C(2p - 1, p) ≡ 1 (mod p)` on concrete small primes.
-/

namespace WeakPolyComposite

/-- At `D = 1`, the `1 × 1` Lucas matrix determinant is `C(r, 0) = 1`. -/
theorem d8_det_D1 (r : ℕ) : Nat.choose r 0 = 1 :=
  Nat.choose_zero_right r

/-- At `D = 2`, the `2 × 2` Lucas matrix determinant `C(2r, r) - C(2r - 1, r)` equals
`C(2r - 1, r)` for every `r ≥ 1`. -/
theorem d8_det_D2_formula (r : ℕ) (hr : 1 ≤ r) :
    Nat.choose (2 * r) r - Nat.choose (2 * r - 1) r = Nat.choose (2 * r - 1) r := by
  have h2r : 2 * r = (2 * r - 1) + 1 := by omega
  have hr1 : r = (r - 1) + 1 := by omega
  have hpascal : Nat.choose (2 * r) r =
      Nat.choose (2 * r - 1) (r - 1) + Nat.choose (2 * r - 1) r := by
    conv_lhs => rw [h2r, hr1]
    exact Nat.choose_succ_succ (2 * r - 1) (r - 1)
  have hsymm : Nat.choose (2 * r - 1) (r - 1) = Nat.choose (2 * r - 1) r := by
    have hsub : (2 * r - 1) - (r - 1) = r := by omega
    rw [← Nat.choose_symm (by omega : r - 1 ≤ 2 * r - 1), hsub]
  rw [hpascal, hsymm]
  omega

/-- Concrete values of `|det A(2, r)| = C(2r - 1, r)` for `r = 2, 3, 4, 5`. -/
theorem d8_det_D2_values :
    Nat.choose 3 2 = 3 ∧ Nat.choose 5 3 = 10 ∧
    Nat.choose 7 4 = 35 ∧ Nat.choose 9 5 = 126 := by
  decide

/-- Strict variable saving of Theorem D8 `(D*r + 1)*p - 1` over rigid `r`-tupling
`(D + 1)*r*(p - 1) + r` in the universal regime `r ≥ 2`, `p ≥ D * r`, `D ≥ 1`. -/
theorem d8_saving_over_rigid (D r p : ℕ) (hD : 1 ≤ D) (hr : 2 ≤ r) (hp : D * r ≤ p) :
    (D * r + 1) * p - 1 < (D + 1) * r * (p - 1) + r := by
  have hp1 : 1 ≤ p := by nlinarith
  have h1 : (D + 1) * r * (p - 1) + r = (D * r + r) * p - D * r := by
    zify [hp1, (by nlinarith : D * r ≤ (D * r + r) * p)]
    ring
  have h2 : (D * r + 1) * p ≤ (D * r + r) * p - D * r := by
    have hgap : D * r ≤ (r - 1) * p := by nlinarith
    zify [hgap, (by nlinarith : D * r ≤ (D * r + r) * p)]
    nlinarith
  have hpos : 1 ≤ (D * r + 1) * p := by nlinarith
  omega

/-- At `D = 1`, Theorem D8 strictly beats rigid `r`-tupling for EVERY `r ≥ 2` and `p ≥ 2`. -/
theorem d8_saving_D1 (r p : ℕ) (hr : 2 ≤ r) (hp : 2 ≤ p) :
    (r + 1) * p - 1 < 2 * r * (p - 1) + r := by
  have hp1 : 1 ≤ p := by omega
  have hmul : (r + 1) * p ≤ 2 * r * (p - 1) + r := by
    zify [hp1]
    nlinarith
  have hpos : 1 ≤ (r + 1) * p := by nlinarith
  omega

/-- Concrete prime-power corollary checks `C(2p - 1, p) % p = 1` for `p ∈ {2, 3, 5, 7}`. -/
theorem d8_prime_power_small :
    Nat.choose 3 2 % 2 = 1 ∧ Nat.choose 5 3 % 3 = 1 ∧
    Nat.choose 9 5 % 5 = 1 ∧ Nat.choose 13 7 % 7 = 1 := by
  decide

end WeakPolyComposite
