import WeakPolyComposite.Basic

/-!
# Theorem D6 (3-Root `B_2`-Symmetry on `mod 10`, Lucas Even-Multiple Slices, and `m = 3q` Coupling)

Formalizes the core arithmetic obstructions of Theorem D6:
1. `D6.1` (`p = 5`, `n = 33`): the simultaneous 3-root `D4` inequalities with `B_2(w_1, w_2)`
   parity (`y_1` odd, `w_2` odd) have no integer solution, proving `N_gen^(10)(2, OR) ≤ 32`.
2. `D6.2` (Lucas Even-Multiple Slice Theorem): `(6 - 3 : ZMod p) ≠ 0` for every prime `p ≥ 5`,
   yielding `n ≤ 7p - 4` whenever `B_p(w_1, ·)` is constant on `U_1`.
3. `D6.3` (`m = 3q`): `19(q - 1) + 6 < 20(q - 1) + 1` for every prime `q ≥ 7`.
-/

namespace WeakPolyComposite

/-- Core arithmetic contradiction of Theorem D6.1 at `p = 5, n = 33`:
if `w1 = 8`, `w2 ≤ 11` odd, `v0 + v1 = 25 - w2`, `y0 + y1 = 8` with `y1` odd,
then the two D4 inequalities on `w2` and `w1 + w2` cannot both hold. -/
theorem d6_mod10_n33_impossible (w2 y0 y1 v0 v1 : ℕ)
    (hw2_le : w2 ≤ 11) (hw2_odd : w2 % 2 = 1)
    (hy_sum : y0 + y1 = 8) (hy1_odd : y1 % 2 = 1)
    (hv_sum : v0 + v1 + w2 = 25)
    (h_d4_w2 : v0 + y0 + (v1 + y1) / 2 ≤ 12)
    (h_d4_w12 : v1 + v0 / 2 ≤ 12) : False := by
  omega

/-- Lucas slice elimination coefficient of Theorem D6.2: for every prime `p ≥ 5`, `(3 : ZMod p) ≠ 0`. -/
theorem d6_lucas_three_ne_zero (p : ℕ) [Fact (Nat.Prime p)] (hp5 : 5 ≤ p) :
    (3 : ZMod p) ≠ 0 := by
  intro heq
  have hdvd : p ∣ 3 := (ZMod.natCast_eq_zero_iff 3 p).mp heq
  have hle : p ≤ 3 := Nat.le_of_dvd (by omega) hdvd
  omega

/-- Theorem D6.3 beats Theorem D3 on `m = 3q` for every prime `q ≥ 7`. -/
theorem d6_m3q_beats_d3 (q : ℕ) (hq : 7 ≤ q) :
    19 * (q - 1) + 6 < 20 * (q - 1) + 1 := by
  omega

end WeakPolyComposite
