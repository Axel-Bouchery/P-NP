import Mathlib.Tactic

namespace WeakPolyComposite

/-!
# Theorem D15 — Stronger Phase-Shifted Cosine, Monotonicity Lemma & Five-Band Crossover

Formalizes the core exact algebraic identities underlying Theorem D15:
1. Stronger Phase-Shifted Cosine Lemma (D15.1):
   For all `alpha in [0, P/2]` and `delta in [-P/4, P/4]` (with `P = pi > 0`),
   `dist(alpha + delta, P * Z)^2 >= (1/4) delta^2 + c(alpha) delta + (1/8) alpha^2`,
   where `c(alpha) = 2 alpha` on `[0, P/4]` and `c(alpha) = P/2 - alpha` on `(P/4, P/2]`.
2. Monotonicity Lemma (D15.3):
   For `0 < A <= B` and `x >= 0`, `(x + A)^2 * B >= A^2 * (x + B)`,
   so `(x + A)^2 / (x + B) >= A^2 / B` (global minimum at `x = 0`).
3. Five-Band Crossover Algebra (D15.3):
   The product of Band 2.75 (`n / (K^2 * Delta^2 * u^2)`) and Band 3 (`n * u^2 / Delta`)
   cancels `u^2` identically and equals `(n / (K * Delta^(3/2)))^2`.
-/

/-- D15.1 Case 1 (`0 <= alpha <= P/4`): exact nonnegative quadratic form. -/
theorem d15_phase_cosine_inner (alpha delta : ℚ) :
    (alpha + delta) ^ 2 - (1 / 4 : ℚ) * delta ^ 2 - (2 * alpha) * delta - (1 / 8 : ℚ) * alpha ^ 2
      = (7 / 8 : ℚ) * alpha ^ 2 + (3 / 4 : ℚ) * delta ^ 2 := by
  ring

/-- D15.1 Subcase 2a endpoint (`delta = -P/4`): exact completion of square with min `(5/448) P^2`. -/
theorem d15_phase_cosine_subcase2a_endpoint (alpha P : ℚ) :
    (alpha - P / 4) ^ 2 - (1 / 4 : ℚ) * (-P / 4) ^ 2 - (P / 2 - alpha) * (-P / 4) - (1 / 8 : ℚ) * alpha ^ 2
      = (7 / 8 : ℚ) * (alpha - (3 / 7 : ℚ) * P) ^ 2 + (5 / 448 : ℚ) * P ^ 2 := by
  ring

/-- D15.1 Subcase 2a vertex positivity at `alpha = P/4` and `alpha = 7P/24`. -/
theorem d15_phase_cosine_subcase2a_vertex_endpoints (P : ℚ) :
    (-(17 / 8 : ℚ) * (P / 4) ^ 2 + P * (P / 4) - (1 / 12 : ℚ) * P ^ 2 = (13 / 384 : ℚ) * P ^ 2) ∧
    (-(17 / 8 : ℚ) * ((7 / 24 : ℚ) * P) ^ 2 + P * ((7 / 24 : ℚ) * P) - (1 / 12 : ℚ) * P ^ 2
      = (127 / 4608 : ℚ) * P ^ 2) := by
  constructor <;> ring

/-- D15.1 Subcase 2b endpoint (`delta = P/4`): exact nonnegative identity with min `(1/64) P^2`. -/
theorem d15_phase_cosine_subcase2b_endpoint (alpha P : ℚ) :
    (P - alpha - P / 4) ^ 2 - (1 / 4 : ℚ) * (P / 4) ^ 2 - (P / 2 - alpha) * (P / 4) - (1 / 8 : ℚ) * alpha ^ 2
      = (7 / 8 : ℚ) * (alpha - P / 2) ^ 2 + (3 / 8 : ℚ) * P * (P / 2 - alpha) + (1 / 64 : ℚ) * P ^ 2 := by
  ring

/-- D15.3 Monotonicity Lemma: for `0 < A <= B` and `x >= 0`, `(x + A)^2 / (x + B) >= A^2 / B`. -/
theorem d15_warnke_monotonicity (A B x : ℚ) (hA : 0 < A) (hAB : A ≤ B) (hx : 0 ≤ x) :
    A ^ 2 * (x + B) ≤ (x + A) ^ 2 * B := by
  have hdiff : (x + A) ^ 2 * B - A ^ 2 * (x + B) = x * (x * B + A * (2 * B - A)) := by ring
  rw [← sub_nonneg, hdiff]
  nlinarith

/-- D15.3 Band 2.75 / Band 3 crossover: `u^2` cancels identically in the product. -/
theorem d15_crossover_cancellation (n K Delta u : ℚ) (hK : K ≠ 0) (hD : Delta ≠ 0) (hu : u ≠ 0) :
    (n / (K ^ 2 * Delta ^ 2 * u ^ 2)) * (n * u ^ 2 / Delta) = n ^ 2 / (K ^ 2 * Delta ^ 3) := by
  field_simp

end WeakPolyComposite
