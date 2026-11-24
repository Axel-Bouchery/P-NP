import Mathlib.Data.Real.Basic
import Mathlib.Tactic

/-!
# Theorem D11 — Log-Quasi-Regular Core, Universal Phase-Shifted Cosine Bound & CRT Rank

Algebraic core lemmas for Theorem D11:
1. Exact quadratic lower-bound identities underlying Lemma D11.2 (Universal Phase-Shifted Cosine Bound):
   - Inner regime `|α| ≤ π/4`: `(α + δ)² - (1/4) δ² - 2 α δ = α² + (3/4) δ² ≥ 0`.
   - Outer regime `α ∈ [π/4, π/2]`, `δ ∈ [-π/4, π/2 - α]`:
     at the boundary `δ = -π/4`,
     `(α - π/4)² - (1/4)(-π/4)² - ((π/2 - α)/2)(-π/4) = (α - 5π/16)² + (3/256) π² > 0`.
   - Outer regime `δ ∈ [π/2 - α, π/4]`: `(1/4) δ² + ((π/2 - α)/2) δ ≤ (3/4) δ² ≤ (3/4) (π/4)² < (π/4)²`.
2. Core peeling density contraction ratio `(1/2) / (11/12)² = 72 / 121 < 2 / 3` (Lemma D11.1).
-/

namespace WeakPolyComposite

/-- Inner-regime quadratic identity of Lemma D11.2. -/
theorem d11_phase_cosine_inner (α δ : ℝ) :
    (α + δ) ^ 2 - (1 / 4 : ℝ) * δ ^ 2 - (2 * α) * δ = α ^ 2 + (3 / 4 : ℝ) * δ ^ 2 := by
  ring

/-- Inner-regime non-negativity of Lemma D11.2. -/
theorem d11_phase_cosine_inner_nonneg (α δ : ℝ) :
    0 ≤ (α + δ) ^ 2 - (1 / 4 : ℝ) * δ ^ 2 - (2 * α) * δ := by
  rw [d11_phase_cosine_inner]
  nlinarith

/-- Outer-regime boundary identity at `δ = -P/4` (`P = π`): exact minimum `(3/256) P²`. -/
theorem d11_phase_cosine_outer_boundary (α P : ℝ) :
    (α - P / 4) ^ 2 - (1 / 4 : ℝ) * (-P / 4) ^ 2 - ((P / 2 - α) / 2) * (-P / 4)
      = (α - 5 * P / 16) ^ 2 + (3 / 256 : ℝ) * P ^ 2 := by
  ring

/-- Outer-regime monotonicity: for `α ≥ P/4` and `δ ≥ -P/4`, the difference grows with `δ`. -/
theorem d11_phase_cosine_outer_left (α δ P : ℝ) (hα : P / 4 ≤ α) (hδ : -P / 4 ≤ δ) :
    (3 / 256 : ℝ) * P ^ 2 ≤
      (α + δ) ^ 2 - (1 / 4 : ℝ) * δ ^ 2 - ((P / 2 - α) / 2) * δ := by
  have h1 : (α + δ) ^ 2 - (1 / 4 : ℝ) * δ ^ 2 - ((P / 2 - α) / 2) * δ
      - ((α - P / 4) ^ 2 - (1 / 4 : ℝ) * (-P / 4) ^ 2 - ((P / 2 - α) / 2) * (-P / 4))
      = (δ + P / 4) * ((3 / 4 : ℝ) * (δ + P / 4) + (5 / 2 : ℝ) * (α - P / 4)) := by
    ring
  have h2 := d11_phase_cosine_outer_boundary α P
  nlinarith

/-- Outer-regime right half (`P/2 - α ≤ δ ≤ P/4`): upper bound `(3/4)(P/4)²` vs `(P/4)²`. -/
theorem d11_phase_cosine_outer_right (α δ P : ℝ)
    (h1 : P / 2 - α ≤ δ) (h2 : δ ≤ P / 4) (hδ : 0 ≤ δ) :
    (1 / 4 : ℝ) * δ ^ 2 + ((P / 2 - α) / 2) * δ ≤ (3 / 4 : ℝ) * (P / 4) ^ 2 := by
  nlinarith

/-- Core peeling density contraction factor for `K ≥ 12`: `(1/2)/(11/12)² = 72/121 < 2/3`. -/
theorem d11_core_density_contraction : (72 : ℚ) / 121 < 2 / 3 := by
  norm_num

end WeakPolyComposite
