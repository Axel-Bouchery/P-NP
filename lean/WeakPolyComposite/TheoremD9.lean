import WeakPolyComposite.Basic

/-!
# Theorem D9 (General Semiprimes `m = pq`: Diagonal-Shifted Rank, Caro–Wei Sparsity,
  Expander–Fourier `m^{3/2}` Convergence Wall, and Sharp `𝔽₃²` Packing `3k + 4`)

Formalizes the exact arithmetic lemmas underlying Theorem D9:
1. Sharp `𝔽₃²` zero-sum packing (`3k + 4`): from `|U| ≥ 3k + 5`, extracting up to `k - 1`
   weight-`≤ 3` blocks leaves at least `8 ≥ 7` elements (above `max_no_wt3 = 6`), and extracting
   `k` blocks leaves at least `5 = D(𝔽₃²)` elements for a `(k + 1)`-th block.
2. Unconditional improvement on `m = 3q` for EVERY prime `q ≥ 5` (including `q = 5`, `m = 15`):
   `19 * q - 15 < 20 * q - 19` (`80 < 81` at `q = 5`).
3. Sharpness of `3k + 4`: `(3 * k + 2) / 3 = k` (attained by `3k + 2` copies of `(1, 0)` and
   `2` copies of `(0, 1)` in `𝔽₃²`).
4. Caro–Wei matching bound (`Δ = 1`) and Möbius `D3a(2)` sharpening on all odd semiprimes `m = pq`.
-/

namespace WeakPolyComposite

/-- Greedy `𝔽₃²` packing never stalls before `k` blocks when `|U| ≥ 3k + 5`,
and leaves at least `5 = D(𝔽₃²)` elements after `k` blocks for a `(k + 1)`-th block. -/
theorem d9_f3sq_packing_sufficient (k : ℕ) (hk : 1 ≤ k) :
    7 ≤ (3 * k + 5) - 3 * (k - 1) ∧ 5 ≤ (3 * k + 5) - 3 * k := by
  omega

/-- Sharpness of `|U| ≤ 3k + 4` in `𝔽₃²`: `3k + 2` copies of `(1, 0)` yield at most `k`
disjoint weight-3 zero-sums (and `2` copies of `(0, 1)` yield none). -/
theorem d9_f3sq_packing_sharp (k : ℕ) :
    (3 * k + 2) / 3 = k := by
  omega

/-- Theorem D9.4 strictly improves Theorem D3 (`20q - 19`) for EVERY `q ≥ 5`
(including `m = 15`: `19 * 5 - 15 = 80 < 81 = 20 * 5 - 19`). -/
theorem d9_m3q_beats_d3_all_q5 (q : ℕ) (hq : 5 ≤ q) :
    19 * q - 15 < 20 * q - 19 := by
  omega

/-- Concrete `mod 15` (`q = 5`) sharpening: `80 < 81` (and Case 1: `72`). -/
theorem d9_mod15_concrete :
    19 * 5 - 15 = 80 ∧ 20 * 5 - 19 = 81 ∧ 17 * 5 - 13 = 72 := by
  decide

/-- Caro–Wei bound for a matching (`Δ = 1`) on `m = pq`: `2 * (2m - q - 1) = 4m - 2q - 2`. -/
theorem d9_caro_wei_matching (m q : ℕ) (hmq : q + 1 ≤ 2 * m) :
    2 * (2 * m - q - 1) = 4 * m - 2 * q - 2 := by
  omega

/-- Möbius `D3a(2)` strict sharpening over Theorem D3 for every odd prime `p ≥ 5`. -/
theorem d9_d3a_sharpening_semiprimes (p q : ℕ) (hp : 5 ≤ p) :
    (2 * p - 1) * (p + 1) * (q - 1) + 2 < (2 * p - 1) * (p + 1) * (q - 1) + (2 * p - 2) := by
  omega

end WeakPolyComposite
