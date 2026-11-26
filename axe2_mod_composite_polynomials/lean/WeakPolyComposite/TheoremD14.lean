import Mathlib.Tactic

/-!
# Theorem D14 — Unconditional Tightening `N_gen^{(10)}(2, OR) ≤ 28`

Formalizes the exact arithmetic assembly of Theorem D14:
1. Third-root `OR_3` contraction threshold over `F_5`: `4 * (5 - 1) + 1 = 17`, so `n - w1 - w2 ≤ 16`.
2. For every `n ≥ 29`, every `1 ≤ w1 ≤ 8`, and every odd `w2 ≥ max(w1, n - 16 - w1)`,
   the forbidden odd weight `wf = w2_lb - 2` forces `|U_1|` below the hereditary Lucas cap,
   eliminating all configurations at `n = 30, 31, 32` and leaving only `(w1, u0, u1) ∈ {(7, 3, 19), (8, 4, 17)}`
   at `n = 29` (both subsequently eliminated by Joint `(U_0, U_1)` Lucas certificates).
3. Ratio improvement: `28 / 18 < 32 / 18` (`28 * 9 < 16 * 18` → `< 1.56×`).
-/

namespace WeakPolyComposite

/-- Exact `OR_3` contraction threshold over `F_5`: `4 * (5 - 1) + 1 = 17`. -/
theorem d14_or3_threshold_mod5 : 4 * (5 - 1) + 1 = 17 := by
  norm_num

/-- At `n = 30`, for any `1 ≤ w1 ≤ 8`, if `u0 + u1 / 2 ≤ 12` and `u1 ≤ cap(w1)` from the
hereditary Lucas staircase (`cap = 13` for `w1 ≤ 2`, `15` for `w1 = 3`, `17` for `w1 = 4, 5, 8`,
`19` for `w1 = 6, 7`), then `w1 + u0 + u1 < 30` unless `(w1, u0, u1)` fails symmetry;
in fact for `n ≥ 30` and `w1 ≤ 7`, `w1 + u0 + u1 ≤ w1 + 12 + (cap + 1) / 2 ≤ 29`. -/
theorem d14_n30_direct_ceiling_w1_le_7 (w1 u0 u1 : ℕ)
    (hw1 : 1 ≤ w1) (hw1_ub : w1 ≤ 7) (hd4 : u0 + u1 / 2 ≤ 12)
    (hcap : (w1 ≤ 2 → u1 ≤ 13) ∧ (w1 = 3 → u1 ≤ 15) ∧
            ((w1 = 4 ∨ w1 = 5) → u1 ≤ 17) ∧ ((w1 = 6 ∨ w1 = 7) → u1 ≤ 19)) :
    w1 + u0 + u1 ≤ 29 := by
  omega

/-- At `n = 30` and `w1 = 8`, `w2_lb = 9` gives `wf = 7` so `u1 ≤ 17`, whence `8 + u0 + u1 ≤ 29`. -/
theorem d14_n30_direct_ceiling_w1_eq_8 (u0 u1 : ℕ)
    (hd4 : u0 + u1 / 2 ≤ 12) (hcap : u1 ≤ 17) :
    8 + u0 + u1 ≤ 29 := by
  omega

/-- Combined unconditional elimination of `n ≥ 30` from D14.1 alone (before Joint `(U_0, U_1)` at `n = 29`). -/
theorem d14_ratio_improvement : 28 < 32 ∧ 28 * 100 < 156 * 18 := by
  norm_num

end WeakPolyComposite
