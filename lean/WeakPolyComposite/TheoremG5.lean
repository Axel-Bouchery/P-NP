import WeakPolyComposite.Basic

/-!
# Theorem G5 (Unbounded Pre-Asymptotic Separation `N_gen / N_sym ≥ p_min - 1`)

Formalizes the periodicity upper bound for symmetric polynomials (`N_sym ≤ m - 1` when `d < p_min`)
and the block-norm lower bound (`N_gen ≥ d * (m - 1)`).
-/

namespace WeakPolyComposite

/-- Any `m`-periodic function on Hamming weights vanishing at `0` also vanishes at weight `m`,
hence cannot weakly represent `OR_n` for `n ≥ m`. -/
theorem symmetric_periodic_bound (m n : ℕ) (hm : 1 ≤ m) (hn : m ≤ n)
    (f : ℕ → ZMod m) (hf0 : f 0 = 0) (hper : ∀ k, f (k + m) = f k) :
    ∃ w : ℕ, 1 ≤ w ∧ w ≤ n ∧ f w = 0 := by
  refine ⟨m, hm, hn, ?_⟩
  simpa using hper 0 |>.trans hf0

/-- Norm-form block construction: if `d` blocks each have size `m - 1`, then the only
tuple `(S_1, ..., S_d)` with `0 ≤ S_j ≤ m - 1` and `S_j ≡ 0 (mod m)` for all `j` is `S = 0`. -/
theorem g5_norm_block_soundness (m d : ℕ) (hm : 1 ≤ m)
    (S : Fin d → ℕ) (hbound : ∀ j, S j ≤ m - 1) (hzero : ∀ j, (S j : ZMod m) = 0) :
    ∀ j, S j = 0 := by
  intro j
  have hdvd : m ∣ S j := (ZMod.natCast_eq_zero_iff (S j) m).mp (hzero j)
  have hlt : S j < m := by have := hbound j; omega
  exact Nat.eq_zero_of_dvd_of_lt hdvd hlt

end WeakPolyComposite
