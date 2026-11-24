import Mathlib.Data.Nat.Choose.Basic
import Mathlib.Tactic

/-!
# Theorem D12 — BBR Margin-Saturation Barrier: Semi-Symmetric Lucas Slice Certificates

Exact finite-field Lucas slice-sum certificates establishing:
1. At `d = 3, n = 28 (mod 6)`:
   - Optimal symmetric `P_3*` (`zeros mod 9 = {0, 1, 5}`):
     `Σ_1 + Σ_10 + Σ_19 + Σ_28 ≡ G(0) (mod 2)` for every `deg(G) ≤ 3`.
   - Single-slice kills (`w = 8` and `w = 16`):
     `C(28 - k, 8 - k) ≡ 1_{k=0} (mod 2)` and `C(28 - k, 16 - k) ≡ 1_{k=0} (mod 2)` for `k = 0..3`.
2. Dual at `d = 3, n = 40 (mod 6)` (`P_2 = I_4(S)`):
   `2 Σ_4 + 2 Σ_8 + Σ_36 + Σ_40 ≡ G(0) (mod 3)` for every `deg(G) ≤ 6`,
   proving `N_gen^{(6)}(3, OR; P_2 = I_4(S)) ≤ 39`.
-/

namespace WeakPolyComposite

/-- Exact `F_2` Lucas slice certificate for `P_3*` (`zeros mod 9 = {0, 1, 5}`) at `d = 3, n = 28`. -/
theorem d12_mod6_d3_n28_opt_cert :
    ∀ k : Fin 4,
      (Nat.choose (28 - k.1) (1 - k.1)
        + Nat.choose (28 - k.1) (10 - k.1)
        + Nat.choose (28 - k.1) (19 - k.1)
        + Nat.choose (28 - k.1) (28 - k.1)) % 2
        = if k.1 = 0 then 1 else 0 := by
  decide

/-- Single-slice `F_2` certificates `w = 8` and `w = 16` at `d = 3, n = 28`. -/
theorem d12_mod6_d3_n28_single_slices :
    (∀ k : Fin 4, Nat.choose (28 - k.1) (8 - k.1) % 2 = if k.1 = 0 then 1 else 0) ∧
    (∀ k : Fin 4, Nat.choose (28 - k.1) (16 - k.1) % 2 = if k.1 = 0 then 1 else 0) := by
  decide

/-- Exact 4-term `F_3` Lucas slice certificate for `P_2 = I_4(S)` at `d = 3, n = 40`. -/
theorem d12_mod6_d3_n40_i4_cert :
    ∀ k : Fin 7,
      (2 * Nat.choose (40 - k.1) (4 - k.1)
        + 2 * Nat.choose (40 - k.1) (8 - k.1)
        + Nat.choose (40 - k.1) (36 - k.1)
        + Nat.choose (40 - k.1) (40 - k.1)) % 3
        = if k.1 = 0 then 1 else 0 := by
  decide

/-- Exact `F_2` Lucas slice certificate for `I_9(S)` at `d = 6, n = 63` (`N_gen(6; P_3 sym) = 62`). -/
theorem d12_mod6_d6_n63_i9_cert :
    ∀ k : Fin 7,
      (Nat.choose (63 - k.1) (9 - k.1)
        + Nat.choose (63 - k.1) (18 - k.1)
        + Nat.choose (63 - k.1) (27 - k.1)
        + Nat.choose (63 - k.1) (36 - k.1)
        + Nat.choose (63 - k.1) (45 - k.1)
        + Nat.choose (63 - k.1) (54 - k.1)
        + Nat.choose (63 - k.1) (63 - k.1)) % 2
        = if k.1 = 0 then 1 else 0 := by
  decide

end WeakPolyComposite
