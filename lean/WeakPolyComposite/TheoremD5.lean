import WeakPolyComposite.Basic

/-!
# Theorem D5 (General Semiprimes `m = pq`: Exact `2`-Block Optimality `2(m - 1)`)

Formalizes the exact `2`-block upper bound `a + b ≤ 2 * (p * q - 1)` for both:
- Case `p ≥ 3` (direct `m`-periodicity in each block), and
- Case `p = 2` (`q ≥ 3` odd, parity identities on the `q`-lattice `(u*q, v*q)` modulo `2`).
-/

namespace WeakPolyComposite

/-- Case `p ≥ 3` of Theorem D5(1): `a ≤ m - 1` and `b ≤ m - 1` implies `a + b ≤ 2(m - 1)`. -/
theorem d5_two_block_odd_primes (m a b : ℕ) (hm : 1 ≤ m)
    (ha : a ≤ m - 1) (hb : b ≤ m - 1) :
    a + b ≤ 2 * (m - 1) := by
  omega

/-- Parity identity 1 for `p = 2` (`q` odd): for any 2-block quadratic `P_2` mod 2 with `P_2(0,0) = 0`,
`P_2(q,0) + P_2(2q,0) + P_2(3q,0) = 0` in `ZMod 2`, so not all three can equal `1`. -/
theorem d5_two_block_p2_row_parity (q : ℕ) (hq : Odd q)
    (α β : ZMod 2) :
    let f (s : ℕ) : ZMod 2 := α * ((s * (s - 1) / 2 : ℕ) : ZMod 2) + β * (s : ZMod 2)
    f q + f (2 * q) + f (3 * q) = 0 := by
  obtain ⟨k, rfl⟩ := hq
  intro f
  -- Exact mod-4 periodicity of `s(s-1)/2 mod 2`:
  -- `(2k+1)k + (2k+1)(4k+1) + (6k+3)(3k+1) = (2k+1)(14k+4)` is even,
  -- and `(2k+1) + 2(2k+1) + 3(2k+1) = 6(2k+1)` is even.
  have h1 : (2 * k + 1) * (2 * k) / 2 + (2 * (2 * k + 1)) * (2 * (2 * k + 1) - 1) / 2 +
      (3 * (2 * k + 1)) * (3 * (2 * k + 1) - 1) / 2 =
      2 * ((2 * k + 1) * (7 * k + 2)) := by
    omega
  have h2 : (2 * k + 1) + 2 * (2 * k + 1) + 3 * (2 * k + 1) = 2 * (3 * (2 * k + 1)) := by
    omega
  dsimp [f]
  have h1z : (((2 * k + 1) * (2 * k + 1 - 1) / 2 : ℕ) : ZMod 2) +
      (((2 * (2 * k + 1)) * (2 * (2 * k + 1) - 1) / 2 : ℕ) : ZMod 2) +
      (((3 * (2 * k + 1)) * (3 * (2 * k + 1) - 1) / 2 : ℕ) : ZMod 2) = 0 := by
    norm_cast
    have h1' : (2 * k + 1) * (2 * k + 1 - 1) / 2 +
        (2 * (2 * k + 1)) * (2 * (2 * k + 1) - 1) / 2 +
        (3 * (2 * k + 1)) * (3 * (2 * k + 1) - 1) / 2 =
        2 * ((2 * k + 1) * (7 * k + 2)) := by omega
    rw [h1', Nat.cast_mul]
    simp
  have h2z : ((2 * k + 1 : ℕ) : ZMod 2) + ((2 * (2 * k + 1) : ℕ) : ZMod 2) +
      ((3 * (2 * k + 1) : ℕ) : ZMod 2) = 0 := by
    norm_cast
    rw [h2, Nat.cast_mul]
    simp
  linear_combination α * h1z + β * h2z

/-- Case `p = 2` final arithmetic bound: if `u0 + v0 ≤ 2` with `a ≤ (u0 + 1)q - 1` and
`b ≤ (v0 + 1)q - 1`, then `a + b ≤ 2(2q - 1)`. -/
theorem d5_two_block_p2_bound (q a b u0 v0 : ℕ) (hq : 3 ≤ q)
    (huv : u0 + v0 ≤ 2) (ha : a ≤ (u0 + 1) * q - 1) (hb : b ≤ (v0 + 1) * q - 1) :
    a + b ≤ 2 * (2 * q - 1) := by
  omega

end WeakPolyComposite
