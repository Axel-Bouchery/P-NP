import Mathlib.Data.ZMod.Basic
import Mathlib.Data.Finset.Basic
import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.NumberTheory.Wilson

/-!
# Core Definitions for Weak Polynomial Representations Modulo Composite Integers

This file defines Boolean inputs `Fin n → Bool`, Hamming weight,
quadratic multilinear polynomials over `ZMod q`, their polar bilinear forms,
and weak representations of `OR_n` modulo `2p` or `pq`.
-/

namespace WeakPolyComposite

/-- Hamming weight of a Boolean vector `x : Fin n → Bool`. -/
def hammingWeight {n : ℕ} (x : Fin n → Bool) : ℕ :=
  (Finset.univ.filter (fun i => x i = true)).card

/-- A Boolean vector is nonzero iff at least one coordinate is `true`. -/
def isNonzero {n : ℕ} (x : Fin n → Bool) : Prop :=
  ∃ i : Fin n, x i = true

/-- Evaluation of a quadratic multilinear polynomial `P_q` with `P_q(0) = 0` over `ZMod q`. -/
structure QuadPolyMod (n q : ℕ) where
  lin : Fin n → ZMod q
  quad : Fin n → Fin n → ZMod q
  quad_symm : ∀ i j, quad i j = quad j i
  quad_diag : ∀ i, quad i i = 0

/-- Evaluate `P` on a Boolean vector `x : Fin n → Bool`. -/
def QuadPolyMod.eval {n q : ℕ} (P : QuadPolyMod n q) (x : Fin n → Bool) : ZMod q :=
  (∑ i : Fin n, if x i then P.lin i else 0) +
  (∑ i : Fin n, ∑ j : Fin n, if i < j ∧ x i ∧ x j then P.quad i j else 0)

/-- Polar bilinear form `B(u, v) = ∑_{i ∈ supp(u), j ∈ supp(v)} quad i j` for disjoint `u, v`. -/
def QuadPolyMod.polar {n q : ℕ} (P : QuadPolyMod n q) (u v : Fin n → Bool) : ZMod q :=
  ∑ i : Fin n, ∑ j : Fin n, if u i ∧ v j then P.quad i j else 0

/-- Weak representation of `OR_n` modulo `p * q` by a pair `(P_p, P_q)`. -/
def WeakRepOR {n p q : ℕ} (Pp : QuadPolyMod n p) (Pq : QuadPolyMod n q) : Prop :=
  ∀ x : Fin n → Bool, isNonzero x → (Pp.eval x ≠ 0 ∨ Pq.eval x ≠ 0)

end WeakPolyComposite
