import Spqr.Planar
import Mathlib.Data.Finset.Card

/-!
# 2-sum gluing of rotation systems

Definitions for PROOF.md §8.5: the 2-sum `G₁ ⊕ G₂` of two edge lists along virtual edges
`e₁`, `e₂`, the spliced rotation system, and an abstract orbit count used to relate the
face/vertex orbits of the splice to those of the summands.
-/

namespace Spqr

/-- `f` permutes the finite set `S` and fixes everything outside it. -/
structure IsPermOn (f : Nat → Nat) (S : Finset Nat) : Prop where
  maps : ∀ q ∈ S, f q ∈ S
  inj : ∀ q ∈ S, ∀ r ∈ S, f q = f r → q = r
  fix : ∀ q, q ∉ S → f q = q

/-- `q` and `r` lie on the same orbit of `f` (equivalence closure of `f a = b`). -/
def SameOrbit (f : Nat → Nat) (q r : Nat) : Prop :=
  Relation.EqvGen (fun a b => f a = b) q r

/-- Forward reachability under `f`. -/
def Reach (f : Nat → Nat) (q r : Nat) : Prop :=
  Relation.ReflTransGen (fun a b => f a = b) q r

open Classical in
/-- The orbit of `q` inside `S`. -/
noncomputable def orbit (f : Nat → Nat) (S : Finset Nat) (q : Nat) : Finset Nat :=
  S.filter (SameOrbit f q)

/-- Number of orbits of `f` on `S`. -/
noncomputable def orbitCount (f : Nat → Nat) (S : Finset Nat) : Nat :=
  (S.image (orbit f S)).card

/-- `f` with `p` skipped: the predecessor of `p` now maps to `f p`, and `p` is fixed. -/
def delete (f : Nat → Nat) (p a : Nat) : Nat :=
  if a = p then p else if f a = p then f p else f a

/-- `f` with the images of `x` and `y` exchanged. -/
def swapImg (f : Nat → Nat) (x y a : Nat) : Nat :=
  if a = x then f y else if a = y then f x else f a

/-- Disjoint union of two steps, `f` on `S` and `g` elsewhere. -/
def unionStep (f g : Nat → Nat) (S : Finset Nat) (q : Nat) : Nat :=
  if q ∈ S then f q else g q

/-- `f` with the points of `P` skipped: each is fixed, and its predecessor jumps over it. -/
def skip (f : Nat → Nat) (P : List Nat) (a : Nat) : Nat :=
  if a ∈ P then a else if f a ∈ P then f (f a) else f a

namespace RotationSystem

variable (rs : RotationSystem)

/-- Totalised step `q ↦ rs (q ^^^ c)`: `c = 3` walks faces, `c = 1` walks vertices. -/
def stepC (c : Nat) : Nat → Nat := stepFn fun q => rs.get (q ^^^ c)

end RotationSystem

/-- Index of edge `f ≠ e` after deleting edge `e`. -/
def delIdx (e f : Nat) : Nat := if f < e then f else f - 1

/-- Original index of the `f`-th edge remaining after deleting edge `e`. -/
def insIdx (e f : Nat) : Nat := if f < e then f else f + 1

/-- Two edge lists with a virtual edge in each: `es₁[e₁] = (u₁, v₁)`, `es₂[e₂] = (u₂, v₂)`.
The 2-sum deletes both virtual edges and identifies `u₂ ↦ u₁`, `v₂ ↦ v₁`. -/
structure TwoSum where
  es₁ : List (Nat × Nat)
  n₁ : Nat
  e₁ : Nat
  u₁ : Nat
  v₁ : Nat
  es₂ : List (Nat × Nat)
  n₂ : Nat
  e₂ : Nat
  u₂ : Nat
  v₂ : Nat

namespace TwoSum

variable (T : TwoSum)

def m₁ : Nat := T.es₁.length
def m₂ : Nat := T.es₂.length

/-- Vertex of `G₂` in the 2-sum: the terminals are identified, the rest shifted by `n₁`. -/
def vert₂ (w : Nat) : Nat :=
  if w = T.u₂ then T.u₁ else if w = T.v₂ then T.v₁ else T.n₁ + w

def edges : List (Nat × Nat) :=
  T.es₁.eraseIdx T.e₁ ++ (T.es₂.eraseIdx T.e₂).map fun p => (T.vert₂ p.1, T.vert₂ p.2)

def nVerts : Nat := T.n₁ + T.n₂

/-- Number of quarter-edges of `G₁ − e₁`. -/
def off : Nat := 4 * (T.m₁ - 1)

/-- Glued index of a quarter-edge of `G₁` (edge `≠ e₁`). -/
def qe₁ (q : Nat) : Nat := 4 * delIdx T.e₁ (q / 4) + q % 4

/-- Glued index of a quarter-edge of `G₂` (edge `≠ e₂`). -/
def qe₂ (q : Nat) : Nat := T.off + 4 * delIdx T.e₂ (q / 4) + q % 4

/-- Original `G₁` quarter-edge of a glued index `< off`. -/
def pre₁ (r : Nat) : Nat := 4 * insIdx T.e₁ (r / 4) + r % 4

/-- Original `G₂` quarter-edge of a glued index `≥ off`. -/
def pre₂ (r : Nat) : Nat := 4 * insIdx T.e₂ ((r - T.off) / 4) + (r - T.off) % 4

/-- The rotation of the splice. A corner `q ↔ 4·e₁ + k` of `G₁` is rewired to the corner
`q ↔ rs₂ (4·e₂ + (k ^^^ 1))` of `G₂`, and symmetrically. -/
def spliceFn (rs₁ rs₂ : RotationSystem) (r : Nat) : Option Nat :=
  if r < T.off then
    (rs₁.get (T.pre₁ r)).bind fun s =>
      if s / 4 = T.e₁ then (rs₂.get (4 * T.e₂ + (s % 4 ^^^ 1))).map T.qe₂
      else some (T.qe₁ s)
  else
    (rs₂.get (T.pre₂ r)).bind fun s =>
      if s / 4 = T.e₂ then (rs₁.get (4 * T.e₁ + (s % 4 ^^^ 1))).map T.qe₁
      else some (T.qe₂ s)

/-- The spliced rotation system of the 2-sum. -/
def splice (rs₁ rs₂ : RotationSystem) : RotationSystem :=
  ⟨((List.range (4 * (T.m₁ + T.m₂ - 2))).map (T.spliceFn rs₁ rs₂)).toArray⟩

/-- Vertex connectivity through an edge list. -/
def _root_.Spqr.EdgesConn (es : List (Nat × Nat)) (a b : Nat) : Prop :=
  Relation.ReflTransGen (fun x y => (x, y) ∈ es ∨ (y, x) ∈ es) a b

/-- Well-formedness of a 2-sum of two planar embeddings. `deg` says the virtual edge is not
glued to itself (its ends have degree `≥ 2`); `face` says one of the virtual edges separates
two distinct faces, `conn` that one of them is not a bridge (in a planar embedding these are
equivalent; both are automatic when a skeleton is 2-connected). -/
structure WF (rs₁ rs₂ : RotationSystem) : Prop where
  e₁ : T.es₁[T.e₁]? = some (T.u₁, T.v₁)
  e₂ : T.es₂[T.e₂]? = some (T.u₂, T.v₂)
  uv₁ : T.u₁ ≠ T.v₁
  uv₂ : T.u₂ ≠ T.v₂
  deg₁ : ∀ k, k < 4 → ∀ s ∈ rs₁.get (4 * T.e₁ + k), s / 4 ≠ T.e₁
  deg₂ : ∀ k, k < 4 → ∀ s ∈ rs₂.get (4 * T.e₂ + k), s / 4 ≠ T.e₂
  emb₁ : IsPlanarEmbedding T.es₁ T.n₁ rs₁
  emb₂ : IsPlanarEmbedding T.es₂ T.n₂ rs₂
  face : ¬rs₁.SameFaceOrbit (4 * T.e₁) (4 * T.e₁ + 2) ∨
    ¬rs₂.SameFaceOrbit (4 * T.e₂) (4 * T.e₂ + 2)
  conn : EdgesConn (T.es₁.eraseIdx T.e₁) T.u₁ T.v₁ ∨ EdgesConn (T.es₂.eraseIdx T.e₂) T.u₂ T.v₂

end TwoSum

end Spqr
