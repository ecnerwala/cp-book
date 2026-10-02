import Spqr.PlanarLayout
import Spqr.PlanarWalk

/-!
# Invariant P and the gluing lemmas

`InvariantP` is the per-entry invariant of the planar walk (PROOF.md §8.2): the piece of a tstack
entry is embedded in the plane with both terminals on the outer face, and the entry's two sides
(`PlSide.bot`, `PlSide.top`) list exactly the exposed ends of the two outer-face boundary walks,
split at the terminals. `StackInv` lifts it to a whole walk state; its preservation by the walk
(`planarWalkOut_stackInv`) is admitted.

`twoSum_planar` / `oneSum_planar` / `disjointUnion_planar` are the graph-level gluing facts used
by `planarEmbed` (S/P/R nodes glue children through twin edges, V items glue blocks at a vertex,
F items collect components). They are stated on `Planar` and admitted here; the explicit splicing
construction with `f = f₁ + f₂ - 2` faces is developed in `Spqr.PlanarGlue`.
-/

namespace Spqr

/-! ### Pieces -/

/-- The piece of a tstack entry: the virtual edges `ves` of `E(t)` with endpoints `ends` in
`[0, nVerts)`, where `top` is the contracted upper DFS stack `T` and `bot` is `vStart`. -/
structure Piece where
  ves : List Nat
  ends : Nat → Nat × Nat
  nVerts : Nat
  bot : Nat
  top : Nat

namespace Piece

variable (P : Piece)

/-- The edge list of the piece, indexed by position in `ves`. -/
def es : List (Nat × Nat) := P.ves.map P.ends

/-- The global quarter-edge `q` (of virtual edge `QE.edge q`) belongs to the piece. -/
def Mem (q : Nat) : Prop := QE.edge q ∈ P.ves

/-- Local index of the global quarter-edge `q` in `P.es`. -/
def loc (q : Nat) : Option Nat := (P.ves.findIdx? (· == QE.edge q)).map fun k => 4 * k + q % 4

end Piece

/-- The exposed ends listed by one side, in boundary order: the outer and inner tree-side ends,
then the outermost and innermost back-edge ends at `T`. -/
def PlSide.ends (s : PlSide) : List Nat :=
  (s.bot.map fun b => [b.1, b.2]).getD [] ++ (s.top.map fun t => [t.ends.1, t.ends.2]).getD []

/-- Side `sd` of the planarity data (`false` = side 0). -/
def Planarity.side (p : Planarity) (sd : Bool) : PlSide := if sd then p.sides.2 else p.sides.1

/-- All exposed ends of both sides. -/
def Planarity.ends (p : Planarity) : List Nat := p.sides.1.ends ++ p.sides.2.ends

/-- Span items of side `sd` of a tstack entry. -/
def spanOf (e : TEntry) (sd : Bool) : List ItemId := if sd then e.spans.2 else e.spans.1

/-- Flip bits of side `sd` of a planar entry (parallel to `spanOf`). -/
def flipsOf (pe : PlEntry) (sd : Bool) : List Bool := if sd then pe.flips.2 else pe.flips.1

/-- The first `k` corners of the face walk of `rs` starting at `q`. -/
def RotationSystem.faceWalk (rs : RotationSystem) (q k : Nat) : List Nat :=
  (List.range k).map fun j => (stepFn rs.faceStep)^[j] q

/-- The side-`sd` boundary walk `w` (a face walk of the piece's embedding) of an entry: it visits
the side's listed ends in order and no other exposed end, and the cap edge of every span item of
the side lies on it, with the side bit of the corner given by the item's flip. -/
structure SideWalk (nv : Nat) (e : TEntry) (pe : PlEntry) (p : Planarity) (P : Piece) (sd : Bool)
    (w : List Nat) : Prop where
  ends_sub : ((p.side sd).ends.filterMap P.loc).Sublist w
  ends_only : ∀ q ∈ p.ends, ∀ l, P.loc q = some l → l ∈ w → q ∈ (p.side sd).ends
  spans : ∀ x ∈ (spanOf e sd).zip (flipsOf pe sd), ∃ l ∈ w, ∃ q, P.loc q = some l ∧
    QE.edge q = x.1 - (1 + nv) ∧ QE.side q = (xor sd x.2).toNat

/-- **Invariant P** (PROOF.md §8.2) for the tstack entry `e` with planarity payload `pe`, data
`p`, piece `P`, embedding `ρ` and outer-face corner `f`, in a walk over `nv` vertices with
quarter-edge matches `qem`:

* *(embedded)* `ρ` is a planar embedding of the piece (with `T` contracted), `qem` on the piece's
  quarter-edges is the restriction of `ρ` to the non-exposed corners, and the unmatched corners
  are exactly the ends listed by the two sides;
* *(outer face)* the face of `f` contains a corner at `vStart` and a corner at `T`, and every
  listed end lies on it; `bot` ends sit at `vStart` and `top` ends at `T`;
* *(boundary, split at the terminals)* the face walk from the first end of each side is that
  side's boundary walk (`SideWalk`);
* *(nesting)* on each side the depths increase inwards and are at least `topDepth`; side 0
  holds a minimal return whenever any back edge is open;
* *(spans)* the span items' caps belong to the piece and the flip lists are parallel to the
  spans. -/
structure InvariantP (nv : Nat) (qem : Qem) (e : TEntry) (pe : PlEntry) (p : Planarity)
    (P : Piece) (ρ : RotationSystem) (f : Nat) : Prop where
  pl : pe.pl = some p
  embedded : IsPlanarEmbedding P.es P.nVerts ρ
  agree : ∀ q r, P.Mem q → qem[q]? = some (some r) →
    ∃ lq lr, P.loc q = some lq ∧ P.loc r = some lr ∧ ρ.get lq = some lr
  exposed : ∀ q, P.Mem q → (qem[q]? = some none ↔ q ∈ p.ends)
  outer_bot : QE.vert P.es f = some P.bot
  outer_top : ∃ q, ρ.SameFaceOrbit f q ∧ QE.vert P.es q = some P.top
  ends_outer : ∀ q ∈ p.ends, ∃ l, P.loc q = some l ∧ ρ.SameFaceOrbit f l
  bot_at : ∀ sd : Bool, ∀ b ∈ (p.side sd).bot, ∃ l, P.loc b.1 = some l ∧ QE.vert P.es l = some P.bot
  top_at : ∀ sd : Bool, ∀ t ∈ (p.side sd).top, ∀ q ∈ [t.ends.1, t.ends.2],
    ∃ l, P.loc q = some l ∧ QE.vert P.es l = some P.top
  side_walk : ∀ sd : Bool, ∀ q₀ ∈ (p.side sd).ends.head?,
    ∃ l₀ k, P.loc q₀ = some l₀ ∧ SideWalk nv e pe p P sd (ρ.faceWalk l₀ k)
  nested : ∀ sd : Bool, ∀ t ∈ (p.side sd).top, e.topDepth ≤ t.depths.1 ∧ t.depths.1 ≤ t.depths.2
  minimal : (p.sides.1.top.isSome ∨ p.sides.2.top.isSome) →
    ∃ t ∈ p.sides.1.top, t.depths.1 = e.topDepth
  spans_mem : ∀ it ∈ e.spans.1 ++ e.spans.2, it - (1 + nv) ∈ P.ves
  flips_len : pe.flips.1.length = e.spans.1.length ∧ pe.flips.2.length = e.spans.2.length

/-- Every tstack entry of a planar walk state that is still flagged planar satisfies Invariant P
for some piece, embedding and outer face. -/
def StackInv (s : PlanarWalkState) : Prop :=
  ∀ x ∈ s.base.tstack.zip s.aux.plStack, ∀ p, x.2.pl = some p →
    ∃ P ρ f, InvariantP s.base.g.nv s.aux.qem x.1 x.2 p P ρ f

/-- The planar walk preserves Invariant P (PROOF.md §8.3–8.4: `makeEdgePlanarity` creates a
planar piece, `mergePlanarity` glues along the shared terminal and its nesting test is the only
obstruction, `finishMatches` / `unwrapPlanarity` replace a piece by its cap and back, `closeSide`
/ `pruneSide` close returns). Admitted. -/
theorem planarWalkOut_stackInv (v d : Nat) (o : DfsOut) (hasVert : Bool) (s : PlanarWalkState)
    (h : StackInv s) : StackInv ((planarWalkOut v d o hasVert).run s).2 := by
  sorry

/-! ### Gluing two embedded graphs -/

/-- Vertex map used to glue a second graph onto a first one on `n₁` vertices: a vertex `w` of
the second graph listed in `ids` as `(v₁, w)` becomes `v₁`, every other vertex is appended after
the first graph (keeping the order of the remaining vertices). -/
def glueVert (n₁ : Nat) (ids : List (Nat × Nat)) (w : Nat) : Nat :=
  match ids.find? (·.2 == w) with
  | some p => p.1
  | none => n₁ + w - (ids.filter (·.2 < w)).length

/-- The union of `es₁` and `es₂` with the vertices of `es₂` mapped through `glueVert`. -/
def glueEdges (n₁ : Nat) (es₁ es₂ : List (Nat × Nat)) (ids : List (Nat × Nat)) : List (Nat × Nat) :=
  es₁ ++ es₂.map fun p => (glueVert n₁ ids p.1, glueVert n₁ ids p.2)

/-- The 2-sum of `es₁` and `es₂` along the twin edges `es₁[e₁] = (u₁, v₁)` and `es₂[e₂] = (u₂, v₂)`:
both twins are deleted and `u₂ ↦ u₁`, `v₂ ↦ v₁`. -/
def twoSumEdges (n₁ : Nat) (es₁ es₂ : List (Nat × Nat)) (e₁ e₂ : Nat) (u₁ v₁ u₂ v₂ : Nat) :
    List (Nat × Nat) :=
  glueEdges n₁ (es₁.eraseIdx e₁) (es₂.eraseIdx e₂) [(u₁, u₂), (v₁, v₂)]

/-- The 1-sum of `es₁` and `es₂` identifying `v₂ ↦ v₁`. -/
def oneSumEdges (n₁ : Nat) (es₁ es₂ : List (Nat × Nat)) (v₁ v₂ : Nat) : List (Nat × Nat) :=
  glueEdges n₁ es₁ es₂ [(v₁, v₂)]

/-- The disjoint union of `es₁` and `es₂`. -/
def disjointUnionEdges (n₁ : Nat) (es₁ es₂ : List (Nat × Nat)) : List (Nat × Nat) :=
  glueEdges n₁ es₁ es₂ []

/-- **2-sum gluing.** Two planar embedded graphs sharing a twin edge glue to a planar graph:
splice the rotations at `u` and at `v` (the corners facing the twin in each embedding become
adjacent) and delete the twins; the face count is `f₁ + f₂ - 2` and Euler's formula follows.
Admitted; the explicit splice is `Spqr.PlanarGlue`. -/
theorem twoSum_planar (n₁ n₂ : Nat) (es₁ es₂ : List (Nat × Nat)) (rs₁ rs₂ : RotationSystem)
    (h₁ : IsPlanarEmbedding es₁ n₁ rs₁) (h₂ : IsPlanarEmbedding es₂ n₂ rs₂) (e₁ e₂ u₁ v₁ u₂ v₂ : Nat)
    (he₁ : es₁[e₁]? = some (u₁, v₁)) (he₂ : es₂[e₂]? = some (u₂, v₂)) (hu₁ : u₁ ≠ v₁) (hu₂ : u₂ ≠ v₂) :
    Planar (twoSumEdges n₁ es₁ es₂ e₁ e₂ u₁ v₁ u₂ v₂) (n₁ + n₂ - 2) := by
  sorry

/-- **1-sum gluing.** Identifying one vertex of two planar embedded graphs keeps planarity:
splice the rotation of `v₂` into that of `v₁` at an outer corner; the face count drops by one
and so does the number of components. Admitted. -/
theorem oneSum_planar (n₁ n₂ : Nat) (es₁ es₂ : List (Nat × Nat)) (rs₁ rs₂ : RotationSystem)
    (h₁ : IsPlanarEmbedding es₁ n₁ rs₁) (h₂ : IsPlanarEmbedding es₂ n₂ rs₂) (v₁ v₂ : Nat)
    (hv₁ : v₁ < n₁) (hv₂ : v₂ < n₂) :
    Planar (oneSumEdges n₁ es₁ es₂ v₁ v₂) (n₁ + n₂ - 1) := by
  sorry

/-- The disjoint union of two planar embedded graphs is planar (components are embedded
separately in `EulerFormula`). Admitted. -/
theorem disjointUnion_planar (n₁ n₂ : Nat) (es₁ es₂ : List (Nat × Nat)) (rs₁ rs₂ : RotationSystem)
    (h₁ : IsPlanarEmbedding es₁ n₁ rs₁) (h₂ : IsPlanarEmbedding es₂ n₂ rs₂) :
    Planar (disjointUnionEdges n₁ es₁ es₂) (n₁ + n₂) := by
  sorry

/-! ### Renumbering a laid-out node's edges from `0` -/

theorem rotS_shift (n neSt ne r : Nat) (hn : 1 ≤ n) : rotS n neSt ne r = 4 * neSt + rotS n 0 ne r := by
  have h : neSt + n - 1 = neSt + (n - 1) := by omega
  simp only [rotS, h, Nat.zero_add]
  split_ifs <;> omega

theorem rotP_shift (k neSt ne r : Nat) : rotP k neSt ne r = 4 * neSt + rotP k 0 ne r := by
  unfold rotP
  dsimp only
  split_ifs <;> omega

/-- An S layout with its node-edges renumbered from `0` is the cycle layout. -/
theorem layoutRot_S_shift (n neSt : Nat) (hn : 2 ≤ n) (ev : List Nat) (mr : Nat → Array (Option Nat)) (cv : Nat) :
    (⟨(layoutRot .S n neSt (neSt + n) ev mr cv).map (·.map (· - 4 * neSt))⟩ : RotationSystem) = cycleRot n := by
  unfold cycleRot
  congr 1
  have hs0 := layoutRot_S_size n 0 hn [] (fun _ => #[]) 0
  rw [Nat.zero_add] at hs0
  apply Array.ext
  · rw [Array.size_map, layoutRot_S_size n neSt hn, hs0]
  · intro i h1 h2
    rw [Array.size_map, layoutRot_S_size n neSt hn] at h1
    obtain ⟨ne, r, hr, rfl⟩ := exists_decomp i
    have hg0 := layoutRot_S_get n 0 hn [] (fun _ => #[]) 0 ne r (by omega) hr
    rw [Nat.zero_add] at hg0
    rw [← Option.some_inj, ← Array.getElem?_eq_getElem, ← Array.getElem?_eq_getElem, Array.getElem?_map,
      layoutRot_S_get n neSt hn ev mr cv ne r (by omega) hr, hg0]
    simp only [Option.map_some, rotS_shift n neSt ne r (by omega), Nat.add_sub_cancel_left]

/-- A P layout with its node-edges renumbered from `0` is the bond layout. -/
theorem layoutRot_P_shift (k neSt : Nat) (ev : List Nat) (mr : Nat → Array (Option Nat)) (cv : Nat) :
    (⟨(layoutRot .P 2 neSt (neSt + k) ev mr cv).map (·.map (· - 4 * neSt))⟩ : RotationSystem) = bondRot k := by
  unfold bondRot
  congr 1
  have hs0 := layoutRot_P_size k 0 [] (fun _ => #[]) 0
  rw [Nat.zero_add] at hs0
  apply Array.ext
  · rw [Array.size_map, layoutRot_P_size k neSt, hs0]
  · intro i h1 h2
    rw [Array.size_map, layoutRot_P_size k neSt] at h1
    obtain ⟨ne, r, hr, rfl⟩ := exists_decomp i
    have hg0 := layoutRot_P_get k 0 [] (fun _ => #[]) 0 ne r (by omega) hr
    rw [Nat.zero_add] at hg0
    rw [← Option.some_inj, ← Array.getElem?_eq_getElem, ← Array.getElem?_eq_getElem, Array.getElem?_map,
      layoutRot_P_get k neSt ev mr cv ne r (by omega) hr, hg0]
    simp only [Option.map_some, rotP_shift k neSt ne r, Nat.add_sub_cancel_left]

end Spqr
