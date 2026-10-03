import Spqr.Proofs.RInvFrame

/-!
# R content at the `finishEdge` and `walkTree` sites

The block-content facts the R-side site admissions need and no record exports, stated at the
site states so that the walk backbone can carry them: `RCloseContent` at the pre-state of a
returning tree-edge `finishEdge` (the settled loop-1 frontier, the `feS₂` top, the R-branch
fields of every `.R` loop-1 iterate, the type-1 vertex close) and `REntryContent` at a child
`walkTree` entry. Each field is mirrored one-to-one by a `rc.<field>` check of
`checks/WalkInvCheck/RContent.lean` (0 violations, seeds 0..3000 × both modes); PROOF.md §4.5.
-/

namespace Spqr
open WalkM
namespace WalkState
variable {s : WalkState} {dfs : DfsData}

/-- The `k`-th loop-1 iterate of the tree edge `o` finished at depth `d` from the site state `s`. -/
abbrev rl1Iter (d : Nat) (o : DfsOut) (s : WalkState) (k : Nat) : WalkState :=
  iter (loop1Body d s.stackDir[d]!) k (ceS₁ o.dest d o.e (feS₀ d o s))

/-- The `RBranch` fields beyond the shape, `mid`, `interior` and `cur_ne`. -/
structure RBranchFields (s : WalkState) (d : Nat) (cur nxt : TEntry) : Prop where
  cur_piece : ∃ i, i ∈ s.entryPieceItems cur
  cur_vs : ∀ i ∈ s.entryPieceItems cur,
    Items.vs s.items i = (some cur.vStart, some s.stackVerts[d]!) ∨
      Items.vs s.items i = (some s.stackVerts[d]!, some cur.vStart)
  nxt_ne : ∃ e, e < s.g.ne ∧ nxt.edges s.g s.items e
  proper : ∃ e, e < s.g.ne ∧ ¬s.rU cur nxt e
  nxt_touch_top : s.g.Touches (nxt.edges s.g s.items) s.stackVerts[d]!
  nxt_touch_bot : s.g.Touches (nxt.edges s.g s.items) nxt.vStart
  nxt_no_cu : ∀ e, nxt.edges s.g s.items e → ¬s.g.Joins e cur.vStart s.stackVerts[d]!

/-- The children of the type-1 vertex close: the spans of `c`, `py`, `vy` as the two
`mergeTstackTops` and the `retarget` concatenate them. -/
def vItems (c py vy : TEntry) : List ItemId :=
  ((c.spans.1 ++ py.spans.1) ++ vy.spans.1) ++ (vy.spans.2 ++ (py.spans.2 ++ c.spans.2))

/-- The edges closed by the type-1 vertex close. -/
def vU (s : WalkState) (c py vy : TEntry) (e : Nat) : Prop :=
  c.edges s.g s.items e ∨ py.edges s.g s.items e ∨ vy.edges s.g s.items e

/-- The sub-pieces of the type-1 vertex close: the merged items other than vertex items. -/
def vPieceItems (s : WalkState) (c py vy : TEntry) : List ItemId :=
  (vItems c py vy).filter fun i => decide (Items.type s.items i ≠ .V)

/-- The state after the P-check of the tree-edge branch. -/
def feP (curV d : Nat) (o : DfsOut) (s : WalkState) : WalkState :=
  after (Spqr.finishP curV (o.cls.lowval d) o.cls.isType1) (feS₂ d o s)

end WalkState

/-- The Hopcroft–Tarjan content of `RCloseShape` (`single`, `maximal`, `type1`, `bond`, `type2`),
the fields a stack shape cannot give. -/
structure RHT (g : Graph) (d : DfsData) (P : Pieces) (U : Nat → Prop) (s t : Nat) : Prop where
  single : ∀ e e', e < g.ne → e' < g.ne → U e → U e' → g.SepClass s t e e'
  maximal : ∀ i, i < P.k → ∀ e e', e < g.ne → e' < g.ne → ¬P.Mem i e → ¬P.Mem i e' →
    g.SepClass (P.x i) (P.y i) e e'
  type1 : ∀ a b, P.SkelPair g U s t a b → d.Anc a b → ∀ o ∈ d.outs b,
    o.cls = .ret (d.depth a) .type1Child →
      P.LaminarWith (fun e => e < g.ne ∧ U e) (d.EndIn o.dest · g)
  bond : ∀ a b e₁ e₂, e₁ ≠ e₂ → g.Joins e₁ a b → g.Joins e₂ a b → U e₁ → U e₂ →
    ∃ i, P.Mem i e₁ ∧ P.Mem i e₂
  type2 : ∀ a b, P.SkelPair g U s t a b → d.Type2Pair a b → ∀ o ∈ d.outs a, o.isTree = true →
    d.Anc o.dest b → P.LaminarWith (fun e => e < g.ne ∧ U e) (g.SepClass a b o.e)

theorem RCloseShape.ht {g : Graph} {d : DfsData} {P : Pieces} {U : Nat → Prop} {s t : Nat}
    (h : RCloseShape g d P U s t) : RHT g d P U s t :=
  ⟨h.single, h.maximal, h.type1, h.bond, h.type2⟩

namespace WalkState
variable {s : WalkState} {dfs : DfsData}

/-- The R content at the pre-state `s` of `finishEdge curV d o origTstack hasVert` for a
returning tree edge `o` (`o.cls.lowval d < d`, `stackVerts[d] = curV`), all with respect to the
DFS `dfs` of `s.g`. -/
structure RCloseContent (dfs : DfsData) (curV d : Nat) (o : DfsOut) (origTstack : Nat)
    (hasVert : Bool) (s : WalkState) : Prop where
  /-- The entries loop 1 leaves above the base (other than the top) are settled, unless exempt. -/
  settled : ∀ t ∈ (feS₁ d o s).tstack.tail.take ((feS₁ d o s).tstack.length - 1 - origTstack),
    d ≤ t.topDepth → t.vStart ≠ curV → (feS₁ d o s).EntryR dfs t
  /-- Without a vertex entry, a `feS₂` top topping out at `d` and not starting at `curV` is settled. -/
  s2_top : hasVert = false → ∀ c ∈ (feS₂ d o s).tstack.head?, c.topDepth = d → c.vStart ≠ curV →
    (feS₂ d o s).EntryR dfs c
  /-- Without a vertex entry, no entry after the P-check owns an edge of `vertItem curV`. -/
  vert_own : hasVert = false → ∀ t ∈ (feP curV d o s).tstack, ∀ e,
    t.edges (feP curV d o s).g (feP curV d o s).items e →
    ¬ Items.EdgeBelow (feP curV d o s).g (feP curV d o s).items (vertItem curV) e
  /-- At every `.R` iterate of loop 1 (head `cur`, next `nxt`): the child is the head's bottom or
  interior to it (`RBranch.mid`). -/
  l1_mid : ∀ k, (∀ j, j ≤ k → result (loop1Cond d) (rl1Iter d o s j) = true) →
    l1Ty d s.stackDir[d]! (rl1Iter d o s k) = .R →
    ∀ cur nxt rest, (rl1Iter d o s k).tstack = cur :: nxt :: rest →
    (rl1Iter d o s k).stackVerts[d + 1]! = cur.vStart ∨
      (rl1Iter d o s k).g.Interior (cur.edges (rl1Iter d o s k).g (rl1Iter d o s k).items)
        (rl1Iter d o s k).stackVerts[d + 1]!
  /-- ... the interval/saturation fields of the R branch. -/
  l1_fields : ∀ k, (∀ j, j ≤ k → result (loop1Cond d) (rl1Iter d o s j) = true) →
    l1Ty d s.stackDir[d]! (rl1Iter d o s k) = .R →
    ∀ cur nxt rest, (rl1Iter d o s k).tstack = cur :: nxt :: rest →
    (rl1Iter d o s k).RBranchFields d cur nxt
  /-- ... `cur` and `nxt` are settled and edge-disjoint. -/
  l1_top : ∀ k, (∀ j, j ≤ k → result (loop1Cond d) (rl1Iter d o s j) = true) →
    l1Ty d s.stackDir[d]! (rl1Iter d o s k) = .R →
    ∀ cur nxt rest, (rl1Iter d o s k).tstack = cur :: nxt :: rest →
    (rl1Iter d o s k).RTop dfs cur nxt
  /-- The type-1 vertex close with loop 2 fired (`feSingle = false`, stack `c :: py :: vy :: _` at
  `feS₂`): `c`, `py` and `vy` are settled, and the union `vU c py vy` with pieces `vPieceItems` has the
  HT content at `(curV, stackVerts[lowval])`. -/
  v_close : hasVert = true → o.cls.isType1 = true → feSingle d o s = false →
    ∀ c py vy rest, (feS₂ d o s).tstack = c :: py :: vy :: rest →
    (feS₂ d o s).EntryR dfs c ∧ (feS₂ d o s).EntryR dfs py ∧ (feS₂ d o s).EntryR dfs vy ∧
    RHT (feS₂ d o s).g dfs
      (Pieces.ofItems (feS₂ d o s).g (feS₂ d o s).items ((feS₂ d o s).vPieceItems c py vy))
      ((feS₂ d o s).vU c py vy) curV (feS₂ d o s).stackVerts[o.cls.lowval d]!

/-- The R content at a child entry `walkTree c (d + 1)` (parent `stackVerts[d]`, before
`stackVerts.set! (d + 1) c`). -/
structure REntryContent (dfs : DfsData) (d c : Nat) (s : WalkState) : Prop where
  /-- A settled entry topping out at `d + 1` stays settled when `stackVerts[d + 1]` is overwritten
  by `c` (the other entries do not read it: `EntryR.set_stackVerts`). -/
  stab : ∀ t ∈ s.tstack, t.topDepth = d + 1 → s.EntryR dfs t →
    ({ s with stackVerts := s.stackVerts.set! (d + 1) c } : WalkState).EntryR dfs t
  /-- An entry starting at the parent tops out at depth `≤ d`. -/
  parent_top : ∀ t ∈ s.tstack, t.vStart = s.stackVerts[d]! → t.topDepth ≤ d

end WalkState
end Spqr
