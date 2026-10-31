import Spqr
/-!
# The R handoff candidate E3 (`buried_vacuous`) is false

PROOF.md §4.5 proposed, at a child-entry site `walkTree c (d + 1)` (parent `stackVerts[d]`), for every
open entry `t` with `d < t.topDepth`:
`(∀ e e', e < ne → e' < ne → t.edges e → t.edges e' → e = e') ∨ g.TwoAttached (t.edges) t.vStart t.vStart`
(a buried entry is a single edge or attached only at its `vStart`), to make `EntryR.single` vacuous
under `stackVerts.set! (d + 1) c` (`rSide_entry_site`). `checks/WalkInvCheck.lean` (field
`e3.buried_vacuous`, seeds 0..3000 × both modes) refutes it: 6 violations (seeds 666, 989, 1966, both
modes), every one a buried P-bond attached at `{t.vStart, stackVerts[t.topDepth]}`, both stale
vertices of a finished sibling subtree. Below, the edge-minimal block (`isBlock`, 8 vertices,
11 edges) of seed 1966's violation, in the relevant case `t.topDepth = d + 1`: the walk enters the
child `7` of `1` (depth 2) after finishing the subtree `1 → 4 → 6`; the open entry
`⟨6, 3, 2, ([], [21])⟩` holds the P item `21` of the two parallel edges `4 – 6` (edges 6, 8), attached
at `4 = stackVerts[3]` (stale) and `6 = vStart`. Kernel-checked: `decide +kernel` runs the library
walk (`walkOuts`/`walkOutPre`) to the site. The restated field is `e3.stab_single` (see
`WalkInvCheck.Extra.e3Check`): for `t.topDepth = d + 1`, `TwoAttached (t.edges) t.vStart c ∨ ∀ e e' ∈
t.edges, SepClass t.vStart c e e'` — exactly `EntryR.single` in the updated state (here the two
parallel edges share the vertex `4 ∉ {6, 7}`, so they are one `SepClass 6 7` class).
-/
namespace WalkInvCheck.E3False
open Spqr WalkM

instance : Inhabited DfsOut := ⟨.back 0 0 .selfLoop⟩
instance : Inhabited DfsTree := ⟨.node 0 []⟩
instance (items : Items) (p c : ItemId) : Decidable (items.IsParent p c) := by
  unfold Items.IsParent; infer_instance
instance (g : Graph) (e v : Nat) : Decidable (g.Inc e v) := by
  unfold Graph.Inc; infer_instance

deriving instance DecidableEq for Spqr.TEntry

def childOf : DfsOut → DfsTree
  | .tree _ _ c => c
  | .back .. => default

/-- One frame of the library walk towards a child-entry site: `walkTree`'s `stackVerts` update at
`t = .node v outs`, `walkOuts` of the first `k` outs and `walkOut`'s prefix (`walkOutPre`,
`firstOccurrence`) of `outs[k]`; returns that out's subtree. After the last frame the state is the
one `walkTree child (d + 1)` starts from (before its own `stackVerts.set!`). -/
def frame (t : DfsTree) (d k : Nat) : WalkM DfsTree := do
  match t with
  | .node v outs =>
    modify fun s => { s with stackVerts := s.stackVerts.set! d v }
    let hv ← walkOuts v d (outs.take k) false
    let o := outs[k]!
    let _ ← walkOutPre v d o hv
    modify fun s => { s with firstOccurrence := s.firstOccurrence.set! d s.g.ne }
    pure (childOf o)

def g : Graph := ⟨8, #[(5, 1), (5, 0), (4, 5), (5, 7), (6, 3), (5, 3), (4, 6), (1, 7), (6, 4), (0, 6), (4, 1)]⟩

deriving instance DecidableEq for Spqr.Graph
deriving instance BEq for Spqr.DfsOut, Spqr.DfsTree

/-- The first tree of `g.dfsForest [] []` (`dfsForest` does not reduce in the kernel; the equality is
checked by `#guard` below). -/
def tree : DfsTree :=
  .node 0 [.tree 1 .component (.node 5 [
    .tree 0 (.ret 0 .type1Child) (.node 1 [
      .tree 10 (.ret 0 .type2Child) (.node 4 [
        .tree 6 (.ret 0 .type2Child) (.node 6 [
          .back 9 0 (.ret 0 .backEdge),
          .tree 4 (.ret 1 .type1Child) (.node 3 [.back 5 5 (.ret 1 .backEdge)]),
          .back 8 4 (.ret 3 .backEdge)]),
        .back 2 5 (.ret 1 .backEdge)]),
      .tree 7 (.ret 1 .type1Child) (.node 7 [.back 3 5 (.ret 1 .backEdge)])])])]

#guard g.dfsForest [] [] == [tree, .node 2 []]

/-- The site: root `0`, first out (to `5`), first out of `5` (to `1`), second out of `1` (the tree
edge `7` to the child `7`, after the subtree `1 → 4 → 6`). -/
def site : WalkState :=
  ((do let t₁ ← frame tree 0 0; let t₂ ← frame t₁ 1 0; let _ ← frame t₂ 2 1 :
    WalkM Unit).run (WalkState.init g false)).2

theorem mem_of_below_closed {items : Items} (S : List ItemId)
    (hS : ∀ p ∈ S, ∀ c ∈ items.ch p, c ∈ S) {a i : ItemId} (h : items.Below a i) (ha : a ∈ S) :
    i ∈ S := by
  induction h with
  | refl => exact ha
  | tail _ hpc ih => exact hS _ ih _ hpc

/-- The entry: `vStart = 6`, `topDepth = 3 = d + 1`, the P item `21`. -/
def tE : TEntry := ⟨6, 3, 2, ([], [21])⟩

theorem site_tstack : site.tstack =
    [⟨1, 2, 4, ([], [2])⟩, ⟨4, 2, 4, ([], [19])⟩, ⟨4, 1, 3, ([], [11])⟩, ⟨4, 3, 3, ([], [5])⟩, tE,
     ⟨6, 1, 1, ([], [20])⟩, ⟨6, 0, 0, ([18], [])⟩, ⟨6, 4, 0, ([], [7])⟩, ⟨5, 1, 0, ([], [6])⟩] := by
  decide +kernel

theorem e3_buried_vacuous_false :
    site.stackVerts[2]! = 1 ∧ site.stackVerts[3]! = 4 ∧
    ∃ t ∈ site.tstack, t.topDepth = 2 + 1 ∧
      ¬ ((∀ e e', e < site.g.ne → e' < site.g.ne → t.edges site.g site.items e →
            t.edges site.g site.items e' → e = e') ∨
        site.g.TwoAttached (t.edges site.g site.items) t.vStart t.vStart) := by
  refine ⟨by decide +kernel, by decide +kernel, tE, by rw [site_tstack]; simp, rfl, ?_⟩
  have hg : site.g = g := by decide +kernel
  have hS : ∀ p ∈ [21, 17, 15], ∀ c ∈ Items.ch site.items p, c ∈ [21, 17, 15] := by decide +kernel
  have hP8 : Items.IsParent site.items 21 (edgeItem g 8) := by decide +kernel
  have hP6 : Items.IsParent site.items 21 (edgeItem g 6) := by decide +kernel
  rw [hg]
  generalize site.items = I at hS hP8 hP6 ⊢
  have h21 : (21 : ItemId) ∈ tE.spans.1 ++ tE.spans.2 := by decide
  have hE8 : tE.edges g I 8 := ⟨21, h21, Relation.ReflTransGen.single hP8⟩
  have hE6 : tE.edges g I 6 := ⟨21, h21, Relation.ReflTransGen.single hP6⟩
  rintro (hs | hA)
  · exact absurd (hs 8 6 (by decide) (by decide) hE8 hE6) (by decide)
  · have hnE : ¬ tE.edges g I 2 := by
      rintro ⟨i, hi, hb⟩
      have hi' : i = 21 := by simpa [tE] using hi
      subst hi'
      exact absurd (mem_of_below_closed [21, 17, 15] hS hb (by decide)) (by decide)
    have hA' : ∀ v e e', e < g.ne → e' < g.ne → tE.edges g I e → ¬ tE.edges g I e' →
        g.Inc e v → g.Inc e' v → v = tE.vStart ∨ v = tE.vStart := hA
    have h46 := hA' 4 8 2 (by decide) (by decide) hE8 hnE (by decide) (by decide)
    have hv : tE.vStart = 6 := rfl
    rw [hv] at h46
    omega

#print axioms e3_buried_vacuous_false

end WalkInvCheck.E3False
