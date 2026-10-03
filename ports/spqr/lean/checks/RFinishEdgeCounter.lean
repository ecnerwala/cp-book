import Spqr.Proofs.RInvFrame
import Mathlib.Tactic.IntervalCases

/-!
# Kernel-checked counterexamples: child-return settling is bounded by `topDepth`

Both states below are reached by the actual walk (`dfsForest_eq` fixes the DFS tree, the state
is the walk's own prefix).

* `Deep`: after the parent's `finishEdge` the frontier entry `(3, 1)` (bottom `3`, top at depth
  `1`, strictly above the parent `2` at depth `2`) still holds the P piece `{3,1}` of the two
  back edges `3 → 1` while the pending path class `1 - 2 - 3` is a second `{3,1}` class. So the
  settled-at-`curV` output `RInvAt dfs curV` of the (old) `finishEdge_rInvAt` is false even
  though the base below the split is settled (`base_entryR`): the parent's loops settle only
  the entries whose top is at depth `≥ d`.
* `Base`: before a child return, a base entry `(1, 0)` of the grandparent `1` (not the parent
  `3`) holds the closed S piece `{6,1}` of the first ear while the tree edge `1 → 3` and the back
  edge `1 → 6` are still outside it; `RReturn.entries` with the `vStart ≠ parent` exemption is
  false, the `d ≤ topDepth` restriction (`RReturn` now) is what the walk keeps.
-/

namespace Spqr.RFinishEdgeCounter

open WalkM WalkState

namespace Deep

def g : Graph := ⟨4, #[(0, 1), (1, 2), (2, 3), (3, 1), (3, 1), (3, 0)]⟩

def leaf : DfsTree := .node 3 [
  .back 5 0 (.ret 0 .backEdge),
  .back 3 1 (.ret 1 .backEdge),
  .back 4 1 (.ret 1 .backEdge)]

def leafOut : DfsOut := .tree 2 (.ret 0 .type2Child) leaf

def midOut : DfsOut := .tree 1 (.ret 0 .type1Child) (.node 2 [leafOut])

def rootOut : DfsOut := .tree 0 .component (.node 1 [midOut])

theorem dfsForest_eq : g.dfsForest [] [] = [.node 0 [rootOut]] := by cbv

/-- The state at vertex `2` (depth `2`) just before the child `3` is walked. -/
def before : WalkState :=
  let s := WalkState.init g false
  let s := after (walkOutPre 0 0 rootOut false) s
  let s := { s with firstOccurrence := s.firstOccurrence.set! 0 g.ne }
  let s := { s with stackVerts := s.stackVerts.set! 1 1 }
  let s := after (walkOutPre 1 1 midOut false) s
  let s := { s with firstOccurrence := s.firstOccurrence.set! 1 g.ne }
  let s := { s with stackVerts := s.stackVerts.set! 2 2 }
  let s := after (walkOutPre 2 2 leafOut false) s
  { s with firstOccurrence := s.firstOccurrence.set! 2 g.ne }

def returned : WalkState := after (walkTree leaf 3) before

/-- After the parent's `finishEdge` of the tree edge `2 → 3`. -/
def finished : WalkState := after (finishEdge 2 2 leafOut before.tstack.length false) returned

theorem before_tstack : before.tstack = [⟨1, 1, 0, ([], [2])⟩] := by cbv

def pEntry : TEntry := ⟨3, 1, 1, ([], [11])⟩

theorem pEntry_mem : pEntry ∈ finished.tstack := by
  cbv
  exact List.mem_cons_of_mem _ (List.mem_cons_of_mem _ (List.mem_cons_self ..))

theorem pEntry_piece : 11 ∈ finished.entryPieceItems pEntry := by
  cbv
  exact List.mem_cons_self ..

theorem pEntry_vs : Items.vs finished.items 11 = (some 3, some 1) := by cbv

theorem pending_no_parent {j : Nat} (hj : j = 6 ∨ j = 10) (p : Nat) :
    ¬Items.IsParent finished.items p j := by
  by_cases hp : p < 12
  · interval_cases p <;> rcases hj with rfl | rfl <;> cbv <;>
      change ¬(_ : Nat) ∈ (_ : List Nat) <;> decide
  · have hsize : finished.items.size = 12 := by cbv
    have hle : finished.items.size ≤ p := by omega
    simp [Items.IsParent, Items.ch_of_le _ _ hle]

theorem pending_outside {e : Nat} (he : e = 1 ∨ e = 5) :
    ¬Items.EdgeBelow g finished.items 11 e := by
  intro h
  have hroot : ∀ p, ¬Items.IsParent finished.items p (edgeItem g e) :=
    pending_no_parent (by rcases he with rfl | rfl <;> decide)
  have heq := Items.Below.eq_of_no_parent hroot h
  rcases he with rfl | rfl <;> cases heq

theorem reach_two {x y : Nat} (h : g.Reach (fun z => z ≠ 3 ∧ z ≠ 1) x y) : x = 2 → y = 2 := by
  induction h with
  | refl _ => exact id
  | tail _ hadj hok ih =>
    intro hx
    have hy := ih hx
    subst hy
    obtain ⟨e, he⟩ := hadj
    by_cases he6 : e < 6
    · interval_cases e <;> simp [Graph.Joins, g] at he <;> omega
    · have hn : g.edges[e]? = none := Array.getElem?_eq_none_iff.2 (by simp [g]; omega)
      simp [Graph.Joins, hn] at he

theorem not_sepClass : ¬g.SepClass 3 1 1 5 := by
  rintro (h | ⟨x, y, ⟨z, hx⟩, ⟨w, hy⟩, hr⟩)
  · cases h
  · have hl := hr.ok_left
    have hr' := hr.ok_right
    simp [Graph.Joins, g] at hx hy
    have hx2 : x = 2 := by rcases hx with ⟨h1, h2⟩ | ⟨h1, h2⟩ <;> omega
    have hy0 : y = 0 := by rcases hy with ⟨h1, h2⟩ | ⟨h1, h2⟩ <;> omega
    subst hx2 hy0
    exact absurd (reach_two hr rfl) (by decide)

/-- The settled-at-`curV` output of the old `finishEdge_rInvAt` fails for every `dfs`. -/
theorem finished_not_rInvAt (dfs : DfsData) : ¬finished.RInvAt dfs 2 := by
  intro h
  have hentry := h.entries pEntry pEntry_mem (by decide)
  have hg : finished.g = g := by cbv
  have hmax := hentry.maximal 11 pEntry_piece 3 1 pEntry_vs
  rw [hg] at hmax
  exact not_sepClass (hmax 1 5 (by decide) (by decide)
    (pending_outside (.inl rfl)) (pending_outside (.inr rfl)))

theorem returned_base :
    returned.tstack.drop (returned.tstack.length - before.tstack.length) = before.tstack := by
  cbv

theorem vert_ch : Items.ch returned.items 2 = [] := by cbv

theorem vert_no_edges (e : Nat) : ¬(⟨1, 1, 0, ([], [2])⟩ : TEntry).edges returned.g returned.items e := by
  rintro ⟨i, hi, hb⟩
  have hi2 : i = 2 := by simpa using hi
  subst hi2
  have hg : returned.g = g := by cbv
  rw [hg] at hb
  rcases Relation.ReflTransGen.cases_head hb with h | ⟨c, hc, -⟩
  · simp [edgeItem, g] at h
    have h' : (2 : Nat) = 5 + e := h
    omega
  · simp [Items.IsParent, vert_ch] at hc

/-- The base below the split is settled: the hypothesis of the revised contract holds. -/
theorem base_entryR (dfs : DfsData) :
    ∀ t ∈ returned.tstack.drop (returned.tstack.length - before.tstack.length),
      t.vStart ≠ 2 → returned.EntryR dfs t := by
  rw [returned_base, before_tstack]
  intro t ht _
  have ht' : t = ⟨1, 1, 0, ([], [2])⟩ := by simpa using ht
  subst ht'
  have hp : returned.entryPieceItems ⟨1, 1, 0, ([], [2])⟩ = [] := by cbv
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_⟩
  · rw [hp]; exact ⟨List.nodup_nil, by simp, by simp, by simp, by simp, by simp⟩
  · rw [hp]; simp
  · exact fun a b e₁ e₂ _ _ _ h1 _ => absurd h1 (vert_no_edges e₁)
  · exact fun _ e e' _ _ h1 _ => absurd h1 (vert_no_edges e)
  · exact fun a b _ _ o _ _ => .inr (.inl fun e _ he => vert_no_edges e he)
  · exact fun a b _ _ o _ _ _ => .inr (.inl fun e _ he => vert_no_edges e he)

end Deep

namespace Base

def g : Graph := ⟨7, #[(6, 1), (3, 2), (4, 6), (3, 1), (6, 1), (1, 4), (2, 6), (1, 4)]⟩

def leaf : DfsTree := .node 2 [.back 6 6 (.ret 0 .backEdge)]

def leafOut : DfsOut := .tree 1 (.ret 0 .type1Child) leaf

def out3 : DfsOut := .tree 3 (.ret 0 .type1Child) (.node 3 [leafOut])

def out5 : DfsOut := .tree 5 (.ret 0 .type1Child)
  (.node 4 [.back 2 6 (.ret 0 .backEdge), .back 7 1 (.ret 1 .backEdge)])

def rootOut : DfsOut := .tree 4 .component
  (.node 1 [out5, out3, .back 0 6 (.ret 0 .backEdge)])

theorem dfsForest_eq :
    g.dfsForest [6, 5, 4, 3, 2, 1, 0] [4, 5, 6, 7, 0, 1, 2, 3] =
      [.node 6 [rootOut], .node 5 [], .node 0 []] := by cbv

/-- The state at vertex `3` (depth `2`) just before the child `2` is walked; the first ear of
`1` (`1 - 4 - 6`) is already closed into the S piece `17` of the entry `(1, 0)`. -/
def before : WalkState :=
  let s := WalkState.init g false
  let s := { s with stackVerts := s.stackVerts.set! 0 6 }
  let s := after (walkOutPre 6 0 rootOut false) s
  let s := { s with firstOccurrence := s.firstOccurrence.set! 0 g.ne }
  let s := { s with stackVerts := s.stackVerts.set! 1 1 }
  let (hv, s) := (walkOut 1 1 out5 false).run s
  let s := after (walkOutPre 1 1 out3 hv) s
  let s := { s with firstOccurrence := s.firstOccurrence.set! 1 g.ne }
  let s := { s with stackVerts := s.stackVerts.set! 2 3 }
  let s := after (walkOutPre 3 2 leafOut false) s
  { s with firstOccurrence := s.firstOccurrence.set! 2 g.ne }

def returned : WalkState := after (walkTree leaf 3) before

def sEntry : TEntry := ⟨1, 0, 0, ([17], [])⟩

theorem returned_base :
    returned.tstack.drop (returned.tstack.length - before.tstack.length) = before.tstack := by
  cbv

theorem sEntry_mem :
    sEntry ∈ returned.tstack.drop (returned.tstack.length - before.tstack.length) := by
  cbv
  exact List.mem_cons_of_mem _ (List.mem_cons_self ..)

theorem sEntry_piece : 17 ∈ returned.entryPieceItems sEntry := by
  cbv
  exact List.mem_cons_self ..

theorem sEntry_vs : Items.vs returned.items 17 = (some 6, some 1) := by cbv

theorem pending_no_parent {j : Nat} (hj : j = 8 ∨ j = 11) (p : Nat) :
    ¬Items.IsParent returned.items p j := by
  by_cases hp : p < 18
  · interval_cases p <;> rcases hj with rfl | rfl <;> cbv <;>
      change ¬(_ : Nat) ∈ (_ : List Nat) <;> decide
  · have hsize : returned.items.size = 18 := by cbv
    have hle : returned.items.size ≤ p := by omega
    simp [Items.IsParent, Items.ch_of_le _ _ hle]

theorem pending_outside {e : Nat} (he : e = 0 ∨ e = 3) :
    ¬Items.EdgeBelow g returned.items 17 e := by
  intro h
  have hroot : ∀ p, ¬Items.IsParent returned.items p (edgeItem g e) :=
    pending_no_parent (by rcases he with rfl | rfl <;> decide)
  have heq := Items.Below.eq_of_no_parent hroot h
  rcases he with rfl | rfl <;> cases heq

/-- The back edge `1 → 6` joins the pair itself: a class of its own. -/
theorem not_sepClass : ¬g.SepClass 6 1 0 3 := by
  rintro (h | ⟨x, y, ⟨z, hx⟩, _, hr⟩)
  · cases h
  · simp [Graph.Joins, g] at hx
    rcases hx with ⟨rfl, _⟩ | ⟨_, rfl⟩
    · exact hr.ok_left.1 rfl
    · exact hr.ok_left.2 rfl

/-- `RReturn.entries` with only the parent exempt is false for every `dfs`. -/
theorem not_base_entries (dfs : DfsData) :
    ¬∀ t ∈ returned.tstack.drop (returned.tstack.length - before.tstack.length),
      t.vStart ≠ 3 → returned.EntryR dfs t := by
  intro h
  have hentry := h sEntry sEntry_mem (by decide)
  have hg : returned.g = g := by cbv
  have hmax := hentry.maximal 17 sEntry_piece 6 1 sEntry_vs
  rw [hg] at hmax
  exact not_sepClass (hmax 0 3 (by decide) (by decide)
    (pending_outside (.inl rfl)) (pending_outside (.inr rfl)))

end Base

#print axioms Deep.dfsForest_eq
#print axioms Deep.finished_not_rInvAt
#print axioms Deep.base_entryR
#print axioms Base.dfsForest_eq
#print axioms Base.not_base_entries

end Spqr.RFinishEdgeCounter
