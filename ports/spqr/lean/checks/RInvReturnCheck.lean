import Spqr.Proofs.RInvFrame
import Mathlib.Tactic.IntervalCases

namespace Spqr.RInvReturnCheck

open WalkM WalkState

def g : Graph := ⟨3, #[(0, 1), (1, 2), (1, 2), (1, 2), (0, 2)]⟩

def child : DfsTree := .node 2 [
  .back 4 0 (.ret 0 .backEdge),
  .back 2 1 (.ret 1 .backEdge),
  .back 3 1 (.ret 1 .backEdge)]

def childOut : DfsOut := .tree 1 (.ret 0 .type1Child) child

def rootOut : DfsOut := .tree 0 .component (.node 1 [childOut])

theorem dfsForest_eq : g.dfsForest [] [] = [.node 0 [rootOut]] := by cbv

def before : WalkState :=
  let s := { WalkState.init g false with stackVerts := #[0, 0, 0] }
  let s := (walkOutPre 0 0 rootOut false).run s |>.2
  let s := { s with firstOccurrence := s.firstOccurrence.set! 0 g.ne }
  let s := { s with stackVerts := s.stackVerts.set! 1 1 }
  let s := (walkOutPre 1 1 childOut false).run s |>.2
  { s with firstOccurrence := s.firstOccurrence.set! 1 g.ne }

def returned : WalkState := after (walkTree child 2) before

def pEntry : TEntry := ⟨2, 1, 1, ([], [9])⟩

theorem pEntry_mem : pEntry ∈ returned.tstack := by
  cbv
  exact List.mem_cons_self ..

theorem pEntry_piece : 9 ∈ returned.entryPieceItems pEntry := by
  cbv
  exact List.mem_cons_self ..

theorem pEntry_vs : Items.vs returned.items 9 = (some 2, some 1) := by cbv

theorem pending_no_parent {j : Nat} (hj : j = 4 ∨ j = 5) (p : Nat) :
    ¬Items.IsParent returned.items p j := by
  by_cases hp : p < 10
  · interval_cases p <;> rcases hj with rfl | rfl <;> cbv <;>
      change ¬(_ : Nat) ∈ (_ : List Nat) <;> decide
  · have hsize : returned.items.size = 10 := by cbv
    have hle : returned.items.size ≤ p := by omega
    simp [Items.IsParent, Items.ch_of_le _ _ hle]

theorem pending_outside {e : Nat} (he : e = 0 ∨ e = 1) :
    ¬Items.EdgeBelow g returned.items 9 e := by
  intro h
  have hroot : ∀ p, ¬Items.IsParent returned.items p (edgeItem g e) :=
    pending_no_parent (by rcases he with rfl | rfl <;> decide)
  have heq := Items.Below.eq_of_no_parent hroot h
  rcases he with rfl | rfl <;> cases heq

theorem not_sepClass : ¬g.SepClass 2 1 1 0 := by
  rintro (h | ⟨x, y, ⟨z, hx⟩, _, hr⟩)
  · cases h
  · simp [Graph.Joins, g] at hx
    rcases hx with ⟨rfl, _⟩ | ⟨_, rfl⟩
    · exact hr.ok_left.2 rfl
    · exact hr.ok_left.1 rfl

theorem returned_not_rInvAt (dfs : DfsData) : ¬returned.RInvAt dfs 1 := by
  intro h
  have hentry := h.entries pEntry pEntry_mem (by decide)
  have hg : returned.g = g := by cbv
  have hmax := hentry.maximal 9 pEntry_piece 2 1 pEntry_vs
  rw [hg] at hmax
  exact not_sepClass (hmax 1 0 (by decide) (by decide)
    (pending_outside (.inr rfl)) (pending_outside (.inl rfl)))

#print axioms pEntry_mem
#print axioms dfsForest_eq
#print axioms pEntry_piece
#print axioms pEntry_vs
#print axioms pending_no_parent
#print axioms pending_outside
#print axioms not_sepClass
#print axioms returned_not_rInvAt

end Spqr.RInvReturnCheck
