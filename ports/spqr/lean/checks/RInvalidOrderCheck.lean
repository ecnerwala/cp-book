import Spqr.Proofs.RItems

namespace Spqr.RInvalidOrderCheck

set_option maxRecDepth 4096
set_option maxHeartbeats 2000000

def k4 : Graph := ⟨4, #[(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]⟩

def output : SpqrTree := k4.spqrTree false [] [6]

def items : Items := (k4.walk false (k4.dfsForest [] [6])).items

theorem wf : k4.WF := by
  simp [Graph.WF, k4]

theorem bad_order : ¬OrderOK k4.ne [6] := by
  simp [OrderOK, Graph.ne, k4]

theorem r_item : 11 < items.size ∧ items.type 11 = .R ∧
    items.vs 11 = (some 0, none) := by
  cbv

theorem not_rSkel3 : ¬Items.RSkel3 k4 items 11 := by
  rintro ⟨s, t, h, -⟩
  rw [r_item.2.2] at h
  simp at h

theorem r_node : 3 < output.size ∧ output.type 3 = .R ∧ output.nVerts 3 = 1 := by
  cbv

theorem not_threeConnected :
    ¬SpqrTree.ThreeConnected (output.nVerts 3)
      ((output.skeleton 3).map fun p =>
        (p.1 - (output.nvRange 3).1, p.2 - (output.nvRange 3).1)) := by
  intro h
  have := h.1
  rw [r_node.2.2] at this
  omega

#print axioms r_node
#print axioms not_threeConnected
#print axioms wf
#print axioms bad_order
#print axioms r_item
#print axioms not_rSkel3

end Spqr.RInvalidOrderCheck
