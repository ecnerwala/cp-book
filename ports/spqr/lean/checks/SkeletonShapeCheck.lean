import Spqr.ItemSpec

namespace Spqr.SkeletonShapeCheck

def k4 : Graph := ⟨4, #[(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]⟩

def k4Items : Items := (k4.walk false (k4.dfsForest [] [])).items

theorem k4_shape : 11 < k4Items.size ∧ k4Items.type 11 = .R ∧
    (k4Items.virtualEdges 11).length = 5 ∧
    ((k4Items.ch 11).filter fun c => k4Items.type c = .V).length = 2 := by
  cbv

def triangle : Graph := ⟨3, #[(0, 1), (1, 2), (0, 2)]⟩

def triangleItems : Items := (triangle.walk false (triangle.dfsForest [] [])).items

theorem triangle_shape : 7 < triangleItems.size ∧ triangleItems.type 7 = .S ∧
    (triangleItems.virtualEdges 7).length = 2 ∧
    ((triangleItems.ch 7).filter fun c => triangleItems.type c = .V).length = 1 := by
  cbv

#print axioms k4_shape
#print axioms triangle_shape

end Spqr.SkeletonShapeCheck
