import Spqr.PlanarEmbed
import Spqr.Spec
import Spqr.PieceSep

/-!
# Per-node specification data of the planar variant

`localSkeleton`, `nodeRot`, `isPlanar`: the skeleton, local rotation system and planarity flag of
one node, shared by `PlanarSpec` (`nodePlanar_sound`) and the gluing fold (`embedItem_step_node_faces`).
-/

namespace Spqr

namespace PlanarSpqrTree

variable (t : PlanarSpqrTree)

/-- The skeleton of node `i` with node-vertices renumbered from `0`. -/
def localSkeleton (i : Nat) : List (Nat × Nat) :=
  (t.toSpqrTree.skeleton i).map fun p => (p.1 - (t.toSpqrTree.nvRange i).1, p.2 - (t.toSpqrTree.nvRange i).1)

/-- The local rotation system of node `i`, on the quarter-edges of its own node-edges. -/
def nodeRot (i : Nat) : RotationSystem :=
  let (neSt, neEn) := t.toSpqrTree.neRange i
  ⟨(t.neRotAdj.extract (4 * neSt) (4 * neEn)).map (·.map (· - 4 * neSt))⟩

def isPlanar (i : Nat) : Bool := t.nodePlanar[i]?.getD false

/-- Every node is flagged planar, so each node in range is. -/
theorem isPlanar_of_all {t : PlanarSpqrTree} (hall : t.nodePlanar.all id = true) {i : Nat}
    (hi : i < t.nodePlanar.size) : t.isPlanar i = true := by
  have h' : ∀ j (h : j < t.nodePlanar.size), t.nodePlanar[j] = true := by
    simpa [Array.all_eq_true] using hall
  unfold isPlanar
  rw [Array.getElem?_eq_getElem hi]
  exact h' i hi

/-- The `neRotAdj` entries of an `R` node stay inside the node's own quarter-edge segment (so
`nodeRot i`, which subtracts `4 * neSt`, records them faithfully). -/
def NodeRotClosed : Prop :=
  ∀ i, i < t.size → t.toSpqrTree.type i = .R → ∀ ta tb,
    4 * (t.toSpqrTree.neRange i).1 ≤ ta → ta < 4 * (t.toSpqrTree.neRange i).2 →
    t.neRotAdj[ta]? = some (some tb) →
    4 * (t.toSpqrTree.neRange i).1 ≤ tb ∧ tb < 4 * (t.toSpqrTree.neRange i).2

/-- `ta` is the quarter-edge `(side 1, dir 0)` of a non-cap edge of node `i` ending at node-vertex
`nv` whose facing quarter-edge is a `(side 0, dir 1)` quarter-edge: the corner where the executable
node step splices the vertex item of `nv`. -/
def CornerAt (i nv ta : Nat) : Prop :=
  4 * ((t.toSpqrTree.neRange i).1 + 1) ≤ ta ∧ ta < 4 * (t.toSpqrTree.neRange i).2 ∧ ta % 4 = 2 ∧
    (∃ d, t.nodeEdges[ta / 4]? = some d ∧ d.nvs.2 = nv) ∧
    ∃ tb, t.neRotAdj[ta]? = some (some tb) ∧ tb % 4 = 1

/-- Corner structure of the local rotation systems used by the node step: cap endpoints have no
splice corner; every inner node-vertex has exactly one, and its facing quarter-edge comes later
(so the node step's `ta < tb` scan visits it). -/
structure NodeCorners : Prop where
  capEnd : ∀ i nv ta, i < t.size →
    (t.toSpqrTree.type i = .S ∨ t.toSpqrTree.type i = .P ∨ t.toSpqrTree.type i = .R) →
    t.toSpqrTree.CapEnd i nv → ¬ t.CornerAt i nv ta
  inner : ∀ i nv, i < t.size →
    (t.toSpqrTree.type i = .S ∨ t.toSpqrTree.type i = .P ∨ t.toSpqrTree.type i = .R) →
    t.toSpqrTree.NvOf i nv → ¬ t.toSpqrTree.CapEnd i nv →
    ∃ ta, t.CornerAt i nv ta ∧ (∀ tb, t.neRotAdj[ta]? = some (some tb) → ta < tb) ∧
      ∀ ta', t.CornerAt i nv ta' → ta' = ta

end PlanarSpqrTree

end Spqr
