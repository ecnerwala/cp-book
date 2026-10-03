import Spqr.PlanarEmbed
import Spqr.Spec

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

end PlanarSpqrTree

end Spqr
