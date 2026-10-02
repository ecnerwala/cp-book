import Spqr.PlanarEmbed
import Spqr.PlanarRelabelProj
import Spqr.Build
import Spqr.Spec

/-!
# Specification of the planar variant

`planarWalk_proj` / `planarRelabel_proj` (proved in `PlanarWalkProj` / `PlanarRelabelProj`): the
planar variant is a conservative extension (its base component is the ordinary algorithm). `nodePlanar_sound` / `nodePlanar_complete`: the per-node
flags are right. `planarEmbed_sound` / `planarEmbed_isSome_iff` / `spqrTree_planar`: the glued
embedding is a planar embedding of `g`, and exists iff `g` is planar.
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

/-- `planarEmbed` is defined exactly when every item is flagged planar. -/
theorem planarEmbed_isSome_iff :
    t.planarEmbed.isSome ↔ ∀ i (h : i < t.nodePlanar.size), t.nodePlanar[i] = true := by
  unfold planarEmbed
  split <;> rename_i h
  · simpa [Array.all_eq_true] using h
  · simp only [Option.isSome_none, Bool.false_eq_true, false_iff]
    simpa [Array.all_eq_true] using h

end PlanarSpqrTree

/-- The planar tree of `g`: the ordinary tree plus `nodePlanar` / `neRotAdj`. -/
abbrev Graph.planarTree (g : Graph) (ternarize : Bool) (vertOrder edgeOrder : List Nat) : PlanarSpqrTree :=
  g.planarSpqrTree ternarize vertOrder edgeOrder

/-- Soundness of the per-node flag: a planar S/P/R node's local rotation system is a planar
embedding of its skeleton. Admitted; plan in PROOF.md §8.4 (S: the cycle layout `layoutRot .S` is
the two-face embedding; P: the bond layout is the `k`-face embedding of `k` parallel edges;
R: the walk's `node_planarity` matches, i.e. the t-stack invariant of §8.2 at the moment the R item
is finished, mapped through `mapRot`). -/
theorem nodePlanar_sound (g : Graph) (ternarize : Bool) (vertOrder edgeOrder : List Nat) (i : Nat)
    (hty : (g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.type i = .S ∨
      (g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.type i = .P ∨
      (g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.type i = .R)
    (h : (g.planarTree ternarize vertOrder edgeOrder).isPlanar i = true) :
    IsPlanarEmbedding ((g.planarTree ternarize vertOrder edgeOrder).localSkeleton i)
      ((g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.nVerts i)
      ((g.planarTree ternarize vertOrder edgeOrder).nodeRot i) := by
  sorry

/-- Completeness of the per-node flag: a node flagged nonplanar has a nonplanar skeleton.
Admitted (Kuratowski-style; plan in PROOF.md §8.3: the `mergePlanarity` nesting obstruction
exhibits a `K₅` / `K₃,₃` subdivision in the skeleton). -/
theorem nodePlanar_complete (g : Graph) (ternarize : Bool) (vertOrder edgeOrder : List Nat) (i : Nat)
    (hi : i < (g.planarTree ternarize vertOrder edgeOrder).size)
    (h : (g.planarTree ternarize vertOrder edgeOrder).isPlanar i = false) :
    ¬ Planar ((g.planarTree ternarize vertOrder edgeOrder).localSkeleton i)
      ((g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.nVerts i) := by
  sorry

/-- Soundness of the glued embedding: it is a planar embedding of `g` (Euler's formula per
component). Admitted; plan in PROOF.md §8.5 (2-sum of planar embeddings along twin edges). -/
theorem planarEmbed_sound (g : Graph) (ternarize : Bool) (vertOrder edgeOrder : List Nat)
    (rs : RotationSystem) (h : (g.planarTree ternarize vertOrder edgeOrder).planarEmbed = some rs) :
    IsPlanarEmbedding g.edges.toList g.nv rs := by
  sorry

/-- The eventual target: the planar SPQR tree yields an embedding iff `g` is planar. `→` is
`planarEmbed_sound`; `←` needs `nodePlanar_complete` plus the fact that a minor (each skeleton is a
minor of `g`) of a planar graph is planar. Admitted. -/
theorem spqrTree_planar (g : Graph) (ternarize : Bool) (vertOrder edgeOrder : List Nat) :
    (g.planarTree ternarize vertOrder edgeOrder).planarEmbed.isSome ↔ Planar g.edges.toList g.nv := by
  sorry

end Spqr
