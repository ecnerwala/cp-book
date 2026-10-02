import Spqr.PlanarEmbed
import Spqr.PlanarRelabelProj
import Spqr.PlanarRotSpec
import Spqr.Build
import Spqr.Correctness
import Spqr.PlanarInv
import Spqr.PlanarEmbedSteps
import Spqr.PlanarEmbedFold
import Spqr.RelabelChildShape
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

/-- Relabel bookkeeping: the segment of `neRotAdj` belonging to node `i` is the `layoutRot` of its
type, computed by `planarRelabel` when the node was laid out (`PlanarRotSpec.lean`, from the
admitted fold characterization `planarRelabel_rot_spec`). -/
theorem neRotAdj_segment (g : Graph) (ternarize : Bool) (vertOrder edgeOrder : List Nat) (i : Nat)
    (hi : i < (g.planarTree ternarize vertOrder edgeOrder).size) :
    ∃ (edgeVes : List Nat) (mapRot : Nat → Array (Option Nat)) (capVe : Nat),
      (g.planarTree ternarize vertOrder edgeOrder).neRotAdj.extract
          (4 * ((g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.neRange i).1)
          (4 * ((g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.neRange i).2) =
        layoutRot ((g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.type i)
          ((g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.nVerts i)
          ((g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.neRange i).1
          ((g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.neRange i).2 edgeVes mapRot capVe :=
  PlanarRot.neRotAdj_segment' g ternarize vertOrder edgeOrder i hi

theorem planarTree_shape (g : Graph) (ternarize : Bool) (vertOrder edgeOrder : List Nat) (i : Nat)
    (hi : i < (g.planarTree ternarize vertOrder edgeOrder).size) :
    (g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.Shape i := by
  have hwf := spqrTree_wf g ternarize vertOrder edgeOrder
  rw [← planarRelabel_proj] at hwf
  exact hwf.shape i hi

theorem skeleton_length (t : SpqrTree) (i : Nat) : (t.skeleton i).length = t.nEdges i := by
  simp [SpqrTree.skeleton, SpqrTree.nodeEdgesOf, SpqrTree.nEdges]

/-- The S case of `nodePlanar_sound`: the skeleton is the cycle (`Shape`), the local rotation is
the cycle layout (`neRotAdj_segment`, `layoutRot_S_shift`), and the cycle layout is a planar
embedding (`cycleRot_isPlanarEmbedding`). -/
theorem nodePlanar_sound_S (g : Graph) (ternarize : Bool) (vertOrder edgeOrder : List Nat) (i : Nat)
    (hi : i < (g.planarTree ternarize vertOrder edgeOrder).size)
    (hS : (g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.type i = .S) :
    IsPlanarEmbedding ((g.planarTree ternarize vertOrder edgeOrder).localSkeleton i)
      ((g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.nVerts i)
      ((g.planarTree ternarize vertOrder edgeOrder).nodeRot i) := by
  have hsh := planarTree_shape g ternarize vertOrder edgeOrder i hi
  simp only [SpqrTree.Shape, hS] at hsh
  obtain ⟨⟨s, e⟩, hnv⟩ : ∃ p, (g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.nvRange i = p := ⟨_, rfl⟩
  rw [hnv] at hsh
  simp only at hsh
  obtain ⟨hn3, hsk⟩ := hsh
  have hnv' : (g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.nVerts i = e - s := by unfold SpqrTree.nVerts; rw [hnv]
  have hsk' : (g.planarTree ternarize vertOrder edgeOrder).localSkeleton i = cycleEdges (e - s) := by
    unfold PlanarSpqrTree.localSkeleton
    rw [hsk, hnv]
    simp only [cycleEdges, List.map_cons, List.map_map]
    refine List.cons_eq_cons.mpr ⟨?_, ?_⟩
    · simp only [Prod.mk.injEq]; omega
    · apply List.map_congr_left
      intro k hk
      simp only [Function.comp, Prod.mk.injEq]
      omega
  have hrot : (g.planarTree ternarize vertOrder edgeOrder).nodeRot i = cycleRot (e - s) := by
    obtain ⟨ev, mr, cv, hseg⟩ := neRotAdj_segment g ternarize vertOrder edgeOrder i hi
    have hlen := skeleton_length (g.planarTree ternarize vertOrder edgeOrder).toSpqrTree i
    rw [hsk] at hlen
    simp only [List.length_cons, List.length_map, List.length_range] at hlen
    unfold PlanarSpqrTree.nodeRot
    obtain ⟨⟨neSt, neEn⟩, hne⟩ : ∃ p, (g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.neRange i = p := ⟨_, rfl⟩
    rw [hne] at hseg ⊢
    unfold SpqrTree.nEdges at hlen
    rw [hne] at hlen
    simp only at hlen
    obtain ⟨k, rfl⟩ : ∃ k, neEn = neSt + k := ⟨neEn - neSt, by omega⟩
    have hk : k = e - s := by omega
    subst hk
    rw [hS, hnv'] at hseg
    simp only [hseg]
    exact layoutRot_S_shift (e - s) neSt (by omega) ev mr cv
  rw [hsk', hrot, hnv']
  exact cycleRot_isPlanarEmbedding _ (by omega)

/-- The P case of `nodePlanar_sound`: the skeleton is a bond (`Shape`), the local rotation is
the bond layout (`neRotAdj_segment`, `layoutRot_P_shift`), and the bond layout is a planar
embedding (`bondRot_isPlanarEmbedding`). -/
theorem nodePlanar_sound_P (g : Graph) (ternarize : Bool) (vertOrder edgeOrder : List Nat) (i : Nat)
    (hi : i < (g.planarTree ternarize vertOrder edgeOrder).size)
    (hP : (g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.type i = .P) :
    IsPlanarEmbedding ((g.planarTree ternarize vertOrder edgeOrder).localSkeleton i)
      ((g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.nVerts i)
      ((g.planarTree ternarize vertOrder edgeOrder).nodeRot i) := by
  have hsh := planarTree_shape g ternarize vertOrder edgeOrder i hi
  simp only [SpqrTree.Shape, hP] at hsh
  obtain ⟨⟨s, e⟩, hnv⟩ : ∃ p, (g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.nvRange i = p := ⟨_, rfl⟩
  rw [hnv] at hsh
  simp only at hsh
  obtain ⟨hn2, hk3, hall⟩ := hsh
  have hnv' : (g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.nVerts i = 2 := by unfold SpqrTree.nVerts; rw [hnv]; exact hn2
  have hlen := skeleton_length (g.planarTree ternarize vertOrder edgeOrder).toSpqrTree i
  have hsk' : (g.planarTree ternarize vertOrder edgeOrder).localSkeleton i = bondEdges ((g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.nEdges i) := by
    unfold PlanarSpqrTree.localSkeleton bondEdges
    rw [List.eq_replicate_iff]
    refine ⟨by rw [List.length_map, hlen], fun b hb => ?_⟩
    rw [List.mem_map] at hb
    obtain ⟨p, hp, rfl⟩ := hb
    rw [hall p hp, hnv]
    simp
  have hrot : (g.planarTree ternarize vertOrder edgeOrder).nodeRot i = bondRot ((g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.nEdges i) := by
    obtain ⟨ev, mr, cv, hseg⟩ := neRotAdj_segment g ternarize vertOrder edgeOrder i hi
    unfold PlanarSpqrTree.nodeRot
    obtain ⟨⟨neSt, neEn⟩, hne⟩ : ∃ p, (g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.neRange i = p := ⟨_, rfl⟩
    rw [hne] at hseg ⊢
    have hk : (g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.nEdges i = neEn - neSt := by unfold SpqrTree.nEdges; rw [hne]
    rw [hk] at hk3 ⊢
    obtain ⟨k, rfl⟩ : ∃ k, neEn = neSt + k := ⟨neEn - neSt, by omega⟩
    rw [Nat.add_sub_cancel_left]
    rw [hP, hnv'] at hseg
    simp only [hseg]
    exact layoutRot_P_shift k neSt ev mr cv
  rw [hsk', hrot, hnv']
  exact bondRot_isPlanarEmbedding _ (by omega)

/-- The R case of `nodePlanar_sound`. Admitted; plan: when the R item is finished, Invariant P
(`InvariantP`, maintained by `planarWalkOut_stackInv`) gives a planar embedding of its piece whose
exposed ends are the four recorded cap matches (`finishMatches`); `planarRelabel` maps it through
`mapRot` into the node's `neRotAdj` segment (`neRotAdj_segment`), and the piece with its cap is
the node's skeleton. -/
theorem nodePlanar_sound_R (g : Graph) (ternarize : Bool) (vertOrder edgeOrder : List Nat) (i : Nat)
    (hi : i < (g.planarTree ternarize vertOrder edgeOrder).size)
    (hR : (g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.type i = .R)
    (h : (g.planarTree ternarize vertOrder edgeOrder).isPlanar i = true) :
    IsPlanarEmbedding ((g.planarTree ternarize vertOrder edgeOrder).localSkeleton i)
      ((g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.nVerts i)
      ((g.planarTree ternarize vertOrder edgeOrder).nodeRot i) := by
  sorry

/-- Soundness of the per-node flag: a planar S/P/R node's local rotation system is a planar
embedding of its skeleton. S and P are `nodePlanar_sound_S` / `nodePlanar_sound_P` (proved modulo
`neRotAdj_segment`); R is `nodePlanar_sound_R` (admitted). -/
theorem nodePlanar_sound (g : Graph) (ternarize : Bool) (vertOrder edgeOrder : List Nat) (i : Nat)
    (hty : (g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.type i = .S ∨
      (g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.type i = .P ∨
      (g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.type i = .R)
    (h : (g.planarTree ternarize vertOrder edgeOrder).isPlanar i = true) :
    IsPlanarEmbedding ((g.planarTree ternarize vertOrder edgeOrder).localSkeleton i)
      ((g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.nVerts i)
      ((g.planarTree ternarize vertOrder edgeOrder).nodeRot i) := by
  have hi : i < (g.planarTree ternarize vertOrder edgeOrder).size := by
    by_contra hge
    have hF : (g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.type i = .F := by
      unfold SpqrTree.type
      rw [Array.getElem?_eq_none (by unfold PlanarSpqrTree.size at hge; omega)]
      rfl
    rcases hty with h' | h' | h' <;> rw [hF] at h' <;> cases h'
  rcases hty with hS | hP | hR
  · exact nodePlanar_sound_S g ternarize vertOrder edgeOrder i hi hS
  · exact nodePlanar_sound_P g ternarize vertOrder edgeOrder i hi hP
  · exact nodePlanar_sound_R g ternarize vertOrder edgeOrder i hi hR h

/-- Completeness of the per-node flag: a node flagged nonplanar has a nonplanar skeleton.
Admitted (Kuratowski-style; plan in PROOF.md §8.3: the `mergePlanarity` nesting obstruction
exhibits a `K₅` / `K₃,₃` subdivision in the skeleton). -/
theorem nodePlanar_complete (g : Graph) (ternarize : Bool) (vertOrder edgeOrder : List Nat) (i : Nat)
    (hi : i < (g.planarTree ternarize vertOrder edgeOrder).size)
    (h : (g.planarTree ternarize vertOrder edgeOrder).isPlanar i = false) :
    ¬ Planar ((g.planarTree ternarize vertOrder edgeOrder).localSkeleton i)
      ((g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.nVerts i) := by
  sorry

/-- At the root (item `0`, the `F` item of `g`, no parent): `edgesBelow 0` is every edge of `g`,
`outerE[0]` is all unset, so `GluedUpTo 0` says the glued rotation is total and agrees with a
planar embedding of `g.edges` under the identity renumbering. Admitted (tree facts about the
root item). -/
theorem glued_root (g : Graph) (ternarize : Bool) (vertOrder edgeOrder : List Nat)
    (hwf : (g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.WF) (s : PlanarSpqrTree.EmbedState)
    (h : (g.planarTree ternarize vertOrder edgeOrder).GluedUpTo g 0 s) :
    IsPlanarEmbedding g.edges.toList g.nv ⟨s.rotAdj⟩ := by
  sorry

/-- Soundness of the glued embedding: it is a planar embedding of `g` (Euler's formula per
component). The reverse-preorder fold is `forM_reverse_range_inv` with the invariant `GluedUpTo`
(`gluedUpTo_init` proved); the per-item steps `embedItem_step_{F,V,Q,leaf,node}` and the root
assembly `glued_root` are admitted (plan in PROOF.md §8.5: `twoSum_planar` for S/P/R with
`nodePlanar_sound`, `oneSum_planar` for V/Q, `disjointUnion_planar` for F). The tree's `WF` comes
from `spqrTree_wf'`, hence the `g.WF` / `OrderOK` hypotheses. -/
theorem planarEmbed_sound (g : Graph) (hg : g.WF) (ternarize : Bool) (vertOrder edgeOrder : List Nat)
    (hvo : OrderOK g.nv vertOrder) (heo : OrderOK g.ne edgeOrder)
    (rs : RotationSystem) (h : (g.planarTree ternarize vertOrder edgeOrder).planarEmbed = some rs) :
    IsPlanarEmbedding g.edges.toList g.nv rs := by
  unfold PlanarSpqrTree.planarEmbed at h
  split at h
  · rename_i hall
    simp only [Option.some.injEq] at h
    subst h
    have hwf : (g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.WF := by
      rw [planarRelabel_proj]; exact spqrTree_wf' g hg ternarize vertOrder edgeOrder hvo heo
    have hsh : (g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.ChildShape := by
      rw [planarRelabel_proj]; exact spqrTree_childShape g ternarize vertOrder edgeOrder
    have hrep : (g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.Represents g := by
      rw [planarRelabel_proj]; exact spqrTree_represents g ternarize vertOrder edgeOrder
    exact glued_root g ternarize vertOrder edgeOrder hwf _
      ((g.planarTree ternarize vertOrder edgeOrder).gluedUpTo_planarEmbed g hwf hsh hrep hall)
  · cases h

/-- The eventual target: the planar SPQR tree yields an embedding iff `g` is planar. `→` is
`planarEmbed_sound`; `←` needs `nodePlanar_complete` plus the fact that a minor (each skeleton is a
minor of `g`) of a planar graph is planar. Admitted. -/
theorem spqrTree_planar (g : Graph) (ternarize : Bool) (vertOrder edgeOrder : List Nat) :
    (g.planarTree ternarize vertOrder edgeOrder).planarEmbed.isSome ↔ Planar g.edges.toList g.nv := by
  sorry

end Spqr
