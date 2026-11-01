import Spqr.PlanarEmbed
import Spqr.PlanarRelabelProj
import Spqr.PlanarRotSpec
import Spqr.Build
import Spqr.Correctness
import Spqr.PlanarInv
import Spqr.PlanarEmbedSteps
import Spqr.PlanarEmbedFold
import Spqr.PlanarEmbedFacesFold
import Spqr.PlanarNodeSpec
import Spqr.PlanarEmbedRoot
import Spqr.WalkPieceSep
import Spqr.PlanarSoundR
import Spqr.PlanarRelabelBridge
import Spqr.PlanarRelabelCornersR
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
theorem neRotAdj_segment (g : Graph) (hg : g.WF) (ternarize : Bool) (vertOrder edgeOrder : List Nat)
    (hvo : OrderOK g.nv vertOrder) (heo : OrderOK g.ne edgeOrder) (i : Nat)
    (hi : i < (g.planarTree ternarize vertOrder edgeOrder).size) :
    ∃ (edgeVes : List Nat) (mapRot : Nat → Array (Option Nat)) (capVe : Nat),
      (g.planarTree ternarize vertOrder edgeOrder).neRotAdj.extract
          (4 * ((g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.neRange i).1)
          (4 * ((g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.neRange i).2) =
        layoutRot ((g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.type i)
          ((g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.nVerts i)
          ((g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.neRange i).1
          ((g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.neRange i).2 edgeVes mapRot capVe :=
  PlanarRot.neRotAdj_segment' g ternarize vertOrder edgeOrder hg hvo heo i hi

theorem planarTree_shape (g : Graph) (hg : g.WF) (ternarize : Bool) (vertOrder edgeOrder : List Nat)
    (hvo : OrderOK g.nv vertOrder) (heo : OrderOK g.ne edgeOrder) (i : Nat)
    (hi : i < (g.planarTree ternarize vertOrder edgeOrder).size) :
    (g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.Shape i := by
  have hwf := spqrTree_wf g hg ternarize vertOrder edgeOrder hvo heo
  rw [← planarRelabel_proj] at hwf
  exact hwf.shape i hi

theorem skeleton_length (t : SpqrTree) (i : Nat) : (t.skeleton i).length = t.nEdges i := by
  simp [SpqrTree.skeleton, SpqrTree.nodeEdgesOf, SpqrTree.nEdges]

/-- The S case of `nodePlanar_sound`: the skeleton is the cycle (`Shape`), the local rotation is
the cycle layout (`neRotAdj_segment`, `layoutRot_S_shift`), and the cycle layout is a planar
embedding (`cycleRot_isPlanarEmbedding`). -/
theorem nodePlanar_sound_S (g : Graph) (hg : g.WF) (ternarize : Bool) (vertOrder edgeOrder : List Nat)
    (hvo : OrderOK g.nv vertOrder) (heo : OrderOK g.ne edgeOrder) (i : Nat)
    (hi : i < (g.planarTree ternarize vertOrder edgeOrder).size)
    (hS : (g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.type i = .S) :
    IsPlanarEmbedding ((g.planarTree ternarize vertOrder edgeOrder).localSkeleton i)
      ((g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.nVerts i)
      ((g.planarTree ternarize vertOrder edgeOrder).nodeRot i) := by
  have hsh := planarTree_shape g hg ternarize vertOrder edgeOrder hvo heo i hi
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
    obtain ⟨ev, mr, cv, hseg⟩ := neRotAdj_segment g hg ternarize vertOrder edgeOrder hvo heo i hi
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
theorem nodePlanar_sound_P (g : Graph) (hg : g.WF) (ternarize : Bool) (vertOrder edgeOrder : List Nat)
    (hvo : OrderOK g.nv vertOrder) (heo : OrderOK g.ne edgeOrder) (i : Nat)
    (hi : i < (g.planarTree ternarize vertOrder edgeOrder).size)
    (hP : (g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.type i = .P) :
    IsPlanarEmbedding ((g.planarTree ternarize vertOrder edgeOrder).localSkeleton i)
      ((g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.nVerts i)
      ((g.planarTree ternarize vertOrder edgeOrder).nodeRot i) := by
  have hsh := planarTree_shape g hg ternarize vertOrder edgeOrder hvo heo i hi
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
    obtain ⟨ev, mr, cv, hseg⟩ := neRotAdj_segment g hg ternarize vertOrder edgeOrder hvo heo i hi
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

/-- The R case of `nodePlanar_sound`: from the walk-side record `PlanarFinish`
(`planarWalk_planarFinish`) and the relabel-side record `RelabelNodeR`
(`planarRelabelTree_relabelNodeR`) by `nodePlanar_sound_R_of` (reorder the piece's edges to
`Items.ordered` order, relabel its vertices to node-vertex positions). -/
theorem nodePlanar_sound_R (g : Graph) (hg : g.WF) (ternarize : Bool) (vertOrder edgeOrder : List Nat)
    (hvo : OrderOK g.nv vertOrder) (heo : OrderOK g.ne edgeOrder) (i : Nat)
    (hi : i < (g.planarTree ternarize vertOrder edgeOrder).size)
    (hR : (g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.type i = .R)
    (h : (g.planarTree ternarize vertOrder edgeOrder).isPlanar i = true) :
    IsPlanarEmbedding ((g.planarTree ternarize vertOrder edgeOrder).localSkeleton i)
      ((g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.nVerts i)
      ((g.planarTree ternarize vertOrder edgeOrder).nodeRot i) := by
  simp only [Graph.planarTree, Graph.planarSpqrTree] at hi hR h ⊢
  have hwf := planarWalk_items_wf' g hg ternarize vertOrder edgeOrder hvo heo
  exact nodePlanar_sound_R_of g _ _ hwf (planarWalk_planarFinish g ternarize _) i
    (planarRelabelTree_relabelNodeR g _ hwf (planarWalk_planarFinish g ternarize _) i hi hR h)

/-- Soundness of the per-node flag: a planar S/P/R node's local rotation system is a planar
embedding of its skeleton. S and P are `nodePlanar_sound_S` / `nodePlanar_sound_P` (proved modulo
`neRotAdj_segment`); R is `nodePlanar_sound_R` (admitted). -/
theorem nodePlanar_sound (g : Graph) (hg : g.WF) (ternarize : Bool) (vertOrder edgeOrder : List Nat)
    (hvo : OrderOK g.nv vertOrder) (heo : OrderOK g.ne edgeOrder) (i : Nat)
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
  · exact nodePlanar_sound_S g hg ternarize vertOrder edgeOrder hvo heo i hi hS
  · exact nodePlanar_sound_P g hg ternarize vertOrder edgeOrder hvo heo i hi hP
  · exact nodePlanar_sound_R g hg ternarize vertOrder edgeOrder hvo heo i hi hR h

/-- Corner structure of the node rotation systems (`NodeCorners`), under the all-planar flag.
S/P: from the explicit `layoutRot` cycle/bond layouts (`PlanarRelabelCornersSP.lean`); R: from the
relabel rows `NodeRow` and `PlanarFinish.corners` (`PlanarRelabelCornersR.lean`). -/
theorem planarTree_nodeCorners (g : Graph) (hg : g.WF) (ternarize : Bool)
    (vertOrder edgeOrder : List Nat) (hvo : OrderOK g.nv vertOrder) (heo : OrderOK g.ne edgeOrder)
    (hall : (g.planarTree ternarize vertOrder edgeOrder).nodePlanar.all id = true) :
    (g.planarTree ternarize vertOrder edgeOrder).NodeCorners := by
  have hwf : (g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.WF := by
    rw [planarRelabel_proj]; exact spqrTree_wf' g hg ternarize vertOrder edgeOrder hvo heo
  exact planarRelabelTree_nodeCorners g _
    (planarWalk_items_wf' g hg ternarize vertOrder edgeOrder hvo heo)
    (planarWalk_planarFinish g ternarize _) hwf
    (fun i hi => neRotAdj_segment g hg ternarize vertOrder edgeOrder hvo heo i hi) hall

/-- The rotation entries of every `R` node stay inside the node's own segment (`NodeRotClosed`),
under the all-planar flag: from the relabel rows `NodeRow` and `PlanarFinish.closed`
(`PlanarRelabelClosed.lean`). -/
theorem planarTree_nodeRotClosed (g : Graph) (hg : g.WF) (ternarize : Bool)
    (vertOrder edgeOrder : List Nat) (hvo : OrderOK g.nv vertOrder) (heo : OrderOK g.ne edgeOrder)
    (hall : (g.planarTree ternarize vertOrder edgeOrder).nodePlanar.all id = true) :
    (g.planarTree ternarize vertOrder edgeOrder).NodeRotClosed :=
  planarRelabelTree_nodeRotClosed g _
    (planarWalk_items_wf' g hg ternarize vertOrder edgeOrder hvo heo)
    (planarWalk_planarFinish g ternarize _) hall

/-- Completeness of the per-node flag: a node flagged nonplanar has a nonplanar skeleton.
Admitted (Kuratowski-style; plan in PROOF.md §8.3: the `mergePlanarity` nesting obstruction
exhibits a `K₅` / `K₃,₃` subdivision in the skeleton). -/
theorem nodePlanar_complete (g : Graph) (ternarize : Bool) (vertOrder edgeOrder : List Nat) (i : Nat)
    (hi : i < (g.planarTree ternarize vertOrder edgeOrder).size)
    (h : (g.planarTree ternarize vertOrder edgeOrder).isPlanar i = false) :
    ¬ Planar ((g.planarTree ternarize vertOrder edgeOrder).localSkeleton i)
      ((g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.nVerts i) := by
  sorry

/-- `planarRelabel` pushes one `nodePlanar` flag per item (next to the `types` push), so the flag
array has one entry per item (`RotInv.np_size`). -/
theorem planarTree_nodePlanar_size (g : Graph) (ternarize : Bool) (vertOrder edgeOrder : List Nat) :
    (g.planarTree ternarize vertOrder edgeOrder).nodePlanar.size =
      (g.planarTree ternarize vertOrder edgeOrder).size :=
  (PlanarRot.finalState_inv g (g.planarWalk ternarize (g.dfsForestFast vertOrder edgeOrder))).np_size

/-- The root piece contains every edge, and its local quarter-edge numbering permutes the
original quarter-edges. Its closed planar embedding transports to the glued rotation. -/
theorem glued_root (g : Graph) (hg : g.WF) (ternarize : Bool) (vertOrder edgeOrder : List Nat)
    (hvo : OrderOK g.nv vertOrder) (heo : OrderOK g.ne edgeOrder)
    (hwf : (g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.WF) (s : PlanarSpqrTree.EmbedState)
    (h : (g.planarTree ternarize vertOrder edgeOrder).GluedUpTo g 0 s) :
    IsPlanarEmbedding g.edges.toList g.nv ⟨s.rotAdj⟩ := by
  apply (g.planarTree ternarize vertOrder edgeOrder).glued_root_of g hwf ?_ ?_ s h
  · rw [planarRelabel_proj]
    exact spqrTree_childShape g hg ternarize vertOrder edgeOrder hvo heo
  · rw [planarRelabel_proj]
    rfl

/-- Soundness of the glued embedding: it is a planar embedding of `g` (Euler's formula per
component). The reverse-preorder fold is `forM_reverse_range_inv` with the invariant `GluedUpTo`
strengthened to `GluedFaces` (same-witness cofacial caps, PROOF.md §8.6); leaf, F, V and Q steps
are proved, the S/P/R step remains the admission `embedItem_step_node_faces`. The tree's `WF` comes
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
      rw [planarRelabel_proj]; exact spqrTree_childShape g hg ternarize vertOrder edgeOrder hvo heo
    have hrep : (g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.Represents g := by
      rw [planarRelabel_proj]; exact spqrTree_represents g hg ternarize vertOrder edgeOrder hvo heo
    have hsep : (g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.PieceSep g := by
      rw [planarRelabel_proj]; exact spqrTree_pieceSep g hg ternarize vertOrder edgeOrder hvo heo
    have hloc : ∀ i, i < (g.planarTree ternarize vertOrder edgeOrder).size →
        (g.planarTree ternarize vertOrder edgeOrder).types[i]! = .S ∨
          (g.planarTree ternarize vertOrder edgeOrder).types[i]! = .P ∨
          (g.planarTree ternarize vertOrder edgeOrder).types[i]! = .R →
        IsPlanarEmbedding ((g.planarTree ternarize vertOrder edgeOrder).localSkeleton i)
          ((g.planarTree ternarize vertOrder edgeOrder).toSpqrTree.nVerts i)
          ((g.planarTree ternarize vertOrder edgeOrder).nodeRot i) := by
      intro i hi hty
      refine nodePlanar_sound g hg ternarize vertOrder edgeOrder hvo heo i ?_
        (PlanarSpqrTree.isPlanar_of_all hall ?_)
      · rw [(g.planarTree ternarize vertOrder edgeOrder).type_eq_of_lt i hi]; exact hty
      · rw [planarTree_nodePlanar_size]; exact hi
    have hlay : ∀ i, i < (g.planarTree ternarize vertOrder edgeOrder).size →
        (g.planarTree ternarize vertOrder edgeOrder).LayoutAt i :=
      fun i hi => neRotAdj_segment g hg ternarize vertOrder edgeOrder hvo heo i hi
    exact glued_root g hg ternarize vertOrder edgeOrder hvo heo hwf _
      ((g.planarTree ternarize vertOrder edgeOrder).gluedFaces_planarEmbed g hg hwf hsh hrep hsep
        hloc hlay (planarTree_nodeRotClosed g hg ternarize vertOrder edgeOrder hvo heo hall)
        (planarTree_nodeCorners g hg ternarize vertOrder edgeOrder hvo heo hall)).toGluedUpTo
  · cases h

/-- The eventual target: the planar SPQR tree yields an embedding iff `g` is planar. `→` is
`planarEmbed_sound`; `←` needs `nodePlanar_complete` plus the fact that a minor (each skeleton is a
minor of `g`) of a planar graph is planar. Admitted. -/
theorem spqrTree_planar (g : Graph) (hg : g.WF) (ternarize : Bool) (vertOrder edgeOrder : List Nat)
    (hvo : OrderOK g.nv vertOrder) (heo : OrderOK g.ne edgeOrder) :
    (g.planarTree ternarize vertOrder edgeOrder).planarEmbed.isSome ↔ Planar g.edges.toList g.nv := by
  sorry

end Spqr
