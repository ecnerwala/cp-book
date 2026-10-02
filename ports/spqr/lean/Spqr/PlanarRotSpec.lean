import Spqr.PlanarRelabelProj
import Spqr.PlanarLayout
import Spqr.Correctness

/-!
# `neRotAdj` segments of the planar relabel

`planarRelabel_rot_spec` (admitted; the `planarRelabel` fold, mirroring `relabel_node_spec`): the
block `planarRelabel` appends to `neRotAdj` for the node numbered `n` is
`layoutRot (type n) (nVerts n) neSt neEn edgeVes mapRot (2 g.ne)`, and `neRotAdj` grows in lockstep
with `nodeEdges` (four quarter-edges per node-edge). `layoutRot_size` (every type has `4 · nEdges`
entries, from `WF.shape`) and `neRotAdj_segment'` turn it into the `extract` form used by
`PlanarSpec.neRotAdj_segment`.
-/

namespace Spqr

namespace PlanarRot

theorem foldl_append_size {α : Type} (mr : Nat → Array α) (hmr : ∀ ve, (mr ve).size = 4)
    (ev : List Nat) (init : Array α) :
    (ev.foldl (fun acc ve => acc ++ mr ve) init).size = init.size + 4 * ev.length := by
  induction ev generalizing init with
  | nil => simp
  | cons v ev ih => simp only [List.foldl_cons, ih, Array.size_append, hmr, List.length_cons]; omega

theorem skeleton_length (t : SpqrTree) (i : Nat) : (t.skeleton i).length = t.nEdges i := by
  simp [SpqrTree.skeleton, SpqrTree.nodeEdgesOf, SpqrTree.nEdges]

theorem neBounds_le (t : SpqrTree) (hwf : t.WF) :
    ∀ k j, j + k ≤ t.size → t.neBounds[j]! ≤ t.neBounds[j + k]! := by
  intro k
  induction k with
  | zero => intros; simp
  | succ k ih =>
    intro j hj
    refine le_trans (ih j (by omega)) ?_
    rw [show j + (k + 1) = j + k + 1 by omega]
    exact hwf.own.ne_bounds_mono (j + k) (by rw [hwf.sizes.neBounds]; omega)

theorem neRange_eq (t : SpqrTree) (i : Nat) :
    t.neRange i = (t.neBounds[i]!, t.neBounds[i + 1]!) := by
  simp [SpqrTree.neRange, getElem!_def]
  constructor <;> (split <;> simp_all)

theorem neRange_mono (t : SpqrTree) (hwf : t.WF) (i : Nat) (hi : i < t.size) :
    (t.neRange i).1 ≤ (t.neRange i).2 := by
  rw [neRange_eq]
  exact hwf.own.ne_bounds_mono i (by rw [hwf.sizes.neBounds]; omega)

theorem neRange_le_size (t : SpqrTree) (hwf : t.WF) (i : Nat) (hi : i < t.size) :
    (t.neRange i).2 ≤ t.nodeEdges.size := by
  rw [neRange_eq, ← hwf.own.ne_last]
  have := neBounds_le t hwf (t.size - (i + 1)) (i + 1) (by omega)
  rwa [show i + 1 + (t.size - (i + 1)) = t.size by omega] at this

/-- `layoutRot` of a well-formed node has four entries per node-edge. -/
theorem layoutRot_size (t : SpqrTree) (hwf : t.WF) (n : Nat) (hn : n < t.size)
    (ev : List Nat) (mr : Nat → Array (Option Nat)) (cv : Nat) (hmr : ∀ ve, (mr ve).size = 4)
    (hR : t.type n = .R → ev.length + 1 = t.nEdges n) :
    (layoutRot (t.type n) (t.nVerts n) (t.neRange n).1 (t.neRange n).2 ev mr cv).size =
      4 * t.nEdges n := by
  have hsh := hwf.shape n hn
  have hlen := skeleton_length t n
  have hmono := neRange_mono t hwf n hn
  have hcap := hwf.twins.cap_none n hn
  obtain ⟨ty, ht⟩ : ∃ ty, t.type n = ty := ⟨_, rfl⟩
  obtain ⟨⟨s, e⟩, hnv⟩ : ∃ p, t.nvRange n = p := ⟨_, rfl⟩
  obtain ⟨⟨neSt, neEn⟩, hne⟩ : ∃ p, t.neRange n = p := ⟨_, rfl⟩
  unfold SpqrTree.Shape at hsh
  unfold SpqrTree.nVerts SpqrTree.nEdges at *
  rw [hnv] at hsh ⊢
  rw [hne] at hsh hlen hR hmono hcap ⊢
  rw [ht] at hsh hR hcap ⊢
  simp only at hsh hlen hR hmono hcap ⊢
  obtain ⟨k, rfl⟩ : ∃ k, neEn = neSt + k := ⟨neEn - neSt, by omega⟩
  simp only [Nat.add_sub_cancel_left] at hsh hlen hR hcap ⊢
  cases ty <;> simp only [NodeType.isNode, Bool.false_eq_true, not_false_eq_true, forall_const] at hsh hcap
  case F => simp [layoutRot, hcap]
  case V => simp [layoutRot, hcap]
  case Q =>
    rcases hsh with ⟨h1, h2⟩ | ⟨h1, h2⟩ <;> rw [h2] at hlen <;> simp at hlen <;>
      simp [layoutRot, h1, ← hlen]
  case I =>
    obtain ⟨h1, h2⟩ := hsh
    rw [h2] at hlen; simp at hlen
    simp [layoutRot, h1, ← hlen]
  case O =>
    obtain ⟨h1, h2⟩ := hsh
    rw [h2] at hlen; simp at hlen
    simp [layoutRot, h1, ← hlen]
  case P =>
    rw [hsh.1]
    exact layoutRot_P_size k neSt ev mr cv
  case S =>
    obtain ⟨h1, h2⟩ := hsh
    rw [h2] at hlen; simp at hlen
    have hk : k = e - s := by omega
    subst hk
    exact layoutRot_S_size (e - s) neSt (by omega) ev mr cv
  case R =>
    have h1 : (e - s == 1) = false := by simp; omega
    have hR' := hR rfl
    simp [layoutRot, h1, foldl_append_size mr hmr, hmr]
    omega

section

variable (g : Graph) (ternarize : Bool) (vertOrder edgeOrder : List Nat)

local notation "PT" => g.planarSpqrTree ternarize vertOrder edgeOrder

/-- Admitted (the `planarRelabel` fold; mirrors `relabel_node_spec`): `neRotAdj` has four entries
per node-edge, and node `n`'s entries `4 neSt + j`, `j < 4 nEdges n`, are those of `layoutRot` with
the `edgeVes`/`mapRot` computed when `n` was laid out (`mapRot` always yields four entries; for an R
node `edgeVes` lists its `nEdges n - 1` non-cap edges). -/
theorem planarRelabel_rot_spec :
    (PT).neRotAdj.size = 4 * (PT).toSpqrTree.nodeEdges.size ∧
    ∀ n, n < (PT).toSpqrTree.size →
      ∃ (edgeVes : List Nat) (mapRot : Nat → Array (Option Nat)),
        (∀ ve, (mapRot ve).size = 4) ∧
        ((PT).toSpqrTree.type n = .R → edgeVes.length + 1 = (PT).toSpqrTree.nEdges n) ∧
        ∀ j, j < 4 * (PT).toSpqrTree.nEdges n →
          (PT).neRotAdj[4 * ((PT).toSpqrTree.neRange n).1 + j]! =
            (layoutRot ((PT).toSpqrTree.type n) ((PT).toSpqrTree.nVerts n) ((PT).toSpqrTree.neRange n).1
              ((PT).toSpqrTree.neRange n).2 edgeVes mapRot (2 * g.ne))[j]! := by
  sorry

/-- `PlanarSpec.neRotAdj_segment` from `planarRelabel_rot_spec`. -/
theorem neRotAdj_segment' (i : Nat) (hi : i < (PT).toSpqrTree.size) :
    ∃ (edgeVes : List Nat) (mapRot : Nat → Array (Option Nat)) (capVe : Nat),
      (PT).neRotAdj.extract (4 * ((PT).toSpqrTree.neRange i).1) (4 * ((PT).toSpqrTree.neRange i).2) =
        layoutRot ((PT).toSpqrTree.type i) ((PT).toSpqrTree.nVerts i) ((PT).toSpqrTree.neRange i).1
          ((PT).toSpqrTree.neRange i).2 edgeVes mapRot capVe := by
  obtain ⟨hsz, hnode⟩ := planarRelabel_rot_spec g ternarize vertOrder edgeOrder
  obtain ⟨ev, mr, hmr, hR, hget⟩ := hnode i hi
  have hwf : (PT).toSpqrTree.WF := by
    rw [planarRelabel_proj]; exact spqrTree_wf g ternarize vertOrder edgeOrder
  refine ⟨ev, mr, 2 * g.ne, ?_⟩
  have hrs := layoutRot_size (PT).toSpqrTree hwf i hi ev mr (2 * g.ne) hmr hR
  have hle := neRange_le_size (PT).toSpqrTree hwf i hi
  have hmono := neRange_mono (PT).toSpqrTree hwf i hi
  have hsz' : ((PT).neRotAdj.extract (4 * ((PT).toSpqrTree.neRange i).1)
      (4 * ((PT).toSpqrTree.neRange i).2)).size = 4 * (PT).toSpqrTree.nEdges i := by
    rw [Array.size_extract, hsz]; unfold SpqrTree.nEdges; omega
  apply Array.ext
  · rw [hsz', hrs]
  · intro j hj1 hj2
    rw [hsz'] at hj1
    have h := hget j hj1
    rw [getElem!_pos _ _ (by rw [hsz]; unfold SpqrTree.nEdges at hj1; omega),
      getElem!_pos _ _ (by rw [hrs]; exact hj1)] at h
    rw [Array.getElem_extract]
    exact h

end

end PlanarRot

end Spqr
