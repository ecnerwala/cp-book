import Spqr.PlanarRelabelProj
import Spqr.PlanarLayout
import Spqr.Correctness
import Spqr.PlanarRotFold

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

theorem nvRange_eq (t : SpqrTree) (i : Nat) :
    t.nvRange i = (t.nvBounds[i]!, t.nvBounds[i + 1]!) := by
  simp [SpqrTree.nvRange, getElem!_def]
  constructor <;> (split <;> simp_all)

theorem nVerts_eq (t : SpqrTree) (i : Nat) : t.nVerts i = t.nvBounds[i + 1]! - t.nvBounds[i]! := by
  simp [SpqrTree.nVerts, nvRange_eq]

theorem nEdges_eq (t : SpqrTree) (i : Nat) : t.nEdges i = t.neBounds[i + 1]! - t.neBounds[i]! := by
  simp [SpqrTree.nEdges, neRange_eq]

theorem type_eq (t : SpqrTree) (i : Nat) (hi : i < t.size) : t.type i = t.types[i]! := by
  have hi' : i < t.types.size := hi
  unfold SpqrTree.type; rw [getElem?_pos t.types i hi', getElem!_pos t.types i hi']; rfl

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

/-- `planarRelabel_rot_spec` for any tree whose fields are read off a `RotInv` state. -/
theorem rot_spec_of_inv (g : Graph) (T : PlanarSpqrTree) (s : PlanarRelabelState)
    (hty : T.toSpqrTree.types = s.base.types) (hne : T.toSpqrTree.neBounds = s.base.neBounds)
    (hnv : T.toSpqrTree.nvBounds = s.base.nvBounds) (hE : T.toSpqrTree.nodeEdges = s.base.nodeEdges)
    (hrot : T.neRotAdj = s.aux.neRotAdj) (hwf : T.toSpqrTree.WF) (hinv : RotInv g s) :
    T.neRotAdj.size = 4 * T.toSpqrTree.nodeEdges.size ∧
    ∀ n, n < T.toSpqrTree.size →
      ∃ (edgeVes : List Nat) (mapRot : Nat → Array (Option Nat)),
        (∀ ve, (mapRot ve).size = 4) ∧
        (T.toSpqrTree.type n = .R → edgeVes.length + 1 = T.toSpqrTree.nEdges n) ∧
        ∀ j, j < 4 * T.toSpqrTree.nEdges n →
          T.neRotAdj[4 * (T.toSpqrTree.neRange n).1 + j]! =
            (layoutRot (T.toSpqrTree.type n) (T.toSpqrTree.nVerts n) (T.toSpqrTree.neRange n).1
              (T.toSpqrTree.neRange n).2 edgeVes mapRot (2 * g.ne))[j]! := by
  obtain ⟨_, h1, h2, h3, h4, h5, ⟨ev, mr, hmr, hR, hb⟩, _⟩ := hinv
  have hsz : T.toSpqrTree.size = s.base.types.size := congrArg Array.size hty
  have hRT : ∀ n, n < s.base.types.size → T.toSpqrTree.type n = .R →
      (ev n).length + 1 = T.toSpqrTree.nEdges n := by
    intro n hn hRn
    rw [nEdges_eq, hne]
    exact hR n hn (by rwa [type_eq _ _ (by omega), hty] at hRn)
  have hblk : ∀ n, n < s.base.types.size →
      rotBlock g.ne s.base.types s.base.nvBounds s.base.neBounds ev mr n =
        layoutRot (T.toSpqrTree.type n) (T.toSpqrTree.nVerts n) (T.toSpqrTree.neRange n).1
          (T.toSpqrTree.neRange n).2 (ev n) (mr n) (2 * g.ne) := by
    intro n hn
    rw [rotBlock, type_eq _ _ (by omega), nVerts_eq, neRange_eq, hty, hnv, hne]
  have hbsz : ∀ n, n < s.base.types.size →
      (rotBlock g.ne s.base.types s.base.nvBounds s.base.neBounds ev mr n).size =
        4 * (s.base.neBounds[n + 1]! - s.base.neBounds[n]!) := by
    intro n hn
    rw [hblk n hn, layoutRot_size T.toSpqrTree hwf n (by omega) (ev n) (mr n) _ (hmr n) (hRT n hn), nEdges_eq, hne]
  have hcs : ∀ n, n ≤ s.base.types.size →
      (concatBlocks (rotBlock g.ne s.base.types s.base.nvBounds s.base.neBounds ev mr) n).size =
        4 * s.base.neBounds[n]! := by
    intro n
    induction n with
    | zero => intro _; rw [concatBlocks.zero, h3]; rfl
    | succ n ih =>
      intro hn
      rw [concatBlocks.succ, Array.size_append, ih (by omega), hbsz n (by omega)]
      have := neRange_mono T.toSpqrTree hwf n (by omega)
      rw [neRange_eq, hne] at this
      simp only at this
      omega
  refine ⟨?_, fun n hn => ⟨ev n, mr n, hmr n, hRT n (by omega), fun j hj => ?_⟩⟩
  · rw [hrot, hb, hcs _ le_rfl, h4, hE]
  · rw [hsz] at hn
    have hidx : 4 * (T.toSpqrTree.neRange n).1 =
        (concatBlocks (rotBlock g.ne s.base.types s.base.nvBounds s.base.neBounds ev mr) n).size := by
      rw [hcs n (by omega), neRange_eq, hne]
    have hj' : j < (rotBlock g.ne s.base.types s.base.nvBounds s.base.neBounds ev mr n).size := by
      rwa [hbsz n hn, ← hne, ← nEdges_eq]
    rw [hrot, hb, hidx, concatBlocks.get hn hj', hblk n hn]

section

variable (g : Graph) (w : PlanarWalkState)

/-- The final `planarRelabel` state behind `planarRelabelTree`. -/
def finalState : PlanarRelabelState :=
  (planarRelabel w.base.items.size rootItem none none none).run (PlanarRelabelState.init g w) |>.2

theorem finalState_inv : RotInv g (finalState g w) :=
  planarRelabel_rotInv g _ _ _ _ _ _ (RotInv.init g w)

theorem planarRelabelTree_rot_spec (hwf : (planarRelabelTree g w).toSpqrTree.WF) :
    (planarRelabelTree g w).neRotAdj.size = 4 * (planarRelabelTree g w).toSpqrTree.nodeEdges.size ∧
    ∀ n, n < (planarRelabelTree g w).toSpqrTree.size →
      ∃ (edgeVes : List Nat) (mapRot : Nat → Array (Option Nat)),
        (∀ ve, (mapRot ve).size = 4) ∧
        ((planarRelabelTree g w).toSpqrTree.type n = .R →
          edgeVes.length + 1 = (planarRelabelTree g w).toSpqrTree.nEdges n) ∧
        ∀ j, j < 4 * (planarRelabelTree g w).toSpqrTree.nEdges n →
          (planarRelabelTree g w).neRotAdj[4 * ((planarRelabelTree g w).toSpqrTree.neRange n).1 + j]! =
            (layoutRot ((planarRelabelTree g w).toSpqrTree.type n) ((planarRelabelTree g w).toSpqrTree.nVerts n)
              ((planarRelabelTree g w).toSpqrTree.neRange n).1 ((planarRelabelTree g w).toSpqrTree.neRange n).2
              edgeVes mapRot (2 * g.ne))[j]! :=
  rot_spec_of_inv g (planarRelabelTree g w) (finalState g w) rfl rfl rfl rfl rfl hwf (finalState_inv g w)

end

section

variable (g : Graph) (ternarize : Bool) (vertOrder edgeOrder : List Nat)

local notation "PT" => g.planarSpqrTree ternarize vertOrder edgeOrder

/-- Admitted (the `planarRelabel` fold; mirrors `relabel_node_spec`): `neRotAdj` has four entries
per node-edge, and node `n`'s entries `4 neSt + j`, `j < 4 nEdges n`, are those of `layoutRot` with
the `edgeVes`/`mapRot` computed when `n` was laid out (`mapRot` always yields four entries; for an R
node `edgeVes` lists its `nEdges n - 1` non-cap edges). -/
theorem planarRelabel_rot_spec (hg : g.WF) (hvo : OrderOK g.nv vertOrder)
    (heo : OrderOK g.ne edgeOrder) :
    (PT).neRotAdj.size = 4 * (PT).toSpqrTree.nodeEdges.size ∧
    ∀ n, n < (PT).toSpqrTree.size →
      ∃ (edgeVes : List Nat) (mapRot : Nat → Array (Option Nat)),
        (∀ ve, (mapRot ve).size = 4) ∧
        ((PT).toSpqrTree.type n = .R → edgeVes.length + 1 = (PT).toSpqrTree.nEdges n) ∧
        ∀ j, j < 4 * (PT).toSpqrTree.nEdges n →
          (PT).neRotAdj[4 * ((PT).toSpqrTree.neRange n).1 + j]! =
            (layoutRot ((PT).toSpqrTree.type n) ((PT).toSpqrTree.nVerts n) ((PT).toSpqrTree.neRange n).1
              ((PT).toSpqrTree.neRange n).2 edgeVes mapRot (2 * g.ne))[j]! := by
  have hwf : (PT).toSpqrTree.WF := by
    rw [planarRelabel_proj]; exact spqrTree_wf g hg ternarize vertOrder edgeOrder hvo heo
  exact planarRelabelTree_rot_spec g _ hwf

/-- `PlanarSpec.neRotAdj_segment` from `planarRelabel_rot_spec`. -/
theorem neRotAdj_segment' (hg : g.WF) (hvo : OrderOK g.nv vertOrder)
    (heo : OrderOK g.ne edgeOrder) (i : Nat) (hi : i < (PT).toSpqrTree.size) :
    ∃ (edgeVes : List Nat) (mapRot : Nat → Array (Option Nat)) (capVe : Nat),
      (PT).neRotAdj.extract (4 * ((PT).toSpqrTree.neRange i).1) (4 * ((PT).toSpqrTree.neRange i).2) =
        layoutRot ((PT).toSpqrTree.type i) ((PT).toSpqrTree.nVerts i) ((PT).toSpqrTree.neRange i).1
          ((PT).toSpqrTree.neRange i).2 edgeVes mapRot capVe := by
  obtain ⟨hsz, hnode⟩ := planarRelabel_rot_spec g ternarize vertOrder edgeOrder hg hvo heo
  obtain ⟨ev, mr, hmr, hR, hget⟩ := hnode i hi
  have hwf : (PT).toSpqrTree.WF := by
    rw [planarRelabel_proj]; exact spqrTree_wf g hg ternarize vertOrder edgeOrder hvo heo
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
