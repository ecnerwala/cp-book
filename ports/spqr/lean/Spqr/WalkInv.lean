import Spqr.WalkSpec
import Spqr.Frame

/-!
# From the tstack guards to `FinishOk`

`finishOk_of_guards` assembles `WalkState.FinishOk` (the per-block hypotheses of `finishEdge_inv`)
from `FinishGuards`, the invariant, bookkeeping facts about the edge being finished, and the ear
facts `ear_*` below, which `FinishGuards` does not provide and are left to the ear invariant
(`EarShape`).
-/

namespace Spqr
open WalkM

namespace WalkState

variable {D : Nat} {s : WalkState}

theorem edgeBelow_vert_nil {v : Nat} (hv : v < s.g.nv) (hch : Items.ch s.items (vertItem v) = []) (e : Nat) :
    ¬ Items.EdgeBelow s.g s.items (vertItem v) e := by
  intro h
  rcases Relation.ReflTransGen.cases_head h with h | ⟨c, hc, -⟩
  · exact absurd h (by show (1 + v : Nat) ≠ 1 + s.g.nv + e; omega)
  · rw [Items.IsParent, hch] at hc; exact List.not_mem_nil hc

theorem finishTailOk_of_nil {curV d : Nat} {hasVert isSingle : Bool} (hv : curV < s.g.nv)
    (hch : hasVert = false → Items.ch s.items (vertItem curV) = [])
    (hm : hasVert = false → isSingle = false → MergeTopOk D (after (pushVertTstack curV d) s)) :
    FinishTailOk D curV d hasVert isSingle s :=
  ⟨fun h => Graph.ConnEdges.empty fun e _ => edgeBelow_vert_nil hv (hch h) e,
   fun h => Graph.TwoAttached.empty fun e _ => edgeBelow_vert_nil hv (hch h) e, hm⟩

section Ear

variable {curV d lv : Nat} {kind : RetKind} {o : DfsOut} {origTstack : Nat} {hasVert : Bool}

/-! Ear facts (Invariant W / `EarShape`), one per `FinishOk` field that `FinishGuards` does not
cover. Each is stated at the state where the block runs. -/

/-- Loop 1: every iteration merges/unwraps/closes a finished sub-ear (`Loop1BodyOk`). -/
theorem ear_loop1 (ho : o.cls = .ret lv kind) (hlow : lv < d) (ht : o.cls.isTree = true)
    (hg : FinishGuards d o origTstack hasVert s) (hi : s.Inv D) (hs : Shape s) (k : Nat)
    (hk : ∀ j, j ≤ k → result (loop1Cond d)
      (iter (Spqr.loop1Body d s.stackDir[d]!) j (ceS₁ o.dest d o.e (feS₀ d o s))) = true) :
    Loop1BodyOk D d s.stackDir[d]! (iter (Spqr.loop1Body d s.stackDir[d]!) k (ceS₁ o.dest d o.e (feS₀ d o s))) := by
  sorry

/-- Loop 2: every late merge joins entries sharing a terminal (`MergeTopOk`). -/
theorem ear_mergeLate (ho : o.cls = .ret lv kind) (hlow : lv < d) (ht : o.cls.isTree = true)
    (hg : FinishGuards d o origTstack hasVert s) (hi : s.Inv D) (hs : Shape s) :
    MergeLateOk D d (feS₁ d o s) := by
  sorry

/-- The vertex close: loop 3 merges, the unwrap, the two merges, the retarget and the type-1 close. -/
theorem ear_closeVert (ho : o.cls = .ret lv kind) (hlow : lv < d) (ht : o.cls.isTree = true)
    (hg : FinishGuards d o origTstack hasVert s) (hi : s.Inv D) (hs : Shape s) (hv : hasVert = true) :
    CloseVertOk D curV s.stackDir[d]! o.cls.isType1 origTstack (feSingle d o s) (feS₂ d o s) := by
  sorry

/-- The P-check after the vertex close. -/
theorem ear_finishP_vert (ho : o.cls = .ret lv kind) (hlow : lv < d) (ht : o.cls.isTree = true)
    (hg : FinishGuards d o origTstack hasVert s) (hi : s.Inv D) (hs : Shape s) (hv : hasVert = true) :
    FinishPOk D curV lv o.cls.isType1 (feS₃ curV d o origTstack s) := by
  sorry

/-- The P-check of a first tree edge (no vertex entry yet). -/
theorem ear_finishP_tree (ho : o.cls = .ret lv kind) (hlow : lv < d) (ht : o.cls.isTree = true)
    (hg : FinishGuards d o origTstack hasVert s) (hi : s.Inv D) (hs : Shape s) (hv : hasVert = false) :
    FinishPOk D curV lv o.cls.isType1 (feS₂ d o s) := by
  sorry

/-- The P-check of a back edge. -/
theorem ear_finishP_back (ho : o.cls = .ret lv kind) (hlow : lv < d) (hb : o.cls.isTree = false)
    (hg : FinishGuards d o origTstack hasVert s) (hi : s.Inv D) (hs : Shape s) :
    FinishPOk D curV lv o.cls.isType1 (feBack curV lv d o s) := by
  sorry

/-- The merge of the vertex entry into the type-2 first-edge entry. -/
theorem ear_tail_tree (ho : o.cls = .ret lv kind) (hlow : lv < d) (ht : o.cls.isTree = true)
    (hg : FinishGuards d o origTstack hasVert s) (hi : s.Inv D) (hs : Shape s) (hv : hasVert = false)
    (hsingle : feSingle d o s = false) :
    MergeTopOk D (after (pushVertTstack curV d) (after (finishP curV lv o.cls.isType1) (feS₂ d o s))) := by
  sorry

/-! Bookkeeping frame facts through the blocks: `g` is constant and the vertex item of `curV` gains
no children before its entry is pushed (only node items' `ch` change). -/

theorem tail_frame_tree (ho : o.cls = .ret lv kind) (hlow : lv < d) (ht : o.cls.isTree = true)
    (hg : FinishGuards d o origTstack hasVert s) (hi : s.Inv D) (hs : Shape s) (hv : hasVert = false) :
    (after (finishP curV lv o.cls.isType1) (feS₂ d o s)).g = s.g ∧
    (Items.ch s.items (vertItem curV) = [] →
      Items.ch (after (finishP curV lv o.cls.isType1) (feS₂ d o s)).items (vertItem curV) = []) := by
  sorry

theorem tail_frame_back (ho : o.cls = .ret lv kind) (hlow : lv < d) (hb : o.cls.isTree = false)
    (hg : FinishGuards d o origTstack hasVert s) (hi : s.Inv D) (hs : Shape s) :
    (after (finishP curV lv o.cls.isType1) (feBack curV lv d o s)).g = s.g ∧
    (Items.ch s.items (vertItem curV) = [] →
      Items.ch (after (finishP curV lv o.cls.isType1) (feBack curV lv d o s)).items (vertItem curV) = []) := by
  sorry

/-- `FinishOk` from the guards, the invariant, the bookkeeping facts of the finished edge
(`he`, `hq`, `hends`, `hch`) and the ear facts. -/
theorem finishOk_of_guards (ho : o.cls = .ret lv kind) (hlow : lv < d)
    (hg : FinishGuards d o origTstack hasVert s) (hi : s.Inv D) (hs : Shape s)
    (hdD : d ≤ D) (hv : curV < s.g.nv) (he : o.e < s.g.ne)
    (hq : Items.ch s.items (edgeItem s.g o.e) = [])
    (hends : Items.PairEq (if o.cls.isTree then (o.dest, s.stackVerts[d]!) else (curV, s.stackVerts[lv]!))
      s.g.edges[o.e]!)
    (hch : hasVert = false → Items.ch s.items (vertItem curV) = []) :
    FinishOk D curV d lv o origTstack hasVert s where
  e_lt := he
  ears ht := ⟨he,
    by show Items.ch (s.items.modify _ _) _ = []; rw [Items.ch_modify_ch_eq (edgeItem s.g o.e)
        (fun it => { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) }) (fun _ => rfl)]; exact hq,
    by show Items.PairEq (o.dest, s.stackVerts[d]!) s.g.edges[o.e]!; simpa [ht] using hends,
    hdD, ear_loop1 ho hlow ht hg hi hs⟩
  late ht := ear_mergeLate ho hlow ht hg hi hs
  vert ht hv' := ear_closeVert ho hlow ht hg hi hs hv'
  rest_vert ht hv' := ⟨ear_finishP_vert ho hlow ht hg hi hs hv',
    fun h => by simp [hv'] at h, fun h => by simp [hv'] at h, fun h => by simp [hv'] at h⟩
  rest_tree ht hv' :=
    have hf := tail_frame_tree (curV := curV) ho hlow ht hg hi hs hv'
    ⟨ear_finishP_tree ho hlow ht hg hi hs hv',
     finishTailOk_of_nil (by rw [hf.1]; exact hv) (fun _ => hf.2 (hch hv')) fun _ => ear_tail_tree ho hlow ht hg hi hs hv'⟩
  q _ := hq
  ends hb := by simpa [hb] using hends
  lv_le _ := by omega
  rest_back hb :=
    have hf := tail_frame_back (curV := curV) ho hlow hb hg hi hs
    ⟨ear_finishP_back ho hlow hb hg hi hs,
     finishTailOk_of_nil (by rw [hf.1]; exact hv) (fun h => hf.2 (hch h)) fun _ h => by cases h⟩

end Ear

end WalkState
end Spqr
