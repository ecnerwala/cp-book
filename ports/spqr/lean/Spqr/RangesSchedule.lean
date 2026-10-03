import Spqr.RangesFrontier

namespace Spqr.WalkState
open WalkM

variable {σ : List Nat} {n D : Nat} {s : WalkState}

/-- Range ownership still supplied by the walk induction, not by the ear frontier. -/
structure FinishCover (σ : List Nat) (n curV d : Nat) (o : DfsOut) (orig : Nat)
    (hasVert : Bool) (s : WalkState) : Prop where
  vert : hasVert = false → PushVertR σ n curV s
  p_vert : o.cls.isTree = true → hasVert = true →
    FinishPCover σ (n + 1) curV (o.cls.lowval d) o.cls.isType1 (feS₃ curV d o orig s)
  p_tree : o.cls.isTree = true → hasVert = false →
    FinishPCover σ (n + 1) curV (o.cls.lowval d) o.cls.isType1 (feS₂ d o s)
  p_back : o.cls.isTree = false →
    FinishPCover σ (n + 1) curV (o.cls.lowval d) o.cls.isType1 (feBack curV (o.cls.lowval d) d o s)

theorem Step.pushVertR {v m : Nat} {s' : WalkState} (st : Step D v s s')
    (hv : v < s.g.nv) (h : PushVertR σ n v s) (hn : n ≤ m) : PushVertR σ m v s' := by
  refine ⟨st.shape.vert v (by rw [st.g]; exact hv), ?_⟩
  intro e he hb
  rw [st.g] at he hb
  exact Nat.lt_of_lt_of_le (h.below e he ((st.below _).1 hb)) hn

theorem finishRestAdj_of_cover {curV d lv : Nat} {isType1 hasVert isSingle : Bool}
    (h : s.RangesInv σ n D) (hs : Shape s) (hv : curV < s.g.nv) (hnd : σ.Nodup)
    (hσ : ∀ e ∈ σ, e < s.g.ne) (hok : FinishRestOk D curV d lv isType1 hasVert isSingle s)
    (hc : FinishPCover σ n curV lv isType1 s) (hr : hasVert = false → PushVertR σ n curV s) :
    FinishRestAdj σ n curV d lv isType1 hasVert isSingle s := by
  have ha := finishPAdj_of_baseCover h hs hnd hσ hok.p hc
  have st := RgStep.finishP h hs hnd hσ hv hok.p ha
  exact ⟨ha, finishTailAdj_of_vert (fun hh => st.step.pushVertR hv (hr hh) (Nat.le_refl _))⟩

theorem finishAdj_of_cover {curV d lv orig : Nat} {kind : RetKind} {o : DfsOut} {hasVert : Bool}
    (he : o.cls = .ret lv kind) (hl : lv < d)
    (h : s.RangesInv σ n D) (hs : Shape s) (hnd : σ.Nodup) (hσ : ∀ e ∈ σ, e < s.g.ne)
    (hb : FinishBook curV d o orig hasVert s) (hf : Frontier (o := o) d orig s)
    (hok : FinishOk D curV d lv o orig hasVert s)
    (hpos : σ[n]? = some o.e) (hblock : o.block <:+: σ) (hc : FinishCover σ n curV d o orig hasVert s) :
    FinishAdj σ n curV d lv o orig hasVert s := by
  have hlv : o.cls.lowval d = lv := by rw [he]; rfl
  have hlow : o.cls.lowval d < d := by rwa [hlv]
  have hj : edgeItem s.g o.e < 1 + s.g.nv + s.g.ne := by
    show 1 + s.g.nv + o.e < _; have := hb.e_lt; omega
  have st₀ : RgStep σ n D curV s (feS₀ d o s) :=
    ⟨Step.modifyVs h.inv hs (edgeItem s.g o.e) _ hj, h.modifyVs (edgeItem s.g o.e) _ hj⟩
  have hv₀ : curV < (feS₀ d o s).g.nv := by rw [st₀.step.g]; exact hb.v_lt
  have hety : Items.type (feS₀ d o s).items (edgeItem (feS₀ d o s).g o.e) ≠ .V := by
    rw [st₀.step.shape.edge o.e (by rw [st₀.step.g]; exact hb.e_lt)]; decide
  have ha₁ := fun ht => closeEarsAdj_of_frontier hf ht hlow st₀.ranges st₀.step.shape hv₀ hnd
    (st₀.hσ hσ) hblock hpos hety (hok.ears ht)
  have st₁ : ∀ ht : o.cls.isTree = true, RgStep σ (n + 1) D curV s (feS₁ d o s) := fun ht =>
    st₀.trans (RgStep.closeEars st₀.ranges st₀.step.shape hnd (st₀.hσ hσ) hv₀ (hok.ears ht) (ha₁ ht))
  have ha₂ := fun ht => mergeLateAdj_of_frontier hf ht hlow (st₁ ht).ranges hnd
    ((st₁ ht).hσ hσ) hblock (hok.late ht)
  have st₂ : ∀ ht : o.cls.isTree = true, RgStep σ (n + 1) D curV s (feS₂ d o s) := fun ht =>
    (st₁ ht).trans (RgStep.mergeLate (st₁ ht).ranges (st₁ ht).step.shape hnd ((st₁ ht).hσ hσ)
      (by rw [(st₁ ht).step.g]; exact hb.v_lt) (hok.late ht) (ha₂ ht))
  have ha₃ := fun ht hv => closeVertAdj_of_frontier hf ht hlow (st₂ ht).ranges (st₂ ht).step.shape
    (by rw [(st₂ ht).step.g]; exact hb.v_lt) hnd ((st₂ ht).hσ hσ) hblock
    (hb.late_length ht hlow) (hok.vert ht hv)
  refine ⟨ha₁, ha₂, ha₃, ?_, ?_, fun _ => hpos, fun _ => hety, ?_⟩
  · intro ht hv
    have st : RgStep σ (n + 1) D curV s (feS₃ curV d o orig s) :=
      (st₂ ht).trans (RgStep.closeVert' (st₂ ht).ranges (st₂ ht).step.shape hnd ((st₂ ht).hσ hσ)
        (by rw [(st₂ ht).step.g]; exact hb.v_lt) (hok.vert ht hv) (ha₃ ht hv))
    exact finishRestAdj_of_cover st.ranges st.step.shape (by rw [st.step.g]; exact hb.v_lt) hnd
      (st.hσ hσ) (hok.rest_vert ht hv) (by simpa only [hlv] using hc.p_vert ht hv)
      (fun hh => by simp [hv] at hh)
  · intro ht hv
    exact finishRestAdj_of_cover (st₂ ht).ranges (st₂ ht).step.shape
      (by rw [(st₂ ht).step.g]; exact hb.v_lt) hnd ((st₂ ht).hσ hσ) (hok.rest_tree ht hv)
      (by simpa only [hlv] using hc.p_tree ht hv)
      (fun hh => (st₂ ht).step.pushVertR hb.v_lt (hc.vert hh) (by omega))
  · intro ht
    have hq : Items.ch (feS₀ d o s).items (edgeItem (feS₀ d o s).g o.e) = [] := by
      show Items.ch (s.items.modify _ _) (edgeItem s.g o.e) = []
      rw [Items.ch_modify_ch_eq (edgeItem s.g o.e)
        (fun it => { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) }) (fun _ => rfl)]
      exact hok.q ht
    have st₁ : RgStep σ (n + 1) D curV _ (after (pushEdgeTstack curV lv o.e) (feS₀ d o s)) :=
      RgStep.pushEdge st₀.ranges st₀.step.shape hnd curV lv o.e hb.e_lt hq hety (hok.ends ht)
        (hok.lv_le ht) hpos
    have st : RgStep σ (n + 1) D curV s (feBack curV lv d o s) := st₀.trans (st₁.trans
      ⟨Step.frame st₁.step.inv st₁.step.shape rfl rfl rfl rfl, st₁.ranges.frame rfl rfl rfl rfl⟩)
    exact finishRestAdj_of_cover st.ranges st.step.shape (by rw [st.step.g]; exact hb.v_lt) hnd
      (st.hσ hσ) (hok.rest_back ht) (by simpa only [hlv] using hc.p_back ht)
      (fun hh => st.step.pushVertR hb.v_lt (hc.vert hh) (by omega))

theorem finishR_of_cover {curV d orig : Nat} {o : DfsOut} {hasVert : Bool}
    (hD : D = if o.cls.isTree then d + 1 else d)
    (h : s.RangesInv σ n D) (hs : Shape s) (hnd : σ.Nodup) (hσ : ∀ e ∈ σ, e < s.g.ne)
    (hg : FinishGuards d o orig hasVert s) (hb : FinishBook curV d o orig hasVert s)
    (hf : Frontier (o := o) d orig s) (hpos : σ[n]? = some o.e) (hblock : o.block <:+: σ)
    (hc : FinishCover σ n curV d o orig hasVert s) : FinishR σ n curV d o orig hasVert s := by
  refine ⟨?_, fun hge => boundaryAdj_of_book hge hb hnd hσ hpos hblock⟩
  intro lv kind he hl
  obtain ⟨sub, base, hlen, hE⟩ := hb.ear
  have hok := finishOk_of_guards he hl hg hE hlen h.inv hs hD hb.v_lt hb.e_lt hb.q (hb.ends lv kind he) hb.vert
  exact finishAdj_of_cover he hl h hs hnd hσ hb hf hok hpos hblock hc

end Spqr.WalkState
