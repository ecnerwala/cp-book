import Spqr.RangesSchedule
import Spqr.RangesClose
import Spqr.EarWalk

namespace Spqr.WalkState
open WalkM

theorem Place.pushVertR {g : Graph} {P X : ItemId → Prop} {s : WalkState} {σ : List Nat} {n v : Nat}
    (h : s.Place g P X) (hv : v < g.nv) (hP : ∀ e, e < g.ne → P (edgeItem g e) → σ.idxOf e < n) :
    PushVertR σ n v s := by
  refine ⟨h.vert v hv, ?_⟩
  intro e he hb
  rw [h.g_eq] at he hb
  by_cases hc : 0 < s.cnt (edgeItem g e)
  · exact hP e he (h.fixed _ (by simp [edgeItem]) (by dsimp [edgeItem]; omega) hc)
  · have hn := noParent_of_cnt_eq_zero (Nat.eq_zero_of_not_pos hc)
    exact (vertItem_ne_edgeItem hv e (Items.Below.eq_of_no_parent hn hb)).elim

def PostAt (σ : List Nat) (n : Nat) (xs : List Nat) : Prop :=
  ∃ pre post, pre.length = n ∧ σ = pre ++ xs ++ post

namespace PostAt
variable {σ xs ys : List Nat} {n : Nat}

theorem block_infix (h : PostAt σ n xs) : xs <:+: σ := by
  obtain ⟨pre, post, _, rfl⟩ := h
  exact ⟨pre, post, rfl⟩

theorem left (h : PostAt σ n (xs ++ ys)) : PostAt σ n xs := by
  obtain ⟨pre, post, hn, rfl⟩ := h
  exact ⟨pre, ys ++ post, hn, by simp [List.append_assoc]⟩

theorem right (h : PostAt σ n (xs ++ ys)) : PostAt σ (n + xs.length) ys := by
  obtain ⟨pre, post, hn, rfl⟩ := h
  exact ⟨pre ++ xs, post, by simp [hn], by simp [List.append_assoc]⟩

theorem singleton {e : Nat} (h : PostAt σ n [e]) : σ[n]? = some e := by
  obtain ⟨pre, post, rfl, rfl⟩ := h
  simp [List.append_assoc]

theorem idx_bounds {e : Nat} (h : PostAt σ n xs) (hnd : σ.Nodup) (he : e ∈ xs) :
    n ≤ σ.idxOf e ∧ σ.idxOf e < n + xs.length := by
  obtain ⟨pre, post, rfl, rfl⟩ := h
  have hdisj := (List.nodup_append.mp (List.nodup_append.mp hnd).1).2.2
  have hpre : e ∉ pre := fun hp => hdisj e hp e he rfl
  simp only [List.idxOf_append, hpre, ite_false, List.mem_append, he, false_or, ite_true]
  have := List.idxOf_lt_length_iff.mpr he
  omega

end PostAt

theorem pushed_past {g : Graph} {σ : List Nat} {n : Nat} {P : ItemId → Prop}
    {vs es block : List Nat} (hP : ∀ e, e < g.ne → P (edgeItem g e) → σ.idxOf e < n)
    (hv : ∀ v ∈ vs, v < g.nv) (he : ∀ e ∈ es, e ∈ block)
    (hp : PostAt σ n block) (hnd : σ.Nodup) :
    ∀ e, e < g.ne → Pushed g P vs es (edgeItem g e) → σ.idxOf e < n + block.length := by
  rintro e helt (h | ⟨v, hv', hve⟩ | ⟨e', he', hee⟩)
  · exact Nat.lt_of_lt_of_le (hP e helt h) (Nat.le_add_right _ _)
  · exact (edgeItem_ne_vertItem (hv v hv') e hve).elim
  · have heq := edgeItem_inj hee
    exact (hp.idx_bounds hnd (heq ▸ he e' he')).2

theorem walkTree_past {g : Graph} {σ : List Nat} {n : Nat} {P X : ItemId → Prop}
    {s : WalkState} (t : DfsTree) (d : Nat) (h : s.Place g P X)
    (hv : ∀ v ∈ t.verts, v < g.nv) (he : ∀ e ∈ t.edges, e < g.ne)
    (hvn : t.verts.Nodup) (hen : t.edges.Nodup)
    (hPv : ∀ v ∈ t.verts, ¬ P (vertItem v)) (hPe : ∀ e ∈ t.edges, ¬ P (edgeItem g e))
    (hP : ∀ e, e < g.ne → P (edgeItem g e) → σ.idxOf e < n)
    (hp : PostAt σ n t.edgePostorder) (hnd : σ.Nodup) :
    wp (walkTree t d) (fun _ s' => ∀ v, v < g.nv → PushVertR σ (n + t.edgePostorder.length) v s') s := by
  apply wp_mono _ ((walk_place_aux g).1 t d P X s h hv he hvn hen hPv hPe)
  intro _ s' h' v hvg
  exact h'.pushVertR hvg (pushed_past hP hv
    (fun e he' => t.edgePostorder_perm_edges.symm.subset he') hp hnd)

theorem walkOutPre_place {g : Graph} {P X : ItemId → Prop} {s : WalkState}
    {v d : Nat} {o : DfsOut} {hasVert : Bool} (h : s.Place g P X) (hv : v < g.nv)
    (hf : hasVert = false → ¬ P (vertItem v)) :
    wp (walkOutPre v d o hasVert) (fun hv' s' =>
      s'.Place g (fun i => P i ∨ (hv' = true ∧ i = vertItem v)) X ∧
      (hasVert = true → hv' = true)) s := by
  unfold walkOutPre
  dsimp only
  rw [bind_stackDir]
  refine bind_spec (setStackDir_spec _ _) ?_
  rintro _ _ rfl
  split
  · next hc =>
    have hh : hasVert = false := by revert hc; cases hasVert <;> simp
    refine bind_spec (pushVertTstack_spec v d) ?_
    rintro _ _ rfl
    rw [wp_pure]
    refine ⟨((h.set_stackDir _).cons_fixed (by show 0 < 1 + v; omega)
      (by show 1 + v < _; omega) (hf hh) v d _ _).mono ?_ (fun _ hi => hi), ?_⟩
    · intro i hi; simpa using hi
    · intro _; rfl
  · rw [wp_pure]
    exact ⟨(h.set_stackDir _).mono (fun _ hi => Or.inl hi) (fun _ hi => hi), fun hi => hi⟩

theorem walkOutPre_past {g : Graph} {P X : ItemId → Prop} {s : WalkState} {σ : List Nat} {n : Nat}
    {v d : Nat} {o : DfsOut} {hasVert : Bool} (h : s.Place g P X) (hv : v < g.nv)
    (hf : hasVert = false → ¬ P (vertItem v))
    (hP : ∀ e, e < g.ne → P (edgeItem g e) → σ.idxOf e < n) :
    wp (walkOutPre v d o hasVert) (fun _ s' => ∀ w, w < g.nv → PushVertR σ n w s') s := by
  apply wp_mono _ (walkOutPre_place h hv hf)
  rintro hv' s' ⟨h', _⟩ w hw
  apply h'.pushVertR hw
  rintro e he (hp | ⟨_, hve⟩)
  · exact hP e he hp
  · exact (edgeItem_ne_vertItem hv e hve).elim

theorem finishEdge_past {g : Graph} {P X : ItemId → Prop} {s : WalkState} {σ : List Nat} {n : Nat}
    (v d : Nat) (o : DfsOut) (orig : Nat) (hasVert : Bool) (h : s.Place g P X)
    (hv : v < g.nv) (he : o.e < g.ne) (hPe : ¬ P (edgeItem g o.e))
    (hPv : hasVert = false → ¬ P (vertItem v))
    (hP : ∀ e, e < g.ne → P (edgeItem g e) → σ.idxOf e < n)
    (hp : σ[n]? = some o.e) (hnd : σ.Nodup) :
    wp (finishEdge v d o orig hasVert) (fun _ s' => ∀ w, w < g.nv → PushVertR σ (n + 1) w s') s := by
  apply wp_mono _ (finishEdge_place v d o orig hasVert h hv he hPe hPv)
  rintro hv' s' ⟨h', _⟩ w hw
  apply h'.pushVertR hw
  rintro e he' (h | h | ⟨_, h⟩)
  · exact Nat.lt_succ_of_lt (hP e he' h)
  · have heq := edgeItem_inj h
    subst e
    have hn := (List.getElem?_eq_some_iff.mp hp).1
    have heq : σ[n]! = o.e := by simp [hp]
    rw [← heq, idxOf_getElem! hnd hn]
    omega
  · exact (edgeItem_ne_vertItem hv e h).elim

theorem RootState.pushVertR {g : Graph} {pre rest : List DfsTree} {s : WalkState}
    (h : RootState g pre s) (hf : ForestOK g (pre ++ rest)) {v : Nat} (hv : v < g.nv) :
    PushVertR (edgePostorderForest (pre ++ rest)) (edgePostorderForest pre).length v s := by
  apply h.place.pushVertR hv
  have hp : PostAt (edgePostorderForest (pre ++ rest)) 0 (edgePostorderForest pre) :=
    ⟨[], edgePostorderForest rest, rfl, by simp [edgePostorderForest]⟩
  have hn := DfsData.edgePostorderForest_perm.nodup_iff.mpr hf.edges_nodup
  have hh := pushed_past (P := fun _ => False) (vs := pre.flatMap DfsTree.verts)
    (es := pre.flatMap DfsTree.edges) (fun _ _ h => h.elim)
    (fun v hv => hf.verts_lt v (by simpa only [List.flatMap_append] using List.mem_append_left (rest.flatMap DfsTree.verts) hv))
    (fun e he => DfsData.edgePostorderForest_perm.symm.subset he) hp hn
  simpa using hh

theorem walkOutPre_closeInv {s : WalkState} (h : s.CloseInv) (v d : Nat) (o : DfsOut)
    (hasVert : Bool) (hv : v < s.g.nv) (hs : s.g.nv < s.items.size) :
    wp (walkOutPre v d o hasVert) (fun _ s' => s'.CloseInv) s := by
  unfold walkOutPre
  dsimp only
  rw [bind_stackDir]
  refine bind_spec (setStackDir_spec _ _) ?_
  rintro _ _ rfl
  have h' := h.frame (s' := { s with stackDir := s.stackDir.set! d (if o.cls.lowval d ≥ d then false else !s.stackDir[o.cls.lowval d]!) }) rfl rfl (fun _ hi => hi)
  split
  · refine bind_spec (pushVertTstack_spec v d) ?_
    rintro _ _ rfl
    rw [wp_pure]
    exact h'.pushVert v d hv hs
  · rw [wp_pure]; exact h'

theorem rootAppend_closeInv {s : WalkState} {P X : ItemId → Prop} (h : s.CloseInv)
    (hp : s.Place s.g P X) (hk : RootOK s) :
    wp (popTstack >>= fun top => modifyItem rootItem fun it => { it with ch := it.ch ++ top.spans.2 })
      (fun _ s' => s'.CloseInv) s := by
  obtain ⟨t, ht, _, hvs⟩ := hk
  simp only [wp_bind, wp_popTstack, wp_modifyItem]
  apply h.pop.root_append
  · exact hp.of_le rfl (Nat.le_refl _) (fun _ _ => rfl) (fun i => by
      have := spansCount_tail_le s.tstack i
      dsimp [cnt]; omega)
  · simpa only [ht, head!_cons] using hvs

mutual
def CoverTree (σ : List Nat) (n : Nat) (t : DfsTree) (d : Nat) (s : WalkState) : Prop :=
  match t with
  | .node v outs => CoverOuts σ n v d outs false { s with stackVerts := s.stackVerts.set! d v }

def CoverOuts (σ : List Nat) (n v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  match outs with
  | [] => (hasVert = false → PushVertR σ n v s) ∧ (hasVert = false → VertFree v s)
  | o :: rest => CoverOut σ n v d o hasVert s ∧
      wp (walkOut v d o hasVert) (fun hasVert' s' => CoverOuts σ (n + o.block.length) v d rest hasVert' s') s

def CoverOut (σ : List Nat) (n v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  (hasVert = false → PushVertR σ n v s) ∧
  (hasVert = false → VertFree v s) ∧
  wp (walkOutPre v d o hasVert) (fun hasVert' s₁ =>
    match o with
    | .tree _ _ child =>
      wp (modify fun s => { s with firstOccurrence := s.firstOccurrence.set! d s.g.ne }) (fun _ s₂ =>
        CoverTree σ n child (d + 1) s₂ ∧
        wp (walkTree child (d + 1)) (fun _ s₃ =>
          FinishCover σ (n + child.edgePostorder.length) v d o s₁.tstack.length hasVert' s₃) s₂) s₁
    | .back .. => FinishCover σ n v d o s₁.tstack.length hasVert' s₁) s
end

theorem walkOut_past {g : Graph} {P X : ItemId → Prop} {s : WalkState} {σ : List Nat} {n : Nat}
    (v d : Nat) (o : DfsOut) (hasVert : Bool) (h : s.Place g P X) (hv : v < g.nv)
    (hvs : ∀ w ∈ DfsOut.vertsList [o], w < g.nv) (hes : ∀ e ∈ DfsOut.edgesList [o], e < g.ne)
    (hvn : (DfsOut.vertsList [o]).Nodup) (hen : (DfsOut.edgesList [o]).Nodup)
    (hvc : v ∉ DfsOut.vertsList [o])
    (hPv : ∀ w ∈ DfsOut.vertsList [o], ¬ P (vertItem w))
    (hPe : ∀ e ∈ DfsOut.edgesList [o], ¬ P (edgeItem g e))
    (hf : hasVert = false → ¬ P (vertItem v))
    (hP : ∀ e, e < g.ne → P (edgeItem g e) → σ.idxOf e < n)
    (hat : PostAt σ n o.block) (hnd : σ.Nodup) :
    wp (walkOut v d o hasVert) (fun _ s' => ∀ w, w < g.nv → PushVertR σ (n + o.block.length) w s') s := by
  apply wp_mono _ ((walk_place_aux g).2.2 v d o hasVert P X s h hv hvs hes hvn hen hvc hPv hPe hf)
  rintro hv' s' ⟨h', _⟩ w hw
  apply h'.pushVertR hw
  apply pushed_past (vs := DfsOut.vertsList [o]) (es := DfsOut.edgesList [o]) _ hvs _ hat hnd
  · rintro e he (hh | ⟨_, heq⟩)
    · exact hP e he hh
    · exact (edgeItem_ne_vertItem hv e heq).elim
  · intro e he
    have hm := (DfsOut.edgePostorderList_perm_edgesList [o]).symm.subset he
    simpa only [DfsOut.edgePostorderList_cons, DfsOut.edgePostorderList, List.append_nil] using hm

theorem coverOut_back {g : Graph} {P X : ItemId → Prop} {s : WalkState} {σ : List Nat} {n : Nat}
    (v d e dest : Nat) (cls : OutClass) (hasVert : Bool) (h : s.Place g P X)
    (hv : v < g.nv) (hf : hasVert = false → ¬ P (vertItem v))
    (hP : ∀ e, e < g.ne → P (edgeItem g e) → σ.idxOf e < n)
    (hc : wp (walkOutPre v d (.back e dest cls) hasVert) (fun hv' s' =>
      FinishPOwnership σ n v d (.back e dest cls) s'.tstack.length hv' s') s) :
    CoverOut σ n v d (.back e dest cls) hasVert s := by
  refine ⟨fun _ => h.pushVertR hv hP, fun h0 => h.vertFree hv (hf h0), ?_⟩
  refine wp_imp (wp_imp (wp_of_forall fun hv' s' hp hc' => ?_)
    (walkOutPre_past h hv hf hP)) hc
  exact { hc' with vert := fun _ => hp v hv }

theorem coverOut_tree {g : Graph} {P X : ItemId → Prop} {s : WalkState} {σ : List Nat} {n : Nat}
    (v d e : Nat) (cls : OutClass) (child : DfsTree) (hasVert : Bool) (h : s.Place g P X)
    (hv : v < g.nv) (hf : hasVert = false → ¬ P (vertItem v))
    (hvs : ∀ w ∈ child.verts, w < g.nv) (hes : ∀ e ∈ child.edges, e < g.ne)
    (hvn : child.verts.Nodup) (hen : child.edges.Nodup) (hvc : v ∉ child.verts)
    (hPv : ∀ w ∈ child.verts, ¬ P (vertItem w)) (hPe : ∀ e ∈ child.edges, ¬ P (edgeItem g e))
    (hP : ∀ e, e < g.ne → P (edgeItem g e) → σ.idxOf e < n)
    (hat : PostAt σ n child.edgePostorder) (hnd : σ.Nodup)
    (hc : wp (walkOutPre v d (.tree e cls child) hasVert) (fun hv' s₁ =>
      wp (modify fun s => { s with firstOccurrence := s.firstOccurrence.set! d s.g.ne }) (fun _ s₂ =>
        CoverTree σ n child (d + 1) s₂ ∧ wp (walkTree child (d + 1)) (fun _ s₃ =>
          FinishPOwnership σ (n + child.edgePostorder.length) v d (.tree e cls child)
            s₁.tstack.length hv' s₃) s₂) s₁) s) :
    CoverOut σ n v d (.tree e cls child) hasVert s := by
  refine ⟨fun _ => h.pushVertR hv hP, fun h0 => h.vertFree hv (hf h0), ?_⟩
  refine wp_imp (wp_imp (wp_of_forall fun hv' s₁ hp hc' => ?_)
    (walkOutPre_place h hv hf)) hc
  simp only [wp_modify] at hc' ⊢
  refine ⟨hc'.1, ?_⟩
  have hp' : ({ s₁ with firstOccurrence := s₁.firstOccurrence.set! d s₁.g.ne } : WalkState).Place
      g (fun i => P i ∨ (hv' = true ∧ i = vertItem v)) X :=
    hp.1.of_le rfl (Nat.le_refl _) (fun _ _ => rfl) (fun _ => Nat.le_refl _)
  have hpv : ∀ w ∈ child.verts, ¬ (P (vertItem w) ∨ (hv' = true ∧ vertItem w = vertItem v)) := by
    rintro w hw (h | ⟨_, heq⟩)
    · exact hPv w hw h
    · have hwv : w = v := by have hh : 1 + w = 1 + v := heq; omega
      exact hvc (hwv ▸ hw)
  have hpe : ∀ e ∈ child.edges, ¬ (P (edgeItem g e) ∨ (hv' = true ∧ edgeItem g e = vertItem v)) := by
    rintro e he (h | ⟨_, heq⟩)
    · exact hPe e he h
    · exact edgeItem_ne_vertItem hv e heq
  have hpast : ∀ e, e < g.ne → (P (edgeItem g e) ∨ (hv' = true ∧ edgeItem g e = vertItem v)) →
      σ.idxOf e < n := by
    rintro e he (h | ⟨_, heq⟩)
    · exact hP e he h
    · exact (edgeItem_ne_vertItem hv e heq).elim
  refine wp_imp (wp_imp (wp_of_forall fun _ s₃ hr hc₃ => ?_)
    (walkTree_past child (d + 1) hp' hvs hes hvn hen hpv hpe hpast hat hnd)) hc'.2
  exact { hc₃ with vert := fun _ => hr v hv }

abbrev ScheduleTree (σ : List Nat) (n : Nat) (t : DfsTree) (d : Nat) (s : WalkState) : Prop :=
  (∀ v outs, t = .node v outs → ({ s with stackVerts := s.stackVerts.set! d v } : WalkState).RangesInv σ n d) →
  Shape s → σ.Nodup → (∀ e ∈ σ, e < s.g.ne) → GuardsTree t d s → BookTree t d s →
    FrontiersTree t d s → CoverTree σ n t d s → PostAt σ n t.edgePostorder → RgTree σ n t d s

abbrev ScheduleOuts (σ : List Nat) (n v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  RgS σ n d s → σ.Nodup → GuardsOuts v d outs hasVert s → BookOuts v d outs hasVert s →
    FrontiersOuts v d outs hasVert s → CoverOuts σ n v d outs hasVert s →
    PostAt σ n (DfsOut.edgePostorderList outs) → RgOuts σ n v d outs hasVert s

abbrev ScheduleOut (σ : List Nat) (n v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  RgS σ n d s → σ.Nodup → GuardsOut v d o hasVert s → BookOut v d o hasVert s →
    FrontiersOut v d o hasVert s → CoverOut σ n v d o hasVert s →
    PostAt σ n o.block → RgOut σ n v d o hasVert s

mutual
theorem scheduleTree : ∀ σ n t d s, ScheduleTree σ n t d s
  | σ, n, .node v outs, d, s => fun hi hs hnd hσ hg hb hf hc hp => by
    unfold RgTree
    exact scheduleOuts σ n v d outs false _ ⟨hi v outs rfl, hs.frame', hσ⟩ hnd hg hb hf hc hp

theorem scheduleOuts : ∀ σ n v d outs hasVert s, ScheduleOuts σ n v d outs hasVert s
  | _, _, _, _, [], _, _ => fun _ _ _ _ _ hc _ => hc.1
  | σ, n, v, d, o :: rest, hasVert, s => fun hrs hnd hg hb hf hc hp => by
    unfold GuardsOuts at hg; unfold BookOuts at hb; unfold FrontiersOuts at hf; unfold CoverOuts at hc
    rw [DfsOut.edgePostorderList_cons] at hp
    have hr := scheduleOut σ n v d o hasVert s hrs hnd hg.1 hb.1 hf.1 hc.1 hp.left
    refine ⟨hr, ?_⟩
    exact wp_imp (wp_imp (wp_imp (wp_imp (wp_imp (wp_of_forall fun hv' s' hrs' hg' hb' hf' hc' =>
      scheduleOuts σ (n + o.block.length) v d rest hv' s' hrs' hnd hg' hb' hf' hc' hp.right)
      (rgOut σ n v d o hasVert s hrs hnd hg.1 hb.1 hr)) hg.2) hb.2) hf.2) hc.2

theorem scheduleOut : ∀ σ n v d o hasVert s, ScheduleOut σ n v d o hasVert s
  | σ, n, v, d, o, hasVert, s => fun ⟨hi, hs, hσ⟩ hnd hg hb hf hc hp => by
    unfold RgOut
    unfold GuardsOut at hg; unfold BookOut at hb; unfold FrontiersOut at hf; unfold CoverOut at hc
    refine ⟨hc.1, ?_⟩
    refine wp_imp (wp_imp (wp_imp (wp_imp (wp_imp (wp_of_forall fun hv' s₁ ⟨hi₁, hs₁, hσ₁⟩ hg₁ hb₁ hf₁ hc₁ => ?_)
      (walkOutPre_ranges hi hs hσ hb.1 hc.1)) hg) hb.2) hf) hc.2.2
    cases o with
    | back e cls dest =>
      exact finishR_of_cover (D := d)
        (by rw [if_neg fun h => by obtain ⟨_, _, _, h⟩ := hb₁.tree.1 h; cases h])
        hi₁ hs₁ hnd hσ₁ hg₁ hb₁ hf₁ hp.singleton hp.block_infix hc₁
    | tree e cls child =>
      try simp only [wp_modify] at hg₁ hb₁ hf₁ hc₁
      simp only [wp_modify]
      have pre : ∀ w outs, child = .node w outs →
          ({ { s₁ with firstOccurrence := s₁.firstOccurrence.set! d s₁.g.ne } with
            stackVerts := s₁.stackVerts.set! (d + 1) w } : WalkState).RangesInv σ n (d + 1) :=
        fun w outs _ => (hi₁.frame' (s' := { s₁ with firstOccurrence := s₁.firstOccurrence.set! d s₁.g.ne })).setSv w
      have hr := scheduleTree σ n child (d + 1) _ pre hs₁.frame' hnd hσ₁ hg₁.1 hb₁.1 hf₁.1 hc₁.1 hp.left
      refine ⟨hr, ?_⟩
      exact wp_imp (wp_imp (wp_of_forall fun _ s₃ ⟨hg₃, hb₃, hf₃, hc₃⟩ ⟨hi₃, hs₃, hσ₃⟩ =>
        finishR_of_cover (D := d + 1) (by rw [if_pos (hb₃.tree.2 ⟨_, _, _, rfl⟩)])
          hi₃ hs₃ hnd hσ₃ hg₃ hb₃ hf₃ hp.right.singleton hp.block_infix hc₃)
        (wp_and hg₁.2 (wp_and hb₁.2 (wp_and hf₁.2 hc₁.2))))
        (rgTree σ n child (d + 1) _ pre hs₁.frame' hnd hσ₁ hg₁.1 hb₁.1 hr)
end

theorem init_rangesInv (g : Graph) (tern : Bool) {σ : List Nat} (hnd : σ.Nodup) :
    (WalkState.init g tern).RangesInv σ 0 0 := by
  refine ⟨init_inv g tern, ?_, ?_, ?_, ?_⟩
  · intro t ht; simp [WalkState.init] at ht
  · intro above t below ht; simp [WalkState.init] at ht
  · intro t ht; simp [WalkState.init] at ht
  · intro i _ _ a b c hab hbc hcl ha hc
    have hae := Items.Below_of_ch_nil ha.edgeBelow (Items.initialItems_ch g i)
    have hce := Items.Below_of_ch_nil hc.edgeBelow (Items.initialItems_ch g i)
    have heq := edgeItem_inj (hae.trans hce.symm)
    have hia := idxOf_getElem! hnd (by omega : a < σ.length)
    have hic := idxOf_getElem! hnd hcl
    rw [heq] at hia
    have hba : b = a := by omega
    subst b
    exact ha.edgeBelow

theorem RangesInv.stackVerts_of_nil {σ : List Nat} {n : Nat} {s : WalkState}
    (h : s.RangesInv σ n 0) (ht : s.tstack = []) (sv : Array Nat) :
    ({ s with stackVerts := sv } : WalkState).RangesInv σ n 0 :=
  ⟨h.inv.stackVerts_of_nil ht sv, h.processed, h.ordered, h.convex, h.closed⟩

theorem RangesInv.root_append {σ : List Nat} {n : Nat} {s : WalkState}
    (h : s.RangesInv σ n 0) (hs : Shape s) (hk : RootOK s)
    (hp : ∀ p, ¬ Items.IsParent s.items p rootItem) :
    wp (popTstack >>= fun top => modifyItem rootItem fun it => { it with ch := it.ch ++ top.spans.2 })
      (fun _ s' => s'.RangesInv σ n 0) s := by
  obtain ⟨t, ht, -⟩ := hk
  have hi : ({ s with tstack := s.tstack.tail } : WalkState).Inv' 0 :=
    ⟨fun _ _ _ hx => by simp [ht] at hx, fun i h1 h2 => ⟨(h.inv.nodes i h1 h2).conn, (h.inv.nodes i h1 h2).attached⟩⟩
  simp only [wp_bind, wp_popTstack, wp_modifyItem]
  exact (h.pop' hi).modifyCh rootItem _ (by change 0 < 1 + s.g.nv + s.g.ne; omega) (fun _ => rfl) hp
    (by simp [ht]) (fun hn => by simp [hs.root] at hn)

def RootsCover (σ : List Nat) (n : Nat) : List DfsTree → WalkState → Prop
  | [], _ => True
  | t :: rest, s => CoverTree σ n t 0 s ∧
      wp (walkTree t 0) (fun _ s₁ =>
        wp (popTstack >>= fun top => modifyItem rootItem fun it => { it with ch := it.ch ++ top.spans.2 })
          (fun _ s₂ => RootsCover σ (n + t.edgePostorder.length) rest s₂) s₁) s

theorem forest_ranges_of_cover {g : Graph} {σ : List Nat} (hnd : σ.Nodup) (hσ : ∀ e ∈ σ, e < g.ne) :
    ∀ forest pre n s, RootState g pre s → s.RangesInv σ n 0 → ForestOK g (pre ++ forest) →
      (∀ t ∈ forest, t.WF []) → (∀ t ∈ forest, t.Ends g) →
      (∀ t ∈ forest, ∀ e, e < g.ne → ∀ x, g.Inc e x → x ∈ t.verts → e ∈ t.edges) →
      RootsCover σ n forest s →
      PostAt σ n (edgePostorderForest forest) →
      wp (walkForest forest) (fun _ s' => s'.RangesInv σ (n + (edgePostorderForest forest).length) 0) s
  | [], _, _, _, _, hr, _, _, _, _, _, _ => by simpa [walkForest, edgePostorderForest, wp_pure] using hr
  | t :: rest, pre, n, s, h, hr, hf, hwf, hends, hcomp, hc, hat => by
    have hb := h.book hf (hwf t (by simp)) (hends t (by simp)) (hcomp t (by simp))
    have hg := gbTree t 0 s hb
    have hi : ∀ v outs, t = .node v outs →
        ({ s with stackVerts := s.stackVerts.set! 0 v } : WalkState).RangesInv σ n 0 :=
      fun v _ _ => hr.stackVerts_of_nil h.tstack _
    have hσ' : ∀ e ∈ σ, e < s.g.ne := by rwa [h.g_eq]
    have hfront := walkTree_frontiers t 0 s (fun v outs ht => (hi v outs ht).inv) h.shape hg hb
    have hat' : PostAt σ n (t.edgePostorder ++ edgePostorderForest rest) := hat
    have hsched := scheduleTree σ n t 0 s hi h.shape hnd hσ' hg hb hfront hc.1 hat'.left
    have hrg := rgTree σ n t 0 s hi h.shape hnd hσ' hg hb hsched
    have hnv : 0 < s.g.nv := by
      obtain ⟨v, outs⟩ := t
      have hv := RootState.hvlt hf v (by simp [DfsTree.verts])
      rw [h.g_eq]; omega
    have hk : wp (walkTree t 0) (fun _ s' => RootOK s') s :=
      walkTree_rootOK t s (fun v outs ht => (hi v outs ht).inv) h.shape hg hb h.tstack
      (by rw [h.sd, ← h.g_eq]; exact hnv)
    have hp := (walk_place_aux g).1 t 0 _ _ s h.place (RootState.hvlt hf) (RootState.helt hf)
      (RootState.hvn hf).1 (RootState.hen hf).1 (RootState.hPv hf) (RootState.hPe hf)
    have hst := h.step hf (hwf t (by simp)) (hends t (by simp)) (hcomp t (by simp))
    show wp ((walkTree t 0 >>= fun _ => popTstack >>= fun top =>
      modifyItem rootItem fun it => { it with ch := it.ch ++ top.spans.2 }) >>= fun _ => walkForest rest) _ s
    rw [show edgePostorderForest (t :: rest) = t.edgePostorder ++ edgePostorderForest rest from rfl,
      List.length_append, ← Nat.add_assoc]
    simp only [wp_bind]
    refine wp_imp (wp_imp (wp_imp (wp_imp (wp_imp (wp_of_forall fun _ s₁ hrs hk hp hst hc => ?_)
      hrg) hk) hp) hst) hc.2
    have hrpop := hrs.1.root_append hrs.2.1 hk (noParent_of_cnt_eq_zero hp.root)
    exact wp_mono (popTstack >>= fun top => modifyItem rootItem fun it => { it with ch := it.ch ++ top.spans.2 })
      (wp_and hrpop (wp_and hst hc)) fun _ s₂ ⟨hr₂, hs₂, hc₂⟩ =>
      forest_ranges_of_cover hnd hσ rest (pre ++ [t]) (n + t.edgePostorder.length) s₂ hs₂ hr₂
        (by simpa using hf) (fun t' ht' => hwf t' (by simp [ht']))
        (fun t' ht' => hends t' (by simp [ht'])) (fun t' ht' => hcomp t' (by simp [ht'])) hc₂ hat'.right

theorem walk_rangesInv_of_cover (g : Graph) (tern : Bool) (forest : List DfsTree)
    (hf : ForestOK g forest) (hwf : ∀ t ∈ forest, t.WF []) (hends : ∀ t ∈ forest, t.Ends g)
    (hecov : ∀ e, e < g.ne → e ∈ forest.flatMap DfsTree.edges)
    (hc : RootsCover (edgePostorderForest forest) 0 forest (WalkState.init g tern)) :
    (g.walk tern forest).RangesInv (edgePostorderForest forest) (edgePostorderForest forest).length 0 := by
  have hnd := DfsData.edgePostorderForest_perm.nodup_iff.2 hf.edges_nodup
  have hσ : ∀ e ∈ edgePostorderForest forest, e < g.ne :=
    fun e he => hf.edges_lt e (DfsData.edgePostorderForest_perm.subset he)
  have h := forest_ranges_of_cover hnd hσ forest [] 0 (WalkState.init g tern)
    (rootState_init g tern) (init_rangesInv g tern hnd) (by simpa using hf) hwf hends
    (fun t ht => comp_of_forest hf hwf hends hecov ht) hc
    ⟨[], [], rfl, by simp⟩
  simpa only [wp, Graph.walk, Nat.zero_add] using h

end Spqr.WalkState
