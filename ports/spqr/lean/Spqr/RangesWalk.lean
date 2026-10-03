import Spqr.RangesSchedule
import Spqr.EarWalk

namespace Spqr.WalkState
open WalkM

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

end PostAt

mutual
def CoverTree (σ : List Nat) (n : Nat) (t : DfsTree) (d : Nat) (s : WalkState) : Prop :=
  match t with
  | .node v outs => CoverOuts σ n v d outs false { s with stackVerts := s.stackVerts.set! d v }

def CoverOuts (σ : List Nat) (n v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  match outs with
  | [] => hasVert = false → PushVertR σ n v s
  | o :: rest => CoverOut σ n v d o hasVert s ∧
      wp (walkOut v d o hasVert) (fun hasVert' s' => CoverOuts σ (n + o.block.length) v d rest hasVert' s') s

def CoverOut (σ : List Nat) (n v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  (hasVert = false → PushVertR σ n v s) ∧
  wp (walkOutPre v d o hasVert) (fun hasVert' s₁ =>
    match o with
    | .tree _ _ child =>
      wp (modify fun s => { s with firstOccurrence := s.firstOccurrence.set! d s.g.ne }) (fun _ s₂ =>
        CoverTree σ n child (d + 1) s₂ ∧
        wp (walkTree child (d + 1)) (fun _ s₃ =>
          FinishCover σ (n + child.edgePostorder.length) v d o s₁.tstack.length hasVert' s₃) s₂) s₁
    | .back .. => FinishCover σ n v d o s₁.tstack.length hasVert' s₁) s
end

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
  | _, _, _, _, [], _, _ => fun _ _ _ _ _ hc _ => hc
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
      (walkOutPre_ranges hi hs hσ hb.1 hc.1)) hg) hb.2) hf) hc.2
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
      (∀ t ∈ forest, t.WF []) → (∀ t ∈ forest, t.Ends g) → RootsCover σ n forest s →
      PostAt σ n (edgePostorderForest forest) →
      wp (walkForest forest) (fun _ s' => s'.RangesInv σ (n + (edgePostorderForest forest).length) 0) s
  | [], _, _, _, _, hr, _, _, _, _, _ => by simpa [walkForest, edgePostorderForest, wp_pure] using hr
  | t :: rest, pre, n, s, h, hr, hf, hwf, hends, hc, hat => by
    have hb := h.book hf (hwf t (by simp)) (hends t (by simp))
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
    have hst := h.step hf (hwf t (by simp)) (hends t (by simp))
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
        (fun t' ht' => hends t' (by simp [ht'])) hc₂ hat'.right

theorem walk_rangesInv_of_cover (g : Graph) (tern : Bool) (forest : List DfsTree)
    (hf : ForestOK g forest) (hwf : ∀ t ∈ forest, t.WF []) (hends : ∀ t ∈ forest, t.Ends g)
    (hc : RootsCover (edgePostorderForest forest) 0 forest (WalkState.init g tern)) :
    (g.walk tern forest).RangesInv (edgePostorderForest forest) (edgePostorderForest forest).length 0 := by
  have hnd := DfsData.edgePostorderForest_perm.nodup_iff.2 hf.edges_nodup
  have hσ : ∀ e ∈ edgePostorderForest forest, e < g.ne :=
    fun e he => hf.edges_lt e (DfsData.edgePostorderForest_perm.subset he)
  have h := forest_ranges_of_cover hnd hσ forest [] 0 (WalkState.init g tern)
    (rootState_init g tern) (init_rangesInv g tern hnd) (by simpa using hf) hwf hends hc
    ⟨[], [], rfl, by simp⟩
  simpa only [wp, Graph.walk, Nat.zero_add] using h

end Spqr.WalkState
