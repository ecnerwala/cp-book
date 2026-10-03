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

end Spqr.WalkState
