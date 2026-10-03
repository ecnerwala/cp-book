import Spqr.RangesOwned
import Spqr.RangesCloseTree

/-! # Root coverage: prefix ownership through the walk

`OwnedD` (`RangesOwned.lean`) is carried through the walk of a DFS tree in the `RgTree` style; at every
`finishEdge` it is instantiated to `Owned` for `finishP_ownership`, which supplies the P-site part of
`CoverOut`, and the vertex part comes from `Place`. -/

namespace Spqr
open WalkM
namespace WalkState

variable {σ sts origs : List Nat} {P : ItemId → Prop} {d n : Nat} {s : WalkState}

theorem OwnedD.mono {P' : ItemId → Prop} (h : OwnedD σ sts origs P d n s) (hP : ∀ i, P i → P' i) :
    OwnedD σ sts origs P' d n s :=
  { h with
    vis := fun t ht => (h.vis t ht).imp_left (hP _)
    fresh := fun w hw hP' hk => h.fresh w hw (fun hp => hP' (hP _ hp)) hk }

theorem OwnedD.frame {s' : WalkState} (h : OwnedD σ sts origs P d n s) (hg : s'.g = s.g)
    (hit : s'.items = s.items) (hts : s'.tstack = s.tstack) (hsv : s'.stackVerts = s.stackVerts) :
    OwnedD σ sts origs P d n s' := by
  obtain ⟨l, c, a, hi, lo, nw, o, vi, fr⟩ := h
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;> simp only [hg, hit, hts, hsv]
  exacts [l, c, a, hi, lo, nw, o, vi, fr]

/-- Pushing the vertex entry of `stackVerts[d]`. -/
theorem OwnedD.pushVert {s' : WalkState} {v fi : Nat} {dir : Bool} (h : OwnedD σ sts origs P d n s)
    (hsv : s.stackVerts[d]! = v) (hsts : ∀ k, k ≤ d → sts[k]! ≤ sts[d]!)
    (hg : s'.g = s.g) (hit : s'.items = s.items) (hsv' : s'.stackVerts = s.stackVerts)
    (hts : s'.tstack = ⟨v, d, fi, setSides dir [vertItem v] []⟩ :: s.tstack) :
    OwnedD σ sts origs P d n s' := by
  obtain ⟨l, c, a, hi, lo, nw, o, vi, fr⟩ := h
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;> simp only [hg, hit, hts, hsv', List.length_cons]
  · intro k hk; have := l k hk; omega
  · intro b hb hb'
    rcases c b hb hb' with ⟨t, ht, he⟩ | hp
    · exact Or.inl ⟨t, List.mem_cons_of_mem _ ht, he⟩
    · exact Or.inr hp
  · exact a
  · exact hi
  · exact lo
  · intro k hk t ht e he hte
    have hl := l k hk
    rw [show s.tstack.length + 1 - origs[k]! = (s.tstack.length - origs[k]!) + 1 by omega,
      List.take_succ_cons, List.mem_cons] at ht
    rcases ht with rfl | ht
    · obtain ⟨i, hi', hb⟩ := hte
      rw [mem_setSides, List.mem_singleton] at hi'
      subst hi'
      have := lo d (Nat.le_refl _) e he (by rw [hsv]; exact hb)
      exact Nat.le_trans (hsts k hk) this
    · exact nw k hk t ht e he hte
  · intro k hk t ht
    have hl := l k hk
    rw [show s.tstack.length + 1 - origs[k]! = (s.tstack.length - origs[k]!) + 1 by omega,
      List.drop_succ_cons] at ht
    exact o k hk t ht
  · intro t ht
    rcases List.mem_cons.1 ht with rfl | ht
    · exact Or.inr ⟨d, Nat.le_refl _, hsv.symm⟩
    · exact vi t ht
  · exact fr

theorem vertCover_pushVert {s' : WalkState} {v fi : Nat} {dir : Bool}
    (hts : s'.tstack = ⟨v, d, fi, setSides dir [vertItem v] []⟩ :: s.tstack) : VertCover v s' := by
  intro e he hb
  rw [hts]
  exact ⟨_, List.mem_cons_self .., vertItem v, (mem_setSides dir _ _).2 (List.mem_singleton_self _), hb⟩

/-- Entering the child `w` at depth `d + 1`: its start is `n`, its entry stack length the current one. -/
theorem OwnedD.entry {s' : WalkState} (h : OwnedD σ sts origs P d n s) (hl : sts.length = d + 1)
    (hol : origs.length = d + 1) {w : Nat} (hw : w < s.g.nv) (hsz : d + 1 < s.stackVerts.size)
    (hPw : ¬ P (vertItem w)) (hwk : ∀ k, k ≤ d → w ≠ s.stackVerts[k]!)
    (hg : s'.g = s.g) (hit : s'.items = s.items) (hts : s'.tstack = s.tstack)
    (hsv' : s'.stackVerts = s.stackVerts.set! (d + 1) w) :
    OwnedD σ (sts ++ [n]) (origs ++ [s.tstack.length]) P (d + 1) n s' := by
  have hsv : ∀ k, k ≤ d → s'.stackVerts[k]! = s.stackVerts[k]! :=
    fun k hk => by rw [hsv']; exact getElem!_set!_ne _ _ _ _ (by omega)
  have hsvd : s'.stackVerts[d + 1]! = w := by rw [hsv']; exact getElem!_set!_self _ _ _ hsz
  have hst : ∀ k, k ≤ d → (sts ++ [n])[k]! = sts[k]! := fun k hk => getElem!_append_left' (by omega)
  have hstd : (sts ++ [n])[d + 1]! = n := by rw [← hl]; exact getElem!_concat_length' _ _
  have hor : ∀ k, k ≤ d → (origs ++ [s.tstack.length])[k]! = origs[k]! :=
    fun k hk => getElem!_append_left' (by omega)
  have hord : (origs ++ [s.tstack.length])[d + 1]! = s.tstack.length := by
    rw [← hol]; exact getElem!_concat_length' _ _
  have hnb : ∀ e, ¬ Items.EdgeBelow s.g s.items (vertItem w) e :=
    edgeBelow_vert_nil hw (h.fresh w hw hPw hwk)
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;> simp only [hg, hit, hts]
  · intro k hk
    rcases Nat.lt_or_ge k (d + 1) with hk' | hk'
    · rw [hor k (by omega)]; exact h.len k (by omega)
    · obtain rfl : k = d + 1 := by omega
      rw [hord]
  · intro b hb hb'
    rw [hst 0 (Nat.zero_le _)] at hb
    rcases h.cover b hb hb' with hs | ⟨k, hk, hb⟩
    · exact Or.inl hs
    · exact Or.inr ⟨k, by omega, by rw [hsv k hk]; exact hb⟩
  · intro k hk e he hb
    rcases Nat.lt_or_ge k d with hk' | hk'
    · rw [hst (k + 1) (by omega)]
      rw [hsv k (by omega)] at hb
      exact h.anc k hk' e he hb
    · obtain rfl : k = d := by omega
      rw [hstd]
      rw [hsv k (Nat.le_refl _)] at hb
      exact h.hi e he hb
  · intro e he hb
    rw [hsvd] at hb
    exact (hnb e hb).elim
  · intro k hk e he hb
    rcases Nat.lt_or_ge k (d + 1) with hk' | hk'
    · rw [hst k (by omega)]
      rw [hsv k (by omega)] at hb
      exact h.lo k (by omega) e he hb
    · obtain rfl : k = d + 1 := by omega
      rw [hsvd] at hb
      exact (hnb e hb).elim
  · intro k hk t ht e he hte
    rcases Nat.lt_or_ge k (d + 1) with hk' | hk'
    · rw [hor k (by omega)] at ht
      rw [hst k (by omega)]
      exact h.new k (by omega) t ht e he hte
    · obtain rfl : k = d + 1 := by omega
      rw [hord, Nat.sub_self, List.take_zero] at ht
      simp at ht
  · intro k hk t ht
    rcases Nat.lt_or_ge k (d + 1) with hk' | hk'
    · rw [hor k (by omega)] at ht
      rw [hsv k (by omega)]
      exact h.old k (by omega) t ht
    · obtain rfl : k = d + 1 := by omega
      rw [hord, Nat.sub_self, List.drop_zero] at ht
      rw [hsvd]
      intro heq
      rcases h.vis t ht with hp | ⟨k', hk', hk''⟩
      · rw [heq] at hp; exact hPw hp
      · rw [heq] at hk''; exact hwk k' hk' hk''
  · intro t ht
    rcases h.vis t ht with hp | ⟨k, hk, hk'⟩
    · exact Or.inl hp
    · exact Or.inr ⟨k, by omega, by rw [hsv k hk]; exact hk'⟩
  · intro x hx hPx hxk
    exact h.fresh x hx hPx fun k hk => by rw [← hsv k hk]; exact hxk k (by omega)

/-- Leaving the child at depth `d + 1` once it is visited (`P`) and its vertex item is on the stack. -/
theorem OwnedD.exit {m orig : Nat} (h : OwnedD σ (sts ++ [m]) (origs ++ [orig]) P (d + 1) n s)
    (hl : sts.length = d + 1) (hol : origs.length = d + 1) (hm : m ≤ n)
    (hσ : ∀ e ∈ σ, e < s.g.ne) (hnσ : n ≤ σ.length)
    (hc : VertCover s.stackVerts[d + 1]! s) (hPw : P (vertItem s.stackVerts[d + 1]!)) :
    OwnedD σ sts origs P d n s := by
  have hst : ∀ k, k ≤ d → (sts ++ [m])[k]! = sts[k]! := fun k hk => getElem!_append_left' (by omega)
  have hstd : (sts ++ [m])[d + 1]! = m := by rw [← hl]; exact getElem!_concat_length' _ _
  have hor : ∀ k, k ≤ d → (origs ++ [orig])[k]! = origs[k]! := fun k hk => getElem!_append_left' (by omega)
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · intro k hk; rw [← hor k hk]; exact h.len k (by omega)
  · intro b hb hb'
    rw [← hst 0 (Nat.zero_le _)] at hb
    rcases h.cover b hb hb' with hs | ⟨k, hk, hbk⟩
    · exact Or.inl hs
    · rcases Nat.lt_or_ge k (d + 1) with hk' | hk'
      · exact Or.inr ⟨k, by omega, hbk⟩
      · obtain rfl : k = d + 1 := by omega
        refine Or.inl (hc _ (hσ _ ?_) hbk)
        rw [getElem!_pos σ b (by omega)]; exact List.getElem_mem _
  · intro k hk e he hb; rw [← hst (k + 1) (by omega)]; exact h.anc k (by omega) e he hb
  · intro e he hb
    have := h.anc d (Nat.lt_succ_self _) e he hb
    rw [hstd] at this; omega
  · intro k hk e he hb; rw [← hst k hk]; exact h.lo k (by omega) e he hb
  · intro k hk t ht e he hte
    rw [← hst k hk]; rw [← hor k hk] at ht
    exact h.new k (by omega) t ht e he hte
  · intro k hk t ht; rw [← hor k hk] at ht; exact h.old k (by omega) t ht
  · intro t ht
    rcases h.vis t ht with hp | ⟨k, hk, hk'⟩
    · exact Or.inl hp
    · rcases Nat.lt_or_ge k (d + 1) with hk'' | hk''
      · exact Or.inr ⟨k, by omega, hk'⟩
      · obtain rfl : k = d + 1 := by omega
        rw [hk']; exact Or.inl hPw
  · intro w hw hPw' hwk
    refine h.fresh w hw hPw' fun k hk => ?_
    rcases Nat.lt_or_ge k (d + 1) with hk' | hk'
    · exact hwk k (by omega)
    · obtain rfl : k = d + 1 := by omega
      intro heq; rw [heq] at hPw'; exact hPw' hPw

theorem OwnedD.toOwned (h : OwnedD σ sts origs P d n s) {v : Nat} (hv : s.stackVerts[d]! = v) :
    Owned σ sts[0]! sts[d]! n v d origs[d]! P sts s where
  sv := hv
  len := h.len d (Nat.le_refl _)
  cover := h.cover
  anc := h.anc
  vertHi := by rw [← hv]; exact h.hi
  vertLo := by rw [← hv]; exact h.lo d (Nat.le_refl _)
  new := h.new d (Nat.le_refl _)
  old := by rw [← hv]; exact h.old d (Nat.le_refl _)
  vis := h.vis
  fresh := h.fresh

theorem walkOutPre_hv (v d : Nat) (o : DfsOut) (hasVert : Bool) :
    wp (walkOutPre v d o hasVert) (fun hv' _ =>
      (o.cls.lowval d < d → o.cls.isType1 = true → hv' = true) ∧ (hasVert = true → hv' = true)) s := by
  unfold walkOutPre
  simp only [wp_bind, wp_stackDir, wp_setStackDir, wp_ite, wp_pushVertTstack, wp_pure]
  split
  · exact ⟨fun _ _ => by trivial, fun _ => by trivial⟩
  · rename_i hc
    refine ⟨fun h1 h2 => ?_, fun _ => by trivial⟩
    cases hasVert
    · exact absurd (by simp [h1, h2]) hc
    · trivial


theorem walkOutPre_ownedD {v : Nat} (d : Nat) (o : DfsOut) (hasVert : Bool)
    (h : OwnedD σ sts origs P d n s) (hsv : s.stackVerts[d]! = v)
    (hsts : ∀ k, k ≤ d → sts[k]! ≤ sts[d]!) (hvc : hasVert = true → VertCover v s) :
    wp (walkOutPre v d o hasVert) (fun hv' s' =>
      OwnedD σ sts origs P d n s' ∧ (hv' = true → VertCover v s')) s := by
  unfold walkOutPre
  simp only [wp_bind, wp_stackDir, wp_setStackDir, wp_ite, wp_pushVertTstack, wp_pure]
  split
  · exact ⟨h.pushVert hsv hsts rfl rfl rfl rfl, fun _ => vertCover_pushVert rfl⟩
  · exact ⟨h.frame rfl rfl rfl rfl, hvc⟩

theorem Pushed.vert {g : Graph} {P : ItemId → Prop} {v : Nat} {vs es : List Nat} {b : Bool} (i : ItemId)
    (h : Pushed g (fun i => P i ∨ (b = true ∧ i = vertItem v)) vs es i) : Pushed g P (v :: vs) es i := by
  rcases h with (h | ⟨-, rfl⟩) | ⟨w, hw, rfl⟩ | ⟨e, he, rfl⟩
  · exact Or.inl h
  · exact Or.inr (Or.inl ⟨v, List.mem_cons_self .., rfl⟩)
  · exact Or.inr (Or.inl ⟨w, List.mem_cons_of_mem _ hw, rfl⟩)
  · exact Or.inr (Or.inr ⟨e, he, rfl⟩)

theorem Pushed.append {g : Graph} {P : ItemId → Prop} {v : Nat} {vs₁ es₁ vs₂ es₂ : List Nat} {b₁ b₂ : Bool}
    (hb : b₁ = true → b₂ = true) (i : ItemId)
    (h : Pushed g (fun i => Pushed g (fun i => P i ∨ (b₁ = true ∧ i = vertItem v)) vs₁ es₁ i ∨
      (b₂ = true ∧ i = vertItem v)) vs₂ es₂ i) :
    Pushed g (fun i => P i ∨ (b₂ = true ∧ i = vertItem v)) (vs₁ ++ vs₂) (es₁ ++ es₂) i := by
  rcases h with (((h | ⟨h1, rfl⟩) | ⟨w, hw, rfl⟩ | ⟨e, he, rfl⟩) | ⟨h2, rfl⟩) | ⟨w, hw, rfl⟩ | ⟨e, he, rfl⟩
  · exact Or.inl (Or.inl h)
  · exact Or.inl (Or.inr ⟨hb h1, rfl⟩)
  · exact Or.inr (Or.inl ⟨w, List.mem_append_left _ hw, rfl⟩)
  · exact Or.inr (Or.inr ⟨e, List.mem_append_left _ he, rfl⟩)
  · exact Or.inl (Or.inr ⟨h2, rfl⟩)
  · exact Or.inr (Or.inl ⟨w, List.mem_append_right _ hw, rfl⟩)
  · exact Or.inr (Or.inr ⟨e, List.mem_append_right _ he, rfl⟩)

theorem Pushed.finish {g : Graph} {P : ItemId → Prop} {v e : Nat} {vs es : List Nat} {b₁ b₂ : Bool}
    (hb : b₁ = true → b₂ = true) (i : ItemId)
    (h : Pushed g (fun i => P i ∨ (b₁ = true ∧ i = vertItem v)) vs es i ∨
      i = edgeItem g e ∨ (b₂ = true ∧ i = vertItem v)) :
    Pushed g (fun i => P i ∨ (b₂ = true ∧ i = vertItem v)) vs (e :: es) i := by
  rcases h with ((h | ⟨h1, rfl⟩) | ⟨w, hw, rfl⟩ | ⟨e', he', rfl⟩) | rfl | ⟨h2, rfl⟩
  · exact Or.inl (Or.inl h)
  · exact Or.inl (Or.inr ⟨hb h1, rfl⟩)
  · exact Or.inr (Or.inl ⟨w, hw, rfl⟩)
  · exact Or.inr (Or.inr ⟨e', List.mem_cons_of_mem _ he', rfl⟩)
  · exact Or.inr (Or.inr ⟨e, List.mem_cons_self .., rfl⟩)
  · exact Or.inl (Or.inr ⟨h2, rfl⟩)

theorem DfsOut.edges_subset_block (o : DfsOut) : ∀ e ∈ o.edges, e ∈ o.block := by
  cases o with
  | back e dest cls => intro e' he'; simpa [DfsOut.edges, DfsOut.block] using he'
  | tree e cls child =>
    intro e' he'
    simp only [DfsOut.edges, List.mem_cons] at he'
    simp only [DfsOut.block, List.mem_append, List.mem_singleton]
    rcases he' with rfl | he'
    · exact Or.inr rfl
    · exact Or.inl (child.edgePostorder_perm_edges.symm.subset he')

/-! ### The induction -/

/-- `CoverTree` for the walk of `t` at depth `d`, in the `RgTree` style, carrying `Place`, the
depth-indexed ownership `OwnedD` (`sts`/`origs` = schedule starts / entry stack lengths of the strict
ancestors, extended by `n`/the current stack length for `t`'s root), the ancestor path `anc`, and the
DFS facts of `t`; afterwards `Place`/`OwnedD` hold for the pushed items and the root's vertex item is
on the stack (`VertCover`). -/
abbrev CvPost (g : Graph) (σ : List Nat) (n v d : Nat) (sts origs : List Nat) (P X : ItemId → Prop)
    (hasVert : Bool) (vs es : List Nat) (hv' : Bool) (s' : WalkState) : Prop :=
  s'.Place g (Pushed g (fun i => P i ∨ (hv' = true ∧ i = vertItem v)) vs es) X ∧
  OwnedD σ sts origs (Pushed g (fun i => P i ∨ (hv' = true ∧ i = vertItem v)) vs es) d n s' ∧
  s'.stackVerts[d]! = v ∧ (hv' = true → VertCover v s') ∧ (hasVert = true → hv' = true)

abbrev CvTree (g : Graph) (σ : List Nat) (n : Nat) (t : DfsTree) (d : Nat) (anc sts origs : List Nat)
    (P X : ItemId → Prop) (s : WalkState) : Prop :=
  (∀ v outs, t = .node v outs → ({ s with stackVerts := s.stackVerts.set! d v } : WalkState).RangesInv σ n d) →
  Shape s → σ.Nodup → (∀ e ∈ σ, e < s.g.ne) → GuardsTree t d s → BookTree t d s → FrontiersTree t d s →
  CsTree σ n t d s →
  s.Place g P X → (∀ e, e < g.ne → P (edgeItem g e) → σ.idxOf e < n) →
  (∀ v outs, t = .node v outs →
    OwnedD σ (sts ++ [n]) (origs ++ [s.tstack.length]) P d n { s with stackVerts := s.stackVerts.set! d v }) →
  sts.length = d → origs.length = d → (∀ k, k < d → sts[k]! ≤ n) → (∀ k, k < d → origs[k]! ≤ s.tstack.length) →
  d = anc.length → (anc ++ t.verts).Nodup → (∀ a ∈ anc, a < g.nv) →
  (∀ v ∈ t.verts, v < g.nv) → (∀ e ∈ t.edges, e < g.ne) → t.edges.Nodup →
  (∀ v ∈ t.verts, ¬ P (vertItem v)) → (∀ e ∈ t.edges, ¬ P (edgeItem g e)) → PostAt σ n t.edgePostorder →
  s.stackVerts.size = g.nv → (∀ k, k < d → anc[k]? = some s.stackVerts[k]!) →
  CoverTree σ n t d s ∧
  wp (walkTree t d) (fun _ s' =>
    s'.Place g (Pushed g P t.verts t.edges) X ∧
    OwnedD σ (sts ++ [n]) (origs ++ [s.tstack.length]) (Pushed g P t.verts t.edges) d
      (n + t.edgePostorder.length) s' ∧
    ∀ v outs, t = .node v outs → s'.stackVerts[d]! = v ∧ VertCover v s') s

abbrev CvOuts (g : Graph) (σ : List Nat) (n v d : Nat) (outs : List DfsOut) (hasVert : Bool)
    (anc sts origs : List Nat) (P X : ItemId → Prop) (s : WalkState) : Prop :=
  RgS σ n d s → σ.Nodup → GuardsOuts v d outs hasVert s → BookOuts v d outs hasVert s →
  FrontiersOuts v d outs hasVert s → CsOuts σ n v d outs hasVert s → s.Place g P X → (∀ e, e < g.ne → P (edgeItem g e) → σ.idxOf e < n) →
  OwnedD σ sts origs P d n s → sts.length = d + 1 → origs.length = d + 1 →
  (∀ k, k ≤ d → sts[k]! ≤ sts[d]!) → sts[d]! ≤ n → (∀ k, k ≤ d → origs[k]! ≤ origs[d]!) →
  d = anc.length → (anc ++ v :: DfsOut.vertsList outs).Nodup → (∀ a ∈ anc, a < g.nv) → v < g.nv →
  (∀ w ∈ DfsOut.vertsList outs, w < g.nv) → (∀ e ∈ DfsOut.edgesList outs, e < g.ne) →
  (DfsOut.edgesList outs).Nodup →
  (∀ w ∈ DfsOut.vertsList outs, ¬ P (vertItem w)) → (∀ e ∈ DfsOut.edgesList outs, ¬ P (edgeItem g e)) →
  (hasVert = false → ¬ P (vertItem v)) → PostAt σ n (DfsOut.edgePostorderList outs) →
  s.stackVerts.size = g.nv → (∀ k, k < d → anc[k]? = some s.stackVerts[k]!) → s.stackVerts[d]! = v →
  (hasVert = true → VertCover v s) →
  CoverOuts σ n v d outs hasVert s ∧
  wp (walkOuts v d outs hasVert) (CvPost g σ (n + (DfsOut.edgePostorderList outs).length) v d sts origs P X hasVert
    (DfsOut.vertsList outs) (DfsOut.edgesList outs)) s

abbrev CvOut (g : Graph) (σ : List Nat) (n v d : Nat) (o : DfsOut) (hasVert : Bool)
    (anc sts origs : List Nat) (P X : ItemId → Prop) (s : WalkState) : Prop :=
  RgS σ n d s → σ.Nodup → GuardsOut v d o hasVert s → BookOut v d o hasVert s →
  FrontiersOut v d o hasVert s → CsOut σ n v d o hasVert s → s.Place g P X → (∀ e, e < g.ne → P (edgeItem g e) → σ.idxOf e < n) →
  OwnedD σ sts origs P d n s → sts.length = d + 1 → origs.length = d + 1 →
  (∀ k, k ≤ d → sts[k]! ≤ sts[d]!) → sts[d]! ≤ n → (∀ k, k ≤ d → origs[k]! ≤ origs[d]!) →
  d = anc.length → (anc ++ v :: o.verts).Nodup → (∀ a ∈ anc, a < g.nv) → v < g.nv →
  (∀ w ∈ o.verts, w < g.nv) → (∀ e ∈ o.edges, e < g.ne) → o.edges.Nodup →
  (∀ w ∈ o.verts, ¬ P (vertItem w)) → (∀ e ∈ o.edges, ¬ P (edgeItem g e)) →
  (hasVert = false → ¬ P (vertItem v)) → PostAt σ n o.block →
  s.stackVerts.size = g.nv → (∀ k, k < d → anc[k]? = some s.stackVerts[k]!) → s.stackVerts[d]! = v →
  (hasVert = true → VertCover v s) →
  CoverOut σ n v d o hasVert s ∧
  wp (walkOut v d o hasVert) (CvPost g σ (n + o.block.length) v d sts origs P X hasVert o.verts o.edges) s

mutual
theorem cvTree : ∀ (g : Graph) (σ : List Nat) (n : Nat) (t : DfsTree) (d : Nat) (anc sts origs : List Nat)
    (P X : ItemId → Prop) (s : WalkState), CvTree g σ n t d anc sts origs P X s
  | g, σ, n, .node v outs, d, anc, sts, origs, P, X, s =>
    fun hi hs hnd hσ hg hb hf hc hpl hP ho hls hlo hsts horigs hd hnodup hanc hw he hen hPv hPe hat
      hsz hsvk => by
    unfold GuardsTree at hg; unfold BookTree at hb; unfold FrontiersTree at hf
    unfold CsTree at hc
    simp only [DfsTree.verts, DfsTree.edges, DfsTree.edgePostorder] at hnodup hw he hen hPv hPe hat ⊢
    unfold CoverTree
    unfold walkTree
    simp only [wp_bind, wp_modify]
    have hv : v < g.nv := hw v (List.mem_cons_self ..)
    have hvo : v ∉ DfsOut.vertsList outs := (List.nodup_cons.1 (List.nodup_append.1 hnodup).2.1).1
    have hdlt : d < g.nv := by
      have hsub : anc ++ [v] ⊆ List.range g.nv := fun a ha => by
        rw [List.mem_range]
        rcases List.mem_append.1 ha with ha | ha
        · exact hanc a ha
        · exact (List.mem_singleton.1 ha) ▸ hv
      have hnd' : (anc ++ [v]).Nodup :=
        hnodup.sublist (List.Sublist.append_left (List.singleton_sublist.2 (List.mem_cons_self ..)) _)
      have := (List.Nodup.subperm hnd' hsub).length_le
      simp at this; omega
    have hsz₀ : ({ s with stackVerts := s.stackVerts.set! d v } : WalkState).stackVerts.size = g.nv := by
      simp [hsz]
    have hsvd₀ : ({ s with stackVerts := s.stackVerts.set! d v } : WalkState).stackVerts[d]! = v :=
      getElem!_set!_self _ _ _ (by rw [hsz]; exact hdlt)
    have hstd : (sts ++ [n])[d]! = n := by rw [← hls]; exact getElem!_concat_length' _ _
    have hord : (origs ++ [s.tstack.length])[d]! = s.tstack.length := by
      rw [← hlo]; exact getElem!_concat_length' _ _
    have hsts' : ∀ k, k ≤ d → (sts ++ [n])[k]! ≤ (sts ++ [n])[d]! := by
      intro k hk
      rw [hstd]
      rcases Nat.lt_or_ge k d with hk' | hk'
      · rw [getElem!_append_left' (by omega)]; exact hsts k hk'
      · obtain rfl : k = d := by omega
        rw [hstd]
    have horigs' : ∀ k, k ≤ d → (origs ++ [s.tstack.length])[k]! ≤ (origs ++ [s.tstack.length])[d]! := by
      intro k hk
      rw [hord]
      rcases Nat.lt_or_ge k d with hk' | hk'
      · rw [getElem!_append_left' (by omega)]; exact horigs k hk'
      · obtain rfl : k = d := by omega
        rw [hord]
    obtain ⟨hcov, hwp⟩ := cvOuts g σ n v d outs false anc (sts ++ [n]) (origs ++ [s.tstack.length]) P X _
      ⟨hi v outs rfl, hs.frame', hσ⟩ hnd hg hb hf hc
      (hpl.of_le rfl (Nat.le_refl _) (fun _ _ => rfl) fun _ => Nat.le_refl _) hP (ho v outs rfl)
      (by simp [hls]) (by simp [hlo]) hsts' (by rw [hstd]) horigs'
      hd hnodup hanc hv (fun w hw' => hw w (List.mem_cons_of_mem _ hw')) he hen
      (fun w hw' => hPv w (List.mem_cons_of_mem _ hw')) hPe (fun _ => hPv v (List.mem_cons_self ..)) hat
      hsz₀ (fun k hk => by rw [getElem!_set!_ne _ _ _ _ (by omega)]; exact hsvk k hk) hsvd₀
      (fun h => by cases h)
    refine ⟨hcov, ?_⟩
    refine wp_mono _ hwp fun hv' s' hpost => ?_
    obtain ⟨hpl', ho', hsv', hvc', -⟩ := hpost
    cases hv'
    · simp only [Bool.false_eq_true, ↓reduceIte, wp_bind, wp_setStackDir, wp_pushVertTstack, wp_pure]
      have hP' : ¬ Pushed g (fun i => P i ∨ (false = true ∧ i = vertItem v))
          (DfsOut.vertsList outs) (DfsOut.edgesList outs) (vertItem v) := by
        rintro ((hi | ⟨hb, _⟩) | ⟨w, hw, hwe⟩ | ⟨e, _, hev⟩)
        · exact hPv v (List.mem_cons_self ..) hi
        · cases hb
        · exact hvo (by rw [vertItem_inj hwe]; exact hw)
        · exact vertItem_ne_edgeItem hv e hev
      refine ⟨?_, OwnedD.mono ?_ fun i h => Pushed.vert (b := false) i h,
        fun _ _ h => by cases h; exact ⟨hsv', vertCover_pushVert rfl⟩⟩
      swap
      · exact ho'.pushVert hsv' hsts' rfl rfl rfl rfl
      refine ((hpl'.set_stackDir (s'.stackDir.set! d true)).cons_fixed (by show 0 < 1 + v; omega)
        (by show 1 + v < _; omega) hP' v d s'.nxtEdgeIdx (s'.stackDir.set! d true)[d]!).mono
        (fun i hi => ?_) fun _ hi => hi
      rcases hi with hi | rfl
      · exact Pushed.vert i hi
      · exact Or.inr (Or.inl ⟨v, List.mem_cons_self .., rfl⟩)
    · exact ⟨hpl'.mono (fun i h => Pushed.vert i h) fun _ h => h, ho'.mono fun i h => Pushed.vert i h,
        fun _ _ h => by cases h; exact ⟨hsv', hvc' rfl⟩⟩

theorem cvOuts : ∀ (g : Graph) (σ : List Nat) (n v d : Nat) (outs : List DfsOut) (hasVert : Bool)
    (anc sts origs : List Nat) (P X : ItemId → Prop) (s : WalkState),
    CvOuts g σ n v d outs hasVert anc sts origs P X s
  | g, σ, n, v, d, [], hasVert, anc, sts, origs, P, X, s =>
    fun _ _ _ _ _ _ hpl hP ho _ _ _ _ _ _ _ _ hv _ _ _ _ _ _ _ _ _ hsvd hvc => by
    unfold walkOuts
    simp only [wp_pure]
    unfold CoverOuts
    rw [show (DfsOut.edgePostorderList []).length = 0 by simp [DfsOut.edgePostorderList]]
    exact ⟨fun _ => hpl.pushVertR hv hP,
      hpl.mono (fun i h => Or.inl (Or.inl h)) fun _ h => h, ho.mono fun i h => Or.inl (Or.inl h), hsvd, hvc,
      fun h => h⟩
  | g, σ, n, v, d, o :: rest, hasVert, anc, sts, origs, P, X, s =>
    fun hrs hnd hg hb hf hc hpl hP ho hls hlo hsts hn horigs hd hnodup hanc hv hw he hen hPv hPe hPf hat
      hsz hsvk hsvd hvc => by
    unfold GuardsOuts at hg; unfold BookOuts at hb; unfold FrontiersOuts at hf
    unfold CsOuts at hc
    rw [DfsOut.vertsList_eq, List.flatMap_cons, ← DfsOut.vertsList_eq] at hnodup hw hPv ⊢
    rw [DfsOut.edgesList_eq, List.flatMap_cons, ← DfsOut.edgesList_eq] at he hen hPe ⊢
    rw [DfsOut.edgePostorderList_cons, List.length_append, ← Nat.add_assoc]
    rw [DfsOut.edgePostorderList_cons] at hat
    unfold walkOuts
    simp only [wp_bind]
    unfold CoverOuts
    have hT : Types g s := hpl.types
    have hK := kOut v d o hasVert s g (d + 1) 0 s hT (Nat.le_refl _) (by omega) (vertItem_ne_zero v)
      (fun w _ => vertItem_ne_zero w) (fun e _ => edgeItem_ne_zero g e) Keep.refl
    have hcons := List.nodup_cons.1 (List.nodup_append.1 hnodup).2.1
    have hdisjV := (List.nodup_append.1 hcons.2).2.2
    have hdisjE := (List.nodup_append.1 hen).2.2
    obtain ⟨hcov, hwp⟩ := cvOut g σ n v d o hasVert anc sts origs P X s hrs hnd hg.1 hb.1 hf.1 hc.1
      hpl hP ho hls hlo hsts hn horigs hd
      (hnodup.sublist (List.Sublist.append_left (List.Sublist.cons_cons v (List.sublist_append_left _ _)) _))
      hanc hv (fun w hw' => hw w (List.mem_append_left _ hw')) (fun e he' => he e (List.mem_append_left _ he'))
      (hen.sublist (List.sublist_append_left _ _))
      (fun w hw' => hPv w (List.mem_append_left _ hw')) (fun e he' => hPe e (List.mem_append_left _ he'))
      hPf hat.left hsz hsvk hsvd hvc
    have hboth : wp (walkOut v d o hasVert) (fun hv₁ s₁ =>
        CoverOuts σ (n + o.block.length) v d rest hv₁ s₁ ∧
        wp (walkOuts v d rest hv₁) (CvPost g σ (n + o.block.length + (DfsOut.edgePostorderList rest).length) v d
          sts origs P X hasVert (o.verts ++ DfsOut.vertsList rest) (o.edges ++ DfsOut.edgesList rest)) s₁) s := by
      refine wp_mono _ (wp_and (rgOut σ n v d o hasVert s hrs hnd hg.1 hb.1
        (scheduleOut σ n v d o hasVert s hrs hnd hg.1 hb.1 hf.1 hcov hat.left)) (wp_and hg.2 (wp_and hb.2
      (wp_and hf.2 (wp_and hc.2 (wp_and hK hwp))))))
        fun hv₁ s₁ ⟨hrs₁, hg₁, hb₁, hf₁, hc₁, hK₁, hpost₁⟩ => ?_
      obtain ⟨hpl₁, ho₁, hsv₁, hvc₁, hhv₁⟩ := hpost₁
      have hhf : hv₁ = false → hasVert = false := fun h => by cases hasVert <;> simp_all
      have hP₁ : ∀ e, e < g.ne →
          Pushed g (fun i => P i ∨ (hv₁ = true ∧ i = vertItem v)) o.verts o.edges (edgeItem g e) →
          σ.idxOf e < n + o.block.length := by
        refine pushed_past (fun e he' h => ?_) (fun w hw' => hw w (List.mem_append_left _ hw'))
          (DfsOut.edges_subset_block o) hat.left hnd
        rcases h with h | ⟨_, h⟩
        · exact hP e he' h
        · exact absurd h.symm (vertItem_ne_edgeItem hv e)
      have hPv₁ : ∀ w ∈ DfsOut.vertsList rest,
          ¬ Pushed g (fun i => P i ∨ (hv₁ = true ∧ i = vertItem v)) o.verts o.edges (vertItem w) := by
        rintro w hw' ((h | ⟨_, h⟩) | ⟨w', hw'', h⟩ | ⟨e, _, h⟩)
        · exact hPv w (List.mem_append_right _ hw') h
        · exact hcons.1 (by rw [← vertItem_inj h]; exact List.mem_append_right _ hw')
        · exact hdisjV w' hw'' w hw' (vertItem_inj h).symm
        · exact vertItem_ne_edgeItem (hw w (List.mem_append_right _ hw')) e h
      have hPe₁ : ∀ e ∈ DfsOut.edgesList rest,
          ¬ Pushed g (fun i => P i ∨ (hv₁ = true ∧ i = vertItem v)) o.verts o.edges (edgeItem g e) := by
        rintro e he' ((h | ⟨_, h⟩) | ⟨w', hw'', h⟩ | ⟨e', he'', h⟩)
        · exact hPe e (List.mem_append_right _ he') h
        · exact vertItem_ne_edgeItem hv e h.symm
        · exact vertItem_ne_edgeItem (hw w' (List.mem_append_left _ hw'')) e h.symm
        · exact hdisjE e' he'' e he' (edgeItem_inj h).symm
      have hPf₁ : hv₁ = false →
          ¬ Pushed g (fun i => P i ∨ (hv₁ = true ∧ i = vertItem v)) o.verts o.edges (vertItem v) := by
        rintro h₀ ((h | ⟨h1, _⟩) | ⟨w', hw'', h⟩ | ⟨e, _, h⟩)
        · exact hPf (hhf h₀) h
        · rw [h₀] at h1; cases h1
        · exact hcons.1 (by rw [vertItem_inj h]; exact List.mem_append_left _ hw'')
        · exact vertItem_ne_edgeItem hv e h
      obtain ⟨hcov₁, hwp₁⟩ := cvOuts g σ (n + o.block.length) v d rest hv₁ anc sts origs
        (Pushed g (fun i => P i ∨ (hv₁ = true ∧ i = vertItem v)) o.verts o.edges) X s₁
        hrs₁ hnd hg₁ hb₁ hf₁ hc₁ hpl₁ hP₁ ho₁ hls hlo hsts (Nat.le_trans hn (Nat.le_add_right _ _)) horigs hd
        (hnodup.sublist (List.Sublist.append_left (List.Sublist.cons_cons v (List.sublist_append_right _ _)) _))
        hanc hv (fun w hw' => hw w (List.mem_append_right _ hw')) (fun e he' => he e (List.mem_append_right _ he'))
        (hen.sublist (List.sublist_append_right _ _)) hPv₁ hPe₁ hPf₁ hat.right
        (hK₁.sv.trans hsz) (fun k hk => by rw [hK₁.svlo k (by omega)]; exact hsvk k hk) hsv₁ hvc₁
      refine ⟨hcov₁, wp_mono _ hwp₁ fun hv' s' hp => ?_⟩
      obtain ⟨hpl', ho', hsv', hvc', hhv'⟩ := hp
      exact ⟨hpl'.mono (fun i h => Pushed.append hhv' i h) fun _ h => h,
        ho'.mono fun i h => Pushed.append hhv' i h, hsv', hvc', fun h => hhv' (hhv₁ h)⟩
    exact ⟨⟨hcov, wp_mono _ hboth fun _ _ h => h.1⟩, wp_mono _ hboth fun _ _ h => h.2⟩

theorem cvOut : ∀ (g : Graph) (σ : List Nat) (n v d : Nat) (o : DfsOut) (hasVert : Bool)
    (anc sts origs : List Nat) (P X : ItemId → Prop) (s : WalkState),
    CvOut g σ n v d o hasVert anc sts origs P X s
  | g, σ, n, v, d, o, hasVert, anc, sts, origs, P, X, s =>
    fun ⟨hi, hs, hσ⟩ hnd hg hb hf hc hpl hP ho hls hlo hsts hn horigs hd hnodup hanc hv hw he hen hPv hPe
      hPf hat hsz hsvk hsvd hvc => by
    unfold GuardsOut at hg; unfold BookOut at hb; unfold FrontiersOut at hf
    unfold CsOut at hc
    have hT : Types g s := hpl.types
    have hdisjA := (List.nodup_append.1 hnodup).2.2
    have hvanc : v ∉ anc := fun hva => hdisjA v hva v (List.mem_cons_self ..) rfl
    have hpre := wp_and (walkOutPre_ranges hi hs hσ hb.1 fun _ => hpl.pushVertR hv hP) (wp_and hg (wp_and hb.2
      (wp_and hf (wp_and hc (wp_and (walkOutPre_place hpl hv hPf)
      (wp_and (walkOutPre_hv v d o hasVert) (wp_and (walkOutPre_ownedD d o hasVert ho hsvd hsts hvc)
      (keep_walkOutPre (D := d + 1) (j := 0) v d o hasVert (Keep.refl (s := s))))))))))
    rw [walkOut_eq, wp_bind]
    cases o with
    | back e dest cls =>
      simp only [DfsOut.edges, List.mem_singleton, forall_eq] at he hPe
      refine ⟨coverOut_back v d e dest cls hasVert hpl hv hPf hP (wp_mono _ hpre
        fun hv₁ s₁ ⟨⟨hi₁, hs₁, hσ₁⟩, hg₁, hb₁, hf₁, hc₁, ⟨hpl₁, hhv₁⟩, ⟨hhv₁', _⟩, ⟨ho₁, hvc₁⟩, hK₁⟩ =>
          ?_), ?_⟩
      · try simp only [wp_pure] at hg₁ hb₁ hf₁ hc₁
        exact finishP_ownership (D := d)
          ⟨by rw [if_neg fun h => by obtain ⟨_, _, _, h⟩ := hb₁.tree.1 h; cases h],
            hnd, ⟨hi₁, hs₁, hσ₁⟩, hg₁, hb₁, hc₁.pos⟩
          (ho₁.toOwned (by rw [hK₁.svlo d (Nat.lt_succ_self _)]; exact hsvd)) (hsts 0 (Nat.zero_le _)) hsts hhv₁'
      · refine wp_mono _ hpre
          fun hv₁ s₁ ⟨⟨hi₁, hs₁, hσ₁⟩, hg₁, hb₁, hf₁, hc₁, ⟨hpl₁, hhv₁⟩, ⟨hhv₁', _⟩, ⟨ho₁, hvc₁⟩, hK₁⟩ =>
            ?_
        try simp only [wp_pure] at hg₁ hb₁ hf₁ hc₁
        have hbase : OwnSite σ n d v d (.back e dest cls) s₁.tstack.length hv₁ s₁ :=
          ⟨by rw [if_neg fun h => by obtain ⟨_, _, _, h⟩ := hb₁.tree.1 h; cases h],
            hnd, ⟨hi₁, hs₁, hσ₁⟩, hg₁, hb₁, hc₁.pos⟩
        have hsv₁ : s₁.stackVerts[d]! = v := by rw [hK₁.svlo d (Nat.lt_succ_self _)]; exact hsvd
        have hhf : hv₁ = false → hasVert = false := fun h => by cases hasVert <;> simp_all
        have hpath : ∀ k, k < d → s₁.stackVerts[k]! ≠ v := fun k hk heq => by
          rw [hK₁.svlo k (by omega)] at heq
          exact hvanc (heq ▸ List.mem_of_getElem? (hsvk k hk))
        unfold walkOutRest
        rw [wp_bind, wp_tstackSize]
        simp only [wp_bind, wp_pure]
        have hfo : wp (finishEdge v d (.back e dest cls) s₁.tstack.length hv₁)
            (fun _ s' => OwnedD σ sts origs P d (n + 1) s') s₁ :=
          finishEdge_ownedD hbase ho₁ hsts hn horigs (ho₁.len d (Nat.le_refl _)) hpath
        refine wp_mono _ (wp_and (finishEdge_place v d (.back e dest cls) _ hv₁ hpl₁ hv he
            (fun h => ?_) (fun h₀ h => ?_)) (wp_and hfo (wp_and
            (keep_finishEdge (D := d + 1) hpl₁.types (by omega) v d (.back e dest cls) _ hv₁ (vertItem_ne_zero v)
              (edgeItem_ne_zero g _) Keep.refl)
            (finishEdge_vertCover hbase))))
          fun hv' s' ⟨⟨hpl', hhv'⟩, ho', hK', hvc'⟩ => ?_
        · rcases h with h | ⟨_, h⟩
          · exact hPe h
          · exact vertItem_ne_edgeItem hv e h.symm
        · rcases h with h | ⟨h1, _⟩
          · exact hPf (hhf h₀) h
          · rw [h₀] at h1; cases h1
        · refine ⟨hpl'.mono (fun i h => ?_) fun _ h => h, ho'.mono fun i h => Or.inl (Or.inl h),
            by rw [hK'.svlo d (Nat.lt_succ_self _)]; exact hsv₁, hvc', fun h => hhv' (hhv₁ h)⟩
          rcases h with (h | ⟨h1, rfl⟩) | rfl | ⟨h2, rfl⟩
          · exact Or.inl (Or.inl h)
          · exact Or.inl (Or.inr ⟨hhv' h1, rfl⟩)
          · exact Or.inr (Or.inr ⟨e, List.mem_singleton_self e, rfl⟩)
          · exact Or.inl (Or.inr ⟨h2, rfl⟩)
    | tree e cls child =>
      obtain ⟨w, outs⟩ := child
      simp only [DfsOut.verts] at hnodup hw hPv
      simp only [DfsOut.edges, List.mem_cons, forall_eq_or_imp] at he hPe
      simp only [DfsOut.edges] at hen
      simp only [DfsOut.block] at hat
      rw [show (DfsOut.tree e cls (.node w outs)).block.length = (DfsTree.node w outs).edgePostorder.length + 1
        by simp [DfsOut.block], ← Nat.add_assoc]
      have hcv := List.nodup_cons.1 (List.nodup_append.1 hnodup).2.1
      have hdisjC := (List.nodup_append.1 hnodup).2.2
      have hce := List.nodup_cons.1 hen
      have hwc : w ∈ (DfsTree.node w outs).verts := by simp [DfsTree.verts]
      have hdlt : d + 1 < g.nv := by
        have hsub : anc ++ [v] ++ [w] ⊆ List.range g.nv := fun a ha => by
          rw [List.mem_range]
          simp only [List.mem_append, List.mem_singleton] at ha
          rcases ha with (ha | rfl) | rfl
          · exact hanc a ha
          · exact hv
          · exact hw _ hwc
        have hnd' : (anc ++ [v] ++ [w]).Nodup := by
          refine hnodup.sublist ?_
          rw [List.append_assoc, List.singleton_append]
          exact List.Sublist.append_left (List.Sublist.cons_cons v (List.singleton_sublist.2 hwc)) _
        have := (List.Nodup.subperm hnd' hsub).length_le
        simp at this; omega
      have hPv' : ∀ w' ∈ (DfsTree.node w outs).verts, ∀ b : Bool,
          ¬ (P (vertItem w') ∨ (b = true ∧ vertItem w' = vertItem v)) := by
        rintro w' hw' b (h | ⟨_, h⟩)
        · exact hPv w' hw' h
        · exact hcv.1 (by rw [← vertItem_inj h]; exact hw')
      have hPe' : ∀ e' ∈ (DfsTree.node w outs).edges, ∀ b : Bool,
          ¬ (P (edgeItem g e') ∨ (b = true ∧ edgeItem g e' = vertItem v)) := by
        rintro e' he' b (h | ⟨_, h⟩)
        · exact hPe.2 e' he' h
        · exact vertItem_ne_edgeItem hv e' h.symm
      have hP' : ∀ b : Bool, ∀ e', e' < g.ne →
          (P (edgeItem g e') ∨ (b = true ∧ edgeItem g e' = vertItem v)) → σ.idxOf e' < n := by
        rintro b e' he' (h | ⟨_, h⟩)
        · exact hP e' he' h
        · exact absurd h.symm (vertItem_ne_edgeItem hv e')
      have hnσ : n + (DfsTree.node w outs).edgePostorder.length ≤ σ.length := by
        obtain ⟨pre', post', hl, hσe⟩ := hat
        rw [hσe]; simp only [List.length_append, List.length_singleton]; omega
      have hchild : wp (walkOutPre v d (.tree e cls (.node w outs)) hasVert) (fun hv₁ s₁ =>
          (hasVert = true → hv₁ = true) ∧
          CoverTree σ n (.node w outs) (d + 1) { s₁ with firstOccurrence := s₁.firstOccurrence.set! d s₁.g.ne } ∧
          wp (walkTree (.node w outs) (d + 1)) (fun _ s₃ =>
            OwnSite σ (n + (DfsTree.node w outs).edgePostorder.length) (d + 1) v d
              (.tree e cls (.node w outs)) s₁.tstack.length hv₁ s₃ ∧
            s₃.Place g (Pushed g (fun i => P i ∨ (hv₁ = true ∧ i = vertItem v))
              (DfsTree.node w outs).verts (DfsTree.node w outs).edges) X ∧
            OwnedD σ sts origs (Pushed g (fun i => P i ∨ (hv₁ = true ∧ i = vertItem v))
              (DfsTree.node w outs).verts (DfsTree.node w outs).edges) d
              (n + (DfsTree.node w outs).edgePostorder.length) s₃ ∧
            s₃.stackVerts[d]! = v ∧ (∀ k, k < d → s₃.stackVerts[k]! ≠ v) ∧
            ((DfsOut.tree e cls (.node w outs)).cls.lowval d < d →
              (DfsOut.tree e cls (.node w outs)).cls.isType1 = true → hv₁ = true) ∧
            origs[d]! ≤ s₁.tstack.length)
            { s₁ with firstOccurrence := s₁.firstOccurrence.set! d s₁.g.ne }) s := by
        refine wp_mono _ hpre
          fun hv₁ s₁ ⟨⟨hi₁, hs₁, hσ₁⟩, hg₁, hb₁, hf₁, hc₁, ⟨hpl₁, hhv₁⟩, ⟨hhv₁', _⟩, ⟨ho₁, hvc₁⟩, hK₁⟩ =>
            ?_
        try simp only [wp_bind, wp_modify] at hg₁ hb₁ hf₁ hc₁
        have hg₁' : s₁.g = g := hK₁.g.trans hT.g_eq
        have hsv₁ : s₁.stackVerts[d]! = v := by rw [hK₁.svlo d (Nat.lt_succ_self _)]; exact hsvd
        have hsvk₁ : ∀ k, k < d → anc[k]? = some s₁.stackVerts[k]! := fun k hk => by
          rw [hK₁.svlo k (by omega)]; exact hsvk k hk
        have hpath₁ : ∀ k, k ≤ d → w ≠ s₁.stackVerts[k]! := by
          intro k hk heq
          rcases Nat.lt_or_ge k d with hk' | hk'
          · exact hdisjC _ (List.mem_of_getElem? (hsvk₁ k hk')) w (List.mem_cons_of_mem _ hwc) heq.symm
          · obtain rfl : k = d := by omega
            rw [hsv₁] at heq
            exact hcv.1 (by rw [← heq]; exact hwc)
        have pre : ∀ w' outs', DfsTree.node w outs = .node w' outs' →
            ({ { s₁ with firstOccurrence := s₁.firstOccurrence.set! d s₁.g.ne } with
              stackVerts := s₁.stackVerts.set! (d + 1) w' } : WalkState).RangesInv σ n (d + 1) :=
          fun w' outs' _ =>
            (hi₁.frame' (s' := { s₁ with firstOccurrence := s₁.firstOccurrence.set! d s₁.g.ne })).setSv w'
        have hpl₂ : ({ s₁ with firstOccurrence := s₁.firstOccurrence.set! d s₁.g.ne } : WalkState).Place g
            (fun i => P i ∨ (hv₁ = true ∧ i = vertItem v)) X :=
          hpl₁.of_le rfl (Nat.le_refl _) (fun _ _ => rfl) fun _ => Nat.le_refl _
        obtain ⟨hcovc, hwpc⟩ := cvTree g σ n (.node w outs) (d + 1) (anc ++ [v]) sts origs
          (fun i => P i ∨ (hv₁ = true ∧ i = vertItem v)) X _ pre hs₁.frame' hnd hσ₁ hg₁.1 hb₁.1 hf₁.1 hc₁.1
          hpl₂ (hP' hv₁)
          (fun w' outs' h => by
            cases h
            exact (ho₁.mono fun i h => Or.inl h).entry hls hlo (by rw [hg₁']; exact hw w hwc)
              (by rw [hK₁.sv, hsz]; exact hdlt) (hPv' w hwc hv₁) hpath₁ rfl rfl rfl rfl)
          hls hlo (fun k hk => Nat.le_trans (hsts k (by omega)) hn) (fun k hk => ho₁.len k (by omega))
          (by simp [hd]) (by simpa using hnodup)
          (fun a ha => by
            rcases List.mem_append.1 ha with ha | ha
            · exact hanc a ha
            · exact (List.mem_singleton.1 ha) ▸ hv)
          hw he.2 hce.2 (fun w' hw' => hPv' w' hw' hv₁) (fun e' he' => hPe' e' he' hv₁) hat.left
          (hK₁.sv.trans hsz)
          (fun k hk => by
            rcases Nat.lt_or_ge k d with hk' | hk'
            · rw [List.getElem?_append_left (by omega)]
              exact hsvk₁ k hk'
            · rw [show k = anc.length by omega, List.getElem?_concat_length, ← hd]
              show some v = some s₁.stackVerts[d]!
              rw [hsv₁])
        refine ⟨hhv₁, hcovc, ?_⟩
        have hT₂ : Types g { s₁ with firstOccurrence := s₁.firstOccurrence.set! d s₁.g.ne } := hpl₂.types
        refine wp_mono _ (wp_and (rgTree σ n (.node w outs) (d + 1) _ pre hs₁.frame' hnd hσ₁ hg₁.1 hb₁.1
            (scheduleTree σ n (.node w outs) (d + 1) _ pre hs₁.frame' hnd hσ₁ hg₁.1 hb₁.1 hf₁.1 hcovc hat.left))
          (wp_and hg₁.2 (wp_and hb₁.2 (wp_and hf₁.2 (wp_and hc₁.2
          (wp_and (kTree (.node w outs) (d + 1) _ g (d + 1) 0 _ hT₂ (Nat.le_refl _) (by omega)
            (fun w' _ => vertItem_ne_zero w') (fun e' _ => edgeItem_ne_zero g e') Keep.refl) hwpc))))))
          fun _ s₃ ⟨⟨hi₃, hs₃, hσ₃⟩, hg₃, hb₃, hf₃, hc₃, hK₃, hpl₃, ho₃, hcc₃⟩ => ?_
        obtain ⟨hsvw, hvcw⟩ := hcc₃ w outs rfl
        have hsv₃ : s₃.stackVerts[d]! = v := by
          rw [hK₃.svlo d (by omega)]; exact hsv₁
        have hpath₃ : ∀ k, k < d → s₃.stackVerts[k]! ≠ v := fun k hk heq => by
          rw [hK₃.svlo k (by omega)] at heq
          exact hvanc (heq ▸ List.mem_of_getElem? (hsvk₁ k hk))
        refine ⟨⟨by rw [if_pos (hb₃.tree.2 ⟨_, _, _, rfl⟩)], hnd, ⟨hi₃, hs₃, hσ₃⟩, hg₃, hb₃, hc₃.pos⟩,
          hpl₃, ?_, hsv₃, hpath₃, hhv₁', ho₁.len d (Nat.le_refl _)⟩
        refine ho₃.exit hls hlo (Nat.le_add_right _ _) hσ₃ hnσ ?_ ?_
        · rw [hsvw]; exact hvcw
        · rw [hsvw]; exact Or.inr (Or.inl ⟨w, hwc, rfl⟩)
      refine ⟨coverOut_tree v d e cls (.node w outs) hasVert hpl hv hPf hw he.2 hcv.2 hce.2 hcv.1 hPv hPe.2 hP
        hat.left hnd (wp_mono _ hchild fun hv₁ s₁ ⟨_, hcovc, hwpc⟩ => ?_), ?_⟩
      · simp only [wp_modify]
        exact ⟨hcovc, wp_mono _ hwpc fun _ s₃ ⟨hbase, _, ho₂, hsv₃, _, hhv₁', _⟩ =>
          finishP_ownership hbase (ho₂.toOwned hsv₃) (hsts 0 (Nat.zero_le _)) hsts hhv₁'⟩
      · refine wp_mono _ hchild fun hv₁ s₁ ⟨hhv₁, _, hwpc⟩ => ?_
        unfold walkOutRest
        rw [wp_bind, wp_tstackSize]
        simp only [wp_bind, wp_modify]
        refine wp_mono _ hwpc fun _ s₃ ⟨hbase, hpl₃, ho₂, hsv₃, hpath₃, _, horig₁⟩ => ?_
        have hhf : hv₁ = false → hasVert = false := fun h => by cases hasVert <;> simp_all
        have hfo : wp (finishEdge v d (.tree e cls (.node w outs)) s₁.tstack.length hv₁)
            (fun _ s' => OwnedD σ sts origs (Pushed g (fun i => P i ∨ (hv₁ = true ∧ i = vertItem v))
              (DfsTree.node w outs).verts (DfsTree.node w outs).edges) d
              (n + (DfsTree.node w outs).edgePostorder.length + 1) s') s₃ :=
          finishEdge_ownedD hbase ho₂ hsts (Nat.le_trans hn (Nat.le_add_right _ _)) horigs
            horig₁ hpath₃
        refine wp_mono _ (wp_and (finishEdge_place v d (.tree e cls (.node w outs)) _ hv₁ hpl₃ hv he.1
            (fun h => ?_) (fun h₀ h => ?_)) (wp_and hfo (wp_and
            (keep_finishEdge (D := d + 1) hpl₃.types (by omega) v d (.tree e cls (.node w outs)) _ hv₁
              (vertItem_ne_zero v) (edgeItem_ne_zero g _) Keep.refl)
            (finishEdge_vertCover hbase))))
          fun hv' s' ⟨⟨hpl', hhv'⟩, ho', hK', hvc'⟩ => ?_
        · rcases h with (h | ⟨_, h⟩) | ⟨w', hw', h⟩ | ⟨e', he', h⟩
          · exact hPe.1 h
          · exact vertItem_ne_edgeItem hv e h.symm
          · exact vertItem_ne_edgeItem (hw w' hw') e h.symm
          · exact hce.1 (by rw [show e = e' from edgeItem_inj h]; exact he')
        · rcases h with (h | ⟨h1, _⟩) | ⟨w', hw', h⟩ | ⟨e', _, h⟩
          · exact hPf (hhf h₀) h
          · rw [h₀] at h1; cases h1
          · exact hcv.1 (by rw [vertItem_inj h]; exact hw')
          · exact vertItem_ne_edgeItem hv e' h
        · refine ⟨hpl'.mono (fun i h => ?_) fun _ h => h, ho'.mono fun i h => Pushed.finish hhv' i (Or.inl h),
            by rw [hK'.svlo d (Nat.lt_succ_self _)]; exact hsv₃, hvc', fun h => hhv' (hhv₁ h)⟩
          exact Pushed.finish hhv' i h
end

theorem EdgeBelow_vert_nil {items : Items} {g : Graph} {v e : Nat} (hv : v < g.nv)
    (hch : Items.ch items (vertItem v) = []) (h : Items.EdgeBelow g items (vertItem v) e) : False := by
  rcases Relation.ReflTransGen.cases_head h with heq | ⟨c, hc, -⟩
  · exact vertItem_ne_edgeItem hv e heq
  · simp [Items.IsParent, hch] at hc

/-- Ownership at the entry of a root tree: nothing is owned yet, the stack is empty and every
unvisited vertex item is childless (`RootState.fresh`). -/
theorem OwnedD.root {g : Graph} {pre : List DfsTree} (h : RootState g pre s) {v : Nat} (hv : v < g.nv)
    (hPv : ¬ Pushed g (fun _ => False) (pre.flatMap DfsTree.verts) (pre.flatMap DfsTree.edges) (vertItem v))
    (σ : List Nat) (n : Nat) :
    OwnedD σ [n] [s.tstack.length]
      (Pushed g (fun _ => False) (pre.flatMap DfsTree.verts) (pre.flatMap DfsTree.edges)) 0 n
      { s with stackVerts := s.stackVerts.set! 0 v } := by
  have hsv : (s.stackVerts.set! 0 v)[0]! = v := by
    simp [show 0 < s.stackVerts.size by rw [h.sv]; omega]
  have hch : ∀ w, w < g.nv →
      ¬ Pushed g (fun _ => False) (pre.flatMap DfsTree.verts) (pre.flatMap DfsTree.edges) (vertItem w) →
      Items.ch s.items (vertItem w) = [] :=
    fun w hw hPw => h.fresh (vertItem w) (by show 0 < 1 + w; omega) (by show 1 + w < _; omega)
      (fun v' hv' heq => hPw (Or.inr (Or.inl ⟨v', hv', heq⟩)))
      (fun e he heq => hPw (Or.inr (Or.inr ⟨e, he, heq⟩)))
  have hvb : ∀ e, Items.EdgeBelow s.g s.items (vertItem v) e → False := fun e hb => by
    rw [h.g_eq] at hb
    exact EdgeBelow_vert_nil hv (hch v hv hPv) hb
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · intro k hk; obtain rfl : k = 0 := Nat.le_zero.1 hk; simp
  · intro b hb hb'; simp at hb; omega
  · intro k hk; exact absurd hk (Nat.not_lt_zero _)
  · intro e _ hb; rw [hsv] at hb; exact (hvb e hb).elim
  · intro k hk e _ hb; obtain rfl : k = 0 := Nat.le_zero.1 hk; rw [hsv] at hb; exact (hvb e hb).elim
  · intro k _ t ht; simp [h.tstack] at ht
  · intro k _ t ht; simp [h.tstack] at ht
  · intro t ht; simp [h.tstack] at ht
  · intro w hw hPw _
    exact hch w (by rw [← h.g_eq]; exact hw) hPw

/-- `RootsCover` over the DFS forest: each root tree is covered by `cvTree`, and the root
append (`RootState.step`, `RangesInv.root_append`) re-establishes the root-state premises. -/
theorem forest_rootsCover {g : Graph} {σ : List Nat} (hnd : σ.Nodup) (hσ : ∀ e ∈ σ, e < g.ne) :
    ∀ forest pre n s, RootState g pre s → s.RangesInv σ n 0 → ForestOK g (pre ++ forest) →
      (∀ t ∈ forest, t.WF []) → (∀ t ∈ forest, t.Ends g) →
      (∀ t ∈ forest, ∀ e, e < g.ne → ∀ x, g.Inc e x → x ∈ t.verts → e ∈ t.edges) →
      PostAt σ n (edgePostorderForest forest) →
      (∀ e ∈ pre.flatMap DfsTree.edges, σ.idxOf e < n) →
      RootsCover σ n forest s
  | [], _, _, _, _, _, _, _, _, _, _, _ => trivial
  | t :: rest, pre, n, s, h, hr, hf, hwf, hends, hcomp, hat, hpre => by
    have hb := h.book hf (hwf t (by simp)) (hends t (by simp)) (hcomp t (by simp))
    have hg := gbTree t 0 s hb
    have hi : ∀ v outs, t = .node v outs →
        ({ s with stackVerts := s.stackVerts.set! 0 v } : WalkState).RangesInv σ n 0 :=
      fun v _ _ => hr.stackVerts_of_nil h.tstack _
    have hσ' : ∀ e ∈ σ, e < s.g.ne := by rwa [h.g_eq]
    have hfront := walkTree_frontiers t 0 s (fun v outs ht => (hi v outs ht).inv) h.shape hg hb
    have hat' : PostAt σ n (t.edgePostorder ++ edgePostorderForest rest) := hat
    have hcs : CsTree σ n t 0 s :=
      csTree_of_dfs (anc := []) t 0 s h.place.types rfl (hwf t (by simp)) (hends t (by simp))
        (hcomp t (by simp)) (by simpa using (RootState.hvn hf).1) (by simp) (RootState.hvlt hf) (RootState.helt hf)
        (RootState.hen hf).1 h.sv (fun k hk => absurd hk (Nat.not_lt_zero _)) hnd hat'.left
        (fun v outs _ k hk => absurd hk (by simp))
    have hP : ∀ e, e < g.ne →
        Pushed g (fun _ => False) (pre.flatMap DfsTree.verts) (pre.flatMap DfsTree.edges) (edgeItem g e) →
        σ.idxOf e < n := by
      rintro e _ (h | ⟨w, hw, hew⟩ | ⟨e', he', hee⟩)
      · exact h.elim
      · refine absurd hew.symm (vertItem_ne_edgeItem (hf.verts_lt w ?_) e)
        rw [List.flatMap_append]; exact List.mem_append_left _ hw
      · rw [edgeItem_inj hee]; exact hpre e' he'
    obtain ⟨hcov, hwp⟩ := cvTree g σ n t 0 [] [] [] _ (fun _ => False) s hi h.shape hnd hσ' hg hb hfront hcs
      h.place hP
      (fun v outs ht =>
        OwnedD.root h (RootState.hvlt hf v (by subst ht; simp [DfsTree.verts]))
          (RootState.hPv hf v (by subst ht; simp [DfsTree.verts])) σ n)
      rfl rfl (fun k hk => absurd hk (Nat.not_lt_zero _)) (fun k hk => absurd hk (Nat.not_lt_zero _)) rfl
      (by simpa using (RootState.hvn hf).1) (by simp)
      (RootState.hvlt hf) (RootState.helt hf) (RootState.hen hf).1 (RootState.hPv hf) (RootState.hPe hf)
      hat'.left h.sv (fun k hk => absurd hk (Nat.not_lt_zero _))
    have hrg := rgTree σ n t 0 s hi h.shape hnd hσ' hg hb
      (scheduleTree σ n t 0 s hi h.shape hnd hσ' hg hb hfront hcov hat'.left)
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
    have hpre' : ∀ e ∈ (pre ++ [t]).flatMap DfsTree.edges, σ.idxOf e < n + t.edgePostorder.length := by
      intro e he
      simp only [List.flatMap_append, List.flatMap_cons, List.flatMap_nil, List.append_nil] at he
      rcases List.mem_append.1 he with he | he
      · exact Nat.lt_of_lt_of_le (hpre e he) (Nat.le_add_right _ _)
      · exact pushed_past hP (RootState.hvlt hf)
          (fun e he => (DfsTree.edgePostorder_perm_edges t).mem_iff.2 he) hat'.left hnd e
          (RootState.helt hf e he) (Or.inr (Or.inr ⟨e, he, rfl⟩))
    refine ⟨hcov, ?_⟩
    refine wp_mono _ (wp_and hrg (wp_and hk (wp_and hp hst))) fun _ s₁ ⟨hrs, hk, hp, hst⟩ => ?_
    have hrpop := hrs.1.root_append hrs.2.1 hk (noParent_of_cnt_eq_zero hp.root)
    exact wp_mono (popTstack >>= fun top => modifyItem rootItem fun it => { it with ch := it.ch ++ top.spans.2 })
      (wp_and hrpop hst) fun _ s₂ ⟨hr₂, hs₂⟩ =>
      forest_rootsCover hnd hσ rest (pre ++ [t]) (n + t.edgePostorder.length) s₂ hs₂ hr₂
        (by simpa using hf) (fun t' ht' => hwf t' (by simp [ht']))
        (fun t' ht' => hends t' (by simp [ht'])) (fun t' ht' => hcomp t' (by simp [ht']))
        hat'.right hpre'

/-- `RootsCover` of the DFS forest walk from the initial state. -/
theorem walk_rootsCover' (g : Graph) (tern : Bool) (forest : List DfsTree)
    (hf : ForestOK g forest) (hwf : ∀ t ∈ forest, t.WF []) (hends : ∀ t ∈ forest, t.Ends g)
    (hecov : ∀ e, e < g.ne → e ∈ forest.flatMap DfsTree.edges) :
    RootsCover (edgePostorderForest forest) 0 forest (WalkState.init g tern) := by
  have hnd := DfsData.edgePostorderForest_perm.nodup_iff.2 hf.edges_nodup
  have hσ : ∀ e ∈ edgePostorderForest forest, e < g.ne :=
    fun e he => hf.edges_lt e (DfsData.edgePostorderForest_perm.subset he)
  exact forest_rootsCover hnd hσ forest [] 0 (WalkState.init g tern) (rootState_init g tern)
    (init_rangesInv g tern hnd) (by simpa using hf) hwf hends
    (fun t ht => comp_of_forest hf hwf hends hecov ht) ⟨[], [], rfl, by simp⟩ (by simp)

end WalkState
end Spqr
