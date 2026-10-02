import Spqr.EarInv

/-!
# Loop 1 from `EarFinish.loop1`

`loop1_bodyOk` derives `Loop1BodyOk` at every iterate of loop 1 (`closeEars`) from `Loop1Spec`
(the per-split merge/unwrap/close facts of `EarFinish.loop1`) and the span/edge disjointness of the
stack, by induction on the iterate with the invariant `L1At`: the top entry is the loop-1 piece
(bottom `l1Bot`, edges `l1Edges`, a single fresh one-sided item), the remaining entries are
untouched (`Below`/`ch`/parents/types as in the original state).
-/

namespace Spqr
open WalkM WalkState

theorem iter_succ' (body : WalkM Unit) (k : Nat) (s : WalkState) :
    iter body (k + 1) s = (body.run (iter body k s)).2 := by
  induction k generalizing s with
  | zero => rfl
  | succ k ih => exact ih _

theorem List.pairwise_split {α : Type} {R : α → α → Prop} {l l₁ l₂ : List α} (h : l.Pairwise R)
    (hl : l = l₁ ++ l₂) : ∀ a ∈ l₁, ∀ b ∈ l₂, R a b := by
  subst hl; exact (List.pairwise_append.1 h).2.2

theorem l1Bot_concat (o : DfsOut) (done : List TEntry) (t : TEntry) :
    WalkState.l1Bot o (done ++ [t]) = t.vStart := by
  simp [WalkState.l1Bot]

theorem l1Edges_concat (o : DfsOut) (s : WalkState) (done : List TEntry) (t : TEntry) (e : Nat) :
    WalkState.l1Edges o s (done ++ [t]) e ↔ WalkState.l1Edges o s done e ∨ t.edges s.g s.items e := by
  simp [WalkState.l1Edges, or_assoc]

namespace Items
variable {items : Items}

/-- A root item on no entry is not below any span item. -/
theorem not_below_of_root {i j : ItemId} (hroot : ∀ p, ¬ items.IsParent p j) (hne : i ≠ j) :
    ¬ items.Below i j := fun h => hne (Below.eq_of_no_parent hroot h)

end Items

namespace TEntry
variable {g : Graph} {items : Items}

/-- The closed entry `setSides dir [item] []` over `items.modify item (ch := getSide t.spans dir)`
holds exactly the edges of the one-sided `t`. -/
theorem edges_finish (dir : Bool) (t : TEntry) (item : ItemId) (vStart topDepth fi : Nat)
    (vsv : Option Nat × Option Nat)
    (hitem : item < items.size) (hroot : ∀ p, ¬ items.IsParent p item)
    (hfree : item ∉ t.spans.1 ++ t.spans.2) (hedge : ∀ e, e < g.ne → edgeItem g e ≠ item)
    (hside : getSide t.spans (!dir) = []) {e : Nat} (he : e < g.ne) :
    TEntry.edges g (items.modify item fun it => { it with vs := vsv, ch := getSide t.spans dir })
      ⟨vStart, topDepth, fi, setSides dir [item] []⟩ e ↔ t.edges g items e := by
  set f : Item → Item := fun it => { it with vs := vsv, ch := getSide t.spans dir }
  have hch : Items.ch (items.modify item f) item = getSide t.spans dir := by
    rw [Items.ch_modify_at item f hitem]
  have hne : ∀ c ∈ getSide t.spans dir, c ≠ item := fun c hc hce => by
    subst hce; exact hfree ((mem_of_getSide_nil dir t.spans hside c).2 hc)
  have hnb : ∀ c, c ≠ item → ¬ items.Below c item := fun c hc => Items.not_below_of_root hroot hc
  simp only [TEntry.edges, Items.EdgeBelow, mem_of_getSide_nil dir t.spans hside, mem_setSides,
    List.mem_singleton, exists_eq_left]
  constructor
  · intro hb
    rcases hb.head_cases with heq | ⟨c, hc, hb⟩
    · exact absurd heq.symm (hedge e he)
    · simp only [Items.IsParent, hch] at hc
      exact ⟨c, hc, (Items.Below_modify_of_not_below item f (hnb c (hne c hc))).1 hb⟩
  · rintro ⟨c, hc, hb⟩
    exact .head (by simpa [Items.IsParent, hch] using hc)
      ((Items.Below_modify_of_not_below item f (hnb c (hne c hc))).2 hb)

theorem getSide_mergeInto (dir : Bool) (cur nxt : TEntry) (h₁ : getSide cur.spans dir = [])
    (h₂ : getSide nxt.spans dir = []) : getSide (mergeInto cur nxt).spans dir = [] := by
  cases dir <;> simp_all [mergeInto, getSide]

theorem mem_mergeInto_spans (cur nxt : TEntry) (i : ItemId) :
    i ∈ (mergeInto cur nxt).spans.1 ++ (mergeInto cur nxt).spans.2 ↔
      i ∈ cur.spans.1 ++ cur.spans.2 ∨ i ∈ nxt.spans.1 ++ nxt.spans.2 := by
  simp only [mergeInto, List.mem_append]; tauto

end TEntry

namespace Graph
theorem Touches.congr {g : Graph} {E₁ E₂ : Nat → Prop} (h : ∀ e, e < g.ne → (E₁ e ↔ E₂ e)) {v : Nat} :
    g.Touches E₁ v ↔ g.Touches E₂ v := by
  constructor <;> rintro ⟨e, he, hE, hv⟩
  · exact ⟨e, he, (h e he).1 hE, hv⟩
  · exact ⟨e, he, (h e he).2 hE, hv⟩
end Graph

namespace WalkState
variable {D d : Nat} {s st : WalkState} {o : DfsOut}

/-- The remaining entries `R` of the stack are untouched since `s`: same descendants, children and
parents of their span items; all original item types are kept. -/
structure L1Keep (s st : WalkState) (R : List TEntry) : Prop where
  size : s.items.size ≤ st.items.size
  type : ∀ j, j < s.items.size → Items.type st.items j = Items.type s.items j
  below : ∀ u ∈ R, ∀ i ∈ u.spans.1 ++ u.spans.2, ∀ j, Items.Below st.items i j ↔ Items.Below s.items i j
  ch : ∀ u ∈ R, ∀ i ∈ u.spans.1 ++ u.spans.2, Items.ch st.items i = Items.ch s.items i
  parent : ∀ u ∈ R, ∀ i ∈ u.spans.1 ++ u.spans.2, ∀ p, Items.IsParent st.items p i → Items.IsParent s.items p i

theorem L1Keep.edges (h : L1Keep s st R) {u : TEntry} (hu : u ∈ R) (e : Nat) :
    u.edges s.g st.items e ↔ u.edges s.g s.items e :=
  TEntry.edges_congr (fun i hi _ => h.below u hu i hi _) e

theorem L1Keep.items_eq (h : L1Keep s st R) {st' : WalkState} (hi : st'.items = st.items) :
    L1Keep s st' R := by
  obtain ⟨size, type, below, ch, parent⟩ := h
  exact ⟨by rw [hi]; exact size, by rw [hi]; exact type, by rw [hi]; exact below, by rw [hi]; exact ch,
    by rw [hi]; exact parent⟩

theorem L1Keep.mono (h : L1Keep s st R) {R' : List TEntry} (hR : ∀ u ∈ R', u ∈ R) : L1Keep s st R' :=
  ⟨h.size, h.type, fun u hu => h.below u (hR u hu), fun u hu => h.ch u (hR u hu),
    fun u hu => h.parent u (hR u hu)⟩

/-- The loop-1 piece on top of the stack: bottom `v`, edges `E`, one single fresh item on side `dir`. -/
structure L1Piece (s st : WalkState) (d : Nat) (dir : Bool) (v : Nat) (E : Nat → Prop) (c : TEntry)
    (R : List TEntry) : Prop where
  tstack : st.tstack = c :: R
  top : c.topDepth = d
  bot : c.vStart = v
  side : getSide c.spans (!dir) = []
  root : ∀ i ∈ c.spans.1 ++ c.spans.2,
    (∀ p, ¬ Items.IsParent s.items p i) ∨ ∃ u ∈ s.tstack, i ∈ u.spans.1 ++ u.spans.2
  fresh : ∀ i ∈ c.spans.1 ++ c.spans.2, ∀ u ∈ R, i ∉ u.spans.1 ++ u.spans.2
  edges : ∀ e, e < s.g.ne → (c.edges s.g st.items e ↔ E e)

/-- The state facts loop 1 keeps relative to the state `s` before `finishEdge`. -/
structure L1Frame (s st : WalkState) (D d : Nat) (dir : Bool) : Prop where
  g : st.g = s.g
  sv : st.stackVerts = s.stackVerts
  dir : st.stackDir[d]! = dir
  inv : st.Inv' D
  shape : Shape st


/-- Closing the entry `t₂` (the unwrapped-or-not `t`) under the loop-1 piece `c`: merge, then finish
into the fresh `item`. -/
theorem l1_close_core {dir : Bool} {v v₀ : Nat} {E : Nat → Prop} {c t t₂ : TEntry} {R : List TEntry}
    {item : ItemId} {s₂ : WalkState}
    (hD : D = d + 1) (hF : L1Frame s s₂ D d dir) (hP : L1Piece s s₂ d dir v E c (t₂ :: R))
    (hK : L1Keep s s₂ R) (hC : L1Close d s v E t)
    (htop : d ≤ t₂.topDepth) (hbot : t₂.vStart = t.vStart) (hside : getSide t₂.spans (!dir) = [])
    (hT : ∀ e, e < s.g.ne → (t₂.edges s.g s₂.items e ↔ t.edges s.g s.items e))
    (hfr : ∀ i ∈ t₂.spans.1 ++ t₂.spans.2, ∀ u ∈ R, i ∉ u.spans.1 ++ u.spans.2)
    (hdisj : ∀ u ∈ R, ∀ e, e < s.g.ne → u.edges s.g s.items e → ¬ (E e ∨ t.edges s.g s.items e))
    (hfree : ItemFree s₂ item) (hiroot : ∀ p, ¬ Items.IsParent s.items p item) (hv : v₀ < s₂.g.nv) :
    CloseTwoOk D s₂ ∧
    ∃ c', L1Piece s ((finishTstackTop item).run (mergeTstackTops.run s₂).2).2 d dir t.vStart
        (fun e => E e ∨ t.edges s.g s.items e) c' R ∧
      L1Keep s ((finishTstackTop item).run (mergeTstackTops.run s₂).2).2 R ∧
      L1Frame s ((finishTstackTop item).run (mergeTstackTops.run s₂).2).2 D d dir := by
  subst hD
  have hts := hP.tstack
  have hg := hF.g
  have hsv := hF.sv
  have hmtop : (TEntry.mergeInto c t₂).topDepth = d := by
    show min t₂.topDepth c.topDepth = d; rw [hP.top]; exact Nat.min_eq_right htop
  have hmv : (TEntry.mergeInto c t₂).vStart = t.vStart := hbot
  have hmside : getSide (TEntry.mergeInto c t₂).spans (!dir) = [] :=
    TEntry.getSide_mergeInto _ _ _ hP.side hside
  have hmE : ∀ e, e < s.g.ne →
      ((TEntry.mergeInto c t₂).edges s.g s₂.items e ↔ E e ∨ t.edges s.g s.items e) := fun e he => by
    rw [TEntry.edges_mergeInto, hP.edges e he, hT e he]
  have hmerge : MergeTopOk (d + 1) s₂ := by
    intro cur nxt rest hr
    rw [hts] at hr
    simp only [List.cons.injEq] at hr
    obtain ⟨rfl, rfl, rfl⟩ := hr
    refine ⟨⟨fun _ hne => ?_, ?_⟩, fun u hu e he hue => ?_⟩
    · obtain ⟨e, he, het⟩ := hne
      rw [hg] at he het
      obtain ⟨x, hx₁, hx₂⟩ := hC.share ⟨e, he, (hT e he).1 het⟩
      rw [hg]
      exact ⟨x, (Graph.Touches.congr fun e he => (hP.edges e he).symm).1 hx₁,
        (Graph.Touches.congr fun e he => (hT e he).symm).1 hx₂⟩
    · rcases hC.bottom with h | ⟨k, hk₁, hk₂, hk⟩ | h
      · left; left; rw [hP.bot, h]; exact hbot.symm
      · left; right; exact ⟨k, by rw [hmtop]; exact hk₁, hk₂, by rw [hP.bot, hsv]; exact hk⟩
      · right; rw [hg, hP.bot]
        exact (Graph.Interior.congr fun e he => by rw [hP.edges e he, hT e he]).2 h
    · rw [hg] at he hue
      have := hdisj u hu e he ((hK.edges hu e).1 hue)
      rintro (h | h) <;> rw [hg] at h
      · exact this (.inl ((hP.edges e he).1 h))
      · exact this (.inr ((hT e he).1 h))
  have hrun₁ := mergeTstackTops_run_eq s₂ c t₂ R hts
  have hfin : FinishTopOk (d + 1) (mergeTstackTops.run s₂).2 := by
    rw [hrun₁]
    have hcur : curE { s₂ with tstack := TEntry.mergeInto c t₂ :: R } = TEntry.mergeInto c t₂ := rfl
    refine ⟨by simp, ?_, ?_⟩
    · rw [hcur]; dsimp only; rw [hmtop, hF.dir]; exact hmside
    · intro k hk₁ hk₂
      rw [hcur] at hk₁ ⊢
      rw [hmtop] at hk₁
      have hk : k = d + 1 := by omega
      subst hk
      dsimp only
      rw [hg, hsv, hmv]
      rcases hC.mid with h | h | h
      · exact .inl h
      · exact .inr (.inl ((Graph.Interior.congr fun e he => (hmE e he).symm).1 h))
      · exact .inr (.inr fun h' => h ((Graph.Touches.congr fun e he => hmE e he).1 h'))
  have hclose : CloseTwoOk (d + 1) s₂ := ⟨hmerge, hfin⟩
  have hstep := Step.closeTwo (v := v₀) hF.inv hF.shape hv item hclose hfree
  refine ⟨hclose, ?_⟩
  have hts₃ : (mergeTstackTops.run s₂).2.tstack = TEntry.mergeInto c t₂ :: R := by rw [hrun₁]
  have hrun₂ := finishTstackTop_run_eq (mergeTstackTops.run s₂).2 item (TEntry.mergeInto c t₂) R hts₃
  have hsd₃ : (mergeTstackTops.run s₂).2.stackDir[(TEntry.mergeInto c t₂).topDepth]! = dir := by
    rw [hrun₁]; dsimp only; rw [hmtop]; exact hF.dir
  rw [hsd₃] at hrun₂
  have hg₃ : (mergeTstackTops.run s₂).2.g = s.g := by rw [hrun₁]; exact hg
  have hsv₃ : (mergeTstackTops.run s₂).2.stackVerts = s.stackVerts := by rw [hrun₁]; exact hsv
  have hsd₃' : (mergeTstackTops.run s₂).2.stackDir = s₂.stackDir := by rw [hrun₁]
  have hit₃ : (mergeTstackTops.run s₂).2.items = s₂.items := by rw [hrun₁]
  set s' := ((finishTstackTop item).run (mergeTstackTops.run s₂).2).2 with hs'
  have hmem₂ : ∀ u ∈ R, u ∈ s₂.tstack := fun u hu => by
    rw [hts]; exact List.mem_cons_of_mem _ (List.mem_cons_of_mem _ hu)
  have hitem_not : ∀ u ∈ R, ∀ i ∈ u.spans.1 ++ u.spans.2, i ≠ item := fun u hu i hi hie => by
    subst hie; exact hfree.free u (hmem₂ u hu) hi
  have hnotm : item ∉ (TEntry.mergeInto c t₂).spans.1 ++ (TEntry.mergeInto c t₂).spans.2 := by
    rw [TEntry.mem_mergeInto_spans]
    rintro (h | h)
    · exact hfree.free c (by rw [hts]; exact List.mem_cons_self) h
    · exact hfree.free t₂ (by rw [hts]; exact List.mem_cons_of_mem _ List.mem_cons_self) h
  have hedge : ∀ e, e < s.g.ne → edgeItem s.g e ≠ item := fun e he => by
    have := hfree.node; rw [hg] at this; show 1 + s.g.nv + e ≠ item; omega
  have hmspan : ∀ i ∈ getSide (TEntry.mergeInto c t₂).spans dir, ∀ u ∈ R,
      i ∉ u.spans.1 ++ u.spans.2 := by
    intro i hi u hu
    have hi' := (mem_of_getSide_nil dir _ hmside i).2 hi
    rw [TEntry.mem_mergeInto_spans] at hi'
    rcases hi' with h | h
    · exact hP.fresh i h u (List.mem_cons_of_mem _ hu)
    · exact hfr i h u hu
  have hsT : s'.tstack = ⟨(TEntry.mergeInto c t₂).vStart, (TEntry.mergeInto c t₂).topDepth,
      (TEntry.mergeInto c t₂).firstIdx, setSides dir [item] []⟩ :: R := by rw [hs', hrun₂]
  have hsG : s'.g = s.g := by rw [hs', hrun₂]; exact hg₃
  have hsSV : s'.stackVerts = s.stackVerts := by rw [hs', hrun₂]; exact hsv₃
  have hsSD : s'.stackDir = s₂.stackDir := by rw [hs', hrun₂]; exact hsd₃'
  obtain ⟨vsv, hsI⟩ : ∃ vsv, s'.items = s₂.items.modify item fun it =>
      { it with vs := vsv, ch := getSide (TEntry.mergeInto c t₂).spans dir } :=
    ⟨_, by rw [hs', hrun₂, hit₃]⟩
  refine ⟨_, ⟨hsT, hmtop, hmv, getSide_setSides_not dir [item], ?_, ?_, ?_⟩, ⟨?_, ?_, ?_, ?_, ?_⟩,
    ⟨hsG, hsSV, by rw [hsSD]; exact hF.dir, hstep.inv, hstep.shape⟩⟩
  · intro i hi
    rw [mem_setSides, List.mem_singleton] at hi
    subst hi; exact .inl hiroot
  · intro i hi u hu
    rw [mem_setSides, List.mem_singleton] at hi
    subst hi; exact hfree.free u (hmem₂ u hu)
  · intro e he
    rw [hsI, TEntry.edges_finish dir (TEntry.mergeInto c t₂) item _ _ _ _ hfree.lt hfree.root hnotm
      hedge hmside he]
    exact hmE e he
  · rw [hsI, Array.size_modify]; exact hK.size
  · intro j hj
    rw [hsI, Items.type_modify_type_eq (f := fun it =>
      { it with vs := vsv, ch := getSide (TEntry.mergeInto c t₂).spans dir }) _ (fun _ => rfl)]
    exact hK.type j hj
  · intro u hu i hi j
    rw [hsI, Items.Below_modify_of_not_below _ _
      (Items.not_below_of_root hfree.root (hitem_not u hu i hi))]
    exact hK.below u hu i hi j
  · intro u hu i hi
    rw [hsI, Items.ch_modify_of_ne _ _ (hitem_not u hu i hi)]; exact hK.ch u hu i hi
  · intro u hu i hi p hp
    rw [hsI] at hp
    rcases Items.IsParent_modify hp with h | ⟨_, _, h⟩
    · exact hK.parent u hu i hi p h
    · exact absurd hi (hmspan i h u hu)

/-- The state after `maybeUnwrapNxt ty`, `mergeTstackTops`, `finishTstackTop` (the tail of
`loop1Body`). -/
def closeAt (ty : NodeType) (st : WalkState) : WalkState :=
  ((finishTstackTop (result (maybeUnwrapNxt ty) st)).run
    (mergeTstackTops.run (after (maybeUnwrapNxt ty) st)).2).2

theorem Items.IsParent_push {items : Items} {x : Item} (hx : x.ch = []) {p i : ItemId}
    (h : Items.IsParent (items.push x) p i) : Items.IsParent items p i := by
  unfold Items.IsParent at h ⊢; rwa [Items.ch_push_nil _ hx] at h

theorem Bool.eq_or_eq_not' (a b : Bool) : a = b ∨ a = !b := by cases a <;> cases b <;> simp

/-- `maybeUnwrapNxt ty` on the loop-1 stack `c :: t :: R`, then close `t` into `c`. -/
theorem l1_close {dir : Bool} {v v₀ : Nat} {E : Nat → Prop} {c t : TEntry} {R : List TEntry}
    {ty : NodeType}
    (hD : D = d + 1) (hty : ty ∉ [NodeType.F, .V, .Q]) (hs0 : Shape s)
    (hF : L1Frame s st D d dir) (hP : L1Piece s st d dir v E c (t :: R)) (hK : L1Keep s st (t :: R))
    (hts : t ∈ s.tstack) (hR : ∀ u ∈ R, u ∈ s.tstack)
    (hside : getSide t.spans (!dir) = []) (htop : d ≤ t.topDepth)
    (hC : L1Close d s v E t) (hU : ty ≠ .R → L1Unwrap s ty t)
    (hsd : ∀ i ∈ t.spans.1 ++ t.spans.2, ∀ u ∈ R, i ∉ u.spans.1 ++ u.spans.2)
    (hdisj : ∀ u ∈ R, ∀ e, e < s.g.ne → u.edges s.g s.items e → ¬ (E e ∨ t.edges s.g s.items e))
    (hv : v₀ < st.g.nv) :
    UnwrapOk ty st ∧ CloseTwoOk D (after (maybeUnwrapNxt ty) st) ∧
    ∃ c', L1Piece s (closeAt ty st) d dir t.vStart (fun e => E e ∨ t.edges s.g s.items e) c' R ∧
      L1Keep s (closeAt ty st) R ∧ L1Frame s (closeAt ty st) D d dir := by
  have hts' := hP.tstack
  have hg := hF.g
  have hnxt : nxtE st = t := by rw [nxtE, hts']; rfl
  have hcur : curE st = c := by rw [curE, hts']; rfl
  have htail : st.tstack.tail.tail = R := by rw [hts']; rfl
  have hnd : nxtDir st = st.stackDir[t.topDepth]! := by rw [nxtDir, hnxt]
  have hh : nxtHead st = (getSide t.spans st.stackDir[t.topDepth]!).head! := by
    rw [nxtHead, hnd, hnxt]
  have hrun := maybeUnwrapNxt_run_eq ty st c t R hts' _ rfl _ rfl
  rw [← hh] at hrun
  have hmemt : t ∈ st.tstack := by rw [hts']; exact List.mem_cons_of_mem _ List.mem_cons_self
  have htyF : ty ≠ .F := fun h => hty (h ▸ List.mem_cons_self)
  have htyQ : ty ≠ .Q := fun h => hty (h ▸ by simp)
  have hKR : L1Keep s st R := hK.mono fun u hu => List.mem_cons_of_mem _ hu
  have hKt := hK.edges List.mem_cons_self
  -- generic alloc branch
  have alloc : (maybeUnwrapNxt ty).run st = (allocItem ty).run st →
      CloseTwoOk D (after (maybeUnwrapNxt ty) st) ∧
      ∃ c', L1Piece s (closeAt ty st) d dir t.vStart (fun e => E e ∨ t.edges s.g s.items e) c' R ∧
        L1Keep s (closeAt ty st) R ∧ L1Frame s (closeAt ty st) D d dir := by
    intro hrun'
    have r := Step.allocRes (v := v₀) hF.inv hF.shape ty
    rw [run_allocItem] at hrun' r
    unfold closeAt after result
    rw [hrun']
    dsimp only
    have hsz := hK.size
    have hB : ∀ i j, Items.Below (st.items.push ⟨ty, (none, none), []⟩) i j ↔ Items.Below st.items i j :=
      fun i j => Items.Below_push_nil _ rfl
    have hE : ∀ (u : TEntry) e, u.edges s.g (st.items.push ⟨ty, (none, none), []⟩) e ↔ u.edges s.g st.items e :=
      fun u e => TEntry.edges_congr (fun i _ e' => hB i (edgeItem s.g e')) e
    refine l1_close_core hD ⟨hg, hF.sv, hF.dir, r.step.inv, r.step.shape⟩
      ⟨hts', hP.top, hP.bot, hP.side, hP.root, hP.fresh, fun e he => by rw [hE]; exact hP.edges e he⟩
      ⟨by rw [Array.size_push]; exact Nat.le_trans hsz (Nat.le_succ _), fun j hj => by
          rw [Items.type_push_of_ne _ (Nat.ne_of_lt (Nat.lt_of_lt_of_le hj hsz))]; exact hK.type j hj,
        fun u hu i hi j => by rw [hB]; exact hKR.below u hu i hi j,
        fun u hu i hi => by rw [Items.ch_push_nil _ rfl]; exact hKR.ch u hu i hi,
        fun u hu i hi p hp => hKR.parent u hu i hi p (Items.IsParent_push rfl hp)⟩
      hC htop rfl hside (fun e he => by rw [hE]; exact hKt e) hsd hdisj r.free
      (fun p hp => Nat.lt_irrefl _ (Nat.lt_of_lt_of_le (hs0.ch_lt p _ hp) hsz)) hv
  by_cases h1 : ty = .R ∨ st.ternarize = true
  · rw [if_pos h1] at hrun
    exact ⟨⟨by rw [hts']; simp, fun h => absurd h1 h⟩, alloc hrun⟩
  rw [if_neg h1] at hrun
  -- the head item
  have hh0 : nxtHead st < st.items.size ∧ nxtHead st < s.items.size ∧
      (getSide t.spans st.stackDir[t.topDepth]! ≠ [] → nxtHead st ∈ t.spans.1 ++ t.spans.2) := by
    by_cases hemp : getSide t.spans st.stackDir[t.topDepth]! = []
    · have h0 : ([] : List ItemId).head! = 0 := rfl
      rw [hh, hemp, h0]
      exact ⟨Nat.lt_of_lt_of_le (by omega) hF.shape.size,
        Nat.lt_of_lt_of_le (by omega) hs0.size, fun h => absurd rfl h⟩
    · have hmem : nxtHead st ∈ t.spans.1 ++ t.spans.2 := by
        rw [hh]
        exact (mem_of_getSide_nil dir t.spans hside _).2
          (by
            rcases Bool.eq_or_eq_not' st.stackDir[t.topDepth]! dir with h | h
            · rw [← h]; exact List.head!_mem_self hemp
            · rw [h, hside] at hemp; exact absurd rfl hemp)
      exact ⟨hF.shape.span t hmemt _ hmem, hs0.span t hts _ hmem, fun _ => hmem⟩
  have htype : st.items[nxtHead st]!.type = Items.type s.items (nxtHead st) := by
    rw [getElem!_pos st.items _ hh0.1, ← Items.type_eq_getElem hh0.1]; exact hK.type _ hh0.2.1
  by_cases h2 : st.items[nxtHead st]!.type = ty
  · rw [if_pos h2] at hrun
    rw [htype] at h2
    obtain ⟨hsingle, hroot, hchild⟩ := hU (fun h => h1 (.inl h)) _ _ hh.symm h2
    have hnd' : st.stackDir[t.topDepth]! = dir := by
      rcases Bool.eq_or_eq_not' st.stackDir[t.topDepth]! dir with h | h
      · exact h
      · rw [h, hside] at hsingle; cases hsingle
    rw [hnd'] at hsingle hrun
    have hmem : nxtHead st ∈ t.spans.1 ++ t.spans.2 :=
      (mem_of_getSide_nil dir t.spans hside _).2 (by rw [hsingle]; exact List.mem_singleton_self _)
    have hUok : UnwrapOk ty st :=
      ⟨by rw [hts']; simp, fun _ _ => ⟨by rw [hnxt, hnd, hnd']; exact hside, by rw [hnxt, hnd, hnd']; exact hsingle,
        fun p hp => hroot p (hK.parent t List.mem_cons_self _ hmem p hp),
        fun u hu => by
          rw [hcur, htail] at hu
          rcases List.mem_cons.1 hu with rfl | hu
          · exact fun hc => hP.fresh _ hc t List.mem_cons_self hmem
          · exact hsd _ hmem u hu⟩⟩
    refine ⟨hUok, ?_⟩
    have r := maybeUnwrapNxt_spec (v := v₀) hF.inv hF.shape hty hUok
    rw [hrun] at r
    have hchS : st.items[nxtHead st]!.ch = Items.ch s.items (nxtHead st) := by
      rw [getElem!_pos st.items _ hh0.1, ← Items.ch_eq_getElem hh0.1]; exact hK.ch t List.mem_cons_self _ hmem
    rw [hchS] at hrun r
    unfold closeAt after result
    rw [hrun]
    dsimp only
    have hedge : ∀ e, e < s.g.ne → edgeItem s.g e ≠ nxtHead st := fun e he heq => by
      have := hs0.edge e he; rw [heq, h2] at this; exact htyQ this
    have hT : ∀ e, e < s.g.ne →
        (TEntry.edges s.g st.items ⟨t.vStart, t.topDepth, t.firstIdx, setSides dir (Items.ch s.items (nxtHead st)) []⟩ e ↔
          t.edges s.g s.items e) := fun e he => by
      rw [← hK.ch t List.mem_cons_self _ hmem, TEntry.edges_unwrap dir _ _ _ _ hedge he, ← hKt e,
        TEntry.edges_single dir _ hside hsingle]
    have hfr : ∀ i ∈ (setSides dir (Items.ch s.items (nxtHead st)) []).1 ++
        (setSides dir (Items.ch s.items (nxtHead st)) []).2,
        ∀ u ∈ s.tstack, i ∉ u.spans.1 ++ u.spans.2 := by
      intro i hi
      rw [mem_setSides] at hi
      exact hchild i hi
    refine l1_close_core hD ⟨hg, hF.sv, hF.dir, r.step.inv, r.step.shape⟩
      ⟨rfl, hP.top, hP.bot, hP.side, hP.root, ?_, hP.edges⟩
      (hKR.items_eq rfl) hC htop rfl (getSide_setSides_not dir _) hT (fun i hi u hu => hfr i hi u (hR u hu)) hdisj
      r.free hroot hv
    intro i hi u hu
    rcases List.mem_cons.1 hu with rfl | hu
    · intro hi'
      rw [mem_setSides] at hi'
      rcases hP.root i hi with hr | ⟨u', hu', hiu⟩
      · exact hr _ hi'
      · exact hchild i hi' u' hu' hiu
    · exact hP.fresh i hi u (List.mem_cons_of_mem _ hu)
  · rw [if_neg h2] at hrun
    exact ⟨⟨by rw [hts']; simp, fun _ h => absurd h h2⟩, alloc hrun⟩

theorem L1Reach.append {hi done rest : List TEntry} (h : L1Reach d hi done rest) : hi = done ++ rest := by
  induction h with
  | nil => rfl
  | close _ _ ih => rw [ih]; simp
  | series _ _ ih => rw [ih]; simp

theorem l1Bot_concat₂ (o : DfsOut) (done : List TEntry) (t t' : TEntry) :
    l1Bot o (done ++ [t, t']) = t'.vStart := by simp [l1Bot]

theorem l1Edges_concat₂ (o : DfsOut) (s : WalkState) (done : List TEntry) (t t' : TEntry) (e : Nat) :
    l1Edges o s (done ++ [t, t']) e ↔ l1Edges o s (done ++ [t]) e ∨ t'.edges s.g s.items e := by
  rw [show done ++ [t, t'] = done ++ [t] ++ [t'] by simp]; exact l1Edges_concat o s _ t' e

theorem result_loop1Cond_iff (d : Nat) (st : WalkState) :
    result (loop1Cond d) st = true ↔ 2 ≤ st.tstack.length ∧ d ≤ st.tstack.tail.head!.topDepth := by
  unfold result; rw [run_loop1Cond]; simp

theorem loop1Body_run_eq (d : Nat) (dir : Bool) (st : WalkState) :
    ((loop1Body d dir).run st).2 = closeAt (l1Ty d dir st) (l1S₁ d dir st) := rfl

/-- The facts about the state `s` before `finishEdge` that loop 1 consumes: the stack splits as
`hi ++ lo ++ base` (`hi` the range), `lo` starts below `d`, and the `EarFinish` ownership facts. -/
structure L1Ctx (D d : Nat) (o : DfsOut) (s : WalkState) (hi lo base : List TEntry) : Prop where
  D_eq : D = d + 1
  shape : Shape s
  tstack : s.tstack = hi ++ lo ++ base
  lo_top : ∀ t ∈ lo.head?, t.topDepth < d
  lo_ne : lo ≠ []
  hi_top : ∀ t ∈ hi, d ≤ t.topDepth
  side : ∀ t ∈ hi, getSide t.spans (!s.stackDir[d]!) = []
  spec : Loop1Spec d o s hi
  q_root : ∀ p, ¬ Items.IsParent s.items p (edgeItem s.g o.e)
  q_free : ∀ t ∈ s.tstack, edgeItem s.g o.e ∉ t.spans.1 ++ t.spans.2
  disj : s.tstack.Pairwise fun t t' => ∀ e, e < s.g.ne → t.edges s.g s.items e → ¬ t'.edges s.g s.items e
  span_disj : s.tstack.Pairwise fun t t' => ∀ i ∈ t.spans.1 ++ t.spans.2, i ∉ t'.spans.1 ++ t'.spans.2

theorem List.pairwise_mid {α : Type} {R : α → α → Prop} {l l₁ l₂ : List α} {a : α} (h : l.Pairwise R)
    (hl : l = l₁ ++ a :: l₂) : ∀ b ∈ l₂, R a b := by
  subst hl
  exact (List.pairwise_cons.1 (List.pairwise_append.1 h).2.1).1

/-- Entries under the next entry `t` are edge-disjoint from the piece and from `t`, and span-disjoint
from `t`. -/
theorem L1Ctx.disj_split {hi lo base : List TEntry} (hc : L1Ctx D d o s hi lo base)
    {done : List TEntry} {t : TEntry} {R : List TEntry} (hsplit : s.tstack = done ++ t :: R) :
    (∀ u ∈ R, ∀ e, e < s.g.ne → u.edges s.g s.items e → ¬ (l1Edges o s done e ∨ t.edges s.g s.items e)) ∧
    (∀ i ∈ t.spans.1 ++ t.spans.2, ∀ u ∈ R, i ∉ u.spans.1 ++ u.spans.2) := by
  refine ⟨fun u hu e he hue h => ?_, fun i hi u hu => List.pairwise_mid hc.span_disj hsplit u hu i hi⟩
  rcases h with (rfl | ⟨t', ht', hte⟩) | h
  · obtain ⟨i, hi, hb⟩ := hue
    have hb' : Items.Below s.items i (edgeItem s.g o.e) := hb
    rw [Items.Below.eq_of_no_parent hc.q_root hb'] at hi
    exact hc.q_free u (by rw [hsplit]; exact List.mem_append_right _ (List.mem_cons_of_mem _ hu)) hi
  · exact List.pairwise_split hc.disj hsplit t' ht' u (List.mem_cons_of_mem _ hu) e he hte hue
  · exact List.pairwise_mid hc.disj hsplit u hu e he h hue

/-- The loop-1 invariant at an iteration boundary: some reached split `done ++ rest` of the range,
the piece on top with bottom `l1Bot o done` and edges `l1Edges o s done`, and the rest kept. -/
def L1Inv (D d : Nat) (o : DfsOut) (s : WalkState) (hi lo base : List TEntry) (st : WalkState) : Prop :=
  ∃ done rest c, L1Reach d hi done rest ∧
    L1Piece s st d s.stackDir[d]! (l1Bot o done) (l1Edges o s done) c (rest ++ lo ++ base) ∧
    L1Keep s st (rest ++ lo ++ base) ∧ L1Frame s st D d s.stackDir[d]!

theorem l1_step {hi lo base : List TEntry} {v₀ : Nat} (hc : L1Ctx D d o s hi lo base)
    (hI : L1Inv D d o s hi lo base st) (hcond : result (loop1Cond d) st = true) (hv : v₀ < st.g.nv) :
    Loop1BodyOk D d s.stackDir[d]! st ∧
      L1Inv D d o s hi lo base ((loop1Body d s.stackDir[d]!).run st).2 := by
  obtain ⟨done, rest, c, hreach, hP, hK, hF⟩ := hI
  have hts := hP.tstack
  rw [result_loop1Cond_iff] at hcond
  obtain ⟨t, rest', rfl⟩ : ∃ t rest', rest = t :: rest' := by
    cases rest with
    | cons t rest' => exact ⟨t, rest', rfl⟩
    | nil =>
      exfalso
      obtain ⟨l, lo', rfl⟩ : ∃ l lo', lo = l :: lo' := by
        cases lo with
        | nil => exact absurd rfl hc.lo_ne
        | cons l lo' => exact ⟨l, lo', rfl⟩
      have h1 := hc.lo_top l rfl
      have h2 := hcond.2
      rw [hts] at h2
      simp only [List.nil_append, List.cons_append, List.tail_cons, List.head!_cons] at h2
      omega
  have hsplit : s.tstack = done ++ t :: (rest' ++ lo ++ base) := by
    rw [hc.tstack, hreach.append]; simp
  obtain ⟨hdisj, hsd⟩ := hc.disj_split hsplit
  have hthi : t ∈ hi := by rw [hreach.append]; simp
  have htmem : t ∈ s.tstack := by rw [hsplit]; simp
  have hR : ∀ u ∈ rest' ++ lo ++ base, u ∈ s.tstack := fun u hu => by
    rw [hsplit]; exact List.mem_append_right _ (List.mem_cons_of_mem _ hu)
  have hside := hc.side t hthi
  have htop := hc.hi_top t hthi
  have hspec := hc.spec done t rest' hreach
  have hnxt : nxtE st = t := by rw [nxtE, hts]; rfl
  have hcur : curE st = c := by rw [curE, hts]; rfl
  have hrunB := loop1Body_run_eq d s.stackDir[d]! st
  have hres := loop1Type_result d s.stackDir[d]! st
  by_cases hgt : d < t.topDepth
  · obtain ⟨hM, t', rest'', rfl, hC', hU'⟩ := hspec.2 hgt
    have hrun := loop1Type_run d s.stackDir[d]! st
    rw [hnxt, if_pos hgt] at hrun
    have hts' : ({ st with stackDir := st.stackDir.set! t.topDepth s.stackDir[d]! } : WalkState).tstack =
        c :: t :: (t' :: rest'' ++ lo ++ base) := hts
    rw [mergeTstackTops_run_eq _ c t _ hts'] at hrun
    have hty : l1Ty d s.stackDir[d]! st = .S := by unfold l1Ty result; rw [hrun]
    have hS₁ : l1S₁ d s.stackDir[d]! st =
        { st with stackDir := st.stackDir.set! t.topDepth s.stackDir[d]!,
                  tstack := TEntry.mergeInto c t :: (t' :: rest'' ++ lo ++ base) } := by
      unfold l1S₁ after; rw [hrun]
    have hmergeS : MergeTopOk D { st with stackDir := st.stackDir.set! t.topDepth s.stackDir[d]! } := by
      intro cur nxt R hr
      rw [hts'] at hr
      simp only [List.cons.injEq] at hr
      obtain ⟨rfl, rfl, rfl⟩ := hr
      have hg := hF.g
      have hKt := hK.edges (List.mem_cons_self (a := t)) 
      refine ⟨⟨fun _ hne => ?_, ?_⟩, fun u hu e he hue => ?_⟩
      · dsimp only at hne ⊢
        rw [hg] at hne ⊢
        obtain ⟨e, he, het⟩ := hne
        obtain ⟨x, hx₁, hx₂⟩ := hM.share ⟨e, he, (hKt e).1 het⟩
        exact ⟨x, (Graph.Touches.congr fun e he => (hP.edges e he).symm).1 hx₁,
          (Graph.Touches.congr fun e he => (hKt e).symm).1 hx₂⟩
      · have hmtop : (TEntry.mergeInto c t).topDepth = d := by
          show min t.topDepth c.topDepth = d; rw [hP.top]; omega
        unfold TEntry.Term
        dsimp only
        rw [hg, hF.sv]
        rcases hM.bottom with h | ⟨k, hk₁, hk₂, hk⟩ | h
        · left; left; exact hP.bot.trans h
        · left; right; exact ⟨k, by rw [hmtop]; exact hk₁, by rw [hc.D_eq]; exact hk₂, by rw [hP.bot]; exact hk⟩
        · right
          rw [hP.bot]
          exact (Graph.Interior.congr fun e he => by rw [hP.edges e he, hKt e]).2 h
      · dsimp only at he hue ⊢
        rw [hg] at he hue
        have := hdisj u hu e he ((hK.edges (List.mem_cons_of_mem _ hu) e).1 hue)
        rintro (h | h) <;> rw [hg] at h
        · exact this (.inl ((hP.edges e he).1 h))
        · exact this (.inr ((hKt e).1 h))
    have st₁ := Step.loop1Type (v := v₀) hF.inv hF.shape (d := d) (edgeDir := s.stackDir[d]!)
      (fun _ => by rw [hnxt]; exact hmergeS)
    have hinv : (l1S₁ d s.stackDir[d]! st).Inv' D := st₁.inv
    have hshape : Shape (l1S₁ d s.stackDir[d]! st) := st₁.shape
    rw [hS₁] at hinv hshape
    have hF₁ : L1Frame s (l1S₁ d s.stackDir[d]! st) D d s.stackDir[d]! := by
      rw [hS₁]
      exact ⟨hF.g, hF.sv, by
        show (st.stackDir.set! t.topDepth s.stackDir[d]!)[d]! = _
        rw [Array.getElem!_set!_ne _ _ _ _ (by omega)]; exact hF.dir, hinv, hshape⟩
    have hsplit₁ : s.tstack = (done ++ [t]) ++ t' :: (rest'' ++ lo ++ base) := by
      rw [hsplit]; simp
    obtain ⟨hdisj₁, hsd₁⟩ := hc.disj_split hsplit₁
    have ht'hi : t' ∈ hi := by rw [hreach.append]; simp
    have hP₁ : L1Piece s (l1S₁ d s.stackDir[d]! st) d s.stackDir[d]! t.vStart (l1Edges o s (done ++ [t]))
        (TEntry.mergeInto c t) (t' :: rest'' ++ lo ++ base) := by
      rw [hS₁]
      refine ⟨rfl, ?_, rfl, TEntry.getSide_mergeInto _ _ _ hP.side hside, ?_, ?_, ?_⟩
      · show min t.topDepth c.topDepth = d; rw [hP.top]; omega
      · intro i hi
        rw [TEntry.mem_mergeInto_spans] at hi
        rcases hi with h | h
        · exact hP.root i h
        · exact .inr ⟨t, htmem, h⟩
      · intro i hi u hu
        rw [TEntry.mem_mergeInto_spans] at hi
        rcases hi with h | h
        · exact hP.fresh i h u (List.mem_cons_of_mem _ hu)
        · exact hsd i h u hu
      · intro e he
        show TEntry.edges s.g st.items _ e ↔ _
        rw [TEntry.edges_mergeInto, hP.edges e he, hK.edges List.mem_cons_self e, l1Edges_concat]
    have hK₁ : L1Keep s (l1S₁ d s.stackDir[d]! st) (t' :: rest'' ++ lo ++ base) := by
      rw [hS₁]; exact (hK.mono fun u hu => List.mem_cons_of_mem _ hu).items_eq rfl
    have hv₁ : v₀ < (l1S₁ d s.stackDir[d]! st).g.nv := by rw [hS₁]; exact hv
    obtain ⟨hUok, hclose, c', hP', hK', hF'⟩ := l1_close (ty := .S) hc.D_eq (by simp) hc.shape hF₁ hP₁ hK₁
      (hR t' (by simp)) (fun u hu => hR u (List.mem_cons_of_mem _ hu)) (hc.side t' ht'hi) (hc.hi_top t' ht'hi)
      hC' (fun _ => hU') hsd₁ hdisj₁ hv₁
    refine ⟨⟨fun _ => by rw [hnxt]; exact hmergeS, by rw [hty]; exact hUok,
      by show CloseTwoOk D (after (maybeUnwrapNxt (l1Ty d _ st)) (l1S₁ d _ st)); rw [hty]; exact hclose⟩, ?_⟩
    rw [hrunB, hty]
    refine ⟨done ++ [t, t'], rest'', c', L1Reach.series hreach hgt, ?_, hK', hF'⟩
    rw [l1Bot_concat₂]
    exact { hP' with edges := fun e he => (hP'.edges e he).trans (l1Edges_concat₂ o s done t t' e).symm }
  · have htd : t.topDepth = d := by omega
    obtain ⟨hC, hUP⟩ := hspec.1 htd
    have hrun := loop1Type_run d s.stackDir[d]! st
    rw [hnxt, hcur, if_neg hgt] at hrun
    have hS₁ : l1S₁ d s.stackDir[d]! st = st := by
      unfold l1S₁ after; rw [hrun]; split <;> rfl
    have hty : l1Ty d s.stackDir[d]! st = if t.vStart = c.vStart then NodeType.P else .R := by
      unfold l1Ty result; rw [hrun]
      by_cases hp : t.vStart = c.vStart
      · rw [if_pos (beq_iff_eq.2 hp), if_pos hp]
      · rw [if_neg (by simpa using hp), if_neg hp]
    have hU : l1Ty d s.stackDir[d]! st ≠ .R → L1Unwrap s (l1Ty d s.stackDir[d]! st) t := by
      rw [hty]
      split
      · rename_i hp; exact fun _ => hUP (hp.trans hP.bot)
      · exact fun h => absurd rfl h
    obtain ⟨hUok, hclose, c', hP', hK', hF'⟩ := l1_close (v₀ := v₀) hc.D_eq hres hc.shape hF hP hK htmem hR
      hside htop hC hU hsd hdisj hv
    refine ⟨⟨fun h => absurd h (by rw [hnxt]; exact hgt), by rw [hS₁]; exact hUok,
      by show CloseTwoOk D (after (maybeUnwrapNxt (l1Ty d _ st)) (l1S₁ d _ st)); rw [hS₁]; exact hclose⟩, ?_⟩
    rw [hrunB, hS₁]
    refine ⟨done ++ [t], rest', c', L1Reach.close hreach htd, ?_, hK', hF'⟩
    rw [l1Bot_concat]
    exact { hP' with edges := fun e he => (hP'.edges e he).trans (l1Edges_concat o s done t e).symm }

theorem l1_iter {hi lo base : List TEntry} {v₀ : Nat} (hc : L1Ctx D d o s hi lo base) {st₀ : WalkState}
    (h0 : L1Inv D d o s hi lo base st₀) (hv : v₀ < s.g.nv) (k : Nat)
    (hk : ∀ j, j < k → result (loop1Cond d) (iter (loop1Body d s.stackDir[d]!) j st₀) = true) :
    L1Inv D d o s hi lo base (iter (loop1Body d s.stackDir[d]!) k st₀) := by
  induction k with
  | zero => exact h0
  | succ k ih =>
    rw [iter_succ']
    have hI := ih fun j hj => hk j (Nat.lt_succ_of_lt hj)
    obtain ⟨_, _, _, _, _, _, hF⟩ := id hI
    exact (l1_step hc hI (hk k (Nat.lt_succ_self k)) (by rw [hF.g]; exact hv)).2

theorem loop1_ok {hi lo base : List TEntry} {v₀ : Nat} (hc : L1Ctx D d o s hi lo base) {st₀ : WalkState}
    (h0 : L1Inv D d o s hi lo base st₀) (hv : v₀ < s.g.nv) (k : Nat)
    (hk : ∀ j, j ≤ k → result (loop1Cond d) (iter (loop1Body d s.stackDir[d]!) j st₀) = true) :
    Loop1BodyOk D d s.stackDir[d]! (iter (loop1Body d s.stackDir[d]!) k st₀) := by
  have hI := l1_iter hc h0 hv k fun j hj => hk j (Nat.le_of_lt hj)
  obtain ⟨_, _, _, _, _, _, hF⟩ := id hI
  exact (l1_step hc hI (hk k (Nat.le_refl k)) (by rw [hF.g]; exact hv)).1

theorem List.mem_takeWhile' {α : Type} {p : α → Bool} : ∀ {l : List α} {a : α}, a ∈ l.takeWhile p → p a = true
  | [], _, h => by simp at h
  | b :: l, a, h => by
    rw [List.takeWhile_cons] at h
    split at h
    · rcases List.mem_cons.1 h with rfl | h
      · assumption
      · exact List.mem_takeWhile' h
    · simp at h

/-- The loop-1 context of a tree edge with `lowval < d`: the range is the longest prefix of `sub` at
depth `≥ d`; it is proper since the `(y, lowval)` piece lies below it. -/
theorem L1Ctx.ofEar {curV : Nat} {hasVert : Bool} {sub base : List TEntry}
    (hE : s.EarFinish curV d o hasVert sub base) (hs : Shape s) (hD : D = d + 1)
    (ht : o.cls.isTree = true) (hlow : o.cls.lowval d < d) :
    ∃ hi lo, sub = hi ++ lo ∧ L1Ctx D d o s hi lo base := by
  refine ⟨sub.takeWhile (fun t => decide (d ≤ t.topDepth)), sub.dropWhile (fun t => decide (d ≤ t.topDepth)),
    List.takeWhile_append_dropWhile.symm, ?_⟩
  have hhi : ∀ t ∈ sub.takeWhile (fun t => decide (d ≤ t.topDepth)), d ≤ t.topDepth := fun t ht => by
    simpa using List.mem_takeWhile' ht
  have hlo : ∀ t ∈ (sub.dropWhile fun t => decide (d ≤ t.topDepth)).head?, t.topDepth < d := by
    intro t ht
    have := List.head?_dropWhile_not (fun t : TEntry => decide (d ≤ t.topDepth)) sub
    rw [Option.mem_def] at ht
    rw [ht] at this
    simpa using this
  have hrange : Loop1Range d sub (sub.takeWhile fun t => decide (d ≤ t.topDepth)) :=
    ⟨_, List.takeWhile_append_dropWhile.symm, hhi, hlo⟩
  refine ⟨hD, hs, by rw [hE.tstack, List.takeWhile_append_dropWhile], hlo, ?_, hhi,
    hE.loop1_side hlow _ hrange, hE.loop1 ht hlow _ hrange, hE.q_root, hE.q_free, hE.disj, hE.span_disj⟩
  intro h
  obtain ⟨mid, py, vy, hsub, hB⟩ := hE.bottom ht hlow
  have hsub' : sub.takeWhile (fun t => decide (d ≤ t.topDepth)) = sub := by
    have := @List.takeWhile_append_dropWhile _ (fun t => decide (d ≤ t.topDepth)) sub
    rwa [h, List.append_nil] at this
  have hpy : py ∈ sub.takeWhile fun t => decide (d ≤ t.topDepth) := by rw [hsub', hsub]; simp
  have := hhi py hpy
  rw [hB.py_top] at this
  omega

/-- The loop-1 invariant at the start of loop 1: the tree-edge entry on top of the untouched stack. -/
theorem l1_init {curV : Nat} {hasVert : Bool} {sub base : List TEntry}
    (hE : s.EarFinish curV d o hasVert sub base) (hi : s.Inv' D) (hs : Shape s)
    (hD : D = d + 1) (he : o.e < s.g.ne) (hq : Items.ch s.items (edgeItem s.g o.e) = [])
    (hends : Items.PairEq (o.dest, s.stackVerts[d]!) s.g.edges[o.e]!)
    {hi' lo : List TEntry} (hrange : sub = hi' ++ lo) :
    L1Inv D d o s hi' lo base (ceS₁ o.dest d o.e (feS₀ d o s)) := by
  have hmem : ∀ u ∈ hi' ++ lo ++ base, u ∈ s.tstack := fun u hu => by rw [hE.tstack, hrange]; exact hu
  have hne : ∀ u ∈ s.tstack, ∀ i ∈ u.spans.1 ++ u.spans.2, i ≠ edgeItem s.g o.e :=
    fun u hu i hi h => hE.q_free u hu (h ▸ hi)
  set f : Item → Item := fun it => { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) }
    with hf
  have hst : ceS₁ o.dest d o.e (feS₀ d o s) =
      { s with items := s.items.modify (edgeItem s.g o.e) f,
               tstack := ⟨o.dest, d, s.nxtEdgeIdx, setSides s.stackDir[d]! [edgeItem s.g o.e] []⟩ :: s.tstack } := rfl
  have hq₀ : Items.ch (s.items.modify (edgeItem s.g o.e) f) (edgeItem s.g o.e) = [] := by
    rw [Items.ch_modify_ch_eq _ f (fun _ => rfl)]; exact hq
  have st₀ : Step D curV s (feS₀ d o s) :=
    Step.modifyVs hi hs (edgeItem s.g o.e) _ (by show 1 + s.g.nv + o.e < _; omega)
  have st₁ : Step D curV (feS₀ d o s) (ceS₁ o.dest d o.e (feS₀ d o s)) :=
    Step.pushEdge st₀.inv st₀.shape o.dest d o.e he hq₀ hends (by omega)
  refine ⟨[], hi', ⟨o.dest, d, s.nxtEdgeIdx, setSides s.stackDir[d]! [edgeItem s.g o.e] []⟩, L1Reach.nil _,
    ⟨?_, rfl, rfl, getSide_setSides_not _ _, ?_, ?_, ?_⟩, ⟨?_, ?_, ?_, ?_, ?_⟩,
    ⟨rfl, rfl, rfl, st₁.inv, st₁.shape⟩⟩
  · rw [hst]; dsimp only; rw [hE.tstack, hrange]
  · intro i hi
    rw [mem_setSides, List.mem_singleton] at hi
    subst hi; exact .inl hE.q_root
  · intro i hi u hu
    rw [mem_setSides, List.mem_singleton] at hi
    subst hi; exact hE.q_free u (hmem u hu)
  · intro e he
    rw [hst]; dsimp only
    rw [TEntry.edges_edgeEntry _ _ _ _ _ hq₀]
    simp [l1Edges]
  · rw [hst]; dsimp only; simp
  · intro j hj
    rw [hst]; dsimp only
    exact Items.type_modify_type_eq _ f (fun _ => rfl) j
  · intro u hu i hi j
    rw [hst]; dsimp only
    exact Items.Below_modify_of_not_below _ f (Items.not_below_of_root hE.q_root (hne u (hmem u hu) i hi))
  · intro u hu i hi
    rw [hst]; dsimp only
    exact Items.ch_modify_of_ne _ f (hne u (hmem u hu) i hi)
  · intro u hu i hi p hp
    rw [hst] at hp; dsimp only at hp
    rcases Items.IsParent_modify hp with h | ⟨_, hj, h⟩
    · exact h
    · simp only [hf] at h
      rw [← Items.ch_eq_getElem hj, hq] at h
      cases h

end WalkState
end Spqr
