import Spqr.RangesInv

/-!
# `RangesInv` across the blocks of `finishEdge`

Companion of `WalkSpec.Step`: the same block structure (`closeEars`, `mergeLate`, `closeVert'`,
`finishRest`), each block taking its `WalkSpec` `*Ok` bundle (stack shape, `Inv'`) plus a
`*Adj` bundle giving the local adjacency condition `MergeAdj` at every merge site and the edge
positions at the push sites. `MergeAdj` is not derivable from `RangesInv` (PROOF.md §4.6: the holes
between non-adjacent pieces are blocks under unpushed path vertex items), so it is a hypothesis
per site, like `MergeTopOk`.
-/

namespace Spqr
open WalkM

namespace WalkState
variable {s s' : WalkState} {σ : List Nat} {n D : Nat}

theorem RangesInv.frame (h : s.RangesInv σ n D) (hg : s'.g = s.g) (hsv : s'.stackVerts = s.stackVerts)
    (hitems : s'.items = s.items) (hts : s'.tstack = s.tstack) : s'.RangesInv σ n D := by
  refine ⟨h.inv.frame hg hsv hitems hts, ?_, ?_, ?_, ?_⟩ <;> simp only [hg, hitems, hts]
  exacts [h.processed, h.ordered, h.convex, h.closed]

theorem RangesInv.modifyVs (j : ItemId) (vsv : Option Nat × Option Nat) (h : s.RangesInv σ n D)
    (hj : j < 1 + s.g.nv + s.g.ne) :
    RangesInv σ n D { s with items := s.items.modify j fun it => { it with vs := vsv } } := by
  have hinv := h.inv.modifyVs j vsv hj
  set f : Item → Item := fun it => { it with vs := vsv } with hf
  have hch : ∀ p, Items.ch (s.items.modify j f) p = Items.ch s.items p :=
    Items.ch_modify_ch_eq j f fun _ => rfl
  have hty : ∀ p, Items.type (s.items.modify j f) p = Items.type s.items p :=
    Items.type_modify_type_eq j f fun _ => rfl
  have hB : ∀ a i, Items.Below (s.items.modify j f) a i ↔ Items.Below s.items a i :=
    fun _ _ => Items.Below_modify_ch_eq j f fun _ => rfl
  have hN : ∀ a i, Items.BelowNoV (s.items.modify j f) a i ↔ Items.BelowNoV s.items a i :=
    fun _ _ => Items.BelowNoV_congr hch (fun _ c _ => hty c)
  refine h.items_congr hinv (fun t _ e => TEntry.edges_congr (fun i _ e => hB i _) e)
    (fun t _ e => TEntry.piece_congr (fun i _ => hty i) (fun i _ e => hN i _) e) ?_
  intro i hi hty' a b c hab hbc hc hpa hpc
  simp only [Array.size_modify] at hi
  rw [hty] at hty'
  exact (hB i _).2 (h.closed i hi hty' a b c hab hbc hc ((hN i _).1 hpa) ((hN i _).1 hpc))

/-- Replacing the top entry by one with the same edges and piece. -/
theorem RangesInv.replaceCur (t t' : TEntry) (rest : List TEntry) (hs : s.tstack = t :: rest)
    (h : s.RangesInv σ n D) (hσ : ∀ e ∈ σ, e < s.g.ne) (hinv : Inv' D { s with tstack := t' :: rest })
    (hE : ∀ e, e < s.g.ne → (t'.edges s.g s.items e ↔ t.edges s.g s.items e))
    (hP : ∀ e, e < s.g.ne → (t'.piece s.g s.items e ↔ t.piece s.g s.items e)) :
    RangesInv σ n D { s with tstack := t' :: rest } := by
  refine ⟨hinv, ?_, ?_, ?_, h.closed⟩
  · intro u hu e he hue
    rcases List.mem_cons.1 hu with rfl | hu
    · exact h.processed t (by simp [hs]) e he ((hE e he).1 hue)
    · exact h.processed u (by simp [hs, hu]) e he hue
  · intro above u below hs' u' hu' e e' he he' hp hp'
    cases above with
    | nil =>
      simp only [List.nil_append, List.cons.injEq] at hs'; obtain ⟨rfl, rfl⟩ := hs'
      exact h.ordered [] t rest hs u' hu' e e' he he' ((hP e he).1 hp) hp'
    | cons a above =>
      simp only [List.cons_append, List.cons.injEq] at hs'; obtain ⟨rfl, hs'⟩ := hs'
      exact h.ordered (t :: above) u below (by rw [hs, hs']; rfl) u' hu' e e' he he' hp hp'
  · intro u hu a b c hab hbc hc hpa hpc
    have hσa := getElem!_lt hσ (s := s) (by omega : a < σ.length)
    have hσb := getElem!_lt hσ (s := s) (by omega : b < σ.length)
    have hσc := getElem!_lt hσ (s := s) hc
    rcases List.mem_cons.1 hu with rfl | hu
    · exact (hE _ hσb).2 (h.convex t (by simp [hs]) a b c hab hbc hc ((hP _ hσa).1 hpa) ((hP _ hσc).1 hpc))
    · exact h.convex u (by simp [hs, hu]) a b c hab hbc hc hpa hpc

/-- Replacing the entry below the top by one with the same edges and piece. -/
theorem RangesInv.replaceNxt (a b b' : TEntry) (rest : List TEntry) (hs : s.tstack = a :: b :: rest)
    (h : s.RangesInv σ n D) (hσ : ∀ e ∈ σ, e < s.g.ne) (hinv : Inv' D { s with tstack := a :: b' :: rest })
    (hE : ∀ e, e < s.g.ne → (b'.edges s.g s.items e ↔ b.edges s.g s.items e))
    (hP : ∀ e, e < s.g.ne → (b'.piece s.g s.items e ↔ b.piece s.g s.items e)) :
    RangesInv σ n D { s with tstack := a :: b' :: rest } := by
  refine ⟨hinv, ?_, ?_, ?_, h.closed⟩
  · intro u hu e he hue
    rcases List.mem_cons.1 hu with rfl | hu
    · exact h.processed u (by simp [hs]) e he hue
    rcases List.mem_cons.1 hu with rfl | hu
    · exact h.processed b (by simp [hs]) e he ((hE e he).1 hue)
    · exact h.processed u (by simp [hs, hu]) e he hue
  · intro above u below hs' u' hu' e e' he he' hp hp'
    match above, hs' with
    | [], hs' =>
      simp only [List.nil_append, List.cons.injEq] at hs'; obtain ⟨rfl, rfl⟩ := hs'
      rcases List.mem_cons.1 hu' with rfl | hu'
      · exact h.ordered [] a (b :: rest) hs b (by simp) e e' he he' hp ((hE e' he').1 hp')
      · exact h.ordered [] a (b :: rest) hs u' (by simp [hu']) e e' he he' hp hp'
    | [_], hs' =>
      simp only [List.cons_append, List.nil_append, List.cons.injEq] at hs'; obtain ⟨rfl, rfl, rfl⟩ := hs'
      exact h.ordered [a] b rest hs u' hu' e e' he he' ((hP e he).1 hp) hp'
    | _ :: _ :: above, hs' =>
      simp only [List.cons_append, List.cons.injEq] at hs'; obtain ⟨rfl, rfl, hs'⟩ := hs'
      exact h.ordered (a :: b :: above) u below (by rw [hs, hs']; rfl) u' hu' e e' he he' hp hp'
  · intro u hu a' b₁ c hab hbc hc hpa hpc
    have hσa := getElem!_lt hσ (s := s) (by omega : a' < σ.length)
    have hσb := getElem!_lt hσ (s := s) (by omega : b₁ < σ.length)
    have hσc := getElem!_lt hσ (s := s) hc
    rcases List.mem_cons.1 hu with rfl | hu
    · exact h.convex u (by simp [hs]) a' b₁ c hab hbc hc hpa hpc
    rcases List.mem_cons.1 hu with rfl | hu
    · exact (hE _ hσb).2 (h.convex b (by simp [hs]) a' b₁ c hab hbc hc ((hP _ hσa).1 hpa) ((hP _ hσc).1 hpc))
    · exact h.convex u (by simp [hs, hu]) a' b₁ c hab hbc hc hpa hpc

/-! ### The primitives under their `WalkSpec` bundles -/

/-- The local adjacency condition of `RangesInv.mergeTop` for the top two entries: every
`σ`-position between a piece edge of `nxt` and one of `cur` is an edge of one of them. -/
def MergeAdj (σ : List Nat) (s : WalkState) : Prop :=
  ∀ cur nxt rest, s.tstack = cur :: nxt :: rest → ∀ a b c, a ≤ b → b ≤ c → c < σ.length →
    nxt.piece s.g s.items σ[a]! → cur.piece s.g s.items σ[c]! →
      cur.edges s.g s.items σ[b]! ∨ nxt.edges s.g s.items σ[b]!

theorem RangesInv.mergeTop' (h : s.RangesInv σ n D) (hnd : σ.Nodup) (hσ : ∀ e ∈ σ, e < s.g.ne)
    (hok : MergeTopOk D s) (hadj : MergeAdj σ s) : RangesInv σ n D (mergeTstackTops.run s).2 := by
  match hts : s.tstack with
  | [] | [_] =>
    rw [run_mergeTstackTops, hts]
    exact ⟨h.inv.setTstack _ Stack_nil, by simp [mergeTops], by simp [mergeTops], by simp [mergeTops], h.closed⟩
  | cur :: nxt :: rest =>
    exact h.mergeTop cur nxt rest hts hnd hσ (hok cur nxt rest hts).1 (hok cur nxt rest hts).2
      (hadj cur nxt rest hts)

theorem RangesInv.finishTop' (h : s.RangesInv σ n D) (hσ : ∀ e ∈ σ, e < s.g.ne) (item : ItemId)
    (hok : FinishTopOk D s) (hf : ItemFree s item) : RangesInv σ n D ((finishTstackTop item).run s).2 := by
  match hts : s.tstack with
  | [] => exact absurd hts hok.nonempty
  | t :: rest =>
    have hc : curE s = t := by rw [curE, hts, List.head!_cons']
    have hside := hok.side
    have hmid := hok.mid
    rw [hc] at hside hmid
    exact h.finishTop item t rest hts hσ hf.lt hf.node hf.root hf.free hside hmid

theorem RangesInv.retarget' (h : s.RangesInv σ n D) (hσ : ∀ e ∈ σ, e < s.g.ne) (curV : Nat) (edgeDir : Bool)
    (hok : RetargetOk D curV s) : RangesInv σ n D ((retarget curV edgeDir).run s).2 := by
  match hts : s.tstack with
  | [] => exact absurd hts hok.nonempty
  | t :: rest =>
    have hc : curE s = t := by rw [curE, hts, List.head!_cons']
    have old := hok.old
    rw [hc] at old
    have hinv := Inv'.retarget curV edgeDir t rest hts h.inv old fun u hu e he hue => by
      have := hok.disj u (by rw [hts]; exact hu) e he hue
      rwa [hc] at this
    rw [retarget_run_eq curV edgeDir s t rest hts] at hinv ⊢
    refine h.replaceCur t _ rest hts hσ hinv (fun e _ => ?_) (fun e _ => ?_)
    · simp only [TEntry.edges, mem_setSides, List.append_nil]
    · simp only [TEntry.piece, mem_setSides, List.append_nil]

theorem RangesInv.alloc' (h : s.RangesInv σ n D) (hσ : ∀ e ∈ σ, e < s.g.ne) (hs : Shape s) (ty : NodeType) :
    RangesInv σ n D ((allocItem ty).run s).2 := by
  rw [run_allocItem]
  exact h.alloc ty hσ hs.size (fun p c hpc => hs.ch_lt p c hpc) (fun t ht i hi => hs.span t ht i hi)

theorem RangesInv.unwrap' {ty : NodeType} (h : s.RangesInv σ n D) (hσ : ∀ e ∈ σ, e < s.g.ne) (hs : Shape s)
    (hty : ty ∉ [NodeType.F, .V, .Q]) (hok : UnwrapOk ty s) :
    RangesInv σ n D ((maybeUnwrapNxt ty).run s).2 := by
  have hinv := (maybeUnwrapNxt_spec (v := 0) h.inv hs hty hok).step.inv
  match hts : s.tstack with
  | [] | [_] => have := hok.two; rw [hts] at this; simp at this
  | a :: b :: rest =>
    have hn : nxtE s = b := by rw [nxtE, hts]; rfl
    have hc : curE s = a := by rw [curE, hts]; rfl
    have hd : nxtDir s = s.stackDir[b.topDepth]! := by rw [nxtDir, hn]
    have hh : nxtHead s = (getSide b.spans s.stackDir[b.topDepth]!).head! := by rw [nxtHead, hn, hd]
    rw [maybeUnwrapNxt_run_eq ty s a b rest hts _ rfl _ rfl] at hinv ⊢
    by_cases h1 : ty = .R ∨ s.ternarize = true
    · simp only [h1, ↓reduceIte]; exact h.alloc' hσ hs ty
    simp only [h1, ↓reduceIte] at hinv ⊢
    by_cases h2 : s.items[(getSide b.spans s.stackDir[b.topDepth]!).head!]!.type = ty
    · simp only [h2, ↓reduceIte] at hinv ⊢
      obtain ⟨side, single, -, -⟩ := hok.unwrap h1 (by rw [hh]; exact h2)
      rw [hn, hd] at side single
      rw [hh] at single
      set dir := s.stackDir[b.topDepth]!
      set i := (getSide b.spans dir).head!
      have hib : i ∈ b.spans.1 ++ b.spans.2 := by
        rw [mem_of_getSide_nil dir b.spans side, single]; exact List.mem_singleton_self i
      have hilt : i < s.items.size := hs.span b (by simp [hts]) i hib
      have hget : s.items[i]! = s.items[i] := getElem!_pos s.items i hilt
      have hity : Items.type s.items i = ty := by rw [Items.type_eq_getElem hilt, ← hget]; exact h2
      have hinode := hs.node_of_type hilt hty hity
      have hie : ∀ e, e < s.g.ne → edgeItem s.g e ≠ i := fun e he => hs.edgeItem_ne he hinode
      have hch : s.items[i]!.ch = Items.ch s.items i := by rw [Items.ch_eq_getElem hilt, hget]
      have hiv : Items.type s.items i ≠ .V := by rw [hity]; rintro rfl; simp at hty
      set b' : TEntry := { b with spans := setSides dir s.items[i]!.ch [] } with hb'
      refine h.replaceNxt a b b' rest hts hσ hinv (fun e he => ?_) (fun e he => ?_)
      · rw [TEntry.edges_single dir i side single, hb', hch]
        exact TEntry.edges_unwrap dir b.vStart b.topDepth b.firstIdx i hie he
      · simp only [TEntry.piece, hb', mem_setSides, hch, mem_of_getSide_nil dir b.spans side, single,
          List.mem_singleton]
        constructor
        · rintro ⟨c, hc, hcv, hp⟩
          exact ⟨i, rfl, hiv, Relation.ReflTransGen.head ⟨hc, hcv⟩ hp⟩
        · rintro ⟨j, rfl, -, hp⟩
          rcases Relation.ReflTransGen.cases_head hp with heq | ⟨c, ⟨hc, hcv⟩, hp⟩
          · exact absurd heq.symm (hie e he)
          · exact ⟨c, hc, hcv, hp⟩
    · simp only [h2, ↓reduceIte]; exact h.alloc' hσ hs ty

/-- The push sites: the vertex item is a `V` item whose hanging blocks are processed, the edge item
is not a `V` item and sits at position `n` of `σ`. -/
structure PushVertR (σ : List Nat) (n v : Nat) (s : WalkState) : Prop where
  vtype : Items.type s.items (vertItem v) = .V
  below : ∀ e, e < s.g.ne → Items.EdgeBelow s.g s.items (vertItem v) e → σ.idxOf e < n

theorem RangesInv.pushVert' (h : s.RangesInv σ n D) (v d : Nat)
    (hc : s.g.ConnEdges (Items.EdgeBelow s.g s.items (vertItem v)))
    (ha : s.g.TwoAttached (Items.EdgeBelow s.g s.items (vertItem v)) v v) (hr : PushVertR σ n v s) :
    RangesInv σ n D ((pushVertTstack v d).run s).2 :=
  h.pushVert v d hc ha hr.vtype hr.below

theorem RangesInv.pushEdge' (h : s.RangesInv σ n D) (hnd : σ.Nodup) (vStart topDepth e : Nat)
    (he : e < s.g.ne) (hq : Items.ch s.items (edgeItem s.g e) = [])
    (htq : Items.type s.items (edgeItem s.g e) ≠ .V)
    (hend : Items.PairEq (vStart, s.stackVerts[topDepth]!) s.g.edges[e]!) (hD : topDepth ≤ D)
    (hn : σ[n]? = some e) : RangesInv σ (n + 1) D ((pushEdgeTstack vStart topDepth e).run s).2 :=
  h.pushEdge vStart topDepth e hnd he hq htq hend hD hn

/-! ### Steps carrying the range invariant -/

/-- A `Step` whose target satisfies the range invariant after `n` edges. -/
structure RgStep (σ : List Nat) (n D v : Nat) (s s' : WalkState) : Prop where
  step : Step D v s s'
  ranges : s'.RangesInv σ n D

variable {v : Nat}

theorem RgStep.trans {n' : Nat} {s₁ s₂ s₃ : WalkState} (h₁ : RgStep σ n D v s₁ s₂) (h₂ : RgStep σ n' D v s₂ s₃) :
    RgStep σ n' D v s₁ s₃ := ⟨h₁.step.trans h₂.step, h₂.ranges⟩

theorem RgStep.refl (h : s.RangesInv σ n D) (hs : Shape s) : RgStep σ n D v s s := ⟨Step.refl h.inv hs, h⟩

theorem RgStep.hσ (st : RgStep σ n D v s s') (hσ : ∀ e ∈ σ, e < s.g.ne) : ∀ e ∈ σ, e < s'.g.ne := by
  rw [st.step.g]; exact hσ

theorem RgStep.mergeTop (h : s.RangesInv σ n D) (hs : Shape s) (hnd : σ.Nodup) (hσ : ∀ e ∈ σ, e < s.g.ne)
    (hok : MergeTopOk D s) (hadj : MergeAdj σ s) : RgStep σ n D v s (mergeTstackTops.run s).2 :=
  ⟨Step.mergeTop h.inv hs hok, h.mergeTop' hnd hσ hok hadj⟩

theorem RgStep.finishTop (h : s.RangesInv σ n D) (hs : Shape s) (hσ : ∀ e ∈ σ, e < s.g.ne) (hv : v < s.g.nv)
    (item : ItemId) (hok : FinishTopOk D s) (hf : ItemFree s item) :
    RgStep σ n D v s ((finishTstackTop item).run s).2 :=
  ⟨Step.finishTop h.inv hs hv item hok hf, h.finishTop' hσ item hok hf⟩

theorem RgStep.closeTwo (h : s.RangesInv σ n D) (hs : Shape s) (hnd : σ.Nodup) (hσ : ∀ e ∈ σ, e < s.g.ne)
    (hv : v < s.g.nv) (item : ItemId) (hok : CloseTwoOk D s) (hadj : MergeAdj σ s) (hf : ItemFree s item) :
    RgStep σ n D v s ((finishTstackTop item).run (mergeTstackTops.run s).2).2 := by
  have st := RgStep.mergeTop (v := v) h hs hnd hσ hok.merge hadj
  exact st.trans (RgStep.finishTop st.ranges st.step.shape (st.hσ hσ) (by rw [st.step.g]; exact hv) item
    hok.finish hf.merge)

theorem RgStep.retarget (h : s.RangesInv σ n D) (hs : Shape s) (hσ : ∀ e ∈ σ, e < s.g.ne) (curV : Nat)
    (edgeDir : Bool) (hok : RetargetOk D curV s) : RgStep σ n D v s ((retarget curV edgeDir).run s).2 :=
  ⟨Step.retarget h.inv hs curV edgeDir hok, h.retarget' hσ curV edgeDir hok⟩

structure RUnwrapRes (σ : List Nat) (n D v : Nat) (s : WalkState) (r : ItemId × WalkState) : Prop where
  step : RgStep σ n D v s r.2
  free : ItemFree r.2 r.1

theorem RgStep.unwrapRes {ty : NodeType} (h : s.RangesInv σ n D) (hs : Shape s) (hσ : ∀ e ∈ σ, e < s.g.ne)
    (hty : ty ∉ [NodeType.F, .V, .Q]) (hok : UnwrapOk ty s) :
    RUnwrapRes σ n D v s ((maybeUnwrapNxt ty).run s) :=
  let r := maybeUnwrapNxt_spec (v := v) h.inv hs hty hok
  ⟨⟨r.step, h.unwrap' hσ hs hty hok⟩, r.free⟩

theorem RgStep.pushVert (h : s.RangesInv σ n D) (hs : Shape s) (d : Nat) (hv : v < s.g.nv)
    (hc : s.g.ConnEdges (Items.EdgeBelow s.g s.items (vertItem v)))
    (ha : s.g.TwoAttached (Items.EdgeBelow s.g s.items (vertItem v)) v v) (hr : PushVertR σ n v s) :
    RgStep σ n D v s ((pushVertTstack v d).run s).2 :=
  ⟨Step.pushVert h.inv hs d hv hc ha, h.pushVert' v d hc ha hr⟩

theorem RgStep.pushEdge (h : s.RangesInv σ n D) (hs : Shape s) (hnd : σ.Nodup) (vStart topDepth e : Nat)
    (he : e < s.g.ne) (hq : Items.ch s.items (edgeItem s.g e) = [])
    (htq : Items.type s.items (edgeItem s.g e) ≠ .V)
    (hend : Items.PairEq (vStart, s.stackVerts[topDepth]!) s.g.edges[e]!) (hD : topDepth ≤ D)
    (hn : σ[n]? = some e) : RgStep σ (n + 1) D v s ((pushEdgeTstack vStart topDepth e).run s).2 :=
  ⟨Step.pushEdge h.inv hs vStart topDepth e he hq hend hD, h.pushEdge' hnd vStart topDepth e he hq htq hend hD hn⟩

theorem RgStep.loop (cond : WalkM Bool) (body : WalkM Unit) (Ok Adj : WalkState → Prop) (fuel : Nat)
    (hcond : ∀ s, (cond.run s).2 = s)
    (hbody : ∀ s, v < s.g.nv → s.RangesInv σ n D → Shape s → (∀ e ∈ σ, e < s.g.ne) →
      (cond.run s).1 = true → Ok s → Adj s → RgStep σ n D v s (body.run s).2)
    (hv : v < s.g.nv) (h : s.RangesInv σ n D) (hs : Shape s) (hσ : ∀ e ∈ σ, e < s.g.ne)
    (hok : ∀ k, (∀ j, j ≤ k → (cond.run (iter body j s)).1 = true) → Ok (iter body k s))
    (hadj : ∀ k, (∀ j, j ≤ k → (cond.run (iter body j s)).1 = true) → Adj (iter body k s)) :
    RgStep σ n D v s ((loop fuel cond body).run s).2 := by
  induction fuel generalizing s with
  | zero => exact RgStep.refl h hs
  | succ fuel ih =>
    rw [loop_succ_run fuel cond body s (hcond s)]
    by_cases hc : (cond.run s).1 = true
    · simp only [hc, ↓reduceIte]
      have st := hbody s hv h hs hσ hc (hok 0 fun j hj => by rw [Nat.le_zero.1 hj]; exact hc)
        (hadj 0 fun j hj => by rw [Nat.le_zero.1 hj]; exact hc)
      refine st.trans (ih (by rw [st.step.g]; exact hv) st.ranges st.step.shape (st.hσ hσ)
        (fun k hk => hok (k + 1) fun j hj => ?_) (fun k hk => hadj (k + 1) fun j hj => ?_)) <;>
      · cases j with
        | zero => exact hc
        | succ j => exact hk j (Nat.le_of_succ_le_succ hj)
    · simp only [hc]; exact RgStep.refl h hs

/-! ### The blocks of `finishEdge` -/

theorem RgStep.loop1Type (h : s.RangesInv σ n D) (hs : Shape s) (hnd : σ.Nodup) (hσ : ∀ e ∈ σ, e < s.g.ne)
    {d : Nat} {edgeDir : Bool}
    (hok : (nxtE s).topDepth > d → MergeTopOk D { s with stackDir := s.stackDir.set! (nxtE s).topDepth edgeDir })
    (hadj : (nxtE s).topDepth > d → MergeAdj σ { s with stackDir := s.stackDir.set! (nxtE s).topDepth edgeDir }) :
    RgStep σ n D v s (l1S₁ d edgeDir s) := by
  unfold l1S₁ after; rw [loop1Type_run]
  by_cases hc : (nxtE s).topDepth > d
  · simp only [hc, ↓reduceIte]
    have st : RgStep σ n D v s { s with stackDir := s.stackDir.set! (nxtE s).topDepth edgeDir } :=
      ⟨Step.frame h.inv hs rfl rfl rfl rfl, h.frame rfl rfl rfl rfl⟩
    exact st.trans (RgStep.mergeTop st.ranges st.step.shape hnd hσ (hok hc) (hadj hc))
  · simp only [hc, ↓reduceIte]
    first | exact RgStep.refl h hs | (split <;> exact RgStep.refl h hs)

/-- Loop-1 body: adjacency at the S-merge and at the closing merge. -/
structure Loop1BodyAdj (σ : List Nat) (d : Nat) (edgeDir : Bool) (s : WalkState) : Prop where
  mergeS : (nxtE s).topDepth > d → MergeAdj σ { s with stackDir := s.stackDir.set! (nxtE s).topDepth edgeDir }
  close : MergeAdj σ (l1S₂ d edgeDir s)

theorem RgStep.loop1Body (h : s.RangesInv σ n D) (hs : Shape s) (hnd : σ.Nodup) (hσ : ∀ e ∈ σ, e < s.g.ne)
    (hv : v < s.g.nv) {d : Nat} {edgeDir : Bool} (hok : Loop1BodyOk D d edgeDir s)
    (hadj : Loop1BodyAdj σ d edgeDir s) : RgStep σ n D v s ((loop1Body d edgeDir).run s).2 := by
  have st₁ : RgStep σ n D v s (l1S₁ d edgeDir s) := RgStep.loop1Type h hs hnd hσ hok.mergeS hadj.mergeS
  have hσ₁ := st₁.hσ hσ
  have hv₁ : v < (l1S₁ d edgeDir s).g.nv := by rw [st₁.step.g]; exact hv
  have r := RgStep.unwrapRes (v := v) st₁.ranges st₁.step.shape hσ₁ (loop1Type_result d edgeDir s) hok.unwrap
  have st₃ := RgStep.closeTwo r.step.ranges r.step.step.shape hnd (r.step.hσ hσ₁)
    (by rw [r.step.step.g]; exact hv₁) _ hok.close hadj.close r.free
  exact st₁.trans (r.step.trans st₃)

/-- `closeEars`: the tree edge sits at position `n` of `σ`, and adjacency at every loop-1 merge. -/
structure CloseEarsAdj (σ : List Nat) (n nxtV d e : Nat) (edgeDir : Bool) (s : WalkState) : Prop where
  pos : σ[n]? = some e
  etype : Items.type s.items (edgeItem s.g e) ≠ .V
  body : ∀ k, (∀ j, j ≤ k → result (loop1Cond d) (iter (loop1Body d edgeDir) j (ceS₁ nxtV d e s)) = true) →
    Loop1BodyAdj σ d edgeDir (iter (loop1Body d edgeDir) k (ceS₁ nxtV d e s))

theorem RgStep.closeEars (h : s.RangesInv σ n D) (hs : Shape s) (hnd : σ.Nodup) (hσ : ∀ e ∈ σ, e < s.g.ne)
    (hv : v < s.g.nv) {nxtV d e : Nat} {edgeDir : Bool} (hok : CloseEarsOk D nxtV d e edgeDir s)
    (hadj : CloseEarsAdj σ n nxtV d e edgeDir s) :
    RgStep σ (n + 1) D v s ((closeEars nxtV d e edgeDir).run s).2 := by
  have st₁ : RgStep σ (n + 1) D v s (ceS₁ nxtV d e s) :=
    RgStep.pushEdge h hs hnd nxtV d e hok.e_lt hok.q hadj.etype hok.ends hok.d_le hadj.pos
  exact st₁.trans (RgStep.loop (loop1Cond d) (Spqr.loop1Body d edgeDir) (Loop1BodyOk D d edgeDir)
    (Loop1BodyAdj σ d edgeDir) _ (fun _ => rfl)
    (fun _ hv h hs hσ _ hok hadj => RgStep.loop1Body h hs hnd hσ hv hok hadj)
    (by rw [st₁.step.g]; exact hv) st₁.ranges st₁.step.shape (st₁.hσ hσ) hok.body hadj.body)

/-- `mergeLate`: adjacency at every late merge. -/
structure MergeLateAdj (σ : List Nat) (d : Nat) (s : WalkState) : Prop where
  body : (curE s).firstIdx > s.firstOccurrence[d]! → ∀ k,
    (∀ j, j ≤ k → result (loop2Cond s.firstOccurrence[d]!) (iter mergeTstackTops j s) = true) →
    MergeAdj σ (iter mergeTstackTops k s)

theorem RgStep.mergeLate (h : s.RangesInv σ n D) (hs : Shape s) (hnd : σ.Nodup) (hσ : ∀ e ∈ σ, e < s.g.ne)
    (hv : v < s.g.nv) {d : Nat} (hok : MergeLateOk D d s) (hadj : MergeLateAdj σ d s) :
    RgStep σ n D v s ((mergeLate d).run s).2 := by
  rw [mergeLate_run]
  by_cases hc : (curE s).firstIdx > s.firstOccurrence[d]!
  · simp only [hc, ↓reduceIte]
    exact RgStep.loop _ mergeTstackTops (MergeTopOk D) (MergeAdj σ) _ (fun _ => rfl)
      (fun _ _ h hs hσ _ hok hadj => RgStep.mergeTop h hs hnd hσ hok hadj) hv h hs hσ (hok.body hc) (hadj.body hc)
  · simp only [hc, ↓reduceIte]; exact RgStep.refl h hs

/-- `closeVert'`: adjacency at the loop-3 merges and the two merges into the vertex entry. -/
structure CloseVertAdj (σ : List Nat) (isType1 : Bool) (origTstack : Nat) (isSingle : Bool) (s : WalkState) :
    Prop where
  loop3 : isType1 = false → ∀ k,
    (∀ j, j ≤ k → result (loop3Cond origTstack) (iter mergeTstackTops j s) = true) →
    MergeAdj σ (iter mergeTstackTops k s)
  merge₁ : MergeAdj σ (cvS₂ isType1 origTstack isSingle s)
  merge₂ : MergeAdj σ (cvS₃ isType1 origTstack isSingle s)

theorem RgStep.vertPre (h : s.RangesInv σ n D) (hs : Shape s) (hnd : σ.Nodup) (hσ : ∀ e ∈ σ, e < s.g.ne)
    (hv : v < s.g.nv) {isType1 : Bool} {origTstack : Nat} {isSingle : Bool}
    (hok : isType1 = false → ∀ k,
      (∀ j, j ≤ k → result (loop3Cond origTstack) (iter mergeTstackTops j s) = true) →
      MergeTopOk D (iter mergeTstackTops k s))
    (hadj : isType1 = false → ∀ k,
      (∀ j, j ≤ k → result (loop3Cond origTstack) (iter mergeTstackTops j s) = true) →
      MergeAdj σ (iter mergeTstackTops k s)) :
    RgStep σ n D v s ((vertPre isType1 origTstack isSingle).run s).2 := by
  cases isType1
  · simp only [WalkState.vertPre, Bool.not_false, ↓reduceIte, WalkM.run_bind, run_tstackSize, WalkM.pure_run]
    exact RgStep.loop _ mergeTstackTops (MergeTopOk D) (MergeAdj σ) _ (fun _ => rfl)
      (fun _ _ h hs hσ _ hok hadj => RgStep.mergeTop h hs hnd hσ hok hadj) hv h hs hσ (hok rfl) (hadj rfl)
  · exact RgStep.refl h hs

theorem RgStep.vertUnwrap (h : s.RangesInv σ n D) (hs : Shape s) (hσ : ∀ e ∈ σ, e < s.g.ne)
    {isType1 isSingle : Bool} (hok : isType1 = true → UnwrapOk (if isSingle then .S else .R) s) :
    RgStep σ n D v s ((vertUnwrap isType1 isSingle).run s).2 ∧
      ∀ i, ((vertUnwrap isType1 isSingle).run s).1 = some i →
        ItemFree ((vertUnwrap isType1 isSingle).run s).2 i := by
  cases isType1
  · exact ⟨RgStep.refl h hs, fun i h => by simp [WalkState.vertUnwrap] at h⟩
  · have r := RgStep.unwrapRes (v := v) h hs hσ (ty := if isSingle then .S else .R)
      (by cases isSingle <;> decide) (hok rfl)
    simp only [WalkState.vertUnwrap, ↓reduceIte, WalkM.map_run, Option.some.injEq]
    exact ⟨r.step, fun i h => h ▸ r.free⟩

theorem RgStep.closeVert' (h : s.RangesInv σ n D) (hs : Shape s) (hnd : σ.Nodup) (hσ : ∀ e ∈ σ, e < s.g.ne)
    (hv : v < s.g.nv) {curV : Nat} {edgeDir isType1 : Bool} {origTstack : Nat} {isSingle : Bool}
    (hok : CloseVertOk D curV edgeDir isType1 origTstack isSingle s)
    (hadj : CloseVertAdj σ isType1 origTstack isSingle s) :
    RgStep σ n D v s ((closeVert' curV edgeDir isType1 origTstack isSingle).run s).2 := by
  have st₁ : RgStep σ n D v s (cvS₁ isType1 origTstack isSingle s) :=
    RgStep.vertPre h hs hnd hσ hv hok.loop3 hadj.loop3
  have hσ₁ := st₁.hσ hσ
  have hv₁ : v < (cvS₁ isType1 origTstack isSingle s).g.nv := by rw [st₁.step.g]; exact hv
  obtain ⟨st₂, hfree⟩ := RgStep.vertUnwrap (v := v) st₁.ranges st₁.step.shape hσ₁ (isType1 := isType1)
    (isSingle := cvB₁ isType1 origTstack isSingle s)
    (fun h => by subst h; exact hok.unwrap rfl)
  have st₂ : RgStep σ n D v (cvS₁ isType1 origTstack isSingle s) (cvS₂ isType1 origTstack isSingle s) := st₂
  have hσ₂ := st₂.hσ hσ₁
  have st₃ : RgStep σ n D v _ (cvS₃ isType1 origTstack isSingle s) :=
    RgStep.mergeTop st₂.ranges st₂.step.shape hnd hσ₂ hok.merge₁ hadj.merge₁
  have hσ₃ := st₃.hσ hσ₂
  have st₄ : RgStep σ n D v _ (cvS₄ isType1 origTstack isSingle s) :=
    RgStep.mergeTop st₃.ranges st₃.step.shape hnd hσ₃ hok.merge₂ hadj.merge₂
  have hσ₄ := st₄.hσ hσ₃
  have st₅ : RgStep σ n D v _ (cvS₅ curV edgeDir isType1 origTstack isSingle s) :=
    RgStep.retarget st₄.ranges st₄.step.shape hσ₄ curV edgeDir hok.retarget
  have hσ₅ := st₅.hσ hσ₄
  have hv₅ : v < (cvS₅ curV edgeDir isType1 origTstack isSingle s).g.nv := by
    rw [st₅.step.g, st₄.step.g, st₃.step.g, st₂.step.g]; exact hv₁
  have st := st₁.trans (st₂.trans (st₃.trans (st₄.trans st₅)))
  cases isType1
  · exact st
  · have hf : ItemFree (cvS₂ true origTstack isSingle s)
        ((maybeUnwrapNxt (if isSingle then .S else .R)).run s).1 := hfree _ rfl
    have hf₅ := ((hf.merge).merge).retarget curV edgeDir
    exact st.trans (RgStep.finishTop st₅.ranges st₅.step.shape hσ₅ hv₅ _ (hok.finish rfl) hf₅)

/-- `finishP`: adjacency at the P-close merge. -/
structure FinishPAdj (σ : List Nat) (curV lowval : Nat) (isType1 : Bool) (s : WalkState) : Prop where
  adj : result (condP curV lowval isType1) s = true → MergeAdj σ (after (maybeUnwrapNxt .P) s)

theorem RgStep.finishP (h : s.RangesInv σ n D) (hs : Shape s) (hnd : σ.Nodup) (hσ : ∀ e ∈ σ, e < s.g.ne)
    (hv : v < s.g.nv) {curV lowval : Nat} {isType1 : Bool} (hok : FinishPOk D curV lowval isType1 s)
    (hadj : FinishPAdj σ curV lowval isType1 s) : RgStep σ n D v s ((finishP curV lowval isType1).run s).2 := by
  simp only [Spqr.finishP, WalkM.run_bind, run_condP]
  by_cases hc : result (condP curV lowval isType1) s = true
  · have h' : (isType1 && decide (s.tstack.length ≥ 2) && (s.tstack.tail.head!.vStart == curV) &&
        (s.tstack.tail.head!.topDepth == lowval)) = true := hc
    simp only [h', ↓reduceIte, WalkM.run_bind]
    obtain ⟨hu, hcl⟩ := hok.ok hc
    have r := RgStep.unwrapRes (v := v) h hs hσ (by decide) hu
    exact r.step.trans (RgStep.closeTwo r.step.ranges r.step.step.shape hnd (r.step.hσ hσ)
      (by rw [r.step.step.g]; exact hv) _ hcl (hadj.adj hc) r.free)
  · have h' : (isType1 && decide (s.tstack.length ≥ 2) && (s.tstack.tail.head!.vStart == curV) &&
        (s.tstack.tail.head!.topDepth == lowval)) = false := Bool.eq_false_iff.2 hc
    simp only [h', Bool.false_eq_true, ↓reduceIte]
    exact RgStep.refl h hs

/-- `finishTail`: the vertex push data and adjacency at the merge into the first ear. -/
structure FinishTailAdj (σ : List Nat) (n curV d : Nat) (hasVert isSingle : Bool) (s : WalkState) : Prop where
  vert : hasVert = false → PushVertR σ n curV s
  merge : hasVert = false → isSingle = false → MergeAdj σ (after (pushVertTstack curV d) s)

theorem RgStep.finishTail (h : s.RangesInv σ n D) (hs : Shape s) (hnd : σ.Nodup) (hσ : ∀ e ∈ σ, e < s.g.ne)
    {curV d : Nat} (hv : curV < s.g.nv) {hasVert isSingle : Bool}
    (hok : FinishTailOk D curV d hasVert isSingle s) (hadj : FinishTailAdj σ n curV d hasVert isSingle s) :
    RgStep σ n D curV s ((finishTail curV d hasVert isSingle).run s).2 := by
  cases hasVert
  · have st₁ := RgStep.pushVert h hs d hv (hok.conn rfl) (hok.att rfl) (hadj.vert rfl)
    cases isSingle
    · exact st₁.trans (RgStep.mergeTop st₁.ranges st₁.step.shape hnd (st₁.hσ hσ) (hok.merge rfl rfl)
        (hadj.merge rfl rfl))
    · exact st₁
  · exact RgStep.refl h hs

structure FinishRestAdj (σ : List Nat) (n curV d lowval : Nat) (isType1 hasVert isSingle : Bool)
    (s : WalkState) : Prop where
  p : FinishPAdj σ curV lowval isType1 s
  tail : FinishTailAdj σ n curV d hasVert isSingle (after (finishP curV lowval isType1) s)

theorem RgStep.finishRest (h : s.RangesInv σ n D) (hs : Shape s) (hnd : σ.Nodup) (hσ : ∀ e ∈ σ, e < s.g.ne)
    {curV d lowval : Nat} (hv : curV < s.g.nv) {isType1 hasVert isSingle : Bool}
    (hok : FinishRestOk D curV d lowval isType1 hasVert isSingle s)
    (hadj : FinishRestAdj σ n curV d lowval isType1 hasVert isSingle s) :
    RgStep σ n D curV s ((finishRest curV d lowval isType1 hasVert isSingle).run s).2 := by
  have st₁ := RgStep.finishP (v := curV) h hs hnd hσ hv hok.p hadj.p
  exact st₁.trans (RgStep.finishTail st₁.ranges st₁.step.shape hnd (st₁.hσ hσ) (by rw [st₁.step.g]; exact hv)
    hok.tail hadj.tail)

/-- The range-side hypotheses of `finishEdge` for an out-edge `o.cls = .ret lv kind`, `lv < d`,
block by block, at the state where the block runs (mirrors `FinishOk`): the edge `o.e` is `σ[n]`,
its item is a `Q` item, and adjacency holds at every merge site. -/
structure FinishAdj (σ : List Nat) (n curV d lv : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool)
    (s : WalkState) : Prop where
  ears : o.cls.isTree = true → CloseEarsAdj σ n o.dest d o.e s.stackDir[d]! (feS₀ d o s)
  late : o.cls.isTree = true → MergeLateAdj σ d (feS₁ d o s)
  vert : o.cls.isTree = true → hasVert = true →
    CloseVertAdj σ o.cls.isType1 origTstack (feSingle d o s) (feS₂ d o s)
  rest_vert : o.cls.isTree = true → hasVert = true →
    FinishRestAdj σ (n + 1) curV d lv o.cls.isType1 hasVert (feB₃ curV d o origTstack s)
      (feS₃ curV d o origTstack s)
  rest_tree : o.cls.isTree = true → hasVert = false →
    FinishRestAdj σ (n + 1) curV d lv o.cls.isType1 hasVert (feSingle d o s) (feS₂ d o s)
  pos : o.cls.isTree = false → σ[n]? = some o.e
  etype : o.cls.isTree = false → Items.type (feS₀ d o s).items (edgeItem (feS₀ d o s).g o.e) ≠ .V
  rest_back : o.cls.isTree = false →
    FinishRestAdj σ (n + 1) curV d lv o.cls.isType1 hasVert true (feBack curV lv d o s)

/-- `finishEdge` preserves the range invariant on a returning edge (`lv < d`), advancing `n`, under
`FinishOk` and `FinishAdj`. -/
theorem finishEdge_rangesInv (curV d lv : Nat) (kind : RetKind) (o : DfsOut) (origTstack : Nat)
    (hasVert : Bool) (ho : o.cls = .ret lv kind) (hlow : lv < d) (hv : curV < s.g.nv)
    (h : s.RangesInv σ n D) (hs : Shape s) (hnd : σ.Nodup) (hσ : ∀ e ∈ σ, e < s.g.ne)
    (hok : FinishOk D curV d lv o origTstack hasVert s)
    (hadj : FinishAdj σ n curV d lv o origTstack hasVert s) :
    RangesInv σ (n + 1) D ((finishEdge curV d o origTstack hasVert).run s).2 := by
  have hlv : o.cls.lowval d = lv := by rw [ho]; rfl
  have hge : ¬ (lv ≥ d) := by omega
  rw [finishEdge_eq]
  simp only [finishEdge', finishTree, finishBack, hlv, hge, ↓reduceIte, WalkM.run_bind, WalkM.get_run,
    run_stackDir, run_makeVs, run_modifyItem]
  have hj : edgeItem s.g o.e < 1 + s.g.nv + s.g.ne := by
    show 1 + s.g.nv + o.e < _; have := hok.e_lt; omega
  have st₀ : RgStep σ n D curV s (feS₀ d o s) :=
    ⟨Step.modifyVs h.inv hs (edgeItem s.g o.e) _ hj, h.modifyVs (edgeItem s.g o.e) _ hj⟩
  have hσ₀ := st₀.hσ hσ
  have hv₀ : curV < (feS₀ d o s).g.nv := by rw [st₀.step.g]; exact hv
  by_cases ht : o.cls.isTree = true
  · simp only [ht, ↓reduceIte, closeVert_eq, WalkM.run_bind]
    have st₁ : RgStep σ (n + 1) D curV _ (feS₁ d o s) :=
      RgStep.closeEars st₀.ranges st₀.step.shape hnd hσ₀ hv₀ (hok.ears ht) (hadj.ears ht)
    have hσ₁ := st₁.hσ hσ₀
    have hv₁ : curV < (feS₁ d o s).g.nv := by rw [st₁.step.g]; exact hv₀
    have st₂ : RgStep σ (n + 1) D curV _ (feS₂ d o s) :=
      RgStep.mergeLate st₁.ranges st₁.step.shape hnd hσ₁ hv₁ (hok.late ht) (hadj.late ht)
    have hσ₂ := st₂.hσ hσ₁
    have hv₂ : curV < (feS₂ d o s).g.nv := by rw [st₂.step.g]; exact hv₁
    cases hasVert
    · simp only [Bool.false_eq_true, ↓reduceIte]
      exact (st₀.trans (st₁.trans (st₂.trans (RgStep.finishRest st₂.ranges st₂.step.shape hnd hσ₂ hv₂
        (hok.rest_tree ht rfl) (hadj.rest_tree ht rfl))))).ranges
    · simp only [↓reduceIte, WalkM.run_bind]
      have st₃ : RgStep σ (n + 1) D curV _ (feS₃ curV d o origTstack s) :=
        RgStep.closeVert' st₂.ranges st₂.step.shape hnd hσ₂ hv₂ (hok.vert ht rfl) (hadj.vert ht rfl)
      have hσ₃ := st₃.hσ hσ₂
      have hv₃ : curV < (feS₃ curV d o origTstack s).g.nv := by rw [st₃.step.g]; exact hv₂
      exact (st₀.trans (st₁.trans (st₂.trans (st₃.trans (RgStep.finishRest st₃.ranges st₃.step.shape hnd hσ₃
        hv₃ (hok.rest_vert ht rfl) (hadj.rest_vert ht rfl)))))).ranges
  · have ht' : o.cls.isTree = false := Bool.eq_false_iff.2 ht
    simp only [ht', Bool.false_eq_true, ↓reduceIte, WalkM.run_bind]
    have hq : Items.ch (feS₀ d o s).items (edgeItem (feS₀ d o s).g o.e) = [] := by
      show Items.ch (s.items.modify _ _) (edgeItem s.g o.e) = []
      rw [Items.ch_modify_ch_eq (edgeItem s.g o.e)
        (fun it => { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) }) (fun _ => rfl)]
      exact hok.q ht'
    have st₁ : RgStep σ (n + 1) D curV _ (after (pushEdgeTstack curV lv o.e) (feS₀ d o s)) :=
      RgStep.pushEdge st₀.ranges st₀.step.shape hnd curV lv o.e hok.e_lt hq (hadj.etype ht') (hok.ends ht')
        (hok.lv_le ht') (hadj.pos ht')
    have st₂ : RgStep σ (n + 1) D curV _ (feBack curV lv d o s) :=
      ⟨Step.frame st₁.step.inv st₁.step.shape rfl rfl rfl rfl, st₁.ranges.frame rfl rfl rfl rfl⟩
    have hσ₂ := st₂.hσ (st₁.hσ hσ₀)
    have hv₂ : curV < (feBack curV lv d o s).g.nv := by rw [st₂.step.g, st₁.step.g]; exact hv₀
    exact (st₀.trans (st₁.trans (st₂.trans (RgStep.finishRest st₂.ranges st₂.step.shape hnd hσ₂ hv₂
      (hok.rest_back ht') (hadj.rest_back ht'))))).ranges

end WalkState

end Spqr
