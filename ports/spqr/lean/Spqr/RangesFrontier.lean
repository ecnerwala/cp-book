import Spqr.RangesTree
import Spqr.EarSides
import Spqr.Proofs.ForestSpec

namespace Spqr.WalkState
open WalkM

variable {σ : List Nat} {n D : Nat} {s : WalkState}

theorem infix_interval (hnd : σ.Nodup) {xs : List Nat} (hx : xs <:+: σ) :
    ∀ a b c, a ≤ b → b ≤ c → c < σ.length → σ[a]! ∈ xs → σ[c]! ∈ xs → σ[b]! ∈ xs := by
  obtain ⟨pre, post, rfl⟩ := hx
  intro a b c hab hbc hc ha hcc
  have hi : ∀ e ∈ xs, (pre ++ xs ++ post).idxOf e = pre.length + xs.idxOf e := by
    intro e he
    have hn : e ∉ pre := by
      have hd := (List.nodup_append.1 (List.nodup_append.1 hnd).1).2.2
      exact fun hp => hd e hp e he rfl
    rw [List.append_assoc, List.idxOf_append_of_notMem hn, List.idxOf_append_of_mem he]
  have hal := List.idxOf_lt_length_iff.2 ha
  have hcl := List.idxOf_lt_length_iff.2 hcc
  have hia := hi _ ha
  have hic := hi _ hcc
  rw [idxOf_getElem! hnd (by omega : a < (pre ++ xs ++ post).length)] at hia
  rw [idxOf_getElem! hnd hc] at hic
  have hb : b < (pre ++ xs ++ post).length := by omega
  have hpre : pre.length ≤ b := by omega
  have hxs : b < (pre ++ xs).length := by simp only [List.length_append]; omega
  rw [getElem!_pos (pre ++ xs ++ post) b hb, List.getElem_append_left hxs, List.getElem_append_right hpre]
  exact List.getElem_mem _

theorem subEdges_iff_mem_block (o : DfsOut) (e : Nat) : subEdges o e ↔ e ∈ o.block := by
  cases o with
  | back => simp [subEdges, DfsOut.block, DfsOut.e]
  | tree e' cls child =>
    simp only [subEdges, DfsOut.e, DfsOut.block, List.mem_append, List.mem_singleton]
    rw [child.edgePostorder_perm_edges.mem_iff]
    exact or_comm

theorem subEdges_interval (hnd : σ.Nodup) {o : DfsOut} (ho : o.block <:+: σ) :
    ∀ a b c, a ≤ b → b ≤ c → c < σ.length →
      subEdges o σ[a]! → subEdges o σ[c]! → subEdges o σ[b]! := by
  simpa only [subEdges_iff_mem_block] using infix_interval hnd ho

/-- The ear's boundary split owns every gap before its pending tree edge. -/
theorem boundaryAdj_of_book {curV d orig : Nat} {o : DfsOut} {hasVert : Bool}
    (hge : d ≤ o.cls.lowval d) (hb : FinishBook curV d o orig hasVert s)
    (hnd : σ.Nodup) (hσ : ∀ e ∈ σ, e < s.g.ne) (hpos : σ[n]? = some o.e)
    (ho : o.block <:+: σ) : BoundaryAdj σ n d o s := by
  obtain ⟨sub, base, hl, he⟩ := hb.ear
  obtain ⟨hn, hnval⟩ := List.getElem?_eq_some_iff.1 hpos
  have hnval : σ[n]! = o.e := by rw [getElem!_pos σ n hn]; exact hnval
  have hfill : ∀ a b, a ≤ b → b < n →
      (∃ t ∈ sub, t.piece s.g s.items σ[a]!) → ∃ t ∈ sub, t.edges s.g s.items σ[b]! := by
    intro a b hab hbn ha
    obtain ⟨t, ht, ha⟩ := ha
    have hsub := he.sub_edges t ht _ (getElem!_lt hσ (by omega : a < σ.length)) ha.edges
    have hnext : subEdges o σ[n]! := by rw [hnval]; exact Or.inl rfl
    have hmid := subEdges_interval hnd ho a b n hab (Nat.le_of_lt hbn) hn hsub hnext
    apply he.sub_cover _ (getElem!_lt hσ (by omega : b < σ.length)) _ hmid
    intro heq
    have hid := idxOf_getElem! hnd (by omega : b < σ.length)
    rw [heq, ← hnval, idxOf_getElem! hnd hn] at hid
    omega
  refine ⟨hpos, ?_, ?_⟩
  · intro ht hlv t rest hts
    obtain ⟨u, hsub, -, -⟩ := he.bd_bridge ht hlv
    have hu := he.tstack
    rw [hsub, List.singleton_append, hts, List.cons.injEq] at hu
    obtain ⟨rfl, rfl⟩ := hu
    have hs := he.bd_side ht hge
    simp only [hlv, ↓reduceIte, hts, List.head?_cons, Option.mem_def, Option.some.injEq] at hs
    refine ⟨hs _ rfl, ?_⟩
    intro a b hab hbn ha
    obtain ⟨u, hu, huE⟩ := hfill a b hab hbn ⟨t, by simp [hsub], ha⟩
    have hu : u = t := by simpa [hsub] using hu
    exact hu ▸ huE
  · intro ht hlv t₁ t₂ rest hts
    obtain ⟨u₁, u₂, hsub, -, -, -, -⟩ := he.bd_comp ht hge hlv
    have hu := he.tstack
    simp only [hsub, List.cons_append, List.nil_append, hts, List.cons.injEq] at hu
    obtain ⟨rfl, rfl, rfl⟩ := hu
    have hs := he.bd_side ht hge
    simp only [hlv, ↓reduceIte, hts, List.head?_cons, List.tail_cons, Option.mem_def, Option.some.injEq] at hs
    refine ⟨hs.1 _ rfl, hs.2 _ rfl, ?_⟩
    intro a b hab hbn ha
    have hx : ∃ t ∈ sub, t.piece s.g s.items σ[a]! := by simpa [hsub] using ha
    simpa [hsub] using hfill a b hab hbn hx

/-- A covered gap between the top two pieces cannot belong to a lower entry. -/
theorem RangesInv.mergeAdj_of_cover (h : s.RangesInv σ n D) (hnd : σ.Nodup)
    (hσ : ∀ e ∈ σ, e < s.g.ne)
    (hcover : ∀ cur nxt rest, s.tstack = cur :: nxt :: rest →
      ∀ a b c, a ≤ b → b ≤ c → c < σ.length →
        nxt.piece s.g s.items σ[a]! → cur.piece s.g s.items σ[c]! →
        ∃ t ∈ s.tstack, t.edges s.g s.items σ[b]!) : MergeAdj σ s := by
  intro cur nxt rest hts a b c hab hbc hc ha hcc
  obtain ⟨t, ht, he⟩ := hcover cur nxt rest hts a b c hab hbc hc ha hcc
  rw [hts] at ht
  rcases List.mem_cons.1 ht with rfl | ht
  · exact Or.inl he
  rcases List.mem_cons.1 ht with rfl | ht
  · exact Or.inr he
  have hlt := h.ordered [cur] nxt rest hts t ht _ _
    (getElem!_lt hσ (by omega : a < σ.length)) (getElem!_lt hσ (by omega : b < σ.length)) ha he
  rw [idxOf_getElem! hnd (by omega : a < σ.length),
    idxOf_getElem! hnd (by omega : b < σ.length)] at hlt
  omega

/-- Exact frontier ownership of an interval supplies the missing merge adjacency. -/
theorem RangesInv.mergeAdj_of_frontier (h : s.RangesInv σ n D) (hnd : σ.Nodup)
    (hσ : ∀ e ∈ σ, e < s.g.ne) {orig : Nat} {base : List TEntry} {E : Nat → Prop}
    (hf : FrontierOwns orig base E s) (hlen : orig + 2 ≤ s.tstack.length)
    (hinterval : ∀ a b c, a ≤ b → b ≤ c → c < σ.length →
      E σ[a]! → E σ[c]! → E σ[b]!) : MergeAdj σ s := by
  apply h.mergeAdj_of_cover hnd hσ
  intro cur nxt rest hts a b c hab hbc hc ha hcc
  have hlen' : 2 ≤ s.tstack.length - orig := by omega
  have hfront : s.tstack.take (s.tstack.length - orig) =
      cur :: nxt :: rest.take (s.tstack.length - orig - 2) := by
    generalize s.tstack.length - orig = k at hlen' ⊢
    obtain ⟨k, rfl⟩ := Nat.exists_eq_add_of_le hlen'
    rw [hts, show 2 + k = k + 1 + 1 by omega]
    simp
  have he : E σ[b]! := hinterval a b c hab hbc hc
    ((hf.2.2 _ (getElem!_lt hσ (by omega : a < σ.length))).1 ⟨nxt, by simp [hfront], ha.edges⟩)
    ((hf.2.2 _ (getElem!_lt hσ hc)).1 ⟨cur, by simp [hfront], hcc.edges⟩)
  obtain ⟨t, ht, he⟩ := (hf.2.2 _ (getElem!_lt hσ (by omega : b < σ.length))).2 he
  exact ⟨t, List.mem_of_mem_take ht, he⟩

theorem RangesInv.iter_merge_ranges (h : s.RangesInv σ n D) (hnd : σ.Nodup)
    (hσ : ∀ e ∈ σ, e < s.g.ne) {cond : WalkM Bool} {orig : Nat} {base : List TEntry} {E : Nat → Prop}
    (hinterval : ∀ a b c, a ≤ b → b ≤ c → c < σ.length → E σ[a]! → E σ[c]! → E σ[b]!)
    (hok : ∀ k, (∀ j, j ≤ k → result cond (iter mergeTstackTops j s) = true) →
      MergeTopOk D (iter mergeTstackTops k s))
    (hf : ∀ k, (∀ j, j ≤ k → result cond (iter mergeTstackTops j s) = true) →
      FrontierOwns orig base E (iter mergeTstackTops k s) ∧ orig + 2 ≤ (iter mergeTstackTops k s).tstack.length) :
    ∀ k, (∀ j, j < k → result cond (iter mergeTstackTops j s) = true) →
      (iter mergeTstackTops k s).RangesInv σ n D ∧ (iter mergeTstackTops k s).g = s.g := by
  intro k
  induction k with
  | zero => exact fun _ => ⟨h, rfl⟩
  | succ k ih =>
    intro hk
    obtain ⟨hr, hg⟩ := ih fun j hj => hk j (by omega)
    have hk' : ∀ j, j ≤ k → result cond (iter mergeTstackTops j s) = true := fun j hj => hk j (by omega)
    have hσ' : ∀ e ∈ σ, e < (iter mergeTstackTops k s).g.ne := by rwa [hg]
    have ha := hr.mergeAdj_of_frontier hnd hσ' (hf k hk').1 (hf k hk').2 hinterval
    rw [iter_succ']
    exact ⟨hr.mergeTop' hnd hσ' (hok k hk') ha, by rw [run_mergeTstackTops]; exact hg⟩

theorem RangesInv.iter_mergeAdj (h : s.RangesInv σ n D) (hnd : σ.Nodup)
    (hσ : ∀ e ∈ σ, e < s.g.ne) {cond : WalkM Bool} {orig : Nat} {base : List TEntry} {E : Nat → Prop}
    (hinterval : ∀ a b c, a ≤ b → b ≤ c → c < σ.length → E σ[a]! → E σ[c]! → E σ[b]!)
    (hok : ∀ k, (∀ j, j ≤ k → result cond (iter mergeTstackTops j s) = true) →
      MergeTopOk D (iter mergeTstackTops k s))
    (hf : ∀ k, (∀ j, j ≤ k → result cond (iter mergeTstackTops j s) = true) →
      FrontierOwns orig base E (iter mergeTstackTops k s) ∧ orig + 2 ≤ (iter mergeTstackTops k s).tstack.length)
    (k : Nat) (hk : ∀ j, j ≤ k → result cond (iter mergeTstackTops j s) = true) :
    MergeAdj σ (iter mergeTstackTops k s) := by
  obtain ⟨hr, hg⟩ := h.iter_merge_ranges hnd hσ hinterval hok hf k (fun j hj => hk j (by omega))
  exact hr.mergeAdj_of_frontier hnd (by rwa [hg]) (hf k hk).1 (hf k hk).2 hinterval

theorem mergeLateAdj_of_frontier {o : DfsOut} {d orig : Nat}
    (hf : Frontier (o := o) d orig s) (ht : o.cls.isTree = true) (hlow : o.cls.lowval d < d)
    (h : (feS₁ d o s).RangesInv σ n D) (hnd : σ.Nodup)
    (hσ : ∀ e ∈ σ, e < (feS₁ d o s).g.ne) (ho : o.block <:+: σ)
    (hok : MergeLateOk D d (feS₁ d o s)) : MergeLateAdj σ d (feS₁ d o s) := by
  refine ⟨fun hc => h.iter_mergeAdj hnd hσ (orig := orig)
    (base := s.tstack.drop (s.tstack.length - orig)) (subEdges_interval hnd ho) (hok.body hc) ?_⟩
  intro k hk
  obtain ⟨hown, hlen⟩ := hf.loop2 ht hlow k (fun j hj => hk j (by omega))
  exact ⟨hown, hlen (hk k (Nat.le_refl _))⟩

theorem loop3Adj_of_frontier {o : DfsOut} {d orig curV : Nat}
    (hf : Frontier (o := o) d orig s) (ht : o.cls.isTree = true) (hlow : o.cls.lowval d < d)
    (h : (feS₂ d o s).RangesInv σ n D) (hnd : σ.Nodup)
    (hσ : ∀ e ∈ σ, e < (feS₂ d o s).g.ne) (ho : o.block <:+: σ)
    (hok : CloseVertOk D curV s.stackDir[d]! o.cls.isType1 orig (feSingle d o s) (feS₂ d o s))
    (hty : o.cls.isType1 = false) :
    ∀ k, (∀ j, j ≤ k → result (loop3Cond orig) (iter mergeTstackTops j (feS₂ d o s)) = true) →
      MergeAdj σ (iter mergeTstackTops k (feS₂ d o s)) := by
  apply h.iter_mergeAdj hnd hσ (orig := orig)
    (base := s.tstack.drop (s.tstack.length - orig)) (subEdges_interval hnd ho) (hok.loop3 hty)
  intro k hk
  obtain ⟨hown, hlen⟩ := hf.loop3 ht hlow hty k (fun j hj => hk j (by omega))
  exact ⟨hown, hlen (hk k (Nat.le_refl _))⟩

theorem FrontierOwns.replaceNxt {orig : Nat} {base : List TEntry} {E : Nat → Prop}
    (hf : FrontierOwns orig base E s) (hlen : orig + 2 ≤ s.tstack.length)
    {a b b' : TEntry} {rest : List TEntry} (hts : s.tstack = a :: b :: rest) {s' : WalkState}
    (hg : s'.g = s.g) (hts' : s'.tstack = a :: b' :: rest)
    (hE : ∀ (t : TEntry) e, e < s.g.ne → (t.edges s'.g s'.items e ↔ t.edges s.g s.items e))
    (hB : ∀ e, e < s.g.ne → (b'.edges s'.g s'.items e ↔ b.edges s.g s.items e)) :
    FrontierOwns orig base E s' := by
  have hr : orig ≤ rest.length := by rw [hts] at hlen; simpa using hlen
  have hlen' : s'.tstack.length = s.tstack.length := by simp [hts, hts']
  have hcalc : s.tstack.length - orig = rest.length - orig + 2 := by
    rw [hts]; simp only [List.length_cons]; omega
  refine ⟨by rw [hlen']; exact hf.1, ?_, ?_⟩
  · have hd := hf.2.1
    rw [hcalc, hts] at hd
    rw [hlen', hcalc, hts']
    simpa only [List.drop_succ_cons] using hd
  · intro e he
    rw [hg] at he
    rw [← hf.2.2 e he, hlen', hcalc, hts, hts']
    simp only [List.take_succ_cons, List.mem_cons, exists_eq_or_imp]
    rw [hB e he]
    simp only [hE _ e he]

theorem FrontierOwns.mergeTop {orig : Nat} {base : List TEntry} {E : Nat → Prop}
    (hf : FrontierOwns orig base E s) (hlen : orig + 2 ≤ s.tstack.length) :
    FrontierOwns orig base E (after mergeTstackTops s) ∧
      (after mergeTstackTops s).tstack.length + 1 = s.tstack.length := by
  match hts : s.tstack with
  | [] => rw [hts] at hlen; simp at hlen
  | [a] => rw [hts] at hlen; simp at hlen
  | a :: b :: rest =>
    have hr : orig ≤ rest.length := by rw [hts] at hlen; simpa using hlen
    have hcalc : s.tstack.length - orig = rest.length - orig + 2 := by
      rw [hts]; simp only [List.length_cons]; omega
    have hcalc' : rest.length + 1 - orig = rest.length - orig + 1 := by omega
    rw [after, run_mergeTstackTops, hts, mergeTops]
    refine ⟨⟨by simp; omega, ?_, ?_⟩, by simp⟩
    · have hd := hf.2.1
      rw [hcalc, hts] at hd
      simpa only [List.length_cons, hcalc', List.drop_succ_cons] using hd
    · intro e he
      rw [← hf.2.2 e he, hcalc, hts]
      simp only [List.length_cons, hcalc', List.take_succ_cons, List.mem_cons, exists_eq_or_imp]
      change ((TEntry.mergeInto a b).edges s.g s.items e ∨ _) ↔ _
      rw [TEntry.edges_mergeInto, or_assoc]

theorem FrontierOwns.congr {orig : Nat} {base : List TEntry} {E : Nat → Prop} {s' : WalkState}
    (hf : FrontierOwns orig base E s) (hg : s'.g = s.g) (hts : s'.tstack = s.tstack)
    (hE : ∀ (t : TEntry) e, t.edges s'.g s'.items e ↔ t.edges s.g s.items e) :
    FrontierOwns orig base E s' := by
  refine ⟨by rw [hts]; exact hf.1, by rw [hts]; exact hf.2.1, ?_⟩
  intro e he
  rw [hg] at he
  simpa only [hts, hE] using hf.2.2 e he

theorem FrontierOwns.unwrap {orig : Nat} {base : List TEntry} {E : Nat → Prop}
    (hf : FrontierOwns orig base E s) (hlen : orig + 2 ≤ s.tstack.length)
    {ty : NodeType} (hs : Shape s) (hty : ty ∉ [NodeType.F, .V, .Q]) (hok : UnwrapOk ty s) :
    FrontierOwns orig base E (after (maybeUnwrapNxt ty) s) ∧
      (after (maybeUnwrapNxt ty) s).tstack.length = s.tstack.length := by
  have halloc : FrontierOwns orig base E (after (allocItem ty) s) ∧
      (after (allocItem ty) s).tstack.length = s.tstack.length := by
    refine ⟨hf.congr rfl rfl ?_, rfl⟩
    exact fun t e => TEntry.edges_congr (fun i _ e => Items.Below_push_nil _ rfl) e
  match hts : s.tstack with
  | [] | [_] => have := hok.two; rw [hts] at this; simp at this
  | a :: b :: rest =>
    have hn : nxtE s = b := by rw [nxtE, hts]; rfl
    have hd : nxtDir s = s.stackDir[b.topDepth]! := by rw [nxtDir, hn]
    have hh : nxtHead s = (getSide b.spans s.stackDir[b.topDepth]!).head! := by rw [nxtHead, hn, hd]
    rw [after, maybeUnwrapNxt_run_eq ty s a b rest hts _ rfl _ rfl]
    by_cases h1 : ty = .R ∨ s.ternarize = true
    · simpa only [h1, ↓reduceIte, after, hts] using halloc
    simp only [h1, ↓reduceIte]
    by_cases h2 : s.items[(getSide b.spans s.stackDir[b.topDepth]!).head!]!.type = ty
    · simp only [h2, ↓reduceIte]
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
      have hie : ∀ e, e < s.g.ne → edgeItem s.g e ≠ i := fun e he =>
        hs.edgeItem_ne he (hs.node_of_type hilt hty hity)
      have hch : s.items[i]!.ch = Items.ch s.items i := by rw [Items.ch_eq_getElem hilt, hget]
      refine ⟨hf.replaceNxt hlen hts rfl rfl (fun _ _ _ => Iff.rfl) ?_, by simp⟩
      intro e he
      rw [TEntry.edges_single dir i side single, hch]
      exact TEntry.edges_unwrap dir b.vStart b.topDepth b.firstIdx i hie he
    · simpa only [h2, ↓reduceIte, after, hts] using halloc

theorem loop1BodyAdj_of_frontier {orig d : Nat} {base : List TEntry} {E : Nat → Prop} {edgeDir : Bool}
    (hf : FrontierOwns orig base E s)
    (hlen : orig + (if d < (nxtE s).topDepth then 3 else 2) ≤ s.tstack.length)
    (h : s.RangesInv σ n D) (hs : Shape s) (hnd : σ.Nodup) (hσ : ∀ e ∈ σ, e < s.g.ne)
    (hok : Loop1BodyOk D d edgeDir s)
    (hinterval : ∀ a b c, a ≤ b → b ≤ c → c < σ.length → E σ[a]! → E σ[c]! → E σ[b]!) :
    Loop1BodyAdj σ d edgeDir s := by
  have hm : orig + 2 ≤ s.tstack.length := by split_ifs at hlen <;> omega
  have ha : ∀ hlt : d < (nxtE s).topDepth,
      MergeAdj σ { s with stackDir := s.stackDir.set! (nxtE s).topDepth edgeDir } := by
    intro _
    exact (h.frame rfl rfl rfl rfl).mergeAdj_of_frontier hnd hσ
      (hf.congr rfl rfl (fun _ _ => Iff.rfl)) hm hinterval
  have st : RgStep σ n D 0 s (l1S₁ d edgeDir s) := RgStep.loop1Type h hs hnd hσ hok.mergeS ha
  have hf₁ : FrontierOwns orig base E (l1S₁ d edgeDir s) ∧ orig + 2 ≤ (l1S₁ d edgeDir s).tstack.length := by
    unfold l1S₁ after
    rw [loop1Type_run]
    by_cases hc : d < (nxtE s).topDepth
    · simp only [hc, ↓reduceIte] at hlen ⊢
      have hx := (hf.congr (s' := { s with stackDir := s.stackDir.set! (nxtE s).topDepth edgeDir })
        rfl rfl (fun _ _ => Iff.rfl)).mergeTop hm
      refine ⟨hx.1, ?_⟩
      have hl := hx.2
      dsimp only [after] at hl
      omega
    · simp only [hc, ↓reduceIte]
      split <;> exact ⟨hf, hm⟩
  have hu := hf₁.1.unwrap hf₁.2 st.step.shape (loop1Type_result d edgeDir s) hok.unwrap
  dsimp only [after] at hu
  have hr := st.ranges.unwrap' (st.hσ hσ) st.step.shape (loop1Type_result d edgeDir s) hok.unwrap
  have hσ₂ := (maybeUnwrapNxt_spec (v := 0) st.ranges.inv st.step.shape
    (loop1Type_result d edgeDir s) hok.unwrap).step.g
  refine ⟨ha, hr.mergeAdj_of_frontier hnd ?_ hu.1 (by rw [hu.2]; exact hf₁.2) hinterval⟩
  rw [hσ₂, st.step.g]
  exact hσ

theorem RgStep.iter_loop1_ranges {v d orig : Nat} {edgeDir : Bool} {base : List TEntry} {E : Nat → Prop}
    (h : s.RangesInv σ n D) (hs : Shape s) (hv : v < s.g.nv)
    (hnd : σ.Nodup) (hσ : ∀ e ∈ σ, e < s.g.ne)
    (hinterval : ∀ a b c, a ≤ b → b ≤ c → c < σ.length → E σ[a]! → E σ[c]! → E σ[b]!)
    (hok : ∀ k, (∀ j, j ≤ k → result (loop1Cond d) (iter (Spqr.loop1Body d edgeDir) j s) = true) →
      Loop1BodyOk D d edgeDir (iter (Spqr.loop1Body d edgeDir) k s))
    (hf : ∀ k, (∀ j, j ≤ k → result (loop1Cond d) (iter (Spqr.loop1Body d edgeDir) j s) = true) →
      let sk := iter (Spqr.loop1Body d edgeDir) k s
      FrontierOwns orig base E sk ∧ orig + (if d < (nxtE sk).topDepth then 3 else 2) ≤ sk.tstack.length) :
    ∀ k, (∀ j, j < k → result (loop1Cond d) (iter (Spqr.loop1Body d edgeDir) j s) = true) →
      RgStep σ n D v s (iter (Spqr.loop1Body d edgeDir) k s) := by
  intro k
  induction k with
  | zero => exact fun _ => RgStep.refl h hs
  | succ k ih =>
    intro hk
    have st := ih fun j hj => hk j (by omega)
    have hk' : ∀ j, j ≤ k → result (loop1Cond d) (iter (Spqr.loop1Body d edgeDir) j s) = true :=
      fun j hj => hk j (by omega)
    have ha := loop1BodyAdj_of_frontier (hf k hk').1 (hf k hk').2 st.ranges st.step.shape hnd
      (st.hσ hσ) (hok k hk') hinterval
    rw [iter_succ']
    exact st.trans (RgStep.loop1Body st.ranges st.step.shape hnd (st.hσ hσ)
      (by rw [st.step.g]; exact hv) (hok k hk') ha)

theorem closeEarsAdj_of_frontier {o : DfsOut} {d orig v : Nat}
    (hf : Frontier (o := o) d orig s) (ht : o.cls.isTree = true) (hlow : o.cls.lowval d < d)
    (h : (feS₀ d o s).RangesInv σ n D) (hs : Shape (feS₀ d o s)) (hv : v < (feS₀ d o s).g.nv)
    (hnd : σ.Nodup) (hσ : ∀ e ∈ σ, e < (feS₀ d o s).g.ne) (ho : o.block <:+: σ)
    (hpos : σ[n]? = some o.e) (hety : Items.type (feS₀ d o s).items (edgeItem (feS₀ d o s).g o.e) ≠ .V)
    (hok : CloseEarsOk D o.dest d o.e s.stackDir[d]! (feS₀ d o s)) :
    CloseEarsAdj σ n o.dest d o.e s.stackDir[d]! (feS₀ d o s) := by
  have st : RgStep σ (n + 1) D v (feS₀ d o s) (ceS₁ o.dest d o.e (feS₀ d o s)) :=
    RgStep.pushEdge h hs hnd _ _ _ hok.e_lt hok.q hety hok.ends hok.d_le hpos
  have hf' : ∀ k, (∀ j, j ≤ k → result (loop1Cond d)
      (iter (loop1Body d s.stackDir[d]!) j (ceS₁ o.dest d o.e (feS₀ d o s))) = true) →
      let sk := iter (loop1Body d s.stackDir[d]!) k (ceS₁ o.dest d o.e (feS₀ d o s))
      FrontierOwns orig (s.tstack.drop (s.tstack.length - orig)) (subEdges o) sk ∧
        orig + (if d < (nxtE sk).topDepth then 3 else 2) ≤ sk.tstack.length := by
    intro k hk
    obtain ⟨hown, hlen⟩ := hf.loop1 ht hlow k (fun j hj => hk j (by omega))
    exact ⟨hown, hlen (hk k (Nat.le_refl _))⟩
  refine ⟨hpos, hety, fun k hk => ?_⟩
  have sk := RgStep.iter_loop1_ranges st.ranges st.step.shape (v := v) (by rw [st.step.g]; exact hv)
    hnd (st.hσ hσ) (subEdges_interval hnd ho) hok.body hf' k (fun j hj => hk j (by omega))
  exact loop1BodyAdj_of_frontier (hf' k hk).1 (hf' k hk).2 sk.ranges sk.step.shape hnd
    (sk.hσ (st.hσ hσ)) (hok.body k hk) (subEdges_interval hnd ho)

theorem finishTailAdj_of_vert {curV d : Nat} {hasVert isSingle : Bool}
    (hv : hasVert = false → PushVertR σ n curV s) : FinishTailAdj σ n curV d hasVert isSingle s := by
  refine ⟨hv, fun hh _ cur nxt rest hts a b c _ _ _ _ hp => ?_⟩
  have ht := (hv hh).vtype
  rw [after, pushVertTstack, run_pushTstack] at hts
  have hc := (List.cons.inj hts).1
  subst cur
  exact (TEntry.piece_vertEntry _ _ _ _ ht _ hp).elim

theorem Frontier.owns_late {o : DfsOut} {d orig : Nat}
    (hf : Frontier (o := o) d orig s) (ht : o.cls.isTree = true) (hlow : o.cls.lowval d < d) :
    FrontierOwns orig (s.tstack.drop (s.tstack.length - orig)) (subEdges o) (feS₂ d o s) := by
  unfold feS₂ after
  rw [mergeLate_run]
  split
  · obtain ⟨k, -, hk, hj, -⟩ := loop_run_iter (feS₁ d o s).tstack.length
      (loop2Cond (feS₁ d o s).firstOccurrence[d]!) mergeTstackTops (feS₁ d o s) (fun _ => rfl)
    rw [hk]
    exact (hf.loop2 ht hlow k hj).1
  · exact (hf.loop2 ht hlow 0 (fun j hj => by omega)).1

theorem closeVertAdj_of_frontier {o : DfsOut} {d orig curV : Nat}
    (hf : Frontier (o := o) d orig s) (ht : o.cls.isTree = true) (hlow : o.cls.lowval d < d)
    (h : (feS₂ d o s).RangesInv σ n D) (hs : Shape (feS₂ d o s)) (hv : curV < (feS₂ d o s).g.nv)
    (hnd : σ.Nodup) (hσ : ∀ e ∈ σ, e < (feS₂ d o s).g.ne) (ho : o.block <:+: σ)
    (hlen : orig + 3 ≤ (feS₂ d o s).tstack.length)
    (hok : CloseVertOk D curV s.stackDir[d]! o.cls.isType1 orig (feSingle d o s) (feS₂ d o s)) :
    CloseVertAdj σ o.cls.isType1 orig (feSingle d o s) (feS₂ d o s) := by
  have ha := loop3Adj_of_frontier hf ht hlow h hnd hσ ho hok
  have st₁ : RgStep σ n D curV (feS₂ d o s)
      (cvS₁ o.cls.isType1 orig (feSingle d o s) (feS₂ d o s)) :=
    RgStep.vertPre h hs hnd hσ hv hok.loop3 ha
  have hf₁ : FrontierOwns orig (s.tstack.drop (s.tstack.length - orig)) (subEdges o)
      (cvS₁ o.cls.isType1 orig (feSingle d o s) (feS₂ d o s)) ∧
      orig + 3 ≤ (cvS₁ o.cls.isType1 orig (feSingle d o s) (feS₂ d o s)).tstack.length := by
    cases hty : o.cls.isType1
    · simp only [cvS₁, after, vertPre, Bool.not_false, ↓reduceIte, WalkM.run_bind, WalkM.pure_run, run_tstackSize]
      obtain ⟨k, -, hk, hj, -⟩ := loop_run_iter (feS₂ d o s).tstack.length (loop3Cond orig)
        mergeTstackTops (feS₂ d o s) (fun _ => rfl)
      rw [hk]
      refine ⟨(hf.loop3 ht hlow hty k hj).1, ?_⟩
      exact iter_merge_length_ge _ (by omega) k _ hlen fun j hj' => by
        have hc := hj j hj'
        rw [run_loop3Cond] at hc
        simpa using hc
    · exact ⟨hf.owns_late ht hlow, hlen⟩
  have hf₂ : FrontierOwns orig (s.tstack.drop (s.tstack.length - orig)) (subEdges o)
      (cvS₂ o.cls.isType1 orig (feSingle d o s) (feS₂ d o s)) ∧
      orig + 3 ≤ (cvS₂ o.cls.isType1 orig (feSingle d o s) (feS₂ d o s)).tstack.length := by
    cases hty : o.cls.isType1
    · rw [hty] at hf₁; exact hf₁
    · rw [hty] at hf₁ st₁ hok
      have hu := hf₁.1.unwrap (by omega) st₁.step.shape
        (ty := if feSingle d o s then .S else .R) (by split <;> decide) (hok.unwrap rfl)
      refine ⟨hu.1, ?_⟩
      change orig + 3 ≤ (after (maybeUnwrapNxt (if feSingle d o s then .S else .R))
        (cvS₁ true orig (feSingle d o s) (feS₂ d o s))).tstack.length
      rw [hu.2]
      exact hf₁.2
  obtain ⟨st₂, -⟩ := RgStep.vertUnwrap (v := curV) st₁.ranges st₁.step.shape (st₁.hσ hσ)
    (isType1 := o.cls.isType1) (isSingle := cvB₁ o.cls.isType1 orig (feSingle d o s) (feS₂ d o s))
    (fun hty => by rw [hty]; exact hok.unwrap hty)
  have st₂ : RgStep σ n D curV (cvS₁ o.cls.isType1 orig (feSingle d o s) (feS₂ d o s))
      (cvS₂ o.cls.isType1 orig (feSingle d o s) (feS₂ d o s)) := st₂
  have ha₁ := st₂.ranges.mergeAdj_of_frontier hnd (st₂.hσ (st₁.hσ hσ)) hf₂.1 (by omega)
    (subEdges_interval hnd ho)
  have st₃ := RgStep.mergeTop (v := curV) st₂.ranges st₂.step.shape hnd (st₂.hσ (st₁.hσ hσ)) hok.merge₁ ha₁
  have hf₃ := hf₂.1.mergeTop (by omega)
  refine ⟨ha, ha₁, st₃.ranges.mergeAdj_of_frontier hnd (st₃.hσ (st₂.hσ (st₁.hσ hσ))) hf₃.1 ?_
    (subEdges_interval hnd ho)⟩
  have hl := hf₃.2
  dsimp only [after] at hl
  omega

theorem FinishBook.late_length {curV d orig : Nat} {o : DfsOut} {hasVert : Bool}
    (hb : FinishBook curV d o orig hasVert s) (ht : o.cls.isTree = true) (hlow : o.cls.lowval d < d) :
    orig + 3 ≤ (feS₂ d o s).tstack.length := by
  obtain ⟨sub, base, hlen, he⟩ := hb.ear
  obtain ⟨c, mid, py, vy, hts, -⟩ := he.loops ht hlow
  rw [hts]
  simp only [List.length_cons, List.length_append, List.length_nil]
  omega

theorem finishPAdj_of_frontier {curV lowval orig : Nat} {isType1 : Bool}
    {base : List TEntry} {E : Nat → Prop}
    (h : s.RangesInv σ n D) (hs : Shape s) (hnd : σ.Nodup) (hσ : ∀ e ∈ σ, e < s.g.ne)
    (hok : FinishPOk D curV lowval isType1 s)
    (hf : result (condP curV lowval isType1) s = true →
      FrontierOwns orig base E s ∧ orig + 2 ≤ s.tstack.length)
    (hinterval : ∀ a b c, a ≤ b → b ≤ c → c < σ.length → E σ[a]! → E σ[c]! → E σ[b]!) :
    FinishPAdj σ curV lowval isType1 s := by
  refine ⟨fun hc => ?_⟩
  have hu := (hok.ok hc).1
  have hf' := (hf hc).1.unwrap (hf hc).2 hs (by decide) hu
  have hr := h.unwrap' hσ hs (by decide) hu
  have hg := (maybeUnwrapNxt_spec (v := curV) h.inv hs (by decide) hu).step.g
  dsimp only [after] at hf'
  exact hr.mergeAdj_of_frontier hnd (by rwa [hg]) hf'.1 (by rw [hf'.2]; exact (hf hc).2) hinterval

end Spqr.WalkState
