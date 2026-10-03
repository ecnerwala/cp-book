import Spqr.Proofs.RInvFrame

/-!
# `RInvTop` across the back-edge branch of `finishEdge`

Per-primitive preservation of `RInvTop dfs v d` for the operations of the back-edge branch
(`finishBack`): the pushed back-edge entry, the type-1 P close (`maybeUnwrapNxt .P`,
`mergeTstackTops`, `finishTstackTop`) and the bookkeeping update. Every entry these operations
create or rewrite starts at the current vertex `v`, so it is exempt from `RInvTop.entries`; the
work is the whole-stack edge-disjointness and the frame for the untouched entries.
-/

namespace Spqr
open WalkM
namespace WalkState
variable {s : WalkState} {dfs : DfsData}

/-- Pushing the single edge entry of `e` (a childless Q item) starting at the exempt vertex `v`:
no open entry owns `e`. -/
theorem RInvTop.pushEdge {v k d e : Nat} (hq : Items.ch s.items (edgeItem s.g e) = [])
    (hown : ∀ t ∈ s.tstack, ¬ t.edges s.g s.items e) (h : s.RInvTop dfs v d) :
    (after (pushEdgeTstack v k e) s).RInvTop dfs v d := by
  show RInvTop { s with tstack := ⟨v, k, s.nxtEdgeIdx, setSides s.stackDir[k]! [edgeItem s.g e] []⟩ :: s.tstack } dfs v d
  have hE := TEntry.edges_edgeEntry (g := s.g) (items := s.items) s.stackDir[k]! v k s.nxtEdgeIdx e hq
  refine ⟨fun t ht hd hne => ?_, List.pairwise_cons.2 ⟨fun t' ht' e' he' => ?_, h.disj⟩⟩
  · rcases List.mem_cons.1 ht with rfl | ht
    · exact absurd rfl hne
    · exact EntryR.congr (s := s) rfl rfl (fun _ _ => rfl) (fun _ _ => rfl) (fun _ _ _ => Iff.rfl)
        (h.entries t ht hd hne)
  · rw [(hE e').1 he']; exact hown t' ht'

/-- Merging the top two entries when the lower one starts at the exempt vertex `v`. -/
theorem RInvTop.mergeTop_exempt {v d : Nat} (cur nxt : TEntry) (rest : List TEntry)
    (hts : s.tstack = cur :: nxt :: rest) (hnxt : nxt.vStart = v) (h : s.RInvTop dfs v d) :
    (after mergeTstackTops s).RInvTop dfs v d := by
  show (mergeTstackTops.run s).2.RInvTop dfs v d
  rw [mergeTstackTops_run_eq s cur nxt rest hts]
  have hd := h.disj; rw [hts] at hd
  obtain ⟨hcur, hd⟩ := List.pairwise_cons.1 hd
  obtain ⟨hnxt', hrest⟩ := List.pairwise_cons.1 hd
  refine ⟨fun t ht hdep hne => ?_, List.pairwise_cons.2 ⟨fun t' ht' e he => ?_, hrest⟩⟩
  · rcases List.mem_cons.1 ht with rfl | ht
    · exact absurd hnxt hne
    · exact EntryR.congr (s := s) rfl rfl (fun _ _ => rfl) (fun _ _ => rfl) (fun _ _ _ => Iff.rfl)
        (h.entries t (by simp [hts, ht]) hdep hne)
  · rcases (TEntry.edges_mergeInto cur nxt e).1 he with he | he
    · exact hcur t' (by simp [ht']) e he
    · exact hnxt' t' ht' e he

/-- Finishing the top entry (starting at the exempt vertex `v`) into a fresh root item: the new
single-item entry owns exactly the old entry's edges, the other entries are untouched. -/
theorem RInvTop.finishTop_exempt {v d : Nat} (item : ItemId) (t : TEntry) (rest : List TEntry)
    (hts : s.tstack = t :: rest) (hv : t.vStart = v) (hitem : item < s.items.size)
    (hroot : ∀ p, ¬ Items.IsParent s.items p item)
    (hfree : ∀ t' ∈ s.tstack, item ∉ t'.spans.1 ++ t'.spans.2)
    (hside : getSide t.spans (!s.stackDir[t.topDepth]!) = [])
    (h : s.RInvTop dfs v d) : (after (finishTstackTop item) s).RInvTop dfs v d := by
  show ((finishTstackTop item).run s).2.RInvTop dfs v d
  rw [finishTstackTop_run_eq s item t rest hts]
  set dir := s.stackDir[t.topDepth]!
  set f : Item → Item := fun it =>
    { it with vs := setSides dir (some s.stackVerts[t.topDepth]!) (some t.vStart),
              ch := getSide t.spans dir }
  have hnb : ∀ i, i ≠ item → ¬ Items.Below s.items i item := fun i hi hb =>
    hi (Items.Below.eq_of_no_parent hroot hb)
  have hch : Items.ch (s.items.modify item f) item = getSide t.spans dir :=
    Items.ch_modify_at item f hitem
  have hne : ∀ c ∈ getSide t.spans dir, c ≠ item := fun c hc hce => by
    subst hce; exact hfree t (by simp [hts]) ((mem_of_getSide_nil dir t.spans hside c).2 hc)
  -- the new top owns a subset of the old top's edges (equality off the junk index `item`)
  have hsub : ∀ e, TEntry.edges s.g (s.items.modify item f) { t with spans := setSides dir [item] [] } e →
      t.edges s.g s.items e ∨ edgeItem s.g e = item := by
    rintro e ⟨i, hi, hb⟩
    rw [List.mem_singleton.1 ((mem_setSides dir [item] i).1 hi)] at hb
    rcases hb.head_cases with heq | ⟨c, hc, hb⟩
    · exact .inr heq.symm
    · simp only [Items.IsParent, hch] at hc
      exact .inl ⟨c, (mem_of_getSide_nil dir t.spans hside c).2 hc,
        (Items.Below_modify_of_not_below item f (hnb c (hne c hc))).1 hb⟩
  have hd := h.disj; rw [hts] at hd
  obtain ⟨htop, hrest⟩ := List.pairwise_cons.1 hd
  have hrestE : ∀ u ∈ rest, ∀ e, u.edges s.g (s.items.modify item f) e ↔ u.edges s.g s.items e :=
    fun u hu e => TEntry.edges_modify_of_not_mem item f hroot (hfree u (by simp [hts, hu])) e
  refine ⟨fun u hu hdep hne' => ?_, List.pairwise_cons.2 ⟨fun u hu e he hue => ?_, ?_⟩⟩
  · rcases List.mem_cons.1 hu with rfl | hu
    · exact absurd hv hne'
    · have hfu := hfree u (by simp [hts, hu])
      refine EntryR.congr (s := s) rfl rfl (fun i hi => ?_) (fun i hi => ?_) (fun i hi e => ?_)
        (h.entries u (by simp [hts, hu]) hdep hne')
      · have hiu : i ≠ item := fun hi' => hfu (hi' ▸ hi)
        show Items.type (s.items.modify item f) i = Items.type s.items i
        simp [Items.type, Array.getElem?_modify, hiu.symm]
      · have hiu : i ≠ item := fun hi' => hfu (hi' ▸ hi)
        exact Items.vs_modify_of_ne item f hiu
      · have hiu : i ≠ item := fun hi' => hfu (hi' ▸ hi)
        exact Items.Below_modify_of_not_below item f (hnb i hiu)
  · rw [hrestE u hu] at hue
    rcases hsub e he with he | he
    · exact htop u hu e he hue
    · obtain ⟨i, hi, hb⟩ := hue
      rw [Items.EdgeBelow, he] at hb
      exact hfree u (by simp [hts, hu]) ((Items.Below.eq_of_no_parent hroot hb) ▸ hi)
  · exact List.Pairwise.imp_of_mem (fun {a b} ha hb hab e he hbe =>
      hab e ((hrestE a ha e).1 he) ((hrestE b hb e).1 hbe)) hrest

/-- Rewriting the terminals of an item that no open entry spans (children and type unchanged). -/
theorem RInvTop.modifyVs_free {v d : Nat} (j : ItemId) (f : Item → Item)
    (hch : ∀ it, (f it).ch = it.ch) (hty : ∀ it, (f it).type = it.type)
    (hfree : ∀ t ∈ s.tstack, j ∉ t.spans.1 ++ t.spans.2)
    (h : s.RInvTop dfs v d) : (after (modifyItem j f) s).RInvTop dfs v d := by
  refine RInvTop.congr (s := s) (s' := after (modifyItem j f) s) rfl rfl rfl
    (fun t ht i hi => ?_) (fun t ht i hi => ?_) (fun t ht i hi e => ?_) h
  · exact Items.type_modify_type_eq j f hty i
  · exact Items.vs_modify_of_ne j f fun hij => hfree t ht (hij ▸ hi)
  · exact Items.Below_modify_ch_eq j f hch

/-- `maybeUnwrapNxt ty` when the entry below the top starts at the exempt vertex `v`: either a fresh
item is allocated (nothing an entry spans changes) or `nxt` is replaced by the children of its
single item (its edge set shrinks). -/
theorem RInvTop.unwrapNxt_exempt {v d : Nat} {ty : NodeType} (hs : Shape s) (hok : UnwrapOk ty s) (hnxt : (nxtE s).vStart = v)
    (h : s.RInvTop dfs v d) : (after (maybeUnwrapNxt ty) s).RInvTop dfs v d := by
  have halloc : ((allocItem ty).run s).2.RInvTop dfs v d := by
    rw [run_allocItem]
    refine RInvTop.congr (s := s) rfl rfl rfl (fun t ht i hi => ?_) (fun t ht i hi => ?_)
      (fun t ht i hi e => ?_) h
    · exact Items.type_push_of_ne _ (Nat.ne_of_lt (hs.span t ht i hi))
    · exact Items.vs_push_of_ne _ (Nat.ne_of_lt (hs.span t ht i hi))
    · exact Items.Below_push_nil _ rfl
  show ((maybeUnwrapNxt ty).run s).2.RInvTop dfs v d
  match hts : s.tstack with
  | [] | [_] => have := hok.two; rw [hts] at this; simp at this
  | a :: b :: rest =>
    have hn : nxtE s = b := by rw [nxtE, hts]; rfl
    have hd : nxtDir s = s.stackDir[b.topDepth]! := by rw [nxtDir, hn]
    have hh : nxtHead s = (getSide b.spans s.stackDir[b.topDepth]!).head! := by rw [nxtHead, hn, hd]
    rw [maybeUnwrapNxt_run_eq ty s a b rest hts _ rfl _ rfl]
    by_cases h1 : ty = .R ∨ s.ternarize = true
    · simp only [h1, ↓reduceIte]; exact halloc
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
      have hch : s.items[i]!.ch = Items.ch s.items i := by rw [Items.ch_eq_getElem hilt, hget]
      set b' : TEntry := { b with spans := setSides dir s.items[i]!.ch [] } with hb'
      have hsub : ∀ e, b'.edges s.g s.items e → b.edges s.g s.items e := by
        rintro e ⟨c, hc, hbe⟩
        rw [hb', mem_setSides, hch] at hc
        exact ⟨i, hib, .head hc hbe⟩
      have hdj := h.disj; rw [hts] at hdj
      obtain ⟨ha, hdj⟩ := List.pairwise_cons.1 hdj
      obtain ⟨hb, hrest⟩ := List.pairwise_cons.1 hdj
      refine ⟨fun t ht hdep hne => ?_, List.pairwise_cons.2 ⟨fun t ht e he hte => ?_,
        List.pairwise_cons.2 ⟨fun t ht e he hte => hb t ht e (hsub e he) hte, hrest⟩⟩⟩
      · rcases List.mem_cons.1 ht with rfl | ht
        · exact EntryR.congr (s := s) rfl rfl (fun _ _ => rfl) (fun _ _ => rfl) (fun _ _ _ => Iff.rfl)
            (h.entries t (by simp [hts]) hdep hne)
        rcases List.mem_cons.1 ht with rfl | ht
        · rw [hn] at hnxt; exact absurd hnxt hne
        · exact EntryR.congr (s := s) rfl rfl (fun _ _ => rfl) (fun _ _ => rfl) (fun _ _ _ => Iff.rfl)
            (h.entries t (by simp [hts, ht]) hdep hne)
      · rcases List.mem_cons.1 ht with rfl | ht
        · exact ha b (by simp) e he (hsub e hte)
        · exact ha t (by simp [ht]) e he hte
    · simp only [h2, ↓reduceIte]; exact halloc

/-- `maybeUnwrapNxt` keeps the stack shape: the entry below the top keeps its start and top. -/
theorem maybeUnwrapNxt_tstack (ty : NodeType) (a b : TEntry) (rest : List TEntry)
    (hts : s.tstack = a :: b :: rest) :
    ∃ b', (after (maybeUnwrapNxt ty) s).tstack = a :: b' :: rest ∧
      b'.vStart = b.vStart ∧ b'.topDepth = b.topDepth := by
  show ∃ b', ((maybeUnwrapNxt ty).run s).2.tstack = a :: b' :: rest ∧ _
  rw [maybeUnwrapNxt_run_eq ty s a b rest hts _ rfl _ rfl]
  by_cases h1 : ty = .R ∨ s.ternarize = true
  · simp only [h1, ↓reduceIte, run_allocItem]; exact ⟨b, hts, rfl, rfl⟩
  simp only [h1, ↓reduceIte]
  by_cases h2 : s.items[(getSide b.spans s.stackDir[b.topDepth]!).head!]!.type = ty
  · simp only [h2, ↓reduceIte]; exact ⟨_, rfl, rfl, rfl⟩
  · simp only [h2, ↓reduceIte, run_allocItem]; exact ⟨b, hts, rfl, rfl⟩

theorem finishTail_true_run (curV d : Nat) (isSingle : Bool) (s : WalkState) :
    (finishTail curV d true isSingle).run s = (true, s) := rfl

/-- The type-1 P-check `finishP curV lowval isType1`: when it fires, the entry below the top starts
at `curV`, so every entry it creates (the unwrapped `nxt`, the merge, the finished P entry) is
exempt from `RInvTop.entries`. -/
theorem RInvTop.finishP {D curV lowval d : Nat} {isType1 : Bool} (hi : s.Inv' D) (hs : Shape s)
    (hok : FinishPOk D curV lowval isType1 s) (h : s.RInvTop dfs curV d) :
    (after (Spqr.finishP curV lowval isType1) s).RInvTop dfs curV d := by
  show ((Spqr.finishP curV lowval isType1).run s).2.RInvTop dfs curV d
  simp only [Spqr.finishP, WalkM.run_bind, run_condP]
  by_cases hc : result (condP curV lowval isType1) s = true
  · have h' : (isType1 && decide (s.tstack.length ≥ 2) && (s.tstack.tail.head!.vStart == curV) &&
        (s.tstack.tail.head!.topDepth == lowval)) = true := hc
    simp only [h', ↓reduceIte, WalkM.run_bind]
    obtain ⟨hu, hc2⟩ := hok.ok hc
    have hnxt : (nxtE s).vStart = curV := by
      simp only [Bool.and_eq_true, beq_iff_eq] at h'
      exact h'.1.2
    have r := maybeUnwrapNxt_spec (v := curV) hi hs (by decide) hu
    have R₁ := h.unwrapNxt_exempt hs hu hnxt
    match hts : s.tstack with
    | [] | [_] => have := hu.two; rw [hts] at this; simp at this
    | a :: b :: rest =>
      obtain ⟨b', hts₁, hb'v, -⟩ := maybeUnwrapNxt_tstack .P a b rest hts
      have hbv : b.vStart = curV := by rw [nxtE, hts] at hnxt; exact hnxt
      set s₁ := ((maybeUnwrapNxt .P).run s).2 with hs₁
      have hts₁' : s₁.tstack = a :: b' :: rest := hts₁
      have R₂ := R₁.mergeTop_exempt a b' rest hts₁' (hb'v.trans hbv)
      have hts₂ : (mergeTstackTops.run s₁).2.tstack = TEntry.mergeInto a b' :: rest := by
        rw [mergeTstackTops_run_eq s₁ a b' rest hts₁']
      have hf := r.free.merge
      have hside := hc2.finish.side
      have hcur : curE (after mergeTstackTops (after (maybeUnwrapNxt .P) s)) = TEntry.mergeInto a b' := by
        show (mergeTstackTops.run s₁).2.tstack.head! = _
        rw [hts₂]; rfl
      rw [hcur] at hside
      exact R₂.finishTop_exempt _ (TEntry.mergeInto a b') rest hts₂ (hb'v.trans hbv) hf.lt hf.root
        hf.free hside
  · have h' : (isType1 && decide (s.tstack.length ≥ 2) && (s.tstack.tail.head!.vStart == curV) &&
        (s.tstack.tail.head!.topDepth == lowval)) = false := Bool.eq_false_iff.2 hc
    simp only [h', Bool.false_eq_true, ↓reduceIte]
    exact h

/-- The back-edge branch of `finishEdge` (`o.cls = .ret lv .backEdge`, `lv < d`) keeps `RInvTop`
at `curV`: the pushed edge entry and the P close start at `curV`, nothing else is touched. The
call sites have `hasVert = true` (the prelude pushes the vertex entry before a type-1 edge) and no
open entry owns the unprocessed edge `o.e` (`Frontier.base_disj` with the whole stack as base). -/
theorem finishEdge_back_rInvTop {D : Nat} (curV d lv : Nat) (o : DfsOut) (origTstack : Nat)
    (ho : o.cls = .ret lv .backEdge) (hlow : lv < d) (hi : s.Inv' D)
    (hs : Shape s) (he : o.e < s.g.ne) (hq : Items.ch s.items (edgeItem s.g o.e) = [])
    (hend : Items.PairEq (curV, s.stackVerts[lv]!) s.g.edges[o.e]!) (hD : lv ≤ D)
    (hrest : FinishRestOk D curV d lv true true true (feBack curV lv d o s))
    (hown : ∀ t ∈ s.tstack, ¬ t.edges s.g s.items o.e) (hR : s.RInvTop dfs curV d) :
    (after (finishEdge curV d o origTstack true) s).RInvTop dfs curV d := by
  have hlv : o.cls.lowval d = lv := by rw [ho]; rfl
  have hge : ¬ (lv ≥ d) := by omega
  have ht' : o.cls.isTree = false := by rw [ho]; rfl
  have h1 : o.cls.isType1 = true := by rw [ho]; rfl
  show ((finishEdge curV d o origTstack true).run s).2.RInvTop dfs curV d
  rw [finishEdge_eq]
  simp only [finishEdge', finishBack, hlv, hge, ↓reduceIte, WalkM.run_bind, WalkM.get_run,
    run_stackDir, run_makeVs, run_modifyItem, ht', Bool.false_eq_true, h1]
  set f : Item → Item := fun it =>
    { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) } with hf
  have hfree₀ : ∀ t ∈ s.tstack, edgeItem s.g o.e ∉ t.spans.1 ++ t.spans.2 := fun t ht hmem =>
    hown t ht ⟨_, hmem, .refl⟩
  have st₀ : Step D curV s (feS₀ d o s) :=
    Step.modifyVs hi hs (edgeItem s.g o.e) _ (by show 1 + s.g.nv + o.e < _; omega)
  have R₀ : (feS₀ d o s).RInvTop dfs curV d :=
    hR.modifyVs_free (edgeItem s.g o.e) f (fun _ => rfl) (fun _ => rfl) hfree₀
  have hq₀ : Items.ch (feS₀ d o s).items (edgeItem (feS₀ d o s).g o.e) = [] := by
    show Items.ch (s.items.modify _ f) (edgeItem s.g o.e) = []
    rw [Items.ch_modify_ch_eq (edgeItem s.g o.e) f (fun _ => rfl)]; exact hq
  have hown₀ : ∀ t ∈ (feS₀ d o s).tstack, ¬ t.edges (feS₀ d o s).g (feS₀ d o s).items o.e := by
    intro t ht hte
    exact hown t ht ((TEntry.edges_congr (fun i _ e =>
      Items.Below_modify_ch_eq (edgeItem s.g o.e) f (fun _ => rfl)) o.e).1 hte)
  have st₁ : Step D curV _ (after (pushEdgeTstack curV lv o.e) (feS₀ d o s)) :=
    Step.pushEdge st₀.inv st₀.shape curV lv o.e he hq₀ hend hD
  have R₁ : (after (pushEdgeTstack curV lv o.e) (feS₀ d o s)).RInvTop dfs curV d :=
    R₀.pushEdge hq₀ hown₀
  have st₂ : Step D curV _ (feBack curV lv d o s) := Step.frame st₁.inv st₁.shape rfl rfl rfl rfl
  have R₂ : (feBack curV lv d o s).RInvTop dfs curV d := R₁.modify _ rfl rfl rfl rfl
  have R₃ := R₂.finishP st₂.inv st₂.shape hrest.p
  exact R₃

end WalkState
end Spqr
