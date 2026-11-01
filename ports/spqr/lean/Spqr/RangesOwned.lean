import Spqr.RangesCloseTree

/-!
# Ownership of the processed prefix (`walk_rootsCover`)

`Owned` is the invariant threaded through the walk of a vertex `v` (depth `d`, subtree starting at
position `n₀` of the schedule `σ`, root tree starting at `nR`, `orig` stack entries at its start,
`sts[k]` the start of the path vertex at depth `k`, `P` the placement predicate): every processed
edge of the root tree is held by a stack entry or hangs under the vertex item of a path vertex; the
vertex items of strict ancestors hold only edges before the next path vertex's start; the entries
above the `orig` old ones hold only edges of `v`'s subtree, the old ones do not bottom at `v`; every
entry bottoms at a visited vertex and unvisited vertex items are empty.
(Checker: `checkOwned`, kinds `own_*`, at every `finishEdge` pre-state, P site and vertex end.)
-/

namespace Spqr
open WalkM
namespace WalkState

variable {σ : List Nat} {n D : Nat} {s : WalkState}

structure Owned (σ : List Nat) (nR n₀ n v d orig : Nat) (P : ItemId → Prop) (sts : List Nat)
    (s : WalkState) : Prop where
  sv : s.stackVerts[d]! = v
  len : orig ≤ s.tstack.length
  cover : ∀ b, nR ≤ b → b < n → (∃ t ∈ s.tstack, t.edges s.g s.items σ[b]!) ∨
    ∃ k, k ≤ d ∧ Items.EdgeBelow s.g s.items (vertItem s.stackVerts[k]!) σ[b]!
  anc : ∀ k, k < d → ∀ e, e < s.g.ne → Items.EdgeBelow s.g s.items (vertItem s.stackVerts[k]!) e →
    σ.idxOf e < sts[k + 1]!
  vertHi : ∀ e, e < s.g.ne → Items.EdgeBelow s.g s.items (vertItem v) e → σ.idxOf e < n
  vertLo : ∀ e, e < s.g.ne → Items.EdgeBelow s.g s.items (vertItem v) e → n₀ ≤ σ.idxOf e
  new : ∀ t ∈ s.tstack.take (s.tstack.length - orig), ∀ e, e < s.g.ne → t.edges s.g s.items e →
    n₀ ≤ σ.idxOf e
  old : ∀ t ∈ s.tstack.drop (s.tstack.length - orig), t.vStart ≠ v
  vis : ∀ t ∈ s.tstack, P (vertItem t.vStart) ∨ ∃ k, k ≤ d ∧ t.vStart = s.stackVerts[k]!
  fresh : ∀ w, w < s.g.nv → ¬ P (vertItem w) → (∀ k, k ≤ d → w ≠ s.stackVerts[k]!) →
    Items.ch s.items (vertItem w) = []

/-! ### Entry edge sets across an unwrap, a close and the type-1 vertex close -/

theorem edges_push_nil' {g : Graph} {items : Items} (x : Item) (hx : x.ch = []) (t : TEntry) (e : Nat) :
    t.edges g (items.push x) e ↔ t.edges g items e :=
  TEntry.edges_congr (fun _ _ _ => by unfold Items.EdgeBelow; exact Items.Below_push_nil x hx) e

theorem edges_retargetEntry {g : Graph} {items : Items} (curV : Nat) (dir : Bool) (t : TEntry) (e : Nat) :
    TEntry.edges g items { t with vStart := curV, spans := setSides (!dir) (t.spans.1 ++ t.spans.2) [] } e ↔
      t.edges g items e := by
  simp only [TEntry.edges, mem_setSides]

/-- `maybeUnwrapNxt` keeps every entry's edge set (the unwrapped `nxt` included). -/
theorem unwrap_edges (hs : Shape s) {ty : NodeType} (hty : ty ∉ [NodeType.F, .V, .Q])
    (hok : UnwrapOk ty s) {a b : TEntry} {rest : List TEntry} (hts : s.tstack = a :: b :: rest) :
    (after (maybeUnwrapNxt ty) s).g = s.g ∧
    ∃ b', (after (maybeUnwrapNxt ty) s).tstack = a :: b' :: rest ∧ b'.vStart = b.vStart ∧
      (∀ e, e < s.g.ne → (b'.edges s.g (after (maybeUnwrapNxt ty) s).items e ↔ b.edges s.g s.items e)) ∧
      ∀ (t : TEntry) e, t.edges s.g (after (maybeUnwrapNxt ty) s).items e ↔ t.edges s.g s.items e := by
  rw [show after (maybeUnwrapNxt ty) s = ((maybeUnwrapNxt ty).run s).2 from rfl,
    maybeUnwrapNxt_run_eq ty s a b rest hts _ rfl _ rfl]
  have halloc : ∀ ty', ((allocItem ty').run s).2.g = s.g ∧
      ∃ b', ((allocItem ty').run s).2.tstack = a :: b' :: rest ∧ b'.vStart = b.vStart ∧
        (∀ e, e < s.g.ne → (b'.edges s.g ((allocItem ty').run s).2.items e ↔ b.edges s.g s.items e)) ∧
        ∀ (t : TEntry) e, t.edges s.g ((allocItem ty').run s).2.items e ↔ t.edges s.g s.items e := by
    intro ty'
    rw [run_allocItem]
    exact ⟨rfl, b, hts, rfl, fun e _ => edges_push_nil' _ rfl b e, fun t e => edges_push_nil' _ rfl t e⟩
  by_cases h1 : ty = .R ∨ s.ternarize = true
  · rw [ite_eq_left_of_eq_true _ _ (eq_true h1)]; exact halloc ty
  rw [ite_eq_right_of_eq_false _ _ (eq_false h1)]
  set dir := s.stackDir[b.topDepth]! with hdir
  set i := (getSide b.spans dir).head! with hi_def
  have hn : nxtE s = b := by rw [nxtE, hts]; rfl
  have hd : nxtDir s = dir := by rw [nxtDir, hn]
  have hh : nxtHead s = i := by rw [nxtHead, hn, hd]
  by_cases h2 : s.items[i]!.type = ty
  · rw [ite_eq_left_of_eq_true _ _ (eq_true h2)]
    have hu : UnwrapAt s := hok.unwrap h1 (by rw [hh]; exact h2)
    have hside := hu.side; rw [hn, hd] at hside
    have hsingle := hu.single; rw [hn, hd, hh] at hsingle
    have hi : i < s.items.size :=
      hs.span b (by rw [hts]; simp) i (mem_of_mem_getSide (by rw [hsingle]; simp))
    have hch : s.items[i]!.ch = Items.ch s.items i := by simp [Items.ch_eq_getElem hi, hi]
    have htyi : s.items[i]!.type = Items.type s.items i := by simp [Items.type_eq_getElem hi, hi]
    have hne : ∀ e, e < s.g.ne → edgeItem s.g e ≠ i := by
      intro e he heq
      have h3 := hs.edge e he
      rw [heq, ← htyi, h2] at h3
      subst h3
      simp at hty
    refine ⟨rfl, _, rfl, rfl, fun e he => ?_, fun t e => Iff.rfl⟩
    rw [hch, TEntry.edges_unwrap dir b.vStart b.topDepth b.firstIdx i hne he,
      TEntry.edges_single dir i hside hsingle]
  · rw [ite_eq_right_of_eq_false _ _ (eq_false h2)]; exact halloc ty

/-- `finishTstackTop x` keeps the top entry's edge set and those of the entries below. -/
theorem finishTop_edges {x : ItemId} (hf : ItemFree s x) {t : TEntry} {rest : List TEntry}
    (hts : s.tstack = t :: rest) (hside : getSide t.spans (!s.stackDir[t.topDepth]!) = []) :
    (after (finishTstackTop x) s).g = s.g ∧
    ∃ t', (after (finishTstackTop x) s).tstack = t' :: rest ∧ t'.vStart = t.vStart ∧
      (∀ e, e < s.g.ne → (t'.edges s.g (after (finishTstackTop x) s).items e ↔ t.edges s.g s.items e)) ∧
      ∀ u ∈ rest, ∀ e, u.edges s.g (after (finishTstackTop x) s).items e ↔ u.edges s.g s.items e := by
  rw [show after (finishTstackTop x) s = ((finishTstackTop x).run s).2 from rfl,
    finishTstackTop_run_eq s x t rest hts]
  set dir := s.stackDir[t.topDepth]! with hdir
  set f : Item → Item := fun it =>
    { it with vs := setSides dir (some s.stackVerts[t.topDepth]!) (some t.vStart),
              ch := getSide t.spans dir } with hf_def
  have hitem_edge : ∀ e, e < s.g.ne → edgeItem s.g e ≠ x := by
    intro e he; show 1 + s.g.nv + e ≠ x; have := hf.node; omega
  have hnb : ∀ i, i ≠ x → ¬ Items.Below s.items i x := fun i hi hb =>
    hi (Items.Below.eq_of_no_parent hf.root hb)
  have hch : Items.ch (s.items.modify x f) x = getSide t.spans dir := by
    rw [Items.ch_modify_at x f hf.lt]
  have hne : ∀ c ∈ getSide t.spans dir, c ≠ x := fun c hc hce => by
    subst hce; exact hf.free t (by simp [hts]) ((mem_of_getSide_nil dir t.spans hside c).2 hc)
  refine ⟨rfl, _, rfl, rfl, fun e he => ?_, fun u hu e => ?_⟩
  · simp only [TEntry.edges, mem_setSides, List.mem_singleton, exists_eq_left, Items.EdgeBelow,
      mem_of_getSide_nil dir t.spans hside]
    constructor
    · intro hb
      rcases hb.head_cases with heq | ⟨c, hc, hb⟩
      · exact absurd heq.symm (hitem_edge e he)
      · simp only [Items.IsParent, hch] at hc
        exact ⟨c, hc, (Items.Below_modify_of_not_below x f (hnb c (hne c hc))).1 hb⟩
    · rintro ⟨c, hc, hb⟩
      exact .head (by simpa [Items.IsParent, hch] using hc)
        ((Items.Below_modify_of_not_below x f (hnb c (hne c hc))).2 hb)
  · exact TEntry.edges_modify_of_not_mem x f hf.root (hf.free u (by simp [hts, hu])) e

/-- The type-1 vertex close from `c :: py :: vy :: base`: one entry bottoming at `curV` that holds
the three edge sets, above `base` with unchanged edge sets. -/
theorem closeVert'_edges {curV : Nat} (hi : s.Inv' D) (hs : Shape s)
    {dir isType1 : Bool} {orig : Nat} {isSingle : Bool} (h1 : isType1 = true)
    (hcv : CloseVertOk D curV dir isType1 orig isSingle s)
    {c py vy : TEntry} {base : List TEntry} (hts : s.tstack = c :: py :: vy :: base) :
    (after (closeVert' curV dir isType1 orig isSingle) s).g = s.g ∧
    ∃ t₃, (after (closeVert' curV dir isType1 orig isSingle) s).tstack = t₃ :: base ∧ t₃.vStart = curV ∧
      (∀ e, e < s.g.ne → (c.edges s.g s.items e ∨ py.edges s.g s.items e ∨ vy.edges s.g s.items e) →
        t₃.edges s.g (after (closeVert' curV dir isType1 orig isSingle) s).items e) ∧
      ∀ u ∈ base, ∀ e, e < s.g.ne →
        (u.edges s.g (after (closeVert' curV dir isType1 orig isSingle) s).items e ↔ u.edges s.g s.items e) := by
  subst h1
  have hS₂ : cvS₂ true orig isSingle s = after (maybeUnwrapNxt (if isSingle then .S else .R)) s := by
    simp only [cvS₂, cvS₁, cvB₁, after, result, vertPre, vertUnwrap, Bool.not_true, Bool.false_eq_true,
      ↓reduceIte, WalkM.pure_run, WalkM.map_run]
  have hS₅ : cvS₅ curV dir true orig isSingle s = after (WalkState.retarget curV dir)
      (after mergeTstackTops (after mergeTstackTops (after (maybeUnwrapNxt (if isSingle then .S else .R)) s))) := by
    simp only [cvS₅, cvS₄, cvS₃, hS₂]
  have hS : after (closeVert' curV dir true orig isSingle) s =
      after (finishTstackTop (result (maybeUnwrapNxt (if isSingle then .S else .R)) s))
        (cvS₅ curV dir true orig isSingle s) := by
    simp only [closeVert', cvS₅, cvS₄, cvS₃, cvS₂, cvS₁, cvB₁, after, result, vertPre, vertUnwrap,
      vertFinish, Bool.not_true, Bool.false_eq_true, ↓reduceIte, WalkM.run_bind, WalkM.pure_run,
      WalkM.map_run]
    rfl
  have hty : (if isSingle then NodeType.S else .R) ∉ [NodeType.F, .V, .Q] := by cases isSingle <;> decide
  have hu : UnwrapOk (if isSingle then .S else .R) s := hcv.unwrap rfl
  have hm₁ := hcv.merge₁; rw [hS₂] at hm₁
  have hm₂ := hcv.merge₂; simp only [cvS₃, hS₂] at hm₂
  have hrt := hcv.retarget; simp only [cvS₄, cvS₃, hS₂] at hrt
  have hfin := hcv.finish rfl; rw [hS₅] at hfin
  rw [hS, hS₅]
  generalize (if isSingle then NodeType.S else NodeType.R) = ty at hty hu hm₁ hm₂ hrt hfin ⊢
  obtain ⟨hg₂, py', hts₂, hvs₂, hpy', hall₂⟩ := unwrap_edges hs hty hu hts
  have r := maybeUnwrapNxt_spec (v := curV) hi hs hty hu
  set s₂ := after (maybeUnwrapNxt ty) s with hs₂def
  have st₂ : Step D curV s s₂ := r.step
  have hfree₂ : ItemFree s₂ (result (maybeUnwrapNxt ty) s) := r.free
  have hS₃ := after_mergeTstackTops_eq hts₂
  have st₃ : Step D curV s₂ (after mergeTstackTops s₂) := Step.mergeTop st₂.inv st₂.shape hm₁
  have hts₃ : (after mergeTstackTops s₂).tstack = TEntry.mergeInto c py' :: vy :: base := by rw [hS₃]
  have hS₄ := after_mergeTstackTops_eq hts₃
  have st₄ : Step D curV _ (after mergeTstackTops (after mergeTstackTops s₂)) :=
    Step.mergeTop st₃.inv st₃.shape hm₂
  set s₄ := after mergeTstackTops (after mergeTstackTops s₂) with hs₄def
  set m := TEntry.mergeInto (TEntry.mergeInto c py') vy with hm_def
  have hts₄ : s₄.tstack = m :: base := by rw [hS₄]
  have hit₄ : s₄.items = s₂.items := by rw [hS₄, hS₃]
  have st₅ : Step D curV s₄ (after (WalkState.retarget curV dir) s₄) :=
    Step.retarget st₄.inv st₄.shape curV dir hrt
  have hS₅' : after (WalkState.retarget curV dir) s₄ =
      { s₄ with tstack := { m with vStart := curV, spans := setSides (!dir) (m.spans.1 ++ m.spans.2) [] } :: base } := by
    show ((WalkState.retarget curV dir).run s₄).2 = _; rw [retarget_run_eq curV dir s₄ m base hts₄]
  set t₅ : TEntry := { m with vStart := curV, spans := setSides (!dir) (m.spans.1 ++ m.spans.2) [] } with ht₅
  set s₅ := after (WalkState.retarget curV dir) s₄ with hs₅def
  have hts₅ : s₅.tstack = t₅ :: base := by rw [hS₅']
  have hit₅ : s₅.items = s₂.items := by rw [hS₅', hit₄]
  have hg₅ : s₅.g = s.g := by rw [st₅.g, st₄.g, st₃.g, st₂.g]
  have hfree₅ : ItemFree s₅ (result (maybeUnwrapNxt ty) s) := (hfree₂.merge.merge).retarget curV dir
  have hside : getSide t₅.spans (!s₅.stackDir[t₅.topDepth]!) = [] := by
    have := hfin.side; simp only [curE, hts₅, List.head!_cons] at this; exact this
  obtain ⟨hg₆, t₆, hts₆, hvs₆, ht₆, hbase₆⟩ := finishTop_edges hfree₅ hts₅ hside
  refine ⟨by rw [hg₆, hg₅], t₆, hts₆, by rw [hvs₆], fun e he hce => ?_, fun u hu e he => ?_⟩
  · have h5 : t₅.edges s₅.g s₅.items e := by
      rw [hit₅, ht₅, edges_retargetEntry, hm_def, TEntry.edges_mergeInto, TEntry.edges_mergeInto, hg₅]
      rcases hce with hc | hp | hv'
      · exact .inl (.inl ((hall₂ c e).2 hc))
      · exact .inl (.inr ((hpy' e he).2 hp))
      · exact .inr ((hall₂ vy e).2 hv')
    have := (ht₆ e (by rw [hg₅]; exact he)).2 h5
    rwa [hg₅] at this
  · have h1 := hbase₆ u hu e
    rw [hit₅, hg₅] at h1
    rw [h1]; exact hall₂ u e

/-- The exports `finishP_ownership`/`finishEdge_ownedD` consume at a `finishEdge` site. -/
structure OwnSite (σ : List Nat) (n D curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool)
    (s : WalkState) : Prop where
  hD : D = if o.cls.isTree then d + 1 else d
  nodup : σ.Nodup
  rgs : RgS σ n D s
  guards : FinishGuards d o origTstack hasVert s
  book : FinishBook curV d o origTstack hasVert s
  pos : σ[n]? = some o.e

/-- P-site coverage at a `finishEdge` call: with the prefix owned (`Owned`), the merge base `nxt`
(bottoming at `curV`, hence above the old entries) holds only edges of `curV`'s subtree, and the
subtree's processed edges `n₀..n` are all on the stack, so `MergeBaseCover σ (min n₀ (n + 1)) (n + 1)`
holds (checker: `checkOwned`, kinds `own_*`). -/
theorem finishP_ownership {curV d : Nat} {o : DfsOut} {origTstack : Nat} {hasVert : Bool}
    (h : OwnSite σ n D curV d o origTstack hasVert s)
    {nR n₀ orig : Nat} {P : ItemId → Prop} {sts : List Nat}
    (ho : Owned σ nR n₀ n curV d orig P sts s) (hnR : nR ≤ n₀)
    (hsts : ∀ k, k ≤ d → sts[k]! ≤ n₀)
    (hhv : o.cls.lowval d < d → o.cls.isType1 = true → hasVert = true) :
    FinishPOwnership σ n curV d o origTstack hasVert s := by
  obtain ⟨sub, base, hlen, hE⟩ := h.book.ear
  have hn : n < σ.length := (List.getElem?_eq_some_iff.1 h.pos).1
  have hσn : σ[n]! = o.e := by
    rw [getElem!_pos σ n hn]; exact (List.getElem?_eq_some_iff.1 h.pos).2
  have hj : edgeItem s.g o.e < 1 + s.g.nv + s.g.ne := by
    show 1 + s.g.nv + o.e < _; have := h.book.e_lt; omega
  have st₀ : Step D curV s (feS₀ d o s) := Step.modifyVs h.rgs.1.inv h.rgs.2.1 (edgeItem s.g o.e) _ hj
  have hv₀ : curV < (feS₀ d o s).g.nv := by rw [st₀.g]; exact h.book.v_lt
  have hbase_mem : ∀ t ∈ base, t ∈ s.tstack := fun t ht => by
    rw [hE.tstack]; exact List.mem_append_right _ ht
  have hnew : ∀ b ∈ s.tstack, b.vStart = curV → b ∈ s.tstack.take (s.tstack.length - orig) := by
    intro b hb hbv
    rw [← List.take_append_drop (s.tstack.length - orig) s.tstack] at hb
    rcases List.mem_append.1 hb with hb | hb
    · exact hb
    · exact absurd hbv (ho.old b hb)
  refine ⟨fun hlow ht hv hp => ?_, fun hlow hb hp => ?_⟩
  · obtain ⟨lv, kind, hoc, hl⟩ := ret_of_lowval_lt hlow
    have hlv : o.cls.lowval d = lv := by rw [hoc]; rfl
    have hok := finishOk_of_guards hoc hl h.guards hE hlen h.rgs.1.inv h.rgs.2.1 h.hD h.book.v_lt
      h.book.e_lt h.book.q (h.book.ends lv kind hoc) h.book.vert
    have hp' := hp
    simp only [result, run_condP, Bool.and_eq_true, decide_eq_true_eq, beq_iff_eq] at hp'
    obtain ⟨⟨⟨ht1, hlen2⟩, hnv⟩, -⟩ := hp'
    have st₁ : Step D curV _ (feS₁ d o s) := Step.closeEars st₀.inv st₀.shape hv₀ (hok.ears ht)
    have hv₁ : curV < (feS₁ d o s).g.nv := by rw [st₁.g]; exact hv₀
    have st₂ : Step D curV _ (feS₂ d o s) := Step.mergeLate st₁.inv st₁.shape hv₁ (hok.late ht)
    have hv₂ : curV < (feS₂ d o s).g.nv := by rw [st₂.g]; exact hv₁
    have hg₂ : (feS₂ d o s).g = s.g := by rw [st₂.g, st₁.g, st₀.g]
    obtain ⟨c, mid, py, vy, hC⟩ := hE.close ht hlow
    have hmid := (hC.type1 ht1).1
    have hts₂ : (feS₂ d o s).tstack = c :: py :: vy :: base := by rw [hC.tstack, hmid]; rfl
    have hcv := hok.vert ht hv
    obtain ⟨hg₃, t₃, hts₃, hvs₃, hsub₃, hbase₃⟩ := closeVert'_edges st₂.inv st₂.shape ht1 hcv hts₂
    have hpo : FinishPOk D curV lv o.cls.isType1 (feS₃ curV d o origTstack s) := (hok.rest_vert ht hv).p
    have st₈ : Step D curV (feS₂ d o s) (feS₃ curV d o origTstack s) :=
      Step.closeVert' st₂.inv st₂.shape hv₂ hcv
    have hS₃ : feS₃ curV d o origTstack s =
        after (closeVert' curV s.stackDir[d]! o.cls.isType1 origTstack (feSingle d o s)) (feS₂ d o s) := rfl
    rw [hS₃] at hp hlen2 hnv hpo st₈ ⊢
    set S := after (closeVert' curV s.stackDir[d]! o.cls.isType1 origTstack (feSingle d o s)) (feS₂ d o s)
      with hSdef
    obtain ⟨b, rest, rfl⟩ : ∃ b rest, base = b :: rest := by
      rcases base with _ | ⟨b, rest⟩
      · rw [hts₃] at hlen2; simp at hlen2
      · exact ⟨b, rest, rfl⟩
    have hnv' : b.vStart = curV := by rw [hts₃] at hnv; exact hnv
    have hu : UnwrapOk .P S := (hpo.ok (by rw [hlv] at hp; exact hp)).1
    obtain ⟨hgP, b', htsP, hvsP, hb', hallP⟩ := unwrap_edges st₈.shape (by decide) hu hts₃
    have hgS : S.g = s.g := by rw [st₈.g, hg₂]
    have hgP' : (after (maybeUnwrapNxt .P) S).g = s.g := by rw [hgP, hgS]
    have hbnew := hnew b (hbase_mem b (by simp)) hnv'
    have hbaseS : ∀ t ∈ b :: rest, ∀ e, e < s.g.ne → (t.edges s.g S.items e ↔ t.edges s.g s.items e) := by
      intro t ht e he
      have h3 := hbase₃ t ht e (by rw [hg₂]; exact he)
      rw [hg₂] at h3
      exact h3.trans (hC.base_edges t ht e)
    have hbaseP : ∀ t ∈ b :: rest, ∀ e, e < s.g.ne → t.edges s.g s.items e →
        ∃ t' ∈ (after (maybeUnwrapNxt .P) S).tstack, t'.edges s.g (after (maybeUnwrapNxt .P) S).items e := by
      intro t ht e he hte
      have h3 : t.edges s.g S.items e := (hbaseS t ht e he).2 hte
      rcases List.mem_cons.1 ht with rfl | ht
      · refine ⟨b', by rw [htsP]; simp, ?_⟩
        have := hb' e (by rw [hgS]; exact he); rw [hgS] at this; exact this.2 h3
      · refine ⟨t, by rw [htsP]; simp [ht], ?_⟩
        have := hallP t e; rw [hgS] at this; exact this.2 h3
    have hsubP : ∀ e, e < s.g.ne → subEdges o e →
        ∃ t' ∈ (after (maybeUnwrapNxt .P) S).tstack, t'.edges s.g (after (maybeUnwrapNxt .P) S).items e := by
      intro e he hsub
      obtain ⟨t, ht, hte⟩ := hC.sub_cover e he hsub
      rw [hmid] at ht
      have h3 : t₃.edges s.g S.items e := by
        have := hsub₃ e (by rw [hg₂]; exact he)
        rw [hg₂] at this
        apply this
        simp at ht
        rcases ht with rfl | rfl | rfl
        · exact .inl hte
        · exact .inr (.inl hte)
        · exact .inr (.inr hte)
      refine ⟨t₃, by rw [htsP]; simp, ?_⟩
      have := hallP t₃ e; rw [hgS] at this; exact this.2 h3
    have hcover : ∀ bb, n₀ ≤ bb → bb < n →
        ∃ t' ∈ (after (maybeUnwrapNxt .P) S).tstack, t'.edges s.g (after (maybeUnwrapNxt .P) S).items σ[bb]! := by
      intro bb hlo hbn
      have hbl : bb < σ.length := by omega
      have he : σ[bb]! < s.g.ne := getElem!_lt h.rgs.2.2 hbl
      rcases ho.cover bb (by omega) hbn with ⟨t, ht, hte⟩ | ⟨k, hk, hbel⟩
      · rw [hE.tstack] at ht
        rcases List.mem_append.1 ht with ht | ht
        · exact hsubP _ he (hE.sub_edges t ht _ he hte)
        · exact hbaseP t ht _ he hte
      · rcases Nat.lt_or_ge k d with hkd | hkd
        · exfalso
          have h1 := ho.anc k hkd _ he hbel
          have h2 := hsts (k + 1) (by omega)
          rw [idxOf_getElem! h.nodup hbl] at h1
          omega
        · have hkd' : k = d := by omega
          subst hkd'
          rw [ho.sv] at hbel
          obtain ⟨t, ht, -, hvi⟩ := hE.vert hv
          exact hbaseP t ht _ he ⟨vertItem curV, hvi, hbel⟩
    refine ⟨min n₀ (n + 1), ⟨?_, ?_⟩⟩
    · intro cur nxt rest' hts' a bb hab hbb hpiece
      exfalso
      rw [htsP] at hts'
      obtain ⟨rfl, rfl, rfl⟩ : t₃ = cur ∧ b' = nxt ∧ rest = rest' := by simpa using hts'
      have hal : a < σ.length := by omega
      have he : σ[a]! < s.g.ne := getElem!_lt h.rgs.2.2 hal
      have h1 : b'.edges s.g (after (maybeUnwrapNxt .P) S).items σ[a]! := by
        have := hpiece.edges; rwa [hgP'] at this
      have h2 : b.edges s.g s.items σ[a]! := by
        have h3 := hb' _ (by rw [hgS]; exact he)
        rw [hgS] at h3
        exact (hbaseS b (by simp) _ he).1 (h3.1 h1)
      have := ho.new b hbnew _ he h2
      rw [idxOf_getElem! h.nodup hal] at this
      omega
    · intro bb hlo hbn
      rw [hgP']
      rcases Nat.lt_or_ge bb n with hbn' | hbn'
      · exact hcover bb (by omega) hbn'
      · have hbe : bb = n := by omega
        subst hbe
        rw [hσn]
        exact hsubP o.e h.book.e_lt (.inl rfl)
  · obtain ⟨lv, kind, hoc, hl⟩ := ret_of_lowval_lt hlow
    have hlv : o.cls.lowval d = lv := by rw [hoc]; rfl
    have hok := finishOk_of_guards hoc hl h.guards hE hlen h.rgs.1.inv h.rgs.2.1 h.hD h.book.v_lt
      h.book.e_lt h.book.q (h.book.ends lv kind hoc) h.book.vert
    rw [hlv] at hp ⊢
    have hp' := hp
    simp only [result, run_condP, Bool.and_eq_true, decide_eq_true_eq, beq_iff_eq] at hp'
    obtain ⟨⟨⟨ht1, hlen2⟩, hnv⟩, -⟩ := hp'
    have hvt : hasVert = true := hhv hlow ht1
    have hq : Items.ch (feS₀ d o s).items (edgeItem (feS₀ d o s).g o.e) = [] := by
      show Items.ch (s.items.modify _ _) (edgeItem s.g o.e) = []
      rw [Items.ch_modify_ch_eq (edgeItem s.g o.e)
        (fun it => { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) }) (fun _ => rfl)]
      exact hok.q hb
    have st₁ : Step D curV _ (after (pushEdgeTstack curV lv o.e) (feS₀ d o s)) :=
      Step.pushEdge st₀.inv st₀.shape curV lv o.e hok.e_lt hq (hok.ends hb) (hok.lv_le hb)
    have st₂ : Step D curV _ (feBack curV lv d o s) := Step.frame st₁.inv st₁.shape rfl rfl rfl rfl
    have hpo : FinishPOk D curV lv o.cls.isType1 (feBack curV lv d o s) := (hok.rest_back hb).p
    set S := feBack curV lv d o s with hSdef
    have htsS : S.tstack =
        ⟨curV, lv, s.nxtEdgeIdx, setSides s.stackDir[lv]! [edgeItem s.g o.e] []⟩ :: s.tstack := rfl
    have hitS : S.items = s.items.modify (edgeItem s.g o.e)
        (fun it => { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) }) := rfl
    have hgS : S.g = s.g := rfl
    obtain ⟨b, rest, hts⟩ : ∃ b rest, s.tstack = b :: rest := by
      rcases hts0 : s.tstack with _ | ⟨b, rest⟩
      · rw [htsS, hts0] at hlen2; simp at hlen2
      · exact ⟨b, rest, rfl⟩
    have hnv' : b.vStart = curV := by rw [htsS, hts] at hnv; exact hnv
    have hu : UnwrapOk .P S := (hpo.ok hp).1
    have htsS' : S.tstack =
        ⟨curV, lv, s.nxtEdgeIdx, setSides s.stackDir[lv]! [edgeItem s.g o.e] []⟩ :: b :: rest := by
      rw [htsS, hts]
    obtain ⟨hgP, b', htsP, hvsP, hb', hallP⟩ := unwrap_edges st₂.shape (by decide) hu htsS'
    have hgP' : (after (maybeUnwrapNxt .P) S).g = s.g := by rw [hgP, hgS]
    have hitemsS : ∀ (t : TEntry) e, t.edges s.g S.items e ↔ t.edges s.g s.items e := fun t e => by
      rw [hitS]
      exact TEntry.edges_congr (fun _ _ _ => by
        unfold Items.EdgeBelow
        exact Items.Below_modify_ch_eq (edgeItem s.g o.e)
          (fun it => { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) })
          (fun _ => rfl)) e
    have hbnew := hnew b (by rw [hts]; simp) hnv'
    have hbaseP : ∀ t ∈ b :: rest, ∀ e, e < s.g.ne → t.edges s.g s.items e →
        ∃ t' ∈ (after (maybeUnwrapNxt .P) S).tstack, t'.edges s.g (after (maybeUnwrapNxt .P) S).items e := by
      intro t ht e he hte
      have h3 : t.edges s.g S.items e := (hitemsS t e).2 hte
      rcases List.mem_cons.1 ht with rfl | ht
      · refine ⟨b', by rw [htsP]; simp, ?_⟩
        have := hb' e (by rw [hgS]; exact he); rw [hgS] at this; exact this.2 h3
      · refine ⟨t, by rw [htsP]; simp [ht], ?_⟩
        have := hallP t e; rw [hgS] at this; exact this.2 h3
    have hcover : ∀ bb, n₀ ≤ bb → bb < n →
        ∃ t' ∈ (after (maybeUnwrapNxt .P) S).tstack, t'.edges s.g (after (maybeUnwrapNxt .P) S).items σ[bb]! := by
      intro bb hlo hbn
      have hbl : bb < σ.length := by omega
      have he : σ[bb]! < s.g.ne := getElem!_lt h.rgs.2.2 hbl
      rcases ho.cover bb (by omega) hbn with ⟨t, ht, hte⟩ | ⟨k, hk, hbel⟩
      · rw [hts] at ht
        exact hbaseP t ht _ he hte
      · rcases Nat.lt_or_ge k d with hkd | hkd
        · exfalso
          have h1 := ho.anc k hkd _ he hbel
          have h2 := hsts (k + 1) (by omega)
          rw [idxOf_getElem! h.nodup hbl] at h1
          omega
        · have hkd' : k = d := by omega
          subst hkd'
          rw [ho.sv] at hbel
          obtain ⟨t, ht, -, hvi⟩ := hE.vert hvt
          have ht' := hbase_mem t ht
          rw [hts] at ht'
          exact hbaseP t ht' _ he ⟨vertItem curV, hvi, hbel⟩
    refine ⟨min n₀ (n + 1), ⟨?_, ?_⟩⟩
    · intro cur nxt rest' hts' a bb hab hbb hpiece
      exfalso
      rw [htsP] at hts'
      obtain ⟨-, rfl, rfl⟩ : (⟨curV, lv, s.nxtEdgeIdx, setSides s.stackDir[lv]! [edgeItem s.g o.e] []⟩ : TEntry) = cur ∧
          b' = nxt ∧ rest = rest' := by simpa using hts'
      have hal : a < σ.length := by omega
      have he : σ[a]! < s.g.ne := getElem!_lt h.rgs.2.2 hal
      have h1 : b'.edges s.g (after (maybeUnwrapNxt .P) S).items σ[a]! := by
        have := hpiece.edges; rwa [hgP'] at this
      have h2 : b.edges s.g s.items σ[a]! := by
        have h3 := hb' _ (by rw [hgS]; exact he)
        rw [hgS] at h3
        exact (hitemsS b _).1 (h3.1 h1)
      have := ho.new b hbnew _ he h2
      rw [idxOf_getElem! h.nodup hal] at this
      omega
    · intro bb hlo hbn
      rw [hgP']
      rcases Nat.lt_or_ge bb n with hbn' | hbn'
      · exact hcover bb (by omega) hbn'
      · have hbe : bb = n := by omega
        subst hbe
        rw [hσn]
        refine ⟨⟨curV, lv, s.nxtEdgeIdx, setSides s.stackDir[lv]! [edgeItem s.g o.e] []⟩, by rw [htsP]; simp, ?_⟩
        have hq' : Items.ch S.items (edgeItem s.g o.e) = [] := hq
        have := hallP ⟨curV, lv, s.nxtEdgeIdx, setSides s.stackDir[lv]! [edgeItem s.g o.e] []⟩ o.e
        rw [hgS] at this
        exact this.2 ((TEntry.edges_edgeEntry s.stackDir[lv]! curV lv s.nxtEdgeIdx o.e hq' o.e).2 rfl)

/-! ### Depth-indexed ownership (`walk_rootsCover`'s induction) -/

/-- `Owned` at every path vertex at once: `sts[k]`/`origs[k]` are the schedule position and stack
length at the entry of the path vertex of depth `k` (`k ≤ d`), `n` the processed prefix.
`cover` ranges over the root tree (`sts[0]`); the vertex item of the current vertex `stackVerts[d]`
holds exactly `sts[d]..n`, that of a strict ancestor `k` holds `sts[k]..sts[k+1]`.
(Checker: `checkOwned`, kinds `own_*`.) -/
structure OwnedD (σ sts origs : List Nat) (P : ItemId → Prop) (d n : Nat) (s : WalkState) : Prop where
  len : ∀ k, k ≤ d → origs[k]! ≤ s.tstack.length
  cover : ∀ b, sts[0]! ≤ b → b < n → (∃ t ∈ s.tstack, t.edges s.g s.items σ[b]!) ∨
    ∃ k, k ≤ d ∧ Items.EdgeBelow s.g s.items (vertItem s.stackVerts[k]!) σ[b]!
  anc : ∀ k, k < d → ∀ e, e < s.g.ne → Items.EdgeBelow s.g s.items (vertItem s.stackVerts[k]!) e →
    σ.idxOf e < sts[k + 1]!
  hi : ∀ e, e < s.g.ne → Items.EdgeBelow s.g s.items (vertItem s.stackVerts[d]!) e → σ.idxOf e < n
  lo : ∀ k, k ≤ d → ∀ e, e < s.g.ne → Items.EdgeBelow s.g s.items (vertItem s.stackVerts[k]!) e →
    sts[k]! ≤ σ.idxOf e
  new : ∀ k, k ≤ d → ∀ t ∈ s.tstack.take (s.tstack.length - origs[k]!), ∀ e, e < s.g.ne →
    t.edges s.g s.items e → sts[k]! ≤ σ.idxOf e
  old : ∀ k, k ≤ d → ∀ t ∈ s.tstack.drop (s.tstack.length - origs[k]!), t.vStart ≠ s.stackVerts[k]!
  vis : ∀ t ∈ s.tstack, P (vertItem t.vStart) ∨ ∃ k, k ≤ d ∧ t.vStart = s.stackVerts[k]!
  fresh : ∀ w, w < s.g.nv → ¬ P (vertItem w) → (∀ k, k ≤ d → w ≠ s.stackVerts[k]!) →
    Items.ch s.items (vertItem w) = []

theorem CloseBase.toOwnSite {curV d : Nat} {o : DfsOut} {origTstack : Nat} {hasVert : Bool}
    (h : CloseBase σ n D curV d o origTstack hasVert s) : OwnSite σ n D curV d o origTstack hasVert s :=
  ⟨h.hD, h.nodup, h.rgs, h.guards, h.book, h.site.pos⟩

/-- Every edge under `vertItem v` is held by an open entry. -/
def VertCover (v : Nat) (s : WalkState) : Prop :=
  ∀ e, e < s.g.ne → Items.EdgeBelow s.g s.items (vertItem v) e → ∃ t ∈ s.tstack, t.edges s.g s.items e

/-- `vertItem v` is spanned by an open entry or lies under one of its spanned items; implies
`VertCover v` and is kept by every primitive of `finishEdge`. -/
def VSpan (v : Nat) (s : WalkState) : Prop :=
  ∃ t ∈ s.tstack, ∃ i ∈ t.spans.1 ++ t.spans.2, Items.Below s.items i (vertItem v)

theorem VSpan.vertCover {v : Nat} (h : VSpan v s) : VertCover v s := by
  obtain ⟨t, ht, i, hi, hb⟩ := h
  intro e _ he
  exact ⟨t, ht, i, hi, Relation.ReflTransGen.trans hb he⟩

theorem VSpan.frame {v : Nat} {s s' : WalkState} (h : VSpan v s)
    (hch : ∀ p, Items.ch s'.items p = Items.ch s.items p)
    (hsp : ∀ t ∈ s.tstack, ∀ i ∈ t.spans.1 ++ t.spans.2,
      ∃ t' ∈ s'.tstack, i ∈ t'.spans.1 ++ t'.spans.2) :
    VSpan v s' := by
  obtain ⟨t, ht, i, hi, hb⟩ := h
  obtain ⟨t', ht', hi'⟩ := hsp t ht i hi
  exact ⟨t', ht', i, hi', (Items.Below_congr hch).2 hb⟩

theorem VSpan.frame' {v : Nat} {s s' : WalkState} (h : VSpan v s) (hi : s'.items = s.items)
    (hts : s'.tstack = s.tstack) : VSpan v s' :=
  h.frame (by rw [hi]; exact fun _ => rfl) (fun t ht i hi' => ⟨t, by rw [hts]; exact ht, hi'⟩)

theorem VSpan.ne_nil {v : Nat} (h : VSpan v s) : s.tstack ≠ [] := by
  obtain ⟨t, ht, -⟩ := h
  exact List.ne_nil_of_mem ht

theorem mem_spans_mergeTop {b a : TEntry} {rest : List TEntry} {t : TEntry} (ht : t ∈ b :: a :: rest)
    {i : ItemId} (hi : i ∈ t.spans.1 ++ t.spans.2) :
    ∃ t' ∈ WalkM.mergeTop (b :: a :: rest), i ∈ t'.spans.1 ++ t'.spans.2 := by
  simp only [WalkM.mergeTop]
  rcases List.mem_cons.1 ht with rfl | ht
  · exact ⟨_, List.mem_cons_self .., by simp only [List.mem_append] at hi ⊢; tauto⟩
  rcases List.mem_cons.1 ht with rfl | ht
  · exact ⟨_, List.mem_cons_self .., by simp only [List.mem_append] at hi ⊢; tauto⟩
  · exact ⟨t, List.mem_cons_of_mem _ ht, hi⟩

theorem two_of_mergeTop_ne_nil {l : List TEntry} (h : WalkM.mergeTop l ≠ []) :
    ∃ a b rest, l = a :: b :: rest := by
  rcases l with _ | ⟨a, _ | ⟨b, rest⟩⟩
  · exact absurd rfl h
  · exact absurd rfl h
  · exact ⟨a, b, rest, rfl⟩

theorem VSpan.mergeTop {v : Nat} (h : VSpan v s) {a b : TEntry} {rest : List TEntry}
    (hts : s.tstack = a :: b :: rest) : VSpan v (after mergeTstackTops s) := by
  rw [after_mergeTstackTops]
  refine h.frame (fun _ => rfl) (fun t ht i hi => ?_)
  show ∃ t' ∈ WalkM.mergeTop s.tstack, _
  rw [hts] at ht ⊢
  exact mem_spans_mergeTop ht hi

/-- A merge whose result is nonempty had two entries. -/
theorem VSpan.mergeTop' {v : Nat} (h : VSpan v s) (hne : (after mergeTstackTops s).tstack ≠ []) :
    VSpan v (after mergeTstackTops s) := by
  rw [after_mergeTstackTops] at hne
  obtain ⟨a, b, rest, hts⟩ := two_of_mergeTop_ne_nil hne
  exact h.mergeTop hts

theorem VSpan.push {v : Nat} (h : VSpan v s) (w d : Nat) (i : ItemId) :
    VSpan v (after (pushTstack w d i) s) :=
  h.frame (fun _ => rfl) (fun t ht _ hj => ⟨t, List.mem_cons_of_mem _ ht, hj⟩)

theorem VSpan.pushEdge {v : Nat} (h : VSpan v s) (w d e : Nat) :
    VSpan v (after (pushEdgeTstack w d e) s) :=
  h.frame (fun _ => rfl) (fun t ht _ hj => ⟨t, List.mem_cons_of_mem _ ht, hj⟩)

theorem VSpan.pushVert (v d : Nat) (s : WalkState) : VSpan v (after (pushVertTstack v d) s) :=
  ⟨_, List.mem_cons_self .., vertItem v, (mem_setSides _ _ _).2 (List.mem_singleton.2 rfl),
    Relation.ReflTransGen.refl⟩

theorem VSpan.alloc {v : Nat} (h : VSpan v s) (ty : NodeType) : VSpan v (after (allocItem ty) s) := by
  show VSpan v ((allocItem ty).run s).2
  rw [run_allocItem]
  exact h.frame (fun p => Items.ch_push_nil _ rfl p) (fun t ht i hi => ⟨t, ht, hi⟩)

theorem VSpan.modifyVs {v : Nat} (h : VSpan v s) (j : ItemId) (vs : Option Nat × Option Nat) :
    VSpan v (after (modifyItem j fun it => { it with vs := vs }) s) := by
  show VSpan v ((modifyItem j fun it => { it with vs := vs }).run s).2
  rw [run_modifyItem]
  exact h.frame (fun p => by
    show Items.ch (s.items.modify j fun it => { it with vs := vs }) p = _
    exact Items.ch_modify_ch_eq (items := s.items) j (fun it => { it with vs := vs }) (fun _ => rfl) p) (fun t ht i hi => ⟨t, ht, hi⟩)

theorem VSpan.retarget {v : Nat} (h : VSpan v s) (curV : Nat) (dir : Bool) :
    VSpan v (after (WalkState.retarget curV dir) s) := by
  show VSpan v ((WalkState.retarget curV dir).run s).2
  rcases hts : s.tstack with _ | ⟨t, rest⟩
  · exact absurd hts h.ne_nil
  rw [retarget_run_eq curV dir s t rest hts]
  refine h.frame (fun _ => rfl) (fun u hu i hi => ?_)
  rw [hts] at hu
  rcases List.mem_cons.1 hu with rfl | hu
  · exact ⟨_, List.mem_cons_self .., by
      show i ∈ (setSides (!dir) (u.spans.1 ++ u.spans.2) []).1 ++ (setSides (!dir) (u.spans.1 ++ u.spans.2) []).2
      rw [mem_setSides]; exact hi⟩
  · exact ⟨u, List.mem_cons_of_mem _ hu, hi⟩

theorem VSpan.maybeUnwrap {v : Nat} (h : VSpan v s) (hs : Shape s) (hv : v < s.g.nv) {ty : NodeType}
    (hty : ty ≠ .V) (hok : UnwrapOk ty s) {a b : TEntry} {rest : List TEntry}
    (hts : s.tstack = a :: b :: rest) : VSpan v (after (maybeUnwrapNxt ty) s) := by
  show VSpan v ((maybeUnwrapNxt ty).run s).2
  rw [maybeUnwrapNxt_run_eq ty s a b rest hts _ rfl _ rfl]
  have halloc : VSpan v ((allocItem ty).run s).2 := h.alloc ty
  by_cases h1 : ty = .R ∨ s.ternarize = true
  · simp only [h1, ↓reduceIte]; exact halloc
  simp only [h1, ↓reduceIte]
  have hn : nxtE s = b := by rw [nxtE, hts]; rfl
  have hd : nxtDir s = s.stackDir[b.topDepth]! := by rw [nxtDir, hn]
  have hh : nxtHead s = (getSide b.spans s.stackDir[b.topDepth]!).head! := by rw [nxtHead, hn, hd]
  by_cases h2 : s.items[(getSide b.spans s.stackDir[b.topDepth]!).head!]!.type = ty
  · simp only [h2, ↓reduceIte]
    have hu : UnwrapAt s := hok.unwrap h1 (by rw [hh]; exact h2)
    have hside := hu.side; rw [hn, hd] at hside
    have hsingle := hu.single; rw [hn, hd, hh] at hsingle
    generalize hi_def : (getSide b.spans s.stackDir[b.topDepth]!).head! = i at h2 hsingle ⊢
    have hi : i < s.items.size :=
      hs.span b (by rw [hts]; simp) i (mem_of_mem_getSide (by rw [hsingle]; simp))
    have hch : s.items[i]!.ch = Items.ch s.items i := by simp [Items.ch_eq_getElem hi, hi]
    have htyi : s.items[i]!.type = Items.type s.items i := by simp [Items.type_eq_getElem hi, hi]
    obtain ⟨t, ht, j, hj, hb⟩ := h
    rw [hts] at ht
    rcases List.mem_cons.1 ht with ht' | ht
    · rw [ht'] at hj; exact ⟨a, List.mem_cons_self .., j, hj, hb⟩
    rcases List.mem_cons.1 ht with ht' | ht
    · rw [ht'] at hj
      have hji : j = i := by
        have := (mem_of_getSide_nil _ b.spans hside j).1 hj
        rw [hsingle] at this; exact List.mem_singleton.1 this
      rw [hji] at hb
      rcases hb.head_cases with heq | ⟨c, hc, hb⟩
      · exfalso
        have := hs.vert v hv
        rw [← heq, ← htyi, h2] at this
        exact hty this
      · refine ⟨_, List.mem_cons_of_mem _ (List.mem_cons_self ..), c, ?_, hb⟩
        show c ∈ (setSides _ s.items[i]!.ch []).1 ++ (setSides _ s.items[i]!.ch []).2
        rw [mem_setSides, hch]; exact hc
    · exact ⟨t, List.mem_cons_of_mem _ (List.mem_cons_of_mem _ ht), j, hj, hb⟩
  · simp only [h2, ↓reduceIte]; exact halloc

theorem VSpan.finishTop {v : Nat} (h : VSpan v s) {x : ItemId} (hf : ItemFree s x) {t : TEntry}
    {rest : List TEntry} (hts : s.tstack = t :: rest)
    (hside : getSide t.spans (!s.stackDir[t.topDepth]!) = []) :
    VSpan v (after (finishTstackTop x) s) := by
  show VSpan v ((finishTstackTop x).run s).2
  rw [finishTstackTop_run_eq s x t rest hts]
  set dir := s.stackDir[t.topDepth]! with hdir
  set f : Item → Item := fun it =>
    { it with vs := setSides dir (some s.stackVerts[t.topDepth]!) (some t.vStart),
              ch := getSide t.spans dir } with hf_def
  have hnb : ∀ i, i ≠ x → ¬ Items.Below s.items i x := fun i hi hb =>
    hi (Items.Below.eq_of_no_parent hf.root hb)
  have hch : Items.ch (s.items.modify x f) x = getSide t.spans dir := by
    rw [Items.ch_modify_at x f hf.lt]
  obtain ⟨u, hu, j, hj, hb⟩ := h
  have hjx : j ≠ x := fun hjx => hf.free u hu (hjx ▸ hj)
  have hb' : Items.Below (s.items.modify x f) j (vertItem v) :=
    (Items.Below_modify_of_not_below x f (hnb j hjx)).2 hb
  rw [hts] at hu
  rcases List.mem_cons.1 hu with hu' | hu
  · rw [hu'] at hj
    refine ⟨_, List.mem_cons_self .., x, ?_, .head ?_ hb'⟩
    · show x ∈ (setSides dir [x] []).1 ++ (setSides dir [x] []).2
      rw [mem_setSides]; exact List.mem_singleton.2 rfl
    · show j ∈ Items.ch (s.items.modify x f) x
      rw [hch]; exact (mem_of_getSide_nil dir t.spans hside j).1 hj
  · exact ⟨u, List.mem_cons_of_mem _ hu, j, hj, hb'⟩

/-- A loop preserves a predicate preserved by each taken iteration. -/
theorem loop_pred (P : WalkState → Prop) (cond : WalkM Bool) (body : WalkM Unit) (fuel : Nat)
    (hcond : ∀ s, (cond.run s).2 = s)
    (hbody : ∀ k, (∀ j, j ≤ k → (cond.run (iter body j s)).1 = true) →
      P (iter body k s) → P (after body (iter body k s)))
    (h : P s) : P (after (WalkM.loop fuel cond body) s) := by
  induction fuel generalizing s with
  | zero => exact h
  | succ fuel ih =>
    show P ((WalkM.loop (fuel + 1) cond body).run s).2
    rw [loop_succ_run fuel cond body s (hcond s)]
    by_cases hc : (cond.run s).1 = true
    · simp only [hc, ↓reduceIte]
      have h₁ := hbody 0 (fun j hj => by rw [Nat.le_zero.1 hj]; exact hc) h
      exact ih (s := (body.run s).2) (fun k hk => hbody (k + 1) fun j hj => by
        cases j with
        | zero => exact hc
        | succ j => exact hk j (Nat.le_of_succ_le_succ hj)) h₁
    · simp only [hc]; exact h

theorem mergeLoop_nil (cond : WalkM Bool) (fuel : Nat) (hcond : ∀ s, (cond.run s).2 = s)
    (hs : s.tstack = []) : (after (WalkM.loop fuel cond mergeTstackTops) s).tstack = [] := by
  induction fuel generalizing s with
  | zero => exact hs
  | succ fuel ih =>
    show ((WalkM.loop (fuel + 1) cond mergeTstackTops).run s).2.tstack = []
    rw [loop_succ_run fuel cond mergeTstackTops s (hcond s)]
    by_cases hc : (cond.run s).1 = true
    · simp only [hc, ↓reduceIte]
      exact ih (by show (after mergeTstackTops s).tstack = []; rw [after_mergeTstackTops, hs]; rfl)
    · simp only [hc]; exact hs

/-- A merge loop ending with a nonempty stack keeps `VSpan`. -/
theorem VSpan.mergeLoop {v : Nat} (h : VSpan v s) (cond : WalkM Bool) (fuel : Nat)
    (hcond : ∀ s, (cond.run s).2 = s)
    (hne : (after (WalkM.loop fuel cond mergeTstackTops) s).tstack ≠ []) :
    VSpan v (after (WalkM.loop fuel cond mergeTstackTops) s) := by
  induction fuel generalizing s with
  | zero => exact h
  | succ fuel ih =>
    show VSpan v ((WalkM.loop (fuel + 1) cond mergeTstackTops).run s).2
    have hl := loop_succ_run fuel cond mergeTstackTops s (hcond s)
    rw [hl]
    by_cases hc : (cond.run s).1 = true
    · simp only [hc, ↓reduceIte]
      have hne₁ : ((WalkM.loop (fuel + 1) cond mergeTstackTops).run s).2.tstack ≠ [] := hne
      rw [hl] at hne₁; simp only [hc, ↓reduceIte] at hne₁
      have hne' : (after mergeTstackTops s).tstack ≠ [] := fun h0 =>
        hne₁ (mergeLoop_nil cond fuel hcond h0)
      exact ih (h.mergeTop' hne') hne₁
    · simp only [hc]; exact h

theorem VSpan.closeTwo {v : Nat} (h : VSpan v s) {x : ItemId} (hf : ItemFree s x)
    (hok : CloseTwoOk D s) :
    VSpan v ((finishTstackTop x).run (mergeTstackTops.run s).2).2 := by
  have hne := hok.finish.nonempty
  have h₁ : VSpan v (after mergeTstackTops s) := h.mergeTop' hne
  obtain ⟨t, rest, hts⟩ := List.exists_cons_of_ne_nil hne
  have hside := hok.finish.side
  rw [curE, hts, List.head!_cons] at hside
  exact h₁.finishTop hf.merge hts hside

theorem VSpan.loop1Body {v : Nat} (h : VSpan v s) (hi : s.Inv' D) (hs : Shape s) (hv : v < s.g.nv)
    {d : Nat} {dir : Bool} (hok : Loop1BodyOk D d dir s) (htwo : 2 ≤ s.tstack.length) :
    VSpan v (after (Spqr.loop1Body d dir) s) := by
  have h₁ : VSpan v (l1S₁ d dir s) := by
    show VSpan v ((loop1Type d dir).run s).2
    rw [loop1Type_run]
    split
    · obtain ⟨a, b, rest, hts⟩ := two_entries_of_le htwo
      exact (h.frame' (s' := { s with stackDir := s.stackDir.set! (nxtE s).topDepth dir }) rfl rfl).mergeTop
        (a := a) (b := b) (rest := rest) hts
    · split <;> exact h
  have st₁ : Step D v s (l1S₁ d dir s) := Step.loop1Type hi hs hok.mergeS
  have hv₁ : v < (l1S₁ d dir s).g.nv := by rw [st₁.g]; exact hv
  have hty := loop1Type_result d dir s
  have hty' : l1Ty d dir s ≠ .V := fun h => hty (by show l1Ty d dir s ∈ _; rw [h]; simp)
  obtain ⟨a, b, rest, hts₁⟩ := two_entries_of_le hok.unwrap.two
  have h₂ : VSpan v (l1S₂ d dir s) := h₁.maybeUnwrap st₁.shape hv₁ hty' hok.unwrap hts₁
  have r := maybeUnwrapNxt_spec (v := v) st₁.inv st₁.shape hty hok.unwrap
  exact h₂.closeTwo r.free hok.close

/-- `Step.loop`, at each iterate. -/
theorem Step.iter' {v : Nat} (cond : WalkM Bool) (body : WalkM Unit) (Ok : WalkState → Prop)
    (hbody : ∀ s, v < s.g.nv → s.Inv' D → Shape s → (cond.run s).1 = true → Ok s →
      Step D v s (body.run s).2)
    (hv : v < s.g.nv) (hi : s.Inv' D) (hs : Shape s)
    (hok : ∀ k, (∀ j, j ≤ k → (cond.run (iter body j s)).1 = true) → Ok (iter body k s)) :
    ∀ k, (∀ j, j < k → (cond.run (iter body j s)).1 = true) → Step D v s (iter body k s)
  | 0, _ => Step.refl hi hs
  | k + 1, hk => by
    have st := Step.iter' cond body Ok hbody hv hi hs hok k fun j hj => hk j (by omega)
    rw [iter_succ']
    exact st.trans (hbody _ (by rw [st.g]; exact hv) st.inv st.shape
      (hk k (Nat.lt_succ_self k)) (hok k fun j hj => hk j (by omega)))

theorem VSpan.closeEars {v : Nat} (h : VSpan v s) (hi : s.Inv' D) (hs : Shape s) (hv : v < s.g.nv)
    {nxtV d e : Nat} {dir : Bool} (hok : CloseEarsOk D nxtV d e dir s) :
    VSpan v (after (closeEars nxtV d e dir) s) := by
  have st₁ : Step D v s (ceS₁ nxtV d e s) := Step.pushEdge hi hs nxtV d e hok.e_lt hok.q hok.ends hok.d_le
  have hv₁ : v < (ceS₁ nxtV d e s).g.nv := by rw [st₁.g]; exact hv
  show VSpan v (after (WalkM.loop _ (loop1Cond d) (Spqr.loop1Body d dir)) (ceS₁ nxtV d e s))
  refine loop_pred (VSpan v) _ _ _ (fun _ => rfl) (fun k hk h => ?_) (h.pushEdge nxtV d e)
  have st : Step D v _ (iter (Spqr.loop1Body d dir) k (ceS₁ nxtV d e s)) :=
    Step.iter' (loop1Cond d) (Spqr.loop1Body d dir) (Loop1BodyOk D d dir)
      (fun _ hv hi hs _ hok => Step.loop1Body hi hs hv hok) hv₁ st₁.inv st₁.shape hok.body k
      (fun j hj => hk j (Nat.le_of_lt hj))
  have htwo : 2 ≤ (iter (Spqr.loop1Body d dir) k (ceS₁ nxtV d e s)).tstack.length := by
    have := hk k (Nat.le_refl k)
    rw [run_loop1Cond] at this
    simp only [Bool.and_eq_true, decide_eq_true_eq] at this
    exact this.1
  exact h.loop1Body st.inv st.shape (by rw [st.g]; exact hv₁) (hok.body k hk) htwo

theorem VSpan.mergeLate {v : Nat} (h : VSpan v s) (d : Nat)
    (hne : (after (Spqr.mergeLate d) s).tstack ≠ []) : VSpan v (after (Spqr.mergeLate d) s) := by
  have hne' : ((Spqr.mergeLate d).run s).2.tstack ≠ [] := hne
  show VSpan v ((Spqr.mergeLate d).run s).2
  rw [mergeLate_run] at hne' ⊢
  by_cases hc : (curE s).firstIdx > s.firstOccurrence[d]!
  · simp only [hc, ↓reduceIte] at hne' ⊢
    exact h.mergeLoop _ _ (fun _ => rfl) hne'
  · simp only [hc, ↓reduceIte]; exact h

theorem VSpan.finishP {v : Nat} (h : VSpan v s) (hi : s.Inv' D) (hs : Shape s) (hv : v < s.g.nv)
    {curV lowval : Nat} {isType1 : Bool} (hok : FinishPOk D curV lowval isType1 s) :
    VSpan v (after (Spqr.finishP curV lowval isType1) s) := by
  show VSpan v ((Spqr.finishP curV lowval isType1).run s).2
  simp only [Spqr.finishP, WalkM.run_bind, run_condP]
  by_cases hc : result (condP curV lowval isType1) s = true
  · have hc' : (isType1 && decide (s.tstack.length ≥ 2) && (s.tstack.tail.head!.vStart == curV) &&
        (s.tstack.tail.head!.topDepth == lowval)) = true := hc
    simp only [hc', ↓reduceIte, WalkM.run_bind]
    obtain ⟨hu, hcl⟩ := hok.ok hc
    obtain ⟨a, b, rest, hts⟩ := two_entries_of_le hu.two
    have r := maybeUnwrapNxt_spec (v := v) hi hs (by decide) hu
    exact (h.maybeUnwrap hs hv (by decide) hu hts).closeTwo r.free hcl
  · have hc' : (isType1 && decide (s.tstack.length ≥ 2) && (s.tstack.tail.head!.vStart == curV) &&
        (s.tstack.tail.head!.topDepth == lowval)) = false := Bool.eq_false_iff.2 hc
    simp only [hc', Bool.false_eq_true, ↓reduceIte]
    exact h

theorem finishP_ne_nil {curV lowval : Nat} {isType1 : Bool} (hok : FinishPOk D curV lowval isType1 s)
    (hne : s.tstack ≠ []) : (after (Spqr.finishP curV lowval isType1) s).tstack ≠ [] := by
  show ((Spqr.finishP curV lowval isType1).run s).2.tstack ≠ []
  simp only [Spqr.finishP, WalkM.run_bind, run_condP]
  by_cases hc : result (condP curV lowval isType1) s = true
  · have hc' : (isType1 && decide (s.tstack.length ≥ 2) && (s.tstack.tail.head!.vStart == curV) &&
        (s.tstack.tail.head!.topDepth == lowval)) = true := hc
    simp only [hc', ↓reduceIte, WalkM.run_bind]
    obtain ⟨hu, hcl⟩ := hok.ok hc
    obtain ⟨t, rest, hts⟩ := List.exists_cons_of_ne_nil hcl.finish.nonempty
    show ((finishTstackTop _).run (after mergeTstackTops (after (maybeUnwrapNxt .P) s))).2.tstack ≠ []
    rw [finishTstackTop_run_eq _ _ t rest hts]
    simp
  · have hc' : (isType1 && decide (s.tstack.length ≥ 2) && (s.tstack.tail.head!.vStart == curV) &&
        (s.tstack.tail.head!.topDepth == lowval)) = false := Bool.eq_false_iff.2 hc
    simp only [hc', Bool.false_eq_true, ↓reduceIte]
    exact hne

/-- `Spqr.finishRest` ends with `vertItem curV` spanned whenever it reports `hasVert`. -/
theorem VSpan.finishRest (hi : s.Inv' D) (hs : Shape s) {curV d lowval : Nat} (hv : curV < s.g.nv)
    {isType1 hasVert isSingle : Bool}
    (hok : FinishRestOk D curV d lowval isType1 hasVert isSingle s)
    (h : hasVert = true → VSpan curV s) (hne : hasVert = false → s.tstack ≠ []) :
    wp (Spqr.finishRest curV d lowval isType1 hasVert isSingle)
      (fun hv' s' => hv' = true → VSpan curV s') s := by
  have st₁ : Step D curV s (after (Spqr.finishP curV lowval isType1) s) := Step.finishP hi hs hv hok.p
  show ((finishTail curV d hasVert isSingle).run (after (Spqr.finishP curV lowval isType1) s)).1 = true →
    VSpan curV ((finishTail curV d hasVert isSingle).run (after (Spqr.finishP curV lowval isType1) s)).2
  cases hasVert
  · intro _
    have hne₁ := finishP_ne_nil hok.p (hne rfl)
    have hp := VSpan.pushVert curV d (after (Spqr.finishP curV lowval isType1) s)
    cases isSingle
    · show VSpan curV (after mergeTstackTops (after (pushVertTstack curV d) _))
      obtain ⟨t, rest, hts⟩ := List.exists_cons_of_ne_nil hne₁
      have hts' : (after (pushVertTstack curV d) (after (Spqr.finishP curV lowval isType1) s)).tstack =
          ⟨curV, d, (after (Spqr.finishP curV lowval isType1) s).nxtEdgeIdx,
            setSides (after (Spqr.finishP curV lowval isType1) s).stackDir[d]! [vertItem curV] []⟩ ::
            t :: rest := by
        show (_ :: (after (Spqr.finishP curV lowval isType1) s).tstack) = _; rw [hts]
      exact hp.mergeTop hts'
    · exact hp
  · intro _
    exact (h rfl).finishP hi hs hv hok.p

theorem VSpan.closeVert' {v : Nat} (h : VSpan v s) (hi : s.Inv' D) (hs : Shape s) (hv : v < s.g.nv)
    {curV : Nat} {dir isType1 : Bool} {orig : Nat} {isSingle : Bool}
    (hok : CloseVertOk D curV dir isType1 orig isSingle s) :
    VSpan v (after (closeVert' curV dir isType1 orig isSingle) s) := by
  have st₁ : Step D v s (cvS₁ isType1 orig isSingle s) := Step.vertPre hi hs hv hok.loop3
  have hv₁ : v < (cvS₁ isType1 orig isSingle s).g.nv := by rw [st₁.g]; exact hv
  obtain ⟨st₂, hfree⟩ := Step.vertUnwrap (v := v) st₁.inv st₁.shape (isType1 := isType1)
    (isSingle := cvB₁ isType1 orig isSingle s) (fun h => by subst h; exact hok.unwrap rfl)
  have st₂ : Step D v (cvS₁ isType1 orig isSingle s) (cvS₂ isType1 orig isSingle s) := st₂
  have st₃ : Step D v _ (cvS₃ isType1 orig isSingle s) := Step.mergeTop st₂.inv st₂.shape hok.merge₁
  have st₄ : Step D v _ (cvS₄ isType1 orig isSingle s) := Step.mergeTop st₃.inv st₃.shape hok.merge₂
  have hne₄ : (cvS₄ isType1 orig isSingle s).tstack ≠ [] := hok.retarget.nonempty
  have hne₃ : (cvS₃ isType1 orig isSingle s).tstack ≠ [] := by
    intro h0; apply hne₄
    show (after mergeTstackTops _).tstack = []
    rw [after_mergeTstackTops, h0]; rfl
  have hne₂ : (cvS₂ isType1 orig isSingle s).tstack ≠ [] := by
    intro h0; apply hne₃
    show (after mergeTstackTops _).tstack = []
    rw [after_mergeTstackTops, h0]; rfl
  show VSpan v (after (vertFinish (result (vertUnwrap isType1 (cvB₁ isType1 orig isSingle s))
    (cvS₁ isType1 orig isSingle s)) (cvB₁ isType1 orig isSingle s)) (cvS₅ curV dir isType1 orig isSingle s))
  cases isType1
  · have hS₂ : cvS₂ false orig isSingle s = cvS₁ false orig isSingle s := rfl
    have hS₁ : cvS₁ false orig isSingle s =
        after (WalkM.loop s.tstack.length (loop3Cond orig) mergeTstackTops) s := rfl
    have h₁ : VSpan v (cvS₁ false orig isSingle s) := by
      rw [hS₁]
      exact h.mergeLoop _ _ (fun _ => rfl) (by rw [← hS₁, ← hS₂]; exact hne₂)
    have h₂ : VSpan v (cvS₂ false orig isSingle s) := by rw [hS₂]; exact h₁
    exact (((h₂.mergeTop' hne₃).mergeTop' hne₄).retarget curV dir)
  · have hS₂ : cvS₂ true orig isSingle s =
        after (maybeUnwrapNxt (if isSingle then .S else .R)) (cvS₁ true orig isSingle s) := by
      simp only [cvS₂, after, WalkState.vertUnwrap, ↓reduceIte, WalkM.map_run]; rfl
    have hu : UnwrapOk (if isSingle then .S else .R) (cvS₁ true orig isSingle s) := hok.unwrap rfl
    obtain ⟨a, b, rest, hts₁⟩ := two_entries_of_le hu.two
    have h₂ : VSpan v (cvS₂ true orig isSingle s) := by
      rw [hS₂]
      exact (show VSpan v (cvS₁ true orig isSingle s) from h).maybeUnwrap st₁.shape hv₁
        (by cases isSingle <;> decide) hu hts₁
    have h₅ := ((h₂.mergeTop' hne₃).mergeTop' hne₄).retarget curV dir
    have hf : ItemFree (cvS₂ true orig isSingle s) ((maybeUnwrapNxt (if isSingle then .S else .R)).run
        (cvS₁ true orig isSingle s)).1 := hfree _ rfl
    have hf₅ := ((hf.merge).merge).retarget curV dir
    have hfin := hok.finish rfl
    obtain ⟨t, rest, hts₅⟩ := List.exists_cons_of_ne_nil hfin.nonempty
    have hside := hfin.side
    rw [curE, hts₅, List.head!_cons] at hside
    exact h₅.finishTop hf₅ hts₅ hside

/-- Once the vertex entry of `curV` is on the stack (`hasVert`), `finishEdge` keeps every edge under
`vertItem curV` on the stack (`EarFinish.vert` at the pre-state: `vertItem curV` is spanned by a `base`
entry; the merges and closes keep it spanned or put it under the pushed node). (Checker: `own_vcover`
at `post`/`end`.) -/
theorem finishEdge_vertCover {curV d : Nat} {o : DfsOut} {origTstack : Nat} {hasVert : Bool}
    (h : OwnSite σ n D curV d o origTstack hasVert s) :
    wp (finishEdge curV d o origTstack hasVert) (fun hv' s' => hv' = true → VertCover curV s') s := by
  obtain ⟨sub, base, hlen, hE⟩ := h.book.ear
  have hi := h.rgs.1.inv
  have hs := h.rgs.2.1
  have hv := h.book.v_lt
  by_cases hge : d ≤ o.cls.lowval d
  · have hnv := hE.bd_noVert hge
    subst hnv
    have hge' : o.cls.lowval d ≥ d := hge
    rw [finishEdge_eq]
    simp only [finishEdge', wp_bind, wp_get, wp_stackDir, hge', ↓reduceIte]
    unfold finishBoundary
    simp only [wp_bind, wp_modifyItem, wp_modify, wp_ite, wp_allocItem, wp_makeVs, wp_popTstack, wp_pure]
    simp
  have hlow : o.cls.lowval d < d := Nat.lt_of_not_le hge
  obtain ⟨lv, kind, ho, hl⟩ := ret_of_lowval_lt hlow
  have hok := finishOk_of_guards ho hl h.guards hE hlen hi hs h.hD hv h.book.e_lt h.book.q
    (h.book.ends lv kind ho) h.book.vert
  have hlv : o.cls.lowval d = lv := by rw [ho]; rfl
  have hge' : ¬ (lv ≥ d) := by omega
  have hj : edgeItem s.g o.e < 1 + s.g.nv + s.g.ne := by
    show 1 + s.g.nv + o.e < _; have := hok.e_lt; omega
  have st₀ : Step D curV s (feS₀ d o s) := Step.modifyVs hi hs (edgeItem s.g o.e) _ hj
  have hv₀ : curV < (feS₀ d o s).g.nv := by rw [st₀.g]; exact hv
  have hpre : hasVert = true → VSpan curV (feS₀ d o s) := fun hhv => by
    obtain ⟨t, ht, -, hmem⟩ := hE.vert hhv
    exact VSpan.modifyVs ⟨t, by rw [hE.tstack]; exact List.mem_append_right _ ht, vertItem curV, hmem,
      Relation.ReflTransGen.refl⟩ _ _
  show ((finishEdge curV d o origTstack hasVert).run s).1 = true →
    VertCover curV ((finishEdge curV d o origTstack hasVert).run s).2
  rw [finishEdge_eq]
  simp only [finishEdge', finishTree, finishBack, hlv, hge', ↓reduceIte, WalkM.run_bind, WalkM.get_run,
    run_stackDir, run_makeVs, run_modifyItem]
  by_cases ht : o.cls.isTree = true
  · simp only [ht, ↓reduceIte, closeVert_eq, WalkM.run_bind]
    have st₁ : Step D curV _ (feS₁ d o s) := Step.closeEars st₀.inv st₀.shape hv₀ (hok.ears ht)
    have hv₁ : curV < (feS₁ d o s).g.nv := by rw [st₁.g]; exact hv₀
    have st₂ : Step D curV _ (feS₂ d o s) := Step.mergeLate st₁.inv st₁.shape hv₁ (hok.late ht)
    have hv₂ : curV < (feS₂ d o s).g.nv := by rw [st₂.g]; exact hv₁
    have hne₂ : (feS₂ d o s).tstack ≠ [] := by
      obtain ⟨c, mid, py, vy, hts, -⟩ := hE.loops ht hlow
      rw [hts]; simp
    have h₂ : hasVert = true → VSpan curV (feS₂ d o s) := fun hhv =>
      ((hpre hhv).closeEars st₀.inv st₀.shape hv₀ (hok.ears ht)).mergeLate d hne₂
    cases hasVert
    · simp only [Bool.false_eq_true, ↓reduceIte]
      have := VSpan.finishRest st₂.inv st₂.shape hv₂ (hok.rest_tree ht rfl) h₂ (fun _ => hne₂)
      exact fun hv' => (this hv').vertCover
    · simp only [↓reduceIte, WalkM.run_bind]
      have hcv := hok.vert ht rfl
      have h₃ : VSpan curV (feS₃ curV d o origTstack s) :=
        (h₂ rfl).closeVert' st₂.inv st₂.shape hv₂ hcv
      have st₃ : Step D curV _ (feS₃ curV d o origTstack s) := Step.closeVert' st₂.inv st₂.shape hv₂ hcv
      have hv₃ : curV < (feS₃ curV d o origTstack s).g.nv := by rw [st₃.g]; exact hv₂
      have := VSpan.finishRest st₃.inv st₃.shape hv₃ (hok.rest_vert ht rfl) (fun _ => h₃)
        (fun h => Bool.noConfusion h)
      exact fun hv' => (this hv').vertCover
  · have ht' : o.cls.isTree = false := Bool.eq_false_iff.2 ht
    simp only [ht', Bool.false_eq_true, ↓reduceIte, WalkM.run_bind]
    have hq : Items.ch (feS₀ d o s).items (edgeItem (feS₀ d o s).g o.e) = [] := by
      show Items.ch (s.items.modify _ _) (edgeItem s.g o.e) = []
      rw [Items.ch_modify_ch_eq (edgeItem s.g o.e)
        (fun it => { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) }) (fun _ => rfl)]
      exact hok.q ht'
    have st₁ : Step D curV _ (after (pushEdgeTstack curV lv o.e) (feS₀ d o s)) :=
      Step.pushEdge st₀.inv st₀.shape curV lv o.e hok.e_lt hq (hok.ends ht') (hok.lv_le ht')
    have st₂ : Step D curV _ (feBack curV lv d o s) := Step.frame st₁.inv st₁.shape rfl rfl rfl rfl
    have hv₂ : curV < (feBack curV lv d o s).g.nv := by rw [st₂.g, st₁.g]; exact hv₀
    have hB : hasVert = true → VSpan curV (feBack curV lv d o s) := fun hhv =>
      ((hpre hhv).pushEdge curV lv o.e).frame' rfl rfl
    have hneB : (feBack curV lv d o s).tstack ≠ [] := by
      show (_ :: _) ≠ []; simp
    have := VSpan.finishRest st₂.inv st₂.shape hv₂ (hok.rest_back ht') hB (fun _ => hneB)
    exact fun hv' => (this hv').vertCover

/-! ### `finishEdge_ownedD`: the processed prefix grows by `o.e` -/

/-- Frame of a returning edge relative to its pre-state `s₀`: graph and path fixed, the vertex and
edge items (`< 1 + nv + ne`) keep their children and subtrees, and every entry bottoms at `curV` or
at the bottom of an old entry. -/
structure OFrame (s₀ : WalkState) (curV : Nat) (s : WalkState) : Prop where
  g : s.g = s₀.g
  sv : s.stackVerts = s₀.stackVerts
  ch : ∀ p, p < 1 + s₀.g.nv + s₀.g.ne → Items.ch s.items p = Items.ch s₀.items p
  below : ∀ a, a < 1 + s₀.g.nv + s₀.g.ne → ∀ i, Items.Below s.items a i ↔ Items.Below s₀.items a i
  vstart : ∀ t ∈ s.tstack, t.vStart = curV ∨ ∃ u ∈ s₀.tstack, u.vStart = t.vStart

namespace OFrame
variable {s₀ : WalkState} {curV : Nat}

theorem refl : OFrame s₀ curV s₀ :=
  ⟨rfl, rfl, fun _ _ => rfl, fun _ _ _ => Iff.rfl, fun t ht => Or.inr ⟨t, ht, rfl⟩⟩

theorem frame {s' : WalkState} (h : OFrame s₀ curV s) (hg : s'.g = s.g) (hsv : s'.stackVerts = s.stackVerts)
    (hch : ∀ p, p < 1 + s₀.g.nv + s₀.g.ne → Items.ch s'.items p = Items.ch s.items p)
    (hbl : ∀ a, a < 1 + s₀.g.nv + s₀.g.ne → ∀ i, Items.Below s'.items a i ↔ Items.Below s.items a i)
    (hvs : ∀ t ∈ s'.tstack, t.vStart = curV ∨ ∃ u ∈ s.tstack, u.vStart = t.vStart) :
    OFrame s₀ curV s' :=
  ⟨hg.trans h.g, hsv.trans h.sv, fun p hp => (hch p hp).trans (h.ch p hp),
   fun a ha i => (hbl a ha i).trans (h.below a ha i),
   fun t ht => by
    rcases hvs t ht with hc | ⟨u, hu, hut⟩
    · exact Or.inl hc
    · rcases h.vstart u hu with hc | ⟨u', hu', hu'u⟩
      · exact Or.inl (hut ▸ hc)
      · exact Or.inr ⟨u', hu', hu'u.trans hut⟩⟩

theorem items_eq {s' : WalkState} (h : OFrame s₀ curV s) (hg : s'.g = s.g) (hsv : s'.stackVerts = s.stackVerts)
    (hit : s'.items = s.items)
    (hvs : ∀ t ∈ s'.tstack, t.vStart = curV ∨ ∃ u ∈ s.tstack, u.vStart = t.vStart) : OFrame s₀ curV s' :=
  h.frame hg hsv (fun _ _ => by rw [hit]) (fun _ _ _ => by rw [hit]) hvs

theorem same {s' : WalkState} (h : OFrame s₀ curV s) (hg : s'.g = s.g) (hsv : s'.stackVerts = s.stackVerts)
    (hit : s'.items = s.items) (hts : s'.tstack = s.tstack) : OFrame s₀ curV s' :=
  h.items_eq hg hsv hit (fun t ht => Or.inr ⟨t, hts ▸ ht, rfl⟩)

theorem mergeTop (h : OFrame s₀ curV s) {a b : TEntry} {rest : List TEntry}
    (hts : s.tstack = a :: b :: rest) : OFrame s₀ curV (after mergeTstackTops s) := by
  rw [after_mergeTstackTops_eq hts]
  refine h.items_eq rfl rfl rfl fun t ht => ?_
  rcases List.mem_cons.1 ht with rfl | ht
  · exact Or.inr ⟨b, by rw [hts]; simp, rfl⟩
  · exact Or.inr ⟨t, by rw [hts]; exact List.mem_cons_of_mem _ (List.mem_cons_of_mem _ ht), rfl⟩

theorem mergeTop' (h : OFrame s₀ curV s) (hne : (after mergeTstackTops s).tstack ≠ []) :
    OFrame s₀ curV (after mergeTstackTops s) := by
  rw [after_mergeTstackTops] at hne
  obtain ⟨a, b, rest, hts⟩ := two_of_mergeTop_ne_nil hne
  exact h.mergeTop hts

theorem mergeLoop (h : OFrame s₀ curV s) (cond : WalkM Bool) (fuel : Nat)
    (hcond : ∀ s, (cond.run s).2 = s)
    (hne : (after (WalkM.loop fuel cond mergeTstackTops) s).tstack ≠ []) :
    OFrame s₀ curV (after (WalkM.loop fuel cond mergeTstackTops) s) := by
  induction fuel generalizing s with
  | zero => exact h
  | succ fuel ih =>
    show OFrame s₀ curV ((WalkM.loop (fuel + 1) cond mergeTstackTops).run s).2
    have hl := loop_succ_run fuel cond mergeTstackTops s (hcond s)
    rw [hl]
    by_cases hc : (cond.run s).1 = true
    · simp only [hc, ↓reduceIte]
      have hne₁ : ((WalkM.loop (fuel + 1) cond mergeTstackTops).run s).2.tstack ≠ [] := hne
      rw [hl] at hne₁; simp only [hc, ↓reduceIte] at hne₁
      have hne' : (after mergeTstackTops s).tstack ≠ [] := fun h0 =>
        hne₁ (mergeLoop_nil cond fuel hcond h0)
      exact ih (h.mergeTop' hne') hne₁
    · simp only [hc]; exact h

theorem push (h : OFrame s₀ curV s) (d : Nat) (i : ItemId) :
    OFrame s₀ curV (after (pushTstack curV d i) s) :=
  h.items_eq rfl rfl rfl fun t ht => by
    rcases List.mem_cons.1 ht with rfl | ht
    · exact Or.inl rfl
    · exact Or.inr ⟨t, ht, rfl⟩

theorem pushVert (h : OFrame s₀ curV s) (d : Nat) : OFrame s₀ curV (after (pushVertTstack curV d) s) :=
  h.push d _

theorem pushEdge (h : OFrame s₀ curV s) (d e : Nat) : OFrame s₀ curV (after (pushEdgeTstack curV d e) s) :=
  h.push d _

theorem alloc (h : OFrame s₀ curV s) (ty : NodeType) : OFrame s₀ curV (after (allocItem ty) s) := by
  show OFrame s₀ curV ((allocItem ty).run s).2
  rw [run_allocItem]
  exact h.frame rfl rfl (fun p _ => Items.ch_push_nil _ rfl p) (fun _ _ _ => Items.Below_push_nil _ rfl)
    (fun t ht => Or.inr ⟨t, ht, rfl⟩)

theorem modifyVs (h : OFrame s₀ curV s) (j : ItemId) (vs : Option Nat × Option Nat) :
    OFrame s₀ curV (after (modifyItem j fun it => { it with vs := vs }) s) := by
  show OFrame s₀ curV ((modifyItem j fun it => { it with vs := vs }).run s).2
  rw [run_modifyItem]
  exact h.frame rfl rfl
    (fun p _ => Items.ch_modify_ch_eq (items := s.items) j (fun it => { it with vs := vs }) (fun _ => rfl) p)
    (fun _ _ _ => Items.Below_modify_ch_eq j (fun it => { it with vs := vs }) (fun _ => rfl))
    (fun t ht => Or.inr ⟨t, ht, rfl⟩)

theorem retarget (h : OFrame s₀ curV s) (dir : Bool) {t : TEntry} {rest : List TEntry}
    (hts : s.tstack = t :: rest) : OFrame s₀ curV (after (WalkState.retarget curV dir) s) := by
  show OFrame s₀ curV ((WalkState.retarget curV dir).run s).2
  rw [retarget_run_eq curV dir s t rest hts]
  refine h.items_eq rfl rfl rfl fun u hu => ?_
  rcases List.mem_cons.1 hu with rfl | hu
  · exact Or.inl rfl
  · exact Or.inr ⟨u, by rw [hts]; exact List.mem_cons_of_mem _ hu, rfl⟩

theorem maybeUnwrap (h : OFrame s₀ curV s) (ty : NodeType) {a b : TEntry} {rest : List TEntry}
    (hts : s.tstack = a :: b :: rest) : OFrame s₀ curV (after (maybeUnwrapNxt ty) s) := by
  show OFrame s₀ curV ((maybeUnwrapNxt ty).run s).2
  rw [maybeUnwrapNxt_run_eq ty s a b rest hts _ rfl _ rfl]
  have halloc : OFrame s₀ curV ((allocItem ty).run s).2 := h.alloc ty
  by_cases h1 : ty = .R ∨ s.ternarize = true
  · simp only [h1, ↓reduceIte]; exact halloc
  simp only [h1, ↓reduceIte]
  by_cases h2 : s.items[(getSide b.spans s.stackDir[b.topDepth]!).head!]!.type = ty
  · simp only [h2, ↓reduceIte]
    refine h.items_eq rfl rfl rfl fun u hu => ?_
    rw [hts]
    simp only [List.mem_cons] at hu
    rcases hu with hu | hu | hu
    · exact Or.inr ⟨a, by simp, hu ▸ rfl⟩
    · exact Or.inr ⟨b, by simp, hu ▸ rfl⟩
    · exact Or.inr ⟨u, by simp [hu], rfl⟩
  · simp only [h2, ↓reduceIte]; exact halloc

theorem finishTop (h : OFrame s₀ curV s) {x : ItemId} (hf : ItemFree s x) {t : TEntry}
    {rest : List TEntry} (hts : s.tstack = t :: rest) : OFrame s₀ curV (after (finishTstackTop x) s) := by
  show OFrame s₀ curV ((finishTstackTop x).run s).2
  rw [finishTstackTop_run_eq s x t rest hts]
  have hx : 1 + s₀.g.nv + s₀.g.ne ≤ x := by have := hf.node; rw [h.g] at this; exact this
  refine h.frame rfl rfl (fun p hp => Items.ch_modify_of_ne (p := p) x _ (Nat.ne_of_lt (Nat.lt_of_lt_of_le hp hx)))
    (fun a ha i => Items.Below_modify_of_not_below (a := a) x _ fun hb => ?_) fun u hu => ?_
  · have : a = x := Items.Below.eq_of_no_parent hf.root hb; omega
  · rcases List.mem_cons.1 hu with rfl | hu
    · exact Or.inr ⟨t, by rw [hts]; simp, rfl⟩
    · exact Or.inr ⟨u, by rw [hts]; exact List.mem_cons_of_mem _ hu, rfl⟩

theorem closeTwo (h : OFrame s₀ curV s) {x : ItemId} (hf : ItemFree s x) (hok : CloseTwoOk D s) :
    OFrame s₀ curV ((finishTstackTop x).run (mergeTstackTops.run s).2).2 := by
  have hne := hok.finish.nonempty
  have h₁ : OFrame s₀ curV (after mergeTstackTops s) := h.mergeTop' hne
  obtain ⟨t, rest, hts⟩ := List.exists_cons_of_ne_nil hne
  exact h₁.finishTop hf.merge hts

theorem loop1Body (h : OFrame s₀ curV s) (hi : s.Inv' D) (hs : Shape s)
    {d : Nat} {dir : Bool} (hok : Loop1BodyOk D d dir s) (htwo : 2 ≤ s.tstack.length) :
    OFrame s₀ curV (after (Spqr.loop1Body d dir) s) := by
  have h₁ : OFrame s₀ curV (l1S₁ d dir s) := by
    show OFrame s₀ curV ((loop1Type d dir).run s).2
    rw [loop1Type_run]
    split
    · obtain ⟨a, b, rest, hts⟩ := two_entries_of_le htwo
      exact (h.same (s' := { s with stackDir := s.stackDir.set! (nxtE s).topDepth dir }) rfl rfl rfl rfl).mergeTop
        (a := a) (b := b) (rest := rest) hts
    · split <;> exact h
  have st₁ : Step D curV s (l1S₁ d dir s) := Step.loop1Type hi hs hok.mergeS
  have hty := loop1Type_result d dir s
  obtain ⟨a, b, rest, hts₁⟩ := two_entries_of_le hok.unwrap.two
  have h₂ : OFrame s₀ curV (l1S₂ d dir s) := h₁.maybeUnwrap _ hts₁
  have r := maybeUnwrapNxt_spec (v := curV) st₁.inv st₁.shape hty hok.unwrap
  exact h₂.closeTwo r.free hok.close

theorem closeEars (h : OFrame s₀ curV s) (hi : s.Inv' D) (hs : Shape s) (hv : curV < s.g.nv)
    {d e : Nat} {dir : Bool} (hok : CloseEarsOk D curV d e dir s) :
    OFrame s₀ curV (after (Spqr.closeEars curV d e dir) s) := by
  have st₁ : Step D curV s (ceS₁ curV d e s) := Step.pushEdge hi hs curV d e hok.e_lt hok.q hok.ends hok.d_le
  have hv₁ : curV < (ceS₁ curV d e s).g.nv := by rw [st₁.g]; exact hv
  show OFrame s₀ curV (after (WalkM.loop _ (loop1Cond d) (Spqr.loop1Body d dir)) (ceS₁ curV d e s))
  refine loop_pred (OFrame s₀ curV) _ _ _ (fun _ => rfl) (fun k hk h => ?_) (h.pushEdge d e)
  have st : Step D curV _ (iter (Spqr.loop1Body d dir) k (ceS₁ curV d e s)) :=
    Step.iter' (loop1Cond d) (Spqr.loop1Body d dir) (Loop1BodyOk D d dir)
      (fun _ hv hi hs _ hok => Step.loop1Body hi hs hv hok) hv₁ st₁.inv st₁.shape hok.body k
      (fun j hj => hk j (Nat.le_of_lt hj))
  have htwo : 2 ≤ (iter (Spqr.loop1Body d dir) k (ceS₁ curV d e s)).tstack.length := by
    have := hk k (Nat.le_refl k)
    rw [run_loop1Cond] at this
    simp only [Bool.and_eq_true, decide_eq_true_eq] at this
    exact this.1
  exact h.loop1Body st.inv st.shape (hok.body k hk) htwo

theorem mergeLate (h : OFrame s₀ curV s) (d : Nat)
    (hne : (after (Spqr.mergeLate d) s).tstack ≠ []) : OFrame s₀ curV (after (Spqr.mergeLate d) s) := by
  have hne' : ((Spqr.mergeLate d).run s).2.tstack ≠ [] := hne
  show OFrame s₀ curV ((Spqr.mergeLate d).run s).2
  rw [mergeLate_run] at hne' ⊢
  by_cases hc : (curE s).firstIdx > s.firstOccurrence[d]!
  · simp only [hc, ↓reduceIte] at hne' ⊢
    exact h.mergeLoop _ _ (fun _ => rfl) hne'
  · simp only [hc, ↓reduceIte]; exact h

end OFrame

theorem OFrame.closeVert' (h : OFrame s₀ curV s) (hi : s.Inv' D) (hs : Shape s) (hv : curV < s.g.nv)
    {dir isType1 : Bool} {orig : Nat} {isSingle : Bool}
    (hok : CloseVertOk D curV dir isType1 orig isSingle s) :
    OFrame s₀ curV (after (closeVert' curV dir isType1 orig isSingle) s) := by
  have st₁ : Step D curV s (cvS₁ isType1 orig isSingle s) := Step.vertPre hi hs hv hok.loop3
  obtain ⟨st₂, hfree⟩ := Step.vertUnwrap (v := curV) st₁.inv st₁.shape (isType1 := isType1)
    (isSingle := cvB₁ isType1 orig isSingle s) (fun h => by subst h; exact hok.unwrap rfl)
  have st₂ : Step D curV (cvS₁ isType1 orig isSingle s) (cvS₂ isType1 orig isSingle s) := st₂
  have hne₄ : (cvS₄ isType1 orig isSingle s).tstack ≠ [] := hok.retarget.nonempty
  have hne₃ : (cvS₃ isType1 orig isSingle s).tstack ≠ [] := by
    intro h0; apply hne₄
    show (after mergeTstackTops _).tstack = []
    rw [after_mergeTstackTops, h0]; rfl
  have hne₂ : (cvS₂ isType1 orig isSingle s).tstack ≠ [] := by
    intro h0; apply hne₃
    show (after mergeTstackTops _).tstack = []
    rw [after_mergeTstackTops, h0]; rfl
  obtain ⟨t₄, rest₄, hts₄⟩ := List.exists_cons_of_ne_nil hne₄
  show OFrame s₀ curV (after (vertFinish (result (vertUnwrap isType1 (cvB₁ isType1 orig isSingle s))
    (cvS₁ isType1 orig isSingle s)) (cvB₁ isType1 orig isSingle s)) (cvS₅ curV dir isType1 orig isSingle s))
  cases isType1
  · have hS₂ : cvS₂ false orig isSingle s = cvS₁ false orig isSingle s := rfl
    have hS₁ : cvS₁ false orig isSingle s =
        after (WalkM.loop s.tstack.length (loop3Cond orig) mergeTstackTops) s := rfl
    have h₁ : OFrame s₀ curV (cvS₁ false orig isSingle s) := by
      rw [hS₁]
      exact h.mergeLoop _ _ (fun _ => rfl) (by rw [← hS₁, ← hS₂]; exact hne₂)
    have h₂ : OFrame s₀ curV (cvS₂ false orig isSingle s) := by rw [hS₂]; exact h₁
    exact ((h₂.mergeTop' hne₃).mergeTop' hne₄).retarget dir hts₄
  · have hS₂ : cvS₂ true orig isSingle s =
        after (maybeUnwrapNxt (if isSingle then .S else .R)) (cvS₁ true orig isSingle s) := by
      simp only [cvS₂, after, WalkState.vertUnwrap, ↓reduceIte, WalkM.map_run]; rfl
    have hu : UnwrapOk (if isSingle then .S else .R) (cvS₁ true orig isSingle s) := hok.unwrap rfl
    obtain ⟨a, b, rest, hts₁⟩ := two_entries_of_le hu.two
    have h₂ : OFrame s₀ curV (cvS₂ true orig isSingle s) := by
      rw [hS₂]
      exact (show OFrame s₀ curV (cvS₁ true orig isSingle s) from h).maybeUnwrap _ hts₁
    have h₅ := ((h₂.mergeTop' hne₃).mergeTop' hne₄).retarget dir hts₄
    have hf : ItemFree (cvS₂ true orig isSingle s) ((maybeUnwrapNxt (if isSingle then .S else .R)).run
        (cvS₁ true orig isSingle s)).1 := hfree _ rfl
    have hf₅ := ((hf.merge).merge).retarget curV dir
    have hfin := hok.finish rfl
    obtain ⟨t, rest, hts₅⟩ := List.exists_cons_of_ne_nil hfin.nonempty
    exact h₅.finishTop hf₅ hts₅

theorem OFrame.finishP (h : OFrame s₀ curV s) (hi : s.Inv' D) (hs : Shape s)
    {lowval : Nat} {isType1 : Bool} (hok : FinishPOk D curV lowval isType1 s) :
    OFrame s₀ curV (after (Spqr.finishP curV lowval isType1) s) := by
  show OFrame s₀ curV ((Spqr.finishP curV lowval isType1).run s).2
  simp only [Spqr.finishP, WalkM.run_bind, run_condP]
  by_cases hc : result (condP curV lowval isType1) s = true
  · have hc' : (isType1 && decide (s.tstack.length ≥ 2) && (s.tstack.tail.head!.vStart == curV) &&
        (s.tstack.tail.head!.topDepth == lowval)) = true := hc
    simp only [hc', ↓reduceIte, WalkM.run_bind]
    obtain ⟨hu, hcl⟩ := hok.ok hc
    obtain ⟨a, b, rest, hts⟩ := two_entries_of_le hu.two
    have r := maybeUnwrapNxt_spec (v := curV) hi hs (by decide) hu
    exact (h.maybeUnwrap _ hts).closeTwo r.free hcl
  · have hc' : (isType1 && decide (s.tstack.length ≥ 2) && (s.tstack.tail.head!.vStart == curV) &&
        (s.tstack.tail.head!.topDepth == lowval)) = false := Bool.eq_false_iff.2 hc
    simp only [hc', Bool.false_eq_true, ↓reduceIte]
    exact h

theorem OFrame.finishRest (h : OFrame s₀ curV s) (hi : s.Inv' D) (hs : Shape s)
    {d lowval : Nat} {isType1 hasVert isSingle : Bool}
    (hok : FinishRestOk D curV d lowval isType1 hasVert isSingle s) :
    OFrame s₀ curV (after (Spqr.finishRest curV d lowval isType1 hasVert isSingle) s) := by
  have h₁ := h.finishP hi hs hok.p
  show OFrame s₀ curV ((finishTail curV d hasVert isSingle).run (after (Spqr.finishP curV lowval isType1) s)).2
  cases hasVert
  · have hp := h₁.pushVert d
    cases isSingle
    · show OFrame s₀ curV (after mergeTstackTops (after (pushVertTstack curV d) _))
      exact hp.items_eq rfl rfl (by rw [after_mergeTstackTops]) fun t ht => by
        rw [after_mergeTstackTops] at ht
        obtain ⟨a, b, rest, hab⟩ := two_of_mergeTop_ne_nil (List.ne_nil_of_mem ht)
        rw [hab] at ht; simp only [WalkM.mergeTop, List.mem_cons] at ht
        rcases ht with rfl | ht
        · exact Or.inr ⟨b, by rw [hab]; simp, rfl⟩
        · exact Or.inr ⟨t, by rw [hab]; simp [ht], rfl⟩
    · exact hp
  · exact h₁

/-- Sources of an edge held by a new entry: the initial new entries `N₀`, the absorbed old entries
`Bh`, or the subtree of `vertItem curV` (all read in `s₂`). -/
def XEdges (s₂ : WalkState) (N₀ Bh : List TEntry) (curV e : Nat) : Prop :=
  (∃ t ∈ N₀, t.edges s₂.g s₂.items e) ∨ (∃ b ∈ Bh, b.edges s₂.g s₂.items e) ∨
    Items.EdgeBelow s₂.g s₂.items (vertItem curV) e

/-- Stack decomposition after the setup state `s₂` (stack `N₀ ++ Bh ++ Bt`): the stack is `N ++ Bt`
with `Bt` untouched, and the new entries `N` hold every `N₀`/`Bh` edge and only `XEdges`. -/
def Dec (s₂ : WalkState) (N₀ Bh Bt : List TEntry) (curV : Nat) (s : WalkState) : Prop :=
  ∃ N : List TEntry, s.tstack = N ++ Bt ∧ N ≠ [] ∧
    (∀ e, e < s₂.g.ne → ((∃ t ∈ N₀, t.edges s₂.g s₂.items e) ∨ ∃ b ∈ Bh, b.edges s₂.g s₂.items e) →
      ∃ t ∈ N, t.edges s₂.g s.items e) ∧
    (∀ t ∈ N, ∀ e, e < s₂.g.ne → t.edges s₂.g s.items e → XEdges s₂ N₀ Bh curV e) ∧
    (∀ t ∈ Bt, ∀ e, e < s₂.g.ne → (t.edges s₂.g s.items e ↔ t.edges s₂.g s₂.items e))

/-- `Dec` for some split `B = Bh ++ Bt` whose absorbed part bottoms at `curV`. -/
def DecE (s₂ : WalkState) (N₀ B : List TEntry) (curV : Nat) (s : WalkState) : Prop :=
  ∃ Bh Bt, B = Bh ++ Bt ∧ (∀ b ∈ Bh, b.vStart = curV) ∧ Dec s₂ N₀ Bh Bt curV s

namespace Dec
variable {s₂ : WalkState} {N₀ Bh Bt : List TEntry} {curV : Nat}

theorem init {B : List TEntry} (hts : s₂.tstack = N₀ ++ B) (hN : N₀ ≠ []) : Dec s₂ N₀ [] B curV s₂ :=
  ⟨N₀, hts, hN,
   fun e _ h => by
    rcases h with h | ⟨b, hb, _⟩
    · exact h
    · simp at hb,
   fun t ht e _ he => Or.inl ⟨t, ht, he⟩, fun _ _ _ _ => Iff.rfl⟩

theorem decE (h : Dec s₂ N₀ Bh Bt curV s) (hBh : ∀ b ∈ Bh, b.vStart = curV) :
    DecE s₂ N₀ (Bh ++ Bt) curV s :=
  ⟨Bh, Bt, rfl, hBh, h⟩

theorem frame {s' : WalkState} (h : Dec s₂ N₀ Bh Bt curV s) (hit : s'.items = s.items)
    (hts : s'.tstack = s.tstack) : Dec s₂ N₀ Bh Bt curV s' := by
  rw [Dec, hit, hts]; exact h

theorem items_congr {s' : WalkState} (h : Dec s₂ N₀ Bh Bt curV s) (hts : s'.tstack = s.tstack)
    (hbl : ∀ a i, Items.Below s'.items a i ↔ Items.Below s.items a i) : Dec s₂ N₀ Bh Bt curV s' := by
  have hed : ∀ (t : TEntry) e, t.edges s₂.g s'.items e ↔ t.edges s₂.g s.items e := fun t e => by
    simp only [TEntry.edges, Items.EdgeBelow, hbl]
  obtain ⟨N, h1, h3, h5, h6, h7⟩ := h
  refine ⟨N, hts.trans h1, h3, fun e he hy => ?_, fun t ht e he hte => ?_,
    fun t ht e he => (hed t e).trans (h7 t ht e he)⟩
  · obtain ⟨t, ht, hte⟩ := h5 e he hy; exact ⟨t, ht, (hed t e).2 hte⟩
  · exact h6 t ht e he ((hed t e).1 hte)

theorem alloc (h : Dec s₂ N₀ Bh Bt curV s) (ty : NodeType) :
    Dec s₂ N₀ Bh Bt curV (after (allocItem ty) s) := by
  show Dec s₂ N₀ Bh Bt curV ((allocItem ty).run s).2
  rw [run_allocItem]
  exact h.items_congr rfl fun _ _ => Items.Below_push_nil _ rfl

theorem modifyVs (h : Dec s₂ N₀ Bh Bt curV s) (j : ItemId) (vs : Option Nat × Option Nat) :
    Dec s₂ N₀ Bh Bt curV (after (modifyItem j fun it => { it with vs := vs }) s) := by
  show Dec s₂ N₀ Bh Bt curV ((modifyItem j fun it => { it with vs := vs }).run s).2
  rw [run_modifyItem]
  exact h.items_congr rfl fun _ _ => Items.Below_modify_ch_eq j (fun it => { it with vs := vs }) (fun _ => rfl)

theorem two (h : Dec s₂ N₀ Bh Bt curV s) (hlen : Bt.length + 2 ≤ s.tstack.length) :
    ∃ a b N', s.tstack = a :: b :: (N' ++ Bt) ∧ Dec s₂ N₀ Bh Bt curV s := by
  obtain ⟨N, h1, h3, h5, h6, h7⟩ := h
  have hl : s.tstack.length = N.length + Bt.length := by rw [h1, List.length_append]
  rcases N with _ | ⟨a, _ | ⟨b, N'⟩⟩
  · exact absurd rfl h3
  · simp only [List.length_singleton] at hl; omega
  · exact ⟨a, b, N', h1, _, h1, h3, h5, h6, h7⟩

/-- Merging two new entries. -/
theorem mergeTop (h : Dec s₂ N₀ Bh Bt curV s) (hlen : Bt.length + 2 ≤ s.tstack.length) :
    Dec s₂ N₀ Bh Bt curV (after mergeTstackTops s) := by
  obtain ⟨a, b, N', hts, N, h1, h3, h5, h6, h7⟩ := h.two hlen
  have hN : N = a :: b :: N' :=
    List.append_cancel_right (show N ++ Bt = (a :: b :: N') ++ Bt from h1.symm.trans hts)
  subst hN
  rw [after_mergeTstackTops_eq hts]
  refine ⟨TEntry.mergeInto a b :: N', rfl, List.cons_ne_nil _ _, fun e he hy => ?_,
    fun t ht e he hte => ?_, fun t ht e he => h7 t ht e he⟩
  · obtain ⟨t, ht, hte⟩ := h5 e he hy
    simp only [List.mem_cons] at ht
    rcases ht with rfl | rfl | ht
    · exact ⟨_, List.mem_cons_self .., (TEntry.edges_mergeInto _ _ _).2 (Or.inl hte)⟩
    · exact ⟨_, List.mem_cons_self .., (TEntry.edges_mergeInto _ _ _).2 (Or.inr hte)⟩
    · exact ⟨t, List.mem_cons_of_mem _ ht, hte⟩
  · rcases List.mem_cons.1 ht with rfl | ht
    · rcases (TEntry.edges_mergeInto _ _ _).1 hte with hte | hte
      · exact h6 a (by simp) e he hte
      · exact h6 b (by simp) e he hte
    · exact h6 t (by simp [ht]) e he hte

/-- Absorbing the old entry `b` (the top of `Bt`) into the single new entry. -/
theorem mergeTop_absorb {b : TEntry} {Bt' : List TEntry} (h : Dec s₂ N₀ Bh (b :: Bt') curV s)
    (hlen : s.tstack.length ≤ Bt'.length + 2) :
    Dec s₂ N₀ (Bh ++ [b]) Bt' curV (after mergeTstackTops s) := by
  obtain ⟨N, h1, h3, h5, h6, h7⟩ := h
  have hl : s.tstack.length = N.length + (Bt'.length + 1) := by rw [h1, List.length_append]; rfl
  obtain ⟨a, N', rfl⟩ := List.exists_cons_of_ne_nil h3
  have hN' : N' = [] := by
    rw [List.length_cons] at hl
    exact List.eq_nil_of_length_eq_zero (by omega)
  subst hN'
  have hts : s.tstack = a :: b :: Bt' := h1
  have hb' := h7 b (List.mem_cons_self ..)
  rw [after_mergeTstackTops_eq hts]
  refine ⟨[TEntry.mergeInto a b], rfl, List.cons_ne_nil _ _, fun e he hy => ?_,
    fun t ht e he hte => ?_, fun t ht e he => h7 t (List.mem_cons_of_mem _ ht) e he⟩
  · refine ⟨_, List.mem_singleton_self _, (TEntry.edges_mergeInto _ _ _).2 ?_⟩
    rcases hy with hy | ⟨u, hu, hue⟩
    · obtain ⟨t, ht, hte⟩ := h5 e he (Or.inl hy)
      rw [List.mem_singleton] at ht; subst ht; exact Or.inl hte
    · rw [List.mem_append, List.mem_singleton] at hu
      rcases hu with hu | rfl
      · obtain ⟨t, ht, hte⟩ := h5 e he (Or.inr ⟨u, hu, hue⟩)
        rw [List.mem_singleton] at ht; subst ht; exact Or.inl hte
      · exact Or.inr ((hb' e he).2 hue)
  · rw [List.mem_singleton] at ht; subst ht
    rcases (TEntry.edges_mergeInto _ _ _).1 hte with hte | hte
    · rcases h6 a (List.mem_singleton_self _) e he hte with hx | ⟨u, hu, hue⟩ | hx
      · exact Or.inl hx
      · exact Or.inr (Or.inl ⟨u, List.mem_append_left _ hu, hue⟩)
      · exact Or.inr (Or.inr hx)
    · exact Or.inr (Or.inl ⟨b, List.mem_append_right _ (List.mem_singleton_self _), (hb' e he).1 hte⟩)

/-- Unwrapping a new `nxt`. -/
theorem unwrap (h : Dec s₂ N₀ Bh Bt curV s) (hg : s.g = s₂.g) (hs : Shape s) {ty : NodeType}
    (hty : ty ∉ [NodeType.F, .V, .Q]) (hok : UnwrapOk ty s) (hlen : Bt.length + 2 ≤ s.tstack.length) :
    Dec s₂ N₀ Bh Bt curV (after (maybeUnwrapNxt ty) s) := by
  obtain ⟨a, b, N', hts, N, h1, h3, h5, h6, h7⟩ := h.two hlen
  have hN : N = a :: b :: N' :=
    List.append_cancel_right (show N ++ Bt = (a :: b :: N') ++ Bt from h1.symm.trans hts)
  subst hN
  obtain ⟨-, b', hts', -, hb', hall⟩ := unwrap_edges hs hty hok hts
  rw [hg] at hb' hall
  refine ⟨a :: b' :: N', hts', List.cons_ne_nil _ _, fun e he hy => ?_, fun t ht e he hte => ?_,
    fun t ht e he => (hall t e).trans (h7 t ht e he)⟩
  · obtain ⟨t, ht, hte⟩ := h5 e he hy
    simp only [List.mem_cons] at ht
    rcases ht with rfl | rfl | ht
    · exact ⟨t, by simp, (hall t e).2 hte⟩
    · exact ⟨b', by simp, (hb' e (hg ▸ he)).2 hte⟩
    · exact ⟨t, by simp [ht], (hall t e).2 hte⟩
  · simp only [List.mem_cons] at ht
    rcases ht with rfl | rfl | ht
    · exact h6 t (by simp) e he ((hall t e).1 hte)
    · exact h6 b (by simp) e he ((hb' e (hg ▸ he)).1 hte)
    · exact h6 t (by simp [ht]) e he ((hall t e).1 hte)

/-- Unwrapping the old `nxt = b` (absorbed into `Bh`). -/
theorem unwrap_absorb {b : TEntry} {Bt' : List TEntry} (h : Dec s₂ N₀ Bh (b :: Bt') curV s)
    (hg : s.g = s₂.g) (hs : Shape s) {ty : NodeType} (hty : ty ∉ [NodeType.F, .V, .Q])
    (hok : UnwrapOk ty s) (hlen : s.tstack.length ≤ Bt'.length + 2) :
    Dec s₂ N₀ (Bh ++ [b]) Bt' curV (after (maybeUnwrapNxt ty) s) := by
  obtain ⟨N, h1, h3, h5, h6, h7⟩ := h
  have hl : s.tstack.length = N.length + (Bt'.length + 1) := by rw [h1, List.length_append]; rfl
  obtain ⟨a, N', rfl⟩ := List.exists_cons_of_ne_nil h3
  have hN' : N' = [] := by
    rw [List.length_cons] at hl
    exact List.eq_nil_of_length_eq_zero (by omega)
  subst hN'
  have hts : s.tstack = a :: b :: Bt' := h1
  obtain ⟨-, b', hts', -, hb', hall⟩ := unwrap_edges hs hty hok hts
  rw [hg] at hb' hall
  have hbe := h7 b (List.mem_cons_self ..)
  refine ⟨[a, b'], hts', List.cons_ne_nil _ _, fun e he hy => ?_, fun t ht e he hte => ?_,
    fun t ht e he => (hall t e).trans (h7 t (List.mem_cons_of_mem _ ht) e he)⟩
  · rcases hy with hy | ⟨u, hu, hue⟩
    · obtain ⟨t, ht, hte⟩ := h5 e he (Or.inl hy)
      rw [List.mem_singleton] at ht; subst ht; exact ⟨t, by simp, (hall t e).2 hte⟩
    · rw [List.mem_append, List.mem_singleton] at hu
      rcases hu with hu | rfl
      · obtain ⟨t, ht, hte⟩ := h5 e he (Or.inr ⟨u, hu, hue⟩)
        rw [List.mem_singleton] at ht; subst ht; exact ⟨t, by simp, (hall t e).2 hte⟩
      · exact ⟨b', by simp, (hb' e (hg ▸ he)).2 ((hbe e he).2 hue)⟩
  · simp only [List.mem_cons, List.not_mem_nil, or_false] at ht
    rcases ht with rfl | rfl
    · rcases h6 t (List.mem_singleton_self _) e he ((hall t e).1 hte) with hx | ⟨u, hu, hue⟩ | hx
      · exact Or.inl hx
      · exact Or.inr (Or.inl ⟨u, List.mem_append_left _ hu, hue⟩)
      · exact Or.inr (Or.inr hx)
    · exact Or.inr (Or.inl ⟨b, List.mem_append_right _ (List.mem_singleton_self _),
        (hbe e he).1 ((hb' e (hg ▸ he)).1 hte)⟩)

theorem finishTop (h : Dec s₂ N₀ Bh Bt curV s) (hg : s.g = s₂.g) {x : ItemId} (hf : ItemFree s x)
    {t : TEntry} {rest : List TEntry} (hts : s.tstack = t :: rest)
    (hside : getSide t.spans (!s.stackDir[t.topDepth]!) = []) :
    Dec s₂ N₀ Bh Bt curV (after (finishTstackTop x) s) := by
  obtain ⟨-, t', hts', -, ht', hrest⟩ := finishTop_edges hf hts hside
  obtain ⟨N, h1, h3, h5, h6, h7⟩ := h
  rw [hg] at ht' hrest
  obtain ⟨t₀, N', rfl⟩ := List.exists_cons_of_ne_nil h3
  have hab : t₀ = t ∧ N' ++ Bt = rest := by
    have := h1.symm.trans hts; simp only [List.cons_append, List.cons.injEq] at this; exact this
  obtain ⟨rfl, hrest'⟩ := hab
  have hmem : ∀ u ∈ N' ++ Bt, u ∈ rest := fun u hu => hrest' ▸ hu
  refine ⟨t' :: N', by rw [hts', ← hrest']; rfl, List.cons_ne_nil _ _, fun e he hy => ?_,
    fun u hu e he hue => ?_,
    fun u hu e he => (hrest u (hmem u (List.mem_append_right _ hu)) e).trans (h7 u hu e he)⟩
  · obtain ⟨u, hu, hue⟩ := h5 e he hy
    rcases List.mem_cons.1 hu with rfl | hu
    · exact ⟨t', List.mem_cons_self .., (ht' e (hg ▸ he)).2 hue⟩
    · exact ⟨u, List.mem_cons_of_mem _ hu, (hrest u (hmem u (List.mem_append_left _ hu)) e).2 hue⟩
  · rcases List.mem_cons.1 hu with rfl | hu
    · exact h6 t₀ (List.mem_cons_self ..) e he ((ht' e (hg ▸ he)).1 hue)
    · exact h6 u (List.mem_cons_of_mem _ hu) e he ((hrest u (hmem u (List.mem_append_left _ hu)) e).1 hue)

theorem retarget (h : Dec s₂ N₀ Bh Bt curV s) (dir : Bool) {t : TEntry} {rest : List TEntry}
    (hts : s.tstack = t :: rest) : Dec s₂ N₀ Bh Bt curV (after (WalkState.retarget curV dir) s) := by
  show Dec s₂ N₀ Bh Bt curV ((WalkState.retarget curV dir).run s).2
  rw [retarget_run_eq curV dir s t rest hts]
  obtain ⟨N, h1, h3, h5, h6, h7⟩ := h
  obtain ⟨t₀, N', rfl⟩ := List.exists_cons_of_ne_nil h3
  have hab : t₀ = t ∧ N' ++ Bt = rest := by
    have := h1.symm.trans hts; simp only [List.cons_append, List.cons.injEq] at this; exact this
  obtain ⟨rfl, hrest'⟩ := hab
  have hed : ∀ e, TEntry.edges s₂.g s.items
      { t₀ with vStart := curV, spans := setSides (!dir) (t₀.spans.1 ++ t₀.spans.2) [] } e ↔
      t₀.edges s₂.g s.items e := fun e => by
    simp only [TEntry.edges, mem_setSides]
  refine ⟨_ :: N', by rw [← hrest']; rfl, List.cons_ne_nil _ _, fun e he hy => ?_,
    fun u hu e he hue => ?_, fun u hu e he => h7 u hu e he⟩
  · obtain ⟨u, hu, hue⟩ := h5 e he hy
    rcases List.mem_cons.1 hu with rfl | hu
    · exact ⟨_, List.mem_cons_self .., (hed e).2 hue⟩
    · exact ⟨u, List.mem_cons_of_mem _ hu, hue⟩
  · rcases List.mem_cons.1 hu with rfl | hu
    · exact h6 t₀ (List.mem_cons_self ..) e he ((hed e).1 hue)
    · exact h6 u (List.mem_cons_of_mem _ hu) e he hue

/-- The vertex push of `finishTail`. -/
theorem pushVert (h : Dec s₂ N₀ Bh Bt curV s) (d : Nat)
    (hbl : ∀ i, Items.Below s.items (vertItem curV) i ↔ Items.Below s₂.items (vertItem curV) i) :
    Dec s₂ N₀ Bh Bt curV (after (pushVertTstack curV d) s) := by
  obtain ⟨N, h1, h3, h5, h6, h7⟩ := h
  refine ⟨⟨curV, d, s.nxtEdgeIdx, setSides s.stackDir[d]! [vertItem curV] []⟩ :: N,
    by show _ :: s.tstack = _; rw [h1]; rfl, List.cons_ne_nil _ _, fun e he hy => ?_,
    fun u hu e he hue => ?_, fun u hu e he => h7 u hu e he⟩
  · obtain ⟨u, hu, hue⟩ := h5 e he hy
    exact ⟨u, List.mem_cons_of_mem _ hu, hue⟩
  · rcases List.mem_cons.1 hu with rfl | hu
    · obtain ⟨i, hi, hb⟩ := hue
      rw [mem_setSides, List.mem_singleton] at hi; subst hi
      exact Or.inr (Or.inr ((hbl _).1 hb))
    · exact h6 u hu e he hue

theorem pushVert_merge (h : Dec s₂ N₀ Bh Bt curV s) (d : Nat)
    (hbl : ∀ i, Items.Below s.items (vertItem curV) i ↔ Items.Below s₂.items (vertItem curV) i) :
    Dec s₂ N₀ Bh Bt curV (after mergeTstackTops (after (pushVertTstack curV d) s)) := by
  obtain ⟨N, h1, h3, -, -, -⟩ := id h
  refine (h.pushVert d hbl).mergeTop ?_
  show Bt.length + 2 ≤ (_ :: s.tstack).length
  rw [List.length_cons, h1, List.length_append]
  have := List.length_pos_of_ne_nil h3
  omega

/-- The P merge: `nxt` is new, or the old `(curV, lowval)` entry on top of `Bt`. -/
theorem closeTwo (h : Dec s₂ N₀ Bh Bt curV s) (hg : s.g = s₂.g) (hs : Shape s) {ty : NodeType}
    (hty : ty ∉ [NodeType.F, .V, .Q]) (hok : UnwrapOk ty s) {x : ItemId}
    (hf : ItemFree (after (maybeUnwrapNxt ty) s) x) (hcl : CloseTwoOk D (after (maybeUnwrapNxt ty) s))
    (hcase : Bt.length + 2 ≤ s.tstack.length ∨ s.tstack.tail.head!.vStart = curV) :
    ∃ Bh' Bt', Bh ++ Bt = Bh' ++ Bt' ∧ (∀ b ∈ Bh', b.vStart = curV ∨ b ∈ Bh) ∧
      Dec s₂ N₀ Bh' Bt' curV
        ((finishTstackTop x).run (mergeTstackTops.run (after (maybeUnwrapNxt ty) s)).2).2 := by
  obtain ⟨a, b, rest, hts⟩ := two_entries_of_le hok.two
  have hg₁ : (after (maybeUnwrapNxt ty) s).g = s.g := (unwrap_edges hs hty hok hts).1
  obtain ⟨-, b', hts', hbv, -, -⟩ := unwrap_edges hs hty hok hts
  have hne := hcl.finish.nonempty
  obtain ⟨t, rest', hts₂⟩ := List.exists_cons_of_ne_nil hne
  have hside := hcl.finish.side
  rw [curE, hts₂, List.head!_cons] at hside
  have hg₂ : (after mergeTstackTops (after (maybeUnwrapNxt ty) s)).g = s₂.g := by
    rw [after_mergeTstackTops]; exact hg₁.trans hg
  by_cases hlen : Bt.length + 2 ≤ s.tstack.length
  · have h₁ := h.unwrap hg hs hty hok hlen
    have h₂ := h₁.mergeTop (by rw [hts']; rw [hts] at hlen; exact hlen)
    exact ⟨Bh, Bt, rfl, fun b hb => Or.inr hb, h₂.finishTop hg₂ hf.merge hts₂ hside⟩
  · have hnxt : b.vStart = curV := by
      rcases hcase with hcase | hcase
      · exact absurd hcase hlen
      · rw [hts] at hcase; exact hcase
    obtain ⟨N, h1, h3, -, -, -⟩ := id h
    have hl : s.tstack.length = N.length + Bt.length := by rw [h1, List.length_append]
    have hN : N.length = 1 := by have := List.length_pos_of_ne_nil h3; omega
    obtain ⟨a', N', rfl⟩ := List.exists_cons_of_ne_nil h3
    have hN' : N' = [] := List.eq_nil_of_length_eq_zero (by simp only [List.length_cons] at hN; omega)
    subst hN'
    rcases Bt with _ | ⟨b₀, Bt'⟩
    · have := congrArg List.length hts
      simp only [List.length_nil, List.length_cons] at hl this; omega
    · have hbr : b₀ = b ∧ Bt' = rest := by
        have := h1.symm.trans hts; simp only [List.singleton_append, List.cons.injEq] at this
        exact this.2
      obtain ⟨rfl, rfl⟩ := hbr
      have h₁ := h.unwrap_absorb hg hs hty hok
        (by rw [h1]; simp only [List.singleton_append, List.length_cons]; omega)
      have h₂ := h₁.mergeTop (by rw [hts']; simp only [List.length_cons]; omega)
      refine ⟨Bh ++ [b₀], Bt', by simp, fun u hu => ?_, h₂.finishTop hg₂ hf.merge hts₂ hside⟩
      rw [List.mem_append, List.mem_singleton] at hu
      rcases hu with hu | rfl
      · exact Or.inr hu
      · exact Or.inl hnxt

end Dec

/-! ### `finishEdge_ownedD`: assembly from a stack decomposition -/

theorem drop_append_len {α : Type} {A B : List α} {m : Nat} (hm : m ≤ B.length) :
    (A ++ B).drop ((A ++ B).length - m) = B.drop (B.length - m) := by
  rw [List.length_append, show A.length + B.length - m = A.length + (B.length - m) by omega,
    List.drop_append, List.drop_eq_nil_of_le (by omega), Nat.add_sub_cancel_left,
    List.nil_append]

theorem take_append_len {α : Type} {A B : List α} {m : Nat} (hm : m ≤ B.length) :
    (A ++ B).take ((A ++ B).length - m) = A ++ B.take (B.length - m) := by
  rw [List.length_append, show A.length + B.length - m = A.length + (B.length - m) by omega,
    List.take_append, List.take_of_length_le (by omega), Nat.add_sub_cancel_left]

theorem le_length_of_disjoint_suffix {α : Type} {Bh Bt : List α} {p : α → Prop} {m : Nat}
    (hp : ∀ b ∈ Bh, p b) (hs : ∀ b ∈ (Bh ++ Bt).drop ((Bh ++ Bt).length - m), ¬ p b)
    (hm : m ≤ (Bh ++ Bt).length) : m ≤ Bt.length := by
  by_contra hlt
  have hl : (Bh ++ Bt).length = Bh.length + Bt.length := List.length_append
  have hi : (Bh ++ Bt).length - m < Bh.length := by omega
  have hiB : (Bh ++ Bt).length - m < (Bh ++ Bt).length := by omega
  have hmem : (Bh ++ Bt)[(Bh ++ Bt).length - m] ∈ (Bh ++ Bt).drop ((Bh ++ Bt).length - m) := by
    rw [List.drop_eq_getElem_cons hiB]; exact List.mem_cons_self
  apply hs _ hmem
  rw [List.getElem_append_left hi]
  exact hp _ (List.getElem_mem hi)

/-- `OwnedD` after a returning edge, from the frame `OFrame`, the stack split `sub ++ (Bh ++ Bt)`
(`sub` popped, `Bh` absorbed into the new entries `N`, `Bt` untouched), the sources of the new
entries' edges and their coverage of the popped/absorbed/returned edges. -/
theorem OwnedD.of_split {σ sts origs : List Nat} {P : ItemId → Prop} {d n : Nat} {s s' : WalkState}
    {curV : Nat} {o : DfsOut} {sub Bh Bt N : List TEntry}
    (ho : OwnedD σ sts origs P d n s)
    (hg : s'.g = s.g) (hsv : s'.stackVerts = s.stackVerts)
    (hch : ∀ p, p < 1 + s.g.nv + s.g.ne → Items.ch s'.items p = Items.ch s.items p)
    (hbl : ∀ a, a < 1 + s.g.nv + s.g.ne → ∀ i, Items.Below s'.items a i ↔ Items.Below s.items a i)
    (hvs : ∀ t ∈ s'.tstack, P (vertItem t.vStart) ∨ t.vStart = curV ∨
      ∃ u ∈ s.tstack, u.vStart = t.vStart)
    (hsvlt : ∀ k, k ≤ d → s.stackVerts[k]! < s.g.nv) (hcur : s.stackVerts[d]! = curV)
    (hsts : ∀ k, k ≤ d → sts[k]! ≤ sts[d]!) (hn : sts[d]! ≤ n)
    (horigs : ∀ k, k ≤ d → origs[k]! ≤ origs[d]!) (horig : origs[d]! ≤ Bh.length + Bt.length)
    (hts : s.tstack = sub ++ (Bh ++ Bt)) (hts' : s'.tstack = N ++ Bt)
    (hBh : ∀ b ∈ Bh, b.vStart = curV)
    (hBt : ∀ t ∈ Bt, ∀ e, e < s.g.ne → (t.edges s'.g s'.items e ↔ t.edges s.g s.items e))
    (hsrc : ∀ t ∈ N, ∀ e, e < s.g.ne → t.edges s'.g s'.items e →
      (∃ u ∈ sub, u.edges s.g s.items e) ∨ e = o.e ∨ (∃ b ∈ Bh, b.edges s.g s.items e) ∨
        Items.EdgeBelow s.g s.items (vertItem curV) e)
    (hcov : ∀ e, e < s.g.ne →
      ((∃ u ∈ sub, u.edges s.g s.items e) ∨ e = o.e ∨ ∃ b ∈ Bh, b.edges s.g s.items e) →
      ∃ t ∈ N, t.edges s'.g s'.items e)
    (hnd : σ.Nodup) (hσ : ∀ e ∈ σ, e < s.g.ne) (hpos : σ[n]? = some o.e) :
    OwnedD σ sts origs P d (n + 1) s' := by
  obtain ⟨hnlen, hσn⟩ := List.getElem?_eq_some_iff.1 hpos
  have hσn! : σ[n]! = o.e := by rw [getElem!_pos σ n hnlen]; exact hσn
  have hidx : σ.idxOf o.e = n := by rw [← hσn!]; exact idxOf_getElem! hnd hnlen
  have hoe : o.e < s.g.ne := hσ _ (hσn ▸ List.getElem_mem hnlen)
  have horig' : origs[d]! ≤ (Bh ++ Bt).length := by rw [List.length_append]; exact horig
  have hBl : origs[d]! ≤ Bt.length :=
    le_length_of_disjoint_suffix (p := fun b => b.vStart = curV) hBh
      (fun b hb => by
        have := ho.old d (Nat.le_refl _) b (by rw [hts, drop_append_len horig']; exact hb)
        rwa [hcur] at this) horig'
  have hlen' : s'.tstack.length = N.length + Bt.length := by rw [hts', List.length_append]
  have hvi : ∀ k, k ≤ d → vertItem s.stackVerts[k]! < 1 + s.g.nv + s.g.ne := fun k hk => by
    have := hsvlt k hk
    show 1 + _ < _
    exact Nat.lt_of_lt_of_le (Nat.add_lt_add_left this 1) (Nat.le_add_right _ _)
  have hbelow : ∀ k, k ≤ d → ∀ e, Items.EdgeBelow s'.g s'.items (vertItem s'.stackVerts[k]!) e ↔
      Items.EdgeBelow s.g s.items (vertItem s.stackVerts[k]!) e := fun k hk e => by
    rw [hsv, hg]; exact hbl _ (hvi k hk) _
  have htake : ∀ k, k ≤ d → s.tstack.take (s.tstack.length - origs[k]!) =
      sub ++ (Bh ++ Bt.take (Bt.length - origs[k]!)) := fun k hk => by
    have h1 : origs[k]! ≤ Bt.length := Nat.le_trans (horigs k hk) hBl
    rw [hts, take_append_len (by rw [List.length_append]; exact Nat.le_trans h1 (Nat.le_add_left _ _)),
      take_append_len h1]
  have hdrop : ∀ k, k ≤ d → s.tstack.drop (s.tstack.length - origs[k]!) =
      Bt.drop (Bt.length - origs[k]!) := fun k hk => by
    have h1 : origs[k]! ≤ Bt.length := Nat.le_trans (horigs k hk) hBl
    rw [hts, drop_append_len (by rw [List.length_append]; exact Nat.le_trans h1 (Nat.le_add_left _ _)),
      drop_append_len h1]
  have htake' : ∀ k, k ≤ d → s'.tstack.take (s'.tstack.length - origs[k]!) =
      N ++ Bt.take (Bt.length - origs[k]!) := fun k hk => by
    rw [hts', take_append_len (Nat.le_trans (horigs k hk) hBl)]
  have hdrop' : ∀ k, k ≤ d → s'.tstack.drop (s'.tstack.length - origs[k]!) =
      Bt.drop (Bt.length - origs[k]!) := fun k hk => by
    rw [hts', drop_append_len (Nat.le_trans (horigs k hk) hBl)]
  have hsrc_lo : ∀ k, k ≤ d → ∀ t ∈ N, ∀ e, e < s.g.ne → t.edges s'.g s'.items e →
      sts[k]! ≤ σ.idxOf e := by
    intro k hk t ht e he hte
    rcases hsrc t ht e he hte with ⟨u, hu, hue⟩ | rfl | ⟨b, hb, hbe⟩ | hv
    · exact ho.new k hk u (by rw [htake k hk]; exact List.mem_append_left _ hu) e he hue
    · rw [hidx]; exact Nat.le_trans (hsts k hk) hn
    · exact ho.new k hk b
        (by rw [htake k hk]; exact List.mem_append_right _ (List.mem_append_left _ hb)) e he hbe
    · exact Nat.le_trans (hsts k hk) (ho.lo d (Nat.le_refl _) e he (hcur ▸ hv))
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · intro k hk
    rw [hlen']
    exact Nat.le_trans (horigs k hk) (Nat.le_trans hBl (Nat.le_add_left _ _))
  · intro b hb0 hb
    rcases Nat.lt_or_ge b n with hbn | hbn
    · have hbe : σ[b]! < s.g.ne := getElem!_lt hσ (Nat.lt_trans hbn hnlen)
      rcases ho.cover b hb0 hbn with ⟨t, ht, hte⟩ | ⟨k, hk, hke⟩
      · rw [hts] at ht
        rcases List.mem_append.1 ht with hu | hb'
        · obtain ⟨t', ht', h'⟩ := hcov _ hbe (Or.inl ⟨t, hu, hte⟩)
          exact Or.inl ⟨t', hts' ▸ List.mem_append_left _ ht', h'⟩
        · rcases List.mem_append.1 hb' with hbh | hbt
          · obtain ⟨t', ht', h'⟩ := hcov _ hbe (Or.inr (Or.inr ⟨t, hbh, hte⟩))
            exact Or.inl ⟨t', hts' ▸ List.mem_append_left _ ht', h'⟩
          · exact Or.inl ⟨t, hts' ▸ List.mem_append_right _ hbt, (hBt t hbt _ hbe).2 hte⟩
      · exact Or.inr ⟨k, hk, (hbelow k hk _).2 hke⟩
    · obtain rfl : n = b := by omega
      obtain ⟨t', ht', h'⟩ := hcov _ hoe (Or.inr (Or.inl rfl))
      exact Or.inl ⟨t', hts' ▸ List.mem_append_left _ ht', hσn! ▸ h'⟩
  · intro k hk e he hke
    exact ho.anc k hk e (hg ▸ he) ((hbelow k (Nat.le_of_lt hk) e).1 hke)
  · intro e he hke
    exact Nat.lt_succ_of_lt (ho.hi e (hg ▸ he) ((hbelow d (Nat.le_refl _) e).1 hke))
  · intro k hk e he hke
    exact ho.lo k hk e (hg ▸ he) ((hbelow k hk e).1 hke)
  · intro k hk t ht e he hte
    rw [htake' k hk] at ht
    have he' : e < s.g.ne := hg ▸ he
    rcases List.mem_append.1 ht with hN | hbt
    · exact hsrc_lo k hk t hN e he' hte
    · have hbt' : t ∈ Bt := List.mem_of_mem_take hbt
      exact ho.new k hk t
        (by rw [htake k hk]; exact List.mem_append_right _ (List.mem_append_right _ hbt)) e he'
        ((hBt t hbt' e he').1 hte)
  · intro k hk t ht
    rw [hdrop' k hk] at ht
    rw [hsv]
    exact ho.old k hk t (by rw [hdrop k hk]; exact ht)
  · intro t ht
    rcases hvs t ht with hp | hc | ⟨u, hu, hut⟩
    · exact Or.inl hp
    · exact Or.inr ⟨d, Nat.le_refl _, by rw [hsv, hcur]; exact hc⟩
    · rcases ho.vis u hu with hp | ⟨k, hk, hku⟩
      · exact Or.inl (hut ▸ hp)
      · exact Or.inr ⟨k, hk, by rw [hsv, ← hut, hku]⟩
  · intro w hw hP hne
    rw [hch _ (by rw [hg] at hw; show 1 + w < _; omega)]
    exact ho.fresh w (hg ▸ hw) hP (fun k hk => by have := hne k hk; rwa [hsv] at this)

/-- `OwnedD` after a boundary return: the popped entries `sub` and the returned edge move under
`vertItem curV`, the base `base` and all other vertex items are untouched. -/
theorem OwnedD.of_boundary {σ sts origs : List Nat} {P : ItemId → Prop} {d n : Nat}
    {s s' : WalkState} {curV : Nat} {o : DfsOut} {sub base : List TEntry}
    (ho : OwnedD σ sts origs P d n s)
    (hg : s'.g = s.g) (hsv : s'.stackVerts = s.stackVerts)
    (hsvlt : ∀ k, k ≤ d → s.stackVerts[k]! < s.g.nv) (hcur : s.stackVerts[d]! = curV)
    (hpath : ∀ k, k < d → s.stackVerts[k]! ≠ curV)
    (hsts : ∀ k, k ≤ d → sts[k]! ≤ sts[d]!) (hn : sts[d]! ≤ n)
    (horigs : ∀ k, k ≤ d → origs[k]! ≤ origs[d]!) (horig : origs[d]! ≤ base.length)
    (hts : s.tstack = sub ++ base) (hts' : s'.tstack = base)
    (hch : ∀ w, w < s.g.nv → w ≠ curV →
      Items.ch s'.items (vertItem w) = Items.ch s.items (vertItem w))
    (hbl : ∀ w, w < s.g.nv → w ≠ curV → ∀ i,
      Items.Below s'.items (vertItem w) i ↔ Items.Below s.items (vertItem w) i)
    (hbase : ∀ t ∈ base, ∀ e, e < s.g.ne → (t.edges s'.g s'.items e ↔ t.edges s.g s.items e))
    (hV : ∀ e, e < s.g.ne → (Items.EdgeBelow s'.g s'.items (vertItem curV) e ↔
      Items.EdgeBelow s.g s.items (vertItem curV) e ∨ e = o.e ∨ ∃ t ∈ sub, t.edges s.g s.items e))
    (hproc : ∀ t ∈ s.tstack, ∀ e, e < s.g.ne → t.edges s.g s.items e → σ.idxOf e < n)
    (hnd : σ.Nodup) (hσ : ∀ e ∈ σ, e < s.g.ne) (hpos : σ[n]? = some o.e) :
    OwnedD σ sts origs P d (n + 1) s' := by
  obtain ⟨hnlen, hσn⟩ := List.getElem?_eq_some_iff.1 hpos
  have hσn! : σ[n]! = o.e := by rw [getElem!_pos σ n hnlen]; exact hσn
  have hidx : σ.idxOf o.e = n := by rw [← hσn!]; exact idxOf_getElem! hnd hnlen
  have hoe : o.e < s.g.ne := hσ _ (hσn ▸ List.getElem_mem hnlen)
  have hbelow : ∀ k, k < d → ∀ e, Items.EdgeBelow s'.g s'.items (vertItem s'.stackVerts[k]!) e ↔
      Items.EdgeBelow s.g s.items (vertItem s.stackVerts[k]!) e := fun k hk e => by
    rw [hsv, hg]; exact hbl _ (hsvlt k (Nat.le_of_lt hk)) (hpath k hk) _
  have hVd : ∀ e, e < s.g.ne → (Items.EdgeBelow s'.g s'.items (vertItem s'.stackVerts[d]!) e ↔
      Items.EdgeBelow s.g s.items (vertItem s.stackVerts[d]!) e ∨ e = o.e ∨
        ∃ t ∈ sub, t.edges s.g s.items e) :=
    fun e he => by rw [hsv, hcur]; exact hV e he
  have htake : ∀ k, k ≤ d → s.tstack.take (s.tstack.length - origs[k]!) =
      sub ++ base.take (base.length - origs[k]!) :=
    fun k hk => by rw [hts, take_append_len (Nat.le_trans (horigs k hk) horig)]
  have hdrop : ∀ k, k ≤ d → s.tstack.drop (s.tstack.length - origs[k]!) =
      base.drop (base.length - origs[k]!) :=
    fun k hk => by rw [hts, drop_append_len (Nat.le_trans (horigs k hk) horig)]
  have hsub_lo : ∀ k, k ≤ d → ∀ t ∈ sub, ∀ e, e < s.g.ne → t.edges s.g s.items e →
      sts[k]! ≤ σ.idxOf e :=
    fun k hk t ht e he hte =>
      ho.new k hk t (by rw [htake k hk]; exact List.mem_append_left _ ht) e he hte
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · intro k hk; rw [hts']; exact Nat.le_trans (horigs k hk) horig
  · intro b hb0 hb
    rcases Nat.lt_or_ge b n with hbn | hbn
    · have hbe : σ[b]! < s.g.ne := getElem!_lt hσ (Nat.lt_trans hbn hnlen)
      rcases ho.cover b hb0 hbn with ⟨t, ht, hte⟩ | ⟨k, hk, hke⟩
      · rw [hts] at ht
        rcases List.mem_append.1 ht with hu | hb'
        · exact Or.inr ⟨d, Nat.le_refl _, (hVd _ hbe).2 (Or.inr (Or.inr ⟨t, hu, hte⟩))⟩
        · exact Or.inl ⟨t, hts' ▸ hb', (hbase t hb' _ hbe).2 hte⟩
      · rcases Nat.lt_or_ge k d with hkd | hkd
        · exact Or.inr ⟨k, hk, (hbelow k hkd _).2 hke⟩
        · obtain rfl : d = k := by omega
          exact Or.inr ⟨d, Nat.le_refl _, (hVd _ hbe).2 (Or.inl hke)⟩
    · obtain rfl : n = b := by omega
      rw [hσn!]
      exact Or.inr ⟨d, Nat.le_refl _, (hVd _ hoe).2 (Or.inr (Or.inl rfl))⟩
  · intro k hk e he hke
    exact ho.anc k hk e (hg ▸ he) ((hbelow k hk e).1 hke)
  · intro e he hke
    have he' : e < s.g.ne := hg ▸ he
    rcases (hVd e he').1 hke with h | rfl | ⟨t, ht, hte⟩
    · exact Nat.lt_succ_of_lt (ho.hi e he' h)
    · rw [hidx]; exact Nat.lt_succ_self _
    · exact Nat.lt_succ_of_lt (hproc t (by rw [hts]; exact List.mem_append_left _ ht) e he' hte)
  · intro k hk e he hke
    have he' : e < s.g.ne := hg ▸ he
    rcases Nat.lt_or_ge k d with hkd | hkd
    · exact ho.lo k hk e he' ((hbelow k hkd e).1 hke)
    · obtain rfl : d = k := by omega
      rcases (hVd e he').1 hke with h | rfl | ⟨t, ht, hte⟩
      · exact ho.lo d (Nat.le_refl _) e he' h
      · rw [hidx]; exact hn
      · exact hsub_lo d (Nat.le_refl _) t ht e he' hte
  · intro k hk t ht e he hte
    rw [hts'] at ht
    have hbt : t ∈ base := List.mem_of_mem_take ht
    have he' : e < s.g.ne := hg ▸ he
    exact ho.new k hk t (by rw [htake k hk]; exact List.mem_append_right _ ht) e he'
      ((hbase t hbt e he').1 hte)
  · intro k hk t ht
    rw [hts'] at ht
    rw [hsv]
    exact ho.old k hk t (by rw [hdrop k hk]; exact ht)
  · intro t ht
    rw [hts'] at ht
    rcases ho.vis t (by rw [hts]; exact List.mem_append_right _ ht) with hp | ⟨k, hk, hk'⟩
    · exact Or.inl hp
    · exact Or.inr ⟨k, hk, by rw [hsv]; exact hk'⟩
  · intro w hw hP hne
    have hw' : w < s.g.nv := hg ▸ hw
    have hwc : w ≠ curV := fun h => hne d (Nat.le_refl _) (by rw [hsv, hcur]; exact h)
    rw [hch w hw' hwc]
    exact ho.fresh w hw' hP (fun k hk => by have := hne k hk; rwa [hsv] at this)

namespace Dec
variable {s₂ : WalkState} {N₀ Bh Bt : List TEntry} {curV : Nat}

theorem finishP (h : Dec s₂ N₀ Bh Bt curV s) (hg : s.g = s₂.g) (hi : s.Inv' D) (hs : Shape s)
    {lowval : Nat} {isType1 : Bool} (hok : FinishPOk D curV lowval isType1 s) :
    ∃ Bh' Bt', Bh ++ Bt = Bh' ++ Bt' ∧ (∀ b ∈ Bh', b.vStart = curV ∨ b ∈ Bh) ∧
      Dec s₂ N₀ Bh' Bt' curV (after (Spqr.finishP curV lowval isType1) s) := by
  show ∃ Bh' Bt', _ ∧ _ ∧ Dec s₂ N₀ Bh' Bt' curV ((Spqr.finishP curV lowval isType1).run s).2
  simp only [Spqr.finishP, WalkM.run_bind, run_condP]
  by_cases hc : result (condP curV lowval isType1) s = true
  · have hc' : (isType1 && decide (s.tstack.length ≥ 2) && (s.tstack.tail.head!.vStart == curV) &&
        (s.tstack.tail.head!.topDepth == lowval)) = true := hc
    simp only [hc', ↓reduceIte, WalkM.run_bind]
    obtain ⟨hu, hcl⟩ := hok.ok hc
    have r := maybeUnwrapNxt_spec (v := curV) hi hs (by decide) hu
    have hnxt : s.tstack.tail.head!.vStart = curV := by
      simp only [Bool.and_eq_true, beq_iff_eq] at hc'
      exact hc'.1.2
    exact h.closeTwo hg hs (by decide) hu r.free hcl (Or.inr hnxt)
  · have hc' : (isType1 && decide (s.tstack.length ≥ 2) && (s.tstack.tail.head!.vStart == curV) &&
        (s.tstack.tail.head!.topDepth == lowval)) = false := Bool.eq_false_iff.2 hc
    simp only [hc', Bool.false_eq_true, ↓reduceIte]
    exact ⟨Bh, Bt, rfl, fun b hb => Or.inr hb, h⟩

theorem finishRest (h : Dec s₂ N₀ Bh Bt curV s) (hg : s.g = s₂.g) (hi : s.Inv' D) (hs : Shape s)
    (hv : curV < s.g.nv)
    (hbl : ∀ i, Items.Below s.items (vertItem curV) i ↔ Items.Below s₂.items (vertItem curV) i)
    {d lowval : Nat} {isType1 hasVert isSingle : Bool}
    (hok : FinishRestOk D curV d lowval isType1 hasVert isSingle s) :
    ∃ Bh' Bt', Bh ++ Bt = Bh' ++ Bt' ∧ (∀ b ∈ Bh', b.vStart = curV ∨ b ∈ Bh) ∧
      Dec s₂ N₀ Bh' Bt' curV (after (Spqr.finishRest curV d lowval isType1 hasVert isSingle) s) := by
  obtain ⟨Bh', Bt', heq, hb, h₁⟩ := h.finishP hg hi hs hok.p
  have st₁ : Step D curV s (after (Spqr.finishP curV lowval isType1) s) := Step.finishP hi hs hv hok.p
  have hbl₁ : ∀ i, Items.Below (after (Spqr.finishP curV lowval isType1) s).items (vertItem curV) i ↔
      Items.Below s₂.items (vertItem curV) i := fun i => (st₁.below _).trans (hbl i)
  refine ⟨Bh', Bt', heq, hb, ?_⟩
  show Dec s₂ N₀ Bh' Bt' curV ((Spqr.finishRest curV d lowval isType1 hasVert isSingle).run s).2
  simp only [Spqr.finishRest, Spqr.finishTail, WalkM.run_bind]
  cases hasVert
  · simp only [Bool.not_false, ↓reduceIte, WalkM.run_bind]
    cases isSingle
    · simp only [Bool.not_false, ↓reduceIte, WalkM.run_bind]
      exact h₁.pushVert_merge d hbl₁
    · simp only [Bool.not_true, Bool.false_eq_true, ↓reduceIte]
      exact h₁.pushVert d hbl₁
  · simp only [Bool.not_true, Bool.false_eq_true, ↓reduceIte]
    exact h₁

end Dec

theorem length_after_mergeTop {a b : TEntry} {rest : List TEntry} (hts : s.tstack = a :: b :: rest) :
    (after mergeTstackTops s).tstack.length = s.tstack.length - 1 := by
  rw [after_mergeTstackTops_eq hts, hts]; simp

theorem loop3_length_le (orig fuel : Nat) (s : WalkState) (hfuel : s.tstack.length ≤ fuel + orig + 3) :
    (after (WalkM.loop fuel (loop3Cond orig) mergeTstackTops) s).tstack.length ≤ orig + 3 := by
  induction fuel generalizing s with
  | zero => show s.tstack.length ≤ orig + 3; omega
  | succ fuel ih =>
    show ((WalkM.loop (fuel + 1) (loop3Cond orig) mergeTstackTops).run s).2.tstack.length ≤ orig + 3
    rw [loop_succ_run fuel _ _ s rfl, run_loop3Cond]
    by_cases hc : s.tstack.length > orig + 3
    · rw [if_pos (by simpa using hc)]
      obtain ⟨a, b, rest, hts⟩ := two_entries_of_le (l := s.tstack) (by omega)
      exact ih (after mergeTstackTops s) (by rw [length_after_mergeTop hts]; omega)
    · rw [if_neg (by simpa using hc)]
      exact Nat.le_of_not_lt hc

namespace Dec
variable {s₂ : WalkState} {N₀ Bh Bt : List TEntry} {curV : Nat}

/-- The type-1/type-2 vertex close on a stack of exactly `Bt.length + 3` entries above the kept
suffix (`mid = []` for type 1; the type-2 loop stops at `orig + 3`). -/
theorem closeVert' (h : Dec s₂ N₀ Bh Bt curV s) (hg : s.g = s₂.g) (hi : s.Inv' D) (hs : Shape s)
    (hv : curV < s.g.nv) {dir isType1 isSingle : Bool}
    (hlen : Bt.length + 3 ≤ s.tstack.length) (h1 : isType1 = true → s.tstack.length = Bt.length + 3)
    (hok : CloseVertOk D curV dir isType1 Bt.length isSingle s) :
    Dec s₂ N₀ Bh Bt curV (after (closeVert' curV dir isType1 Bt.length isSingle) s) := by
  have st₁ : Step D curV s (cvS₁ isType1 Bt.length isSingle s) := Step.vertPre hi hs hv hok.loop3
  obtain ⟨st₂, hfree⟩ := Step.vertUnwrap (v := curV) st₁.inv st₁.shape (isType1 := isType1)
    (isSingle := cvB₁ isType1 Bt.length isSingle s) (fun h => by subst h; exact hok.unwrap rfl)
  have st₂ : Step D curV (cvS₁ isType1 Bt.length isSingle s) (cvS₂ isType1 Bt.length isSingle s) := st₂
  have key : Dec s₂ N₀ Bh Bt curV (cvS₂ isType1 Bt.length isSingle s) →
      (cvS₂ isType1 Bt.length isSingle s).tstack.length = Bt.length + 3 →
      Dec s₂ N₀ Bh Bt curV (cvS₅ curV dir isType1 Bt.length isSingle s) ∧
        (cvS₅ curV dir isType1 Bt.length isSingle s).g = s₂.g := by
    intro hd hl₂
    obtain ⟨a, b, rest, hts₂⟩ := two_entries_of_le (l := (cvS₂ isType1 Bt.length isSingle s).tstack) (by omega)
    have h₃ := hd.mergeTop (by omega)
    have hl₃ : (cvS₃ isType1 Bt.length isSingle s).tstack.length = Bt.length + 2 := by
      show (after mergeTstackTops _).tstack.length = _
      rw [length_after_mergeTop hts₂, hl₂]; rfl
    obtain ⟨a', b', rest', hts₃⟩ := two_entries_of_le (l := (cvS₃ isType1 Bt.length isSingle s).tstack) (by omega)
    have h₄ := h₃.mergeTop (show Bt.length + 2 ≤ (cvS₃ isType1 Bt.length isSingle s).tstack.length by omega)
    have hl₄ : (cvS₄ isType1 Bt.length isSingle s).tstack.length = Bt.length + 1 := by
      show (after mergeTstackTops _).tstack.length = _
      rw [length_after_mergeTop hts₃, hl₃]; rfl
    obtain ⟨t₄, rest₄, hts₄⟩ := List.exists_cons_of_ne_nil (l := (cvS₄ isType1 Bt.length isSingle s).tstack)
      (by intro h0; rw [h0] at hl₄; simp at hl₄)
    refine ⟨h₄.retarget dir hts₄, ?_⟩
    show ((WalkState.retarget curV dir).run _).2.g = _
    rw [retarget_run_eq curV dir _ t₄ rest₄ hts₄]
    show (after mergeTstackTops (after mergeTstackTops _)).g = _
    rw [after_mergeTstackTops, after_mergeTstackTops]
    exact st₂.g.trans (st₁.g.trans hg)
  show Dec s₂ N₀ Bh Bt curV (after (vertFinish (result (vertUnwrap isType1 (cvB₁ isType1 Bt.length isSingle s))
    (cvS₁ isType1 Bt.length isSingle s)) (cvB₁ isType1 Bt.length isSingle s))
    (cvS₅ curV dir isType1 Bt.length isSingle s))
  cases isType1
  · have hS₂ : cvS₂ false Bt.length isSingle s = cvS₁ false Bt.length isSingle s := rfl
    have hS₁ : cvS₁ false Bt.length isSingle s =
        after (WalkM.loop s.tstack.length (loop3Cond Bt.length) mergeTstackTops) s := rfl
    have h₁ : Dec s₂ N₀ Bh Bt curV (cvS₁ false Bt.length isSingle s) ∧
        Bt.length + 3 ≤ (cvS₁ false Bt.length isSingle s).tstack.length := by
      rw [hS₁]
      refine loop_pred (fun s' => Dec s₂ N₀ Bh Bt curV s' ∧ Bt.length + 3 ≤ s'.tstack.length) _ _ _
        (fun _ => rfl) (fun k hk hd => ?_) ⟨h, hlen⟩
      have hc := hk k (Nat.le_refl _)
      rw [run_loop3Cond] at hc
      simp only [decide_eq_true_eq] at hc
      obtain ⟨a, b, rest, hts⟩ := two_entries_of_le (l := (iter mergeTstackTops k s).tstack) (by omega)
      exact ⟨hd.1.mergeTop (by omega), by rw [length_after_mergeTop hts]; omega⟩
    have hle := loop3_length_le Bt.length s.tstack.length s (by omega)
    have hl₂ : (cvS₂ false Bt.length isSingle s).tstack.length = Bt.length + 3 := by
      rw [hS₂, hS₁]; rw [hS₁] at h₁; omega
    exact (key (hS₂ ▸ h₁.1) hl₂).1
  · have hS₂ : cvS₂ true Bt.length isSingle s =
        after (maybeUnwrapNxt (if isSingle then .S else .R)) (cvS₁ true Bt.length isSingle s) := by
      simp only [cvS₂, after, WalkState.vertUnwrap, ↓reduceIte, WalkM.map_run]; rfl
    have hu : UnwrapOk (if isSingle then .S else .R) (cvS₁ true Bt.length isSingle s) := hok.unwrap rfl
    have hty : (if isSingle then NodeType.S else .R) ∉ [NodeType.F, .V, .Q] := by cases isSingle <;> decide
    obtain ⟨a, b, rest, hts₁⟩ := two_entries_of_le hu.two
    have hl₁ := h1 rfl
    have hd₂ : Dec s₂ N₀ Bh Bt curV (cvS₂ true Bt.length isSingle s) := by
      rw [hS₂]
      exact (show Dec s₂ N₀ Bh Bt curV (cvS₁ true Bt.length isSingle s) from h).unwrap hg hs hty hu
        (show Bt.length + 2 ≤ s.tstack.length by omega)
    obtain ⟨-, b', hts₂, -, -, -⟩ := unwrap_edges hs hty hu hts₁
    have hl₂ : (cvS₂ true Bt.length isSingle s).tstack.length = Bt.length + 3 := by
      rw [hS₂]
      show (after (maybeUnwrapNxt _) s).tstack.length = _
      rw [hts₂]
      have hl₁' : (a :: b :: rest).length = Bt.length + 3 := by rw [← hts₁]; exact hl₁
      simpa using hl₁'
    obtain ⟨h₅, hg₅⟩ := key hd₂ hl₂
    have hf : ItemFree (cvS₂ true Bt.length isSingle s) ((maybeUnwrapNxt (if isSingle then .S else .R)).run
        (cvS₁ true Bt.length isSingle s)).1 := hfree _ rfl
    have hf₅ := ((hf.merge).merge).retarget curV dir
    have hfin := hok.finish rfl
    obtain ⟨t, rest, hts₅⟩ := List.exists_cons_of_ne_nil hfin.nonempty
    have hside := hfin.side
    rw [curE, hts₅] at hside
    exact h₅.finishTop hg₅ hf₅ hts₅ hside

end Dec

/-! ### `finishEdge_ownedD`: the step -/

theorem ch_of_size_le {items : Items} {p : ItemId} (h : items.size ≤ p) : Items.ch items p = [] := by
  simp [Items.ch, Array.getElem?_eq_none h]

theorem Below_of_ch_nil {items : Items} {j i : ItemId} (hj : Items.ch items j = []) :
    Items.Below items j i ↔ i = j := by
  constructor
  · intro h
    rcases h.head_cases with h | ⟨c, hc, _⟩
    · exact h.symm
    · simp [Items.IsParent, hj] at hc
  · rintro rfl; exact .refl

theorem Below_of_root {items : Items} {a j : ItemId} (hj : ∀ p, ¬ Items.IsParent items p j) :
    Items.Below items a j ↔ a = j :=
  ⟨Items.Below.eq_of_no_parent hj, fun h => h ▸ .refl⟩

theorem Below_trans {items : Items} {x y z : ItemId} (h₁ : Items.Below items x y)
    (h₂ : Items.Below items y z) : Items.Below items x z := by
  induction h₂ with
  | refl => exact h₁
  | tail _ h ih => exact ih.tail h

theorem Below_modify_ch_append {items : Items} (j : ItemId) (L : List ItemId) (hj : j < items.size)
    {a i : ItemId} :
    Items.Below (items.modify j fun it => { it with ch := it.ch ++ L }) a i ↔
      Items.Below items a i ∨ (Items.Below items a j ∧ ∃ q ∈ L, Items.Below items q i) := by
  have hchj : Items.ch items j = items[j].ch := by simp [Items.ch, Array.getElem?_eq_getElem hj]
  have hmono : ∀ p c, Items.IsParent items p c →
      Items.IsParent (items.modify j fun it => { it with ch := it.ch ++ L }) p c := by
    intro p c h
    by_cases hp : p = j
    · subst hp
      rw [Items.IsParent, Items.ch_modify_at p _ hj]
      rw [Items.IsParent, hchj] at h
      exact List.mem_append_left _ h
    · rwa [Items.IsParent, Items.ch_modify_of_ne j _ hp]
  have hmono' : ∀ a i, Items.Below items a i →
      Items.Below (items.modify j fun it => { it with ch := it.ch ++ L }) a i := by
    intro a i h
    induction h with
    | refl => exact .refl
    | @tail b c _ hbc ih => exact ih.tail (hmono b c hbc)
  constructor
  · intro h
    induction h with
    | refl => exact Or.inl .refl
    | @tail b c _ hbc ih =>
      rcases Items.IsParent_modify hbc with hbc | ⟨hbj, _, hc⟩
      · rcases ih with ih | ⟨haj, q, hq, hqb⟩
        · exact Or.inl (ih.tail hbc)
        · exact Or.inr ⟨haj, q, hq, hqb.tail hbc⟩
      · subst hbj
        simp only [List.mem_append] at hc
        have hab : Items.Below items a b := by
          rcases ih with ih | ⟨haj, _⟩
          · exact ih
          · exact haj
        rcases hc with hc | hc
        · exact Or.inl (hab.tail (by rw [Items.IsParent, hchj]; exact hc))
        · exact Or.inr ⟨hab, c, hc, .refl⟩
  · rintro (h | ⟨haj, q, hq, hqi⟩)
    · exact hmono' _ _ h
    · refine Below_trans (Below_trans (hmono' _ _ haj) (.tail .refl ?_)) (hmono' _ _ hqi)
      rw [Items.IsParent, Items.ch_modify_at j _ hj]
      exact List.mem_append_right _ hq

theorem Below_modify_ch_set {items : Items} (j : ItemId) (L : List ItemId) (hj : j < items.size)
    (h0 : Items.ch items j = []) {a i : ItemId} :
    Items.Below (items.modify j fun it => { it with ch := L }) a i ↔
      Items.Below items a i ∨ (Items.Below items a j ∧ ∃ q ∈ L, Items.Below items q i) := by
  rw [← Below_modify_ch_append j L hj]
  apply Items.Below_congr
  intro p
  by_cases hp : p = j
  · subst hp
    rw [Items.ch_modify_at p _ hj, Items.ch_modify_at p _ hj]
    have : items[p].ch = [] := by
      rw [← h0]; simp [Items.ch, Array.getElem?_eq_getElem hj]
    simp [this]
  · rw [Items.ch_modify_of_ne j _ hp, Items.ch_modify_of_ne j _ hp]

/-- `OwnedD.of_boundary` for the item writes of `finishBoundary` on top of `M` (the items after
the `vs` writes and the leaf allocation): `Q` gets the children `L`, then goes under
`vertItem curV`. -/
theorem OwnedD.of_boundary' {σ sts origs : List Nat} {P : ItemId → Prop} {d n : Nat}
    {s s' : WalkState} {curV : Nat} {o : DfsOut} {sub base : List TEntry} {M : Items}
    {L : List ItemId}
    (ho : OwnedD σ sts origs P d n s)
    (hg : s'.g = s.g) (hsv : s'.stackVerts = s.stackVerts)
    (hsvlt : ∀ k, k ≤ d → s.stackVerts[k]! < s.g.nv) (hcur : s.stackVerts[d]! = curV)
    (hpath : ∀ k, k < d → s.stackVerts[k]! ≠ curV)
    (hsts : ∀ k, k ≤ d → sts[k]! ≤ sts[d]!) (hn : sts[d]! ≤ n)
    (horigs : ∀ k, k ≤ d → origs[k]! ≤ origs[d]!) (horig : origs[d]! ≤ base.length)
    (hts : s.tstack = sub ++ base) (hts' : s'.tstack = base)
    (hit : s'.items = (M.modify (edgeItem s.g o.e) fun it => { it with ch := L }).modify
      (vertItem curV) fun it => { it with ch := it.ch ++ [edgeItem s.g o.e] })
    (hsz : s.items.size ≤ M.size)
    (hMch : ∀ p, p < s.items.size → Items.ch M p = Items.ch s.items p)
    (hMch' : ∀ p, s.items.size ≤ p → Items.ch M p = [])
    (hMbl : ∀ a i, Items.Below M a i ↔ Items.Below s.items a i)
    (hshape : Shape s) (hv : curV < s.g.nv) (he : o.e < s.g.ne)
    (hq : Items.ch s.items (edgeItem s.g o.e) = [])
    (hqr : ∀ p, ¬ Items.IsParent s.items p (edgeItem s.g o.e))
    (hvr : ∀ p, ¬ Items.IsParent s.items p (vertItem curV))
    (hqf : ∀ t ∈ base, edgeItem s.g o.e ∉ t.spans.1 ++ t.spans.2)
    (hvf : ∀ t ∈ base, vertItem curV ∉ t.spans.1 ++ t.spans.2)
    (hVL : vertItem curV ∉ L)
    (hL : ∀ e, e < s.g.ne → ((∃ q ∈ L, Items.Below s.items q (edgeItem s.g e)) ↔
      ∃ t ∈ sub, t.edges s.g s.items e))
    (hproc : ∀ t ∈ s.tstack, ∀ e, e < s.g.ne → t.edges s.g s.items e → σ.idxOf e < n)
    (hnd : σ.Nodup) (hσ : ∀ e ∈ σ, e < s.g.ne) (hpos : σ[n]? = some o.e) :
    OwnedD σ sts origs P d (n + 1) s' := by
  set Q := edgeItem s.g o.e with hQ
  set V := vertItem curV with hV
  have hQsz : Q < s.items.size :=
    Nat.lt_of_lt_of_le (by show 1 + s.g.nv + o.e < _; omega) hshape.size
  have hVsz : V < s.items.size := Nat.lt_of_lt_of_le (by show 1 + curV < _; omega) hshape.size
  have hQV : V ≠ Q := by show 1 + curV ≠ 1 + s.g.nv + o.e; omega
  have hMq : Items.ch M Q = [] := by rw [hMch Q hQsz]; exact hq
  have hMqr : ∀ p, ¬ Items.IsParent M p Q := fun p h => by
    rcases Nat.lt_or_ge p s.items.size with hp | hp
    · exact hqr p (by rwa [Items.IsParent, hMch p hp] at h)
    · simp [Items.IsParent, hMch' p hp] at h
  have hb₁ : ∀ a i, Items.Below (M.modify Q fun it => { it with ch := L }) a i ↔
      Items.Below s.items a i ∨ (a = Q ∧ ∃ q ∈ L, Items.Below s.items q i) := by
    intro a i
    rw [Below_modify_ch_set Q L (Nat.lt_of_lt_of_le hQsz hsz) hMq, hMbl, Below_of_root hMqr]
    simp only [hMbl]
  have hb₂ : ∀ a i, Items.Below s'.items a i ↔
      Items.Below (M.modify Q fun it => { it with ch := L }) a i ∨
      (Items.Below (M.modify Q fun it => { it with ch := L }) a V ∧
        Items.Below (M.modify Q fun it => { it with ch := L }) Q i) := by
    intro a i
    rw [hit, Below_modify_ch_append V [Q]
      (by rw [Array.size_modify]; exact Nat.lt_of_lt_of_le hVsz hsz)]
    simp only [List.mem_singleton, exists_eq_left]
  have hb₁V : ∀ a, Items.Below (M.modify Q fun it => { it with ch := L }) a V ↔ a = V := by
    intro a
    rw [hb₁, Below_of_root hvr]
    constructor
    · rintro (h | ⟨_, q, hq, hqV⟩)
      · exact h
      · rw [(Below_of_root hvr).1 hqV] at hq
        exact absurd hq hVL
    · intro h; exact Or.inl h
  have hb₁Q : ∀ i, Items.Below (M.modify Q fun it => { it with ch := L }) Q i ↔
      i = Q ∨ ∃ q ∈ L, Items.Below s.items q i := by
    intro i
    rw [hb₁, Below_of_ch_nil hq]
    simp
  have hb : ∀ a i, Items.Below s'.items a i ↔
      (Items.Below s.items a i ∨ (a = Q ∧ ∃ q ∈ L, Items.Below s.items q i)) ∨
      (a = V ∧ (i = Q ∨ ∃ q ∈ L, Items.Below s.items q i)) := by
    intro a i; rw [hb₂, hb₁, hb₁V, hb₁Q]
  refine OwnedD.of_boundary ho hg hsv hsvlt hcur hpath hsts hn horigs horig hts hts'
    ?_ ?_ ?_ ?_ hproc hnd hσ hpos
  · intro w hw hwc
    have hwV : vertItem w ≠ V := fun h => hwc (vertItem_inj h)
    have hwQ : vertItem w ≠ Q := by show 1 + w ≠ 1 + s.g.nv + o.e; omega
    rw [hit, Items.ch_modify_of_ne _ _ hwV, Items.ch_modify_of_ne _ _ hwQ]
    exact hMch _ (Nat.lt_of_lt_of_le (by show 1 + w < _; omega) hshape.size)
  · intro w hw hwc i
    have hwV : vertItem w ≠ V := fun h => hwc (vertItem_inj h)
    have hwQ : vertItem w ≠ Q := by show 1 + w ≠ 1 + s.g.nv + o.e; omega
    rw [hb]
    simp only [hwV, hwQ, false_and, or_false]
  · intro t ht e he
    rw [hg]
    simp only [TEntry.edges, Items.EdgeBelow]
    constructor
    · rintro ⟨i, hi, hb'⟩
      rw [hb] at hb'
      rcases hb' with (hb' | ⟨rfl, -⟩) | ⟨rfl, -⟩
      · exact ⟨i, hi, hb'⟩
      · exact absurd hi (hqf t ht)
      · exact absurd hi (hvf t ht)
    · rintro ⟨i, hi, hb'⟩
      exact ⟨i, hi, (hb _ _).2 (Or.inl (Or.inl hb'))⟩
  · intro e he
    rw [hg]
    simp only [Items.EdgeBelow]
    have hxQ : edgeItem s.g e = Q ↔ e = o.e := ⟨edgeItem_inj, fun h => by rw [h]⟩
    rw [hb, hL e he, hxQ]
    have hVV : vertItem curV = V := rfl
    simp only [hVV, hQV, false_and, or_false, true_and, eq_self_iff_true]

/-- `OwnedD` is preserved by `finishEdge` (the processed prefix grows by `o.e`), given monotone
starts/stack lengths along the path, that the path below `d` avoids `curV`, in-range path
vertices, and that the child of a returning tree edge is visited.
(Checker: `own_*` at `post`.) -/
theorem finishEdge_ownedD {curV d : Nat} {o : DfsOut} {origTstack : Nat} {hasVert : Bool}
    (h : OwnSite σ n D curV d o origTstack hasVert s)
    {sts origs : List Nat} {P : ItemId → Prop} (ho : OwnedD σ sts origs P d n s)
    (hsts : ∀ k, k ≤ d → sts[k]! ≤ sts[d]!) (hn : sts[d]! ≤ n)
    (horigs : ∀ k, k ≤ d → origs[k]! ≤ origs[d]!) (horig : origs[d]! ≤ origTstack)
    (hpath : ∀ k, k < d → s.stackVerts[k]! ≠ curV)
    (hsvlt : ∀ k, k ≤ d → s.stackVerts[k]! < s.g.nv)
    (hdest : o.cls.isTree = true → P (vertItem o.dest)) :
    OwnedD σ sts origs P d (n + 1) (after (finishEdge curV d o origTstack hasVert) s) := by
  obtain ⟨sub, base, hlen, hE⟩ := h.book.ear
  subst hlen
  have hi := h.rgs.1.inv
  have hs := h.rgs.2.1
  have hσ := h.rgs.2.2
  have hv := h.book.v_lt
  have hoe := h.book.e_lt
  have hnd := h.nodup
  have hpos := h.pos
  have hcur := hE.sv_d
  have hproc := h.rgs.1.processed
  have hts := hE.tstack
  have hVsz : vertItem curV < s.items.size :=
    Nat.lt_of_lt_of_le (by show 1 + curV < _; omega) hs.size
  have hchB : ∀ (vs : Option Nat × Option Nat) p,
      Items.ch (s.items.modify (edgeItem s.g o.e) fun it => { it with vs := vs }) p =
        Items.ch s.items p := fun vs p =>
    Items.ch_modify_ch_eq (edgeItem s.g o.e) (fun it => { it with vs := vs }) (fun _ => rfl) p
  have hblB : ∀ (vs : Option Nat × Option Nat) a i,
      Items.Below (s.items.modify (edgeItem s.g o.e) fun it => { it with vs := vs }) a i ↔
        Items.Below s.items a i := fun vs a i =>
    Items.Below_modify_ch_eq (edgeItem s.g o.e) (fun it => { it with vs := vs }) (fun _ => rfl)
  have hch₀ : ∀ p, Items.ch (feS₀ d o s).items p = Items.ch s.items p := fun p => hchB _ p
  have hbl₀ : ∀ a i, Items.Below (feS₀ d o s).items a i ↔ Items.Below s.items a i :=
    fun a i => hblB _ a i
  by_cases hge : d ≤ o.cls.lowval d
  · have hnv := hE.bd_noVert hge
    subst hnv
    have hge' : o.cls.lowval d ≥ d := hge
    have hvf : ∀ t ∈ base, vertItem curV ∉ t.spans.1 ++ t.spans.2 := fun t ht hm =>
      Bool.noConfusion (hE.vert_free t (by rw [hts]; exact List.mem_append_right _ ht) hm).1
    have hqf : ∀ t ∈ base, edgeItem s.g o.e ∉ t.spans.1 ++ t.spans.2 := fun t ht =>
      hE.q_free t (by rw [hts]; exact List.mem_append_right _ ht)
    have hvsub : ∀ t ∈ sub, vertItem curV ∉ t.spans.1 ++ t.spans.2 := fun t ht hm =>
      Bool.noConfusion (hE.vert_free t (by rw [hts]; exact List.mem_append_left _ ht) hm).1
    have hI : ∀ x, x < s.items.size → ¬ Items.Below s.items s.items.size x := fun x hx hb => by
      rcases hb.head_cases with h | ⟨c, hc, _⟩
      · exact absurd h (Nat.ne_of_gt hx)
      · rw [Items.IsParent, ch_of_size_le (Nat.le_refl _)] at hc
        simp at hc
    show wp (finishEdge curV d o base.length false) (fun _ s' => OwnedD σ sts origs P d (n + 1) s') s
    rw [finishEdge_eq]
    simp only [finishEdge', wp_bind, wp_get, wp_stackDir, hge', ↓reduceIte]
    unfold finishBoundary
    simp only [wp_bind, wp_modifyItem, wp_modify, wp_ite, wp_allocItem, wp_makeVs, wp_popTstack,
      wp_pure]
    split
    · rename_i hT
      split
      · rename_i hL
        have hL' : o.cls.lowval d = d + 1 := by simpa using hL
        obtain ⟨t, hsub, -⟩ := hE.bd_bridge hT hL'
        have hts' : s.tstack = t :: base := by rw [hts, hsub]; rfl
        have hside : t.spans.1 = [] := by
          have := hE.bd_side hT hge
          rw [if_pos hL', hts'] at this
          exact this t rfl
        rw [hts']
        simp only [List.head!_cons, List.tail_cons]
        refine OwnedD.of_boundary' (sub := [t]) ho rfl rfl hsvlt hcur hpath hsts hn horigs horig
          hts' rfl rfl ?_ ?_ ?_ ?_ hs hv hoe h.book.q hE.q_root hE.v_root hqf hvf ?_ ?_ hproc hnd
          hσ hpos
        · simp only [Array.size_modify, Array.size_push]; omega
        · intro p hp
          refine (Items.ch_modify_ch_eq _ _ ?_ p).trans ((Items.ch_push_nil _ ?_ p).trans (hchB _ p))
          · exact fun _ => rfl
          · rfl
        · intro p hp
          refine ((Items.ch_modify_ch_eq _ _ ?_ p).trans
            ((Items.ch_push_nil _ ?_ p).trans (hchB _ p))).trans (ch_of_size_le hp)
          · exact fun _ => rfl
          · rfl
        · intro a i
          refine (Items.Below_modify_ch_eq _ _ ?_).trans ((Items.Below_push_nil _ ?_).trans (hblB _ a i))
          · exact fun _ => rfl
          · rfl
        · intro hm
          rcases List.mem_cons.1 hm with hm | hm
          · rw [Array.size_modify] at hm
            exact absurd hm (Nat.ne_of_lt hVsz)
          · exact hvsub t (by rw [hsub]; simp) (List.mem_append_right _ hm)
        · intro e he
          have hx : edgeItem s.g e < s.items.size :=
            Nat.lt_of_lt_of_le (by show 1 + s.g.nv + e < _; omega) hs.size
          constructor
          · rintro ⟨q, hq, hb⟩
            refine ⟨t, List.mem_singleton_self _, ?_⟩
            rcases List.mem_cons.1 hq with rfl | hq
            · rw [Array.size_modify] at hb
              exact absurd hb (hI _ hx)
            · exact ⟨q, by rw [hside]; exact hq, hb⟩
          · rintro ⟨t', ht', q, hq, hb⟩
            rw [List.mem_singleton] at ht'
            subst ht'
            rw [hside, List.nil_append] at hq
            exact ⟨q, List.mem_cons_of_mem _ hq, hb⟩
      · rename_i hL
        have hL' : o.cls.lowval d ≠ d + 1 := by simpa using hL
        obtain ⟨t₁, t₂, hsub, -⟩ := hE.bd_comp hT hge hL'
        have hts' : s.tstack = t₁ :: t₂ :: base := by rw [hts, hsub]; rfl
        obtain ⟨hside₁, hside₂⟩ : t₁.spans.2 = [] ∧ t₂.spans.1 = [] := by
          have := hE.bd_side hT hge
          rw [if_neg hL', hts'] at this
          exact ⟨this.1 t₁ rfl, this.2 t₂ rfl⟩
        rw [hts']
        simp only [List.head!_cons, List.tail_cons]
        refine OwnedD.of_boundary' (sub := [t₁, t₂]) ho rfl rfl hsvlt hcur hpath hsts hn horigs
          horig hts' rfl rfl (le_of_eq (Array.size_modify ..).symm) (fun p _ => hchB _ p)
          (fun p hp => (hchB _ p).trans (ch_of_size_le hp)) (fun a i => hblB _ a i) hs hv hoe
          h.book.q hE.q_root hE.v_root hqf hvf ?_ ?_ hproc hnd hσ hpos
        · intro hm
          rcases List.mem_append.1 hm with hm | hm
          · exact hvsub t₁ (by rw [hsub]; simp) (List.mem_append_left _ hm)
          · exact hvsub t₂ (by rw [hsub]; simp) (List.mem_append_right _ hm)
        · intro e he
          simp only [TEntry.edges, Items.EdgeBelow]
          constructor
          · rintro ⟨q, hq, hb⟩
            rcases List.mem_append.1 hq with hq | hq
            · exact ⟨t₁, by simp, q, by rw [hside₁, List.append_nil]; exact hq, hb⟩
            · exact ⟨t₂, by simp, q, by rw [hside₂, List.nil_append]; exact hq, hb⟩
          · rintro ⟨t', ht', q, hq, hb⟩
            simp only [List.mem_cons, List.not_mem_nil, or_false] at ht'
            rcases ht' with rfl | rfl
            · rw [hside₁, List.append_nil] at hq
              exact ⟨q, List.mem_append_left _ hq, hb⟩
            · rw [hside₂, List.nil_append] at hq
              exact ⟨q, List.mem_append_right _ hq, hb⟩
    · rename_i hT
      have hT' : o.cls.isTree = false := Bool.eq_false_iff.2 hT
      have hsub := hE.back_nil hT'
      subst hsub
      have hts' : s.tstack = base := hts
      refine OwnedD.of_boundary' (sub := []) ho rfl rfl hsvlt hcur hpath hsts hn horigs horig
        hts hts' rfl ?_ ?_ ?_ ?_ hs hv hoe h.book.q hE.q_root hE.v_root hqf hvf ?_ ?_ hproc hnd
        hσ hpos
      · simp only [Array.size_modify, Array.size_push]; omega
      · intro p hp
        refine (Items.ch_modify_ch_eq _ _ ?_ p).trans ((Items.ch_push_nil _ ?_ p).trans (hchB _ p))
        · exact fun _ => rfl
        · rfl
      · intro p hp
        refine ((Items.ch_modify_ch_eq _ _ ?_ p).trans
          ((Items.ch_push_nil _ ?_ p).trans (hchB _ p))).trans (ch_of_size_le hp)
        · exact fun _ => rfl
        · rfl
      · intro a i
        refine (Items.Below_modify_ch_eq _ _ ?_).trans ((Items.Below_push_nil _ ?_).trans (hblB _ a i))
        · exact fun _ => rfl
        · rfl
      · simp only [List.mem_singleton, Array.size_modify]
        exact Nat.ne_of_lt hVsz
      · intro e he
        have hx : edgeItem s.g e < s.items.size :=
          Nat.lt_of_lt_of_le (by show 1 + s.g.nv + e < _; omega) hs.size
        simp only [List.mem_singleton, exists_eq_left, List.not_mem_nil, false_and, exists_false,
          iff_false, Array.size_modify]
        exact hI _ hx
  · have hlow : o.cls.lowval d < d := Nat.lt_of_not_le hge
    obtain ⟨lv, kind, hocls, hl⟩ := ret_of_lowval_lt hlow
    have hok := finishOk_of_guards hocls hl h.guards hE rfl hi hs h.hD hv h.book.e_lt h.book.q
      (h.book.ends lv kind hocls) h.book.vert
    have hlv : o.cls.lowval d = lv := by rw [hocls]; rfl
    have hge' : ¬ (lv ≥ d) := by omega
    have hj : edgeItem s.g o.e < 1 + s.g.nv + s.g.ne := by show 1 + s.g.nv + o.e < _; omega
    have st₀ : Step D curV s (feS₀ d o s) := Step.modifyVs hi hs (edgeItem s.g o.e) _ hj
    have hv₀ : curV < (feS₀ d o s).g.nv := by rw [st₀.g]; exact hv
    show OwnedD σ sts origs P d (n + 1) ((finishEdge curV d o base.length hasVert).run s).2
    rw [finishEdge_eq]
    simp only [finishEdge', finishTree, finishBack, hlv, hge', ↓reduceIte, WalkM.run_bind,
      WalkM.get_run, run_stackDir, run_makeVs, run_modifyItem]
    by_cases ht : o.cls.isTree = true
    · simp only [ht, ↓reduceIte, closeVert_eq, WalkM.run_bind]
      have hdP := hdest ht
      have hce := hok.ears ht
      have st₁ : Step D curV _ (feS₁ d o s) := Step.closeEars st₀.inv st₀.shape hv₀ hce
      have hv₁ : curV < (feS₁ d o s).g.nv := by rw [st₁.g]; exact hv₀
      have st₂ : Step D curV _ (feS₂ d o s) := Step.mergeLate st₁.inv st₁.shape hv₁ (hok.late ht)
      have hv₂ : curV < (feS₂ d o s).g.nv := by rw [st₂.g]; exact hv₁
      have st₀₂ : Step D curV s (feS₂ d o s) := (st₀.trans st₁).trans st₂
      have hg₂ : (feS₂ d o s).g = s.g := st₀₂.g
      obtain ⟨c, mid, py, vy, hC⟩ := hE.close ht hlow
      have hts₂ : (feS₂ d o s).tstack = (c :: mid ++ [py, vy]) ++ base := by
        rw [hC.tstack]
      have hne₂ : (feS₂ d o s).tstack ≠ [] := by rw [hts₂]; simp
      have hD₂ : Dec (feS₂ d o s) (c :: mid ++ [py, vy]) [] base curV (feS₂ d o s) :=
        Dec.init hts₂ (by simp)
      have st₁' : Step D curV (feS₀ d o s) (ceS₁ o.dest d o.e (feS₀ d o s)) :=
        Step.pushEdge st₀.inv st₀.shape o.dest d o.e hce.e_lt hce.q hce.ends hce.d_le
      have hv₁' : curV < (ceS₁ o.dest d o.e (feS₀ d o s)).g.nv := by rw [st₁'.g]; exact hv₀
      have hg₁ : (ceS₁ o.dest d o.e (feS₀ d o s)).g = s.g := by rw [st₁'.g, st₀.g]
      have hit₁ : (ceS₁ o.dest d o.e (feS₀ d o s)).items = (feS₀ d o s).items := rfl
      have hts₁ : (ceS₁ o.dest d o.e (feS₀ d o s)).tstack =
          ⟨o.dest, d, (feS₀ d o s).nxtEdgeIdx,
            setSides (feS₀ d o s).stackDir[d]! [edgeItem (feS₀ d o s).g o.e] []⟩ :: s.tstack := rfl
      have fr₁ : OFrame (ceS₁ o.dest d o.e (feS₀ d o s)) curV (feS₁ d o s) := by
        show OFrame _ curV (after (WalkM.loop _ (loop1Cond d) (Spqr.loop1Body d s.stackDir[d]!))
          (ceS₁ o.dest d o.e (feS₀ d o s)))
        refine loop_pred (OFrame (ceS₁ o.dest d o.e (feS₀ d o s)) curV) _ _ _ (fun _ => rfl)
          (fun k hk h => ?_) OFrame.refl
        have st : Step D curV _
            (iter (Spqr.loop1Body d s.stackDir[d]!) k (ceS₁ o.dest d o.e (feS₀ d o s))) :=
          Step.iter' (loop1Cond d) (Spqr.loop1Body d s.stackDir[d]!) (Loop1BodyOk D d s.stackDir[d]!)
            (fun _ hv hi hs _ hok => Step.loop1Body hi hs hv hok) hv₁' st₁'.inv st₁'.shape hce.body k
            (fun j hj => hk j (Nat.le_of_lt hj))
        have htwo : 2 ≤ (iter (Spqr.loop1Body d s.stackDir[d]!) k
            (ceS₁ o.dest d o.e (feS₀ d o s))).tstack.length := by
          have := hk k (Nat.le_refl k)
          rw [run_loop1Cond] at this
          simp only [Bool.and_eq_true, decide_eq_true_eq] at this
          exact this.1
        exact h.loop1Body st.inv st.shape (hce.body k hk) htwo
      have fr₂ : OFrame (ceS₁ o.dest d o.e (feS₀ d o s)) curV (feS₂ d o s) := fr₁.mergeLate d hne₂
      have hbase₂ : ∀ t ∈ base, ∀ e, t.edges s.g (feS₂ d o s).items e ↔ t.edges s.g s.items e :=
        hC.base_edges
      have hsubE : ∀ e, e < s.g.ne → (∃ t ∈ c :: mid ++ [py, vy], t.edges s.g (feS₂ d o s).items e) →
          (∃ u ∈ sub, u.edges s.g s.items e) ∨ e = o.e := fun e he ⟨t, ht, hte⟩ => by
        by_cases heo : e = o.e
        · exact Or.inr heo
        · exact Or.inl (hE.sub_cover e he heo (hC.sub_edges t ht e he hte))
      have key : ∀ s' : WalkState, OFrame (ceS₁ o.dest d o.e (feS₀ d o s)) curV s' →
          (∃ Bh Bt, ([] : List TEntry) ++ base = Bh ++ Bt ∧
            (∀ b ∈ Bh, b.vStart = curV ∨ b ∈ ([] : List TEntry)) ∧
            Dec (feS₂ d o s) (c :: mid ++ [py, vy]) Bh Bt curV s') →
          OwnedD σ sts origs P d (n + 1) s' := by
        rintro s' fr ⟨Bh, Bt, hBB, hBh, N, hts', -, hcov', hsrc', hBt'⟩
        rw [List.nil_append] at hBB
        have hg' : s'.g = s.g := fr.g.trans hg₁
        have hmemB : ∀ b ∈ Bh, b ∈ base := fun b hb => by
          rw [hBB]; exact List.mem_append_left _ hb
        have hmemT : ∀ b ∈ Bt, b ∈ base := fun b hb => by
          rw [hBB]; exact List.mem_append_right _ hb
        refine OwnedD.of_split (sub := sub) (Bh := Bh) (Bt := Bt) (N := N) (o := o) ho hg'
          (fr.sv.trans st₁'.sv |>.trans st₀.sv) ?_ ?_ ?_ hsvlt hcur hsts hn horigs ?_ ?_ hts'
          ?_ ?_ ?_ ?_ hnd hσ hpos
        · intro p hp
          rw [fr.ch p (by rw [hg₁]; exact hp), hit₁]
          exact hch₀ p
        · intro a ha i
          rw [fr.below a (by rw [hg₁]; exact ha) i, hit₁]
          exact hbl₀ a i
        · intro t ht
          rcases fr.vstart t ht with hc | ⟨u, hu, hut⟩
          · exact Or.inr (Or.inl hc)
          · rw [hts₁] at hu
            rcases List.mem_cons.1 hu with rfl | hu
            · exact Or.inl (by rw [← hut]; exact hdP)
            · exact Or.inr (Or.inr ⟨u, hu, hut⟩)
        · rw [← List.length_append, ← hBB]; exact horig
        · rw [hts, hBB]
        · exact fun b hb => (hBh b hb).resolve_right (by simp)
        · intro t ht e he
          have h1 := hBt' t ht e (by rw [hg₂]; exact he)
          rw [hg₂] at h1
          rw [hg']
          exact h1.trans (hbase₂ t (hmemT t ht) e)
        · intro t ht e he hte
          rw [hg'] at hte
          have hx := hsrc' t ht e (by rw [hg₂]; exact he) (by rw [hg₂]; exact hte)
          unfold XEdges at hx
          rcases hx with hx | ⟨b, hb, hbe⟩ | hx
          · rw [hg₂] at hx
            rcases hsubE e he hx with hx | hx
            · exact Or.inl hx
            · exact Or.inr (Or.inl hx)
          · rw [hg₂] at hbe
            exact Or.inr (Or.inr (Or.inl ⟨b, hb, (hbase₂ b (hmemB b hb) e).1 hbe⟩))
          · rw [hg₂] at hx
            exact Or.inr (Or.inr (Or.inr ((st₀₂.below _).1 hx)))
        · intro e he hx
          rw [hg']
          have hcv : (∃ t ∈ c :: mid ++ [py, vy], t.edges (feS₂ d o s).g (feS₂ d o s).items e) ∨
              ∃ b ∈ Bh, b.edges (feS₂ d o s).g (feS₂ d o s).items e := by
            rw [hg₂]
            rcases hx with ⟨u, hu, hue⟩ | rfl | ⟨b, hb, hbe⟩
            · exact Or.inl (hC.sub_cover e he (hE.sub_edges u hu e he hue))
            · exact Or.inl ⟨c, by simp, hC.c_edge⟩
            · exact Or.inr ⟨b, hb, (hbase₂ b (hmemB b hb) e).2 hbe⟩
          have := hcov' e (by rw [hg₂]; exact he) hcv
          rw [hg₂] at this
          exact this
      cases hasVert
      · simp only [Bool.false_eq_true, ↓reduceIte]
        exact key _ (fr₂.finishRest st₂.inv st₂.shape (hok.rest_tree ht rfl))
          (hD₂.finishRest rfl st₂.inv st₂.shape hv₂ (fun _ => Iff.rfl) (hok.rest_tree ht rfl))
      · simp only [↓reduceIte, WalkM.run_bind]
        have hcv := hok.vert ht rfl
        have hlen3 : base.length + 3 ≤ (feS₂ d o s).tstack.length := by
          rw [hC.tstack]
          simp only [List.length_cons, List.length_append, List.length_nil]; omega
        have h1 : o.cls.isType1 = true → (feS₂ d o s).tstack.length = base.length + 3 := fun h1 => by
          rw [hC.tstack, (hC.type1 h1).1]
          simp only [List.length_cons, List.length_append, List.length_nil]; omega
        have hD₃ := hD₂.closeVert' rfl st₂.inv st₂.shape hv₂ hlen3 h1 hcv
        have st₃ : Step D curV _ (feS₃ curV d o base.length s) :=
          Step.closeVert' st₂.inv st₂.shape hv₂ hcv
        have hv₃ : curV < (feS₃ curV d o base.length s).g.nv := by rw [st₃.g]; exact hv₂
        have fr₃ : OFrame (ceS₁ o.dest d o.e (feS₀ d o s)) curV (feS₃ curV d o base.length s) :=
          fr₂.closeVert' st₂.inv st₂.shape hv₂ hcv
        exact key _ (fr₃.finishRest st₃.inv st₃.shape (hok.rest_vert ht rfl))
          (hD₃.finishRest st₃.g st₃.inv st₃.shape hv₃ st₃.below (hok.rest_vert ht rfl))
    · have ht' : o.cls.isTree = false := Bool.eq_false_iff.2 ht
      simp only [ht', Bool.false_eq_true, ↓reduceIte, WalkM.run_bind]
      have hq : Items.ch (feS₀ d o s).items (edgeItem (feS₀ d o s).g o.e) = [] := by
        rw [hch₀]; exact hok.q ht'
      have st₁ : Step D curV _ (after (pushEdgeTstack curV lv o.e) (feS₀ d o s)) :=
        Step.pushEdge st₀.inv st₀.shape curV lv o.e hok.e_lt hq (hok.ends ht') (hok.lv_le ht')
      have st₂ : Step D curV _ (feBack curV lv d o s) := Step.frame st₁.inv st₁.shape rfl rfl rfl rfl
      have st₀₂ : Step D curV s (feBack curV lv d o s) := (st₀.trans st₁).trans st₂
      have hv₂ : curV < (feBack curV lv d o s).g.nv := by rw [st₀₂.g]; exact hv
      have hgB : (feBack curV lv d o s).g = s.g := st₀₂.g
      have hsub := hE.back_nil ht'
      subst hsub
      show OwnedD σ sts origs P d (n + 1)
        (after (finishRest curV d lv o.cls.isType1 hasVert true) (feBack curV lv d o s))
      have htsB : (feBack curV lv d o s).tstack =
          [⟨curV, lv, (feS₀ d o s).nxtEdgeIdx,
            setSides (feS₀ d o s).stackDir[lv]! [edgeItem (feS₀ d o s).g o.e] []⟩] ++ base := by
        show _ :: s.tstack = _
        rw [hts]; rfl
      have fr : OFrame s curV (feBack curV lv d o s) :=
        ((OFrame.refl.modifyVs (edgeItem s.g o.e)
          (setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest))).pushEdge lv o.e).same
          rfl rfl rfl rfl
      have hqeE : ∀ e, TEntry.edges s.g (feBack curV lv d o s).items
          ⟨curV, lv, (feS₀ d o s).nxtEdgeIdx,
            setSides (feS₀ d o s).stackDir[lv]! [edgeItem (feS₀ d o s).g o.e] []⟩ e ↔ e = o.e :=
        fun e => TEntry.edges_edgeEntry (g := s.g) (items := (feBack curV lv d o s).items)
          (feS₀ d o s).stackDir[lv]! curV lv (feS₀ d o s).nxtEdgeIdx o.e hq e
      have hbaseB : ∀ t ∈ base, ∀ e, t.edges s.g (feBack curV lv d o s).items e ↔
          t.edges s.g s.items e := fun t ht e =>
        TEntry.edges_modify_of_not_mem (edgeItem s.g o.e)
          (fun it => { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) })
          hE.q_root (hE.q_free t (by rw [hts]; exact ht)) e
      obtain ⟨Bh, Bt, hBB, hBh, N, hts', -, hcov', hsrc', hBt'⟩ :=
        (Dec.init (s₂ := feBack curV lv d o s) (curV := curV) htsB (by simp)).finishRest rfl
          st₀₂.inv st₀₂.shape hv₂ (fun _ => Iff.rfl) (hok.rest_back ht')
      rw [List.nil_append] at hBB
      have frR := fr.finishRest st₀₂.inv st₀₂.shape (hok.rest_back ht')
      have hmemB : ∀ b ∈ Bh, b ∈ base := fun b hb => by
        rw [hBB]; exact List.mem_append_left _ hb
      have hmemT : ∀ b ∈ Bt, b ∈ base := fun b hb => by
        rw [hBB]; exact List.mem_append_right _ hb
      refine OwnedD.of_split (sub := []) (Bh := Bh) (Bt := Bt) (N := N) (o := o) ho frR.g frR.sv
        frR.ch frR.below ?_ hsvlt hcur hsts hn horigs ?_ ?_ hts' ?_ ?_ ?_ ?_ hnd hσ hpos
      · intro t ht
        rcases frR.vstart t ht with hc | ⟨u, hu, hut⟩
        · exact Or.inr (Or.inl hc)
        · exact Or.inr (Or.inr ⟨u, hu, hut⟩)
      · rw [← List.length_append, ← hBB]; exact horig
      · rw [hts, hBB]
      · exact fun b hb => (hBh b hb).resolve_right (by simp)
      · intro t ht e he
        have h1 := hBt' t ht e (by rw [hgB]; exact he)
        rw [hgB] at h1
        rw [frR.g]
        exact h1.trans (hbaseB t (hmemT t ht) e)
      · intro t ht e he hte
        rw [frR.g] at hte
        have hx := hsrc' t ht e (by rw [hgB]; exact he) (by rw [hgB]; exact hte)
        unfold XEdges at hx
        rcases hx with ⟨u, hu, hue⟩ | ⟨b, hb, hbe⟩ | hx
        · rw [List.mem_singleton] at hu
          subst hu
          rw [hgB] at hue
          exact Or.inr (Or.inl ((hqeE e).1 hue))
        · rw [hgB] at hbe
          exact Or.inr (Or.inr (Or.inl ⟨b, hb, (hbaseB b (hmemB b hb) e).1 hbe⟩))
        · rw [hgB] at hx
          exact Or.inr (Or.inr (Or.inr ((st₀₂.below _).1 hx)))
      · intro e he hx
        rw [frR.g]
        have hcv : (∃ t ∈ [(⟨curV, lv, (feS₀ d o s).nxtEdgeIdx,
              setSides (feS₀ d o s).stackDir[lv]! [edgeItem (feS₀ d o s).g o.e] []⟩ : TEntry)],
              t.edges (feBack curV lv d o s).g (feBack curV lv d o s).items e) ∨
            ∃ b ∈ Bh, b.edges (feBack curV lv d o s).g (feBack curV lv d o s).items e := by
          rw [hgB]
          rcases hx with ⟨u, hu, _⟩ | rfl | ⟨b, hb, hbe⟩
          · simp at hu
          · exact Or.inl ⟨_, List.mem_singleton_self _, (hqeE _).2 rfl⟩
          · exact Or.inr ⟨b, hb, (hbaseB b (hmemB b hb) e).2 hbe⟩
        have := hcov' e (by rw [hgB]; exact he) hcv
        rw [hgB] at this
        exact this

end WalkState
end Spqr
