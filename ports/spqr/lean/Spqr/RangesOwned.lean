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

theorem after_mergeTstackTops_eq {a b : TEntry} {rest : List TEntry} (hts : s.tstack = a :: b :: rest) :
    after mergeTstackTops s = { s with tstack := TEntry.mergeInto a b :: rest } := by
  show (mergeTstackTops.run s).2 = _; rw [mergeTstackTops_run_eq s a b rest hts]

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

/-- `OwnedD` is preserved by `finishEdge` (the processed prefix grows by `o.e`), given monotone
starts/stack lengths along the path and that the path below `d` avoids `curV`.
(Checker: `own_*` at `post`.) -/
theorem finishEdge_ownedD {curV d : Nat} {o : DfsOut} {origTstack : Nat} {hasVert : Bool}
    (h : OwnSite σ n D curV d o origTstack hasVert s)
    {sts origs : List Nat} {P : ItemId → Prop} (ho : OwnedD σ sts origs P d n s)
    (hsts : ∀ k, k ≤ d → sts[k]! ≤ sts[d]!) (hn : sts[d]! ≤ n)
    (horigs : ∀ k, k ≤ d → origs[k]! ≤ origs[d]!) (horig : origs[d]! ≤ origTstack)
    (hpath : ∀ k, k < d → s.stackVerts[k]! ≠ curV) :
    OwnedD σ sts origs P d (n + 1) (after (finishEdge curV d o origTstack hasVert) s) := by
  sorry

/-- Every edge under `vertItem v` is held by an open entry. -/
def VertCover (v : Nat) (s : WalkState) : Prop :=
  ∀ e, e < s.g.ne → Items.EdgeBelow s.g s.items (vertItem v) e → ∃ t ∈ s.tstack, t.edges s.g s.items e

/-- Once the vertex entry of `curV` is on the stack (`hasVert`), `finishEdge` keeps every edge under
`vertItem curV` on the stack (`EarFinish.vert` at the pre-state: `vertItem curV` is spanned by a `base`
entry; the merges and closes keep it spanned or put it under the pushed node). (Checker: `own_vcover`
at `post`/`end`.) -/
theorem finishEdge_vertCover {curV d : Nat} {o : DfsOut} {origTstack : Nat} {hasVert : Bool}
    (h : OwnSite σ n D curV d o origTstack hasVert s) :
    wp (finishEdge curV d o origTstack hasVert) (fun hv' s' => hv' = true → VertCover curV s') s := by
  sorry

end WalkState
end Spqr
