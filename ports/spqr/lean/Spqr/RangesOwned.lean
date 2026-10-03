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

end WalkState
end Spqr
