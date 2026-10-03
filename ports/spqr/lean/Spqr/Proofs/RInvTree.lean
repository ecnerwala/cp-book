import Spqr.Proofs.RInvBack

/-!
# `RInvTop` across the tree-edge branch of `finishEdge`

The tree-edge branch (`finishTree`) rewrites the stack above the split `origTstack` through
loop 1 (positionally: entries below the split are untouched), then loop 2, `closeVert` and
`finishRest` merge the frontier into entries starting at `curV` or topping out strictly above
`d` — all exempt from `RInvTop.entries`. `RInvG` is the one invariant covering both phases: the
settled bottom `n` of the stack (positional, as in `RInvFront`) with one hole `c`, the single
entry whose settledness is the R content proper (`finishEdge_tree_cur_entryR`); `RInvH` is
`RInvG` over the whole stack.
-/

namespace Spqr
open WalkM
namespace WalkState
variable {s : WalkState} {dfs : DfsData}

/-- An entry exempt from `RInvTop.entries`: tops out strictly above `d` or starts at `v`. -/
def Exempt (v d : Nat) (t : TEntry) : Prop := d ≤ t.topDepth → t.vStart = v

theorem Exempt.mergeInto {v d : Nat} (cur : TEntry) {nxt : TEntry} (h : Exempt v d nxt) :
    Exempt v d (TEntry.mergeInto cur nxt) := fun hd =>
  h (le_trans hd (Nat.min_le_left _ _))

theorem Exempt.of_vStart {v d : Nat} {t : TEntry} (h : t.vStart = v) : Exempt v d t := fun _ => h

theorem Exempt.of_lt {v d : Nat} {t : TEntry} (h : t.topDepth < d) : Exempt v d t := fun hd =>
  absurd hd (Nat.not_le.2 h)

/-- `EntryR` for the non-exempt entries of the bottom `n` of the stack other than the hole `c`
(all entries when `s.tstack.length ≤ n`), and whole-stack edge-disjointness. -/
structure RInvG (s : WalkState) (dfs : DfsData) (v d n : Nat) : Prop where
  entries : ∀ t ∈ s.tstack.drop (s.tstack.length - n), d ≤ t.topDepth → t.vStart ≠ v →
    s.EntryR dfs t
  disj : s.tstack.Pairwise fun t t' => ∀ e, t.edges s.g s.items e → ¬t'.edges s.g s.items e

/-- Everything below the top entry is settled: the top is the entry under construction. -/
def RInvH (s : WalkState) (dfs : DfsData) (v d : Nat) : Prop :=
  s.RInvG dfs v d (s.tstack.length - 1)

/-- The top entry is settled. -/
def SettledTop (s : WalkState) (dfs : DfsData) (v d : Nat) : Prop :=
  s.tstack ≠ [] → d ≤ (curE s).topDepth → (curE s).vStart ≠ v → s.EntryR dfs (curE s)

theorem RInvG.of_front {v d n : Nat} (h : s.RInvFront dfs v d n) : s.RInvG dfs v d n :=
  ⟨fun t ht hd hv => h.base t ht hd hv, h.disj⟩

theorem RInvG.mono {v d n n' : Nat} (hn : n' ≤ n) (h : s.RInvG dfs v d n) : s.RInvG dfs v d n' :=
  ⟨fun t ht => h.entries t (by
    have : s.tstack.drop (s.tstack.length - n') =
        (s.tstack.drop (s.tstack.length - n)).drop
          ((s.tstack.length - n') - (s.tstack.length - n)) := by
      rw [List.drop_drop]; congr 1; omega
    rw [this] at ht; exact List.mem_of_mem_drop ht), h.disj⟩

theorem RInvTop.toG {v d : Nat} (n : Nat) (h : s.RInvTop dfs v d) : s.RInvG dfs v d n :=
  ⟨fun t ht => h.entries t (List.mem_of_mem_drop ht), h.disj⟩

theorem RInvTop.toH {v d : Nat} (h : s.RInvTop dfs v d) : s.RInvH dfs v d := h.toG _

theorem mem_tail_of_mem_drop_pred {α : Type*} {t : α} {l : List α}
    (h : t ∈ l.drop (l.length - (l.length - 1))) : t ∈ l.tail := by
  match l, h with
  | [], h => simp at h
  | _ :: l, h => simpa using h

theorem mem_drop_pred_of_mem_tail {α : Type*} {t : α} {l : List α} (h : t ∈ l.tail) :
    t ∈ l.drop (l.length - (l.length - 1)) := by
  match l, h with
  | [], h => simp at h
  | _ :: l, h => simpa using h

theorem RInvH.entries {v d : Nat} (h : s.RInvH dfs v d) :
    ∀ t ∈ s.tstack.tail, d ≤ t.topDepth → t.vStart ≠ v → s.EntryR dfs t :=
  fun t ht => RInvG.entries h t (mem_drop_pred_of_mem_tail ht)

theorem RInvH.disj {v d : Nat} (h : s.RInvH dfs v d) :
    s.tstack.Pairwise fun t t' => ∀ e, t.edges s.g s.items e → ¬t'.edges s.g s.items e :=
  RInvG.disj h

theorem RInvH.of_tail {v d : Nat}
    (hent : ∀ t ∈ s.tstack.tail, d ≤ t.topDepth → t.vStart ≠ v → s.EntryR dfs t)
    (hdisj : s.tstack.Pairwise fun t t' => ∀ e, t.edges s.g s.items e → ¬t'.edges s.g s.items e) :
    s.RInvH dfs v d :=
  ⟨fun t ht => hent t (mem_tail_of_mem_drop_pred ht), hdisj⟩

theorem RInvH.toTop {v d : Nat} (hc : s.SettledTop dfs v d) (h : s.RInvH dfs v d) :
    s.RInvTop dfs v d := by
  refine ⟨fun t ht hd hne => ?_, h.disj⟩
  match hts : s.tstack with
  | [] => rw [hts] at ht; simp at ht
  | a :: l =>
    rw [hts] at ht
    rcases List.mem_cons.1 ht with rfl | ht
    · have hc' := hc (by rw [hts]; simp)
      rw [curE, hts] at hc'
      exact hc' hd hne
    · exact h.entries t (by rw [hts]; exact ht) hd hne

theorem mem_drop_cons {α : Type*} {x t : α} {l : List α} {k : Nat} (h : t ∈ (x :: l).drop k) :
    (k = 0 ∧ t = x) ∨ t ∈ l.drop (k - 1) := by
  cases k with
  | zero =>
    rw [List.drop_zero] at h
    rcases List.mem_cons.1 h with rfl | h
    · exact .inl ⟨rfl, rfl⟩
    · exact .inr (by rw [Nat.zero_sub, List.drop_zero]; exact h)
  | succ k => exact .inr (by rw [List.drop_succ_cons] at h; rw [Nat.add_sub_cancel]; exact h)

theorem mem_drop_cons_of {α : Type*} {x t : α} {l : List α} {k : Nat} (h : t ∈ l.drop (k - 1)) :
    t ∈ (x :: l).drop k := by
  cases k with
  | zero => rw [Nat.zero_sub, List.drop_zero] at h; exact List.mem_cons_of_mem x h
  | succ k => rw [List.drop_succ_cons]; rw [Nat.add_sub_cancel] at h; exact h

/-! ### Frame -/

theorem RInvG.congr {s' : WalkState} {v d n : Nat} (hg : s'.g = s.g)
    (hsv : s'.stackVerts = s.stackVerts) (hts : s'.tstack = s.tstack)
    (hty : ∀ t ∈ s.tstack, ∀ i ∈ t.spans.1 ++ t.spans.2, Items.type s'.items i = Items.type s.items i)
    (hvs : ∀ t ∈ s.tstack, ∀ i ∈ t.spans.1 ++ t.spans.2, Items.vs s'.items i = Items.vs s.items i)
    (hE : ∀ t ∈ s.tstack, ∀ i ∈ t.spans.1 ++ t.spans.2, ∀ e,
      Items.EdgeBelow s.g s'.items i e ↔ Items.EdgeBelow s.g s.items i e)
    (h : s.RInvG dfs v d n) : s'.RInvG dfs v d n := by
  refine ⟨fun t ht hd hv => ?_, ?_⟩
  · rw [hts] at ht
    exact EntryR.congr (s := s) hg hsv (hty t (List.mem_of_mem_drop ht))
      (hvs t (List.mem_of_mem_drop ht)) (hE t (List.mem_of_mem_drop ht)) (h.entries t ht hd hv)
  · rw [hts, hg]
    exact List.Pairwise.imp_of_mem (fun {a b} ha hb hab e he hbe =>
      hab e ((TEntry.edges_congr (hE a ha) e).1 he) ((TEntry.edges_congr (hE b hb) e).1 hbe)) h.disj

theorem RInvG.of_eq {s' : WalkState} {v d n : Nat} (hg : s'.g = s.g)
    (hi : s'.items = s.items) (hsv : s'.stackVerts = s.stackVerts) (hts : s'.tstack = s.tstack)
    (h : s.RInvG dfs v d n) : s'.RInvG dfs v d n :=
  h.congr hg hsv hts (fun _ _ _ _ => by rw [hi]) (fun _ _ _ _ => by rw [hi]) (fun _ _ _ _ _ => by rw [hi])

theorem RInvG.setStackDir {v d n k : Nat} {b : Bool} (h : s.RInvG dfs v d n) :
    (after (setStackDir k b) s).RInvG dfs v d n :=
  RInvG.of_eq (s := s) (s' := after (WalkM.setStackDir k b) s) rfl rfl rfl rfl h

theorem RInvG.modify {v d n : Nat} (f : WalkState → WalkState) (hg : (f s).g = s.g)
    (hi : (f s).items = s.items) (hsv : (f s).stackVerts = s.stackVerts)
    (hts : (f s).tstack = s.tstack) (h : s.RInvG dfs v d n) :
    (after (modify f : WalkM Unit) s).RInvG dfs v d n :=
  RInvG.of_eq (s := s) (s' := after (_root_.modify f : WalkM Unit) s) hg hi hsv hts h

theorem RInvG.modifyVs_free {v d n : Nat} (j : ItemId) (f : Item → Item)
    (hch : ∀ it, (f it).ch = it.ch) (hty : ∀ it, (f it).type = it.type)
    (hfree : ∀ t ∈ s.tstack, j ∉ t.spans.1 ++ t.spans.2)
    (h : s.RInvG dfs v d n) : (after (modifyItem j f) s).RInvG dfs v d n := by
  refine RInvG.congr (s := s) (s' := after (modifyItem j f) s) rfl rfl rfl
    (fun t ht i hi => ?_) (fun t ht i hi => ?_) (fun t ht i hi e => ?_) h
  · exact Items.type_modify_type_eq j f hty i
  · exact Items.vs_modify_of_ne j f fun hij => hfree t ht (hij ▸ hi)
  · exact Items.Below_modify_ch_eq j f hch

theorem RInvG.alloc {v d n : Nat} (ty : NodeType) (hs : Shape s) (h : s.RInvG dfs v d n) :
    ((allocItem ty).run s).2.RInvG dfs v d n := by
  rw [run_allocItem]
  refine RInvG.congr (s := s) rfl rfl rfl (fun t ht i hi => ?_) (fun t ht i hi => ?_)
    (fun t ht i hi e => ?_) h
  · exact Items.type_push_of_ne _ (Nat.ne_of_lt (hs.span t ht i hi))
  · exact Items.vs_push_of_ne _ (Nat.ne_of_lt (hs.span t ht i hi))
  · exact Items.Below_push_nil _ rfl

/-! ### Stack primitives -/

/-- Pushing the single edge entry of `e` (a childless Q item) above the settled bottom `n`. -/
theorem RInvG.pushEdge {v d n : Nat} (u k e : Nat) (hn : n ≤ s.tstack.length)
    (hq : Items.ch s.items (edgeItem s.g e) = [])
    (hown : ∀ t ∈ s.tstack, ¬ t.edges s.g s.items e) (h : s.RInvG dfs v d n) :
    (after (pushEdgeTstack u k e) s).RInvG dfs v d n := by
  show RInvG { s with tstack := ⟨u, k, s.nxtEdgeIdx, setSides s.stackDir[k]! [edgeItem s.g e] []⟩ :: s.tstack } dfs v d n
  have hE := TEntry.edges_edgeEntry (g := s.g) (items := s.items) s.stackDir[k]! u k s.nxtEdgeIdx e hq
  refine ⟨fun t ht hd hne => ?_, List.pairwise_cons.2 ⟨fun t' ht' e' he' => ?_, h.disj⟩⟩
  · simp only [List.length_cons] at ht
    rw [Nat.succ_sub hn, List.drop_succ_cons] at ht
    exact EntryR.congr (s := s) rfl rfl (fun _ _ => rfl) (fun _ _ => rfl) (fun _ _ _ => Iff.rfl)
      (h.entries t ht hd hne)
  · rw [(hE e').1 he']; exact hown t' ht'

/-- Pushing the vertex entry of `v` (exempt). -/
theorem RInvG.pushVert {v d n k : Nat}
    (hown : ∀ t ∈ s.tstack, ∀ e, t.edges s.g s.items e → ¬ Items.EdgeBelow s.g s.items (vertItem v) e)
    (h : s.RInvG dfs v d n) : (after (pushVertTstack v k) s).RInvG dfs v d n := by
  show RInvG { s with tstack := ⟨v, k, s.nxtEdgeIdx, setSides s.stackDir[k]! [vertItem v] []⟩ :: s.tstack } dfs v d n
  have hed : ∀ e, (⟨v, k, s.nxtEdgeIdx, setSides s.stackDir[k]! [vertItem v] []⟩ : TEntry).edges
      s.g s.items e → Items.EdgeBelow s.g s.items (vertItem v) e := by
    rintro e ⟨i, hi, hb⟩
    rw [List.mem_singleton.1 ((mem_setSides _ _ i).1 hi)] at hb
    exact hb
  refine ⟨fun t ht hd hne => ?_, List.pairwise_cons.2 ⟨fun t' ht' e he hte => hown t' ht' e hte (hed e he), h.disj⟩⟩
  rcases mem_drop_cons ht with ⟨-, rfl⟩ | ht
  · exact absurd rfl hne
  · simp only [List.length_cons] at ht
    exact EntryR.congr (s := s) rfl rfl (fun _ _ => rfl) (fun _ _ => rfl) (fun _ _ _ => Iff.rfl)
      (h.entries t (by
        rcases Nat.lt_or_ge n (s.tstack.length + 1) with hlt | hge
        · have : s.tstack.length + 1 - n - 1 = s.tstack.length - n := by omega
          rw [this] at ht; exact ht
        · rw [Nat.sub_eq_zero_of_le (by omega), List.drop_zero]; exact List.mem_of_mem_drop ht) hd hne)

/-- Merging the top two entries; the merged entry is exempt whenever it lands in the bottom `n`. -/
theorem RInvG.mergeTop {v d n : Nat} (cur nxt : TEntry) (rest : List TEntry)
    (hts : s.tstack = cur :: nxt :: rest)
    (hx : s.tstack.length ≤ n + 1 → Exempt v d (TEntry.mergeInto cur nxt))
    (h : s.RInvG dfs v d n) : (after mergeTstackTops s).RInvG dfs v d n := by
  show (mergeTstackTops.run s).2.RInvG dfs v d n
  rw [mergeTstackTops_run_eq s cur nxt rest hts]
  have hd := h.disj; rw [hts] at hd
  obtain ⟨hcur, hd⟩ := List.pairwise_cons.1 hd
  obtain ⟨hnxt', hrest⟩ := List.pairwise_cons.1 hd
  have hlen : s.tstack.length = rest.length + 2 := by rw [hts]; rfl
  refine ⟨fun t ht hdep hne => ?_, List.pairwise_cons.2 ⟨fun t' ht' e he => ?_, hrest⟩⟩
  · simp only [List.length_cons] at ht
    rcases mem_drop_cons ht with ⟨hk, rfl⟩ | ht
    · exact absurd (hx (by omega) hdep) hne
    · refine EntryR.congr (s := s) rfl rfl (fun _ _ => rfl) (fun _ _ => rfl) (fun _ _ _ => Iff.rfl)
        (h.entries t ?_ hdep hne)
      rw [hts]; simp only [List.length_cons]
      apply mem_drop_cons_of; apply mem_drop_cons_of
      have : rest.length + 1 + 1 - n - 1 - 1 = rest.length + 1 - n - 1 := by omega
      rw [this]; exact ht
  · rcases (TEntry.edges_mergeInto cur nxt e).1 he with he | he
    · exact hcur t' (by simp [ht']) e he
    · exact hnxt' t' ht' e he

/-- Finishing the top entry into a fresh root item; the top is exempt whenever it lies in the
bottom `n`. -/
theorem RInvG.finishTop {v d n : Nat} (item : ItemId) (t : TEntry) (rest : List TEntry)
    (hts : s.tstack = t :: rest) (hx : s.tstack.length ≤ n → Exempt v d t) (hitem : item < s.items.size)
    (hroot : ∀ p, ¬ Items.IsParent s.items p item)
    (hfree : ∀ t' ∈ s.tstack, item ∉ t'.spans.1 ++ t'.spans.2)
    (hside : getSide t.spans (!s.stackDir[t.topDepth]!) = [])
    (h : s.RInvG dfs v d n) : (after (finishTstackTop item) s).RInvG dfs v d n := by
  show ((finishTstackTop item).run s).2.RInvG dfs v d n
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
  have hlen : s.tstack.length = rest.length + 1 := by rw [hts]; rfl
  refine ⟨fun u hu hdep hne' => ?_, List.pairwise_cons.2 ⟨fun u hu e he hue => ?_, ?_⟩⟩
  · simp only [List.length_cons] at hu
    rcases mem_drop_cons hu with ⟨hk, rfl⟩ | hu
    · exact absurd (hx (by omega) hdep) hne'
    · have hu' : u ∈ s.tstack.drop (s.tstack.length - n) := by
        rw [hts]; simp only [List.length_cons]; exact mem_drop_cons_of hu
      have hfu := hfree u (by simp [hts]; exact .inr (List.mem_of_mem_drop hu))
      refine EntryR.congr (s := s) rfl rfl (fun i hi => ?_) (fun i hi => ?_) (fun i hi e => ?_)
        (h.entries u hu' hdep hne')
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

/-- `maybeUnwrapNxt ty`; the entry below the top is exempt whenever it lies in the bottom `n`. -/
theorem RInvG.unwrapNxt {v d n : Nat} {ty : NodeType} (hs : Shape s) (hok : UnwrapOk ty s)
    (hx : s.tstack.length ≤ n + 1 → Exempt v d (nxtE s))
    (h : s.RInvG dfs v d n) : (after (maybeUnwrapNxt ty) s).RInvG dfs v d n := by
  have halloc := h.alloc ty hs
  show ((maybeUnwrapNxt ty).run s).2.RInvG dfs v d n
  match hts : s.tstack with
  | [] | [_] => have := hok.two; rw [hts] at this; simp at this
  | a :: b :: rest =>
    have hn : nxtE s = b := by rw [nxtE, hts]; rfl
    have hd : nxtDir s = s.stackDir[b.topDepth]! := by rw [nxtDir, hn]
    have hh : nxtHead s = (getSide b.spans s.stackDir[b.topDepth]!).head! := by rw [nxtHead, hn, hd]
    have hlen : s.tstack.length = rest.length + 2 := by rw [hts]; rfl
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
      · simp only [List.length_cons] at ht
        rcases mem_drop_cons ht with ⟨hk, rfl⟩ | ht
        · exact EntryR.congr (s := s) rfl rfl (fun _ _ => rfl) (fun _ _ => rfl) (fun _ _ _ => Iff.rfl)
            (h.entries t (by rw [hts]; simp only [List.length_cons]; rw [hk]; simp) hdep hne)
        rcases mem_drop_cons ht with ⟨hk, rfl⟩ | ht
        · rw [hn] at hx
          exact absurd (hx (by omega) hdep) hne
        · exact EntryR.congr (s := s) rfl rfl (fun _ _ => rfl) (fun _ _ => rfl) (fun _ _ _ => Iff.rfl)
            (h.entries t (by
              rw [hts]; simp only [List.length_cons]
              exact mem_drop_cons_of (mem_drop_cons_of ht)) hdep hne)
      · rcases List.mem_cons.1 ht with rfl | ht
        · exact ha b (by simp) e he (hsub e hte)
        · exact ha t (by simp [ht]) e he hte
    · simp only [h2, ↓reduceIte]; exact halloc

/-- Re-targeting the top entry to `v` (exempt; its edge set is unchanged). -/
theorem RInvG.retarget {v d n : Nat} (edgeDir : Bool) (t : TEntry) (rest : List TEntry)
    (hts : s.tstack = t :: rest) (h : s.RInvG dfs v d n) :
    (after (WalkState.retarget v edgeDir) s).RInvG dfs v d n := by
  show ((WalkState.retarget v edgeDir).run s).2.RInvG dfs v d n
  rw [retarget_run_eq v edgeDir s t rest hts]
  have hE : ∀ e, TEntry.edges s.g s.items { t with vStart := v, spans := setSides (!edgeDir) (t.spans.1 ++ t.spans.2) [] } e ↔
      t.edges s.g s.items e := by
    intro e; constructor
    · rintro ⟨i, hi, hb⟩; exact ⟨i, (mem_setSides _ _ i).1 hi, hb⟩
    · rintro ⟨i, hi, hb⟩; exact ⟨i, (mem_setSides _ _ i).2 hi, hb⟩
  have hd := h.disj; rw [hts] at hd
  obtain ⟨htop, hrest⟩ := List.pairwise_cons.1 hd
  refine ⟨fun u hu hdep hne => ?_, List.pairwise_cons.2 ⟨fun u hu e he => htop u hu e ((hE e).1 he), hrest⟩⟩
  simp only [List.length_cons] at hu
  rcases mem_drop_cons hu with ⟨-, rfl⟩ | hu
  · exact absurd rfl hne
  · exact EntryR.congr (s := s) rfl rfl (fun _ _ => rfl) (fun _ _ => rfl) (fun _ _ _ => Iff.rfl)
      (h.entries u (by rw [hts]; simp only [List.length_cons]; exact mem_drop_cons_of hu) hdep hne)

/-! ### Loop 1 (positional) -/

/-- One iteration of loop 1 above the settled bottom `n`: the three entries it touches lie in the
frontier (`hlen`, from `Frontier.loop1`), so every entry it creates is in the frontier too. -/
theorem RInvG.loop1Body {D v d d₀ n : Nat} {edgeDir : Bool} (hi : s.Inv' D) (hs : Shape s)
    (hok : Loop1BodyOk D d edgeDir s) (hc : result (loop1Cond d) s = true)
    (hlen : n + (if d < (nxtE s).topDepth then 3 else 2) ≤ s.tstack.length)
    (h : s.RInvG dfs v d₀ n) :
    (after (Spqr.loop1Body d edgeDir) s).RInvG dfs v d₀ n ∧
      (after (Spqr.loop1Body d edgeDir) s).tstack.length + 1 ≤ s.tstack.length := by
  have hc' : 2 ≤ s.tstack.length := by
    have : (decide (s.tstack.length ≥ 2) && decide (s.tstack.tail.head!.topDepth ≥ d)) = true := hc
    simp only [Bool.and_eq_true, decide_eq_true_eq] at this; exact this.1
  have st₁ : Step D v s (l1S₁ d edgeDir s) := Step.loop1Type hi hs hok.mergeS
  have key₁ : (l1S₁ d edgeDir s).RInvG dfs v d₀ n ∧ n + 2 ≤ (l1S₁ d edgeDir s).tstack.length ∧
      (l1S₁ d edgeDir s).tstack.length ≤ s.tstack.length := by
    unfold l1S₁ after; rw [loop1Type_run]
    by_cases hgt : (nxtE s).topDepth > d
    · simp only [hgt, ↓reduceIte]; simp only [hgt, ↓reduceIte] at hlen
      obtain ⟨a, b, rest, hts⟩ : ∃ a b rest, s.tstack = a :: b :: rest := by
        match hts : s.tstack with
        | [] | [_] => rw [hts] at hc'; simp at hc'
        | a :: b :: rest => exact ⟨a, b, rest, rfl⟩
      have h' : ({ s with stackDir := s.stackDir.set! (nxtE s).topDepth edgeDir } : WalkState).RInvG dfs v d₀ n :=
        h.setStackDir
      have hl : s.tstack.length = rest.length + 2 := by rw [hts]; rfl
      refine ⟨h'.mergeTop a b rest hts (fun hl' => absurd hl' (by show ¬ s.tstack.length ≤ n + 1; omega)), ?_⟩
      rw [mergeTstackTops_run_eq { s with stackDir := s.stackDir.set! (nxtE s).topDepth edgeDir } a b rest hts]
      simp only [List.length_cons]; omega
    · simp only [hgt, ↓reduceIte]; simp only [hgt, ↓reduceIte] at hlen
      split <;> exact ⟨h, hlen, le_rfl⟩
  obtain ⟨R₁, hlen₁, hle₁⟩ := key₁
  have r := maybeUnwrapNxt_spec (v := v) st₁.inv st₁.shape (loop1Type_result d edgeDir s) hok.unwrap
  have R₂ := R₁.unwrapNxt st₁.shape hok.unwrap (fun hl => absurd hl (by omega))
  match hts₁ : (l1S₁ d edgeDir s).tstack with
  | [] | [_] => have := hok.unwrap.two; rw [hts₁] at this; simp at this
  | a :: b :: rest =>
    have hl₁ : (l1S₁ d edgeDir s).tstack.length = rest.length + 2 := by rw [hts₁]; rfl
    obtain ⟨b', hts₂, -, -⟩ := maybeUnwrapNxt_tstack (l1Ty d edgeDir s) a b rest hts₁
    set s₂ := after (maybeUnwrapNxt (l1Ty d edgeDir s)) (l1S₁ d edgeDir s) with hs₂
    have hts₂' : s₂.tstack = a :: b' :: rest := hts₂
    have R₃ := R₂.mergeTop a b' rest hts₂' (fun hl => absurd hl (by rw [hts₂']; simp only [List.length_cons]; omega))
    have hts₃ : (mergeTstackTops.run s₂).2.tstack = TEntry.mergeInto a b' :: rest := by
      rw [mergeTstackTops_run_eq s₂ a b' rest hts₂']
    have hf := r.free.merge
    have hside := hok.close.finish.side
    have hcur : curE (after mergeTstackTops s₂) = TEntry.mergeInto a b' := by
      show (mergeTstackTops.run s₂).2.tstack.head! = _
      rw [hts₃]; rfl
    have hcl : l1S₂ d edgeDir s = s₂ := rfl
    rw [hcl] at hside
    rw [hcur] at hside
    have R₄ := R₃.finishTop _ (TEntry.mergeInto a b') rest hts₃
      (fun hl => absurd hl (by rw [hts₃]; simp only [List.length_cons]; omega)) hf.lt hf.root hf.free hside
    refine ⟨R₄, ?_⟩
    show ((finishTstackTop _).run (mergeTstackTops.run s₂).2).2.tstack.length + 1 ≤ s.tstack.length
    rw [finishTstackTop_run_eq _ _ (TEntry.mergeInto a b') rest hts₃]
    simp only [List.length_cons]; omega

/-- Loop 1 above the settled bottom `n`, with its exit: the final state is some iterate at which
the condition fails, and `Step` carries `Inv'`/`Shape`. -/
theorem RInvG.loop1 {D v w d d₀ n : Nat} {edgeDir : Bool} (fuel : Nat)
    (hfuel : s.tstack.length ≤ fuel) (hv : w < s.g.nv) (hi : s.Inv' D) (hs : Shape s)
    (hok : ∀ k, (∀ j, j ≤ k → result (loop1Cond d) (iter (Spqr.loop1Body d edgeDir) j s) = true) →
      Loop1BodyOk D d edgeDir (iter (Spqr.loop1Body d edgeDir) k s))
    (hlen : ∀ k, (∀ j, j < k → result (loop1Cond d) (iter (Spqr.loop1Body d edgeDir) j s) = true) →
      result (loop1Cond d) (iter (Spqr.loop1Body d edgeDir) k s) = true →
      n + (if d < (nxtE (iter (Spqr.loop1Body d edgeDir) k s)).topDepth then 3 else 2) ≤
        (iter (Spqr.loop1Body d edgeDir) k s).tstack.length)
    (h : s.RInvG dfs v d₀ n) :
    ∃ k, after (loop fuel (loop1Cond d) (Spqr.loop1Body d edgeDir)) s = iter (Spqr.loop1Body d edgeDir) k s ∧
      (∀ j, j < k → result (loop1Cond d) (iter (Spqr.loop1Body d edgeDir) j s) = true) ∧
      result (loop1Cond d) (iter (Spqr.loop1Body d edgeDir) k s) = false ∧
      (iter (Spqr.loop1Body d edgeDir) k s).RInvG dfs v d₀ n ∧ Step D w s (iter (Spqr.loop1Body d edgeDir) k s) := by
  induction fuel generalizing s with
  | zero =>
    refine ⟨0, rfl, fun j hj => absurd hj (Nat.not_lt_zero j), ?_, h, Step.refl hi hs⟩
    have h0 : s.tstack.length = 0 := Nat.le_zero.1 hfuel
    show (decide (s.tstack.length ≥ 2) && decide (s.tstack.tail.head!.topDepth ≥ d)) = false
    rw [h0]; rfl
  | succ fuel ih =>
    unfold after
    rw [loop_succ_run fuel _ _ s rfl]
    by_cases hc : ((loop1Cond d).run s).1 = true
    · simp only [hc, ↓reduceIte]
      have h0 : ∀ j, j ≤ 0 → result (loop1Cond d) (iter (Spqr.loop1Body d edgeDir) j s) = true :=
        fun j hj => by rw [Nat.le_zero.1 hj]; exact hc
      obtain ⟨R, hl⟩ := h.loop1Body hi hs (hok 0 h0) hc (hlen 0 (fun j hj => absurd hj (Nat.not_lt_zero j)) hc)
      have st := Step.loop1Body hi hs hv (hok 0 h0)
      obtain ⟨k, heq, hj, hexit, R', st'⟩ := ih (by show (after (Spqr.loop1Body d edgeDir) s).tstack.length ≤ fuel; omega)
        (by rw [st.g]; exact hv) st.inv st.shape
        (fun k hk => hok (k + 1) fun j hj => by
          cases j with
          | zero => exact hc
          | succ j => exact hk j (Nat.le_of_succ_le_succ hj))
        (fun k hk hck => hlen (k + 1) (fun j hj => by
          cases j with
          | zero => exact hc
          | succ j => exact hk j (Nat.lt_of_succ_lt_succ hj)) hck) R
      refine ⟨k + 1, heq, fun j hjk => ?_, hexit, R', st.trans st'⟩
      cases j with
      | zero => exact hc
      | succ j => exact hj j (Nat.lt_of_succ_lt_succ hjk)
    · simp only [hc, Bool.false_eq_true, ↓reduceIte]
      exact ⟨0, rfl, fun j hj => absurd hj (Nat.not_lt_zero j), Bool.eq_false_iff.2 hc, h, Step.refl hi hs⟩

/-- `closeEars` above the settled bottom `n`: the pushed tree edge and loop 1. -/
theorem RInvG.closeEars {D v w d d₀ n : Nat} {nxtV e : Nat} {edgeDir : Bool} (hi : s.Inv' D)
    (hs : Shape s) (hv : w < s.g.nv) (hn : n ≤ s.tstack.length) (hok : CloseEarsOk D nxtV d e edgeDir s)
    (hown : ∀ t ∈ s.tstack, ¬ t.edges s.g s.items e)
    (hlen : ∀ k, (∀ j, j < k → result (loop1Cond d) (iter (Spqr.loop1Body d edgeDir) j (ceS₁ nxtV d e s)) = true) →
      result (loop1Cond d) (iter (Spqr.loop1Body d edgeDir) k (ceS₁ nxtV d e s)) = true →
      n + (if d < (nxtE (iter (Spqr.loop1Body d edgeDir) k (ceS₁ nxtV d e s))).topDepth then 3 else 2) ≤
        (iter (Spqr.loop1Body d edgeDir) k (ceS₁ nxtV d e s)).tstack.length)
    (h : s.RInvG dfs v d₀ n) :
    ∃ k, after (Spqr.closeEars nxtV d e edgeDir) s = iter (Spqr.loop1Body d edgeDir) k (ceS₁ nxtV d e s) ∧
      (∀ j, j < k → result (loop1Cond d) (iter (Spqr.loop1Body d edgeDir) j (ceS₁ nxtV d e s)) = true) ∧
      result (loop1Cond d) (iter (Spqr.loop1Body d edgeDir) k (ceS₁ nxtV d e s)) = false ∧
      (iter (Spqr.loop1Body d edgeDir) k (ceS₁ nxtV d e s)).RInvG dfs v d₀ n ∧
      Step D w s (iter (Spqr.loop1Body d edgeDir) k (ceS₁ nxtV d e s)) := by
  have st₁ : Step D w s (ceS₁ nxtV d e s) := Step.pushEdge hi hs nxtV d e hok.e_lt hok.q hok.ends hok.d_le
  have R₁ : (ceS₁ nxtV d e s).RInvG dfs v d₀ n := h.pushEdge nxtV d e hn hok.q hown
  obtain ⟨k, heq, hj, hexit, R, st⟩ := R₁.loop1 (ceS₁ nxtV d e s).tstack.length le_rfl
    (by rw [st₁.g]; exact hv) st₁.inv st₁.shape hok.body hlen
  exact ⟨k, heq, hj, hexit, R, st₁.trans st⟩

/-- After loop 1: the base is settled positionally and the frontier entries below the top are
settled by hypothesis, so everything below the top is settled. -/
theorem RInvG.widen {v d n : Nat}
    (hset : ∀ t ∈ s.tstack.tail.take (s.tstack.length - 1 - n), d ≤ t.topDepth → t.vStart ≠ v →
      s.EntryR dfs t)
    (h : s.RInvG dfs v d n) : s.RInvH dfs v d := by
  refine RInvH.of_tail (fun t ht hd hne => ?_) h.disj
  rw [← List.take_append_drop (s.tstack.length - 1 - n) s.tstack.tail] at ht
  rcases List.mem_append.1 ht with ht | ht
  · exact hset t ht hd hne
  · refine h.entries t ?_ hd hne
    rw [← List.drop_one, List.drop_drop] at ht
    have : s.tstack.drop (1 + (s.tstack.length - 1 - n)) =
        (s.tstack.drop (s.tstack.length - n)).drop
          ((1 + (s.tstack.length - 1 - n)) - (s.tstack.length - n)) := by
      rw [List.drop_drop]; congr 1; omega
    rw [this] at ht
    exact List.mem_of_mem_drop ht

/-! ### Whole-stack phase (`RInvH`: everything below the top) -/

theorem RInvG.of_nil {v d n : Nat} (hts : s.tstack = []) : s.RInvG dfs v d n :=
  ⟨fun t ht => by rw [hts] at ht; simp at ht, by rw [hts]; exact List.Pairwise.nil⟩

theorem mergeTstackTops_run_short (h : s.tstack.length ≤ 1) : (mergeTstackTops.run s).2.tstack = [] := by
  obtain ⟨g, tern, items, sv, sd, nei, fo, ts, tb, tsl⟩ := s
  match ts, h with
  | [], _ => rfl
  | [_], _ => rfl
  | _ :: _ :: _, h => simp at h

/-- Merging the top into the entry below it: both are above everything settled. -/
theorem RInvH.mergeTop {v d : Nat} (h : s.RInvH dfs v d) : (after mergeTstackTops s).RInvH dfs v d := by
  match hts : s.tstack with
  | [] | [_] =>
    exact RInvG.of_nil (mergeTstackTops_run_short (by rw [hts]; simp))
  | a :: b :: rest =>
    have hlen : s.tstack.length = rest.length + 2 := by rw [hts]; rfl
    have R := (h.mono (n' := rest.length) (by omega)).mergeTop a b rest hts
      (fun h => absurd h (by omega))
    have hl : (after mergeTstackTops s).tstack.length - 1 = rest.length := by
      show ((mergeTstackTops.run s).2).tstack.length - 1 = _
      rw [mergeTstackTops_run_eq s a b rest hts]; simp
    show RInvG _ dfs v d _
    rw [hl]; exact R

theorem RInvH.unwrapNxt {v d : Nat} {ty : NodeType} (hs : Shape s) (hok : UnwrapOk ty s)
    (hx : Exempt v d (nxtE s)) (h : s.RInvH dfs v d) : (after (maybeUnwrapNxt ty) s).RInvH dfs v d := by
  have R := RInvG.unwrapNxt hs hok (fun _ => hx) h
  match hts : s.tstack with
  | [] | [_] => have := hok.two; rw [hts] at this; simp at this
  | a :: b :: rest =>
    obtain ⟨b', hts', -, -⟩ := maybeUnwrapNxt_tstack ty a b rest hts
    have hl : (after (maybeUnwrapNxt ty) s).tstack.length - 1 = s.tstack.length - 1 := by
      rw [hts', hts]; simp
    show RInvG _ dfs v d _
    rw [hl]; exact R

theorem RInvH.finishTop {v d : Nat} (item : ItemId) (t : TEntry) (rest : List TEntry)
    (hts : s.tstack = t :: rest) (hf : ItemFree s item)
    (hside : getSide t.spans (!s.stackDir[t.topDepth]!) = [])
    (h : s.RInvH dfs v d) : (after (finishTstackTop item) s).RInvH dfs v d := by
  have hlen : s.tstack.length = rest.length + 1 := by rw [hts]; rfl
  have R := (h.mono (n' := rest.length) (by omega)).finishTop item t rest hts
    (fun h => absurd h (by omega)) hf.lt hf.root hf.free hside
  have hl : (after (finishTstackTop item) s).tstack.length - 1 = rest.length := by
    show ((finishTstackTop item).run s).2.tstack.length - 1 = _
    rw [finishTstackTop_run_eq s item t rest hts]; simp
  show RInvG _ dfs v d _
  rw [hl]; exact R

theorem RInvH.retarget {v d : Nat} (edgeDir : Bool) (t : TEntry) (rest : List TEntry)
    (hts : s.tstack = t :: rest) (h : s.RInvH dfs v d) :
    (after (WalkState.retarget v edgeDir) s).RInvH dfs v d := by
  have hlen : s.tstack.length = rest.length + 1 := by rw [hts]; rfl
  have R := (h.mono (n' := rest.length) (by omega)).retarget edgeDir t rest hts
  have hl : (after (WalkState.retarget v edgeDir) s).tstack.length - 1 = rest.length := by
    show ((WalkState.retarget v edgeDir).run s).2.tstack.length - 1 = _
    rw [retarget_run_eq v edgeDir s t rest hts]; simp
  show RInvG _ dfs v d _
  rw [hl]; exact R

/-- Pushing the vertex entry of `v` on a settled stack: the old top goes below the top. -/
theorem RInvTop.pushVert_own {v d k : Nat}
    (hown : ∀ t ∈ s.tstack, ∀ e, t.edges s.g s.items e → ¬ Items.EdgeBelow s.g s.items (vertItem v) e)
    (h : s.RInvTop dfs v d) : (after (pushVertTstack v k) s).RInvH dfs v d := by
  have R := (h.toG s.tstack.length).pushVert (k := k) hown
  show RInvG _ dfs v d _
  have hl : (after (pushVertTstack v k) s).tstack.length - 1 = s.tstack.length := by
    show (_ :: s.tstack).length - 1 = _; simp
  rw [hl]; exact R

/-- A merge loop (loops 2 and 3). -/
theorem RInvH.mergeLoop {D v d : Nat} (cond : WalkM Bool) (fuel : Nat)
    (hcond : ∀ s, (cond.run s).2 = s) (hv : v < s.g.nv) (hi : s.Inv' D) (hs : Shape s)
    (hok : ∀ k, (∀ j, j ≤ k → (cond.run (iter mergeTstackTops j s)).1 = true) →
      MergeTopOk D (iter mergeTstackTops k s))
    (h : s.RInvH dfs v d) :
    (after (loop fuel cond mergeTstackTops) s).RInvH dfs v d ∧
      Step D v s (after (loop fuel cond mergeTstackTops) s) := by
  induction fuel generalizing s with
  | zero => exact ⟨h, Step.refl hi hs⟩
  | succ fuel ih =>
    show ((loop (fuel + 1) cond mergeTstackTops).run s).2.RInvH dfs v d ∧
      Step D v s ((loop (fuel + 1) cond mergeTstackTops).run s).2
    rw [loop_succ_run fuel cond mergeTstackTops s (hcond s)]
    by_cases hc : (cond.run s).1 = true
    · simp only [hc, ↓reduceIte]
      have h0 : ∀ j, j ≤ 0 → (cond.run (iter mergeTstackTops j s)).1 = true :=
        fun j hj => by rw [Nat.le_zero.1 hj]; exact hc
      have st := Step.mergeTop (v := v) hi hs (hok 0 h0)
      obtain ⟨R, st'⟩ := ih (by rw [st.g]; exact hv) st.inv st.shape
        (fun k hk => hok (k + 1) fun j hj => by
          cases j with
          | zero => exact hc
          | succ j => exact hk j (Nat.le_of_succ_le_succ hj)) h.mergeTop
      exact ⟨R, st.trans st'⟩
    · simp only [hc, Bool.false_eq_true, ↓reduceIte]; exact ⟨h, Step.refl hi hs⟩

/-- Loop 2. -/
theorem RInvH.mergeLate {D v d : Nat} (hv : v < s.g.nv) (hi : s.Inv' D) (hs : Shape s)
    (hok : MergeLateOk D d s) (h : s.RInvH dfs v d) : (after (Spqr.mergeLate d) s).RInvH dfs v d := by
  show ((Spqr.mergeLate d).run s).2.RInvH dfs v d
  rw [mergeLate_run]
  by_cases hf : (curE s).firstIdx > s.firstOccurrence[d]!
  · simp only [hf, ↓reduceIte]; exact (h.mergeLoop _ _ (fun _ => rfl) hv hi hs (hok.body hf)).1
  · simp only [hf, ↓reduceIte]; exact h

/-- `closeVert`: loop 3, the unwrap (of an exempt entry), the two merges, the re-targeting to `v`
and the type-1 close only touch the top, except the unwrap. -/
theorem RInvH.closeVert' {D v d : Nat} {edgeDir isType1 : Bool} {origTstack : Nat}
    {isSingle : Bool} (hv : v < s.g.nv) (hi : s.Inv' D) (hs : Shape s)
    (hok : CloseVertOk D v edgeDir isType1 origTstack isSingle s)
    (hx1 : isType1 = true → Exempt v d (nxtE s))
    (h : s.RInvH dfs v d) :
    (after (WalkState.closeVert' v edgeDir isType1 origTstack isSingle) s).RInvH dfs v d := by
  have st₁ : Step D v s (cvS₁ isType1 origTstack isSingle s) := Step.vertPre hi hs hv hok.loop3
  have R₁ : (cvS₁ isType1 origTstack isSingle s).RInvH dfs v d := by
    cases isType1
    · show ((vertPre false origTstack isSingle).run s).2.RInvH dfs v d
      simp only [WalkState.vertPre, Bool.not_false, ↓reduceIte, WalkM.run_bind, run_tstackSize, WalkM.pure_run]
      exact (h.mergeLoop _ _ (fun _ => rfl) hv hi hs (hok.loop3 rfl)).1
    · exact h
  obtain ⟨st₂, hfree⟩ := Step.vertUnwrap (v := v) st₁.inv st₁.shape (isType1 := isType1)
    (isSingle := cvB₁ isType1 origTstack isSingle s) (fun h => by subst h; exact hok.unwrap rfl)
  have st₂ : Step D v (cvS₁ isType1 origTstack isSingle s) (cvS₂ isType1 origTstack isSingle s) := st₂
  have R₂ : (cvS₂ isType1 origTstack isSingle s).RInvH dfs v d := by
    cases isType1
    · exact R₁
    · show ((some <$> maybeUnwrapNxt (if isSingle then NodeType.S else .R)).run s).2.RInvH dfs v d
      rw [WalkM.map_run]
      exact R₁.unwrapNxt st₁.shape (hok.unwrap rfl) (hx1 rfl)
  have st₃ : Step D v _ (cvS₃ isType1 origTstack isSingle s) := Step.mergeTop st₂.inv st₂.shape hok.merge₁
  have R₃ : (cvS₃ isType1 origTstack isSingle s).RInvH dfs v d := R₂.mergeTop
  have st₄ : Step D v _ (cvS₄ isType1 origTstack isSingle s) := Step.mergeTop st₃.inv st₃.shape hok.merge₂
  have R₄ : (cvS₄ isType1 origTstack isSingle s).RInvH dfs v d := R₃.mergeTop
  match hts₄ : (cvS₄ isType1 origTstack isSingle s).tstack with
  | [] => exact absurd hts₄ hok.retarget.nonempty
  | t :: rest =>
    have R₅ : (cvS₅ v edgeDir isType1 origTstack isSingle s).RInvH dfs v d := R₄.retarget edgeDir t rest hts₄
    have hts₅ : (cvS₅ v edgeDir isType1 origTstack isSingle s).tstack =
        { t with vStart := v, spans := setSides (!edgeDir) (t.spans.1 ++ t.spans.2) [] } :: rest := by
      show ((WalkState.retarget v edgeDir).run _).2.tstack = _
      rw [retarget_run_eq v edgeDir _ t rest hts₄]
    cases isType1
    · exact R₅
    · have hf : ItemFree (cvS₂ true origTstack isSingle s) ((maybeUnwrapNxt (if isSingle then NodeType.S else .R)).run s).1 :=
        hfree _ rfl
      have hf₅ := ((hf.merge).merge).retarget v edgeDir
      have hside := (hok.finish rfl).side
      have hcur : curE (cvS₅ v edgeDir true origTstack isSingle s) =
          { t with vStart := v, spans := setSides (!edgeDir) (t.spans.1 ++ t.spans.2) [] } := by
        rw [curE, hts₅]; rfl
      rw [hcur] at hside
      exact R₅.finishTop _ _ rest hts₅ hf₅ hside

/-- The type-1 P-check: when it fires, the entry below the top starts at `v`. -/
theorem RInvH.finishP {D v lowval d : Nat} {isType1 : Bool} (hi : s.Inv' D) (hs : Shape s)
    (hok : FinishPOk D v lowval isType1 s) (h : s.RInvH dfs v d) :
    (after (Spqr.finishP v lowval isType1) s).RInvH dfs v d := by
  show ((Spqr.finishP v lowval isType1).run s).2.RInvH dfs v d
  simp only [Spqr.finishP, WalkM.run_bind, run_condP]
  by_cases hc : result (condP v lowval isType1) s = true
  · have h' : (isType1 && decide (s.tstack.length ≥ 2) && (s.tstack.tail.head!.vStart == v) &&
        (s.tstack.tail.head!.topDepth == lowval)) = true := hc
    simp only [h', ↓reduceIte, WalkM.run_bind]
    obtain ⟨hu, hc2⟩ := hok.ok hc
    have hnxt : (nxtE s).vStart = v := by
      simp only [Bool.and_eq_true, beq_iff_eq] at h'
      exact h'.1.2
    have r := maybeUnwrapNxt_spec (v := v) hi hs (by decide) hu
    have R₁ := h.unwrapNxt hs hu (Exempt.of_vStart hnxt)
    match hts : s.tstack with
    | [] | [_] => have := hu.two; rw [hts] at this; simp at this
    | a :: b :: rest =>
      obtain ⟨b', hts₁, -, -⟩ := maybeUnwrapNxt_tstack .P a b rest hts
      set s₁ := ((maybeUnwrapNxt .P).run s).2 with hs₁
      have hts₁' : s₁.tstack = a :: b' :: rest := hts₁
      have R₂ := R₁.mergeTop
      have hts₂ : (mergeTstackTops.run s₁).2.tstack = TEntry.mergeInto a b' :: rest := by
        rw [mergeTstackTops_run_eq s₁ a b' rest hts₁']
      have hf := r.free.merge
      have hside := hc2.finish.side
      have hcur : curE (after mergeTstackTops (after (maybeUnwrapNxt .P) s)) = TEntry.mergeInto a b' := by
        show (mergeTstackTops.run s₁).2.tstack.head! = _
        rw [hts₂]; rfl
      rw [hcur] at hside
      exact R₂.finishTop _ (TEntry.mergeInto a b') rest hts₂ hf hside
  · have h' : (isType1 && decide (s.tstack.length ≥ 2) && (s.tstack.tail.head!.vStart == v) &&
        (s.tstack.tail.head!.topDepth == lowval)) = false := Bool.eq_false_iff.2 hc
    simp only [h', Bool.false_eq_true, ↓reduceIte]
    exact h

/-- The first-edge vertex push: the old top goes below the vertex entry, so it must be settled
(`htop`); the merge (loop 2 fired) then only touches the top. -/
theorem RInvH.finishTail {v d k : Nat} {hasVert isSingle : Bool}
    (hown : hasVert = false → ∀ t ∈ s.tstack, ∀ e, t.edges s.g s.items e →
      ¬ Items.EdgeBelow s.g s.items (vertItem v) e)
    (htop : hasVert = false → s.SettledTop dfs v d)
    (h : s.RInvH dfs v d) : (after (Spqr.finishTail v k hasVert isSingle) s).RInvH dfs v d := by
  cases hasVert
  · have R₁ := (h.toTop (htop rfl)).pushVert_own (k := k) (hown rfl)
    cases isSingle
    · exact R₁.mergeTop
    · exact R₁
  · exact h

theorem RInvH.finishRest {D v d k lowval : Nat} {isType1 hasVert isSingle : Bool}
    (hi : s.Inv' D) (hs : Shape s)
    (hok : FinishRestOk D v k lowval isType1 hasVert isSingle s)
    (hown : hasVert = false → ∀ t ∈ (after (Spqr.finishP v lowval isType1) s).tstack, ∀ e,
      t.edges (after (Spqr.finishP v lowval isType1) s).g (after (Spqr.finishP v lowval isType1) s).items e →
      ¬ Items.EdgeBelow (after (Spqr.finishP v lowval isType1) s).g (after (Spqr.finishP v lowval isType1) s).items (vertItem v) e)
    (htop : hasVert = false → (after (Spqr.finishP v lowval isType1) s).SettledTop dfs v d)
    (h : s.RInvH dfs v d) :
    (after (Spqr.finishRest v k lowval isType1 hasVert isSingle) s).RInvH dfs v d :=
  (h.finishP hi hs hok.p).finishTail hown htop

/-! ### The tree-edge branch -/

/-- The state after the P-check of the tree-edge branch. -/
def feP (curV d : Nat) (o : DfsOut) (s : WalkState) : WalkState :=
  after (Spqr.finishP curV (o.cls.lowval d) o.cls.isType1) (feS₂ d o s)

/-- What the tree-edge branch of `finishEdge` needs beyond `Frontier`: the tree edge is not yet
owned; after loop 1 every frontier entry below the top is settled (loop 1 runs while
`nxt.topDepth ≥ d`, so these top out above `d` — vertex entries of the child, `EntryR.vert`);
the type-1 `closeVert` unwraps an exempt entry (the one below the top starts at `curV` or tops
out above `d`); and the first-edge vertex entry of `curV` takes edges owned by nobody. -/
structure FinishRShape (dfs : DfsData) (curV d : Nat) (o : DfsOut) (origTstack : Nat)
    (hasVert : Bool) (s : WalkState) : Prop where
  pend : ∀ t ∈ s.tstack, ¬ t.edges s.g s.items o.e
  settled : ∀ t ∈ (feS₁ d o s).tstack.tail.take ((feS₁ d o s).tstack.length - 1 - origTstack),
    d ≤ t.topDepth → t.vStart ≠ curV → (feS₁ d o s).EntryR dfs t
  unwrap : hasVert = true → o.cls.isType1 = true → Exempt curV d (nxtE (feS₂ d o s))
  vert_own : hasVert = false → ∀ t ∈ (feP curV d o s).tstack, ∀ e,
    t.edges (feP curV d o s).g (feP curV d o s).items e →
    ¬ Items.EdgeBelow (feP curV d o s).g (feP curV d o s).items (vertItem curV) e

/-- The top entry (if any) starts at `v`. -/
def TopStart (v : Nat) (s : WalkState) : Prop := ∀ t rest, s.tstack = t :: rest → t.vStart = v

theorem TopStart.settledTop {v d : Nat} (h : TopStart v s) : s.SettledTop dfs v d := by
  intro hne _ hv
  match hts : s.tstack with
  | [] => exact absurd hts hne
  | t :: rest => exact absurd (by rw [curE, hts]; exact h t rest hts) hv

theorem topStart_modifyCur {v : Nat} (f : TEntry → TEntry) (hf : ∀ t, (f t).vStart = t.vStart)
    (h : TopStart v s) : TopStart v (after (modifyCur f) s) := by
  intro t rest hts
  change (match s.tstack with | a :: rest => f a :: rest | [] => []) = t :: rest at hts
  match hs : s.tstack with
  | [] => rw [hs] at hts; cases hts
  | a :: r => rw [hs] at hts; cases hts; rw [hf]; exact h _ _ hs

theorem topStart_cvS₅ {v : Nat} (edgeDir isType1 : Bool) (origTstack : Nat) (isSingle : Bool) :
    TopStart v (cvS₅ v edgeDir isType1 origTstack isSingle s) := by
  intro t rest hts
  change (match (cvS₄ isType1 origTstack isSingle s).tstack with
    | a :: rest => { a with vStart := v, spans := setSides (!edgeDir) (a.spans.1 ++ a.spans.2) [] } :: rest
    | [] => []) = t :: rest at hts
  match hs : (cvS₄ isType1 origTstack isSingle s).tstack with
  | [] => rw [hs] at hts; cases hts
  | a :: r => rw [hs] at hts; cases hts; rfl

theorem topStart_finishTstackTop {v : Nat} (i : ItemId) (h : TopStart v s) :
    TopStart v (after (finishTstackTop i) s) := by
  intro t rest hts
  match hs : s.tstack with
  | [] =>
    change (match s.tstack with | a :: rest => _ :: rest | [] => []) = t :: rest at hts
    rw [hs] at hts; cases hts
  | a :: r =>
    obtain ⟨a', it, he, hv, -, -⟩ := finishTstackTop_run i s hs
    change (((finishTstackTop i).run s).2).tstack = _ at hts
    rw [he] at hts; cases hts; rw [hv]; exact h _ _ hs

theorem topStart_vertFinish {v : Nat} (item : Option ItemId) (b : Bool) (h : TopStart v s) :
    TopStart v (after (vertFinish item b) s) := by
  cases item with
  | none => exact h
  | some i => exact topStart_finishTstackTop i h

theorem topStart_closeVert' {v : Nat} (edgeDir isType1 : Bool) (origTstack : Nat) (isSingle : Bool) :
    TopStart v (after (closeVert' v edgeDir isType1 origTstack isSingle) s) :=
  topStart_vertFinish _ _ (topStart_cvS₅ edgeDir isType1 origTstack isSingle)

theorem topStart_finishP {v lowval : Nat} (isType1 : Bool) (h : TopStart v s) :
    TopStart v (after (Spqr.finishP v lowval isType1) s) := by
  show TopStart v ((Spqr.finishP v lowval isType1).run s).2
  simp only [Spqr.finishP, WalkM.run_bind, run_condP]
  by_cases hc : (isType1 && decide (s.tstack.length ≥ 2) && (s.tstack.tail.head!.vStart == v) &&
      (s.tstack.tail.head!.topDepth == lowval)) = true
  · simp only [hc, ↓reduceIte, WalkM.run_bind]
    simp only [Bool.and_eq_true, decide_eq_true_eq, beq_iff_eq] at hc
    match hts : s.tstack with
    | [] | [_] => rw [hts] at hc; simp at hc
    | a :: b :: rest =>
      have hbv : b.vStart = v := by rw [hts] at hc; exact hc.1.2
      obtain ⟨b', hts₁, hb'v, -⟩ := maybeUnwrapNxt_tstack .P a b rest hts
      have hts₁' : ((maybeUnwrapNxt .P).run s).2.tstack = a :: b' :: rest := hts₁
      have hts₂ : (mergeTstackTops.run ((maybeUnwrapNxt .P).run s).2).2.tstack =
          TEntry.mergeInto a b' :: rest := by
        rw [mergeTstackTops_run_eq _ a b' rest hts₁']
      refine topStart_finishTstackTop _ ?_
      intro t r ht
      rw [hts₂] at ht; cases ht
      show b'.vStart = v
      rw [hb'v, hbv]
  · have h' : (isType1 && decide (s.tstack.length ≥ 2) && (s.tstack.tail.head!.vStart == v) &&
        (s.tstack.tail.head!.topDepth == lowval)) = false := Bool.eq_false_iff.2 hc
    simp only [h', Bool.false_eq_true, ↓reduceIte]
    exact h

theorem topStart_finishRest_vert {v k lowval : Nat} (isType1 isSingle : Bool) (h : TopStart v s) :
    TopStart v (after (Spqr.finishRest v k lowval isType1 true isSingle) s) := by
  show TopStart v ((Spqr.finishRest v k lowval isType1 true isSingle).run s).2
  simp only [Spqr.finishRest, Spqr.finishTail, Bool.not_true, Bool.false_eq_true, ↓reduceIte,
    WalkM.run_bind, WalkM.pure_run]
  exact topStart_finishP isType1 h

theorem finishEdge_tree_vert_topStart (curV d lv : Nat) (kind : RetKind) (o : DfsOut)
    (origTstack : Nat) (ho : o.cls = .ret lv kind) (hk : kind ≠ .backEdge) (hlow : lv < d) :
    TopStart curV (after (finishEdge curV d o origTstack true) s) := by
  have ht : o.cls.isTree = true := by
    rw [ho]; cases kind with
    | backEdge => exact absurd rfl hk
    | type1Child => rfl
    | type2Child => rfl
  have hlv : o.cls.lowval d = lv := by rw [ho]; rfl
  have hge : ¬ (lv ≥ d) := by omega
  rw [finishEdge_eq]
  simp only [finishEdge', finishTree, hlv, hge, ↓reduceIte, ht, closeVert_eq]
  show TopStart curV (after (Spqr.finishRest curV d lv o.cls.isType1 true (feB₃ curV d o origTstack s))
    (feS₃ curV d o origTstack s))
  exact topStart_finishRest_vert (k := d) o.cls.isType1 (feB₃ curV d o origTstack s)
    (topStart_closeVert' _ _ _ _)

/-- Admitted: (Lemma 4.3, R-maximality content proper, first-edge case) the entry the tree-edge branch builds on
top of the stack is `EntryR` whenever it tops out at depth `≥ d` and does not start at `curV`:
after the P-check (`feP`; the first-edge vertex entry is pushed on top of it) and in the output
(the vertex entry, or the old top with the vertex entry merged in). With `hasVert = true` the top
starts at `curV` after `closeVert` and is exempt (`finishEdge_tree_vert_topStart`), so only the
first-edge case is content. Its pieces are the S/P/R items closed by loop 1 (the complement of a closed
`(nxt.vStart, d)` set is one class because all `(nxt.vStart, d)` classes were P-merged when
`nxt.vStart` was finished and the path class is merged here), merged by loop 2 and `closeVert`;
the R close is `loop1_rBranch`. Everything else the branch does keeps the entries below the top
settled (`finishEdge_tree_rInvTop_of_top`). Checked on seeds 0..300 and 6000 extra multigraphs,
both ternarize modes (`checks/RFinishEdgeCheck.lean`, contract B, 6414 sites; `ptop` lines, 278
non-exempt first-edge tops). -/
theorem finishEdge_tree_top_settled_first {D : Nat} (curV d lv : Nat) (kind : RetKind) (o : DfsOut)
    (origTstack : Nat) (ho : o.cls = .ret lv kind) (hk : kind ≠ .backEdge)
    (hlow : lv < d) (hv : curV < s.g.nv)
    (hi : s.Inv' D) (hs : Shape s) (hok : FinishOk D curV d lv o origTstack false s)
    (hg : FinishGuards d o origTstack false s)
    (hfront : Frontier (o := o) d origTstack s)
    (hshape : FinishRShape dfs curV d o origTstack false s)
    (h2 : s.g.TwoConnected) (hsp : dfs.Spec s.g) (hrt : dfs.Rooted s.g)
    (hd : dfs.depth curV = d) (hcur : s.stackVerts[d]! = curV)
    (hanc : ∀ k, k ≤ d → dfs.Anc s.stackVerts[k]! curV ∧ dfs.depth s.stackVerts[k]! = k)
    (hR : s.RInvFront dfs curV d origTstack) :
    (feP curV d o s).SettledTop dfs curV d ∧
      (after (finishEdge curV d o origTstack false) s).SettledTop dfs curV d := by
  sorry

/-- The top of the tree-edge branch's output is settled: for `hasVert = true` the vertex close
re-targets it to `curV` (`finishEdge_tree_vert_topStart`, exempt); the first-edge case is the
admitted `finishEdge_tree_top_settled_first`. -/
theorem finishEdge_tree_top_settled {D : Nat} (curV d lv : Nat) (kind : RetKind) (o : DfsOut)
    (origTstack : Nat) (hasVert : Bool) (ho : o.cls = .ret lv kind) (hk : kind ≠ .backEdge)
    (hlow : lv < d) (hv : curV < s.g.nv)
    (hi : s.Inv' D) (hs : Shape s) (hok : FinishOk D curV d lv o origTstack hasVert s)
    (hg : FinishGuards d o origTstack hasVert s)
    (hfront : Frontier (o := o) d origTstack s)
    (hshape : FinishRShape dfs curV d o origTstack hasVert s)
    (h2 : s.g.TwoConnected) (hsp : dfs.Spec s.g) (hrt : dfs.Rooted s.g)
    (hd : dfs.depth curV = d) (hcur : s.stackVerts[d]! = curV)
    (hanc : ∀ k, k ≤ d → dfs.Anc s.stackVerts[k]! curV ∧ dfs.depth s.stackVerts[k]! = k)
    (hR : s.RInvFront dfs curV d origTstack) :
    (hasVert = false → (feP curV d o s).SettledTop dfs curV d) ∧
      (after (finishEdge curV d o origTstack hasVert) s).SettledTop dfs curV d := by
  cases hasVert
  · have h := finishEdge_tree_top_settled_first curV d lv kind o origTstack ho hk hlow hv hi hs hok
      hg hfront hshape h2 hsp hrt hd hcur hanc hR
    exact ⟨fun _ => h.1, h.2⟩
  · exact ⟨(fun h => nomatch h),
      (finishEdge_tree_vert_topStart curV d lv kind o origTstack ho hk hlow).settledTop⟩

/-- The tree-edge branch modulo `finishEdge_tree_top_settled`: loop 1 touches only the frontier
(`Frontier.loop1`), so the settled base survives positionally (`RInvG`); after it the frontier
below the top is settled (`FinishRShape.settled`, `RInvG.widen`), and loop 2, `closeVert` and
`finishRest` only touch the top (`RInvH`), except the unwrap of an exempt entry and the first-edge
vertex push, which needs the top settled. -/
theorem finishEdge_tree_rInvTop_of_top {D : Nat} (curV d lv : Nat) (kind : RetKind) (o : DfsOut)
    (origTstack : Nat) (hasVert : Bool) (ho : o.cls = .ret lv kind) (hk : kind ≠ .backEdge)
    (hlow : lv < d) (hv : curV < s.g.nv) (hi : s.Inv' D) (hs : Shape s)
    (hok : FinishOk D curV d lv o origTstack hasVert s)
    (hfront : Frontier (o := o) d origTstack s)
    (hshape : FinishRShape dfs curV d o origTstack hasVert s)
    (hR : s.RInvFront dfs curV d origTstack)
    (hP : hasVert = false → (feP curV d o s).SettledTop dfs curV d)
    (hc : (after (finishEdge curV d o origTstack hasVert) s).SettledTop dfs curV d) :
    (after (finishEdge curV d o origTstack hasVert) s).RInvTop dfs curV d := by
  have ht : o.cls.isTree = true := by
    rw [ho]; cases kind with
    | backEdge => exact absurd rfl hk
    | type1Child => rfl
    | type2Child => rfl
  have hlv : o.cls.lowval d = lv := by rw [ho]; rfl
  have hge : ¬ (lv ≥ d) := by omega
  have hlow' : o.cls.lowval d < d := by rw [hlv]; exact hlow
  have st₀ : Step D curV s (feS₀ d o s) :=
    Step.modifyVs hi hs (edgeItem s.g o.e) _ (by show 1 + s.g.nv + o.e < _; have := hok.e_lt; omega)
  have hv₀ : curV < (feS₀ d o s).g.nv := by rw [st₀.g]; exact hv
  set f : Item → Item := fun it =>
    { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) } with hf
  have hfree₀ : ∀ t ∈ s.tstack, edgeItem s.g o.e ∉ t.spans.1 ++ t.spans.2 := fun t ht hmem =>
    hshape.pend t ht ⟨_, hmem, .refl⟩
  have R₀ : (feS₀ d o s).RInvG dfs curV d origTstack :=
    (RInvG.of_front hR).modifyVs_free (edgeItem s.g o.e) f (fun _ => rfl) (fun _ => rfl) hfree₀
  have hown₀ : ∀ t ∈ (feS₀ d o s).tstack, ¬ t.edges (feS₀ d o s).g (feS₀ d o s).items o.e :=
    fun t ht hte => hshape.pend t ht ((TEntry.edges_congr (fun i _ e =>
      Items.Below_modify_ch_eq (edgeItem s.g o.e) f (fun _ => rfl)) o.e).1 hte)
  have hn₀ : origTstack ≤ (feS₀ d o s).tstack.length := hfront.size
  obtain ⟨k, heq, -, -, R₁, -⟩ := R₀.closeEars st₀.inv st₀.shape hv₀ hn₀ (hok.ears ht) hown₀
    (fun k hk hck => (hfront.loop1 ht hlow' k hk).2 hck)
  have st₁ : Step D curV _ (feS₁ d o s) := Step.closeEars st₀.inv st₀.shape hv₀ (hok.ears ht)
  have hv₁ : curV < (feS₁ d o s).g.nv := by rw [st₁.g]; exact hv₀
  have R₁' : (feS₁ d o s).RInvG dfs curV d origTstack := by
    have heq' : feS₁ d o s =
        iter (Spqr.loop1Body d s.stackDir[d]!) k (ceS₁ o.dest d o.e (feS₀ d o s)) := heq
    rw [heq']; exact R₁
  have H₁ : (feS₁ d o s).RInvH dfs curV d := R₁'.widen hshape.settled
  have H₂ : (feS₂ d o s).RInvH dfs curV d := H₁.mergeLate hv₁ st₁.inv st₁.shape (hok.late ht)
  have st₂ : Step D curV _ (feS₂ d o s) := Step.mergeLate st₁.inv st₁.shape hv₁ (hok.late ht)
  have hv₂ : curV < (feS₂ d o s).g.nv := by rw [st₂.g]; exact hv₁
  rw [finishEdge_eq] at hc ⊢
  simp only [finishEdge', finishTree, hlv, hge, ↓reduceIte, ht, closeVert_eq] at hc ⊢
  cases hasVert
  · simp only [Bool.false_eq_true, ↓reduceIte] at hc ⊢
    have hown := hshape.vert_own rfl
    have hP' := hP rfl
    simp only [feP, hlv] at hown hP'
    exact (H₂.finishRest st₂.inv st₂.shape (hok.rest_tree ht rfl) (fun _ => hown) (fun _ => hP')).toTop hc
  · simp only [↓reduceIte] at hc ⊢
    have st₃ : Step D curV _ (feS₃ curV d o origTstack s) :=
      Step.closeVert' st₂.inv st₂.shape hv₂ (hok.vert ht rfl)
    have H₃ : (feS₃ curV d o origTstack s).RInvH dfs curV d :=
      H₂.closeVert' hv₂ st₂.inv st₂.shape (hok.vert ht rfl) (hshape.unwrap rfl)
    exact (H₃.finishRest st₃.inv st₃.shape (hok.rest_vert ht rfl) (fun h => nomatch h)
      (fun h => nomatch h)).toTop hc

/-- Lemma 4.3 per `finishEdge`, tree-edge branch: everything but the `EntryR` of the entry built
on top (`finishEdge_tree_top_settled`). -/
theorem finishEdge_tree_rInvTop {D : Nat} (curV d lv : Nat) (kind : RetKind) (o : DfsOut)
    (origTstack : Nat) (hasVert : Bool) (ho : o.cls = .ret lv kind) (hk : kind ≠ .backEdge)
    (hlow : lv < d) (hv : curV < s.g.nv)
    (hi : s.Inv' D) (hs : Shape s) (hok : FinishOk D curV d lv o origTstack hasVert s)
    (hg : FinishGuards d o origTstack hasVert s)
    (hfront : Frontier (o := o) d origTstack s)
    (hshape : FinishRShape dfs curV d o origTstack hasVert s)
    (h2 : s.g.TwoConnected) (hsp : dfs.Spec s.g) (hrt : dfs.Rooted s.g)
    (hd : dfs.depth curV = d) (hcur : s.stackVerts[d]! = curV)
    (hanc : ∀ k, k ≤ d → dfs.Anc s.stackVerts[k]! curV ∧ dfs.depth s.stackVerts[k]! = k)
    (hR : s.RInvFront dfs curV d origTstack) :
    (after (finishEdge curV d o origTstack hasVert) s).RInvTop dfs curV d :=
  have h := finishEdge_tree_top_settled curV d lv kind o origTstack hasVert ho hk hlow hv hi hs hok
    hg hfront hshape h2 hsp hrt hd hcur hanc hR
  finishEdge_tree_rInvTop_of_top curV d lv kind o origTstack hasVert ho hk hlow hv hi hs hok hfront
    hshape hR h.1 h.2

/-- Lemma 4.3 per `finishEdge`: the back-edge branch is `finishEdge_back_rInvTop`, the tree-edge
branch `finishEdge_tree_rInvTop` (admitting only `finishEdge_tree_top_settled`). `hback` records
the `walkOut` call site of a back edge: the vertex entry was pushed by `walkOutPre` (a back edge
is type 1 with `lv < d`) and `origTstack` was read with nothing pushed since, so the base is the
whole stack. -/
theorem finishEdge_rInvTop {D : Nat} (curV d lv : Nat) (kind : RetKind) (o : DfsOut) (origTstack : Nat)
    (hasVert : Bool) (ho : o.cls = .ret lv kind) (hlow : lv < d) (hv : curV < s.g.nv)
    (hi : s.Inv' D) (hs : Shape s) (hok : FinishOk D curV d lv o origTstack hasVert s)
    (hg : FinishGuards d o origTstack hasVert s)
    (hfront : Frontier (o := o) d origTstack s)
    (hshape : kind ≠ .backEdge → FinishRShape dfs curV d o origTstack hasVert s)
    (hback : kind = .backEdge → hasVert = true ∧ s.tstack.length ≤ origTstack)
    (h2 : s.g.TwoConnected) (hsp : dfs.Spec s.g) (hrt : dfs.Rooted s.g)
    (hd : dfs.depth curV = d) (hcur : s.stackVerts[d]! = curV)
    (hanc : ∀ k, k ≤ d → dfs.Anc s.stackVerts[k]! curV ∧ dfs.depth s.stackVerts[k]! = k)
    (hR : s.RInvFront dfs curV d origTstack) :
    (after (finishEdge curV d o origTstack hasVert) s).RInvTop dfs curV d := by
  by_cases hk : kind = .backEdge
  · subst hk
    obtain ⟨hhv, hlen⟩ := hback rfl
    subst hhv
    have ht' : o.cls.isTree = false := by rw [ho]; rfl
    have h1 : o.cls.isType1 = true := by rw [ho]; rfl
    have hdrop : s.tstack.drop (s.tstack.length - origTstack) = s.tstack := by
      rw [Nat.sub_eq_zero_of_le hlen, List.drop_zero]
    have hRt : s.RInvTop dfs curV d :=
      ⟨fun t ht => hR.base t (by rw [hdrop]; exact ht), hR.disj⟩
    have hown : ∀ t ∈ s.tstack, ¬ t.edges s.g s.items o.e := fun t ht hte =>
      hfront.base_disj t (by rw [hdrop]; exact ht) o.e hok.e_lt hte (Or.inl rfl)
    have hrest := hok.rest_back ht'
    rw [h1] at hrest
    exact finishEdge_back_rInvTop curV d lv o origTstack ho hlow hi hs hok.e_lt (hok.q ht')
      (hok.ends ht') (hok.lv_le ht') hrest hown hRt
  · exact finishEdge_tree_rInvTop curV d lv kind o origTstack hasVert ho hk hlow hv hi hs hok hg
      hfront (hshape hk) h2 hsp hrt hd hcur hanc hR

end WalkState
end Spqr
