import Spqr.Frame
import Spqr.WalkItemsWF

/-! # `I` items hang under `Q` items (`PROOF.md` §7)

`IOk`: every id in a child list or a stack span is allocated and not an `I` item, and every `I`
child has a `Q` parent. The only `I` allocation (`finishBoundary`) puts the new item straight into
the block root's child list; everything else moves non-`I` ids between spans and child lists. -/

namespace Spqr
open WalkM

/-- Every id in the spans of `ts` satisfies `P`. -/
def SpanAll (P : ItemId → Prop) (ts : List TEntry) : Prop :=
  ∀ t ∈ ts, ∀ i ∈ t.spans.1 ++ t.spans.2, P i

theorem mem_setSides_append {i : ItemId} (dir : Bool) (a b : List ItemId) :
    i ∈ (setSides dir a b).1 ++ (setSides dir a b).2 ↔ i ∈ a ++ b := by
  unfold setSides; split <;> simp [or_comm]

theorem mem_getSide {i : ItemId} (p : List ItemId × List ItemId) (dir : Bool) (h : i ∈ getSide p dir) :
    i ∈ p.1 ++ p.2 := by
  unfold getSide at h; split at h <;> simp [h]

namespace SpanAll
variable {P : ItemId → Prop} {ts : List TEntry}

theorem mono {P' : ItemId → Prop} (h : SpanAll P ts) (hP : ∀ i, P i → P' i) : SpanAll P' ts :=
  fun t ht i hi => hP i (h t ht i hi)

theorem tail (h : SpanAll P ts) : SpanAll P ts.tail :=
  fun t ht => h t (List.mem_of_mem_tail ht)

theorem head (h : SpanAll P ts) : ∀ i ∈ ts.head!.spans.1 ++ ts.head!.spans.2, P i := by
  cases ts with
  | nil => intro i hi; simp [List.head!_nil, default_spans] at hi
  | cons t ts => exact h t (by simp)

theorem cons {e : TEntry} (he : ∀ i ∈ e.spans.1 ++ e.spans.2, P i) (h : SpanAll P ts) :
    SpanAll P (e :: ts) := by
  intro t ht
  rcases List.mem_cons.mp ht with rfl | ht
  · exact he
  · exact h t ht

theorem cons_single {a b c : Nat} {dir : Bool} {item : ItemId} (hi : P item) (h : SpanAll P ts) :
    SpanAll P (⟨a, b, c, setSides dir [item] []⟩ :: ts) :=
  cons (fun i hi' => by rw [mem_setSides_append] at hi'; simp at hi'; subst hi'; exact hi) h

theorem merge (h : SpanAll P ts) : SpanAll P (mergeTop ts) := by
  match ts with
  | [] => simp [SpanAll, mergeTop]
  | [_] => simp [SpanAll, mergeTop]
  | b :: a :: rest =>
    intro t ht
    simp only [mergeTop, List.mem_cons] at ht
    rcases ht with rfl | ht
    · intro i hi
      simp only [List.mem_append] at hi
      rcases hi with (hi | hi) | (hi | hi)
      · exact h b (by simp) i (by simp [hi])
      · exact h a (by simp) i (by simp [hi])
      · exact h a (by simp) i (by simp [hi])
      · exact h b (by simp) i (by simp [hi])
    · exact h t (by simp [ht])

theorem modifyCur (f : TEntry → TEntry) (h : SpanAll P ts)
    (hf : ∀ t ∈ ts, ∀ i ∈ (f t).spans.1 ++ (f t).spans.2, P i) :
    SpanAll P (match (generalizing := false) ts with | a :: rest => f a :: rest | [] => []) := by
  cases ts with
  | nil => simp [SpanAll]
  | cons a rest =>
    intro t ht
    rcases List.mem_cons.mp ht with rfl | ht
    · exact hf _ (List.mem_cons_self ..)
    · exact h _ (List.mem_cons_of_mem _ ht)

theorem modifyNxt (f : TEntry → TEntry) (h : SpanAll P ts)
    (hf : ∀ t ∈ ts, ∀ i ∈ (f t).spans.1 ++ (f t).spans.2, P i) :
    SpanAll P (match (generalizing := false) ts with | a :: b :: rest => a :: f b :: rest | l => l) := by
  match ts with
  | [] => exact h
  | [_] => exact h
  | a :: b :: rest =>
    intro t ht
    rcases List.mem_cons.mp ht with rfl | ht
    · exact h _ (List.mem_cons_self ..)
    rcases List.mem_cons.mp ht with rfl | ht
    · exact hf b (by simp)
    · exact h _ (by simp [ht])

end SpanAll

/-- `I` items hang under `Q` items, plus the bookkeeping that keeps it so. -/
structure IOk (g : Graph) (s : WalkState) : Prop where
  g_eq : s.g = g
  size : 1 + g.nv + g.ne ≤ s.items.size
  vert : ∀ v, v < g.nv → Items.type s.items (vertItem v) = .V
  edge : ∀ e, e < g.ne → Items.type s.items (edgeItem g e) = .Q
  ch_lt : ∀ p c, c ∈ Items.ch s.items p → c < s.items.size
  spans : SpanAll (fun i => i < s.items.size ∧ Items.type s.items i ≠ .I) s.tstack
  parent : ∀ p c, c ∈ Items.ch s.items p → Items.type s.items c = .I → Items.type s.items p = .Q

theorem type_modify_vs_ch (items : Items) (i j : Nat) (vs : Option Nat × Option Nat) (ch : List ItemId) :
    Items.type (items.modify i fun it => { it with vs := vs, ch := ch }) j = Items.type items j :=
  Items.type_modify _ _ _ _ (fun _ => rfl)

theorem ch_getElem_eq (items : Items) {p : Nat} (hp : p < items.size) : items[p].ch = Items.ch items p := by
  simp [Items.ch, Array.getElem?_eq_getElem hp]

namespace IOk
variable {g : Graph} {s s' : WalkState}

theorem of_eq (h : IOk g s) (hg : s'.g = s.g) (hi : s'.items = s.items) (ht : s'.tstack = s.tstack) :
    IOk g s' := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;> simp only [hg, hi, ht]
  exacts [h.g_eq, h.size, h.vert, h.edge, h.ch_lt, h.spans, h.parent]

theorem tstack (h : IOk g s) (ts : List TEntry)
    (hs : SpanAll (fun i => i < s.items.size ∧ Items.type s.items i ≠ .I) ts) :
    IOk g { s with tstack := ts } :=
  ⟨h.g_eq, h.size, h.vert, h.edge, h.ch_lt, hs, h.parent⟩

theorem good_vert (h : IOk g s) {v : Nat} (hv : v < g.nv) :
    vertItem v < s.items.size ∧ Items.type s.items (vertItem v) ≠ .I :=
  ⟨by have := h.size; show 1 + v < _; omega, by rw [h.vert v hv]; decide⟩

theorem good_edge (h : IOk g s) {e : Nat} (he : e < g.ne) :
    edgeItem s.g e < s.items.size ∧ Items.type s.items (edgeItem s.g e) ≠ .I := by
  rw [h.g_eq]
  exact ⟨by have := h.size; show 1 + g.nv + e < _; omega, by rw [h.edge e he]; decide⟩

theorem edge_lt (h : IOk g s) {e : Nat} (he : e < g.ne) : edgeItem s.g e < s.items.size :=
  (h.good_edge he).1

theorem edge_type (h : IOk g s) {e : Nat} (he : e < g.ne) : Items.type s.items (edgeItem s.g e) = .Q := by
  rw [h.g_eq]; exact h.edge e he

theorem push_vert (h : IOk g s) {v : Nat} (hv : v < g.nv) (d : Nat) (dir : Bool) :
    IOk g { s with tstack := ⟨v, d, s.nxtEdgeIdx, setSides dir [vertItem v] []⟩ :: s.tstack } :=
  h.tstack _ (SpanAll.cons_single (h.good_vert hv) h.spans)

theorem push_edge (h : IOk g s) {e : Nat} (he : e < g.ne) (v d : Nat) (dir : Bool) :
    IOk g { s with tstack := ⟨v, d, s.nxtEdgeIdx, setSides dir [edgeItem s.g e] []⟩ :: s.tstack } :=
  h.tstack _ (SpanAll.cons_single (h.good_edge he) h.spans)

theorem merge (h : IOk g s) : IOk g { s with tstack := mergeTop s.tstack } :=
  h.tstack _ (SpanAll.merge h.spans)

theorem pop (h : IOk g s) : IOk g { s with tstack := s.tstack.tail } :=
  h.tstack _ (SpanAll.tail h.spans)

theorem stackVerts (h : IOk g s) (a : Array Nat) : IOk g { s with stackVerts := a } :=
  ⟨h.g_eq, h.size, h.vert, h.edge, h.ch_lt, h.spans, h.parent⟩

theorem stackDir (h : IOk g s) (a : Array Bool) : IOk g { s with stackDir := a } :=
  ⟨h.g_eq, h.size, h.vert, h.edge, h.ch_lt, h.spans, h.parent⟩

theorem firstOcc (h : IOk g s) (a : Array Nat) : IOk g { s with firstOccurrence := a } :=
  ⟨h.g_eq, h.size, h.vert, h.edge, h.ch_lt, h.spans, h.parent⟩

theorem firstOccNxt (h : IOk g s) (a : Array Nat) (b : Nat) :
    IOk g { s with firstOccurrence := a, nxtEdgeIdx := b } :=
  ⟨h.g_eq, h.size, h.vert, h.edge, h.ch_lt, h.spans, h.parent⟩

theorem totBlocks (h : IOk g s) : IOk g { s with totBlocks := s.totBlocks + 1 } :=
  ⟨h.g_eq, h.size, h.vert, h.edge, h.ch_lt, h.spans, h.parent⟩

theorem totSelfLoops (h : IOk g s) : IOk g { s with totSelfLoops := s.totSelfLoops + 1 } :=
  ⟨h.g_eq, h.size, h.vert, h.edge, h.ch_lt, h.spans, h.parent⟩

theorem modifyCur (h : IOk g s) (f : TEntry → TEntry)
    (hf : ∀ t ∈ s.tstack, ∀ i ∈ (f t).spans.1 ++ (f t).spans.2,
      i < s.items.size ∧ Items.type s.items i ≠ .I) :
    IOk g { s with tstack := match s.tstack with | a :: rest => f a :: rest | [] => [] } :=
  h.tstack _ (SpanAll.modifyCur f h.spans hf)

theorem modifyNxt (h : IOk g s) (f : TEntry → TEntry)
    (hf : ∀ t ∈ s.tstack, ∀ i ∈ (f t).spans.1 ++ (f t).spans.2,
      i < s.items.size ∧ Items.type s.items i ≠ .I) :
    IOk g { s with tstack := match s.tstack with | a :: b :: rest => a :: f b :: rest | l => l } :=
  h.tstack _ (SpanAll.modifyNxt f h.spans hf)

theorem modify_vs (h : IOk g s) (i : Nat) (f : Item → Item) (ht : ∀ it, (f it).type = it.type)
    (hc : ∀ it, (f it).ch = it.ch) : IOk g { s with items := s.items.modify i f } := by
  have e1 : ∀ j, Items.type (s.items.modify i f) j = Items.type s.items j :=
    fun j => Items.type_modify _ _ _ _ ht
  have e2 : ∀ j, Items.ch (s.items.modify i f) j = Items.ch s.items j :=
    fun j => Items.ch_modify_of_ch _ _ _ _ hc
  refine ⟨h.g_eq, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;> simp only [e1, e2, Array.size_modify]
  exacts [h.size, h.vert, h.edge, h.ch_lt, h.spans, h.parent]

theorem modify_ch (h : IOk g s) (p : Nat) (f : Item → Item) (ht : ∀ it, (f it).type = it.type)
    (hp : p < s.items.size)
    (hnew : ∀ c ∈ (f s.items[p]).ch,
      c < s.items.size ∧ (Items.type s.items c = .I → Items.type s.items p = .Q)) :
    IOk g { s with items := s.items.modify p f } := by
  have e1 : ∀ j, Items.type (s.items.modify p f) j = Items.type s.items j :=
    fun j => Items.type_modify _ _ _ _ ht
  have e2 : ∀ j, Items.ch (s.items.modify p f) j =
      if j = p then (f s.items[p]).ch else Items.ch s.items j := by
    intro j; split
    · subst j; exact Items.ch_modify_self _ _ _ hp
    · exact Items.ch_modify_ne _ _ _ _ (Ne.symm ‹_›)
  refine ⟨h.g_eq, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;> simp only [e1, e2, Array.size_modify]
  · exact h.size
  · exact h.vert
  · exact h.edge
  · intro q c hc; split at hc
    · exact (hnew c hc).1
    · exact h.ch_lt q c hc
  · exact h.spans
  · intro q c hc hI; split at hc
    · subst ‹q = p›; exact (hnew c hc).2 hI
    · exact h.parent q c hc hI

theorem modify_vert_ch (h : IOk g s) {v : Nat} (hv : v < g.nv) (f : Item → Item)
    (ht : ∀ it, (f it).type = it.type)
    (hnew : IOk g s → ∀ c ∈ (f (s.items[vertItem v]'(h.good_vert hv).1)).ch,
      c < s.items.size ∧ Items.type s.items c ≠ .I) :
    IOk g { s with items := s.items.modify (vertItem v) f } :=
  h.modify_ch _ f ht (h.good_vert hv).1 fun c hc =>
    ⟨(hnew h c hc).1, fun hI => absurd hI (hnew h c hc).2⟩

theorem modify_edge_ch (h : IOk g s) {e : Nat} (he : e < g.ne) (f : Item → Item)
    (ht : ∀ it, (f it).type = it.type)
    (hnew : IOk g s → ∀ c ∈ (f (s.items[edgeItem s.g e]'(h.edge_lt he))).ch, c < s.items.size) :
    IOk g { s with items := s.items.modify (edgeItem s.g e) f } :=
  h.modify_ch _ f ht (h.edge_lt he) fun c hc => ⟨hnew h c hc, fun _ => h.edge_type he⟩

theorem vert_append_ok (h : IOk g s) {v e : Nat} (hv : v < g.nv) (he : e < g.ne) :
    ∀ c ∈ (s.items[vertItem v]'(h.good_vert hv).1).ch ++ [edgeItem s.g e],
      c < s.items.size ∧ Items.type s.items c ≠ .I := by
  intro c hc
  rw [ch_getElem_eq] at hc
  simp only [List.mem_append, List.mem_singleton] at hc
  rcases hc with hc | rfl
  · exact ⟨h.ch_lt _ c hc, fun hI => by
      have := h.parent _ c hc hI; rw [h.vert v hv] at this; cases this⟩
  · exact h.good_edge he

theorem push (h : IOk g s) (x : Item) (hx : x.ch = []) : IOk g { s with items := s.items.push x } := by
  have e1 : ∀ j, j < s.items.size → Items.type (s.items.push x) j = Items.type s.items j := by
    intro j hj; rw [Items.type_push, ite_eq_right (Nat.ne_of_lt hj)]
  have e2 : ∀ j, Items.ch (s.items.push x) j = if j = s.items.size then [] else Items.ch s.items j := by
    intro j; rw [Items.ch_push, hx]
  have hsz : s.items.size ≤ (s.items.push x).size := by simp
  refine ⟨h.g_eq, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · exact Nat.le_trans h.size hsz
  · intro v hv
    show Items.type (s.items.push x) (vertItem v) = .V
    rw [e1 _ (h.good_vert hv).1]; exact h.vert v hv
  · intro e he
    show Items.type (s.items.push x) (edgeItem g e) = .Q
    rw [e1 _ (by have := h.size; show 1 + g.nv + e < _; omega)]; exact h.edge e he
  · intro q c hc
    change c ∈ Items.ch (s.items.push x) q at hc
    rw [e2] at hc; split at hc
    · simp at hc
    · exact Nat.lt_of_lt_of_le (h.ch_lt q c hc) hsz
  · intro t ht i hi
    obtain ⟨hlt, hne⟩ := h.spans t ht i hi
    exact ⟨Nat.lt_of_lt_of_le hlt hsz, by
      show Items.type (s.items.push x) i ≠ .I
      rw [e1 _ hlt]; exact hne⟩
  · intro q c hc hI
    change c ∈ Items.ch (s.items.push x) q at hc
    change Items.type (s.items.push x) c = .I at hI
    show Items.type (s.items.push x) q = .Q
    rw [e2] at hc; split at hc
    · simp at hc
    · have hq : q < s.items.size :=
        Nat.lt_of_not_le fun hle => by rw [Items.ch_of_le _ _ hle] at hc; simp at hc
      rw [e1 _ (h.ch_lt q c hc)] at hI
      rw [e1 _ hq]
      exact h.parent q c hc hI

end IOk

section
variable {g : Graph} {s : WalkState}

theorem iOk_mergeTstackTops (h : IOk g s) : wp mergeTstackTops (fun _ s' => IOk g s') s := by
  simp only [wp_mergeTstackTops]; exact h.merge

theorem iOk_finishTstackTop (item : ItemId) (hi : item < s.items.size)
    (hI : Items.type s.items item ≠ .I) (h : IOk g s) :
    wp (finishTstackTop item) (fun _ s' => IOk g s') s := by
  unfold finishTstackTop
  simp (config := { proj := false }) only [wp_bind, wp_cur, wp_stackDir, wp_makeVs, wp_modifyItem,
    wp_modifyCur]
  have hhd := SpanAll.head h.spans
  refine (h.modify_ch item _ (by intros; rfl) hi ?_).modifyCur _ ?_
  · intro c hc
    have := hhd c (mem_getSide _ _ hc)
    exact ⟨this.1, fun hI' => absurd hI' this.2⟩
  · intro t _ i hi'
    simp only [mem_setSides_append, List.append_nil, List.mem_singleton] at hi'
    rw [hi']
    exact ⟨by simpa using hi, by rw [type_modify_vs_ch]; exact hI⟩

theorem iOk_maybeUnwrapNxt (ty : NodeType) (hty : ty = .S ∨ ty = .P ∨ ty = .R) (h : IOk g s) :
    wp (maybeUnwrapNxt ty)
      (fun item s' => IOk g s' ∧ item < s'.items.size ∧ Items.type s'.items item ≠ .I) s := by
  have hF : ty ≠ .F := by rcases hty with rfl | rfl | rfl <;> decide
  have hQ : ty ≠ .Q := by rcases hty with rfl | rfl | rfl <;> decide
  have hI : ty ≠ .I := by rcases hty with rfl | rfl | rfl <;> decide
  have alloc : IOk g { s with items := s.items.push ⟨ty, (none, none), []⟩ } ∧
      s.items.size < (s.items.push ⟨ty, (none, none), []⟩).size ∧
      Items.type (s.items.push ⟨ty, (none, none), []⟩) s.items.size ≠ .I :=
    ⟨h.push _ rfl, by simp, by rw [Items.type_push, ite_eq_left rfl]; exact hI⟩
  unfold maybeUnwrapNxt
  simp only [wp_bind, wp_ite, wp_get, wp_allocItem, wp_pure, wp_nxt, wp_stackDir, wp_getItem,
    wp_modifyNxt]
  split
  · exact alloc
  · split
    · rename_i heq
      have heq' := beq_iff_eq.mp heq
      rw [Items.type_getElem!'] at heq'
      have hlt : _ < s.items.size := Nat.lt_of_not_le fun hle => hF (heq' ▸ Items.type_of_le _ _ hle)
      refine ⟨h.modifyNxt _ fun t _ i hi => ?_, hlt, heq' ▸ hI⟩
      rw [mem_setSides_append] at hi; simp only [List.append_nil] at hi
      rw [Items.ch_getElem! _ _ hlt] at hi
      exact ⟨h.ch_lt _ i hi, fun hI' => hQ (heq' ▸ h.parent _ i hi hI')⟩
    · exact alloc

theorem iOk_loop1Body (d : Nat) (edgeDir : Bool) (h : IOk g s) :
    wp (loop1Body d edgeDir) (fun _ s' => IOk g s') s := by
  unfold loop1Body loop1Type
  simp only [wp_bind, wp_ite, wp_nxt, wp_cur, wp_setStackDir, wp_mergeTstackTops, wp_pure]
  have fin : ∀ (ty : NodeType), ty = .S ∨ ty = .P ∨ ty = .R → ∀ s₁, IOk g s₁ →
      wp (maybeUnwrapNxt ty) (fun item s₂ => wp (finishTstackTop item) (fun _ s' => IOk g s')
        { s₂ with tstack := mergeTop s₂.tstack }) s₁ := by
    intro ty hty s₁ h₁
    exact wp_mono _ (iOk_maybeUnwrapNxt ty hty h₁) fun item s₂ ⟨h₂, hlt, hI⟩ =>
      iOk_finishTstackTop item hlt hI h₂.merge
  split
  · exact fin .S (by simp) _ (h.stackDir (s.stackDir.set! s.tstack.tail.head!.topDepth edgeDir)).merge
  · split
    · exact fin .P (by simp) _ h
    · exact fin .R (by simp) _ h

theorem iOk_loop (n : Nat) (cond : WalkM Bool) (body : WalkM Unit) (hc : ∀ s, (cond.run s).2 = s)
    (hb : ∀ s, IOk g s → wp body (fun _ s' => IOk g s') s) (h : IOk g s) :
    wp (loop n cond body) (fun _ s' => IOk g s') s :=
  wp_loop (fun s' => IOk g s') n cond body
    (fun s₁ h₁ => by show IOk g (cond.run s₁).2; rw [hc]; exact h₁) hb h fun _ h => h

theorem iOk_finishTail (curV d : Nat) (hasVert isSingle : Bool) (hv : curV < g.nv) (h : IOk g s) :
    wp (finishTail curV d hasVert isSingle) (fun _ s' => IOk g s') s := by
  unfold finishTail
  simp only [wp_bind, wp_ite, wp_pure, wp_pushVertTstack, wp_mergeTstackTops]
  split
  · split
    · exact (h.push_vert hv d _).merge
    · exact h.push_vert hv d _
  · exact h

theorem iOk_finishRest (curV d lowval : Nat) (isType1 hasVert isSingle : Bool) (hv : curV < g.nv)
    (h : IOk g s) :
    wp (finishRest curV d lowval isType1 hasVert isSingle) (fun _ s' => IOk g s') s := by
  unfold finishRest finishP condP
  simp only [wp_bind, wp_ite, wp_tstackSize, wp_nxt, wp_pure, wp_mergeTstackTops]
  split
  · refine wp_mono _ (iOk_maybeUnwrapNxt .P (by simp) h) fun item s₁ ⟨h₁, hlt, hI⟩ => ?_
    exact wp_mono _ (iOk_finishTstackTop item hlt hI h₁.merge) fun _ s₂ h₂ =>
      iOk_finishTail _ _ _ _ hv h₂
  · exact iOk_finishTail _ _ _ _ hv h

theorem iOk_closeVertTail (curV : Nat) (edgeDir isSingle : Bool) (k : Bool → WalkM Bool)
    (item : Option ItemId)
    (hitem : ∀ i, item = some i → i < s.items.size ∧ Items.type s.items i ≠ .I)
    (hk : ∀ b s, IOk g s → wp (k b) (fun _ s' => IOk g s') s) (h : IOk g s) :
    wp (closeVertTail curV edgeDir isSingle k item) (fun _ s' => IOk g s') s := by
  unfold closeVertTail
  have h₃ := h.merge.merge.modifyCur
    (fun t => { t with vStart := curV, spans := setSides (!edgeDir) (t.spans.1 ++ t.spans.2) [] })
    fun t ht i hi => by
      simp only [mem_setSides_append, List.append_nil] at hi
      exact h.merge.merge.spans t ht i hi
  cases item with
  | none =>
    simp only [wp_bind, wp_mergeTstackTops, wp_modifyCur]
    exact hk _ _ h₃
  | some item =>
    simp only [wp_bind, wp_mergeTstackTops, wp_modifyCur]
    obtain ⟨hlt, hI⟩ := hitem item rfl
    exact wp_mono _ (iOk_finishTstackTop item hlt hI h₃) fun _ s₂ h₂ => hk _ _ h₂

theorem iOk_closeVert (curV : Nat) (edgeDir isType1 : Bool) (origTstack : Nat) (isSingle : Bool)
    (k : Bool → WalkM Bool)
    (hk : ∀ b s, IOk g s → wp (k b) (fun _ s' => IOk g s') s) (h : IOk g s) :
    wp (closeVert curV edgeDir isType1 origTstack isSingle k) (fun _ s' => IOk g s') s := by
  unfold closeVert
  simp only [wp_bind, wp_ite, wp_tstackSize, wp_map]
  split
  · refine wp_mono _ (iOk_loop _ _ _ (fun s => by rw [run_loop3Cond])
      (fun s h => iOk_mergeTstackTops h) h) fun _ s₃ h₃ => ?_
    exact iOk_closeVertTail _ _ _ _ _ (fun _ h => by cases h) hk h₃
  · refine wp_mono _ (iOk_maybeUnwrapNxt _ (by split <;> simp) h) fun item s₃ ⟨h₃, hlt, hI⟩ => ?_
    exact iOk_closeVertTail _ _ _ _ _ (fun i hi => by cases hi; exact ⟨hlt, hI⟩) hk h₃

theorem iOk_finishTree (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert edgeDir : Bool)
    (hv : curV < g.nv) (he : o.e < g.ne) (h : IOk g s) :
    wp (finishTree curV d o origTstack hasVert edgeDir) (fun _ s' => IOk g s') s := by
  unfold finishTree closeEars mergeLate
  simp only [wp_bind, wp_pushEdgeTstack, wp_tstackSize, wp_get, wp_cur, wp_ite, wp_pure]
  refine wp_mono _ (iOk_loop _ _ _ (fun s => by rw [run_loop1Cond])
    (fun s h => iOk_loop1Body d _ h) (h.push_edge he _ _ _)) fun _ s₁ h₁ => ?_
  have fin : ∀ (b : Bool) (s₂ : WalkState), IOk g s₂ →
      if hasVert = true then
        wp (closeVert curV edgeDir o.cls.isType1 origTstack b
          (finishRest curV d (o.cls.lowval d) o.cls.isType1 hasVert)) (fun _ s' => IOk g s') s₂
      else wp (finishRest curV d (o.cls.lowval d) o.cls.isType1 hasVert b)
        (fun _ s' => IOk g s') s₂ := by
    intro b s₂ h₂
    split
    · exact iOk_closeVert _ _ _ _ _ _ (fun b s h => iOk_finishRest _ _ _ _ _ _ hv h) h₂
    · exact iOk_finishRest _ _ _ _ _ _ hv h₂
  split
  · exact wp_mono _ (iOk_loop _ _ _ (fun s => by rw [run_loop2Cond])
      (fun s h => iOk_mergeTstackTops h) h₁) fun _ s₂ h₂ => fin false s₂ h₂
  · exact fin true s₁ h₁

theorem iOk_finishBoundary (curV d : Nat) (o : DfsOut) (hasVert : Bool)
    (hv : curV < g.nv) (he : o.e < g.ne) (h : IOk g s) :
    wp (finishBoundary curV d o (edgeItem s.g o.e) hasVert) (fun _ s' => IOk g s') s := by
  unfold finishBoundary
  simp (config := { proj := false }) only [wp_bind, wp_modifyItem, wp_modify, wp_ite, wp_allocItem,
    wp_makeVs, wp_popTstack, wp_pure]
  have h₁ := (h.modify_vs (edgeItem s.g o.e) (fun it => { it with vs := (some curV, none) })
    (by intros; rfl) (by intros; rfl)).totBlocks
  split
  · split
    · -- bridge: the new `I` item goes straight under the block root `Q`
      have h₂ := h₁.push ⟨.I, (none, none), []⟩ rfl
      refine IOk.modify_vert_ch ?_ hv _ (by intros; rfl) (fun h₄ => h₄.vert_append_ok hv he)
      refine IOk.modify_edge_ch ((h₂.modify_vs _ _ (by intros; rfl) (by intros; rfl)).pop) he _
        (by intros; rfl) ?_
      intro _ c hc
      simp only [List.mem_cons] at hc
      rcases hc with rfl | hc
      · simp
      · have := SpanAll.head h₂.spans c (by simp [hc])
        simpa using this.1
    · refine IOk.modify_vert_ch ?_ hv _ (by intros; rfl) (fun h₄ => h₄.vert_append_ok hv he)
      refine IOk.modify_edge_ch h₁.pop.pop he _ (by intros; rfl) ?_
      intro _ c hc
      simp only [List.mem_append] at hc
      rcases hc with hc | hc
      · exact (SpanAll.head h₁.spans c (by simp [hc])).1
      · exact (SpanAll.head (SpanAll.tail h₁.spans) c (by simp [hc])).1
  · refine IOk.modify_vert_ch ?_ hv _ (by intros; rfl) (fun h₄ => h₄.vert_append_ok hv he)
    refine IOk.modify_edge_ch ((h₁.totSelfLoops.push ⟨.O, (none, none), []⟩ rfl).modify_vs _ _
      (by intros; rfl) (by intros; rfl)) he _ (by intros; rfl) ?_
    intro _ c hc
    simp only [List.mem_singleton] at hc
    subst hc
    simp

theorem iOk_finishEdge (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool)
    (hv : curV < g.nv) (he : o.e < g.ne) (h : IOk g s) :
    wp (finishEdge curV d o origTstack hasVert) (fun _ s' => IOk g s') s := by
  by_cases hge : o.cls.lowval d ≥ d
  · simp only [finishEdge_eq, finishEdge', wp_bind, wp_get, wp_stackDir, hge, ↓reduceIte]
    exact iOk_finishBoundary _ _ _ _ hv he h
  · by_cases ht : o.cls.isTree = true
    · simp only [finishEdge_eq, finishEdge', wp_bind, wp_get, wp_stackDir, hge, ht, ↓reduceIte,
        wp_makeVs, wp_modifyItem]
      exact iOk_finishTree _ _ _ _ _ _ hv he (h.modify_vs _ _ (by intros; rfl) (by intros; rfl))
    · have ht' : o.cls.isTree = false := Bool.eq_false_iff.2 ht
      simp only [finishEdge_eq, finishEdge', finishBack, wp_bind, wp_get, wp_stackDir, hge, ht',
        Bool.false_eq_true, ↓reduceIte, wp_makeVs, wp_modifyItem, wp_pushEdgeTstack, wp_modify]
      refine iOk_finishRest _ _ _ _ _ _ hv ?_
      exact IOk.firstOccNxt (IOk.push_edge (IOk.modify_vs h (edgeItem s.g o.e)
        (fun it => { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) })
        (by intros; rfl) (by intros; rfl)) he _ _ _) _ _

theorem iOk_walk_aux (g : Graph) :
    (∀ (t : DfsTree) (d : Nat) (s : WalkState), IOk g s → t.Bounded g.nv g.ne →
      wp (walkTree t d) (fun _ s' => IOk g s') s) ∧
    (∀ (v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState),
      IOk g s → v < g.nv → DfsOut.BoundedList g.nv g.ne outs →
      wp (walkOuts v d outs hasVert) (fun _ s' => IOk g s') s) ∧
    (∀ (v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState),
      IOk g s → v < g.nv → o.Bounded g.nv g.ne →
      wp (walkOut v d o hasVert) (fun _ s' => IOk g s') s) := by
  refine walkTree.mutual_induct _ _ _ ?_ ?_ ?_ ?_
  · intro d v outs ih s h hb
    obtain ⟨hv, hb⟩ := hb
    simp only [walkTree, wp_bind, wp_modify]
    refine wp_mono _ (ih _ (h.stackVerts (s.stackVerts.set! d v)) hv hb) fun hasVert s₁ h₁ => ?_
    split
    · exact h₁
    · simp only [wp_bind, wp_setStackDir, wp_pushVertTstack]
      exact (h₁.stackDir (s₁.stackDir.set! d true)).push_vert hv d _
  · intro o v d hasVert ih s h hv hb
    have hdir := h.stackDir (s.stackDir.set! d
      (if o.cls.lowval d ≥ d then false else !s.stackDir[o.cls.lowval d]!))
    cases o with
    | back e dest cls =>
      unfold walkOut
      simp only [wp_bind, wp_stackDir, wp_setStackDir, wp_ite, wp_pushVertTstack, wp_pure,
        wp_tstackSize]
      split
      · exact iOk_finishEdge _ _ _ _ _ hv hb.1 (hdir.push_vert hv d _)
      · exact iOk_finishEdge _ _ _ _ _ hv hb.1 hdir
    | tree e cls child =>
      dsimp only at ih
      unfold walkOut
      simp only [wp_bind, wp_stackDir, wp_setStackDir, wp_ite, wp_pushVertTstack, wp_pure,
        wp_tstackSize, wp_modify]
      split
      · exact wp_mono _ (ih _ ((hdir.push_vert hv d _).firstOcc (s.firstOccurrence.set! d s.g.ne)) hb.2)
          fun _ s₃ h₃ => iOk_finishEdge _ _ _ _ _ hv hb.1 h₃
      · exact wp_mono _ (ih _ (hdir.firstOcc (s.firstOccurrence.set! d s.g.ne)) hb.2)
          fun _ s₃ h₃ => iOk_finishEdge _ _ _ _ _ hv hb.1 h₃
  · intro v d hasVert s h hv hb
    exact h
  · intro v d hasVert o rest ih₁ ih₂ s h hv hb
    simp only [walkOuts, wp_bind]
    exact wp_mono _ (ih₁ s h hv hb.1) fun hasVert₁ s₁ h₁ => ih₂ hasVert₁ s₁ h₁ hv hb.2

theorem iOk_walkForest (forest : List DfsTree) (hb : ∀ t ∈ forest, t.Bounded g.nv g.ne) (h : IOk g s) :
    wp (walkForest forest) (fun _ s' => IOk g s') s := by
  induction forest generalizing s with
  | nil => exact h
  | cons t rest ih =>
    show wp ((walkTree t 0 >>= fun _ => popTstack >>= fun top =>
      modifyItem rootItem fun it => { it with ch := it.ch ++ top.spans.2 }) >>= fun _ => walkForest rest) _ s
    simp only [wp_bind, wp_popTstack, wp_modifyItem]
    refine wp_mono _ ((iOk_walk_aux g).1 t 0 s h (hb t (by simp))) fun _ s₁ h₁ => ?_
    have h₂ := h₁.pop
    have hroot : rootItem < s₁.items.size := by
      have := h₁.size; show 0 < s₁.items.size; omega
    refine ih (fun t ht => hb t (by simp [ht])) (h₂.modify_ch rootItem _ (by intros; rfl) hroot ?_)
    intro c hc
    simp only [List.mem_append, ch_getElem_eq] at hc
    rcases hc with hc | hc
    · exact ⟨h₂.ch_lt _ c hc, h₂.parent _ c hc⟩
    · have := SpanAll.head h₁.spans c (by simp [hc])
      exact ⟨this.1, fun hI => absurd hI this.2⟩

theorem iOk_init (g : Graph) (tern : Bool) : IOk g (WalkState.init g tern) where
  g_eq := rfl
  size := Nat.le_of_eq (Items.initialItems_size g).symm
  vert v hv := by
    have h1 : (1 + v : Nat) < 1 + g.nv := by omega
    show Items.type (initialItems g) (vertItem v) = .V
    rw [Items.initialItems_type]; simp [vertItem, h1]
  edge e he := by
    have h1 : ¬ (1 + g.nv + e : Nat) < 1 + g.nv := by omega
    have h2 : (1 + g.nv + e : Nat) < 1 + g.nv + g.ne := by omega
    show Items.type (initialItems g) (edgeItem g e) = .Q
    rw [Items.initialItems_type]; simp [edgeItem, h1, h2]
  ch_lt p c hc := by
    change c ∈ Items.ch (initialItems g) p at hc
    rw [Items.initialItems_ch] at hc; simp at hc
  spans t ht := by simp [WalkState.init] at ht
  parent p c hc := by
    change c ∈ Items.ch (initialItems g) p at hc
    rw [Items.initialItems_ch] at hc; simp at hc

end

/-- An `I` item hangs under a `Q` item. -/
theorem walk_i_parent' (g : Graph) (tern : Bool) (forest : List DfsTree)
    (hb : ∀ t ∈ forest, t.Bounded g.nv g.ne) :
    ∀ p c, Items.IsParent (g.walk tern forest).items p c →
      Items.type (g.walk tern forest).items c = .I →
      Items.type (g.walk tern forest).items p = .Q :=
  (iOk_walkForest forest hb (iOk_init g tern)).parent

end Spqr
