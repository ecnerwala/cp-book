import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Spqr.WalkTyping

/-!
# Span discipline of the walk

Every item id is in at most one place: counted over all tstack span lists and all `ch` lists of
allocated items, each id occurs at most once, everything that occurs is allocated, the root never
occurs, and a `vertItem`/`edgeItem` occurs only once its vertex has been pushed / its edge
finished. This is the "at most once" half of the tree shape of `Items.WF` (`ch_lt`, `ch_nodup`,
`root_no_parent`, uniqueness of the parent); it holds regardless of the tstack semantics.
-/

namespace Spqr

/-- Occurrences of `i` in the span lists of the tstack. -/
def spansCount (ts : List TEntry) (i : ItemId) : Nat :=
  (ts.map fun t => (t.spans.1 ++ t.spans.2).count i).sum

/-- Occurrences of `i` in the child lists of the allocated items. -/
def chCount (items : Items) (i : ItemId) : Nat :=
  ∑ j ∈ Finset.range items.size, (items.ch j).count i

/-- Occurrences of `i` anywhere (spans or child lists). -/
def WalkState.cnt (s : WalkState) (i : ItemId) : Nat := spansCount s.tstack i + chCount s.items i

theorem spansCount_nil (i : ItemId) : spansCount [] i = 0 := rfl
theorem spansCount_cons (t : TEntry) (ts : List TEntry) (i : ItemId) :
    spansCount (t :: ts) i = (t.spans.1 ++ t.spans.2).count i + spansCount ts i := by
  simp [spansCount]
theorem spansCount_tail_le (ts : List TEntry) (i : ItemId) : spansCount ts.tail i ≤ spansCount ts i := by
  cases ts with
  | nil => exact Nat.le_refl _
  | cons t ts => rw [List.tail_cons, spansCount_cons]; omega
theorem spansCount_mergeTop_le (ts : List TEntry) (i : ItemId) :
    spansCount (WalkM.mergeTop ts) i ≤ spansCount ts i := by
  match ts with
  | b :: a :: rest =>
    simp only [WalkM.mergeTop, spansCount_cons, List.count_append]; omega
  | [] => exact Nat.le_refl _
  | [_] => simp [WalkM.mergeTop, spansCount_nil]

theorem count_setSides (b : Bool) (l₁ l₂ : List ItemId) (i : ItemId) :
    ((setSides b l₁ l₂).1 ++ (setSides b l₁ l₂).2).count i = (l₁ ++ l₂).count i := by
  cases b <;> simp [setSides, List.count_append]; omega
theorem count_getSide_le (p : List ItemId × List ItemId) (b : Bool) (i : ItemId) :
    (getSide p b).count i ≤ (p.1 ++ p.2).count i := by
  cases b <;> simp [getSide, List.count_append]

theorem default_spans : (default : TEntry).spans = ([], []) := rfl
theorem head!_cons (t : TEntry) (ts : List TEntry) : (t :: ts).head! = t := rfl
theorem spansCount_head_tail (ts : List TEntry) (i : ItemId) :
    spansCount ts i = (ts.head!.spans.1 ++ ts.head!.spans.2).count i + spansCount ts.tail i := by
  cases ts with
  | nil => rfl
  | cons t ts => simp [spansCount_cons]

namespace Items

variable (items : Items)

theorem count_le_chCount {j : Nat} (hj : j < items.size) (i : ItemId) :
    (items.ch j).count i ≤ chCount items i :=
  Finset.single_le_sum (f := fun j => (items.ch j).count i) (fun _ _ => Nat.zero_le _)
    (Finset.mem_range.mpr hj)

theorem chCount_push (x : Item) (hx : x.ch = []) (i : ItemId) :
    chCount (items.push x) i = chCount items i := by
  simp only [chCount, Array.size_push, Finset.sum_range_succ, ch_push, ite_true, hx, List.count_nil,
    Nat.add_zero]
  exact Finset.sum_congr rfl fun j hj => by
    simp [Nat.ne_of_lt (Finset.mem_range.mp hj)]

theorem chCount_modify {j : Nat} (hj : j < items.size) (f : Item → Item) (i : ItemId) :
    chCount (items.modify j f) i + (items.ch j).count i =
      chCount items i + ((f items[j]).ch).count i := by
  have hm : j ∈ Finset.range (items.modify j f).size := by rw [Array.size_modify]; exact Finset.mem_range.mpr hj
  have hm' : j ∈ Finset.range items.size := Finset.mem_range.mpr hj
  unfold chCount
  rw [← Finset.sum_erase_add _ _ hm, ← Finset.sum_erase_add _ _ hm', ch_modify_self _ _ _ hj,
    Array.size_modify]
  have : ∑ k ∈ (Finset.range items.size).erase j, (Items.ch (items.modify j f) k).count i =
      ∑ k ∈ (Finset.range items.size).erase j, (items.ch k).count i :=
    Finset.sum_congr rfl fun k hk => by rw [ch_modify_ne items j k f (Finset.mem_erase.mp hk).1.symm]
  rw [this]; omega

theorem chCount_modify_of_ch (j : Nat) (f : Item → Item) (hf : ∀ it, (f it).ch = it.ch) (i : ItemId) :
    chCount (items.modify j f) i = chCount items i := by
  unfold chCount
  rw [Array.size_modify]
  exact Finset.sum_congr rfl fun k _ => by rw [ch_modify_of_ch _ _ _ _ hf]

theorem chCount_modify_of_le {j : Nat} (hj : items.size ≤ j) (f : Item → Item) (i : ItemId) :
    chCount (items.modify j f) i = chCount items i := by
  unfold chCount
  rw [Array.size_modify]
  exact Finset.sum_congr rfl fun k hk => by
    rw [ch_modify_ne _ _ _ _ (Nat.ne_of_gt (Nat.lt_of_lt_of_le (Finset.mem_range.mp hk) hj))]

theorem chCount_pos {i : ItemId} (h : 0 < chCount items i) : ∃ j, j < items.size ∧ i ∈ items.ch j := by
  by_contra hc
  simp only [not_exists, not_and] at hc
  refine Nat.lt_irrefl 0 (Nat.lt_of_lt_of_eq h (Finset.sum_eq_zero fun j hj => ?_))
  exact List.count_eq_zero.mpr (hc j (Finset.mem_range.mp hj))

theorem two_le_chCount {i j j' : ItemId} (hj : j < items.size) (hj' : j' < items.size) (hne : j ≠ j')
    (hi : i ∈ items.ch j) (hi' : i ∈ items.ch j') : 2 ≤ chCount items i := by
  have hsub : ({j, j'} : Finset Nat) ⊆ Finset.range items.size := by
    intro k hk
    simp only [Finset.mem_insert, Finset.mem_singleton] at hk
    rcases hk with rfl | rfl <;> exact Finset.mem_range.mpr ‹_›
  have key : (items.ch j).count i + (items.ch j').count i ≤ chCount items i := by
    have := Finset.sum_le_sum_of_subset (f := fun j => (items.ch j).count i) hsub
    rw [Finset.sum_pair hne] at this
    exact this
  have h1 := List.count_pos_iff.mpr hi
  have h2 := List.count_pos_iff.mpr hi'
  omega

theorem ch_nodup_of_chCount (h : ∀ i, chCount items i ≤ 1) {j : Nat} (hj : j < items.size) :
    (items.ch j).Nodup := by
  rw [List.nodup_iff_count_le_one]
  exact fun a => Nat.le_trans (count_le_chCount items hj a) (h a)

end Items

namespace WalkState

/-- `s` with the child list of `x` cleared (what `finishTstackTop x` is about to overwrite). -/
def dropCh (s : WalkState) (x : ItemId) : WalkState :=
  { s with items := s.items.modify x fun it => { it with ch := [] } }

theorem dropCh_tstack (s : WalkState) (x : ItemId) (ts : List TEntry) :
    ({ s with tstack := ts }).dropCh x = { s.dropCh x with tstack := ts } := rfl

/-- Span discipline: every id is in at most one place (a tstack span or a `ch` list), everything
placed is allocated, the root is never placed, a fixed item (`vertItem`/`edgeItem`) is placed only
if `P` says it has been pushed / finished, and the "loose" nodes `X` (allocated, not yet placed) are
nowhere. -/
structure Place (g : Graph) (P X : ItemId → Prop) (s : WalkState) : Prop where
  g_eq : s.g = g
  size : 1 + g.nv + g.ne ≤ s.items.size
  root_type : Items.type s.items rootItem = .F
  vert : ∀ v, v < g.nv → Items.type s.items (vertItem v) = .V
  edge : ∀ e, e < g.ne → Items.type s.items (edgeItem g e) = .Q
  le : ∀ i, s.cnt i ≤ 1
  alloc : ∀ i, 0 < s.cnt i → i < s.items.size
  root : s.cnt rootItem = 0
  fixed : ∀ i : Nat, 0 < i → i < 1 + g.nv + g.ne → 0 < s.cnt i → P i
  loose : ∀ i, X i → s.cnt i = 0
  loose_node : ∀ i, X i → 1 + g.nv + g.ne ≤ i ∧ i < s.items.size

namespace Place

variable {g : Graph} {P P' X X' : ItemId → Prop} {s s' : WalkState}

theorem mono (h : s.Place g P X) (hP : ∀ i, P i → P' i) (hX : ∀ i, X' i → X i) : s.Place g P' X' :=
  { h with
    fixed := fun i h0 hi hc => hP i (h.fixed i h0 hi hc)
    loose := fun i hi => h.loose i (hX i hi)
    loose_node := fun i hi => h.loose_node i (hX i hi) }

theorem edgeItem_g (h : s.Place g P X) (e : Nat) : edgeItem s.g e = edgeItem g e := by rw [h.g_eq]

theorem cnt_eq_zero (h : s.Place g P X) {i : Nat} (h0 : 0 < i) (hi : i < 1 + g.nv + g.ne) (hP : ¬ P i) :
    s.cnt i = 0 :=
  Nat.eq_zero_of_not_pos fun hc => hP (h.fixed i h0 hi hc)

theorem cnt_size (h : s.Place g P X) : s.cnt s.items.size = 0 :=
  Nat.eq_zero_of_not_pos fun hc => Nat.lt_irrefl _ (h.alloc _ hc)

theorem node_ge (h : s.Place g P X) {x : Nat} (hty : Items.type s.items x ∈ [NodeType.S, .P, .R]) :
    1 + g.nv + g.ne ≤ x := by
  by_contra hlt
  rcases Nat.lt_or_ge x (1 + g.nv) with hx | hx
  · rcases Nat.eq_zero_or_pos x with rfl | hx0
    · have hr : Items.type s.items 0 = .F := h.root_type
      rw [hr] at hty; simp at hty
    · have := h.vert (x - 1) (by omega)
      rw [show vertItem (x - 1) = x from Nat.add_sub_of_le hx0] at this
      rw [this] at hty; simp at hty
  · have := h.edge (x - (1 + g.nv)) (by omega)
    rw [show edgeItem g (x - (1 + g.nv)) = x from Nat.add_sub_of_le hx] at this
    rw [this] at hty; simp at hty

/-- Generic step: the counts of `s'` are bounded by those of `s`, plus one for the so far unplaced
`x`, which `s'` places. -/
theorem step (h : s.Place g P X) {x : Nat} (hg : s'.g = s.g) (hsz : s.items.size ≤ s'.items.size)
    (hty : ∀ j, j < s.items.size → Items.type s'.items j = Items.type s.items j)
    (hx : s.cnt x = 0) (hxlt : x < s'.items.size) (hx0 : x ≠ 0)
    (hle : ∀ i, s'.cnt i ≤ s.cnt i + if i = x then 1 else 0) :
    s'.Place g (fun i => P i ∨ i = x) (fun i => X i ∧ i ≠ x) where
  g_eq := hg.trans h.g_eq
  size := Nat.le_trans h.size hsz
  root_type := by rw [hty rootItem (by have := h.size; show 0 < _; omega)]; exact h.root_type
  vert v hv := by rw [hty _ (by have := h.size; show 1 + v < _; omega)]; exact h.vert v hv
  edge e he := by rw [hty _ (by have := h.size; show 1 + g.nv + e < _; omega)]; exact h.edge e he
  le i := by
    by_cases hix : i = x
    · subst hix; have h1 := hle i; simp only [ite_true] at h1; omega
    · have h1 := hle i; have h2 := h.le i; simp only [hix, ite_false] at h1; omega
  alloc i hi := by
    by_cases hix : i = x
    · exact hix ▸ hxlt
    · have h1 := hle i; simp only [hix, ite_false] at h1
      exact Nat.lt_of_lt_of_le (h.alloc i (by omega)) hsz
  root := by
    have h1 := hle rootItem; have h2 := h.root
    rw [h2] at h1; split_ifs at h1 with h'
    · exact absurd h'.symm hx0
    · omega
  fixed i h0 hi hc := by
    by_cases hix : i = x
    · exact Or.inr hix
    · have h1 := hle i; simp only [hix, ite_false] at h1
      exact Or.inl (h.fixed i h0 hi (by omega))
  loose i hi := by
    have h1 := hle i; have h2 := h.loose i hi.1
    simp only [hi.2, ite_false] at h1; omega
  loose_node i hi := ⟨(h.loose_node i hi.1).1, Nat.lt_of_lt_of_le (h.loose_node i hi.1).2 hsz⟩

/-- A step that places nothing new. -/
theorem of_le (h : s.Place g P X) (hg : s'.g = s.g) (hsz : s.items.size ≤ s'.items.size)
    (hty : ∀ j, j < s.items.size → Items.type s'.items j = Items.type s.items j)
    (hle : ∀ i, s'.cnt i ≤ s.cnt i) : s'.Place g P X :=
  { h with
    g_eq := hg.trans h.g_eq
    size := Nat.le_trans h.size hsz
    root_type := by rw [hty rootItem (by have := h.size; show 0 < _; omega)]; exact h.root_type
    vert := fun v hv => by rw [hty _ (by have := h.size; show 1 + v < _; omega)]; exact h.vert v hv
    edge := fun e he => by rw [hty _ (by have := h.size; show 1 + g.nv + e < _; omega)]; exact h.edge e he
    le := fun i => Nat.le_trans (hle i) (h.le i)
    alloc := fun i hi => Nat.lt_of_lt_of_le (h.alloc i (Nat.lt_of_lt_of_le hi (hle i))) hsz
    root := Nat.eq_zero_of_le_zero (h.root ▸ hle rootItem)
    fixed := fun i h0 hi hc => h.fixed i h0 hi (Nat.lt_of_lt_of_le hc (hle i))
    loose := fun i hi => Nat.eq_zero_of_le_zero (h.loose i hi ▸ hle i)
    loose_node := fun i hi => ⟨(h.loose_node i hi).1, Nat.lt_of_lt_of_le (h.loose_node i hi).2 hsz⟩ }

/-- Placing a node `x` (loose) strips it from `X` and leaves `P` alone. -/
theorem node_of {x : Nat} (h : s.Place g (fun i => P i ∨ i = x) X') (hx : 1 + g.nv + g.ne ≤ x) :
    s.Place g P X' :=
  { h with fixed := fun i h0 hi hc =>
      (h.fixed i h0 hi hc).resolve_right fun h' => by have : i = x := h'; omega }

/-- Placing a fixed item `x` leaves `X` alone. -/
theorem fixed_of {x : Nat} (h : s.Place g P' (fun i => X i ∧ i ≠ x)) (hx : ¬ X x) : s.Place g P' X :=
  h.mono (fun _ hi => hi) fun _ hi => ⟨hi, fun h' => hx (h' ▸ hi)⟩

theorem addLoose (h : s.Place g P X) {x : Nat} (hx : s.cnt x = 0) (hn : 1 + g.nv + g.ne ≤ x)
    (hlt : x < s.items.size) : s.Place g P (fun i => X i ∨ i = x) :=
  { h with
    loose := fun i hi => hi.elim (h.loose i) fun h' => h' ▸ hx
    loose_node := fun i hi => hi.elim (h.loose_node i) fun h' => h' ▸ ⟨hn, hlt⟩ }

theorem of_eq (h : s.Place g P X) (hg : s'.g = s.g) (hi : s'.items = s.items) (ht : s'.tstack = s.tstack) :
    s'.Place g P X :=
  h.of_le hg (hi ▸ Nat.le_refl _) (fun _ _ => by rw [hi]) fun _ => by simp [cnt, hi, ht]

theorem tstack_le (h : s.Place g P X) {ts : List TEntry}
    (hts : ∀ i, spansCount ts i ≤ spansCount s.tstack i) : ({ s with tstack := ts }).Place g P X :=
  h.of_le rfl (Nat.le_refl _) (fun _ _ => rfl) fun i => by simp only [cnt]; have := hts i; omega

theorem tail (h : s.Place g P X) : ({ s with tstack := s.tstack.tail }).Place g P X :=
  h.tstack_le fun _ => spansCount_tail_le _ _
theorem mergeTop (h : s.Place g P X) : ({ s with tstack := WalkM.mergeTop s.tstack }).Place g P X :=
  h.tstack_le fun _ => spansCount_mergeTop_le _ _
theorem modifyCur_reorder (h : s.Place g P X) (b : Bool) :
    ({ s with tstack := match s.tstack with
      | a :: rest => { a with spans := setSides b (a.spans.1 ++ a.spans.2) [] } :: rest
      | [] => [] }).Place g P X := by
  refine h.tstack_le fun i => ?_
  cases s.tstack with
  | nil => exact Nat.le_refl _
  | cons a rest => simp only [spansCount_cons, count_setSides, List.append_nil]; exact Nat.le_refl _

theorem set_stackDir (h : s.Place g P X) (a : Array Bool) : ({ s with stackDir := a }).Place g P X :=
  h.of_eq rfl rfl rfl
theorem set_stackVerts (h : s.Place g P X) (a : Array Nat) : ({ s with stackVerts := a }).Place g P X :=
  h.of_eq rfl rfl rfl

theorem push (h : s.Place g P X) (ty : NodeType) :
    ({ s with items := s.items.push ⟨ty, (none, none), []⟩ }).Place g P
      (fun i => X i ∨ i = s.items.size) := by
  have h' : ({ s with items := s.items.push ⟨ty, (none, none), []⟩ }).Place g P X :=
    h.of_le rfl (by simp) (fun j hj => by rw [Items.type_push]; simp [Nat.ne_of_lt hj])
      fun i => by simp only [cnt]; rw [Items.chCount_push _ _ rfl]
  refine h'.addLoose ?_ h.size (by simp)
  show spansCount s.tstack s.items.size + chCount (s.items.push _) s.items.size = 0
  rw [Items.chCount_push _ _ rfl]; exact h.cnt_size

theorem modify_vs (h : s.Place g P X) (x : ItemId) (vs : Option Nat × Option Nat) :
    ({ s with items := s.items.modify x fun it => { it with vs := vs } }).Place g P X :=
  h.of_le rfl (by simp) (fun j _ => Items.type_modify _ _ _ _ fun _ => rfl)
    fun i => by
      simp only [cnt]
      rw [Items.chCount_modify_of_ch _ x (fun it => { it with vs := vs }) (fun _ => rfl)]

theorem dropCh (h : s.Place g P X) (x : ItemId) : (s.dropCh x).Place g P X := by
  refine h.of_le rfl (by simp [WalkState.dropCh]) (fun j _ => Items.type_modify _ _ _ _ fun _ => rfl)
    fun i => ?_
  simp only [cnt, WalkState.dropCh]
  rcases Nat.lt_or_ge x s.items.size with hx | hx
  · have := Items.chCount_modify s.items hx (fun it => { it with ch := [] }) i
    simp only [List.count_nil, Nat.add_zero] at this; omega
  · rw [Items.chCount_modify_of_le _ hx]

theorem cons (h : s.Place g P X) {x : Nat} (hx : s.cnt x = 0) (hxlt : x < s.items.size) (hx0 : x ≠ 0)
    (vS tD fI : Nat) (dir : Bool) :
    ({ s with tstack := ⟨vS, tD, fI, setSides dir [x] []⟩ :: s.tstack }).Place g (fun i => P i ∨ i = x)
      (fun i => X i ∧ i ≠ x) :=
  h.step rfl (Nat.le_refl _) (fun _ _ => rfl) hx hxlt hx0 fun i => by
    simp only [cnt, spansCount_cons, count_setSides, List.append_nil]
    by_cases hix : i = x
    · subst hix; simp only [List.count_singleton, beq_self_eq_true, ite_true]; omega
    · simp [Ne.symm hix, hix]

/-- Pushing a fixed item `x` (its vertex / edge is now "done"). -/
theorem cons_fixed (h : s.Place g P X) {x : Nat} (hx0 : 0 < x) (hxlt : x < 1 + g.nv + g.ne) (hP : ¬ P x)
    (vS tD fI : Nat) (dir : Bool) :
    ({ s with tstack := ⟨vS, tD, fI, setSides dir [x] []⟩ :: s.tstack }).Place g (fun i => P i ∨ i = x) X :=
  (h.cons (h.cnt_eq_zero hx0 hxlt hP) (Nat.lt_of_lt_of_le hxlt h.size) (Nat.ne_of_gt hx0) vS tD fI dir).fixed_of
    fun hx => by have := (h.loose_node x hx).1; omega

/-- `finishTstackTop x`: `x`'s children become one side of the top ear, and `x` replaces the ear. -/
theorem finishTop {x : Nat} (h : (s.dropCh x).Place g P X) (hx : X x) (dir : Bool) (vs : Option Nat × Option Nat) :
    ({ s with
        items := s.items.modify x fun it => { it with vs := vs, ch := getSide s.tstack.head!.spans dir },
        tstack := match s.tstack with
          | a :: rest => { a with spans := setSides dir [x] [] } :: rest
          | [] => [] }).Place g P (fun i => X i ∧ i ≠ x) := by
  have hn := h.loose_node x hx
  have hlt : x < s.items.size := by simpa [WalkState.dropCh] using hn.2
  refine (h.step (x := x) (by rfl) (by simp [WalkState.dropCh]) ?_ (h.loose x hx) (by simpa [WalkState.dropCh] using hlt)
    (Nat.ne_of_gt (Nat.lt_of_lt_of_le (by omega) hn.1)) fun i => ?_).node_of hn.1
  · intro j _
    simp only [WalkState.dropCh]
    rw [Items.type_modify (items := s.items) x j (fun it => { it with vs := vs, ch := getSide s.tstack.head!.spans dir })
      fun _ => rfl, Items.type_modify (items := s.items) x j (fun it => { it with ch := [] }) fun _ => rfl]
  · simp only [cnt, WalkState.dropCh]
    have h1 := Items.chCount_modify s.items hlt
      (fun it => { it with vs := vs, ch := getSide s.tstack.head!.spans dir }) i
    have h2 := Items.chCount_modify s.items hlt (fun it => { it with ch := [] }) i
    simp only [List.count_nil, Nat.add_zero] at h1 h2
    have h3 := count_getSide_le s.tstack.head!.spans dir i
    cases hts : s.tstack with
    | nil =>
      have h4 : getSide ([] : List TEntry).head!.spans dir = [] := by cases dir <;> rfl
      simp only [hts, spansCount_nil, h4, List.count_nil] at h1 ⊢
      split_ifs <;> omega
    | cons a rest =>
      simp only [hts, head!_cons, spansCount_cons, count_setSides, List.append_nil] at h1 h3 ⊢
      by_cases hix : i = x
      · subst hix; simp only [List.count_singleton, beq_self_eq_true, ite_true]; omega
      · simp only [List.count_singleton, hix, Ne.symm hix, beq_iff_eq, ite_false]; omega

/-- `maybeUnwrapNxt` reusing the node `x` at the head of a side of `nxt`: its children move into
that side and `x` becomes loose. -/
theorem unwrap (h : s.Place g P X) {a t : TEntry} {rest : List TEntry} (hts : s.tstack = a :: t :: rest)
    {x : Nat} {dir : Bool} (hx : x ∈ getSide t.spans dir) (hn : 1 + g.nv + g.ne ≤ x)
    (hlt : x < s.items.size) :
    (({ s with tstack := a :: { t with spans := setSides dir (Items.ch s.items x) [] } :: rest }).dropCh x).Place
      g P (fun i => X i ∨ i = x) := by
  have hmem : x ∈ t.spans.1 ++ t.spans.2 := by
    cases dir <;> simp [getSide] at hx <;> simp [hx]
  have hpos := List.count_pos_iff.mpr hmem
  have key : ∀ i, (({ s with tstack := a :: { t with spans := setSides dir (Items.ch s.items x) [] } :: rest }).dropCh x).cnt i
      + (t.spans.1 ++ t.spans.2).count i = s.cnt i := by
    intro i
    simp only [cnt, WalkState.dropCh, hts, spansCount_cons, count_setSides, List.append_nil]
    have h2 := Items.chCount_modify s.items hlt (fun it => { it with ch := [] }) i
    simp only [List.count_nil, Nat.add_zero] at h2
    omega
  refine Place.addLoose ?_ ?_ hn (by simpa [WalkState.dropCh] using hlt)
  · exact h.of_le rfl (by simp [WalkState.dropCh]) (fun j _ => Items.type_modify _ _ _ _ fun _ => rfl)
      fun i => by have := key i; omega
  · have := key x; have := h.le x; omega

theorem ch_eq_getElem (items : Items) {j : Nat} (hj : j < items.size) : Items.ch items j = items[j].ch := by
  simp [Items.ch, Array.getElem?_eq_getElem hj]

theorem dropCh_tstack_le (h : (s.dropCh x).Place g P X) {ts : List TEntry}
    (hts : ∀ i, spansCount ts i ≤ spansCount s.tstack i) : (({ s with tstack := ts }).dropCh x).Place g P X := by
  rw [dropCh_tstack]; exact h.tstack_le hts
theorem dropCh_mergeTop (h : (s.dropCh x).Place g P X) :
    (({ s with tstack := WalkM.mergeTop s.tstack }).dropCh x).Place g P X :=
  h.dropCh_tstack_le fun _ => spansCount_mergeTop_le _ _
theorem spansCount_reorder_le (ts : List TEntry) (b : Bool) (i : ItemId) :
    spansCount (match ts with
      | a :: rest => { a with spans := setSides b (a.spans.1 ++ a.spans.2) [] } :: rest
      | [] => []) i ≤ spansCount ts i := by
  cases ts with
  | nil => exact Nat.le_refl _
  | cons a rest => simp only [spansCount_cons, count_setSides, List.append_nil]; exact Nat.le_refl _
theorem spansCount_modifyCur_le (f : TEntry → TEntry)
    (hf : ∀ t i, ((f t).spans.1 ++ (f t).spans.2).count i ≤ (t.spans.1 ++ t.spans.2).count i)
    (ts : List TEntry) (i : ItemId) :
    spansCount (match ts with | a :: rest => f a :: rest | [] => []) i ≤ spansCount ts i := by
  cases ts with
  | nil => exact Nat.le_refl _
  | cons a rest => simp only [spansCount_cons]; have := hf a i; omega
theorem modifyCur (h : s.Place g P X) (f : TEntry → TEntry)
    (hf : ∀ t i, ((f t).spans.1 ++ (f t).spans.2).count i ≤ (t.spans.1 ++ t.spans.2).count i) :
    ({ s with tstack := match s.tstack with | a :: rest => f a :: rest | [] => [] }).Place g P X :=
  h.tstack_le fun _ => spansCount_modifyCur_le f hf _ _
theorem dropCh_modifyCur (h : (s.dropCh x).Place g P X) (f : TEntry → TEntry)
    (hf : ∀ t i, ((f t).spans.1 ++ (f t).spans.2).count i ≤ (t.spans.1 ++ t.spans.2).count i) :
    (({ s with tstack := match s.tstack with | a :: rest => f a :: rest | [] => [] }).dropCh x).Place g P X :=
  h.dropCh_tstack_le fun _ => spansCount_modifyCur_le f hf _ _
theorem dropCh_reorder (h : (s.dropCh x).Place g P X) (b : Bool) :
    (({ s with tstack := match s.tstack with
      | a :: rest => { a with spans := setSides b (a.spans.1 ++ a.spans.2) [] } :: rest
      | [] => [] }).dropCh x).Place g P X :=
  h.dropCh_tstack_le fun _ => spansCount_reorder_le _ _ _

/-- Pop the top ear and write `x :: top.spans.2` as the children of the allocated `q`. -/
theorem q_cons (h : s.Place g P X) {x q : Nat} (hx : X x) (hq : q < s.items.size) :
    ({ s with
        items := s.items.modify q fun it => { it with ch := x :: s.tstack.head!.spans.2 },
        tstack := s.tstack.tail }).Place g P (fun i => X i ∧ i ≠ x) := by
  have hn := h.loose_node x hx
  refine Place.node_of (x := x) ?_ hn.1
  refine h.step (by rfl) (by simp) (fun j _ => Items.type_modify _ _ _ _ fun _ => rfl) (h.loose x hx)
    (by simpa using hn.2) (Nat.ne_of_gt (Nat.lt_of_lt_of_le (by omega) hn.1)) fun i => ?_
  simp only [cnt]
  have h1 := Items.chCount_modify s.items hq (fun it => { it with ch := x :: s.tstack.head!.spans.2 }) i
  have h2 := spansCount_head_tail s.tstack i
  simp only [List.count_cons, List.count_append, beq_iff_eq] at h1 h2
  by_cases hix : i = x
  · subst hix; simp only [ite_true] at h1 ⊢; omega
  · simp only [Ne.symm hix, hix, ite_false] at h1 ⊢; omega

/-- Pop two ears and write `backedge.spans.1 ++ t.spans.2` as the children of `q`. -/
theorem q_merge (h : s.Place g P X) {q : Nat} (hq : q < s.items.size) :
    ({ s with
        items := s.items.modify q fun it => { it with ch := s.tstack.head!.spans.1 ++ s.tstack.tail.head!.spans.2 },
        tstack := s.tstack.tail.tail }).Place g P X := by
  refine h.of_le rfl (by simp) (fun j _ => Items.type_modify _ _ _ _ fun _ => rfl) fun i => ?_
  simp only [cnt]
  have h1 := Items.chCount_modify s.items hq
    (fun it => { it with ch := s.tstack.head!.spans.1 ++ s.tstack.tail.head!.spans.2 }) i
  have h2 := spansCount_head_tail s.tstack i
  have h3 := spansCount_head_tail s.tstack.tail i
  simp only [List.count_append] at h1 h2 h3
  omega

/-- Write `[x]` as the children of `q`. -/
theorem q_single (h : s.Place g P X) {x q : Nat} (hx : X x) (hq : q < s.items.size) :
    ({ s with items := s.items.modify q fun it => { it with ch := [x] } }).Place g P
      (fun i => X i ∧ i ≠ x) := by
  have hn := h.loose_node x hx
  refine Place.node_of (x := x) ?_ hn.1
  refine h.step (by rfl) (by simp) (fun j _ => Items.type_modify _ _ _ _ fun _ => rfl) (h.loose x hx)
    (by simpa using hn.2) (Nat.ne_of_gt (Nat.lt_of_lt_of_le (by omega) hn.1)) fun i => ?_
  simp only [cnt]
  have h1 := Items.chCount_modify s.items hq (fun it => { it with ch := [x] }) i
  simp only [List.count_singleton, beq_iff_eq] at h1
  by_cases hix : i = x
  · subst hix; simp only [ite_true] at h1 ⊢; omega
  · simp only [Ne.symm hix, hix, ite_false] at h1 ⊢; omega

/-- Append the fixed item `x` to the children of `j`. -/
theorem append_fixed (h : s.Place g P X) {j x : Nat} (hj : j < s.items.size) (hx0 : 0 < x)
    (hxlt : x < 1 + g.nv + g.ne) (hP : ¬ P x) :
    ({ s with items := s.items.modify j fun it => { it with ch := it.ch ++ [x] } }).Place g
      (fun i => P i ∨ i = x) X := by
  refine Place.fixed_of (x := x) ?_ fun hx => by have := (h.loose_node x hx).1; omega
  refine h.step (by rfl) (by simp) (fun j _ => Items.type_modify _ _ _ _ fun _ => rfl)
    (h.cnt_eq_zero hx0 hxlt hP) (by simpa using Nat.lt_of_lt_of_le hxlt h.size) (Nat.ne_of_gt hx0) fun i => ?_
  simp only [cnt]
  have h1 := Items.chCount_modify s.items hj (fun it => { it with ch := it.ch ++ [x] }) i
  simp only [List.count_append, List.count_singleton, beq_iff_eq, ← ch_eq_getElem s.items hj] at h1
  by_cases hix : i = x
  · subst hix; simp only [ite_true] at h1 ⊢; omega
  · simp only [Ne.symm hix, hix, ite_false] at h1 ⊢; omega

/-- Pop the top ear and append its second side to the root's children. -/
theorem root_append (h : s.Place g P X) :
    ({ s with
        items := s.items.modify rootItem fun it => { it with ch := it.ch ++ s.tstack.head!.spans.2 },
        tstack := s.tstack.tail }).Place g P X := by
  have h0 : rootItem < s.items.size := Nat.lt_of_lt_of_le (show (0 : Nat) < 1 + g.nv + g.ne by omega) h.size
  refine h.of_le rfl (by simp) (fun j _ => Items.type_modify _ _ _ _ fun _ => rfl) fun i => ?_
  simp only [cnt]
  have h1 := Items.chCount_modify s.items h0 (fun it => { it with ch := it.ch ++ s.tstack.head!.spans.2 }) i
  have h2 := spansCount_head_tail s.tstack i
  simp only [List.count_append, ← ch_eq_getElem s.items h0] at h1 h2
  omega

end Place

end WalkState

namespace WalkM

open WalkState

variable {g : Graph} {P X : ItemId → Prop} {s : WalkState}

theorem modifyItem_place {i : Nat} {f : Item → Item} {P' X' : ItemId → Prop}
    (hf : ({ s with items := s.items.modify i f }).Place g P' X') :
    wp (modifyItem i f) (fun _ s' => s'.Place g P' X') s := hf
theorem modify_place {f : WalkState → WalkState} {P' X' : ItemId → Prop} (hf : (f s).Place g P' X') :
    wp (modify f) (fun _ s' => s'.Place g P' X') s := hf
theorem allocItem_place (h : s.Place g P X) (ty : NodeType) :
    wp (allocItem ty) (fun item s' => s'.Place g P (fun i => X i ∨ i = item) ∧ ¬ X item) s :=
  ⟨h.push ty, fun hx => Nat.lt_irrefl _ (h.loose_node _ hx).2⟩
theorem makeVs_pure (a b : Nat) : wp (makeVs a b) (fun _ s' => s' = s) s := rfl
theorem popTstack_spec' :
    wp popTstack (fun t s' => t = s.tstack.head! ∧ s' = { s with tstack := s.tstack.tail }) s := ⟨rfl, rfl⟩

theorem loop_place {n : Nat} {cond : WalkM Bool} {body : WalkM Unit}
    (hc : ∀ s, s.Place g P X → wp cond (fun _ s' => s'.Place g P X) s)
    (hb : ∀ s, s.Place g P X → wp body (fun _ s' => s'.Place g P X) s) (h : s.Place g P X) :
    wp (loop n cond body) (fun _ s' => s'.Place g P X) s :=
  wp_loop (fun s => s.Place g P X) n cond body hc hb h fun _ h => h

theorem maybeUnwrapNxt_place (h : s.Place g P X) (ty : NodeType) (hty : ty ∈ [NodeType.S, .P, .R]) :
    wp (maybeUnwrapNxt ty)
      (fun item s' => (s'.dropCh item).Place g P (fun i => X i ∨ i = item) ∧ ¬ X item) s := by
  have halloc : wp (allocItem ty)
      (fun item s' => (s'.dropCh item).Place g P (fun i => X i ∨ i = item) ∧ ¬ X item) s :=
    ⟨(h.push ty).dropCh _, fun hx => Nat.lt_irrefl _ (h.loose_node _ hx).2⟩
  have hroot : Items.type s.items 0 ≠ ty := by
    have hr : Items.type s.items 0 = .F := h.root_type
    rw [hr]; rintro rfl; simp at hty
  unfold maybeUnwrapNxt
  simp only [wp_bind, wp_get, wp_ite, wp_pure, wp_nxt, wp_stackDir, wp_getItem, wp_modifyNxt]
  split
  · exact halloc
  · split
    · next heq =>
      rw [beq_iff_eq, Items.type_getElem!'] at heq
      rcases hts : s.tstack with _ | ⟨a, _ | ⟨t, rest⟩⟩
      · rw [hts] at heq
        exact absurd heq (by cases s.stackDir[([] : List TEntry).tail.head!.topDepth]! <;> exact hroot)
      · rw [hts] at heq
        exact absurd heq (by cases s.stackDir[[a].tail.head!.topDepth]! <;> exact hroot)
      · rw [hts] at heq
        simp only [List.tail_cons, head!_cons] at heq ⊢
        generalize s.stackDir[t.topDepth]! = dir at heq ⊢
        have hlt : (getSide t.spans dir).head! < s.items.size := by
          by_contra hge
          rw [Items.type_of_le _ _ (Nat.le_of_not_lt hge)] at heq
          rw [← heq] at hty; simp at hty
        have hn := h.node_ge (x := (getSide t.spans dir).head!) (heq ▸ hty)
        have hmem : (getSide t.spans dir).head! ∈ getSide t.spans dir := by
          cases hl : getSide t.spans dir with
          | nil =>
            rw [hl] at heq
            exact absurd heq hroot
          | cons y ys => exact List.mem_cons_self
        rw [Items.ch_getElem! _ _ hlt]
        exact ⟨h.unwrap hts hmem hn hlt, fun hx => by
          have := h.loose _ hx
          have h1 := List.count_pos_iff.mpr hmem
          have h2 := count_getSide_le t.spans dir (getSide t.spans dir).head!
          have h3 : spansCount s.tstack (getSide t.spans dir).head! ≤ s.cnt (getSide t.spans dir).head! :=
            Nat.le_add_right _ _
          rw [hts] at h3
          simp only [spansCount_cons] at h3
          omega⟩
    · exact halloc

theorem finishTstackTop_place {x : Nat} (h : (s.dropCh x).Place g P X) (hx : X x) :
    wp (finishTstackTop x) (fun _ s' => s'.Place g P (fun i => X i ∧ i ≠ x)) s := by
  unfold finishTstackTop
  simp only [wp_bind, wp_cur, wp_stackDir, wp_makeVs, wp_modifyItem, wp_modifyCur]
  exact h.finishTop hx _ _

theorem vertItem_ne_edgeItem {v : Nat} (hv : v < g.nv) (e : Nat) : vertItem v ≠ edgeItem g e := by
  show (1 + v : Nat) ≠ 1 + g.nv + e; omega

theorem finishEdge_place (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool)
    (h : s.Place g P X) (hv : curV < g.nv) (he : o.e < g.ne)
    (hPe : ¬ P (edgeItem g o.e)) (hPv : hasVert = false → ¬ P (vertItem curV)) :
    wp (finishEdge curV d o origTstack hasVert)
      (fun hv' s' => s'.Place g (fun i => P i ∨ i = edgeItem g o.e ∨ (hv' = true ∧ i = vertItem curV)) X ∧
        (hasVert = true → hv' = true)) s := by
  unfold finishEdge
  rw [bind_get]
  extract_lets g' nxtV e lowval isTree isType1 qItem jpB isT jp4 jpN jpS isF
  rw [bind_stackDir]
  extract_lets jp2 jpT
  have hq : qItem = edgeItem g o.e := congrArg (edgeItem · o.e) h.g_eq
  have hqlt : edgeItem g o.e < s.items.size := Nat.lt_of_lt_of_le (by show 1 + g.nv + o.e < _; omega) h.size
  have hqpos : 0 < edgeItem g o.e := by show 0 < 1 + g.nv + o.e; omega
  have hqlt' : edgeItem g o.e < 1 + g.nv + g.ne := by show 1 + g.nv + o.e < _; omega
  have hvpos : 0 < vertItem curV := by show 0 < 1 + curV; omega
  have hvlt : vertItem curV < 1 + g.nv + g.ne := by show 1 + curV < _; omega
  have hspr : ∀ b : Bool, (if b = true then NodeType.S else NodeType.R) ∈ [NodeType.S, .P, .R] := by
    intro b; cases b <;> simp
  -- Invariants before / after the edge item has been placed.
  let I₁ : WalkState → Prop := fun s' => s'.Place g (fun i => P i ∨ i = edgeItem g o.e) X
  let Post : Bool → WalkState → Prop := fun hv' s' =>
    s'.Place g (fun i => P i ∨ i = edgeItem g o.e ∨ (hv' = true ∧ i = vertItem curV)) X ∧
      (hasVert = true → hv' = true)
  have hI₁ : ∀ s', I₁ s' → Post hasVert s' := fun s' h' =>
    ⟨h'.mono (fun i hi => hi.elim Or.inl fun h => Or.inr (Or.inl h)) fun _ hi => hi, fun h => h⟩
  have hjpB : ∀ s', s'.Place g P X → wp (jpB ()) Post s' := by
    intro s' h'
    dsimp -zeta only [jpB]
    rw [hq]
    refine bind_spec (modifyItem_place (h'.append_fixed (x := edgeItem g o.e)
      (Nat.lt_of_lt_of_le (by show 1 + curV < _; omega) h'.size) hqpos hqlt' hPe)) fun _ s'' h'' => ?_
    exact hI₁ _ h''
  have hjpN : ∀ r b s', I₁ s' → wp (jpN r b) Post s' := by
    intro r b s' h'
    dsimp -zeta only [jpN]
    rw [bind_tstackSize, bind_nxt, bind_nxt]
    extract_lets jp5
    have hjp5 : ∀ s'', I₁ s'' → wp (jp5 ()) Post s'' := by
      intro s'' h''
      dsimp -zeta only [jp5]
      split
      · next hnv =>
        have hnv' : hasVert = false := by simpa using hnv
        have hP' : ¬ (P (vertItem curV) ∨ vertItem curV = edgeItem g o.e) :=
          fun h => h.elim (hPv hnv') (vertItem_ne_edgeItem hv o.e)
        have hc := h''.cons_fixed hvpos hvlt hP' curV d s''.nxtEdgeIdx s''.stackDir[d]!
        refine bind_spec (pushVertTstack_spec curV d) ?_
        rintro _ _ rfl
        have hpost : ∀ s₃, s₃.Place g (fun i => (P i ∨ i = edgeItem g o.e) ∨ i = vertItem curV) X →
            Post true s₃ := fun s₃ h₃ =>
          ⟨h₃.mono (fun i hi => hi.elim (fun h => h.elim Or.inl fun h => Or.inr (Or.inl h))
            fun h => Or.inr (Or.inr ⟨rfl, h⟩)) fun _ hi => hi, fun _ => rfl⟩
        split
        · refine bind_spec mergeTstackTops_spec ?_
          rintro _ _ rfl
          exact hpost _ hc.mergeTop
        · exact hpost _ hc
      · exact hI₁ _ h''
    split
    · refine bind_spec (maybeUnwrapNxt_place h' .P (by simp)) ?_
      rintro item s₁ ⟨h₁, hX₁⟩
      refine bind_spec mergeTstackTops_spec ?_
      rintro _ _ rfl
      refine bind_spec (finishTstackTop_place h₁.dropCh_mergeTop (Or.inr rfl)) fun _ s₂ h₂ => ?_
      exact hjp5 _ (h₂.mono (fun _ hi => hi) fun i hi => ⟨Or.inl hi, fun h => hX₁ (h ▸ hi)⟩)
    · exact hjp5 _ h'
  have hjpS : ∀ ty, ty ∈ [NodeType.S, .P, .R] → ∀ s', I₁ s' → wp (jpS ty) (fun _ => I₁) s' := by
    intro ty hty s' h'
    dsimp -zeta only [jpS]
    refine bind_spec (maybeUnwrapNxt_place h' ty hty) ?_
    rintro item s₁ ⟨h₁, hX₁⟩
    refine bind_spec mergeTstackTops_spec ?_
    rintro _ _ rfl
    exact wp_mono _ (finishTstackTop_place h₁.dropCh_mergeTop (Or.inr rfl)) fun _ s₂ h₂ =>
      h₂.mono (fun _ hi => hi) fun i hi => ⟨Or.inl hi, fun h => hX₁ (h ▸ hi)⟩
  have hjp2 : ∀ r b s', I₁ s' → wp (jp2 r b) Post s' := by
    intro r b s' h'
    dsimp -zeta only [jp2]
    extract_lets jp3
    have hjp3 : ∀ item s'', (match item with
        | some i => (s''.dropCh i).Place g (fun j => P j ∨ j = edgeItem g o.e) (fun j => X j ∨ j = i) ∧ ¬ X i
        | none => I₁ s'') →
        wp (jp3 item) Post s'' := by
      intro item s'' hi
      dsimp -zeta only [jp3]
      refine bind_spec mergeTstackTops_spec ?_
      rintro _ _ rfl
      refine bind_spec mergeTstackTops_spec ?_
      rintro _ _ rfl
      refine bind_spec (modifyCur_spec _) ?_
      rintro _ _ rfl
      cases item with
      | some i =>
        obtain ⟨hi, hX⟩ := hi
        refine bind_spec (finishTstackTop_place (hi.dropCh_mergeTop.dropCh_mergeTop.dropCh_modifyCur _
          fun t i => by simp only [count_setSides, List.append_nil]; exact Nat.le_refl _) (by exact Or.inr rfl))
          fun _ s₃ h₃ => ?_
        exact hjpN () _ _ (h₃.mono (fun _ hj => hj) fun j hj => ⟨Or.inl hj, fun h => hX (h ▸ hj)⟩)
      | none =>
        exact hjpN () _ _ (hi.mergeTop.mergeTop.modifyCur _ fun t i => by
          simp only [count_setSides, List.append_nil]; exact Nat.le_refl _)
    split
    · refine bind_spec (map_spec some (maybeUnwrapNxt_place h' _ (hspr b))) ?_
      rintro _ s₁ ⟨i, rfl, h₁⟩
      exact hjp3 (some i) s₁ h₁
    · exact hjp3 none s' h'
  have hjpT : ∀ r b s', I₁ s' → wp (jpT r b) Post s' := by
    intro r b s' h'
    dsimp -zeta only [jpT]
    split
    · split
      · rw [bind_tstackSize]
        refine bind_spec (loop_place (fun _ h => h)
          (fun s h => by rw [wp_mergeTstackTops]; exact h.mergeTop) h') fun _ s₁ h₁ => ?_
        exact hjp2 () _ _ h₁
      · exact hjp2 () _ _ h'
    · exact hjpN () _ _ h'
  split
  · rw [hq]
    refine bind_spec (modifyItem_place (h.modify_vs _ (some curV, none))) fun _ s₁ h₁ => ?_
    refine bind_spec (modify_place (h₁.of_eq rfl rfl rfl)) fun _ s₂ h₂ => ?_
    split
    · split
      · refine bind_spec (allocItem_place h₂ .I) ?_
        rintro item s₃ ⟨h₃, hX₃⟩
        refine bind_spec (makeVs_pure _ _) ?_
        rintro vs _ rfl
        refine bind_spec (modifyItem_place (h₃.modify_vs item vs)) fun _ s₄ h₄ => ?_
        refine bind_spec popTstack_spec' ?_
        rintro _ _ ⟨rfl, rfl⟩
        refine bind_spec (modifyItem_place (g := g) (P' := P) (X' := fun i => (X i ∨ i = item) ∧ i ≠ item) ?_)
          fun _ s₅ h₅ => ?_
        · exact h₄.q_cons (x := item) (Or.inr rfl) (Nat.lt_of_lt_of_le hqlt' h₄.size)
        · exact hjpB _ (h₅.mono (fun _ hi => hi) fun i hi => ⟨Or.inl hi, fun h => hX₃ (h ▸ hi)⟩)
      · refine bind_spec popTstack_spec' ?_
        rintro _ _ ⟨rfl, rfl⟩
        refine bind_spec popTstack_spec' ?_
        rintro _ _ ⟨rfl, rfl⟩
        refine bind_spec (modifyItem_place (g := g) (P' := P) (X' := X) ?_) fun _ s₃ h₃ => ?_
        · exact h₂.q_merge (Nat.lt_of_lt_of_le hqlt' h₂.size)
        · exact hjpB _ h₃
    · refine bind_spec (modify_place (h₂.of_eq rfl rfl rfl)) fun _ s₃ h₃ => ?_
      refine bind_spec (allocItem_place h₃ .O) ?_
      rintro item s₄ ⟨h₄, hX₄⟩
      refine bind_spec (modifyItem_place (h₄.modify_vs item (some curV, none))) fun _ s₅ h₅ => ?_
      refine bind_spec (modifyItem_place (h₅.q_single (x := item) (Or.inr rfl)
        (Nat.lt_of_lt_of_le hqlt' h₅.size))) fun _ s₆ h₆ => ?_
      exact hjpB _ (h₆.mono (fun _ hi => hi) fun i hi => ⟨Or.inl hi, fun h => hX₄ (h ▸ hi)⟩)
  · refine bind_spec (makeVs_pure _ _) ?_
    rintro vs _ rfl
    rw [hq]
    refine bind_spec (modifyItem_place (h.modify_vs _ vs)) fun _ s₁ h₁ => ?_
    split
    · refine bind_spec (pushEdgeTstack_spec nxtV d e) ?_
      rintro _ _ rfl
      rw [h₁.edgeItem_g]
      have h₁' : I₁ _ := h₁.cons_fixed hqpos hqlt' hPe nxtV d s₁.nxtEdgeIdx s₁.stackDir[d]!
      rw [bind_tstackSize]
      refine bind_spec (loop_place (fun _ h => h) ?_ h₁') fun _ s₂ h₂ => ?_
      · intro s₂ h₂
        rw [bind_nxt]
        split
        · rw [bind_nxt]
          refine bind_spec (setStackDir_spec _ _) ?_
          rintro _ _ rfl
          refine bind_spec mergeTstackTops_spec ?_
          rintro _ _ rfl
          exact hjpS .S (by simp) _ (h₂.set_stackDir _).mergeTop
        · rw [bind_nxt, bind_cur]
          split
          · exact hjpS .P (by simp) _ h₂
          · exact hjpS .R (by simp) _ h₂
      rw [bind_get]
      extract_lets fo
      rw [bind_cur]
      split
      · rw [bind_tstackSize]
        refine bind_spec (loop_place (fun _ h => h)
          (fun s h => by rw [wp_mergeTstackTops]; exact h.mergeTop) h₂) fun _ s₃ h₃ => ?_
        exact hjpT () _ _ h₃
      · exact hjpT () _ _ h₂
    · refine bind_spec (pushEdgeTstack_spec curV lowval e) ?_
      rintro _ _ rfl
      rw [h₁.edgeItem_g]
      have h₁' : I₁ _ := h₁.cons_fixed hqpos hqlt' hPe curV lowval s₁.nxtEdgeIdx s₁.stackDir[lowval]!
      refine bind_spec (modify_place (h₁'.of_eq rfl rfl rfl)) fun _ s₂ h₂ => ?_
      exact hjpN () _ _ h₂

theorem DfsOut.vertsList_cons (o : DfsOut) (rest : List DfsOut) :
    DfsOut.vertsList (o :: rest) = DfsOut.vertsList [o] ++ DfsOut.vertsList rest := by
  cases o <;> simp [DfsOut.vertsList]

/-- Items that may be placed after walking vertices `vs` and edges `es`. -/
def Pushed (g : Graph) (P : ItemId → Prop) (vs es : List Nat) (i : ItemId) : Prop :=
  P i ∨ (∃ v ∈ vs, i = vertItem v) ∨ ∃ e ∈ es, i = edgeItem g e

theorem vertItem_inj {v w : Nat} (h : vertItem v = vertItem w) : v = w := by
  have : (1 + v : Nat) = 1 + w := h; omega
theorem edgeItem_inj {e e' : Nat} (h : edgeItem g e = edgeItem g e') : e = e' := by
  have : (1 + g.nv + e : Nat) = 1 + g.nv + e' := h; omega
theorem edgeItem_ne_vertItem {v : Nat} (hv : v < g.nv) (e : Nat) : edgeItem g e ≠ vertItem v :=
  (vertItem_ne_edgeItem hv e).symm

theorem walk_place_aux (g : Graph) :
    (∀ (t : DfsTree) (d : Nat) (P X : ItemId → Prop) (s : WalkState), s.Place g P X →
      (∀ v ∈ t.verts, v < g.nv) → (∀ e ∈ t.edges, e < g.ne) → t.verts.Nodup → t.edges.Nodup →
      (∀ v ∈ t.verts, ¬ P (vertItem v)) → (∀ e ∈ t.edges, ¬ P (edgeItem g e)) →
      wp (walkTree t d) (fun _ s' => s'.Place g (Pushed g P t.verts t.edges) X) s) ∧
    (∀ (v d : Nat) (outs : List DfsOut) (hasVert : Bool) (P X : ItemId → Prop) (s : WalkState),
      s.Place g P X → v < g.nv →
      (∀ w ∈ DfsOut.vertsList outs, w < g.nv) → (∀ e ∈ DfsOut.edgesList outs, e < g.ne) →
      (DfsOut.vertsList outs).Nodup → (DfsOut.edgesList outs).Nodup → v ∉ DfsOut.vertsList outs →
      (∀ w ∈ DfsOut.vertsList outs, ¬ P (vertItem w)) → (∀ e ∈ DfsOut.edgesList outs, ¬ P (edgeItem g e)) →
      (hasVert = false → ¬ P (vertItem v)) →
      wp (walkOuts v d outs hasVert) (fun hv' s' =>
        s'.Place g (Pushed g (fun i => P i ∨ (hv' = true ∧ i = vertItem v))
          (DfsOut.vertsList outs) (DfsOut.edgesList outs)) X ∧ (hasVert = true → hv' = true)) s) ∧
    (∀ (v d : Nat) (o : DfsOut) (hasVert : Bool) (P X : ItemId → Prop) (s : WalkState),
      s.Place g P X → v < g.nv →
      (∀ w ∈ DfsOut.vertsList [o], w < g.nv) → (∀ e ∈ DfsOut.edgesList [o], e < g.ne) →
      (DfsOut.vertsList [o]).Nodup → (DfsOut.edgesList [o]).Nodup → v ∉ DfsOut.vertsList [o] →
      (∀ w ∈ DfsOut.vertsList [o], ¬ P (vertItem w)) → (∀ e ∈ DfsOut.edgesList [o], ¬ P (edgeItem g e)) →
      (hasVert = false → ¬ P (vertItem v)) →
      wp (walkOut v d o hasVert) (fun hv' s' =>
        s'.Place g (Pushed g (fun i => P i ∨ (hv' = true ∧ i = vertItem v))
          (DfsOut.vertsList [o]) (DfsOut.edgesList [o])) X ∧ (hasVert = true → hv' = true)) s) := by
  refine walkTree.mutual_induct _ _ _ ?_ ?_ ?_ ?_
  · intro d v outs ih P X s h hvlt helt hvn hen hPv hPe
    simp only [DfsTree.verts, DfsTree.edges] at *
    have hv : v < g.nv := hvlt v List.mem_cons_self
    obtain ⟨hvo, hvn⟩ := List.nodup_cons.1 hvn
    simp only [walkTree, wp_bind, wp_modify]
    refine wp_mono _ (ih P X _ (h.set_stackVerts _) hv (fun w hw => hvlt w (List.mem_cons_of_mem _ hw))
      helt hvn hen hvo (fun w hw => hPv w (List.mem_cons_of_mem _ hw)) hPe
      (fun _ => hPv v List.mem_cons_self)) fun hasVert s₁ ⟨h₁, _⟩ => ?_
    split
    · next hhv =>
      exact h₁.mono (fun i hi => by
        rcases hi with (hi | ⟨_, rfl⟩) | ⟨w, hw, rfl⟩ | ⟨e, he, rfl⟩
        · exact Or.inl hi
        · exact Or.inr (Or.inl ⟨v, List.mem_cons_self, rfl⟩)
        · exact Or.inr (Or.inl ⟨w, List.mem_cons_of_mem _ hw, rfl⟩)
        · exact Or.inr (Or.inr ⟨e, he, rfl⟩)) fun _ hi => hi
    · next hhv =>
      have hhv : hasVert = false := by simpa using hhv
      refine bind_spec (setStackDir_spec _ _) ?_
      rintro _ _ rfl
      refine wp_mono _ (pushVertTstack_spec v d) ?_
      rintro _ _ rfl
      have hP' : ¬ Pushed g (fun i => P i ∨ (hasVert = true ∧ i = vertItem v))
          (DfsOut.vertsList outs) (DfsOut.edgesList outs) (vertItem v) := by
        rintro ((hi | ⟨hb, _⟩) | ⟨w, hw, hwe⟩ | ⟨e, _, hev⟩)
        · exact hPv v List.mem_cons_self hi
        · simp [hhv] at hb
        · exact hvo (vertItem_inj hwe ▸ hw)
        · exact vertItem_ne_edgeItem hv e hev
      refine ((h₁.set_stackDir _).cons_fixed (by show 0 < 1 + v; omega) (by show 1 + v < _; omega) hP'
        v d _ _).mono (fun i hi => ?_) fun _ hi => hi
      rcases hi with ((hi | ⟨hb, _⟩) | ⟨w, hw, rfl⟩ | ⟨e, he, rfl⟩) | rfl
      · exact Or.inl hi
      · simp [hhv] at hb
      · exact Or.inr (Or.inl ⟨w, List.mem_cons_of_mem _ hw, rfl⟩)
      · exact Or.inr (Or.inr ⟨e, he, rfl⟩)
      · exact Or.inr (Or.inl ⟨v, List.mem_cons_self, rfl⟩)
  · intro o v d hasVert ih P X s h hv hvlt helt hvn hen hvo hPv hPe hPcur
    have hvpos : 0 < vertItem v := by show 0 < 1 + v; omega
    have hvlt' : vertItem v < 1 + g.nv + g.ne := by show 1 + v < _; omega
    cases o with
    | back e dest cls =>
      simp only [DfsOut.vertsList, DfsOut.edgesList, List.mem_singleton, List.not_mem_nil,
        forall_eq, false_implies, implies_true] at hvlt helt hvn hen hvo hPv hPe ⊢
      have hfin : ∀ (P₁ : ItemId → Prop) (hasVert₁ : Bool) (origTstack : Nat) (s₁ : WalkState),
          s₁.Place g P₁ X → ¬ P₁ (edgeItem g e) → (hasVert₁ = false → ¬ P₁ (vertItem v)) →
          (∀ i, P₁ i → P i ∨ (hasVert₁ = true ∧ i = vertItem v)) → (hasVert = true → hasVert₁ = true) →
          wp (finishEdge v d (.back e dest cls) origTstack hasVert₁)
            (fun hv' s' => s'.Place g (Pushed g (fun i => P i ∨ (hv' = true ∧ i = vertItem v))
              [] [e]) X ∧ (hasVert = true → hv' = true)) s₁ := by
        intro P₁ hasVert₁ origTstack s₁ h₁ hPe₁ hPv₁ hP₁ hhv₁
        refine wp_mono _ (finishEdge_place v d (.back e dest cls) origTstack hasVert₁ h₁ hv helt hPe₁ hPv₁)
          fun hv' s₂ ⟨h₂, hhv₂⟩ => ⟨h₂.mono (fun i hi => ?_) fun _ hi => hi, fun hh => hhv₂ (hhv₁ hh)⟩
        rcases hi with hi | rfl | ⟨rfl, rfl⟩
        · rcases hP₁ i hi with hi | ⟨hb, rfl⟩
          · exact Or.inl (Or.inl hi)
          · exact Or.inl (Or.inr ⟨hhv₂ hb, rfl⟩)
        · exact Or.inr (Or.inr ⟨e, List.mem_singleton_self e, rfl⟩)
        · exact Or.inl (Or.inr ⟨rfl, rfl⟩)
      unfold walkOut
      dsimp only
      rw [bind_stackDir]
      refine bind_spec (setStackDir_spec _ _) ?_
      rintro _ _ rfl
      split
      · next hc =>
        have hhv : hasVert = false := by revert hc; cases hasVert <;> simp
        refine bind_spec (pushVertTstack_spec v d) ?_
        rintro _ _ rfl
        rw [bind_pure, bind_tstackSize]
        refine hfin _ true _ _ ((h.set_stackDir _).cons_fixed hvpos hvlt' (hPcur hhv) v d _ _)
          (fun hh => hh.elim hPe (edgeItem_ne_vertItem hv e)) (fun hh => nomatch hh)
          (fun i hi => hi.elim Or.inl fun hi => Or.inr ⟨rfl, hi⟩) fun _ => rfl
      · rw [bind_pure, bind_tstackSize]
        exact hfin _ hasVert _ _ (h.set_stackDir _) hPe hPcur (fun _ hi => Or.inl hi) fun hh => hh
    | tree e cls child =>
      obtain ⟨cv, couts⟩ := child
      simp only [DfsOut.vertsList, DfsOut.edgesList, List.append_nil, DfsTree.verts, DfsTree.edges]
        at hvlt helt hvn hen hvo hPv hPe ih ⊢
      obtain ⟨he, helt⟩ := List.forall_mem_cons.1 helt
      obtain ⟨heo, hen⟩ := List.nodup_cons.1 hen
      obtain ⟨hPe, hPe'⟩ := List.forall_mem_cons.1 hPe
      have hfin : ∀ (P₁ : ItemId → Prop) (hasVert₁ : Bool) (origTstack : Nat) (s₁ : WalkState),
          s₁.Place g P₁ X → (∀ w ∈ cv :: DfsOut.vertsList couts, ¬ P₁ (vertItem w)) →
          (∀ e' ∈ DfsOut.edgesList couts, ¬ P₁ (edgeItem g e')) →
          ¬ P₁ (edgeItem g e) → (hasVert₁ = false → ¬ P₁ (vertItem v)) →
          (∀ i, P₁ i → P i ∨ (hasVert₁ = true ∧ i = vertItem v)) → (hasVert = true → hasVert₁ = true) →
          wp ((modify fun s => { s with firstOccurrence := s.firstOccurrence.set! d s.g.ne }) >>= fun _ =>
              walkTree (.node cv couts) (d + 1) >>= fun _ =>
              finishEdge v d (.tree e cls (.node cv couts)) origTstack hasVert₁)
            (fun hv' s' => s'.Place g (Pushed g (fun i => P i ∨ (hv' = true ∧ i = vertItem v))
              (cv :: DfsOut.vertsList couts) (e :: DfsOut.edgesList couts)) X ∧
              (hasVert = true → hv' = true)) s₁ := by
        intro P₁ hasVert₁ origTstack s₁ h₁ hPv₁ hPe₁ hPe₁' hPcur₁ hP₁ hhv₁
        refine bind_spec (modify_place (h₁.of_eq rfl rfl rfl)) fun _ s₂ h₂ => ?_
        refine bind_spec (ih P₁ X _ h₂ hvlt helt hvn hen hPv₁ hPe₁) fun _ s₃ h₃ => ?_
        have hv₂ : ∀ w ∈ cv :: DfsOut.vertsList couts, w < g.nv := hvlt
        have hPe₃ : ¬ Pushed g P₁ (cv :: DfsOut.vertsList couts) (DfsOut.edgesList couts) (edgeItem g e) := by
          rintro (hh | ⟨w, hw, hwe⟩ | ⟨e', he', hee⟩)
          · exact hPe₁' hh
          · exact edgeItem_ne_vertItem (hv₂ w hw) e hwe
          · exact heo (edgeItem_inj hee ▸ he')
        have hPv₃ : hasVert₁ = false →
            ¬ Pushed g P₁ (cv :: DfsOut.vertsList couts) (DfsOut.edgesList couts) (vertItem v) := by
          rintro hh (hh' | ⟨w, hw, hwe⟩ | ⟨e', _, hee⟩)
          · exact hPcur₁ hh hh'
          · exact hvo (vertItem_inj hwe ▸ hw)
          · exact vertItem_ne_edgeItem hv e' hee
        refine wp_mono _ (finishEdge_place v d (.tree e cls (.node cv couts)) origTstack hasVert₁ h₃ hv he hPe₃ hPv₃)
          fun hv' s₄ ⟨h₄, hhv₄⟩ => ⟨h₄.mono (fun i hi => ?_) fun _ hi => hi, fun hh => hhv₄ (hhv₁ hh)⟩
        rcases hi with (hi | ⟨w, hw, rfl⟩ | ⟨e', he', rfl⟩) | rfl | ⟨rfl, rfl⟩
        · rcases hP₁ i hi with hi | ⟨hb, rfl⟩
          · exact Or.inl (Or.inl hi)
          · exact Or.inl (Or.inr ⟨hhv₄ hb, rfl⟩)
        · exact Or.inr (Or.inl ⟨w, hw, rfl⟩)
        · exact Or.inr (Or.inr ⟨e', List.mem_cons_of_mem _ he', rfl⟩)
        · exact Or.inr (Or.inr ⟨e, List.mem_cons_self, rfl⟩)
        · exact Or.inl (Or.inr ⟨rfl, rfl⟩)
      unfold walkOut
      dsimp only
      rw [bind_stackDir]
      refine bind_spec (setStackDir_spec _ _) ?_
      rintro _ _ rfl
      split
      · next hc =>
        have hhv : hasVert = false := by revert hc; cases hasVert <;> simp
        refine bind_spec (pushVertTstack_spec v d) ?_
        rintro _ _ rfl
        rw [bind_pure, bind_tstackSize]
        refine hfin _ true _ _ ((h.set_stackDir _).cons_fixed hvpos hvlt' (hPcur hhv) v d _ _)
          (fun w hw hh => hh.elim (hPv w hw) fun hh => hvo (vertItem_inj hh ▸ hw))
          (fun e' he' hh => hh.elim (hPe' e' he') (edgeItem_ne_vertItem hv e'))
          (fun hh => hh.elim hPe (edgeItem_ne_vertItem hv e)) (fun hh => nomatch hh)
          (fun i hi => hi.elim Or.inl fun hi => Or.inr ⟨rfl, hi⟩) fun _ => rfl
      · rw [bind_pure, bind_tstackSize]
        exact hfin _ hasVert _ _ (h.set_stackDir _) hPv hPe' hPe hPcur (fun _ hi => Or.inl hi) fun hh => hh
  · intro v d hasVert P X s h hv _ _ _ _ _ _ _ hPcur
    simp only [walkOuts, wp_pure]
    exact ⟨h.mono (fun i hi => Or.inl (Or.inl hi)) fun _ hi => hi, fun hh => hh⟩
  · intro v d hasVert o rest ih₁ ih₂ P X s h hv hvlt helt hvn hen hvo hPv hPe hPcur
    rw [DfsOut.vertsList_cons] at hvlt hvn hvo hPv ⊢
    rw [DfsOut.edgesList_cons] at helt hen hPe ⊢
    obtain ⟨hvlt₁, hvlt₂⟩ := List.forall_mem_append.1 hvlt
    obtain ⟨helt₁, helt₂⟩ := List.forall_mem_append.1 helt
    obtain ⟨hPv₁, hPv₂⟩ := List.forall_mem_append.1 hPv
    obtain ⟨hPe₁, hPe₂⟩ := List.forall_mem_append.1 hPe
    rw [List.mem_append, not_or] at hvo
    obtain ⟨hvn₁, hvn₂, hvd⟩ := List.nodup_append.1 hvn
    obtain ⟨hen₁, hen₂, hed⟩ := List.nodup_append.1 hen
    simp only [walkOuts]
    refine bind_spec (ih₁ P X s h hv hvlt₁ helt₁ hvn₁ hen₁ hvo.1 hPv₁ hPe₁ hPcur) ?_
    rintro hv₁ s₁ ⟨h₁, hhv₁⟩
    refine wp_mono _ (ih₂ hv₁ _ X s₁ h₁ hv hvlt₂ helt₂ hvn₂ hen₂ hvo.2 ?_ ?_ ?_) ?_
    · rintro w hw ((hh | ⟨_, hwv⟩) | ⟨w', hw', hww⟩ | ⟨e', _, hwe⟩)
      · exact hPv₂ w hw hh
      · exact hvo.2 (vertItem_inj hwv ▸ hw)
      · exact hvd w' hw' w hw (vertItem_inj hww).symm
      · exact vertItem_ne_edgeItem (hvlt₂ w hw) e' hwe
    · rintro e' he' ((hh | ⟨_, hev⟩) | ⟨w', hw', hew⟩ | ⟨e'', he'', hee⟩)
      · exact hPe₂ e' he' hh
      · exact edgeItem_ne_vertItem hv e' hev
      · exact edgeItem_ne_vertItem (hvlt₁ w' hw') e' hew
      · exact hed e'' he'' e' he' (edgeItem_inj hee).symm
    · rintro hh ((hh' | ⟨hb, _⟩) | ⟨w', hw', hvw⟩ | ⟨e', _, hve⟩)
      · exact hPcur (by cases hasVert <;> simp_all) hh'
      · simp [hh] at hb
      · exact hvo.1 (vertItem_inj hvw ▸ hw')
      · exact vertItem_ne_edgeItem hv e' hve
    · rintro hv₂ s₂ ⟨h₂, hhv₂⟩
      refine ⟨h₂.mono (fun i hi => ?_) fun _ hi => hi, fun hh => hhv₂ (hhv₁ hh)⟩
      rcases hi with (((hi | ⟨hb, rfl⟩) | ⟨w, hw, rfl⟩ | ⟨e', he', rfl⟩) | ⟨hb, rfl⟩) | ⟨w, hw, rfl⟩ | ⟨e', he', rfl⟩
      · exact Or.inl (Or.inl hi)
      · exact Or.inl (Or.inr ⟨hhv₂ hb, rfl⟩)
      · exact Or.inr (Or.inl ⟨w, List.mem_append_left _ hw, rfl⟩)
      · exact Or.inr (Or.inr ⟨e', List.mem_append_left _ he', rfl⟩)
      · exact Or.inl (Or.inr ⟨hb, rfl⟩)
      · exact Or.inr (Or.inl ⟨w, List.mem_append_right _ hw, rfl⟩)
      · exact Or.inr (Or.inr ⟨e', List.mem_append_right _ he', rfl⟩)

theorem walkForest_place (forest : List DfsTree) (h : s.Place g P X)
    (hvlt : ∀ v ∈ forest.flatMap DfsTree.verts, v < g.nv)
    (helt : ∀ e ∈ forest.flatMap DfsTree.edges, e < g.ne)
    (hvn : (forest.flatMap DfsTree.verts).Nodup) (hen : (forest.flatMap DfsTree.edges).Nodup)
    (hPv : ∀ v ∈ forest.flatMap DfsTree.verts, ¬ P (vertItem v))
    (hPe : ∀ e ∈ forest.flatMap DfsTree.edges, ¬ P (edgeItem g e)) :
    wp (walkForest forest)
      (fun _ s' => s'.Place g (Pushed g P (forest.flatMap DfsTree.verts) (forest.flatMap DfsTree.edges)) X) s := by
  induction forest generalizing P s with
  | nil => exact h.mono (fun _ hi => Or.inl hi) fun _ hi => hi
  | cons t rest ih =>
    rw [List.flatMap_cons] at hvlt helt hvn hen hPv hPe ⊢
    obtain ⟨hvlt₁, hvlt₂⟩ := List.forall_mem_append.1 hvlt
    obtain ⟨helt₁, helt₂⟩ := List.forall_mem_append.1 helt
    obtain ⟨hPv₁, hPv₂⟩ := List.forall_mem_append.1 hPv
    obtain ⟨hPe₁, hPe₂⟩ := List.forall_mem_append.1 hPe
    obtain ⟨hvn₁, hvn₂, hvd⟩ := List.nodup_append.1 hvn
    obtain ⟨hen₁, hen₂, hed⟩ := List.nodup_append.1 hen
    show wp ((walkTree t 0 >>= fun _ => popTstack >>= fun top =>
      modifyItem rootItem fun it => { it with ch := it.ch ++ top.spans.2 }) >>= fun _ => walkForest rest) _ s
    refine bind_spec (P := fun _ s' => s'.Place g (Pushed g P t.verts t.edges) X) ?_ fun _ s₁ h₁ => ?_
    · refine bind_spec ((walk_place_aux g).1 t 0 P X s h hvlt₁ helt₁ hvn₁ hen₁ hPv₁ hPe₁) fun _ s₁ h₁ => ?_
      refine bind_spec popTstack_spec' ?_
      rintro _ _ ⟨rfl, rfl⟩
      exact h₁.root_append
    · refine wp_mono _ (ih h₁ hvlt₂ helt₂ hvn₂ hen₂ ?_ ?_) fun _ s₂ h₂ => h₂.mono (fun i hi => ?_) fun _ hi => hi
      · rintro v hv (hh | ⟨w, hw, hvw⟩ | ⟨e, _, hve⟩)
        · exact hPv₂ v hv hh
        · exact hvd w hw v hv (vertItem_inj hvw).symm
        · exact vertItem_ne_edgeItem (hvlt₂ v hv) e hve
      · rintro e he (hh | ⟨w, hw, hew⟩ | ⟨e', he', hee⟩)
        · exact hPe₂ e he hh
        · exact edgeItem_ne_vertItem (hvlt₁ w hw) e hew
        · exact hed e' he' e he (edgeItem_inj hee).symm
      rcases hi with (hi | ⟨w, hw, rfl⟩ | ⟨e, he, rfl⟩) | ⟨w, hw, rfl⟩ | ⟨e, he, rfl⟩
      · exact Or.inl hi
      · exact Or.inr (Or.inl ⟨w, List.mem_append_left _ hw, rfl⟩)
      · exact Or.inr (Or.inr ⟨e, List.mem_append_left _ he, rfl⟩)
      · exact Or.inr (Or.inl ⟨w, List.mem_append_right _ hw, rfl⟩)
      · exact Or.inr (Or.inr ⟨e, List.mem_append_right _ he, rfl⟩)

end WalkM

theorem WalkState.init_cnt (g : Graph) (ternarize : Bool) (i : ItemId) :
    (WalkState.init g ternarize).cnt i = 0 := by
  simp only [WalkState.cnt, WalkState.init, spansCount_nil, Nat.zero_add, chCount]
  exact Finset.sum_eq_zero fun j _ => by rw [Items.initialItems_ch]; rfl

theorem WalkState.init_place (g : Graph) (ternarize : Bool) :
    (WalkState.init g ternarize).Place g (fun _ => False) (fun _ => False) where
  g_eq := rfl
  size := Nat.le_of_eq (Items.initialItems_size g).symm
  root_type := by show Items.type (initialItems g) 0 = .F; rw [Items.initialItems_type]; rfl
  vert v hv := by
    have h1 : (1 + v : Nat) < 1 + g.nv := by omega
    show Items.type (initialItems g) (vertItem v) = .V
    rw [Items.initialItems_type]; simp [vertItem, h1]
  edge e he := by
    have h1 : ¬ (1 + g.nv + e : Nat) < 1 + g.nv := by omega
    have h2 : (1 + g.nv + e : Nat) < 1 + g.nv + g.ne := by omega
    show Items.type (initialItems g) (edgeItem g e) = .Q
    rw [Items.initialItems_type]; simp [edgeItem, h1, h2]
  le i := by rw [WalkState.init_cnt]; exact Nat.zero_le _
  alloc i hi := by rw [WalkState.init_cnt] at hi; exact absurd hi (Nat.lt_irrefl 0)
  root := WalkState.init_cnt g ternarize _
  fixed i _ _ hi := by rw [WalkState.init_cnt] at hi; exact absurd hi (Nat.lt_irrefl 0)
  loose _ h := h.elim
  loose_node _ h := h.elim

/-- Hypotheses on the DFS forest used by the placement invariant: vertices and edges in range and
visited once (all provided by `dfsForest_spanning`). -/
structure ForestOK (g : Graph) (forest : List DfsTree) : Prop where
  verts_lt : ∀ v ∈ forest.flatMap DfsTree.verts, v < g.nv
  edges_lt : ∀ e ∈ forest.flatMap DfsTree.edges, e < g.ne
  verts_nodup : (forest.flatMap DfsTree.verts).Nodup
  edges_nodup : (forest.flatMap DfsTree.edges).Nodup

theorem ForestOK.of_perm {g : Graph} {forest : List DfsTree}
    (hv : (forest.flatMap DfsTree.verts).Perm (List.range g.nv))
    (he : (forest.flatMap DfsTree.edges).Perm (List.range g.ne)) : ForestOK g forest where
  verts_lt _ hv' := List.mem_range.1 (hv.subset hv')
  edges_lt _ he' := List.mem_range.1 (he.subset he')
  verts_nodup := hv.nodup_iff.2 List.nodup_range
  edges_nodup := he.nodup_iff.2 List.nodup_range

/-- The span discipline at the end of the walk: every id in a span or `ch` list is allocated,
is not the root, occurs at most once over all spans and `ch` lists, and fixed items occur only if
their vertex / edge was walked. -/
theorem walk_place (g : Graph) (ternarize : Bool) (forest : List DfsTree) (hf : ForestOK g forest) :
    (g.walk ternarize forest).Place g
      (WalkM.Pushed g (fun _ => False) (forest.flatMap DfsTree.verts) (forest.flatMap DfsTree.edges))
      (fun _ => False) :=
  WalkM.walkForest_place forest (WalkState.init_place g ternarize) hf.verts_lt hf.edges_lt
    hf.verts_nodup hf.edges_nodup (fun _ _ h => h) fun _ _ h => h

section Consequences

variable (g : Graph) (ternarize : Bool) (forest : List DfsTree) (hf : ForestOK g forest)
include hf

theorem walk_chCount_le : ∀ i, chCount (g.walk ternarize forest).items i ≤ 1 := fun i =>
  Nat.le_trans (Nat.le_add_left _ _) ((walk_place g ternarize forest hf).le i)

theorem walk_ch_nodup : ∀ p, (Items.ch (g.walk ternarize forest).items p).Nodup := fun p => by
  by_cases hp : p < (g.walk ternarize forest).items.size
  · exact Items.ch_nodup_of_chCount _ (walk_chCount_le g ternarize forest hf) hp
  · rw [Items.ch_of_le _ _ (Nat.le_of_not_lt hp)]; exact List.nodup_nil

theorem walk_ch_lt : ∀ p c, Items.IsParent (g.walk ternarize forest).items p c →
    c < (g.walk ternarize forest).items.size := fun p c hpc => by
  have hp : p < (g.walk ternarize forest).items.size := by
    by_contra hp
    rw [Items.IsParent, Items.ch_of_le _ _ (Nat.le_of_not_lt hp)] at hpc
    exact List.not_mem_nil hpc
  refine (walk_place g ternarize forest hf).alloc c (Nat.lt_of_lt_of_le ?_ (Nat.le_add_left _ _))
  exact Nat.lt_of_lt_of_le (List.count_pos_iff.2 hpc) (Items.count_le_chCount _ hp c)

theorem walk_root_no_parent : ∀ p, ¬ Items.IsParent (g.walk ternarize forest).items p rootItem := fun p hp => by
  have hlt := walk_ch_lt g ternarize forest hf p rootItem hp
  have hp' : p < (g.walk ternarize forest).items.size := by
    by_contra hp'
    rw [Items.IsParent, Items.ch_of_le _ _ (Nat.le_of_not_lt hp')] at hp
    exact List.not_mem_nil hp
  have h1 := Nat.lt_of_lt_of_le (List.count_pos_iff.2 hp) (Items.count_le_chCount _ hp' rootItem)
  have h2 := (walk_place g ternarize forest hf).root
  simp only [WalkState.cnt] at h2
  omega

theorem walk_parent_unique : ∀ c p p', Items.IsParent (g.walk ternarize forest).items p c →
    Items.IsParent (g.walk ternarize forest).items p' c → p' = p := fun c p p' hp hp' => by
  by_contra hne
  have hlt : ∀ q, Items.IsParent (g.walk ternarize forest).items q c →
      q < (g.walk ternarize forest).items.size := fun q hq => by
    by_contra hq'
    rw [Items.IsParent, Items.ch_of_le _ _ (Nat.le_of_not_lt hq')] at hq
    exact List.not_mem_nil hq
  have := Items.two_le_chCount _ (hlt p hp) (hlt p' hp') (Ne.symm hne) hp hp'
  have := walk_chCount_le g ternarize forest hf c
  omega

/-- A vertex item has a parent only if its vertex was visited. -/
theorem walk_vert_parent_mem : ∀ v p, v < g.nv →
    Items.IsParent (g.walk ternarize forest).items p (vertItem v) → v ∈ forest.flatMap DfsTree.verts :=
  fun v p hv hp => by
  have hp' : p < (g.walk ternarize forest).items.size := by
    by_contra hp'
    rw [Items.IsParent, Items.ch_of_le _ _ (Nat.le_of_not_lt hp')] at hp
    exact List.not_mem_nil hp
  have h1 := Nat.lt_of_lt_of_le (List.count_pos_iff.2 hp) (Items.count_le_chCount _ hp' (vertItem v))
  rcases (walk_place g ternarize forest hf).fixed (vertItem v) (by show 0 < 1 + v; omega)
    (by show 1 + v < _; omega) (Nat.lt_of_lt_of_le h1 (Nat.le_add_left _ _)) with
    h | ⟨w, hw, hvw⟩ | ⟨e, _, hve⟩
  · exact h.elim
  · exact WalkM.vertItem_inj hvw ▸ hw
  · exact absurd hve (WalkM.vertItem_ne_edgeItem hv e)

/-- An edge item has a parent only if its edge was visited. -/
theorem walk_edge_parent_mem : ∀ e p, e < g.ne →
    Items.IsParent (g.walk ternarize forest).items p (edgeItem g e) → e ∈ forest.flatMap DfsTree.edges :=
  fun e p he hp => by
  have hp' : p < (g.walk ternarize forest).items.size := by
    by_contra hp'
    rw [Items.IsParent, Items.ch_of_le _ _ (Nat.le_of_not_lt hp')] at hp
    exact List.not_mem_nil hp
  have h1 := Nat.lt_of_lt_of_le (List.count_pos_iff.2 hp) (Items.count_le_chCount _ hp' (edgeItem g e))
  rcases (walk_place g ternarize forest hf).fixed (edgeItem g e) (by show 0 < 1 + g.nv + e; omega)
    (by show 1 + g.nv + e < _; omega) (Nat.lt_of_lt_of_le h1 (Nat.le_add_left _ _)) with
    h | ⟨w, hw, hew⟩ | ⟨e', he', hee⟩
  · exact h.elim
  · exact absurd hew (WalkM.edgeItem_ne_vertItem (hf.verts_lt w hw) e)
  · exact WalkM.edgeItem_inj hee ▸ he'

end Consequences

end Spqr
