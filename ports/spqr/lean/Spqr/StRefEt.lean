import Spqr.StRef
import Spqr.Proofs.ForestSpec

/-!
# Even–Tarjan on the reference order

Every reference block (`refBlocks`) is an st-numbering (`StBlock.St`): the proof is the
standard Even–Tarjan induction over `refTree`, phrased on the pieces of the open ear. The
invariant `VInv` says that every vertex item of the pieces other than the current vertex has a
neighbour on each side, where a back edge to the path vertex at depth `l` counts as a
neighbour on side `dirs[l]` (`Nb`), the current vertex has a neighbour on the side of its first
return, and every other vertex of the pieces lies on the side `dirs[l]` of the current vertex
for some return depth `l` (`root_side`). No tstack reasoning here.
-/

namespace Spqr

/-! Lowval / return-depth facts about a DFS out-edge (moved from `StEar.lean`, which now imports this
file: `StEar` pulls in `Spqr.Ear`, whose `Frame` would clash with `PlanarRelabelM.Frame` downstream). -/

theorem lowval_classify_tree {d : Nat} {n : Lowvals} (h : n.1 ≤ d + 1) :
    (classify d true n).lowval d = n.1 := by
  unfold classify
  split
  · next h1 =>
    split
    · next h2 => simp at h2; simp [h2, OutClass.lowval]
    · next h2 => simp at h1 h2; simp [OutClass.lowval]; omega
  · rfl

theorem lowval_classify_lt_iff_ret {d : Nat} {isTree : Bool} {n : Lowvals} :
    (classify d isTree n).lowval d < d ↔ ∃ l k, classify d isTree n = .ret l k := by
  unfold classify
  split
  · next h1 =>
    split
    · cases isTree <;> simp [OutClass.lowval]
    · simp [OutClass.lowval]
  · next h1 => simp [OutClass.lowval]; omega

namespace DfsOut

/-- The lowval of a well-formed out-edge at depth `d` is the least depth its piece returns to
(`d + 1` for a bridge). -/
theorem lowval_eq_lmin {anc : List Nat} {v : Nat} {o : DfsOut} (hwf : o.WF anc v) :
    o.cls.lowval anc.length = lmin (anc.length + 1) (o.retDepths anc.length) := by
  cases o with
  | back e dest cls =>
    rw [DfsOut.WF] at hwf
    obtain ⟨i, hi, rfl⟩ := hwf
    have hi' : i ≤ anc.length := by
      have := (List.getElem?_eq_some_iff.1 hi).1; simp at this; omega
    simp only [DfsOut.cls, retDepths, lmin_cons, lmin_nil, lowval_classify_back hi']
    omega
  | tree e cls child =>
    rw [DfsOut.WF] at hwf
    obtain ⟨-, rfl⟩ := hwf
    simp only [DfsOut.cls, retDepths]
    exact lowval_classify_tree (lmin_le _ _)

theorem lowval_mem_retDepths {anc : List Nat} {v : Nat} {o : DfsOut} (hwf : o.WF anc v)
    (hret : o.cls.lowval anc.length < anc.length) :
    o.cls.lowval anc.length ∈ o.retDepths anc.length := by
  rw [lowval_eq_lmin hwf] at hret ⊢
  rcases lmin_eq_or_mem (anc.length + 1) (o.retDepths anc.length) with h | h
  · omega
  · exact h

theorem lowval_le_retDepths {anc : List Nat} {v : Nat} {o : DfsOut} (hwf : o.WF anc v) :
    ∀ y ∈ o.retDepths anc.length, o.cls.lowval anc.length ≤ y := by
  intro y hy; rw [lowval_eq_lmin hwf]; exact lmin_le_of_mem hy

theorem ret_of_lowval_lt {anc : List Nat} {v : Nat} {o : DfsOut} (hwf : o.WF anc v)
    (hret : o.cls.lowval anc.length < anc.length) : ∃ l k, o.cls = .ret l k := by
  cases o with
  | back e dest cls =>
    rw [DfsOut.WF] at hwf
    obtain ⟨i, -, rfl⟩ := hwf
    exact lowval_classify_lt_iff_ret.1 hret
  | tree e cls child =>
    rw [DfsOut.WF] at hwf
    obtain ⟨-, rfl⟩ := hwf
    exact lowval_classify_lt_iff_ret.1 hret

end DfsOut

end Spqr

namespace Spqr.StRefEt

/-! ### Order in a list -/

/-- `a` occurs strictly before `b` in `L`. -/
def Before {α : Type} (L : List α) (a b : α) : Prop := ∃ l₁ l₂ l₃, L = l₁ ++ a :: l₂ ++ b :: l₃

namespace Before

variable {α : Type} {L L' : List α} {a b : α}

theorem mem_left (h : Before L a b) : a ∈ L := by
  obtain ⟨l₁, l₂, l₃, rfl⟩ := h; simp

theorem mem_right (h : Before L a b) : b ∈ L := by
  obtain ⟨l₁, l₂, l₃, rfl⟩ := h; simp

theorem sublist_iff : List.Sublist [a, b] L ↔ Before L a b := by
  constructor
  · intro h
    obtain ⟨r₁, r₂, rfl, ha, hb⟩ := List.cons_sublist_iff.mp h
    obtain ⟨s, t, rfl⟩ := List.append_of_mem ha
    obtain ⟨s', t', rfl⟩ := List.append_of_mem (List.singleton_sublist.mp hb)
    exact ⟨s, t ++ s', t', by simp⟩
  · rintro ⟨l₁, l₂, l₃, rfl⟩
    exact List.cons_sublist_iff.mpr ⟨l₁ ++ [a], l₂ ++ b :: l₃, by simp, by simp,
      List.singleton_sublist.mpr (by simp)⟩

theorem sub (h : Before L a b) (hL : List.Sublist L L') : Before L' a b :=
  sublist_iff.mp ((sublist_iff.mpr h).trans hL)

theorem of_mem_mem {L₁ L₂ : List α} (ha : a ∈ L₁) (hb : b ∈ L₂) : Before (L₁ ++ L₂) a b := by
  obtain ⟨s, t, rfl⟩ := List.append_of_mem ha
  obtain ⟨s', t', rfl⟩ := List.append_of_mem hb
  exact ⟨s, t ++ s', t', by simp⟩

theorem filterMap {β : Type} {f : α → Option β} {a' b' : β} (h : Before L a b)
    (ha : f a = some a') (hb : f b = some b') : Before (L.filterMap f) a' b' := by
  obtain ⟨l₁, l₂, l₃, rfl⟩ := h
  exact ⟨l₁.filterMap f, l₂.filterMap f, l₃.filterMap f, by
    simp [List.filterMap_append, ha, hb]⟩

theorem idxOf_lt [DecidableEq α] (hnd : L.Nodup) (h : Before L a b) : L.idxOf a < L.idxOf b := by
  obtain ⟨l₁, l₂, l₃, rfl⟩ := h
  have hb : b ∉ l₁ ++ a :: l₂ := fun hb =>
    (List.nodup_cons.mp (List.nodup_middle.mp hnd)).1 (List.mem_append_left _ hb)
  have ha : a ∈ l₁ ++ a :: l₂ := by simp
  rw [List.idxOf_append_of_mem ha, List.idxOf_append_of_notMem hb]
  have := List.idxOf_lt_length_of_mem ha
  rw [List.idxOf_cons_self]; omega

theorem asymm [DecidableEq α] (hnd : L.Nodup) (h₁ : Before L a b) (h₂ : Before L b a) : False := by
  have := h₁.idxOf_lt hnd; have := h₂.idxOf_lt hnd; omega

theorem ne [DecidableEq α] (hnd : L.Nodup) (h : Before L a b) : a ≠ b := by
  rintro rfl; exact h.asymm hnd h

end Before

/-- If every other element precedes `w`, `w` is last. -/
theorem getLast?_eq_of_before {α : Type} [DecidableEq α] {L : List α} {w : α} (hnd : L.Nodup)
    (hw : w ∈ L) (h : ∀ y ∈ L, y ≠ w → Before L y w) : L.getLast? = some w := by
  obtain ⟨s, t, rfl⟩ := List.append_of_mem hw
  cases t with
  | nil => simp
  | cons y t =>
    have hyw : y ≠ w := by
      intro hyw; subst hyw
      exact (List.nodup_cons.mp (List.nodup_middle.mp hnd)).1 (by simp)
    have h₁ : Before (s ++ w :: y :: t) w y := ⟨s, [], t, by simp⟩
    exact (h₁.asymm hnd (h y (by simp) hyw)).elim

/-! ### Items of a block -/

/-- The vertex of a vertex item (`StBlock.seq` is `filterMap (vertOf g)` behind the root). -/
def vertOf (g : Graph) (x : ItemId) : Option Nat :=
  if 1 ≤ x ∧ x < 1 + g.nv then some (x - 1) else none

/-- The edge of an edge item (`StBlock.edges` is `filterMap (edgeOf g)` behind the root). -/
def edgeOf (g : Graph) (x : ItemId) : Option (Nat × Nat) :=
  if 1 + g.nv ≤ x ∧ x < 1 + g.nv + g.ne then some g.edges[x - (1 + g.nv)]! else none

theorem _root_.Spqr.StBlock.seq_eq (g : Graph) (b : StBlock) :
    b.seq g = (b.root.map (·.1)).toList ++ b.items.filterMap (vertOf g) := rfl

theorem _root_.Spqr.StBlock.edges_eq (g : Graph) (b : StBlock) :
    b.edges g = b.root.toList ++ b.items.filterMap (edgeOf g) := rfl

theorem vertOf_vertItem {g : Graph} {x : Nat} (hx : x < g.nv) : vertOf g (vertItem x) = some x := by
  simp [vertOf, vertItem]; omega

theorem vertOf_edgeItem (g : Graph) (e : Nat) : vertOf g (edgeItem g e) = none := by
  have : ¬ (1 ≤ edgeItem g e ∧ edgeItem g e < 1 + g.nv) := by
    show ¬ (1 ≤ 1 + g.nv + e ∧ 1 + g.nv + e < 1 + g.nv); omega
  simp [vertOf, this]

theorem edgeOf_vertItem (g : Graph) {x : Nat} (hx : x < g.nv) : edgeOf g (vertItem x) = none := by
  have : ¬ (1 + g.nv ≤ vertItem x ∧ vertItem x < 1 + g.nv + g.ne) := by
    show ¬ (1 + g.nv ≤ 1 + x ∧ 1 + x < 1 + g.nv + g.ne); omega
  simp [edgeOf, this]

theorem edgeOf_edgeItem (g : Graph) {e : Nat} (he : e < g.ne) :
    edgeOf g (edgeItem g e) = some g.edges[e]! := by
  simp [edgeOf, edgeItem]; omega

theorem vertItem_inj {x y : Nat} (h : vertItem x = vertItem y) : x = y := by
  simp [vertItem] at h; exact h

theorem vertItem_ne_edgeItem {g : Graph} {x e : Nat} (hx : x < g.nv) : vertItem x ≠ edgeItem g e := by
  intro h; have h' : 1 + x = 1 + g.nv + e := h; omega

/-- Every item is a vertex item of a vertex or an edge item of an edge. -/
def ItemsShape (g : Graph) (L : List ItemId) : Prop :=
  ∀ i ∈ L, (∃ x, x < g.nv ∧ i = vertItem x) ∨ ∃ e, e < g.ne ∧ i = edgeItem g e

theorem mem_filterMap_vertOf {g : Graph} {L : List ItemId} (hs : ItemsShape g L) {x : Nat}
    (h : x ∈ L.filterMap (vertOf g)) : vertItem x ∈ L ∧ x < g.nv := by
  obtain ⟨i, hi, hx⟩ := List.mem_filterMap.mp h
  rcases hs i hi with ⟨y, hy, rfl⟩ | ⟨e, he, rfl⟩
  · rw [vertOf_vertItem hy] at hx; cases hx; exact ⟨hi, hy⟩
  · rw [vertOf_edgeItem] at hx; cases hx

theorem mem_filterMap_vertOf_of {g : Graph} {L : List ItemId} {x : Nat} (hx : x < g.nv)
    (h : vertItem x ∈ L) : x ∈ L.filterMap (vertOf g) :=
  List.mem_filterMap.mpr ⟨_, h, vertOf_vertItem hx⟩

theorem mem_filterMap_edgeOf {g : Graph} {L : List ItemId} (hs : ItemsShape g L) (p : Nat × Nat) :
    p ∈ L.filterMap (edgeOf g) ↔ ∃ e, e < g.ne ∧ edgeItem g e ∈ L ∧ g.edges[e]! = p := by
  rw [List.mem_filterMap]
  constructor
  · rintro ⟨i, hi, hp⟩
    rcases hs i hi with ⟨y, hy, rfl⟩ | ⟨e, he, rfl⟩
    · rw [edgeOf_vertItem g hy] at hp; cases hp
    · rw [edgeOf_edgeItem g he] at hp; cases hp; exact ⟨e, he, hi, rfl⟩
  · rintro ⟨e, he, hi, rfl⟩
    exact ⟨_, hi, edgeOf_edgeItem g he⟩

theorem ItemsShape.sub {g : Graph} {L L' : List ItemId} (hs : ItemsShape g L') (h : List.Sublist L L') :
    ItemsShape g L := fun i hi => hs i (h.subset hi)

theorem ItemsShape.append {g : Graph} {L L' : List ItemId} (h : ItemsShape g L)
    (h' : ItemsShape g L') : ItemsShape g (L ++ L') := by
  intro i hi
  rcases List.mem_append.mp hi with hi | hi
  · exact h i hi
  · exact h' i hi

/-! ### Neighbours on a side -/

/-- Edge `e` joins `x` and `y`. -/
def Joins (g : Graph) (e x y : Nat) : Prop := g.edges[e]! = (x, y) ∨ g.edges[e]! = (y, x)

theorem Joins.symm {g : Graph} {e x y : Nat} (h : Joins g e x y) : Joins g e y x := Or.symm h

/-- `y` lies on side `s` of `x` among the vertex items of `L` (`false`: before, `true`: after). -/
def Side (L : List ItemId) (s : Bool) (x y : Nat) : Prop :=
  if s then Before L (vertItem x) (vertItem y) else Before L (vertItem y) (vertItem x)

theorem Side.flip {L : List ItemId} {s : Bool} {x y : Nat} (h : Side L s x y) : Side L (!s) y x := by
  cases s <;> simpa [Side] using h

theorem Side.sub {L L' : List ItemId} {s : Bool} {x y : Nat} (h : Side L s x y) (hL : List.Sublist L L') :
    Side L' s x y := by
  cases s <;> simp only [Side] at h ⊢ <;> exact h.sub hL

theorem Side.mem_right {L : List ItemId} {s : Bool} {x y : Nat} (h : Side L s x y) :
    vertItem y ∈ L := by
  cases s <;> simp only [Side] at h
  · exact h.mem_left
  · exact h.mem_right

/-- `x` has a neighbour on side `s`: an edge item of `L` joining `x` to a vertex item of `L` on
that side, or to the path vertex `anc[l]` of a return depth `l`, which lies on side `dirs[l]`. -/
def Nb (g : Graph) (anc : List Nat) (dirs : List Bool) (rets : List Nat) (L : List ItemId)
    (x : Nat) (s : Bool) : Prop :=
  ∃ e y, edgeItem g e ∈ L ∧ Joins g e x y ∧
    (Side L s x y ∨ ∃ l, anc[l]? = some y ∧ dirs[l]? = some s ∧ l ∈ rets)

theorem Nb.mono {g : Graph} {anc : List Nat} {dirs : List Bool} {rets rets' : List Nat}
    {L L' : List ItemId} {x : Nat} {s : Bool} (h : Nb g anc dirs rets L x s) (hL : List.Sublist L L')
    (hr : ∀ l ∈ rets, l ∈ rets') : Nb g anc dirs rets' L' x s := by
  obtain ⟨e, y, he, hj, h⟩ := h
  refine ⟨e, y, hL.subset he, hj, ?_⟩
  rcases h with h | ⟨l, h1, h2, h3⟩
  · exact .inl (h.sub hL)
  · exact .inr ⟨l, h1, h2, hr l h3⟩

/-- Leaving the vertex `v`: a neighbour "at the path vertex `v`" becomes an in-list neighbour. -/
theorem Nb.pop {g : Graph} {anc : List Nat} {dirs : List Bool} {v : Nat} {sd : Bool}
    {rets rets' : List Nat} {L L' : List ItemId} {x : Nat} {s : Bool}
    (hlen : dirs.length = anc.length) (h : Nb g (anc ++ [v]) (dirs ++ [sd]) rets L x s)
    (hL : List.Sublist L L') (hr : ∀ l ∈ rets, l ∈ rets') (hv : Side L' sd x v) :
    Nb g anc dirs rets' L' x s := by
  obtain ⟨e, y, he, hj, h⟩ := h
  refine ⟨e, y, hL.subset he, hj, ?_⟩
  rcases h with h | ⟨l, h1, h2, h3⟩
  · exact .inl (h.sub hL)
  · by_cases hl : l < anc.length
    · rw [List.getElem?_append_left hl] at h1
      rw [List.getElem?_append_left (hlen ▸ hl)] at h2
      exact .inr ⟨l, h1, h2, hr l h3⟩
    · have hl' : l = anc.length := by
        have := (List.getElem?_eq_some_iff.mp h1).1
        simp only [List.length_append, List.length_singleton] at this; omega
      subst hl'
      rw [← hlen, List.getElem?_append_right (Nat.le_refl _)] at h2
      simp at h1 h2
      subst h1; subst h2
      exact .inl hv

/-- The edge endpoints of a well-formed graph are vertices. -/
theorem _root_.Spqr.Graph.WF.ends_lt {g : Graph} (hg : g.WF) {e : Nat} (he : e < g.ne) :
    (g.edges[e]!).1 < g.nv ∧ (g.edges[e]!).2 < g.nv := by
  rw [getElem!_pos g.edges e he]
  exact hg _ (Array.getElem_mem he)

theorem Joins.lt {g : Graph} (hg : g.WF) {e x y : Nat} (he : e < g.ne) (h : Joins g e x y) :
    x < g.nv ∧ y < g.nv := by
  have := hg.ends_lt he
  rcases h with h | h <;> rw [h] at this <;> simp at this <;> omega

/-! ### The invariant of the open ear -/

/-- The pieces `L` (the open ear at the vertex `v` of depth `anc.length`, below the path `anc`
with directions `dirs`) after processing the out-edges whose subtree vertices are `vs` and whose
return depths are `rets`: every vertex item other than `v` has a neighbour on each side (`nb`), `v`
has one on the side `dirs[l]` of its first return `l` (`root_nb`), and every other vertex lies on
the side `dirs[l]` of `v` for a return depth `l` (`root_side`). `hasVert = false` means nothing has
been pushed yet (all out-edges so far were block boundaries). -/
structure VInv (g : Graph) (anc : List Nat) (dirs : List Bool) (v : Nat) (vs rets : List Nat)
    (hasVert : Bool) (L : List ItemId) : Prop where
  shape : ItemsShape g L
  vert_mem : ∀ x, x < g.nv → vertItem x ∈ L → x = v ∨ x ∈ vs
  vert_nodup : (L.filterMap (vertOf g)).Nodup
  edge_mem : ∀ e, edgeItem g e ∈ L → ∃ x y, Joins g e x y ∧ x ≠ y ∧ vertItem x ∈ L ∧
    (vertItem y ∈ L ∨ ∃ l, anc[l]? = some y ∧ l ∈ rets)
  nb : ∀ x, x < g.nv → vertItem x ∈ L → x ≠ v →
    Nb g anc dirs rets L x false ∧ Nb g anc dirs rets L x true
  root_side : ∀ y, y < g.nv → vertItem y ∈ L → y ≠ v →
    ∃ l s, l ∈ rets ∧ dirs[l]? = some s ∧ Side L s v y
  root_nb : hasVert = true → vertItem v ∈ L ∧ ((∃ l ∈ rets, l < anc.length) →
    ∃ l s, l ∈ rets ∧ (∀ l' ∈ rets, l ≤ l') ∧ dirs[l]? = some s ∧ Nb g anc dirs rets L v s)
  empty : hasVert = false → L = [] ∧ ∀ l ∈ rets, anc.length ≤ l

theorem VInv.nil (g : Graph) (anc : List Nat) (dirs : List Bool) (v : Nat) :
    VInv g anc dirs v [] [] false [] where
  shape := fun _ h => absurd h List.not_mem_nil
  vert_mem := fun _ _ h => absurd h List.not_mem_nil
  vert_nodup := List.nodup_nil
  edge_mem := fun _ h => absurd h List.not_mem_nil
  nb := fun _ _ h => absurd h List.not_mem_nil
  root_side := fun _ _ h => absurd h List.not_mem_nil
  root_nb := fun h => by cases h
  empty := fun _ => ⟨rfl, fun _ h => absurd h List.not_mem_nil⟩

/-- Growing `vs` and adding return depths `≥ anc.length` (a block boundary) keeps the invariant. -/
theorem VInv.extend {g : Graph} {anc : List Nat} {dirs : List Bool} {v : Nat} {vs vs' rets rets' : List Nat}
    {hasVert : Bool} {L : List ItemId} (hlen : dirs.length = anc.length)
    (h : VInv g anc dirs v vs rets hasVert L)
    (hvs : ∀ x ∈ vs, x ∈ vs') (hr : ∀ l ∈ rets, l ∈ rets')
    (hr' : ∀ l ∈ rets', l ∈ rets ∨ anc.length ≤ l) : VInv g anc dirs v vs' rets' hasVert L where
  shape := h.shape
  vert_mem := fun x hx hxL => (h.vert_mem x hx hxL).imp_right (hvs x)
  vert_nodup := h.vert_nodup
  edge_mem := fun e he => by
    obtain ⟨x, y, hj, hne, hx, hy⟩ := h.edge_mem e he
    exact ⟨x, y, hj, hne, hx, hy.imp_right fun ⟨l, h1, h2⟩ => ⟨l, h1, hr l h2⟩⟩
  nb := fun x hx hxL hne =>
    ⟨(h.nb x hx hxL hne).1.mono (List.Sublist.refl _) hr,
      (h.nb x hx hxL hne).2.mono (List.Sublist.refl _) hr⟩
  root_side := fun y hy hyL hne => by
    obtain ⟨l, s, h1, h2, h3⟩ := h.root_side y hy hyL hne
    exact ⟨l, s, hr l h1, h2, h3⟩
  root_nb := fun hv => by
    obtain ⟨hvL, hret⟩ := h.root_nb hv
    refine ⟨hvL, fun ⟨l, hl, hld⟩ => ?_⟩
    have hl' : l ∈ rets := by
      rcases hr' l hl with h | h
      · exact h
      · omega
    obtain ⟨l₀, s, h1, h2, h3, h4⟩ := hret ⟨l, hl', hld⟩
    have hl₀ : l₀ < anc.length := by
      have := (List.getElem?_eq_some_iff.mp h3).1; omega
    refine ⟨l₀, s, hr l₀ h1, fun l' hl' => ?_, h3, h4.mono (List.Sublist.refl _) hr⟩
    rcases hr' l' hl' with h | h
    · exact h2 l' h
    · omega
  empty := fun hv => by
    obtain ⟨hL, hret⟩ := h.empty hv
    refine ⟨hL, fun l hl => ?_⟩
    rcases hr' l hl with h | h
    · exact hret l h
    · exact h

/-! ### One returning out-edge -/

/-- Neighbours inside the material `M` pushed for one returning out-edge of `v`: as `Nb`, but `v`
itself counts as a neighbour on side `sd` (`v` will lie on side `sd` of all of `M`). -/
def NbV (g : Graph) (anc : List Nat) (dirs : List Bool) (rets : List Nat) (v : Nat) (sd : Bool)
    (M : List ItemId) (x : Nat) (s : Bool) : Prop :=
  ∃ e y, edgeItem g e ∈ M ∧ Joins g e x y ∧
    (Side M s x y ∨ (y = v ∧ s = sd) ∨ ∃ l, anc[l]? = some y ∧ dirs[l]? = some s ∧ l ∈ rets)

theorem NbV.toNb {g : Graph} {anc : List Nat} {dirs : List Bool} {rets rets' : List Nat} {v : Nat}
    {sd : Bool} {M L' : List ItemId} {x : Nat} {s : Bool} (h : NbV g anc dirs rets v sd M x s)
    (hM : List.Sublist M L') (hr : ∀ l ∈ rets, l ∈ rets') (hv : Side L' sd x v) :
    Nb g anc dirs rets' L' x s := by
  obtain ⟨e, y, he, hj, h⟩ := h
  refine ⟨e, y, hM.subset he, hj, ?_⟩
  rcases h with h | ⟨rfl, rfl⟩ | ⟨l, h1, h2, h3⟩
  · exact .inl (h.sub hM)
  · exact .inl hv
  · exact .inr ⟨l, h1, h2, hr l h3⟩

/-- A neighbour fact below `v` (path `anc ++ [v]`, directions `dirs ++ [sd]`) read at `v`. -/
theorem Nb.popV {g : Graph} {anc : List Nat} {dirs : List Bool} {v : Nat} {sd : Bool}
    {rets rets' : List Nat} {L M : List ItemId} {x : Nat} {s : Bool}
    (hlen : dirs.length = anc.length) (h : Nb g (anc ++ [v]) (dirs ++ [sd]) rets L x s)
    (hL : List.Sublist L M) (hr : ∀ l ∈ rets, l ∈ rets') : NbV g anc dirs rets' v sd M x s := by
  obtain ⟨e, y, he, hj, h⟩ := h
  refine ⟨e, y, hL.subset he, hj, ?_⟩
  rcases h with h | ⟨l, h1, h2, h3⟩
  · exact .inl (h.sub hL)
  · by_cases hl : l < anc.length
    · rw [List.getElem?_append_left hl] at h1
      rw [List.getElem?_append_left (hlen ▸ hl)] at h2
      exact .inr (.inr ⟨l, h1, h2, hr l h3⟩)
    · have hl' : l = anc.length := by
        have := (List.getElem?_eq_some_iff.mp h1).1
        simp only [List.length_append, List.length_singleton] at this; omega
      subst hl'
      rw [← hlen, List.getElem?_append_right (Nat.le_refl _)] at h2
      simp at h1 h2
      exact .inr (.inl ⟨h1.symm, h2.symm⟩)

/-- Processing one returning out-edge of `v` with lowval `l`: the material `M` is placed on side
`lowDir = dirs[l]` of the previous pieces (or of `v`'s own vertex item if nothing was pushed yet). -/
theorem VInv.step {g : Graph} {anc : List Nat} {dirs : List Bool} {v : Nat} {vs rets : List Nat}
    {hasVert : Bool} {L : List ItemId} (hlen : dirs.length = anc.length)
    (hinv : VInv g anc dirs v vs rets hasVert L) (hv : v < g.nv)
    {l : Nat} {lowDir : Bool} (hl : dirs[l]? = some lowDir)
    {cvs crets : List Nat} {M L' : List ItemId}
    (hL' : L' = if lowDir then (if hasVert then L else [vertItem v]) ++ M
      else M ++ (if hasVert then L else [vertItem v]))
    (hshape : ItemsShape g M)
    (hMv : ∀ x, x < g.nv → vertItem x ∈ M → x ∈ cvs)
    (hMnd : (M.filterMap (vertOf g)).Nodup)
    (hdisj : ∀ x ∈ cvs, x ≠ v ∧ x ∉ vs)
    (hMe : ∀ e, edgeItem g e ∈ M → ∃ x y, Joins g e x y ∧ x ≠ y ∧ (vertItem x ∈ M ∨ x = v) ∧
      (vertItem y ∈ M ∨ y = v ∨ ∃ l', anc[l']? = some y ∧ l' ∈ crets))
    (hMnb : ∀ x, x < g.nv → vertItem x ∈ M →
      NbV g anc dirs crets v (!lowDir) M x false ∧ NbV g anc dirs crets v (!lowDir) M x true)
    (hroot : ∃ e y, edgeItem g e ∈ M ∧ Joins g e v y ∧ (vertItem y ∈ M ∨ anc[l]? = some y))
    (hlc : l ∈ crets) (hlc' : ∀ l' ∈ crets, l ≤ l')
    (hsorted : (∃ l₀ ∈ rets, l₀ < anc.length) → ∃ l₁ ∈ rets, l₁ ≤ l) :
    VInv g anc dirs v (vs ++ cvs) (rets ++ crets) true L' := by
  have hd : l < anc.length := by have := (List.getElem?_eq_some_iff.mp hl).1; omega
  set L₀ : List ItemId := if hasVert then L else [vertItem v] with hL₀
  have hL₀L : ∀ i ∈ L₀, (hasVert = true ∧ i ∈ L) ∨ i = vertItem v := by
    intro i hi; cases hasVert <;> simp_all
  have hvL₀ : vertItem v ∈ L₀ := by
    cases hasVert
    · simp [hL₀]
    · simpa [hL₀] using (hinv.root_nb rfl).1
  have hLL₀ : hasVert = true → L₀ = L := fun h => by simp [hL₀, h]
  have hL₀sub : List.Sublist L₀ L' := by
    cases lowDir <;> simp [hL']
  have hMsub : List.Sublist M L' := by
    cases lowDir <;> simp [hL']
  have hmem : ∀ i ∈ L', i ∈ L₀ ∨ i ∈ M := by
    intro i hi; cases lowDir <;> simp [hL'] at hi <;> tauto
  have hsideM : ∀ x, vertItem x ∈ M → Side L' (!lowDir) x v := by
    intro x hx
    cases lowDir
    · show Before L' (vertItem x) (vertItem v)
      rw [hL']; simpa using Before.of_mem_mem hx hvL₀
    · show Before L' (vertItem v) (vertItem x)
      rw [hL']; simpa using Before.of_mem_mem hvL₀ hx
  have hshape₀ : ItemsShape g L₀ := by
    intro i hi
    rcases hL₀L i hi with ⟨h, hi⟩ | rfl
    · exact hinv.shape i hi
    · exact .inl ⟨v, hv, rfl⟩
  have hvert₀ : ∀ x, x < g.nv → vertItem x ∈ L₀ → x = v ∨ x ∈ vs := by
    intro x hxlt hx
    rcases hL₀L _ hx with ⟨h, hx⟩ | hx
    · exact hinv.vert_mem x hxlt hx
    · exact .inl (vertItem_inj hx)
  have hnd₀ : (L₀.filterMap (vertOf g)).Nodup := by
    cases hasVert
    · simp [hL₀, vertOf_vertItem hv]
    · simpa [hL₀] using hinv.vert_nodup
  have hedge₀ : ∀ e, edgeItem g e ∈ L₀ → hasVert = true ∧ edgeItem g e ∈ L := by
    intro e he
    rcases hL₀L _ he with h | h
    · exact h
    · exact absurd h.symm (vertItem_ne_edgeItem hv)
  have hrets : ∀ l' ∈ rets, l' ∈ rets ++ crets := fun l' h => List.mem_append_left _ h
  have hcrets : ∀ l' ∈ crets, l' ∈ rets ++ crets := fun l' h => List.mem_append_right _ h
  have hshape' : ItemsShape g L' := by
    cases lowDir <;> simp only [hL', ite_true, ite_false, Bool.false_eq_true] <;>
      first | exact hshape₀.append hshape | exact hshape.append hshape₀
  refine
    { shape := hshape'
      vert_mem := ?_
      vert_nodup := ?_
      edge_mem := ?_
      nb := ?_
      root_side := ?_
      root_nb := ?_
      empty := fun h => nomatch h }
  · intro x hxlt hx
    rcases hmem _ hx with hx | hx
    · exact (hvert₀ x hxlt hx).imp_right (List.mem_append_left _)
    · exact .inr (List.mem_append_right _ (hMv x hxlt hx))
  · have hdisj' : (L₀.filterMap (vertOf g)).Disjoint (M.filterMap (vertOf g)) := by
      intro x h1 h2
      obtain ⟨h1, hxlt⟩ := mem_filterMap_vertOf hshape₀ h1
      have h2 := (mem_filterMap_vertOf hshape h2).1
      have := hdisj x (hMv x hxlt h2)
      rcases hvert₀ x hxlt h1 with h | h
      · exact this.1 h
      · exact this.2 h
    cases lowDir
    · simp only [hL', Bool.false_eq_true, ite_false, List.filterMap_append]
      exact List.Nodup.append hMnd hnd₀ (fun x h1 h2 => hdisj' h2 h1)
    · simp only [hL', ite_true, List.filterMap_append]
      exact List.Nodup.append hnd₀ hMnd hdisj'
  · intro e he
    rcases hmem _ he with he | he
    · obtain ⟨hhv, he⟩ := hedge₀ e he
      obtain ⟨x, y, hj, hne, hx, hy⟩ := hinv.edge_mem e he
      have hLsub : List.Sublist L L' := (hLL₀ hhv) ▸ hL₀sub
      refine ⟨x, y, hj, hne, hLsub.subset hx, ?_⟩
      rcases hy with hy | ⟨l', h1, h2⟩
      · exact .inl (hLsub.subset hy)
      · exact .inr ⟨l', h1, hrets _ h2⟩
    · obtain ⟨x, y, hj, hne, hx, hy⟩ := hMe e he
      refine ⟨x, y, hj, hne, ?_, ?_⟩
      · rcases hx with hx | rfl
        · exact hMsub.subset hx
        · exact hL₀sub.subset hvL₀
      · rcases hy with hy | rfl | ⟨l', h1, h2⟩
        · exact .inl (hMsub.subset hy)
        · exact .inl (hL₀sub.subset hvL₀)
        · exact .inr ⟨l', h1, hcrets _ h2⟩
  · intro x hxlt hx hne
    rcases hmem _ hx with hx | hx
    · rcases hL₀L _ hx with ⟨hhv, hx⟩ | hx
      · have hLsub : List.Sublist L L' := (hLL₀ hhv) ▸ hL₀sub
        exact ⟨(hinv.nb x hxlt hx hne).1.mono hLsub hrets,
          (hinv.nb x hxlt hx hne).2.mono hLsub hrets⟩
      · exact absurd (vertItem_inj hx) hne
    · exact ⟨(hMnb x hxlt hx).1.toNb hMsub hcrets (hsideM x hx),
        (hMnb x hxlt hx).2.toNb hMsub hcrets (hsideM x hx)⟩
  · intro y hylt hy hne
    rcases hmem _ hy with hy | hy
    · rcases hL₀L _ hy with ⟨hhv, hy⟩ | hy
      · have hLsub : List.Sublist L L' := (hLL₀ hhv) ▸ hL₀sub
        obtain ⟨l', s, h1, h2, h3⟩ := hinv.root_side y hylt hy hne
        exact ⟨l', s, hrets _ h1, h2, h3.sub hLsub⟩
      · exact absurd (vertItem_inj hy) hne
    · refine ⟨l, lowDir, hcrets _ hlc, hl, ?_⟩
      simpa using (hsideM y hy).flip
  · intro _
    refine ⟨hL₀sub.subset hvL₀, fun _ => ?_⟩
    by_cases hex : ∃ l₀ ∈ rets, l₀ < anc.length
    · have hhv : hasVert = true := by
        cases h : hasVert
        · obtain ⟨l₀, h1, h2⟩ := hex
          have := (hinv.empty h).2 l₀ h1; omega
        · rfl
      obtain ⟨l₀, s, h1, h2, h3, h4⟩ := (hinv.root_nb hhv).2 hex
      have hl₀ : l₀ < anc.length := by
        have := (List.getElem?_eq_some_iff.mp h3).1; omega
      have hLsub : List.Sublist L L' := (hLL₀ hhv) ▸ hL₀sub
      refine ⟨l₀, s, hrets _ h1, fun l' hl' => ?_, h3, h4.mono hLsub hrets⟩
      rcases List.mem_append.mp hl' with hl' | hl'
      · exact h2 l' hl'
      · obtain ⟨l₁, hl₁, hle⟩ := hsorted hex
        exact Nat.le_trans (h2 l₁ hl₁) (Nat.le_trans hle (hlc' l' hl'))
    · refine ⟨l, lowDir, hcrets _ hlc, fun l' hl' => ?_, hl, ?_⟩
      · rcases List.mem_append.mp hl' with hl' | hl'
        · by_contra hlt
          exact hex ⟨l', hl', by omega⟩
        · exact hlc' l' hl'
      · obtain ⟨e, y, he, hj, hy⟩ := hroot
        refine ⟨e, y, hMsub.subset he, hj, ?_⟩
        rcases hy with hy | hy
        · exact .inl (by simpa using (hsideM y hy).flip)
        · exact .inr ⟨l, hy, hl, hcrets _ hlc⟩

/-- The path clause of `Nb` / `edge_mem` read at the block boundary: `anc ++ [v]` indexed by a
return depth `≥ anc.length` is `v`. -/
theorem path_clause_boundary {anc : List Nat} {v : Nat} {rets : List Nat}
    (hrets : ∀ l ∈ rets, anc.length ≤ l) {l y : Nat} (h1 : (anc ++ [v])[l]? = some y) (h3 : l ∈ rets) :
    y = v ∧ l = anc.length := by
  have hl := (List.getElem?_eq_some_iff.mp h1).1
  simp only [List.length_append, List.length_singleton] at hl
  have hl' : l = anc.length := by have := hrets l h3; omega
  subst hl'
  rw [List.getElem?_append_right (Nat.le_refl _)] at h1
  simp at h1
  exact ⟨h1.symm, rfl⟩

/-- The back edge `e` of `v` to the path vertex `dest = anc[l]`. -/
theorem VInv.step_back {g : Graph} {anc : List Nat} {dirs : List Bool} {v : Nat} {vs rets : List Nat}
    {hasVert : Bool} {L : List ItemId} (hlen : dirs.length = anc.length)
    (hinv : VInv g anc dirs v vs rets hasVert L) (hv : v < g.nv)
    {l : Nat} {lowDir : Bool} (hl : dirs[l]? = some lowDir) {dest e : Nat}
    (hdest : anc[l]? = some dest) (hne : dest ≠ v) (he : e < g.ne) (hj : Joins g e v dest)
    {L' : List ItemId}
    (hL' : L' = if lowDir then (if hasVert then L else [vertItem v]) ++ [edgeItem g e]
      else [edgeItem g e] ++ (if hasVert then L else [vertItem v]))
    (hsorted : (∃ l₀ ∈ rets, l₀ < anc.length) → ∃ l₁ ∈ rets, l₁ ≤ l) :
    VInv g anc dirs v (vs ++ []) (rets ++ [l]) true L' := by
  refine hinv.step hlen hv hl hL' ?_ ?_ ?_ ?_ ?_ ?_ ?_ ?_ ?_ hsorted
  · intro i hi; simp at hi; exact .inr ⟨e, he, hi⟩
  · intro x hx h; simp at h; exact absurd h (vertItem_ne_edgeItem hx)
  · simp [vertOf_edgeItem]
  · intro x hx; simp at hx
  · intro e' he'
    simp at he'
    have : e' = e := by simp [edgeItem] at he'; omega
    subst this
    exact ⟨v, dest, hj, hne.symm, .inr rfl, .inr (.inr ⟨l, hdest, by simp⟩)⟩
  · intro x hx h; simp at h; exact absurd h (vertItem_ne_edgeItem hx)
  · exact ⟨e, dest, by simp, hj, .inr hdest⟩
  · simp
  · intro l' hl'; simp at hl'; omega

/-- The returning tree edge `e` of `v` to the child `c`, whose pieces `Lc` satisfy the invariant
below `v` (path `anc ++ [v]`, directions `dirs ++ [!lowDir]`). -/
theorem VInv.step_tree {g : Graph} {anc : List Nat} {dirs : List Bool} {v : Nat} {vs rets : List Nat}
    {hasVert : Bool} {L : List ItemId} (hlen : dirs.length = anc.length)
    (hinv : VInv g anc dirs v vs rets hasVert L) (hv : v < g.nv)
    {l : Nat} {lowDir : Bool} (hl : dirs[l]? = some lowDir) {c e : Nat} {cvs crets : List Nat}
    {Lc : List ItemId} (hc : VInv g (anc ++ [v]) (dirs ++ [!lowDir]) c cvs crets true Lc)
    (hdisj : ∀ x ∈ c :: cvs, x ≠ v ∧ x ∉ vs)
    (he : e < g.ne) (hj : Joins g e v c) (hlc : l ∈ crets) (hlc' : ∀ l' ∈ crets, l ≤ l')
    {L' : List ItemId}
    (hL' : L' = if lowDir then (if hasVert then L else [vertItem v]) ++ ([edgeItem g e] ++ Lc)
      else (Lc ++ [edgeItem g e]) ++ (if hasVert then L else [vertItem v]))
    (hsorted : (∃ l₀ ∈ rets, l₀ < anc.length) → ∃ l₁ ∈ rets, l₁ ≤ l) :
    VInv g anc dirs v (vs ++ c :: cvs) (rets ++ crets) true L' := by
  have hd : l < anc.length := by have := (List.getElem?_eq_some_iff.mp hl).1; omega
  set M : List ItemId := if lowDir then [edgeItem g e] ++ Lc else Lc ++ [edgeItem g e] with hM
  have hLcM : List.Sublist Lc M := by cases lowDir <;> simp [hM]
  have hQM : edgeItem g e ∈ M := by cases lowDir <;> simp [hM]
  have hmemM : ∀ i ∈ M, i = edgeItem g e ∨ i ∈ Lc := by
    intro i hi; cases lowDir <;> simp [hM] at hi <;> tauto
  have hVM : ∀ x, x < g.nv → vertItem x ∈ M → vertItem x ∈ Lc := by
    intro x hx h
    rcases hmemM _ h with h | h
    · exact absurd h (vertItem_ne_edgeItem hx)
    · exact h
  have hcL : vertItem c ∈ Lc := (hc.root_nb rfl).1
  have hL'' : L' = if lowDir then (if hasVert then L else [vertItem v]) ++ M
      else M ++ (if hasVert then L else [vertItem v]) := by
    cases lowDir <;> simpa [hM] using hL'
  refine hinv.step hlen hv hl hL'' ?_ ?_ ?_ hdisj ?_ ?_ ?_ hlc hlc' hsorted
  · intro i hi
    rcases hmemM _ hi with rfl | hi
    · exact .inr ⟨e, he, rfl⟩
    · exact hc.shape i hi
  · intro x hx h
    rcases hc.vert_mem x hx (hVM x hx h) with h | h
    · exact h ▸ List.mem_cons_self ..
    · exact List.mem_cons_of_mem _ h
  · cases lowDir <;> simpa [hM, vertOf_edgeItem] using hc.vert_nodup
  · intro e' he'
    rcases hmemM _ he' with he' | he'
    · have : e' = e := by simp [edgeItem] at he'; omega
      subst this
      exact ⟨v, c, hj, (hdisj c (List.mem_cons_self ..)).1.symm, .inr rfl, .inl (hLcM.subset hcL)⟩
    · obtain ⟨x, y, hj', hne, hx, hy⟩ := hc.edge_mem e' he'
      refine ⟨x, y, hj', hne, .inl (hLcM.subset hx), ?_⟩
      rcases hy with hy | ⟨l', h1, h2⟩
      · exact .inl (hLcM.subset hy)
      · by_cases hl' : l' < anc.length
        · rw [List.getElem?_append_left hl'] at h1
          exact .inr (.inr ⟨l', h1, h2⟩)
        · have := path_clause_boundary (rets := [l']) (fun _ h => by simp at h; omega) h1 (by simp)
          exact .inr (.inl this.1)
  · intro x hx h
    have hxL := hVM x hx h
    by_cases hxc : x = c
    · rw [hxc]
      obtain ⟨_, h2⟩ := hc.root_nb rfl
      obtain ⟨l₀, s, h1, hmin, h3, h4⟩ := h2 ⟨l, hlc, by simp; omega⟩
      have hl₀ : l₀ = l := Nat.le_antisymm (hmin l hlc) (hlc' l₀ h1)
      subst hl₀
      rw [List.getElem?_append_left (by omega), hl] at h3
      obtain rfl := Option.some.inj h3
      have hside : NbV g anc dirs crets v (!lowDir) M c lowDir :=
        h4.popV hlen hLcM (fun _ h => h)
      have hother : NbV g anc dirs crets v (!lowDir) M c (!lowDir) :=
        ⟨e, v, hQM, hj.symm, .inr (.inl ⟨rfl, rfl⟩)⟩
      cases lowDir <;> exact ⟨by assumption, by assumption⟩
    · exact ⟨(hc.nb x hx hxL hxc).1.popV hlen hLcM (fun _ h => h),
        (hc.nb x hx hxL hxc).2.popV hlen hLcM (fun _ h => h)⟩
  · exact ⟨e, c, hQM, hj, .inl (hLcM.subset hcL)⟩

/-! ### A completed block -/

/-- A completed block below `v` (boundary tree edge `(v, c)`, direction `false` at `v`): `v`
followed by the child's pieces is an st-order of the block. -/
theorem block_st {g : Graph} (hg : g.WF) {anc : List Nat} {dirs : List Bool} {v c : Nat}
    {vs rets : List Nat} {L : List ItemId} (hlen : dirs.length = anc.length) (hcv : c < g.nv)
    (hc : VInv g (anc ++ [v]) (dirs ++ [false]) c vs rets true L)
    (hrets : ∀ l ∈ rets, anc.length ≤ l) (hvL : ∀ x, x < g.nv → vertItem x ∈ L → x ≠ v) :
    StBlock.St g ⟨some (v, c), L⟩ := by
  have hcL : vertItem c ∈ L := (hc.root_nb rfl).1
  have hclause : ∀ {l y}, l ∈ rets → (anc ++ [v])[l]? = some y → y = v ∧ l = anc.length :=
    fun hl h => path_clause_boundary hrets h hl
  have hdirs : ∀ {l s}, l ∈ rets → (dirs ++ [false])[l]? = some s → s = false := by
    intro l s hl h
    have hl' : l = anc.length := by
      have := (List.getElem?_eq_some_iff.mp h).1
      simp only [List.length_append, List.length_singleton] at this
      have := hrets l hl; omega
    subst hl'
    rw [← hlen, List.getElem?_append_right (Nat.le_refl _)] at h
    simpa using h.symm
  have hedge : ∀ {e}, edgeItem g e ∈ L → e < g.ne := by
    intro e he
    rcases hc.shape _ he with ⟨y, hy, hxy⟩ | ⟨e', he', hxy⟩
    · exact absurd hxy.symm (vertItem_ne_edgeItem hy)
    · have : e = e' := by simp [edgeItem] at hxy; omega
      subst this; exact he'
  show Items.StList _ _
  rw [StBlock.seq_eq, StBlock.edges_eq]
  simp only [Option.map, Option.toList, List.singleton_append]
  set S := L.filterMap (vertOf g) with hS
  set E := L.filterMap (edgeOf g) with hE
  have hSmem : ∀ {x}, x ∈ S → vertItem x ∈ L ∧ x < g.nv := fun h => mem_filterMap_vertOf hc.shape h
  have hSmem' : ∀ {x}, x < g.nv → vertItem x ∈ L → x ∈ S := fun hx h => mem_filterMap_vertOf_of hx h
  have hvS : v ∉ S := fun h => hvL v (hSmem h).2 (hSmem h).1 rfl
  have hnd : (v :: S).Nodup := List.nodup_cons.mpr ⟨hvS, hc.vert_nodup⟩
  have hbef : ∀ {x y}, x < g.nv → y < g.nv → Before L (vertItem x) (vertItem y) →
      Before (v :: S) x y := fun hx hy h =>
    (h.filterMap (vertOf_vertItem hx) (vertOf_vertItem hy)).sub (List.sublist_cons_self _ _)
  have hvbef : ∀ {x}, x ∈ S → Before (v :: S) v x := fun h =>
    Before.of_mem_mem (List.mem_singleton_self v) h
  have hidx : ∀ {x y}, Before (v :: S) x y → (v :: S).idxOf x < (v :: S).idxOf y :=
    fun h => h.idxOf_lt hnd
  have hedgeE : ∀ {e}, edgeItem g e ∈ L → g.edges[e]! ∈ (v, c) :: E := fun he =>
    List.mem_cons_of_mem _ ((mem_filterMap_edgeOf hc.shape _).mpr ⟨_, hedge he, he, rfl⟩)
  have hpos : ∀ {x y s}, x < g.nv → y < g.nv → vertItem x ∈ L →
      (Side L s x y ∨ ∃ l, (anc ++ [v])[l]? = some y ∧ (dirs ++ [false])[l]? = some s ∧ l ∈ rets) →
      if s then Before (v :: S) x y else Before (v :: S) y x := by
    intro x y s hx hy hxL h
    rcases h with h | ⟨l, h1, h2, h3⟩
    · cases s
      · exact hbef hy hx h
      · exact hbef hx hy h
    · obtain ⟨rfl, -⟩ := hclause h3 h1
      obtain rfl := hdirs h3 h2
      exact hvbef (hSmem' hx hxL)
  have hclF : ∀ {x}, x < g.nv → vertItem x ∈ L →
      Nb g (anc ++ [v]) (dirs ++ [false]) rets L x false →
      ∃ p ∈ (v, c) :: E, (p.1 = x ∧ (v :: S).idxOf p.2 < (v :: S).idxOf x) ∨
        (p.2 = x ∧ (v :: S).idxOf p.1 < (v :: S).idxOf x) := by
    intro x hx hxL ⟨e, y, he, hj, h⟩
    have hy : y < g.nv := (hj.lt hg (hedge he)).2
    have hb : Before (v :: S) y x := hpos hx hy hxL h
    refine ⟨g.edges[e]!, hedgeE he, ?_⟩
    rcases hj with hj | hj <;> rw [hj]
    · exact .inl ⟨rfl, hidx hb⟩
    · exact .inr ⟨rfl, hidx hb⟩
  have hclT : ∀ {x}, x < g.nv → vertItem x ∈ L →
      Nb g (anc ++ [v]) (dirs ++ [false]) rets L x true →
      ∃ p ∈ (v, c) :: E, (p.1 = x ∧ (v :: S).idxOf x < (v :: S).idxOf p.2) ∨
        (p.2 = x ∧ (v :: S).idxOf x < (v :: S).idxOf p.1) := by
    intro x hx hxL ⟨e, y, he, hj, h⟩
    have hy : y < g.nv := (hj.lt hg (hedge he)).2
    have hb : Before (v :: S) x y := hpos hx hy hxL h
    refine ⟨g.edges[e]!, hedgeE he, ?_⟩
    rcases hj with hj | hj <;> rw [hj]
    · exact .inl ⟨rfl, hidx hb⟩
    · exact .inr ⟨rfl, hidx hb⟩
  refine ⟨hnd, ?_, ?_⟩
  · intro p hp
    rcases List.mem_cons.mp hp with rfl | hp
    · exact ⟨List.mem_cons_self .., List.mem_cons_of_mem _ (hSmem' hcv hcL), (hvL c hcv hcL).symm⟩
    · obtain ⟨e, he, heL, rfl⟩ := (mem_filterMap_edgeOf hc.shape p).mp hp
      obtain ⟨x, y, hj, hne, hxL, hy⟩ := hc.edge_mem e heL
      obtain ⟨hx, hy'⟩ := hj.lt hg he
      have hxS : x ∈ v :: S := List.mem_cons_of_mem _ (hSmem' hx hxL)
      have hyS : y ∈ v :: S := by
        rcases hy with hy | ⟨l, h1, h3⟩
        · exact List.mem_cons_of_mem _ (hSmem' hy' hy)
        · obtain ⟨rfl, -⟩ := hclause h3 h1
          exact List.mem_cons_self ..
      rcases hj with hj | hj <;> rw [hj]
      · exact ⟨hxS, hyS, hne⟩
      · exact ⟨hyS, hxS, hne.symm⟩
  · intro x hx hhead hlast
    rcases List.mem_cons.mp hx with rfl | hxS
    · simp at hhead
    obtain ⟨hxL, hxlt⟩ := hSmem hxS
    by_cases hxc : x = c
    · subst hxc
      refine absurd (getLast?_eq_of_before hnd hx fun y hy hyc => ?_) hlast
      rcases List.mem_cons.mp hy with rfl | hyS
      · exact hvbef hxS
      · obtain ⟨hyL, hylt⟩ := hSmem hyS
        obtain ⟨l, s, h1, h2, h3⟩ := hc.root_side y hylt hyL hyc
        obtain rfl := hdirs h1 h2
        exact hbef hylt hxlt h3
    · obtain ⟨nb0, nb1⟩ := hc.nb x hxlt hxL hxc
      exact ⟨hclF hxlt hxL nb0, hclT hxlt hxL nb1⟩

/-! ### `refOut` unfolded -/

theorem refOut_boundary_back {g : Graph} {v d : Nat} {dirs : List Bool} {e dest : Nat}
    {cls : OutClass} {hv : Bool} (h : d ≤ cls.lowval d) :
    refOut g v d dirs (.back e dest cls) hv = ([], [], hv) := by
  simp [refOut, DfsOut.cls, h]

theorem refOut_boundary_tree {g : Graph} {v d : Nat} {dirs : List Bool} {e : Nat} {cls : OutClass}
    {child : DfsTree} {hv : Bool} (h : d ≤ cls.lowval d) :
    refOut g v d dirs (.tree e cls child) hv =
      ([], (refTree g child (d + 1) (dirs ++ [false])).2 ++
        [⟨some (v, child.v), stNest (refTree g child (d + 1) (dirs ++ [false])).1⟩], hv) := by
  simp [refOut, DfsOut.cls, DfsOut.dest, h]

theorem refOut_ret_back {g : Graph} {v d : Nat} {dirs : List Bool} {e dest l : Nat} {hv : Bool}
    (h : l < d) :
    refOut g v d dirs (.back e dest (.ret l .backEdge)) hv =
      (if hv then [⟨dirs.getD l false, [edgeItem g e]⟩]
        else [⟨!dirs.getD l false, [vertItem v]⟩, ⟨dirs.getD l false, [edgeItem g e]⟩], [], true) := by
  cases hv <;> simp [refOut, DfsOut.cls, OutClass.lowval, OutClass.isType1, Nat.not_le.mpr h]

theorem refOut_ret_tree {g : Graph} {v d : Nat} {dirs : List Bool} {e l : Nat} {k : RetKind}
    {child : DfsTree} {hv : Bool} (h : l < d) :
    refOut g v d dirs (.tree e (.ret l k) child) hv =
      ((if hv then [⟨dirs.getD l false, stNest ((refTree g child (d + 1) (dirs ++ [!dirs.getD l false])).1 ++
            [⟨!dirs.getD l false, [edgeItem g e]⟩])⟩]
        else if k = .type2Child then
          (refTree g child (d + 1) (dirs ++ [!dirs.getD l false])).1 ++
            [⟨!dirs.getD l false, [edgeItem g e]⟩] ++ [⟨!dirs.getD l false, [vertItem v]⟩]
        else [⟨!dirs.getD l false, [vertItem v]⟩, ⟨dirs.getD l false, stNest
          ((refTree g child (d + 1) (dirs ++ [!dirs.getD l false])).1 ++
            [⟨!dirs.getD l false, [edgeItem g e]⟩])⟩]),
        (refTree g child (d + 1) (dirs ++ [!dirs.getD l false])).2, true) := by
  cases hv <;> cases k <;>
    simp [refOut, DfsOut.cls, OutClass.lowval, OutClass.isType1, Nat.not_le.mpr h]

theorem refOuts_nil {g : Graph} {v d : Nat} {dirs : List Bool} {hv : Bool} :
    refOuts g v d dirs [] hv = ([], [], hv) := by
  rw [refOuts]

theorem refOuts_cons {g : Graph} {v d : Nat} {dirs : List Bool} {o : DfsOut} {rest : List DfsOut}
    {hv : Bool} :
    refOuts g v d dirs (o :: rest) hv =
      ((refOut g v d dirs o hv).1 ++ (refOuts g v d dirs rest (refOut g v d dirs o hv).2.2).1,
        (refOut g v d dirs o hv).2.1 ++ (refOuts g v d dirs rest (refOut g v d dirs o hv).2.2).2.1,
        (refOuts g v d dirs rest (refOut g v d dirs o hv).2.2).2.2) := by
  rw [refOuts]

theorem refTree_node {g : Graph} {v : Nat} {outs : List DfsOut} {d : Nat} {dirs : List Bool} :
    refTree g (.node v outs) d dirs =
      (if (refOuts g v d dirs outs false).2.2 then (refOuts g v d dirs outs false).1
        else (refOuts g v d dirs outs false).1 ++ [StPiece.mk true [vertItem v]],
        (refOuts g v d dirs outs false).2.1) := by
  rw [refTree]

/-- At depth `0` every out-edge is a block boundary: nothing is pushed. -/
theorem refOuts_zero {g : Graph} {v : Nat} {dirs : List Bool} :
    ∀ (outs : List DfsOut) (hv : Bool),
      (refOuts g v 0 dirs outs hv).1 = [] ∧ (refOuts g v 0 dirs outs hv).2.2 = hv
  | [], hv => by simp [refOuts_nil]
  | o :: rest, hv => by
    rw [refOuts_cons]
    have h : (refOut g v 0 dirs o hv).1 = [] ∧ (refOut g v 0 dirs o hv).2.2 = hv := by
      cases o
      · rw [refOut_boundary_back (Nat.zero_le _)]; exact ⟨rfl, rfl⟩
      · rw [refOut_boundary_tree (Nat.zero_le _)]; exact ⟨rfl, rfl⟩
    rw [h.1, h.2]
    exact refOuts_zero rest hv

/-- A vertex with no returning out-edge: its own vertex item alone. -/
theorem VInv.single {g : Graph} {anc : List Nat} {dirs : List Bool} {v : Nat} {vs rets : List Nat}
    (hv : v < g.nv) (h : VInv g anc dirs v vs rets false []) :
    VInv g anc dirs v vs rets true [vertItem v] where
  shape := fun i hi => by simp at hi; exact .inl ⟨v, hv, hi⟩
  vert_mem := fun x _ hx => by simp at hx; exact .inl (vertItem_inj hx)
  vert_nodup := by simp [vertOf_vertItem hv]
  edge_mem := fun e he => by simp at he; exact absurd he.symm (vertItem_ne_edgeItem hv)
  nb := fun x _ hx hne => by simp at hx; exact absurd (vertItem_inj hx) hne
  root_side := fun y _ hy hne => by simp at hy; exact absurd (vertItem_inj hy) hne
  root_nb := fun _ => ⟨by simp, fun ⟨l, hl, hlt⟩ => absurd ((h.empty rfl).2 l hl) (by omega)⟩
  empty := fun h => nomatch h

/-- The block of a DFS root: just the root vertex. -/
theorem root_block_st (g : Graph) {v : Nat} (hv : v < g.nv) : StBlock.St g ⟨none, [vertItem v]⟩ := by
  show Items.StList _ _
  rw [StBlock.seq_eq, StBlock.edges_eq]
  simp [Option.map, Option.toList, vertOf_vertItem hv, edgeOf_vertItem g hv, Items.StList]

theorem rank_ret_le {l : Nat} {k : RetKind} {c : OutClass}
    (h : (OutClass.ret l k).rank ≤ c.rank) : ∃ l' k', c = .ret l' k' ∧ l ≤ l' := by
  cases c with
  | ret l' k' =>
    exact ⟨l', k', rfl, by cases k <;> cases k' <;> simp [OutClass.rank, RetKind.rank] at h <;> omega⟩
  | bridge => exfalso; cases k <;> simp only [OutClass.rank, RetKind.rank] at h <;> omega
  | component => exfalso; cases k <;> simp only [OutClass.rank, RetKind.rank] at h <;> omega
  | selfLoop => exfalso; cases k <;> simp only [OutClass.rank, RetKind.rank] at h <;> omega

theorem classify_back_lt {d i : Nat} (h : i < d) : classify d false (i, d) = .ret i .backEdge := by
  unfold classify; simp [Nat.not_le.mpr h]

theorem _root_.Spqr.DfsTree.verts_eq (t : DfsTree) : t.verts = t.v :: t.verts.tail := by
  cases t; rfl

theorem getElem?_eq_some_getD {dirs : List Bool} {l : Nat} (h : l < dirs.length) :
    dirs[l]? = some (dirs.getD l false) := by
  rw [List.getD_eq_getElem?_getD, List.getElem?_eq_getElem h]; rfl

theorem joins_of_Joins {g : Graph} {e x y : Nat} (h : g.Joins e x y) : e < g.ne ∧ Joins g e x y := by
  rcases h with h | h
  · obtain ⟨he, h'⟩ := Array.getElem?_eq_some_iff.mp h
    exact ⟨he, .inl (by rw [getElem!_pos g.edges e he]; exact h')⟩
  · obtain ⟨he, h'⟩ := Array.getElem?_eq_some_iff.mp h
    exact ⟨he, .inr (by rw [getElem!_pos g.edges e he]; exact h')⟩

/-- Every out-edge of the tree joins its vertex to its destination (`DfsData.Spec.joins`). -/
def TreeJoins (g : Graph) (t : DfsTree) : Prop := ∀ p ∈ t.allOuts, g.Joins p.2.e p.1 p.2.dest

def OutsJoins (g : Graph) (v : Nat) (outs : List DfsOut) : Prop :=
  ∀ p ∈ DfsOut.allOutsList v outs, g.Joins p.2.e p.1 p.2.dest

/-! ### The induction -/

mutual
/-- The pieces of a well-formed subtree at depth `anc.length` satisfy the invariant with
`hasVert = true`, and every block completed inside is st-numbered. -/
theorem refTree_inv {g : Graph} (hg : g.WF) : ∀ (t : DfsTree) (anc : List Nat) (dirs : List Bool),
    dirs.length = anc.length → t.WF anc → TreeJoins g t → (anc ++ t.verts).Nodup →
    (∀ x ∈ t.verts, x < g.nv) →
    (∀ b ∈ (refTree g t anc.length dirs).2, b.St g) ∧
    VInv g anc dirs t.v t.verts.tail (t.retDepths anc.length) true
      (stNest (refTree g t anc.length dirs).1)
  | .node v outs, anc, dirs, hlen, hwf, hJ, hnd, hB => by
    rw [DfsTree.WF] at hwf
    obtain ⟨hsort, hwf⟩ := hwf
    have h := refOuts_inv hg outs v anc dirs [] [] false [] hlen hsort hwf hJ
      (by simpa [DfsTree.verts] using hnd) (by simpa [DfsTree.verts] using hB)
      (fun ⟨_, h, _⟩ => absurd h List.not_mem_nil)
      (by simpa [stNest, stNestL, stNestR] using VInv.nil g anc dirs v)
    simp only [refTree_node, DfsTree.v, DfsTree.verts, List.tail_cons, DfsTree.retDepths]
    refine ⟨h.1, ?_⟩
    have h2 := h.2
    simp only [List.nil_append] at h2
    cases hv : (refOuts g v anc.length dirs outs false).2.2
    · rw [hv] at h2
      have hP : stNest (refOuts g v anc.length dirs outs false).1 = [] := (h2.empty rfl).1
      simp only [Bool.false_eq_true, ite_false]
      rw [stNest_append, hP]
      rw [hP] at h2
      simpa [stNest, stNestL, stNestR] using h2.single (hB v (by simp [DfsTree.verts]))
    · rw [hv] at h2
      simpa using h2
theorem refOuts_inv {g : Graph} (hg : g.WF) : ∀ (outs : List DfsOut) (v : Nat) (anc : List Nat)
    (dirs : List Bool) (vs rets : List Nat) (hv : Bool) (P : List StPiece),
    dirs.length = anc.length → outs.Pairwise (fun a b => a.cls.rank ≤ b.cls.rank) →
    (∀ o ∈ outs, o.WF anc v) → OutsJoins g v outs →
    (anc ++ v :: (vs ++ DfsOut.vertsList outs)).Nodup →
    (∀ x ∈ v :: DfsOut.vertsList outs, x < g.nv) →
    ((∃ l₀ ∈ rets, l₀ < anc.length) →
      ∀ o ∈ outs, ∃ l k, o.cls = .ret l k ∧ ∃ l₁ ∈ rets, l₁ ≤ l) →
    VInv g anc dirs v vs rets hv (stNest P) →
    (∀ b ∈ (refOuts g v anc.length dirs outs hv).2.1, b.St g) ∧
    VInv g anc dirs v (vs ++ DfsOut.vertsList outs) (rets ++ DfsOut.retDepthsList anc.length outs)
      (refOuts g v anc.length dirs outs hv).2.2
      (stNest (P ++ (refOuts g v anc.length dirs outs hv).1))
  | [], v, anc, dirs, vs, rets, hv, P, _, _, _, _, _, _, _, hinv => by
    simp only [refOuts_nil, DfsOut.vertsList, DfsOut.retDepthsList, List.append_nil]
    exact ⟨fun _ h => absurd h List.not_mem_nil, hinv⟩
  | o :: rest, v, anc, dirs, vs, rets, hv, P, hlen, hsort, hwf, hJ, hnd, hB, hrs, hinv => by
    have hvlt : v < g.nv := hB v (List.mem_cons_self ..)
    have hvanc : v ∉ anc := fun h => List.disjoint_of_nodup_append hnd h (List.mem_cons_self ..)
    have hwfo := hwf o (List.mem_cons_self ..)
    have hdirs : ∀ {l : Nat}, l < anc.length → dirs[l]? = some (dirs.getD l false) :=
      fun h => getElem?_eq_some_getD (by omega)
    have key : ∃ cvs, DfsOut.vertsList (o :: rest) = cvs ++ DfsOut.vertsList rest ∧
        (∀ b ∈ (refOut g v anc.length dirs o hv).2.1, b.St g) ∧
        VInv g anc dirs v (vs ++ cvs) (rets ++ o.retDepths anc.length)
          (refOut g v anc.length dirs o hv).2.2
          (stNest (P ++ (refOut g v anc.length dirs o hv).1)) := by
      cases o with
      | back e dest cls =>
        have hJ' : g.Joins e v dest := hJ (v, .back e dest cls) (by simp [DfsOut.allOutsList])
        obtain ⟨hJe, hJj⟩ := joins_of_Joins hJ'
        rw [DfsOut.WF] at hwfo
        obtain ⟨i, hi, hcls⟩ := hwfo
        refine ⟨[], by simp [DfsOut.vertsList], ?_⟩
        by_cases hb : anc.length ≤ cls.lowval anc.length
        · rw [refOut_boundary_back hb]
          refine ⟨fun _ h => absurd h List.not_mem_nil, ?_⟩
          simp only [List.append_nil, DfsOut.retDepths]
          exact hinv.extend hlen (fun x hx => hx)
            (fun l hl => List.mem_append_left _ hl)
            (fun l hl => by
              rcases List.mem_append.mp hl with hl | hl
              · exact .inl hl
              · simp at hl; subst hl; exact .inr hb)
        · have hi' : i < anc.length + 1 := by
            have := (List.getElem?_eq_some_iff.mp hi).1; simpa using this
          have hlow : cls.lowval anc.length = i := by
            rw [hcls]; exact lowval_classify_back (by omega)
          have hilt : i < anc.length := by omega
          rw [classify_back_lt hilt] at hcls
          subst hcls
          rw [List.getElem?_append_left hilt] at hi
          have hne : dest ≠ v := fun h => hvanc (h ▸ List.mem_of_getElem? hi)
          rw [refOut_ret_back hilt]
          refine ⟨fun _ h => absurd h List.not_mem_nil, ?_⟩
          simp only [DfsOut.retDepths, OutClass.lowval]
          refine hinv.step_back hlen hvlt (hdirs hilt) hi hne hJe hJj ?_ ?_
          · cases hv
            · have hP := (hinv.empty rfl).1
              rw [stNest_append, hP]
              cases dirs.getD i false <;> simp [stNestL, stNestR]
            · rw [stNest_append]
              generalize stNest P = L0
              cases dirs.getD i false <;> simp [stNestL, stNestR]
          · intro hex
            obtain ⟨l, k, h1, l₁, h2, h3⟩ := hrs hex _ (List.mem_cons_self ..)
            simp only [DfsOut.cls, OutClass.ret.injEq] at h1
            exact ⟨l₁, h2, h1.1 ▸ h3⟩
      | tree e cls child =>
        have hJ' : g.Joins e v child.v := hJ (v, .tree e cls child) (by simp [DfsOut.allOutsList])
        obtain ⟨hJe, hJj⟩ := joins_of_Joins hJ'
        rw [DfsOut.WF] at hwfo
        obtain ⟨hcwf, hcls⟩ := hwfo
        have hcJ : TreeJoins g child := fun p hp => hJ p (by simp [DfsOut.allOutsList, hp])
        have hcB : ∀ x ∈ child.verts, x < g.nv := fun x hx => hB x (by simp [DfsOut.vertsList, hx])
        have hnd2 : (anc ++ v :: (vs ++ (child.verts ++ DfsOut.vertsList rest))).Nodup := by
          simpa [DfsOut.vertsList] using hnd
        have hcnd : (anc ++ [v] ++ child.verts).Nodup := by
          refine hnd2.sublist ?_
          rw [List.append_assoc, List.singleton_append]
          exact (List.Sublist.refl anc).append (List.Sublist.cons_cons v
            ((List.sublist_append_left _ _).trans (List.sublist_append_right _ _)))
        have hcdisj : ∀ x ∈ child.verts, x ≠ v ∧ x ∉ vs := by
          intro x hx
          have h1 := hnd2.of_append_right
          rw [List.nodup_cons] at h1
          exact ⟨fun h => h1.1 (h ▸ List.mem_append_right _ (List.mem_append_left _ hx)),
            fun h => List.disjoint_of_nodup_append h1.2 h (List.mem_append_left _ hx)⟩
        have hcverts : child.verts = child.v :: child.verts.tail := child.verts_eq
        have hwfo' := hwf _ (List.mem_cons_self ..)
        refine ⟨child.verts, by simp [DfsOut.vertsList], ?_⟩
        by_cases hb : anc.length ≤ cls.lowval anc.length
        · have ih := refTree_inv hg child (anc ++ [v]) (dirs ++ [false]) (by simp [hlen]) hcwf hcJ
            hcnd hcB
          simp only [List.length_append, List.length_singleton] at ih
          rw [refOut_boundary_tree hb]
          refine ⟨fun b hb' => ?_, ?_⟩
          · rcases List.mem_append.mp hb' with hb' | hb'
            · exact ih.1 b hb'
            · simp at hb'; subst hb'
              refine block_st hg hlen (hcB _ (by rw [hcverts]; exact List.mem_cons_self ..)) ih.2
                ?_ ?_
              · intro l hl
                have : cls.lowval anc.length ≤ l := DfsOut.lowval_le_retDepths hwfo' l hl
                omega
              · intro x hx hxL
                rcases ih.2.vert_mem x hx hxL with h | h
                · exact (hcdisj x (by rw [hcverts, h]; exact List.mem_cons_self ..)).1
                · exact (hcdisj x (by rw [hcverts]; exact List.mem_cons_of_mem _ h)).1
          · simp only [List.append_nil, DfsOut.retDepths]
            exact hinv.extend hlen (fun x hx => List.mem_append_left _ hx)
              (fun l hl => List.mem_append_left _ hl)
              (fun l hl => by
                rcases List.mem_append.mp hl with hl | hl
                · exact .inl hl
                · have : cls.lowval anc.length ≤ l := DfsOut.lowval_le_retDepths hwfo' l hl
                  exact .inr (by omega))
        · have hret : cls.lowval anc.length < anc.length := Nat.lt_of_not_le hb
          obtain ⟨l, k, hcls'⟩ := DfsOut.ret_of_lowval_lt hwfo' hret
          simp only [DfsOut.cls] at hcls'
          subst hcls'
          have hlt : l < anc.length := hret
          obtain ⟨lowDir, hlowDir⟩ : ∃ b, dirs.getD l false = b := ⟨_, rfl⟩
          have hl : dirs[l]? = some lowDir := hlowDir ▸ hdirs hlt
          have hlc : l ∈ child.retDepths (anc.length + 1) := DfsOut.lowval_mem_retDepths hwfo' hret
          have hlc' : ∀ l' ∈ child.retDepths (anc.length + 1), l ≤ l' :=
            DfsOut.lowval_le_retDepths hwfo'
          have ih := refTree_inv hg child (anc ++ [v]) (dirs ++ [!lowDir]) (by simp [hlen]) hcwf hcJ
            hcnd hcB
          simp only [List.length_append, List.length_singleton] at ih
          rw [refOut_ret_tree hlt, hlowDir]
          refine ⟨ih.1, ?_⟩
          simp only [DfsOut.retDepths]
          rw [hcverts]
          refine hinv.step_tree hlen hvlt hl ih.2 (fun x hx => hcdisj x (by rw [hcverts]; exact hx))
            hJe hJj hlc hlc' ?_ ?_
          · generalize (refTree g child (anc.length + 1) (dirs ++ [!lowDir])).1 = ps
            cases hv
            · have hP := (hinv.empty rfl).1
              rw [stNest_append, hP]
              cases k <;> cases lowDir <;> simp [stNestL_append, stNestR_append, stNest, stNestL, stNestR]
            · rw [stNest_append]
              generalize stNest P = L0
              cases lowDir <;> simp [stNestL_append, stNestR_append, stNest, stNestL, stNestR]
          · intro hex
            obtain ⟨l', k', h1, l₁, h2, h3⟩ := hrs hex _ (List.mem_cons_self ..)
            simp only [DfsOut.cls, OutClass.ret.injEq] at h1
            exact ⟨l₁, h2, h1.1 ▸ h3⟩
    obtain ⟨cvs, hcv, hbl, hinv'⟩ := key
    have hcr : DfsOut.retDepthsList anc.length (o :: rest) =
        o.retDepths anc.length ++ DfsOut.retDepthsList anc.length rest := by
      cases o <;> simp [DfsOut.retDepths]
    have hrs' : (∃ l₀ ∈ rets ++ o.retDepths anc.length, l₀ < anc.length) →
        ∀ o' ∈ rest, ∃ l k, o'.cls = .ret l k ∧ ∃ l₁ ∈ rets ++ o.retDepths anc.length, l₁ ≤ l := by
      rintro ⟨l₀, hl₀, hlt⟩ o' ho'
      rcases List.mem_append.mp hl₀ with hl₀ | hl₀
      · obtain ⟨l, k, h1, l₁, h2, h3⟩ := hrs ⟨l₀, hl₀, hlt⟩ o' (List.mem_cons_of_mem _ ho')
        exact ⟨l, k, h1, l₁, List.mem_append_left _ h2, h3⟩
      · have hret : o.cls.lowval anc.length < anc.length := by
          rw [DfsOut.lowval_eq_lmin hwfo]; exact Nat.lt_of_le_of_lt (lmin_le_of_mem hl₀) hlt
        obtain ⟨l, k, hcls⟩ := DfsOut.ret_of_lowval_lt hwfo hret
        have hrank := List.rel_of_pairwise_cons hsort ho'
        rw [hcls] at hrank
        obtain ⟨l', k', h1, h2⟩ := rank_ret_le hrank
        refine ⟨l', k', h1, l, List.mem_append_right _ ?_, h2⟩
        have := DfsOut.lowval_mem_retDepths hwfo hret
        rwa [hcls] at this
    have ih := refOuts_inv hg rest v anc dirs (vs ++ cvs) (rets ++ o.retDepths anc.length) _
      (P ++ (refOut g v anc.length dirs o hv).1) hlen hsort.of_cons
      (fun o' ho' => hwf o' (List.mem_cons_of_mem _ ho'))
      (fun p hp => hJ p (by
        cases o <;> simp only [DfsOut.allOutsList]
        · exact List.mem_cons_of_mem _ hp
        · exact List.mem_cons_of_mem _ (List.mem_append_right _ hp)))
      (by rw [hcv] at hnd; simpa [List.append_assoc] using hnd)
      (fun x hx => hB x (by
        rw [hcv]
        rcases List.mem_cons.mp hx with rfl | hx
        · exact List.mem_cons_self ..
        · exact List.mem_cons_of_mem _ (List.mem_append_right _ hx))) hrs' hinv'
    rw [refOuts_cons]
    refine ⟨fun b hb => ?_, ?_⟩
    · rcases List.mem_append.mp hb with hb | hb
      · exact hbl b hb
      · exact ih.1 b hb
    · rw [hcv, hcr]
      simpa [List.append_assoc] using ih.2
end

/-- Even–Tarjan on the reference: every block of `refBlocks` is st-numbered by its sequence
(`refTree_inv`: every piece spliced at depth `l` joins the open path at depth `l` on side
`dirs[l]`, so every vertex other than the block's terminals has a neighbour on each side). The
DFS facts used (`DfsTree.WF`, `DfsForestSpec.joins`, distinct vertices below `g.nv`) need `g.WF`
and the order hypotheses, as `dfsForest_spanning` does; without `g.WF` the statement fails
(`PROOF.md` §7.6). -/
theorem refBlocks_st {g : Graph} (hg : g.WF) {vo eo : List Nat} (hvo : OrderOK g.nv vo)
    (heo : OrderOK g.ne eo) : ∀ b ∈ refBlocks g (g.dfsForest vo eo), b.St g := by
  intro b hb
  have spec := dfsForestSpec_of_dfsForest hg hvo heo
  have hwf := dfsForest_wf hg hvo heo
  have hnd : ((g.dfsForest vo eo).flatMap DfsTree.verts).Nodup := spec.verts_nodup
  have hB : ∀ x ∈ (g.dfsForest vo eo).flatMap DfsTree.verts, x < g.nv := fun x hx =>
    List.mem_range.mp ((dfsForest_spanning' hg hvo heo).1.subset hx)
  unfold refBlocks at hb
  obtain ⟨t, ht, hb⟩ := List.mem_flatMap.mp hb
  have hJ : TreeJoins g t := fun p hp => by
    obtain ⟨x, o⟩ := p
    exact spec.joins x o (DfsData.mem_forestAllOuts.mp (List.mem_flatMap.mpr ⟨t, ht, hp⟩))
  have h := refTree_inv hg t [] [] rfl (hwf t ht) hJ
    (by simpa using (List.nodup_flatMap.mp hnd).1 t ht)
    (fun x hx => hB x (List.mem_flatMap.mpr ⟨t, ht, hx⟩))
  simp only [List.length_nil] at h
  rcases List.mem_append.mp hb with hb | hb
  · exact h.1 b hb
  · simp at hb; subst hb
    obtain ⟨v, outs⟩ := t
    rw [refTree_node, (refOuts_zero outs false).1, (refOuts_zero outs false).2]
    simp only [Bool.false_eq_true, ite_false, List.nil_append]
    rw [stNest_single]
    exact root_block_st g (hB v (List.mem_flatMap.mpr ⟨_, ht, by simp [DfsTree.verts]⟩))

end Spqr.StRefEt

namespace Spqr
export StRefEt (refBlocks_st)
end Spqr
