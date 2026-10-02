import Spqr.StRef

/-!
# The simulation relation for `walk_st'` / `walk_vsOriented` (PROOF.md §7.6)

The walk's tstack above the base of the current frame reads (`readStack`) as `stNest` of the
reference's pieces for that frame, up to expanding the S / P / R items closed inside the ear back to
their leaves (`Expands`); every S / P / R item is either live (below an item of the open stack) or
finished, i.e. a contiguous segment of a reference block with the `VsOriented` clauses already
satisfied; the blocks used are those of the *truncated* DFS tree (`truncTree`: along the open path
every vertex keeps its finished out-edges and the current tree edge), whose st-orders only grow
outwards as the walk proceeds. Checked at every `finishEdge` boundary by `check_stsim`.
-/

namespace Spqr

/-! ### Leaf expansion -/

/-- `ExpandsList items xs L`: `L` is the concatenated leaf expansion of `xs` (V / Q items are their
own leaf, any other item is replaced by its children). -/
inductive ExpandsList (items : Items) : List ItemId → List ItemId → Prop
  | nil : ExpandsList items [] []
  | leaf {x : ItemId} {xs L : List ItemId} (h : Items.type items x = .V ∨ Items.type items x = .Q)
      (hxs : ExpandsList items xs L) : ExpandsList items (x :: xs) (x :: L)
  | node {x : ItemId} {xs L : List ItemId} (h : ¬ (Items.type items x = .V ∨ Items.type items x = .Q))
      (hxs : ExpandsList items (Items.ch items x ++ xs) L) : ExpandsList items (x :: xs) L

/-- The leaves of one item. -/
def Expands (items : Items) (x : ItemId) (L : List ItemId) : Prop := ExpandsList items [x] L

namespace ExpandsList

variable {items : Items}

theorem append {xs ys L M : List ItemId} (h₁ : ExpandsList items xs L) (h₂ : ExpandsList items ys M) :
    ExpandsList items (xs ++ ys) (L ++ M) := by
  induction h₁ with
  | nil => simpa using h₂
  | leaf h _ ih => exact .leaf h ih
  | node h _ ih =>
    rw [List.append_assoc] at ih
    exact .node h ih

theorem append_inv {xs ys L : List ItemId} (h : ExpandsList items (xs ++ ys) L) :
    ∃ L₁ L₂, L = L₁ ++ L₂ ∧ ExpandsList items xs L₁ ∧ ExpandsList items ys L₂ := by
  generalize hz : xs ++ ys = zs at h
  induction h generalizing xs with
  | nil =>
    rcases List.append_eq_nil_iff.1 hz with ⟨rfl, rfl⟩
    exact ⟨[], [], rfl, .nil, .nil⟩
  | @leaf x zs L h hzs ih =>
    cases xs with
    | nil =>
      simp only [List.nil_append] at hz
      subst hz
      exact ⟨[], x :: L, rfl, .nil, .leaf h hzs⟩
    | cons x' xs =>
      simp only [List.cons_append, List.cons.injEq] at hz
      obtain ⟨rfl, hz⟩ := hz
      obtain ⟨L₁, L₂, rfl, h₁, h₂⟩ := ih hz
      exact ⟨_ :: L₁, L₂, rfl, .leaf h h₁, h₂⟩
  | @node x zs L h hzs ih =>
    cases xs with
    | nil =>
      simp only [List.nil_append] at hz
      subst hz
      exact ⟨[], L, rfl, .nil, .node h hzs⟩
    | cons x' xs =>
      simp only [List.cons_append, List.cons.injEq] at hz
      obtain ⟨rfl, hz⟩ := hz
      obtain ⟨L₁, L₂, rfl, h₁, h₂⟩ := ih (by rw [← hz, List.append_assoc])
      exact ⟨L₁, L₂, rfl, .node h h₁, h₂⟩

theorem cons_iff {x : ItemId} {xs L : List ItemId} :
    ExpandsList items (x :: xs) L ↔ ∃ L₁ L₂, L = L₁ ++ L₂ ∧ Expands items x L₁ ∧ ExpandsList items xs L₂ :=
  ⟨fun h => append_inv (xs := [x]) h, fun ⟨_, _, h, h₁, h₂⟩ => h ▸ append h₁ h₂⟩

theorem det {xs L L' : List ItemId} (h : ExpandsList items xs L) (h' : ExpandsList items xs L') : L = L' := by
  induction h generalizing L' with
  | nil => cases h'; rfl
  | leaf h _ ih =>
    cases h' with
    | leaf _ h' => rw [ih h']
    | node h' => exact absurd h h'
  | node h _ ih =>
    cases h' with
    | leaf h' => exact absurd h' h
    | node _ h' => exact ih h'

/-- Framing: the expansion only looks at the types and children below `xs`. -/
theorem congr {items' : Items} {xs L : List ItemId}
    (H : ∀ x ∈ xs, ∀ y, Items.Below items x y →
      Items.type items' y = Items.type items y ∧ Items.ch items' y = Items.ch items y)
    (h : ExpandsList items xs L) : ExpandsList items' xs L := by
  induction h with
  | nil => exact .nil
  | @leaf x xs L h _ ih =>
    have hx := H x (by simp) x .refl
    exact .leaf (hx.1 ▸ h) (ih fun x' hx' => H x' (by simp [hx']))
  | @node x xs L h _ ih =>
    have hx := H x (by simp) x .refl
    refine .node (hx.1 ▸ h) (hx.2 ▸ ih fun x' hx' y hy => ?_)
    rcases List.mem_append.1 hx' with hx' | hx'
    · exact H x (by simp) y (.head hx' hy)
    · exact H x' (by simp [hx']) y hy

/-- Expanding `i` (not a leaf) to its own children in the list does not change the expansion. -/
theorem expandItem_self_iff {i : ItemId} (hi : ¬ (Items.type items i = .V ∨ Items.type items i = .Q))
    {l M : List ItemId} :
    ExpandsList items (expandItem i (Items.ch items i) l) M ↔ ExpandsList items l M := by
  induction l generalizing M with
  | nil => simp [expandItem]
  | cons x l ih =>
    by_cases hx : x = i
    · subst hx
      simp only [expandItem, List.flatMap_cons, ite_true]
      constructor
      · intro h
        obtain ⟨L₁, L₂, rfl, h₁, h₂⟩ := append_inv h
        exact .node hi (append h₁ ((ih (M := L₂)).1 h₂))
      · intro h
        cases h with
        | leaf h => exact absurd h hi
        | node _ h =>
          obtain ⟨L₁, L₂, rfl, h₁, h₂⟩ := append_inv h
          exact append h₁ ((ih (M := L₂)).2 h₂)
    · simp only [expandItem, List.flatMap_cons, ite_eq_right_iff.2 (fun h => absurd h hx), List.singleton_append]
      change ExpandsList items (x :: expandItem i (Items.ch items i) l) M ↔ _
      rw [cons_iff, cons_iff]
      exact ⟨fun ⟨L₁, L₂, e, h₁, h₂⟩ => ⟨L₁, L₂, e, h₁, (ih (M := L₂)).1 h₂⟩,
        fun ⟨L₁, L₂, e, h₁, h₂⟩ => ⟨L₁, L₂, e, h₁, (ih (M := L₂)).2 h₂⟩⟩

end ExpandsList

/-- A root is below no other item. -/
theorem Items.not_below_root_of_ne {items : Items} {i x : ItemId} (hr : ∀ p, ¬ Items.IsParent items p i)
    (hx : x ≠ i) : ¬ Items.Below items x i := by
  intro h
  rcases Relation.ReflTransGen.cases_tail h with h | ⟨p, _, hp⟩
  · exact hx h.symm
  · exact hr p hp

/-- Modifying `i` does not change the expansion of lists `i` is not below. -/
theorem ExpandsList.modify_of_not_below {items : Items} {i : ItemId} (f : Item → Item)
    {xs L : List ItemId} (hi : ∀ x ∈ xs, ¬ Items.Below items x i)
    (h : ExpandsList items xs L) : ExpandsList (items.modify i f) xs L :=
  h.congr fun x hx y hy => by
    have hyi : y ≠ i := fun e => hi x hx (e ▸ hy)
    exact ⟨by simp [Items.type, Array.getElem?_modify, hyi.symm], Items.ch_modify_ne _ _ _ _ (Ne.symm hyi)⟩

theorem ExpandsList.modify_root {items : Items} {i : ItemId} (f : Item → Item)
    (hr : ∀ p, ¬ Items.IsParent items p i) {xs L : List ItemId} (hi : i ∉ xs)
    (h : ExpandsList items xs L) : ExpandsList (items.modify i f) xs L :=
  h.modify_of_not_below f fun x hx => Items.not_below_root_of_ne hr fun e => hi (e ▸ hx)

/-- Closing: `i` (not a leaf, below nothing in `l`) gets the children `ch`; a list `l'` that expands
through `i` to `l` reads as `l` did. -/
theorem ExpandsList.close {items : Items} {i : ItemId} (f : Item → Item)
    (hf : ∀ it, (f it).type = it.type) (hi : i < items.size)
    (hty : ¬ (Items.type items i = .V ∨ Items.type items i = .Q))
    {ch : List ItemId} (hch : (f items[i]).ch = ch)
    {l l' M : List ItemId} (hl : expandItem i ch l' = l) (hnb : ∀ x ∈ l, ¬ Items.Below items x i)
    (h : ExpandsList items l M) : ExpandsList (items.modify i f) l' M := by
  have hty' : ¬ (Items.type (items.modify i f) i = .V ∨ Items.type (items.modify i f) i = .Q) := by
    rwa [Items.type_modify _ _ _ _ hf]
  refine (ExpandsList.expandItem_self_iff hty').1 ?_
  rw [Items.ch_modify_self _ _ _ hi, hch, hl]
  exact h.modify_of_not_below f hnb

/-- Pushing a new item does not change the expansion of lists whose subtrees are in range. -/
theorem ExpandsList.push {items : Items} (it : Item) {xs L : List ItemId}
    (hb : ∀ x ∈ xs, ∀ y, Items.Below items x y → y < items.size)
    (h : ExpandsList items xs L) : ExpandsList (items.push it) xs L :=
  h.congr fun x hx y hy => by
    have := Nat.ne_of_lt (hb x hx y hy)
    rw [Items.type_push, Items.ch_push, ite_eq_right_iff.2 (fun e => absurd e this), ite_eq_right_iff.2 (fun e => absurd e this)]
    exact ⟨rfl, rfl⟩

/-- The reference's `dirs` are the walk's `stackDir` below `d`. -/
def DirsOf (s : WalkState) (d : Nat) : List Bool := (List.range d).map fun k => s.stackDir[k]!

/-- The stack segment `new` reads, up to expansion, as the pieces `ps`. -/
def StRead (items : Items) (new : List TEntry) (ps : List StPiece) : Prop :=
  ExpandsList items (readL new) (stNestL ps) ∧ ExpandsList items (readR new) (stNestR ps)

theorem StRead.flat {items : Items} {new : List TEntry} {ps : List StPiece} (h : StRead items new ps) :
    ExpandsList items (readStack new) (stNest ps) :=
  h.1.append h.2

theorem StRead.congr {items items' : Items} {new : List TEntry} {ps : List StPiece}
    (H : ∀ x ∈ readStack new, ∀ y, Items.Below items x y →
      Items.type items' y = Items.type items y ∧ Items.ch items' y = Items.ch items y)
    (h : StRead items new ps) : StRead items' new ps :=
  ⟨h.1.congr fun x hx => H x (List.mem_append_left _ hx),
    h.2.congr fun x hx => H x (List.mem_append_right _ hx)⟩

theorem mem_readStack_of_readL {x : ItemId} {ts : List TEntry} (h : x ∈ readL ts) : x ∈ readStack ts :=
  List.mem_append_left _ h

theorem mem_readStack_of_readR {x : ItemId} {ts : List TEntry} (h : x ∈ readR ts) : x ∈ readStack ts :=
  List.mem_append_right _ h

theorem readL_cons_setSides (v d idx : Nat) (dir : Bool) (items : List ItemId) (ts : List TEntry) :
    readL (⟨v, d, idx, setSides dir items []⟩ :: ts) = stNestL [⟨dir, items⟩] ++ readL ts := by
  cases dir <;> simp [readL, setSides, stNestL]

theorem readR_cons_setSides (v d idx : Nat) (dir : Bool) (items : List ItemId) (ts : List TEntry) :
    readR (⟨v, d, idx, setSides dir items []⟩ :: ts) = readR ts ++ stNestR [⟨dir, items⟩] := by
  cases dir <;> simp [readR, setSides, stNestR]

theorem expandItem_cons_self (i : ItemId) (ch l : List ItemId) :
    expandItem i ch (i :: l) = ch ++ expandItem i ch l := by
  simp [expandItem]

theorem expandItem_nil (i : ItemId) (ch : List ItemId) : expandItem i ch [] = [] := rfl

theorem readL_close (dir : Bool) (item : ItemId) (t : TEntry) (rest : List TEntry)
    (hside : getSide t.spans (!dir) = []) (hnew : item ∉ readL rest) :
    expandItem item (getSide t.spans dir) (readL ({ t with spans := setSides dir [item] [] } :: rest)) =
      readL (t :: rest) := by
  cases dir <;> simp [getSide] at hside ⊢ <;>
    simp [readL, setSides, expandItem_append, expandItem_cons_self, expandItem_nil, expandItem_of_not_mem _ _ _ hnew, hside]

theorem readR_close (dir : Bool) (item : ItemId) (t : TEntry) (rest : List TEntry)
    (hside : getSide t.spans (!dir) = []) (hnew : item ∉ readR rest) :
    expandItem item (getSide t.spans dir) (readR ({ t with spans := setSides dir [item] [] } :: rest)) =
      readR (t :: rest) := by
  cases dir <;> simp [getSide] at hside ⊢ <;>
    simp [readR, setSides, expandItem_append, expandItem_cons_self, expandItem_nil, expandItem_of_not_mem _ _ _ hnew, hside]

theorem readL_reopen (a b : TEntry) (rest : List TEntry) (dir : Bool) (i : ItemId) (ch : List ItemId)
    (hb : b.spans = setSides dir [i] []) (ha : i ∉ a.spans.1 ++ a.spans.2) (hrest : i ∉ readStack rest) :
    readL (a :: { b with spans := setSides dir ch [] } :: rest) = expandItem i ch (readL (a :: b :: rest)) := by
  have h1 : i ∉ a.spans.1 := fun h => ha (List.mem_append_left _ h)
  have h2 : i ∉ readL rest := fun h => hrest (mem_readStack_of_readL h)
  cases dir <;> simp [setSides] at hb <;>
    simp [readL, setSides, hb, expandItem_append, expandItem_cons_self, expandItem_of_not_mem _ _ _ h1,
      expandItem_of_not_mem _ _ _ h2]

theorem readR_reopen (a b : TEntry) (rest : List TEntry) (dir : Bool) (i : ItemId) (ch : List ItemId)
    (hb : b.spans = setSides dir [i] []) (ha : i ∉ a.spans.1 ++ a.spans.2) (hrest : i ∉ readStack rest) :
    readR (a :: { b with spans := setSides dir ch [] } :: rest) = expandItem i ch (readR (a :: b :: rest)) := by
  have h1 : i ∉ a.spans.2 := fun h => ha (List.mem_append_right _ h)
  have h2 : i ∉ readR rest := fun h => hrest (mem_readStack_of_readR h)
  cases dir <;> simp [setSides] at hb <;>
    simp [readR, setSides, hb, expandItem_append, expandItem_cons_self, expandItem_of_not_mem _ _ _ h1,
      expandItem_of_not_mem _ _ _ h2]

/-! ### Items against blocks -/

/-- The body of `VsOriented` for one item `i` with leaves `L` and one block `b`. -/
def VsOrientedAt (g : Graph) (items : Items) (b : StBlock) (i : ItemId) (L : List ItemId) : Prop :=
  (∀ x ∈ L, x ∈ b.items) ∧
  (∀ e, e < g.ne → (edgeItem g e ∈ b.items ∨ ∃ r ∈ b.root, Items.PairEq g.edges[e]! r) →
    Items.Below items i (edgeItem g e) → edgeItem g e ∈ L) ∧
  Oriented (b.seq g) (Items.vs items i) ∧
  (∀ s t, Items.vs items i = (some s, some t) → ∀ c ∈ Items.ch items i, Items.type items c = .V →
    Precedes (b.seq g) s (c - 1) ∧ Precedes (b.seq g) (c - 1) t) ∧
  (∀ c ∈ Items.ch items i, Items.type items c ≠ .V → Oriented (b.seq g) (Items.vs items c)) ∧
  ∀ c ∈ Items.ch items i, Items.type items c ≠ .V → ∀ u v, Items.vs items c = (some u, some v) →
    ∀ e, e < g.ne → (edgeItem g e ∈ b.items ∨ ∃ r ∈ b.root, Items.PairEq g.edges[e]! r) →
    Items.Below items c (edgeItem g e) →
    ∀ y, (g.edges[e]!).1 = y ∨ (g.edges[e]!).2 = y →
      (u = y ∨ Precedes (b.seq g) u y) ∧ (y = v ∨ Precedes (b.seq g) y v)

theorem vsOriented_iff (g : Graph) (items : Items) (blocks : List StBlock) :
    VsOriented g items blocks ↔
      ∀ i, i < items.size → Items.type items i = .S ∨ Items.type items i = .P ∨ Items.type items i = .R →
        ∃ b ∈ blocks, VsOrientedAt g items b i (Items.leaves items items.size i) := Iff.rfl

/-- `i` is finished in the block `b`: its leaves `L` are a contiguous segment of `b.items` and the
`VsOriented` clauses hold for them. -/
def InBlock (g : Graph) (items : Items) (b : StBlock) (i : ItemId) : Prop :=
  ∃ L, Expands items i L ∧ (∃ A B, b.items = A ++ L ++ B) ∧ VsOrientedAt g items b i L

/-- The items relative to the stack and the blocks: stack items are roots and distinct; every
S / P / R item is live (below a stack item) or finished in one of `blocks`; the stack's subtrees and
all children are in range, child lists have no repeats. -/
structure StItems (g : Graph) (s : WalkState) (blocks : List StBlock) : Prop where
  roots : ∀ x ∈ readStack s.tstack, ∀ p, ¬ Items.IsParent s.items p x
  nodup : (readStack s.tstack).Nodup
  bounded : ∀ x ∈ readStack s.tstack, ∀ y, Items.Below s.items x y → y < s.items.size
  chLt : ∀ p c, Items.IsParent s.items p c → c < s.items.size
  chNodup : ∀ p, (Items.ch s.items p).Nodup
  closed : ∀ i, i < s.items.size →
    Items.type s.items i = .S ∨ Items.type s.items i = .P ∨ Items.type s.items i = .R →
    (∃ x ∈ readStack s.tstack, Items.Below s.items x i) ∨ ∃ b ∈ blocks, InBlock g s.items b i

/-! ### The truncated tree and the relation at the boundaries -/

/-- A frame of the open path: vertex `v`, its finished out-edges `done`, the current tree edge
`o` (whose child is the next frame). -/
structure PathFrame where
  v : Nat
  done : List DfsOut
  o : DfsOut

/-- The DFS tree truncated along the open path, with the subtree `t` at its bottom. -/
def truncTree : List PathFrame → DfsTree → DfsTree
  | [], t => t
  | f :: fs, t => .node f.v (f.done ++ [.tree f.o.e f.o.cls (truncTree fs t)])

/-- The relation at the end of `walkTree t d` (`d = fs.length`) under the path `fs`, with `prev`
the finished trees of the forest and `base` the tstack when the walk of `t` started. -/
structure StSim (g : Graph) (prev : List DfsTree) (fs : List PathFrame) (t : DfsTree)
    (base : List TEntry) (s : WalkState) : Prop where
  read : ∃ new, s.tstack = new ++ base ∧ StRead s.items new (refTree g t fs.length (DirsOf s fs.length)).1
  items : StItems g s (refBlocks g (prev ++ [truncTree fs t]))

/-- The relation at the start of an out-edge of `v` (`d = fs.length`): the out-edges `done` are
finished, `hasVert` is the walk's flag. -/
structure StSimOuts (g : Graph) (prev : List DfsTree) (fs : List PathFrame) (v : Nat)
    (done : List DfsOut) (hasVert : Bool) (base : List TEntry) (s : WalkState) : Prop where
  read : ∃ new, s.tstack = new ++ base ∧
    StRead s.items new (refOuts g v fs.length (DirsOf s fs.length) done false).1
  hasVert : (refOuts g v fs.length (DirsOf s fs.length) done false).2.2 = hasVert
  items : StItems g s (refBlocks g (prev ++ [truncTree fs (.node v done)]))

end Spqr
