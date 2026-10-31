/-!
# Join lists

Lists with `O(1)` append, used for the ear spans and item children of the walk. Each one is
flattened by `toList` at most once, when it becomes a node's child list.
-/

namespace Spqr

inductive CatList (α : Type u) where
  | nil
  | single (a : α)
  | append (l r : CatList α)
deriving Repr, Inhabited

namespace CatList

variable {α : Type u}

/-- The reference semantics. Compiled code uses `toListTR` (see `toList_eq_toListTR`). -/
def toList : CatList α → List α
  | nil => []
  | single a => [a]
  | append l r => l.toList ++ r.toList

def toListTR.go : CatList α → List α → List α
  | nil, acc => acc
  | single a, acc => a :: acc
  | append l r, acc => go l (go r acc)

def toListTR (c : CatList α) : List α := toListTR.go c []

theorem toListTR.go_eq (c : CatList α) (acc : List α) : go c acc = c.toList ++ acc := by
  induction c generalizing acc with
  | nil => rfl
  | single a => rfl
  | append l r ihl ihr => simp [go, toList, ihl, ihr]

@[csimp] theorem toList_eq_toListTR : @toList = @toListTR := by
  funext α c
  simp [toListTR, toListTR.go_eq]

/-- Append, skipping empty operands so that every `append` node has a non-empty list on each
side. -/
def cat : CatList α → CatList α → CatList α
  | nil, r => r
  | l, nil => l
  | l, r => append l r

instance : Append (CatList α) := ⟨cat⟩

@[simp] theorem toList_nil : (nil : CatList α).toList = [] := rfl
@[simp] theorem toList_single (a : α) : (single a).toList = [a] := rfl
@[simp] theorem toList_append' (l r : CatList α) : (append l r).toList = l.toList ++ r.toList := rfl

@[simp] theorem toList_cat (l r : CatList α) : (l ++ r).toList = l.toList ++ r.toList := by
  show (cat l r).toList = _
  cases l <;> cases r <;> simp [cat, toList]

def head? : CatList α → Option α
  | nil => none
  | single a => some a
  | append l r => l.head? <|> r.head?

@[simp] theorem head?_eq_toList_head? (c : CatList α) : c.head? = c.toList.head? := by
  induction c with
  | nil => rfl
  | single a => rfl
  | append l r ihl ihr =>
    simp only [head?, toList, ihl, ihr]
    cases l.toList <;> simp

def isEmpty (c : CatList α) : Bool := c.head?.isNone

theorem isEmpty_iff (c : CatList α) : c.isEmpty = true ↔ c.toList = [] := by
  simp [isEmpty, Option.isNone_iff_eq_none, List.head?_eq_none_iff]

end CatList

/-- A join list with its first element cached, so that `head!` is `O(1)`. -/
structure Span (α : Type u) where
  cat : CatList α
  head? : Option α
  valid : head? = cat.toList.head?

namespace Span

variable {α : Type u}

def nil : Span α := ⟨.nil, none, rfl⟩
def single (a : α) : Span α := ⟨.single a, some a, rfl⟩
def append (l r : Span α) : Span α :=
  ⟨l.cat ++ r.cat, l.head?.or r.head?, by simp [l.valid, r.valid]⟩
def toList (s : Span α) : List α := s.cat.toList
def head! [Inhabited α] (s : Span α) : α := s.head?.getD default

instance : Append (Span α) := ⟨append⟩
instance : Inhabited (Span α) := ⟨nil⟩
instance [Repr α] : Repr (Span α) := ⟨fun s n => reprPrec s.toList n⟩

@[simp] theorem toList_nil : (nil : Span α).toList = [] := rfl
@[simp] theorem toList_single (a : α) : (single a).toList = [a] := rfl
@[simp] theorem toList_append (l r : Span α) : (l ++ r).toList = l.toList ++ r.toList := by
  show (append l r).toList = _
  simp [append, toList]

theorem head!_eq [Inhabited α] (s : Span α) : s.head! = s.toList.head! := by
  show s.head?.getD default = s.cat.toList.head!
  rw [s.valid]
  cases s.cat.toList <;> rfl

theorem head!_toList [Inhabited α] (s : Span α) : s.toList.head! = s.head! := (head!_eq s).symm

end Span

end Spqr
