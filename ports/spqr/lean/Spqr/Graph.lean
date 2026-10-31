namespace Spqr

/-- An undirected multigraph on vertices `0, …, nv-1`; `edges[e] = (u, v)`. Self-loops and
parallel edges are allowed. -/
structure Graph where
  nv : Nat
  edges : Array (Nat × Nat)
deriving Repr

namespace Graph

def ne (g : Graph) : Nat := g.edges.size

end Graph

/-- The indices listed in `order` (in that order), followed by the remaining indices in `[0, n)`
in increasing order. Compiled code uses `inOrderFast` (see `inOrder_eq_inOrderFast`). -/
def inOrder (n : Nat) (order : List Nat) : List Nat :=
  order ++ (List.range n).filter (fun i => !order.contains i)

/-- `marks n order` marks the indices in `order` that are `< n`. -/
def marks (n : Nat) (order : List Nat) : Array Bool :=
  order.foldl (init := Array.replicate n false) fun m i =>
    if h : i < m.size then m.set i true h else m

def inOrderFast (n : Nat) (order : List Nat) : List Nat :=
  let m := marks n order
  order ++ (List.range n).filter (fun i => !m.getD i false)

theorem size_marks_foldl (init : Array Bool) (order : List Nat) :
    (order.foldl (fun m i => if h : i < m.size then m.set i true h else m) init).size = init.size := by
  induction order generalizing init with
  | nil => rfl
  | cons j order ih =>
    simp only [List.foldl]
    split <;> simp [ih]

theorem getElem_marks_foldl (init : Array Bool) (order : List Nat) (i : Nat) (hi : i < init.size) :
    (order.foldl (fun m i => if h : i < m.size then m.set i true h else m) init)[i]'(by
        rw [size_marks_foldl]; exact hi) = (init[i] || order.contains i) := by
  induction order generalizing init with
  | nil => simp
  | cons j order ih =>
    simp only [List.foldl]
    split
    · rw [ih _ (by simpa using hi)]
      by_cases hij : i = j
      · subst hij; simp
      · simp [Array.getElem_set, Ne.symm hij, hij]
    · rw [ih _ hi]
      have hne : i ≠ j := by omega
      simp [hne]

theorem size_marks (n : Nat) (order : List Nat) : (marks n order).size = n := by
  simp [marks, size_marks_foldl]

theorem marks_getD (n : Nat) (order : List Nat) (i : Nat) (hi : i < n) :
    (marks n order).getD i false = order.contains i := by
  rw [Array.getD_eq_getD_getElem?, Array.getElem?_eq_getElem (by rw [size_marks]; exact hi)]
  simp only [Option.getD_some, marks]
  rw [getElem_marks_foldl _ _ _ (by simpa using hi)]
  simp

@[csimp] theorem inOrder_eq_inOrderFast : @inOrder = @inOrderFast := by
  funext n order
  simp only [inOrder, inOrderFast]
  congr 1
  apply List.filter_congr
  intro i hi
  rw [marks_getD n order i (List.mem_range.mp hi)]

/-- Adjacency lists: `adj[v]` lists `(dest, e)` for every edge `e` incident to `v`, in `edgeOrder`
order. A self-loop appears once in its vertex's list. -/
def Graph.adjacency (g : Graph) (edgeOrder : List Nat) : Array (List (Nat × Nat)) :=
  let adj := (inOrder g.ne edgeOrder).foldl (init := Array.replicate g.nv [])
    fun adj e =>
      let (u, v) := g.edges[e]!
      let adj := adj.modify u ((v, e) :: ·)
      if u != v then adj.modify v ((u, e) :: ·) else adj
  adj.map List.reverse

end Spqr
