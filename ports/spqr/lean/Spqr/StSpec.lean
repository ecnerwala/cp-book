import Mathlib.Data.List.Nodup
import Spqr.Correctness

/-!
# The st-ordering layer of the output

On top of the decomposition, the node-vertex order of every S/P/R node is an st-numbering of its
skeleton (`StNumbered`), its node-edges are listed in a dominance order (`EdgeDominance`), and the
adjacency rows are the "center is longest" bracket order (`AdjBracket`). `spqrTree_st` states this
for the output; it is split into the walk-level `walk_st` (the items' children lists are
st-orderings, `Items.StNumbered`) and the relabel-level `relabel_st`, following the route of
`Spqr.Correctness`. The walk invariant behind `walk_st` is in `Spqr.StWalk`; see `PROOF.md` §7.
-/

namespace Spqr

namespace SpqrTree

variable (t : SpqrTree)

/-- Row `r` of the adjacency CSR (`2 nv` and `2 nv + 1` are the two rows of node-vertex `nv`). -/
def adjRow (r : Nat) : List NodeAdj :=
  (List.range (t.adjBounds[r + 1]! - t.adjBounds[r]!)).map fun k => t.adjDat[t.adjBounds[r]! + k]!

/-- The node-edges of `i` other than its cap. -/
def ownEdges (i : Nat) : List NodeEdge := (t.nodeEdgesOf i).drop (if t.hasCap i then 1 else 0)

/-- The node-vertex order of node `i` is an st-numbering of its skeleton: every edge goes from a
lower to a higher node-vertex, and every interior node-vertex has an edge to a lower and an edge
to a higher node-vertex. (The first and last node-vertices are the cap endpoints.) -/
def StNumbered (i : Nat) : Prop :=
  (∀ p ∈ t.skeleton i, p.1 < p.2) ∧
  ∀ nv, (t.nvRange i).1 < nv → nv + 1 < (t.nvRange i).2 →
    (∃ p ∈ t.skeleton i, p.2 = nv) ∧ (∃ p ∈ t.skeleton i, p.1 = nv)

/-- Among the non-cap node-edges of `i`, a later edge is never dominated componentwise by an
earlier edge with different endpoints (the cap `(s, e-1)` dominates everything and comes first). -/
def EdgeDominance (i : Nat) : Prop :=
  ((t.ownEdges i).map (·.nvs)).Pairwise fun p q => p ≠ q → ¬ (q.1 ≤ p.1 ∧ q.2 ≤ p.2)

/-- Bracket order of the adjacency of node-vertex `nv`: row `2 nv` holds the incidences to lower
node-vertices, row `2 nv + 1` those to higher ones, each with non-increasing destination, so that
reading `row 2nv ++ row (2nv+1)` the center is the longest edge. -/
def AdjBracket (nv : Nat) : Prop :=
  (∀ a ∈ t.adjRow (2 * nv), a.destNv < nv) ∧ (∀ a ∈ t.adjRow (2 * nv + 1), nv < a.destNv) ∧
  ((t.adjRow (2 * nv)).map (·.destNv)).Pairwise (· ≥ ·) ∧
  ((t.adjRow (2 * nv + 1)).map (·.destNv)).Pairwise (· ≥ ·)

/-- The st-ordering guarantees of the output. -/
structure StOrder : Prop where
  st : ∀ i, i < t.size → t.type i = .S ∨ t.type i = .P ∨ t.type i = .R → t.StNumbered i
  dom : ∀ i, i < t.size → t.EdgeDominance i
  /-- Nodes with a single node-vertex (loops) are excepted. -/
  adj : ∀ nv i, nv < t.nodeVerts.size → t.nodeOfNv nv = some i → 1 < t.nVerts i → t.AdjBracket nv

/-- The V children of a node occur in `children i` in increasing node-vertex order: the `k`-th V
child is node-vertex `nvSt + |a| + k` of `nv_layout`. -/
theorem vchildren_nv_increasing (h : t.WF) (i : Nat) (hi : i < t.size) :
    (((t.children i).filter fun c => t.type c = .V).map fun c => (t.nodeVertsOf i).idxOf ⟨i, c⟩).Pairwise (· < ·) := by
  obtain ⟨a, b, hl, -, -, -⟩ := h.own.nv_layout i hi
  have hnd : (t.nodeVertsOf i).Nodup := List.Nodup.of_map _ (h.own.nv_distinct i hi)
  generalize hm0 : ((t.children i).filter fun c => t.type c = .V) = m0 at hl ⊢
  have hidx : ∀ j (hj : j < m0.length), (t.nodeVertsOf i).idxOf ⟨i, m0[j]⟩ = a.length + j := by
    intro j hj
    have hlt : a.length + j < (t.nodeVertsOf i).length := by
      rw [hl]; simp only [List.length_append, List.length_map]; omega
    have hget : (t.nodeVertsOf i)[a.length + j] = ⟨i, m0[j]⟩ := by
      simp only [hl]
      rw [List.getElem_append_left (by simp only [List.length_append, List.length_map]; omega),
        List.getElem_append_right (Nat.le_add_right _ _), List.getElem_map]
      simp
    rw [← hget]; exact hnd.idxOf_getElem _ _
  rw [List.pairwise_map, List.pairwise_iff_getElem]
  intro j k hj hk hjk
  rw [hidx j hj, hidx k hk]; omega

end SpqrTree

namespace Items

variable (items : Items)

/-- The vertex list of an S/P/R item: first endpoint, the V children in order, second endpoint. -/
def vertList (i : ItemId) : List Nat :=
  (items.vs i).1.toList ++ (((items.ch i).filter fun c => items.type c = .V).map fun c => c - 1) ++
    (items.vs i).2.toList

/-- `xs` is an st-numbering of the edge list `es`: `xs` has no repeats, the edges join vertices of
`xs`, and every vertex of `xs` other than the first and the last has a neighbour earlier and a
neighbour later in `xs`. -/
def StList (xs : List Nat) (es : List (Nat × Nat)) : Prop :=
  xs.Nodup ∧ (∀ p ∈ es, p.1 ∈ xs ∧ p.2 ∈ xs ∧ p.1 ≠ p.2) ∧
  ∀ x ∈ xs, xs.head? ≠ some x → xs.getLast? ≠ some x →
    (∃ p ∈ es, (p.1 = x ∧ xs.idxOf p.2 < xs.idxOf x) ∨ (p.2 = x ∧ xs.idxOf p.1 < xs.idxOf x)) ∧
    (∃ p ∈ es, (p.1 = x ∧ xs.idxOf x < xs.idxOf p.2) ∨ (p.2 = x ∧ xs.idxOf x < xs.idxOf p.1))

/-- The children list of item `i` is in s-t order: with endpoints `(s, t)`, `vertList i` is an
st-numbering of the skeleton (the children's `vs` plus `(s, t)`), and every child's `vs` is
oriented along it. -/
def StItem (i : ItemId) : Prop :=
  ∃ s t, items.vs i = (some s, some t) ∧
    StList (items.vertList i) ((s, t) :: items.virtualEdges i) ∧
    ∀ p ∈ items.virtualEdges i, (items.vertList i).idxOf p.1 < (items.vertList i).idxOf p.2

/-- Item-level st-ordering: every S/P/R item is in s-t order. -/
def StNumbered : Prop :=
  ∀ i, i < items.size → items.type i = .S ∨ items.type i = .P ∨ items.type i = .R → items.StItem i

end Items

/-! ### Relabel-level facts -/

/-- A list sorted by endpoint-position sum is in dominance order. -/
theorem pairwise_dominance_of_sorted_sum (l : List (Nat × Nat))
    (h : l.Pairwise fun p q => p.1 + p.2 ≤ q.1 + q.2) :
    l.Pairwise fun p q => p ≠ q → ¬ (q.1 ≤ p.1 ∧ q.2 ≤ p.2) := by
  refine h.imp ?_
  intro p q hle hne ⟨h1, h2⟩
  apply hne
  ext <;> omega

/-- `orderedChildren` returns a permutation of the children sorted (stably) by endpoint-position
sum; for non-R items it is the children list itself. -/
theorem orderedChildren_sorted (it : Item) (nvSt : Nat) (s : RelabelState) (hR : it.type = .R) :
    let loc (c : ItemId) : Nat :=
      if c < 1 + s.g.nv then 2 * (s.vertPos[c - 1]! - nvSt)
      else (s.vertPos[(s.items[c]!.vs).1.getD 0]! - nvSt) + (s.vertPos[(s.items[c]!.vs).2.getD 0]! - nvSt)
    let out := ((RelabelM.orderedChildren it nvSt).run s).1
    out.Perm it.ch ∧ out.Pairwise fun a b => loc a ≤ loc b := by
  intro loc out
  have hout : out = it.ch.mergeSort fun a b => decide (loc a ≤ loc b) := by
    simp only [out, RelabelM.orderedChildren, hR]; rfl
  rw [hout]
  refine ⟨List.mergeSort_perm _ _, ?_⟩
  have := List.pairwise_mergeSort (le := fun a b => decide (loc a ≤ loc b))
    (fun a b c => by simp only [decide_eq_true_eq]; omega)
    (fun a b => by simp only [Bool.or_eq_true, decide_eq_true_eq]; omega) it.ch
  exact this.imp fun h => by simpa using h

theorem orderedChildren_eq_of_ne_R (it : Item) (nvSt : Nat) (s : RelabelState) (hR : it.type ≠ .R) :
    ((RelabelM.orderedChildren it nvSt).run s).1 = it.ch := by
  have : (it.type != NodeType.R) = true := by simpa using hR
  simp only [RelabelM.orderedChildren, this]; rfl

/-- Edge children of an R node, stably sorted by endpoint-position sum and oriented along the
node-vertex order, are in dominance order. -/
theorem edgeChildren_dominance (l : List (Nat × Nat))
    (h : l.Pairwise fun p q => p.1 + p.2 ≤ q.1 + q.2) :
    l.Pairwise fun p q => p ≠ q → ¬ (q.1 ≤ p.1 ∧ q.2 ≤ p.2) :=
  pairwise_dominance_of_sorted_sum l h

/-- The reverse fill of `layoutNode` for an R node: with `edgeChildren` in dominance order and each
`(a, c)` oriented `a < c`, every row `2 nv` of the layout holds the lower incidences of `nv`, row
`2 nv + 1` the higher ones, with non-increasing destination. -/
theorem layoutNode_r_bracket (node nvSt nvEn neSt neEn : Nat) (edgeChildren : List (Nat × Nat))
    (hor : ∀ p ∈ edgeChildren, nvSt ≤ p.1 ∧ p.1 < p.2 ∧ p.2 < nvEn)
    (hdom : edgeChildren.Pairwise fun p q => p ≠ q → ¬ (q.1 ≤ p.1 ∧ q.2 ≤ p.2))
    (hnd : edgeChildren.Nodup) (hne : neEn = neSt + edgeChildren.length + 1) :
    let l := layoutNode .R node nvSt nvEn neSt neEn edgeChildren
    ∀ nv, nvSt ≤ nv → nv < nvEn →
      let row (r : Nat) : List NodeAdj :=
        (List.range (l.adjBounds[r + 1 - 2 * nvSt]! - l.adjBounds[r - 2 * nvSt]!)).map fun k =>
          l.adjDat[l.adjBounds[r - 2 * nvSt]! + k - 2 * neSt]!
      (∀ a ∈ row (2 * nv), a.destNv < nv) ∧ (∀ a ∈ row (2 * nv + 1), nv < a.destNv) ∧
      ((row (2 * nv)).map (·.destNv)).Pairwise (· ≥ ·) ∧
      ((row (2 * nv + 1)).map (·.destNv)).Pairwise (· ≥ ·) := by
  sorry

/-! ### The theorems -/

/-- Phase 2: the walk's children lists are in s-t order (`PROOF.md` §7, via `WalkState.StInv`). -/
theorem walk_st (g : Graph) (tern : Bool) (vo eo : List Nat) :
    Items.StNumbered (g.walk tern (g.dfsForest vo eo)).items := by
  sorry

/-- Phase 3: relabelling well-formed items in s-t order gives st-ordered output. -/
theorem relabel_st (g : Graph) (items : Items) (hst : items.StNumbered) (h : items.WF g) :
    (relabelTree g items).StOrder := by
  sorry

theorem spqrTree_st (g : Graph) (tern : Bool) (vo eo : List Nat) : (g.spqrTree tern vo eo).StOrder := by
  rw [spqrTree_eq]; exact relabel_st g _ (walk_st g tern vo eo) (walk_items_wf g tern vo eo)

end Spqr
