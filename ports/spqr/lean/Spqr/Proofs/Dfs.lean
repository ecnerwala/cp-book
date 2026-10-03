import Mathlib.Data.List.Nodup
import Mathlib.Data.List.Perm.Subperm
import Spqr.Dfs

/-!
# Phase 1 proofs: the lowpoint DFS

`dfsVisit_spec` is the per-visit invariant; `dfsForest_spanning'` is Lemma 1.1 of `PROOF.md`,
the lowvals / `classify` lemmas at the end are Lemma 1.2.
-/

namespace Spqr

/-! ### Two smallest distinct values -/

/-- `min` over `d :: xs`. -/
def lmin (d : Nat) (xs : List Nat) : Nat := xs.foldr min d

/-- The two smallest distinct values of `d :: xs`; the second is `d` when there is only one.
This is what `mergeLowvals` computes. -/
def low2 (d : Nat) (xs : List Nat) : Lowvals :=
  (lmin d xs, lmin d (xs.filter (· ≠ lmin d xs)))

@[simp] theorem lmin_nil (d : Nat) : lmin d [] = d := rfl
@[simp] theorem lmin_cons (d x : Nat) (xs : List Nat) : lmin d (x :: xs) = min x (lmin d xs) := rfl

theorem lmin_le (d : Nat) (xs : List Nat) : lmin d xs ≤ d := by
  induction xs with
  | nil => simp
  | cons x xs ih => simp; omega

theorem lmin_le_of_mem {d x : Nat} {xs : List Nat} (h : x ∈ xs) : lmin d xs ≤ x := by
  induction xs with
  | nil => simp at h
  | cons y xs ih =>
    simp at h ⊢
    rcases h with rfl | h
    · omega
    · have := ih h; omega

theorem lmin_eq_or_mem (d : Nat) (xs : List Nat) : lmin d xs = d ∨ lmin d xs ∈ xs := by
  induction xs with
  | nil => simp
  | cons x xs ih =>
    simp only [lmin_cons, List.mem_cons]
    rcases Nat.le_total x (lmin d xs) with hx | hx
    · rw [Nat.min_eq_left hx]; simp
    · rw [Nat.min_eq_right hx]
      rcases ih with h | h
      · exact .inl h
      · exact .inr (.inr h)

theorem le_lmin {d m : Nat} {xs : List Nat} (hd : m ≤ d) (h : ∀ x ∈ xs, m ≤ x) : m ≤ lmin d xs := by
  induction xs with
  | nil => simpa
  | cons x xs ih =>
    simp at h
    exact Nat.le_min.2 ⟨h.1, ih h.2⟩

theorem lmin_eq_iff {d m : Nat} {xs : List Nat} :
    lmin d xs = m ↔ m ≤ d ∧ (∀ x ∈ xs, m ≤ x) ∧ (m = d ∨ m ∈ xs) := by
  constructor
  · rintro rfl
    exact ⟨lmin_le d xs, fun x hx => lmin_le_of_mem hx, lmin_eq_or_mem d xs⟩
  · rintro ⟨h1, h2, h3⟩
    have := le_lmin h1 h2
    rcases h3 with rfl | h3
    · have := lmin_le m xs; omega
    · have := lmin_le_of_mem (d := d) h3; omega

theorem lmin_congr {d : Nat} {xs ys : List Nat} (h : ∀ x, x ∈ xs ↔ x ∈ ys) : lmin d xs = lmin d ys := by
  rw [lmin_eq_iff]
  refine ⟨lmin_le d ys, fun x hx => lmin_le_of_mem ((h x).1 hx), ?_⟩
  rcases lmin_eq_or_mem d ys with h' | h'
  · exact .inl h'
  · exact .inr ((h _).2 h')

theorem lmin_append (d : Nat) (xs ys : List Nat) : lmin d (xs ++ ys) = min (lmin d xs) (lmin d ys) := by
  induction xs with
  | nil => have := lmin_le d ys; simp; omega
  | cons x xs ih => simp only [List.cons_append, lmin_cons, ih]; omega

theorem min_lmin (a b : Nat) (xs : List Nat) : min a (lmin b xs) = lmin (min a b) xs := by
  induction xs with
  | nil => simp
  | cons x xs ih => simp [← ih]; omega

theorem low2_congr {d : Nat} {xs ys : List Nat} (h : ∀ x, x ∈ xs ↔ x ∈ ys) : low2 d xs = low2 d ys := by
  unfold low2
  rw [lmin_congr h]
  congr 1
  exact lmin_congr fun x => by simp only [List.mem_filter, h x]

theorem low2_perm {d : Nat} {xs ys : List Nat} (h : xs.Perm ys) : low2 d xs = low2 d ys :=
  low2_congr fun _ => h.mem_iff

theorem filter_ne_eq_self (xs : List Nat) (m : Nat) (h : ∀ x ∈ xs, m < x) :
    xs.filter (· ≠ m) = xs :=
  List.filter_eq_self.2 fun x hx => by have := h x hx; simp; omega

/-- `mergeLowvals` merges the lowvals of a subtree hanging at depth `s ≥ d` into those at depth `d`. -/
theorem mergeLowvals_low2 {d s : Nat} (hs : d ≤ s) (xs ys : List Nat) :
    mergeLowvals (low2 d xs) (low2 s ys) = low2 d (xs ++ ys) := by
  have hdy : ∀ zs, lmin d zs = min d (lmin s zs) := fun zs => by
    rw [min_lmin, Nat.min_eq_left hs]
  have ha := lmin_le d xs
  have hb := lmin_le d (xs.filter (· ≠ lmin d xs))
  unfold mergeLowvals low2
  simp only [lmin_append, List.filter_append]
  split
  · next hlt =>
    have hmin : min (lmin d xs) (lmin d ys) = lmin s ys := by rw [hdy ys]; omega
    rw [hmin, filter_ne_eq_self xs (lmin s ys) fun x hx => by
      have := lmin_le_of_mem (d := d) hx; omega, hdy (ys.filter _)]
    simp only [Prod.mk.injEq, true_and]
    omega
  · next hlt =>
    have hmin : min (lmin d xs) (lmin d ys) = lmin d xs := by rw [hdy ys]; omega
    rw [hmin, hdy (ys.filter _)]
    split
    · next heq =>
      rw [beq_iff_eq] at heq
      rw [heq]
      simp only [Prod.mk.injEq, true_and]
      omega
    · next hne =>
      rw [beq_iff_eq] at hne
      rw [filter_ne_eq_self ys (lmin d xs) fun x hx => by
        have := lmin_le_of_mem (d := s) hx; omega]
      simp only [Prod.mk.injEq, true_and]
      omega

theorem low2_singleton {d x : Nat} (h : x ≤ d) : low2 d [x] = (x, d) := by
  unfold low2
  simp [lmin, h]

/-! ### Trees -/

namespace DfsOut

def verts : DfsOut → List Nat
  | back .. => []
  | tree _ _ c => c.verts

def edges : DfsOut → List Nat
  | back e .. => [e]
  | tree e _ c => e :: c.edges

theorem vertsList_eq (l : List DfsOut) : vertsList l = l.flatMap verts := by
  induction l with
  | nil => rfl
  | cons o rest ih => cases o <;> simp [vertsList, verts, ih]

theorem edgesList_eq (l : List DfsOut) : edgesList l = l.flatMap edges := by
  induction l with
  | nil => rfl
  | cons o rest ih => cases o <;> simp [edgesList, edges, ih]

end DfsOut

mutual
/-- Depths returned to by the back edges in a tree whose root is at depth `d`. -/
def DfsTree.retDepths (d : Nat) : DfsTree → List Nat
  | .node _ outs => DfsOut.retDepthsList d outs
def DfsOut.retDepthsList (d : Nat) : List DfsOut → List Nat
  | [] => []
  | .back _ _ cls :: rest => cls.lowval d :: DfsOut.retDepthsList d rest
  | .tree _ _ child :: rest => child.retDepths (d + 1) ++ DfsOut.retDepthsList d rest
end

namespace DfsOut

def retDepths (d : Nat) : DfsOut → List Nat
  | back _ _ cls => [cls.lowval d]
  | tree _ _ c => c.retDepths (d + 1)

theorem retDepthsList_eq (d : Nat) (l : List DfsOut) : retDepthsList d l = l.flatMap (retDepths d) := by
  induction l with
  | nil => rfl
  | cons o rest ih => cases o <;> simp [retDepthsList, retDepths, ih]

end DfsOut

mutual
/-- Out-edge `o` of the vertex `v` whose proper ancestors are `anc` (`anc[i]` at depth `i`): a back
edge goes to an ancestor or to `v` and is classified by that depth; a tree edge is classified by
the child's lowvals, and the child is itself well formed. -/
def DfsOut.WF (anc : List Nat) (v : Nat) : DfsOut → Prop
  | .back _ dest cls =>
    ∃ i, (anc ++ [v])[i]? = some dest ∧ cls = classify anc.length false (i, anc.length)
  | .tree _ cls child =>
    child.WF (anc ++ [v]) ∧
      cls = classify anc.length true (low2 (anc.length + 1) (child.retDepths (anc.length + 1)))
/-- A DFS tree hanging below the ancestor path `anc`: every out-list is sorted by `OutClass.rank`
and every out-edge is `DfsOut.WF`. -/
def DfsTree.WF (anc : List Nat) : DfsTree → Prop
  | .node v outs => outs.Pairwise (fun a b => a.cls.rank ≤ b.cls.rank) ∧ ∀ o ∈ outs, o.WF anc v
end

/-! ### `dfsVisit` unfolded -/

/-- The body of the fold in `dfsVisit`. -/
def dfsStep (adj : Array (List (Nat × Nat))) (fuel d : Nat) (prvE : Option Nat) :
    List DfsOut × Lowvals × Array (Option Nat) → Nat × Nat →
      List DfsOut × Lowvals × Array (Option Nat)
  | (outs, lv, depth), (nxt, e) =>
    if some e == prvE || (depth[nxt]!).any (· > d) then (outs, lv, depth)
    else match depth[nxt]! with
      | none =>
        let (child, n, depth) := dfsVisit adj fuel nxt (d + 1) (some e) depth
        (.tree e (classify d true n) child :: outs, mergeLowvals lv n, depth)
      | some nd =>
        (.back e nxt (classify d false (nd, d)) :: outs, mergeLowvals lv (nd, d), depth)

theorem dfsVisit_succ (adj : Array (List (Nat × Nat))) (fuel v d : Nat) (prvE : Option Nat)
    (depth : Array (Option Nat)) :
    dfsVisit adj (fuel + 1) v d prvE depth =
      let r := adj[v]!.foldl (dfsStep adj fuel d prvE) ([], (d, d), depth.set! v (some d))
      (.node v (r.1.reverse.mergeSort fun a b => a.cls.rank ≤ b.cls.rank), r.2.1, r.2.2) :=
  rfl

/-! ### Array helpers -/

theorem getElem!_modify_self {α : Type} [Inhabited α] (A : Array α) (i : Nat) (f : α → α)
    (hi : i < A.size) : (A.modify i f)[i]! = f A[i]! := by
  rw [getElem!_pos A i hi, getElem!_pos _ i (by simpa using hi), Array.getElem_modify_self]

theorem getElem!_modify_ne {α : Type} [Inhabited α] (A : Array α) (i j : Nat) (f : α → α)
    (h : i ≠ j) : (A.modify i f)[j]! = A[j]! := by
  by_cases hj : j < A.size
  · rw [getElem!_pos A j hj, getElem!_pos _ j (by simpa using hj),
      Array.getElem_modify_of_ne h]
  · rw [getElem!_neg A j hj, getElem!_neg _ j (by simpa using hj)]

theorem getElem!_map' {α β : Type} [Inhabited α] [Inhabited β] (f : α → β) (hf : f default = default)
    (A : Array α) (i : Nat) : (A.map f)[i]! = f A[i]! := by
  by_cases hi : i < A.size
  · rw [getElem!_pos A i hi, getElem!_pos _ i (by simpa using hi), Array.getElem_map]
  · rw [getElem!_neg A i hi, getElem!_neg _ i (by simpa using hi), hf]

theorem getElem!_replicate' {α : Type} [Inhabited α] (n : Nat) (x : α) (i : Nat) (hi : i < n) :
    (Array.replicate n x)[i]! = x := by
  rw [getElem!_pos _ i (by simpa using hi), Array.getElem_replicate]

theorem lt_size_of_getElem!_ne {α : Type} [Inhabited α] {A : Array α} {i : Nat}
    (h : A[i]! ≠ default) : i < A.size := by
  by_contra hi
  exact h (getElem!_neg A i hi)

/-! ### Adjacency lists -/

/-- All edge endpoints are vertices. -/
def Graph.WF (g : Graph) : Prop := ∀ p ∈ g.edges, p.1 < g.nv ∧ p.2 < g.nv

/-- A (possibly partial) order: distinct indices below `n`. -/
def OrderOK (n : Nat) (order : List Nat) : Prop := order.Nodup ∧ ∀ x ∈ order, x < n

theorem inOrder_perm {n : Nat} {order : List Nat} (h : OrderOK n order) :
    (inOrder n order).Perm (List.range n) := by
  have hnd : (inOrder n order).Nodup := by
    unfold inOrder
    refine List.Nodup.append h.1 (List.Nodup.filter _ List.nodup_range) ?_
    intro x hx hx'
    simp at hx'
    exact hx'.2 hx
  refine (List.perm_ext_iff_of_nodup hnd List.nodup_range).2 fun x => ?_
  simp only [inOrder, List.mem_append, List.mem_filter, List.mem_range, Bool.not_eq_eq_eq_not,
    Bool.not_true, List.contains_eq_mem, decide_eq_false_iff_not]
  constructor
  · rintro (hx | ⟨hx, _⟩)
    · exact h.2 x hx
    · exact hx
  · intro hx
    by_cases hm : x ∈ order
    · exact .inl hm
    · exact .inr ⟨hx, hm⟩

/-- What the DFS needs from the adjacency lists: entries are in range, each edge index occurs at
most once per list, in the lists of both of its endpoints and nowhere else. -/
structure AdjOK (g : Graph) (adj : Array (List (Nat × Nat))) : Prop where
  bounds : ∀ x y e : Nat, (y, e) ∈ adj[x]! → x < g.nv ∧ y < g.nv ∧ e < g.ne
  nodup : ∀ x : Nat, (List.map (Prod.snd : Nat × Nat → Nat) adj[x]!).Nodup
  symm : ∀ x y e : Nat, (y, e) ∈ adj[x]! → (x, e) ∈ adj[y]!
  ends : ∀ x y x' y' e : Nat, (y, e) ∈ adj[x]! → (y', e) ∈ adj[x']! → x' = x ∨ x' = y
  cover : ∀ e, e < g.ne → ∃ x y : Nat, (y, e) ∈ adj[x]!

theorem AdjOK.snd_inj {g : Graph} {adj : Array (List (Nat × Nat))} (h : AdjOK g adj) {x y y' e : Nat}
    (h1 : (y, e) ∈ adj[x]!) (h2 : (y', e) ∈ adj[x]!) : y = y' :=
  (Prod.mk.inj (List.inj_on_of_nodup_map (h.nodup x) h1 h2 rfl)).1

/-- The body of the fold in `Graph.adjacency`. -/
def adjStep (g : Graph) (adj : Array (List (Nat × Nat))) (e : Nat) : Array (List (Nat × Nat)) :=
  let (u, v) := g.edges[e]!
  let adj := adj.modify u ((v, e) :: ·)
  if u != v then adj.modify v ((u, e) :: ·) else adj

theorem adjacency_eq (g : Graph) (eo : List Nat) :
    g.adjacency eo = ((inOrder g.ne eo).foldl (adjStep g) (Array.replicate g.nv [])).map List.reverse :=
  rfl

/-- Invariant of the fold in `Graph.adjacency` after processing the edges `es`. -/
def AdjInv (g : Graph) (es : List Nat) (A : Array (List (Nat × Nat))) : Prop :=
  A.size = g.nv ∧
  (∀ x y e : Nat, (y, e) ∈ A[x]! ↔
    e ∈ es ∧ e < g.ne ∧ x < g.nv ∧ (g.edges[e]! = (x, y) ∨ g.edges[e]! = (y, x))) ∧
  ∀ x : Nat, (List.map (Prod.snd : Nat × Nat → Nat) A[x]!).Nodup

theorem adjInv_nil (g : Graph) : AdjInv g [] (Array.replicate g.nv []) := by
  have h : ∀ x : Nat, (Array.replicate g.nv ([] : List (Nat × Nat)))[x]! = [] := by
    intro x
    by_cases hx : x < g.nv
    · exact getElem!_replicate' _ _ _ hx
    · exact getElem!_neg _ _ (by simpa using hx)
  refine ⟨by simp, fun x y e => ?_, fun x => ?_⟩ <;> simp [h]

theorem mem_modify_cons (B : Array (List (Nat × Nat))) {z : Nat} (hz : z < B.size) (q : Nat × Nat)
    (x : Nat) (p : Nat × Nat) : p ∈ (B.modify z (q :: ·))[x]! ↔ p ∈ B[x]! ∨ (x = z ∧ p = q) := by
  by_cases hx : x = z
  · subst hx
    rw [getElem!_modify_self _ _ _ hz]
    simp only [List.mem_cons, true_and]
    exact or_comm
  · rw [getElem!_modify_ne _ _ _ _ (Ne.symm hx)]
    simp [hx]

theorem nodup_modify_cons (B : Array (List (Nat × Nat))) {z : Nat} (hz : z < B.size) (q : Nat × Nat)
    (hnd : ∀ x : Nat, (List.map (Prod.snd : Nat × Nat → Nat) B[x]!).Nodup)
    (hq : q.2 ∉ List.map (Prod.snd : Nat × Nat → Nat) B[z]!) (x : Nat) :
    (List.map (Prod.snd : Nat × Nat → Nat) (B.modify z (q :: ·))[x]!).Nodup := by
  by_cases hx : x = z
  · subst hx
    rw [getElem!_modify_self _ _ _ hz]
    simp only [List.map_cons, List.nodup_cons]
    exact ⟨hq, hnd x⟩
  · rw [getElem!_modify_ne _ _ _ _ (Ne.symm hx)]
    exact hnd x

theorem Graph.WF.getElem! {g : Graph} (hg : g.WF) {e : Nat} (he : e < g.ne) :
    (g.edges[e]!).1 < g.nv ∧ (g.edges[e]!).2 < g.nv := by
  apply hg
  rw [getElem!_pos g.edges e (show e < g.edges.size from he)]
  exact Array.getElem_mem he

theorem adjInv_step {g : Graph} (hg : g.WF) {es : List Nat} {A : Array (List (Nat × Nat))}
    (hA : AdjInv g es A) {e : Nat} (he : e < g.ne) (hes : e ∉ es) :
    AdjInv g (es ++ [e]) (adjStep g A e) := by
  obtain ⟨hsize, hmem, hnd⟩ := hA
  have huv := hg.getElem! he
  rcases hE : g.edges[e]! with ⟨u, v⟩
  rw [hE] at huv
  have hnew : ∀ x : Nat, e ∉ List.map (Prod.snd : Nat × Nat → Nat) A[x]! := by
    intro x hx
    simp only [List.mem_map] at hx
    obtain ⟨⟨y, e'⟩, hy, rfl⟩ := hx
    exact hes ((hmem x y _).1 hy).1
  have hu : u < A.size := hsize ▸ huv.1
  have hv : v < A.size := hsize ▸ huv.2
  have hv' : v < (A.modify u ((v, e) :: ·)).size := by simpa using hv
  simp only [adjStep, hE]
  split
  · next hne =>
    have hne : u ≠ v := by simpa using hne
    refine ⟨by simpa using hsize, fun x y e' => ?_, ?_⟩
    · rw [mem_modify_cons _ hv', mem_modify_cons _ hu, hmem]
      simp only [List.mem_append, List.mem_singleton, Prod.mk.injEq]
      constructor
      · rintro ((⟨h1, h2, h3, h4⟩ | ⟨rfl, rfl, rfl⟩) | ⟨rfl, rfl, rfl⟩)
        · exact ⟨.inl h1, h2, h3, h4⟩
        · exact ⟨.inr rfl, he, huv.1, .inl hE⟩
        · exact ⟨.inr rfl, he, huv.2, .inr hE⟩
      · rintro ⟨h1 | rfl, h2, h3, h4⟩
        · exact .inl (.inl ⟨h1, h2, h3, h4⟩)
        · rw [hE] at h4
          rcases h4 with h4 | h4 <;> simp only [Prod.mk.injEq] at h4
          · exact .inl (.inr ⟨h4.1.symm, h4.2.symm, rfl⟩)
          · exact .inr ⟨h4.2.symm, h4.1.symm, rfl⟩
    · refine nodup_modify_cons _ hv' _ (nodup_modify_cons _ hu _ hnd (hnew u)) ?_
      rw [getElem!_modify_ne _ _ _ _ hne]
      exact hnew v
  · next hne =>
    have hne : u = v := by simpa using hne
    subst hne
    refine ⟨by simpa using hsize, fun x y e' => ?_, nodup_modify_cons _ hu _ hnd (hnew u)⟩
    rw [mem_modify_cons _ hu, hmem]
    simp only [List.mem_append, List.mem_singleton, Prod.mk.injEq]
    constructor
    · rintro (⟨h1, h2, h3, h4⟩ | ⟨rfl, rfl, rfl⟩)
      · exact ⟨.inl h1, h2, h3, h4⟩
      · exact ⟨.inr rfl, he, huv.1, .inl hE⟩
    · rintro ⟨h1 | rfl, h2, h3, h4⟩
      · exact .inl ⟨h1, h2, h3, h4⟩
      · rw [hE] at h4
        rcases h4 with h4 | h4 <;> simp only [Prod.mk.injEq] at h4
        · exact .inr ⟨h4.1.symm, h4.2.symm, rfl⟩
        · exact .inr ⟨h4.2.symm, h4.1.symm, rfl⟩

theorem adjInv_foldl {g : Graph} (hg : g.WF) :
    ∀ (rest done : List Nat) (A : Array (List (Nat × Nat))), AdjInv g done A →
      (done ++ rest).Nodup → (∀ e ∈ rest, e < g.ne) →
      AdjInv g (done ++ rest) (rest.foldl (adjStep g) A)
  | [], done, A, hA, _, _ => by simpa using hA
  | e :: rest, done, A, hA, hnd, hlt => by
    have hes : e ∉ done := fun h =>
      (List.nodup_cons.1 (List.nodup_middle.1 hnd)).1 (List.mem_append_left _ h)
    have h := adjInv_foldl hg rest (done ++ [e]) (adjStep g A e)
      (adjInv_step hg hA (hlt e (by simp)) hes) (by simpa using hnd)
      (fun e' h' => hlt e' (by simp [h']))
    simpa using h

/-- The adjacency entry `(y, e)` of `x` is the edge `e = {x, y}` of `g`. -/
theorem adjacency_ends {g : Graph} (hg : g.WF) {eo : List Nat} (heo : OrderOK g.ne eo)
    (x y e : Nat) (h : (y, e) ∈ (g.adjacency eo)[x]!) :
    g.edges[e]! = (x, y) ∨ g.edges[e]! = (y, x) := by
  have hperm := inOrder_perm heo
  have hinv := adjInv_foldl hg (inOrder g.ne eo) [] _ (adjInv_nil g)
    (by simpa using hperm.nodup_iff.2 List.nodup_range) (fun e he => by simpa using hperm.subset he)
  simp only [List.nil_append] at hinv
  obtain ⟨_, hmem, -⟩ := hinv
  rw [adjacency_eq, getElem!_map' List.reverse rfl _ _, List.mem_reverse, hmem] at h
  exact h.2.2.2

theorem adjacency_ok {g : Graph} (hg : g.WF) {eo : List Nat} (heo : OrderOK g.ne eo) :
    AdjOK g (g.adjacency eo) := by
  have hperm := inOrder_perm heo
  have hinv := adjInv_foldl hg (inOrder g.ne eo) [] _ (adjInv_nil g)
    (by simpa using hperm.nodup_iff.2 List.nodup_range) (fun e he => by simpa using hperm.subset he)
  simp only [List.nil_append] at hinv
  obtain ⟨_, hmem, hnd⟩ := hinv
  have hrev : ∀ x : Nat, (g.adjacency eo)[x]! =
      ((inOrder g.ne eo).foldl (adjStep g) (Array.replicate g.nv []))[x]!.reverse := fun x => by
    rw [adjacency_eq]; exact getElem!_map' List.reverse rfl _ _
  have hmem' : ∀ x y e : Nat, (y, e) ∈ (g.adjacency eo)[x]! ↔
      e < g.ne ∧ x < g.nv ∧ (g.edges[e]! = (x, y) ∨ g.edges[e]! = (y, x)) := by
    intro x y e
    rw [hrev, List.mem_reverse, hmem]
    constructor
    · rintro ⟨_, h⟩; exact h
    · intro h; exact ⟨hperm.mem_iff.2 (List.mem_range.2 h.1), h⟩
  refine ⟨?_, ?_, ?_, ?_, ?_⟩
  · intro x y e h
    obtain ⟨he, hx, hxy⟩ := (hmem' x y e).1 h
    have := hg.getElem! he
    rcases hxy with hxy | hxy <;> rw [hxy] at this
    · exact ⟨hx, this.2, he⟩
    · exact ⟨hx, this.1, he⟩
  · intro x
    rw [hrev, List.map_reverse]
    exact List.nodup_reverse.2 (hnd x)
  · intro x y e h
    obtain ⟨he, _, hxy⟩ := (hmem' x y e).1 h
    have := hg.getElem! he
    refine (hmem' y x e).2 ⟨he, ?_, hxy.symm⟩
    rcases hxy with hxy | hxy <;> rw [hxy] at this
    · exact this.2
    · exact this.1
  · intro x y x' y' e h h'
    obtain ⟨_, _, hxy⟩ := (hmem' x y e).1 h
    obtain ⟨_, _, hxy'⟩ := (hmem' x' y' e).1 h'
    rcases hxy with hxy | hxy <;> rcases hxy' with hxy' | hxy' <;> rw [hxy] at hxy' <;>
      simp only [Prod.mk.injEq] at hxy' <;> omega
  · intro e he
    have := hg.getElem! he
    exact ⟨_, _, (hmem' _ _ e).2 ⟨he, this.1, .inl rfl⟩⟩

/-! ### Tree lemmas -/

namespace DfsTree

theorem v_mem_verts (t : DfsTree) : t.v ∈ t.verts := by
  cases t; simp [verts, v]

theorem verts_perm_of_outs_perm {v : Nat} {outs outs' : List DfsOut} (h : outs.Perm outs') :
    (node v outs).verts.Perm (node v outs').verts := by
  simp only [verts, DfsOut.vertsList_eq]
  exact (h.flatMap_right _).cons v

theorem edges_perm_of_outs_perm {v : Nat} {outs outs' : List DfsOut} (h : outs.Perm outs') :
    (node v outs).edges.Perm (node v outs').edges := by
  simp only [edges, DfsOut.edgesList_eq]
  exact h.flatMap_right _

theorem retDepths_perm_of_outs_perm {v d : Nat} {outs outs' : List DfsOut} (h : outs.Perm outs') :
    ((node v outs).retDepths d).Perm ((node v outs').retDepths d) := by
  simp only [retDepths, DfsOut.retDepthsList_eq]
  exact h.flatMap_right _

end DfsTree

theorem lowval_classify_back {d nd : Nat} (h : nd ≤ d) : (classify d false (nd, d)).lowval d = nd := by
  unfold classify
  split
  · next h1 =>
    have : nd = d := by simp at h1; omega
    simp [this, OutClass.lowval]
  · rfl

/-! ### The `dfsVisit` invariant -/

/-- Precondition of `dfsVisit adj _ v d prvE depth`: the vertices on the active path `anc`
(`anc[i]` at depth `i`, `d = anc.length`) are the only visited vertices adjacent to an unvisited
one; `v` is unvisited; the parent edge `prvE` leads from `v` to a visited vertex. -/
structure Pre (g : Graph) (adj : Array (List (Nat × Nat))) (anc : List Nat) (v d : Nat)
    (prvE : Option Nat) (depth : Array (Option Nat)) : Prop where
  size : depth.size = g.nv
  hd : anc.length = d
  path : ∀ i w : Nat, anc[i]? = some w → depth[w]! = some i
  hv : v < g.nv
  unvisited : depth[v]! = none
  closed : ∀ w : Nat, w ∉ anc → depth[w]! ≠ none → ∀ p ∈ adj[w]!, depth[p.1]! ≠ none
  parent : ∀ pe : Nat, prvE = some pe → ∃ y : Nat, (y, pe) ∈ adj[v]! ∧ depth[y]! ≠ none

/-- Postcondition of `dfsVisit adj _ v d prvE depth = (t, lv, depth')`: `t` is a well-formed tree
rooted at `v` whose vertices are exactly the newly visited ones, whose edges are exactly the edges
at those vertices other than `prvE`, and `lv` are its lowvals. -/
structure Post (g : Graph) (adj : Array (List (Nat × Nat))) (anc : List Nat) (v d : Nat)
    (prvE : Option Nat) (depth : Array (Option Nat)) (t : DfsTree) (lv : Lowvals)
    (depth' : Array (Option Nat)) : Prop where
  size : depth'.size = g.nv
  mono : ∀ w k : Nat, depth[w]! = some k → depth'[w]! = some k
  root : t.v = v
  verts_iff : ∀ w : Nat, w ∈ t.verts ↔ depth[w]! = none ∧ depth'[w]! ≠ none
  verts_depth : ∀ w ∈ t.verts, ∃ k, d ≤ k ∧ depth'[w]! = some k
  verts_nodup : t.verts.Nodup
  closed : ∀ w : Nat, w ∉ anc → depth'[w]! ≠ none → ∀ p ∈ adj[w]!, depth'[p.1]! ≠ none
  edges_iff : ∀ e : Nat, e ∈ t.edges ↔ some e ≠ prvE ∧ ∃ x ∈ t.verts, ∃ y : Nat, (y, e) ∈ adj[x]!
  edges_nodup : t.edges.Nodup
  wf : t.WF anc
  lv : lv = low2 d (t.retDepths d)

/-- Invariant of the fold in `dfsVisit` after processing the adjacency entries `done` of `v`. -/
structure FoldInv (g : Graph) (adj : Array (List (Nat × Nat))) (anc : List Nat) (v d : Nat)
    (prvE : Option Nat) (depth : Array (Option Nat)) (done : List (Nat × Nat)) (outs : List DfsOut)
    (lv : Lowvals) (cur : Array (Option Nat)) : Prop where
  size : cur.size = g.nv
  mono : ∀ w k : Nat, depth[w]! = some k → cur[w]! = some k
  hv : cur[v]! = some d
  new : ∀ w : Nat, w ∈ DfsOut.vertsList outs ↔ w ≠ v ∧ depth[w]! = none ∧ cur[w]! ≠ none
  new_depth : ∀ w ∈ DfsOut.vertsList outs, ∃ k, d < k ∧ cur[w]! = some k
  nodup : (DfsOut.vertsList outs).Nodup
  closed : ∀ w : Nat, w ∉ anc → w ≠ v → cur[w]! ≠ none → ∀ p ∈ adj[w]!, cur[p.1]! ≠ none
  done_vis : ∀ p ∈ done, cur[p.1]! ≠ none
  edges_iff : ∀ e : Nat, e ∈ DfsOut.edgesList outs ↔ some e ≠ prvE ∧
    ((∃ y : Nat, (y, e) ∈ done) ∨ ∃ x ∈ DfsOut.vertsList outs, ∃ y : Nat, (y, e) ∈ adj[x]!)
  edges_nodup : (DfsOut.edgesList outs).Nodup
  wf : ∀ o ∈ outs, o.WF anc v
  lv : lv = low2 d (DfsOut.retDepthsList d outs)

variable {g : Graph} {adj : Array (List (Nat × Nat))} {anc : List Nat} {v d : Nat} {prvE : Option Nat}
  {depth : Array (Option Nat)}

namespace Pre

variable (h : Pre g adj anc v d prvE depth)
include h

theorem anc_depth {w : Nat} (hw : w ∈ anc) : ∃ i, i < d ∧ anc[i]? = some w ∧ depth[w]! = some i := by
  obtain ⟨i, hi⟩ := List.mem_iff_getElem?.1 hw
  exact ⟨i, h.hd ▸ (List.getElem?_eq_some_iff.1 hi).1, hi, h.path i w hi⟩

theorem not_mem_anc : v ∉ anc := fun hv => by
  obtain ⟨_, _, _, h'⟩ := h.anc_depth hv
  rw [h.unvisited] at h'
  exact nomatch h'

theorem anc_nodup : anc.Nodup := by
  rw [List.nodup_iff_injective_getElem]
  rintro ⟨i, hi⟩ ⟨j, hj⟩ hij
  simp only at hij
  have h1 := h.path i anc[i] (List.getElem?_eq_getElem hi)
  have h2 := h.path j anc[j] (List.getElem?_eq_getElem hj)
  rw [← hij, h1] at h2
  exact Fin.ext (Option.some.inj h2)

theorem anc_lt {w : Nat} (hw : w ∈ anc) : w < g.nv := by
  obtain ⟨_, _, _, h'⟩ := h.anc_depth hw
  exact h.size ▸ lt_size_of_getElem!_ne (by rw [h']; simp)

theorem succ_le : d + 1 ≤ g.nv := by
  have hnd : (anc ++ [v]).Nodup :=
    List.Nodup.append h.anc_nodup (List.nodup_singleton v) (fun w hw hw' => by
      rw [List.mem_singleton] at hw'
      exact h.not_mem_anc (hw' ▸ hw))
  have hsub : anc ++ [v] ⊆ List.range g.nv := fun w hw => by
    rw [List.mem_range]
    rcases List.mem_append.1 hw with hw | hw
    · exact h.anc_lt hw
    · exact (List.mem_singleton.1 hw) ▸ h.hv
  have := (List.Nodup.subperm hnd hsub).length_le
  simpa [h.hd] using this

theorem nbr_in_anc (hadj : AdjOK g adj) {y e : Nat} (hy : (y, e) ∈ adj[v]!)
    (hvis : depth[y]! ≠ none) : y ∈ anc := by
  by_contra hn
  exact h.closed y hn hvis (v, e) (hadj.symm _ _ _ hy) h.unvisited

end Pre

theorem foldInv_init (hpre : Pre g adj anc v d prvE depth) :
    FoldInv g adj anc v d prvE depth [] [] (d, d) (depth.set! v (some d)) := by
  have hsize : v < depth.size := hpre.size ▸ hpre.hv
  have hget : ∀ w, (depth.set! v (some d))[w]! = if w = v then some d else depth[w]! := by
    intro w
    by_cases hw : w = v
    · subst hw; rw [Array.getElem!_set!_self _ _ _ hsize]; simp
    · rw [Array.getElem!_set!_ne _ _ _ _ (Ne.symm hw)]; simp [hw]
  refine ⟨by simpa using hpre.size, ?_, by rw [hget]; simp, ?_, by simp [DfsOut.vertsList_eq],
    List.nodup_nil, ?_, by simp, by simp [DfsOut.edgesList_eq, DfsOut.vertsList_eq], List.nodup_nil,
    by simp, by simp [low2, lmin, DfsOut.retDepthsList_eq]⟩
  · intro w k hw
    rw [hget]
    split
    · next h => subst h; rw [hpre.unvisited] at hw; exact nomatch hw
    · exact hw
  · intro w
    simp only [DfsOut.vertsList_eq, List.flatMap_nil, List.not_mem_nil, false_iff, not_and, hget]
    intro hwv h1
    simp [hwv, h1]
  · intro w hw hwv hvis p hp
    rw [hget] at hvis ⊢
    simp only [hwv, ↓reduceIte] at hvis
    split
    · simp
    · exact hpre.closed w hw hvis p hp

theorem post_of_foldInv (hpre : Pre g adj anc v d prvE depth)
    {outs : List DfsOut} {lv : Lowvals} {cur : Array (Option Nat)}
    (hinv : FoldInv g adj anc v d prvE depth adj[v]! outs lv cur) :
    Post g adj anc v d prvE depth
      (.node v (outs.reverse.mergeSort fun a b => a.cls.rank ≤ b.cls.rank)) lv cur := by
  set sorted := outs.reverse.mergeSort fun a b => a.cls.rank ≤ b.cls.rank with hsorted
  have hperm : sorted.Perm outs := (List.mergeSort_perm _ _).trans (List.reverse_perm _)
  have hV := (DfsTree.verts_perm_of_outs_perm (v := v) hperm)
  have hE := (DfsTree.edges_perm_of_outs_perm (v := v) hperm)
  have hR := (DfsTree.retDepths_perm_of_outs_perm (v := v) (d := d) hperm)
  have hVmem : ∀ w, w ∈ (DfsTree.node v sorted).verts ↔ w = v ∨ w ∈ DfsOut.vertsList outs := by
    intro w; rw [hV.mem_iff]; simp [DfsTree.verts]
  have hEmem : ∀ e, e ∈ (DfsTree.node v sorted).edges ↔ e ∈ DfsOut.edgesList outs := by
    intro e; rw [hE.mem_iff]; simp [DfsTree.edges]
  have hvV : v ∉ DfsOut.vertsList outs := fun h => ((hinv.new v).1 h).1 rfl
  refine ⟨hinv.size, hinv.mono, rfl, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · intro w
    rw [hVmem, hinv.new]
    constructor
    · rintro (rfl | ⟨_, h1, h2⟩)
      · exact ⟨hpre.unvisited, by simp [hinv.hv]⟩
      · exact ⟨h1, h2⟩
    · rintro ⟨h1, h2⟩
      by_cases hw : w = v
      · exact .inl hw
      · exact .inr ⟨hw, h1, h2⟩
  · intro w hw
    rcases (hVmem w).1 hw with rfl | hw
    · exact ⟨d, Nat.le_refl d, hinv.hv⟩
    · obtain ⟨k, hk, hk'⟩ := hinv.new_depth w hw
      exact ⟨k, Nat.le_of_lt hk, hk'⟩
  · refine hV.symm.nodup_iff.1 ?_
    simp only [DfsTree.verts, List.nodup_cons]
    exact ⟨hvV, hinv.nodup⟩
  · intro w hw hvis p hp
    by_cases hwv : w = v
    · subst hwv; exact hinv.done_vis p hp
    · exact hinv.closed w hw hwv hvis p hp
  · intro e
    rw [hEmem, hinv.edges_iff]
    constructor
    · rintro ⟨h1, ⟨y, hy⟩ | ⟨x, hx, y, hy⟩⟩
      · exact ⟨h1, v, (hVmem v).2 (.inl rfl), y, hy⟩
      · exact ⟨h1, x, (hVmem x).2 (.inr hx), y, hy⟩
    · rintro ⟨h1, x, hx, y, hy⟩
      rcases (hVmem x).1 hx with rfl | hx
      · exact ⟨h1, .inl ⟨y, hy⟩⟩
      · exact ⟨h1, .inr ⟨x, hx, y, hy⟩⟩
  · exact hE.symm.nodup_iff.1 (by simpa [DfsTree.edges] using hinv.edges_nodup)
  · simp only [DfsTree.WF]
    refine ⟨?_, fun o ho => hinv.wf o (hperm.mem_iff.1 ho)⟩
    have := List.pairwise_mergeSort (le := fun a b : DfsOut => decide (a.cls.rank ≤ b.cls.rank))
      (fun a b c h1 h2 => by simp at h1 h2 ⊢; omega) (fun a b => by simp; omega) outs.reverse
    exact this.imp fun h => by simpa using h
  · rw [hinv.lv]
    exact low2_perm (by simpa [DfsTree.retDepths] using hR.symm)

namespace DfsOut

@[simp] theorem vertsList_back (e dest : Nat) (cls : OutClass) (rest : List DfsOut) :
    vertsList (.back e dest cls :: rest) = vertsList rest := rfl
@[simp] theorem vertsList_tree (e : Nat) (cls : OutClass) (child : DfsTree) (rest : List DfsOut) :
    vertsList (.tree e cls child :: rest) = child.verts ++ vertsList rest := rfl
@[simp] theorem edgesList_back (e dest : Nat) (cls : OutClass) (rest : List DfsOut) :
    edgesList (.back e dest cls :: rest) = e :: edgesList rest := rfl
@[simp] theorem edgesList_tree (e : Nat) (cls : OutClass) (child : DfsTree) (rest : List DfsOut) :
    edgesList (.tree e cls child :: rest) = e :: child.edges ++ edgesList rest := rfl
@[simp] theorem retDepthsList_back (d e dest : Nat) (cls : OutClass) (rest : List DfsOut) :
    retDepthsList d (.back e dest cls :: rest) = cls.lowval d :: retDepthsList d rest := rfl
@[simp] theorem retDepthsList_tree (d e : Nat) (cls : OutClass) (child : DfsTree) (rest : List DfsOut) :
    retDepthsList d (.tree e cls child :: rest) = child.retDepths (d + 1) ++ retDepthsList d rest := rfl

end DfsOut

/-- `dfsVisit` with `fuel` satisfies its spec whenever `fuel + d ≥ nv`. -/
def VisitSpec (g : Graph) (adj : Array (List (Nat × Nat))) (fuel : Nat) : Prop :=
  ∀ (anc : List Nat) (v d : Nat) (prvE : Option Nat) (depth : Array (Option Nat)),
    g.nv ≤ fuel + d → Pre g adj anc v d prvE depth →
    let r := dfsVisit adj fuel v d prvE depth
    Post g adj anc v d prvE depth r.1 r.2.1 r.2.2

theorem foldInv_step (hadj : AdjOK g adj) (hpre : Pre g adj anc v d prvE depth) {fuel : Nat}
    (hfuel : g.nv ≤ fuel + 1 + d) (ih : VisitSpec g adj fuel)
    {done : List (Nat × Nat)} {outs : List DfsOut} {lv : Lowvals} {cur : Array (Option Nat)}
    (hinv : FoldInv g adj anc v d prvE depth done outs lv cur) {nxt e : Nat}
    (hmem : (nxt, e) ∈ adj[v]!) (hdone : ∀ y, (y, e) ∉ done) (hsub : ∀ p ∈ done, p ∈ adj[v]!) :
    let r := dfsStep adj fuel d prvE (outs, lv, cur) (nxt, e)
    FoldInv g adj anc v d prvE depth (done ++ [(nxt, e)]) r.1 r.2.1 r.2.2 := by
  have hvV : v ∉ DfsOut.vertsList outs := fun h => ((hinv.new v).1 h).1 rfl
  have hancV : ∀ w ∈ anc, w ∉ DfsOut.vertsList outs := fun w hw hwV => by
    obtain ⟨_, _, _, hdw⟩ := hpre.anc_depth hw
    have := ((hinv.new w).1 hwV).2.1
    rw [hdw] at this
    exact nomatch this
  have heE : nxt ∉ DfsOut.vertsList outs → e ∉ DfsOut.edgesList outs := fun hnV he => by
    rcases ((hinv.edges_iff e).1 he).2 with ⟨y, hy⟩ | ⟨x, hx, y, hy⟩
    · exact hdone y hy
    · rcases hadj.ends _ _ _ _ _ hmem hy with h | h <;> subst h
      · exact hvV hx
      · exact hnV hx
  have hnxt_anc : depth[nxt]! ≠ none →
      ∃ i, i < d ∧ (anc ++ [v])[i]? = some nxt ∧ cur[nxt]! = some i := fun hvis => by
    obtain ⟨i, hi, hanc, hdi⟩ := hpre.anc_depth (hpre.nbr_in_anc hadj hmem hvis)
    exact ⟨i, hi, by rw [List.getElem?_append_left (hpre.hd ▸ hi)]; exact hanc, hinv.mono _ _ hdi⟩
  simp only [dfsStep]
  split
  · next hskip =>
    simp only [Bool.or_eq_true, beq_iff_eq, Option.any_eq_true, decide_eq_true_eq] at hskip
    have hvis : cur[nxt]! ≠ none := by
      rcases hskip with hp | ⟨k, hk, _⟩
      · obtain ⟨y, hy, hyv⟩ := hpre.parent e hp.symm
        rw [hadj.snd_inj hy hmem] at hyv
        obtain ⟨_, _, _, h⟩ := hnxt_anc hyv
        rw [h]; simp
      · rw [hk]; simp
    refine ⟨hinv.size, hinv.mono, hinv.hv, hinv.new, hinv.new_depth, hinv.nodup, hinv.closed, ?_, ?_,
      hinv.edges_nodup, hinv.wf, hinv.lv⟩
    · intro p hp
      rcases List.mem_append.1 hp with hp | hp
      · exact hinv.done_vis p hp
      · rw [List.mem_singleton.1 hp]; exact hvis
    · intro e'
      rw [hinv.edges_iff]
      constructor
      · rintro ⟨h1, ⟨y, hy⟩ | h2⟩
        · exact ⟨h1, .inl ⟨y, List.mem_append_left _ hy⟩⟩
        · exact ⟨h1, .inr h2⟩
      · rintro ⟨h1, ⟨y, hy⟩ | h2⟩
        · rcases List.mem_append.1 hy with hy | hy
          · exact ⟨h1, .inl ⟨y, hy⟩⟩
          · obtain ⟨hy1, hy2⟩ := Prod.mk.inj (List.mem_singleton.1 hy)
            rw [hy2] at h1 ⊢
            clear hy hy1 hy2
            rcases hskip with hp | ⟨k, hk, hkd⟩
            · exact absurd hp h1
            · have hnV : nxt ∈ DfsOut.vertsList outs := by
                refine (hinv.new nxt).2 ⟨?_, ?_, hvis⟩
                · rintro rfl
                  rw [hinv.hv] at hk
                  have := Option.some.inj hk
                  omega
                · by_contra hn
                  obtain ⟨i, hi, _, h⟩ := hnxt_anc hn
                  rw [h] at hk
                  have := Option.some.inj hk
                  omega
              exact ⟨h1, .inr ⟨nxt, hnV, v, hadj.symm _ _ _ hmem⟩⟩
        · exact ⟨h1, .inr h2⟩
  · next hskip =>
    simp only [Bool.or_eq_true, beq_iff_eq, Option.any_eq_true, decide_eq_true_eq, not_or, not_exists,
      not_and] at hskip
    obtain ⟨hne, hgt⟩ := hskip
    split
    · next hnone =>
      have hpre' : Pre g adj (anc ++ [v]) nxt (d + 1) (some e) cur := by
        refine ⟨hinv.size, by simp [hpre.hd], ?_, (hadj.bounds _ _ _ hmem).2.1, hnone, ?_, ?_⟩
        · intro i w hi
          by_cases hid : i < anc.length
          · rw [List.getElem?_append_left hid] at hi
            exact hinv.mono _ _ (hpre.path i w hi)
          · rw [List.getElem?_append_right (Nat.le_of_not_lt hid), List.getElem?_singleton] at hi
            split at hi
            · have hid' : i = d := by have := hpre.hd; omega
              rw [← Option.some.inj hi, hid']
              exact hinv.hv
            · exact nomatch hi
        · intro w hw hvis p hp
          simp only [List.mem_append, List.mem_singleton, not_or] at hw
          exact hinv.closed w hw.1 hw.2 hvis p hp
        · intro pe hpe
          cases hpe
          exact ⟨v, hadj.symm _ _ _ hmem, by rw [hinv.hv]; simp⟩
      have hpost := ih (anc ++ [v]) nxt (d + 1) (some e) cur (by omega) hpre'
      obtain ⟨⟨child, n, cur'⟩, hr⟩ : ∃ r, dfsVisit adj fuel nxt (d + 1) (some e) cur = r := ⟨_, rfl⟩
      rw [hr] at hpost ⊢
      simp only at hpost
      have hcm : ∀ w k : Nat, cur[w]! = some k → cur'[w]! = some k := hpost.mono
      have hcv : ∀ w : Nat, w ∈ child.verts ↔ cur[w]! = none ∧ cur'[w]! ≠ none := hpost.verts_iff
      have hchild_v : ∀ w ∈ child.verts, cur[w]! = none := fun w hw => ((hcv w).1 hw).1
      have hV_cur : ∀ w ∈ DfsOut.vertsList outs, cur[w]! ≠ none := fun w hw => ((hinv.new w).1 hw).2.2
      have hdisj : ∀ w, w ∈ child.verts → w ∉ DfsOut.vertsList outs := fun w hw hw' =>
        hV_cur w hw' (hchild_v w hw)
      have hvc : v ∉ child.verts := fun h => by
        have := hchild_v v h
        rw [hinv.hv] at this
        exact nomatch this
      have hnxtc : nxt ∈ child.verts := hpost.root ▸ child.v_mem_verts
      refine ⟨hpost.size, fun w k h => hcm w k (hinv.mono w k h), hcm v d hinv.hv, ?_, ?_, ?_, ?_, ?_, ?_,
        ?_, ?_, ?_⟩
      · intro w
        simp only [DfsOut.vertsList_tree, List.mem_append]
        constructor
        · rintro (hw | hw)
          · have ⟨h1, h2⟩ := (hcv w).1 hw
            refine ⟨?_, ?_, h2⟩
            · rintro rfl
              rw [hinv.hv] at h1
              exact nomatch h1
            · by_contra hn
              obtain ⟨k, hk⟩ := Option.ne_none_iff_exists'.1 hn
              rw [hinv.mono w k hk] at h1
              exact nomatch h1
          · obtain ⟨h1, h2, h3⟩ := (hinv.new w).1 hw
            obtain ⟨k, hk⟩ := Option.ne_none_iff_exists'.1 h3
            exact ⟨h1, h2, by rw [hcm w k hk]; simp⟩
        · rintro ⟨h1, h2, h3⟩
          by_cases hc : cur[w]! = none
          · exact .inl ((hcv w).2 ⟨hc, h3⟩)
          · exact .inr ((hinv.new w).2 ⟨h1, h2, hc⟩)
      · intro w hw
        simp only [DfsOut.vertsList_tree, List.mem_append] at hw
        rcases hw with hw | hw
        · obtain ⟨k, hk, hk'⟩ := hpost.verts_depth w hw
          exact ⟨k, by omega, hk'⟩
        · obtain ⟨k, hk, hk'⟩ := hinv.new_depth w hw
          exact ⟨k, hk, hcm w k hk'⟩
      · simp only [DfsOut.vertsList_tree]
        exact List.Nodup.append hpost.verts_nodup hinv.nodup hdisj
      · intro w hw hwv hvis p hp
        exact hpost.closed w (by simp [hw, hwv]) hvis p hp
      · intro p hp
        rcases List.mem_append.1 hp with hp | hp
        · obtain ⟨k, hk⟩ := Option.ne_none_iff_exists'.1 (hinv.done_vis p hp)
          rw [hcm _ k hk]
          simp
        · rw [List.mem_singleton.1 hp]
          exact ((hcv nxt).1 hnxtc).2
      · intro e'
        simp only [DfsOut.edgesList_tree, List.cons_append, List.mem_cons, List.mem_append,
          List.not_mem_nil, or_false, DfsOut.vertsList_tree]
        rw [hinv.edges_iff, hpost.edges_iff]
        constructor
        · rintro (rfl | ⟨_, x, hx, y, hy⟩ | ⟨h1, ⟨y, hy⟩ | ⟨x, hx, y, hy⟩⟩)
          · exact ⟨hne, .inl ⟨nxt, .inr rfl⟩⟩
          · refine ⟨?_, .inr ⟨x, .inl hx, y, hy⟩⟩
            intro hp
            obtain ⟨y0, hy0, hy0v⟩ := hpre.parent e' hp.symm
            rcases hadj.ends _ _ _ _ _ hy0 hy with h | h <;> subst h
            · exact hvc hx
            · obtain ⟨k, hk⟩ := Option.ne_none_iff_exists'.1 hy0v
              have := hchild_v _ hx
              rw [hinv.mono _ k hk] at this
              exact nomatch this
          · exact ⟨h1, .inl ⟨y, .inl hy⟩⟩
          · exact ⟨h1, .inr ⟨x, .inr hx, y, hy⟩⟩
        · rintro ⟨h1, ⟨y, hy | hy⟩ | ⟨x, hx | hx, y, hy⟩⟩
          · exact .inr (.inr ⟨h1, .inl ⟨y, hy⟩⟩)
          · exact .inl (Prod.mk.inj hy).2
          · by_cases hee : e' = e
            · exact .inl hee
            · exact .inr (.inl ⟨by simpa using hee, x, hx, y, hy⟩)
          · exact .inr (.inr ⟨h1, .inr ⟨x, hx, y, hy⟩⟩)
      · simp only [DfsOut.edgesList_tree, List.cons_append, List.nodup_cons, List.mem_append, not_or]
        refine ⟨⟨fun h => ((hpost.edges_iff e).1 h).1 rfl, heE (hdisj nxt hnxtc)⟩,
          List.Nodup.append hpost.edges_nodup hinv.edges_nodup ?_⟩
        intro e' h1 h2
        obtain ⟨_, x, hx, y, hy⟩ := (hpost.edges_iff e').1 h1
        rcases ((hinv.edges_iff e').1 h2).2 with ⟨z, hz⟩ | ⟨x', hx', y', hy'⟩
        · rcases hadj.ends _ _ _ _ _ (hsub _ hz) hy with h | h <;> subst h
          · exact hvc hx
          · exact hinv.done_vis _ hz (hchild_v _ hx)
        · rcases hadj.ends _ _ _ _ _ hy' hy with h | h <;> subst h
          · exact hV_cur _ hx' (hchild_v _ hx)
          · exact hinv.closed _ (fun ha => hancV _ ha hx') (fun hxv => hvV (hxv ▸ hx')) (hV_cur _ hx')
              _ hy' (hchild_v _ hx)
      · intro o ho
        rcases List.mem_cons.1 ho with rfl | ho
        · simp only [DfsOut.WF]
          refine ⟨hpost.wf, ?_⟩
          rw [hpost.lv, hpre.hd]
        · exact hinv.wf o ho
      · rw [DfsOut.retDepthsList_tree, hinv.lv, hpost.lv, mergeLowvals_low2 (Nat.le_succ d)]
        exact low2_perm List.perm_append_comm
    · next nd hsome =>
      have hnd : nd ≤ d := Nat.le_of_not_lt (hgt nd hsome)
      have hnV : nxt ∉ DfsOut.vertsList outs := fun h => by
        obtain ⟨k, hk, hk'⟩ := hinv.new_depth nxt h
        rw [hsome] at hk'
        have := Option.some.inj hk'
        omega
      have hidx : (anc ++ [v])[nd]? = some nxt := by
        by_cases hnv : nxt = v
        · subst hnv
          rw [hinv.hv] at hsome
          rw [← Option.some.inj hsome, ← hpre.hd]
          exact List.getElem?_concat_length
        · have hvis : depth[nxt]! ≠ none := fun h =>
            hnV ((hinv.new nxt).2 ⟨hnv, h, by rw [hsome]; simp⟩)
          obtain ⟨i, _, hi, h⟩ := hnxt_anc hvis
          rw [hsome] at h
          rw [Option.some.inj h]
          exact hi
      have hlow : (classify d false (nd, d)).lowval d = nd := lowval_classify_back hnd
      refine ⟨hinv.size, hinv.mono, hinv.hv, hinv.new, hinv.new_depth, hinv.nodup, hinv.closed, ?_, ?_,
        ?_, ?_, ?_⟩
      · intro p hp
        rcases List.mem_append.1 hp with hp | hp
        · exact hinv.done_vis p hp
        · rw [List.mem_singleton.1 hp, hsome]; simp
      · intro e'
        simp only [DfsOut.edgesList_back, List.mem_cons]
        rw [hinv.edges_iff]
        constructor
        · rintro (rfl | ⟨h1, ⟨y, hy⟩ | h2⟩)
          · exact ⟨hne, .inl ⟨nxt, List.mem_append_right _ (List.mem_singleton_self _)⟩⟩
          · exact ⟨h1, .inl ⟨y, List.mem_append_left _ hy⟩⟩
          · exact ⟨h1, .inr h2⟩
        · rintro ⟨h1, ⟨y, hy⟩ | h2⟩
          · rcases List.mem_append.1 hy with hy | hy
            · exact .inr ⟨h1, .inl ⟨y, hy⟩⟩
            · exact .inl (Prod.mk.inj (List.mem_singleton.1 hy)).2
          · exact .inr ⟨h1, .inr h2⟩
      · simp only [DfsOut.edgesList_back, List.nodup_cons]
        exact ⟨heE hnV, hinv.edges_nodup⟩
      · intro o ho
        rcases List.mem_cons.1 ho with rfl | ho
        · simp only [DfsOut.WF]
          exact ⟨nd, hidx, by rw [hpre.hd]⟩
        · exact hinv.wf o ho
      · rw [DfsOut.retDepthsList_back, hlow, hinv.lv]
        have h1 : (nd, d) = low2 d [nd] := by simp [low2, lmin, Nat.min_eq_left hnd]
        rw [h1, mergeLowvals_low2 (Nat.le_refl d)]
        exact low2_perm List.perm_append_comm

theorem foldInv_foldl (hadj : AdjOK g adj) (hpre : Pre g adj anc v d prvE depth) {fuel : Nat}
    (hfuel : g.nv ≤ fuel + 1 + d) (ih : VisitSpec g adj fuel) :
    ∀ (rest done : List (Nat × Nat)) (outs : List DfsOut) (lv : Lowvals) (cur : Array (Option Nat)),
      adj[v]! = done ++ rest → FoldInv g adj anc v d prvE depth done outs lv cur →
      let r := rest.foldl (dfsStep adj fuel d prvE) (outs, lv, cur)
      FoldInv g adj anc v d prvE depth (done ++ rest) r.1 r.2.1 r.2.2
  | [], done, outs, lv, cur, _, hinv => by simpa using hinv
  | (nxt, e) :: rest, done, outs, lv, cur, hsplit, hinv => by
    have hmem : (nxt, e) ∈ adj[v]! := by rw [hsplit]; simp
    have hdone : ∀ y, (y, e) ∉ done := fun y hy => by
      have hnd := hadj.nodup v
      rw [hsplit, List.map_append, List.map_cons] at hnd
      exact (List.nodup_cons.1 (List.nodup_middle.1 hnd)).1
        (List.mem_append_left _ (List.mem_map.2 ⟨(y, e), hy, rfl⟩))
    have hsub : ∀ p ∈ done, p ∈ adj[v]! := fun p hp => by
      rw [hsplit]; exact List.mem_append_left _ hp
    have hstep := foldInv_step hadj hpre hfuel ih hinv hmem hdone hsub
    have := foldInv_foldl hadj hpre hfuel ih rest (done ++ [(nxt, e)]) _ _ _ (by simpa using hsplit) hstep
    simpa using this

theorem dfsVisit_spec (hadj : AdjOK g adj) : ∀ fuel, VisitSpec g adj fuel
  | 0 => fun _ _ _ _ _ hfuel hpre => absurd hpre.succ_le (by omega)
  | fuel + 1 => fun anc v d prvE depth hfuel hpre => by
    show Post _ _ _ _ _ _ _ (dfsVisit adj (fuel + 1) v d prvE depth).1 _ _
    rw [dfsVisit_succ]
    exact post_of_foldInv hpre (foldInv_foldl hadj hpre hfuel (dfsVisit_spec hadj fuel) adj[v]! [] [] (d, d)
      _ rfl (foldInv_init hpre))

/-! ### The forest -/

/-- The body of the fold in `Graph.dfsForest`. -/
def forestStep (adj : Array (List (Nat × Nat))) (nv : Nat) :
    List DfsTree × Array (Option Nat) → Nat → List DfsTree × Array (Option Nat)
  | (roots, depth), rt =>
    if depth[rt]!.isSome then (roots, depth)
    else
      let (t, _, depth) := dfsVisit adj nv rt 0 none depth
      (t :: roots, depth)

theorem dfsForest_eq (g : Graph) (vo eo : List Nat) :
    g.dfsForest vo eo =
      ((inOrder g.nv vo).foldl (forestStep (g.adjacency eo) g.nv)
        ([], Array.replicate g.nv none)).1.reverse :=
  rfl

/-- Invariant of the fold in `Graph.dfsForest` after the root candidates `done`. -/
structure ForestInv (g : Graph) (adj : Array (List (Nat × Nat))) (done : List Nat)
    (roots : List DfsTree) (depth : Array (Option Nat)) : Prop where
  size : depth.size = g.nv
  done_vis : ∀ w ∈ done, depth[w]! ≠ none
  closed : ∀ w : Nat, depth[w]! ≠ none → ∀ p ∈ adj[w]!, depth[p.1]! ≠ none
  verts_iff : ∀ w : Nat, w ∈ roots.flatMap DfsTree.verts ↔ depth[w]! ≠ none
  verts_nodup : (roots.flatMap DfsTree.verts).Nodup
  edges_iff : ∀ e : Nat, e ∈ roots.flatMap DfsTree.edges ↔
    ∃ x ∈ roots.flatMap DfsTree.verts, ∃ y : Nat, (y, e) ∈ adj[x]!
  edges_nodup : (roots.flatMap DfsTree.edges).Nodup
  wf : ∀ t ∈ roots, t.WF []

theorem forestInv_nil (g : Graph) (adj : Array (List (Nat × Nat))) :
    ForestInv g adj [] [] (Array.replicate g.nv none) := by
  have h : ∀ w : Nat, (Array.replicate g.nv (none : Option Nat))[w]! = none := by
    intro w
    by_cases hw : w < g.nv
    · exact getElem!_replicate' _ _ _ hw
    · exact getElem!_neg _ _ (by simpa using hw)
  exact ⟨by simp, by simp, fun w hw => absurd (h w) hw, by simp [h], List.nodup_nil, by simp,
    List.nodup_nil, by simp⟩

theorem forestInv_step (hadj : AdjOK g adj) {done : List Nat} {roots : List DfsTree}
    {depth : Array (Option Nat)} (hinv : ForestInv g adj done roots depth) {rt : Nat} (hrt : rt < g.nv) :
    let r := forestStep adj g.nv (roots, depth) rt
    ForestInv g adj (done ++ [rt]) r.1 r.2 := by
  simp only [forestStep]
  split
  · next hsome =>
    refine ⟨hinv.size, ?_, hinv.closed, hinv.verts_iff, hinv.verts_nodup, hinv.edges_iff,
      hinv.edges_nodup, hinv.wf⟩
    intro w hw
    rcases List.mem_append.1 hw with hw | hw
    · exact hinv.done_vis w hw
    · rw [List.mem_singleton.1 hw]; exact Option.isSome_iff_ne_none.1 hsome
  · next hnone =>
    have hnone : depth[rt]! = none := Option.not_isSome_iff_eq_none.1 hnone
    have hpre : Pre g adj [] rt 0 none depth :=
      ⟨hinv.size, rfl, fun i w h => by simp at h, hrt, hnone, fun w _ => hinv.closed w,
        fun pe h => nomatch h⟩
    have hpost := dfsVisit_spec hadj g.nv [] rt 0 none depth (by simp) hpre
    obtain ⟨⟨t, n, depth'⟩, hr⟩ : ∃ r, dfsVisit adj g.nv rt 0 none depth = r := ⟨_, rfl⟩
    rw [hr] at hpost ⊢
    simp only at hpost
    have hmono : ∀ w : Nat, depth[w]! ≠ none → depth'[w]! ≠ none := fun w hw => by
      obtain ⟨k, hk⟩ := Option.ne_none_iff_exists'.1 hw
      rw [hpost.mono w k hk]; simp
    have htv : ∀ w ∈ t.verts, depth[w]! = none := fun w hw => ((hpost.verts_iff w).1 hw).1
    have hdisj : ∀ w, w ∈ t.verts → w ∉ roots.flatMap DfsTree.verts := fun w hw hw' =>
      (hinv.verts_iff w).1 hw' (htv w hw)
    refine ⟨hpost.size, ?_, fun w hw => hpost.closed w (by simp) hw, ?_, ?_, ?_, ?_, ?_⟩
    · intro w hw
      rcases List.mem_append.1 hw with hw | hw
      · exact hmono w (hinv.done_vis w hw)
      · rw [List.mem_singleton.1 hw]
        exact ((hpost.verts_iff rt).1 (hpost.root ▸ t.v_mem_verts)).2
    · intro w
      simp only [List.flatMap_cons, List.mem_append]
      constructor
      · rintro (hw | hw)
        · exact ((hpost.verts_iff w).1 hw).2
        · exact hmono w ((hinv.verts_iff w).1 hw)
      · intro hw
        by_cases hd : depth[w]! = none
        · exact .inl ((hpost.verts_iff w).2 ⟨hd, hw⟩)
        · exact .inr ((hinv.verts_iff w).2 hd)
    · simp only [List.flatMap_cons]
      exact List.Nodup.append hpost.verts_nodup hinv.verts_nodup hdisj
    · intro e
      simp only [List.flatMap_cons, List.mem_append]
      rw [hinv.edges_iff, hpost.edges_iff]
      constructor
      · rintro (⟨_, x, hx, y, hy⟩ | ⟨x, hx, y, hy⟩)
        · exact ⟨x, .inl hx, y, hy⟩
        · exact ⟨x, .inr hx, y, hy⟩
      · rintro ⟨x, hx | hx, y, hy⟩
        · exact .inl ⟨by simp, x, hx, y, hy⟩
        · exact .inr ⟨x, hx, y, hy⟩
    · simp only [List.flatMap_cons]
      refine List.Nodup.append hpost.edges_nodup hinv.edges_nodup ?_
      intro e h1 h2
      obtain ⟨_, x, hx, y, hy⟩ := (hpost.edges_iff e).1 h1
      obtain ⟨x', hx', y', hy'⟩ := (hinv.edges_iff e).1 h2
      rcases hadj.ends _ _ _ _ _ hy' hy with h | h <;> subst h
      · exact (hinv.verts_iff _).1 hx' (htv _ hx)
      · exact hinv.closed _ ((hinv.verts_iff _).1 hx') _ hy' (htv _ hx)
    · intro t' ht'
      rcases List.mem_cons.1 ht' with rfl | ht'
      · exact hpost.wf
      · exact hinv.wf t' ht'

theorem forestInv_foldl (hadj : AdjOK g adj) :
    ∀ (rest done : List Nat) (roots : List DfsTree) (depth : Array (Option Nat)),
      (∀ x ∈ rest, x < g.nv) → ForestInv g adj done roots depth →
      let r := rest.foldl (forestStep adj g.nv) (roots, depth)
      ForestInv g adj (done ++ rest) r.1 r.2
  | [], done, roots, depth, _, hinv => by simpa using hinv
  | rt :: rest, done, roots, depth, hlt, hinv => by
    have := forestInv_foldl hadj rest (done ++ [rt]) _ _ (fun x hx => hlt x (by simp [hx]))
      (forestInv_step hadj hinv (hlt rt (by simp)))
    simpa using this

/-- Lemma 1.1: every vertex appears exactly once in the forest, and every edge index exactly once
as an out-edge. -/
theorem dfsForest_spanning' (hg : g.WF) {vo eo : List Nat} (hvo : OrderOK g.nv vo)
    (heo : OrderOK g.ne eo) :
    ((g.dfsForest vo eo).flatMap DfsTree.verts).Perm (List.range g.nv) ∧
    ((g.dfsForest vo eo).flatMap DfsTree.edges).Perm (List.range g.ne) := by
  have hadj := adjacency_ok hg heo
  have hperm := inOrder_perm hvo
  have hinv := forestInv_foldl hadj (inOrder g.nv vo) [] [] _
    (fun x hx => List.mem_range.1 (hperm.subset hx)) (forestInv_nil g _)
  simp only [List.nil_append] at hinv
  rw [dfsForest_eq]
  refine ⟨((List.reverse_perm _).flatMap_right _).trans ?_,
    ((List.reverse_perm _).flatMap_right _).trans ?_⟩
  · refine (List.perm_ext_iff_of_nodup hinv.verts_nodup List.nodup_range).2 fun w => ?_
    rw [hinv.verts_iff, List.mem_range]
    constructor
    · intro h; exact hinv.size ▸ lt_size_of_getElem!_ne h
    · intro h; exact hinv.done_vis w (hperm.mem_iff.2 (List.mem_range.2 h))
  · refine (List.perm_ext_iff_of_nodup hinv.edges_nodup List.nodup_range).2 fun e => ?_
    rw [hinv.edges_iff, List.mem_range]
    constructor
    · rintro ⟨x, _, y, hy⟩; exact (hadj.bounds _ _ _ hy).2.2
    · intro h
      obtain ⟨x, y, hxy⟩ := hadj.cover e h
      exact ⟨x, (hinv.verts_iff x).2 (hinv.done_vis x
        (hperm.mem_iff.2 (List.mem_range.2 (hadj.bounds _ _ _ hxy).1))), y, hxy⟩

/-- `dfsForest`'s trees are well-formed: out-lists sorted by rank, back edges go to ancestors with
the right `classify`, tree children classified by their lowvals. -/
theorem dfsForest_wf (hg : g.WF) {vo eo : List Nat} (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) :
    ∀ t ∈ g.dfsForest vo eo, t.WF [] := by
  have hinv := forestInv_foldl (adjacency_ok hg heo) (inOrder g.nv vo) [] [] _
    (fun x hx => List.mem_range.1 ((inOrder_perm hvo).subset hx)) (forestInv_nil g _)
  rw [dfsForest_eq]
  intro t ht
  exact hinv.wf t (List.mem_reverse.1 ht)

/-! ### Lemma 1.2: lowvals and the `classify` table -/

theorem lmin_lt_iff {d k : Nat} {xs : List Nat} : lmin d xs < k ↔ d < k ∨ ∃ x ∈ xs, x < k := by
  constructor
  · intro h
    rcases lmin_eq_or_mem d xs with h' | h'
    · exact .inl (h' ▸ h)
    · exact .inr ⟨_, h', h⟩
  · rintro (h | ⟨x, hx, h⟩)
    · exact Nat.lt_of_le_of_lt (lmin_le d xs) h
    · exact Nat.lt_of_le_of_lt (lmin_le_of_mem hx) h

theorem le_lmin_iff {d k : Nat} {xs : List Nat} : k ≤ lmin d xs ↔ k ≤ d ∧ ∀ x ∈ xs, k ≤ x :=
  ⟨fun h => ⟨Nat.le_trans h (lmin_le d xs), fun _ hx => Nat.le_trans h (lmin_le_of_mem hx)⟩,
    fun h => le_lmin h.1 h.2⟩

@[simp] theorem low2_fst (d : Nat) (xs : List Nat) : (low2 d xs).1 = lmin d xs := rfl
theorem low2_snd (d : Nat) (xs : List Nat) : (low2 d xs).2 = lmin d (xs.filter (· ≠ lmin d xs)) := rfl

theorem classify_eq_bridge_iff {d : Nat} {isTree : Bool} {n : Lowvals} :
    classify d isTree n = .bridge ↔ d < n.1 := by
  unfold classify
  split
  · next h1 =>
    split
    · next h2 =>
      rw [beq_iff_eq] at h2
      cases isTree <;> simp <;> omega
    · next h2 =>
      rw [beq_iff_eq] at h2
      simp; omega
  · next h1 => simp; omega

theorem classify_eq_component_iff {d : Nat} {n : Lowvals} :
    classify d true n = .component ↔ n.1 = d := by
  unfold classify
  split
  · next h1 =>
    split
    · next h2 => rw [beq_iff_eq] at h2; simp [h2]
    · next h2 => rw [beq_iff_eq] at h2; simp [h2]
  · next h1 => simp; omega

theorem classify_eq_selfLoop_iff {d : Nat} {n : Lowvals} :
    classify d false n = .selfLoop ↔ n.1 = d := by
  unfold classify
  split
  · next h1 =>
    split
    · next h2 => rw [beq_iff_eq] at h2; simp [h2]
    · next h2 => rw [beq_iff_eq] at h2; simp [h2]
  · next h1 => simp; omega

theorem classify_eq_ret_iff {d l : Nat} {isTree : Bool} {k : RetKind} {n : Lowvals} :
    classify d isTree n = .ret l k ↔
      n.1 = l ∧ l < d ∧
        k = if isTree then (if n.2 < d then .type2Child else .type1Child) else .backEdge := by
  unfold classify
  split
  · next h1 =>
    split
    · cases isTree <;> simp <;> omega
    · simp; omega
  · next h1 =>
    simp only [OutClass.ret.injEq]
    constructor
    · rintro ⟨h3, h4⟩
      exact ⟨h3, by omega, h4.symm⟩
    · rintro ⟨h3, _, h4⟩
      exact ⟨h3, h4.symm⟩

theorem classify_eq_ret_iff_tree {d l : Nat} {k : RetKind} {n : Lowvals} :
    classify d true n = .ret l k ↔
      n.1 = l ∧ l < d ∧ k = if n.2 < d then .type2Child else .type1Child := by
  rw [classify_eq_ret_iff]; simp only [↓reduceIte]

theorem classify_eq_ret_iff_back {d l : Nat} {k : RetKind} {n : Lowvals} :
    classify d false n = .ret l k ↔ n.1 = l ∧ l < d ∧ k = .backEdge := by
  rw [classify_eq_ret_iff]; simp only [Bool.false_eq_true, ↓reduceIte]

/-! A tree child at depth `d + 1` whose back edges return to depths `ys`. -/

theorem classify_child_bridge_iff {d : Nat} {ys : List Nat} :
    classify d true (low2 (d + 1) ys) = .bridge ↔ ∀ y ∈ ys, d < y := by
  rw [classify_eq_bridge_iff, low2_fst, Nat.lt_iff_add_one_le, le_lmin_iff]
  constructor
  · intro h y hy; have := h.2 y hy; omega
  · intro h; exact ⟨Nat.le_refl _, fun y hy => h y hy⟩

theorem classify_child_component_iff {d : Nat} {ys : List Nat} :
    classify d true (low2 (d + 1) ys) = .component ↔ d ∈ ys ∧ ∀ y ∈ ys, d ≤ y := by
  rw [classify_eq_component_iff, low2_fst, lmin_eq_iff]
  constructor
  · rintro ⟨_, h2, h3⟩
    rcases h3 with h3 | h3
    · omega
    · exact ⟨h3, h2⟩
  · rintro ⟨h1, h2⟩
    exact ⟨Nat.le_succ d, h2, .inr h1⟩

theorem classify_child_type1_iff {d l : Nat} {ys : List Nat} :
    classify d true (low2 (d + 1) ys) = .ret l .type1Child ↔
      l ∈ ys ∧ l < d ∧ ∀ y ∈ ys, y = l ∨ d ≤ y := by
  rw [classify_eq_ret_iff_tree, low2_fst, low2_snd]
  constructor
  · rintro ⟨hl, hld, hk⟩
    split at hk
    · exact absurd hk (by decide)
    · next h2 =>
      refine ⟨?_, hld, ?_⟩
      · rcases lmin_eq_or_mem (d + 1) ys with h | h
        · omega
        · exact hl ▸ h
      · intro y hy
        by_cases hyl : y = l
        · exact .inl hyl
        · right
          have hmem : y ∈ ys.filter (· ≠ lmin (d + 1) ys) :=
            List.mem_filter.2 ⟨hy, by simpa [hl] using hyl⟩
          have := lmin_le_of_mem (d := d + 1) hmem
          omega
  · rintro ⟨h1, h2, h3⟩
    have hl : lmin (d + 1) ys = l := lmin_eq_iff.2 ⟨by omega, fun y hy => by
      rcases h3 y hy with rfl | h <;> omega, .inr h1⟩
    have hc : ¬ lmin (d + 1) (ys.filter (· ≠ lmin (d + 1) ys)) < d := by
      rw [Nat.not_lt, le_lmin_iff]
      refine ⟨by omega, fun y hy => ?_⟩
      rw [List.mem_filter, hl] at hy
      rcases h3 y hy.1 with h | h
      · simp [h] at hy
      · exact h
    exact ⟨hl, h2, by simp only [hc, ↓reduceIte]⟩

theorem classify_child_type2_iff {d l : Nat} {ys : List Nat} :
    classify d true (low2 (d + 1) ys) = .ret l .type2Child ↔
      l ∈ ys ∧ l < d ∧ (∀ y ∈ ys, l ≤ y) ∧ ∃ y ∈ ys, y ≠ l ∧ y < d := by
  rw [classify_eq_ret_iff_tree, low2_fst, low2_snd]
  constructor
  · rintro ⟨hl, hld, hk⟩
    split at hk
    · next h2 =>
      refine ⟨?_, hld, fun y hy => hl ▸ lmin_le_of_mem hy, ?_⟩
      · rcases lmin_eq_or_mem (d + 1) ys with h | h
        · omega
        · exact hl ▸ h
      · rcases lmin_lt_iff.1 h2 with h | ⟨y, hy, hyd⟩
        · omega
        · rw [List.mem_filter, hl] at hy
          exact ⟨y, hy.1, by simpa using hy.2, hyd⟩
    · exact absurd hk (by decide)
  · rintro ⟨h1, h2, h3, y, hy, hyl, hyd⟩
    have hl : lmin (d + 1) ys = l := lmin_eq_iff.2 ⟨by omega, h3, .inr h1⟩
    have hmem : y ∈ ys.filter (· ≠ lmin (d + 1) ys) :=
      List.mem_filter.2 ⟨hy, by simpa [hl] using hyl⟩
    have hc : lmin (d + 1) (ys.filter (· ≠ lmin (d + 1) ys)) < d :=
      lmin_lt_iff.2 (.inr ⟨y, hmem, hyd⟩)
    exact ⟨hl, h2, by simp only [hc, ↓reduceIte]⟩

/-- A back edge from depth `d` to depth `nd ≤ d`. -/
theorem classify_back_eq {d nd : Nat} (h : nd ≤ d) :
    classify d false (nd, d) = if nd = d then .selfLoop else .ret nd .backEdge := by
  by_cases hnd : nd = d
  · subst hnd; simp [classify]
  · simp [classify, hnd, Nat.not_le.2 (Nat.lt_of_le_of_ne h hnd)]

end Spqr
