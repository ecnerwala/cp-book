import Spqr.Dfs
import Spqr.BucketSort

/-!
# Phase 1, linear time

`dfsVisit` sorts the out-edges of every vertex separately. `dfsForestFast` computes the same forest
with a single global stable bucket sort of all out-edges by rank, followed by a stable distribution
back to their nodes.
-/

namespace Spqr

abbrev DfsState := List DfsOut × Lowvals × Array (Option Nat)

/-- `dfsVisit` without the per-vertex sort: out-edges are in adjacency order. -/
def dfsVisitRaw (adj : Array (List (Nat × Nat))) :
    Nat → (v d : Nat) → Option Nat → Array (Option Nat) → DfsTree × Lowvals × Array (Option Nat)
  | 0, v, _, _, depth => (.node v [], (0, 0), depth)
  | fuel + 1, v, d, prvE, depth =>
    let depth := depth.set! v (some d)
    let (outs, lv, depth) := adj[v]!.foldl (init := (([] : List DfsOut), ((d, d) : Lowvals), depth))
      (dfsStepRaw adj fuel d prvE)
    (.node v outs.reverse, lv, depth)
where
  dfsStepRaw (adj : Array (List (Nat × Nat))) (fuel d : Nat) (prvE : Option Nat) :
      DfsState → Nat × Nat → DfsState
    | (outs, lv, depth), (nxt, e) =>
      if some e == prvE || (depth[nxt]!).any (· > d) then (outs, lv, depth)
      else match depth[nxt]! with
        | none =>
          let (child, n, depth) := dfsVisitRaw adj fuel nxt (d + 1) (some e) depth
          (.tree e (classify d true n) child :: outs, mergeLowvals lv n, depth)
        | some nd =>
          (.back e nxt (classify d false (nd, d)) :: outs, mergeLowvals lv (nd, d), depth)

def Graph.dfsForestRaw (g : Graph) (vertOrder edgeOrder : List Nat) : List DfsTree :=
  let adj := g.adjacency edgeOrder
  let (roots, _) := (inOrder g.nv vertOrder).foldl (init := (([] : List DfsTree), Array.replicate g.nv (none : Option Nat)))
    (rootStep adj g.nv)
  roots.reverse
where
  rootStep (adj : Array (List (Nat × Nat))) (fuel : Nat) : List DfsTree × Array (Option Nat) → Nat → List DfsTree × Array (Option Nat)
    | (roots, depth), rt =>
      if depth[rt]!.isSome then (roots, depth)
      else
        let (t, _, depth) := dfsVisitRaw adj fuel rt 0 none depth
        (t :: roots, depth)

mutual
/-- Sort the out-edges of every vertex by rank. -/
def DfsTree.sortOuts : DfsTree → DfsTree
  | .node v outs => .node v ((sortOutsList outs).mergeSort fun a b => a.cls.rank ≤ b.cls.rank)
def sortOutsList : List DfsOut → List DfsOut
  | [] => []
  | o :: rest => o.sortOuts :: sortOutsList rest
def DfsOut.sortOuts : DfsOut → DfsOut
  | .back e dest cls => .back e dest cls
  | .tree e cls child => .tree e cls child.sortOuts
end

theorem sortOutsList_eq_map (l : List DfsOut) : sortOutsList l = l.map DfsOut.sortOuts := by
  induction l with
  | nil => rfl
  | cons o rest ih => rw [sortOutsList, ih, List.map_cons]

private def liftSort (s : DfsState) : DfsState := (s.1.map DfsOut.sortOuts, s.2.1, s.2.2)

theorem dfsVisit_eq_raw (adj : Array (List (Nat × Nat))) (fuel : Nat) :
    ∀ (v d : Nat) (prvE : Option Nat) (depth : Array (Option Nat)),
      dfsVisit adj fuel v d prvE depth =
        ((dfsVisitRaw adj fuel v d prvE depth).1.sortOuts, (dfsVisitRaw adj fuel v d prvE depth).2.1,
          (dfsVisitRaw adj fuel v d prvE depth).2.2) := by
  induction fuel with
  | zero => intros; simp [dfsVisit, dfsVisitRaw, DfsTree.sortOuts, sortOutsList]
  | succ fuel ih =>
    intro v d prvE depth
    rw [dfsVisit, dfsVisitRaw]
    rw (occs := .pos [1]) [show (([] : List DfsOut), ((d, d) : Lowvals), depth.set! v (some d)) =
      liftSort ([], (d, d), depth.set! v (some d)) from rfl]
    rw [List.foldl_hom liftSort (g₁ := dfsVisitRaw.dfsStepRaw adj fuel d prvE)]
    · simp [liftSort, DfsTree.sortOuts, sortOutsList_eq_map, List.map_reverse]
    · intro s x
      obtain ⟨outs, lv, dep⟩ := s
      obtain ⟨nxt, e⟩ := x
      unfold liftSort dfsVisitRaw.dfsStepRaw
      dsimp only
      by_cases hc : (some e == prvE || (dep[nxt]!).any (· > d)) = true
      · simp only [hc, ↓reduceIte]
      · simp only [hc]
        cases dep[nxt]! with
        | none => simp [ih, DfsOut.sortOuts]
        | some nd => simp [DfsOut.sortOuts]

private def liftForest (s : List DfsTree × Array (Option Nat)) : List DfsTree × Array (Option Nat) :=
  (s.1.map DfsTree.sortOuts, s.2)

theorem Graph.dfsForest_eq_raw (g : Graph) (vertOrder edgeOrder : List Nat) :
    g.dfsForest vertOrder edgeOrder = (g.dfsForestRaw vertOrder edgeOrder).map DfsTree.sortOuts := by
  rw [Graph.dfsForest, Graph.dfsForestRaw]
  rw (occs := .pos [1]) [show (([] : List DfsTree), Array.replicate g.nv (none : Option Nat)) =
    liftForest ([], Array.replicate g.nv none) from rfl]
  rw [List.foldl_hom liftForest (g₁ := Graph.dfsForestRaw.rootStep (g.adjacency edgeOrder) g.nv)]
  · simp [liftForest, List.map_reverse]
  · intro s rt
    obtain ⟨roots, depth⟩ := s
    unfold liftForest Graph.dfsForestRaw.rootStep
    dsimp only
    by_cases hc : depth[rt]!.isSome = true
    · simp only [hc, ↓reduceIte]
    · simp [hc, dfsVisit_eq_raw]

/-- `(parent id, out-edge, child id)` in a preorder numbering of the forest; the child id is
unused for back edges. -/
abbrev Pair := Nat × DfsOut × Nat

def Pair.rank (p : Pair) : Nat := p.2.1.cls.rank

mutual
/-- Number the nodes in preorder starting from `n`; emit each node's out-edges, newest first. -/
def DfsTree.flat : DfsTree → Nat → List Pair → Nat × List Pair
  | .node _ outs, n, acc => flatOuts outs n (n + 1) acc
def flatOuts : List DfsOut → Nat → Nat → List Pair → Nat × List Pair
  | [], _, n, acc => (n, acc)
  | o :: rest, p, n, acc =>
    let (n', acc') := flatOut o p n acc
    flatOuts rest p n' acc'
def flatOut : DfsOut → Nat → Nat → List Pair → Nat × List Pair
  | .back e dest cls, p, n, acc => (n, (p, .back e dest cls, 0) :: acc)
  | .tree e cls child, p, n, acc => child.flat n ((p, .tree e cls child, n) :: acc)
end

/-- Flatten a forest, returning each root with its id. -/
def flatForest : List DfsTree → Nat → List Pair → List (DfsTree × Nat) × Nat × List Pair
  | [], n, acc => ([], n, acc)
  | t :: ts, n, acc =>
    let (n', acc') := t.flat n acc
    let (tagged, n'', acc'') := flatForest ts n' acc'
    ((t, n) :: tagged, n'', acc'')

/-- Rebuild node `n`'s subtree with the out-edge lists in `tbl`; `fuel` bounds the height. -/
def DfsTree.reorder (tbl : Array (List (DfsOut × Nat))) : Nat → Nat → DfsTree → DfsTree
  | 0, _, t => t
  | fuel + 1, n, .node v _ => .node v (tbl[n]!.map (reorderOut tbl fuel))
where
  reorderOut (tbl : Array (List (DfsOut × Nat))) (fuel : Nat) : DfsOut × Nat → DfsOut
    | (.tree e cls child, cid) => .tree e cls (child.reorder tbl fuel cid)
    | (o, _) => o

mutual
def DfsTree.height : DfsTree → Nat
  | .node _ outs => heightOuts outs + 1
def heightOuts : List DfsOut → Nat
  | [] => 0
  | o :: rest => max o.height (heightOuts rest)
def DfsOut.height : DfsOut → Nat
  | .back .. => 0
  | .tree _ _ child => child.height
end

def maxKey (key : α → Nat) (l : List α) : Nat := l.foldl (fun b x => max b (key x + 1)) 0

/-- `dfsForest`, with one global bucket sort of all out-edges by rank instead of a sort per vertex. -/
def Graph.dfsForestFast (g : Graph) (vertOrder edgeOrder : List Nat) : List DfsTree :=
  let (tagged, cnt, acc) := flatForest (g.dfsForestRaw vertOrder edgeOrder) 0 []
  let pairs := acc.reverse
  let sorted := bucketSort Pair.rank (maxKey Pair.rank pairs) pairs
  let tbl := (bucketLists (·.1) cnt sorted).map (·.map (·.2))
  tagged.map fun x => x.1.reorder tbl (g.nv + 1) x.2

section flatSpec

mutual
def DfsTree.fnext : DfsTree → Nat → Nat
  | .node _ outs, n => fnextOuts outs (n + 1)
def fnextOuts : List DfsOut → Nat → Nat
  | [], n => n
  | o :: rest, n => fnextOuts rest (o.fnext n)
def DfsOut.fnext : DfsOut → Nat → Nat
  | .back .., n => n
  | .tree _ _ child, n => child.fnext n
end

mutual
def DfsTree.flist : DfsTree → Nat → List Pair
  | .node _ outs, n => flistOuts outs n (n + 1)
def flistOuts : List DfsOut → Nat → Nat → List Pair
  | [], _, _ => []
  | o :: rest, p, n => flistOut o p n ++ flistOuts rest p (o.fnext n)
def flistOut : DfsOut → Nat → Nat → List Pair
  | .back e dest cls, p, _ => [(p, .back e dest cls, 0)]
  | .tree e cls child, p, n => (p, .tree e cls child, n) :: child.flist n
end

mutual
theorem DfsTree.flat_eq : ∀ (t : DfsTree) (n : Nat) (acc : List Pair),
    t.flat n acc = (t.fnext n, (t.flist n).reverse ++ acc)
  | .node _ outs, n, acc => by rw [DfsTree.flat, flatOuts_eq, DfsTree.fnext, DfsTree.flist]
theorem flatOuts_eq : ∀ (outs : List DfsOut) (p n : Nat) (acc : List Pair),
    flatOuts outs p n acc = (fnextOuts outs n, (flistOuts outs p n).reverse ++ acc)
  | [], _, _, _ => rfl
  | o :: rest, p, n, acc => by
    rw [flatOuts, flatOut_eq]
    dsimp only
    rw [flatOuts_eq, fnextOuts, flistOuts, List.reverse_append, List.append_assoc]
theorem flatOut_eq : ∀ (o : DfsOut) (p n : Nat) (acc : List Pair),
    flatOut o p n acc = (o.fnext n, (flistOut o p n).reverse ++ acc)
  | .back .., _, _, _ => rfl
  | .tree _ _ child, p, n, acc => by
    rw [flatOut, DfsTree.flat_eq, DfsOut.fnext, flistOut, List.reverse_cons, List.append_assoc,
      List.singleton_append]
end

def fnextForest : List DfsTree → Nat → Nat
  | [], n => n
  | t :: ts, n => fnextForest ts (t.fnext n)
def flistForest : List DfsTree → Nat → List Pair
  | [], _ => []
  | t :: ts, n => t.flist n ++ flistForest ts (t.fnext n)
def tagRoots : List DfsTree → Nat → List (DfsTree × Nat)
  | [], _ => []
  | t :: ts, n => (t, n) :: tagRoots ts (t.fnext n)

theorem flatForest_eq : ∀ (ts : List DfsTree) (n : Nat) (acc : List Pair),
    flatForest ts n acc = (tagRoots ts n, fnextForest ts n, (flistForest ts n).reverse ++ acc)
  | [], _, _ => rfl
  | t :: ts, n, acc => by
    rw [flatForest, DfsTree.flat_eq]
    dsimp only
    rw [flatForest_eq]
    dsimp only
    rw [tagRoots, fnextForest, flistForest, List.reverse_append, List.append_assoc]

mutual
theorem DfsTree.le_fnext : ∀ (t : DfsTree) (n : Nat), n < t.fnext n
  | .node _ outs, n => by rw [DfsTree.fnext]; exact Nat.lt_of_lt_of_le (Nat.lt_succ_self n) (le_fnextOuts outs _)
theorem le_fnextOuts : ∀ (outs : List DfsOut) (n : Nat), n ≤ fnextOuts outs n
  | [], _ => Nat.le_refl _
  | o :: rest, n => by rw [fnextOuts]; exact Nat.le_trans (DfsOut.le_fnext o n) (le_fnextOuts rest _)
theorem DfsOut.le_fnext : ∀ (o : DfsOut) (n : Nat), n ≤ o.fnext n
  | .back .., _ => Nat.le_refl _
  | .tree _ _ child, n => Nat.le_of_lt (child.le_fnext n)
end

theorem le_fnextForest : ∀ (ts : List DfsTree) (n : Nat), n ≤ fnextForest ts n
  | [], _ => Nat.le_refl _
  | t :: ts, n => Nat.le_trans (Nat.le_of_lt (t.le_fnext n)) (le_fnextForest ts _)

mutual
theorem DfsTree.id_mem_flist : ∀ (t : DfsTree) (n : Nat), ∀ q ∈ t.flist n, n ≤ q.1 ∧ q.1 < t.fnext n
  | .node _ outs, n, q, hq => by
    rw [DfsTree.flist] at hq
    rw [DfsTree.fnext]
    rcases id_mem_flistOuts outs n (n + 1) q hq with h | h
    · exact ⟨Nat.le_of_eq h.symm, Nat.lt_of_lt_of_le (h ▸ Nat.lt_succ_self n) (le_fnextOuts outs _)⟩
    · exact ⟨Nat.le_of_lt h.1, h.2⟩
theorem id_mem_flistOuts : ∀ (outs : List DfsOut) (p n : Nat), ∀ q ∈ flistOuts outs p n,
    q.1 = p ∨ (n ≤ q.1 ∧ q.1 < fnextOuts outs n)
  | [], _, _, _, hq => by simp [flistOuts] at hq
  | o :: rest, p, n, q, hq => by
    rw [flistOuts, List.mem_append] at hq
    rw [fnextOuts]
    rcases hq with hq | hq
    · rcases id_mem_flistOut o p n q hq with h | h
      · exact Or.inl h
      · exact Or.inr ⟨h.1, Nat.lt_of_lt_of_le h.2 (le_fnextOuts rest _)⟩
    · rcases id_mem_flistOuts rest p _ q hq with h | h
      · exact Or.inl h
      · exact Or.inr ⟨Nat.le_trans (o.le_fnext n) h.1, h.2⟩
theorem id_mem_flistOut : ∀ (o : DfsOut) (p n : Nat), ∀ q ∈ flistOut o p n,
    q.1 = p ∨ (n ≤ q.1 ∧ q.1 < o.fnext n)
  | .back .., p, n, q, hq => by
    rw [flistOut, List.mem_singleton] at hq
    exact Or.inl (hq ▸ rfl)
  | .tree _ _ child, p, n, q, hq => by
    rw [flistOut, List.mem_cons] at hq
    rw [DfsOut.fnext]
    rcases hq with hq | hq
    · exact Or.inl (hq ▸ rfl)
    · exact Or.inr (child.id_mem_flist n q hq)
end

theorem id_mem_flistForest : ∀ (ts : List DfsTree) (n : Nat), ∀ q ∈ flistForest ts n,
    n ≤ q.1 ∧ q.1 < fnextForest ts n
  | [], _, _, hq => by simp [flistForest] at hq
  | t :: ts, n, q, hq => by
    rw [flistForest, List.mem_append] at hq
    rw [fnextForest]
    rcases hq with hq | hq
    · have := t.id_mem_flist n q hq
      exact ⟨this.1, Nat.lt_of_lt_of_le this.2 (le_fnextForest ts _)⟩
    · have := id_mem_flistForest ts _ q hq
      exact ⟨Nat.le_trans (Nat.le_of_lt (t.le_fnext n)) this.1, this.2⟩

end flatSpec

section heightSpec

theorem height_le_heightOuts : ∀ (outs : List DfsOut), ∀ o ∈ outs, o.height ≤ heightOuts outs
  | o' :: rest, o, h => by
    rw [heightOuts]
    rcases List.mem_cons.1 h with rfl | h
    · exact Nat.le_max_left _ _
    · exact Nat.le_trans (height_le_heightOuts rest o h) (Nat.le_max_right _ _)

theorem heightOuts_le : ∀ (outs : List DfsOut) (k : Nat), (∀ o ∈ outs, o.height ≤ k) → heightOuts outs ≤ k
  | [], _, _ => Nat.zero_le _
  | o :: rest, k, h => by
    rw [heightOuts]
    exact Nat.max_le.2 ⟨h o (List.mem_cons_self ..), heightOuts_le rest k fun o' h' => h o' (List.mem_cons_of_mem _ h')⟩

theorem height_dfsVisitRaw (adj : Array (List (Nat × Nat))) (fuel : Nat) :
    ∀ (v d : Nat) (prvE : Option Nat) (depth : Array (Option Nat)),
      (dfsVisitRaw adj fuel v d prvE depth).1.height ≤ fuel + 1 := by
  induction fuel with
  | zero => intros; simp [dfsVisitRaw, DfsTree.height, heightOuts]
  | succ fuel ih =>
    intro v d prvE depth
    rw [dfsVisitRaw]
    dsimp only
    rw [DfsTree.height]
    refine Nat.succ_le_succ (heightOuts_le _ _ fun o ho => ?_)
    rw [List.mem_reverse] at ho
    have key : ∀ (l : List (Nat × Nat)) (init : DfsState), (∀ o ∈ init.1, o.height ≤ fuel + 1) →
        ∀ o ∈ (l.foldl (dfsVisitRaw.dfsStepRaw adj fuel d prvE) init).1, o.height ≤ fuel + 1 := by
      intro l
      induction l with
      | nil => intro init h; exact h
      | cons x l ihl =>
        intro init h
        rw [List.foldl_cons]
        refine ihl _ ?_
        obtain ⟨outs, lv, dep⟩ := init
        obtain ⟨nxt, e⟩ := x
        unfold dfsVisitRaw.dfsStepRaw
        dsimp only
        split
        · exact h
        · split
          · intro o ho
            rcases List.mem_cons.1 ho with rfl | ho
            · exact ih _ _ _ _
            · exact h o ho
          · intro o ho
            rcases List.mem_cons.1 ho with rfl | ho
            · exact Nat.zero_le _
            · exact h o ho
    exact key _ _ (by simp) o ho

end heightSpec

section reorderSpec

def pairSorted (P : List Pair) (m : Nat) : List (DfsOut × Nat) :=
  ((P.filter (·.1 == m)).mergeSort fun a b => Pair.rank a ≤ Pair.rank b).map (·.2)

/-- `tbl[m]` holds node `m`'s out-edges in rank order, for every id `m ∈ [lo, hi)`. -/
def TblOk (tbl : Array (List (DfsOut × Nat))) (P : List Pair) (lo hi : Nat) : Prop :=
  ∀ m, lo ≤ m → m < hi → tbl[m]! = pairSorted P m

theorem TblOk.mono {tbl : Array (List (DfsOut × Nat))} {P : List Pair} {lo hi lo' hi' : Nat}
    (h : TblOk tbl P lo hi) (hlo : lo ≤ lo') (hhi : hi' ≤ hi) : TblOk tbl P lo' hi' :=
  fun m h1 h2 => h m (Nat.le_trans hlo h1) (Nat.lt_of_lt_of_le h2 hhi)

theorem TblOk.congr {tbl : Array (List (DfsOut × Nat))} {P P' : List Pair} {lo hi : Nat} (h : TblOk tbl P lo hi)
    (hf : ∀ m, lo ≤ m → m < hi → P.filter (·.1 == m) = P'.filter (·.1 == m)) : TblOk tbl P' lo hi :=
  fun m h1 h2 => by rw [h m h1 h2]; unfold pairSorted; rw [hf m h1 h2]

theorem filter_id_eq_nil {P : List Pair} {m : Nat} (h : ∀ q ∈ P, q.1 ≠ m) : P.filter (·.1 == m) = [] :=
  List.filter_eq_nil_iff.2 fun q hq => by simpa using h q hq

theorem cls_reorderOut (tbl : Array (List (DfsOut × Nat))) (fuel : Nat) (q : DfsOut × Nat) :
    (DfsTree.reorder.reorderOut tbl fuel q).cls = q.1.cls := by
  obtain ⟨o, cid⟩ := q
  cases o <;> simp [DfsTree.reorder.reorderOut, DfsOut.cls]

theorem reorderOuts_eq (tbl : Array (List (DfsOut × Nat))) (fuel : Nat)
    (ih : ∀ (t : DfsTree) (n : Nat), t.height ≤ fuel → TblOk tbl (t.flist n) n (t.fnext n) →
      t.reorder tbl fuel n = t.sortOuts) :
    ∀ (outs : List DfsOut) (p n : Nat), p < n → (∀ o ∈ outs, o.height ≤ fuel) →
      TblOk tbl (flistOuts outs p n) n (fnextOuts outs n) →
      ((flistOuts outs p n).filter (·.1 == p)).map (fun q => DfsTree.reorder.reorderOut tbl fuel q.2) =
        outs.map DfsOut.sortOuts
  | [], _, _, _, _, _ => rfl
  | o :: rest, p, n, hp, hh, htbl => by
    rw [flistOuts, fnextOuts] at htbl
    rw [flistOuts, List.filter_append, List.map_append, List.map_cons]
    rw [reorderOuts_eq tbl fuel ih rest p (o.fnext n) (Nat.lt_of_lt_of_le hp (o.le_fnext n))
      (fun o' h' => hh o' (List.mem_cons_of_mem _ h'))
      ((htbl.mono (o.le_fnext n) (Nat.le_refl _)).congr fun m h1 _ => by
        rw [List.filter_append, filter_id_eq_nil (P := flistOut o p n) fun q hq => ?_, List.nil_append]
        have := o.le_fnext n
        rcases id_mem_flistOut o p n q hq with h | h <;> omega)]
    congr 1
    cases o with
    | back e dest cls => simp [flistOut, DfsTree.reorder.reorderOut, DfsOut.sortOuts]
    | tree e cls child =>
      rw [flistOut, List.filter_cons_of_pos (by simp), filter_id_eq_nil fun q hq => ?_, List.map_singleton]
      · simp only [DfsTree.reorder.reorderOut, DfsOut.sortOuts]
        rw [ih child n (hh _ (List.mem_cons_self ..))
          ((htbl.mono (Nat.le_refl n) (le_fnextOuts rest _)).congr fun m h1 h2 => by
            rw [flistOut, List.filter_append, List.filter_cons_of_neg (by simp; omega),
              filter_id_eq_nil (P := flistOuts rest p _) fun q hq => ?_, List.append_nil]
            rcases id_mem_flistOuts rest p _ q hq with h | h
            · omega
            · have := h.1; simp only [DfsOut.fnext] at this; omega)]
        exact List.singleton_append
      · have := (child.id_mem_flist n q hq).1; omega

theorem reorder_eq_sortOuts (tbl : Array (List (DfsOut × Nat))) :
    ∀ (fuel : Nat) (t : DfsTree) (n : Nat), t.height ≤ fuel → TblOk tbl (t.flist n) n (t.fnext n) →
      t.reorder tbl fuel n = t.sortOuts := by
  intro fuel
  induction fuel with
  | zero => intro t n h; cases t; simp [DfsTree.height] at h
  | succ fuel ih =>
    intro t n hh htbl
    obtain ⟨v, outs⟩ := t
    rw [DfsTree.reorder, DfsTree.sortOuts, sortOutsList_eq_map, htbl n (Nat.le_refl n) (DfsTree.le_fnext _ n)]
    unfold pairSorted
    rw [List.map_map, List.map_mergeSort (s := fun a b => a.cls.rank ≤ b.cls.rank)]
    · congr 2
      rw [DfsTree.flist]
      rw [DfsTree.flist, DfsTree.fnext] at htbl
      rw [DfsTree.height] at hh
      exact reorderOuts_eq tbl fuel ih outs n (n + 1) (Nat.lt_succ_self n)
        (fun o ho => Nat.le_trans (height_le_heightOuts outs o ho) (Nat.le_of_succ_le_succ hh))
        (htbl.mono (Nat.le_succ n) (Nat.le_refl _))
    · intro a _ b _
      simp only [Function.comp, cls_reorderOut]
      rfl

theorem tagRoots_reorder (tbl : Array (List (DfsOut × Nat))) (fuel : Nat) :
    ∀ (ts : List DfsTree) (n : Nat), (∀ t ∈ ts, t.height ≤ fuel) →
      TblOk tbl (flistForest ts n) n (fnextForest ts n) →
      (tagRoots ts n).map (fun x => x.1.reorder tbl fuel x.2) = ts.map DfsTree.sortOuts
  | [], _, _, _ => rfl
  | t :: ts, n, hh, htbl => by
    rw [flistForest, fnextForest] at htbl
    rw [tagRoots, List.map_cons, List.map_cons]
    congr 1
    · exact reorder_eq_sortOuts tbl fuel t n (hh t (List.mem_cons_self ..))
        ((htbl.mono (Nat.le_refl _) (le_fnextForest ts _)).congr fun m _ h2 => by
          rw [List.filter_append, filter_id_eq_nil (P := flistForest ts _) fun q hq => ?_, List.append_nil]
          have := (id_mem_flistForest ts _ q hq).1; omega)
    · exact tagRoots_reorder tbl fuel ts _ (fun t' h' => hh t' (List.mem_cons_of_mem _ h'))
        ((htbl.mono (Nat.le_of_lt (t.le_fnext n)) (Nat.le_refl _)).congr fun m h1 _ => by
          rw [List.filter_append, filter_id_eq_nil (P := t.flist n) fun q hq => ?_, List.nil_append]
          have := (t.id_mem_flist n q hq).2; omega)

theorem height_dfsForestRaw (g : Graph) (vertOrder edgeOrder : List Nat) :
    ∀ t ∈ g.dfsForestRaw vertOrder edgeOrder, t.height ≤ g.nv + 1 := by
  rw [Graph.dfsForestRaw]
  dsimp only
  intro t ht
  rw [List.mem_reverse] at ht
  have key : ∀ (l : List Nat) (init : List DfsTree × Array (Option Nat)), (∀ t ∈ init.1, t.height ≤ g.nv + 1) →
      ∀ t ∈ (l.foldl (Graph.dfsForestRaw.rootStep (g.adjacency edgeOrder) g.nv) init).1, t.height ≤ g.nv + 1 := by
    intro l
    induction l with
    | nil => intro init h; exact h
    | cons x l ihl =>
      intro init h
      rw [List.foldl_cons]
      refine ihl _ ?_
      obtain ⟨roots, depth⟩ := init
      unfold Graph.dfsForestRaw.rootStep
      dsimp only
      split
      · exact h
      · intro t ht
        rcases List.mem_cons.1 ht with rfl | ht
        · exact height_dfsVisitRaw _ _ _ _ _ _
        · exact h t ht
  exact key _ _ (by simp) t ht

theorem le_foldl_max (key : α → Nat) : ∀ (l : List α) (init : Nat), init ≤ l.foldl (fun b x => max b (key x + 1)) init
  | [], _ => Nat.le_refl _
  | _ :: l, _ => Nat.le_trans (Nat.le_max_left _ _) (le_foldl_max key l _)

theorem lt_maxKey (key : α → Nat) (l : List α) : ∀ x ∈ l, key x < maxKey key l := by
  unfold maxKey
  suffices ∀ init, ∀ x ∈ l, key x < l.foldl (fun b x => max b (key x + 1)) init from this 0
  induction l with
  | nil => simp
  | cons y l ih =>
    intro init x hx
    rw [List.foldl_cons]
    rcases List.mem_cons.1 hx with rfl | hx
    · exact Nat.lt_of_lt_of_le (Nat.lt_succ_self _) (Nat.le_trans (Nat.le_max_right _ _) (le_foldl_max key l _))
    · exact ih _ x hx

theorem filter_mergeSort_rank (P : List Pair) (m : Nat) :
    (P.mergeSort fun a b => Pair.rank a ≤ Pair.rank b).filter (·.1 == m) =
      (P.filter (·.1 == m)).mergeSort fun a b => Pair.rank a ≤ Pair.rank b := by
  refine eq_of_pairwise_of_filter_eq (key := Pair.rank) ?_ ?_ fun k => ?_
  · exact ((List.pairwise_mergeSort (decide_le_trans _) (decide_le_total _) P).filter _).imp fun h => by simpa using h
  · exact (List.pairwise_mergeSort (decide_le_trans _) (decide_le_total _) _).imp fun h => by simpa using h
  · rw [List.filter_filter, List.mergeSort_filter_of_pairwise _ _ (fun x y hx hy => ?_) (decide_le_trans _) (decide_le_total _),
      List.mergeSort_filter_of_pairwise _ _ (fun x y hx hy => ?_) (decide_le_trans _) (decide_le_total _),
      List.filter_filter]
    all_goals simp only [Bool.and_eq_true, beq_iff_eq] at hx hy; simp; omega

theorem Graph.dfsForestFast_eq (g : Graph) (vertOrder edgeOrder : List Nat) :
    g.dfsForestFast vertOrder edgeOrder = g.dfsForest vertOrder edgeOrder := by
  rw [Graph.dfsForest_eq_raw, Graph.dfsForestFast, flatForest_eq]
  dsimp only
  rw [List.append_nil, List.reverse_reverse, bucketSort_eq_mergeSort _ (lt_maxKey _ _)]
  apply tagRoots_reorder _ _ _ 0 (height_dfsForestRaw g vertOrder edgeOrder)
  intro m _ hm
  have hsize : m < ((bucketLists (·.1) (fnextForest (g.dfsForestRaw vertOrder edgeOrder) 0)
      ((flistForest (g.dfsForestRaw vertOrder edgeOrder) 0).mergeSort fun a b => Pair.rank a ≤ Pair.rank b)).map
        (·.map (·.2))).size := by
    simp [size_bucketLists]; omega
  rw [getElem!_pos (cont := Array (List (DfsOut × Nat))) _ _ hsize, Array.getElem_map, getElem_bucketLists]
  unfold pairSorted
  rw [← filter_mergeSort_rank]
  congr 1
  apply List.filter_congr
  intro q hq
  have := (id_mem_flistForest _ 0 q ((List.mergeSort_perm _ _).mem_iff.1 hq)).2
  simp [Nat.min_eq_left (Nat.le_of_lt this)]

end reorderSpec

end Spqr
