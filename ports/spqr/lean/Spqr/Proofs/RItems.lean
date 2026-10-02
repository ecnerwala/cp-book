import Spqr.Proofs.RInvFrame
import Spqr.Correctness

/-!
# R items are 3-connected (PROOF.md §4.5, item level)

`Items.RSkel3 g items i`: the skeleton of the R item `i` (its non-`V` children as pieces, the
complement of its edge set as the parent piece at `i`'s terminals, contracted) is 3-connected.
`RBranch.rSkel3` establishes it for the item Loop 1's R branch closes (`rCloseItems`: `allocItem .R`,
then `finishTstackTop` of the merged entry), from `RBranch.threeConnected`. `items_r_three_connected`
(admitted) is the statement for every R item of the walk's output on a block.
-/

namespace Spqr

theorem findIdx?_congr_mem {α : Type} {l : List α} {p q : α → Bool} (h : ∀ x ∈ l, p x = q x) :
    l.findIdx? p = l.findIdx? q := by
  induction l with
  | nil => rfl
  | cons a l ih =>
    rw [List.findIdx?_cons, List.findIdx?_cons, h a (.head _), ih fun x hx => h x (.tail _ hx)]

theorem Items.type_modify_of_ne {items : Items} (j : ItemId) (f : Item → Item) {p : ItemId}
    (h : p ≠ j) : Items.type (items.modify j f) p = items.type p := by
  simp [Items.type, Array.getElem?_modify, h.symm]

theorem Pieces.ofItems_piece_none (g : Graph) (items : Items) (L : List ItemId) {e : Nat}
    (he : ¬ e < g.ne) : (Pieces.ofItems g items L).piece e = none := by
  simp [Pieces.ofItems, he]

open Classical in
theorem Pieces.ofItems_congr {g : Graph} {items items' : Items} {L : List ItemId}
    (hE : ∀ i ∈ L, ∀ e, items'.EdgeBelow g i e ↔ items.EdgeBelow g i e)
    (hvs : ∀ i : Nat, Items.vs items' L[i]! = Items.vs items L[i]!) :
    Pieces.ofItems g items' L = Pieces.ofItems g items L := by
  have hp : (fun e => if e < g.ne then L.findIdx? fun i => decide (items'.EdgeBelow g i e) else none) =
      fun e => if e < g.ne then L.findIdx? fun i => decide (items.EdgeBelow g i e) else none := by
    funext e
    split_ifs with he
    · exact findIdx?_congr_mem fun i hi => by rw [decide_eq_decide]; exact hE i hi e
    · rfl
  have hx : (fun i => (Items.vs items' L[i]!).1.getD 0) = fun i => (Items.vs items L[i]!).1.getD 0 :=
    funext fun i => by rw [hvs i]
  have hy : (fun i => (Items.vs items' L[i]!).2.getD 0) = fun i => (Items.vs items L[i]!).2.getD 0 :=
    funext fun i => by rw [hvs i]
  unfold Pieces.ofItems
  rw [hp, hx, hy]

open Classical in
theorem Pieces.addParent_congr (P : Pieces) (g : Graph) {U U' : Nat → Prop} (s t : Nat)
    (hU : ∀ e, e < g.ne → (U e ↔ U' e)) (hP : ∀ e, ¬ e < g.ne → P.piece e = none) :
    P.addParent g U s t = P.addParent g U' s t := by
  have hpiece : (fun e => if U e then P.piece e else if e < g.ne then some P.k else none) =
      fun e => if U' e then P.piece e else if e < g.ne then some P.k else none := by
    funext e
    by_cases he : e < g.ne
    · rw [hU e he]
    · simp [he, hP e he]
  unfold Pieces.addParent
  rw [hpiece]

/-- The skeleton of the R item `i` in the block `g` is 3-connected: its non-`V` children as
pieces, the complement of its edge set as the parent piece at `i`'s terminals, contracted. -/
def Items.RSkel3 (g : Graph) (items : Items) (i : ItemId) : Prop :=
  ∃ s t, (items.vs i = (some s, some t) ∨ items.vs i = (some t, some s)) ∧
    (((Pieces.ofItems g items ((items.ch i).filter fun c => decide (items.type c ≠ .V))).addParent g
      (items.EdgeBelow g i) s t).contract g).ThreeConnected

namespace WalkState

variable {s : WalkState} {d : Nat} {cur nxt : TEntry} {rest : List TEntry} {dfs : DfsData}

/-- The items after Loop 1's R close of `cur`, `nxt` into the fresh item `s.items.size`
(`allocItem .R`, `mergeTstackTops`, then `finishTstackTop` with the merged spans on side `dir`). -/
def rCloseItems (s : WalkState) (d : Nat) (cur nxt : TEntry) (dir : Bool) : Items :=
  (s.items.push ⟨.R, (none, none), []⟩).modify s.items.size fun it =>
    { it with vs := setSides dir (some s.stackVerts[d]!) (some nxt.vStart),
              ch := getSide (TEntry.mergeInto cur nxt).spans dir }

theorem rU_iff_rItems {e : Nat} :
    s.rU cur nxt e ↔ ∃ i ∈ rItems cur nxt, Items.EdgeBelow s.g s.items i e := by
  constructor
  · rintro (⟨i, hi, hb⟩ | ⟨i, hi, hb⟩)
    · exact ⟨i, by simp only [rItems, List.mem_append] at hi ⊢; tauto, hb⟩
    · exact ⟨i, by simp only [rItems, List.mem_append] at hi ⊢; tauto, hb⟩
  · rintro ⟨i, hi, hb⟩
    exact rU_of_mem_rItems hi hb

/-- The item closed at Loop 1's R branch has a 3-connected skeleton. -/
theorem RBranch.rSkel3 (hsh : Shape s) (hb : s.RBranch d cur nxt rest) (hR : s.RTop dfs cur nxt)
    (h : s.Inv (d + 1)) (h2 : s.g.TwoConnected) (hsp : dfs.Spec s.g) (hrt : dfs.Rooted s.g)
    (dir : Bool) (hside : getSide (TEntry.mergeInto cur nxt).spans (!dir) = []) :
    Items.RSkel3 s.g (s.rCloseItems d cur nxt dir) s.items.size := by
  unfold rCloseItems
  set n := s.items.size with hn
  set x0 : Item := ⟨.R, (none, none), []⟩ with hx0
  set f : Item → Item := fun it =>
    { it with vs := setSides dir (some s.stackVerts[d]!) (some nxt.vStart),
              ch := getSide (TEntry.mergeInto cur nxt).spans dir } with hf
  set items' := (s.items.push x0).modify n f with hitems'
  have hL : getSide (TEntry.mergeInto cur nxt).spans dir = rItems cur nxt := by
    cases dir
    · simp only [Bool.not_false, getSide, TEntry.mergeInto, Bool.false_eq_true, ↓reduceIte] at hside ⊢
      rw [rItems, hside, List.append_nil]
    · simp only [Bool.not_true, getSide, TEntry.mergeInto, Bool.false_eq_true, ↓reduceIte] at hside ⊢
      rw [rItems, hside, List.nil_append]
  have hmem : ∀ i ∈ rItems cur nxt, i < n := by
    intro i hi
    simp only [rItems, List.mem_append] at hi
    have hc := hsh.span cur (by simp [hb.tstack])
    have hn' := hsh.span nxt (by simp [hb.tstack])
    simp only [List.mem_append] at hc hn'
    rcases hi with (hi | hi) | (hi | hi)
    · exact hc i (.inl hi)
    · exact hn' i (.inl hi)
    · exact hn' i (.inr hi)
    · exact hc i (.inr hi)
  have hroot : ∀ p, ¬ Items.IsParent s.items p n := fun p hp => Nat.lt_irrefl _ (hsh.ch_lt p n hp)
  have hnlt : n < (s.items.push x0).size := by simp [hn]
  have hty : ∀ i, i < n → Items.type items' i = Items.type s.items i := fun i hi => by
    rw [hitems', Items.type_modify_of_ne n f (Nat.ne_of_lt hi),
      Items.type_push_of_ne x0 (Nat.ne_of_lt hi)]
  have hvs : ∀ i, i ≠ n → Items.vs items' i = Items.vs s.items i := fun i hi => by
    rw [hitems', Items.vs_modify_of_ne n f hi, Items.vs_push_of_ne x0 hi]
  have hE : ∀ i, i < n → ∀ e, Items.EdgeBelow s.g items' i e ↔ Items.EdgeBelow s.g s.items i e := by
    intro i hi e
    show Items.Below ((s.items.push x0).modify n f) i (edgeItem s.g e) ↔
      Items.Below s.items i (edgeItem s.g e)
    rw [Items.Below_modify_of_not_below n f, Items.Below_push_nil x0 rfl]
    intro hb
    rw [Items.Below_push_nil x0 rfl] at hb
    exact Nat.ne_of_lt hi (Items.Below.eq_of_no_parent hroot hb)
  have hvs_n : Items.vs items' n = setSides dir (some s.stackVerts[d]!) (some nxt.vStart) := by
    rw [hitems', Items.vs_modify_at n f hnlt]
  have hch_n : Items.ch items' n = rItems cur nxt := by
    rw [hitems', Items.ch_modify_at n f hnlt, ← hL]
  have hEn : ∀ e, e < s.g.ne → (Items.EdgeBelow s.g items' n e ↔ s.rU cur nxt e) := by
    intro e he
    rw [rU_iff_rItems]
    constructor
    · intro hb
      rcases Items.Below.head_cases hb with hb | ⟨c, hc, hb⟩
      · exfalso
        have h1 := hsh.size
        have h2 : (n : Nat) = 1 + s.g.nv + e := hb
        omega
      · rw [Items.IsParent, hch_n] at hc
        exact ⟨c, hc, (hE c (hmem c hc) e).1 hb⟩
    · rintro ⟨c, hc, hb⟩
      exact Relation.ReflTransGen.head (by rw [Items.IsParent, hch_n]; exact hc) ((hE c (hmem c hc) e).2 hb)
  have hlist : (Items.ch items' n).filter (fun c => decide (Items.type items' c ≠ .V)) =
      s.rPieceItems cur nxt := by
    rw [hch_n]
    exact List.filter_congr fun c hc => by rw [hty c (hmem c hc)]
  have hvs' : ∀ i : Nat, Items.vs items' (s.rPieceItems cur nxt)[i]! =
      Items.vs s.items (s.rPieceItems cur nxt)[i]! := by
    intro i
    refine hvs _ ?_
    by_cases hi : i < (s.rPieceItems cur nxt).length
    · rw [getElem!_pos (s.rPieceItems cur nxt) i hi]
      exact Nat.ne_of_lt (hmem _ (List.mem_of_mem_filter (List.getElem_mem hi)))
    · rw [getElem!_neg (s.rPieceItems cur nxt) i hi]
      have h1 := hsh.size
      show (0 : Nat) ≠ (n : Nat)
      omega
  refine ⟨nxt.vStart, s.stackVerts[d]!, ?_, ?_⟩
  · rw [hvs_n]; cases dir <;> simp [setSides]
  · rw [hlist, Pieces.ofItems_congr (L := s.rPieceItems cur nxt)
      (fun i hi e => hE i (hmem i (List.mem_of_mem_filter hi)) e) hvs',
      Pieces.addParent_congr _ _ _ _ hEn (fun e he => Pieces.ofItems_piece_none _ _ _ he)]
    exact hb.threeConnected h h2 hsp hrt hR

end WalkState

/-- Admitted (PROOF.md §4.5, item level): on a block, every R item of the walk's output has a
3-connected skeleton. Every R item is closed either at Loop 1's R branch (`RBranch.rSkel3` under
`loop1_rBranch`) or as the type-1 R close of `finishEdge` (the `isSingle = false` case, the same
argument with `cur` the type-1 entry); the skeleton is then untouched by later steps (items are
only modified when allocated or closed). -/
theorem items_r_three_connected (g : Graph) (tern : Bool) (vo eo : List Nat) (h2 : g.TwoConnected) :
    let items := (g.walk tern (g.dfsForest vo eo)).items
    ∀ i, i < items.size → Items.type items i = NodeType.R → Items.RSkel3 g items i := by
  sorry

end Spqr
