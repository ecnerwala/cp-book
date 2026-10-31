import Spqr.StTree
import Spqr.StRefEt

/-! Reference-side and bookkeeping lemmas for the walk simulation (`StInduct.lean`): `refOuts` over a
snoc, the pieces of a vertex-less prefix, `DirsOf` under depth / frame changes, `simBlocks` over a
snoc of path frames, and the segmented base stack `segsStack` / `SegRead`. -/

namespace Spqr

open WalkM WalkState StRefEt

theorem refOuts_snoc {g : Graph} {v d : Nat} {dirs : List Bool} :
    ∀ (done : List DfsOut) (o : DfsOut) (hv : Bool),
    refOuts g v d dirs (done ++ [o]) hv =
      ((refOuts g v d dirs done hv).1 ++ (refOut g v d dirs o (refOuts g v d dirs done hv).2.2).1,
       (refOuts g v d dirs done hv).2.1 ++ (refOut g v d dirs o (refOuts g v d dirs done hv).2.2).2.1,
       (refOut g v d dirs o (refOuts g v d dirs done hv).2.2).2.2)
  | [], o, hv => by simp [refOuts_cons, refOuts_nil]
  | o' :: done, o, hv => by
    simp only [List.cons_append, refOuts_cons, refOuts_snoc done o, List.append_assoc]

theorem refOut_hv_false {g : Graph} {v d : Nat} {dirs : List Bool} {o : DfsOut} {hv : Bool}
    (h : (refOut g v d dirs o hv).2.2 = false) : hv = false ∧ (refOut g v d dirs o hv).1 = [] := by
  cases o <;> cases hv <;> simp only [refOut] at h ⊢ <;> split_ifs at h <;> simp_all

theorem refOuts_hv_false {g : Graph} {v d : Nat} {dirs : List Bool} :
    ∀ (outs : List DfsOut) (hv : Bool), (refOuts g v d dirs outs hv).2.2 = false →
      hv = false ∧ (refOuts g v d dirs outs hv).1 = []
  | [], hv, h => by simpa [refOuts_nil] using h
  | o :: rest, hv, h => by
    rw [refOuts_cons] at h ⊢
    obtain ⟨h₁, h₂⟩ := refOuts_hv_false rest _ h
    obtain ⟨h₃, h₄⟩ := refOut_hv_false h₁
    subst h₃
    exact ⟨rfl, by simp [h₂, h₄]⟩

theorem DirsOf_length (s : WalkState) (d : Nat) : (DirsOf s d).length = d := by simp [DirsOf]

theorem DirsOf_succ (s : WalkState) (d : Nat) : DirsOf s (d + 1) = DirsOf s d ++ [s.stackDir[d]!] := by
  simp [DirsOf, List.range_succ]

theorem DirsOf_take (s : WalkState) {k d : Nat} (h : k ≤ d) : (DirsOf s d).take k = DirsOf s k := by
  simp [DirsOf, ← List.map_take, List.take_range, Nat.min_eq_left h]

theorem DirsOf_congr {s s' : WalkState} {d : Nat} (h : ∀ k, k < d → s'.stackDir[k]! = s.stackDir[k]!) :
    DirsOf s' d = DirsOf s d := by
  simp only [DirsOf]
  exact List.map_congr_left fun k hk => h k (List.mem_range.1 hk)

theorem frameBlocks_append (g : Graph) (dirs : List Bool) :
    ∀ (k : Nat) (fs fs' : List PathFrame),
    frameBlocks g dirs k (fs ++ fs') = frameBlocks g dirs k fs ++ frameBlocks g dirs (k + fs.length) fs'
  | k, [], fs' => by simp [frameBlocks]
  | k, f :: fs, fs' => by
    have e : k + 1 + fs.length = k + (fs.length + 1) := by omega
    simp only [List.cons_append, frameBlocks, frameBlocks_append g dirs (k + 1) fs fs', List.append_assoc,
      List.length_cons, e]

theorem frameBlocks_congr {g : Graph} {dirs dirs' : List Bool} :
    ∀ (k : Nat) (fs : List PathFrame), (∀ j, k ≤ j → j < k + fs.length → dirs.take j = dirs'.take j) →
    frameBlocks g dirs k fs = frameBlocks g dirs' k fs
  | _, [], _ => rfl
  | k, f :: fs, h => by
    simp only [frameBlocks]
    rw [h k le_rfl (by simp), frameBlocks_congr (k + 1) fs fun j hj hj' => h j (by omega) (by simp; omega)]

theorem simBlocks_snoc {g : Graph} {prev : List DfsTree} {fs : List PathFrame} {f : PathFrame}
    {dirs : List Bool} :
    simBlocks g prev (fs ++ [f]) dirs =
      simBlocks g prev fs dirs ++ (refOuts g f.v fs.length (dirs.take fs.length) f.done false).2.1 := by
  simp [simBlocks, frameBlocks_append, frameBlocks]

theorem simBlocks_congr {g : Graph} {prev : List DfsTree} {fs : List PathFrame} {dirs dirs' : List Bool}
    (h : ∀ j, j < fs.length → dirs.take j = dirs'.take j) :
    simBlocks g prev fs dirs = simBlocks g prev fs dirs' := by
  simp only [simBlocks]
  rw [frameBlocks_congr 0 fs fun j _ hj => h j (by simpa using hj)]

/-- The frames' blocks only look at the directions strictly below `fs.length`. -/
theorem simBlocks_DirsOf {g : Graph} {prev : List DfsTree} {fs : List PathFrame} {s s' : WalkState}
    {d d' : Nat} (hlen : d = fs.length) (hd : d ≤ d')
    (h : ∀ k, k < d → s'.stackDir[k]! = s.stackDir[k]!) :
    simBlocks g prev fs (DirsOf s' d') = simBlocks g prev fs (DirsOf s d) := by
  apply simBlocks_congr
  intro j hj
  rw [DirsOf_take s' (by omega), DirsOf_take s (by omega)]
  exact DirsOf_congr fun k hk => h k (by omega)

theorem DfsOut.vertsList_append :
    ∀ (a b : List DfsOut), DfsOut.vertsList (a ++ b) = DfsOut.vertsList a ++ DfsOut.vertsList b
  | [], _ => rfl
  | .back .. :: a, b => by simp [DfsOut.vertsList, vertsList_append a b]
  | .tree _ _ _ :: a, b => by simp [DfsOut.vertsList, vertsList_append a b]

theorem DfsOut.edgesList_append :
    ∀ (a b : List DfsOut), DfsOut.edgesList (a ++ b) = DfsOut.edgesList a ++ DfsOut.edgesList b
  | [], _ => rfl
  | .back .. :: a, b => by simp [DfsOut.edgesList, edgesList_append a b]
  | .tree _ _ _ :: a, b => by simp [DfsOut.edgesList, edgesList_append a b]

/-- The base stack below the current frame, as the segments pushed by the enclosing frames (top
first), each with the pieces it reads. -/
def segsStack (segs : List (List TEntry × List StPiece)) : List TEntry := (segs.map Prod.fst).flatten

def SegRead (items : Items) (segs : List (List TEntry × List StPiece)) : Prop :=
  ∀ sg ∈ segs, StRead items sg.1 sg.2

theorem segsStack_cons (sg : List TEntry × List StPiece) (segs : List (List TEntry × List StPiece)) :
    segsStack (sg :: segs) = sg.1 ++ segsStack segs := by simp [segsStack]

theorem mem_readStack_append {x : ItemId} {l l' : List TEntry} :
    x ∈ readStack (l ++ l') ↔ x ∈ readStack l ∨ x ∈ readStack l' := by
  simp only [readStack, readL_append, readR_append, List.mem_append]
  tauto

theorem mem_readStack_segsStack {x : ItemId} :
    ∀ {segs : List (List TEntry × List StPiece)},
    x ∈ readStack (segsStack segs) → ∃ sg ∈ segs, x ∈ readStack sg.1
  | [], h => by simp [segsStack, readStack, readL, readR] at h
  | sg :: segs, h => by
    rw [segsStack_cons, mem_readStack_append] at h
    rcases h with h | h
    · exact ⟨sg, List.mem_cons_self, h⟩
    · obtain ⟨sg', h₁, h₂⟩ := mem_readStack_segsStack h
      exact ⟨sg', List.mem_cons_of_mem _ h₁, h₂⟩

theorem SegRead.congr {items items' : Items} {segs : List (List TEntry × List StPiece)}
    (H : ∀ x ∈ readStack (segsStack segs), ∀ y, Items.Below items x y →
      Items.type items' y = Items.type items y ∧ Items.ch items' y = Items.ch items y)
    (h : SegRead items segs) : SegRead items' segs := by
  intro sg hsg
  refine (h sg hsg).congr fun x hx => H x ?_
  induction segs with
  | nil => exact (List.not_mem_nil hsg).elim
  | cons sg' segs ih =>
    rw [segsStack_cons, mem_readStack_append]
    rcases List.mem_cons.1 hsg with rfl | hsg'
    · exact Or.inl hx
    · exact Or.inr (ih (fun x hx => H x (by rw [segsStack_cons, mem_readStack_append]; exact Or.inr hx))
        (fun sg hs => h sg (List.mem_cons_of_mem _ hs)) hsg')

theorem SegRead.cons {items : Items} {new : List TEntry} {ps : List StPiece}
    {segs : List (List TEntry × List StPiece)} (h₁ : StRead items new ps) (h : SegRead items segs) :
    SegRead items ((new, ps) :: segs) := by
  intro sg hsg
  rcases List.mem_cons.1 hsg with rfl | hsg
  · exact h₁
  · exact h sg hsg

theorem StRead.nil {items : Items} : StRead items [] [] := ⟨.nil, .nil⟩

end Spqr
