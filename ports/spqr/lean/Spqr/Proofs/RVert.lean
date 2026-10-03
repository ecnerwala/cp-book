import Spqr.Proofs.RSkel3
import Spqr.Proofs.RInvTree

/-!
# The type-1 vertex close builds an `RSkel3` item

`closeVert' curV edgeDir true origTstack false` on a stack `c :: py :: vy :: rest` allocates an `R`
item, merges the three entries, retargets the result to `curV` and finishes it. The new item's
children are `vItems c py vy` (the concatenated spans), its edge set is `vU c py vy`, and its
terminals are `curV` and `stackVerts[topDepth]`. Given an `RCloseShape` for these data
(`Spqr.RCloseShape`, the same interface the loop-1 R branch uses), the item is `Items.RSkel3`
(`vClose_rSkel3`); this mirrors `RBranch.rSkel3`.
-/

namespace Spqr
open WalkM

namespace WalkState

variable {s : WalkState} {dfs : DfsData}

/-- The children of the type-1 vertex close: the spans of `c`, `py`, `vy` as the two
`mergeTstackTops` and the `retarget` concatenate them. -/
def vItems (c py vy : TEntry) : List ItemId :=
  ((c.spans.1 ++ py.spans.1) ++ vy.spans.1) ++ (vy.spans.2 ++ (py.spans.2 ++ c.spans.2))

/-- The edges closed by the type-1 vertex close. -/
def vU (s : WalkState) (c py vy : TEntry) (e : Nat) : Prop :=
  c.edges s.g s.items e ∨ py.edges s.g s.items e ∨ vy.edges s.g s.items e

/-- The sub-pieces of the type-1 vertex close: the merged items other than vertex items. -/
def vPieceItems (s : WalkState) (c py vy : TEntry) : List ItemId :=
  (vItems c py vy).filter fun i => decide (Items.type s.items i ≠ .V)

/-- The merged, retargeted top entry of the type-1 vertex close (`cvS₅`). -/
def vMerged (curV : Nat) (edgeDir : Bool) (c py vy : TEntry) : TEntry :=
  { TEntry.mergeInto (TEntry.mergeInto c py) vy with
    vStart := curV, spans := setSides (!edgeDir) (vItems c py vy) [] }

theorem vMerged_topDepth (curV : Nat) (edgeDir : Bool) (c py vy : TEntry) :
    (vMerged curV edgeDir c py vy).topDepth = min vy.topDepth (min py.topDepth c.topDepth) := rfl

theorem vU_iff_vItems {c py vy : TEntry} {e : Nat} :
    s.vU c py vy e ↔ ∃ i ∈ vItems c py vy, Items.EdgeBelow s.g s.items i e := by
  simp only [vU, TEntry.edges, vItems, List.mem_append]
  constructor
  · rintro (⟨i, hi, hb⟩ | ⟨i, hi, hb⟩ | ⟨i, hi, hb⟩) <;> exact ⟨i, by tauto, hb⟩
  · rintro ⟨i, hi, hb⟩
    rcases hi with ((h | h) | h) | (h | (h | h))
    · exact .inl ⟨i, .inl h, hb⟩
    · exact .inr (.inl ⟨i, .inl h, hb⟩)
    · exact .inr (.inr ⟨i, .inl h, hb⟩)
    · exact .inr (.inr ⟨i, .inr h, hb⟩)
    · exact .inr (.inl ⟨i, .inr h, hb⟩)
    · exact .inl ⟨i, .inr h, hb⟩

theorem getSide_vMerged {curV : Nat} {edgeDir dir : Bool} {c py vy : TEntry}
    (hside : getSide (vMerged curV edgeDir c py vy).spans (!dir) = []) :
    getSide (vMerged curV edgeDir c py vy).spans dir = vItems c py vy := by
  cases edgeDir <;> cases dir <;>
    simp only [vMerged, getSide, setSides, Bool.not_true, Bool.not_false, ↓reduceIte] at hside ⊢ <;>
    first | rfl | exact hside.symm

/-- The state after the merges and the retarget of a type-1, non-single vertex close. -/
theorem cvS₅_type1 (curV : Nat) (edgeDir : Bool) (origTstack : Nat) {c py vy : TEntry}
    {rest : List TEntry} (hts : s.tstack = c :: py :: vy :: rest) :
    cvS₅ curV edgeDir true origTstack false s =
      { s with items := s.items.push ⟨.R, (none, none), []⟩,
               tstack := vMerged curV edgeDir c py vy :: rest } := by
  have e2 := maybeUnwrapNxt_run_eq .R s c py (vy :: rest) hts _ rfl _ rfl
  simp only [true_or, ↓reduceIte, run_allocItem] at e2
  have h2 : cvS₂ true origTstack false s = { s with items := s.items.push ⟨.R, (none, none), []⟩ } := by
    unfold cvS₂ after
    rw [show cvB₁ true origTstack false s = false from rfl, show cvS₁ true origTstack false s = s from rfl]
    show ((maybeUnwrapNxt .R).run s).2 = _
    rw [e2]
  have h3 : cvS₃ true origTstack false s =
      { s with items := s.items.push ⟨.R, (none, none), []⟩,
               tstack := TEntry.mergeInto c py :: vy :: rest } := by
    unfold cvS₃ after
    rw [h2, mergeTstackTops_run_eq { s with items := s.items.push ⟨.R, (none, none), []⟩ } c py (vy :: rest) hts]
  have h4 : cvS₄ true origTstack false s =
      { s with items := s.items.push ⟨.R, (none, none), []⟩,
               tstack := TEntry.mergeInto (TEntry.mergeInto c py) vy :: rest } := by
    unfold cvS₄ after
    rw [h3, mergeTstackTops_run_eq _ (TEntry.mergeInto c py) vy rest rfl]
  unfold cvS₅ after
  rw [h4, retarget_run_eq curV edgeDir _ _ rest rfl]
  rfl

/-- `vClose_rSkel3`: the item built by a type-1, non-single vertex close on `c :: py :: vy :: rest`
is `Items.RSkel3`, given the `RCloseShape` of its data (children `vPieceItems`, edge set `vU`,
terminals `curV` and `stackVerts[topDepth]`) and the top-hole fact that the retargeted entry's
spans lie on the side `stackDir[topDepth]` (`FinishTopOk.side`). -/
theorem vClose_rSkel3 (hsh : Shape s) {curV : Nat} {edgeDir : Bool} {origTstack : Nat}
    {c py vy : TEntry} {rest : List TEntry} (hts : s.tstack = c :: py :: vy :: rest)
    (hside : getSide (vMerged curV edgeDir c py vy).spans
      (!s.stackDir[(vMerged curV edgeDir c py vy).topDepth]!) = [])
    (hR : RCloseShape s.g dfs (Pieces.ofItems s.g s.items (s.vPieceItems c py vy)) (s.vU c py vy)
      curV s.stackVerts[(vMerged curV edgeDir c py vy).topDepth]!)
    (hsp : dfs.Spec s.g) (hrt : dfs.Rooted s.g) (h2 : s.g.TwoConnected) :
    Items.RSkel3 s.g ((closeVert' curV edgeDir true origTstack false).run s).2.items
      s.items.size := by
  have key : ((closeVert' curV edgeDir true origTstack false).run s).2 =
      ((finishTstackTop ((maybeUnwrapNxt .R).run (cvS₁ true origTstack false s)).1).run
        (cvS₅ curV edgeDir true origTstack false s)).2 := rfl
  have e2 := maybeUnwrapNxt_run_eq .R s c py (vy :: rest) hts _ rfl _ rfl
  simp only [true_or, ↓reduceIte, run_allocItem] at e2
  have hitem : ((maybeUnwrapNxt .R).run (cvS₁ true origTstack false s)).1 = s.items.size := by
    rw [show cvS₁ true origTstack false s = s from rfl, e2]
  rw [key, hitem, cvS₅_type1 curV edgeDir origTstack hts,
    finishTstackTop_run_eq _ _ (vMerged curV edgeDir c py vy) rest rfl]
  dsimp only
  set n := s.items.size with hn
  set x0 : Item := ⟨.R, (none, none), []⟩ with hx0
  set top := (vMerged curV edgeDir c py vy).topDepth with htop
  have hL := getSide_vMerged hside
  set f : Item → Item := fun it =>
    { it with vs := setSides s.stackDir[top]! (some s.stackVerts[top]!) (some curV),
              ch := getSide (vMerged curV edgeDir c py vy).spans s.stackDir[top]! } with hf
  show Items.RSkel3 s.g ((s.items.push x0).modify n f) n
  set items' := (s.items.push x0).modify n f with hitems'
  have hmem : ∀ i ∈ vItems c py vy, i < n := by
    intro i hi
    have hc := hsh.span c (by simp [hts])
    have hp := hsh.span py (by simp [hts])
    have hv := hsh.span vy (by simp [hts])
    simp only [vItems, List.mem_append] at hi
    simp only [List.mem_append] at hc hp hv
    rcases hi with ((h | h) | h) | (h | (h | h))
    exacts [hc i (.inl h), hp i (.inl h), hv i (.inl h), hv i (.inr h), hp i (.inr h), hc i (.inr h)]
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
  have hvs_n : Items.vs items' n = setSides s.stackDir[top]! (some s.stackVerts[top]!) (some curV) := by
    rw [hitems', Items.vs_modify_at n f hnlt]
  have hch_n : Items.ch items' n = vItems c py vy := by
    rw [hitems', Items.ch_modify_at n f hnlt]
    exact hL
  have hEn : ∀ e, e < s.g.ne → (Items.EdgeBelow s.g items' n e ↔ s.vU c py vy e) := by
    intro e he
    rw [vU_iff_vItems]
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
      exact Relation.ReflTransGen.head (by rw [Items.IsParent, hch_n]; exact hc)
        ((hE c (hmem c hc) e).2 hb)
  have hlist : (Items.ch items' n).filter (fun c => decide (Items.type items' c ≠ .V)) =
      s.vPieceItems c py vy := by
    rw [hch_n]
    exact List.filter_congr fun c hc => by rw [hty c (hmem c hc)]
  have hvs' : ∀ i : Nat, Items.vs items' (s.vPieceItems c py vy)[i]! =
      Items.vs s.items (s.vPieceItems c py vy)[i]! := by
    intro i
    refine hvs _ ?_
    by_cases hi : i < (s.vPieceItems c py vy).length
    · rw [getElem!_pos (s.vPieceItems c py vy) i hi]
      exact Nat.ne_of_lt (hmem _ (List.mem_of_mem_filter (List.getElem_mem hi)))
    · rw [getElem!_neg (s.vPieceItems c py vy) i hi]
      have h1 := hsh.size
      show (0 : Nat) ≠ (n : Nat)
      omega
  refine ⟨curV, s.stackVerts[top]!, ?_, ?_⟩
  · rw [hvs_n]; cases s.stackDir[top]! <;> simp [setSides]
  · rw [hlist, Pieces.ofItems_congr (L := s.vPieceItems c py vy)
      (fun i hi e => hE i (hmem i (List.mem_of_mem_filter hi)) e) hvs',
      Pieces.addParent_congr _ _ _ _ hEn (fun e he => Pieces.ofItems_piece_none _ _ _ he)]
    exact hR.threeConnected hsp hrt h2

/-- Admitted (PROOF.md §4.5): the `RCloseShape` of the type-1 vertex close at a returning tree
edge with `feSingle = false`. At `feS₂` the stack is `c :: py :: vy :: base` (`EarClose`, `mid = []`
for type 1): `c` is the loop-1/2 output holding the child subtree's ear, `py` the single root item
of the `(o.dest, lv)` bottom, `vy` the vertex entry of `o.dest`; their union `vU` is the whole
sub-ear (`EarClose.sub_edges`/`sub_cover`), touching only `curV`, `stackVerts[lv]` and interior
vertices (`EarClose.type1`). Obligations: the non-V items are pairwise edge-disjoint 2-attached
pieces (`EntryR.pieces` of the settled entries, `RInvFront` at `feS₂`), the union is a proper
2-attached connected edge set at `(curV, stackVerts[lv])`, and the HT content (`single`,
`maximal`, `type1`, `bond`, `type2`) from run saturation and `RangesInv` interval ownership. -/
theorem closeVert_type1_rCloseShape {D : Nat} (curV d lv : Nat) (kind : RetKind) (o : DfsOut)
    (origTstack : Nat) (ho : o.cls = .ret lv kind) (hk : kind ≠ .backEdge) (hlow : lv < d)
    (hv : curV < s.g.nv) (hi : s.Inv' D) (hs : Shape s)
    (hok : FinishOk D curV d lv o origTstack true s)
    (hfront : Frontier (o := o) d origTstack s)
    (hshape : FinishRShape dfs curV d o origTstack true s)
    (h2 : s.g.TwoConnected) (hsp : dfs.Spec s.g) (hrt : dfs.Rooted s.g)
    (hd : dfs.depth curV = d) (hcur : s.stackVerts[d]! = curV)
    (hanc : ∀ k, k ≤ d → dfs.Anc s.stackVerts[k]! curV ∧ dfs.depth s.stackVerts[k]! = k)
    (hR : s.RInvFront dfs curV d origTstack)
    (h1 : o.cls.isType1 = true) (hsingle : feSingle d o s = false)
    {sub base : List TEntry} {c py vy : TEntry} (hbl : base.length = origTstack)
    (hE : s.EarFinish curV d o true sub base)
    (hC : EarClose curV d lv o true base s (feS₂ d o s) c [] py vy) :
    RCloseShape (feS₂ d o s).g dfs
      (Pieces.ofItems (feS₂ d o s).g (feS₂ d o s).items ((feS₂ d o s).vPieceItems c py vy))
      ((feS₂ d o s).vU c py vy) curV (feS₂ d o s).stackVerts[lv]! := by
  sorry

end WalkState

end Spqr
