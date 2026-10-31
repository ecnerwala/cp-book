import Spqr.Proofs.RSkel3
import Spqr.Proofs.RInvTree
import Spqr.RangesCloseTree

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

theorem vMerged_edges {g : Graph} {items : Items} {curV : Nat} {edgeDir : Bool} {c py vy : TEntry}
    {e : Nat} :
    (vMerged curV edgeDir c py vy).edges g items e ↔ ∃ i ∈ vItems c py vy, items.EdgeBelow g i e := by
  cases edgeDir <;> simp [TEntry.edges, vMerged, setSides]

theorem mem_vPieceItems_iff {c py vy : TEntry} {i : ItemId} :
    i ∈ s.vPieceItems c py vy ↔
      i ∈ s.entryPieceItems c ∨ i ∈ s.entryPieceItems py ∨ i ∈ s.entryPieceItems vy := by
  simp only [vPieceItems, vItems, entryPieceItems, List.mem_filter, List.mem_append]
  tauto

theorem vPieceItems_perm {c py vy : TEntry} :
    (s.vPieceItems c py vy).Perm
      (s.entryPieceItems c ++ s.entryPieceItems py ++ s.entryPieceItems vy) := by
  simp only [vPieceItems, vItems, entryPieceItems, ← List.filter_append]
  apply List.Perm.filter
  rw [List.perm_iff_count]; intro a; simp only [List.count_append]; omega

/-- The merged items of three pairwise edge-disjoint settled entries are `PieceItems`. -/
theorem vPieceItems_pieceItems {c py vy : TEntry} {rest : List TEntry}
    (hts : s.tstack = c :: py :: vy :: rest)
    (hdisj : s.tstack.Pairwise fun t t' => ∀ e, e < s.g.ne →
      t.edges s.g s.items e → ¬t'.edges s.g s.items e)
    (hc : s.EntryR dfs c) (hp : s.EntryR dfs py) (hy : s.EntryR dfs vy) :
    s.PieceItems (s.vPieceItems c py vy) := by
  rw [hts] at hdisj
  obtain ⟨h1, h2⟩ := List.pairwise_cons.1 hdisj
  obtain ⟨h3, -⟩ := List.pairwise_cons.1 h2
  have dcp := h1 py (by simp)
  have dcy := h1 vy (by simp)
  have dpy := h3 vy (by simp)
  have hc' := hc.pieces
  have hp' := hp.pieces
  have hy' := hy.pieces
  have both : ∀ {t t' : TEntry},
      (∀ e, e < s.g.ne → t.edges s.g s.items e → ¬t'.edges s.g s.items e) →
      s.PieceItems (s.entryPieceItems t) →
      ∀ i, i ∈ s.entryPieceItems t → i ∈ s.entryPieceItems t' → False := by
    intro t t' hd ht i hi hi'
    obtain ⟨e, he, hb⟩ := ht.ne i hi
    exact hd e he (edges_of_entryPieceItems hi hb) (edges_of_entryPieceItems hi' hb)
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_⟩
  · refine vPieceItems_perm.nodup_iff.2 (List.nodup_append.2
      ⟨List.nodup_append.2 ⟨hc'.nodup, hp'.nodup, ?_⟩, hy'.nodup, ?_⟩)
    · intro i hi j hj hij; subst hij; exact both dcp hc' i hi hj
    · intro i hi j hj hij; subst hij
      rcases List.mem_append.1 hi with hi | hi
      · exact both dcy hc' i hi hj
      · exact both dpy hp' i hi hj
  · intro i hi
    rcases mem_vPieceItems_iff.1 hi with hi | hi | hi
    exacts [hc'.vs i hi, hp'.vs i hi, hy'.vs i hi]
  · intro i hi
    rcases mem_vPieceItems_iff.1 hi with hi | hi | hi
    exacts [hc'.ne i hi, hp'.ne i hi, hy'.ne i hi]
  · intro i hi
    rcases mem_vPieceItems_iff.1 hi with hi | hi | hi
    exacts [hc'.conn i hi, hp'.conn i hi, hy'.conn i hi]
  · intro i hi
    rcases mem_vPieceItems_iff.1 hi with hi | hi | hi
    exacts [hc'.attached i hi, hp'.attached i hi, hy'.attached i hi]
  · intro i hi j hj hij e he hei hej
    have cross : ∀ {t t' : TEntry},
        (∀ e, e < s.g.ne → t.edges s.g s.items e → ¬t'.edges s.g s.items e) →
        i ∈ s.entryPieceItems t → j ∈ s.entryPieceItems t' → False := fun hd hi hj =>
      hd e he (edges_of_entryPieceItems hi hei) (edges_of_entryPieceItems hj hej)
    have cross' : ∀ {t t' : TEntry},
        (∀ e, e < s.g.ne → t.edges s.g s.items e → ¬t'.edges s.g s.items e) →
        j ∈ s.entryPieceItems t → i ∈ s.entryPieceItems t' → False := fun hd hj hi =>
      hd e he (edges_of_entryPieceItems hj hej) (edges_of_entryPieceItems hi hei)
    rcases mem_vPieceItems_iff.1 hi with hi | hi | hi <;>
      rcases mem_vPieceItems_iff.1 hj with hj | hj | hj
    · exact hc'.disj i hi j hj hij e he hei hej
    · exact cross dcp hi hj
    · exact cross dcy hi hj
    · exact cross' dcp hj hi
    · exact hp'.disj i hi j hj hij e he hei hej
    · exact cross dpy hi hj
    · exact cross' dcy hj hi
    · exact cross' dpy hj hi
    · exact hy'.disj i hi j hj hij e he hei hej

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
    (hshape : FinishRShape dfs curV d o origTstack true s)
    (h2 : s.g.TwoConnected) (hd : dfs.depth curV = d)
    (h1 : o.cls.isType1 = true) (hsingle : feSingle d o s = false)
    {base : List TEntry} {c py vy : TEntry}
    (hC : EarClose curV d lv o true base s (feS₂ d o s) c [] py vy)
    {σ : List Nat} {n : Nat} (hcb : CloseBase σ n D curV d o origTstack true s)
    (hc : CloseContent curV d o origTstack true s) :
    RCloseShape (feS₂ d o s).g dfs
      (Pieces.ofItems (feS₂ d o s).g (feS₂ d o s).items ((feS₂ d o s).vPieceItems c py vy))
      ((feS₂ d o s).vU c py vy) curV (feS₂ d o s).stackVerts[lv]! := by
  have ht : o.cls.isTree = true := by
    rw [ho]; cases kind <;> first | rfl | exact absurd rfl hk
  have hlv : o.cls.lowval d = lv := by rw [ho]; rfl
  have hlow' : o.cls.lowval d < d := by omega
  have hts : (feS₂ d o s).tstack = c :: py :: vy :: base := by simpa using hC.tstack
  obtain ⟨hEc, hEpy, hEvy, hht⟩ := hshape.content.v_close rfl h1 hsingle c py vy base hts
  rw [hlv] at hht
  have hdisj := hC.disj
  rw [← hC.g] at hdisj
  have hpi := vPieceItems_pieceItems hts hdisj hEc hEpy hEvy
  obtain ⟨t, vc⟩ := closeCtx_v_site hcb hc ht hlow' rfl h1
  rw [hsingle, cvS₅_type1 curV s.stackDir[d]! origTstack hts] at vc
  obtain ⟨rest', hts'⟩ := vc.stack
  obtain rfl : t = vMerged curV s.stackDir[d]! c py vy := (List.cons.inj hts').1.symm
  have htop : (vMerged curV s.stackDir[d]! c py vy).topDepth = lv := by
    rw [vMerged_topDepth]
    have := hC.c_top; have := hC.py_top; have := hC.vy_top
    omega
  have hedges : ∀ e, (vMerged curV s.stackDir[d]! c py vy).edges (feS₂ d o s).g
      ((feS₂ d o s).items.push ⟨.R, (none, none), []⟩) e ↔ (feS₂ d o s).vU c py vy e := by
    intro e
    rw [TEntry.edges_congr (fun _ _ _ => Items.Below_push_nil _ rfl) e, vU_iff_vItems, vMerged_edges]
  have hproper : ∃ e, e < (feS₂ d o s).g.ne ∧ ¬(feS₂ d o s).vU c py vy e := by
    obtain ⟨e, he, -, hnot⟩ := vc.pend curV (.inl rfl)
    exact ⟨e, he, fun hU => hnot ((hedges e).2 hU)⟩
  have hvok := (hcb.finishOk ho hlow).vert ht rfl
  have hvadj := (hcb.finishR.1 lv kind ho hlow).vert ht rfl
  rw [h1] at hvok hvadj
  have st₂ := rangesInv_feS₂ hcb ht hlow'
  have st₅ := RgStep.cvS₅ st₂.ranges st₂.step.shape hcb.nodup (st₂.hσ hcb.rgs.2.2)
    (by rw [st₂.step.g]; exact hcb.book.v_lt) hvok hvadj
  have hinv := st₅.step.inv
  rw [hsingle, cvS₅_type1 curV s.stackDir[d]! origTstack hts] at hinv
  have hconn := (hinv.entries [] _ base rfl).conn
  have hatt := (hC.type1 h1).2
  rw [← hC.g, ← hC.sv] at hatt
  have h2' : (feS₂ d o s).g.TwoConnected := by rw [hC.g]; exact h2
  obtain ⟨hwf, hsub⟩ := hpi.wf'
    (fun i hi e hb => vU_iff_vItems.2 ⟨i, List.mem_of_mem_filter hi, hb⟩) h2' hproper
  refine ⟨hwf, hsub, (Graph.ConnEdges.congr fun e _ => hedges e).1 hconn, ?_, ?_, ?_, ?_, hproper,
    hht.single, hht.maximal, hht.type1, hht.bond, hht.type2⟩
  · intro w e e' he he' hU hU' hw hw'
    rcases hatt w ⟨e, he, hU, hw⟩ with h | h | h
    · exact .inl h
    · exact .inr h
    · exact absurd (h e' he' hw') hU'
  · obtain ⟨e, he, hE, hw⟩ := vc.touch.1
    exact ⟨e, he, (hedges e).1 hE, hw⟩
  · have h := vc.touch.2
    rw [htop] at h
    obtain ⟨e, he, hE, hw⟩ := h
    exact ⟨e, he, (hedges e).1 hE, hw⟩
  · have h := vc.ne
    rw [htop] at h
    exact h

end WalkState

end Spqr
