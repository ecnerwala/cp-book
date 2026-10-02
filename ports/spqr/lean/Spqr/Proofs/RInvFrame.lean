import Spqr.Proofs.RInv
import Spqr.WalkInv
import Spqr.EarShape

/-!
# Frame lemmas for the R-maximality invariant, and its admitted preservation

`EntryR`, `RInv`, `RInvAt` read the state only through `g`, `stackVerts`, the tstack, and the
types, terminals and edge sets of the entries' items (`EntryR.congr`, `RInvAt.congr`). This
transports them across the bookkeeping steps of `finishEdge`: `modifyItem` of a free item
(`feS₀`), `setStackDir`, the `firstOccurrence`/`nxtEdgeIdx` update, and `pushVertTstack` (a vertex
entry has no edges and no pieces: `EntryR.vert`).

The content-carrying blocks are admitted: `finishEdge_rInvAt` (Loop 1's closes, Loop 2,
`closeVert`, the P-check, the back-edge push keep the invariant settled at `curV`),
`walkTree_rInvAt` (finishing a child settles its entries), and `loop1_rBranch` (the stack shape at
Loop 1's R branch). `loop1_r_threeConnected` combines the last with `RBranch.threeConnected`.
-/

namespace Spqr
open WalkM

namespace WalkState

variable {s s' : WalkState} {dfs : DfsData} {t : TEntry}

theorem mem_spans_of_entryPieceItems {i : ItemId} (hi : i ∈ s.entryPieceItems t) :
    i ∈ t.spans.1 ++ t.spans.2 :=
  (List.mem_filter.1 hi).1

theorem entryPieceItems_congr
    (hty : ∀ i ∈ t.spans.1 ++ t.spans.2, Items.type s'.items i = Items.type s.items i) :
    s'.entryPieceItems t = s.entryPieceItems t := by
  unfold entryPieceItems
  exact List.filter_congr fun i hi => by rw [hty i hi]

/-- `EntryR` depends on the state only through `g`, `stackVerts`, and the types, terminals and
edge sets of `t`'s items. -/
theorem EntryR.congr (hg : s'.g = s.g) (hsv : s'.stackVerts = s.stackVerts)
    (hty : ∀ i ∈ t.spans.1 ++ t.spans.2, Items.type s'.items i = Items.type s.items i)
    (hvs : ∀ i ∈ t.spans.1 ++ t.spans.2, Items.vs s'.items i = Items.vs s.items i)
    (hE : ∀ i ∈ t.spans.1 ++ t.spans.2, ∀ e,
      Items.EdgeBelow s.g s'.items i e ↔ Items.EdgeBelow s.g s.items i e)
    (h : s.EntryR dfs t) : s'.EntryR dfs t := by
  have hp : s'.entryPieceItems t = s.entryPieceItems t := entryPieceItems_congr hty
  have hed : t.edges s.g s'.items = t.edges s.g s.items :=
    funext fun e => propext (TEntry.edges_congr hE e)
  have hEb : ∀ i ∈ s.entryPieceItems t,
      Items.EdgeBelow s.g s'.items i = Items.EdgeBelow s.g s.items i :=
    fun i hi => funext fun e => propext (hE i (mem_spans_of_entryPieceItems hi) e)
  have hvs' : ∀ i ∈ s.entryPieceItems t, Items.vs s'.items i = Items.vs s.items i :=
    fun i hi => hvs i (mem_spans_of_entryPieceItems hi)
  have hskp : ∀ a b, s'.EntrySkelPair t a b ↔ s.EntrySkelPair t a b := by
    intro a b
    unfold EntrySkelPair
    rw [hp, hg, hsv]
    refine and_congr (forall₂_congr fun i hi => ?_) (and_congr (forall₂_congr fun i hi => ?_) Iff.rfl)
    · rw [hvs' i hi, hEb i hi]
    · rw [hvs' i hi]
  have hlam : ∀ K, s'.EntryLaminar t K ↔ s.EntryLaminar t K := by
    intro K
    unfold EntryLaminar
    rw [hp, hg, hed]
    refine or_congr (exists_congr fun i => ?_) Iff.rfl
    constructor
    · rintro ⟨hi, hK⟩; exact ⟨hi, by rw [← hEb i hi]; exact hK⟩
    · rintro ⟨hi, hK⟩; exact ⟨hi, by rw [hEb i hi]; exact hK⟩
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_⟩
  · rw [hp]
    obtain ⟨nodup, vs, ne, conn, att, disj⟩ := h.pieces
    refine ⟨nodup, fun i hi => ?_, fun i hi => ?_, fun i hi => ?_, fun i hi x y hxy => ?_,
      fun i hi j hj hij e => ?_⟩
    · rw [hvs' i hi]; exact vs i hi
    · rw [hg, hEb i hi]; exact ne i hi
    · rw [hg, hEb i hi]; exact conn i hi
    · rw [hvs' i hi] at hxy; rw [hg, hEb i hi]; exact att i hi x y hxy
    · rw [hg, hEb i hi, hEb j hj]; exact disj i hi j hj hij e
  · rw [hp]
    intro i hi x y hxy e e' he he' hn hn'
    rw [hvs' i hi] at hxy
    rw [hg] at he he' hn hn' ⊢
    rw [hEb i hi] at hn hn'
    exact h.maximal i hi x y hxy e e' he he' hn hn'
  · rw [hp, hg, hed]
    intro a b e₁ e₂ hne h1 h2 t1 t2
    obtain ⟨i, hi, b1, b2⟩ := h.bond a b e₁ e₂ hne h1 h2 t1 t2
    exact ⟨i, hi, by rw [hEb i hi]; exact b1, by rw [hEb i hi]; exact b2⟩
  · rw [hg, hsv, hed]; exact h.single
  · rw [hg]
    intro a b hsk hanc o ho hcls
    exact (hlam _).2 (h.type1 a b ((hskp a b).1 hsk) hanc o ho hcls)
  · rw [hg]
    intro a b hsk ht2 o ho hot hanc
    exact (hlam _).2 (h.type2 a b ((hskp a b).1 hsk) ht2 o ho hot hanc)

theorem RInvAt.congr {v : Nat} (hg : s'.g = s.g) (hsv : s'.stackVerts = s.stackVerts)
    (hts : s'.tstack = s.tstack)
    (hty : ∀ t ∈ s.tstack, ∀ i ∈ t.spans.1 ++ t.spans.2,
      Items.type s'.items i = Items.type s.items i)
    (hvs : ∀ t ∈ s.tstack, ∀ i ∈ t.spans.1 ++ t.spans.2, Items.vs s'.items i = Items.vs s.items i)
    (hE : ∀ t ∈ s.tstack, ∀ i ∈ t.spans.1 ++ t.spans.2, ∀ e,
      Items.EdgeBelow s.g s'.items i e ↔ Items.EdgeBelow s.g s.items i e)
    (h : s.RInvAt dfs v) : s'.RInvAt dfs v := by
  refine ⟨fun t ht hv => ?_, ?_⟩
  · rw [hts] at ht
    exact EntryR.congr hg hsv (hty t ht) (hvs t ht) (hE t ht) (h.entries t ht hv)
  · rw [hts, hg]
    exact h.disj.imp_of_mem fun {t t'} ht ht' hd e he => by
      rw [TEntry.edges_congr (hE t' ht')]
      rw [TEntry.edges_congr (hE t ht)] at he
      exact hd e he

theorem RInv.congr (hg : s'.g = s.g) (hsv : s'.stackVerts = s.stackVerts)
    (hts : s'.tstack = s.tstack)
    (hty : ∀ t ∈ s.tstack, ∀ i ∈ t.spans.1 ++ t.spans.2,
      Items.type s'.items i = Items.type s.items i)
    (hvs : ∀ t ∈ s.tstack, ∀ i ∈ t.spans.1 ++ t.spans.2, Items.vs s'.items i = Items.vs s.items i)
    (hE : ∀ t ∈ s.tstack, ∀ i ∈ t.spans.1 ++ t.spans.2, ∀ e,
      Items.EdgeBelow s.g s'.items i e ↔ Items.EdgeBelow s.g s.items i e)
    (h : s.RInv dfs) : s'.RInv dfs := by
  refine ⟨fun t ht => ?_, ?_⟩
  · rw [hts] at ht
    exact EntryR.congr hg hsv (hty t ht) (hvs t ht) (hE t ht) (h.entries t ht)
  · rw [hts, hg]
    exact h.disj.imp_of_mem fun {t t'} ht ht' hd e he => by
      rw [TEntry.edges_congr (hE t' ht')]
      rw [TEntry.edges_congr (hE t ht)] at he
      exact hd e he

theorem RInvAt.of_eq {v : Nat} (hg : s'.g = s.g) (hi : s'.items = s.items)
    (hsv : s'.stackVerts = s.stackVerts) (hts : s'.tstack = s.tstack) (h : s.RInvAt dfs v) :
    s'.RInvAt dfs v :=
  h.congr hg hsv hts (fun _ _ _ _ => by rw [hi]) (fun _ _ _ _ => by rw [hi])
    (fun _ _ _ _ _ => by rw [hi])

theorem RInv.of_eq (hg : s'.g = s.g) (hi : s'.items = s.items)
    (hsv : s'.stackVerts = s.stackVerts) (hts : s'.tstack = s.tstack) (h : s.RInv dfs) :
    s'.RInv dfs :=
  h.congr hg hsv hts (fun _ _ _ _ => by rw [hi]) (fun _ _ _ _ => by rw [hi])
    (fun _ _ _ _ _ => by rw [hi])

/-! ### Bookkeeping steps of `finishEdge` -/

theorem RInvAt.setStackDir {v d : Nat} {b : Bool} (h : s.RInvAt dfs v) :
    (after (setStackDir d b) s).RInvAt dfs v :=
  RInvAt.of_eq (s := s) (s' := after (WalkM.setStackDir d b) s) rfl rfl rfl rfl h

/-- `modify f` for an `f` that touches none of `g`, `items`, `stackVerts`, `tstack` (the
`firstOccurrence`/`nxtEdgeIdx` update of the back-edge branch). -/
theorem RInvAt.modify {v : Nat} (f : WalkState → WalkState) (hg : (f s).g = s.g)
    (hi : (f s).items = s.items) (hsv : (f s).stackVerts = s.stackVerts)
    (hts : (f s).tstack = s.tstack) (h : s.RInvAt dfs v) :
    (after (modify f : WalkM Unit) s).RInvAt dfs v :=
  RInvAt.of_eq (s := s) (s' := after (_root_.modify f : WalkM Unit) s) hg hi hsv hts h

/-- `modifyItem` of an item in no entry's spans and under no parent (`feS₀`: the tree edge's `Q`
item gets its terminals before it is pushed). -/
theorem RInvAt.modifyItem_free {v : Nat} (j : ItemId) (f : Item → Item)
    (hroot : ∀ p, ¬ Items.IsParent s.items p j)
    (hfree : ∀ t ∈ s.tstack, j ∉ t.spans.1 ++ t.spans.2)
    (h : s.RInvAt dfs v) : (after (modifyItem j f) s).RInvAt dfs v := by
  refine RInvAt.congr (s := s) (s' := after (modifyItem j f) s) rfl rfl rfl
    (fun t ht i hi => ?_) (fun t ht i hi => ?_) (fun t ht i hi e => ?_) h
  · have hne : i ≠ j := fun hij => hfree t ht (hij ▸ hi)
    show Items.type (s.items.modify j f) i = Items.type s.items i
    simp [Items.type, Array.getElem?_modify, hne.symm]
  · exact Items.vs_modify_of_ne j f fun hij => hfree t ht (hij ▸ hi)
  · exact Items.Below_modify_of_not_below j f fun hb =>
      hfree t ht ((Items.Below.eq_of_no_parent hroot hb) ▸ hi)

/-- A vertex entry (`pushVertTstack`) has no pieces and no edges, so `EntryR` holds vacuously. -/
theorem EntryR.vert {v d : Nat} (hv : v < s.g.nv) (hsh : Shape s)
    (hch : Items.ch s.items (vertItem v) = []) (t : TEntry)
    (hsp : t.spans = setSides s.stackDir[d]! [vertItem v] []) : s.EntryR dfs t := by
  have hmem : ∀ i ∈ t.spans.1 ++ t.spans.2, i = vertItem v := fun i hi => by
    rw [hsp] at hi; exact List.mem_singleton.1 ((mem_setSides _ _ _).1 hi)
  have hp : s.entryPieceItems t = [] := List.filter_eq_nil_iff.2 fun i hi => by
    rw [hmem i hi]; simp [hsh.vert v hv]
  have hed : ∀ e, ¬ t.edges s.g s.items e := by
    rintro e ⟨i, hi, hb⟩; rw [hmem i hi] at hb; exact edgeBelow_vert_nil hv hch e hb
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_⟩
  · rw [hp]; exact ⟨List.nodup_nil, by simp, by simp, by simp, by simp, by simp⟩
  · rw [hp]; simp
  · exact fun a b e₁ e₂ _ _ _ h1 _ => absurd h1 (hed e₁)
  · exact fun _ e e' _ _ h1 _ => absurd h1 (hed e)
  · exact fun a b _ _ o _ _ => .inr (.inl fun e _ he => hed e he)
  · exact fun a b _ _ o _ _ _ => .inr (.inl fun e _ he => hed e he)

theorem RInvAt.pushVert {v d w : Nat} (hv : v < s.g.nv) (hsh : Shape s)
    (hch : Items.ch s.items (vertItem v) = []) (h : s.RInvAt dfs w) :
    (after (pushVertTstack v d) s).RInvAt dfs w := by
  show RInvAt { s with tstack := ⟨v, d, s.nxtEdgeIdx, setSides s.stackDir[d]! [vertItem v] []⟩ :: s.tstack } dfs w
  have hcongr : ∀ t, s.EntryR dfs t →
      EntryR { s with tstack := ⟨v, d, s.nxtEdgeIdx, setSides s.stackDir[d]! [vertItem v] []⟩ :: s.tstack } dfs t :=
    fun t => EntryR.congr rfl rfl (fun _ _ => rfl) (fun _ _ => rfl) (fun _ _ _ => Iff.rfl)
  have hed : ∀ e, ¬ (⟨v, d, s.nxtEdgeIdx, setSides s.stackDir[d]! [vertItem v] []⟩ : TEntry).edges
      s.g s.items e := by
    rintro e ⟨i, hi, hb⟩
    rw [List.mem_singleton.1 ((mem_setSides _ _ _).1 hi)] at hb
    exact edgeBelow_vert_nil hv hch e hb
  refine ⟨fun t ht hne => ?_, List.pairwise_cons.2 ⟨fun t' _ e he => absurd he (hed e), h.disj⟩⟩
  rcases List.mem_cons.1 ht with rfl | ht
  · exact hcongr _ (EntryR.vert hv hsh hch _ rfl)
  · exact hcongr t (h.entries t ht hne)

/-! ### Iterates of Loop 1 -/

theorem iter_succ' (body : WalkM Unit) (k : Nat) (s : WalkState) :
    iter body (k + 1) s = (body.run (iter body k s)).2 := by
  induction k generalizing s with
  | zero => rfl
  | succ k ih => exact ih _

/-- `Inv`/`Shape` at every iterate of Loop 1 whose guards `CloseEarsOk.body` supplies. -/
theorem closeEars_iter_step {D v nxtV d e : Nat} {edgeDir : Bool} (hi : s.Inv D) (hs : Shape s)
    (hv : v < s.g.nv) (hok : CloseEarsOk D nxtV d e edgeDir s) (k : Nat)
    (hk : ∀ j, j ≤ k → result (loop1Cond d) (iter (loop1Body d edgeDir) j (ceS₁ nxtV d e s)) = true) :
    Step D v s (iter (loop1Body d edgeDir) k (ceS₁ nxtV d e s)) := by
  induction k with
  | zero => exact Step.pushEdge hi hs nxtV d e hok.e_lt hok.q hok.ends hok.d_le
  | succ k ih =>
    have st := ih fun j hj => hk j (Nat.le_succ_of_le hj)
    rw [iter_succ']
    exact st.trans (Step.loop1Body st.inv st.shape (by rw [st.g]; exact hv)
      (hok.body k fun j hj => hk j (Nat.le_succ_of_le hj)))

/-! ### Admitted content lemmas (PROOF.md §4.5, walk side) -/

/-- Admitted (Lemma 4.3, R-maximality content, per `finishEdge`): with the invariant settled at
`curV` (`stackVerts[d]`), `finishEdge` keeps it settled at `curV`. Loop 1's S/P/R closes build
`EntryR` entries from `EntryR` entries (the merged items are the new entry's pieces; the complement
of a closed `(nxt.vStart, d)` set is one class because all `(nxt.vStart, d)` classes were P-merged
when `nxt.vStart` was finished); Loop 2 and `closeVert` glue the classes returning to `curV` into
one `(curV, ·)` entry; the P-check bonds a type-1 class with the previous `(curV, lv)` entry; the
back-edge branch pushes a `(curV, lv)` entry — all `(curV, ·)` entries are exempt. The statement
takes `FinishOk`/`FinishGuards` (stack shape) as hypotheses; `dfs` is the sorted DFS tree of the
block `s.g` and `stackVerts[0..d]` its ancestor chain of `curV`. -/
theorem finishEdge_rInvAt {D : Nat} (curV d lv : Nat) (kind : RetKind) (o : DfsOut) (origTstack : Nat)
    (hasVert : Bool) (ho : o.cls = .ret lv kind) (hlow : lv < d) (hv : curV < s.g.nv)
    (hi : s.Inv D) (hs : Shape s) (hok : FinishOk D curV d lv o origTstack hasVert s)
    (hg : FinishGuards d o origTstack hasVert s)
    (h2 : s.g.TwoConnected) (hsp : dfs.Spec s.g) (hrt : dfs.Rooted s.g)
    (hd : dfs.depth curV = d) (hcur : s.stackVerts[d]! = curV)
    (hanc : ∀ k, k ≤ d → dfs.Anc s.stackVerts[k]! curV ∧ dfs.depth s.stackVerts[k]! = k)
    (hR : s.RInvAt dfs curV) :
    (after (finishEdge curV d o origTstack hasVert) s).RInvAt dfs curV := by
  sorry

/-- Admitted (Lemma 4.3, ascend): walking the subtree of the child `c = stackVerts[d+1]` from a
state settled at `c` (no entry has bottom `c` yet) leaves a state settled at its parent
`stackVerts[d]`: once `c`'s out-edges are done every `(c, l)` class has been P-merged into the single
`(c, l)` entry, which is then `maximal`. -/
theorem walkTree_rInvAt {D : Nat} (d : Nat) (c : Nat) (outs : List DfsOut) (hi : s.Inv D) (hs : Shape s)
    (hg : GuardsTree (.node c outs) (d + 1) s)
    (h2 : s.g.TwoConnected) (hsp : dfs.Spec s.g) (hrt : dfs.Rooted s.g)
    (hc : dfs.IsParent s.stackVerts[d]! c) (hR : s.RInvAt dfs c) :
    (after (walkTree (.node c outs) (d + 1)) s).RInvAt dfs s.stackVerts[d]! := by
  sorry

/-- Admitted (Lemma 4.3 at the R branch): during `closeEars` of `finishEdge` at depth `d` for the
tree edge `e` to the child `nxtV = stackVerts[d+1]`, every iterate of Loop 1 at which `loop1Type`
answers `.R` has the shape `RBranch`, and its two top entries are `EntryR` and edge-disjoint
(`RTop`). `cur`'s bottom is the child, finished, so `cur`'s `(nxtV, d)` classes are P-merged and
`cur` is settled; `nxt`'s bottom is a finished vertex strictly below the child. The hypotheses
below do not determine which edges the entries hold, so `interior`, `proper`, `nxt_ne`, `cur_c`
need the Loop-1 ear content (entries with `topDepth > d` have `vStart = nxtV`; the entries with
`topDepth ≥ d` hold exactly the tree edge and the child's subtree edges) as a further hypothesis
or from `EarShape`; only `tstack`, `cur_top`, `nxt_top`, `ne` follow from `run_loop1Cond`,
`loop1Type_run` and the head-`topDepth` induction. -/
theorem loop1_rBranch {D nxtV d e : Nat} {edgeDir : Bool} (hi : s.Inv D) (hs : Shape s)
    (hok : CloseEarsOk D nxtV d e edgeDir s)
    (h2 : s.g.TwoConnected) (hsp : dfs.Spec s.g) (hrt : dfs.Rooted s.g)
    (hc : nxtV = s.stackVerts[d + 1]!) (hR : s.RInvAt dfs s.stackVerts[d]!) (k : Nat)
    (hk : ∀ j, j ≤ k → result (loop1Cond d) (iter (loop1Body d edgeDir) j (ceS₁ nxtV d e s)) = true)
    (hty : l1Ty d edgeDir (iter (loop1Body d edgeDir) k (ceS₁ nxtV d e s)) = .R) :
    ∃ cur nxt rest, (iter (loop1Body d edgeDir) k (ceS₁ nxtV d e s)).RBranch d cur nxt rest ∧
      (iter (loop1Body d edgeDir) k (ceS₁ nxtV d e s)).RTop dfs cur nxt := by
  sorry

/-- The R skeleton closed at any R branch of Loop 1 is 3-connected (modulo `loop1_rBranch`): the
walk at depth `d` runs under `Inv (d+1)` (via `Step`; `RBranch.threeConnected` itself only needs
`Inv' (d+1)`), and `Step` keeps `g` and `stackVerts`. -/
theorem loop1_r_threeConnected {nxtV d e : Nat} {edgeDir : Bool} (hi : s.Inv (d + 1)) (hs : Shape s)
    (hok : CloseEarsOk (d + 1) nxtV d e edgeDir s) (hv : nxtV < s.g.nv)
    (h2 : s.g.TwoConnected) (hsp : dfs.Spec s.g) (hrt : dfs.Rooted s.g)
    (hc : nxtV = s.stackVerts[d + 1]!) (hR : s.RInvAt dfs s.stackVerts[d]!) (k : Nat)
    (hk : ∀ j, j ≤ k → result (loop1Cond d) (iter (loop1Body d edgeDir) j (ceS₁ nxtV d e s)) = true)
    (hty : l1Ty d edgeDir (iter (loop1Body d edgeDir) k (ceS₁ nxtV d e s)) = .R) :
    let sk := iter (loop1Body d edgeDir) k (ceS₁ nxtV d e s)
    ∃ cur nxt rest, sk.RBranch d cur nxt rest ∧
      (((Pieces.ofItems sk.g sk.items (sk.rPieceItems cur nxt)).addParent sk.g (sk.rU cur nxt)
        nxt.vStart sk.stackVerts[d]!).contract sk.g).ThreeConnected := by
  intro sk
  obtain ⟨cur, nxt, rest, hb, hR'⟩ := loop1_rBranch hi hs hok h2 hsp hrt hc hR k hk hty
  have st := closeEars_iter_step (v := nxtV) hi hs hv hok k hk
  have hg : sk.g = s.g := st.g
  exact ⟨cur, nxt, rest, hb,
    hb.threeConnected (.of_inv st.inv) (hg ▸ h2) (hg ▸ hsp) (hg ▸ hrt) hR'⟩

end WalkState

end Spqr
