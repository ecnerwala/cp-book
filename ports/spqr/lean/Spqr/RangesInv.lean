import Spqr.Ranges
import Spqr.WalkSpec
import Spqr.Sim

/-!
# The walk-time range invariant (PROOF.md §4.6)

`WalkState.RangesInv σ n D s`: after processing the first `n` edges of the edge order `σ`
(`edgePostorderForest forest`), the open tstack entries own only processed edges, their piece
edges are strictly ordered along `σ` top-down (laminar/consecutive intervals), each entry's piece
is convex in `σ` with its holes inside the entry (blocks hanging from its interior vertices), and
every allocated node item is convex the same way. The attachment half is `Inv' D`
(`WalkSpec.lean`: connected pieces attached at `Term'`, closed items 2-attached at their `vs`).

It is schedule-agnostic: nothing refers to ears or to the merge schedule. The preservation lemmas
for the primitives (`alloc`, `pushVert`, `pushEdge`, `mergeTop`, `finishTop`) take only local
hypotheses on the entries involved (`mergeTop`: the two merged entries are adjacent in `σ`, `hadj`).
Saturation relative to the frontier is not a field here (see §4.6 for the counterexamples to the
attachment-count formulations); the R layer states it as a hypothesis over the entries below the
`Frontier` split (`Proofs/RInvFrame.lean`). Checked empirically at every `finishEdge` by
`checks/RangesInvCheck.lean`.
-/

namespace Spqr
open WalkM

namespace Items
variable {g : Graph} {items items' : Items}

theorem BelowNoV.below {a i : ItemId} (h : items.BelowNoV a i) : items.Below a i := by
  induction h with
  | refl => exact .refl
  | tail _ hstep ih => exact ih.tail hstep.1

theorem PieceEdge.edgeBelow {i e : Nat} (h : items.PieceEdge g i e) : items.EdgeBelow g i e :=
  h.below

theorem BelowNoV_congr (hch : ∀ p, items'.ch p = items.ch p)
    (hty : ∀ p c, items.IsParent p c → items'.type c = items.type c) {a i : ItemId} :
    items'.BelowNoV a i ↔ items.BelowNoV a i := by
  have : (fun p c => items'.IsParent p c ∧ items'.type c ≠ .V) =
      (fun p c => items.IsParent p c ∧ items.type c ≠ .V) := by
    funext p c; apply propext
    rw [IsParent_congr hch]
    constructor <;> rintro ⟨h1, h2⟩
    · exact ⟨h1, by rw [← hty p c h1]; exact h2⟩
    · exact ⟨h1, by rw [hty p c h1]; exact h2⟩
  unfold BelowNoV; rw [this]

theorem BelowNoV_modify_of_not_below (j : ItemId) (f : Item → Item) (hf : ∀ it, (f it).type = it.type)
    {a i : ItemId} (hj : ¬ items.Below a j) :
    Items.BelowNoV (items.modify j f) a i ↔ items.BelowNoV a i := by
  have hty : ∀ c, Items.type (items.modify j f) c = items.type c := type_modify_type_eq j f hf
  constructor <;> intro h
  · induction h with
    | refl => exact .refl
    | @tail b c _ hstep ih =>
      have hb : b ≠ j := fun hb => hj (hb ▸ ih.below)
      rcases IsParent_modify hstep.1 with hp | ⟨hp, _⟩
      · exact ih.tail ⟨hp, by rw [← hty]; exact hstep.2⟩
      · exact absurd hp hb
  · induction h with
    | refl => exact .refl
    | @tail b c hab hstep ih =>
      have hb : b ≠ j := fun hb => hj (hb ▸ Items.BelowNoV.below hab)
      refine ih.tail ⟨?_, by rw [hty]; exact hstep.2⟩
      show c ∈ Items.ch (items.modify j f) b
      rw [ch_modify_of_ne j f hb]; exact hstep.1

end Items

namespace TEntry
variable {g : Graph} {items items' : Items}

/-- The piece edges of an entry: below a non-V span item without passing through a V item. -/
def piece (g : Graph) (items : Items) (t : TEntry) (e : Nat) : Prop :=
  ∃ i ∈ t.spans.1 ++ t.spans.2, items.type i ≠ .V ∧ items.PieceEdge g i e

theorem piece.edges {t : TEntry} {e : Nat} (h : t.piece g items e) : t.edges g items e :=
  let ⟨i, hi, _, hp⟩ := h; ⟨i, hi, hp.edgeBelow⟩

theorem piece_congr {t : TEntry} (hty : ∀ i ∈ t.spans.1 ++ t.spans.2, items'.type i = items.type i)
    (h : ∀ i ∈ t.spans.1 ++ t.spans.2, ∀ e, items'.PieceEdge g i e ↔ items.PieceEdge g i e) (e : Nat) :
    t.piece g items' e ↔ t.piece g items e := by
  simp only [piece]
  constructor <;> rintro ⟨i, hi, hv, hp⟩
  · exact ⟨i, hi, by rwa [hty i hi] at hv, (h i hi e).1 hp⟩
  · exact ⟨i, hi, by rwa [hty i hi], (h i hi e).2 hp⟩

theorem piece_mergeInto (cur nxt : TEntry) (e : Nat) :
    (mergeInto cur nxt).piece g items e ↔ cur.piece g items e ∨ nxt.piece g items e := by
  simp only [piece, mergeInto, List.mem_append]
  constructor
  · rintro ⟨i, hi, hv, hp⟩
    rcases hi with (hi | hi) | (hi | hi)
    exacts [.inl ⟨i, .inl hi, hv, hp⟩, .inr ⟨i, .inl hi, hv, hp⟩, .inr ⟨i, .inr hi, hv, hp⟩,
      .inl ⟨i, .inr hi, hv, hp⟩]
  · rintro (⟨i, hi | hi, hv, hp⟩ | ⟨i, hi | hi, hv, hp⟩)
    exacts [⟨i, .inl (.inl hi), hv, hp⟩, ⟨i, .inr (.inr hi), hv, hp⟩, ⟨i, .inl (.inr hi), hv, hp⟩,
      ⟨i, .inr (.inl hi), hv, hp⟩]

theorem piece_vertEntry (dir : Bool) (v d fi : Nat) (hv : items.type (vertItem v) = .V) (e : Nat) :
    ¬ TEntry.piece g items ⟨v, d, fi, setSides dir [vertItem v] []⟩ e := by
  rintro ⟨i, hi, hne, _⟩
  rw [mem_setSides, List.mem_singleton] at hi; subst hi; exact hne hv

theorem piece_edgeEntry (dir : Bool) (vStart topDepth fi e : Nat)
    (hq : Items.ch items (edgeItem g e) = []) (ht : items.type (edgeItem g e) ≠ .V) (e' : Nat) :
    TEntry.piece g items ⟨vStart, topDepth, fi, setSides dir [edgeItem g e] []⟩ e' ↔ e' = e := by
  constructor
  · intro h; exact (edges_edgeEntry dir vStart topDepth fi e hq e').1 h.edges
  · intro h; rw [h]; exact ⟨edgeItem g e, (mem_setSides dir _ _).2 (List.mem_singleton.2 rfl), ht, .refl⟩

end TEntry

namespace WalkState
variable {s : WalkState} {σ : List Nat} {n D : Nat}

/-- The range invariant after the first `n` edges of `σ`, at depth `D` (see the module doc). -/
structure RangesInv (σ : List Nat) (n D : Nat) (s : WalkState) : Prop where
  /-- Attachments: `Inv' D`. -/
  inv : s.Inv' D
  /-- Open entries own processed edges only. -/
  processed : ∀ t ∈ s.tstack, ∀ e, e < s.g.ne → t.edges s.g s.items e → σ.idxOf e < n
  /-- Every piece edge lies strictly after all edges owned by lower entries. -/
  ordered : ∀ above t below, s.tstack = above ++ t :: below → ∀ t' ∈ below,
    ∀ e e', e < s.g.ne → e' < s.g.ne → t.piece s.g s.items e → t'.edges s.g s.items e' →
      σ.idxOf e' < σ.idxOf e
  /-- An entry's piece is an interval of `σ` whose holes are edges of the entry. -/
  convex : ∀ t ∈ s.tstack, ∀ a b c, a ≤ b → b ≤ c → c < σ.length →
    t.piece s.g s.items σ[a]! → t.piece s.g s.items σ[c]! → t.edges s.g s.items σ[b]!
  /-- A node item's piece is an interval of `σ` whose holes are below the item (`Ranges.convex`). -/
  closed : ∀ i, i < s.items.size → Items.type s.items i ∉ [NodeType.F, .V] →
    ∀ a b c, a ≤ b → b ≤ c → c < σ.length →
      Items.PieceEdge s.g s.items i σ[a]! → Items.PieceEdge s.g s.items i σ[c]! →
        Items.EdgeBelow s.g s.items i σ[b]!

theorem idxOf_getElem! (hnd : σ.Nodup) {b : Nat} (hb : b < σ.length) : σ.idxOf σ[b]! = b := by
  rw [getElem!_pos σ b hb]; exact hnd.idxOf_getElem b hb

theorem getElem!_lt (hσ : ∀ e ∈ σ, e < s.g.ne) {b : Nat} (hb : b < σ.length) : σ[b]! < s.g.ne := by
  rw [getElem!_pos σ b hb]; exact hσ _ (List.getElem_mem hb)

theorem RangesInv.items_congr (h : s.RangesInv σ n D) {items' : Items}
    (hinv : Inv' D { s with items := items' })
    (hE : ∀ t ∈ s.tstack, ∀ e, t.edges s.g items' e ↔ t.edges s.g s.items e)
    (hP : ∀ t ∈ s.tstack, ∀ e, t.piece s.g items' e ↔ t.piece s.g s.items e)
    (hC : ∀ i, i < items'.size → Items.type items' i ∉ [NodeType.F, .V] →
      ∀ a b c, a ≤ b → b ≤ c → c < σ.length →
        Items.PieceEdge s.g items' i σ[a]! → Items.PieceEdge s.g items' i σ[c]! →
          Items.EdgeBelow s.g items' i σ[b]!) :
    RangesInv σ n D { s with items := items' } :=
  ⟨hinv, fun t ht e he hte => h.processed t ht e he ((hE t ht e).1 hte),
   fun above t below hs t' ht' e e' he he' hp hp' =>
     h.ordered above t below hs t' ht' e e' he he' ((hP t (by rw [hs]; simp) e).1 hp)
       ((hE t' (by rw [hs]; simp [ht']) e').1 hp'),
   fun t ht a b c hab hbc hc hpa hpc =>
     (hE t ht _).2 (h.convex t ht a b c hab hbc hc ((hP t ht _).1 hpa) ((hP t ht _).1 hpc)),
   hC⟩

theorem RangesInv.alloc (ty : NodeType) (h : s.RangesInv σ n D) (hσ : ∀ e ∈ σ, e < s.g.ne)
    (hsize : 1 + s.g.nv + s.g.ne ≤ s.items.size)
    (hch : ∀ p c, Items.IsParent s.items p c → c < s.items.size)
    (hsp : ∀ t ∈ s.tstack, ∀ i ∈ t.spans.1 ++ t.spans.2, i < s.items.size) :
    RangesInv σ n D { s with items := s.items.push ⟨ty, (none, none), []⟩ } := by
  have hB : ∀ a i, Items.Below (s.items.push ⟨ty, (none, none), []⟩) a i ↔ Items.Below s.items a i :=
    fun _ _ => Items.Below_push_nil _ rfl
  have hN : ∀ a i, Items.BelowNoV (s.items.push ⟨ty, (none, none), []⟩) a i ↔ Items.BelowNoV s.items a i :=
    fun _ _ => Items.BelowNoV_congr (Items.ch_push_nil _ rfl)
      (fun p c hpc => Items.type_push_of_ne _ (Nat.ne_of_lt (hch p c hpc)))
  refine h.items_congr (h.inv.alloc ty hsize) (fun t _ e => TEntry.edges_congr (fun i _ e => hB i _) e)
    (fun t ht e => TEntry.piece_congr (fun i hi => Items.type_push_of_ne _ (Nat.ne_of_lt (hsp t ht i hi)))
      (fun i _ e => hN i _) e) ?_
  intro i hi hty a b c hab hbc hc hpa hpc
  have hi' : i < s.items.size + 1 := by simpa using hi
  by_cases hi'' : i = s.items.size
  · subst hi''
    exfalso
    rcases Relation.ReflTransGen.cases_head hpa with heq | ⟨c', hc', _⟩
    · have := getElem!_lt hσ (s := s) (by omega : a < σ.length)
      have : s.items.size = 1 + s.g.nv + σ[a]! := heq
      omega
    · simp [Items.IsParent, Items.ch_push_size ⟨ty, (none, none), []⟩ rfl] at hc'
  · rw [Items.type_push_of_ne _ hi''] at hty
    exact (hB i _).2 (h.closed i (by omega) hty a b c hab hbc hc ((hN i _).1 hpa) ((hN i _).1 hpc))

theorem RangesInv.pushVert (v d : Nat) (h : s.RangesInv σ n D)
    (hc : s.g.ConnEdges (Items.EdgeBelow s.g s.items (vertItem v)))
    (ha : s.g.TwoAttached (Items.EdgeBelow s.g s.items (vertItem v)) v v)
    (hv : Items.type s.items (vertItem v) = .V)
    (hn : ∀ e, e < s.g.ne → Items.EdgeBelow s.g s.items (vertItem v) e → σ.idxOf e < n) :
    RangesInv σ n D ((pushVertTstack v d).run s).2 := by
  have hinv := h.inv.pushVert v d hc ha
  rw [pushVertTstack, run_pushTstack] at hinv ⊢
  have hE : ∀ e, TEntry.edges s.g s.items ⟨v, d, s.nxtEdgeIdx, setSides s.stackDir[d]! [vertItem v] []⟩ e ↔
      Items.EdgeBelow s.g s.items (vertItem v) e := by
    intro e; simp only [TEntry.edges, mem_setSides, List.mem_singleton, exists_eq_left]
  have hP := TEntry.piece_vertEntry (g := s.g) s.stackDir[d]! v d s.nxtEdgeIdx hv
  refine ⟨hinv, ?_, ?_, ?_, h.closed⟩
  · intro t ht e he hte
    rcases List.mem_cons.1 ht with rfl | ht
    · exact hn e he ((hE e).1 hte)
    · exact h.processed t ht e he hte
  · intro above t below hs t' ht' e e' he he' hp hp'
    cases above with
    | nil =>
      simp only [List.nil_append, List.cons.injEq] at hs; obtain ⟨rfl, rfl⟩ := hs
      exact absurd hp (hP e)
    | cons a above =>
      simp only [List.cons_append, List.cons.injEq] at hs; obtain ⟨rfl, hs⟩ := hs
      exact h.ordered above t below hs t' ht' e e' he he' hp hp'
  · intro t ht a b c hab hbc hc hpa hpc
    rcases List.mem_cons.1 ht with rfl | ht
    · exact absurd hpa (hP _)
    · exact h.convex t ht a b c hab hbc hc hpa hpc

theorem RangesInv.pushEdge (vStart topDepth e : Nat) (h : s.RangesInv σ n D) (hnd : σ.Nodup)
    (he : e < s.g.ne) (hq : Items.ch s.items (edgeItem s.g e) = [])
    (htq : Items.type s.items (edgeItem s.g e) ≠ .V)
    (hend : Items.PairEq (vStart, s.stackVerts[topDepth]!) s.g.edges[e]!) (hD : topDepth ≤ D)
    (hn : σ[n]? = some e) :
    RangesInv σ (n + 1) D ((pushEdgeTstack vStart topDepth e).run s).2 := by
  have hinv := h.inv.pushEdge vStart topDepth e he hq hend hD
  rw [run_pushEdgeTstack] at hinv ⊢
  have hE := TEntry.edges_edgeEntry (items := s.items) s.stackDir[topDepth]! vStart topDepth s.nxtEdgeIdx e hq
  have hP := TEntry.piece_edgeEntry (items := s.items) s.stackDir[topDepth]! vStart topDepth s.nxtEdgeIdx e hq htq
  obtain ⟨hlt, hget⟩ := List.getElem?_eq_some_iff.1 hn
  have hidx : σ.idxOf e = n := by rw [← hget]; exact hnd.idxOf_getElem n hlt
  refine ⟨hinv, ?_, ?_, ?_, h.closed⟩
  · intro t ht e' he' hte
    rcases List.mem_cons.1 ht with rfl | ht
    · rw [(hE e').1 hte, hidx]; omega
    · exact Nat.lt_succ_of_lt (h.processed t ht e' he' hte)
  · intro above t below hs t' ht' e₁ e₂ he₁ he₂ hp hp'
    cases above with
    | nil =>
      simp only [List.nil_append, List.cons.injEq] at hs; obtain ⟨rfl, rfl⟩ := hs
      rw [(hP e₁).1 hp, hidx]; exact h.processed t' ht' e₂ he₂ hp'
    | cons a above =>
      simp only [List.cons_append, List.cons.injEq] at hs; obtain ⟨rfl, hs⟩ := hs
      exact h.ordered above t below hs t' ht' e₁ e₂ he₁ he₂ hp hp'
  · intro t ht a b c hab hbc hc hpa hpc
    rcases List.mem_cons.1 ht with rfl | ht
    · have ha : a = n := by
        have := idxOf_getElem! hnd (σ := σ) (by omega : a < σ.length)
        rw [(hP _).1 hpa, hidx] at this; omega
      have hc' : c = n := by
        have := idxOf_getElem! hnd (σ := σ) hc
        rw [(hP _).1 hpc, hidx] at this; omega
      have hb : b = n := by omega
      subst hb
      rw [hE, getElem!_pos σ b hlt]; exact hget
    · exact h.convex t ht a b c hab hbc hc hpa hpc

/-- `mergeTstackTops` under `MergeOk` (as `Inv'.mergeTop`) and the local range condition `hadj`:
the two merged entries are adjacent in `σ` (every edge between a piece edge of `nxt` and one of
`cur` belongs to one of them). -/
theorem RangesInv.mergeTop (cur nxt : TEntry) (rest : List TEntry) (hs : s.tstack = cur :: nxt :: rest)
    (h : s.RangesInv σ n D) (hnd : σ.Nodup) (hσ : ∀ e ∈ σ, e < s.g.ne) (hok : MergeOk D s cur nxt)
    (hdisj : ∀ t ∈ rest, ∀ e, e < s.g.ne → t.edges s.g s.items e →
      ¬ (cur.edges s.g s.items e ∨ nxt.edges s.g s.items e))
    (hadj : ∀ a b c, a ≤ b → b ≤ c → c < σ.length →
      nxt.piece s.g s.items σ[a]! → cur.piece s.g s.items σ[c]! →
        cur.edges s.g s.items σ[b]! ∨ nxt.edges s.g s.items σ[b]!) :
    RangesInv σ n D (mergeTstackTops.run s).2 := by
  have hinv := h.inv.mergeTop cur nxt rest hs hok hdisj
  rw [mergeTstackTops_run_eq s cur nxt rest hs] at hinv ⊢
  have hcur : cur ∈ s.tstack := by simp [hs]
  have hnxt : nxt ∈ s.tstack := by simp [hs]
  have hE := TEntry.edges_mergeInto (g := s.g) (items := s.items) cur nxt
  have hP := TEntry.piece_mergeInto (g := s.g) (items := s.items) cur nxt
  refine ⟨hinv, ?_, ?_, ?_, h.closed⟩
  · intro t ht e he hte
    rcases List.mem_cons.1 ht with rfl | ht
    · rcases (hE e).1 hte with h1 | h1
      exacts [h.processed cur hcur e he h1, h.processed nxt hnxt e he h1]
    · exact h.processed t (by simp [hs, ht]) e he hte
  · intro above t below hs' t' ht' e e' he he' hp hp'
    cases above with
    | nil =>
      simp only [List.nil_append, List.cons.injEq] at hs'; obtain ⟨rfl, rfl⟩ := hs'
      rcases (hP e).1 hp with h1 | h1
      · exact h.ordered [] cur (nxt :: rest) hs t' (List.mem_cons_of_mem _ ht') e e' he he' h1 hp'
      · exact h.ordered [cur] nxt rest hs t' ht' e e' he he' h1 hp'
    | cons a above =>
      simp only [List.cons_append, List.cons.injEq] at hs'; obtain ⟨rfl, hs'⟩ := hs'
      exact h.ordered (cur :: nxt :: above) t below (by rw [hs, hs']; rfl) t' ht' e e' he he' hp hp'
  · intro t ht a b c hab hbc hc hpa hpc
    rcases List.mem_cons.1 ht with rfl | ht
    · rw [hE]
      rcases (hP _).1 hpa with ha | ha <;> rcases (hP _).1 hpc with hc' | hc'
      · exact .inl (h.convex cur hcur a b c hab hbc hc ha hc')
      · have := h.ordered [] cur (nxt :: rest) hs nxt (List.mem_cons_self ..) _ _
          (getElem!_lt hσ (s := s) (by omega : a < σ.length)) (getElem!_lt hσ (s := s) hc) ha hc'.edges
        rw [idxOf_getElem! hnd (by omega : a < σ.length), idxOf_getElem! hnd hc] at this; omega
      · exact hadj a b c hab hbc hc ha hc'
      · exact .inr (h.convex nxt hnxt a b c hab hbc hc ha hc')
    · exact h.convex t (by simp [hs, ht]) a b c hab hbc hc hpa hpc

/-- `finishTstackTop` under the hypotheses of `Inv'.finishTop`: the closed item takes over the
top entry's edges and a subset of its pieces, so its convexity is the entry's. -/
theorem RangesInv.finishTop (item : ItemId) (t : TEntry) (rest : List TEntry)
    (hs : s.tstack = t :: rest) (h : s.RangesInv σ n D) (hσ : ∀ e ∈ σ, e < s.g.ne)
    (hitem : item < s.items.size) (hnode : 1 + s.g.nv + s.g.ne ≤ item)
    (hroot : ∀ p, ¬ Items.IsParent s.items p item)
    (hfree : ∀ t' ∈ s.tstack, item ∉ t'.spans.1 ++ t'.spans.2)
    (hside : getSide t.spans (!s.stackDir[t.topDepth]!) = [])
    (hmid : ∀ k, t.topDepth < k → k ≤ D → s.stackVerts[k]! = t.vStart ∨
      s.g.Interior (t.edges s.g s.items) s.stackVerts[k]! ∨
      ¬ s.g.Touches (t.edges s.g s.items) s.stackVerts[k]!) :
    RangesInv σ n D ((finishTstackTop item).run s).2 := by
  have hinv := h.inv.finishTop item t rest hs hitem hnode hroot hfree hside hmid
  rw [finishTstackTop_run_eq s item t rest hs] at hinv ⊢
  set dir := s.stackDir[t.topDepth]!
  set f : Item → Item := fun it =>
    { it with vs := setSides dir (some s.stackVerts[t.topDepth]!) (some t.vStart),
              ch := getSide t.spans dir }
  have hf : ∀ it, (f it).type = it.type := fun _ => rfl
  have hty : ∀ c, Items.type (s.items.modify item f) c = Items.type s.items c :=
    Items.type_modify_type_eq item f hf
  have hitem_edge : ∀ e, e < s.g.ne → edgeItem s.g e ≠ item := by
    intro e he; show 1 + s.g.nv + e ≠ item; omega
  have hnb : ∀ i, i ≠ item → ¬ Items.Below s.items i item := fun i hi hb =>
    hi (Items.Below.eq_of_no_parent hroot hb)
  have hch : Items.ch (s.items.modify item f) item = getSide t.spans dir := by
    rw [Items.ch_modify_at item f hitem]
  have hne : ∀ c ∈ getSide t.spans dir, c ≠ item := fun c hc hce => by
    subst hce; exact hfree t (by simp [hs]) ((mem_of_getSide_nil dir t.spans hside c).2 hc)
  have ht : t ∈ s.tstack := by simp [hs]
  have hsub : ∀ e, e < s.g.ne →
      (Items.EdgeBelow s.g (s.items.modify item f) item e ↔ t.edges s.g s.items e) := by
    intro e he
    simp only [TEntry.edges, Items.EdgeBelow, mem_of_getSide_nil dir t.spans hside]
    constructor
    · intro hb
      rcases hb.head_cases with heq | ⟨c, hc, hb⟩
      · exact absurd heq.symm (hitem_edge e he)
      · simp only [Items.IsParent, hch] at hc
        exact ⟨c, hc, (Items.Below_modify_of_not_below item f (hnb c (hne c hc))).1 hb⟩
    · rintro ⟨c, hc, hb⟩
      exact .head (by simpa [Items.IsParent, hch] using hc)
        ((Items.Below_modify_of_not_below item f (hnb c (hne c hc))).2 hb)
  have hpiece : ∀ e, e < s.g.ne → Items.PieceEdge s.g (s.items.modify item f) item e →
      t.piece s.g s.items e := by
    intro e he hp
    rcases Relation.ReflTransGen.cases_head hp with heq | ⟨c, ⟨hc, hcv⟩, hb⟩
    · exact absurd heq.symm (hitem_edge e he)
    · simp only [Items.IsParent, hch] at hc
      refine ⟨c, (mem_of_getSide_nil dir t.spans hside c).2 hc, by rwa [hty] at hcv, ?_⟩
      exact (Items.BelowNoV_modify_of_not_below item f hf (hnb c (hne c hc))).1 hb
  have hE' : ∀ e, e < s.g.ne →
      (TEntry.edges s.g (s.items.modify item f) { t with spans := setSides dir [item] [] } e ↔
        t.edges s.g s.items e) := by
    intro e he
    rw [← hsub e he]
    simp only [TEntry.edges]
    constructor
    · rintro ⟨i, hi, hb⟩
      rwa [List.mem_singleton.1 ((mem_setSides dir [item] i).1 hi)] at hb
    · exact fun hb => ⟨item, (mem_setSides dir [item] item).2 (List.mem_singleton.2 rfl), hb⟩
  have hP' : ∀ e, e < s.g.ne →
      TEntry.piece s.g (s.items.modify item f) { t with spans := setSides dir [item] [] } e →
        t.piece s.g s.items e := by
    rintro e he ⟨i, hi, _, hp⟩
    rw [List.mem_singleton.1 ((mem_setSides dir [item] i).1 hi)] at hp
    exact hpiece e he hp
  have hEu : ∀ u ∈ rest, ∀ e,
      TEntry.edges s.g (s.items.modify item f) u e ↔ u.edges s.g s.items e :=
    fun u hu e => TEntry.edges_modify_of_not_mem item f hroot (hfree u (by simp [hs, hu])) e
  have hPu : ∀ u ∈ rest, ∀ e,
      TEntry.piece s.g (s.items.modify item f) u e ↔ u.piece s.g s.items e := by
    intro u hu e
    refine TEntry.piece_congr (fun i _ => hty i) (fun i hi e => ?_) e
    have hi' : i ≠ item := fun h => hfree u (by simp [hs, hu]) (h ▸ hi)
    exact Items.BelowNoV_modify_of_not_below item f hf (hnb i hi')
  refine ⟨hinv, ?_, ?_, ?_, ?_⟩
  · intro u hu e he hue
    rcases List.mem_cons.1 hu with rfl | hu
    · exact h.processed t ht e he ((hE' e he).1 hue)
    · exact h.processed u (by simp [hs, hu]) e he ((hEu u hu e).1 hue)
  · intro above u below hs' u' hu' e e' he he' hp hp'
    cases above with
    | nil =>
      simp only [List.nil_append, List.cons.injEq] at hs'; obtain ⟨rfl, rfl⟩ := hs'
      exact h.ordered [] t rest hs u' hu' e e' he he' (hP' e he hp) ((hEu u' hu' e').1 hp')
    | cons a above =>
      simp only [List.cons_append, List.cons.injEq] at hs'; obtain ⟨rfl, hs'⟩ := hs'
      have hu : u ∈ rest := by rw [hs']; simp
      have hu'' : u' ∈ rest := by rw [hs']; simp [hu']
      exact h.ordered (t :: above) u below (by rw [hs, hs']; rfl) u' hu' e e' he he'
        ((hPu u hu e).1 hp) ((hEu u' hu'' e').1 hp')
  · intro u hu a b c hab hbc hc hpa hpc
    have hσa := getElem!_lt hσ (s := s) (by omega : a < σ.length)
    have hσb := getElem!_lt hσ (s := s) (by omega : b < σ.length)
    have hσc := getElem!_lt hσ (s := s) hc
    rcases List.mem_cons.1 hu with rfl | hu
    · exact (hE' _ hσb).2 (h.convex t ht a b c hab hbc hc (hP' _ hσa hpa) (hP' _ hσc hpc))
    · exact (hEu u hu _).2 (h.convex u (by simp [hs, hu]) a b c hab hbc hc
        ((hPu u hu _).1 hpa) ((hPu u hu _).1 hpc))
  · intro i hi hty' a b c hab hbc hc hpa hpc
    have hσa := getElem!_lt hσ (s := s) (by omega : a < σ.length)
    have hσb := getElem!_lt hσ (s := s) (by omega : b < σ.length)
    have hσc := getElem!_lt hσ (s := s) hc
    simp only [Array.size_modify] at hi
    by_cases hi' : i = item
    · rw [hi'] at hpa hpc ⊢
      exact (hsub _ hσb).2 (h.convex t ht a b c hab hbc hc (hpiece _ hσa hpa) (hpiece _ hσc hpc))
    · rw [hty] at hty'
      exact (Items.Below_modify_of_not_below item f (hnb i hi')).2
        (h.closed i hi hty' a b c hab hbc hc
          ((Items.BelowNoV_modify_of_not_below item f hf (hnb i hi')).1 hpa)
          ((Items.BelowNoV_modify_of_not_below item f hf (hnb i hi')).1 hpc))

end WalkState
end Spqr
