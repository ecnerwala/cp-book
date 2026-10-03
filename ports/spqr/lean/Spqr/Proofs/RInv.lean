import Spqr.RInv
import Spqr.Proofs.RClose
import Spqr.Proofs.Type2

/-!
# `RStep` and `RContent` from the R-maximality invariant (PROOF.md §4.5, walk side)

At Loop 1's R branch, `RTop` (`EntryR` of the two top entries, their edge-disjointness) together with the
branch shape `RBranch` and `Inv' (d+1)` give `WalkState.RStep` and `WalkState.RContent`:

* `pieces`: the merged items are the two entries' pieces, disjoint across entries;
* `single`: `nxt` is attached at the child too, so `EntryR.single` makes it one class, and a
  `{nxt.vStart, stackVerts[d]}`-class of a `cur` edge avoiding `nxt` would be attached at
  `stackVerts[d]` alone (`nxt.vStart` is not a vertex of `cur`);
* `maximal`: per piece, from the entry owning it;
* `bond`: parallel edges in different entries would join `cur`'s terminals (`RBranch.nxt_no_cu`);
* `type1`/`type2`: the two entry-level laminarities combine, since both entries touch a vertex
  outside the pair (the child, or `stackVerts[d]`), and the class is closed under `SepClass`
  (`type1_class` for the type-1 class).
-/

namespace Spqr

open DfsData

namespace WalkState

variable {s : WalkState} {d : Nat} {cur nxt : TEntry} {rest : List TEntry} {dfs : DfsData}

theorem edges_of_entryPieceItems {t : TEntry} {i e : Nat} (hi : i ∈ s.entryPieceItems t)
    (hb : Items.EdgeBelow s.g s.items i e) : t.edges s.g s.items e :=
  ⟨i, List.mem_of_mem_filter hi, hb⟩

theorem mem_rPieceItems_iff {i : ItemId} :
    i ∈ s.rPieceItems cur nxt ↔ i ∈ s.entryPieceItems cur ∨ i ∈ s.entryPieceItems nxt := by
  simp only [rPieceItems, rItems, entryPieceItems, List.mem_filter, List.mem_append]
  tauto

theorem rPieceItems_perm :
    (s.rPieceItems cur nxt).Perm (s.entryPieceItems cur ++ s.entryPieceItems nxt) := by
  simp only [rPieceItems, rItems, entryPieceItems, ← List.filter_append]
  apply List.Perm.filter
  rw [List.perm_iff_count]; intro a; simp only [List.count_append]; omega

theorem sepClass_of_common {a b v e e' : Nat} (hv : v ≠ a) (hv' : v ≠ b) (he : s.g.IsEnd e v)
    (he' : s.g.IsEnd e' v) : s.g.SepClass a b e e' :=
  Graph.EdgeConn.of_reach he he' (.refl ⟨hv, hv'⟩)

theorem RInv.toRTop (hR : s.RInv dfs) (hs : s.tstack = cur :: nxt :: rest) : s.RTop dfs cur nxt := by
  have hd := hR.disj
  rw [hs, List.pairwise_cons] at hd
  exact ⟨hR.entries cur (by simp [hs]), hR.entries nxt (by simp [hs]), hd.1 nxt (List.mem_cons_self ..)⟩

theorem RInvAt.toRTop {v : Nat} (hR : s.RInvAt dfs v) (hs : s.tstack = cur :: nxt :: rest)
    (hc : cur.vStart ≠ v) (hn : nxt.vStart ≠ v) : s.RTop dfs cur nxt := by
  have hd := hR.disj
  rw [hs, List.pairwise_cons] at hd
  exact ⟨hR.entries cur (by simp [hs]) hc, hR.entries nxt (by simp [hs]) hn,
    hd.1 nxt (List.mem_cons_self ..)⟩

theorem RInvTop.toRTop {v d : Nat} (hR : s.RInvTop dfs v d) (hs : s.tstack = cur :: nxt :: rest)
    (hcd : d ≤ cur.topDepth) (hnd : d ≤ nxt.topDepth) (hc : cur.vStart ≠ v) (hn : nxt.vStart ≠ v) :
    s.RTop dfs cur nxt := by
  have hd := hR.disj
  rw [hs, List.pairwise_cons] at hd
  exact ⟨hR.entries cur (by simp [hs]) hcd hc, hR.entries nxt (by simp [hs]) hnd hn,
    hd.1 nxt (List.mem_cons_self ..)⟩

theorem RInvAt.toTop {v : Nat} (hR : s.RInvAt dfs v) (d : Nat) : s.RInvTop dfs v d :=
  ⟨fun t ht _ hv => hR.entries t ht hv, hR.disj⟩

theorem RInvTop.toFront {v d : Nat} (hR : s.RInvTop dfs v d) (origTstack : Nat) :
    s.RInvFront dfs v d origTstack :=
  ⟨fun t ht hd hv => hR.entries t (List.mem_of_mem_drop ht) hd hv, hR.disj⟩

theorem RBranch.pieceItems (hR : s.RTop dfs cur nxt) (hb : s.RBranch d cur nxt rest) :
    s.PieceItems (s.rPieceItems cur nxt) := by
  have hc := (hR.entry_cur).pieces
  have hn := (hR.entry_nxt).pieces
  have both : ∀ i, i ∈ s.entryPieceItems cur → i ∈ s.entryPieceItems nxt → False := by
    intro i hic hin
    obtain ⟨e, -, he⟩ := hc.ne i hic
    exact hR.disj _ (edges_of_entryPieceItems hic he) (edges_of_entryPieceItems hin he)
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_⟩
  · exact rPieceItems_perm.nodup_iff.2
      (List.nodup_append.2 ⟨hc.nodup, hn.nodup, by intro i hic j hjn hij; subst hij; exact both i hic hjn⟩)
  · intro i hi
    rcases mem_rPieceItems_iff.1 hi with hi | hi
    · exact hc.vs i hi
    · exact hn.vs i hi
  · intro i hi
    rcases mem_rPieceItems_iff.1 hi with hi | hi
    · exact hc.ne i hi
    · exact hn.ne i hi
  · intro i hi
    rcases mem_rPieceItems_iff.1 hi with hi | hi
    · exact hc.conn i hi
    · exact hn.conn i hi
  · intro i hi
    rcases mem_rPieceItems_iff.1 hi with hi | hi
    · exact hc.attached i hi
    · exact hn.attached i hi
  · intro i hi j hj hij e hei hej
    rcases mem_rPieceItems_iff.1 hi with hi | hi <;> rcases mem_rPieceItems_iff.1 hj with hj | hj
    · exact hc.disj i hi j hj hij e hei hej
    · exact hR.disj _ (edges_of_entryPieceItems hi hei) (edges_of_entryPieceItems hj hej)
    · exact hR.disj _ (edges_of_entryPieceItems hj hej) (edges_of_entryPieceItems hi hei)
    · exact hn.disj i hi j hj hij e hei hej

theorem RBranch.rStep (hR : s.RTop dfs cur nxt) (hb : s.RBranch d cur nxt rest) :
    s.RStep d cur nxt rest :=
  ⟨hb.tstack, hb.cur_top, hb.nxt_top, hb.ne, hb.interior, hb.cur_ne, hb.nxt_ne, hb.proper,
    hb.pieceItems hR⟩

theorem RBranch.hmid (hb : s.RBranch d cur nxt rest) :
    ∀ k, d < k → k ≤ d + 1 → s.stackVerts[k]! = cur.vStart := by
  intro k h1 h2
  obtain rfl : k = d + 1 := by omega
  exact hb.cur_c.symm

section Content

variable (hR : s.RTop dfs cur nxt) (hb : s.RBranch d cur nxt rest)
include hR hb

local notation "L" => s.rPieceItems cur nxt
local notation "P" => Pieces.ofItems s.g s.items (s.rPieceItems cur nxt)

theorem RBranch.mem_iff {i e : Nat} :
    (P).Mem i e ↔ e < s.g.ne ∧ ∃ h : i < (L).length, Items.EdgeBelow s.g s.items (L)[i] e :=
  Pieces.ofItems_mem_iff (hb.pieceItems hR).disj (hb.pieceItems hR).nodup

/-- A set inside one piece of an entry is inside a piece of `P`. -/
theorem RBranch.piece_of_entry {i : ItemId} (hi : i ∈ L) {K : Nat → Prop}
    (hKlt : ∀ e, K e → e < s.g.ne) (hK : ∀ e, K e → Items.EdgeBelow s.g s.items i e) :
    ∃ j, j < (P).k ∧ ∀ e, K e → (P).Mem j e := by
  obtain ⟨j, hj, rfl⟩ := List.mem_iff_getElem.1 hi
  exact ⟨j, hj, fun e he => (hb.mem_iff hR).2 ⟨hKlt e he, hj, hK e he⟩⟩

omit hR in
theorem RBranch.termPair_cu :
    (P).TermPair cur.vStart s.stackVerts[d]! := by
  obtain ⟨i, hi⟩ := hb.cur_piece
  have hiL : i ∈ L := mem_rPieceItems_iff.2 (.inl hi)
  obtain ⟨j, hj, rfl⟩ := List.mem_iff_getElem.1 hiL
  rcases hb.cur_vs _ hi with hvs | hvs
  · exact ⟨j, hj, .inl ⟨(Pieces.ofItems_x hj hvs).symm, (Pieces.ofItems_y hj hvs).symm⟩⟩
  · exact ⟨j, hj, .inr ⟨(Pieces.ofItems_y hj hvs).symm, (Pieces.ofItems_x hj hvs).symm⟩⟩

/-- A skeleton pair of `U` is an entry-level skeleton pair of both entries. -/
theorem RBranch.entrySkelPair {t : TEntry} (ht : t = cur ∨ t = nxt) {a b : Nat}
    (hsk : (P).SkelPair s.g (s.rU cur nxt) nxt.vStart s.stackVerts[d]! a b) :
    s.EntrySkelPair t a b := by
  obtain ⟨-, -, -, hsb, hnt, hst⟩ := hsk
  have htL : ∀ i ∈ s.entryPieceItems t, i ∈ L := fun i hi =>
    mem_rPieceItems_iff.2 (ht.elim (fun h => .inl (h ▸ hi)) fun h => .inr (h ▸ hi))
  refine ⟨?_, ?_, ?_⟩
  · intro i hi x y hvs hTi
    obtain ⟨j, hj, rfl⟩ := List.mem_iff_getElem.1 (htL i hi)
    by_contra hxy
    push Not at hxy
    obtain ⟨e, he, hE, hv⟩ := hTi
    refine hsb j ⟨hj, ⟨e, he, (hb.mem_iff hR).2 ⟨he, hj, hE⟩, hv⟩, ?_, ?_⟩
    · rw [Pieces.ofItems_x hj hvs]; exact hxy.1
    · rw [Pieces.ofItems_y hj hvs]; exact hxy.2
  · intro i hi x y hvs hab
    obtain ⟨j, hj, rfl⟩ := List.mem_iff_getElem.1 (htL i hi)
    refine hnt ⟨j, hj, ?_⟩
    rw [Pieces.ofItems_x hj hvs, Pieces.ofItems_y hj hvs]; exact hab
  · rcases ht with rfl | rfl
    · rw [hb.cur_top]
      rintro (⟨rfl, rfl⟩ | ⟨rfl, rfl⟩)
      · exact hnt (hb.termPair_cu)
      · obtain ⟨j, hj, h⟩ := hb.termPair_cu
        exact hnt ⟨j, hj, h.elim (fun h => .inr ⟨h.2, h.1⟩) fun h => .inl ⟨h.2, h.1⟩⟩
    · rw [hb.nxt_top]; exact hst

/-- Two entry-level laminarities combine to a laminarity with `U`, for a class `K` closed under
`SepClass a b`, when both entries touch a common vertex outside `{a, b}`. -/
theorem RBranch.laminar_union {K : Nat → Prop} {a b : Nat} (hKlt : ∀ e, K e → e < s.g.ne)
    (hcl : ∀ e e', K e → s.g.SepClass a b e e' → K e')
    (hv : ∃ v, v ≠ a ∧ v ≠ b ∧ s.g.Touches (cur.edges s.g s.items) v ∧
      s.g.Touches (nxt.edges s.g s.items) v)
    (hc : s.EntryLaminar cur K) (hn : s.EntryLaminar nxt K) :
    (P).LaminarWith (fun e => e < s.g.ne ∧ s.rU cur nxt e) K := by
  obtain ⟨v, hva, hvb, ⟨e₁, he₁, hE₁, hv₁⟩, ⟨e₂, he₂, hE₂, hv₂⟩⟩ := hv
  have h12 : s.g.SepClass a b e₁ e₂ :=
    sepClass_of_common hva hvb (Graph.isEnd_iff.2 ⟨he₁, hv₁⟩) (Graph.isEnd_iff.2 ⟨he₂, hv₂⟩)
  rcases hc with ⟨i, hi, hK⟩ | hc | hc
  · exact .inl (hb.piece_of_entry hR (mem_rPieceItems_iff.2 (.inl hi)) hKlt hK)
  · rcases hn with ⟨i, hi, hK⟩ | hn | hn
    · exact .inl (hb.piece_of_entry hR (mem_rPieceItems_iff.2 (.inr hi)) hKlt hK)
    · exact .inr (.inl fun e he => fun | ⟨_, .inl h⟩ => hc e he h | ⟨_, .inr h⟩ => hn e he h)
    · exact absurd hE₁ (hc e₁ (hcl _ _ (hn e₂ he₂ hE₂) h12.symm))
  · rcases hn with ⟨i, hi, hK⟩ | hn | hn
    · exact .inl (hb.piece_of_entry hR (mem_rPieceItems_iff.2 (.inr hi)) hKlt hK)
    · exact absurd hE₂ (hn e₂ (hcl _ _ (hc e₁ he₁ hE₁) h12))
    · exact .inr (.inr fun e => fun | ⟨he, .inl h⟩ => hc e he h | ⟨he, .inr h⟩ => hn e he h)

end Content

/-- `RContent` at the R branch, from `RTop`, the branch shape, and `Inv' (d+1)`. -/
theorem RBranch.rContent (h : s.Inv' (d + 1)) (h2 : s.g.TwoConnected) (hs : dfs.Spec s.g)
    (hR : s.RTop dfs cur nxt) (hb : s.RBranch d cur nxt rest) : s.RContent dfs d cur nxt := by
  have hr := hb.rStep hR
  have hcur := hR.entry_cur
  have hnxt := hR.entry_nxt
  obtain ⟨e₀, he₀, hU₀⟩ := hb.proper
  obtain ⟨e₁, he₁, hE₁⟩ := hb.cur_ne
  have hcurA := hr.cur_attached h hb.hmid
  have att := hr.union_attached h hb.hmid
  obtain ⟨hcurT, hcurU⟩ := hcurA.touches h2 he₁ hE₁ he₀ (fun h => hU₀ (.inl h))
  have hnxtC := hr.nxt_touches h2 hcurA
  have hcu : cur.vStart ≠ s.stackVerts[d]! := hcurA.ne h2 he₁ hE₁ he₀ (fun h => hU₀ (.inl h))
  obtain ⟨-, hne, -, -, -⟩ :=
    Graph.twoAttached_union_classes h2 att ⟨e₁, he₁, .inl hE₁⟩ ⟨e₀, he₀, hU₀⟩
  have hp := hb.pieceItems hR
  have memE : ∀ {i} (hi : i < (s.rPieceItems cur nxt).length), ∀ e, e < s.g.ne →
      ((Pieces.ofItems s.g s.items (s.rPieceItems cur nxt)).Mem i e ↔
        Items.EdgeBelow s.g s.items (s.rPieceItems cur nxt)[i] e) := by
    intro i hi e he
    rw [hb.mem_iff hR]; exact ⟨fun h => h.2.2, fun h => ⟨he, hi, h⟩⟩
  -- `nxt.vStart` is not a vertex of `cur`.
  have hbot : ¬s.g.Touches (cur.edges s.g s.items) nxt.vStart := by
    rintro ⟨e₄, he₄, hE₄, hv₄⟩
    obtain ⟨e₅, he₅, hE₅, hv₅⟩ := hb.nxt_touch_bot
    rcases hcurA _ e₄ e₅ he₄ he₅ hE₄ (fun h => hR.disj _ h hE₅) hv₄ hv₅ with h | h
    · exact hb.ne h
    · exact hne h
  -- the shared vertex outside a skeleton pair
  have shared : ∀ {a b}, (Pieces.ofItems s.g s.items (s.rPieceItems cur nxt)).SkelPair s.g
      (s.rU cur nxt) nxt.vStart s.stackVerts[d]! a b →
      ∃ v, v ≠ a ∧ v ≠ b ∧ s.g.Touches (cur.edges s.g s.items) v ∧
        s.g.Touches (nxt.edges s.g s.items) v := by
    intro a b hsk
    have hnt := hsk.2.2.2.2.1
    by_cases hca : cur.vStart = a
    · subst hca
      refine ⟨s.stackVerts[d]!, hcu.symm, fun hbu => hnt ?_, hcurU, hb.nxt_touch_top⟩
      rw [← hbu]; exact hb.termPair_cu
    by_cases hcb : cur.vStart = b
    · subst hcb
      refine ⟨s.stackVerts[d]!, fun hau => hnt ?_, hcu.symm, hcurU, hb.nxt_touch_top⟩
      obtain ⟨j, hj, hh⟩ := hb.termPair_cu
      rw [← hau]; exact ⟨j, hj, hh.elim (fun h => .inr ⟨h.2, h.1⟩) fun h => .inl ⟨h.2, h.1⟩⟩
    exact ⟨cur.vStart, hca, hcb, hcurT, hnxtC⟩
  refine ⟨?single, ?maximal, ?type1, ?bond, ?type2⟩
  case single =>
    obtain ⟨e₂, he₂, hE₂, hv₂⟩ := hnxtC
    obtain ⟨e₃, he₃, hE₃, hv₃⟩ := hcurT
    have hnA : ¬s.g.TwoAttached (nxt.edges s.g s.items) nxt.vStart s.stackVerts[nxt.topDepth]! := by
      rw [hb.nxt_top]
      intro hA
      rcases hA _ e₂ e₃ he₂ he₃ hE₂ (fun h => hR.disj _ hE₃ h) hv₂ hv₃ with h | h
      · exact hb.ne h.symm
      · exact hcu h
    have hn : ∀ e e', e < s.g.ne → e' < s.g.ne → nxt.edges s.g s.items e →
        nxt.edges s.g s.items e' → s.g.SepClass nxt.vStart s.stackVerts[d]! e e' := by
      have := hnxt.single hnA
      rwa [hb.nxt_top] at this
    -- every `cur` edge is in the class of `e₂`
    have hc : ∀ e, e < s.g.ne → cur.edges s.g s.items e →
        s.g.SepClass nxt.vStart s.stackVerts[d]! e e₂ := by
      intro e he hE
      by_contra hno
      have hKsub : ∀ f, s.g.SepClass nxt.vStart s.stackVerts[d]! e f → cur.edges s.g s.items f := by
        intro f hf
        rcases (att.sepClass_mem hf).1 (.inl hE) with h | h
        · exact h
        · exfalso
          have hflt : f < s.g.ne := by
            rcases hf with rfl | ⟨x, y, -, hy, -⟩
            · exact he
            · exact hy.lt
          exact hno (hf.trans (hn f e₂ hflt he₂ h hE₂))
      have hKA : s.g.TwoAttached (s.g.SepClass nxt.vStart s.stackVerts[d]! e)
          s.stackVerts[d]! s.stackVerts[d]! := by
        intro v f f' hf hf' hF hF' hv hv'
        rcases Graph.sepClass_twoAttached v f f' hf hf' hF hF' hv hv' with rfl | h
        · exact absurd ⟨f, hf, hKsub f hF, hv⟩ hbot
        · exact .inl h
      exact hKA.ne h2 he (Graph.EdgeConn.refl _) he₀ (fun hf => hU₀ (.inl (hKsub _ hf))) rfl
    intro e e' he he' hE hE'
    have h1 : s.g.SepClass nxt.vStart s.stackVerts[d]! e e₂ := by
      rcases hE with hE | hE
      · exact hc e he hE
      · exact hn e e₂ he he₂ hE hE₂
    have h2' : s.g.SepClass nxt.vStart s.stackVerts[d]! e' e₂ := by
      rcases hE' with hE' | hE'
      · exact hc e' he' hE'
      · exact hn e' e₂ he' he₂ hE' hE₂
    exact h1.trans h2'.symm
  case maximal =>
    intro i hi e e' he he' hne hne'
    have hi' : i < (s.rPieceItems cur nxt).length := hi
    obtain ⟨x, y, hvs⟩ := hp.vs _ (List.getElem_mem hi')
    rw [memE hi' e he] at hne
    rw [memE hi' e' he'] at hne'
    rw [Pieces.ofItems_x hi' hvs, Pieces.ofItems_y hi' hvs]
    rcases mem_rPieceItems_iff.1 (List.getElem_mem hi') with hm | hm
    · exact hcur.maximal _ hm x y hvs e e' he he' hne hne'
    · exact hnxt.maximal _ hm x y hvs e e' he he' hne hne'
  case bond =>
    intro a b f₁ f₂ hf hj₁ hj₂ hU₁ hU₂
    have same : ∀ {t : TEntry}, (t = cur ∨ t = nxt) → s.EntryR dfs t →
        t.edges s.g s.items f₁ → t.edges s.g s.items f₂ →
        ∃ i, (Pieces.ofItems s.g s.items (s.rPieceItems cur nxt)).Mem i f₁ ∧
          (Pieces.ofItems s.g s.items (s.rPieceItems cur nxt)).Mem i f₂ := by
      intro t ht hE hE₁ hE₂
      obtain ⟨i, hi, hb₁, hb₂⟩ := hE.bond a b f₁ f₂ hf hj₁ hj₂ hE₁ hE₂
      have hiL : i ∈ s.rPieceItems cur nxt :=
        mem_rPieceItems_iff.2 (ht.elim (fun h => .inl (h ▸ hi)) fun h => .inr (h ▸ hi))
      obtain ⟨j, hj, rfl⟩ := List.mem_iff_getElem.1 hiL
      exact ⟨j, (hb.mem_iff hR).2 ⟨hj₁.lt, hj, hb₁⟩, (hb.mem_iff hR).2 ⟨hj₂.lt, hj, hb₂⟩⟩
    -- a parallel pair across the two entries would join `cur`'s terminals
    have cross : ∀ {f₁ f₂}, f₁ ≠ f₂ → s.g.Joins f₁ a b → s.g.Joins f₂ a b →
        cur.edges s.g s.items f₁ → nxt.edges s.g s.items f₂ → False := by
      intro f₁ f₂ hf hj₁ hj₂ hE₁ hE₂
      have hab : a ≠ b := by
        rintro rfl
        rcases h2 a f₁ e₀ hj₁.lt he₀ with rfl | ⟨x, y, hx, -, hr⟩
        · exact hU₀ (.inl hE₁)
        · have hok : ∀ {x y}, s.g.Reach (· ≠ a) x y → x ≠ a := fun hr => by
            induction hr with
            | refl h => exact h
            | tail _ _ _ ih => exact ih
          exact hok hr (by rcases hx.eq_or hj₁ with h | h <;> exact h)
      have hend : ∀ v, s.g.IsEnd f₁ v → v = cur.vStart ∨ v = s.stackVerts[d]! := by
        intro v hv
        have hv₂ : s.g.IsEnd f₂ v := by
          rcases hv.eq_or hj₁ with rfl | rfl
          · exact hj₂.isEnd
          · exact hj₂.symm.isEnd
        exact hcurA v f₁ f₂ hj₁.lt hj₂.lt hE₁ (fun h => hR.disj _ h hE₂)
          (Graph.isEnd_iff.1 hv).2 (Graph.isEnd_iff.1 hv₂).2
      rcases hend a hj₁.isEnd with rfl | rfl <;> rcases hend b hj₁.symm.isEnd with rfl | rfl
      · exact hab rfl
      · exact hb.nxt_no_cu f₂ hE₂ hj₂
      · exact hb.nxt_no_cu f₂ hE₂ hj₂.symm
      · exact hab rfl
    rcases hU₁ with hE₁ | hE₁ <;> rcases hU₂ with hE₂ | hE₂
    · exact same (.inl rfl) hcur hE₁ hE₂
    · exact absurd (cross hf hj₁ hj₂ hE₁ hE₂) not_false
    · exact absurd (cross hf.symm hj₂ hj₁ hE₂ hE₁) not_false
    · exact same (.inr rfl) hnxt hE₁ hE₂
  case type1 =>
    intro a b hsk hanc o ho hcls
    have hKlt : ∀ e, dfs.EndIn o.dest e s.g → e < s.g.ne := fun _ ⟨_, hx, _⟩ => hx.lt
    have hcl : ∀ e e', dfs.EndIn o.dest e s.g → s.g.SepClass a b e e' → dfs.EndIn o.dest e' s.g :=
      fun e e' hK hc => (type1_class hs hanc ho hcls hK).1 hc
    exact hb.laminar_union hR hKlt hcl (shared hsk)
      (hcur.type1 a b (hb.entrySkelPair hR (.inl rfl) hsk) hanc o ho hcls)
      (hnxt.type1 a b (hb.entrySkelPair hR (.inr rfl) hsk) hanc o ho hcls)
  case type2 =>
    intro a b hsk ht2 o ho hot hanc
    have hKlt : ∀ e, s.g.SepClass a b o.e e → e < s.g.ne := by
      rintro e (rfl | ⟨x, y, -, hy, -⟩)
      · exact (hs.joins a o ho).lt
      · exact hy.lt
    have hcl : ∀ e e', s.g.SepClass a b o.e e → s.g.SepClass a b e e' → s.g.SepClass a b o.e e' :=
      fun e e' hK hc => hK.trans hc
    exact hb.laminar_union hR hKlt hcl (shared hsk)
      (hcur.type2 a b (hb.entrySkelPair hR (.inl rfl) hsk) ht2 o ho hot hanc)
      (hnxt.type2 a b (hb.entrySkelPair hR (.inr rfl) hsk) ht2 o ho hot hanc)

/-- The R skeleton closed at Loop 1's R branch is 3-connected, from the walk invariants
`RTop`, `RBranch`, and `Inv' (d+1)` (no `RContent`/`RStep` hypothesis). -/
theorem RBranch.threeConnected (h : s.Inv' (d + 1)) (h2 : s.g.TwoConnected) (hs : dfs.Spec s.g)
    (hrt : dfs.Rooted s.g) (hR : s.RTop dfs cur nxt) (hb : s.RBranch d cur nxt rest) :
    (((Pieces.ofItems s.g s.items (s.rPieceItems cur nxt)).addParent s.g (s.rU cur nxt)
      nxt.vStart s.stackVerts[d]!).contract s.g).ThreeConnected :=
  (hb.rStep hR).threeConnected' h hb.hmid h2 hs hrt (hb.rContent h h2 hs hR)

end WalkState

end Spqr
