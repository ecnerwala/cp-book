import Spqr.RangesBoundary

set_option maxHeartbeats 2000000

/-!
# `walkTree_rangesInv`: the range invariant along the walk

Mirrors `walkTree_inv'` (`WalkInv.lean`): the same mutual induction over `walkTree`/`walkOuts`/
`walkOut`, with the ear guards (`GuardsTree`), the bookkeeping (`BookTree`) and, new here, the
range-side hypotheses `RgTree σ n t d s`: at every `finishEdge` the `FinishAdj`/`BoundaryAdj`
bundle (edge position in `σ`, adjacency at each merge), and at every vertex push `PushVertR`. The
edge counter `n` advances by one per `finishEdge`, i.e. by `o.block.length` per out-edge and by
`t.edgePostorder.length` per subtree.

`finishBoundary_rangesInv` handles the block-closing edge (`d ≤ o.cls.lowval d`) using the
ear-side `BoundaryOk` and the range-side `BoundaryAdj`.
-/

namespace Spqr
open WalkM

theorem DfsOut.edgePostorderList_cons (o : DfsOut) (rest : List DfsOut) :
    DfsOut.edgePostorderList (o :: rest) = o.block ++ DfsOut.edgePostorderList rest := by
  cases o <;> simp [DfsOut.edgePostorderList, DfsOut.block]

namespace WalkState
variable {s : WalkState} {σ : List Nat} {n D : Nat}

theorem RangesInv.withInv {D' : Nat} (h : s.RangesInv σ n D) (hi : s.Inv' D') : s.RangesInv σ n D' :=
  ⟨hi, h.processed, h.ordered, h.convex, h.closed⟩

theorem RangesInv.setSv {d : Nat} (x : Nat) (h : s.RangesInv σ n d) :
    ({ s with stackVerts := s.stackVerts.set! (d + 1) x } : WalkState).RangesInv σ n (d + 1) :=
  ⟨h.inv.setSv x, h.processed, h.ordered, h.convex, h.closed⟩

theorem RangesInv.frame' {s' : WalkState} (h : s.RangesInv σ n D) (hg : s'.g = s.g := by rfl)
    (hsv : s'.stackVerts = s.stackVerts := by rfl) (hi : s'.items = s.items := by rfl)
    (hts : s'.tstack = s.tstack := by rfl) : s'.RangesInv σ n D :=
  h.frame hg hsv hi hts

/-- Range-side hypotheses of `finishBoundary` (the block-closing edge `o.e`, `d ≤ o.cls.lowval d`):
the edge sits at position `n` of `σ`, the entries popped into its `Q` item carry the block on the
side that is moved (`side`), and their edges fill `σ` up to `n` from the first piece edge (`fill`),
so the closed `Q` item is convex. -/
structure BoundaryAdj (σ : List Nat) (n d : Nat) (o : DfsOut) (s : WalkState) : Prop where
  pos : σ[n]? = some o.e
  bridge : o.cls.isTree = true → o.cls.lowval d = d + 1 → ∀ t rest, s.tstack = t :: rest →
    t.spans.1 = [] ∧ ∀ a b, a ≤ b → b < n → t.piece s.g s.items σ[a]! → t.edges s.g s.items σ[b]!
  component : o.cls.isTree = true → o.cls.lowval d ≠ d + 1 → ∀ t₁ t₂ rest, s.tstack = t₁ :: t₂ :: rest →
    t₁.spans.2 = [] ∧ t₂.spans.1 = [] ∧ ∀ a b, a ≤ b → b < n →
      (t₁.piece s.g s.items σ[a]! ∨ t₂.piece s.g s.items σ[a]!) →
      (t₁.edges s.g s.items σ[b]! ∨ t₂.edges s.g s.items σ[b]!)

/-- `finishBoundary` preserves the range invariant, advancing `n`, under `BoundaryOk`
(ear layer) and `BoundaryAdj`. -/
theorem finishBoundary_rangesInv {curV d : Nat} {o : DfsOut} {origTstack : Nat} {hasVert : Bool}
    (hge : d ≤ o.cls.lowval d) (h : s.RangesInv σ n D) (hs : Shape s) (hnd : σ.Nodup)
    (hσ : ∀ e ∈ σ, e < s.g.ne) (hb : FinishBook curV d o origTstack hasVert s)
    (hok : BoundaryOk D curV d o s) (hadj : BoundaryAdj σ n d o s) :
    (after (finishEdge curV d o origTstack hasVert) s).RangesInv σ (n + 1) D ∧
      (after (finishEdge curV d o origTstack hasVert) s).g = s.g := by
  have hi := h.inv
  have hge' : o.cls.lowval d ≥ d := hge
  have hqlt : edgeItem s.g o.e < 1 + s.g.nv + s.g.ne := by
    show 1 + s.g.nv + o.e < _; have := hb.e_lt; omega
  have hvlt : vertItem curV < 1 + s.g.nv + s.g.ne := by
    show 1 + curV < _; have := hb.v_lt; omega
  have hvq : vertItem curV ≠ edgeItem s.g o.e := by
    show 1 + curV ≠ 1 + s.g.nv + o.e; have := hb.v_lt; omega
  have hvsz : vertItem curV < s.items.size := by
    show 1 + curV < _; have := hb.v_lt; have := hs.size; omega
  -- the Q item's `vs`, the block counter
  have b₀ : BStep D s { s with items := s.items.modify (edgeItem s.g o.e) fun it => { it with vs := (some curV, none) } } :=
    BStep.modifyVs hi hs _ _ hqlt
  set s₀ := { s with items := s.items.modify (edgeItem s.g o.e) fun it => { it with vs := (some curV, none) } }
    with hs₀
  have b₁ : BStep D s { s₀ with totBlocks := s₀.totBlocks + 1 } := b₀.trans (BStep.frame b₀.inv b₀.shape)
  set s₁ := { s₀ with totBlocks := s₀.totBlocks + 1 } with hs₁
  have r₀ := h.modifyVs (edgeItem s.g o.e) (some curV, none) hqlt
  have r₁ : s₁.RangesInv σ n D := r₀.frame'
  have hn : n < σ.length := (List.getElem?_eq_some_iff.1 hadj.pos).1
  have hq_root₁ : ∀ p, ¬ Items.IsParent s₁.items p (edgeItem s.g o.e) := fun p h =>
    hok.q_root p ((isParent_modifyVs_iff _ _ _ _ _).1 h)
  have hv_root₁ : ∀ p, ¬ Items.IsParent s₁.items p (vertItem curV) := fun p h =>
    hok.v_root p ((isParent_modifyVs_iff _ _ _ _ _).1 h)
  have hts₁ : s₁.tstack = s.tstack := rfl
  have hsz₁ : s₁.items.size = s.items.size := by simp [hs₁, hs₀]
  have fin : ∀ s₂ : WalkState, BStep D s s₂ → s₂.RangesInv σ n D → s₂.g = s.g →
      (∀ p, ¬ Items.IsParent s₂.items p (vertItem curV)) →
      (∀ t ∈ s₂.tstack, vertItem curV ∉ t.spans.1 ++ t.spans.2) →
      let s₃ := { s₂ with items := s₂.items.modify (vertItem curV) fun it => { it with ch := it.ch ++ [edgeItem s.g o.e] } }
      s₃.RangesInv σ (n + 1) D := by
    intro s₂ b r hg hroot hfree
    refine (r.modifyCh (vertItem curV)
      (fun it => { it with ch := it.ch ++ [edgeItem s.g o.e] }) (by rw [hg]; exact hvlt)
      (fun _ => rfl) hroot hfree ?_).advance (Nat.le_succ n)
    intro ht
    exact (ht (by simp [b.shape.vert curV (by rw [hg]; exact hb.v_lt)])).elim
  show wp (finishEdge curV d o origTstack hasVert)
    (fun _ s' => s'.RangesInv σ (n + 1) D ∧ s'.g = s.g) s
  rw [finishEdge_eq]
  simp only [finishEdge', wp_bind, wp_get, wp_stackDir, hge', ↓reduceIte]
  unfold finishBoundary
  simp only [wp_bind, wp_modifyItem, wp_modify, wp_ite, wp_allocItem, wp_makeVs, wp_popTstack, wp_pure, and_true]
  split
  · rename_i hT
    have hpops := hok.pops hT
    split
    · -- bridge
      rename_i hL
      simp only [hL, ↓reduceIte] at hpops
      obtain ⟨t, rest, hts⟩ : ∃ t rest, s.tstack = t :: rest := by
        match h : s.tstack, hpops with
        | t :: rest, _ => exact ⟨t, rest, rfl⟩
      have ht : t ∈ s.tstack := by rw [hts]; exact List.mem_cons_self ..
      have b₂ := b₁.trans (BStep.alloc b₁.inv b₁.shape .I)
      have b₃ := b₂.trans (BStep.modifyVs_leaf b₂.inv b₂.shape s₁.items.size
        (setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest)) (Items.ch_push_size _ rfl))
      have b₄ := b₃.trans (BStep.pop' b₃.inv b₃.shape (s₀ := s) t rest hts rfl rfl
        (fun a i => (Items.Below_modify_ch_eq s₁.items.size
            (fun it => { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) }) fun _ => rfl).trans
          ((Items.Below_push_nil ⟨.I, (none, none), []⟩ rfl).trans
            (Items.Below_modify_ch_eq (edgeItem s.g o.e) (fun it => { it with vs := (some curV, none) }) fun _ => rfl)))
        (hok.gone hT t rest hts))
      have r₂ := r₁.alloc .I hσ b₁.shape.size b₁.shape.ch_lt b₁.shape.span
      have r₃ := r₂.modifyVs_leaf s₁.items.size
        (setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest)) (Items.ch_push_size _ rfl)
      have r₄ := r₃.pop' b₄.inv
      obtain ⟨hside, hfill⟩ := hadj.bridge hT (by simpa using hL) t rest hts
      have cap := h.cap_of_entries hnd hσ [t] t.spans.2
        (by simpa using ht) (by intro i; simp [hside])
        (by simpa using hfill)
      have cap₁ := cap.modifyVs (edgeItem s.g o.e) (some curV, none)
      have cap₂ := cap₁.alloc .I
        (fun i hi => by rw [hsz₁]; exact hs.span t ht i (by simpa [hside] using hi)) b₁.shape.ch_lt
      have cap₃ := cap₂.modifyVs s₁.items.size
        (setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest))
      have cap₄ := cap₃.cons_leaf s₁.items.size (by have := hs.size; rw [hsz₁]; omega)
        (by rw [Items.ch_modify_ch_eq _ (fun it => { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) }) (fun _ => rfl)]; exact Items.ch_push_size _ rfl) hσ (Nat.le_of_lt hn)
      have hsz₄ : ((s₁.items.push ⟨.I, (none, none), []⟩).modify s₁.items.size fun it =>
          { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) }).size = s.items.size + 1 := by
        simp [hs₁, hs₀]
      have hroot₄ : ∀ p, ¬ Items.IsParent ((s₁.items.push ⟨.I, (none, none), []⟩).modify s₁.items.size fun it =>
          { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) }) p (edgeItem s.g o.e) :=
        fun p h => hq_root₁ p ((isParent_push_iff _ _ _ _).1 ((isParent_modifyVs_iff _ _ _ _ _).1 h))
      have hvroot₄ : ∀ p, ¬ Items.IsParent ((s₁.items.push ⟨.I, (none, none), []⟩).modify s₁.items.size fun it =>
          { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) }) p (vertItem curV) :=
        fun p h => hv_root₁ p ((isParent_push_iff _ _ _ _).1 ((isParent_modifyVs_iff _ _ _ _ _).1 h))
      have b₅ := b₄.trans (BStep.modifyCh b₄.inv b₄.shape (edgeItem s.g o.e)
        (fun it => { it with ch := s₁.items.size :: s.tstack.head!.spans.2 }) hqlt (fun _ => rfl) hroot₄
        (fun t' ht' => hok.q_free t' (List.mem_of_mem_tail ht'))
        (fun hj c hc => by
          show c < ((s₁.items.push _).modify _ _).size
          rw [hsz₄]
          rcases List.mem_cons.1 hc with hc | hc
          · rw [hc, hsz₁]; exact Nat.lt_succ_self _
          · exact Nat.lt_succ_of_lt (hs.span t ht c (List.mem_append_right _ (by simpa [hts] using hc)))))
      have r₅ := r₄.modifyCh (edgeItem s.g o.e)
        (fun it => { it with ch := s₁.items.size :: s.tstack.head!.spans.2 }) hqlt (fun _ => rfl) hroot₄
        (fun t' ht' => hok.q_free t' (List.mem_of_mem_tail ht')) (fun _ => by
          apply Items.cap_convex (by simpa [hts] using cap₄) hnd hadj.pos
            (by rw [hsz₄]; exact Nat.lt_succ_of_lt (Nat.lt_of_lt_of_le hqlt hs.size)) hroot₄
          simp only [List.mem_cons, hts, List.head!_cons]
          rintro (heq | hmem)
          · exact (Nat.ne_of_lt (Nat.lt_of_lt_of_le hqlt hs.size)) (heq.trans hsz₁)
          · exact hok.q_free t ht (List.mem_append_right _ hmem))
      refine fin _ b₅ r₅ rfl (fun p h => ?_) (fun t' ht' => hok.v_free t' (List.mem_of_mem_tail ht'))
      · rcases Items.IsParent_modify h with h | ⟨rfl, hj, hc⟩
        · exact hvroot₄ p h
        · rcases List.mem_cons.1 hc with hc | hc
          · exact absurd (hc.trans hsz₁) (Nat.ne_of_lt hvsz)
          · exact hok.v_free t ht (List.mem_append_right _ (by simpa [hts] using hc))
    · -- component
      rename_i hL
      simp only [hL, Bool.false_eq_true, ↓reduceIte] at hpops
      obtain ⟨t₁, t₂, rest, hts⟩ : ∃ t₁ t₂ rest, s.tstack = t₁ :: t₂ :: rest := by
        match h : s.tstack, hpops with
        | t₁ :: t₂ :: rest, _ => exact ⟨t₁, t₂, rest, rfl⟩
      have ht₁ : t₁ ∈ s.tstack := by rw [hts]; exact List.mem_cons_self ..
      have ht₂ : t₂ ∈ s.tstack := by rw [hts]; exact List.mem_cons_of_mem _ (List.mem_cons_self ..)
      have hB₁ : ∀ a i, Items.Below s₁.items a i ↔ Items.Below s.items a i :=
        fun a i => Items.Below_modify_ch_eq (edgeItem s.g o.e) (fun it => { it with vs := (some curV, none) }) fun _ => rfl
      have b₂ := b₁.trans (BStep.pop' b₁.inv b₁.shape (s₀ := s) t₁ (t₂ :: rest) hts rfl rfl hB₁
        (hok.gone hT t₁ (t₂ :: rest) hts))
      have b₃ := b₂.trans (BStep.pop' b₂.inv b₂.shape (s₀ := s) t₂ rest
        (by show s.tstack.tail = _; rw [hts, List.tail_cons]) rfl rfl hB₁ (hok.gone₂ hT (fun h => hL (by simp [h])) t₁ t₂ rest hts))
      have r₂ := r₁.pop' b₂.inv
      have r₃ := r₂.pop' b₃.inv
      obtain ⟨hside₁, hside₂, hfill⟩ := hadj.component hT (by simpa using hL) t₁ t₂ rest hts
      have cap := h.cap_of_entries hnd hσ [t₁, t₂] (t₁.spans.1 ++ t₂.spans.2)
        (by intro t ht; simp at ht; rcases ht with rfl | rfl <;> assumption)
        (by intro i; simp [hside₁, hside₂])
        (by simpa using hfill)
      have cap₁ := cap.modifyVs (edgeItem s.g o.e) (some curV, none)
      have r₄ := r₃.modifyCh (edgeItem s.g o.e)
        (fun it => { it with ch := s.tstack.head!.spans.1 ++ s.tstack.tail.head!.spans.2 }) hqlt
        (fun _ => rfl) hq_root₁
        (fun t' ht' => hok.q_free t' (List.mem_of_mem_tail (List.mem_of_mem_tail ht'))) (fun _ => by
          apply Items.cap_convex (by simpa [hts] using cap₁) hnd hadj.pos
            (by rw [hsz₁]; exact Nat.lt_of_lt_of_le hqlt hs.size) hq_root₁
          simp only [hts, List.head!_cons, List.tail_cons, List.mem_append]
          rintro (hm | hm)
          · exact hok.q_free t₁ ht₁ (List.mem_append_left _ hm)
          · exact hok.q_free t₂ ht₂ (List.mem_append_right _ hm))
      have b₄ := b₃.trans (BStep.modifyCh b₃.inv b₃.shape (edgeItem s.g o.e)
        (fun it => { it with ch := s.tstack.head!.spans.1 ++ s.tstack.tail.head!.spans.2 }) hqlt (fun _ => rfl) hq_root₁
        (fun t' ht' => hok.q_free t' (List.mem_of_mem_tail (List.mem_of_mem_tail ht')))
        (fun hj c hc => by
          show c < s₁.items.size
          rw [hsz₁]
          rcases List.mem_append.1 hc with hc | hc
          · exact hs.span t₁ ht₁ c (List.mem_append_left _ (by simpa [hts] using hc))
          · exact hs.span t₂ ht₂ c (List.mem_append_right _ (by simpa [hts] using hc))))
      refine fin _ b₄ r₄ rfl (fun p h => ?_)
        (fun t' ht' => hok.v_free t' (List.mem_of_mem_tail (List.mem_of_mem_tail ht')))
      · rcases Items.IsParent_modify h with h | ⟨rfl, hj, hc⟩
        · exact hv_root₁ p h
        · rcases List.mem_append.1 hc with hc | hc
          · exact hok.v_free t₁ ht₁ (List.mem_append_left _ (by simpa [hts] using hc))
          · exact hok.v_free t₂ ht₂ (List.mem_append_right _ (by simpa [hts] using hc))
  · -- self-loop
    have b₂ := b₁.trans (BStep.frame (s' := { s₁ with totSelfLoops := s₁.totSelfLoops + 1 }) b₁.inv b₁.shape)
    set s₂ := { s₁ with totSelfLoops := s₁.totSelfLoops + 1 } with hs₂
    have hsz₂ : s₂.items.size = s.items.size := by simp [hs₂, hs₁, hs₀]
    have b₃ := b₂.trans (BStep.alloc b₂.inv b₂.shape .O)
    have b₄ := b₃.trans (BStep.modifyVs_leaf b₃.inv b₃.shape s₂.items.size (some curV, none) (Items.ch_push_size _ rfl))
    have r₂ : s₂.RangesInv σ n D := r₁.frame'
    have r₃ := r₂.alloc .O hσ b₂.shape.size b₂.shape.ch_lt b₂.shape.span
    have r₄ := r₃.modifyVs_leaf s₂.items.size (some curV, none) (Items.ch_push_size _ rfl)
    have cap : Items.CapRange σ n s.g
        ((s₂.items.push ⟨.O, (none, none), []⟩).modify s₂.items.size fun it => { it with vs := (some curV, none) }) [] :=
      ⟨by simp, by simp⟩
    have cap₄ := cap.cons_leaf s₂.items.size (by have := hs.size; rw [hsz₂]; omega)
      (by rw [Items.ch_modify_ch_eq _ (fun it => { it with vs := (some curV, none) }) (fun _ => rfl)]; exact Items.ch_push_size _ rfl) hσ (Nat.le_of_lt hn)
    have hsz₄ : ((s₂.items.push ⟨.O, (none, none), []⟩).modify s₂.items.size fun it =>
        { it with vs := (some curV, none) }).size = s.items.size + 1 := by simp [hs₂, hs₁, hs₀]
    have hroot₄ : ∀ p, ¬ Items.IsParent ((s₂.items.push ⟨.O, (none, none), []⟩).modify s₂.items.size fun it =>
        { it with vs := (some curV, none) }) p (edgeItem s.g o.e) :=
      fun p h => hq_root₁ p ((isParent_push_iff _ _ _ _).1 ((isParent_modifyVs_iff _ _ _ _ _).1 h))
    have hvroot₄ : ∀ p, ¬ Items.IsParent ((s₂.items.push ⟨.O, (none, none), []⟩).modify s₂.items.size fun it =>
        { it with vs := (some curV, none) }) p (vertItem curV) :=
      fun p h => hv_root₁ p ((isParent_push_iff _ _ _ _).1 ((isParent_modifyVs_iff _ _ _ _ _).1 h))
    have b₅ := b₄.trans (BStep.modifyCh b₄.inv b₄.shape (edgeItem s.g o.e)
      (fun it => { it with ch := [s₂.items.size] }) hqlt (fun _ => rfl) hroot₄
      (fun t' ht' => hok.q_free t' ht')
      (fun hj c hc => by
        show c < ((s₂.items.push _).modify _ _).size
        rw [hsz₄, List.mem_singleton.1 hc, hsz₂]; exact Nat.lt_succ_self _))
    have r₅ := r₄.modifyCh (edgeItem s.g o.e) (fun it => { it with ch := [s₂.items.size] })
      hqlt (fun _ => rfl) hroot₄ (fun t' ht' => hok.q_free t' ht') (fun _ =>
        Items.cap_convex cap₄ hnd hadj.pos (by rw [hsz₄]; exact Nat.lt_succ_of_lt (Nat.lt_of_lt_of_le hqlt hs.size)) hroot₄
          (by simp only [List.mem_singleton]; rw [hsz₂]; exact Nat.ne_of_lt (Nat.lt_of_lt_of_le hqlt hs.size)))
    refine fin _ b₅ r₅ rfl (fun p h => ?_) (fun t' ht' => hok.v_free t' ht')
    · rcases Items.IsParent_modify h with h | ⟨rfl, hj, hc⟩
      · exact hvroot₄ p h
      · rw [List.mem_singleton] at hc
        exact absurd (hc.trans hsz₂) (Nat.ne_of_lt hvsz)

/-- The range hypotheses of one `finishEdge` call, at the state where it runs. -/
def FinishR (σ : List Nat) (n curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool) (s : WalkState) :
    Prop :=
  (∀ lv kind, o.cls = .ret lv kind → lv < d → FinishAdj σ n curV d lv o origTstack hasVert s) ∧
  (d ≤ o.cls.lowval d → BoundaryAdj σ n d o s)

/-- `finishEdge` preserves the range invariant (at depth `d`, like `finishEdge_step`) and advances
`n`, under the guards, the bookkeeping and `FinishR`. -/
theorem finishEdge_ranges {v d : Nat} {o : DfsOut} {m : Nat} {hasVert : Bool} {D : Nat}
    (hD : D = if o.cls.isTree then d + 1 else d) (h : s.RangesInv σ n D) (hs : Shape s) (hnd : σ.Nodup)
    (hσ : ∀ e ∈ σ, e < s.g.ne) (hg : FinishGuards d o m hasVert s) (hb : FinishBook v d o m hasVert s)
    (hr : FinishR σ n v d o m hasVert s) :
    (after (finishEdge v d o m hasVert) s).RangesInv σ (n + 1) d ∧ Shape (after (finishEdge v d o m hasVert) s) ∧
      (after (finishEdge v d o m hasVert) s).g = s.g := by
  obtain ⟨hi', hs'⟩ := finishEdge_step hD h.inv hs hg hb
  by_cases hge : d ≤ o.cls.lowval d
  · obtain ⟨hr', hgg⟩ := finishBoundary_rangesInv hge h hs hnd hσ hb (ear_boundary hge hg h.inv hs hb hD) (hr.2 hge)
    exact ⟨hr'.withInv hi', hs', hgg⟩
  · obtain ⟨lv, kind, ho, hl⟩ := ret_of_lowval_lt (Nat.lt_of_not_le hge)
    obtain ⟨sub, base, hlen, hE⟩ := hb.ear
    have hok := finishOk_of_guards ho hl hg hE hlen h.inv hs hD hb.v_lt hb.e_lt hb.q (hb.ends lv kind ho) hb.vert
    have hr' := finishEdge_rangesInv v d lv kind o m hasVert ho hl hb.v_lt h hs hnd hσ hok (hr.1 lv kind ho hl)
    have st := finishEdge_inv v d lv kind o m hasVert ho hl hb.v_lt h.inv hs hok
    exact ⟨hr'.withInv hi', hs', st.g⟩

mutual
/-- The range-side hypotheses of the walk, in the style of `GuardsTree`. -/
def RgTree (σ : List Nat) (n : Nat) (t : DfsTree) (d : Nat) (s : WalkState) : Prop :=
  match t with
  | .node v outs => RgOuts σ n v d outs false { s with stackVerts := s.stackVerts.set! d v }

def RgOuts (σ : List Nat) (n v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  match outs with
  | [] => hasVert = false → PushVertR σ n v s
  | o :: rest => RgOut σ n v d o hasVert s ∧
      wp (walkOut v d o hasVert) (fun hasVert' s' => RgOuts σ (n + o.block.length) v d rest hasVert' s') s

def RgOut (σ : List Nat) (n v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  (hasVert = false → PushVertR σ n v s) ∧
  wp (walkOutPre v d o hasVert) (fun hasVert' s₁ =>
    match o with
    | .tree _ _ child =>
      wp (modify fun s => { s with firstOccurrence := s.firstOccurrence.set! d s.g.ne }) (fun _ s₂ =>
        RgTree σ n child (d + 1) s₂ ∧
        wp (walkTree child (d + 1)) (fun _ s₃ =>
          FinishR σ (n + child.edgePostorder.length) v d o s₁.tstack.length hasVert' s₃) s₂) s₁
    | .back .. => FinishR σ n v d o s₁.tstack.length hasVert' s₁) s
end

theorem lt_ne_of_g {s s' : WalkState} (hg : s'.g = s.g) (hσ : ∀ e ∈ σ, e < s.g.ne) :
    ∀ e ∈ σ, e < s'.g.ne := by
  rw [hg]; exact hσ

theorem walkOutPre_ranges {v d : Nat} {o : DfsOut} {hasVert : Bool} (hi : s.RangesInv σ n d) (hs : Shape s)
    (hσ : ∀ e ∈ σ, e < s.g.ne) (hb : VertBook v hasVert s) (hr : hasVert = false → PushVertR σ n v s) :
    wp (walkOutPre v d o hasVert) (fun _ s₁ => s₁.RangesInv σ n d ∧ Shape s₁ ∧ ∀ e ∈ σ, e < s₁.g.ne) s := by
  have hi' : ({ s with stackDir := s.stackDir.set! d (if o.cls.lowval d ≥ d then false else !s.stackDir[o.cls.lowval d]!) } : WalkState).RangesInv σ n d :=
    hi.frame'
  have hs' : Shape ({ s with stackDir := s.stackDir.set! d (if o.cls.lowval d ≥ d then false else !s.stackDir[o.cls.lowval d]!) } : WalkState) :=
    hs.frame'
  unfold walkOutPre
  simp only [wp_bind, wp_stackDir, wp_setStackDir, wp_ite]
  split
  · rename_i hc
    have hf : hasVert = false := by cases hasVert <;> simp_all
    obtain ⟨hv, hc, ha⟩ := hb hf
    have st := RgStep.pushVert (v := v) hi' hs' d hv hc ha ⟨(hr hf).vtype, (hr hf).below⟩
    exact ⟨st.ranges, st.step.shape, lt_ne_of_g st.step.g hσ⟩
  · exact ⟨hi', hs', hσ⟩

abbrev RgS (σ : List Nat) (n d : Nat) (s : WalkState) : Prop :=
  s.RangesInv σ n d ∧ Shape s ∧ ∀ e ∈ σ, e < s.g.ne

abbrev RgInvTree (σ : List Nat) (n : Nat) (t : DfsTree) (d : Nat) (s : WalkState) : Prop :=
  (∀ v outs, t = .node v outs → ({ s with stackVerts := s.stackVerts.set! d v } : WalkState).RangesInv σ n d) →
  Shape s → σ.Nodup → (∀ e ∈ σ, e < s.g.ne) → GuardsTree t d s → BookTree t d s → RgTree σ n t d s →
  wp (walkTree t d) (fun _ s' => RgS σ (n + t.edgePostorder.length) d s') s

abbrev RgInvOuts (σ : List Nat) (n v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  RgS σ n d s → σ.Nodup → GuardsOuts v d outs hasVert s → BookOuts v d outs hasVert s →
  RgOuts σ n v d outs hasVert s →
  wp (walkOuts v d outs hasVert) (fun hasVert' s' =>
    RgS σ (n + (DfsOut.edgePostorderList outs).length) d s' ∧ VertBook v hasVert' s' ∧
      (hasVert' = false → PushVertR σ (n + (DfsOut.edgePostorderList outs).length) v s')) s

abbrev RgInvOut (σ : List Nat) (n v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  RgS σ n d s → σ.Nodup → GuardsOut v d o hasVert s → BookOut v d o hasVert s → RgOut σ n v d o hasVert s →
  wp (walkOut v d o hasVert) (fun _ s' => RgS σ (n + o.block.length) d s') s

mutual
theorem rgTree : ∀ (σ : List Nat) (n : Nat) (t : DfsTree) (d : Nat) (s : WalkState), RgInvTree σ n t d s
  | σ, n, .node v outs, d, s => fun hi hs hnd hσ hg hb hr => by
    unfold walkTree
    simp only [wp_bind, wp_modify]
    unfold GuardsTree at hg; unfold BookTree at hb; unfold RgTree at hr
    refine wp_imp (wp_of_forall fun hv s' ⟨⟨hi', hs', hσ'⟩, hvb, hpr⟩ => ?_)
      (rgOuts σ n v d outs false _ ⟨hi v outs rfl, hs.frame', hσ⟩ hnd hg hb hr)
    cases hv
    · simp only [Bool.false_eq_true, ↓reduceIte, wp_bind, wp_setStackDir]
      obtain ⟨hv, hc, ha⟩ := hvb rfl
      have hi₂ : ({ s' with stackDir := s'.stackDir.set! d true } : WalkState).RangesInv σ _ d := hi'.frame'
      have hs₂ : Shape ({ s' with stackDir := s'.stackDir.set! d true } : WalkState) := hs'.frame'
      have st := RgStep.pushVert (v := v) hi₂ hs₂ d hv hc ha ⟨(hpr rfl).vtype, (hpr rfl).below⟩
      exact ⟨st.ranges, st.step.shape, lt_ne_of_g st.step.g hσ'⟩
    · exact ⟨hi', hs', hσ'⟩

theorem rgOuts : ∀ (σ : List Nat) (n v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState),
    RgInvOuts σ n v d outs hasVert s
  | σ, n, v, d, [], hasVert, s => fun ⟨hi, hs, hσ⟩ _ _ hb hr => by
    unfold BookOuts at hb; unfold RgOuts at hr
    unfold walkOuts
    simp only [wp_pure, DfsOut.edgePostorderList, List.length_nil, Nat.add_zero]
    exact ⟨⟨hi, hs, hσ⟩, hb, hr⟩
  | σ, n, v, d, o :: rest, hasVert, s => fun hrs hnd hg hb hr => by
    unfold GuardsOuts at hg; unfold BookOuts at hb; unfold RgOuts at hr
    unfold walkOuts
    simp only [wp_bind]
    rw [DfsOut.edgePostorderList_cons, List.length_append, ← Nat.add_assoc]
    exact wp_imp (wp_imp (wp_imp (wp_imp (wp_of_forall fun hv' s' hrs' hg' hb' hr' =>
      rgOuts σ (n + o.block.length) v d rest hv' s' hrs' hnd hg' hb' hr') (rgOut σ n v d o hasVert s hrs hnd hg.1 hb.1 hr.1)) hg.2) hb.2) hr.2

theorem rgOut : ∀ (σ : List Nat) (n v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState), RgInvOut σ n v d o hasVert s
  | σ, n, v, d, o, hasVert, s => fun ⟨hi, hs, hσ⟩ hnd hg hb hr => by
    rw [walkOut_eq, wp_bind]
    unfold GuardsOut at hg; unfold BookOut at hb; unfold RgOut at hr
    refine wp_imp (wp_imp (wp_imp (wp_imp (wp_of_forall fun hv' s₁ ⟨hi₁, hs₁, hσ₁⟩ hg₁ hb₁ hr₁ => ?_)
      (walkOutPre_ranges hi hs hσ hb.1 hr.1)) hg) hb.2) hr.2
    unfold walkOutRest
    rw [wp_bind, wp_tstackSize]
    cases o with
    | back e cls dest =>
      try simp only [wp_pure] at hg₁ hb₁ hr₁
      simp only [wp_bind, wp_pure, DfsOut.block, List.length_singleton]
      obtain ⟨h₁, h₂, h₃⟩ := finishEdge_ranges (D := d)
        (by rw [if_neg fun h => by obtain ⟨_, _, _, h⟩ := hb₁.tree.1 h; cases h]) hi₁ hs₁ hnd hσ₁ hg₁ hb₁ hr₁
      exact ⟨h₁, h₂, lt_ne_of_g h₃ hσ₁⟩
    | tree e cls child =>
      try simp only [wp_bind, wp_modify] at hg₁ hb₁ hr₁
      simp only [wp_bind, wp_modify, DfsOut.block, List.length_append, List.length_singleton, ← Nat.add_assoc]
      refine wp_imp (wp_imp (wp_of_forall fun _ s₃ ⟨hg₃, hb₃, hr₃⟩ ⟨hi₃, hs₃, hσ₃⟩ => ?body)
        (wp_and hg₁.2 (wp_and hb₁.2 hr₁.2)))
        (rgTree σ n child (d + 1) _ ?pre hs₁.frame' hnd hσ₁ hg₁.1 hb₁.1 hr₁.1)
      case body =>
        obtain ⟨h₁, h₂, h₃⟩ := finishEdge_ranges (D := d + 1)
          (by rw [if_pos (hb₃.tree.2 ⟨_, _, _, rfl⟩)]) hi₃ hs₃ hnd hσ₃ hg₃ hb₃ hr₃
        exact ⟨h₁, h₂, lt_ne_of_g h₃ hσ₃⟩
      case pre =>
        exact fun w outs _ =>
          (hi₁.frame' (s' := { s₁ with firstOccurrence := s₁.firstOccurrence.set! d s₁.g.ne })).setSv w
end

/-- `walkTree` preserves the range invariant and advances `n` by the number of edges below `t`,
given the ear guards, the bookkeeping facts and the range-side hypotheses `RgTree`. -/
theorem walkTree_rangesInv (t : DfsTree) (d : Nat) (s : WalkState)
    (hi : ∀ v outs, t = .node v outs → ({ s with stackVerts := s.stackVerts.set! d v } : WalkState).RangesInv σ n d)
    (hs : Shape s) (hnd : σ.Nodup) (hσ : ∀ e ∈ σ, e < s.g.ne) (hg : GuardsTree t d s) (hb : BookTree t d s)
    (hr : RgTree σ n t d s) :
    ((walkTree t d).run s).2.RangesInv σ (n + t.edgePostorder.length) d ∧ Shape ((walkTree t d).run s).2 :=
  let r := rgTree σ n t d s hi hs hnd hσ hg hb hr
  ⟨r.1, r.2.1⟩

end WalkState
end Spqr
