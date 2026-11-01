import Spqr.StVert
import Spqr.StEars

/-!
# St-reading of a finished tree edge

`finishTree_st` composes the per-block preservation lemmas (`closeEars_st`, `L1StInv.mergeLate`,
`closeVert_st`, `finishP_st`, `finishTail_st`) along `finishTree` for a tree edge with
`lowval < d`: the child's entries `sub` (reading as `ps`) above the finished-out-edge entries `pre`
(reading as `qs`) above `B` become, with the vertex entry present, one folded entry reading as
`qs ++ [⟨!stackDir[d], stNest (ps ++ [Q e])⟩]`, and otherwise the ear entries plus the vertex
entry reading as `ps ++ [Q e, V curV]`.
-/

namespace Spqr

open WalkM WalkState

theorem mem_readStack_of_mem {ts : List TEntry} {u : TEntry} {x : ItemId} (hu : u ∈ ts)
    (hx : x ∈ u.spans.1 ++ u.spans.2) : x ∈ readStack ts := by
  induction ts with
  | nil => exact absurd hu List.not_mem_nil
  | cons t ts ih =>
    rcases List.mem_cons.1 hu with rfl | hu
    · exact mem_readStack_cons.2 (Or.inl hx)
    · exact mem_readStack_cons.2 (Or.inr (ih hu))

theorem mem_spans_setSides_single (dir : Bool) (i : ItemId) :
    i ∈ (setSides dir [i] []).1 ++ (setSides dir [i] []).2 := by
  cases dir <;> simp [setSides]

theorem finishTree_st {D d lv : Nat} {kind : RetKind} {o : DfsOut} {s : WalkState} {curV : Nat}
    {hasVert : Bool} {sub pre B : List TEntry} {g : Graph} {ps qs : List StPiece} {blocks : List StBlock}
    (hE : s.EarFinish curV d o hasVert sub (pre ++ B)) (hi : s.Inv' D) (hs : Shape s) (hD : D = d + 1)
    (hok : FinishOk D curV d lv o (pre ++ B).length hasVert s)
    (ho : o.cls = .ret lv kind) (hlow : lv < d) (ht : o.cls.isTree = true)
    (hv : curV < s.g.nv) (he : o.e < s.g.ne) (hq : Items.ch s.items (edgeItem s.g o.e) = [])
    (hends : Items.PairEq (o.dest, s.stackVerts[d]!) s.g.edges[o.e]!)
    (hqty : Items.type s.items (edgeItem s.g o.e) = .Q)
    (hsd : s.stackDir[d]! = !s.stackDir[lv]!)
    (hB : ∀ t ∈ B, t.vStart ≠ curV)
    (hpre : hasVert = false → pre = [] ∧ qs = [])
    (hvf : hasVert = false →
      (∀ p, ¬ Items.IsParent (after (finishP curV lv o.cls.isType1) (feS₂ d o s)).items p (vertItem curV)) ∧
      vertItem curV ∉ readStack (after (finishP curV lv o.cls.isType1) (feS₂ d o s)).tstack)
    (hR : StRead s.items sub ps) (hRq : StRead s.items pre qs) (hI : StItems g s blocks) :
    let r := (finishTree curV d o (pre ++ B).length hasVert s.stackDir[d]!).run (feS₀ d o s)
    r.1 = true ∧ r.2.g = s.g ∧ (∀ k, k ≤ d → r.2.stackDir[k]! = s.stackDir[k]!) ∧
    ∃ new', r.2.tstack = new' ++ B ∧
      StRead r.2.items new' (qs ++
        if hasVert then [⟨!s.stackDir[d]!, stNest (ps ++ [⟨s.stackDir[d]!, [edgeItem s.g o.e]⟩])⟩]
        else ps ++ [⟨s.stackDir[d]!, [edgeItem s.g o.e]⟩, ⟨s.stackDir[d]!, [vertItem curV]⟩]) ∧
      StItems g r.2 blocks := by
  dsimp only
  have hlv : o.cls.lowval d = lv := by rw [ho]; rfl
  have hlow' : o.cls.lowval d < d := by rw [hlv]; exact hlow
  set Qp : StPiece := ⟨s.stackDir[d]!, [edgeItem s.g o.e]⟩ with hQp
  have hJ₁ : L1StInv g s d (ps ++ [Qp]) blocks (pre ++ B) (feS₁ d o s) :=
    closeEars_st hE hi hs hD he hq hends hv ht hlow' hqty (pre := []) (B := pre ++ B) rfl
      (by rw [List.append_nil]; exact hR) hI
  have hJ₁' : L1StInv g s d (qs ++ ps ++ [Qp]) blocks B (feS₁ d o s) :=
    closeEars_st hE hi hs hD he hq hends hv ht hlow' hqty rfl (StRead.append hR hRq) hI
  obtain ⟨c, mid, py, vy, hts₂, ⟨mid₀, hsub⟩, hcl, -, h1⟩ := hE.loops ht hlow'
  obtain ⟨mid₁, py₁, vy₁, hsub₁, hbot⟩ := hE.bottom ht hlow'
  have hpv : [py₁, vy₁] = [py, vy] := List.append_inj_right' (hsub₁.symm.trans hsub) rfl
  rw [(List.cons.inj hpv).1, (List.cons.inj (List.cons.inj hpv).2).1] at hbot
  have hlen₂ : B.length + 1 ≤ (feS₂ d o s).tstack.length := by
    rw [hts₂]; simp only [List.length_cons, List.length_append]; omega
  have hlen₂' : (pre ++ B).length + 1 ≤ (feS₂ d o s).tstack.length := by
    rw [hts₂]; simp only [List.length_cons, List.length_append]; omega
  have hJ₂ : L1StInv g s d (ps ++ [Qp]) blocks (pre ++ B) (feS₂ d o s) :=
    hJ₁.mergeLate d hlen₂'
  have hJ₂' : L1StInv g s d (qs ++ (ps ++ [Qp])) blocks B (feS₂ d o s) := by
    rw [← List.append_assoc]; exact hJ₁'.mergeLate d hlen₂
  have st₀ : Step D curV s (feS₀ d o s) :=
    Step.modifyVs hi hs (edgeItem s.g o.e) _ (by show 1 + s.g.nv + o.e < _; omega)
  have hv₀ : curV < (feS₀ d o s).g.nv := by rw [st₀.g]; exact hv
  have st₁ : Step D curV _ (feS₁ d o s) := Step.closeEars st₀.inv st₀.shape hv₀ (hok.ears ht)
  have hv₁ : curV < (feS₁ d o s).g.nv := by rw [st₁.g]; exact hv₀
  have st₂ : Step D curV _ (feS₂ d o s) := Step.mergeLate st₁.inv st₁.shape hv₁ (hok.late ht)
  have hv₂ : curV < (feS₂ d o s).g.nv := by rw [st₂.g]; exact hv₁
  have hg₂ : (feS₂ d o s).g = s.g := by rw [st₂.g, st₁.g, st₀.g]
  cases hasVert
  · obtain ⟨rfl, rfl⟩ := hpre rfl
    have hr : (finishTree curV d o ([] ++ B).length false s.stackDir[d]!).run (feS₀ d o s) =
        (finishTail curV d false (feSingle d o s)).run
          ((finishP curV lv o.cls.isType1).run (feS₂ d o s)).2 := by
      simp only [finishTree, finishRest, WalkM.run_bind, Bool.false_eq_true, ↓reduceIte, hlv]; rfl
    rw [hr]
    obtain ⟨new₂, hts₂', hR₂⟩ := hJ₂.read
    have hts₂c : (feS₂ d o s).tstack = c :: ((mid ++ [py, vy]) ++ ([] ++ B)) := hts₂
    have hnew₂ : new₂ = c :: (mid ++ [py, vy]) := List.append_cancel_right (hts₂'.symm.trans hts₂c)
    subst hnew₂
    have hP₂ := finishP_st (feS₂ d o s) curV lv o.cls.isType1 c (mid ++ [py, vy]) ([] ++ B) (ps ++ [Qp])
      blocks hts₂c (fun t ht' => hB t (by simpa using ht'))
      (fun h1' t htm hvs _ => by
        obtain ⟨hmid, -⟩ := h1 h1'
        subst hmid
        exact absurd hvs (hE.sub_bot t (by rw [hsub]; exact List.mem_append_right _ (by simpa using htm))))
      hR₂ hJ₂.items
    dsimp only at hP₂
    obtain ⟨new₃, hne₃, hts₃, hsd₃, hg₃, -, hR₃, hI₃⟩ := hP₂
    set s₃ := ((finishP curV lv o.cls.isType1).run (feS₂ d o s)).2 with hs₃
    have st₃ : Step D curV _ s₃ := Step.finishP (v := curV) st₂.inv st₂.shape hv₂ (hok.rest_tree ht rfl).p
    have hv₃ : curV < s₃.g.nv := by rw [st₃.g]; exact hv₂
    have hT := finishTail_st s₃ curV d false (feSingle d o s) new₃ ([] ++ B) (ps ++ [Qp]) blocks hts₃
      (fun _ _ => hne₃)
      (fun _ => ⟨by have := st₃.shape.size; show 1 + curV < _; omega, st₃.shape.vert curV hv₃,
        (hvf rfl).1, (hvf rfl).2⟩)
      hR₃ hI₃
    dsimp only at hT
    obtain ⟨hr1, hsd₄, hg₄, -, new₄, hts₄, hR₄, hI₄⟩ := hT
    have hd₃ : s₃.stackDir[d]! = s.stackDir[d]! := by rw [hsd₃]; exact hJ₂.dirs d (Nat.le_refl d)
    refine ⟨hr1, by rw [hg₄, hg₃, hg₂], fun k hk => by rw [hsd₄, hsd₃]; exact hJ₂.dirs k hk,
      new₄, by simpa using hts₄, ?_, hI₄⟩
    simp only [Bool.false_eq_true, ↓reduceIte, hd₃, List.append_assoc, List.singleton_append,
      List.nil_append] at hR₄ ⊢
    exact hR₄
  · have hr : (finishTree curV d o (pre ++ B).length true s.stackDir[d]!).run (feS₀ d o s) =
        (finishTail curV d true (feB₃ curV d o (pre ++ B).length s)).run
          ((finishP curV lv o.cls.isType1).run
            ((closeVert' curV s.stackDir[d]! o.cls.isType1 (pre ++ B).length (feSingle d o s)).run
              (feS₂ d o s)).2).2 := by
      simp only [finishTree, finishRest, WalkM.run_bind, ↓reduceIte, closeVert_eq, hlv]; rfl
    rw [hr]
    have hcv := closeVert_st (g := g) (s := s) (st := feS₂ d o s) (d := d) (lv := lv) (ps := ps ++ [Qp])
      (blocks := blocks) (base := pre ++ B) curV (pre ++ B).length s.stackDir[d]! o.cls.isType1
      (feSingle d o s) c mid py vy rfl hJ₂ hJ₂' hts₂ (Nat.le_of_lt hlow) rfl hsd
      (by rw [← hlv]; exact hcl) (by have := hbot.vy_top; omega) (fun h1' => (h1 h1').1)
      (fun _ => ⟨by rw [hbot.py_top, hlv], by
        obtain ⟨i, hi, -⟩ := hbot.py_item; exact ⟨i, by rw [hi, hlv]⟩⟩)
    dsimp only at hcv
    set s₃ := ((closeVert' curV s.stackDir[d]! o.cls.isType1 (pre ++ B).length (feSingle d o s)).run
      (feS₂ d o s)).2 with hs₃
    obtain ⟨m, hts₃, -, hsd₃, hg₃, hR₃, hI₃, hT1⟩ := hcv
    have st₃ : Step D curV _ s₃ := Step.closeVert' st₂.inv st₂.shape hv₂ (hok.vert ht rfl)
    have hv₃ : curV < s₃.g.nv := by rw [st₃.g]; exact hv₂
    have hdlv : s₃.stackDir[lv]! = s.stackDir[lv]! := by
      rw [hsd₃]; exact hJ₂.dirs lv (Nat.le_of_lt hlow)
    have hP₃ := finishP_st s₃ curV lv o.cls.isType1 m pre B (qs ++ [⟨!s.stackDir[d]!, stNest (ps ++ [Qp])⟩])
      blocks hts₃ hB
      (fun h1' t htm hvs htd => by
        obtain ⟨hmt, hms⟩ := hT1 h1'
        obtain ⟨⟨i, hspans, -⟩, -⟩ :=
          hE.p_entry hlow' h1' t (List.mem_append_left B htm) hvs (by rw [hlv]; exact htd)
        rw [htd] at hspans
        have hdlv₂ : (feS₂ d o s).stackDir[lv]! = s.stackDir[lv]! := hJ₂.dirs lv (Nat.le_of_lt hlow)
        refine ⟨by rw [hdlv, ← hdlv₂]; exact hms, by rw [hmt], by rw [hdlv, hspans, getSide_setSides_other], ?_⟩
        have hmem : i ∈ readStack s₃.tstack :=
          mem_readStack_of_mem (by rw [hts₃]; exact List.mem_cons_of_mem _ (List.mem_append_left _ htm))
            (by rw [hspans]; exact mem_spans_setSides_single _ i)
        intro dir h hh hty
        by_cases hdir : dir = s.stackDir[lv]!
        · subst hdir
          rw [hspans, getSide_setSides] at hh ⊢
          simp only [List.head!_cons] at hh
          subst hh
          exact ⟨rfl, hI₃.roots i hmem, fun c hc u hu hcu => hI₃.roots c (mem_readStack_of_mem hu hcu) i hc⟩
        · have hdir' : dir = !s.stackDir[lv]! := by
            cases dir <;> cases hsl : s.stackDir[lv]! <;> simp_all
          rw [hspans, hdir', getSide_setSides_other] at hh
          exfalso
          have hh' : h = rootItem := by rw [← hh]; rfl
          rw [hh', st₃.shape.root] at hty
          cases hty)
      hR₃ hI₃
    dsimp only at hP₃
    obtain ⟨new₄, -, hts₄, hsd₄, hg₄, -, hR₄, hI₄⟩ := hP₃
    have hT := finishTail_st _ curV d true (feB₃ curV d o (pre ++ B).length s) new₄ B _ blocks hts₄
      (fun h => by cases h) (fun h => by cases h) hR₄ hI₄
    dsimp only at hT
    obtain ⟨hr1, hsd₅, hg₅, -, new₅, hts₅, hR₅, hI₅⟩ := hT
    refine ⟨hr1, by rw [hg₅, hg₄, hg₃, hg₂], fun k hk => by rw [hsd₅, hsd₄, hsd₃]; exact hJ₂.dirs k hk,
      new₅, hts₅, ?_, hI₅⟩
    simp only [↓reduceIte] at hR₅ ⊢
    exact hR₅

/-- A back edge with `lowval < d`: the entry `⟨curV, lv, [Q e]⟩` is pushed on `pre` (reading as
`qs`), then `finishP`/`finishTail` run as for a tree edge. -/
theorem finishBack_st {D d lv : Nat} {kind : RetKind} {o : DfsOut} {s : WalkState} {curV : Nat}
    {hasVert : Bool} {sub pre B : List TEntry} {g : Graph} {qs : List StPiece} {blocks : List StBlock}
    (hE : s.EarFinish curV d o hasVert sub (pre ++ B)) (hi : s.Inv' D) (hs : Shape s) (hD : D = d)
    (hok : FinishOk D curV d lv o (pre ++ B).length hasVert s)
    (ho : o.cls = .ret lv kind) (hlow : lv < d) (hb : o.cls.isTree = false)
    (hv : curV < s.g.nv) (he : o.e < s.g.ne) (hq : Items.ch s.items (edgeItem s.g o.e) = [])
    (hends : Items.PairEq (curV, s.stackVerts[lv]!) s.g.edges[o.e]!)
    (hqty : Items.type s.items (edgeItem s.g o.e) = .Q)
    (hB : ∀ t ∈ B, t.vStart ≠ curV)
    (hvf : hasVert = false →
      (∀ p, ¬ Items.IsParent (after (finishP curV lv o.cls.isType1) (feBack curV lv d o s)).items p
        (vertItem curV)) ∧
      vertItem curV ∉ readStack (after (finishP curV lv o.cls.isType1) (feBack curV lv d o s)).tstack)
    (hRq : StRead s.items pre qs) (hI : StItems g s blocks) :
    let r := (finishBack curV d o hasVert).run (feS₀ d o s)
    r.1 = true ∧ r.2.g = s.g ∧ r.2.stackDir = s.stackDir ∧
    ∃ new', r.2.tstack = new' ++ B ∧
      StRead r.2.items new' ((qs ++ [⟨s.stackDir[lv]!, [edgeItem s.g o.e]⟩]) ++
        if hasVert then [] else [⟨s.stackDir[d]!, [vertItem curV]⟩]) ∧
      StItems g r.2 blocks := by
  dsimp only
  have hlv : o.cls.lowval d = lv := by rw [ho]; rfl
  have hts : s.tstack = pre ++ B := by rw [hE.tstack, hE.back_nil hb]; rfl
  have hr : (finishBack curV d o hasVert).run (feS₀ d o s) =
      (finishTail curV d hasVert true).run
        ((finishP curV lv o.cls.isType1).run (feBack curV lv d o s)).2 := by
    simp only [finishBack, finishRest, WalkM.run_bind, hlv]; rfl
  rw [hr]
  set f : Item → Item :=
    fun it => { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) } with hf
  set E : TEntry := ⟨curV, lv, s.nxtEdgeIdx, setSides s.stackDir[lv]! [edgeItem s.g o.e] []⟩ with hEdef
  have hqS : ¬ (Items.type s.items (edgeItem s.g o.e) = .S ∨ Items.type s.items (edgeItem s.g o.e) = .P ∨
      Items.type s.items (edgeItem s.g o.e) = .R) := by rw [hqty]; decide
  have hI₀ : StItems g { s with items := s.items.modify (edgeItem s.g o.e) f } blocks :=
    StItems.modify_root _ f (fun _ => rfl) (fun _ => rfl) hE.q_root hqS hI
  have hty₀ : Items.type (s.items.modify (edgeItem s.g o.e) f) (edgeItem s.g o.e) = .Q := by
    rw [Items.type_modify s.items _ _ f (fun _ => rfl)]; exact hqty
  have hch₀ : Items.ch (s.items.modify (edgeItem s.g o.e) f) (edgeItem s.g o.e) = [] := by
    rw [Items.ch_modify_ch_eq _ f (fun _ => rfl)]; exact hq
  have hroot₀ : ∀ p, ¬ Items.IsParent (s.items.modify (edgeItem s.g o.e) f) p (edgeItem s.g o.e) := by
    intro p hp
    unfold Items.IsParent at hp
    rw [Items.ch_modify_ch_eq _ f (fun _ => rfl)] at hp
    exact hE.q_root p hp
  have hfree : edgeItem s.g o.e ∉ readStack s.tstack := fun hm => by
    obtain ⟨u, hu, hmu⟩ := mem_readStack_exists hm
    exact hE.q_free u hu hmu
  have hqlt : edgeItem s.g o.e < (s.items.modify (edgeItem s.g o.e) f).size := by
    rw [Array.size_modify]; have := hs.size; show 1 + s.g.nv + o.e < _; omega
  have hR₀ : StRead (s.items.modify (edgeItem s.g o.e) f) pre qs :=
    hRq.congr fun x _ y _ =>
      ⟨Items.type_modify s.items _ _ f (fun _ => rfl), Items.ch_modify_ch_eq _ f (fun _ => rfl) y⟩
  have hR₁ : StRead (feBack curV lv d o s).items (E :: pre) (qs ++ [⟨s.stackDir[lv]!, [edgeItem s.g o.e]⟩]) :=
    StRead.pushEntry _ _ _ _ _ (Or.inr hty₀) hR₀
  have hI₁ : StItems g (feBack curV lv d o s) blocks :=
    StItems.perm (StItems.pushEntry (s := { s with items := s.items.modify (edgeItem s.g o.e) f }) curV lv
      s.nxtEdgeIdx s.stackDir[lv]! (edgeItem s.g o.e) hqlt hroot₀ hfree hI₀) (List.Perm.refl _) rfl
  have hts₁ : (feBack curV lv d o s).tstack = E :: (pre ++ B) := by show E :: s.tstack = _; rw [hts]
  have st₀ : Step D curV s (feS₀ d o s) :=
    Step.modifyVs hi hs (edgeItem s.g o.e) _ (by show 1 + s.g.nv + o.e < _; omega)
  have st₁ : Step D curV _ (after (pushEdgeTstack curV lv o.e) (feS₀ d o s)) :=
    Step.pushEdge st₀.inv st₀.shape curV lv o.e he hch₀ hends (by rw [hD]; omega)
  have st₂ : Step D curV _ (feBack curV lv d o s) := st₁.trans (Step.frame st₁.inv st₁.shape rfl rfl rfl rfl)
  have hv₂ : curV < (feBack curV lv d o s).g.nv := hv
  have hP := finishP_st (feBack curV lv d o s) curV lv o.cls.isType1 E pre B
    (qs ++ [⟨s.stackDir[lv]!, [edgeItem s.g o.e]⟩]) blocks hts₁ hB
    (fun h1' t htm hvs htd => by
      obtain ⟨⟨i, hspans, -⟩, -⟩ :=
        hE.p_entry (by rw [hlv]; exact hlow) h1' t (List.mem_append_left B htm) hvs (by rw [hlv]; exact htd)
      rw [htd] at hspans
      refine ⟨getSide_setSides_other _ _ _, Nat.le_refl _, by rw [hspans]; exact getSide_setSides_other _ _ _, ?_⟩
      have hmem : i ∈ readStack (feBack curV lv d o s).tstack :=
        mem_readStack_of_mem (by rw [hts₁]; exact List.mem_cons_of_mem _ (List.mem_append_left _ htm))
          (by rw [hspans]; exact mem_spans_setSides_single _ i)
      intro dir h hh hty
      by_cases hdir : dir = s.stackDir[lv]!
      · subst hdir
        rw [hspans, getSide_setSides] at hh ⊢
        simp only [List.head!_cons] at hh
        subst hh
        exact ⟨rfl, hI₁.roots i hmem, fun c hc u hu hcu => hI₁.roots c (mem_readStack_of_mem hu hcu) i hc⟩
      · have hdir' : dir = !s.stackDir[lv]! := by
          cases dir <;> cases hsl : s.stackDir[lv]! <;> simp_all
        rw [hspans, hdir', getSide_setSides_other] at hh
        exfalso
        have hh' : h = rootItem := by rw [← hh]; rfl
        rw [hh', st₂.shape.root] at hty
        cases hty)
    hR₁ hI₁
  dsimp only at hP
  obtain ⟨new₄, -, hts₄, hsd₄, hg₄, -, hR₄, hI₄⟩ := hP
  set s₄ := ((finishP curV lv o.cls.isType1).run (feBack curV lv d o s)).2 with hs₄
  have st₄ : Step D curV _ s₄ := Step.finishP (v := curV) st₂.inv st₂.shape hv₂ (hok.rest_back hb).p
  have hv₄ : curV < s₄.g.nv := by rw [st₄.g]; exact hv₂
  have hT := finishTail_st s₄ curV d hasVert true new₄ B _ blocks hts₄ (fun _ h => by cases h)
    (fun hh => ⟨by have := st₄.shape.size; show 1 + curV < _; omega, st₄.shape.vert curV hv₄,
      (hvf hh).1, (hvf hh).2⟩) hR₄ hI₄
  dsimp only at hT
  obtain ⟨hr1, hsd₅, hg₅, -, new₅, hts₅, hR₅, hI₅⟩ := hT
  have hd₄ : s₄.stackDir[d]! = s.stackDir[d]! := by rw [hsd₄]; rfl
  refine ⟨hr1, by rw [hg₅, hg₄]; rfl, by rw [hsd₅, hsd₄]; rfl, new₅, hts₅, ?_, hI₅⟩
  cases hasVert <;> simp only [Bool.false_eq_true, ↓reduceIte, hd₄, List.append_nil] at hR₅ ⊢ <;> exact hR₅

/-- `walkOutPre`: `setStackDir d` leaves `DirsOf s d` and the readings fixed; the type-1 vertex push
(`!hasVert`, `lowval < d`, type 1) appends the reference's `pre` piece `⟨!stackDir[lowval], [V v]⟩`. -/
theorem walkOutPre_st {g : Graph} {blocks : List StBlock} {base new : List TEntry} {ps : List StPiece}
    (s : WalkState) (v d : Nat) (o : DfsOut) (hasVert : Bool) (hd : d < s.stackDir.size)
    (hts : s.tstack = new ++ base) (hR : StRead s.items new ps) (hI : StItems g s blocks)
    (hvert : hasVert = false →
      vertItem v < s.items.size ∧ Items.type s.items (vertItem v) = .V ∧
      (∀ p, ¬ Items.IsParent s.items p (vertItem v)) ∧ vertItem v ∉ readStack s.tstack) :
    let r := (walkOutPre v d o hasVert).run s
    let x : Bool := if d ≤ o.cls.lowval d then false else !s.stackDir[o.cls.lowval d]!
    let push : Bool := !hasVert && decide (o.cls.lowval d < d) && o.cls.isType1
    r.1 = (push || hasVert) ∧
    r.2.items = s.items ∧ r.2.g = s.g ∧ r.2.stackVerts = s.stackVerts ∧
    r.2.stackDir = s.stackDir.set! d x ∧ DirsOf r.2 d = DirsOf s d ∧
    ∃ new', r.2.tstack = new' ++ base ∧
      StRead r.2.items new' (ps ++ if push then [⟨x, [vertItem v]⟩] else []) ∧
      StItems g r.2 blocks := by
  dsimp only
  set x : Bool := if d ≤ o.cls.lowval d then false else !s.stackDir[o.cls.lowval d]! with hx
  set s₁ : WalkState := { s with stackDir := s.stackDir.set! d x } with hs₁
  have hdirs : DirsOf s₁ d = DirsOf s d := by
    unfold DirsOf
    apply List.map_congr_left
    intro k hk
    rw [List.mem_range] at hk
    show (s.stackDir.set! d x)[k]! = _
    exact Array.getElem!_set!_ne _ _ _ _ (Nat.ne_of_gt hk)
  have hxd : (s.stackDir.set! d x)[d]! = x := Array.getElem!_set!_self _ _ _ hd
  have hR₁ : StRead s₁.items new ps := hR
  have hI₁ : StItems g s₁ blocks := StItems.perm hI (List.Perm.refl _) rfl
  by_cases hp : (!hasVert && decide (o.cls.lowval d < d) && o.cls.isType1) = true
  · have hr : (walkOutPre v d o hasVert).run s =
        (true, { s₁ with tstack :=
          ⟨v, d, s.nxtEdgeIdx, setSides (s.stackDir.set! d x)[d]! [vertItem v] []⟩ :: s.tstack }) := by
      unfold walkOutPre
      simp only [hp, ↓reduceIte]
      rfl
    rw [hr]
    have hhv : hasVert = false := by cases hasVert <;> simp_all
    obtain ⟨hlt, hty, hroot, hfree⟩ := hvert hhv
    refine ⟨by simp [hp], rfl, rfl, rfl, rfl, hdirs,
      ⟨v, d, s.nxtEdgeIdx, setSides x [vertItem v] []⟩ :: new, ?_, ?_, ?_⟩
    · show ⟨v, d, s.nxtEdgeIdx, setSides (s.stackDir.set! d x)[d]! [vertItem v] []⟩ :: s.tstack = _
      rw [hxd, hts]; rfl
    · simp only [hp, ↓reduceIte]
      exact StRead.pushEntry _ _ _ _ _ (Or.inl hty) hR₁
    · have := StItems.pushEntry (s := s₁) v d s.nxtEdgeIdx x (vertItem v) hlt hroot hfree hI₁
      show StItems g { s₁ with tstack :=
        ⟨v, d, s.nxtEdgeIdx, setSides (s.stackDir.set! d x)[d]! [vertItem v] []⟩ :: s.tstack } blocks
      rw [hxd]; exact this
  · have hr : (walkOutPre v d o hasVert).run s = (hasVert, s₁) := by
      unfold walkOutPre
      simp only [hp, ↓reduceIte]
      rfl
    rw [hr]
    simp only [Bool.not_eq_true] at hp
    refine ⟨by simp [hp], rfl, rfl, rfl, rfl, hdirs, new, hts, ?_, hI₁⟩
    simp only [hp, Bool.false_eq_true, ↓reduceIte, List.append_nil]
    exact hR₁

theorem DirsOf_getD (s : WalkState) {k d : Nat} (h : k < d) : (DirsOf s d).getD k false = s.stackDir[k]! := by
  simp [DirsOf, List.getD_eq_getElem?_getD, List.getElem?_range h]

/-- The state after `finishP` on a returning edge: `feS₂` (tree) resp. `feBack` (back edge). -/
def fePState (curV lv d : Nat) (o : DfsOut) (s : WalkState) : WalkState :=
  if o.cls.isTree then after (finishP curV lv o.cls.isType1) (feS₂ d o s)
  else after (finishP curV lv o.cls.isType1) (feBack curV lv d o s)

/-- `finishEdge` for a returning edge (`lowval = lv < d`): the reference's `mid ++ post` pieces for
`o` (with `lowDir = !stackDir[d]`, `sd = stackDir[d]`) are appended to `qs`, above the untouched `B`. -/
theorem finishEdge_st {D d lv : Nat} {kind : RetKind} {o : DfsOut} {s : WalkState} {curV : Nat}
    {hasVert : Bool} {sub pre B : List TEntry} {g : Graph} {ps qs : List StPiece}
    {blocks : List StBlock}
    (hE : s.EarFinish curV d o hasVert sub (pre ++ B)) (hi : s.Inv' D) (hs : Shape s)
    (hD : D = if o.cls.isTree then d + 1 else d)
    (hok : FinishOk D curV d lv o (pre ++ B).length hasVert s)
    (ho : o.cls = .ret lv kind) (hlow : lv < d)
    (hv : curV < s.g.nv) (he : o.e < s.g.ne) (hq : Items.ch s.items (edgeItem s.g o.e) = [])
    (hends : Items.PairEq (if o.cls.isTree then (o.dest, s.stackVerts[d]!) else (curV, s.stackVerts[lv]!))
      s.g.edges[o.e]!)
    (hqty : Items.type s.items (edgeItem s.g o.e) = .Q)
    (hsd : s.stackDir[d]! = !s.stackDir[lv]!)
    (hB : ∀ t ∈ B, t.vStart ≠ curV)
    (hpre : hasVert = false → pre = [] ∧ qs = [])
    (hvf : hasVert = false →
      (∀ p, ¬ Items.IsParent (fePState curV lv d o s).items p (vertItem curV)) ∧
      vertItem curV ∉ readStack (fePState curV lv d o s).tstack)
    (hR : StRead s.items sub ps) (hRq : StRead s.items pre qs) (hI : StItems g s blocks) :
    let r := (finishEdge curV d o (pre ++ B).length hasVert).run s
    r.1 = true ∧ r.2.g = s.g ∧
    (∀ k, k ≤ d → r.2.stackDir[k]! = s.stackDir[k]!) ∧
    ∃ new', r.2.tstack = new' ++ B ∧
      StRead r.2.items new'
        (qs ++
          (if o.cls.isTree then
            (if hasVert then
              [⟨!s.stackDir[d]!, stNest (ps ++ [⟨s.stackDir[d]!, [edgeItem s.g o.e]⟩])⟩]
            else ps ++ [⟨s.stackDir[d]!, [edgeItem s.g o.e]⟩])
          else [⟨!s.stackDir[d]!, [edgeItem s.g o.e]⟩]) ++
          (if hasVert then [] else [⟨s.stackDir[d]!, [vertItem curV]⟩])) ∧
      StItems g r.2 blocks := by
  dsimp only
  have hlv : o.cls.lowval d = lv := by rw [ho]; rfl
  have hge : ¬ (lv ≥ d) := by omega
  have hr : (finishEdge curV d o (pre ++ B).length hasVert).run s =
      if o.cls.isTree = true then
        (finishTree curV d o (pre ++ B).length hasVert s.stackDir[d]!).run (feS₀ d o s)
      else (finishBack curV d o hasVert).run (feS₀ d o s) := by
    rw [finishEdge_eq]
    simp only [finishEdge', hlv, hge, ↓reduceIte, WalkM.run_bind, WalkM.get_run, run_stackDir, run_makeVs,
      run_modifyItem]
    by_cases ht : o.cls.isTree = true <;> simp only [ht, ↓reduceIte, Bool.false_eq_true] <;> rfl
  rw [hr]
  by_cases ht : o.cls.isTree = true
  · rw [if_pos ht]
    rw [if_pos ht] at hD hends
    simp only [fePState, ht, ↓reduceIte] at hvf
    have hT := finishTree_st hE hi hs hD hok ho hlow ht hv he hq hends hqty hsd hB hpre hvf hR hRq hI
    dsimp only at hT
    obtain ⟨h1, hg, hdirs, new', hts, hR', hI'⟩ := hT
    refine ⟨h1, hg, hdirs, new', hts, ?_, hI'⟩
    simp only [ht, ↓reduceIte]
    cases hasVert <;>
      simp only [Bool.false_eq_true, ↓reduceIte, List.append_assoc, List.append_nil, List.singleton_append,
        List.cons_append, List.nil_append] at hR' ⊢ <;> exact hR'
  · have ht' : o.cls.isTree = false := by simpa using ht
    rw [if_neg ht]
    rw [if_neg ht] at hD hends
    simp only [fePState, ht, Bool.false_eq_true, ↓reduceIte] at hvf
    have hBk := finishBack_st hE hi hs hD hok ho hlow ht' hv he hq hends hqty hB hvf hRq hI
    dsimp only at hBk
    obtain ⟨h1, hg, hsd', new', hts, hR', hI'⟩ := hBk
    refine ⟨h1, hg, fun k _ => by rw [hsd'], new', hts, ?_, hI'⟩
    simp only [ht, Bool.false_eq_true, ↓reduceIte]
    simp only [hsd, Bool.not_not, List.append_assoc] at hR' ⊢
    exact hR'

end Spqr
