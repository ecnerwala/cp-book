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

end Spqr
