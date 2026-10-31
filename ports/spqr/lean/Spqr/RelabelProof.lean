import Spqr.RelabelLoop

/-!
# `relabel`: per-call specification

`relabel_spec`: from a `CallPre` state, a call `relabel fuel cur p pn ct` (with enough fuel) reaches a
`CallPost` state. The `do` block is never unfolded as a whole: each primitive is abstracted by its
`Step*` relation (`wp_abs` / `wp_jpR`), and the children loop by `LoopInv` (`wp_forIn_inv`).
-/

namespace Spqr.Ghost

open RelabelM

/-- Discharge one arm of a join point whose relation is a conjunction of field equations. -/
macro "jp_armR" : tactic =>
  `(tactic| first
      | exact ⟨_, by constructor <;> first | rfl | simp only [*, ↓reduceIte, reduceCtorEq, Items.nvList] | (split_ifs <;> first | rfl | contradiction), arm_modify _ _ _⟩
      | exact ⟨_, by constructor <;> first | rfl | simp only [*, ↓reduceIte, reduceCtorEq, Items.nvList] | (split_ifs <;> first | rfl | contradiction), arm_pure _ _⟩
      | exact ⟨_, by constructor <;> first | rfl | simp only [*, ↓reduceIte, reduceCtorEq, Items.nvList] | (split_ifs <;> first | rfl | contradiction), arm_id _ _⟩)

set_option pp.deepTerms false
set_option pp.deepTerms.threshold 80
set_option pp.proofs.withType false

theorem relabel_spec {g : Graph} {items : Items} (hwf : items.WF g) :
    ∀ (fuel : Nat) (cur : ItemId) (p pn ct : Option Nat) (s : RelabelState),
      CallPre g items cur s → (items.desc cur).card ≤ fuel →
      wp (relabel fuel cur p pn ct) (fun _ s' => CallPost g items cur p pn ct s s') s
  | 0, cur, p, pn, ct, s, hpre, hf => by
    have := Items.one_le_card_desc hpre.lt
    omega
  | fuel + 1, cur, p, pn, ct, s, hpre, hf => by
    have ih := relabel_spec hwf fuel
    have hc := hpre.cons
    have hcur := hpre.lt
    unfold relabel
    simp only [wp_bind, wp_get, wp_item, wp_modify, hc.items_eq, hc.g_eq, Items.getElem!_ch hcur,
      Items.getElem!_type hcur, Items.getElem!_vs hcur]
    try dsimp only
    -- number `cur`
    refine wp_abs _ (StepNum cur (items.type cur) p pn s) _
      (by constructor <;> first | rfl | exact hc.g_eq.symm | exact hc.items_eq.symm) fun s1 h1 => ?_
    -- V/Q bookkeeping
    refine wp_jpR _ _ (StepVQ g items cur s.types.size) _ (by split <;> jp_armR) fun s2 h2 => ?_
    try simp only [wp_bind, wp_get, wp_modify]
    -- node-verts
    refine wp_abs _ (StepNV ((items.nvList g cur).map (⟨s.types.size, ·⟩)) s2) _ (by constructor <;> rfl)
      fun s3 h3 => ?_
    try dsimp only
    -- R: vertex positions
    refine wp_jpR _ _ (StepPos (items.type cur) ((items.nvList g cur).map (⟨s.types.size, ·⟩)) s2.nodeVerts.size) _
      (by split <;> jp_armR) fun s4 h4 => ?_
    try simp only [wp_bind, wp_get, wp_modify]
    refine wp_orderedChildren' _ _ fun children hchildren => ?_
    try simp only [wp_bind, wp_get, wp_modify]
    -- child slots, skeleton
    refine wp_abs _ (StepLay children _ _ _ _ _ s4) _ (by constructor <;> rfl) fun s6 h6 => ?_
    try simp only [wp_bind, wp_get, wp_modify]
    -- cap twin, then the children loop
    split
    · rename_i hcap
      try simp only [wp_bind, wp_modify]
      refine wp_abs _ (StepCap true s4.nodeEdges.size ct s6) _ (by constructor <;> rfl) fun s7 h7 => ?_
      have he := entry_of_steps hwf hpre ct h1 h2 h3 h4 hchildren h6 h7 ⟨fun _ => hcap, fun _ => rfl⟩
      try simp only [wp_bind]
      refine wp_forIn_inv _ _ _ _
        (fun rest b σ => ∃ done, LoopInv g items cur ct s.types.size s2.chDat.size s2.nodeVerts.size
          s4.nodeEdges.size s s7 children done rest b σ) s7 ⟨[], loop_init hpre he rfl (by rw [he.hasCap_eq, if_pos hcap])⟩ ?_ ?_
      · intro c hc rest b σ ⟨done, hi⟩
        have hcc := hi.mem_ch he
        simp only [wp_bind, wp_get, wp_modify, wp_pure, wp_ite]
        split_ifs with hv hn
        · refine wp_abs _ (StepSlotT (s2.chDat.size + b.2.2) σ.types.size none σ) _ (by constructor <;> rfl)
            fun σ1 hs1 => ?_
          refine wp_mono (ih c _ _ _ σ1 (hi.pre hwf hpre he (hs1.cons hi.cons) hs1.order)
            (loop_fuel hwf.tree hcur hcc hf)) fun _ σ3 hpost => ⟨_, rfl, done ++ [c], ?_⟩
          exact hi.step hwf hpre he (Or.inl ⟨hv, rfl, rfl, rfl, rfl⟩) hs1 hpost
        · refine wp_abs _ (StepSlotT (s2.chDat.size + b.2.2) σ.types.size (some b.2.1) σ) _ (by constructor <;> rfl)
            fun σ1 hs1 => ?_
          refine wp_mono (ih c _ _ _ σ1 (hi.pre hwf hpre he (hs1.cons hi.cons) hs1.order)
            (loop_fuel hwf.tree hcur hcc hf)) fun _ σ3 hpost => ⟨_, rfl, done ++ [c], ?_⟩
          exact hi.step hwf hpre he (Or.inr (Or.inl ⟨hv, hn, rfl, rfl, rfl, rfl⟩)) hs1 hpost
        · refine wp_abs _ (StepSlotT (s2.chDat.size + b.2.2) σ.types.size none σ) _ (by constructor <;> rfl)
            fun σ1 hs1 => ?_
          refine wp_mono (ih c _ _ _ σ1 (hi.pre hwf hpre he (hs1.cons hi.cons) hs1.order)
            (loop_fuel hwf.tree hcur hcc hf)) fun _ σ3 hpost => ⟨_, rfl, done ++ [c], ?_⟩
          exact hi.step hwf hpre he (Or.inr (Or.inr ⟨hv, hn, rfl, rfl, rfl, rfl⟩)) hs1 hpost
      · intro b σ ⟨done, hi⟩
        try simp only [wp_modify]
        exact hi.fin hwf hpre he (by constructor <;> rfl)
    · rename_i hcap
      refine wp_abs _ (StepCap false s4.nodeEdges.size ct s6) _ (by constructor <;> rfl) fun s7 h7 => ?_
      have he := entry_of_steps hwf hpre ct h1 h2 h3 h4 hchildren h6 h7
        ⟨fun h => absurd h Bool.false_ne_true, fun h => absurd h hcap⟩
      try simp only [wp_bind]
      refine wp_forIn_inv _ _ _ _
        (fun rest b σ => ∃ done, LoopInv g items cur ct s.types.size s2.chDat.size s2.nodeVerts.size
          s4.nodeEdges.size s s7 children done rest b σ) s7 ⟨[], loop_init hpre he rfl (by rw [he.hasCap_eq, if_neg hcap])⟩ ?_ ?_
      · intro c hc rest b σ ⟨done, hi⟩
        have hcc := hi.mem_ch he
        simp only [wp_bind, wp_get, wp_modify, wp_pure, wp_ite]
        split_ifs with hv hn
        · refine wp_abs _ (StepSlotT (s2.chDat.size + b.2.2) σ.types.size none σ) _ (by constructor <;> rfl)
            fun σ1 hs1 => ?_
          refine wp_mono (ih c _ _ _ σ1 (hi.pre hwf hpre he (hs1.cons hi.cons) hs1.order)
            (loop_fuel hwf.tree hcur hcc hf)) fun _ σ3 hpost => ⟨_, rfl, done ++ [c], ?_⟩
          exact hi.step hwf hpre he (Or.inl ⟨hv, rfl, rfl, rfl, rfl⟩) hs1 hpost
        · refine wp_abs _ (StepSlotT (s2.chDat.size + b.2.2) σ.types.size (some b.2.1) σ) _ (by constructor <;> rfl)
            fun σ1 hs1 => ?_
          refine wp_mono (ih c _ _ _ σ1 (hi.pre hwf hpre he (hs1.cons hi.cons) hs1.order)
            (loop_fuel hwf.tree hcur hcc hf)) fun _ σ3 hpost => ⟨_, rfl, done ++ [c], ?_⟩
          exact hi.step hwf hpre he (Or.inr (Or.inl ⟨hv, hn, rfl, rfl, rfl, rfl⟩)) hs1 hpost
        · refine wp_abs _ (StepSlotT (s2.chDat.size + b.2.2) σ.types.size none σ) _ (by constructor <;> rfl)
            fun σ1 hs1 => ?_
          refine wp_mono (ih c _ _ _ σ1 (hi.pre hwf hpre he (hs1.cons hi.cons) hs1.order)
            (loop_fuel hwf.tree hcur hcc hf)) fun _ σ3 hpost => ⟨_, rfl, done ++ [c], ?_⟩
          exact hi.step hwf hpre he (Or.inr (Or.inr ⟨hv, hn, rfl, rfl, rfl, rfl⟩)) hs1 hpost
      · intro b σ ⟨done, hi⟩
        try simp only [wp_modify]
        exact hi.fin hwf hpre he (by constructor <;> rfl)

end Spqr.Ghost
