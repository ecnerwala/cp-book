import Spqr.PlanarRotInv

/-!
# `planarRelabel` preserves `RotInv`

`RotInv.step` is the per-node bookkeeping (the new node's `layoutRot` block lands at
`4 · neBounds[n]`); `planarRelabel_rotInv` threads it through the fold with the abstract-state
`wp`/`PlanarRelabelM.Frame` calculus of `PlanarRotInv.lean` (states are never unfolded past one primitive).
-/

namespace Spqr

open PlanarRelabelM

theorem RotInv.step {g : Graph} {s s' : PlanarRelabelState} (h : RotInv g s) (ty : NodeType)
    (nvEn neEn : Nat) (edgeVes : List Nat) (mapRot : Nat → Array (Option Nat))
    (hmr : ∀ ve, (mapRot ve).size = 4) (hR : ty = .R → edgeVes.length + 1 = neEn - s.base.nodeEdges.size)
    (hg : s'.base.g = g) (ht : s'.base.types = s.base.types.push ty)
    (hne : s'.base.neBounds = s.base.neBounds.push neEn) (hnv : s'.base.nvBounds = s.base.nvBounds.push nvEn)
    (hE : s'.base.nodeEdges.size = neEn) (hV : s'.base.nodeVerts.size = nvEn)
    (hrot : s'.aux.neRotAdj = s.aux.neRotAdj ++
      layoutRot ty (nvEn - s.base.nodeVerts.size) s.base.nodeEdges.size neEn edgeVes mapRot (2 * g.ne)) :
    RotInv g s' := by
  obtain ⟨hg0, h1, h2, h3, h4, h5, ev, mr, hmr0, hR0, hb⟩ := h
  have hn1 : s.base.types.size + 1 = s.base.neBounds.size := h1.symm
  have hn2 : s.base.types.size + 1 = s.base.nvBounds.size := h2.symm
  have e1 : (s.base.nvBounds.push nvEn)[s.base.types.size + 1]! = nvEn := by rw [hn2, Array.getElem!_push_eq']
  have e2 : (s.base.neBounds.push neEn)[s.base.types.size + 1]! = neEn := by rw [hn1, Array.getElem!_push_eq']
  refine ⟨hg, by rw [hne, ht]; simp [h1], by rw [hnv, ht]; simp [h2], ?_, ?_, ?_,
    fun m => if m = s.base.types.size then edgeVes else ev m,
    fun m => if m = s.base.types.size then mapRot else mr m, ?_, ?_, ?_⟩
  · rw [hne, Array.getElem!_push_lt' _ _ _ (by omega)]; exact h3
  · rw [ht, hne, Array.size_push, hE, e2]
  · rw [ht, hnv, Array.size_push, hV, e1]
  · intro m ve; by_cases hm : m = s.base.types.size <;> simp [hm, hmr, hmr0]
  · intro m hm hty
    rw [ht, Array.size_push] at hm
    dsimp only
    rcases Nat.lt_succ_iff_lt_or_eq.mp hm with hm' | hm'
    · rw [ht, Array.getElem!_push_lt' s.base.types _ m hm'] at hty
      rw [hne, Array.getElem!_push_lt' s.base.neBounds _ (m + 1) (by omega),
        Array.getElem!_push_lt' s.base.neBounds _ m (by omega), if_neg (by omega)]
      exact hR0 m hm' hty
    · subst hm'
      rw [ht, Array.getElem!_push_eq'] at hty
      rw [hne, e2, Array.getElem!_push_lt' s.base.neBounds _ s.base.types.size (by omega), h4, if_pos rfl]
      exact hR hty
  · rw [hrot, hb, ht, Array.size_push, concatBlocks.succ]
    congr 1
    · refine concatBlocks.congr fun m hm => ?_
      simp only [rotBlock, hne, hnv]
      rw [Array.getElem!_push_lt' s.base.types _ m hm, Array.getElem!_push_lt' s.base.nvBounds _ (m + 1) (by omega),
        Array.getElem!_push_lt' s.base.nvBounds _ m (by omega), Array.getElem!_push_lt' s.base.neBounds _ m (by omega),
        Array.getElem!_push_lt' s.base.neBounds _ (m + 1) (by omega), if_neg (show m ≠ _ by omega),
        if_neg (show m ≠ _ by omega)]
    · simp only [rotBlock, hne, hnv, ↓reduceIte]
      rw [Array.getElem!_push_eq' s.base.types, e1, e2,
        Array.getElem!_push_lt' s.base.nvBounds _ s.base.types.size (by omega),
        Array.getElem!_push_lt' s.base.neBounds _ s.base.types.size (by omega), h4, h5]

theorem RotInv.step' {g : Graph} {s s' : PlanarRelabelState} (h : RotInv g s) {ty : NodeType}
    {nV neSt neEn m : Nat} {edgeVes : List Nat} {mapRot : Nat → Array (Option Nat)} {X : Array (Option Nat)}
    (hrot : s'.aux.neRotAdj = X ++ layoutRot ty nV neSt neEn edgeVes mapRot m)
    (hX : X = s.aux.neRotAdj) (hnV : nV = s'.base.nodeVerts.size - s.base.nodeVerts.size)
    (hSt : neSt = s.base.nodeEdges.size) (hEn : neEn = s'.base.nodeEdges.size) (hm : m = 2 * g.ne)
    (hmr : ∀ ve, (mapRot ve).size = 4) (hR : ty = .R → edgeVes.length + 1 = neEn - neSt)
    (hg : s'.base.g = g) (ht : s'.base.types = s.base.types.push ty)
    (hne : s'.base.neBounds = s.base.neBounds.push neEn)
    (hnv : s'.base.nvBounds = s.base.nvBounds.push s'.base.nodeVerts.size) : RotInv g s' := by
  subst hX hnV hSt hEn hm
  exact RotInv.step h ty _ _ edgeVes mapRot hmr hR hg ht hne hnv rfl rfl hrot

theorem liftR_modify_snd (f : RelabelState → RelabelState) (s : PlanarRelabelState) :
    (liftR (modify f) s).2 = { s with base := f s.base } := by rw [liftR_eq]; rfl
theorem modifyAux_snd (f : PlanarRelabelAux → PlanarRelabelAux) (s : PlanarRelabelState) :
    (modifyAux f s).2 = { s with aux := f s.aux } := rfl

theorem wp_jp_match1 {Q : Unit → PlanarRelabelState → Prop} (ty : NodeType) (m₁ m₂ k : PlanarRelabelM Unit)
    (s : PlanarRelabelState) (h₁ : ∀ s, PlanarRelabelM.Frame s (m₁ s).2) (h₂ : ∀ s, PlanarRelabelM.Frame s (m₂ s).2)
    (hk : ∀ s', PlanarRelabelM.Frame s s' → wp k Q s') :
    wp (planarRelabel.match_1 (fun _ => PlanarRelabelM Unit) ty (fun _ => m₁ >>= fun _ => k)
      (fun _ => m₂ >>= fun _ => k) (fun _ => k)) Q s := by
  cases ty
  all_goals first
    | exact hk s (PlanarRelabelM.Frame.refl s)
    | exact wp_frame m₁ h₁ fun _ s' hf => hk s' hf
    | exact wp_frame m₂ h₂ fun _ s' hf => hk s' hf

theorem wp_jp_ite {Q : Unit → PlanarRelabelState → Prop} (c : Prop) [Decidable c] (m k : PlanarRelabelM Unit)
    (s : PlanarRelabelState) (h : ∀ s, PlanarRelabelM.Frame s (m s).2) (hk : ∀ s', PlanarRelabelM.Frame s s' → wp k Q s') :
    wp (if c then m >>= fun _ => k else k) Q s := by
  split
  · exact wp_frame m h fun _ s' hf => hk s' hf
  · exact hk s (PlanarRelabelM.Frame.refl s)

theorem wp_jp_ite3 {Q : Unit → PlanarRelabelState → Prop} (c : Prop) [Decidable c] (m₁ : PlanarRelabelM Unit)
    (m₂ : PlanarRelabelAux → PlanarRelabelM Unit) (k : PlanarRelabelM Unit) (s : PlanarRelabelState)
    (h₁ : ∀ s, PlanarRelabelM.Frame s (m₁ s).2) (h₂ : ∀ a s, PlanarRelabelM.Frame s (m₂ a s).2) (hk : ∀ s', PlanarRelabelM.Frame s s' → wp k Q s') :
    wp (if c then m₁ >>= fun _ => getAux >>= fun a => m₂ a >>= fun _ => k else k) Q s := by
  split
  · rw [wp_bind]
    refine wp_frame m₁ h₁ fun _ s₁ hf₁ => ?_
    simp only [wp_bind, wp_getAux]
    exact wp_frame (m₂ _) (h₂ _) fun _ s₂ hf₂ => hk s₂ (hf₁.trans hf₂)
  · exact hk s (PlanarRelabelM.Frame.refl s)

theorem planarRelabel_rotInv (g : Graph) :
    ∀ (fuel : Nat) (cur : ItemId) (p pn ct : Option Nat) (s : PlanarRelabelState),
      RotInv g s → RotInv g ((planarRelabel fuel cur p pn ct) s).2
  | 0, _, _, _, _, _, h => h
  | fuel + 1, cur, p, pn, ct, s, h => by
    have ih := planarRelabel_rotInv g fuel
    show wp (planarRelabel (fuel + 1) cur p pn ct) (fun _ s' => RotInv g s') s
    unfold planarRelabel
    simp only [wp_bind, wp_liftR_get, wp_liftR_item]
    -- number `cur`
    refine wp_abs _ fun s₁ hg₁ ht₁ hne₁ hnv₁ hE₁ hV₁ hR₁ => ?_
    simp only [liftR_modify_snd] at hg₁ ht₁ hne₁ hnv₁ hE₁ hV₁ hR₁
    try dsimp only
    -- V/Q bookkeeping
    refine wp_jp_match1 _ _ _ _ _ (fun s => PlanarRelabelM.Frame.liftR_modify _ s ⟨rfl, rfl, rfl, rfl, rfl, rfl, rfl⟩)
      (fun s => PlanarRelabelM.Frame.liftR_modify _ s ⟨rfl, rfl, rfl, rfl, rfl, rfl, rfl⟩) fun s₂ hf₂ => ?_
    try simp only [wp_bind]
    refine wp_frame _ (fun s => PlanarRelabelM.Frame.setupNode _ _ _ s) fun planar s₃ hf₃ => ?_
    try simp only [wp_bind]
    refine wp_frame _ (fun s => PlanarRelabelM.Frame.modifyAux _ s rfl) fun _ s₄ hf₄ => ?_
    try simp only [wp_liftR_get]
    -- node-verts
    refine wp_abs _ fun s₅ hg₅ ht₅ hne₅ hnv₅ hE₅ hV₅ hR₅ => ?_
    simp only [liftR_modify_snd, Array.size_append, List.size_toArray] at hg₅ ht₅ hne₅ hnv₅ hE₅ hV₅ hR₅
    try dsimp only
    -- R: vertex positions + flips
    refine wp_jp_ite3 _ _ _ _ _ (fun s => PlanarRelabelM.Frame.liftR_modify _ s ⟨rfl, rfl, rfl, rfl, rfl, rfl, rfl⟩)
      (fun _ s => PlanarRelabelM.Frame.applyFlips _ _ _ s) fun s₆ hf₆ => ?_
    try simp only [wp_bind]
    refine wp_liftR_orderedChildren _ _ fun children hch => ?_
    -- child slots
    refine wp_frame _ (fun s => PlanarRelabelM.Frame.liftR_modify _ s ⟨rfl, rfl, rfl, rfl, rfl, rfl, rfl⟩)
      fun _ s₇ hf₇ => ?_
    try simp only [wp_liftR_get]
    try dsimp only
    -- R: rotEdgeNe
    refine wp_jp_ite _ _ _ _ (fun s => PlanarRelabelM.Frame.modifyAux _ s rfl) fun s₈ hf₈ => ?_
    try simp only [wp_bind, wp_getAux]
    try dsimp only
    -- node edges / bounds
    refine wp_abs _ fun s₉ hg₉ ht₉ hne₉ hnv₉ hE₉ hV₉ hR₉ => ?_
    simp only [liftR_modify_snd, Array.size_append] at hg₉ ht₉ hne₉ hnv₉ hE₉ hV₉ hR₉
    -- rotation block
    refine wp_abs _ fun s₁₀ hg₁₀ ht₁₀ hne₁₀ hnv₁₀ hE₁₀ hV₁₀ hR₁₀ => ?_
    simp only [modifyAux_snd] at hg₁₀ ht₁₀ hne₁₀ hnv₁₀ hE₁₀ hV₁₀ hR₁₀
    try dsimp only
    have hinv₁₀ : RotInv g s₁₀ := by
      refine RotInv.step' h hR₁₀ ?hX ?hnV ?hSt ?hEn ?hm ?hmr ?hR ?hg ?ht ?hne ?hnv
      case hX =>
        simp only [hR₉, hf₈.neRotAdj, hf₇.neRotAdj, hf₆.neRotAdj, hR₅, hf₄.neRotAdj, hf₃.neRotAdj, hf₂.neRotAdj, hR₁]
      case hnV =>
        simp only [hV₁₀, hV₉, hf₈.nodeVerts, hf₇.nodeVerts, hf₆.nodeVerts, hV₅, hf₄.nodeVerts, hf₃.nodeVerts,
          hf₂.nodeVerts, hV₁]
      case hSt =>
        simp only [hf₇.nodeEdges, hf₆.nodeEdges, hE₅, hf₄.nodeEdges, hf₃.nodeEdges, hf₂.nodeEdges, hE₁]
      case hEn =>
        simp only [hE₁₀, hE₉, hf₈.nodeEdges, hf₇.nodeEdges, hf₆.nodeEdges, hE₅, hf₄.nodeEdges, hf₃.nodeEdges,
          hf₂.nodeEdges, hE₁]
        rw [(layoutNode_sz _ _ _ _ _ _ _).1]
        omega
      case hm => rw [h.g_eq]
      case hmr =>
        intro ve
        try dsimp only
        split <;> simp
      case hR =>
        intro hty
        rw [Nat.add_sub_cancel_left]
        simp [hty, NodeType.isNode, List.countP_eq_length_filter]
      case hg =>
        simp only [hg₁₀, hg₉, hf₈.g, hf₇.g, hf₆.g, hg₅, hf₄.g, hf₃.g, hf₂.g, hg₁, h.g_eq]
      case ht =>
        simp only [ht₁₀, ht₉, hf₈.types, hf₇.types, hf₆.types, ht₅, hf₄.types, hf₃.types, hf₂.types, ht₁]
      case hne =>
        simp only [hne₁₀, hne₉, hf₈.neBounds, hf₇.neBounds, hf₆.neBounds, hne₅, hf₄.neBounds, hf₃.neBounds,
          hf₂.neBounds, hne₁, hf₇.nodeEdges, hf₆.nodeEdges, hE₅, hf₄.nodeEdges,
          hf₃.nodeEdges, hf₂.nodeEdges, hE₁]
      case hnv =>
        simp only [hnv₁₀, hnv₉, hf₈.nvBounds, hf₇.nvBounds, hf₆.nvBounds, hnv₅, hf₄.nvBounds, hf₃.nvBounds,
          hf₂.nvBounds, hnv₁, hV₁₀, hV₉, hf₈.nodeVerts, hf₇.nodeVerts, hf₆.nodeVerts, hV₅, hf₄.nodeVerts,
          hf₃.nodeVerts, hf₂.nodeVerts, hV₁]
    -- cap twin
    refine wp_jp_ite _ _ _ _
      (fun s => PlanarRelabelM.Frame.liftR_modify _ s (by constructor <;> simp only [Array.size_modify])) fun s₁₁ hf₁₁ => ?_
    have hinv₁₁ := RotInv.frame hf₁₁ hinv₁₀
    try simp only [wp_bind]
    -- children
    refine wp_forIn_inv _ _ _ _ (fun _ _ s' => RotInv g s') _ hinv₁₁ ?_ ?_
    · intro c hc rest b s' hinv'
      simp only [wp_bind, wp_liftR_get, wp_liftR_modify, wp_pure, wp_ite]
      split_ifs
      all_goals
        unfold wp
        try dsimp only
        exact ⟨_, rfl, ih c _ _ _ _ (RotInv.frame (by constructor <;> simp only [Array.size_modify]) hinv')⟩
    · intro b s' hinv'
      try dsimp only
      exact wp_frame _ (fun s => PlanarRelabelM.Frame.liftR_modify _ s ⟨rfl, rfl, rfl, rfl, rfl, rfl, rfl⟩)
        fun _ s'' hf => RotInv.frame hf hinv'

end Spqr
