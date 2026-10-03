import Spqr.PlanarRelabelRows

/-!
# `RowInv` is preserved by `planarRelabel`

The relabel-side half of `planarRelabelTree_relabelNodeR`: every planar R node laid out by
`planarRelabel` satisfies `NodeRow`, and a subtree only writes `qem` slots below it (`Unch`), so the
slots of later siblings are still the walk's (`Fresh`).
-/

namespace Spqr

open Ghost PlanarRelabelM

theorem RFr.liftR_modify (f : RelabelState → RelabelState) (s : PlanarRelabelState)
    (hg : (f s.base).g = s.base.g) (hi : (f s.base).items = s.base.items) (ht : (f s.base).types = s.base.types)
    (hne : (f s.base).neBounds = s.base.neBounds) (hnv : (f s.base).nvBounds = s.base.nvBounds)
    (hE : (f s.base).nodeEdges = s.base.nodeEdges) (hV : (f s.base).nodeVerts.size = s.base.nodeVerts.size)
    (hvp : (f s.base).vertPos = s.base.vertPos) : RFr s (liftR (modify f) s).2 := by
  rw [liftR_modify_snd]
  exact ⟨hg, hi, ht, hne, hnv, by unfold nvsArr; simp only [hE], hV, hvp, rfl, rfl, rfl, rfl, rfl, rfl⟩

theorem foldl_vertPos_of (f : Array Nat → NodeVert × Nat → Array Nat)
    (hf : ∀ a nv pos, f a (nv, pos) = a.set! nv.vert pos) (vl : List Nat) (idx n : Nat) (a : Array Nat) :
    (((vl.map fun v => (⟨idx, v⟩ : NodeVert)).zipIdx n).foldl f a).size = a.size ∧
    ((∀ v ∈ vl, v < a.size) →
      Items.PosOK n vl fun v => (((vl.map fun v => (⟨idx, v⟩ : NodeVert)).zipIdx n).foldl f a)[v]!) := by
  have e : f = fun a x => a.set! x.1.vert x.2 := by funext a x; exact hf a x.1 x.2
  rw [e, Ghost.foldl_set!_map]
  have e2 : (vl.map fun v => (⟨idx, v⟩ : NodeVert)).map (·.vert) = vl := by simp [Function.comp_def]
  rw [e2]
  exact ⟨Ghost.foldl_set!_size _ _ _, fun h => Ghost.foldl_set!_posOK vl n a h⟩

/-- One child of the loop: the recursive call keeps `RowInv`, writes only below `cur`, and leaves the
later siblings' slots fresh. -/
theorem rowInv_child_step (g : Graph) (w : PlanarWalkState) (hwf : Items.WF g w.base.items) (fuel : Nat)
    (ih : ∀ (cur : ItemId) (p pn ct : Option Nat) (s : PlanarRelabelState),
      cur < w.base.items.size → RowInv g w s → Fresh g w cur s →
      RowInv g w ((planarRelabel fuel cur p pn ct) s).2 ∧
        Unch g w.base.items cur s ((planarRelabel fuel cur p pn ct) s).2)
    {cur c : ItemId} {children rest : List ItemId} (hperm : children.Perm (Items.ch w.base.items cur))
    (hsub : (c :: rest).Sublist children) {s s' X : PlanarRelabelState} (hinv' : RowInv g w s')
    (hu' : Unch g w.base.items cur s s') (hfr' : ∀ c' ∈ c :: rest, Fresh g w c' s') (hfX : RFr s' X)
    (p pn ct : Option Nat) :
    RowInv g w ((planarRelabel fuel c p pn ct) X).2 ∧
      Unch g w.base.items cur s ((planarRelabel fuel c p pn ct) X).2 ∧
      (∀ c' ∈ rest, Fresh g w c' ((planarRelabel fuel c p pn ct) X).2) ∧ rest.Sublist children := by
  have ht := hwf.tree
  have hcp : c ∈ Items.ch w.base.items cur := hperm.mem_iff.1 (hsub.subset (List.mem_cons_self ..))
  have hnd : (c :: rest).Nodup := hsub.nodup (hperm.nodup_iff.2 (ht.ch_nodup cur))
  obtain ⟨hinv'', hu''⟩ := ih c _ _ _ _ (ht.ch_lt _ _ hcp) (hinv'.fr hfX) ((hfr' c (List.mem_cons_self ..)).fr hfX)
  refine ⟨hinv'', hu'.trans ((Unch.rfr_left hfX hu'').of_child hcp), fun c' hc' => ?_,
    (List.sublist_cons_self c rest).trans hsub⟩
  have hne : c ≠ c' := fun e => (List.nodup_cons.1 hnd).1 (e ▸ hc')
  exact Fresh.of_unch ht hcp (hperm.mem_iff.1 (hsub.subset (List.mem_cons_of_mem _ hc'))) hne
    ((hfr' c' (List.mem_cons_of_mem _ hc')).fr hfX) hu''


theorem planarFinish_node (g : Graph) (w : PlanarWalkState) (hwf : Items.WF g w.base.items) (hpf : PlanarFinish g w)
    (cur : ItemId) (hcur : cur < w.base.items.size)
    (hnode : Items.type w.base.items cur = .S ∨ Items.type w.base.items cur = .P ∨ Items.type w.base.items cur = .R)
    (m : Array Nat) (hm : w.aux.nodePlanarity[cur - (1 + g.nv + g.ne)]! = .planar m) :
    1 + g.nv + g.ne ≤ cur ∧ w.aux.nodePlanarity[cur - (1 + g.nv + g.ne)]? = some (.planar m) ∧
      (m.size = 4 ∧ (∀ s, s < 4 → 1 + g.nv + QE.edge m[s]! ∈ Items.ch w.base.items cur) ∧ m.toList.Nodup) ∧
      (∀ c ∈ Items.ch w.base.items cur, c < 1 + g.nv + 2 * g.ne) := by
  have ht := hwf.tree
  have hge : 1 + g.nv + g.ne ≤ cur := ht.node_ge_of (by rcases hnode with h | h | h <;> rw [h] <;> decide)
    (by rcases hnode with h | h | h <;> rw [h] <;> decide) (by rcases hnode with h | h | h <;> rw [h] <;> decide)
  have hmk := nodePlanarity_getElem?_of _ _ _ hm
  have ecur : 1 + g.nv + g.ne + (cur - (1 + g.nv + g.ne)) = cur := by omega
  have hpfit : w.base.items[1 + g.nv + g.ne + (cur - (1 + g.nv + g.ne))]!.type ∈ [NodeType.S, .P, .R] := by
    rw [ecur, Items.getElem!_type hcur]; rcases hnode with h | h | h <;> rw [h] <;> simp
  have h1 := hpf.match_ends _ m hmk hpfit
  have h2 := hpf.ch_lt _ m hmk hpfit
  simp only [ecur, Items.getElem!_ch hcur] at h1 h2
  exact ⟨hge, hmk, h1, h2⟩

theorem setupSpec_keep (g : Graph) (w : PlanarWalkState) (hwf : Items.WF g w.base.items) (hpf : PlanarFinish g w)
    (cur : ItemId) (hcur : cur < w.base.items.size) (a : PlanarRelabelAux) (ha : a.nodePlanarity = w.aux.nodePlanarity)
    (planar : Bool) (q₃ : Qem) (hsp : SetupSpec g (Items.type w.base.items cur) cur a planar q₃) :
    q₃.size = a.qem.size ∧ ∀ q, q < 8 * g.ne →
      (∀ d ∈ Items.ch w.base.items cur, 1 + g.nv ≤ d → ¬ (4 * (d - (1 + g.nv)) ≤ q ∧ q < 4 * (d - (1 + g.nv)) + 4)) →
      q₃[q]! = a.qem[q]! := by
  obtain ⟨hsp1, hsp2, hsp3⟩ := hsp
  by_cases hnode : Items.type w.base.items cur = .S ∨ Items.type w.base.items cur = .P ∨
      Items.type w.base.items cur = .R
  · cases planar with
    | false => rw [hsp2 rfl]; exact ⟨rfl, fun _ _ _ => rfl⟩
    | true =>
      obtain ⟨m, hm, hq⟩ := hsp1 rfl hnode
      rw [ha] at hm
      obtain ⟨_, _, ⟨_, hm2, _⟩, _⟩ := planarFinish_node g w hwf hpf cur hcur hnode m hm
      rw [hq]
      refine ⟨size_capLinked _ _ _, fun q hq hd => capLinked_get!_of g _ m q hq fun s hs heq => ?_⟩
      exact hd _ (hm2 s hs) (by omega) (heq ▸ edge_slot_of _ s)
  · rw [(hsp3 hnode).2]; exact ⟨rfl, fun _ _ _ => rfl⟩

theorem edgeVes_facts (g : Graph) (items : Items) (ht : Items.Tree g items) (cur : ItemId) (children : List ItemId)
    (hperm : children.Perm (Items.ch items cur)) (hcl : ∀ c ∈ Items.ch items cur, c < 1 + g.nv + 2 * g.ne) :
    (∀ v ∈ (children.filter (· ≥ 1 + g.nv)).map (· - (1 + g.nv)),
      ∃ c ∈ Items.ch items cur, 1 + g.nv ≤ c ∧ c < 1 + g.nv + 2 * g.ne ∧ v = c - (1 + g.nv)) ∧
    ((children.filter (· ≥ 1 + g.nv)).map (· - (1 + g.nv))).Nodup := by
  refine ⟨fun v hv => ?_, ?_⟩
  · obtain ⟨c, hc, rfl⟩ := List.mem_map.1 hv
    obtain ⟨hc, hge'⟩ := List.mem_filter.1 hc
    have hc' := hperm.mem_iff.1 hc
    exact ⟨c, hc', by simpa using hge', hcl c hc', rfl⟩
  · refine List.Nodup.map_on ?_ (List.Nodup.filter _ (hperm.nodup_iff.2 (ht.ch_nodup cur)))
    intro x hx y hy hxy
    have hx' : 1 + g.nv ≤ x := by simpa using (List.mem_filter.1 hx).2
    have hy' : 1 + g.nv ≤ y := by simpa using (List.mem_filter.1 hy).2
    (try simp only [ItemId] at *); omega

theorem rotEdgeNe_node (g : Graph) (edgeVes : List Nat) (hnd : edgeVes.Nodup) (hlt : ∀ v ∈ edgeVes, v < 2 * g.ne)
    (r : Array Nat) (hr : r.size = 2 * g.ne + 1) (neSt : Nat) :
    ∀ v ∈ capVe g :: edgeVes,
      (((edgeVes.zipIdx (neSt + 1)).foldl (init := r) (fun r (ve, ne) => r.set! ve ne)).set! (2 * g.ne) neSt)[v]! =
        neSt + (capVe g :: edgeVes).idxOf v := by
  have hcap : 2 * g.ne ∉ edgeVes := fun hv => by have := hlt _ hv; omega
  obtain ⟨hc0, hrest⟩ := rotEdgeNe_init edgeVes neSt r g.ne hnd (fun v hv => by rw [hr]; have := hlt v hv; omega)
    (by omega) hcap
  intro v hv
  rcases List.mem_cons.1 hv with rfl | hv'
  · simp only [capVe, List.idxOf_cons_self, Nat.add_zero]; exact hc0
  · have hne : capVe g ≠ v := fun e => hcap (by rw [show 2 * g.ne = capVe g from rfl, e]; exact hv')
    rw [List.idxOf_cons_ne _ hne, hrest v hv']; omega

theorem layoutNode_nvs_R (ty : NodeType) (hty : ty = .R) (idx nvSt nvEn neSt neEn : Nat) (ec : List (Nat × Nat))
    (hnv : nvSt + 2 ≤ nvEn) (hne : neEn = neSt + ec.length + 1) (k : Nat) (hk : k < neEn - neSt) :
    ((layoutNode ty idx nvSt nvEn neSt neEn ec).edges.map (·.nvs))[k]! = ((nvSt, nvEn - 1) :: ec)[k]! := by
  subst hty
  have hsz := (layoutNode_sz .R idx nvSt nvEn neSt neEn ec).1
  rw [Array.getElem!_map' _ _ _ (by rw [hsz]; exact hk)]
  rw [layoutNode_edges _ _ _ _ _ _ _ rfl (by omega) (by rintro (h | h) <;> cases h) (fun _ => ⟨hnv, hne⟩)
    (by intro h; cases h) k hk]
  show layoutNvs _ _ _ _ k = _
  unfold layoutNvs
  rw [if_neg (by omega), if_neg (by rintro (h | h | h) <;> cases h), if_neg (by intro h; cases h)]
  cases k with
  | zero => simp
  | succ j =>
    rw [if_neg (Nat.succ_ne_zero _), List.getElem!_cons_succ, Nat.add_sub_cancel]
    by_cases hj : j < ec.length
    · rw [List.getElem?_eq_getElem hj, getElem!_pos ec j hj]; rfl
    · rw [List.getElem?_eq_none_iff.2 (Nat.le_of_not_lt hj), getElem!_neg ec j hj]; rfl

theorem mapRot_R_get (nV neSt neEn : Nat) (edgeVes : List Nat) (q : Qem) (r : Array Nat) (cv : Nat) (planar : Bool)
    (hpl : planar = true) (hnV : nV ≠ 1) (l : Nat) (hl : l < 4 * (edgeVes.length + 1)) :
    (layoutRot .R nV neSt neEn edgeVes (fun ve => if (!planar) = true then Array.replicate 4 none else
      Array.map (fun z => Option.map (fun o => 4 * r[QE.edge o]! + (o &&& 2) + (1 - z % 2)) q[4 * ve + z]!)
        (Array.range 4)) cv)[l]! =
    (q[4 * (cv :: edgeVes)[l / 4]! + l % 4]!).map fun o => 4 * r[QE.edge o]! + (o &&& 2) + (1 - l % 4 % 2) := by
  subst hpl
  rw [layoutRot_R_get _ _ _ _ _ _ hnV (fun ve => by simp) l hl]
  simp only [Bool.not_true, Bool.false_eq_true, ↓reduceIte]
  rw [Array.getElem!_map' _ _ _ (by simp; omega), range4_get! _ (Nat.mod_lt _ (by omega))]

theorem qemAgree_node (g : Graph) (w : PlanarWalkState) (cur : ItemId) (s : PlanarRelabelState) (hfr : Fresh g w cur s)
    (hsz : s.aux.qem.size = w.aux.qem.size) (hcl : ∀ c ∈ Items.ch w.base.items cur, c < 1 + g.nv + 2 * g.ne)
    (m : Array Nat) (fl : List Bool) :
    QemAgree (fun q => (∃ c ∈ Items.ch w.base.items cur, 1 + g.nv ≤ c ∧
        4 * (c - (1 + g.nv)) ≤ q ∧ q < 4 * (c - (1 + g.nv)) + 4) ∨ ∃ s, s < 4 ∧ q = 4 * capVe g + s)
      (flipped g (Items.ch w.base.items cur) fl (capLinked g s.aux.qem m))
      (flipped g (Items.ch w.base.items cur) fl (capLinked g w.aux.qem m)) := by
  refine qemAgree_flipped g _ _ (fun c hc hge z hz => Or.inl ⟨c, hc, hge, ?_, ?_⟩) ?_
  · (try simp only [ItemId] at *); omega
  · (try simp only [ItemId] at *); omega
  refine qemAgree_capLinked g m ⟨hsz, fun q hq => ?_⟩
  obtain ⟨c, hc, hge', h1, h2⟩ := hq
  obtain ⟨z, hz, rfl⟩ : ∃ z, z < 4 ∧ q = 4 * (c - (1 + g.nv)) + z :=
    ⟨q - 4 * (c - (1 + g.nv)), by (try simp only [ItemId] at *); omega, by (try simp only [ItemId] at *); omega⟩
  exact hfr cur c Relation.ReflTransGen.refl hc hge' (hcl c hc) z hz

theorem planarRelabel_rowInv (g : Graph) (w : PlanarWalkState) (hwf : Items.WF g w.base.items)
    (hpf : PlanarFinish g w) :
    ∀ (fuel : Nat) (cur : ItemId) (p pn ct : Option Nat) (s : PlanarRelabelState),
      cur < w.base.items.size → RowInv g w s → Fresh g w cur s →
      RowInv g w ((planarRelabel fuel cur p pn ct) s).2 ∧
        Unch g w.base.items cur s ((planarRelabel fuel cur p pn ct) s).2
  | 0, _, _, _, _, s, _, h, _ => ⟨h, Unch.refl s⟩
  | fuel + 1, cur, p, pn, ct, s, hcur, h, hfr => by
    have ih := planarRelabel_rowInv g w hwf hpf fuel
    have ht := hwf.tree
    show wp (planarRelabel (fuel + 1) cur p pn ct) (fun _ s' => RowInv g w s' ∧ Unch g w.base.items cur s s') s
    unfold planarRelabel
    simp only [wp_bind, wp_liftR_get, wp_liftR_item]
    rw [h.g_eq, h.items]
    simp only [Items.getElem!_type hcur, Items.getElem!_ch hcur, Items.getElem!_vs_all]
    -- number `cur`
    refine wp_rabs _ fun s₁ hg₁ hi₁ ht₁ hne₁ hnv₁ hE₁ hV₁ hvp₁ hq₁ hnp₁ hfl₁ hrne₁ hR₁ hN₁ => ?_
    simp only [liftR_modify_snd] at hg₁ hi₁ ht₁ hne₁ hnv₁ hE₁ hV₁ hvp₁ hq₁ hnp₁ hfl₁ hrne₁ hR₁ hN₁
    try dsimp only
    -- V/Q bookkeeping
    refine wp_jp_match1_rfr _ _ _ _ _ (fun s => RFr.liftR_modify _ s rfl rfl rfl rfl rfl rfl rfl rfl)
      (fun s => RFr.liftR_modify _ s rfl rfl rfl rfl rfl rfl rfl rfl) fun s₂ hf₂ => ?_
    try simp only [wp_bind]
    -- setupNode
    refine wp_setupNode _ _ _ _ fun planar q₃ hsp => ?_
    try simp only [wp_bind]
    -- nodePlanar
    refine wp_rabs _ fun s₄ hg₄ hi₄ ht₄ hne₄ hnv₄ hE₄ hV₄ hvp₄ hq₄ hnp₄ hfl₄ hrne₄ hR₄ hN₄ => ?_
    simp only [modifyAux_snd] at hg₄ hi₄ ht₄ hne₄ hnv₄ hE₄ hV₄ hvp₄ hq₄ hnp₄ hfl₄ hrne₄ hR₄ hN₄
    try simp only [wp_liftR_get]
    -- node-verts
    refine wp_rabs _ fun s₅ hg₅ hi₅ ht₅ hne₅ hnv₅ hE₅ hV₅ hvp₅ hq₅ hnp₅ hfl₅ hrne₅ hR₅ hN₅ => ?_
    simp only [liftR_modify_snd, Array.size_append, List.size_toArray] at hg₅ hi₅ ht₅ hne₅ hnv₅ hE₅ hV₅ hvp₅ hq₅ hnp₅ hfl₅ hrne₅ hR₅ hN₅
    try dsimp only
    -- R: vertex positions + flips
    refine wp_ite_jp3 _ _ _ _ _ fun s₆ h₆ => ?_
    try simp only [wp_bind]
    refine wp_liftR_orderedChildren _ _ fun children hch => ?_
    -- child slots
    refine wp_rfr _ (fun s => RFr.liftR_modify _ s rfl rfl rfl rfl rfl rfl rfl rfl) fun _ s₇ hf₇ => ?_
    try simp only [wp_liftR_get]
    try dsimp only
    -- R: rotEdgeNe
    refine wp_ite_jp _ _ _ _ fun s₈ h₈ => ?_
    try simp only [wp_bind, wp_getAux]
    try dsimp only
    -- node edges / bounds
    refine wp_rabs _ fun s₉ hg₉ hi₉ ht₉ hne₉ hnv₉ hE₉ hV₉ hvp₉ hq₉ hnp₉ hfl₉ hrne₉ hR₉ hN₉ => ?_
    simp only [liftR_modify_snd] at hg₉ hi₉ ht₉ hne₉ hnv₉ hE₉ hV₉ hvp₉ hq₉ hnp₉ hfl₉ hrne₉ hR₉ hN₉
    -- rotation block
    refine wp_rabs _ fun s₁₀ hg₁₀ hi₁₀ ht₁₀ hne₁₀ hnv₁₀ hE₁₀ hV₁₀ hvp₁₀ hq₁₀ hnp₁₀ hfl₁₀ hrne₁₀ hR₁₀ hN₁₀ => ?_
    simp only [modifyAux_snd] at hg₁₀ hi₁₀ ht₁₀ hne₁₀ hnv₁₀ hE₁₀ hV₁₀ hvp₁₀ hq₁₀ hnp₁₀ hfl₁₀ hrne₁₀ hR₁₀ hN₁₀
    try dsimp only
    -- cap twin
    refine wp_ite_jp _ _ _ _ fun s₁₁ h₁₁ => ?_
    have hf₁₁ : RFr s₁₀ s₁₁ := by
      rcases h₁₁ with ⟨_, rfl⟩ | ⟨_, rfl⟩
      · rw [liftR_modify_snd]
        exact ⟨rfl, rfl, rfl, rfl, rfl, nvsArr_modify_twin _ _ _, by simp only [Array.size_modify], rfl, rfl,
          rfl, rfl, rfl, rfl, rfl⟩
      · exact RFr.refl _
    have hperm : children.Perm (Items.ch w.base.items cur) := by
      rw [hch, ← Items.getElem!_ch hcur]
      by_cases hR : w.base.items[cur]!.type = .R
      · have e : (w.base.items[cur]!.type != NodeType.R) = false := by simpa using hR
        simp only [RelabelM.orderedChildren, e]
        exact List.mergeSort_perm _ _
      · have e : (w.base.items[cur]!.type != NodeType.R) = true := by simpa using hR
        simp only [RelabelM.orderedChildren, e]
        exact List.Perm.refl _
    -- fold the node-vertex list
    have hnvl : ((Items.vs w.base.items cur).1.toList ++ ((Items.ch w.base.items cur).filter (· < 1 + g.nv)).map (· - 1) ++
        (Items.vs w.base.items cur).2.toList) = Items.nvList g w.base.items cur := rfl
    rw [hnvl] at hV₅ hnv₉ hE₉ hR₁₀ h₆
    simp only [List.length_map, Nat.add_sub_cancel_left] at hV₅ hnv₉ hE₉ hR₁₀
    rw [hf₇.vertPos] at hE₉
    have e7i : s₇.base.items = w.base.items := by
      rw [hf₇.items]
      rcases h₆ with ⟨_, rfl⟩ | ⟨_, rfl⟩
      · simp only [applyFlips_eq, liftR_modify_snd]; rw [hi₅, hi₄, hf₂.items, hi₁, h.items]
      · rw [hi₅, hi₄, hf₂.items, hi₁, h.items]
    rw [e7i] at hE₉
    -- the R-only steps
    obtain ⟨h6g, h6i, h6t, h6ne, h6nv, h6E, h6V, h6vpe, h6np, h6fl, h6rne, h6R, h6N, h6q, h6pos⟩ :
        s₆.base.g = s₅.base.g ∧ s₆.base.items = s₅.base.items ∧ s₆.base.types = s₅.base.types ∧
        s₆.base.neBounds = s₅.base.neBounds ∧ s₆.base.nvBounds = s₅.base.nvBounds ∧
        s₆.base.nodeEdges = s₅.base.nodeEdges ∧ s₆.base.nodeVerts.size = s₅.base.nodeVerts.size ∧
        s₆.base.vertPos.size = s₅.base.vertPos.size ∧ s₆.aux.nodePlanarity = s₅.aux.nodePlanarity ∧
        s₆.aux.itemFlips = s₅.aux.itemFlips ∧ s₆.aux.rotEdgeNe = s₅.aux.rotEdgeNe ∧
        s₆.aux.neRotAdj = s₅.aux.neRotAdj ∧ s₆.aux.nodePlanar = s₅.aux.nodePlanar ∧
        s₆.aux.qem = (if (Items.type w.base.items cur == NodeType.R) = true then
          flipped g (Items.ch w.base.items cur) s₅.aux.itemFlips[cur]! s₅.aux.qem else s₅.aux.qem) ∧
        ((Items.type w.base.items cur == NodeType.R) = true →
          Items.PosOK s₄.base.nodeVerts.size (Items.nvList g w.base.items cur) fun v => s₆.base.vertPos[v]!) := by
      have hlt : ∀ v ∈ Items.nvList g w.base.items cur, v < s₅.base.vertPos.size := by
        rw [hvp₅, hvp₄, hf₂.vertPos, hvp₁, h.vp_size]; exact hwf.nvList_lt cur
      rcases h₆ with ⟨hc, rfl⟩ | ⟨hc, rfl⟩
      · simp only [applyFlips_eq, liftR_modify_snd, hc, ↓reduceIte, true_implies, Items.getElem!_ch hcur]
        refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
        all_goals first
          | trivial
          | rfl
          | exact (foldl_vertPos_of _ (fun _ _ _ => rfl) _ _ _ _).1
          | exact (foldl_vertPos_of _ (fun _ _ _ => rfl) _ _ _ _).2 hlt
      · simp only [hc, Bool.false_eq_true, ↓reduceIte, false_implies]
        refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
        all_goals first | trivial | rfl
    obtain ⟨h8b, h8q, h8np, h8fl, h8R, h8N, h8rsz, h8rne⟩ :
        s₈.base = s₇.base ∧ s₈.aux.qem = s₇.aux.qem ∧ s₈.aux.nodePlanarity = s₇.aux.nodePlanarity ∧
        s₈.aux.itemFlips = s₇.aux.itemFlips ∧ s₈.aux.neRotAdj = s₇.aux.neRotAdj ∧
        s₈.aux.nodePlanar = s₇.aux.nodePlanar ∧ s₈.aux.rotEdgeNe.size = s₇.aux.rotEdgeNe.size ∧
        ((Items.type w.base.items cur == NodeType.R) = true → s₈.aux.rotEdgeNe =
          ((((children.filter (· ≥ 1 + g.nv)).map (· - (1 + g.nv))).zipIdx (s₇.base.nodeEdges.size + 1)).foldl
            (init := s₇.aux.rotEdgeNe) (fun r (ve, ne) => r.set! ve ne)).set! (2 * g.ne) s₇.base.nodeEdges.size) := by
      rcases h₈ with ⟨hc, rfl⟩ | ⟨hc, rfl⟩
      · simp only [modifyAux_snd, hc, true_implies]
        refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
        all_goals first | trivial | rfl | exact rotEdgeNe_init_size _ _ _ _
      · simp only [hc, Bool.false_eq_true, false_implies]
        refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
        all_goals first | trivial | rfl
    -- the fields of `s₁₁`
    have e2q : s₂.aux.qem = s.aux.qem := by rw [hf₂.qem, hq₁]
    have e2np : s₂.aux.nodePlanarity = w.aux.nodePlanarity := by rw [hf₂.np, hnp₁, h.np]
    have e5q : s₅.aux.qem = q₃ := by rw [hq₅, hq₄]
    have e5fl : s₅.aux.itemFlips = w.aux.itemFlips := by rw [hfl₅, hfl₄, hf₂.fl, hfl₁, h.fl]
    have e11q : s₁₁.aux.qem = s₆.aux.qem := by rw [hf₁₁.qem, hq₁₀, hq₉, h8q, hf₇.qem]
    have e11g : s₁₁.base.g = g := by
      rw [hf₁₁.g, hg₁₀, hg₉, h8b, hf₇.g, h6g, hg₅, hg₄, hf₂.g, hg₁, h.g_eq]
    have e11i : s₁₁.base.items = w.base.items := by
      rw [hf₁₁.items, hi₁₀, hi₉, h8b, hf₇.items, h6i, hi₅, hi₄, hf₂.items, hi₁, h.items]
    have e6g : s₆.base.g = g := by rw [h6g, hg₅, hg₄, hf₂.g, hg₁, h.g_eq]
    have e6i : s₆.base.items = w.base.items := by rw [h6i, hi₅, hi₄, hf₂.items, hi₁, h.items]
    have e11np : s₁₁.aux.nodePlanarity = w.aux.nodePlanarity := by
      rw [hf₁₁.np, hnp₁₀, hnp₉, h8np, hf₇.np, h6np, hnp₅, hnp₄, hf₂.np, hnp₁, h.np]
    have e11fl : s₁₁.aux.itemFlips = w.aux.itemFlips := by
      rw [hf₁₁.fl, hfl₁₀, hfl₉, h8fl, hf₇.fl, h6fl, hfl₅, hfl₄, hf₂.fl, hfl₁, h.fl]
    have e7rne : s₇.aux.rotEdgeNe = s.aux.rotEdgeNe := by rw [hf₇.rne, h6rne, hrne₅, hrne₄, hf₂.rne, hrne₁]
    have e11rne : s₁₁.aux.rotEdgeNe.size = 2 * g.ne + 1 := by
      rw [hf₁₁.rne, hrne₁₀, hrne₉, h8rsz, e7rne, h.rne_size]
    have e11T : s₁₁.base.types = s.base.types.push (Items.type w.base.items cur) := by
      rw [hf₁₁.types, ht₁₀, ht₉, h8b, hf₇.types, h6t, ht₅, ht₄, hf₂.types, ht₁]
    have e11N : s₁₁.aux.nodePlanar = s.aux.nodePlanar.push planar := by
      rw [hf₁₁.nodePlanar, hN₁₀, hN₉, h8N, hf₇.nodePlanar, h6N, hN₅, hN₄, hf₂.nodePlanar, hN₁]
    have envSt : s₄.base.nodeVerts.size = s.base.nodeVerts.size := by rw [hV₄, hf₂.nodeVerts, hV₁]
    have eneSt : s₇.base.nodeEdges.size = s.base.nodeEdges.size := by
      rw [← nvsArr_size, hf₇.nvs, nvsArr_size, h6E, hE₅, hE₄, ← nvsArr_size, hf₂.nvs, nvsArr_size, hE₁]
    have e11nvB : s₁₁.base.nvBounds = s.base.nvBounds.push
        (s₄.base.nodeVerts.size + (Items.nvList g w.base.items cur).length) := by
      rw [hf₁₁.nvBounds, hnv₁₀, hnv₉, h8b, hf₇.nvBounds, h6nv, hnv₅, hnv₄, hf₂.nvBounds, hnv₁]
    have e11V : s₁₁.base.nodeVerts.size = s.base.nodeVerts.size + (Items.nvList g w.base.items cur).length := by
      rw [hf₁₁.nodeVerts, hV₁₀, hV₉, h8b, hf₇.nodeVerts, h6V, hV₅, envSt]
    have e11vp : s₁₁.base.vertPos.size = g.nv := by
      rw [hf₁₁.vertPos, hvp₁₀, hvp₉, h8b, hf₇.vertPos, h6vpe, hvp₅, hvp₄, hf₂.vertPos, hvp₁, h.vp_size]
    have nE_eq : (if (Items.type w.base.items cur).isNode = true then children.countP (· ≥ 1 + g.nv) else 0) +
        (if ((Items.type w.base.items cur).isNode &&
          !(Items.type w.base.items cur == NodeType.Q && !children.isEmpty)) = true then 1 else 0) =
        Items.nEdges g w.base.items cur := by
      rw [Items.nEdges, Items.capCount, Items.hasCap, hperm.countP_eq, hperm.isEmpty_eq]
    rw [nE_eq] at hne₉ hE₉ hR₁₀
    have e11neB : s₁₁.base.neBounds = s.base.neBounds.push
        (s₇.base.nodeEdges.size + Items.nEdges g w.base.items cur) := by
      rw [hf₁₁.neBounds, hne₁₀, hne₉, h8b, hf₇.neBounds, h6ne, hne₅, hne₄, hf₂.neBounds, hne₁]
    have e8nvs : nvsArr s₈ = nvsArr s := by
      rw [show nvsArr s₈ = nvsArr s₇ by unfold nvsArr; rw [h8b], hf₇.nvs]
      unfold nvsArr; rw [h6E, hE₅, hE₄]
      exact hf₂.nvs.trans (by unfold nvsArr; rw [hE₁])
    have hev : ((children.filter (· ≥ 1 + g.nv)).map (· - (1 + g.nv))).length =
        (Items.ch w.base.items cur).countP (· ≥ 1 + g.nv) := by
      rw [List.length_map, ← List.countP_eq_length_filter, hperm.countP_eq]
    rw [show (List.map (fun c => (s₆.base.vertPos[(Items.vs w.base.items c).1.getD 0]!,
        s₆.base.vertPos[(Items.vs w.base.items c).2.getD 0]!)) (List.filter (fun x => decide (x ≥ 1 + g.nv)) children)) =
        Items.edgeChildren g w.base.items (fun v => s₆.base.vertPos[v]!) children from rfl] at hE₉
    have e11nvs : nvsArr s₁₁ = nvsArr s ++ (layoutNode (Items.type w.base.items cur) s.base.types.size
        s₄.base.nodeVerts.size (s₄.base.nodeVerts.size + (Items.nvList g w.base.items cur).length)
        s₇.base.nodeEdges.size (s₇.base.nodeEdges.size + Items.nEdges g w.base.items cur)
        (Items.edgeChildren g w.base.items (fun v => s₆.base.vertPos[v]!) children)).edges.map (·.nvs) := by
      rw [hf₁₁.nvs]
      show s₁₀.base.nodeEdges.map (·.nvs) = _
      rw [hE₁₀, hE₉, Array.map_append]
      show nvsArr s₈ ++ _ = nvsArr s ++ _
      rw [e8nvs]
    have hEsz : ((layoutNode (Items.type w.base.items cur) s.base.types.size
        s₄.base.nodeVerts.size (s₄.base.nodeVerts.size + (Items.nvList g w.base.items cur).length)
        s₇.base.nodeEdges.size (s₇.base.nodeEdges.size + Items.nEdges g w.base.items cur)
        (Items.edgeChildren g w.base.items (fun v => s₆.base.vertPos[v]!) children)).edges.map (·.nvs)).size =
        Items.nEdges g w.base.items cur := by
      rw [Array.size_map, (layoutNode_sz _ _ _ _ _ _ _).1, Nat.add_sub_cancel_left]
    have e11E : s₁₁.base.nodeEdges.size = s.base.nodeEdges.size + Items.nEdges g w.base.items cur := by
      rw [← nvsArr_size, e11nvs, Array.size_append, nvsArr_size, hEsz]
    have e11rot : s₁₁.aux.neRotAdj = s.aux.neRotAdj ++ layoutRot (Items.type w.base.items cur)
        (Items.nvList g w.base.items cur).length s₇.base.nodeEdges.size
        (s₇.base.nodeEdges.size + Items.nEdges g w.base.items cur)
        (List.map (fun x => x - (1 + g.nv)) (List.filter (fun x => decide (x ≥ 1 + g.nv)) children))
        (fun ve => if (!planar) = true then Array.replicate 4 none else
          Array.map (fun z => Option.map (fun o => 4 * s₈.aux.rotEdgeNe[QE.edge o]! + (o &&& 2) + (1 - z % 2))
            s₈.aux.qem[4 * ve + z]!) (Array.range 4)) (2 * g.ne) := by
      rw [hf₁₁.neRotAdj, hR₁₀, hR₉, h8R, hf₇.neRotAdj, h6R, hR₅, hR₄, hf₂.neRotAdj, hR₁]
    have hmr : ∀ ve, ((fun ve => if (!planar) = true then Array.replicate 4 none else
        Array.map (fun z => Option.map (fun o => 4 * s₈.aux.rotEdgeNe[QE.edge o]! + (o &&& 2) + (1 - z % 2))
          s₈.aux.qem[4 * ve + z]!) (Array.range 4)) ve).size = 4 := by
      intro ve; dsimp only; split <;> simp
    have hRsz : (layoutRot (Items.type w.base.items cur)
        (Items.nvList g w.base.items cur).length s₇.base.nodeEdges.size
        (s₇.base.nodeEdges.size + Items.nEdges g w.base.items cur)
        (List.map (fun x => x - (1 + g.nv)) (List.filter (fun x => decide (x ≥ 1 + g.nv)) children))
        (fun ve => if (!planar) = true then Array.replicate 4 none else
          Array.map (fun z => Option.map (fun o => 4 * s₈.aux.rotEdgeNe[QE.edge o]! + (o &&& 2) + (1 - z % 2))
            s₈.aux.qem[4 * ve + z]!) (Array.range 4)) (2 * g.ne)).size = 4 * Items.nEdges g w.base.items cur := by
      by_cases hn : (Items.type w.base.items cur).isNode = true
      · exact layoutRot_size_node hwf hcur hn _ _ _ _ hmr hev
      · have hFV : Items.type w.base.items cur = .F ∨ Items.type w.base.items cur = .V := by
          revert hn; cases Items.type w.base.items cur <;> simp [NodeType.isNode]
        have h0 : Items.nEdges g w.base.items cur = 0 := by
          rw [← nE_eq]; rcases hFV with h | h <;> rw [h] <;> rfl
        rw [h0]
        rcases hFV with h | h <;> rw [h] <;> simp [layoutRot]
    -- the writes to `qem`
    have hq3 := setupSpec_keep g w hwf hpf cur hcur s₂.aux e2np planar q₃ hsp
    rw [e2q] at hq3
    have hqsz : s₁₁.aux.qem.size = s.aux.qem.size := by
      rw [e11q, h6q]; split
      · rw [size_flipped, e5q, hq3.1]
      · rw [e5q, hq3.1]
    have hkeep : ∀ q, q < 8 * g.ne →
        (∀ d ∈ Items.ch w.base.items cur, 1 + g.nv ≤ d → ¬ (4 * (d - (1 + g.nv)) ≤ q ∧ q < 4 * (d - (1 + g.nv)) + 4)) →
        s₁₁.aux.qem[q]! = s.aux.qem[q]! := by
      intro q hq hd
      rw [e11q, h6q]; split
      · rw [flipped_get!_of g _ _ _ q hd, e5q, hq3.2 q hq hd]
      · rw [e5q, hq3.2 q hq hd]
    have hunch : Unch g w.base.items cur s s₁₁ := by
      refine ⟨hqsz, fun q hq => hkeep q ?_ ?_⟩
      · by_contra h'; exact hq (Or.inl (by omega))
      · intro d hd hge hin; exact hq (Or.inr ⟨cur, d, Relation.ReflTransGen.refl, hd, hge, hin.1, hin.2⟩)
    have hfresh : ∀ c ∈ children, Fresh g w c s₁₁ := by
      intro c hc j c' hj hc' hge hlt z hz
      have hcc : c ∈ Items.ch w.base.items cur := hperm.mem_iff.1 hc
      rw [← hfr j c' (Relation.ReflTransGen.head hcc hj) hc' hge hlt z hz]
      apply hkeep _ (by (try simp only [ItemId] at hlt hge hz ⊢); omega)
      intro d hd hdge hin
      exact ht.child_ne_of_below hcc hj hc' hd (slot_eq_of hge hdge hz hin)
    -- `RowInv`
    have hrow : RowInv g w s₁₁ := by
      have gnv0 : s₁₁.base.nvBounds[s.base.types.size]! = s₄.base.nodeVerts.size := by
        rw [e11nvB, Array.getElem!_push_lt' _ _ _ (by rw [h.nv_size]; omega), h.nv_last, envSt]
      have gnv1 : s₁₁.base.nvBounds[s.base.types.size + 1]! =
          s₄.base.nodeVerts.size + (Items.nvList g w.base.items cur).length := by
        rw [e11nvB]
        have := Array.getElem!_push_eq' s.base.nvBounds
          (s₄.base.nodeVerts.size + (Items.nvList g w.base.items cur).length)
        rwa [h.nv_size] at this
      have gne0 : s₁₁.base.neBounds[s.base.types.size]! = s₇.base.nodeEdges.size := by
        rw [e11neB, Array.getElem!_push_lt' _ _ _ (by rw [h.ne_size]; omega), h.ne_last, eneSt]
      have gne1 : s₁₁.base.neBounds[s.base.types.size + 1]! =
          s₇.base.nodeEdges.size + Items.nEdges g w.base.items cur := by
        rw [e11neB]
        have := Array.getElem!_push_eq' s.base.neBounds (s₇.base.nodeEdges.size + Items.nEdges g w.base.items cur)
        rwa [h.ne_size] at this
      refine ⟨e11g, e11i, e11np, e11fl, e11rne, ?_, ?_, ?_, ?_, hqsz.trans h.qem_size, e11vp, ?_, ?_, ?_⟩
      · rw [e11neB, e11T, Array.size_push, Array.size_push, h.ne_size]
      · rw [e11nvB, e11T, Array.size_push, Array.size_push, h.nv_size]
      · rw [e11N, e11T, Array.size_push, Array.size_push, h.np_size]
      · rw [e11rot, Array.size_append, hRsz, h.rot_size, e11E]; omega
      · rw [e11T, Array.size_push, gne1, e11E, eneSt]
      · rw [e11T, Array.size_push, gnv1, e11V, envSt]
      intro n hn hR hpl
      rw [e11T, Array.size_push] at hn
      rw [e11T] at hR
      rw [e11N] at hpl
      rcases Nat.lt_succ_iff_lt_or_eq.1 hn with hlt | rfl
      · rw [Array.getElem!_push_lt' _ _ _ hlt] at hR
        rw [Array.getElem!_push_lt' _ _ _ (by rw [h.np_size]; exact hlt)] at hpl
        exact (h.rows n hlt hR hpl).push e11nvB e11neB e11nvs e11rot (by rw [h.nv_size]; omega)
          (by rw [h.ne_size]; omega)
      · rw [Array.getElem!_push_eq'] at hR
        have := Array.getElem!_push_eq' s.aux.nodePlanar planar
        rw [h.np_size] at this
        rw [this] at hpl
        -- the new R row
        have hRb : (Items.type w.base.items cur == NodeType.R) = true := by rw [hR]; rfl
        obtain ⟨hsp1, -, -⟩ := hsp
        obtain ⟨m, hm, hq₃⟩ := hsp1 hpl (Or.inr (Or.inr hR))
        rw [e2np] at hm
        obtain ⟨hge, hmk, ⟨hm1, hm2, hm3⟩, hcl⟩ := planarFinish_node g w hwf hpf cur hcur (Or.inr (Or.inr hR)) m hm
        have hn : (Items.type w.base.items cur).isNode = true := by rw [hR]; rfl
        have hL := (hwf.layout_hyps hcur hn).2.2.2.1 hR
        have hnER : Items.nEdges g w.base.items cur =
            ((children.filter (· ≥ 1 + g.nv)).map (· - (1 + g.nv))).length + 1 := by
          rw [← nE_eq, hR]
          simp only [NodeType.isNode, List.length_map, List.countP_eq_length_filter]; rfl
        have hord : children = Items.ordered g w.base.items cur s₄.base.nodeVerts.size
            fun v => s₆.base.vertPos[v]! := by
          rw [hch, orderedChildren_eq_ordered g w.base.items cur hcur _ s₆.base e6g e6i]
        obtain ⟨hves_mem, hnd⟩ := edgeVes_facts g w.base.items ht cur children hperm hcl
        have hrne := rotEdgeNe_node g _ hnd (fun v hv => by
            obtain ⟨c, _, hge', hc, rfl⟩ := hves_mem v hv; (try simp only [ItemId] at hc hge' ⊢); omega)
          s.aux.rotEdgeNe h.rne_size s₇.base.nodeEdges.size
        rw [← e7rne, ← h8rne hRb] at hrne
        have hagree := qemAgree_node g w cur s hfr h.qem_size hcl m w.aux.itemFlips[cur]!
        have e8q : s₈.aux.qem = flipped g (Items.ch w.base.items cur) w.aux.itemFlips[cur]!
            (capLinked g s.aux.qem m) := by
          rw [h8q, hf₇.qem, h6q, if_pos hRb, e5fl, e5q, hq₃, e2q]
        unfold NodeRow
        refine ⟨cur, m, children, fun v => s₆.base.vertPos[v]!, ?_⟩
        dsimp only
        rw [gnv0, gnv1, gne0, gne1, Items.getElem!_type hcur, Items.getElem!_ch hcur]
        refine ⟨hge, hcur, hR, hmk, hord, h6pos hRb, rfl, by omega, ?_, ?_, ?_, ?_⟩
        · rw [e11nvs, Array.size_append, nvsArr_size, hEsz, ← eneSt]
        · intro k hk
          have hk' : k < ((layoutNode (Items.type w.base.items cur) s.base.types.size
              s₄.base.nodeVerts.size (s₄.base.nodeVerts.size + (Items.nvList g w.base.items cur).length)
              s₇.base.nodeEdges.size (s₇.base.nodeEdges.size + Items.nEdges g w.base.items cur)
              (Items.edgeChildren g w.base.items (fun v => s₆.base.vertPos[v]!) children)).edges.map (·.nvs)).size := by
            rw [hEsz]; omega
          rw [e11nvs, show s₇.base.nodeEdges.size + k = (nvsArr s).size + k by rw [nvsArr_size, eneSt],
            Array.getElem!_append_right' _ _ _ hk']
          exact layoutNode_nvs_R _ hR _ _ _ _ _ _ (by omega)
            (by rw [hnER]; simp only [Items.edgeChildren, List.length_map]; omega) k hk
        · rw [e11rot, Array.size_append, hRsz, h.rot_size, eneSt]; omega
        · intro l hl
          have hrot : s₁₁.aux.neRotAdj[4 * s₇.base.nodeEdges.size + l]! = (layoutRot (Items.type w.base.items cur)
              (Items.nvList g w.base.items cur).length s₇.base.nodeEdges.size
              (s₇.base.nodeEdges.size + Items.nEdges g w.base.items cur)
              (List.map (fun x => x - (1 + g.nv)) (List.filter (fun x => decide (x ≥ 1 + g.nv)) children))
              (fun ve => if (!planar) = true then Array.replicate 4 none else
                Array.map (fun z => Option.map (fun o => 4 * s₈.aux.rotEdgeNe[QE.edge o]! + (o &&& 2) + (1 - z % 2))
                  s₈.aux.qem[4 * ve + z]!) (Array.range 4)) (2 * g.ne))[l]! := by
            rw [e11rot, show 4 * s₇.base.nodeEdges.size = s.aux.neRotAdj.size by rw [h.rot_size, eneSt],
              Array.getElem!_append_right' _ _ _ (by rw [hRsz, hnER]; simpa using hl)]
          rw [hrot, hR, mapRot_R_get _ _ _ _ _ _ _ planar hpl (by omega) l (by simpa using hl), e8q]
          have hmem : (capVe g :: (children.filter (· ≥ 1 + g.nv)).map (· - (1 + g.nv)))[l / 4]! ∈
              capVe g :: (children.filter (· ≥ 1 + g.nv)).map (· - (1 + g.nv)) := by
            rw [List.getElem!_eq_getElem?_getD, List.getElem?_eq_getElem (by simp at hl ⊢; omega)]
            exact List.getElem_mem _
          have hS : (∃ c ∈ Items.ch w.base.items cur, 1 + g.nv ≤ c ∧
              4 * (c - (1 + g.nv)) ≤ 4 * (capVe g :: (children.filter (· ≥ 1 + g.nv)).map (· - (1 + g.nv)))[l / 4]! + l % 4 ∧
              4 * (capVe g :: (children.filter (· ≥ 1 + g.nv)).map (· - (1 + g.nv)))[l / 4]! + l % 4 < 4 * (c - (1 + g.nv)) + 4) ∨
              ∃ s, s < 4 ∧ 4 * (capVe g :: (children.filter (· ≥ 1 + g.nv)).map (· - (1 + g.nv)))[l / 4]! + l % 4 =
                4 * capVe g + s := by
            rcases List.mem_cons.1 hmem with hcap | hv
            · exact Or.inr ⟨l % 4, Nat.mod_lt _ (by omega), by rw [hcap]⟩
            · obtain ⟨c, hc, hge', _, he⟩ := hves_mem _ hv
              exact Or.inl ⟨c, hc, hge', by rw [he]; (try simp only [ItemId]); omega,
                by rw [he]; (try simp only [ItemId]); omega⟩
          rw [show (2 * g.ne : Nat) = capVe g from rfl, hagree.2 _ hS]
          refine ⟨fun h0 => by rw [h0]; rfl, fun o ho hmo => ?_⟩
          rw [ho, Option.map_some, hrne _ hmo]
          congr 1; omega
    try simp only [wp_bind]
    -- children
    refine wp_forIn_inv _ _ _ _
      (fun rest _ s' => RowInv g w s' ∧ Unch g w.base.items cur s s' ∧ (∀ c ∈ rest, Fresh g w c s') ∧
        rest.Sublist children) _ ⟨hrow, hunch, hfresh, List.Sublist.refl _⟩ ?_ ?_
    · intro c hc rest b s' ⟨hinv', hu', hfr', hsub⟩
      simp only [wp_bind, wp_liftR_get, wp_liftR_modify, wp_pure, wp_ite]
      split_ifs
      all_goals
        unfold wp
        try dsimp only
        refine ⟨_, rfl, rowInv_child_step g w hwf fuel ih hperm hsub hinv' hu' hfr' ?_ _ _ _⟩
        exact ⟨rfl, rfl, rfl, rfl, rfl, (by first | rfl | (unfold nvsArr; exact nvsArr_modify_twin _ _ _)), rfl, rfl,
          rfl, rfl, rfl, rfl, rfl, rfl⟩
    · intro b s' ⟨hinv', hu', _, _⟩
      try dsimp only
      exact wp_rfr _ (fun s => RFr.liftR_modify _ s rfl rfl rfl rfl rfl rfl rfl rfl)
        fun _ s'' hf => ⟨hinv'.fr hf, hu'.rfr hf⟩

end Spqr
