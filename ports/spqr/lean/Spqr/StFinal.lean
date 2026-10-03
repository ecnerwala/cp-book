import Spqr.StInduct
import Spqr.StRestrict
import Spqr.WalkItemsWF

/-!
# `walk_st'` / `walk_vsOriented` from the simulation (PROOF.md §7.6)

`walk_sim` runs `stForest` on the real forest: at the end of the walk every S / P / R item is
`InBlock` of a reference block. `walk_vsOriented` is then the orientation clause of `InBlock` once
`Items.leaves items items.size i` is identified with the `Expands` leaves (the item tree is acyclic
with unique parents, so every parent chain is shorter than `items.size`); `walk_st'` is the
restriction computation: the leaves of `i` form a segment of the `Nodup` reference order, each
child's leaves are a nonempty constant run of the restriction, and `collapseRuns` recovers `ch i`.
-/

namespace Spqr

/-! ## Acyclicity: no parent chain repeats an item -/

/-- With unique parents and a root above every item, no item is below one of its children. -/
theorem no_cycle {items : Items} (hlt : ∀ p c, Items.IsParent items p c → c < items.size)
    (huniq : ∀ x p c, x < items.size → c ∈ Items.ch items x → Items.IsParent items p c → p = x)
    (hacyc : Items.Acyc items) : ∀ y z, Items.IsParent items z y → ¬ Items.Below items y z := by
  intro y z hzy hyz
  obtain ⟨r, hr, hry⟩ := hacyc y (hlt z y hzy)
  revert z
  induction hry with
  | refl => exact fun z hz _ => hr z hz
  | tail hrp hpy ih =>
    rename_i p y
    intro z hzy hyz
    have hplt : p < items.size := by
      by_contra hp; simp [Items.IsParent, Items.ch_of_le _ _ (Nat.le_of_not_lt hp)] at hpy
    have hzp : z = p := huniq p z y hplt hpy hzy
    rw [hzp] at hzy hyz
    rcases Relation.ReflTransGen.cases_head hyz with heq | ⟨w, hyw, hwz⟩
    · exact ih p (heq ▸ hzy) Relation.ReflTransGen.refl
    · rcases Relation.ReflTransGen.cases_tail hwz with heq | ⟨q, hwq, hqz⟩
      · exact ih y (heq ▸ hyw) (Relation.ReflTransGen.single hzy)
      · exact ih q hqz (((Relation.ReflTransGen.single hzy).tail hyw).trans hwq)

/-- A parent chain from `x` is `Nodup`, hence shorter than `items.size`. -/
theorem chain_lt_size {items : Items} (hlt : ∀ p c, Items.IsParent items p c → c < items.size)
    (hnc : ∀ y z, Items.IsParent items z y → ¬ Items.Below items y z) {x : ItemId} (hx : x < items.size) :
    ∀ c : List ItemId, List.IsChain (Items.IsParent items) (x :: c) → c.length < items.size := by
  intro c hc
  have hpw : List.Pairwise (Relation.TransGen (Items.IsParent items)) (x :: c) :=
    List.isChain_iff_pairwise.1 (hc.imp fun _ _ h => Relation.TransGen.single h)
  have hnd : (x :: c).Nodup := hpw.imp fun {a b} h e => by
    subst e
    obtain ⟨w, haw, hwa⟩ := Relation.TransGen.tail'_iff.1 h
    exact hnc a w hwa haw
  have hsub : x :: c ⊆ List.range items.size := by
    intro y hy
    rw [List.mem_range]
    rcases List.mem_cons.1 hy with rfl | hy
    · exact hx
    · obtain ⟨w, -, hwy⟩ := Relation.TransGen.tail'_iff.1 (List.rel_of_pairwise_cons hpw hy)
      exact hlt w y hwy
  have := (List.subperm_of_subset hnd hsub).length_le
  simp at this; omega

/-! ## `Items.leaves` agrees with `Expands` -/

theorem ExpandsList.nil_inv {items : Items} {L : List ItemId} (h : ExpandsList items [] L) : L = [] := by
  cases h; rfl

theorem Expands.leaf_inv {items : Items} {x : ItemId} {L : List ItemId} (h : Expands items x L)
    (hx : Items.type items x = .V ∨ Items.type items x = .Q) : L = [x] := by
  cases h with
  | leaf _ h' => rw [ExpandsList.nil_inv h']
  | node h' _ => exact absurd hx h'

theorem Expands.node_inv {items : Items} {x : ItemId} {L : List ItemId} (h : Expands items x L)
    (hx : ¬ (Items.type items x = .V ∨ Items.type items x = .Q)) :
    ExpandsList items (Items.ch items x) L := by
  cases h with
  | leaf h' _ => exact absurd h' hx
  | node _ h' => simpa using h'

theorem leaves_eq_of_expands {items : Items} : ∀ (fuel : Nat) {xs L : List ItemId},
    ExpandsList items xs L →
    (∀ x ∈ xs, ∀ c, List.IsChain (Items.IsParent items) (x :: c) → c.length < fuel) →
    xs.flatMap (Items.leaves items fuel) = L := by
  intro fuel
  induction fuel with
  | zero =>
    intro xs L h hb
    cases h with
    | nil => rfl
    | @leaf x _ _ _ _ => exact absurd (hb x (by simp) [] (by simp)) (by omega)
    | @node x _ _ _ _ => exact absurd (hb x (by simp) [] (by simp)) (by omega)
  | succ fuel ih =>
    intro xs
    induction xs with
    | nil => intro L h _; rw [ExpandsList.nil_inv h]; rfl
    | cons x xs ihx =>
      intro L h hb
      obtain ⟨L₁, L₂, rfl, h₁, h₂⟩ := ExpandsList.cons_iff.1 h
      rw [List.flatMap_cons, ihx h₂ (fun y hy => hb y (by simp [hy]))]
      congr 1
      by_cases hx : Items.type items x = .V ∨ Items.type items x = .Q
      · rw [h₁.leaf_inv hx]
        rcases hx with hx | hx <;> simp [Items.leaves, hx]
      · have e₁ : Items.leaves items (fuel + 1) x =
            (Items.ch items x).flatMap (Items.leaves items fuel) := by
          rcases ht : Items.type items x <;> simp [Items.leaves, ht] <;>
            first | exact absurd (Or.inl ht) hx | exact absurd (Or.inr ht) hx
        rw [e₁]
        refine ih (h₁.node_inv hx) fun c hc cs hcs => ?_
        have := hb x (by simp) (c :: cs) (List.isChain_cons_cons.2 ⟨hc, hcs⟩)
        simp at this; omega

theorem leaves_size_eq {items : Items} (hlt : ∀ p c, Items.IsParent items p c → c < items.size)
    (hnc : ∀ y z, Items.IsParent items z y → ¬ Items.Below items y z) {x : ItemId} {L : List ItemId}
    (hx : x < items.size) (h : Expands items x L) : Items.leaves items items.size x = L := by
  have := leaves_eq_of_expands items.size (xs := [x]) h fun y hy c hc => by
    simp at hy; subst hy; exact chain_lt_size hlt hnc hx c hc
  simpa using this

/-! ## Nonempty expansions -/

/-- Reachability through non-leaf items only. -/
def ExpBelow (items : Items) : ItemId → ItemId → Prop :=
  Relation.ReflTransGen fun p c =>
    Items.IsParent items p c ∧ ¬ (Items.type items p = .V ∨ Items.type items p = .Q)

theorem expandsList_ne_nil {items : Items} : ∀ {xs L : List ItemId}, ExpandsList items xs L → xs ≠ [] →
    (∀ x ∈ xs, ∀ z, ExpBelow items x z → ¬ (Items.type items z = .V ∨ Items.type items z = .Q) →
      Items.ch items z ≠ []) →
    L ≠ [] := by
  intro xs L h
  induction h with
  | nil => exact fun h _ => absurd rfl h
  | leaf _ _ _ => exact fun _ _ => List.cons_ne_nil _ _
  | node hx hxs ih =>
    rename_i x xs L
    intro _ hb
    have hch : Items.ch items x ≠ [] := hb x (by simp) x Relation.ReflTransGen.refl hx
    refine ih (fun e => hch (List.append_eq_nil_iff.1 e).1) ?_
    intro y hy z hyz hz
    rcases List.mem_append.1 hy with hy | hy
    · exact hb x (by simp) z (Relation.ReflTransGen.head ⟨hy, hx⟩ hyz) hz
    · exact hb y (by simp [hy]) z hyz hz

/-- Admitted (walk shape, dump-checked by `check_stsim`): an `I` item hangs under a `Q` item. -/
theorem walk_i_parent (g : Graph) (tern : Bool) (vo eo : List Nat) :
    ∀ p c, Items.IsParent (g.walk tern (g.dfsForest vo eo)).items p c →
      Items.type (g.walk tern (g.dfsForest vo eo)).items c = .I →
      Items.type (g.walk tern (g.dfsForest vo eo)).items p = .Q := by
  sorry

theorem spr_of_expBelow {g : Graph} {items : Items} (hwf : Items.WF g items)
    (hI : ∀ p c, Items.IsParent items p c → Items.type items c = .I → Items.type items p = .Q)
    {i : ItemId} (ht : Items.type items i = .S ∨ Items.type items i = .P ∨ Items.type items i = .R) :
    ∀ z, ExpBelow items i z → ¬ (Items.type items z = .V ∨ Items.type items z = .Q) →
      Items.type items z = .S ∨ Items.type items z = .P ∨ Items.type items z = .R := by
  intro z hz
  unfold ExpBelow at hz
  induction hz with
  | refl => exact fun _ => ht
  | @tail p z _ hpz ih =>
    intro hz
    obtain ⟨hpz, hpl⟩ := hpz
    have hpS := ih hpl
    have hzlt : z < items.size := hwf.tree.ch_lt p z hpz
    have hF : Items.type items z ≠ .F := by
      intro hF
      cases z with
      | zero => exact hwf.tree.root_no_parent p hpz
      | succ m =>
        by_cases h1 : m < g.nv
        · have := hwf.tree.vert m h1
          rw [vertItem, Nat.add_comm] at this
          simp [this] at hF
        by_cases h2 : m + 1 < 1 + g.nv + g.ne
        · have := hwf.tree.edge (m + 1 - (1 + g.nv)) (by omega)
          rw [edgeItem, Nat.add_sub_cancel' (by omega)] at this
          simp [this] at hF
        · have := hwf.tree.node (m + 1) (by omega) hzlt
          simp [hF] at this
    cases htz : Items.type items z
    · exact absurd htz hF
    · exact absurd (Or.inl htz) hz
    · exact absurd (Or.inr htz) hz
    · have := hI p z hpz htz; simp [this] at hpS
    · have := (hwf.shapes.o_parent p z hpz htz).1; simp [this] at hpS
    · simp
    · simp
    · simp

theorem spr_ch_ne_nil {g : Graph} {items : Items} (hwf : Items.WF g items) {i : ItemId}
    (hi : i < items.size)
    (ht : Items.type items i = .S ∨ Items.type items i = .P ∨ Items.type items i = .R) :
    Items.ch items i ≠ [] := by
  intro hnil
  rcases ht with ht | ht | ht
  · obtain ⟨u, v, xs, -, hf, hlen, -⟩ := hwf.shapes.s_shape i hi ht
    rw [hnil] at hf
    simp at hf
    simp [hf] at hlen
  · have := (hwf.shapes.p_shape i hi ht).1
    simp [Items.virtualEdges, hnil] at this
  · have := (hwf.shapes.r_shape i hi ht).1
    simp [hnil] at this

theorem expands_ne_nil {g : Graph} {items : Items} (hwf : Items.WF g items)
    (hI : ∀ p c, Items.IsParent items p c → Items.type items c = .I → Items.type items p = .Q)
    {i : ItemId} (hi : i < items.size)
    (ht : Items.type items i = .S ∨ Items.type items i = .P ∨ Items.type items i = .R)
    {L : List ItemId} (h : Expands items i L) : L ≠ [] := by
  refine expandsList_ne_nil h (by simp) fun x hx z hz hzl => ?_
  simp at hx; subst hx
  have hzS := spr_of_expBelow hwf hI ht z hz hzl
  have hzlt : z < items.size := by
    unfold ExpBelow at hz
    rcases Relation.ReflTransGen.cases_tail hz with rfl | ⟨p, -, hpz, -⟩
    · exact hi
    · exact hwf.tree.ch_lt p z hpz
  exact spr_ch_ne_nil hwf hzlt hzS

/-! ## The restriction recovers `ch i` -/

theorem collapseRuns_cons_ne {x : ItemId} {R : List ItemId} (h : ∀ y ∈ R.head?, y ≠ x) :
    collapseRuns (x :: R) = x :: collapseRuns R := by
  cases R with
  | nil => rfl
  | cons y R => simp [collapseRuns, Ne.symm (h y rfl)]

theorem collapseRuns_replicate_append {x : ItemId} {R : List ItemId} :
    ∀ n, 0 < n → (∀ y ∈ R.head?, y ≠ x) →
      collapseRuns (List.replicate n x ++ R) = x :: collapseRuns R := by
  intro n hn h
  induction n with
  | zero => omega
  | succ n ih =>
    rcases Nat.eq_zero_or_pos n with rfl | hn'
    · simpa using collapseRuns_cons_ne h
    · rw [← ih hn']
      obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by omega⟩
      simp [List.replicate_succ, collapseRuns]

theorem ExpandsList.expands_of_mem {items : Items} : ∀ {cs L : List ItemId}, ExpandsList items cs L →
    ∀ c ∈ cs, ∃ Lc, Expands items c Lc ∧ Lc ⊆ L := by
  intro cs
  induction cs with
  | nil => intro L _ c hc; simp at hc
  | cons c' cs ih =>
    intro L h c hc
    obtain ⟨L₁, L₂, rfl, h₁, h₂⟩ := ExpandsList.cons_iff.1 h
    rcases List.mem_cons.1 hc with rfl | hc
    · exact ⟨L₁, h₁, List.subset_append_left _ _⟩
    · obtain ⟨Lc, hLc, hsub⟩ := ih h₂ c hc
      exact ⟨Lc, hLc, fun {_} ha => List.mem_append_right _ (hsub ha)⟩

theorem restrict_mid {items : Items} (fuel : Nat) : ∀ {cs L : List ItemId}, ExpandsList items cs L →
    (∀ c ∈ cs, ∀ Lc, Expands items c Lc → Items.leaves items fuel c = Lc ∧ Lc ≠ []) →
    L.Nodup → cs.Nodup →
    collapseRuns (L.filterMap fun x => cs.find? fun c => x ∈ Items.leaves items fuel c) = cs := by
  intro cs
  induction cs with
  | nil => intro L h _ _ _; rw [ExpandsList.nil_inv h]; rfl
  | cons c cs ih =>
    intro L h hl hnd hcs
    obtain ⟨L₁, L₂, rfl, h₁, h₂⟩ := ExpandsList.cons_iff.1 h
    obtain ⟨hleq, hne⟩ := hl c (by simp) L₁ h₁
    rw [List.filterMap_append]
    have e₁ : (L₁.filterMap fun x => (c :: cs).find? fun c => x ∈ Items.leaves items fuel c) =
        List.replicate L₁.length c := by
      rw [← List.map_const']
      rw [List.filterMap_eq_map_iff_forall_eq_some]
      intro x hx
      exact List.find?_cons_of_pos (by simp [hleq, hx])
    have e₂ : (L₂.filterMap fun x => (c :: cs).find? fun c => x ∈ Items.leaves items fuel c) =
        L₂.filterMap fun x => cs.find? fun c => x ∈ Items.leaves items fuel c := by
      refine List.filterMap_congr fun x hx => ?_
      refine List.find?_cons_of_neg ?_
      have : x ∉ L₁ := fun hx' => (List.nodup_append.1 hnd).2.2 x hx' x hx rfl
      simp [hleq, this]
    rw [e₁, e₂, collapseRuns_replicate_append _ (List.length_pos_iff_ne_nil.2 hne) ?_,
      ih h₂ (fun c' hc' => hl c' (by simp [hc'])) (List.nodup_append.1 hnd).2.1 (List.nodup_cons.1 hcs).2]
    intro y hy
    have hy' := List.mem_of_mem_head? hy
    obtain ⟨x, -, hx⟩ := List.mem_filterMap.1 hy'
    have := List.mem_of_find?_eq_some hx
    exact fun e => (List.nodup_cons.1 hcs).1 (e ▸ this)

theorem restrictCh_eq {items : Items} {fuel i : Nat} {A L B : List ItemId}
    (hL : ExpandsList items (Items.ch items i) L)
    (hl : ∀ c ∈ Items.ch items i, ∀ Lc, Expands items c Lc →
      Items.leaves items fuel c = Lc ∧ Lc ≠ [])
    (hnd : (A ++ L ++ B).Nodup) (hcs : (Items.ch items i).Nodup) :
    restrictCh items fuel (A ++ L ++ B) i = Items.ch items i := by
  have hnone : ∀ x, x ∉ L →
      ((Items.ch items i).find? fun c => x ∈ Items.leaves items fuel c) = none := by
    intro x hx
    rw [List.find?_eq_none]
    intro c hc
    obtain ⟨Lc, hLc, hsub⟩ := hL.expands_of_mem c hc
    rw [(hl c hc Lc hLc).1]
    simpa using fun h => hx (hsub h)
  have hA : (A.filterMap fun x => (Items.ch items i).find? fun c => x ∈ Items.leaves items fuel c) = [] :=
    List.filterMap_eq_nil_iff.2 fun x hx => hnone x fun hx' =>
      (List.nodup_append.1 (List.nodup_append.1 hnd).1).2.2 x hx x hx' rfl
  have hB : (B.filterMap fun x => (Items.ch items i).find? fun c => x ∈ Items.leaves items fuel c) = [] :=
    List.filterMap_eq_nil_iff.2 fun x hx => hnone x fun hx' =>
      (List.nodup_append.1 hnd).2.2 x (List.mem_append_right _ hx') x hx rfl
  unfold restrictCh
  rw [List.filterMap_append, List.filterMap_append, hA, hB, List.nil_append, List.append_nil]
  exact restrict_mid fuel hL hl (List.nodup_append.1 (List.nodup_append.1 hnd).1).2.1 hcs

/-! ## The final state of the walk -/

theorem Items.WF.no_cycle {g : Graph} {items : Items} (hwf : Items.WF g items) :
    ∀ y z, Items.IsParent items z y → ¬ Items.Below items y z :=
  _root_.Spqr.no_cycle hwf.tree.ch_lt
    (fun x p c _ hc hp => by
      have hc0 : 0 < c := by
        rcases Nat.eq_zero_or_pos c with rfl | h
        · exact absurd hc (hwf.tree.root_no_parent x)
        · exact h
      obtain ⟨p₀, -, huniq⟩ := hwf.tree.unique_parent c hc0 (hwf.tree.ch_lt x c hc)
      rw [huniq p hp, huniq x hc])
    fun i hi => ⟨rootItem, hwf.tree.root_no_parent, hwf.tree.reach i hi⟩

theorem stItems_init (g : Graph) (tern : Bool) :
    StItems g (WalkState.init g tern) (refBlocks g []) where
  roots := by simp [WalkState.init, readStack, readL, readR]
  nodup := by simp [WalkState.init, readStack, readL, readR]
  bounded := by simp [WalkState.init, readStack, readL, readR]
  chLt p c h := by simp [Items.IsParent, WalkState.init, Items.initialItems_ch] at h
  chNodup p := by simp [WalkState.init, Items.initialItems_ch]
  closed i hi ht := by
    simp only [WalkState.init] at hi ht
    rw [Items.initialItems_type] at ht
    rw [Items.initialItems_size] at hi
    split_ifs at ht <;> simp at ht

theorem walk_sim (g : Graph) (tern : Bool) (forest : List DfsTree) (hf : ForestOK g forest)
    (hwf : ∀ t ∈ forest, t.WF []) (hends : ∀ t ∈ forest, t.Ends g) :
    StItems g (g.walk tern forest) (refBlocks g forest) := by
  have hht : ∀ t ∈ forest, t.height ≤ g.nv := by
    intro t ht
    have hnd : t.verts.Nodup := by
      obtain ⟨pre, post, rfl⟩ := List.append_of_mem ht
      have h := hf.verts_nodup
      rw [List.flatMap_append, List.flatMap_cons, List.nodup_append, List.nodup_append] at h
      exact h.2.1.1
    have hsub : t.verts ⊆ List.range g.nv := fun v hv =>
      List.mem_range.2 (hf.verts_lt v (List.mem_flatMap.2 ⟨t, ht, hv⟩))
    have := (List.subperm_of_subset hnd hsub).length_le
    simp at this
    exact (DfsTree.height_le_verts_length t).trans this
  have h := stForest g forest [] (WalkState.init g tern) (fun _ => False) (fun _ => False)
    (WalkState.rootState_init g tern) (WalkState.init_full g tern) (fun _ _ h => h) (fun _ _ h => h)
    (by simpa using hf) hwf hends hht (stItems_init g tern)
  simpa [WalkM.wp, Graph.walk] using h

theorem walk_inBlock (g : Graph) (hg : g.WF) (tern : Bool) (vo eo : List Nat)
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) (i : ItemId)
    (hi : i < (g.walk tern (g.dfsForest vo eo)).items.size)
    (ht : Items.type (g.walk tern (g.dfsForest vo eo)).items i = .S ∨
      Items.type (g.walk tern (g.dfsForest vo eo)).items i = .P ∨
      Items.type (g.walk tern (g.dfsForest vo eo)).items i = .R) :
    ∃ b ∈ refBlocks g (g.dfsForest vo eo), InBlock g (g.walk tern (g.dfsForest vo eo)).items b i := by
  obtain ⟨hvp, hep⟩ := dfsForest_spanning' hg hvo heo
  have hf : ForestOK g (g.dfsForest vo eo) := ForestOK.of_perm hvp hep
  have hS := walk_sim g tern _ hf (dfsForest_wf hg hvo heo) (dfsForest_ends g hg hvo heo)
  have hts := (walk_full g tern _ hf (dfsForest_wf hg hvo heo) (dfsForest_ends g hg hvo heo)).2
  rcases hS.closed i hi ht with ⟨x, hx, -⟩ | h
  · rw [hts] at hx; simp [readStack, readL, readR] at hx
  · exact h

theorem walk_vsOriented (g : Graph) (hg : g.WF) (tern : Bool) (vo eo : List Nat)
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) :
    VsOriented g (g.walk tern (g.dfsForest vo eo)).items (refBlocks g (g.dfsForest vo eo)) := by
  have hwf := walk_items_wf g hg tern vo eo hvo heo
  rw [vsOriented_iff]
  intro i hi ht
  obtain ⟨b, hb, L, hL, -, hV⟩ := walk_inBlock g hg tern vo eo hvo heo i hi ht
  rw [leaves_size_eq hwf.tree.ch_lt hwf.no_cycle hi hL]
  exact ⟨b, hb, hV⟩

/-- Admitted (reference side, dump-checked by `check_stref`): the reference order lists every
V / Q item at most once. -/
theorem refOrder_nodup (g : Graph) (hg : g.WF) (vo eo : List Nat) (hvo : OrderOK g.nv vo)
    (heo : OrderOK g.ne eo) : (refOrder g (g.dfsForest vo eo)).Nodup := by
  sorry

theorem walk_st' (g : Graph) (hg : g.WF) (tern : Bool) (vo eo : List Nat)
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) (i : ItemId)
    (hi : i < (g.walk tern (g.dfsForest vo eo)).items.size)
    (ht : Items.type (g.walk tern (g.dfsForest vo eo)).items i = .S ∨
      Items.type (g.walk tern (g.dfsForest vo eo)).items i = .P ∨
      Items.type (g.walk tern (g.dfsForest vo eo)).items i = .R) :
    Items.ch (g.walk tern (g.dfsForest vo eo)).items i =
      restrictCh (g.walk tern (g.dfsForest vo eo)).items
        (g.walk tern (g.dfsForest vo eo)).items.size
        (refOrder g (g.dfsForest vo eo)) i := by
  have hwf := walk_items_wf g hg tern vo eo hvo heo
  have hI := walk_i_parent g tern vo eo
  obtain ⟨b, hb, L, hL, ⟨A, B, hAB⟩, -⟩ := walk_inBlock g hg tern vo eo hvo heo i hi ht
  obtain ⟨pre, post, hbl⟩ := List.append_of_mem hb
  have horder : refOrder g (g.dfsForest vo eo) =
      (pre.flatMap StBlock.items ++ A) ++ L ++ (B ++ post.flatMap StBlock.items) := by
    simp [refOrder, hbl, List.flatMap_append, List.flatMap_cons, hAB]
  have hnd := refOrder_nodup g hg vo eo hvo heo
  rw [horder] at hnd ⊢
  have hxnl : ¬ (Items.type (g.walk tern (g.dfsForest vo eo)).items i = .V ∨
      Items.type (g.walk tern (g.dfsForest vo eo)).items i = .Q) := by
    rcases ht with h | h | h <;> simp [h]
  refine (restrictCh_eq (hL.node_inv hxnl) ?_ hnd (hwf.tree.ch_nodup i)).symm
  intro c hc Lc hLc
  have hclt := hwf.tree.ch_lt i c hc
  refine ⟨leaves_size_eq hwf.tree.ch_lt hwf.no_cycle hclt hLc, ?_⟩
  by_cases hcl : Items.type (g.walk tern (g.dfsForest vo eo)).items c = .V ∨
      Items.type (g.walk tern (g.dfsForest vo eo)).items c = .Q
  · rw [hLc.leaf_inv hcl]; simp
  · exact expands_ne_nil hwf hI hclt
      (spr_of_expBelow hwf hI ht c (Relation.ReflTransGen.single ⟨hc, hxnl⟩) hcl) hLc

/-- Phase 2: the walk's children lists are in s-t order (`PROOF.md` §7.6): they are the reference
order (`walk_st'`) with oriented `vs` (`walk_vsOriented`), which is an st-order
(`stItem_of_refOrder`, on the st-numbered reference blocks `refBlocks_st`, which needs `g.WF`
and the order hypotheses like `dfsForest_spanning`). -/
theorem walk_st (g : Graph) (hg : g.WF) (tern : Bool) (vo eo : List Nat) (hvo : OrderOK g.nv vo)
    (heo : OrderOK g.ne eo) (hwf : Items.WF g (g.walk tern (g.dfsForest vo eo)).items) :
    Items.StNumbered (g.walk tern (g.dfsForest vo eo)).items :=
  fun i hi ht => stItem_of_refOrder g tern vo eo i hi ht hwf
    (walk_st' g hg tern vo eo hvo heo i hi ht) (walk_vsOriented g hg tern vo eo hvo heo)
    (refBlocks_st hg hvo heo)
    (refBlocks_root_edge hg hvo heo)

theorem spqrTree_st (g : Graph) (hg : g.WF) (tern : Bool) (vo eo : List Nat) (hvo : OrderOK g.nv vo)
    (heo : OrderOK g.ne eo) (hwf : Items.WF g (g.walk tern (g.dfsForest vo eo)).items) : (g.spqrTree tern vo eo).StOrder := by
  rw [spqrTree_eq]; exact relabel_st g _ (walk_st g hg tern vo eo hvo heo hwf) hwf

end Spqr
