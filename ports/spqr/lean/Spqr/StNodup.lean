import Spqr.StSimLemmas
import Spqr.WalkPlace
import Spqr.WalkInv

/-! # The reference order repeats no item (`PROOF.md` §7)

`refTree`/`refOrder` list each vertex and edge item at most once: the pieces and blocks of a
subtree contain only the vertex items of its vertices and the edge items of its edges, each once
(`refTree_items`), and the DFS forest visits every vertex and edge once (`ForestOK`). -/

namespace Spqr
open StRefEt

def pieceItems (ps : List StPiece) : List ItemId := ps.flatMap StPiece.items

theorem pieceItems_append (ps qs : List StPiece) :
    pieceItems (ps ++ qs) = pieceItems ps ++ pieceItems qs := List.flatMap_append ..

theorem stNest_perm : ∀ ps : List StPiece, (stNest ps).Perm (pieceItems ps)
  | [] => by simp [stNest, stNestL, stNestR, pieceItems]
  | p :: ps => by
    refine List.perm_iff_count.2 fun a => ?_
    have ih := (stNest_perm ps).count_eq a
    simp only [stNest, List.count_append, pieceItems, List.flatMap_cons] at ih ⊢
    simp only [stNestL, stNestR, List.count_append]
    by_cases hs : p.side = true <;> simp [hs] <;> omega

theorem edgeItem_inj' {g : Graph} {e e' : Nat} (h : edgeItem g e = edgeItem g e') : e = e' := by
  have : (1 + g.nv + e : Nat) = 1 + g.nv + e' := h; omega

/-- Provenance of an item of the pieces/blocks of out-edges of `v`: the vertex item of a vertex
in `vs`, the edge item of an edge in `es`, or `vertItem v` when it was first pushed here
(`hv = false`, `hv' = true`). -/
def Prov (g : Graph) (v : Nat) (vs es : List Nat) (hv hv' : Bool) (i : ItemId) : Prop :=
  (∃ x ∈ vs, i = vertItem x) ∨ (∃ e ∈ es, i = edgeItem g e) ∨
    (hv = false ∧ hv' = true ∧ i = vertItem v)

theorem Prov.ne {g : Graph} {v : Nat} {vs₁ es₁ vs₂ es₂ : List Nat} {hv hv₁ hv₂ : Bool} {i : ItemId}
    (hvd : ∀ x ∈ vs₁, x ∉ vs₂) (hed : ∀ e ∈ es₁, e ∉ es₂) (hv1 : v ∉ vs₁) (hv2 : v ∉ vs₂)
    (hlt₁ : ∀ x ∈ vs₁, x < g.nv) (hlt₂ : ∀ x ∈ vs₂, x < g.nv) (hvlt : v < g.nv)
    (h₁ : Prov g v vs₁ es₁ hv hv₁ i) (h₂ : Prov g v vs₂ es₂ hv₁ hv₂ i) : False := by
  rcases h₁ with ⟨x, hx, rfl⟩ | ⟨e, he, rfl⟩ | ⟨-, h1, rfl⟩ <;>
    rcases h₂ with ⟨y, hy, h⟩ | ⟨e', he', h⟩ | ⟨h2, -, h⟩
  · exact hvd x hx (by rw [vertItem_inj h]; exact hy)
  · exact vertItem_ne_edgeItem (hlt₁ x hx) h
  · exact hv1 (by rw [← vertItem_inj h]; exact hx)
  · exact vertItem_ne_edgeItem (hlt₂ y hy) h.symm
  · exact hed e he (by rw [edgeItem_inj' h]; exact he')
  · exact vertItem_ne_edgeItem hvlt h.symm
  · exact hv2 (by rw [vertItem_inj h]; exact hy)
  · exact vertItem_ne_edgeItem hvlt h
  · rw [h1] at h2; cases h2

def outsItems (r : List StPiece × List StBlock × Bool) : List ItemId :=
  r.2.1.flatMap StBlock.items ++ pieceItems r.1

def treeItems (r : List StPiece × List StBlock) : List ItemId :=
  r.2.flatMap StBlock.items ++ pieceItems r.1

/-- The items of the pieces of a returning out-edge other than the vertex item of `v`. -/
def outMid (g : Graph) (d : Nat) (dirs : List Bool) : DfsOut → List ItemId
  | .tree e cls child =>
    pieceItems (refTree g child (d + 1) (dirs ++ [!dirs.getD (cls.lowval d) false])).1 ++ [edgeItem g e]
  | .back e _ _ => [edgeItem g e]

def outBlocks (g : Graph) (d : Nat) (dirs : List Bool) : DfsOut → List StBlock
  | .tree _ cls child => (refTree g child (d + 1) (dirs ++ [!dirs.getD (cls.lowval d) false])).2
  | .back .. => []

theorem refOut_ret_items {g : Graph} {v d : Nat} {dirs : List Bool} {o : DfsOut} {hv : Bool}
    (hlt : o.cls.lowval d < d) :
    (refOut g v d dirs o hv).2.2 = true ∧
    (pieceItems (refOut g v d dirs o hv).1).Perm
      ((if hv then [] else [vertItem v]) ++ outMid g d dirs o) ∧
    (refOut g v d dirs o hv).2.1 = outBlocks g d dirs o := by
  obtain ⟨lv, kind, ho, hlv⟩ := WalkState.ret_of_lowval_lt hlt
  cases o with
  | back e dest cls =>
    simp only [DfsOut.cls] at ho; subst ho
    cases hv <;> cases kind <;>
      simp [refOut, DfsOut.cls, OutClass.lowval, OutClass.isType1, Nat.not_le.mpr hlv, pieceItems,
        outMid, outBlocks]
    exact List.Perm.swap _ _ _
  | tree e cls child =>
    simp only [DfsOut.cls] at ho; subst ho
    cases hv <;> cases kind <;>
      simp [refOut, DfsOut.cls, OutClass.lowval, OutClass.isType1, Nat.not_le.mpr hlv, pieceItems,
        outMid, outBlocks] <;>
      first
      | exact (stNest_perm _).trans (by simp [pieceItems])
      | exact List.Perm.cons _ ((stNest_perm _).trans (by simp [pieceItems]))
      | exact List.perm_append_comm
      | exact List.Perm.swap _ _ _
      | exact ((List.Perm.swap _ _ _).append_left _).trans List.perm_middle

theorem refOut_items_of {g : Graph} {v d : Nat} {dirs : List Bool} {o : DfsOut} {hv : Bool}
    (hvn : (v :: DfsOut.vertsList [o]).Nodup) (hlt : ∀ x ∈ v :: DfsOut.vertsList [o], x < g.nv)
    (hen : (DfsOut.edgesList [o]).Nodup)
    (ih : ∀ e cls child, o = .tree e cls child → ∀ dirs' : List Bool,
      (treeItems (refTree g child (d + 1) dirs')).Nodup ∧
      ∀ i ∈ treeItems (refTree g child (d + 1) dirs'),
        (∃ x ∈ child.verts, i = vertItem x) ∨ ∃ e ∈ child.edges, i = edgeItem g e) :
    (outsItems (refOut g v d dirs o hv)).Nodup ∧
    (hv = true → (refOut g v d dirs o hv).2.2 = true) ∧
    ∀ i ∈ outsItems (refOut g v d dirs o hv),
      Prov g v (DfsOut.vertsList [o]) (DfsOut.edgesList [o]) hv (refOut g v d dirs o hv).2.2 i := by
  have hvlt : v < g.nv := hlt v (List.mem_cons_self ..)
  by_cases hb : d ≤ o.cls.lowval d
  · cases o with
    | back e dest cls =>
      simp only [DfsOut.cls] at hb
      rw [refOut_boundary_back hb]; simp [outsItems, pieceItems]
    | tree e cls child =>
      simp only [DfsOut.cls] at hb
      rw [refOut_boundary_tree hb]
      obtain ⟨hnd, hprov⟩ := ih e cls child rfl (dirs ++ [false])
      have hperm : (outsItems ([], (refTree g child (d + 1) (dirs ++ [false])).2 ++
          [⟨some (v, child.v), stNest (refTree g child (d + 1) (dirs ++ [false])).1⟩], hv)).Perm
          (treeItems (refTree g child (d + 1) (dirs ++ [false]))) := by
        simp only [outsItems, treeItems, pieceItems, List.flatMap_append, List.flatMap_cons,
          List.flatMap_nil, List.append_nil]
        exact (stNest_perm _).append_left _
      refine ⟨hperm.nodup_iff.2 hnd, fun h => h, fun i hi => ?_⟩
      rcases hprov i (hperm.subset hi) with ⟨x, hx, h⟩ | ⟨e', he', h⟩
      · exact Or.inl ⟨x, by simp [DfsOut.vertsList, hx], h⟩
      · exact Or.inr (Or.inl ⟨e', by simp [DfsOut.edgesList, he'], h⟩)
  · have hlt' := Nat.lt_of_not_le hb
    obtain ⟨h22, hperm, hbl⟩ := refOut_ret_items (v := v) (dirs := dirs) (hv := hv) hlt'
    have hitems : (outsItems (refOut g v d dirs o hv)).Perm
        ((if hv then [] else [vertItem v]) ++ (outBlocks g d dirs o).flatMap StBlock.items ++
          outMid g d dirs o) := by
      simp only [outsItems]; rw [hbl]
      refine (hperm.append_left _).trans ?_
      rw [← List.append_assoc]; exact List.perm_append_comm.append_right _
    rw [h22]
    cases o with
    | back e dest cls =>
      simp only [outBlocks, outMid, List.flatMap_nil, List.append_nil] at hitems
      refine ⟨hitems.nodup_iff.2 ?_, fun _ => rfl, fun i hi => ?_⟩
      · cases hv <;> simp [vertItem_ne_edgeItem hvlt]
      · have := hitems.subset hi
        cases hv <;> simp at this
        · rcases this with rfl | rfl
          · exact Or.inr (Or.inr (by simp))
          · exact Or.inr (Or.inl ⟨e, by simp [DfsOut.edgesList], rfl⟩)
        · subst this; exact Or.inr (Or.inl ⟨e, by simp [DfsOut.edgesList], rfl⟩)
    | tree e cls child =>
      simp only [outBlocks, outMid] at hitems
      obtain ⟨hnd, hprov⟩ := ih e cls child rfl (dirs ++ [!dirs.getD (cls.lowval d) false])
      have hitems' : (outsItems (refOut g v d dirs (.tree e cls child) hv)).Perm
          ((if hv then [] else [vertItem v]) ++
            (treeItems (refTree g child (d + 1) (dirs ++ [!dirs.getD (cls.lowval d) false])) ++
              [edgeItem g e])) := by
        refine hitems.trans ?_
        simp only [treeItems, List.append_assoc]
        exact List.Perm.refl _
      simp only [DfsOut.vertsList, DfsOut.edgesList, List.append_nil] at hvn hlt hen ⊢
      have hcv : child.verts.Nodup := (List.nodup_cons.1 hvn).2
      have hvc : v ∉ child.verts := (List.nodup_cons.1 hvn).1
      have hltc : ∀ x ∈ child.verts, x < g.nv := fun x hx => hlt x (List.mem_cons_of_mem _ hx)
      have hec : e ∉ child.edges := (List.nodup_cons.1 hen).1
      have hmid : (treeItems (refTree g child (d + 1) (dirs ++ [!dirs.getD (cls.lowval d) false])) ++
          [edgeItem g e]).Nodup := by
        refine List.Nodup.append hnd (List.nodup_singleton _) fun i hi hi' => ?_
        simp at hi'; subst hi'
        rcases hprov _ hi with ⟨x, hx, h⟩ | ⟨e', he', h⟩
        · exact vertItem_ne_edgeItem (hltc x hx) h.symm
        · exact hec (by rw [edgeItem_inj' h]; exact he')
      refine ⟨hitems'.nodup_iff.2 ?_, by simp, fun i hi => ?_⟩
      · cases hv
        · simp only [Bool.false_eq_true, ↓reduceIte, List.singleton_append]
          refine List.nodup_cons.2 ⟨fun hi => ?_, hmid⟩
          rcases List.mem_append.1 hi with hi | hi
          · rcases hprov _ hi with ⟨x, hx, h⟩ | ⟨e', he', h⟩
            · exact hvc (by rw [vertItem_inj h]; exact hx)
            · exact vertItem_ne_edgeItem hvlt h
          · simp at hi; exact vertItem_ne_edgeItem hvlt hi
        · simpa using hmid
      · have hi' := hitems'.subset hi
        rcases List.mem_append.1 hi' with hi' | hi'
        · cases hv <;> simp at hi'
          subst hi'; exact Or.inr (Or.inr (by simp))
        rcases List.mem_append.1 hi' with hi' | hi'
        · rcases hprov _ hi' with ⟨x, hx, h⟩ | ⟨e', he', h⟩
          · exact Or.inl ⟨x, hx, h⟩
          · exact Or.inr (Or.inl ⟨e', List.mem_cons_of_mem _ he', h⟩)
        · simp at hi'; subst hi'; exact Or.inr (Or.inl ⟨e, List.mem_cons_self .., rfl⟩)

mutual
theorem refTree_items {g : Graph} : ∀ (t : DfsTree) (d : Nat) (dirs : List Bool),
    t.verts.Nodup → t.edges.Nodup → (∀ x ∈ t.verts, x < g.nv) →
    (treeItems (refTree g t d dirs)).Nodup ∧
    ∀ i ∈ treeItems (refTree g t d dirs),
      (∃ x ∈ t.verts, i = vertItem x) ∨ ∃ e ∈ t.edges, i = edgeItem g e
  | .node v outs, d, dirs, hvn, hen, hlt => by
    simp only [DfsTree.verts, DfsTree.edges] at hvn hen hlt ⊢
    obtain ⟨hnd, -, hprov⟩ := refOuts_items outs v d dirs false hvn hen hlt
    simp only [outsItems] at hnd hprov
    simp only [refTree_node, treeItems]
    have lift : ∀ i ∈ (refOuts g v d dirs outs false).2.1.flatMap StBlock.items ++
        pieceItems (refOuts g v d dirs outs false).1,
        (∃ x ∈ v :: DfsOut.vertsList outs, i = vertItem x) ∨
          ∃ e ∈ DfsOut.edgesList outs, i = edgeItem g e := fun i hi => by
      rcases hprov i hi with ⟨x, hx, h⟩ | ⟨e, he, h⟩ | ⟨-, -, h⟩
      · exact Or.inl ⟨x, List.mem_cons_of_mem _ hx, h⟩
      · exact Or.inr ⟨e, he, h⟩
      · exact Or.inl ⟨v, List.mem_cons_self .., h⟩
    split
    · exact ⟨hnd, lift⟩
    · rename_i hv'
      have hv'' : (refOuts g v d dirs outs false).2.2 = false := Bool.eq_false_iff.2 hv'
      rw [pieceItems_append, ← List.append_assoc]
      refine ⟨List.Nodup.append hnd (by simp [pieceItems]) fun i hi hi' => ?_, fun i hi => ?_⟩
      · simp [pieceItems] at hi'; subst hi'
        rcases hprov _ hi with ⟨x, hx, h⟩ | ⟨e, he, h⟩ | ⟨-, h, -⟩
        · exact (List.nodup_cons.1 hvn).1 (by rw [vertItem_inj h]; exact hx)
        · exact vertItem_ne_edgeItem (hlt v (List.mem_cons_self ..)) h
        · rw [hv''] at h; cases h
      · rcases List.mem_append.1 hi with hi | hi
        · exact lift i hi
        · simp [pieceItems] at hi; subst hi; exact Or.inl ⟨v, List.mem_cons_self .., rfl⟩

theorem refOuts_items {g : Graph} : ∀ (outs : List DfsOut) (v d : Nat) (dirs : List Bool) (hv : Bool),
    (v :: DfsOut.vertsList outs).Nodup → (DfsOut.edgesList outs).Nodup →
    (∀ x ∈ v :: DfsOut.vertsList outs, x < g.nv) →
    (outsItems (refOuts g v d dirs outs hv)).Nodup ∧
    (hv = true → (refOuts g v d dirs outs hv).2.2 = true) ∧
    ∀ i ∈ outsItems (refOuts g v d dirs outs hv),
      Prov g v (DfsOut.vertsList outs) (DfsOut.edgesList outs) hv (refOuts g v d dirs outs hv).2.2 i
  | [], v, d, dirs, hv, _, _, _ => by
    simp [refOuts_nil, outsItems, pieceItems]
  | o :: rest, v, d, dirs, hv, hvn, hen, hlt => by
    rw [DfsOut.vertsList_cons] at hvn hlt ⊢
    rw [DfsOut.edgesList_cons] at hen ⊢
    have hvlt : v < g.nv := hlt v (List.mem_cons_self ..)
    have hvn' := (List.nodup_cons.1 hvn).2
    have hvnot : v ∉ DfsOut.vertsList [o] ++ DfsOut.vertsList rest := (List.nodup_cons.1 hvn).1
    have hv1 : v ∉ DfsOut.vertsList [o] := fun h => hvnot (List.mem_append_left _ h)
    have hv2 : v ∉ DfsOut.vertsList rest := fun h => hvnot (List.mem_append_right _ h)
    have hvn₁ : (v :: DfsOut.vertsList [o]).Nodup :=
      List.nodup_cons.2 ⟨hv1, hvn'.sublist (List.sublist_append_left _ _)⟩
    have hvn₂ : (v :: DfsOut.vertsList rest).Nodup :=
      List.nodup_cons.2 ⟨hv2, hvn'.sublist (List.sublist_append_right _ _)⟩
    have hdisj : ∀ x ∈ DfsOut.vertsList [o], x ∉ DfsOut.vertsList rest :=
      fun x hx hx' => List.disjoint_of_nodup_append hvn' hx hx'
    have hedisj : ∀ e ∈ DfsOut.edgesList [o], e ∉ DfsOut.edgesList rest :=
      fun e he he' => List.disjoint_of_nodup_append hen he he'
    have hen₁ := hen.sublist (List.sublist_append_left _ _)
    have hen₂ := hen.sublist (List.sublist_append_right _ _)
    have hlt₁ : ∀ x ∈ DfsOut.vertsList [o], x < g.nv :=
      fun x hx => hlt x (List.mem_cons_of_mem _ (List.mem_append_left _ hx))
    have hlt₂ : ∀ x ∈ DfsOut.vertsList rest, x < g.nv :=
      fun x hx => hlt x (List.mem_cons_of_mem _ (List.mem_append_right _ hx))
    have ih₁ := refOut_items_of (dirs := dirs) (hv := hv) hvn₁
      (fun x hx => (List.mem_cons.1 hx).elim (fun h => h ▸ hvlt) (hlt₁ x)) hen₁
      (fun e cls child ho dirs' => by
        subst ho
        simp only [DfsOut.vertsList, DfsOut.edgesList, List.append_nil] at hvn₁ hlt₁ hen₁
        exact refTree_items child (d + 1) dirs' (List.nodup_cons.1 hvn₁).2 (List.nodup_cons.1 hen₁).2
          hlt₁)
    have ih₂ := refOuts_items rest v d dirs (refOut g v d dirs o hv).2.2 hvn₂ hen₂
      (fun x hx => (List.mem_cons.1 hx).elim (fun h => h ▸ hvlt) (hlt₂ x))
    obtain ⟨hnd₁, hhv₁, hprov₁⟩ := ih₁
    obtain ⟨hnd₂, hhv₂, hprov₂⟩ := ih₂
    rw [refOuts_cons]
    have hperm : (outsItems ((refOut g v d dirs o hv).1 ++
        (refOuts g v d dirs rest (refOut g v d dirs o hv).2.2).1,
        (refOut g v d dirs o hv).2.1 ++ (refOuts g v d dirs rest (refOut g v d dirs o hv).2.2).2.1,
        (refOuts g v d dirs rest (refOut g v d dirs o hv).2.2).2.2)).Perm
        (outsItems (refOut g v d dirs o hv) ++
          outsItems (refOuts g v d dirs rest (refOut g v d dirs o hv).2.2)) := by
      simp only [outsItems, List.flatMap_append, pieceItems_append]
      refine List.perm_iff_count.2 fun a => ?_
      simp only [List.count_append]; omega
    refine ⟨hperm.nodup_iff.2 (List.Nodup.append hnd₁ hnd₂ fun i hi₁ hi₂ =>
      Prov.ne hdisj hedisj hv1 hv2 hlt₁ hlt₂ hvlt (hprov₁ i hi₁) (hprov₂ i hi₂)),
      fun h => hhv₂ (hhv₁ h), fun i hi => ?_⟩
    rcases List.mem_append.1 (hperm.subset hi) with hi | hi
    · rcases hprov₁ i hi with ⟨x, hx, h⟩ | ⟨e, he, h⟩ | ⟨h1, h2, h⟩
      · exact Or.inl ⟨x, List.mem_append_left _ hx, h⟩
      · exact Or.inr (Or.inl ⟨e, List.mem_append_left _ he, h⟩)
      · exact Or.inr (Or.inr ⟨h1, hhv₂ h2, h⟩)
    · rcases hprov₂ i hi with ⟨x, hx, h⟩ | ⟨e, he, h⟩ | ⟨h1, h2, h⟩
      · exact Or.inl ⟨x, List.mem_append_right _ hx, h⟩
      · exact Or.inr (Or.inl ⟨e, List.mem_append_right _ he, h⟩)
      · refine Or.inr (Or.inr ⟨?_, h2, h⟩)
        cases hv
        · rfl
        · rw [hhv₁ rfl] at h1; cases h1
end

/-- The reference order of a forest visiting every vertex and edge once repeats no item. -/
theorem refOrder_nodup_of_forestOK {g : Graph} : ∀ (forest : List DfsTree),
    (forest.flatMap DfsTree.verts).Nodup → (forest.flatMap DfsTree.edges).Nodup →
    (∀ x ∈ forest.flatMap DfsTree.verts, x < g.nv) →
    (refOrder g forest).Nodup ∧
    ∀ i ∈ refOrder g forest, (∃ x ∈ forest.flatMap DfsTree.verts, i = vertItem x) ∨
      ∃ e ∈ forest.flatMap DfsTree.edges, i = edgeItem g e
  | [], _, _, _ => by simp [refOrder, refBlocks]
  | t :: rest, hvn, hen, hlt => by
    simp only [List.flatMap_cons] at hvn hen hlt ⊢
    have hvn₁ := hvn.sublist (List.sublist_append_left _ _)
    have hvn₂ := hvn.sublist (List.sublist_append_right _ _)
    have hen₁ := hen.sublist (List.sublist_append_left _ _)
    have hen₂ := hen.sublist (List.sublist_append_right _ _)
    have hlt₁ : ∀ x ∈ t.verts, x < g.nv := fun x hx => hlt x (List.mem_append_left _ hx)
    have hlt₂ : ∀ x ∈ rest.flatMap DfsTree.verts, x < g.nv :=
      fun x hx => hlt x (List.mem_append_right _ hx)
    obtain ⟨hnd₁, hprov₁⟩ := refTree_items (g := g) t 0 [] hvn₁ hen₁ hlt₁
    obtain ⟨hnd₂, hprov₂⟩ := refOrder_nodup_of_forestOK rest hvn₂ hen₂ hlt₂
    have hsplit : refOrder g (t :: rest) = ((refTree g t 0 []).2.flatMap StBlock.items ++
        stNest (refTree g t 0 []).1) ++ refOrder g rest := by
      simp [refOrder, refBlocks]
    have hperm : (refOrder g (t :: rest)).Perm (treeItems (refTree g t 0 []) ++ refOrder g rest) := by
      rw [hsplit]; exact ((stNest_perm _).append_left _).append_right _
    refine ⟨hperm.nodup_iff.2 (List.Nodup.append hnd₁ hnd₂ fun i hi₁ hi₂ => ?_), fun i hi => ?_⟩
    · have h₂ := hprov₂ i hi₂
      rcases hprov₁ i hi₁ with ⟨x, hx, rfl⟩ | ⟨e, he, rfl⟩ <;>
        rcases h₂ with ⟨y, hy, h⟩ | ⟨e', he', h⟩
      · exact List.disjoint_of_nodup_append hvn hx (by rw [vertItem_inj h]; exact hy)
      · exact vertItem_ne_edgeItem (hlt₁ x hx) h
      · exact vertItem_ne_edgeItem (hlt₂ y hy) h.symm
      · exact List.disjoint_of_nodup_append hen he (by rw [edgeItem_inj' h]; exact he')
    · rcases List.mem_append.1 (hperm.subset hi) with hi | hi
      · rcases hprov₁ i hi with ⟨x, hx, h⟩ | ⟨e, he, h⟩
        · exact Or.inl ⟨x, List.mem_append_left _ hx, h⟩
        · exact Or.inr ⟨e, List.mem_append_left _ he, h⟩
      · rcases hprov₂ i hi with ⟨x, hx, h⟩ | ⟨e, he, h⟩
        · exact Or.inl ⟨x, List.mem_append_right _ hx, h⟩
        · exact Or.inr ⟨e, List.mem_append_right _ he, h⟩

theorem refOrder_nodup_of_perm {g : Graph} {forest : List DfsTree}
    (hv : (forest.flatMap DfsTree.verts).Perm (List.range g.nv))
    (he : (forest.flatMap DfsTree.edges).Perm (List.range g.ne)) : (refOrder g forest).Nodup :=
  (refOrder_nodup_of_forestOK forest (hv.nodup_iff.2 List.nodup_range)
    (he.nodup_iff.2 List.nodup_range) fun _ hx => List.mem_range.1 (hv.subset hx)).1

end Spqr
