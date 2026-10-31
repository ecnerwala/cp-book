import Spqr.PlanarEmbedTree
import Spqr.PlanarEmbedLeaf
import Spqr.Proofs.PlanarReindex

namespace Spqr.PlanarSpqrTree

variable (t : PlanarSpqrTree)

theorem edgesBelow_subset_root (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape)
    {j : Nat} (hj : j < t.size) : t.edgesBelow j ⊆ t.edgesBelow 0 := by
  induction j using Nat.strong_induction_on with
  | h j ih =>
    by_cases hj0 : j = 0
    · subst j; exact List.Subset.refl _
    obtain ⟨p, hp, hpj⟩ := hwf.preorder.par_lt j (by omega) hj
    have hp' := lt_trans hpj hj
    apply List.Subset.trans ?_ (ih p hpj hp')
    intro e he
    rw [t.edgesBelow_eq hwf hsh p hp', List.mem_append]
    refine Or.inr (List.mem_flatMap.2 ⟨j, ?_, he⟩)
    rw [← t.children_eq, hwf.preorder.ch_eq p hp', List.mem_filter]
    exact ⟨List.mem_range.2 hj, by simpa using hp⟩

theorem mem_edgesBelow_root_iff (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape)
    (e : Nat) : e ∈ t.edgesBelow 0 ↔ e < t.ne := by
  refine ⟨t.mem_edgesBelow_lt hwf, ?_⟩
  intro he
  obtain ⟨j, _, ht, ho⟩ := hwf.bij.edge_index e he
  have hj := t.orig_some_lt hwf ho
  apply t.edgesBelow_subset_root hwf hsh hj
  have hsub := hwf.preorder.subtree_eq j hj
  refine List.mem_filterMap.2 ⟨j, List.mem_range'.2 ⟨0, by omega, by omega⟩, ?_⟩
  rw [t.type_eq_of_lt j hj] at ht
  simp [ht, ho]

theorem edgesBelow_root_perm (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape) :
    (t.edgesBelow 0).Perm (List.range t.ne) := by
  apply (List.perm_ext_iff_of_nodup (t.edgesBelow_nodup hwf 0) List.nodup_range).2
  simp only [List.mem_range, t.mem_edgesBelow_root_iff hwf hsh, implies_true]

theorem pieceBelow_root_perm (g : Graph) (hwf : t.toSpqrTree.WF)
    (hsh : t.toSpqrTree.ChildShape) (hne : t.ne = g.ne) :
    (t.pieceBelow g 0).es.Perm g.edges.toList := by
  have h := (t.edgesBelow_root_perm hwf hsh).map (fun e => g.edges[e]!)
  change ((t.edgesBelow 0).map fun e => g.edges[e]!).Perm g.edges.toList
  rw [hne] at h
  have heq : (List.range g.ne).map (fun e => g.edges[e]!) = g.edges.toList := by
    apply List.ext_getElem
    · simp [Graph.ne]
    · intro i hi hi'
      simp only [List.length_map, List.length_range] at hi
      simp [getElem!_pos g.edges i hi]
      rfl
  rwa [heq] at h

theorem glued_root_of (g : Graph) (hwf : t.toSpqrTree.WF)
    (hsh : t.toSpqrTree.ChildShape) (hne : t.ne = g.ne) (s : EmbedState)
    (h : t.GluedUpTo g 0 s) : IsPlanarEmbedding g.edges.toList g.nv ⟨s.rotAdj⟩ := by
  by_cases hi : 0 < t.size
  · obtain ⟨ρ, hρ, ha, ho, _⟩ := h.piece 0 ⟨by omega, hi, by
      intro p hp
      have hp' := (t.parent_some_iff 0 p).2 hp
      rw [hwf.preorder.root_par] at hp'
      cases hp'⟩
    let P := t.pieceBelow g 0
    let σ : RotationSystem := ⟨s.rotAdj⟩
    let φ := fun q => (P.loc q).getD 0
    have hsize : σ.size = 4 * g.edges.toList.length := by
      change s.rotAdj.size = 4 * g.ne
      rw [h.rot_size, hne]
    have hm : ∀ q, q < σ.size → P.Mem q := by
      intro q hq
      apply (t.mem_edgesBelow_root_iff hwf hsh _).2
      rw [hsize] at hq
      unfold QE.edge
      change q / 4 < t.ne
      rw [hne]; change q / 4 < g.edges.toList.length
      omega
    have hφ : ∀ q, q < σ.size → P.loc q = some (φ q) := by
      intro q hq
      obtain ⟨l, hl⟩ := Piece.loc_exists (hm q hq)
      simp only [φ, hl, Option.getD_some]
    have hφlt : ∀ q, q < σ.size → φ q < ρ.size := by
      intro q hq
      rw [hρ.size]
      simpa only [Piece.es, List.length_map] using Piece.loc_lt (hφ q hq)
    have hglobal : ∀ {q l}, P.loc q = some l → q < σ.size := by
      intro q l hl
      have he := t.mem_pieceBelow_bound hwf g (Piece.mem_of_loc hl)
      change q < s.rotAdj.size
      rwa [h.rot_size]
    have hnexp : ∀ q, ¬s.exposedAt 0 q := by
      intro q ⟨k, hk⟩
      exact (h.outer_slots 0 k q hk).2.1 hwf.preorder.root_type
    have hpairs : ∀ q, q < σ.size → ∃ r, s.rotAdj[q]? = some (some r) := by
      intro q hq
      have hqs : q < s.rotAdj.size := hq
      rw [Array.getElem?_eq_getElem hqs]
      cases hr : s.rotAdj[q] with
      | none =>
        exact (hnexp q ((ho q (hm q hq)).1 (by
          rw [Array.getElem?_eq_getElem hqs, hr]))).elim
      | some r => exact ⟨r, rfl⟩
    have hget : ∀ q, q < σ.size → ρ.get (φ q) = (σ.get q).map φ := by
      intro q hq
      obtain ⟨r, hr⟩ := hpairs q hq
      obtain ⟨lq, lr, hq', hr', hqr⟩ := ha q r (hm q hq) hr
      have heq : lq = φ q := Option.some.inj (hq'.symm.trans (hφ q hq))
      have her : lr = φ r := Option.some.inj (hr'.symm.trans (hφ r (hglobal hr')))
      subst lq; subst lr
      simpa only [RotationSystem.get, σ, hr, Option.bind_some, id_eq, Option.map_some] using hqr
    have hbound : ∀ q, q < σ.size → ∀ r ∈ σ.get q, r < σ.size := by
      intro q hq r hr
      obtain ⟨r', hr'⟩ := hpairs q hq
      have heq : r' = r := by
        simpa only [Option.mem_def, RotationSystem.get, σ, hr', Option.bind_some, id_eq,
          Option.some.injEq] using hr
      subst r'
      obtain ⟨_, _, _, hl, _⟩ := ha q r (hm q hq) hr'
      exact hglobal hl
    have hperm := t.pieceBelow_root_perm g hwf hsh hne
    apply hρ.reindex φ hsize hperm.length_eq (fun p => hperm.mem_iff)
      hφlt (fun q r hq hr he => Piece.loc_injective (hφ q hq) (he ▸ hφ r hr)) hget hbound
    · intro q hq
      exact Piece.loc_mod_two (hφ q hq)
    · intro q hq
      exact t.pieceBelow_loc_vert hwf hne (hφ q hq)
    · intro q hq c hc
      have hqc : q ^^^ c < σ.size := by
        rw [hsize] at hq ⊢
        exact xor_lt_mul4 hq hc
      exact Option.some.inj ((hφ _ hqc).symm.trans (Piece.loc_xor (hφ q hq) hc))
  · have hne0 : t.ne = 0 := by
      by_contra hn
      obtain ⟨j, _, _, hj⟩ := hwf.bij.edge_index 0 (by omega)
      have := t.orig_some_lt hwf hj
      omega
    have hes : g.edges.toList = [] := by
      apply List.length_eq_zero_iff.1
      change g.ne = 0
      omega
    have hs : s.rotAdj = #[] := by
      apply Array.ext
      · simpa [hne0] using h.rot_size
      · intro j _ hj; simp at hj
    simpa only [hes, hs] using isPlanarEmbedding_nil g.nv

end Spqr.PlanarSpqrTree
