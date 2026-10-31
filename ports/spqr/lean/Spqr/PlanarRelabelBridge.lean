import Spqr.PlanarRelabelRowsFold

/-!
# The relabel bridge: `planarRelabelTree` leaves `RelabelNodeR` at every planar R node

`planarRelabel_rowInv` run from `PlanarRelabelState.init` gives `RowInv` of the final state; the
`NodeRow` of node `n` is `RelabelNodeR` once the tree's `nvRange`/`neRange`/`skeleton`/`neRotAdj`
are read off `SpqrTree.ofState`.
-/

namespace Spqr

theorem getD_zero_eq_get! (a : Array Nat) (i : Nat) : a[i]?.getD 0 = a[i]! := by
  rw [getElem!_def]; rfl

/-- `planarRelabel` leaves `RelabelNodeR` at every planar R node (needs `PlanarFinish`: the
`setupNode` writes of other nodes stay inside their own children's slots). -/
theorem planarRelabelTree_relabelNodeR (g : Graph) (w : PlanarWalkState)
    (hwf : Items.WF g w.base.items) (hpf : PlanarFinish g w) (n : Nat)
    (hn : n < (planarRelabelTree g w).size)
    (hR : (planarRelabelTree g w).toSpqrTree.type n = .R)
    (hpl : (planarRelabelTree g w).isPlanar n = true) :
    RelabelNodeR g w (planarRelabelTree g w) n := by
  have hroot : rootItem < w.base.items.size :=
    Nat.lt_of_lt_of_le (by omega : 0 < 1 + g.nv + g.ne) hwf.tree.size
  have hfr : Fresh g w rootItem (PlanarRelabelState.init g w) := fun _ _ _ _ _ _ _ _ => rfl
  obtain ⟨hinv, -⟩ := planarRelabel_rowInv g w hwf hpf w.base.items.size rootItem none none none _ hroot
    (RowInv.init g w) hfr
  generalize hs : ((planarRelabel w.base.items.size rootItem none none none) (PlanarRelabelState.init g w)).2 = s
    at hinv
  have hT : planarRelabelTree g w =
      { SpqrTree.ofState g s.base with nodePlanar := s.aux.nodePlanar, neRotAdj := s.aux.neRotAdj } := by
    rw [← hs]; rfl
  rw [hT] at hn hR hpl ⊢
  have hn' : n < s.base.types.size := hn
  have hR' : s.base.types[n]! = .R := by
    rw [getElem!_pos s.base.types n hn']
    have : s.base.types[n]?.getD .F = .R := hR
    rwa [Array.getElem?_eq_getElem hn'] at this
  have hpl' : s.aux.nodePlanar[n]! = true := by
    have hn2 : n < s.aux.nodePlanar.size := by rw [hinv.np_size]; exact hn'
    rw [getElem!_pos s.aux.nodePlanar n hn2]
    have : s.aux.nodePlanar[n]?.getD false = true := hpl
    rwa [Array.getElem?_eq_getElem hn2] at this
  have hrow := hinv.rows n hn' hR' hpl'
  unfold NodeRow at hrow
  obtain ⟨it, m, children, pos, h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11, h12⟩ := hrow
  have env0 : ((SpqrTree.ofState g s.base).nvRange n).1 = s.base.nvBounds[n]! := getD_zero_eq_get! _ _
  have env1 : ((SpqrTree.ofState g s.base).nvRange n).2 = s.base.nvBounds[n + 1]! := getD_zero_eq_get! _ _
  have ene0 : ((SpqrTree.ofState g s.base).neRange n).1 = s.base.neBounds[n]! := getD_zero_eq_get! _ _
  have ene1 : ((SpqrTree.ofState g s.base).neRange n).2 = s.base.neBounds[n + 1]! := getD_zero_eq_get! _ _
  have hsk : (SpqrTree.ofState g s.base).skeleton n =
      (s.base.nvBounds[n]!, s.base.nvBounds[n + 1]! - 1) :: Items.edgeChildren g w.base.items pos children := by
    have hsz : (nvsArr s).size = s.base.nodeEdges.size := Array.size_map ..
    unfold SpqrTree.skeleton SpqrTree.nodeEdgesOf
    rw [ene0, ene1, h8]
    have hlen : s.base.neBounds[n]! + ((children.filter (· ≥ 1 + g.nv)).map (· - (1 + g.nv))).length + 1 -
        s.base.neBounds[n]! = ((children.filter (· ≥ 1 + g.nv)).map (· - (1 + g.nv))).length + 1 := by omega
    rw [hlen]
    apply List.ext_getElem
    · simp only [List.length_map, List.length_range, List.length_cons, Items.edgeChildren]
    · intro k hk1 hk2
      have hk : k < ((children.filter (· ≥ 1 + g.nv)).map (· - (1 + g.nv))).length + 1 := by simpa using hk1
      have e1 := h10 k (by rw [h8]; omega)
      rw [getElem!_pos ((s.base.nvBounds[n]!, s.base.nvBounds[n + 1]! - 1) ::
        Items.edgeChildren g w.base.items pos children) k hk2] at e1
      rw [← e1]
      unfold nvsArr
      have hb : s.base.neBounds[n]! + k < s.base.nodeEdges.size := by rw [← hsz]; rw [h8] at h9; omega
      rw [Array.getElem!_map' _ _ _ hb, List.getElem_map, List.getElem_map, List.getElem_range]
      show (s.base.nodeEdges[s.base.neBounds[n]! + k]?.getD default).nvs = _
      rw [Array.getElem?_eq_getElem hb, getElem!_pos s.base.nodeEdges _ hb]
      rfl
  unfold RelabelNodeR
  refine ⟨it, m, children, pos, ?_⟩
  dsimp only
  rw [env0, env1, ene0, ene1, hsk]
  exact ⟨h1, h2, h3, h4, h5, h6, h7, h8, rfl, h11, h12⟩

end Spqr
