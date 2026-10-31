import Spqr.PlanarEmbedTree
import Spqr.Proofs.PieceSplice

namespace Spqr.PlanarSpqrTree

variable (t : PlanarSpqrTree)

theorem v_child_boundary (g : Graph) (hwf : t.toSpqrTree.WF) (hne : t.ne = g.ne)
    (hsep : t.toSpqrTree.PieceSep g) {i v j : Nat} (hi : i < t.size)
    (ht : t.toSpqrTree.type i = .V) (hv : t.origId[i]! = some v)
    (s : EmbedState) (h : t.GluedUpTo g (i + 1) s) (hj : j ∈ t.children i) :
    ∃ a b ρ, s.outerE[j]![0]! = some a ∧ s.outerE[j]![1]! = some b ∧
      (t.pieceBelow g j).OpenEmbedding s.rotAdj ρ a b v := by
  have hm := t.child_maximal hwf hi hj
  obtain ⟨hjs, hp, hij⟩ := t.child_data hwf hi hj
  obtain ⟨ρ, hρ, ha, ho, hpair⟩ := h.piece j hm
  have hs := fun k q hk => (h.outer_slots j k q hk).2.2 i hp (Or.inr ht)
  have hpr := hpair 0 (by omega)
  obtain ⟨q, k, hk⟩ := h.outer_present j (by omega) hjs i hp ht
    (t.v_children_nonempty hwf hsep hi ht hj)
  have hkl : k < 2 := hs k q hk
  have hex : ∃ a, s.outerE[j]?.bind (fun o => o[0]?) = some (some a) := by
    interval_cases k
    · exact ⟨q, hk⟩
    · exact hpr.2 q hk
  obtain ⟨a, ha0⟩ := hex
  obtain ⟨b, hb1, la, lb, hla, hlb, hab⟩ := hpr.1 a ha0
  refine ⟨a, b, ρ, (outer_some_iff s j 0 a).2 ha0, (outer_some_iff s j 1 b).2 hb1,
    hρ, ha, ⟨la, lb, hla, hlb, hab, ?_⟩, ?_, h.outer_dir j 0 a ha0, h.outer_dir j 1 b hb1⟩
  · rw [t.pieceBelow_loc_vert hwf hne hla]
    exact h.outer_at_vertex j i v a hp ht hv ⟨0, ha0⟩
  · intro q hq
    rw [ho q hq]
    constructor
    · rintro ⟨k, hk⟩
      have hk' := hs k q hk
      interval_cases k
      · exact Or.inl (Option.some.inj (Option.some.inj (hk.symm.trans ha0)))
      · exact Or.inr (Option.some.inj (Option.some.inj (hk.symm.trans hb1)))
    · rintro (rfl | rfl)
      · exact ⟨0, ha0⟩
      · exact ⟨1, hb1⟩

end Spqr.PlanarSpqrTree
