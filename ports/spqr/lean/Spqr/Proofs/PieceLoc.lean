import Spqr.PlanarInv

namespace Spqr.Piece

variable {P : Piece} {q r l : Nat}

theorem loc_data (h : P.loc q = some l) :
    ∃ k, k < P.ves.length ∧ P.ves[k]? = some (QE.edge q) ∧ l = 4 * k + q % 4 := by
  obtain ⟨k, hk, hl⟩ := Option.map_eq_some_iff.1 h
  obtain ⟨hkl, heq, _⟩ := List.findIdx?_eq_some_iff_getElem.1 hk
  exact ⟨k, hkl, by simpa [hkl] using heq, hl.symm⟩

theorem mem_of_loc (h : P.loc q = some l) : P.Mem q := by
  obtain ⟨k, _, hk, _⟩ := loc_data h
  exact List.mem_of_getElem? hk

theorem loc_exists (h : P.Mem q) : ∃ l, P.loc q = some l := by
  have hf : P.ves.findIdx? (· == QE.edge q) = some (P.ves.findIdx (· == QE.edge q)) :=
    List.findIdx?_eq_some_of_exists ⟨QE.edge q, h, by simp⟩
  exact ⟨4 * P.ves.findIdx (· == QE.edge q) + q % 4,
    by simp only [loc, hf, Option.map_some]⟩

theorem loc_lt (h : P.loc q = some l) : l < 4 * P.ves.length := by
  obtain ⟨k, hk, _, hl⟩ := loc_data h
  omega

theorem loc_mod_two (h : P.loc q = some l) : l % 2 = q % 2 := by
  obtain ⟨k, _, _, hl⟩ := loc_data h
  omega

theorem vert_loc (h : P.loc q = some l) :
    QE.vert P.es l = some (if QE.side q = 0 then (P.ends (QE.edge q)).1
      else (P.ends (QE.edge q)).2) := by
  obtain ⟨k, _, hk, rfl⟩ := loc_data h
  have he : QE.edge (4 * k + q % 4) = k := by unfold QE.edge; omega
  have hs : QE.side (4 * k + q % 4) = QE.side q := by unfold QE.side; omega
  simp only [QE.vert, he, hs, es, List.getElem?_map, hk, Option.map_some]

theorem loc_injective (hq : P.loc q = some l) (hr : P.loc r = some l) : q = r := by
  obtain ⟨k, _, hk, hl⟩ := loc_data hq
  obtain ⟨k', _, hk', hl'⟩ := loc_data hr
  have hkk : k = k' := by omega
  rw [← hkk, hk] at hk'
  have he : q / 4 = r / 4 := Option.some.inj hk'
  omega

theorem loc_append_left (E : List Nat) (h : P.loc q = some l) :
    ({P with ves := P.ves ++ E} : Piece).loc q = some l := by
  obtain ⟨k, hk, rfl⟩ := Option.map_eq_some_iff.1 h
  simp [loc, List.findIdx?_append, hk]

theorem loc_append_right (E : List Nat) (hne : QE.edge q ∉ E) (h : P.loc q = some l) :
    ({P with ves := E ++ P.ves} : Piece).loc q = some (4 * E.length + l) := by
  have hE : E.findIdx? (· == QE.edge q) = none := by
    simp only [List.findIdx?_eq_none_iff, beq_eq_false_iff_ne]
    intro e he heq
    exact hne (heq ▸ he)
  obtain ⟨k, hk, rfl⟩ := Option.map_eq_some_iff.1 h
  simp [loc, List.findIdx?_append, hE, hk, Nat.mul_add, Nat.add_comm, Nat.add_assoc]

end Spqr.Piece
