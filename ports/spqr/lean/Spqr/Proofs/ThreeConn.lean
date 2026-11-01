import Spqr.Spec
import Spqr.Proofs.CompCount

/-!
# 3-connectivity and bridges

`SpqrTree.ThreeConnected.edgesConn_eraseIdx`: no edge of a 3-connected graph is a bridge, i.e.
the ends of every edge stay connected after deleting it.  Used for the R skeleton: the cap edge
(and every virtual edge) of an R node is not a bridge of the skeleton.
-/

namespace Spqr

theorem exists_ne_two {n a b : Nat} (hn : 3 ≤ n) : ∃ z, z < n ∧ z ≠ a ∧ z ≠ b := by
  by_cases h0 : 0 ≠ a ∧ 0 ≠ b
  · exact ⟨0, by omega, h0.1, h0.2⟩
  by_cases h1 : 1 ≠ a ∧ 1 ≠ b
  · exact ⟨1, by omega, h1.1, h1.2⟩
  exact ⟨2, by omega, by omega, by omega⟩

theorem exists_ne_three {n a b c : Nat} (hn : 4 ≤ n) : ∃ z, z < n ∧ z ≠ a ∧ z ≠ b ∧ z ≠ c := by
  by_cases h0 : 0 ≠ a ∧ 0 ≠ b ∧ 0 ≠ c
  · exact ⟨0, by omega, h0⟩
  by_cases h1 : 1 ≠ a ∧ 1 ≠ b ∧ 1 ≠ c
  · exact ⟨1, by omega, h1⟩
  by_cases h2 : 2 ≠ a ∧ 2 ≠ b ∧ 2 ≠ c
  · exact ⟨2, by omega, h2⟩
  exact ⟨3, by omega, by omega, by omega, by omega⟩

theorem mem_eraseIdx_of_ne {es : List (Nat × Nat)} {e : Nat} {p : Nat × Nat}
    (hp : p ∈ es) (hpe : es[e]? ≠ some p) : p ∈ es.eraseIdx e := by
  obtain ⟨j, hj⟩ := List.getElem?_of_mem hp
  exact List.mem_eraseIdx_iff_getElem?.2 ⟨j, fun h => hpe (h ▸ hj), hj⟩

/-- A path of a 3-connectivity witness avoiding an end of `es[e] = (u, v)` lives in `es − e`. -/
theorem edgesConn_eraseIdx_of_avoid {es : List (Nat × Nat)} {n e u v a b x y : Nat}
    (he : es[e]? = some (u, v)) (hab : a = u ∨ a = v ∨ b = u ∨ b = v)
    (h : Relation.ReflTransGen
      (fun x y => (x < n ∧ x ≠ a ∧ x ≠ b) ∧ (y < n ∧ y ≠ a ∧ y ≠ b) ∧ ((x, y) ∈ es ∨ (y, x) ∈ es))
      x y) :
    EdgesConn (es.eraseIdx e) x y := by
  induction h with
  | refl => exact Relation.ReflTransGen.refl
  | @tail y z _ hstep ih =>
    obtain ⟨hy, hz, hyz⟩ := hstep
    refine ih.tail ?_
    have key : ∀ p : Nat × Nat, p ∈ es → (p = (y, z) ∨ p = (z, y)) → p ∈ es.eraseIdx e := by
      intro p hp hpyz
      refine mem_eraseIdx_of_ne hp fun hpe => ?_
      rw [he] at hpe
      cases hpe
      rcases hpyz with h | h <;> cases h <;> omega
    rcases hyz with h | h
    · exact Or.inl (key _ h (Or.inl rfl))
    · exact Or.inr (key _ h (Or.inr rfl))

/-- No edge of a 3-connected graph is a bridge. -/
theorem SpqrTree.ThreeConnected.edgesConn_eraseIdx {n : Nat} {es : List (Nat × Nat)}
    (h3 : SpqrTree.ThreeConnected n es) (hbound : ∀ p ∈ es, p.1 < n ∧ p.2 < n) {e u v : Nat}
    (he : es[e]? = some (u, v)) (huv : u ≠ v) : EdgesConn (es.eraseIdx e) u v := by
  obtain ⟨hn, hconn⟩ := h3
  refine Classical.byContradiction fun hnot => ?_
  have hu : u < n := (hbound _ (List.mem_of_getElem? he)).1
  have hv : v < n := (hbound _ (List.mem_of_getElem? he)).2
  have hpath : ∀ a b x y, (a = u ∨ a = v ∨ b = u ∨ b = v) → a < n → b < n →
      x < n → x ≠ a → x ≠ b → y < n → y ≠ a → y ≠ b → EdgesConn (es.eraseIdx e) x y := by
    intro a b x y hab ha hb hx hxa hxb hy hya hyb
    exact edgesConn_eraseIdx_of_avoid he hab (hconn a b ha hb x y ⟨hx, hxa, hxb⟩ ⟨hy, hya, hyb⟩)
  have hA : ∀ x, EdgesConn (es.eraseIdx e) u x → x < n := by
    intro x hx
    rcases Relation.ReflTransGen.cases_tail hx with rfl | ⟨y, -, hyx⟩
    · exact hu
    · rcases hyx with h | h
      · exact (hbound _ (List.mem_of_mem_eraseIdx h)).2
      · exact (hbound _ (List.mem_of_mem_eraseIdx h)).1
  obtain ⟨z, hzn, hzu, hzv⟩ := exists_ne_two (n := n) (a := u) (b := v) (by omega)
  by_cases hex : ∃ x, x ≠ u ∧ EdgesConn (es.eraseIdx e) u x
  · obtain ⟨x, hxu, hx⟩ := hex
    have hxv : x ≠ v := fun h => hnot (h ▸ hx)
    have hxn := hA x hx
    by_cases hey : ∃ y, y < n ∧ y ≠ v ∧ ¬EdgesConn (es.eraseIdx e) u y
    · obtain ⟨y, hyn, hyv, hy⟩ := hey
      have hyu : y ≠ u := fun h => hy (h ▸ Relation.ReflTransGen.refl)
      exact hy (hx.trans (hpath u v x y (Or.inl rfl) hu hv hxn hxu hxv hyn hyu hyv))
    · obtain ⟨w, hwn, hwv, hwu, hwz⟩ := exists_ne_three (n := n) (a := v) (b := u) (c := z) hn
      have hw : EdgesConn (es.eraseIdx e) u w :=
        Classical.byContradiction fun h => hey ⟨w, hwn, hwv, h⟩
      exact hnot (hw.trans (edgesConn_symm
        (hpath u z v w (Or.inl rfl) hu hzn hv (Ne.symm huv) (Ne.symm hzv) hwn hwu hwz)))
  · obtain ⟨w, hwn, hwv, hwz, hwu⟩ := exists_ne_three (n := n) (a := v) (b := z) (c := u) hn
    exact hex ⟨w, hwu,
      hpath v z u w (Or.inr (Or.inl rfl)) hv hzn hu huv (Ne.symm hzu) hwn hwv hwz⟩

/-- Every vertex of a 3-connected graph keeps an edge after deleting any one edge. -/
theorem SpqrTree.ThreeConnected.hasEdge_eraseIdx {n : Nat} {es : List (Nat × Nat)}
    (h3 : SpqrTree.ThreeConnected n es) (hbound : ∀ p ∈ es, p.1 < n ∧ p.2 < n) {e u v k : Nat}
    (he : es[e]? = some (u, v)) (huv : u ≠ v) (hk : k < n) : HasEdge (es.eraseIdx e) k := by
  obtain ⟨hn, hconn⟩ := h3
  have hu : u < n := (hbound _ (List.mem_of_getElem? he)).1
  have hv : v < n := (hbound _ (List.mem_of_getElem? he)).2
  obtain ⟨a, ha, hka, hab⟩ : ∃ a, a < n ∧ k ≠ a ∧ (a = u ∨ a = v) := by
    by_cases hku : k = u
    · exact ⟨v, hv, hku ▸ huv, Or.inr rfl⟩
    · exact ⟨u, hu, hku, Or.inl rfl⟩
  obtain ⟨z, hzn, hza, hzk⟩ := exists_ne_two (n := n) (a := a) (b := k) (by omega)
  obtain ⟨w, hwn, hwa, hwz, hwk⟩ := exists_ne_three (n := n) (a := a) (b := z) (c := k) hn
  have hpath : EdgesConn (es.eraseIdx e) k w :=
    edgesConn_eraseIdx_of_avoid he (by rcases hab with h | h <;> simp [h])
      (hconn a z ha hzn k w ⟨hk, hka, Ne.symm hzk⟩ ⟨hwn, hwa, hwz⟩)
  rcases Relation.ReflTransGen.cases_head hpath with h | ⟨y, hky, -⟩
  · exact absurd h.symm hwk
  · rcases hky with h | h
    · exact ⟨_, h, Or.inl rfl⟩
    · exact ⟨_, h, Or.inr rfl⟩

end Spqr
