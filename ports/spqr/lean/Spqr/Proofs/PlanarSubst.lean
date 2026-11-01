import Spqr.Proofs.PlanarGlue
import Spqr.Proofs.PlanarDegree

/-!
# Substitution of an edge by a planar piece on the same vertex space

`IsPlanarEmbedding.subst`: in a planar embedding of `es₁` on `[0, n)` replace the edge
`es₁[e] = (u, v)` by a planar embedding of `es₂` on the same vertex space whose edge `0` is
`(u, v)` and which meets `es₁ − e` only in `u`, `v`.  The result is `TwoSum.splice`
(PROOF.md §8) with the terminals literally shared, read back on `[0, n)` by collapsing the
shifted copy of the second vertex space.
-/

namespace Spqr

/-- The 2-sum of a substitution: both terminals are the same vertices, `e₂ = 0`. -/
def substSum (es₁ es₂ : List (Nat × Nat)) (n e u v : Nat) : TwoSum :=
  ⟨es₁, n, e, u, v, es₂, n, 0, u, v⟩

theorem hasEdge_of_edgesConn_of_ne {es : List (Nat × Nat)} {a b : Nat} (h : EdgesConn es a b)
    (hab : a ≠ b) : HasEdge es a := by
  rcases Relation.ReflTransGen.cases_head h with rfl | ⟨c, hac, -⟩
  · exact absurd rfl hab
  · rcases hac with hac | hac
    · exact ⟨_, hac, Or.inl rfl⟩
    · exact ⟨_, hac, Or.inr rfl⟩

theorem hasEdge_eraseIdx_of_ne {es : List (Nat × Nat)} {e u v x : Nat} (he : es[e]? = some (u, v))
    (h : HasEdge es x) (hu : x ≠ u) (hv : x ≠ v) : HasEdge (es.eraseIdx e) x := by
  obtain ⟨p, hp, hpx⟩ := h
  obtain ⟨j, hj⟩ := List.getElem?_of_mem hp
  refine ⟨p, List.mem_eraseIdx_iff_getElem?.2 ⟨j, ?_, hj⟩, hpx⟩
  rintro rfl
  rw [he] at hj
  cases hj
  rcases hpx with h | h
  · exact hu h.symm
  · exact hv h.symm

/-- The ends of a non-loop, non-bridge edge are never paired with the edge itself. -/
theorem deg_of_edgesConn {es : List (Nat × Nat)} {n e u v : Nat} {ρ : RotationSystem}
    (h : IsPlanarEmbedding es n ρ) (he : es[e]? = some (u, v)) (huv : u ≠ v)
    (hconn : EdgesConn (es.eraseIdx e) u v) :
    ∀ k, k < 4 → ∀ s ∈ ρ.get (4 * e + k), s / 4 ≠ e := by
  intro k hk s hs
  have he' : e < es.length := (List.getElem?_eq_some_iff.1 he).1
  have hq : 4 * e + k < ρ.size := by rw [h.size]; omega
  have hdiv : (4 * e + k) / 4 = e := by omega
  rw [← RotationSystem.rot_eq_of_get hs]
  have hvert : QE.vert es (4 * e + k) = some (if (4 * e + k) / 2 % 2 = 0 then u else v) := by
    unfold QE.vert QE.edge QE.side; rw [hdiv, he]; rfl
  have hx : HasEdge (es.eraseIdx e) (if (4 * e + k) / 2 % 2 = 0 then u else v) := by
    split
    · exact hasEdge_of_edgesConn_of_ne hconn huv
    · exact hasEdge_of_edgesConn_of_ne (edgesConn_symm hconn) (Ne.symm huv)
  obtain ⟨p, hp, hpx⟩ := hx
  obtain ⟨j, hje, hj⟩ := List.mem_eraseIdx_iff_getElem?.1 hp
  have := IsEmbedding.rot_div_ne h.toIsEmbedding hq (by rw [hdiv]; exact he) huv (k := j)
    (by rw [hdiv]; exact hje) ⟨p, hj, by rw [hvert]; rcases hpx with h | h <;> rw [h] <;> simp⟩
  rwa [hdiv] at this

/-- Well-formedness of the substitution 2-sum. -/
theorem substSum_wf {es₁ es₂ : List (Nat × Nat)} {n e u v : Nat}
    {ρ₁ ρ₂ : RotationSystem} (h₁ : IsPlanarEmbedding es₁ n ρ₁) (h₂ : IsPlanarEmbedding es₂ n ρ₂)
    (he : es₁[e]? = some (u, v)) (he₂ : es₂[0]? = some (u, v)) (huv : u ≠ v)
    (hdeg : ∀ k, k < 4 → ∀ s ∈ ρ₁.get (4 * e + k), s / 4 ≠ e)
    (hface : ¬ρ₂.SameFaceOrbit 0 2) (hconn : EdgesConn (es₂.eraseIdx 0) u v) :
    (substSum es₁ es₂ n e u v).WF ρ₁ ρ₂ :=
  { e₁ := he, e₂ := he₂, uv₁ := huv, uv₂ := huv, deg₁ := hdeg
    deg₂ := deg_of_edgesConn h₂ he₂ huv hconn
    emb₁ := h₁, emb₂ := h₂, face := Or.inr hface, conn := Or.inr hconn }

theorem IsPlanarEmbedding.subst {es₁ es₂ : List (Nat × Nat)} {n e u v : Nat}
    {ρ₁ ρ₂ : RotationSystem} (h₁ : IsPlanarEmbedding es₁ n ρ₁) (h₂ : IsPlanarEmbedding es₂ n ρ₂)
    (he : es₁[e]? = some (u, v)) (he₂ : es₂[0]? = some (u, v)) (huv : u ≠ v)
    (hdeg : ∀ k, k < 4 → ∀ s ∈ ρ₁.get (4 * e + k), s / 4 ≠ e)
    (hface : ¬ρ₂.SameFaceOrbit 0 2) (hconn : EdgesConn (es₂.eraseIdx 0) u v)
    (hsep : ∀ w, HasEdge (es₁.eraseIdx e) w → HasEdge (es₂.eraseIdx 0) w → w = u ∨ w = v) :
    IsPlanarEmbedding (es₁.eraseIdx e ++ es₂.eraseIdx 0) n
      ((substSum es₁ es₂ n e u v).splice ρ₁ ρ₂) := by
  set T := substSum es₁ es₂ n e u v with hT
  have W : T.WF ρ₁ ρ₂ := substSum_wf h₁ h₂ he he₂ huv hdeg hface hconn
  have hsp := T.splice_isPlanarEmbedding W
  have hu : u < n := (h₂.verts _ (List.mem_of_getElem? he₂)).1
  have hv : v < n := (h₂.verts _ (List.mem_of_getElem? he₂)).2
  let f : Nat → Nat := fun x => if x < n then x else x - n
  have hf₁ : ∀ x, x < n → f x = x := fun x hx => ite_eq_left hx
  have hf₂ : ∀ x, f (n + x) = x := fun x => by
    show (if n + x < n then n + x else n + x - n) = x
    rw [ite_eq_right (by omega)]; omega
  have hfv : ∀ w, w < n → f (T.vert₂ w) = w := by
    intro w hw
    show f (if w = u then u else if w = v then v else n + w) = w
    split <;> rename_i h1
    · rw [hf₁ u hu, h1]
    · split <;> rename_i h2
      · rw [hf₁ v hv, h2]
      · exact hf₂ w
  have hmap : mapEdges f T.edges = es₁.eraseIdx e ++ es₂.eraseIdx 0 := by
    show (T.es₁.eraseIdx T.e₁ ++ (T.es₂.eraseIdx T.e₂).map fun p => (T.vert₂ p.1, T.vert₂ p.2)).map
      (fun p => (f p.1, f p.2)) = es₁.eraseIdx e ++ es₂.eraseIdx 0
    rw [List.map_append, List.map_map]
    congr 1
    · refine (List.map_congr_left ?_).trans (List.map_id _)
      intro p hp
      have hpm := h₁.verts p (List.mem_of_mem_eraseIdx hp)
      show (f p.1, f p.2) = p
      rw [hf₁ _ hpm.1, hf₁ _ hpm.2]
    · refine (List.map_congr_left ?_).trans (List.map_id _)
      intro p hp
      have hpm := h₂.verts p (List.mem_of_mem_eraseIdx hp)
      show (f (T.vert₂ p.1), f (T.vert₂ p.2)) = p
      rw [hfv _ hpm.1, hfv _ hpm.2]
  have hvT : ∀ p ∈ T.edges, p.1 < T.nVerts ∧ p.2 < T.nVerts := T.edges_verts W
  have hfT : ∀ p ∈ T.edges, f p.1 < n ∧ f p.2 < n := by
    intro p hp
    have : (f p.1, f p.2) ∈ mapEdges f T.edges := List.mem_map_of_mem hp
    rw [hmap, List.mem_append] at this
    rcases this with h | h
    · exact h₁.verts _ (List.mem_of_mem_eraseIdx h)
    · exact h₂.verts _ (List.mem_of_mem_eraseIdx h)
  have hnn : T.nVerts = n + n := rfl
  have hinj : ∀ x y, HasEdge T.edges x → HasEdge T.edges y → f x = f y → x = y := by
    have key : ∀ x y, HasEdge T.edges x → HasEdge T.edges y → x < n → n ≤ y → f x = f y → False := by
      intro x y hx hy hxn hyn hxy
      have hyn' : y < n + n := by
        obtain ⟨p, hp, hpy⟩ := hy
        have := hvT p hp
        rcases hpy with rfl | rfl <;> omega
      obtain ⟨w, rfl⟩ : ∃ w, y = n + w := ⟨y - n, by omega⟩
      rw [hf₁ x hxn, hf₂] at hxy
      subst hxy
      have h2 := (T.hasEdge_edges_add W (show x < T.n₂ by exact hxn)).1 hy
      have h1 := (T.hasEdge_edges_lt W (show x < T.n₁ by exact hxn)).1 hx
      have h1' : HasEdge (es₁.eraseIdx e) x := hasEdge_eraseIdx_of_ne he h1 h2.2.1 h2.2.2
      have h2' : HasEdge (es₂.eraseIdx 0) x := hasEdge_eraseIdx_of_ne he₂ h2.1 h2.2.1 h2.2.2
      rcases hsep x h1' h2' with h | h
      · exact h2.2.1 h
      · exact h2.2.2 h
    intro x y hx hy hxy
    by_cases hxn : x < n <;> by_cases hyn : y < n
    · rwa [hf₁ x hxn, hf₁ y hyn] at hxy
    · exact (key x y hx hy hxn (not_lt.1 hyn) hxy).elim
    · exact (key y x hy hx hyn (not_lt.1 hxn) hxy.symm).elim
    · have : f x = x - n := ite_eq_right hxn
      have : f y = y - n := ite_eq_right hyn
      omega
  have := IsPlanarEmbedding.map hvT hfT hinj hsp
  rwa [hmap] at this

/-- The rotation of the substitution, as the total function `TwoSum.glue`. -/
theorem substSum_rot {es₁ es₂ : List (Nat × Nat)} {n e u v : Nat} {ρ₁ ρ₂ : RotationSystem}
    (W : (substSum es₁ es₂ n e u v).WF ρ₁ ρ₂) {r : Nat}
    (hr : r < 4 * (es₁.length + es₂.length - 2)) :
    ((substSum es₁ es₂ n e u v).splice ρ₁ ρ₂).rot r = (substSum es₁ es₂ n e u v).glue ρ₁ ρ₂ r :=
  (substSum es₁ es₂ n e u v).splice_rot W hr

theorem substSum_size {es₁ es₂ : List (Nat × Nat)} {n e u v : Nat} {ρ₁ ρ₂ : RotationSystem} :
    ((substSum es₁ es₂ n e u v).splice ρ₁ ρ₂).size = 4 * (es₁.length + es₂.length - 2) :=
  (substSum es₁ es₂ n e u v).splice_size

/-- Connectivity survives the substitution of edge `e = (u, v)` by an edge list connecting `u`
to `v`. -/
theorem edgesConn_subst {es₁ es₂ : List (Nat × Nat)} {e u v a b : Nat}
    (he : es₁[e]? = some (u, v)) (hconn : EdgesConn es₂ u v) (h : EdgesConn es₁ a b) :
    EdgesConn (es₁.eraseIdx e ++ es₂) a b := by
  have hconn' : EdgesConn (es₁.eraseIdx e ++ es₂) u v :=
    edgesConn_mono (fun p hp => List.mem_append_right _ hp) hconn
  have step : ∀ x y, (x, y) ∈ es₁ → EdgesConn (es₁.eraseIdx e ++ es₂) x y := by
    intro x y hxy
    by_cases hxy' : (x, y) = (u, v)
    · cases hxy'; exact hconn'
    · obtain ⟨j, hj⟩ := List.getElem?_of_mem hxy
      have hje : j ≠ e := fun h => hxy' (by rw [h, he] at hj; exact (Option.some.inj hj).symm)
      exact Relation.ReflTransGen.single (Or.inl (List.mem_append_left _
        (List.mem_eraseIdx_iff_getElem?.2 ⟨j, hje, hj⟩)))
  induction h with
  | refl => exact Relation.ReflTransGen.refl
  | tail _ hyz ih =>
    refine ih.trans ?_
    rcases hyz with hyz | hyz
    · exact step _ _ hyz
    · exact edgesConn_symm (step _ _ hyz)

end Spqr
