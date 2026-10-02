import Spqr.Proofs.Contract
import Spqr.Spec

namespace Spqr

theorem SpqrTree.ThreeConnected_congr_undirected {n : Nat} {es es' : List (Nat × Nat)}
    (h : ∀ u v, ((u, v) ∈ es ∨ (v, u) ∈ es) ↔ ((u, v) ∈ es' ∨ (v, u) ∈ es')) :
    SpqrTree.ThreeConnected n es ↔ SpqrTree.ThreeConnected n es' := by
  simp only [SpqrTree.ThreeConnected, h]

namespace Graph

variable {g : Graph} {ok : Nat → Prop}

theorem exists_joins_of_mem {u v : Nat} (h : (u, v) ∈ g.edges.toList) :
    ∃ e, g.Joins e u v := by
  obtain ⟨e, he, hval⟩ := List.getElem_of_mem h
  exact ⟨e, .inl (Array.getElem?_eq_some_iff.2 ⟨by simpa using he, by simpa using hval⟩)⟩

theorem IsEnd.reach {e u w : Nat} (hu : g.IsEnd e u) (hw : g.IsEnd e w)
    (hou : ok u) (how : ok w) : g.Reach ok u w := by
  obtain ⟨v, hj⟩ := hu
  rcases hw.eq_or hj with rfl | rfl
  · exact .refl hou
  · exact (Reach.refl hou).tail hj.adj how

theorem EdgeConn.reach {e f u w : Nat} (h : g.EdgeConn ok e f)
    (hu : g.IsEnd e u) (hw : g.IsEnd f w) (hou : ok u) (how : ok w) :
    g.Reach ok u w := by
  rcases h with rfl | ⟨x, y, hx, hy, hr⟩
  · exact hu.reach hw hou how
  · exact (hu.reach hx hou hr.ok_left).trans
      (hr.trans (hy.reach hw hr.ok_right how))

theorem TwoConnected.two_incident (h2 : g.TwoConnected) {e u v f w : Nat}
    (he : g.Joins e u v) (hf : g.IsEnd f w) (hwu : w ≠ u) (hwv : w ≠ v) :
    ∃ f, f ≠ e ∧ g.IsEnd f u := by
  have hfe : f ≠ e := by
    rintro rfl
    exact (hf.eq_or he).elim hwu hwv
  rcases h2 v e f he.lt hf.lt with hef | ⟨x, y, hx, hy, hr⟩
  · exact (hfe hef.symm).elim
  · have hxu : x = u := (hx.eq_or he).resolve_right hr.ok_left
    subst x
    by_cases hyu : y = u
    · exact ⟨f, hfe, hyu ▸ hy⟩
    · rcases hr.exit (ok' := (· = u)) rfl with hin | ⟨z, z', hz, ⟨f', hj⟩, hz', -⟩
      · exact (hyu hin.ok_right.2).elim
      · have hzu : z = u := hz.ok_right.2
        subst z
        refine ⟨f', ?_, hj.isEnd⟩
        rintro rfl
        rcases he.eq_or hj with ⟨-, h⟩ | ⟨-, h⟩
        · exact hz' h.symm
        · exact hz.ok_right.1 h.symm

theorem ThreeConnected.reach (h3 : g.ThreeConnected) (h2 : g.TwoConnected)
    (hdeg : ∀ v, (∃ e, g.IsEnd e v) → ∃ e f, e ≠ f ∧ g.IsEnd e v ∧ g.IsEnd f v)
    {a b u w : Nat} (hu : ∃ e, g.IsEnd e u) (hw : ∃ e, g.IsEnd e w)
    (hua : u ≠ a) (hub : u ≠ b) (hwa : w ≠ a) (hwb : w ≠ b) :
    g.Reach (fun v => v ≠ a ∧ v ≠ b) u w := by
  obtain ⟨e, e', hee', he, he'⟩ := hdeg u hu
  obtain ⟨f, f', hff', hf, hf'⟩ := hdeg w hw
  by_cases hab : a = b
  · subst b
    exact ((h2 a e f he.lt hf.lt).reach he hf hua hwa).mono fun _ h => ⟨h, h⟩
  · by_contra hn
    apply h3 a b
    refine ⟨hab, .inr ⟨e, e', f, f', he.lt, he'.lt, hf.lt, hf'.lt, hee', ?_, hff', ?_, ?_⟩⟩
    · exact .of_reach he he' (.refl ⟨hua, hub⟩)
    · exact .of_reach hf hf' (.refl ⟨hwa, hwb⟩)
    · exact fun h => hn (h.reach he hf ⟨hua, hub⟩ ⟨hwa, hwb⟩)

theorem TwoConnected.two_incident_of_vertices (h2 : g.TwoConnected) {vs : List Nat}
    (hnd : vs.Nodup) (hlen : 3 ≤ vs.length)
    (hactive : ∀ v ∈ vs, ∃ e, g.IsEnd e v) {u : Nat} (hu : ∃ e, g.IsEnd e u) :
    ∃ e f, e ≠ f ∧ g.IsEnd e u ∧ g.IsEnd f u := by
  obtain ⟨e, v, he⟩ := hu
  have third : ∃ w ∈ vs, w ≠ u ∧ w ≠ v := by
    by_contra! hn
    have hsub : vs ⊆ [u, v] := by
      intro w hw
      by_cases hwu : w = u
      · simp [hwu]
      · simp [hn w hw hwu]
    have := hnd.length_le_of_subset hsub
    simp only [List.length_cons, List.length_nil] at this
    omega
  obtain ⟨w, hw, hwu, hwv⟩ := third
  obtain ⟨f, hf⟩ := hactive w hw
  obtain ⟨f', hfe, hf'⟩ := h2.two_incident he hf hwu hwv
  exact ⟨e, f', hfe.symm, he.isEnd, hf'⟩

theorem ThreeConnected.relabel (h3 : g.ThreeConnected) (h2 : g.TwoConnected) {vs : List Nat}
    (hnd : vs.Nodup) (hlen : 4 ≤ vs.length)
    (hactive : ∀ v ∈ vs, ∃ e, g.IsEnd e v)
    (hends : ∀ e u w, g.Joins e u w → u ∈ vs ∧ w ∈ vs) :
    SpqrTree.ThreeConnected vs.length
      (g.edges.toList.map fun p => (vs.idxOf p.1, vs.idxOf p.2)) := by
  refine ⟨hlen, ?_⟩
  intro a b ha hb
  dsimp only
  intro u w hu hw
  have hdeg : ∀ v, (∃ e, g.IsEnd e v) → ∃ e f, e ≠ f ∧ g.IsEnd e v ∧ g.IsEnd f v :=
    fun _ hv => h2.two_incident_of_vertices hnd (by omega) hactive hv
  have hu' : vs[u] ≠ vs[a] ∧ vs[u] ≠ vs[b] := by
    constructor
    · exact fun h => hu.2.1 (hnd.getElem_inj_iff.1 h)
    · exact fun h => hu.2.2 (hnd.getElem_inj_iff.1 h)
  have hw' : vs[w] ≠ vs[a] ∧ vs[w] ≠ vs[b] := by
    constructor
    · exact fun h => hw.2.1 (hnd.getElem_inj_iff.1 h)
    · exact fun h => hw.2.2 (hnd.getElem_inj_iff.1 h)
  have hr := h3.reach h2 hdeg
    (hactive _ (List.getElem_mem hu.1)) (hactive _ (List.getElem_mem hw.1))
    hu'.1 hu'.2 hw'.1 hw'.2
  have alive : ∀ v ∈ vs, v ≠ vs[a] → v ≠ vs[b] →
      vs.idxOf v < vs.length ∧ vs.idxOf v ≠ a ∧ vs.idxOf v ≠ b := by
    intro v hv hva hvb
    refine ⟨List.idxOf_lt_length_iff.2 hv, ?_, ?_⟩
    · intro h
      apply hva
      simpa [h] using (List.getElem_idxOf (List.idxOf_lt_length_iff.2 hv)).symm
    · intro h
      apply hvb
      simpa [h] using (List.getElem_idxOf (List.idxOf_lt_length_iff.2 hv)).symm
  have tr : ∀ {x y}, g.Reach (fun v => v ≠ vs[a] ∧ v ≠ vs[b]) x y →
      Relation.ReflTransGen
        (fun x y => (x < vs.length ∧ x ≠ a ∧ x ≠ b) ∧ (y < vs.length ∧ y ≠ a ∧ y ≠ b) ∧
          ((x, y) ∈ g.edges.toList.map (fun p => (vs.idxOf p.1, vs.idxOf p.2)) ∨
            (y, x) ∈ g.edges.toList.map (fun p => (vs.idxOf p.1, vs.idxOf p.2))))
        (vs.idxOf x) (vs.idxOf y) := by
    intro x y hr
    induction hr with
    | refl => exact .refl
    | @tail y z hr hadj hz ih =>
      obtain ⟨e, hj⟩ := hadj
      obtain ⟨hy, hz'⟩ := hends e y z hj
      refine ih.tail ⟨alive y hy hr.ok_right.1 hr.ok_right.2, alive z hz' hz.1 hz.2, ?_⟩
      rcases hj with hj | hj
      · exact .inl (List.mem_map.2 ⟨(y, z), by simpa using Array.mem_of_getElem? hj, rfl⟩)
      · exact .inr (List.mem_map.2 ⟨(z, y), by simpa using Array.mem_of_getElem? hj, rfl⟩)
  simpa only [hnd.idxOf_getElem u hu.1, hnd.idxOf_getElem w hw.1] using tr hr

end Graph

end Spqr
