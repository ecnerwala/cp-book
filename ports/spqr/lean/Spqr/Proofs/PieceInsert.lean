import Spqr.Proofs.PlanarInsert
import Spqr.Proofs.PieceSplice
import Spqr.PlanarEmbedLink

/-!
# Hanging a `Q` edge into a capped piece

`Piece.Capped` is the certificate of a capped maximal piece (four exposed ends, two facing pairs,
the two outer ends cofacial in the same witness `ρ`); `Piece.Open2` is the certificate of the
`Q` edge glued onto it (two facing pairs, one at each endpoint). `Capped.insert` performs the two
`link`s of the executable `Q` step and produces the witness `ρ.insert`.
-/

namespace Spqr

open PlanarSpqrTree PlanarSpqrTree.EmbedM

namespace Piece

/-- A capped piece: exposed ends `c0 c1 c2 c3` (slots `0..3`), facing pairs `(c0, c1)` at `u`
and `(c2, c3)` at `v`, and `c0`, `c2` on a common face of the certificate `ρ` that agrees with
the glued rotation `A`. -/
structure Capped (P : Piece) (A : Array (Option Nat)) (ρ : RotationSystem)
    (c0 c1 c2 c3 u v : Nat) : Prop where
  planar : IsPlanarEmbedding P.es P.nVerts ρ
  agrees : P.Agrees A ρ
  unset : ∀ q, P.Mem q → (A[q]? = some none ↔ q = c0 ∨ q = c1 ∨ q = c2 ∨ q = c3)
  pair0 : ∃ l0 l1, P.loc c0 = some l0 ∧ P.loc c1 = some l1 ∧ ρ.get l0 = some l1 ∧
    QE.vert P.es l0 = some u
  pair2 : ∃ l2 l3, P.loc c2 = some l2 ∧ P.loc c3 = some l3 ∧ ρ.get l2 = some l3 ∧
    QE.vert P.es l2 = some v
  face : ∃ l0 l2, P.loc c0 = some l0 ∧ P.loc c2 = some l2 ∧ ρ.SameFaceOrbit l0 l2
  dir0 : c0 % 2 = 0
  dir1 : c1 % 2 = 1
  dir2 : c2 % 2 = 0
  dir3 : c3 % 2 = 1

/-- A piece with two exposed facing pairs `(a, b)` at `u` and `(a', b')` at `v`. -/
structure Open2 (P : Piece) (A : Array (Option Nat)) (ρ : RotationSystem)
    (a b a' b' u v : Nat) : Prop where
  planar : IsPlanarEmbedding P.es P.nVerts ρ
  agrees : P.Agrees A ρ
  pair : ∃ la lb, P.loc a = some la ∧ P.loc b = some lb ∧ ρ.get la = some lb ∧
    QE.vert P.es la = some u
  pair' : ∃ la lb, P.loc a' = some la ∧ P.loc b' = some lb ∧ ρ.get la = some lb ∧
    QE.vert P.es la = some v
  unset : ∀ q, P.Mem q → (A[q]? = some none ↔ q = a ∨ q = b ∨ q = a' ∨ q = b')
  dir_a : a % 2 = 0
  dir_b : b % 2 = 1
  dir_a' : a' % 2 = 0
  dir_b' : b' % 2 = 1

theorem mem_cons (P : Piece) (e q : Nat) :
    ({P with ves := e :: P.ves} : Piece).Mem q ↔ q / 4 = e ∨ P.Mem q := by
  simp [Mem, QE.edge]

theorem loc_cons_self (P : Piece) (e k : Nat) (hk : k < 4) :
    ({P with ves := e :: P.ves} : Piece).loc (4 * e + k) = some k := by
  simp [loc, QE.edge, List.findIdx?_cons, show (4 * e + k) / 4 = e by omega,
    show (4 * e + k) % 4 = k by omega]

theorem loc_cons_of_mem (P : Piece) {e q l : Nat} (he : e ∉ P.ves) (h : P.loc q = some l) :
    ({P with ves := e :: P.ves} : Piece).loc q = some (4 + l) := by
  have := loc_append_right (P := P) [e]
    (by simp only [List.mem_singleton]; intro h'; exact he (h' ▸ mem_of_loc h)) h
  simpa using this

theorem edge_ne_of_mem {P : Piece} {e q : Nat} (he : e ∉ P.ves) (hq : P.Mem q) : q / 4 ≠ e := by
  intro h
  unfold Mem QE.edge at hq
  exact he (h ▸ hq)

theorem es_cons (P : Piece) (e : Nat) :
    ({P with ves := e :: P.ves} : Piece).es = P.ends e :: P.es := by
  simp [es]

theorem x_even {e x : Nat} (hx : x = 0 ∨ x = 2) : (4 * e + x) % 2 = 0 := by omega

theorem x_odd {e x : Nat} (hx : x = 0 ∨ x = 2) : (4 * e + (3 - x)) % 2 = 1 := by omega

theorem four_mul_add_mod {q e : Nat} (hqe : q / 4 = e) : 4 * e + q % 4 = q := by omega

theorem qk_ne {x : Nat} (hx : x = 0 ∨ x = 2) :
    (2 - x ≠ x) ∧ (2 - x ≠ 3 - x) ∧ (x + 1 ≠ x) ∧ (x + 1 ≠ 3 - x) := by omega

theorem x_lt4 {x : Nat} (hx : x = 0 ∨ x = 2) : x < 4 := by omega

theorem x3_lt4 {x : Nat} (hx : x = 0 ∨ x = 2) : 3 - x < 4 := by omega

theorem ne_of_div {q e k : Nat} (hq : q / 4 ≠ e) (hk : k < 4) : 4 * e + k ≠ q :=
  fun h => hq (by omega)

theorem four_add_ne {e k k' : Nat} (_hk : k < 4) (_hk' : k' < 4) (h : k ≠ k') :
    4 * e + k ≠ 4 * e + k' := by omega

theorem qe_cases {q e x : Nat} (hx : x = 0 ∨ x = 2) (hqe : q / 4 = e)
    (h1 : q ≠ 4 * e + (2 - x)) (h3 : q ≠ 4 * e + (x + 1)) :
    q = 4 * e + x ∨ q = 4 * e + (3 - x) := by omega

theorem mod4_cases {e k x : Nat} (hx : x = 0 ∨ x = 2) (hk : k < 4) :
    4 * e + k = 4 * e + x ∨ 4 * e + k = 4 * e + (x + 1) ∨ 4 * e + k = 4 * e + (2 - x) ∨
      4 * e + k = 4 * e + (3 - x) := by omega

theorem Capped.frame {P : Piece} {A B : Array (Option Nat)} {ρ : RotationSystem}
    {c0 c1 c2 c3 u v : Nat} (h : P.Capped A ρ c0 c1 c2 c3 u v)
    (hframe : ∀ q, P.Mem q → B[q]? = A[q]?) : P.Capped B ρ c0 c1 c2 c3 u v where
  planar := h.planar
  agrees := fun q r hq hr => h.agrees q r hq (by rwa [hframe q hq] at hr)
  unset := fun q hq => by rw [hframe q hq]; exact h.unset q hq
  pair0 := h.pair0
  pair2 := h.pair2
  face := h.face
  dir0 := h.dir0
  dir1 := h.dir1
  dir2 := h.dir2
  dir3 := h.dir3

/-- The two `link`s of the `Q` step onto its capped child: `4e + x + 1 ↔ c0`, `4e + 2 - x ↔ c3`.
The new ends are `(4e + x, c1)` at `u` and `(c2, 4e + 3 - x)` at `v`, with witness `ρ.insert`. -/
theorem Capped.insert {P : Piece} {s : EmbedState} {ρ : RotationSystem} {c0 c1 c2 c3 u v e x : Nat}
    (h : P.Capped s.rotAdj ρ c0 c1 c2 c3 u v) (hx : x = 0 ∨ x = 2)
    (he : e ∉ P.ves) (hends : P.ends e = if x = 0 then (u, v) else (v, u)) (huv : u ≠ v)
    (hA : ∀ k, k < 4 → s.rotAdj[4 * e + k]? = some none)
    (hbound : ∀ q, ({P with ves := e :: P.ves} : Piece).Mem q → q < s.rotAdj.size) :
    ∃ ρ', ({P with ves := e :: P.ves} : Piece).Open2
      ((link (some (4 * e + (2 - x))) (some c3)).run
        ((link (some (4 * e + (x + 1))) (some c0)).run s).2).2.rotAdj
      ρ' (4 * e + x) c1 c2 (4 * e + (3 - x)) u v := by
  obtain ⟨l0, l1, hl0, hl1, h01, hvu⟩ := h.pair0
  obtain ⟨l2, l3, hl2, hl3, h23, hvv⟩ := h.pair2
  obtain ⟨l0', l2', hl0', hl2', hface⟩ := h.face
  rw [hl0] at hl0'; rw [hl2] at hl2'
  cases hl0'; cases hl2'
  have hpl := h.planar
  have hsz : ρ.size = 4 * P.ves.length := by rw [hpl.size, es, List.length_map]
  have hl0lt : l0 < ρ.size := by rw [hsz]; exact loc_lt hl0
  have hl2lt : l2 < ρ.size := by rw [hsz]; exact loc_lt hl2
  have hl02 : l0 ≠ l2 := fun h => huv (by rw [h] at hvu; exact Option.some.inj (hvu.symm.trans hvv))
  have hl0e : l0 % 2 = 0 := by rw [loc_mod_two hl0]; exact h.dir0
  have hl2e : l2 % 2 = 0 := by rw [loc_mod_two hl2]; exact h.dir2
  have hface' : SameOrbit (ρ.stepC 3) l0 l2 :=
    (RotationSystem.sameFaceOrbit_iff hpl.total hpl.involution hsz hl0lt).1 hface
  have hr0 : ρ.rot l0 = l1 := RotationSystem.rot_eq_of_get h01
  have hr2 : ρ.rot l2 = l3 := RotationSystem.rot_eq_of_get h23
  have hc13 : c1 ≠ c3 := by
    intro hc
    rw [hc] at hl1
    have : l1 = l3 := Option.some.inj (hl1.symm.trans hl3)
    exact hl02 (RotationSystem.rot_inj hpl.total hpl.involution hl0lt hl2lt (by rw [hr0, hr2, this]))
  have hc02 : c0 ≠ c2 := fun hc => hl02 (by rw [hc] at hl0; exact Option.some.inj (hl0.symm.trans hl2))
  have hce := fun q hq => edge_ne_of_mem he (P := P) (q := q) hq
  have hc0e := hce c0 (mem_of_loc hl0)
  have hc1e := hce c1 (mem_of_loc hl1)
  have hc2e := hce c2 (mem_of_loc hl2)
  have hc3e := hce c3 (mem_of_loc hl3)
  have hd0 := h.dir0; have hd1 := h.dir1; have hd2 := h.dir2; have hd3 := h.dir3
  set P' : Piece := {P with ves := e :: P.ves} with hP'
  have hes' : P'.es = P.ends e :: P.es := es_cons P e
  have hvx : QE.vert (P.ends e :: P.es) x = some u := by
    rw [vert_cons_lt _ _ (by omega)]
    rcases hx with rfl | rfl <;> simp [hends]
  have hvx2 : QE.vert (P.ends e :: P.es) (2 - x) = some v := by
    rw [vert_cons_lt _ _ (by omega)]
    rcases hx with rfl | rfl <;> simp [hends]
  have H : RotationSystem.InsertSetting ρ x l0 l2 :=
    ⟨hpl.total, hpl.involution, hpl.opposite_dir, by rw [hsz]; omega, hx, hl0lt, hl2lt, hl0e,
      hl2e, hl02⟩
  have hL0 : P'.loc (4 * e + x) = some x := loc_cons_self P e x (by omega)
  have hL1 : P'.loc (4 * e + (x + 1)) = some (x + 1) := loc_cons_self P e (x + 1) (by omega)
  have hL2 : P'.loc (4 * e + (2 - x)) = some (2 - x) := loc_cons_self P e (2 - x) (by omega)
  have hL3 : P'.loc (4 * e + (3 - x)) = some (3 - x) := loc_cons_self P e (3 - x) (by omega)
  have hC0 : P'.loc c0 = some (4 + l0) := loc_cons_of_mem P he hl0
  have hC1 : P'.loc c1 = some (4 + l1) := loc_cons_of_mem P he hl1
  have hC2 : P'.loc c2 = some (4 + l2) := loc_cons_of_mem P he hl2
  have hC3 : P'.loc c3 = some (4 + l3) := loc_cons_of_mem P he hl3
  have hmemE : ∀ k, k < 4 → P'.Mem (4 * e + k) := fun k hk =>
    (mem_cons P e _).2 (Or.inl (by omega))
  have hmemP : ∀ q, P.Mem q → P'.Mem q := fun q hq => (mem_cons P e q).2 (Or.inr hq)
  have hb1 : 4 * e + (x + 1) < s.rotAdj.size := hbound _ (hmemE _ (by omega))
  have hb2 : 4 * e + (2 - x) < s.rotAdj.size := hbound _ (hmemE _ (by omega))
  have hbc0 : c0 < s.rotAdj.size := hbound _ (hmemP _ (mem_of_loc hl0))
  have hbc3 : c3 < s.rotAdj.size := hbound _ (hmemP _ (mem_of_loc hl3))
  have hget : ∀ q, ((link (some (4 * e + (2 - x))) (some c3)).run
      ((link (some (4 * e + (x + 1))) (some c0)).run s).2).2.rotAdj[q]? =
      if q = 4 * e + (2 - x) then some (some c3) else if q = c3 then some (some (4 * e + (2 - x)))
      else if q = 4 * e + (x + 1) then some (some c0) else if q = c0 then some (some (4 * e + (x + 1)))
      else s.rotAdj[q]? := by
    intro q
    rw [link_rotAdj_get _ _ _ (by rw [link_rotAdj_size]; exact hb2) (by rw [link_rotAdj_size]; exact hbc3)
      (ne_of_div hc3e (by omega)), link_rotAdj_get _ _ _ hb1 hbc0 (ne_of_div hc0e (by omega))]
  refine ⟨ρ.insert x l0 l2, ?_, ?_, ?_, ?_, ?_, x_even hx, hd1, hd2, x_odd hx⟩
  · rw [hes']
    exact hpl.insert hx hl0lt hl2lt hl0e hl2e hl02 hvu hvv huv hface' hvx hvx2
  · intro q r hq hr
    rw [hget] at hr
    split_ifs at hr with h1 h2 h3 h4
    · subst h1; cases hr
      exact ⟨2 - x, 4 + l3, hL2, hC3, by rw [H.insert_get_x2, hr2]⟩
    · subst h2; cases hr
      exact ⟨4 + l3, 2 - x, hC3, hL2, by rw [← hr2]; exact H.insert_get_a3⟩
    · subst h3; cases hr
      exact ⟨x + 1, 4 + l0, hL1, hC0, H.insert_get_x1⟩
    · subst h4; cases hr
      exact ⟨4 + l0, x + 1, hC0, hL1, H.insert_get_a0⟩
    · rcases (mem_cons P e q).1 hq with hqe | hqP
      · have := hA (q % 4) (Nat.mod_lt _ (by decide))
        rw [four_mul_add_mod hqe, hr] at this
        cases this
      · obtain ⟨lq, lr, hlq, hlr, hqr⟩ := h.agrees q r hqP hr
        have hnot : ¬(q = c0 ∨ q = c1 ∨ q = c2 ∨ q = c3) := by
          rw [← h.unset q hqP, hr]; exact fun h => by cases h
        have hlqlt : lq < ρ.size := by rw [hsz]; exact loc_lt hlq
        have hne : ∀ c l, P.loc c = some l → q ≠ c → lq ≠ l := fun c l hc hqc hl =>
          hqc (loc_injective hlq (hl ▸ hc))
        refine ⟨4 + lq, 4 + lr, loc_cons_of_mem P he hlq, loc_cons_of_mem P he hlr, ?_⟩
        rw [H.insert_get_add hlqlt (hne c0 l0 hl0 (fun h => hnot (Or.inl h)))
          (by rw [hr0]; exact hne c1 l1 hl1 (fun h => hnot (Or.inr (Or.inl h))))
          (hne c2 l2 hl2 (fun h => hnot (Or.inr (Or.inr (Or.inl h)))))
          (by rw [hr2]; exact hne c3 l3 hl3 (fun h => hnot (Or.inr (Or.inr (Or.inr h))))),
          RotationSystem.rot_eq_of_get hqr]
  · exact ⟨x, 4 + l1, hL0, hC1, by rw [H.insert_get_x, hr0], by rw [hes']; exact hvx⟩
  · exact ⟨4 + l2, 3 - x, hC2, hL3, H.insert_get_a2, by rw [hes', vert_cons_add]; exact hvv⟩
  · intro q hq
    rw [hget]
    split_ifs with h1 h2 h3 h4
    · subst h1
      refine iff_of_false (by simp) ?_
      rintro (h | h | h | h)
      · exact four_add_ne (by omega) (by omega) (qk_ne hx).1 h
      · exact ne_of_div hc1e (by omega) h
      · exact ne_of_div hc2e (by omega) h
      · exact four_add_ne (by omega) (by omega) (qk_ne hx).2.1 h
    · subst h2
      refine iff_of_false (by simp) ?_
      rintro (h | h | h | h)
      · exact ne_of_div hc3e (by omega) h.symm
      · exact hc13 h.symm
      · omega
      · exact ne_of_div hc3e (by omega) h.symm
    · subst h3
      refine iff_of_false (by simp) ?_
      rintro (h | h | h | h)
      · exact four_add_ne (by omega) (by omega) (qk_ne hx).2.2.1 h
      · exact ne_of_div hc1e (by omega) h
      · exact ne_of_div hc2e (by omega) h
      · exact four_add_ne (by omega) (by omega) (qk_ne hx).2.2.2 h
    · subst h4
      refine iff_of_false (by simp) ?_
      rintro (h | h | h | h)
      · exact ne_of_div hc0e (by omega) h.symm
      · omega
      · exact hc02 h
      · exact ne_of_div hc0e (by omega) h.symm
    · rcases (mem_cons P e q).1 hq with hqe | hqP
      · have := hA (q % 4) (Nat.mod_lt _ (by decide))
        rw [four_mul_add_mod hqe] at this
        exact iff_of_true this ((qe_cases hx hqe h1 h3).imp id fun h => Or.inr (Or.inr h))
      · rw [h.unset q hqP]
        have hqe := hce q hqP
        constructor
        · rintro (h | h | h | h)
          · exact absurd h h4
          · exact Or.inr (Or.inl h)
          · exact Or.inr (Or.inr (Or.inl h))
          · exact absurd h h2
        · rintro (h | h | h | h)
          · exact (ne_of_div hqe (x_lt4 hx) h.symm).elim
          · exact Or.inr (Or.inl h)
          · exact Or.inr (Or.inr (Or.inl h))
          · exact (ne_of_div hqe (x3_lt4 hx) h.symm).elim

/-- A `Q` item with an `I` child: the one-edge piece with witness `single`. -/
theorem single_open2 (P : Piece) {s : EmbedState} {u v e x : Nat} (hx : x = 0 ∨ x = 2)
    (hends : P.ends e = if x = 0 then (u, v) else (v, u)) (huv : u ≠ v)
    (hu : u < P.nVerts) (hv : v < P.nVerts)
    (hA : ∀ k, k < 4 → s.rotAdj[4 * e + k]? = some none) :
    ({P with ves := [e]} : Piece).Open2 s.rotAdj RotationSystem.single
      (4 * e + x) (4 * e + (x + 1)) (4 * e + (2 - x)) (4 * e + (3 - x)) u v := by
  have hloc : ∀ k, k < 4 → ({P with ves := [e]} : Piece).loc (4 * e + k) = some k :=
    fun k hk => loc_cons_self {P with ves := []} e k hk
  have hmem : ∀ q, ({P with ves := [e]} : Piece).Mem q → q = 4 * e + q % 4 ∧ q % 4 < 4 := by
    intro q hq
    simp only [Mem, QE.edge, List.mem_singleton] at hq
    exact ⟨by omega, Nat.mod_lt _ (by decide)⟩
  have hes : ({P with ves := [e]} : Piece).es = [P.ends e] := by simp [es]
  refine ⟨?_, ?_, ?_, ?_, ?_, by omega, by omega, by omega, by omega⟩
  · rw [hes, hends]
    split
    · exact isPlanarEmbedding_single hu hv huv
    · exact isPlanarEmbedding_single hv hu huv.symm
  · intro q r hq hr
    obtain ⟨hq', hk⟩ := hmem q hq
    have := hA _ hk
    rw [← hq', hr] at this
    cases this
  · refine ⟨x, x + 1, hloc x (by omega), hloc (x + 1) (by omega), ?_, ?_⟩
    · rcases hx with rfl | rfl <;> rfl
    · rw [hes, vert_cons_lt _ _ (by omega), hends]
      rcases hx with rfl | rfl <;> simp
  · refine ⟨2 - x, 3 - x, hloc (2 - x) (by omega), hloc (3 - x) (by omega), ?_, ?_⟩
    · rcases hx with rfl | rfl <;> rfl
    · rw [hes, vert_cons_lt _ _ (by omega), hends]
      rcases hx with rfl | rfl <;> simp
  · intro q hq
    obtain ⟨hq', hk⟩ := hmem q hq
    have := hA _ hk
    rw [← hq'] at this
    rw [this]
    refine iff_of_true rfl ?_
    rw [hq']
    exact mod4_cases hx hk

theorem Open2.size {P : Piece} {A : Array (Option Nat)} {ρ : RotationSystem}
    {a b a' b' u v : Nat} (h : P.Open2 A ρ a b a' b' u v) : ρ.size = 4 * P.ves.length := by
  rw [h.planar.size, es, List.length_map]

/-- The two pairs of an `Open2` piece are at distinct vertices, so they are disjoint. -/
theorem Open2.ne {P : Piece} {A : Array (Option Nat)} {ρ : RotationSystem}
    {a b a' b' u v : Nat} (h : P.Open2 A ρ a b a' b' u v) (huv : u ≠ v) : a ≠ a' ∧ b ≠ b' := by
  obtain ⟨la, lb, hla, hlb, hab, hva⟩ := h.pair
  obtain ⟨la', lb', hla', hlb', hab', hva'⟩ := h.pair'
  have hlalt : la < ρ.size := by rw [h.size]; exact loc_lt hla
  have hla'lt : la' < ρ.size := by rw [h.size]; exact loc_lt hla'
  constructor
  · rintro rfl
    rw [hla] at hla'; cases hla'
    exact huv (Option.some.inj (hva.symm.trans hva'))
  · rintro rfl
    rw [hlb] at hlb'; cases hlb'
    have h1 := (h.planar.involution la hlalt lb hab).2
    have h2 := (h.planar.involution la' hla'lt lb hab').2
    rw [h1] at h2; cases h2
    exact huv (Option.some.inj (hva.symm.trans hva'))

/-- Closing the `v` pair of an `Open2` piece on itself (the lower `V` item has no edges). -/
theorem Open2.close {P : Piece} {s : EmbedState} {ρ : RotationSystem} {a b a' b' u v : Nat}
    (h : P.Open2 s.rotAdj ρ a b a' b' u v) (huv : u ≠ v)
    (hbound : ∀ q, P.Mem q → q < s.rotAdj.size) :
    P.OpenEmbedding ((link (some b') (some a')).run s).2.rotAdj ρ a b u := by
  obtain ⟨hne, hne'⟩ := h.ne huv
  obtain ⟨la, lb, hla, hlb, hab, hva⟩ := h.pair
  obtain ⟨la', lb', hla', hlb', hab', hva'⟩ := h.pair'
  have hla'lt : la' < ρ.size := by rw [h.size]; exact loc_lt hla'
  have hba' : ρ.get lb' = some la' := (h.planar.involution la' hla'lt lb' hab').2
  have hda := h.dir_a; have hdb := h.dir_b; have hda' := h.dir_a'; have hdb' := h.dir_b'
  have hget : ∀ q, ((link (some b') (some a')).run s).2.rotAdj[q]? =
      if q = b' then some (some a') else if q = a' then some (some b') else s.rotAdj[q]? :=
    fun q => link_rotAdj_get _ _ _ (hbound _ (mem_of_loc hlb')) (hbound _ (mem_of_loc hla'))
      (by omega) q
  refine ⟨h.planar, ?_, ⟨la, lb, hla, hlb, hab, hva⟩, ?_, hda, hdb⟩
  · intro q r hq hr
    rw [hget] at hr
    split_ifs at hr with h1 h2
    · subst h1; cases hr; exact ⟨lb', la', hlb', hla', hba'⟩
    · subst h2; cases hr; exact ⟨la', lb', hla', hlb', hab'⟩
    · exact h.agrees q r hq hr
  · intro q hq
    rw [hget]
    split_ifs with h1 h2
    · subst h1
      exact iff_of_false (by simp) (by
        rintro (h | h)
        · omega
        · exact hne' h.symm)
    · subst h2
      exact iff_of_false (by simp) (by
        rintro (h | h)
        · exact hne h.symm
        · omega)
    · rw [h.unset q hq]
      constructor
      · rintro (h | h | h | h)
        · exact Or.inl h
        · exact Or.inr h
        · exact absurd h h2
        · exact absurd h h1
      · rintro (h | h)
        · exact Or.inl h
        · exact Or.inr (Or.inl h)

/-- Splicing the `v` pair of an `Open2` piece onto the open piece `E` at `v` (the lower `V`
item): `b' ↔ w0` then `a' ↔ w1`, as the executable `Q` step does. -/
theorem Open2.splice {P : Piece} {E : List Nat} {s : EmbedState} {ρ ρ₂ : RotationSystem}
    {a b a' b' w0 w1 u v : Nat}
    (h : P.Open2 s.rotAdj ρ a b a' b' u v) (huv : u ≠ v)
    (hw : ({P with ves := E} : Piece).OpenEmbedding s.rotAdj ρ₂ w0 w1 v)
    (hdis : List.Disjoint P.ves E)
    (hsep : ∀ w, HasEdge P.es w → HasEdge ({P with ves := E} : Piece).es w → w = v)
    (hbound : ∀ q, ({P with ves := P.ves ++ E} : Piece).Mem q → q < s.rotAdj.size) :
    ∃ ρ', ({P with ves := P.ves ++ E} : Piece).OpenEmbedding
      ((link (some a') (some w1)).run ((link (some b') (some w0)).run s).2).2.rotAdj ρ' a b u := by
  obtain ⟨hne, hne'⟩ := h.ne huv
  obtain ⟨la, lb, hla, hlb, hab, hva⟩ := h.pair
  obtain ⟨la', lb', hla', hlb', hab', hva'⟩ := h.pair'
  obtain ⟨lw0, lw1, hlw0, hlw1, hw01, hvw⟩ := hw.boundary
  have hlalt : la < ρ.size := by rw [h.size]; exact loc_lt hla
  have hla'lt : la' < ρ.size := by rw [h.size]; exact loc_lt hla'
  have hba : ρ.get lb = some la := (h.planar.involution la hlalt lb hab).2
  have hda := h.dir_a; have hdb := h.dir_b; have hda' := h.dir_a'; have hdb' := h.dir_b'
  have hdw0 := hw.left_dir; have hdw1 := hw.right_dir
  have hmemP : ∀ q, P.Mem q → ({P with ves := P.ves ++ E} : Piece).Mem q := fun q hq =>
    List.mem_append_left _ hq
  have hmemE : ∀ q, ({P with ves := E} : Piece).Mem q →
      ({P with ves := P.ves ++ E} : Piece).Mem q := fun q hq => List.mem_append_right _ hq
  have hPE : ∀ q, P.Mem q → ∀ q', ({P with ves := E} : Piece).Mem q' → q ≠ q' :=
    fun q hq q' hq' heq => hdis hq (heq ▸ hq')
  have hma := mem_of_loc hla; have hmb := mem_of_loc hlb
  have hma' := mem_of_loc hla'; have hmb' := mem_of_loc hlb'
  have hmw0 := mem_of_loc hlw0; have hmw1 := mem_of_loc hlw1
  have hba_ : s.rotAdj[a]? = some none := (h.unset a hma).2 (Or.inl rfl)
  have hbb_ : s.rotAdj[b]? = some none := (h.unset b hmb).2 (Or.inr (Or.inl rfl))
  have hget₀ : ∀ q, ((link (some a) (some b)).run s).2.rotAdj[q]? =
      if q = a then some (some b) else if q = b then some (some a) else s.rotAdj[q]? :=
    fun q => link_rotAdj_get _ _ _ (hbound _ (hmemP _ hma)) (hbound _ (hmemP _ hmb)) (by omega) q
  have h₁ : P.OpenEmbedding ((link (some a) (some b)).run s).2.rotAdj ρ a' b' v := by
    refine ⟨h.planar, ?_, ⟨la', lb', hla', hlb', hab', hva'⟩, ?_, hda', hdb'⟩
    · intro q r hq hr
      rw [hget₀] at hr
      split_ifs at hr with h1 h2
      · subst h1; cases hr; exact ⟨la, lb, hla, hlb, hab⟩
      · subst h2; cases hr; exact ⟨lb, la, hlb, hla, hba⟩
      · exact h.agrees q r hq hr
    · intro q hq
      rw [hget₀]
      split_ifs with h1 h2
      · subst h1
        exact iff_of_false (by simp) (by
        rintro (h | h)
        · exact hne h
        · omega)
      · subst h2
        exact iff_of_false (by simp) (by
        rintro (h | h)
        · omega
        · exact hne' h)
      · rw [h.unset q hq]
        constructor
        · rintro (h | h | h | h)
          · exact absurd h h1
          · exact absurd h h2
          · exact Or.inl h
          · exact Or.inr h
        · rintro (h | h)
          · exact Or.inr (Or.inr (Or.inl h))
          · exact Or.inr (Or.inr (Or.inr h))
  have h₂ : ({P with ves := E} : Piece).OpenEmbedding ((link (some a) (some b)).run s).2.rotAdj
      ρ₂ w0 w1 v := by
    refine ⟨hw.planar, ?_, hw.boundary, ?_, hdw0, hdw1⟩
    · intro q r hq hr
      rw [hget₀, ite_eq_right (hPE a hma q hq).symm, ite_eq_right (hPE b hmb q hq).symm] at hr
      exact hw.agrees q r hq hr
    · intro q hq
      rw [hget₀, ite_eq_right (hPE a hma q hq).symm, ite_eq_right (hPE b hmb q hq).symm]
      exact hw.unset q hq
  have hbound' : ∀ q, ({P with ves := P.ves ++ E} : Piece).Mem q →
      q < ((link (some a) (some b)).run s).2.rotAdj.size := by
    rw [link_rotAdj_size]; exact hbound
  obtain ⟨ρ', hρ'⟩ := OpenEmbedding.splice h₁ h₂ hdis hsep hbound'
  have hgetS : ∀ q, ((link (some b') (some w0)).run ((link (some a) (some b)).run s).2).2.rotAdj[q]? =
      if q = b' then some (some w0) else if q = w0 then some (some b') else
      if q = a then some (some b) else if q = b then some (some a) else s.rotAdj[q]? := by
    intro q
    rw [link_rotAdj_get _ _ _ (by rw [link_rotAdj_size]; exact hbound _ (hmemP _ hmb'))
      (by rw [link_rotAdj_size]; exact hbound _ (hmemE _ hmw0)) (by omega), hget₀]
  have hgetT : ∀ q, ((link (some a') (some w1)).run
      ((link (some b') (some w0)).run s).2).2.rotAdj[q]? =
      if q = a' then some (some w1) else if q = w1 then some (some a') else
      if q = b' then some (some w0) else if q = w0 then some (some b') else s.rotAdj[q]? := by
    intro q
    rw [link_rotAdj_get _ _ _ (by rw [link_rotAdj_size]; exact hbound _ (hmemP _ hma'))
      (by rw [link_rotAdj_size]; exact hbound _ (hmemE _ hmw1)) (by omega),
      link_rotAdj_get _ _ _ (hbound _ (hmemP _ hmb')) (hbound _ (hmemE _ hmw0)) (by omega)]
  obtain ⟨la'', lw1', hla'', hlw1', ha'w1, _⟩ := hρ'.boundary
  have hsz' : ρ'.size = 4 * (P.ves ++ E).length := by rw [hρ'.planar.size, es, List.length_map]
  have hla''lt : la'' < ρ'.size := by rw [hsz']; exact loc_lt hla''
  have hw1a' : ρ'.get lw1' = some la'' := (hρ'.planar.involution la'' hla''lt lw1' ha'w1).2
  have hnea' : a' ≠ w0 := hPE a' hma' w0 hmw0
  have hneb'w1 : b' ≠ w1 := hPE b' hmb' w1 hmw1
  have hnea := hPE a hma; have hneb := hPE b hmb
  refine ⟨ρ', hρ'.planar, ?_, ?_, ?_, hda, hdb⟩
  · intro q r hq hr
    rw [hgetT] at hr
    split_ifs at hr with h1 h2 h3 h4
    · subst h1; cases hr; exact ⟨la'', lw1', hla'', hlw1', ha'w1⟩
    · subst h2; cases hr; exact ⟨lw1', la'', hlw1', hla'', hw1a'⟩
    · subst h3; cases hr
      exact hρ'.agrees _ _ hq (by rw [hgetS, ite_eq_left rfl])
    · subst h4; cases hr
      exact hρ'.agrees _ _ hq (by rw [hgetS, ite_eq_right (by omega), ite_eq_left rfl])
    · have hqa : q ≠ a := fun e => by rw [e, hba_] at hr; cases hr
      have hqb : q ≠ b := fun e => by rw [e, hbb_] at hr; cases hr
      exact hρ'.agrees _ _ hq (by rw [hgetS, ite_eq_right h3, ite_eq_right h4, ite_eq_right hqa, ite_eq_right hqb]; exact hr)
  · obtain ⟨la₁, lb₁, hla₁, hlb₁, hab₁⟩ := hρ'.agrees a b (hmemP _ hma)
      (by rw [hgetS, ite_eq_right (by omega), ite_eq_right (hnea w0 hmw0), ite_eq_left rfl])
    refine ⟨la₁, lb₁, hla₁, hlb₁, hab₁, ?_⟩
    rw [loc_append_left E hla] at hla₁
    cases hla₁
    have hlt : la / 4 < (P.ves.map P.ends).length := by
      rw [List.length_map]; have := loc_lt hla; omega
    show QE.vert ((P.ves ++ E).map P.ends) la = some u
    rw [List.map_append, vert_append_left hlt]
    exact hva
  · intro q hq
    rw [hgetT]
    split_ifs with h1 h2 h3 h4
    · subst h1
      exact iff_of_false (by simp) (by
        rintro (h | h)
        · exact hne h.symm
        · omega)
    · subst h2
      exact iff_of_false (by simp) (by
        rintro (h | h)
        · exact hnea _ hmw1 h.symm
        · exact hneb _ hmw1 h.symm)
    · subst h3
      exact iff_of_false (by simp) (by
        rintro (h | h)
        · omega
        · exact hne' h.symm)
    · subst h4
      exact iff_of_false (by simp) (by
        rintro (h | h)
        · exact hnea _ hmw0 h.symm
        · exact hneb _ hmw0 h.symm)
    · rcases List.mem_append.1 hq with hqP | hqE
      · rw [h.unset q hqP]
        constructor
        · rintro (h | h | h | h)
          · exact Or.inl h
          · exact Or.inr h
          · exact absurd h h1
          · exact absurd h h3
        · rintro (h | h)
          · exact Or.inl h
          · exact Or.inr (Or.inl h)
      · rw [hw.unset q hqE]
        constructor
        · rintro (h | h)
          · exact absurd h h4
          · exact absurd h h2
        · rintro (h | h)
          · exact absurd h.symm (hnea q hqE)
          · exact absurd h.symm (hneb q hqE)

end Piece

end Spqr
