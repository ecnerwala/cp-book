import Spqr.Proofs.PlanarSplice
import Spqr.Proofs.OrbitSplit

/-!
# Inserting an edge between two cofacial exposed slots

`RotationSystem.single` / `loop` are the rotation systems of a one-edge piece.
`RotationSystem.insert` hangs a new edge into two exposed slot pairs of a planar embedding `ρ`
whose slot quarter-edges lie on a common face, by two conjugations; it is again a planar
embedding (the face is split: `F + 2`, `V` unchanged).
-/

namespace Spqr

open Classical

namespace RotationSystem

/-- A single non-loop edge with both ends exposed: `0 ↔ 1`, `2 ↔ 3`. -/
def single : RotationSystem := ⟨#[some 1, some 0, some 3, some 2]⟩

/-- A single loop with its inner ends glued: `1 ↔ 2`, exposed `0 ↔ 3`. -/
def loop : RotationSystem := ⟨#[some 3, some 2, some 1, some 0]⟩

theorem single_size : single.size = 4 := rfl
theorem loop_size : loop.size = 4 := rfl

end RotationSystem

theorem hasEdge_single_iff {u v x : Nat} : HasEdge [(u, v)] x ↔ x = u ∨ x = v := by
  simp only [HasEdge, List.mem_singleton, exists_eq_left]
  exact or_congr eq_comm eq_comm

theorem numNonIsolated_single {u v n : Nat} (hu : u < n) (hv : v < n) :
    numNonIsolated [(u, v)] n = if u = v then 1 else 2 := by
  rw [numNonIsolated_eq_card]
  have : (Finset.range n).filter (HasEdge [(u, v)]) = {u, v} := by
    ext x
    simp only [Finset.mem_filter, Finset.mem_range, hasEdge_single_iff, Finset.mem_insert,
      Finset.mem_singleton]
    constructor
    · exact And.right
    · rintro (rfl | rfl) <;> simp [hu, hv]
  rw [this]
  split_ifs with h
  · subst h; simp
  · exact Finset.card_pair h

theorem edgesConn_single (u v : Nat) : EdgesConn [(u, v)] u v :=
  Relation.ReflTransGen.single (Or.inl (List.mem_singleton.2 rfl))

theorem ccCount_single {u v n : Nat} (hu : u < n) (hv : v < n) : ccCount [(u, v)] n = 1 := by
  unfold ccCount
  rw [Finset.card_eq_one]
  refine ⟨(Finset.range n).filter (EdgesConn [(u, v)] u), ?_⟩
  ext S
  simp only [Finset.mem_image, Finset.mem_filter, Finset.mem_range, Finset.mem_singleton,
    hasEdge_single_iff]
  constructor
  · rintro ⟨x, ⟨hx, rfl | rfl⟩, rfl⟩
    · rfl
    · exact (cls_eq_of_edgesConn (edgesConn_single u x)).symm
  · rintro rfl
    exact ⟨u, ⟨hu, Or.inl rfl⟩, rfl⟩

theorem numComponents_single {u v n : Nat} (hu : u < n) (hv : v < n) :
    numComponents [(u, v)] n = 1 := by
  rw [numComponents_eq_ccCount (by simp [hu, hv]), ccCount_single hu hv]

theorem vert_single (u v : Nat) {q : Nat} (hq : q < 4) :
    QE.vert [(u, v)] q = some (if q / 2 % 2 = 0 then u else v) := by
  have : q / 4 = 0 := by omega
  simp only [QE.vert, QE.edge, QE.side, this]
  rfl

theorem isPlanarEmbedding_single {u v n : Nat} (hu : u < n) (hv : v < n) (huv : u ≠ v) :
    IsPlanarEmbedding [(u, v)] n RotationSystem.single where
  size := rfl
  verts := by simp [hu, hv]
  total := by decide
  involution := by decide
  opposite_dir := by decide
  same_vertex := by
    intro q hq r hr
    have hq : q < 4 := hq
    interval_cases q <;> simp [RotationSystem.single, RotationSystem.get] at hr <;> subst hr <;>
      simp [vert_single]
  vertex_orbits := by rw [numNonIsolated_single hu hv, ite_eq_right huv]; decide
  euler := by
    unfold EulerFormula
    rw [numNonIsolated_single hu hv, ite_eq_right huv, numComponents_single hu hv]
    show RotationSystem.single.numFaceOrbits + 2 * 2 = 2 * (2 * 1 + 1)
    decide

theorem isPlanarEmbedding_loop {u n : Nat} (hu : u < n) :
    IsPlanarEmbedding [(u, u)] n RotationSystem.loop where
  size := rfl
  verts := by simp [hu]
  total := by decide
  involution := by decide
  opposite_dir := by decide
  same_vertex := by
    intro q hq r hr
    have hq : q < 4 := hq
    interval_cases q <;> simp [RotationSystem.loop, RotationSystem.get] at hr <;> subst hr <;>
      simp [vert_single]
  vertex_orbits := by rw [numNonIsolated_single hu hu, ite_eq_left rfl]; decide
  euler := by
    unfold EulerFormula
    rw [numNonIsolated_single hu hu, ite_eq_left rfl, numComponents_single hu hu]
    show RotationSystem.loop.numFaceOrbits + 2 * 1 = 2 * (2 * 1 + 1)
    decide

theorem add4_xor (y : Nat) {c : Nat} (hc : c < 4) : (4 + y) ^^^ c = 4 + (y ^^^ c) := by
  rw [xor_eq4 (4 + y) hc, xor_eq4 y hc]
  have h1 : (4 + y) / 4 = y / 4 + 1 := by omega
  have h2 : (4 + y) % 4 = y % 4 := by omega
  rw [h1, h2]; omega

theorem xor3_div2 (a : Nat) : (a ^^^ 3) / 2 % 2 = 1 - a / 2 % 2 := by
  have h : ∀ j, j < 4 → j ^^^ 3 = 3 - j := by decide
  rw [xor_eq4 a (by decide), h (a % 4) (Nat.mod_lt _ (by decide))]
  omega

theorem vert_across (es : List (Nat × Nat)) (a : Nat) :
    QE.vert es (a ^^^ 3) = (es[a / 4]?).map fun p => if a / 2 % 2 = 0 then p.2 else p.1 := by
  unfold QE.vert QE.edge QE.side
  rw [xor_div4 a (by decide), xor3_div2 a]
  congr 1; funext p
  rcases Nat.mod_two_eq_zero_or_one (a / 2) with h | h <;> simp [h]

theorem vert_cons_add (p : Nat × Nat) (es : List (Nat × Nat)) (y : Nat) :
    QE.vert (p :: es) (4 + y) = QE.vert es y := by
  unfold QE.vert QE.edge QE.side
  have h1 : (4 + y) / 4 = y / 4 + 1 := by omega
  have h2 : (4 + y) / 2 % 2 = y / 2 % 2 := by omega
  rw [h1, h2, List.getElem?_cons_succ]

theorem vert_cons_lt (p : Nat × Nat) (es : List (Nat × Nat)) {q : Nat} (hq : q < 4) :
    QE.vert (p :: es) q = some (if q / 2 % 2 = 0 then p.1 else p.2) := by
  unfold QE.vert QE.edge QE.side
  have : q / 4 = 0 := by omega
  rw [this]; rfl

theorem vert_some_of_lt {es : List (Nat × Nat)} {q : Nat} (hq : q < 4 * es.length) :
    ∃ v, QE.vert es q = some v := by
  unfold QE.vert QE.edge
  rw [List.getElem?_eq_getElem (by omega)]
  exact ⟨_, rfl⟩

theorem identEdges_self (a : Nat) (es : List (Nat × Nat)) : identEdges a a es = es := by
  unfold identEdges
  conv_rhs => rw [← List.map_id es]
  apply List.map_congr_left
  intro p _
  have h1 : ∀ y, ident a a y = y := fun y => by
    unfold ident; split_ifs with h
    · exact h.symm
    · rfl
  simp [h1]

namespace RotationSystem

variable {rs : RotationSystem}

theorem edgesConn_vert_stepC3 {es : List (Nat × Nat)} (ht : rs.Total) (hsv : rs.SameVertex es)
    {m : Nat} (hs : rs.size = 4 * m) {a va vb : Nat} (ha : a < rs.size)
    (hva : QE.vert es a = some va) (hvb : QE.vert es (rs.stepC 3 a) = some vb) :
    EdgesConn es va vb := by
  have ha3 : a ^^^ 3 < rs.size := hs ▸ xor_lt_mul4 (hs ▸ ha) (by decide)
  rw [stepC_eq_rot ht hs (by decide) ha, vert_rot ht hsv ha3, vert_across] at hvb
  unfold QE.vert QE.edge QE.side at hva
  obtain ⟨p, hp⟩ : ∃ p, es[a / 4]? = some p := by
    cases h : es[a / 4]? with
    | none => simp [h] at hva
    | some p => exact ⟨p, rfl⟩
  have hmem : p ∈ es := List.mem_of_getElem? hp
  rw [hp] at hva hvb
  simp only [Option.map_some, Option.some.injEq] at hva hvb
  subst hva hvb
  split_ifs
  · exact Relation.ReflTransGen.single (Or.inl hmem)
  · exact Relation.ReflTransGen.single (Or.inr hmem)

theorem edgesConn_of_sameOrbit {es : List (Nat × Nat)} (ht : rs.Total) (hi : rs.Involution)
    (hsv : rs.SameVertex es) (hlen : rs.size = 4 * es.length) {q r vq vr : Nat}
    (h : SameOrbit (rs.stepC 3) q r)
    (hq : QE.vert es q = some vq) (hr : QE.vert es r = some vr) : EdgesConn es vq vr := by
  let ℓ : Nat → Prop := fun z => ∀ vz, QE.vert es z = some vz → EdgesConn es vq vz
  have hf : ∀ a, ℓ (rs.stepC 3 a) = ℓ a := by
    intro a
    by_cases ha : a < rs.size
    · have hfa : rs.stepC 3 a < rs.size := stepC_lt ht hi hlen (by decide) ha
      obtain ⟨va, hva⟩ := vert_some_of_lt (es := es) (q := a) (by omega)
      obtain ⟨vb, hvb⟩ := vert_some_of_lt (es := es) (q := rs.stepC 3 a) (by omega)
      have hab := edgesConn_vert_stepC3 ht hsv hlen ha hva hvb
      apply propext
      constructor
      · intro h1 vz hz
        rw [hva] at hz; cases hz
        exact (h1 vb hvb).trans (edgesConn_symm hab)
      · intro h1 vz hz
        rw [hvb] at hz; cases hz
        exact (h1 va hva).trans hab
    · rw [stepC_of_ge hlen (by decide) (not_lt.1 ha)]
  have := sameOrbit_invariant ℓ hf h
  have hq' : ℓ q := fun vz hz => by rw [hq] at hz; cases hz; exact Relation.ReflTransGen.refl
  rw [this] at hq'
  exact hq' vr hr

end RotationSystem

theorem numNonIsolated_cons {es : List (Nat × Nat)} {n : Nat} {p : Nat × Nat}
    (h1 : HasEdge es p.1) (h2 : HasEdge es p.2) :
    numNonIsolated (p :: es) n = numNonIsolated es n := by
  rw [numNonIsolated_eq_card, numNonIsolated_eq_card]
  congr 1
  apply Finset.filter_congr
  intro z _
  constructor
  · rintro ⟨q, hq, hz⟩
    rcases List.mem_cons.1 hq with rfl | hq
    · rcases hz with rfl | rfl <;> assumption
    · exact ⟨q, hq, hz⟩
  · rintro ⟨q, hq, hz⟩
    exact ⟨q, List.mem_cons_of_mem _ hq, hz⟩

theorem numComponents_cons {es : List (Nat × Nat)} {n : Nat} {p : Nat × Nat}
    (hes : ∀ q ∈ es, q.1 < n ∧ q.2 < n) (hp : p.1 < n ∧ p.2 < n)
    (hne : p.1 ≠ p.2) (hconn : EdgesConn es p.1 p.2) :
    numComponents (p :: es) n = numComponents es n := by
  rw [numComponents_eq_ccCount hes, numComponents_eq_ccCount (fun q hq => by
    rcases List.mem_cons.1 hq with rfl | hq
    · exact hp
    · exact hes q hq)]
  symm
  apply ccCount_congr_of_conn n (fun q hq => List.mem_cons_of_mem _ hq)
  intro q hq
  rcases List.mem_cons.1 hq with rfl | hq
  · exact Or.inr ⟨hne, hconn⟩
  · exact Or.inl hq

namespace RotationSystem

theorem single_total : single.Total := by decide
theorem single_involution : single.Involution := by decide
theorem single_oppositeDir : single.OppositeDir := by decide
theorem single_numFaceOrbits : single.numFaceOrbits = 2 := by decide
theorem single_numVertexOrbits : single.numVertexOrbits = 4 := by decide

theorem sameOrbit_of_eq {f : Nat → Nat} {q r : Nat} (h : f q = r) : SameOrbit f q r :=
  h ▸ SameOrbit.step q

theorem rot_eq_of_get {rs : RotationSystem} {q r : Nat} (h : rs.get q = some r) : rs.rot q = r := by
  simp [rot, h]

variable (ρ : RotationSystem) (x a0 a2 : Nat)

/-- `single.union ρ` with the new edge's inner quarter-edge `x + 1` hung into the slot `a0`. -/
def U₁ : RotationSystem := (single.union ρ).conj (x + 1) (4 + ρ.rot a0)

/-- The new edge (quarter-edges `0..3`, its `x`-side at `a0`'s slot, its other side at `a2`'s
slot) hung into `ρ`, as two conjugations. -/
def insert : RotationSystem := (U₁ ρ x a0).conj (2 - x) (4 + a2)

theorem U_size : (single.union ρ).size = 4 + ρ.size := by rw [union_size]; rfl
theorem U_get_lt {q : Nat} (hq : q < 4) : (single.union ρ).get q = single.get q :=
  union_get_lt _ _ hq
theorem U_get_add (y : Nat) : (single.union ρ).get (4 + y) = (ρ.get y).map (4 + ·) := by
  rw [union_get_ge _ _ (show single.size ≤ 4 + y by show 4 ≤ 4 + y; omega)]
  show (ρ.get (4 + y - 4)).map (4 + ·) = _
  rw [Nat.add_sub_cancel_left]
theorem U₁_size : (U₁ ρ x a0).size = 4 + ρ.size := by rw [U₁, conj_size, U_size]
theorem insert_size : (ρ.insert x a0 a2).size = 4 + ρ.size := by rw [insert, conj_size, U₁_size]
theorem U₁_get {q : Nat} (hq : q < 4 + ρ.size) : (U₁ ρ x a0).get q =
    ((single.union ρ).get (Equiv.swap (x + 1) (4 + ρ.rot a0) q)).map
      (Equiv.swap (x + 1) (4 + ρ.rot a0)) :=
  conj_get_lt _ _ _ (by rw [U_size]; exact hq)
theorem insert_get {q : Nat} (hq : q < 4 + ρ.size) : (ρ.insert x a0 a2).get q =
    ((U₁ ρ x a0).get (Equiv.swap (2 - x) (4 + a2) q)).map (Equiv.swap (2 - x) (4 + a2)) :=
  conj_get_lt _ _ _ (by rw [U₁_size]; exact hq)

/-- Hypotheses for hanging a new edge into `ρ` at the (even, distinct) slots `a0`, `a2`. -/
structure InsertSetting : Prop where
  total : ρ.Total
  involution : ρ.Involution
  oppositeDir : ρ.OppositeDir
  size4 : ρ.size % 4 = 0
  hx : x = 0 ∨ x = 2
  ha0 : a0 < ρ.size
  ha2 : a2 < ρ.size
  h0 : a0 % 2 = 0
  h2 : a2 % 2 = 0
  ha02 : a0 ≠ a2

namespace InsertSetting

variable {ρ x a0 a2} (H : InsertSetting ρ x a0 a2)
include H

theorem hs : ρ.size = 4 * (ρ.size / 4) := by have := H.size4; omega
theorem hs' : 4 + ρ.size = 4 * (1 + ρ.size / 4) := by have := H.size4; omega
theorem xle : x ≤ 2 := by have := H.hx; omega
theorem xor_facts : x ^^^ 3 = 3 - x ∧ (x + 1) ^^^ 3 = 2 - x ∧ (2 - x) ^^^ 3 = x + 1 ∧
    (3 - x) ^^^ 3 = x ∧ x ^^^ 1 = x + 1 ∧ (x + 1) ^^^ 1 = x ∧ (2 - x) ^^^ 1 = 3 - x ∧
    (3 - x) ^^^ 1 = 2 - x := by
  rcases H.hx with rfl | rfl <;> decide
theorem single_gets : single.get x = some (x + 1) ∧ single.get (x + 1) = some x ∧
    single.get (2 - x) = some (3 - x) ∧ single.get (3 - x) = some (2 - x) := by
  rcases H.hx with rfl | rfl <;> decide

theorem a1_lt : ρ.rot a0 < ρ.size := rot_lt H.total H.involution H.ha0
theorem a3_lt : ρ.rot a2 < ρ.size := rot_lt H.total H.involution H.ha2
theorem a1_odd : ρ.rot a0 % 2 = 1 := by
  have := rot_mod2 H.total H.oppositeDir H.ha0; have := H.h0; omega
theorem a3_odd : ρ.rot a2 % 2 = 1 := by
  have := rot_mod2 H.total H.oppositeDir H.ha2; have := H.h2; omega
theorem a1_ne_a0 : ρ.rot a0 ≠ a0 := rot_ne H.total H.oppositeDir H.ha0
theorem a3_ne_a2 : ρ.rot a2 ≠ a2 := rot_ne H.total H.oppositeDir H.ha2
theorem a1_ne_a2 : ρ.rot a0 ≠ a2 := by have := H.a1_odd; have := H.h2; omega
theorem a3_ne_a0 : ρ.rot a2 ≠ a0 := by have := H.a3_odd; have := H.h0; omega
theorem a1_ne_a3 : ρ.rot a0 ≠ ρ.rot a2 :=
  fun h => H.ha02 (rot_inj H.total H.involution H.ha0 H.ha2 h)
theorem rot_a1 : ρ.rot (ρ.rot a0) = a0 := rot_rot H.total H.involution H.ha0
theorem rot_a3 : ρ.rot (ρ.rot a2) = a2 := rot_rot H.total H.involution H.ha2
theorem rot_ne_a1 {y : Nat} (hy : y < ρ.size) (h : y ≠ a0) : ρ.rot y ≠ ρ.rot a0 :=
  fun e => h (rot_inj H.total H.involution hy H.ha0 e)

theorem U_total : (single.union ρ).Total := union_total _ _ single_total H.total
theorem U_involution : (single.union ρ).Involution :=
  union_involution _ _ single_involution H.involution
theorem U_oppositeDir : (single.union ρ).OppositeDir :=
  union_oppositeDir _ _ (m₁ := 1) rfl single_oppositeDir H.oppositeDir
theorem U_rot_add {y : Nat} (hy : y < ρ.size) : (single.union ρ).rot (4 + y) = 4 + ρ.rot y :=
  rot_eq_of_get (by rw [U_get_add, get_eq_rot H.total hy]; rfl)
theorem U_rot_x1 : (single.union ρ).rot (x + 1) = x :=
  rot_eq_of_get (by rw [U_get_lt _ (by have := H.xle; omega), H.single_gets.2.1])
theorem U_rot_a1 : (single.union ρ).rot (4 + ρ.rot a0) = 4 + a0 := by
  rw [H.U_rot_add H.a1_lt, H.rot_a1]
theorem U_numFaceOrbits : (single.union ρ).numFaceOrbits = 2 + ρ.numFaceOrbits := by
  rw [union_numFaceOrbits _ _ single_total single_involution H.total H.involution (m₁ := 1) rfl
    H.hs, single_numFaceOrbits]
theorem U_numVertexOrbits : (single.union ρ).numVertexOrbits = 4 + ρ.numVertexOrbits := by
  rw [union_numVertexOrbits _ _ single_total single_involution H.total H.involution (m₁ := 1) rfl
    H.hs, single_numVertexOrbits]

omit H in
theorem U_stepC_lt4 {c : Nat} (hc : c < 4) (z : Nat) :
    z < 4 ↔ (single.union ρ).stepC c z < 4 := by
  rw [union_stepC _ _ (m₁ := 1) rfl hc]
  by_cases hz : z < 4
  · have h1 : z ∈ Finset.range single.size := by simp [single_size, hz]
    simp only [unionStep, h1, ↓reduceIte]
    exact iff_of_true hz (stepC_lt single_total single_involution (m := 1) rfl hc hz)
  · refine iff_of_false hz ?_
    have h1 : z ∉ Finset.range single.size := fun h => hz (by simpa [single_size] using h)
    have h2 : ¬z < single.size := hz
    simp only [unionStep, shift]
    rw [ite_eq_right h1, ite_eq_right h2]
    show ¬4 + _ < 4
    omega

theorem U₁_total : (U₁ ρ x a0).Total :=
  conj_total _ _ _ H.U_total (by rw [U_size]; have := H.xle; omega)
    (by rw [U_size]; have := H.a1_lt; omega)
theorem U₁_involution : (U₁ ρ x a0).Involution :=
  conj_involution _ _ _ H.U_total H.U_involution (by rw [U_size]; have := H.xle; omega)
    (by rw [U_size]; have := H.a1_lt; omega)
theorem U₁_oppositeDir : (U₁ ρ x a0).OppositeDir :=
  conj_oppositeDir _ _ _ H.U_total H.U_oppositeDir (by rw [U_size]; have := H.xle; omega)
    (by rw [U_size]; have := H.a1_lt; omega) (by have := H.a1_odd; have := H.hx; omega)
theorem U₁_hs : (U₁ ρ x a0).size = 4 * (1 + ρ.size / 4) := by rw [U₁_size]; exact H.hs'

theorem U₁_get_x : (U₁ ρ x a0).get x = some (4 + ρ.rot a0) := by
  have hx := H.xle
  have e1 : Equiv.swap (x + 1) (4 + ρ.rot a0) x = x :=
    Equiv.swap_apply_of_ne_of_ne (by omega) (by omega)
  rw [U₁_get _ _ _ (by omega), e1, U_get_lt _ (by omega), H.single_gets.1]
  simp only [Option.map_some, Equiv.swap_apply_left]
theorem U₁_get_x1 : (U₁ ρ x a0).get (x + 1) = some (4 + a0) := by
  have hx := H.xle; have h1 := H.a1_ne_a0
  have e2 : Equiv.swap (x + 1) (4 + ρ.rot a0) (4 + a0) = 4 + a0 :=
    Equiv.swap_apply_of_ne_of_ne (by omega) (by omega)
  rw [U₁_get _ _ _ (by omega), Equiv.swap_apply_left, U_get_add, get_eq_rot H.total H.a1_lt,
    H.rot_a1]
  simp only [Option.map_some, e2]
theorem U₁_get_x2 : (U₁ ρ x a0).get (2 - x) = some (3 - x) := by
  have hx := H.hx
  have e1 : Equiv.swap (x + 1) (4 + ρ.rot a0) (2 - x) = 2 - x :=
    Equiv.swap_apply_of_ne_of_ne (by omega) (by omega)
  have e2 : Equiv.swap (x + 1) (4 + ρ.rot a0) (3 - x) = 3 - x :=
    Equiv.swap_apply_of_ne_of_ne (by omega) (by omega)
  rw [U₁_get _ _ _ (by omega), e1, U_get_lt _ (by omega), H.single_gets.2.2.1]
  simp only [Option.map_some, e2]
theorem U₁_get_x3 : (U₁ ρ x a0).get (3 - x) = some (2 - x) := by
  have hx := H.hx
  have e1 : Equiv.swap (x + 1) (4 + ρ.rot a0) (3 - x) = 3 - x :=
    Equiv.swap_apply_of_ne_of_ne (by omega) (by omega)
  have e2 : Equiv.swap (x + 1) (4 + ρ.rot a0) (2 - x) = 2 - x :=
    Equiv.swap_apply_of_ne_of_ne (by omega) (by omega)
  rw [U₁_get _ _ _ (by omega), e1, U_get_lt _ (by omega), H.single_gets.2.2.2]
  simp only [Option.map_some, e2]
theorem U₁_get_a0 : (U₁ ρ x a0).get (4 + a0) = some (x + 1) := by
  have hx := H.xle; have h1 := H.a1_ne_a0; have ha0 := H.ha0
  have e1 : Equiv.swap (x + 1) (4 + ρ.rot a0) (4 + a0) = 4 + a0 :=
    Equiv.swap_apply_of_ne_of_ne (by omega) (by omega)
  rw [U₁_get _ _ _ (by omega), e1, U_get_add, get_eq_rot H.total H.ha0]
  simp only [Option.map_some, Equiv.swap_apply_right]
theorem U₁_get_a1 : (U₁ ρ x a0).get (4 + ρ.rot a0) = some x := by
  have hx := H.xle; have ha1 := H.a1_lt
  have e2 : Equiv.swap (x + 1) (4 + ρ.rot a0) x = x :=
    Equiv.swap_apply_of_ne_of_ne (by omega) (by omega)
  rw [U₁_get _ _ _ (by omega), Equiv.swap_apply_right, U_get_lt _ (by omega), H.single_gets.2.1]
  simp only [Option.map_some, e2]
theorem U₁_get_add {y : Nat} (hy : y < ρ.size) (h0 : y ≠ a0) (h1 : y ≠ ρ.rot a0) :
    (U₁ ρ x a0).get (4 + y) = some (4 + ρ.rot y) := by
  have hx := H.xle; have h2 := H.rot_ne_a1 hy h0
  have e1 : Equiv.swap (x + 1) (4 + ρ.rot a0) (4 + y) = 4 + y :=
    Equiv.swap_apply_of_ne_of_ne (by omega) (by omega)
  have e2 : Equiv.swap (x + 1) (4 + ρ.rot a0) (4 + ρ.rot y) = 4 + ρ.rot y :=
    Equiv.swap_apply_of_ne_of_ne (by omega) (by omega)
  rw [U₁_get _ _ _ (by omega), e1, U_get_add, get_eq_rot H.total hy]
  simp only [Option.map_some, e2]

theorem U₁_rot_x : (U₁ ρ x a0).rot x = 4 + ρ.rot a0 := rot_eq_of_get H.U₁_get_x
theorem U₁_rot_x1 : (U₁ ρ x a0).rot (x + 1) = 4 + a0 := rot_eq_of_get H.U₁_get_x1
theorem U₁_rot_x2 : (U₁ ρ x a0).rot (2 - x) = 3 - x := rot_eq_of_get H.U₁_get_x2
theorem U₁_rot_x3 : (U₁ ρ x a0).rot (3 - x) = 2 - x := rot_eq_of_get H.U₁_get_x3
theorem U₁_rot_a0 : (U₁ ρ x a0).rot (4 + a0) = x + 1 := rot_eq_of_get H.U₁_get_a0
theorem U₁_rot_a1 : (U₁ ρ x a0).rot (4 + ρ.rot a0) = x := rot_eq_of_get H.U₁_get_a1
theorem U₁_rot_a2 : (U₁ ρ x a0).rot (4 + a2) = 4 + ρ.rot a2 :=
  rot_eq_of_get (H.U₁_get_add H.ha2 H.ha02.symm H.a1_ne_a2.symm)
theorem U₁_rot_a3 : (U₁ ρ x a0).rot (4 + ρ.rot a2) = 4 + a2 := by
  rw [rot_eq_of_get (H.U₁_get_add H.a3_lt H.a3_ne_a0 H.a1_ne_a3.symm), H.rot_a3]

theorem U₁_stepC {c : Nat} (hc : c < 4) {q : Nat} (hq : q < 4 + ρ.size) :
    (U₁ ρ x a0).stepC c q = (U₁ ρ x a0).rot (q ^^^ c) :=
  stepC_eq_rot H.U₁_total H.U₁_hs hc (by rw [U₁_size]; exact hq)

theorem U₁_numFaceOrbits : (U₁ ρ x a0).numFaceOrbits = ρ.numFaceOrbits := by
  have hx := H.xle; have ha1 := H.a1_lt
  have := conj_numFaceOrbits (single.union ρ) (x + 1) (4 + ρ.rot a0) H.U_total H.U_involution
    H.U_oppositeDir (by rw [U_size]; omega) (by rw [U_size]; omega) (by omega)
    (by rw [H.U_rot_x1]; omega) (m := 1 + ρ.size / 4) (by rw [U_size]; exact H.hs')
    (· < 4) (U_stepC_lt4 (by decide)) (xor_lt4 (by omega) (by decide))
    (by have := ge_of_xor_ge (m := 1) (c := 3) (q := 4 + ρ.rot a0) (by omega) (by decide); omega)
    (by rw [H.U_rot_x1]; exact xor_lt4 (by omega) (by decide))
    (by rw [H.U_rot_a1]
        have := ge_of_xor_ge (m := 1) (c := 3) (q := 4 + a0) (by omega) (by decide); omega)
  rw [H.U_numFaceOrbits] at this
  show ((single.union ρ).conj (x + 1) (4 + ρ.rot a0)).numFaceOrbits = _
  omega

theorem U₁_numVertexOrbits : (U₁ ρ x a0).numVertexOrbits = 2 + ρ.numVertexOrbits := by
  have hx := H.xle; have ha1 := H.a1_lt
  have := conj_numVertexOrbits (single.union ρ) (x + 1) (4 + ρ.rot a0) H.U_total H.U_involution
    H.U_oppositeDir (by rw [U_size]; omega) (by rw [U_size]; omega) (by omega)
    (by rw [H.U_rot_x1]; omega) (m := 1 + ρ.size / 4) (by rw [U_size]; exact H.hs')
    (· < 4) (U_stepC_lt4 (by decide)) (xor_lt4 (by omega) (by decide))
    (by have := ge_of_xor_ge (m := 1) (c := 1) (q := 4 + ρ.rot a0) (by omega) (by decide); omega)
    (by rw [H.U_rot_x1]; exact xor_lt4 (by omega) (by decide))
    (by rw [H.U_rot_a1]
        have := ge_of_xor_ge (m := 1) (c := 1) (q := 4 + a0) (by omega) (by decide); omega)
  rw [H.U_numVertexOrbits] at this
  show ((single.union ρ).conj (x + 1) (4 + ρ.rot a0)).numVertexOrbits = _
  omega

theorem U₂_vertexC (z : Nat) : (z = 2 - x ∨ z = 3 - x) ↔
    ((U₁ ρ x a0).stepC 1 z = 2 - x ∨ (U₁ ρ x a0).stepC 1 z = 3 - x) := by
  have hx := H.xle
  have hf := H.xor_facts
  have f2 : (U₁ ρ x a0).stepC 1 (2 - x) = 2 - x := by
    rw [H.U₁_stepC (by decide) (by omega), hf.2.2.2.2.2.2.1, H.U₁_rot_x3]
  have f3 : (U₁ ρ x a0).stepC 1 (3 - x) = 3 - x := by
    rw [H.U₁_stepC (by decide) (by omega), hf.2.2.2.2.2.2.2, H.U₁_rot_x2]
  constructor
  · rintro (rfl | rfl)
    · exact Or.inl f2
    · exact Or.inr f3
  · intro h
    by_cases hz : z < 4 + ρ.size
    · have hperm := isPermOn_stepC H.U₁_total H.U₁_involution H.U₁_hs (c := 1) (by decide)
      rw [U₁_size] at hperm
      rcases h with h | h
      · exact Or.inl (hperm.inj z (by simp [hz]) (2 - x) (by simp; omega) (h.trans f2.symm))
      · exact Or.inr (hperm.inj z (by simp [hz]) (3 - x) (by simp; omega) (h.trans f3.symm))
    · rwa [stepC_of_ge H.U₁_hs (by decide) (by rw [U₁_size]; omega)] at h

theorem insert_numVertexOrbits : (ρ.insert x a0 a2).numVertexOrbits = ρ.numVertexOrbits := by
  have hx := H.hx; have ha2 := H.ha2; have ha3 := H.a3_lt
  have hf := H.xor_facts
  have := conj_numVertexOrbits (U₁ ρ x a0) (2 - x) (4 + a2) H.U₁_total H.U₁_involution
    H.U₁_oppositeDir (by rw [U₁_size]; omega) (by rw [U₁_size]; omega) (by omega)
    (by rw [H.U₁_rot_x2]; omega) H.U₁_hs (fun z => z = 2 - x ∨ z = 3 - x) H.U₂_vertexC
    (by rw [hf.2.2.2.2.2.2.1]; exact Or.inr rfl)
    (by rw [add4_xor _ (by decide)]; omega)
    (by rw [H.U₁_rot_x2, hf.2.2.2.2.2.2.2]; exact Or.inl rfl)
    (by rw [H.U₁_rot_a2, add4_xor _ (by decide)]; omega)
  rw [H.U₁_numVertexOrbits] at this
  show ((U₁ ρ x a0).conj (2 - x) (4 + a2)).numVertexOrbits = _
  omega

theorem U₁_transport_step (y : Nat) :
    SameOrbit ((U₁ ρ x a0).stepC 3) (4 + y) (4 + ρ.stepC 3 y) := by
  have hx := H.xle; have hf := H.xor_facts; have ha0 := H.ha0; have ha1 := H.a1_lt
  by_cases hy : y < ρ.size
  · have hy3 : y ^^^ 3 < ρ.size := H.hs ▸ xor_lt_mul4 (H.hs ▸ hy) (by decide)
    rw [stepC_eq_rot H.total H.hs (by decide) hy]
    have hstep : (U₁ ρ x a0).stepC 3 (4 + y) = (U₁ ρ x a0).rot (4 + (y ^^^ 3)) := by
      rw [H.U₁_stepC (by decide) (by omega), add4_xor _ (by decide)]
    have s_x : (U₁ ρ x a0).stepC 3 x = 2 - x := by
      rw [H.U₁_stepC (by decide) (by omega), hf.1, H.U₁_rot_x3]
    have s_x2 : (U₁ ρ x a0).stepC 3 (2 - x) = 4 + a0 := by
      rw [H.U₁_stepC (by decide) (by omega), hf.2.2.1, H.U₁_rot_x1]
    have s_x1 : (U₁ ρ x a0).stepC 3 (x + 1) = 3 - x := by
      rw [H.U₁_stepC (by decide) (by omega), hf.2.1, H.U₁_rot_x2]
    have s_x3 : (U₁ ρ x a0).stepC 3 (3 - x) = 4 + ρ.rot a0 := by
      rw [H.U₁_stepC (by decide) (by omega), hf.2.2.2.1, H.U₁_rot_x]
    by_cases h1 : y ^^^ 3 = ρ.rot a0
    · rw [h1, H.rot_a1]
      have s1 : (U₁ ρ x a0).stepC 3 (4 + y) = x := by rw [hstep, h1, H.U₁_rot_a1]
      exact ((sameOrbit_of_eq s1).trans (sameOrbit_of_eq s_x)).trans (sameOrbit_of_eq s_x2)
    by_cases h0 : y ^^^ 3 = a0
    · rw [h0]
      have s1 : (U₁ ρ x a0).stepC 3 (4 + y) = x + 1 := by rw [hstep, h0, H.U₁_rot_a0]
      exact ((sameOrbit_of_eq s1).trans (sameOrbit_of_eq s_x1)).trans (sameOrbit_of_eq s_x3)
    · have s1 : (U₁ ρ x a0).stepC 3 (4 + y) = 4 + ρ.rot (y ^^^ 3) := by
        rw [hstep, rot_eq_of_get (H.U₁_get_add hy3 h0 h1)]
      exact sameOrbit_of_eq s1
  · rw [stepC_of_ge H.hs (by decide) (not_lt.1 hy)]
    exact SameOrbit.refl _ _

theorem U₁_transport {y y' : Nat} (h : SameOrbit (ρ.stepC 3) y y') :
    SameOrbit ((U₁ ρ x a0).stepC 3) (4 + y) (4 + y') := by
  induction h with
  | rel a b hab => exact hab ▸ H.U₁_transport_step a
  | refl => exact SameOrbit.refl _ _
  | symm a b _ ih => exact ih.symm
  | trans a b c _ _ ih1 ih2 => exact ih1.trans ih2

theorem insert_numFaceOrbits (hface : SameOrbit (ρ.stepC 3) a0 a2) :
    (ρ.insert x a0 a2).numFaceOrbits = ρ.numFaceOrbits + 2 := by
  have hx := H.xle; have hf := H.xor_facts; have ha0 := H.ha0; have ha1 := H.a1_lt
  have ha2 := H.ha2; have ha3 := H.a3_lt
  have s_x : (U₁ ρ x a0).stepC 3 x = 2 - x := by
    rw [H.U₁_stepC (by decide) (by omega), hf.1, H.U₁_rot_x3]
  have s_x2 : (U₁ ρ x a0).stepC 3 (2 - x) = 4 + a0 := by
    rw [H.U₁_stepC (by decide) (by omega), hf.2.2.1, H.U₁_rot_x1]
  have s_x1 : (U₁ ρ x a0).stepC 3 (x + 1) = 3 - x := by
    rw [H.U₁_stepC (by decide) (by omega), hf.2.1, H.U₁_rot_x2]
  have s_x3 : (U₁ ρ x a0).stepC 3 (3 - x) = 4 + ρ.rot a0 := by
    rw [H.U₁_stepC (by decide) (by omega), hf.2.2.2.1, H.U₁_rot_x]
  have r1 : SameOrbit (ρ.stepC 3) (a0 ^^^ 3) (ρ.rot a0) := by
    have : ρ.stepC 3 (a0 ^^^ 3) = ρ.rot a0 := by
      rw [stepC_eq_rot H.total H.hs (by decide) (H.hs ▸ xor_lt_mul4 (H.hs ▸ ha0) (by decide)),
        xor_xor_self]
    exact sameOrbit_of_eq this
  have r3 : SameOrbit (ρ.stepC 3) (ρ.rot a2 ^^^ 3) a2 := by
    have : ρ.stepC 3 (ρ.rot a2 ^^^ 3) = a2 := by
      rw [stepC_eq_rot H.total H.hs (by decide) (H.hs ▸ xor_lt_mul4 (H.hs ▸ ha3) (by decide)),
        xor_xor_self, H.rot_a3]
    exact sameOrbit_of_eq this
  have h1 : SameOrbit ((U₁ ρ x a0).stepC 3) ((2 - x) ^^^ 3) ((4 + a2) ^^^ 3) := by
    rw [hf.2.2.1, add4_xor _ (by decide)]
    refine ((sameOrbit_of_eq s_x1).trans (sameOrbit_of_eq s_x3)).trans ?_
    exact H.U₁_transport (r1.symm.trans (sameOrbit_stepC_xor H.total H.involution H.hs (by decide) hface))
  have h2 : SameOrbit ((U₁ ρ x a0).stepC 3) ((U₁ ρ x a0).rot (2 - x) ^^^ 3)
      ((U₁ ρ x a0).rot (4 + a2) ^^^ 3) := by
    rw [H.U₁_rot_x2, H.U₁_rot_a2, hf.2.2.2.1, add4_xor _ (by decide)]
    refine ((sameOrbit_of_eq s_x).trans (sameOrbit_of_eq s_x2)).trans ?_
    exact H.U₁_transport (hface.trans r3.symm)
  have := conj_numFaceOrbits_split (rs := U₁ ρ x a0) (a := 2 - x) (b := 4 + a2) H.U₁_total
    H.U₁_involution H.U₁_oppositeDir (by rw [U₁_size]; omega) (by rw [U₁_size]; omega) (by omega)
    (by rw [H.U₁_rot_x2]; omega) H.U₁_hs h1 h2
  rw [H.U₁_numFaceOrbits] at this
  show ((U₁ ρ x a0).conj (2 - x) (4 + a2)).numFaceOrbits = _
  omega

omit H in
theorem U_sameVertex {es : List (Nat × Nat)} (p : Nat × Nat) (hsv : ρ.SameVertex es) :
    (single.union ρ).SameVertex (p :: es) := by
  intro q hq r hr
  rw [U_size] at hq
  by_cases h4 : q < 4
  · rw [U_get_lt _ h4] at hr
    interval_cases q <;> simp [single, get] at hr <;> subst hr <;> simp [vert_cons_lt]
  · obtain ⟨y, rfl⟩ : ∃ y, q = 4 + y := ⟨q - 4, by omega⟩
    rw [U_get_add, Option.mem_def, Option.map_eq_some_iff] at hr
    obtain ⟨r', hr', rfl⟩ := hr
    rw [vert_cons_add, vert_cons_add]
    exact hsv y (by omega) r' hr'

theorem U₁_sameVertex {es : List (Nat × Nat)} {p : Nat × Nat} {u : Nat} (hsv : ρ.SameVertex es)
    (hu : QE.vert (p :: es) (x + 1) = some u) (hu0 : QE.vert es a0 = some u) :
    (U₁ ρ x a0).SameVertex (p :: es) := by
  have := conj_sameVertex (single.union ρ) (x + 1) (4 + ρ.rot a0) H.U_total
    (by rw [U_size]; have := H.xle; omega) (by rw [U_size]; have := H.a1_lt; omega)
    (U_sameVertex p hsv) hu (by rw [vert_cons_add, vert_rot H.total hsv H.ha0, hu0])
  rwa [identEdges_self] at this

theorem insert_sameVertex {es : List (Nat × Nat)} {p : Nat × Nat} {u v : Nat}
    (hsv : ρ.SameVertex es)
    (hu : QE.vert (p :: es) (x + 1) = some u) (hu0 : QE.vert es a0 = some u)
    (hv : QE.vert (p :: es) (2 - x) = some v) (hv2 : QE.vert es a2 = some v) :
    (ρ.insert x a0 a2).SameVertex (p :: es) := by
  have := conj_sameVertex (U₁ ρ x a0) (2 - x) (4 + a2) H.U₁_total
    (by rw [U₁_size]; have := H.xle; omega) (by rw [U₁_size]; have := H.ha2; omega)
    (H.U₁_sameVertex hsv hu hu0) hv (by rw [vert_cons_add, hv2])
  rwa [identEdges_self] at this

theorem insert_total : (ρ.insert x a0 a2).Total :=
  conj_total _ _ _ H.U₁_total (by rw [U₁_size]; have := H.xle; omega)
    (by rw [U₁_size]; have := H.ha2; omega)
theorem insert_involution : (ρ.insert x a0 a2).Involution :=
  conj_involution _ _ _ H.U₁_total H.U₁_involution (by rw [U₁_size]; have := H.xle; omega)
    (by rw [U₁_size]; have := H.ha2; omega)
theorem insert_oppositeDir : (ρ.insert x a0 a2).OppositeDir :=
  conj_oppositeDir _ _ _ H.U₁_total H.U₁_oppositeDir (by rw [U₁_size]; have := H.xle; omega)
    (by rw [U₁_size]; have := H.ha2; omega) (by have := H.h2; have := H.hx; omega)

theorem a3_ne_a1 : ρ.rot a2 ≠ ρ.rot a0 := H.a1_ne_a3.symm
theorem x_ne_x2 : x ≠ 2 - x := by have := H.hx; omega
theorem x1_ne_x2 : x + 1 ≠ 2 - x := by have := H.hx; omega
theorem x3_ne_x2 : 3 - x ≠ 2 - x := by have := H.hx; omega

theorem insert_get_x : (ρ.insert x a0 a2).get x = some (4 + ρ.rot a0) := by
  have hx := H.hx; have h12 := H.a1_ne_a2
  rw [insert_get _ _ _ _ (by omega), Equiv.swap_apply_of_ne_of_ne H.x_ne_x2 (by omega),
    H.U₁_get_x]
  simp only [Option.map_some, Equiv.swap_apply_of_ne_of_ne (show 4 + ρ.rot a0 ≠ 2 - x by omega)
    (show 4 + ρ.rot a0 ≠ 4 + a2 by omega)]
theorem insert_get_x1 : (ρ.insert x a0 a2).get (x + 1) = some (4 + a0) := by
  have hx := H.hx; have h02 := H.ha02
  rw [insert_get _ _ _ _ (by omega), Equiv.swap_apply_of_ne_of_ne H.x1_ne_x2 (by omega),
    H.U₁_get_x1]
  simp only [Option.map_some, Equiv.swap_apply_of_ne_of_ne (show 4 + a0 ≠ 2 - x by omega)
    (show 4 + a0 ≠ 4 + a2 by omega)]
theorem insert_get_x2 : (ρ.insert x a0 a2).get (2 - x) = some (4 + ρ.rot a2) := by
  have hx := H.hx; have h32 := H.a3_ne_a2
  rw [insert_get _ _ _ _ (by omega), Equiv.swap_apply_left,
    H.U₁_get_add H.ha2 H.ha02.symm H.a1_ne_a2.symm]
  simp only [Option.map_some, Equiv.swap_apply_of_ne_of_ne (show 4 + ρ.rot a2 ≠ 2 - x by omega)
    (show 4 + ρ.rot a2 ≠ 4 + a2 by omega)]
theorem insert_get_x3 : (ρ.insert x a0 a2).get (3 - x) = some (4 + a2) := by
  have hx := H.hx
  rw [insert_get _ _ _ _ (by omega), Equiv.swap_apply_of_ne_of_ne H.x3_ne_x2 (by omega),
    H.U₁_get_x3]
  simp only [Option.map_some, Equiv.swap_apply_left]
theorem insert_get_a0 : (ρ.insert x a0 a2).get (4 + a0) = some (x + 1) := by
  have hx := H.hx; have h02 := H.ha02; have ha0 := H.ha0
  rw [insert_get _ _ _ _ (by omega), Equiv.swap_apply_of_ne_of_ne (by omega) (by omega),
    H.U₁_get_a0]
  simp only [Option.map_some, Equiv.swap_apply_of_ne_of_ne H.x1_ne_x2
    (show x + 1 ≠ 4 + a2 by omega)]
theorem insert_get_a1 : (ρ.insert x a0 a2).get (4 + ρ.rot a0) = some x := by
  have hx := H.hx; have h12 := H.a1_ne_a2; have ha1 := H.a1_lt
  rw [insert_get _ _ _ _ (by omega), Equiv.swap_apply_of_ne_of_ne (by omega) (by omega),
    H.U₁_get_a1]
  simp only [Option.map_some, Equiv.swap_apply_of_ne_of_ne H.x_ne_x2 (show x ≠ 4 + a2 by omega)]
theorem insert_get_a2 : (ρ.insert x a0 a2).get (4 + a2) = some (3 - x) := by
  have hx := H.hx; have ha2 := H.ha2
  rw [insert_get _ _ _ _ (by omega), Equiv.swap_apply_right, H.U₁_get_x2]
  simp only [Option.map_some, Equiv.swap_apply_of_ne_of_ne H.x3_ne_x2
    (show 3 - x ≠ 4 + a2 by omega)]
theorem insert_get_a3 : (ρ.insert x a0 a2).get (4 + ρ.rot a2) = some (2 - x) := by
  have hx := H.hx; have h32 := H.a3_ne_a2; have ha3 := H.a3_lt
  rw [insert_get _ _ _ _ (by omega), Equiv.swap_apply_of_ne_of_ne (by omega) (by omega),
    H.U₁_get_add H.a3_lt H.a3_ne_a0 H.a3_ne_a1, H.rot_a3]
  simp only [Option.map_some, Equiv.swap_apply_right]
theorem insert_get_add {y : Nat} (hy : y < ρ.size) (h0 : y ≠ a0) (h1 : y ≠ ρ.rot a0)
    (h2 : y ≠ a2) (h3 : y ≠ ρ.rot a2) :
    (ρ.insert x a0 a2).get (4 + y) = some (4 + ρ.rot y) := by
  have hx := H.hx
  have h3' : ρ.rot y ≠ a2 := fun e => h3 (by rw [← e, rot_rot H.total H.involution hy])
  rw [insert_get _ _ _ _ (by omega), Equiv.swap_apply_of_ne_of_ne (by omega) (by omega),
    H.U₁_get_add hy h0 h1]
  simp only [Option.map_some, Equiv.swap_apply_of_ne_of_ne (show 4 + ρ.rot y ≠ 2 - x by omega)
    (show 4 + ρ.rot y ≠ 4 + a2 by omega)]

/-- The two sides of the inserted edge lie on distinct faces. -/
theorem insert_not_sameFaceOrbit (hface : SameOrbit (ρ.stepC 3) a0 a2) :
    ¬(ρ.insert x a0 a2).SameFaceOrbit 0 2 := by
  have hx := H.xle; have hf := H.xor_facts; have ha0 := H.ha0; have ha1 := H.a1_lt
  have ha2 := H.ha2; have ha3 := H.a3_lt
  have s_x : (U₁ ρ x a0).stepC 3 x = 2 - x := by
    rw [H.U₁_stepC (by decide) (by omega), hf.1, H.U₁_rot_x3]
  have s_x2 : (U₁ ρ x a0).stepC 3 (2 - x) = 4 + a0 := by
    rw [H.U₁_stepC (by decide) (by omega), hf.2.2.1, H.U₁_rot_x1]
  have s_x1 : (U₁ ρ x a0).stepC 3 (x + 1) = 3 - x := by
    rw [H.U₁_stepC (by decide) (by omega), hf.2.1, H.U₁_rot_x2]
  have s_x3 : (U₁ ρ x a0).stepC 3 (3 - x) = 4 + ρ.rot a0 := by
    rw [H.U₁_stepC (by decide) (by omega), hf.2.2.2.1, H.U₁_rot_x]
  have r1 : SameOrbit (ρ.stepC 3) (a0 ^^^ 3) (ρ.rot a0) := by
    have : ρ.stepC 3 (a0 ^^^ 3) = ρ.rot a0 := by
      rw [stepC_eq_rot H.total H.hs (by decide) (H.hs ▸ xor_lt_mul4 (H.hs ▸ ha0) (by decide)),
        xor_xor_self]
    exact sameOrbit_of_eq this
  have r3 : SameOrbit (ρ.stepC 3) (ρ.rot a2 ^^^ 3) a2 := by
    have : ρ.stepC 3 (ρ.rot a2 ^^^ 3) = a2 := by
      rw [stepC_eq_rot H.total H.hs (by decide) (H.hs ▸ xor_lt_mul4 (H.hs ▸ ha3) (by decide)),
        xor_xor_self, H.rot_a3]
    exact sameOrbit_of_eq this
  have h1 : SameOrbit ((U₁ ρ x a0).stepC 3) ((2 - x) ^^^ 3) ((4 + a2) ^^^ 3) := by
    rw [hf.2.2.1, add4_xor _ (by decide)]
    refine ((sameOrbit_of_eq s_x1).trans (sameOrbit_of_eq s_x3)).trans ?_
    exact H.U₁_transport (r1.symm.trans (sameOrbit_stepC_xor H.total H.involution H.hs (by decide) hface))
  have h2 : SameOrbit ((U₁ ρ x a0).stepC 3) ((U₁ ρ x a0).rot (2 - x) ^^^ 3)
      ((U₁ ρ x a0).rot (4 + a2) ^^^ 3) := by
    rw [H.U₁_rot_x2, H.U₁_rot_a2, hf.2.2.2.1, add4_xor _ (by decide)]
    refine ((sameOrbit_of_eq s_x).trans (sameOrbit_of_eq s_x2)).trans ?_
    exact H.U₁_transport (hface.trans r3.symm)
  have hns : ¬SameOrbit ((ρ.insert x a0 a2).stepC 3) ((2 - x) ^^^ 3) ((4 + a2) ^^^ 3) :=
    conj_not_sameOrbit_split (rs := U₁ ρ x a0) (a := 2 - x) (b := 4 + a2) H.U₁_total
      H.U₁_involution H.U₁_oppositeDir (by rw [U₁_size]; omega) (by rw [U₁_size]; omega) (by omega)
      (by rw [H.U₁_rot_x2]; omega) H.U₁_hs (by decide) (by decide) h1 h2
  have hsz : (ρ.insert x a0 a2).size = 4 * (1 + ρ.size / 4) := by
    rw [insert_size]; have := H.size4; omega
  intro h
  rw [sameFaceOrbit_iff H.insert_total H.insert_involution hsz (by rw [insert_size]; omega)] at h
  have h' := sameOrbit_stepC_xor H.insert_total H.insert_involution hsz (by decide) h
  have e0 : (0 : Nat) ^^^ 3 = 3 := by decide
  have e2 : (2 : Nat) ^^^ 3 = 1 := by decide
  rw [e0, e2] at h'
  have hr : (ρ.insert x a0 a2).stepC 3 ((4 + a2) ^^^ 3) = 3 - x := by
    rw [stepC_eq_rot H.insert_total hsz (by decide)
      (hsz ▸ xor_lt_mul4 (hsz ▸ (by rw [insert_size]; omega)) (by decide)), xor_xor_self]
    exact rot_eq_of_get H.insert_get_a2
  apply hns
  rcases H.hx with rfl | rfl
  · have e : (2 - 0) ^^^ 3 = 1 := by decide
    rw [e]
    exact ((sameOrbit_of_eq hr).trans h').symm
  · have e : (2 - 2) ^^^ 3 = 3 := by decide
    rw [e]
    exact h'.trans (sameOrbit_of_eq hr).symm

end InsertSetting

end RotationSystem

/-- Hanging a new edge `p` between two exposed slots `a0`, `a2` of a planar embedding `ρ` that
lie on a common face yields a planar embedding of `p :: es`. -/
theorem IsPlanarEmbedding.insert {es : List (Nat × Nat)} {n : Nat} {ρ : RotationSystem}
    (h : IsPlanarEmbedding es n ρ) {x a0 a2 u v : Nat} (hx : x = 0 ∨ x = 2)
    (ha0 : a0 < ρ.size) (ha2 : a2 < ρ.size) (h0 : a0 % 2 = 0) (h2 : a2 % 2 = 0) (ha02 : a0 ≠ a2)
    (hu : QE.vert es a0 = some u) (hv : QE.vert es a2 = some v) (huv : u ≠ v)
    (hface : SameOrbit (ρ.stepC 3) a0 a2) {p : Nat × Nat}
    (hp : QE.vert (p :: es) x = some u) (hp2 : QE.vert (p :: es) (2 - x) = some v) :
    IsPlanarEmbedding (p :: es) n (ρ.insert x a0 a2) := by
  have H : RotationSystem.InsertSetting ρ x a0 a2 :=
    ⟨h.total, h.involution, h.opposite_dir, by rw [h.size]; omega, hx, ha0, ha2, h0, h2, ha02⟩
  have hp1 : QE.vert (p :: es) (x + 1) = some u := by
    have e : (x + 1) / 2 % 2 = x / 2 % 2 := by rcases hx with rfl | rfl <;> rfl
    rw [vert_cons_lt _ _ (by omega), e, ← vert_cons_lt _ _ (by omega), hp]
  have hun : u < n := (hasEdge_of_vert hu).lt_of h.verts
  have hvn : v < n := (hasEdge_of_vert hv).lt_of h.verts
  have hpuv : (p.1 = u ∧ p.2 = v) ∨ (p.1 = v ∧ p.2 = u) := by
    rw [vert_cons_lt _ _ (by omega)] at hp hp2
    rcases hx with rfl | rfl <;> simp at hp hp2
    · exact Or.inl ⟨hp, hp2⟩
    · exact Or.inr ⟨hp2, hp⟩
  have hp1e : HasEdge es p.1 := by
    rcases hpuv with ⟨e, _⟩ | ⟨e, _⟩ <;> rw [e]
    · exact hasEdge_of_vert hu
    · exact hasEdge_of_vert hv
  have hp2e : HasEdge es p.2 := by
    rcases hpuv with ⟨_, e⟩ | ⟨_, e⟩ <;> rw [e]
    · exact hasEdge_of_vert hv
    · exact hasEdge_of_vert hu
  have hconn : EdgesConn es p.1 p.2 := by
    have := RotationSystem.edgesConn_of_sameOrbit h.total h.involution h.same_vertex h.size hface hu hv
    rcases hpuv with ⟨e1, e2⟩ | ⟨e1, e2⟩ <;> rw [e1, e2]
    · exact this
    · exact edgesConn_symm this
  have hpne : p.1 ≠ p.2 := by
    rcases hpuv with ⟨e1, e2⟩ | ⟨e1, e2⟩ <;> rw [e1, e2]
    · exact huv
    · exact huv.symm
  exact {
    size := by rw [RotationSystem.insert_size, h.size, List.length_cons]; omega
    verts := fun q hq => by
      rcases List.mem_cons.1 hq with rfl | hq
      · rcases hpuv with ⟨e1, e2⟩ | ⟨e1, e2⟩ <;> rw [e1, e2] <;> exact ⟨by assumption, by assumption⟩
      · exact h.verts q hq
    total := H.insert_total
    involution := H.insert_involution
    opposite_dir := H.insert_oppositeDir
    same_vertex := H.insert_sameVertex h.same_vertex hp1 hu hp2 hv
    vertex_orbits := by
      rw [H.insert_numVertexOrbits, h.vertex_orbits, numNonIsolated_cons hp1e hp2e]
    euler := by
      have := h.euler
      unfold EulerFormula at this ⊢
      rw [H.insert_numFaceOrbits hface, numNonIsolated_cons hp1e hp2e,
        numComponents_cons h.verts ⟨by rcases hpuv with ⟨e, _⟩ | ⟨e, _⟩ <;> rw [e] <;> assumption,
          by rcases hpuv with ⟨_, e⟩ | ⟨_, e⟩ <;> rw [e] <;> assumption⟩ hpne hconn,
        List.length_cons]
      omega }

end Spqr
