import Spqr.Proofs.PlanarVStep
import Spqr.Proofs.PlanarSubstStep
import Spqr.Proofs.PlanarUninsert

/-!
# Gluing an `R` node: skeleton + corner pieces + substituted children

Generic data `RG`: a planar embedding `ρ₀` of the skeleton `cap :: sk` (`m = sk.length` virtual
edges, edge `j` of `sk` at quarter-edges `4 * (j + 1) + r`), `nV` open pieces (`ves v`, `vρ v`,
boundary pair `vw0 v ↔ vw1 v` at the vertex `vx v`) hung at the skeleton corners `vta v ↔ ρ₀.rot
(vta v)` and `m` capped pieces (`ces j`, `cρ j`, exposed slots `cl j 0..3`) substituted for the
virtual edges. Quarter-edges of the final system are addressed by *codes*: `0..3` the cap,
`4 * (j + 1) + r` skeleton edge `j`, `vIdx v + r` piece `v`, `cIdx j + r` child `j`; after `k`
substitutions code `c ≥ 4` sits at position `c - 4 * k` (`pos`).

* `vstage`: the `nV` one-sums (`vstep`), rotation `tgtV`;
* `cstage`: the `m` substitutions (`subst_step`), rotation through `res`/`skt`;
* `glue`: the cap removed (`IsPlanarEmbedding.uninsert`): a planar embedding of
  `vflat nV ++ cflat m` with explicit rotation and the cap's two sides cofacial.
-/

namespace Spqr

structure RG where
  n : Nat
  cap : Nat × Nat
  sk : List (Nat × Nat)
  ρ₀ : RotationSystem
  nV : Nat
  ves : Nat → List (Nat × Nat)
  vρ : Nat → RotationSystem
  vw0 : Nat → Nat
  vw1 : Nat → Nat
  vta : Nat → Nat
  vx : Nat → Nat
  ces : Nat → List (Nat × Nat)
  cρ : Nat → RotationSystem
  cl : Nat → Nat → Nat

namespace RG

open RotationSystem

variable (G : RG)

def m : Nat := G.sk.length
def Sk : List (Nat × Nat) := G.cap :: G.sk
def vflat (k : Nat) : List (Nat × Nat) := (List.range k).flatMap G.ves
def cflat (k : Nat) : List (Nat × Nat) := (List.range k).flatMap G.ces
def vIdx (v : Nat) : Nat := 4 * (G.m + 1) + 4 * (G.vflat v).length
def cIdx (j : Nat) : Nat := G.vIdx G.nV + 4 * (G.cflat j).length

/-- Target of the skeleton quarter-edge `l` after the first `k` corner pieces are hung. -/
def tgtV : Nat → Nat → Nat
  | 0, l => G.ρ₀.rot l
  | k + 1, l =>
    if l = G.vta k then G.vIdx k + G.vw1 k
    else if l = G.ρ₀.rot (G.vta k) then G.vIdx k + G.vw0 k
    else tgtV k l

/-- Code `t` (a skeleton quarter-edge) after the first `k` substitutions: on a substituted edge
it has become the child's exposed slot. -/
def res (k t : Nat) : Nat :=
  if 4 ≤ t ∧ t < 4 * (k + 1) then G.cIdx (t / 4 - 1) + G.cl (t / 4 - 1) (t % 4) else t

/-- Target of the skeleton quarter-edge `l` after all corners and `k` substitutions. -/
def skt (k l : Nat) : Nat := G.res k (G.tgtV G.nV l)

/-- Position of code `c` after `k` substitutions. -/
def pos (k c : Nat) : Nat := if c < 4 then c else c - 4 * k

def skE (k : Nat) : Nat × Nat := G.sk[k]?.getD (0, 0)

/-- Edge list after all corners and `k` substitutions. -/
def E (k : Nat) : List (Nat × Nat) := G.cap :: (G.sk.drop k ++ (G.vflat G.nV ++ G.cflat k))

theorem vflat_succ (k : Nat) : G.vflat (k + 1) = G.vflat k ++ G.ves k := by
  simp [vflat, List.range_succ]

theorem cflat_succ (k : Nat) : G.cflat (k + 1) = G.cflat k ++ G.ces k := by
  simp [cflat, List.range_succ]

theorem vIdx_succ (k : Nat) : G.vIdx (k + 1) = G.vIdx k + 4 * (G.ves k).length := by
  simp only [vIdx, vflat_succ, List.length_append]; omega

theorem cIdx_succ (k : Nat) : G.cIdx (k + 1) = G.cIdx k + 4 * (G.ces k).length := by
  simp only [cIdx, cflat_succ, List.length_append]; omega

theorem vIdx_mono {a b : Nat} (h : a ≤ b) : G.vIdx a ≤ G.vIdx b := by
  induction b with
  | zero => cases Nat.le_zero.1 h; exact le_rfl
  | succ b ih =>
    rcases Nat.lt_or_eq_of_le h with h | h
    · exact (ih (by omega)).trans (by rw [vIdx_succ]; omega)
    · rw [h]

theorem cIdx_mono {a b : Nat} (h : a ≤ b) : G.cIdx a ≤ G.cIdx b := by
  induction b with
  | zero => cases Nat.le_zero.1 h; exact le_rfl
  | succ b ih =>
    rcases Nat.lt_or_eq_of_le h with h | h
    · exact (ih (by omega)).trans (by rw [cIdx_succ]; omega)
    · rw [h]

theorem Sk_length : G.Sk.length = G.m + 1 := by simp [Sk, m]

theorem cIdx_ge (j : Nat) : 4 * (G.m + 1) ≤ G.cIdx j := by unfold cIdx vIdx; omega

theorem pos_zero (c : Nat) : pos 0 c = c := by unfold pos; split_ifs <;> omega

theorem res_zero (t : Nat) : G.res 0 t = t := by unfold res; rw [ite_eq_right (by omega)]

theorem sh1_pos {k c : Nat} (h : c < 4 ∨ 4 * (k + 2) ≤ c) : sh1 (pos k c) = pos (k + 1) c := by
  unfold sh1 pos; split_ifs <;> omega

theorem pre_pos {k c : Nat} (h : c < 4 ∨ 4 * (k + 2) ≤ c) :
    (if pos (k + 1) c < 4 then pos (k + 1) c else pos (k + 1) c + 4) = pos k c := by
  unfold pos; split_ifs <;> omega

theorem pos_div {k t : Nat} (_ht : t < 4 ∨ 4 * (k + 1) ≤ t) :
    pos k t / 4 = 1 ↔ 4 * (k + 1) ≤ t ∧ t < 4 * (k + 2) := by
  unfold pos; split_ifs <;> omega

theorem pick4_cl (k x : Nat) : pick4 (G.cl k 0) (G.cl k 1) (G.cl k 2) (G.cl k 3) x = G.cl k (x % 4) := by
  unfold pick4
  split_ifs with h0 h1 h2
  · rw [h0]
  · rw [h1]
  · rw [h2]
  · have : x % 4 = 3 := by omega
    rw [this]

theorem E_length {k : Nat} (_hk : k ≤ G.m) :
    (G.E k).length = 1 + ((G.m - k) + ((G.vflat G.nV).length + (G.cflat k).length)) := by
  simp [E, m, List.length_drop]; omega

theorem hasEdge_append {es es' : List (Nat × Nat)} {w : Nat} :
    HasEdge (es ++ es') w ↔ HasEdge es w ∨ HasEdge es' w := by
  unfold HasEdge
  constructor
  · rintro ⟨p, hp, hw⟩
    rcases List.mem_append.1 hp with h | h
    · exact Or.inl ⟨p, h, hw⟩
    · exact Or.inr ⟨p, h, hw⟩
  · rintro (⟨p, hp, hw⟩ | ⟨p, hp, hw⟩)
    · exact ⟨p, List.mem_append_left _ hp, hw⟩
    · exact ⟨p, List.mem_append_right _ hp, hw⟩

theorem hasEdge_vflat {k w : Nat} (h : HasEdge (G.vflat k) w) : ∃ v, v < k ∧ HasEdge (G.ves v) w := by
  obtain ⟨p, hp, hw⟩ := h
  obtain ⟨v, hv, hpv⟩ := List.mem_flatMap.1 hp
  exact ⟨v, List.mem_range.1 hv, p, hpv, hw⟩

theorem hasEdge_cflat {k w : Nat} (h : HasEdge (G.cflat k) w) : ∃ j, j < k ∧ HasEdge (G.ces j) w := by
  obtain ⟨p, hp, hw⟩ := h
  obtain ⟨j, hj, hpj⟩ := List.mem_flatMap.1 hp
  exact ⟨j, List.mem_range.1 hj, p, hpj, hw⟩

/-- Hypotheses of the gluing. -/
structure Hyp : Prop where
  planar₀ : IsPlanarEmbedding G.Sk G.n G.ρ₀
  v_planar : ∀ v, v < G.nV → IsPlanarEmbedding (G.ves v) G.n (G.vρ v)
  v_get : ∀ v, v < G.nV → (G.vρ v).get (G.vw0 v) = some (G.vw1 v)
  v_lt : ∀ v, v < G.nV → G.vw0 v < (G.vρ v).size
  v_even : ∀ v, v < G.nV → G.vw0 v % 2 = 0
  ta_ge : ∀ v, v < G.nV → 4 ≤ G.vta v
  ta_lt : ∀ v, v < G.nV → G.vta v < 4 * (G.m + 1)
  ta_even : ∀ v, v < G.nV → G.vta v % 2 = 0
  ta_vert : ∀ v, v < G.nV → QE.vert G.Sk (G.vta v) = some (G.vx v)
  w_vert : ∀ v, v < G.nV → QE.vert (G.ves v) (G.vw0 v) = some (G.vx v)
  x_inj : ∀ v v', v < G.nV → v' < G.nV → G.vx v = G.vx v' → v = v'
  x_cap : ∀ v, v < G.nV → G.vx v ≠ G.cap.1 ∧ G.vx v ≠ G.cap.2
  v_sep : ∀ v w, v < G.nV → HasEdge G.Sk w → HasEdge (G.ves v) w → w = G.vx v
  vv_sep : ∀ v v' w, v < G.nV → v' < G.nV → v ≠ v' → HasEdge (G.ves v) w →
    HasEdge (G.ves v') w → HasEdge G.Sk w
  c_planar : ∀ j, j < G.m → IsPlanarEmbedding (G.ces j) G.n (G.cρ j)
  c_get0 : ∀ j, j < G.m → (G.cρ j).get (G.cl j 0) = some (G.cl j 1)
  c_get2 : ∀ j, j < G.m → (G.cρ j).get (G.cl j 2) = some (G.cl j 3)
  c_even0 : ∀ j, j < G.m → G.cl j 0 % 2 = 0
  c_even2 : ∀ j, j < G.m → G.cl j 2 % 2 = 0
  c_lt0 : ∀ j, j < G.m → G.cl j 0 < (G.cρ j).size
  c_lt2 : ∀ j, j < G.m → G.cl j 2 < (G.cρ j).size
  c_vert : ∀ j p, G.sk[j]? = some p → QE.vert (G.ces j) (G.cl j 0) = some p.1 ∧
    QE.vert (G.ces j) (G.cl j 2) = some p.2 ∧ p.1 ≠ p.2
  c_face : ∀ j, j < G.m → SameOrbit ((G.cρ j).stepC 3) (G.cl j 0) (G.cl j 2)
  c_sep : ∀ j p w, G.sk[j]? = some p → HasEdge (G.ces j) w → HasEdge G.Sk w → w = p.1 ∨ w = p.2
  cc_sep : ∀ j j' w, j < G.m → j' < G.m → j ≠ j' → HasEdge (G.ces j) w → HasEdge (G.ces j') w →
    HasEdge G.Sk w
  cv_sep : ∀ j v w, j < G.m → v < G.nV → HasEdge (G.ces j) w → HasEdge (G.ves v) w → HasEdge G.Sk w
  deg₀ : ∀ l, l < 4 * (G.m + 1) → G.ρ₀.rot l / 4 ≠ l / 4
  cap_ne : G.cap.1 ≠ G.cap.2
  cap_conn : EdgesConn G.sk G.cap.1 G.cap.2

variable {G}

/-- Invariant after hanging the first `k` corner pieces. -/
structure InvV (G : RG) (k : Nat) (ρ : RotationSystem) : Prop where
  planar : IsPlanarEmbedding (G.Sk ++ G.vflat k) G.n ρ
  size : ρ.size = G.vIdx k
  skel : ∀ l, l < 4 * (G.m + 1) → ρ.rot l = G.tgtV k l
  piece : ∀ v, v < k → ∀ r, r < (G.vρ v).size → ρ.rot (G.vIdx v + r) =
    if r = G.vw0 v then G.ρ₀.rot (G.vta v) else if r = G.vw1 v then G.vta v
    else G.vIdx v + (G.vρ v).rot r

/-- Invariant after `k` substitutions. -/
structure InvC (G : RG) (k : Nat) (ρ : RotationSystem) : Prop where
  planar : IsPlanarEmbedding (G.E k) G.n ρ
  size : ρ.size = pos k (G.cIdx k)
  skel : ∀ l, (l < 4 ∨ 4 * (k + 1) ≤ l) → l < 4 * (G.m + 1) →
    ρ.rot (pos k l) = pos k (G.skt k l)
  piece : ∀ v, v < G.nV → ∀ r, r < (G.vρ v).size → ρ.rot (pos k (G.vIdx v + r)) =
    pos k (if r = G.vw0 v then G.res k (G.ρ₀.rot (G.vta v)) else if r = G.vw1 v then G.res k (G.vta v)
      else G.vIdx v + (G.vρ v).rot r)
  child : ∀ j, j < k → ∀ r, r < (G.cρ j).size → ρ.rot (pos k (G.cIdx j + r)) =
    pos k (if r = G.cl j 0 then G.skt k (4 * (j + 1)) else if r = G.cl j 1 then G.skt k (4 * (j + 1) + 1)
      else if r = G.cl j 2 then G.skt k (4 * (j + 1) + 2) else if r = G.cl j 3 then G.skt k (4 * (j + 1) + 3)
      else G.cIdx j + (G.cρ j).rot r)

/-- Partner of target `t` once the cap is removed: targets on the cap are rewired across it. -/
def fin (G : RG) (t : Nat) : Nat :=
  if t = 0 then G.skt G.m 1 else if t = 1 then G.skt G.m 0 else if t = 2 then G.skt G.m 3
  else if t = 3 then G.skt G.m 2 else t

namespace Hyp

variable (H : G.Hyp)
include H

theorem ρ₀_size : G.ρ₀.size = 4 * (G.m + 1) := by rw [H.planar₀.size, Sk_length]

theorem tb_lt {v : Nat} (hv : v < G.nV) : G.ρ₀.rot (G.vta v) < 4 * (G.m + 1) := by
  rw [← H.ρ₀_size]
  exact rot_lt H.planar₀.total H.planar₀.involution (by rw [H.ρ₀_size]; exact H.ta_lt v hv)

theorem tb_vert {v : Nat} (hv : v < G.nV) : QE.vert G.Sk (G.ρ₀.rot (G.vta v)) = some (G.vx v) := by
  rw [← H.ta_vert v hv]
  exact (H.planar₀.same_vertex _ (by rw [H.ρ₀_size]; exact H.ta_lt v hv) _
    (get_eq_rot H.planar₀.total (by rw [H.ρ₀_size]; exact H.ta_lt v hv))).symm

theorem tb_ne_ta {v : Nat} (hv : v < G.nV) : G.ρ₀.rot (G.vta v) ≠ G.vta v :=
  rot_ne H.planar₀.total H.planar₀.opposite_dir (by rw [H.ρ₀_size]; exact H.ta_lt v hv)

/-- Corners of distinct pieces are disjoint. -/
theorem corner_ne {v v' : Nat} (hv : v < G.nV) (hv' : v' < G.nV) (hne : v ≠ v') :
    G.vta v ≠ G.vta v' ∧ G.vta v ≠ G.ρ₀.rot (G.vta v') := by
  constructor
  · intro h
    apply hne
    apply H.x_inj v v' hv hv'
    have := H.ta_vert v hv; rw [h, H.ta_vert v' hv'] at this
    exact (Option.some.inj this).symm
  · intro h
    apply hne
    apply H.x_inj v v' hv hv'
    have := H.ta_vert v hv; rw [h, H.tb_vert hv'] at this
    exact (Option.some.inj this).symm

theorem tgtV_ta {k v : Nat} (hk : k ≤ v) (hv : v < G.nV) : G.tgtV k (G.vta v) = G.ρ₀.rot (G.vta v) := by
  induction k with
  | zero => rfl
  | succ k ih =>
    have hne := H.corner_ne hv (by omega) (by omega : v ≠ k)
    simp only [tgtV, ite_eq_right hne.1, ite_eq_right hne.2]
    exact ih (by omega)

theorem tgtV_lt {k l : Nat} (hk : k ≤ G.nV) (hl : l < 4 * (G.m + 1)) : G.tgtV k l < G.vIdx G.nV := by
  induction k with
  | zero =>
    show G.ρ₀.rot l < _
    have := rot_lt H.planar₀.total H.planar₀.involution (q := l) (by rw [H.ρ₀_size]; exact hl)
    rw [H.ρ₀_size] at this
    unfold vIdx; omega
  | succ k ih =>
    have hw1 : G.vw1 k < (G.vρ k).size := by
      rw [← rot_eq_of_get (H.v_get k (by omega))]
      exact rot_lt (H.v_planar k (by omega)).total (H.v_planar k (by omega)).involution
        (H.v_lt k (by omega))
    have hmono := G.vIdx_mono (a := k + 1) (b := G.nV) hk
    rw [vIdx_succ, (H.v_planar k (by omega)).size] at *
    simp only [tgtV]
    split_ifs
    · omega
    · have := H.v_lt k (by omega); rw [(H.v_planar k (by omega)).size] at this; omega
    · exact ih (by omega)

theorem tgtV_cap {k s : Nat} (hk : k ≤ G.nV) (hs : s < 4) : G.tgtV k s = G.ρ₀.rot s := by
  induction k with
  | zero => rfl
  | succ k ih =>
    have h1 : s ≠ G.vta k := by have := H.ta_ge k (by omega); omega
    have h2 : s ≠ G.ρ₀.rot (G.vta k) := by
      intro h
      have hv := H.tb_vert (v := k) (by omega)
      rw [← h, Sk, vert_cons_lt _ _ hs] at hv
      have hx := H.x_cap k (by omega)
      split_ifs at hv <;> [exact hx.1 (Option.some.inj hv).symm; exact hx.2 (Option.some.inj hv).symm]
    simp only [tgtV, ite_eq_right h1, ite_eq_right h2]
    exact ih (by omega)

theorem invV_zero : G.InvV 0 G.ρ₀ where
  planar := by simpa [vflat] using H.planar₀
  size := by rw [H.ρ₀_size]; simp [vIdx, vflat]
  skel := fun _ _ => rfl
  piece := fun _ h => absurd h (Nat.not_lt_zero _)

theorem invV_step {k : Nat} {ρ : RotationSystem} (hk : k < G.nV) (h : G.InvV k ρ) :
    G.InvV (k + 1) ((ρ.union (G.vρ k)).conj (G.vta k) (ρ.size + G.vw0 k)) := by
  have hvp := H.v_planar k hk
  have hta := H.ta_lt k hk
  have hta_lt : G.vta k < ρ.size := by rw [h.size]; unfold vIdx; omega
  have hva : QE.vert (G.Sk ++ G.vflat k) (G.vta k) = some (G.vx k) := by
    rw [vert_append_left (by rw [Sk_length]; omega)]; exact H.ta_vert k hk
  have hsep : ∀ w, HasEdge (G.Sk ++ G.vflat k) w → HasEdge (G.ves k) w → w = G.vx k := by
    intro w hw hwk
    rcases hasEdge_append.1 hw with hw | hw
    · exact H.v_sep k w hk hw hwk
    · obtain ⟨v, hv, hwv⟩ := G.hasEdge_vflat hw
      exact H.v_sep k w hk (H.vv_sep v k w (by omega) hk (by omega) hwv hwk) hwk
  obtain ⟨hpl, hsz, h1, h2⟩ := vstep h.planar hvp hta_lt (H.v_get k hk) (H.v_lt k hk)
    (by rw [H.ta_even k hk, H.v_even k hk]) hva (H.w_vert k hk) hsep
  have hrta : ρ.rot (G.vta k) = G.ρ₀.rot (G.vta k) := by
    rw [h.skel _ hta, H.tgtV_ta le_rfl hk]
  have htb := H.tb_lt hk
  refine ⟨?_, ?_, ?_, ?_⟩
  · rw [vflat_succ, ← List.append_assoc]; exact hpl
  · rw [hsz, h.size, vIdx_succ, hvp.size]
  · intro l hl
    rw [h1 l (by rw [h.size]; unfold vIdx; omega), hrta, h.size, h.skel l hl]
    rfl
  · intro v hv r hr
    rcases Nat.lt_or_ge v k with hvk | hvk
    · have hlt : G.vIdx v + r < ρ.size := by
        rw [h.size]
        have := G.vIdx_mono (a := v + 1) (b := k) hvk
        rw [vIdx_succ, ← (H.v_planar v (by omega)).size] at this; omega
      rw [h1 _ hlt, hrta, ite_eq_right (by unfold vIdx; omega), ite_eq_right (by unfold vIdx; omega)]
      exact h.piece v hvk r hr
    · have hvk' : v = k := by omega
      subst hvk'
      rw [← h.size, h2 r hr, hrta]

/-- The corner pieces hung one after the other. -/
def vsys (G : RG) : Nat → RotationSystem
  | 0 => G.ρ₀
  | k + 1 => ((vsys G k).union (G.vρ k)).conj (G.vta k) ((vsys G k).size + G.vw0 k)

theorem invV_all {k : Nat} (hk : k ≤ G.nV) : G.InvV k (vsys G k) := by
  induction k with
  | zero => exact H.invV_zero
  | succ k ih => exact H.invV_step (by omega) (ih (by omega))

theorem cl_lt {j s : Nat} (hj : j < G.m) (hs : s < 4) : G.cl j s < (G.cρ j).size := by
  have hc := H.c_planar j hj
  have h1 : G.cl j 1 < (G.cρ j).size := by
    rw [← rot_eq_of_get (H.c_get0 j hj)]; exact rot_lt hc.total hc.involution (H.c_lt0 j hj)
  have h3 : G.cl j 3 < (G.cρ j).size := by
    rw [← rot_eq_of_get (H.c_get2 j hj)]; exact rot_lt hc.total hc.involution (H.c_lt2 j hj)
  rcases s with _ | _ | _ | _ | s
  · exact H.c_lt0 j hj
  · exact h1
  · exact H.c_lt2 j hj
  · exact h3
  · omega

theorem code_lt {j k s : Nat} (hj : j < k) (hk : k ≤ G.m) (hs : s < 4) :
    G.cIdx j + G.cl j s < G.cIdx k := by
  have := H.cl_lt (by omega : j < G.m) hs
  rw [(H.c_planar j (by omega)).size] at this
  have hm := G.cIdx_mono (a := j + 1) (b := k) hj
  rw [cIdx_succ] at hm; omega

omit H in
theorem tgtV_cases {k l : Nat} (hk : k ≤ G.nV) : G.tgtV k l = G.ρ₀.rot l ∨ 4 * (G.m + 1) ≤ G.tgtV k l := by
  induction k with
  | zero => exact Or.inl rfl
  | succ k ih =>
    simp only [tgtV]
    split_ifs
    · right; unfold vIdx; omega
    · right; unfold vIdx; omega
    · exact ih (by omega)

omit H in
theorem res_succ {k : Nat} (hk : k < G.m) (x : Nat) :
    G.res (k + 1) x = if 4 * (k + 1) ≤ G.res k x ∧ G.res k x < 4 * (k + 2) then
      G.cIdx k + G.cl k (G.res k x % 4) else G.res k x := by
  by_cases hA : 4 ≤ x ∧ x < 4 * (k + 1)
  · have hge := G.cIdx_ge (x / 4 - 1)
    simp only [res, ite_eq_left hA, ite_eq_left (show 4 ≤ x ∧ x < 4 * (k + 1 + 1) by omega)]
    rw [ite_eq_right (by omega)]
  · by_cases hB : 4 ≤ x ∧ x < 4 * (k + 1 + 1)
    · have hx : x / 4 - 1 = k := by omega
      simp only [res, ite_eq_right hA, ite_eq_left hB, hx]
      rw [ite_eq_left (by omega)]
    · simp only [res, ite_eq_right hA, ite_eq_right hB]
      rw [ite_eq_right (by omega)]

omit H in
theorem skt_succ {k : Nat} (hk : k < G.m) (l : Nat) :
    G.skt (k + 1) l = if 4 * (k + 1) ≤ G.skt k l ∧ G.skt k l < 4 * (k + 2) then
      G.cIdx k + G.cl k (G.skt k l % 4) else G.skt k l :=
  res_succ hk _

theorem res_live {k x : Nat} (hk : k ≤ G.m) (hx : x < 4 * (G.m + 1)) :
    (G.res k x < 4 ∨ 4 * (k + 1) ≤ G.res k x) ∧ G.res k x < G.cIdx k ∧
    (G.res k x = x ∨ 4 * (G.m + 1) ≤ G.res k x) := by
  unfold res
  split_ifs with hA
  · have := H.code_lt (j := x / 4 - 1) (k := k) (s := x % 4) (by omega) hk (Nat.mod_lt _ (by omega))
    have := G.cIdx_ge (x / 4 - 1)
    omega
  · have := G.cIdx_ge k
    omega

theorem skt_live {k l : Nat} (hk : k ≤ G.m) (hl : l < 4 * (G.m + 1)) :
    (G.skt k l < 4 ∨ 4 * (k + 1) ≤ G.skt k l) ∧ G.skt k l < G.cIdx k ∧ G.skt k l / 4 ≠ l / 4 := by
  unfold skt
  rcases tgtV_cases (G := G) (k := G.nV) (l := l) le_rfl with ht | ht
  · rw [ht]
    have hr := rot_lt H.planar₀.total H.planar₀.involution (q := l) (by rw [H.ρ₀_size]; exact hl)
    rw [H.ρ₀_size] at hr
    obtain ⟨h1, h2, h3⟩ := H.res_live hk hr
    refine ⟨h1, h2, ?_⟩
    rcases h3 with h3 | h3
    · rw [h3]; exact H.deg₀ l hl
    · omega
  · have hlt := H.tgtV_lt (k := G.nV) le_rfl hl
    have hres : G.res k (G.tgtV G.nV l) = G.tgtV G.nV l := by
      unfold res; rw [ite_eq_right (by omega)]
    rw [hres]
    have := G.cIdx_mono (a := 0) (b := k) (Nat.zero_le _)
    have h0 : G.cIdx 0 = G.vIdx G.nV := by simp [cIdx, cflat]
    omega

theorem invC_zero : G.InvC 0 (vsys G G.nV) := by
  have h := H.invV_all (k := G.nV) le_rfl
  refine ⟨?_, ?_, ?_, ?_, ?_⟩
  · have := h.planar
    simpa [E, Sk, cflat] using this
  · rw [h.size, pos_zero]; simp [cIdx, cflat]
  · intro l _ hl
    rw [pos_zero, pos_zero, h.skel l hl, skt, res_zero]
  · intro v hv r hr
    rw [pos_zero, pos_zero, res_zero, res_zero]
    exact h.piece v hv r hr
  · intro j hj; omega

/-- The substitutions one after the other. -/
def csys (G : RG) : Nat → RotationSystem
  | 0 => vsys G G.nV
  | k + 1 => (substSum (G.E k) (G.skE k :: G.ces k) G.n 1 (G.skE k).1 (G.skE k).2).splice
      (csys G k) ((G.cρ k).insert 0 (G.cl k 0) (G.cl k 2))

theorem invC_step {k : Nat} (hk : k < G.m) (h : G.InvC k (csys G k)) : G.InvC (k + 1) (csys G (k + 1)) := by
  have hk' : k < G.sk.length := hk
  obtain ⟨⟨u, v⟩, hp⟩ : ∃ p, G.sk[k]? = some p := ⟨_, List.getElem?_eq_getElem hk'⟩
  have hskE : G.skE k = (u, v) := by simp [skE, hp]
  have hdrop : G.sk.drop k = (u, v) :: G.sk.drop (k + 1) := by
    rw [List.drop_eq_getElem_cons hk']
    rw [List.getElem?_eq_getElem hk'] at hp
    rw [Option.some.inj hp]
  have he : (G.E k)[1]? = some (u, v) := by simp [E, hdrop]
  have hE1 : (G.E k).eraseIdx 1 = G.cap :: (G.sk.drop (k + 1) ++ (G.vflat G.nV ++ G.cflat k)) := by
    simp [E, hdrop]
  have hE2 : (G.E k).eraseIdx 1 ++ G.ces k = G.E (k + 1) := by
    rw [hE1, E, cflat_succ]; simp [List.append_assoc]
  have hlen := G.E_length (k := k) (by omega)
  have hoff : 4 * ((G.E k).length - 1) = pos (k + 1) (G.cIdx k) := by
    have := G.cIdx_ge k
    unfold pos cIdx vIdx; rw [hlen]; split_ifs <;> omega
  obtain ⟨hvu, hvv, huv⟩ := H.c_vert k (u, v) hp
  set ρ := csys G k with hρ
  have hdeg : ∀ s, s < 4 → ρ.rot (4 + s) / 4 ≠ 1 := by
    intro s hs
    have hl : 4 * (k + 1) + s < 4 * (G.m + 1) := by omega
    have := h.skel (4 * (k + 1) + s) (Or.inr (by omega)) hl
    have hpos : pos k (4 * (k + 1) + s) = 4 + s := by unfold pos; split_ifs <;> omega
    rw [hpos] at this
    rw [this]
    obtain ⟨h1, h2, h3⟩ := H.skt_live (k := k) (l := 4 * (k + 1) + s) (by omega) hl
    rw [Ne, pos_div h1]
    omega
  have hsep : ∀ w, HasEdge ((G.E k).eraseIdx 1) w → HasEdge (G.ces k) w → w = u ∨ w = v := by
    intro w hw hwc
    rw [hE1] at hw
    have hSk : HasEdge G.Sk w → w = u ∨ w = v := fun h => H.c_sep k (u, v) w hp hwc h
    obtain ⟨q, hq, hqw⟩ := hw
    rcases List.mem_cons.1 hq with hq | hq
    · exact hSk ⟨q, by rw [hq]; exact List.mem_cons_self, hqw⟩
    rcases List.mem_append.1 hq with hq | hq
    · exact hSk ⟨q, List.mem_cons_of_mem _ (List.mem_of_mem_drop hq), hqw⟩
    rcases List.mem_append.1 hq with hq | hq
    · obtain ⟨v', hv', hwv⟩ := G.hasEdge_vflat ⟨q, hq, hqw⟩
      exact hSk (H.cv_sep k v' w hk hv' hwc hwv)
    · obtain ⟨j, hj, hwj⟩ := G.hasEdge_cflat ⟨q, hq, hqw⟩
      exact hSk (H.cc_sep k j w hk (by omega) (by omega) hwc hwj)
  obtain ⟨hpl, hsz, h1, h2⟩ := subst_step h.planar he huv hdeg (H.c_planar k hk) (H.c_get0 k hk)
    (H.c_get2 k hk) (H.c_even0 k hk) (H.c_even2 k hk) (H.c_lt0 k hk) (H.c_lt2 k hk) hvu hvv
    (H.c_face k hk) hsep
  have hσ : csys G (k + 1) = (substSum (G.E k) ((u, v) :: G.ces k) G.n 1 u v).splice ρ
      ((G.cρ k).insert 0 (G.cl k 0) (G.cl k 2)) := by
    simp only [csys, hskE, hρ]
  rw [hσ]
  set σ := (substSum (G.E k) ((u, v) :: G.ces k) G.n 1 u v).splice ρ
      ((G.cρ k).insert 0 (G.cl k 0) (G.cl k 2)) with hσdef
  have hcge := G.cIdx_ge k
  have step : ∀ c t, (c < 4 ∨ 4 * (k + 2) ≤ c) → c < G.cIdx k → (t < 4 ∨ 4 * (k + 1) ≤ t) →
      t < G.cIdx k → ρ.rot (pos k c) = pos k t →
      σ.rot (pos (k + 1) c) =
        pos (k + 1) (if 4 * (k + 1) ≤ t ∧ t < 4 * (k + 2) then G.cIdx k + G.cl k (t % 4) else t) := by
    intro c t hc hclt ht htlt hρt
    rw [h1 _ (by rw [hoff]; unfold pos; split_ifs <;> omega), pre_pos hc, hρt]
    by_cases hb : 4 * (k + 1) ≤ t ∧ t < 4 * (k + 2)
    · rw [ite_eq_left ((pos_div ht).2 hb), ite_eq_left hb, pick4_cl, hoff]
      have : pos k t % 4 = t % 4 := by unfold pos; split_ifs <;> omega
      rw [this]
      unfold pos; split_ifs <;> omega
    · rw [ite_eq_right (fun h => hb ((pos_div ht).1 h)), ite_eq_right hb]
      exact sh1_pos (by omega)
  refine ⟨?_, ?_, ?_, ?_, ?_⟩
  · rw [← hE2]; exact hpl
  · rw [hsz, hoff, (H.c_planar k hk).size, cIdx_succ]
    unfold pos; split_ifs <;> omega
  · intro l hl hlm
    obtain ⟨ht, htlt, -⟩ := H.skt_live (k := k) (l := l) (by omega) hlm
    rw [step l (G.skt k l) hl (by omega) ht htlt (h.skel l (by omega) hlm), skt_succ hk]
  · intro v' hv' r hr
    have hvlt : G.vIdx v' + r < G.cIdx k := by
      have := G.vIdx_mono (a := v' + 1) (b := G.nV) hv'
      rw [vIdx_succ, ← (H.v_planar v' hv').size] at this
      have := G.cIdx_mono (a := 0) (b := k) (Nat.zero_le _)
      have h0 : G.cIdx 0 = G.vIdx G.nV := by simp [cIdx, cflat]
      omega
    have hvge : 4 * (G.m + 1) ≤ G.vIdx v' := by unfold vIdx; omega
    have hta := H.ta_lt v' hv'
    have htb := H.tb_lt hv'
    have hrot : G.vIdx v' + (G.vρ v').rot r < G.cIdx k := by
      have := G.vIdx_mono (a := v' + 1) (b := G.nV) hv'
      rw [vIdx_succ, ← (H.v_planar v' hv').size] at this
      have := rot_lt (H.v_planar v' hv').total (H.v_planar v' hv').involution hr
      have := G.cIdx_mono (a := 0) (b := k) (Nat.zero_le _)
      have h0 : G.cIdx 0 = G.vIdx G.nV := by simp [cIdx, cflat]
      omega
    split_ifs with e0 e1
    · obtain ⟨ht, htlt, -⟩ := H.res_live (k := k) (x := G.ρ₀.rot (G.vta v')) (by omega) htb
      rw [step _ _ (Or.inr (by omega)) hvlt ht htlt (by rw [h.piece v' hv' r hr, ite_eq_left e0]), res_succ hk]
    · obtain ⟨ht, htlt, -⟩ := H.res_live (k := k) (x := G.vta v') (by omega) hta
      rw [step _ _ (Or.inr (by omega)) hvlt ht htlt
        (by rw [h.piece v' hv' r hr, ite_eq_right e0, ite_eq_left e1]), res_succ hk]
    · rw [step _ _ (Or.inr (by omega)) hvlt (Or.inr (by omega)) hrot
        (by rw [h.piece v' hv' r hr, ite_eq_right e0, ite_eq_right e1]), ite_eq_right (by omega)]
  · intro j hj r hr
    rcases Nat.lt_or_ge j k with hjk | hjk
    · have hjlt : G.cIdx j + r < G.cIdx k := by
        have := G.cIdx_mono (a := j + 1) (b := k) hjk
        rw [cIdx_succ, ← (H.c_planar j (by omega)).size] at this; omega
      have hjge := G.cIdx_ge j
      have hrot : G.cIdx j + (G.cρ j).rot r < G.cIdx k := by
        have := G.cIdx_mono (a := j + 1) (b := k) hjk
        rw [cIdx_succ, ← (H.c_planar j (by omega)).size] at this
        have := rot_lt (H.c_planar j (by omega)).total (H.c_planar j (by omega)).involution hr
        omega
      have hc := h.child j hjk r hr
      have key : ∀ s, s < 4 → r = G.cl j s →
          ρ.rot (pos k (G.cIdx j + r)) = pos k (G.skt k (4 * (j + 1) + s)) → 
          σ.rot (pos (k + 1) (G.cIdx j + r)) = pos (k + 1) (G.skt (k + 1) (4 * (j + 1) + s)) := by
        intro s hs _ hρs
        obtain ⟨ht, htlt, -⟩ := H.skt_live (k := k) (l := 4 * (j + 1) + s) (by omega) (by omega)
        rw [step _ _ (Or.inr (by omega)) hjlt ht htlt hρs, skt_succ hk]
      split_ifs with e0 e1 e2 e3
      · exact key 0 (by omega) e0 (by simp only [Nat.add_zero]; rw [hc, ite_eq_left e0])
      · exact key 1 (by omega) e1 (by rw [hc, ite_eq_right e0, ite_eq_left e1])
      · exact key 2 (by omega) e2 (by rw [hc, ite_eq_right e0, ite_eq_right e1, ite_eq_left e2])
      · exact key 3 (by omega) e3 (by
          rw [hc, ite_eq_right e0, ite_eq_right e1, ite_eq_right e2, ite_eq_left e3])
      · rw [step _ _ (Or.inr (by omega)) hjlt (Or.inr (by omega)) hrot
          (by rw [hc, ite_eq_right e0, ite_eq_right e1, ite_eq_right e2, ite_eq_right e3]),
          ite_eq_right (by omega)]
    · have hjk' : j = k := by omega
      subst hjk'
      have hpos : pos (j + 1) (G.cIdx j + r) = 4 * ((G.E j).length - 1) + r := by
        rw [hoff]; unfold pos; split_ifs <;> omega
      rw [hpos, h2 r hr]
      have key : ∀ s, s < 4 →
          sh1 (ρ.rot (4 + s)) = pos (j + 1) (G.skt (j + 1) (4 * (j + 1) + s)) := by
        intro s hs
        have hl : 4 * (j + 1) + s < 4 * (G.m + 1) := by omega
        have hpos' : pos j (4 * (j + 1) + s) = 4 + s := by unfold pos; split_ifs <;> omega
        have := h.skel (4 * (j + 1) + s) (Or.inr (by omega)) hl
        rw [hpos'] at this
        obtain ⟨ht, htlt, h3⟩ := H.skt_live (k := j) (l := 4 * (j + 1) + s) (by omega) hl
        have hnot : ¬(4 * (j + 1) ≤ G.skt j (4 * (j + 1) + s) ∧ G.skt j (4 * (j + 1) + s) < 4 * (j + 2)) := by
          omega
        rw [this, skt_succ hk, ite_eq_right hnot]
        exact sh1_pos (by omega)
      split_ifs with e0 e1 e2 e3
      · exact key 0 (by omega)
      · exact key 1 (by omega)
      · exact key 2 (by omega)
      · exact key 3 (by omega)
      · rw [hoff]; unfold pos; split_ifs <;> omega

theorem invC_all {k : Nat} (hk : k ≤ G.m) : G.InvC k (csys G k) := by
  induction k with
  | zero => exact H.invC_zero
  | succ k ih => exact H.invC_step (k := k) (by omega) (ih (by omega))

theorem cap_conn_children : EdgesConn (G.cflat G.m) G.cap.1 G.cap.2 := by
  have hstep : ∀ x y, (x, y) ∈ G.sk → EdgesConn (G.cflat G.m) x y := by
    intro x y hxy
    obtain ⟨j, hj⟩ := List.getElem?_of_mem hxy
    have hjm : j < G.m := (List.getElem?_eq_some_iff.1 hj).1
    obtain ⟨hvu, hvv, -⟩ := H.c_vert j (x, y) hj
    have hc := H.c_planar j hjm
    have := edgesConn_of_sameOrbit hc.total hc.involution hc.same_vertex hc.size (H.c_face j hjm) hvu hvv
    exact edgesConn_mono (fun p hp => List.mem_flatMap.2 ⟨j, List.mem_range.2 hjm, hp⟩) this
  have h := H.cap_conn
  generalize G.cap.2 = b at h ⊢
  induction h with
  | refl => exact Relation.ReflTransGen.refl
  | tail _ hbc ih =>
    refine Relation.ReflTransGen.trans ih ?_
    rcases hbc with h | h
    · exact hstep _ _ h
    · exact edgesConn_symm (hstep _ _ h)

theorem skt_cap {s : Nat} (hs : s < 4) : 4 * (G.m + 1) ≤ G.skt G.m s ∧ G.skt G.m s < G.cIdx G.m := by
  unfold skt
  rw [H.tgtV_cap le_rfl hs]
  have hr := rot_lt H.planar₀.total H.planar₀.involution (q := s) (by rw [H.ρ₀_size]; omega)
  rw [H.ρ₀_size] at hr
  have hd := H.deg₀ s (by omega)
  have h4 : 4 ≤ G.ρ₀.rot s := by omega
  unfold res
  rw [ite_eq_left ⟨h4, hr⟩]
  exact ⟨by have := G.cIdx_ge (G.ρ₀.rot s / 4 - 1); omega, H.code_lt (by omega) le_rfl (Nat.mod_lt _ (by omega))⟩

theorem glue : ∃ σ : RotationSystem,
    IsPlanarEmbedding (G.vflat G.nV ++ G.cflat G.m) G.n σ ∧
    σ.size + 4 * (G.m + 1) = G.cIdx G.m ∧
    (∀ v, v < G.nV → ∀ r, r < (G.vρ v).size → σ.rot (G.vIdx v + r - 4 * (G.m + 1)) =
      G.fin (if r = G.vw0 v then G.res G.m (G.ρ₀.rot (G.vta v)) else if r = G.vw1 v then G.res G.m (G.vta v)
        else G.vIdx v + (G.vρ v).rot r) - 4 * (G.m + 1)) ∧
    (∀ j, j < G.m → ∀ r, r < (G.cρ j).size → σ.rot (G.cIdx j + r - 4 * (G.m + 1)) =
      G.fin (if r = G.cl j 0 then G.skt G.m (4 * (j + 1)) else if r = G.cl j 1 then G.skt G.m (4 * (j + 1) + 1)
        else if r = G.cl j 2 then G.skt G.m (4 * (j + 1) + 2) else if r = G.cl j 3 then G.skt G.m (4 * (j + 1) + 3)
        else G.cIdx j + (G.cρ j).rot r) - 4 * (G.m + 1)) ∧
    σ.rot (G.skt G.m 0 - 4 * (G.m + 1)) = G.skt G.m 1 - 4 * (G.m + 1) ∧
    σ.rot (G.skt G.m 2 - 4 * (G.m + 1)) = G.skt G.m 3 - 4 * (G.m + 1) ∧
    SameOrbit (σ.stepC 3) (G.skt G.m 1 - 4 * (G.m + 1)) (G.skt G.m 3 - 4 * (G.m + 1)) := by
  have h := H.invC_all (k := G.m) le_rfl
  set rs := csys G G.m with hrs
  have hE : G.E G.m = G.cap :: (G.vflat G.nV ++ G.cflat G.m) := by
    simp [E, m, List.drop_length]
  have hpl : IsPlanarEmbedding (G.cap :: (G.vflat G.nV ++ G.cflat G.m)) G.n rs := hE ▸ h.planar
  have hconn : EdgesConn (G.vflat G.nV ++ G.cflat G.m) G.cap.1 G.cap.2 :=
    edgesConn_mono (fun p hp => List.mem_append_right _ hp) H.cap_conn_children
  obtain ⟨-, -, -, -, -, σ, hσ, hsz, hget, hface⟩ := hpl.uninsert H.cap_ne hconn
  have hM := G.cIdx_ge G.m
  have hrot : ∀ s, s < 4 → rs.rot s = G.skt G.m s - 4 * G.m := by
    intro s hs
    have := h.skel s (Or.inl hs) (by omega)
    rw [show pos G.m s = s by unfold pos; rw [ite_eq_left hs]] at this
    rw [this]; unfold pos; rw [ite_eq_right (by have := (H.skt_cap hs).1; omega)]
  have hsize : rs.size = G.cIdx G.m - 4 * G.m := by
    rw [h.size]; unfold pos; rw [ite_eq_right (by omega)]
  have hfin : ∀ c t, 4 * (G.m + 1) ≤ c → c < G.cIdx G.m → (t < 4 ∨ 4 * (G.m + 1) ≤ t) →
      t < G.cIdx G.m → rs.rot (pos G.m c) = pos G.m t →
      σ.rot (c - 4 * (G.m + 1)) = G.fin t - 4 * (G.m + 1) := by
    intro c t hc hclt ht htlt hρ
    have hq : c - 4 * (G.m + 1) < σ.size := by omega
    have hpc : pos G.m c = c - 4 * G.m := by unfold pos; rw [ite_eq_right (by omega)]
    rw [hpc] at hρ
    have h4q : 4 + (c - 4 * (G.m + 1)) = c - 4 * G.m := by omega
    apply rot_eq_of_get
    rw [hget _ hq, h4q]
    congr 1
    have hiff : ∀ s, s < 4 → (c - 4 * G.m = rs.rot s ↔ t = s) := by
      intro s hs
      constructor
      · intro he
        have := rot_rot hpl.total hpl.involution (q := s) (by rw [hsize]; omega)
        rw [← he, hρ] at this
        unfold pos at this; split_ifs at this <;> omega
      · intro he
        subst he
        have := rot_rot hpl.total hpl.involution (q := c - 4 * G.m) (by rw [hsize]; omega)
        rw [hρ, show pos G.m t = t by unfold pos; rw [ite_eq_left hs]] at this
        exact this.symm
    have hneg : ∀ s, s < 4 → t ≠ s → ¬(c - 4 * G.m = rs.rot s) := fun s hs hne he => hne ((hiff s hs).1 he)
    unfold fin
    rcases ht with ht | ht
    · obtain rfl | rfl | rfl | rfl : t = 0 ∨ t = 1 ∨ t = 2 ∨ t = 3 := by omega
      · rw [ite_eq_left ((hiff 0 (by omega)).2 rfl), hrot 1 (by omega), ite_eq_left rfl]; omega
      · rw [ite_eq_right (hneg 0 (by omega) (by omega)), ite_eq_left ((hiff 1 (by omega)).2 rfl),
          hrot 0 (by omega), ite_eq_right (by omega : (1 : Nat) ≠ 0), ite_eq_left rfl]
        omega
      · rw [ite_eq_right (hneg 0 (by omega) (by omega)), ite_eq_right (hneg 1 (by omega) (by omega)),
          ite_eq_left ((hiff 2 (by omega)).2 rfl), hrot 3 (by omega),
          ite_eq_right (by omega : (2 : Nat) ≠ 0), ite_eq_right (by omega : (2 : Nat) ≠ 1), ite_eq_left rfl]
        omega
      · rw [ite_eq_right (hneg 0 (by omega) (by omega)), ite_eq_right (hneg 1 (by omega) (by omega)),
          ite_eq_right (hneg 2 (by omega) (by omega)), ite_eq_left ((hiff 3 (by omega)).2 rfl),
          hrot 2 (by omega), ite_eq_right (by omega : (3 : Nat) ≠ 0), ite_eq_right (by omega : (3 : Nat) ≠ 1),
          ite_eq_right (by omega : (3 : Nat) ≠ 2), ite_eq_left rfl]
        omega
    · rw [ite_eq_right (hneg 0 (by omega) (by omega)), ite_eq_right (hneg 1 (by omega) (by omega)),
        ite_eq_right (hneg 2 (by omega) (by omega)), ite_eq_right (hneg 3 (by omega) (by omega)), hρ,
        ite_eq_right (by omega), ite_eq_right (by omega), ite_eq_right (by omega), ite_eq_right (by omega)]
      unfold pos; rw [ite_eq_right (by omega)]; omega
  have h0 : G.cIdx 0 = G.vIdx G.nV := by simp [cIdx, cflat]
  have hcap : ∀ s, s < 4 → rs.rot (pos G.m (G.skt G.m s)) = pos G.m s := by
    intro s hs
    have := rot_rot hpl.total hpl.involution (q := s) (by rw [hsize]; omega)
    rw [hrot s hs] at this
    rw [show pos G.m (G.skt G.m s) = G.skt G.m s - 4 * G.m by
      unfold pos; rw [ite_eq_right (by have := (H.skt_cap hs).1; omega)], this]
    unfold pos; rw [ite_eq_left hs]
  refine ⟨σ, hσ, by omega, ?_, ?_, ?_, ?_, ?_⟩
  · intro v hv r hr
    have hvge : 4 * (G.m + 1) ≤ G.vIdx v + r := by unfold vIdx; omega
    have hvlt : G.vIdx v + r < G.cIdx G.m := by
      have := G.vIdx_mono (a := v + 1) (b := G.nV) hv
      rw [vIdx_succ, ← (H.v_planar v hv).size] at this
      have := G.cIdx_mono (a := 0) (b := G.m) (Nat.zero_le _)
      omega
    have hrot' : G.vIdx v + (G.vρ v).rot r < G.cIdx G.m := by
      have := G.vIdx_mono (a := v + 1) (b := G.nV) hv
      rw [vIdx_succ, ← (H.v_planar v hv).size] at this
      have := rot_lt (H.v_planar v hv).total (H.v_planar v hv).involution hr
      have := G.cIdx_mono (a := 0) (b := G.m) (Nat.zero_le _)
      omega
    apply hfin _ _ hvge hvlt _ _ (h.piece v hv r hr)
    · split_ifs
      · exact (H.res_live le_rfl (H.tb_lt hv)).1
      · exact (H.res_live le_rfl (H.ta_lt v hv)).1
      · right; unfold vIdx; omega
    · split_ifs
      · exact (H.res_live le_rfl (H.tb_lt hv)).2.1
      · exact (H.res_live le_rfl (H.ta_lt v hv)).2.1
      · exact hrot'
  · intro j hj r hr
    have hjge : 4 * (G.m + 1) ≤ G.cIdx j + r := by have := G.cIdx_ge j; omega
    have hjlt : G.cIdx j + r < G.cIdx G.m := by
      have := G.cIdx_mono (a := j + 1) (b := G.m) hj
      rw [cIdx_succ, ← (H.c_planar j hj).size] at this; omega
    have hrot' : G.cIdx j + (G.cρ j).rot r < G.cIdx G.m := by
      have := G.cIdx_mono (a := j + 1) (b := G.m) hj
      rw [cIdx_succ, ← (H.c_planar j hj).size] at this
      have := rot_lt (H.c_planar j hj).total (H.c_planar j hj).involution hr
      omega
    apply hfin _ _ hjge hjlt _ _ (h.child j hj r hr)
    · split_ifs
      · exact (H.skt_live le_rfl (by omega)).1
      · exact (H.skt_live le_rfl (by omega)).1
      · exact (H.skt_live le_rfl (by omega)).1
      · exact (H.skt_live le_rfl (by omega)).1
      · right; have := G.cIdx_ge j; omega
    · split_ifs
      · exact (H.skt_live le_rfl (by omega)).2.1
      · exact (H.skt_live le_rfl (by omega)).2.1
      · exact (H.skt_live le_rfl (by omega)).2.1
      · exact (H.skt_live le_rfl (by omega)).2.1
      · exact hrot'
  · have := hfin (G.skt G.m 0) 0 (H.skt_cap (by omega)).1 (H.skt_cap (by omega)).2 (Or.inl (by omega))
      (by omega) (hcap 0 (by omega))
    rw [this]; rfl
  · have := hfin (G.skt G.m 2) 2 (H.skt_cap (by omega)).1 (H.skt_cap (by omega)).2 (Or.inl (by omega))
      (by omega) (hcap 2 (by omega))
    rw [this]; rfl
  · rw [hrot 1 (by omega), hrot 3 (by omega)] at hface
    have e1 : G.skt G.m 1 - 4 * G.m - 4 = G.skt G.m 1 - 4 * (G.m + 1) := by omega
    have e3 : G.skt G.m 3 - 4 * G.m - 4 = G.skt G.m 3 - 4 * (G.m + 1) := by omega
    rw [e1, e3] at hface
    exact hface

end Hyp

end RG

end Spqr
