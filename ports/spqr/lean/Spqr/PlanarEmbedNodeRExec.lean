import Spqr.PlanarEmbedNodeR

/-!
# The `R` node: pointwise description of the node fold

The node loop visits the node's quarter-edges `4 * neSt + l`, `l < 4 * nE`, in order and processes
each rotation pair `{l, T l}` once, at its smaller member: a cap quarter-edge records the facing
child slot in the node's row, a corner (`l % 4 = 2`, `T l % 4 = 1`, `V` item present) links the two
child slots through the `V` item's boundary pair, every other pair links the two child slots. Since
distinct pairs touch distinct slots, the final `rotAdj` is `partner` on every child slot that does
not face the cap, the `V` pairs are linked back, and everything else is unchanged (`nodeFoldR`).
-/

namespace Spqr.PlanarSpqrTree

open EmbedM

variable (t : PlanarSpqrTree)

/-- The `V` item read by the node step at a corner quarter-edge `ta`. -/
def vOf (ta : Nat) : Nat := (t.nodeVerts[(t.nodeEdges[QE.edge ta]!).nvs.2]!).vert

section RFold

variable (i neSt nE : Nat) (s : EmbedState) (T Q : Nat → Nat)

/-- `l` is a corner of the node step where a `V` item is attached. -/
def Corner (l : Nat) : Prop :=
  4 ≤ l ∧ l < 4 * nE ∧ l % 4 = 2 ∧ T l % 4 = 1 ∧ (s.outerE[t.vOf (4 * neSt + l)]![0]!).isSome = true

instance (l : Nat) : Decidable (t.Corner neSt nE s T l) := by unfold Corner; infer_instance

/-- The two boundary slots of the `V` item attached at corner `l`. -/
def w0 (l : Nat) : Nat := (s.outerE[t.vOf (4 * neSt + l)]![0]!).getD 0
def w1 (l : Nat) : Nat := (s.outerE[t.vOf (4 * neSt + l)]![1]!).getD 0

/-- Final partner of the child slot `Q l` of a non-cap quarter-edge `l` whose rotation `T l` is not
on the cap. -/
def partner (l : Nat) : Nat :=
  if t.Corner neSt nE s T l then t.w1 neSt s l
  else if t.Corner neSt nE s T (T l) then t.w0 neSt s (T l)
  else Q (T l)

/-- Data of the `R` node fold: `T` is the node's own rotation (`neRotAdj` relative to `4 * neSt`),
`Q l` the child slot that quarter-edge `l` reads through `treeQe`, and the `V` items at corners have
two distinct boundary slots, all distinct from each other and from the child slots. -/
structure RFold : Prop where
  rot : ∀ l, l < 4 * nE → t.neRotAdj[4 * neSt + l]! = some (4 * neSt + T l)
  T_lt : ∀ l, l < 4 * nE → T l < 4 * nE
  T_inv : ∀ l, l < 4 * nE → T (T l) = l
  T_edge : ∀ l, l < 4 * nE → T l / 4 ≠ l / 4
  qe : ∀ l, 4 ≤ l → l < 4 * nE → ∀ s' : EmbedState,
    (∀ j, j ≠ i → s'.outerE[j]? = s.outerE[j]?) → t.treeQe s' (4 * neSt + l) = some (Q l)
  Q_lt : ∀ l, 4 ≤ l → l < 4 * nE → Q l < s.rotAdj.size
  Q_inj : ∀ l l', 4 ≤ l → l < 4 * nE → 4 ≤ l' → l' < 4 * nE → Q l = Q l' → l = l'
  v_ne : ∀ l, 4 ≤ l → l < 4 * nE → t.vOf (4 * neSt + l) ≠ i
  v_lt : ∀ l, t.Corner neSt nE s T l → l < T l
  v_some1 : ∀ l, t.Corner neSt nE s T l → (s.outerE[t.vOf (4 * neSt + l)]![1]!).isSome = true
  w_ne : ∀ l, t.Corner neSt nE s T l → t.w0 neSt s l ≠ t.w1 neSt s l
  w_lt : ∀ l, t.Corner neSt nE s T l → t.w0 neSt s l < s.rotAdj.size ∧ t.w1 neSt s l < s.rotAdj.size
  wQ : ∀ l, t.Corner neSt nE s T l → ∀ l', 4 ≤ l' → l' < 4 * nE →
    Q l' ≠ t.w0 neSt s l ∧ Q l' ≠ t.w1 neSt s l
  ww : ∀ l l', t.Corner neSt nE s T l → t.Corner neSt nE s T l' → l ≠ l' →
    t.w0 neSt s l ≠ t.w0 neSt s l' ∧ t.w0 neSt s l ≠ t.w1 neSt s l' ∧
    t.w1 neSt s l ≠ t.w0 neSt s l' ∧ t.w1 neSt s l ≠ t.w1 neSt s l'
  row : ∃ r : Array (Option Nat), s.outerE[i]? = some r ∧ r.size = 4
  nE_pos : 0 < nE

/-- The state after the first `x` quarter-edges of the node loop. -/
def rFold (x : Nat) : EmbedState := (List.range' (4 * neSt) x).foldl (t.nodeStep i neSt) s

/-- Invariant after `x` steps: the pairs `{l, T l}` with `min l (T l) < x` are processed. -/
structure InvR (x : Nat) (sx : EmbedState) : Prop where
  rsize : sx.rotAdj.size = s.rotAdj.size
  rows : ∀ j, j ≠ i → sx.outerE[j]? = s.outerE[j]?
  slot : ∀ l, 4 ≤ l → l < 4 * nE → 4 ≤ T l →
    sx.rotAdj[Q l]? =
      if min l (T l) < x then some (some (t.partner neSt nE s T Q l)) else s.rotAdj[Q l]?
  capslot : ∀ l, 4 ≤ l → l < 4 * nE → T l < 4 → sx.rotAdj[Q l]? = s.rotAdj[Q l]?
  w1 : ∀ l, t.Corner neSt nE s T l →
    sx.rotAdj[t.w1 neSt s l]? =
      if l < x then some (some (Q l)) else s.rotAdj[t.w1 neSt s l]?
  w0 : ∀ l, t.Corner neSt nE s T l →
    sx.rotAdj[t.w0 neSt s l]? =
      if l < x then some (some (Q (T l))) else s.rotAdj[t.w0 neSt s l]?
  other : ∀ q, (∀ l, 4 ≤ l → l < 4 * nE → q ≠ Q l) →
    (∀ l, t.Corner neSt nE s T l → q ≠ t.w0 neSt s l ∧ q ≠ t.w1 neSt s l) →
    sx.rotAdj[q]? = s.rotAdj[q]?
  row : ∃ r : Array (Option Nat), sx.outerE[i]? = some r ∧ r.size = 4 ∧
    ∀ l, l < 4 → r[2 * QE.side l + (1 - QE.dir l)]? =
      if l < x then some (some (Q (T l))) else s.outerE[i]![2 * QE.side l + (1 - QE.dir l)]?

variable {i neSt nE s T Q}

theorem rFold_succ (x : Nat) :
    t.rFold i neSt s (x + 1) = t.nodeStep i neSt (t.rFold i neSt s x) (4 * neSt + x) := by
  simp [rFold, List.range'_1_concat]

theorem invR_zero (H : t.RFold i neSt nE s T Q) : t.InvR i neSt nE s T Q 0 (t.rFold i neSt s 0) := by
  obtain ⟨r, hr, hr4⟩ := H.row
  obtain ⟨hi, hri⟩ := Array.getElem?_eq_some_iff.1 hr
  refine ⟨rfl, fun _ _ => rfl, ?_, fun _ _ _ _ => rfl, ?_, ?_, fun _ _ _ => rfl, r, hr, hr4, ?_⟩
  · intro l _ _ _; simp [rFold]
  · intro l _; simp [rFold]
  · intro l _; simp [rFold]
  · intro l _; rw [getElem!_pos s.outerE i hi, hri]; simp

theorem side_dir_inj {a b : Nat} (ha : a < 4) (hb : b < 4)
    (h : 2 * QE.side a + (1 - QE.dir a) = 2 * QE.side b + (1 - QE.dir b)) : a = b := by
  unfold QE.side QE.dir at h; omega

theorem invR_step (H : t.RFold i neSt nE s T Q) {x : Nat} (hx : x < 4 * nE)
    (h : t.InvR i neSt nE s T Q x (t.rFold i neSt s x)) :
    t.InvR i neSt nE s T Q (x + 1) (t.rFold i neSt s (x + 1)) := by
  rw [rFold_succ]
  set sx := t.rFold i neSt s x with hsx
  have hrot := H.rot x hx
  have hTx := H.T_lt x hx
  have hTT := H.T_inv x hx
  have hTe := H.T_edge x hx
  have hqe : ∀ l, 4 ≤ l → l < 4 * nE → t.treeQe sx (4 * neSt + l) = some (Q l) :=
    fun l h4 hl => H.qe l h4 hl sx h.rows
  -- bookkeeping: `min l (T l) = x` iff `l = x` or `l = T x`
  have hmin : ∀ l, l < 4 * nE → (min l (T l) < x + 1 ↔ min l (T l) < x ∨ l = x ∨ l = T x) := by
    intro l hl
    constructor
    · intro hlt
      by_cases hc : min l (T l) < x
      · exact Or.inl hc
      · have : min l (T l) = x := by omega
        by_cases hlx : l = x
        · exact Or.inr (Or.inl hlx)
        · right; right
          have : T l = x := by omega
          rw [← this, H.T_inv l hl]
    · rintro (hc | rfl | rfl)
      · omega
      · omega
      · rw [hTT]; omega
  by_cases hlt : T x < x
  · -- skip: the pair was processed at `T x`
    rw [t.nodeStep_skip i neSt sx hrot (by omega)]
    refine ⟨h.rsize, h.rows, ?_, h.capslot, ?_, ?_, h.other, ?_⟩
    · intro l h4 hl hT
      rw [h.slot l h4 hl hT]
      have : (min l (T l) < x + 1) ↔ (min l (T l) < x) := by
        rw [hmin l hl]
        constructor
        · rintro (hc | rfl | rfl)
          · exact hc
          · omega
          · rw [hTT]; omega
        · exact Or.inl
      simp only [this]
    · intro l hc
      rw [h.w1 l hc]
      have : l ≠ x := fun e => by subst e; have := H.v_lt l hc; omega
      have : (l < x + 1) ↔ (l < x) := by omega
      simp only [this]
    · intro l hc
      rw [h.w0 l hc]
      have : l ≠ x := fun e => by subst e; have := H.v_lt l hc; omega
      have : (l < x + 1) ↔ (l < x) := by omega
      simp only [this]
    · obtain ⟨r, hr, hr4, hrk⟩ := h.row
      refine ⟨r, hr, hr4, fun l hl => ?_⟩
      rw [hrk l hl]
      have : l ≠ x := by
        intro e; subst e
        have : T l / 4 = l / 4 := by omega
        exact hTe this
      have : (l < x + 1) ↔ (l < x) := by omega
      simp only [this]
  by_cases hcap : x < 4
  · -- cap quarter-edge: record the facing child slot in the node's row
    have hTx4 : 4 ≤ T x := by
      by_contra hh
      exact hTe (by omega)
    rw [t.nodeStep_cap i neSt sx hrot (by omega) (by omega)]
    have hside : QE.side (4 * neSt + x) = QE.side x := by unfold QE.side; omega
    have hdir : QE.dir (4 * neSt + x) = QE.dir x := by unfold QE.dir; omega
    rw [hside, hdir, hqe (T x) hTx4 hTx]
    refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
    · rw [setOuter_rotAdj]; exact h.rsize
    · intro j hj; rw [setOuter_outerE_ne _ _ _ _ _ hj]; exact h.rows j hj
    · intro l h4 hl hT
      rw [setOuter_rotAdj, h.slot l h4 hl hT]
      have : (min l (T l) < x + 1) ↔ (min l (T l) < x) := by
        rw [hmin l hl]
        constructor
        · rintro (hc | rfl | rfl)
          · exact hc
          · omega
          · omega
        · exact Or.inl
      simp only [this]
    · intro l h4 hl hT; rw [setOuter_rotAdj]; exact h.capslot l h4 hl hT
    · intro l hc
      rw [setOuter_rotAdj, h.w1 l hc]
      have : l ≠ x := fun e => by subst e; have := hc.1; omega
      have : (l < x + 1) ↔ (l < x) := by omega
      simp only [this]
    · intro l hc
      rw [setOuter_rotAdj, h.w0 l hc]
      have : l ≠ x := fun e => by subst e; have := hc.1; omega
      have : (l < x + 1) ↔ (l < x) := by omega
      simp only [this]
    · intro q hq hw; rw [setOuter_rotAdj]; exact h.other q hq hw
    · obtain ⟨r, hr, hr4, hrk⟩ := h.row
      refine ⟨r.set! (2 * QE.side x + (1 - QE.dir x)) (some (Q (T x))),
        setOuter_outerE_eq' _ _ _ _ hr, by simp [hr4], fun l hl => ?_⟩
      have hidx : 2 * QE.side x + (1 - QE.dir x) < r.size := by
        rw [hr4]; unfold QE.side QE.dir; omega
      simp only [Array.set!, Array.getElem?_setIfInBounds, hidx, ↓reduceIte]
      by_cases hlx : l = x
      · subst hlx; simp
      · have hne : 2 * QE.side x + (1 - QE.dir x) ≠ 2 * QE.side l + (1 - QE.dir l) :=
          fun e => hlx (side_dir_inj hcap hl e).symm
        rw [if_neg hne, hrk l hl]
        have : (l < x + 1) ↔ (l < x) := by omega
        simp only [this]
  -- a non-cap quarter-edge with `x < T x`
  have hx4 : 4 ≤ x := by omega
  have hxT : x < T x := by omega
  have hTx4 : 4 ≤ T x := by omega
  have hQne : Q x ≠ Q (T x) := fun e => by
    have := H.Q_inj x (T x) hx4 hx hTx4 hTx e; omega
  have hQx := H.Q_lt x hx4 hx
  have hQTx := H.Q_lt (T x) hTx4 hTx
  have hcornTx : ¬ t.Corner neSt nE s T (T x) := fun hc => by
    have := H.v_lt (T x) hc; rw [hTT] at this; omega
  -- the common part of a plain link `Q x ↔ Q (T x)`
  have plain : t.InvR i neSt nE s T Q (x + 1) ((link (some (Q x)) (some (Q (T x)))).run sx).2 →
      True := fun _ => trivial
  have hplain : ¬ t.Corner neSt nE s T x →
      t.InvR i neSt nE s T Q (x + 1) ((link (some (Q x)) (some (Q (T x)))).run sx).2 := by
    intro hcornx
    have hget := fun q => link_rotAdj_get (Q x) (Q (T x)) sx (by rw [h.rsize]; exact hQx)
      (by rw [h.rsize]; exact hQTx) hQne q
    have hpx : t.partner neSt nE s T Q x = Q (T x) := by
      simp [partner, hcornx, hcornTx]
    have hpTx : t.partner neSt nE s T Q (T x) = Q x := by
      simp [partner, hcornx, hcornTx, hTT]
    refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
    · rw [link_rotAdj_size]; exact h.rsize
    · intro j hj; rw [link_outerE]; exact h.rows j hj
    · intro l h4 hl hT
      rw [hget]
      by_cases hlx : l = x
      · subst hlx; simp [hpx]
      by_cases hlT : l = T x
      · subst hlT; rw [if_neg (Ne.symm hQne), if_pos rfl, hpTx, if_pos (by rw [hTT]; omega)]
      have h1 : Q l ≠ Q x := fun e => hlx (H.Q_inj l x h4 hl hx4 hx e)
      have h2 : Q l ≠ Q (T x) := fun e => hlT (H.Q_inj l (T x) h4 hl hTx4 hTx e)
      rw [if_neg h1, if_neg h2, h.slot l h4 hl hT]
      have : (min l (T l) < x + 1) ↔ (min l (T l) < x) := by
        rw [hmin l hl]
        constructor
        · rintro (hc | rfl | rfl)
          · exact hc
          · exact absurd rfl hlx
          · exact absurd rfl hlT
        · exact Or.inl
      simp only [this]
    · intro l h4 hl hT
      rw [hget]
      have h1 : Q l ≠ Q x := fun e => by have := H.Q_inj l x h4 hl hx4 hx e; subst this; omega
      have h2 : Q l ≠ Q (T x) := fun e => by
        have := H.Q_inj l (T x) h4 hl hTx4 hTx e; subst this; rw [hTT] at hT; omega
      rw [if_neg h1, if_neg h2]; exact h.capslot l h4 hl hT
    · intro l hc
      rw [hget]
      have := H.wQ l hc
      rw [if_neg (Ne.symm (this x hx4 hx).2), if_neg (Ne.symm (this (T x) hTx4 hTx).2), h.w1 l hc]
      have : l ≠ x := fun e => hcornx (e ▸ hc)
      have : (l < x + 1) ↔ (l < x) := by omega
      simp only [this]
    · intro l hc
      rw [hget]
      have := H.wQ l hc
      rw [if_neg (Ne.symm (this x hx4 hx).1), if_neg (Ne.symm (this (T x) hTx4 hTx).1), h.w0 l hc]
      have : l ≠ x := fun e => hcornx (e ▸ hc)
      have : (l < x + 1) ↔ (l < x) := by omega
      simp only [this]
    · intro q hq hw
      rw [hget, if_neg (hq x hx4 hx), if_neg (hq (T x) hTx4 hTx)]
      exact h.other q hq hw
    · obtain ⟨r, hr, hr4, hrk⟩ := h.row
      refine ⟨r, by rw [link_outerE]; exact hr, hr4, fun l hl => ?_⟩
      rw [hrk l hl]
      have : (l < x + 1) ↔ (l < x) := by omega
      simp only [this]
  by_cases hshape : x % 4 = 2 ∧ T x % 4 = 1
  · rw [t.nodeStep_corner i neSt sx hrot (by omega) (by omega) (by omega) (by omega)]
    have hv : QE.edge (4 * neSt + x) = QE.edge (4 * neSt + x) := rfl
    show t.InvR i neSt nE s T Q (x + 1)
      (if (sx.outerE[t.vOf (4 * neSt + x)]![0]!).isSome then
        ((link (t.treeQe sx (4 * neSt + T x)) sx.outerE[t.vOf (4 * neSt + x)]![0]!).run
          ((link (t.treeQe sx (4 * neSt + x)) sx.outerE[t.vOf (4 * neSt + x)]![1]!).run sx).2).2
      else ((link (t.treeQe sx (4 * neSt + x)) (t.treeQe sx (4 * neSt + T x))).run sx).2)
    have hrowv : sx.outerE[t.vOf (4 * neSt + x)]! = s.outerE[t.vOf (4 * neSt + x)]! := by
      have := h.rows _ (H.v_ne x hx4 hx)
      simp only [getElem!_def, this]
    rw [hrowv, hqe x hx4 hx, hqe (T x) hTx4 hTx]
    by_cases hsome : (s.outerE[t.vOf (4 * neSt + x)]![0]!).isSome = true
    · rw [if_pos hsome]
      have hcornx : t.Corner neSt nE s T x := ⟨hx4, hx, hshape.1, hshape.2, hsome⟩
      have hs1 := H.v_some1 x hcornx
      obtain ⟨a0, ha0⟩ := Option.isSome_iff_exists.1 hsome
      obtain ⟨a1, ha1⟩ := Option.isSome_iff_exists.1 hs1
      have hw0 : t.w0 neSt s x = a0 := by simp [w0, ha0]
      have hw1 : t.w1 neSt s x = a1 := by simp [w1, ha1]
      rw [ha0, ha1, ← hw0, ← hw1]
      obtain ⟨hw0lt, hw1lt⟩ := H.w_lt x hcornx
      have hwne := H.w_ne x hcornx
      have hwQ := H.wQ x hcornx
      have hQxw1 : Q x ≠ t.w1 neSt s x := (hwQ x hx4 hx).2
      have hQTw0 : Q (T x) ≠ t.w0 neSt s x := (hwQ (T x) hTx4 hTx).1
      set s1 := ((link (some (Q x)) (some (t.w1 neSt s x))).run sx).2 with hs1def
      have hget1 := fun q => link_rotAdj_get (Q x) (t.w1 neSt s x) sx
        (by rw [h.rsize]; exact hQx) (by rw [h.rsize]; exact hw1lt) hQxw1 q
      have hsz1 : s1.rotAdj.size = s.rotAdj.size := by rw [hs1def, link_rotAdj_size]; exact h.rsize
      have hget2 := fun q => link_rotAdj_get (Q (T x)) (t.w0 neSt s x) s1
        (by rw [hsz1]; exact hQTx) (by rw [hsz1]; exact hw0lt) hQTw0 q
      have hpx : t.partner neSt nE s T Q x = t.w1 neSt s x := by simp [partner, hcornx]
      have hpTx : t.partner neSt nE s T Q (T x) = t.w0 neSt s x := by
        simp [partner, hcornTx, hTT, hcornx]
      have hw0Qx : t.w0 neSt s x ≠ Q x := Ne.symm (hwQ x hx4 hx).1
      have hw0w1 : t.w0 neSt s x ≠ t.w1 neSt s x := hwne
      have hw1QT : t.w1 neSt s x ≠ Q (T x) := Ne.symm (hwQ (T x) hTx4 hTx).2
      refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
      · rw [link_rotAdj_size]; exact hsz1
      · intro j hj; rw [link_outerE, hs1def, link_outerE]; exact h.rows j hj
      · intro l h4 hl hT
        rw [hget2]
        by_cases hlT : l = T x
        · subst hlT; rw [if_pos rfl, hpTx, if_pos (by rw [hTT]; omega)]
        have h2 : Q l ≠ Q (T x) := fun e => hlT (H.Q_inj l (T x) h4 hl hTx4 hTx e)
        rw [if_neg h2, if_neg (hwQ l h4 hl).1, hget1]
        by_cases hlx : l = x
        · subst hlx; rw [if_pos rfl, hpx, if_pos (by omega)]
        have h1 : Q l ≠ Q x := fun e => hlx (H.Q_inj l x h4 hl hx4 hx e)
        rw [if_neg h1, if_neg (hwQ l h4 hl).2, h.slot l h4 hl hT]
        have : (min l (T l) < x + 1) ↔ (min l (T l) < x) := by
          rw [hmin l hl]
          constructor
          · rintro (hc | rfl | rfl)
            · exact hc
            · exact absurd rfl hlx
            · exact absurd rfl hlT
          · exact Or.inl
        simp only [this]
      · intro l h4 hl hT
        have h1 : Q l ≠ Q x := fun e => by have := H.Q_inj l x h4 hl hx4 hx e; subst this; omega
        have h2 : Q l ≠ Q (T x) := fun e => by
          have := H.Q_inj l (T x) h4 hl hTx4 hTx e; subst this; rw [hTT] at hT; omega
        rw [hget2, if_neg h2, if_neg (hwQ l h4 hl).1, hget1, if_neg h1, if_neg (hwQ l h4 hl).2]
        exact h.capslot l h4 hl hT
      · intro l hc
        rw [hget2]
        by_cases hlx : l = x
        · subst hlx
          rw [if_neg hw1QT, if_neg (Ne.symm hw0w1), hget1, if_neg (Ne.symm hQxw1),
            if_pos rfl, if_pos (by omega)]
        have hww := H.ww l x hc hcornx hlx
        rw [if_neg (Ne.symm (H.wQ l hc (T x) hTx4 hTx).2), if_neg hww.2.2.1, hget1,
          if_neg (Ne.symm (H.wQ l hc x hx4 hx).2), if_neg hww.2.2.2, h.w1 l hc]
        have : (l < x + 1) ↔ (l < x) := by omega
        simp only [this]
      · intro l hc
        rw [hget2]
        by_cases hlx : l = x
        · subst hlx
          rw [if_neg (Ne.symm hQTw0), if_pos rfl, if_pos (by omega)]
        have hww := H.ww l x hc hcornx hlx
        rw [if_neg (Ne.symm (H.wQ l hc (T x) hTx4 hTx).1), if_neg hww.1, hget1,
          if_neg (Ne.symm (H.wQ l hc x hx4 hx).1), if_neg hww.2.1, h.w0 l hc]
        have : (l < x + 1) ↔ (l < x) := by omega
        simp only [this]
      · intro q hq hw
        rw [hget2, if_neg (hq (T x) hTx4 hTx), if_neg (hw x hcornx).1, hget1, if_neg (hq x hx4 hx),
          if_neg (hw x hcornx).2]
        exact h.other q hq hw
      · obtain ⟨r, hr, hr4, hrk⟩ := h.row
        refine ⟨r, by rw [link_outerE, hs1def, link_outerE]; exact hr, hr4, fun l hl => ?_⟩
        rw [hrk l hl]
        have : (l < x + 1) ↔ (l < x) := by omega
        simp only [this]
    · rw [if_neg hsome]
      exact hplain fun hc => hsome hc.2.2.2.2
  · rw [t.nodeStep_link i neSt sx hrot (by omega) (by omega) (by
      intro hh; apply hshape; constructor <;> omega)]
    rw [hqe x hx4 hx, hqe (T x) hTx4 hTx]
    exact hplain fun hc => hshape ⟨hc.2.2.1, hc.2.2.2.1⟩

theorem invR_all (H : t.RFold i neSt nE s T Q) {x : Nat} (hx : x ≤ 4 * nE) :
    t.InvR i neSt nE s T Q x (t.rFold i neSt s x) := by
  induction x with
  | zero => exact t.invR_zero H
  | succ x ih => exact t.invR_step H (by omega) (ih (by omega))

/-- The node fold of an `R` node, pointwise. -/
theorem nodeFoldR (H : t.RFold i neSt nE s T Q) :
    let s' := (List.range' (4 * neSt) (4 * nE)).foldl (t.nodeStep i neSt) s
    s'.rotAdj.size = s.rotAdj.size ∧
    (∀ j, j ≠ i → s'.outerE[j]? = s.outerE[j]?) ∧
    (∀ l, 4 ≤ l → l < 4 * nE → 4 ≤ T l →
      s'.rotAdj[Q l]? = some (some (t.partner neSt nE s T Q l))) ∧
    (∀ l, 4 ≤ l → l < 4 * nE → T l < 4 → s'.rotAdj[Q l]? = s.rotAdj[Q l]?) ∧
    (∀ l, t.Corner neSt nE s T l →
      s'.rotAdj[t.w1 neSt s l]? = some (some (Q l)) ∧
      s'.rotAdj[t.w0 neSt s l]? = some (some (Q (T l)))) ∧
    (∀ q, (∀ l, 4 ≤ l → l < 4 * nE → q ≠ Q l) →
      (∀ l, t.Corner neSt nE s T l → q ≠ t.w0 neSt s l ∧ q ≠ t.w1 neSt s l) →
      s'.rotAdj[q]? = s.rotAdj[q]?) ∧
    (∃ r : Array (Option Nat), s'.outerE[i]? = some r ∧ r.size = 4 ∧
      ∀ l, l < 4 → r[2 * QE.side l + (1 - QE.dir l)]? = some (some (Q (T l)))) := by
  intro s'
  have h := t.invR_all H (x := 4 * nE) le_rfl
  have hs' : s' = t.rFold i neSt s (4 * nE) := rfl
  rw [hs']
  obtain ⟨r, hr, hr4, hrk⟩ := h.row
  refine ⟨h.rsize, h.rows, ?_, h.capslot, ?_, h.other, r, hr, hr4, ?_⟩
  · intro l h4 hl hT
    rw [h.slot l h4 hl hT, if_pos (by have := H.T_lt l hl; omega)]
  · intro l hc
    exact ⟨by rw [h.w1 l hc, if_pos hc.2.1], by rw [h.w0 l hc, if_pos hc.2.1]⟩
  · intro l hl
    rw [hrk l hl, if_pos (by have := H.nE_pos; omega)]

end RFold

end Spqr.PlanarSpqrTree
