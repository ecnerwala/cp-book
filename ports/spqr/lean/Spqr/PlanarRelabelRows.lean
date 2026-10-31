import Spqr.PlanarRelabelSem
import Spqr.PlanarRotSpec
import Spqr.PlanarRelabelNodeR
import Spqr.PlanarLayout
import Spqr.PlanarRelabelQem
import Spqr.PlanarRotFold

/-!
# `planarRelabel`: the rows of planar R nodes

`RowInv` is carried through `planarRelabel`: every planar R node laid out so far has its
`RelabelNodeR` data (`NodeRow`) readable off the state's arrays. The `qem` read by `mapRot` at an R
node is the walk's `qem` after that node's own `capLinked`/`flipped` only, because the layout of a
subtree writes only the cap slots and the slots of virtual edges below it (`SubSlot`), and
distinct subtrees have distinct virtual edges (`Fresh`).
-/

namespace Spqr
open PlanarRelabelM

/-- Quarter-edge slots the layout of the subtree of `cur` may write: the cap's, and those of the
virtual edges of children of items below `cur`. -/
def SubSlot (g : Graph) (items : Items) (cur : ItemId) (q : Nat) : Prop :=
  8 * g.ne ≤ q ∨ ∃ j c, Items.Below items cur j ∧ c ∈ Items.ch items j ∧ 1 + g.nv ≤ c ∧
    4 * (c - (1 + g.nv)) ≤ q ∧ q < 4 * (c - (1 + g.nv)) + 4

/-- `qem` of `s` agrees with the walk's on the slots of the virtual edges below `cur`. -/
def Fresh (g : Graph) (w : PlanarWalkState) (cur : ItemId) (s : PlanarRelabelState) : Prop :=
  ∀ j c, Items.Below w.base.items cur j → c ∈ Items.ch w.base.items j → 1 + g.nv ≤ c →
    c < 1 + g.nv + 2 * g.ne →
    ∀ z, z < 4 → s.aux.qem[4 * (c - (1 + g.nv)) + z]! = w.aux.qem[4 * (c - (1 + g.nv)) + z]!

/-- Endpoint pairs of the node-edges laid out so far. -/
def nvsArr (s : PlanarRelabelState) : Array (Nat × Nat) := s.base.nodeEdges.map (·.nvs)

/-- `RelabelNodeR` for node `n`, read off the state's arrays. -/
def NodeRow (g : Graph) (w : PlanarWalkState) (s : PlanarRelabelState) (n : Nat) : Prop :=
  ∃ (it : ItemId) (m : Array Nat) (children : List ItemId) (pos : Nat → Nat),
    let items := w.base.items
    let nvSt := s.base.nvBounds[n]!
    let nvEn := s.base.nvBounds[n + 1]!
    let neSt := s.base.neBounds[n]!
    let neEn := s.base.neBounds[n + 1]!
    let edgeVes := (children.filter (· ≥ 1 + g.nv)).map (· - (1 + g.nv))
    let ves := capVe g :: edgeVes
    let Q := flipped g items[it]!.ch w.aux.itemFlips[it]! (capLinked g w.aux.qem m)
    1 + g.nv + g.ne ≤ it ∧ it < items.size ∧ items[it]!.type = .R ∧
    w.aux.nodePlanarity[it - (1 + g.nv + g.ne)]? = some (.planar m) ∧
    children = Items.ordered g items it nvSt pos ∧
    Items.PosOK nvSt (Items.nvList g items it) pos ∧
    nvEn = nvSt + (Items.nvList g items it).length ∧
    neEn = neSt + edgeVes.length + 1 ∧
    neEn ≤ (nvsArr s).size ∧
    (∀ k, k < neEn - neSt →
      (nvsArr s)[neSt + k]! = ((nvSt, nvEn - 1) :: Items.edgeChildren g items pos children)[k]!) ∧
    4 * neEn ≤ s.aux.neRotAdj.size ∧
    ∀ l, l < 4 * ves.length →
      (Q[4 * ves[l / 4]! + l % 4]! = none → s.aux.neRotAdj[4 * neSt + l]! = none) ∧
      ∀ o, Q[4 * ves[l / 4]! + l % 4]! = some o → QE.edge o ∈ ves →
        s.aux.neRotAdj[4 * neSt + l]! =
          some (4 * (neSt + ves.idxOf (QE.edge o)) + (o &&& 2) + (1 - l % 2))

/-- Invariant of `planarRelabel`: bookkeeping sizes, the constant fields, and `NodeRow` at every
planar R node laid out so far. -/
structure RowInv (g : Graph) (w : PlanarWalkState) (s : PlanarRelabelState) : Prop where
  g_eq : s.base.g = g
  items : s.base.items = w.base.items
  np : s.aux.nodePlanarity = w.aux.nodePlanarity
  fl : s.aux.itemFlips = w.aux.itemFlips
  rne_size : s.aux.rotEdgeNe.size = 2 * g.ne + 1
  ne_size : s.base.neBounds.size = s.base.types.size + 1
  nv_size : s.base.nvBounds.size = s.base.types.size + 1
  np_size : s.aux.nodePlanar.size = s.base.types.size
  rot_size : s.aux.neRotAdj.size = 4 * s.base.nodeEdges.size
  qem_size : s.aux.qem.size = w.aux.qem.size
  vp_size : s.base.vertPos.size = g.nv
  ne_last : s.base.neBounds[s.base.types.size]! = s.base.nodeEdges.size
  nv_last : s.base.nvBounds[s.base.types.size]! = s.base.nodeVerts.size
  rows : ∀ n, n < s.base.types.size → s.base.types[n]! = .R → s.aux.nodePlanar[n]! = true →
    NodeRow g w s n

/-- The fields `RowInv`/`NodeRow`/`Fresh` read, unchanged. -/
structure RFr (s s' : PlanarRelabelState) : Prop where
  g : s'.base.g = s.base.g
  items : s'.base.items = s.base.items
  types : s'.base.types = s.base.types
  neBounds : s'.base.neBounds = s.base.neBounds
  nvBounds : s'.base.nvBounds = s.base.nvBounds
  nvs : nvsArr s' = nvsArr s
  nodeVerts : s'.base.nodeVerts.size = s.base.nodeVerts.size
  vertPos : s'.base.vertPos = s.base.vertPos
  qem : s'.aux.qem = s.aux.qem
  np : s'.aux.nodePlanarity = s.aux.nodePlanarity
  fl : s'.aux.itemFlips = s.aux.itemFlips
  rne : s'.aux.rotEdgeNe = s.aux.rotEdgeNe
  neRotAdj : s'.aux.neRotAdj = s.aux.neRotAdj
  nodePlanar : s'.aux.nodePlanar = s.aux.nodePlanar

theorem RFr.refl (s : PlanarRelabelState) : RFr s s :=
  ⟨rfl, rfl, rfl, rfl, rfl, rfl, rfl, rfl, rfl, rfl, rfl, rfl, rfl, rfl⟩
theorem RFr.trans {s₁ s₂ s₃ : PlanarRelabelState} (h₁ : RFr s₁ s₂) (h₂ : RFr s₂ s₃) : RFr s₁ s₃ :=
  ⟨h₂.g.trans h₁.g, h₂.items.trans h₁.items, h₂.types.trans h₁.types, h₂.neBounds.trans h₁.neBounds,
    h₂.nvBounds.trans h₁.nvBounds, h₂.nvs.trans h₁.nvs, h₂.nodeVerts.trans h₁.nodeVerts,
    h₂.vertPos.trans h₁.vertPos, h₂.qem.trans h₁.qem, h₂.np.trans h₁.np, h₂.fl.trans h₁.fl,
    h₂.rne.trans h₁.rne, h₂.neRotAdj.trans h₁.neRotAdj, h₂.nodePlanar.trans h₁.nodePlanar⟩

theorem nvsArr_size (s : PlanarRelabelState) : (nvsArr s).size = s.base.nodeEdges.size := by
  simp [nvsArr]

theorem NodeRow.fr {g : Graph} {w : PlanarWalkState} {s s' : PlanarRelabelState} (hf : RFr s s') {n : Nat}
    (h : NodeRow g w s n) : NodeRow g w s' n := by
  unfold NodeRow at h ⊢
  rw [hf.nvBounds, hf.neBounds, hf.nvs, hf.neRotAdj]
  exact h

theorem RowInv.fr {g : Graph} {w : PlanarWalkState} {s s' : PlanarRelabelState} (hf : RFr s s')
    (h : RowInv g w s) : RowInv g w s' := by
  obtain ⟨h1, h2, h3, h4, h5, h6, h7, h8, h9, h9q, h9v, h9ne, h9nv, h10⟩ := h
  refine ⟨by rw [hf.g, h1], by rw [hf.items, h2], by rw [hf.np, h3], by rw [hf.fl, h4],
    by rw [hf.rne, h5], by rw [hf.neBounds, hf.types, h6], by rw [hf.nvBounds, hf.types, h7],
    by rw [hf.nodePlanar, hf.types, h8], ?_, by rw [hf.qem, h9q], by rw [hf.vertPos, h9v],
    by rw [hf.neBounds, hf.types, h9ne, ← nvsArr_size, ← nvsArr_size, hf.nvs],
    by rw [hf.nvBounds, hf.types, h9nv, hf.nodeVerts], ?_⟩
  · rw [hf.neRotAdj, h9, ← nvsArr_size, ← nvsArr_size, hf.nvs]
  · intro n hn hty hpl
    rw [hf.types] at hn hty
    rw [hf.nodePlanar] at hpl
    exact (h10 n hn hty hpl).fr hf

theorem Fresh.fr {g : Graph} {w : PlanarWalkState} {cur : ItemId} {s s' : PlanarRelabelState} (hf : RFr s s')
    (h : Fresh g w cur s) : Fresh g w cur s' := by
  unfold Fresh at h ⊢; rw [hf.qem]; exact h

/-! ### Array helpers -/

theorem Array.getElem!_append_left' {α : Type} [Inhabited α] (a b : Array α) (i : Nat) (h : i < a.size) :
    (a ++ b)[i]! = a[i]! := by
  rw [getElem!_pos (a ++ b) i (by simp; omega), getElem!_pos a i h, Array.getElem_append_left h]

theorem Array.getElem!_append_right' {α : Type} [Inhabited α] (a b : Array α) (i : Nat) (h : i < b.size) :
    (a ++ b)[a.size + i]! = b[i]! := by
  rw [getElem!_pos (a ++ b) (a.size + i) (by simp; omega), getElem!_pos b i h,
    Array.getElem_append_right (by omega)]
  congr 1; omega

theorem Array.getElem!_map' {α β : Type} [Inhabited α] [Inhabited β] (f : α → β) (a : Array α) (i : Nat)
    (h : i < a.size) : (a.map f)[i]! = f a[i]! := by
  rw [getElem!_pos (a.map f) i (by simpa using h), getElem!_pos a i h, Array.getElem_map]

theorem nvsArr_modify_twin (s : PlanarRelabelState) (k : Nat) (t : Option Nat) :
    (s.base.nodeEdges.modify k fun ne => { ne with twin := t }).map (·.nvs) = nvsArr s := by
  unfold nvsArr
  refine Array.ext (by simp) fun i h1 h2 => ?_
  simp only [Array.getElem_map, Array.getElem_modify]
  split <;> rfl

/-- Entry `i` of a fold of 4-blocks `f` appended to `acc`. -/
theorem foldl_append_get {α : Type} [Inhabited α] (f : Nat → Array α) (hf : ∀ v, (f v).size = 4)
    (l : List Nat) (acc : Array α) (i : Nat) (hi : i < acc.size + 4 * l.length) :
    (l.foldl (fun acc ve => acc ++ f ve) acc)[i]! =
      if i < acc.size then acc[i]! else (f l[(i - acc.size) / 4]!)[(i - acc.size) % 4]! := by
  induction l generalizing acc with
  | nil => simp at hi; simp [hi]
  | cons v l ih =>
    rw [List.foldl_cons, ih _ (by simp only [Array.size_append, hf, List.length_cons] at hi ⊢; omega)]
    simp only [Array.size_append, hf]
    by_cases h1 : i < acc.size
    · rw [if_pos (by omega), if_pos h1, Array.getElem!_append_left' _ _ _ h1]
    · by_cases h2 : i < acc.size + 4
      · rw [if_pos h2, if_neg h1]
        obtain ⟨j, rfl⟩ : ∃ j, i = acc.size + j := ⟨i - acc.size, by omega⟩
        rw [Array.getElem!_append_right' _ _ _ (by rw [hf]; omega)]
        have : (acc.size + j - acc.size) / 4 = 0 := by omega
        rw [this, Nat.add_sub_cancel_left, Nat.mod_eq_of_lt (by omega)]
        rfl
      · rw [if_neg h2, if_neg h1]
        have e1 : (i - (acc.size + 4)) / 4 + 1 = (i - acc.size) / 4 := by omega
        have e2 : (i - (acc.size + 4)) % 4 = (i - acc.size) % 4 := by omega
        rw [← e1, ← e2, List.getElem!_cons_succ]

/-- The rows of an R node's `layoutRot`: block `l / 4` is `mapRot` of the `l / 4`-th edge of
`capVe :: edgeVes`. -/
theorem layoutRot_R_get (nV neSt neEn : Nat) (edgeVes : List Nat) (mapRot : Nat → Array (Option Nat))
    (cv : Nat) (hnV : nV ≠ 1) (hmr : ∀ ve, (mapRot ve).size = 4) (l : Nat) (hl : l < 4 * (edgeVes.length + 1)) :
    (layoutRot .R nV neSt neEn edgeVes mapRot cv)[l]! = (mapRot (cv :: edgeVes)[l / 4]!)[l % 4]! := by
  have h1 : (nV == 1) = false := by simpa using hnV
  simp only [layoutRot, h1, Bool.false_eq_true, ↓reduceIte,
    show (NodeType.R == NodeType.Q || NodeType.R == NodeType.I) = false from rfl,
    show (NodeType.R == NodeType.P) = false from rfl, show (NodeType.R == NodeType.S) = false from rfl]
  rw [foldl_append_get mapRot hmr edgeVes (mapRot cv) l (by rw [hmr]; omega), hmr]
  by_cases h4 : l < 4
  · rw [if_pos h4, Nat.div_eq_of_lt h4, Nat.mod_eq_of_lt h4]; rfl
  · rw [if_neg h4]
    have e1 : (l - 4) / 4 + 1 = l / 4 := by omega
    have e2 : (l - 4) % 4 = l % 4 := by omega
    rw [← e1, ← e2, List.getElem!_cons_succ]

/-- Size of `layoutRot`, from the per-type edge counts. -/
theorem layoutRot_size_of (ty : NodeType) (nV neSt neEn : Nat) (edgeVes : List Nat)
    (mapRot : Nat → Array (Option Nat)) (cv : Nat) (hmr : ∀ ve, (mapRot ve).size = 4)
    (hle : neSt ≤ neEn) (hFV : ty = .F ∨ ty = .V → neEn = neSt)
    (h1 : ty ≠ .F → ty ≠ .V → nV = 1 → neEn = neSt + 1)
    (hQI : ty = .Q ∨ ty = .I → neEn = neSt + 1)
    (hP : ty = .P → nV ≠ 1 → nV = 2)
    (hS : ty = .S → nV ≠ 1 → neEn = neSt + nV ∧ 2 ≤ nV)
    (hRO : ty = .R ∨ ty = .O → neEn = neSt + edgeVes.length + 1) :
    (layoutRot ty nV neSt neEn edgeVes mapRot cv).size = 4 * (neEn - neSt) := by
  have one : ty ≠ NodeType.F → ty ≠ .V → nV = 1 →
      (layoutRot ty nV neSt neEn edgeVes mapRot cv).size = 4 * (neEn - neSt) := by
    intro hF hV hv
    rw [h1 hF hV hv, Nat.add_sub_cancel_left]
    subst hv
    cases ty <;> simp [layoutRot] at hF hV ⊢
  have fold : ty = NodeType.R ∨ ty = .O → nV ≠ 1 →
      (layoutRot ty nV neSt neEn edgeVes mapRot cv).size = 4 * (neEn - neSt) := by
    intro hty hv
    have hv' : (nV == 1) = false := by simpa using hv
    rw [hRO hty]
    rcases hty with rfl | rfl <;> simp [layoutRot, hv', PlanarRot.foldl_append_size mapRot hmr, hmr] <;> omega
  cases ty
  case F => rw [hFV (Or.inl rfl)]; simp [layoutRot]
  case V => rw [hFV (Or.inr rfl)]; simp [layoutRot]
  case R => by_cases hv : nV = 1
            · exact one (by simp) (by simp) hv
            · exact fold (Or.inl rfl) hv
  case O => by_cases hv : nV = 1
            · exact one (by simp) (by simp) hv
            · exact fold (Or.inr rfl) hv
  case Q =>
    by_cases hv : nV = 1
    · exact one (by simp) (by simp) hv
    · have hv' : (nV == 1) = false := by simpa using hv
      rw [hQI (Or.inl rfl)]; simp [layoutRot, hv']
  case I =>
    by_cases hv : nV = 1
    · exact one (by simp) (by simp) hv
    · have hv' : (nV == 1) = false := by simpa using hv
      rw [hQI (Or.inr rfl)]; simp [layoutRot, hv']
  case P =>
    by_cases hv : nV = 1
    · exact one (by simp) (by simp) hv
    · obtain ⟨k, rfl⟩ : ∃ k, neEn = neSt + k := ⟨neEn - neSt, by omega⟩
      rw [hP rfl hv, Nat.add_sub_cancel_left]
      exact layoutRot_P_size k neSt edgeVes mapRot cv
  case S =>
    by_cases hv : nV = 1
    · exact one (by simp) (by simp) hv
    · obtain ⟨h2, h3⟩ := hS rfl hv
      rw [h2, Nat.add_sub_cancel_left]
      exact layoutRot_S_size nV neSt h3 edgeVes mapRot cv

/-! ### `RFr` without `qem`/`vertPos` -/

structure RFrX (s s' : PlanarRelabelState) : Prop where
  g : s'.base.g = s.base.g
  items : s'.base.items = s.base.items
  types : s'.base.types = s.base.types
  neBounds : s'.base.neBounds = s.base.neBounds
  nvBounds : s'.base.nvBounds = s.base.nvBounds
  nvs : nvsArr s' = nvsArr s
  nodeVerts : s'.base.nodeVerts.size = s.base.nodeVerts.size
  np : s'.aux.nodePlanarity = s.aux.nodePlanarity
  fl : s'.aux.itemFlips = s.aux.itemFlips
  rne : s'.aux.rotEdgeNe = s.aux.rotEdgeNe
  neRotAdj : s'.aux.neRotAdj = s.aux.neRotAdj
  nodePlanar : s'.aux.nodePlanar = s.aux.nodePlanar

theorem RFrX.refl (s : PlanarRelabelState) : RFrX s s := ⟨rfl, rfl, rfl, rfl, rfl, rfl, rfl, rfl, rfl, rfl, rfl, rfl⟩
theorem RFr.toX {s s' : PlanarRelabelState} (h : RFr s s') : RFrX s s' :=
  ⟨h.g, h.items, h.types, h.neBounds, h.nvBounds, h.nvs, h.nodeVerts, h.np, h.fl, h.rne, h.neRotAdj, h.nodePlanar⟩
theorem RFrX.toFr {s s' : PlanarRelabelState} (h : RFrX s s') (hq : s'.aux.qem = s.aux.qem)
    (hv : s'.base.vertPos = s.base.vertPos) : RFr s s' :=
  ⟨h.g, h.items, h.types, h.neBounds, h.nvBounds, h.nvs, h.nodeVerts, hv, hq, h.np, h.fl, h.rne, h.neRotAdj,
    h.nodePlanar⟩

/-! ### `wp` helpers -/

theorem wp_rfr {α : Type} {Q : α → PlanarRelabelState → Prop} (m : PlanarRelabelM α) {s : PlanarRelabelState}
    (hm : ∀ s, RFr s (m s).2) (hk : ∀ a s', RFr s s' → Q a s') : wp m Q s :=
  hk _ _ (hm s)

theorem wp_jp_match1_rfr {Q : Unit → PlanarRelabelState → Prop} (ty : NodeType) (m₁ m₂ k : PlanarRelabelM Unit)
    (s : PlanarRelabelState) (h₁ : ∀ s, RFr s (m₁ s).2) (h₂ : ∀ s, RFr s (m₂ s).2)
    (hk : ∀ s', RFr s s' → wp k Q s') :
    wp (planarRelabel.match_1 (fun _ => PlanarRelabelM Unit) ty (fun _ => m₁ >>= fun _ => k)
      (fun _ => m₂ >>= fun _ => k) (fun _ => k)) Q s := by
  cases ty
  all_goals first
    | exact hk s (RFr.refl s)
    | exact wp_rfr m₁ h₁ fun _ s' hf => hk s' hf
    | exact wp_rfr m₂ h₂ fun _ s' hf => hk s' hf

theorem wp_ite_jp {Q : Unit → PlanarRelabelState → Prop} (c : Prop) [Decidable c] (m k : PlanarRelabelM Unit)
    (s : PlanarRelabelState) (hk : ∀ s', (c ∧ s' = (m s).2) ∨ (¬ c ∧ s' = s) → wp k Q s') :
    wp (if c then m >>= fun _ => k else k) Q s := by
  split
  · exact hk _ (Or.inl ⟨‹_›, rfl⟩)
  · exact hk _ (Or.inr ⟨‹_›, rfl⟩)

theorem wp_ite_jp3 {Q : Unit → PlanarRelabelState → Prop} (c : Prop) [Decidable c] (m₁ : PlanarRelabelM Unit)
    (m₂ : PlanarRelabelAux → PlanarRelabelM Unit) (k : PlanarRelabelM Unit) (s : PlanarRelabelState)
    (hk : ∀ s', (c ∧ s' = (m₂ (m₁ s).2.aux (m₁ s).2).2) ∨ (¬ c ∧ s' = s) → wp k Q s') :
    wp (if c then m₁ >>= fun _ => getAux >>= fun a => m₂ a >>= fun _ => k else k) Q s := by
  split
  · exact hk _ (Or.inl ⟨‹_›, rfl⟩)
  · exact hk _ (Or.inr ⟨‹_›, rfl⟩)

/-- Abstract the state after `m` behind the fields `RowInv` reads. -/
theorem wp_rabs {Q : Unit → PlanarRelabelState → Prop} (m : PlanarRelabelM Unit) {s : PlanarRelabelState}
    (hk : ∀ s', s'.base.g = (m s).2.base.g → s'.base.items = (m s).2.base.items →
      s'.base.types = (m s).2.base.types → s'.base.neBounds = (m s).2.base.neBounds →
      s'.base.nvBounds = (m s).2.base.nvBounds → s'.base.nodeEdges = (m s).2.base.nodeEdges →
      s'.base.nodeVerts.size = (m s).2.base.nodeVerts.size → s'.base.vertPos = (m s).2.base.vertPos →
      s'.aux.qem = (m s).2.aux.qem → s'.aux.nodePlanarity = (m s).2.aux.nodePlanarity →
      s'.aux.itemFlips = (m s).2.aux.itemFlips → s'.aux.rotEdgeNe = (m s).2.aux.rotEdgeNe →
      s'.aux.neRotAdj = (m s).2.aux.neRotAdj → s'.aux.nodePlanar = (m s).2.aux.nodePlanar → Q () s') :
    wp m Q s :=
  hk _ rfl rfl rfl rfl rfl rfl rfl rfl rfl rfl rfl rfl rfl rfl

/-- What `setupNode` does to `qem`, by node type and planarity flag. -/
def SetupSpec (g : Graph) (ty : NodeType) (cur : ItemId) (a : PlanarRelabelAux) (b : Bool) (q' : Qem) : Prop :=
  (b = true → (ty = .S ∨ ty = .P ∨ ty = .R) →
    ∃ m, a.nodePlanarity[cur - (1 + g.nv + g.ne)]! = .planar m ∧ q' = capLinked g a.qem m) ∧
  (b = false → q' = a.qem) ∧
  (¬ (ty = .S ∨ ty = .P ∨ ty = .R) → b = true ∧ q' = a.qem)

theorem setupNode_spec (g : Graph) (ty : NodeType) (cur : ItemId) (s : PlanarRelabelState) :
    ∃ b q', setupNode g ty cur s = (b, { s with aux := { s.aux with qem := q' } }) ∧ SetupSpec g ty cur s.aux b q' := by
  unfold SetupSpec
  by_cases hty : ty = .S ∨ ty = .P ∨ ty = .R
  · rcases hm : s.aux.nodePlanarity[cur - (1 + g.nv + g.ne)]! with _ | m | _
    · exact ⟨false, s.aux.qem, setupNode_nonplanar g ty cur s hty (fun m h => by rw [hm] at h; cases h),
        fun h => Bool.noConfusion h, fun _ => rfl, fun h => absurd hty h⟩
    · exact ⟨true, capLinked g s.aux.qem m, setupNode_planar g ty cur s hty m hm,
        fun _ _ => ⟨m, rfl, rfl⟩, fun h => Bool.noConfusion h, fun h => absurd hty h⟩
    · exact ⟨false, s.aux.qem, setupNode_nonplanar g ty cur s hty (fun m h => by rw [hm] at h; cases h),
        fun h => Bool.noConfusion h, fun _ => rfl, fun h => absurd hty h⟩
  · exact ⟨true, s.aux.qem, setupNode_other g ty cur s
      ⟨fun h => hty (Or.inl h), fun h => hty (Or.inr (Or.inl h)), fun h => hty (Or.inr (Or.inr h))⟩,
      fun _ h => absurd h hty, fun h => Bool.noConfusion h, fun _ => ⟨rfl, rfl⟩⟩

theorem wp_setupNode {Q : Bool → PlanarRelabelState → Prop} (g : Graph) (ty : NodeType) (cur : ItemId)
    (s : PlanarRelabelState)
    (hk : ∀ b q', SetupSpec g ty cur s.aux b q' → Q b { s with aux := { s.aux with qem := q' } }) :
    wp (setupNode g ty cur) Q s := by
  obtain ⟨b, q', he, hs⟩ := setupNode_spec g ty cur s
  unfold wp; rw [he]; exact hk b q' hs

theorem nodePlanarity_getElem?_of (a : Array NodePlanarity) (i : Nat) (m : Array Nat) (h : a[i]! = .planar m) :
    a[i]? = some (.planar m) := by
  by_cases hi : i < a.size
  · rw [getElem!_pos a i hi] at h; rw [getElem?_pos a i hi, h]
  · rw [getElem!_neg a i hi] at h
    have e : (default : NodePlanarity) = .unset := rfl
    rw [e] at h; cases h

theorem rotEdgeNe_init_size (edgeVes : List Nat) (neSt : Nat) (r : Array Nat) (ne : Nat) :
    (((edgeVes.zipIdx (neSt + 1)).foldl (init := r) (fun r (ve, ne) => r.set! ve ne)).set! (2 * ne) neSt).size =
      r.size := by
  show (((edgeVes.zipIdx (neSt + 1)).foldl (fun r x => r.set! x.1 x.2) r).set! (2 * ne) neSt).size = r.size
  rw [Array.set!_eq_setIfInBounds, Array.size_setIfInBounds]; exact Ghost.foldl_set!_size _ _ _

/-! ### Transport of `NodeRow` across a node's appends -/

theorem NodeRow.push {g : Graph} {w : PlanarWalkState} {s s' : PlanarRelabelState} {n : Nat} (h : NodeRow g w s n)
    {a b : Nat} {E : Array (Nat × Nat)} {R : Array (Option Nat)}
    (hnvB : s'.base.nvBounds = s.base.nvBounds.push a) (hneB : s'.base.neBounds = s.base.neBounds.push b)
    (hnvs : nvsArr s' = nvsArr s ++ E) (hrot : s'.aux.neRotAdj = s.aux.neRotAdj ++ R)
    (h1 : n + 1 < s.base.nvBounds.size) (h2 : n + 1 < s.base.neBounds.size) : NodeRow g w s' n := by
  unfold NodeRow at h ⊢
  rw [hnvB, hneB, hnvs, hrot, Array.getElem!_push_lt' _ _ _ (by omega), Array.getElem!_push_lt' _ _ _ h1,
    Array.getElem!_push_lt' _ _ _ (by omega), Array.getElem!_push_lt' _ _ _ h2]
  obtain ⟨it, m, children, pos, h⟩ := h
  refine ⟨it, m, children, pos, ?_⟩
  dsimp only at h ⊢
  obtain ⟨h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11, h12⟩ := h
  refine ⟨h1, h2, h3, h4, h5, h6, h7, h8, ?_, ?_, ?_, ?_⟩
  · simp only [Array.size_append]; omega
  · intro k hk; rw [Array.getElem!_append_left' _ _ _ (by omega)]; exact h10 k hk
  · simp only [Array.size_append]; omega
  · intro l hl
    have hl' := hl; simp only [List.length_cons] at hl'
    rw [Array.getElem!_append_left' _ _ _ (by omega)]; exact h12 l hl

/-! ### `Unch`: `qem` outside a subtree's slots -/

/-- `qem` of `s'` agrees with `s` outside the slots the subtree of `cur` may write. -/
def Unch (g : Graph) (items : Items) (cur : ItemId) (s s' : PlanarRelabelState) : Prop :=
  s'.aux.qem.size = s.aux.qem.size ∧ ∀ q, ¬ SubSlot g items cur q → s'.aux.qem[q]! = s.aux.qem[q]!

theorem Unch.refl {g : Graph} {items : Items} {cur : ItemId} (s : PlanarRelabelState) : Unch g items cur s s :=
  ⟨rfl, fun _ _ => rfl⟩
theorem Unch.trans {g : Graph} {items : Items} {cur : ItemId} {s₁ s₂ s₃ : PlanarRelabelState}
    (h₁ : Unch g items cur s₁ s₂) (h₂ : Unch g items cur s₂ s₃) : Unch g items cur s₁ s₃ :=
  ⟨h₂.1.trans h₁.1, fun q hq => (h₂.2 q hq).trans (h₁.2 q hq)⟩
theorem Unch.rfr {g : Graph} {items : Items} {cur : ItemId} {s₁ s₂ s₃ : PlanarRelabelState} (hf : RFr s₂ s₃)
    (h : Unch g items cur s₁ s₂) : Unch g items cur s₁ s₃ := by
  unfold Unch at *; rw [hf.qem]; exact h
theorem Unch.rfr_left {g : Graph} {items : Items} {cur : ItemId} {s₁ s₂ s₃ : PlanarRelabelState} (hf : RFr s₁ s₂)
    (h : Unch g items cur s₂ s₃) : Unch g items cur s₁ s₃ := by
  unfold Unch at *; rw [← hf.qem]; exact h

theorem SubSlot.mono {g : Graph} {items : Items} {cur c : ItemId} (hc : c ∈ Items.ch items cur) {q : Nat}
    (h : SubSlot g items c q) : SubSlot g items cur q := by
  rcases h with h | ⟨j, c', hj, hc', h1, h2, h3⟩
  · exact Or.inl h
  · exact Or.inr ⟨j, c', Relation.ReflTransGen.head hc hj, hc', h1, h2, h3⟩

theorem Unch.of_child {g : Graph} {items : Items} {cur c : ItemId} (hc : c ∈ Items.ch items cur)
    {s s' : PlanarRelabelState} (h : Unch g items c s s') : Unch g items cur s s' :=
  ⟨h.1, fun q hq => h.2 q fun h' => hq (SubSlot.mono hc h')⟩

theorem Items.Tree.below_lt {g : Graph} {items : Items} (ht : items.Tree g) {a x : ItemId} (ha : a < items.size)
    (h : items.Below a x) : x < items.size := by
  rcases Relation.ReflTransGen.cases_tail h with rfl | ⟨b, _, hbx⟩
  · exact ha
  · exact ht.ch_lt b x hbx

theorem Items.Tree.node_ge {g : Graph} {items : Items} (ht : items.Tree g) {i : ItemId} (h : items.type i = .R) :
    1 + g.nv + g.ne ≤ i := by
  by_contra hlt
  by_cases h0 : i = 0
  · subst h0
    have := ht.root; unfold rootItem at this; rw [this] at h; cases h
  · by_cases hv : i < 1 + g.nv
    · have := ht.vert (i - 1) (by (try simp only [ItemId] at *); omega)
      unfold vertItem at this
      rw [show 1 + (i - 1) = i by (try simp only [ItemId] at *); omega, h] at this; cases this
    · have := ht.edge (i - (1 + g.nv)) (by (try simp only [ItemId] at *); omega)
      unfold edgeItem at this
      rw [show 1 + g.nv + (i - (1 + g.nv)) = i by (try simp only [ItemId] at *); omega, h] at this; cases this

/-- A sibling's subtree leaves the slots below `c'` untouched. -/
theorem Fresh.of_unch {g : Graph} {w : PlanarWalkState} (ht : Items.Tree g w.base.items) {p c c' : ItemId}
    (hc : c ∈ Items.ch w.base.items p) (hc' : c' ∈ Items.ch w.base.items p) (hne : c ≠ c')
    {s s' : PlanarRelabelState} (hfr : Fresh g w c' s) (hu : Unch g w.base.items c s s') : Fresh g w c' s' := by
  intro j c'' hj hc'' hge hlt z hz
  rw [← hfr j c'' hj hc'' hge hlt z hz]
  apply hu.2
  rintro (h | ⟨j', d, hj', hd, hge', h1, h2⟩)
  · (try simp only [ItemId] at *); omega
  · have hdc : d = c'' := by (try simp only [ItemId] at *); omega
    subst hdc
    have hsz : d < w.base.items.size := ht.ch_lt _ _ hd
    obtain ⟨par, _, huniq⟩ := ht.unique_parent d (ht.child_pos hd) hsz
    have hjj : j' = j := (huniq j' hd).trans (huniq j hc'').symm
    subst hjj
    have hd1 : j' ∈ Items.desc w.base.items c := Items.mem_desc.2 ⟨ht.below_lt (ht.ch_lt _ _ hc) hj', hj'⟩
    have hd2 : j' ∈ Items.desc w.base.items c' := Items.mem_desc.2 ⟨ht.below_lt (ht.ch_lt _ _ hc') hj, hj⟩
    exact Finset.disjoint_left.1 (ht.desc_disjoint hc hc' hne) hd1 hd2

/-! ### Node-edge counts -/

theorem Items.WF.one_le_nEdges {g : Graph} {items : Items} (hw : items.WF g) {i : ItemId} (hi : i < items.size)
    (hn : (items.type i).isNode = true) : 1 ≤ items.nEdges g i := by
  by_cases hQ : items.type i = .Q
  · obtain ⟨h1, h2⟩ := hw.tree.Q_range hQ
    obtain ⟨e, rfl⟩ : ∃ e, i = edgeItem g e :=
      ⟨i - (1 + g.nv), by unfold edgeItem; (try simp only [ItemId] at *); omega⟩
    rcases hw.shapes.q_children e (by unfold edgeItem at h2; (try simp only [ItemId] at *); omega)
      with h | ⟨c, hcFV, _, hch⟩
    · exact Ghost.Items.one_le_nEdges_of_hasCap (by simp [Items.hasCap, hn, h])
    · rw [Items.nEdges_eq hn]
      have hc : c ∈ items.ch (edgeItem g e) := by rcases hch with h | ⟨v, _, h⟩ <;> simp [h]
      have hcn : (items.type c).isNode = true := by revert hcFV; cases items.type c <;> simp [NodeType.isNode]
      have hge : 1 + g.nv ≤ c := (hw.tree.isNode_iff (hw.tree.child_lt hc)).1 hcn
      have : 0 < (items.ch (edgeItem g e)).countP (· ≥ 1 + g.nv) :=
        List.countP_pos_iff.2 ⟨c, hc, by simpa using hge⟩
      omega
  · exact Ghost.Items.one_le_nEdges_of_hasCap (Items.hasCap_of_ne_Q hn hQ)

theorem layoutRot_size_node {g : Graph} {items : Items} (hw : items.WF g) {cur : ItemId} (hcur : cur < items.size)
    (hn : (items.type cur).isNode = true) (neSt : Nat) (edgeVes : List Nat) (mapRot : Nat → Array (Option Nat))
    (cv : Nat) (hmr : ∀ ve, (mapRot ve).size = 4) (hev : edgeVes.length = (items.ch cur).countP (· ≥ 1 + g.nv)) :
    (layoutRot (items.type cur) (items.nvList g cur).length neSt (neSt + items.nEdges g cur) edgeVes mapRot cv).size =
      4 * items.nEdges g cur := by
  obtain ⟨hL1, hLQI, hLS, hLR, hLcap, hLO⟩ := hw.layout_hyps hcur hn
  have h1n := hw.one_le_nEdges hcur hn
  have h1v := hw.one_le_nvList hcur hn
  have := layoutRot_size_of (items.type cur) (items.nvList g cur).length neSt (neSt + items.nEdges g cur) edgeVes
    mapRot cv hmr (by omega) ?_ ?_ ?_ ?_ ?_ ?_
  · rwa [Nat.add_sub_cancel_left] at this
  · intro h; rcases h with h | h <;> rw [h] at hn <;> simp [NodeType.isNode] at hn
  · intro _ _ hv; have := hL1 hv; omega
  · intro h; have := hLQI h; omega
  · intro hP hv
    have := hLcap (Items.hasCap_of_ne_Q hn (by rw [hP]; decide)) (Or.inr (Or.inr hP)); omega
  · intro hS hv; have := hLS hS; omega
  · intro h
    have hcap : items.hasCap cur = true := Items.hasCap_of_ne_Q hn (by rcases h with h | h <;> rw [h] <;> decide)
    rw [Items.nEdges_eq hn, Items.capCount, if_pos hcap, hev]; omega


theorem RowInv.init (g : Graph) (w : PlanarWalkState) : RowInv g w (PlanarRelabelState.init g w) :=
  ⟨rfl, rfl, rfl, rfl, by simp [PlanarRelabelState.init], rfl, rfl, rfl, rfl, rfl,
    by simp [PlanarRelabelState.init, RelabelState.init], rfl, rfl,
    fun n hn => by simp [PlanarRelabelState.init, RelabelState.init] at hn⟩

theorem Items.Tree.node_ge_of {g : Graph} {items : Items} (ht : items.Tree g) {i : ItemId}
    (hF : items.type i ≠ .F) (hV : items.type i ≠ .V) (hQ : items.type i ≠ .Q) : 1 + g.nv + g.ne ≤ i := by
  by_contra hlt
  by_cases h0 : i = 0
  · subst h0
    have := ht.root; unfold rootItem at this; exact hF this
  · by_cases hv : i < 1 + g.nv
    · have := ht.vert (i - 1) (by (try simp only [ItemId] at *); omega)
      unfold vertItem at this
      rw [show 1 + (i - 1) = i by (try simp only [ItemId] at *); omega] at this; exact hV this
    · have := ht.edge (i - (1 + g.nv)) (by (try simp only [ItemId] at *); omega)
      unfold edgeItem at this
      rw [show 1 + g.nv + (i - (1 + g.nv)) = i by (try simp only [ItemId] at *); omega] at this; exact hQ this

/-- A child of a descendant of the child `c` of `cur` is not a child of `cur`. -/
theorem Items.Tree.child_ne_of_below {g : Graph} {items : Items} (ht : items.Tree g) {cur c j c' d : ItemId}
    (hc : c ∈ items.ch cur) (hj : items.Below c j) (hc' : c' ∈ items.ch j) (hd : d ∈ items.ch cur) : c' ≠ d := by
  rintro rfl
  obtain ⟨par, _, huniq⟩ := ht.unique_parent c' (ht.child_pos hd) (ht.ch_lt _ _ hd)
  have : j = cur := (huniq j hc').trans (huniq cur hd).symm
  subst this
  exact ht.acyclic j (Relation.TransGen.head' hc hj)

theorem slot_eq_of {g : Graph} {c d : ItemId} {z : Nat} (hc : 1 + g.nv ≤ c) (hd : 1 + g.nv ≤ d) (hz : z < 4)
    (h : 4 * (d - (1 + g.nv)) ≤ 4 * (c - (1 + g.nv)) + z ∧ 4 * (c - (1 + g.nv)) + z < 4 * (d - (1 + g.nv)) + 4) :
    c = d := by
  (try simp only [ItemId] at *); omega

theorem edge_slot_of {g : Graph} (o z : Nat) :
    4 * (1 + g.nv + QE.edge o - (1 + g.nv)) ≤ o ∧ o < 4 * (1 + g.nv + QE.edge o - (1 + g.nv)) + 4 := by
  unfold QE.edge; omega

/-- The executable `orderedChildren` is `Items.ordered` with the state's vertex positions. -/
theorem orderedChildren_eq_ordered (g : Graph) (items : Items) (cur : ItemId) (hcur : cur < items.size) (n : Nat)
    (st : RelabelState) (hg : st.g = g) (hi : st.items = items) :
    ((RelabelM.orderedChildren items[cur]! n).run st).1 = Items.ordered g items cur n (fun v => st.vertPos[v]!) := by
  subst hg hi
  unfold RelabelM.orderedChildren Items.ordered Items.loc
  rw [Items.getElem!_type hcur, Items.getElem!_ch hcur]
  by_cases hR : Items.type st.items cur = NodeType.R
  · simp only [hR, bne_self_eq_false, Bool.false_eq_true, ↓reduceIte, ne_eq, not_true_eq_false]
    simp only [bind, pure, StateT.run, StateT.bind, StateT.pure, get, getThe, MonadStateOf.get, StateT.get,
      Items.getElem!_vs_all]
  · have e : (Items.type st.items cur != NodeType.R) = true := by simpa using hR
    rw [if_pos e, if_pos hR]; rfl

end Spqr
