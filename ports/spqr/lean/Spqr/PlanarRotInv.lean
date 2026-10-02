import Spqr.PlanarRelabel
import Spqr.PlanarRelabelProj
import Spqr.RelabelCost
import Spqr.LayoutSize

/-!
# The `neRotAdj` invariant of `planarRelabel`

`RotInv g s`: the planar relabel state `s` has `neRotAdj` equal to the concatenation of one
`layoutRot` block per numbered node (`rotBlock`), with `neBounds`/`nvBounds` one entry longer than
`types` and ending at the current `nodeEdges`/`nodeVerts` sizes. `planarRelabel_rotInv` (the fold
induction, `PlanarRotFold.lean`) preserves it; `PlanarRotSpec.planarRelabel_rot_spec` reads it off
the final state.

The `wp` calculus below is `RelabelM.wp` (`RelabelCost.lean`) for `PlanarRelabelM`, with `Frame`
tracking exactly the fields `RotInv` reads.
-/

namespace Spqr

namespace PlanarRelabelM

variable {α β : Type}

def wp (m : PlanarRelabelM α) (Q : α → PlanarRelabelState → Prop) (s : PlanarRelabelState) : Prop :=
  Q (m s).1 (m s).2

variable {s : PlanarRelabelState}

theorem wp_pure {Q : α → PlanarRelabelState → Prop} (a : α) :
    wp (pure a : PlanarRelabelM α) Q s = Q a s := rfl
theorem wp_bind {Q : α → PlanarRelabelState → Prop} (m : PlanarRelabelM β) (f : β → PlanarRelabelM α) :
    wp (m >>= f) Q s = wp m (fun b s' => wp (f b) Q s') s := rfl
theorem wp_map {Q : α → PlanarRelabelState → Prop} (f : β → α) (m : PlanarRelabelM β) :
    wp (f <$> m) Q s = wp m (fun b s' => Q (f b) s') s := rfl
theorem wp_ite {Q : α → PlanarRelabelState → Prop} (c : Prop) [Decidable c] (m₁ m₂ : PlanarRelabelM α) :
    wp (if c then m₁ else m₂) Q s = if c then wp m₁ Q s else wp m₂ Q s := by split <;> rfl

theorem liftR_eq (m : RelabelM α) (s : PlanarRelabelState) :
    liftR m s = ((m s.base).1, { s with base := (m s.base).2 }) := by
  unfold liftR; rcases hm : m s.base with ⟨a, b⟩; rfl
theorem liftAux_eq (m : StateM PlanarRelabelAux α) (s : PlanarRelabelState) :
    liftAux m s = ((m s.aux).1, { s with aux := (m s.aux).2 }) := by
  unfold liftAux; rcases hm : m s.aux with ⟨a, b⟩; rfl

theorem wp_liftR {Q : α → PlanarRelabelState → Prop} (m : RelabelM α) :
    wp (liftR m) Q s = Q (m s.base).1 { s with base := (m s.base).2 } := by
  unfold wp; rw [liftR_eq]
theorem wp_liftR_get {Q : RelabelState → PlanarRelabelState → Prop} :
    wp (liftR get) Q s = Q s.base s := by rw [wp_liftR]; rfl
theorem wp_liftR_item {Q : Item → PlanarRelabelState → Prop} (i : ItemId) :
    wp (liftR (RelabelM.item i)) Q s = Q s.base.items[i]! s := by rw [wp_liftR]; rfl
theorem wp_liftR_modify {Q : Unit → PlanarRelabelState → Prop} (f : RelabelState → RelabelState) :
    wp (liftR (modify f)) Q s = Q () { s with base := f s.base } := by rw [wp_liftR]; rfl
theorem wp_getAux {Q : PlanarRelabelAux → PlanarRelabelState → Prop} : wp getAux Q s = Q s.aux s := rfl
theorem wp_modifyAux {Q : Unit → PlanarRelabelState → Prop} (f : PlanarRelabelAux → PlanarRelabelAux) :
    wp (modifyAux f) Q s = Q () { s with aux := f s.aux } := rfl
theorem wp_liftR_orderedChildren {Q : List ItemId → PlanarRelabelState → Prop} (it : Item) (n : Nat)
    (h : ∀ children, children = ((RelabelM.orderedChildren it n).run s.base).1 → Q children s) :
    wp (liftR (RelabelM.orderedChildren it n)) Q s := by
  have hs : (RelabelM.orderedChildren it n s.base).2 = s.base := by
    unfold RelabelM.orderedChildren; split <;> rfl
  rw [wp_liftR, hs]; exact h _ rfl

theorem wp_forIn_inv {α β : Type} (l : List α) (init : β) (f : α → β → PlanarRelabelM (ForInStep β))
    (Q : β → PlanarRelabelState → Prop) (I : List α → β → PlanarRelabelState → Prop) (s : PlanarRelabelState)
    (hinit : I l init s)
    (hstep : ∀ a ∈ l, ∀ rest b s', I (a :: rest) b s' →
      wp (f a b) (fun r s'' => ∃ b', r = .yield b' ∧ I rest b' s'') s')
    (hfin : ∀ b s', I [] b s' → Q b s') : wp (forIn l init f) Q s := by
  induction l generalizing init s with
  | nil => exact hfin _ _ hinit
  | cons a l ih =>
    obtain ⟨b', hr, hI⟩ := hstep a (List.mem_cons_self ..) l init s hinit
    have e : (forIn (a :: l) init f) s = (forIn l b' f) ((f a init) s).2 := by
      rw [List.forIn_cons]
      show (match ((f a init) s).1 with
        | ForInStep.done b => pure b | ForInStep.yield b => forIn l b f) ((f a init) s).2 = _
      rw [hr]
    unfold wp; rw [e]
    exact ih b' _ hI fun a' ha' => hstep a' (List.mem_cons_of_mem _ ha')

/-- The fields `RotInv` reads. -/
structure Frame (s s' : PlanarRelabelState) : Prop where
  g : s'.base.g = s.base.g
  types : s'.base.types = s.base.types
  neBounds : s'.base.neBounds = s.base.neBounds
  nvBounds : s'.base.nvBounds = s.base.nvBounds
  nodeEdges : s'.base.nodeEdges.size = s.base.nodeEdges.size
  nodeVerts : s'.base.nodeVerts.size = s.base.nodeVerts.size
  neRotAdj : s'.aux.neRotAdj = s.aux.neRotAdj

theorem Frame.refl (s : PlanarRelabelState) : Frame s s := ⟨rfl, rfl, rfl, rfl, rfl, rfl, rfl⟩
theorem Frame.trans {s₁ s₂ s₃ : PlanarRelabelState} (h₁ : Frame s₁ s₂) (h₂ : Frame s₂ s₃) : Frame s₁ s₃ :=
  ⟨h₂.g.trans h₁.g, h₂.types.trans h₁.types, h₂.neBounds.trans h₁.neBounds, h₂.nvBounds.trans h₁.nvBounds,
    h₂.nodeEdges.trans h₁.nodeEdges, h₂.nodeVerts.trans h₁.nodeVerts, h₂.neRotAdj.trans h₁.neRotAdj⟩

theorem wp_frame {Q : α → PlanarRelabelState → Prop} (m : PlanarRelabelM α)
    (hm : ∀ s, Frame s (m s).2) (hk : ∀ a s', Frame s s' → Q a s') : wp m Q s :=
  hk _ _ (hm s)

/-- Abstract the state after `m` behind the fields `RotInv` reads. -/
theorem wp_abs {Q : Unit → PlanarRelabelState → Prop} (m : PlanarRelabelM Unit)
    (hk : ∀ s', s'.base.g = (m s).2.base.g → s'.base.types = (m s).2.base.types →
      s'.base.neBounds = (m s).2.base.neBounds → s'.base.nvBounds = (m s).2.base.nvBounds →
      s'.base.nodeEdges.size = (m s).2.base.nodeEdges.size →
      s'.base.nodeVerts.size = (m s).2.base.nodeVerts.size →
      s'.aux.neRotAdj = (m s).2.aux.neRotAdj → Q () s') : wp m Q s :=
  hk _ rfl rfl rfl rfl rfl rfl rfl

theorem wp_jp {Q : Unit → PlanarRelabelState → Prop} (m k : PlanarRelabelM Unit) (s : PlanarRelabelState)
    (h : ∀ s, ∃ s', Frame s s' ∧ (m s).2 = (k s').2)
    (hk : ∀ s', Frame s s' → wp k Q s') : wp m Q s := by
  obtain ⟨s', hf, he⟩ := h s
  have := hk s' hf
  unfold wp at this ⊢
  exact he ▸ this

theorem arm_frame (m k : PlanarRelabelM Unit) (s : PlanarRelabelState) (hm : Frame s (m s).2) :
    ∃ s', Frame s s' ∧ ((m >>= fun _ => k) s).2 = (k s').2 :=
  ⟨(m s).2, hm, rfl⟩

theorem arm_id (k : PlanarRelabelM Unit) (s : PlanarRelabelState) :
    ∃ s', Frame s s' ∧ (k s).2 = (k s').2 :=
  ⟨s, Frame.refl s, rfl⟩

theorem Frame.liftR_modify (f : RelabelState → RelabelState) (s : PlanarRelabelState)
    (h : Frame s { s with base := f s.base }) : Frame s (liftR (modify f) s).2 := by
  rw [liftR_eq]; exact h

theorem Frame.bind {m : PlanarRelabelM α} {k : α → PlanarRelabelM β} (s : PlanarRelabelState)
    (hm : Frame s (m s).2) (hk : ∀ a s', Frame s s' → Frame s' (k a s').2) :
    Frame s ((m >>= k) s).2 :=
  hm.trans (hk _ _ hm)

theorem Frame.pure (a : α) (s : PlanarRelabelState) : Frame s ((pure a : PlanarRelabelM α) s).2 := Frame.refl s
theorem Frame.getAux (s : PlanarRelabelState) : Frame s (getAux s).2 := Frame.refl s

/-- A `StateM` loop body that preserves `P` on the state preserves it over the loop. -/
theorem forIn_list_state_pres {σ β γ : Type} (P : σ → Prop) (l : List γ) (init : β)
    (f : γ → β → StateM σ (ForInStep β)) (hf : ∀ a b st, P st → P ((f a b) st).2) :
    ∀ st, P st → P ((forIn l init f) st).2 := by
  induction l generalizing init with
  | nil => intro st h; exact h
  | cons a l ih =>
    intro st h
    rw [List.forIn_cons]
    show P ((match ((f a init) st).1 with
      | ForInStep.done b => pure b | ForInStep.yield b => forIn l b f) ((f a init) st).2).2
    have := hf a init st h
    cases hx : ((f a init) st).1 with
    | done b => exact this
    | yield b => exact ih _ _ this

theorem StateT_bind_pure_snd {σ α β : Type} (m : StateM σ α) (b : β) (st : σ) :
    ((m >>= fun _ => pure b) st).2 = (m st).2 := by
  show (match m st with | (_, s') => (b, s')).2 = _
  rcases m st with ⟨x, y⟩; rfl

theorem setupNode_neRotAdj (g : Graph) (type : NodeType) (cur : ItemId) (s : PlanarRelabelState) :
    (setupNode g type cur s).2.aux.neRotAdj = s.aux.neRotAdj := by
  unfold setupNode; rw [liftAux_eq]
  generalize s.aux = a
  show ((match type with
    | .S | .P | .R => _
    | _ => _ : StateM PlanarRelabelAux Bool) a).2.neRotAdj = a.neRotAdj
  split
  all_goals first
    | rfl
    | (show ((match a.nodePlanarity[cur - (1 + g.nv + g.ne)]! with
          | .planar m => _ | _ => _ : StateM PlanarRelabelAux Bool) a).2.neRotAdj = a.neRotAdj
       generalize a.nodePlanarity[cur - (1 + g.nv + g.ne)]! = p
       cases p
       all_goals first
         | rfl
         | (rw [StateT_bind_pure_snd, Std.Legacy.Range.forIn_eq_forIn_range']
            refine forIn_list_state_pres (fun a' => a'.neRotAdj = a.neRotAdj) _ _ _ ?_ a rfl
            intro _ _ st h; exact h))

theorem applyFlips_neRotAdj (g : Graph) (it : Item) (flips : List Bool) (s : PlanarRelabelState) :
    (applyFlips g it flips s).2.aux.neRotAdj = s.aux.neRotAdj := by
  unfold applyFlips; rw [liftAux_eq]
  generalize s.aux = a
  try rw [StateT_bind_pure_snd]
  refine forIn_list_state_pres (fun a' => a'.neRotAdj = a.neRotAdj) _ _ _ ?_ a rfl
  intro p b st h
  obtain ⟨c, flip⟩ := p
  dsimp only
  first | (split <;> exact h) | (split_ifs <;> exact h)

theorem Frame.setupNode (g : Graph) (type : NodeType) (cur : ItemId) (s : PlanarRelabelState) :
    Frame s (setupNode g type cur s).2 := by
  have hb := AuxOnlyR.setupNode g type cur s
  have ha := setupNode_neRotAdj g type cur s
  exact ⟨by rw [hb], by rw [hb], by rw [hb], by rw [hb], by rw [hb], by rw [hb], ha⟩

theorem Frame.applyFlips (g : Graph) (it : Item) (flips : List Bool) (s : PlanarRelabelState) :
    Frame s (applyFlips g it flips s).2 := by
  have hb := AuxOnlyR.applyFlips g it flips s
  have ha := applyFlips_neRotAdj g it flips s
  exact ⟨by rw [hb], by rw [hb], by rw [hb], by rw [hb], by rw [hb], by rw [hb], ha⟩

theorem Frame.modifyAux (f : PlanarRelabelAux → PlanarRelabelAux) (s : PlanarRelabelState)
    (h : (f s.aux).neRotAdj = s.aux.neRotAdj) : Frame s (modifyAux f s).2 :=
  ⟨rfl, rfl, rfl, rfl, rfl, rfl, h⟩

end PlanarRelabelM

/-! ## Blocks -/

/-- `f 0 ++ f 1 ++ ⋯ ++ f (k-1)`. -/
def concatBlocks {α : Type} (f : Nat → Array α) (k : Nat) : Array α :=
  (List.range k).foldl (fun acc n => acc ++ f n) #[]

namespace concatBlocks

variable {α : Type} {f f' : Nat → Array α}

theorem zero : concatBlocks f 0 = #[] := rfl
theorem succ (k : Nat) : concatBlocks f (k + 1) = concatBlocks f k ++ f k := by
  simp [concatBlocks, List.range_succ]

theorem congr {k : Nat} (h : ∀ n, n < k → f n = f' n) : concatBlocks f k = concatBlocks f' k := by
  induction k with
  | zero => rfl
  | succ k ih => rw [succ, succ, ih (fun n hn => h n (by omega)), h k (by omega)]

theorem size_le {n k : Nat} (h : n ≤ k) : (concatBlocks f n).size ≤ (concatBlocks f k).size := by
  induction k with
  | zero => cases Nat.le_zero.mp h; exact le_rfl
  | succ k ih =>
    rcases Nat.lt_succ_iff_lt_or_eq.mp (Nat.lt_succ_of_le h) with h' | h'
    · rw [succ, Array.size_append]; exact le_trans (ih (by omega)) (Nat.le_add_right _ _)
    · subst h'; exact le_rfl

theorem get [Inhabited α] {n k j : Nat} (hn : n < k) (hj : j < (f n).size) :
    (concatBlocks f k)[(concatBlocks f n).size + j]! = (f n)[j]! := by
  induction k with
  | zero => omega
  | succ k ih =>
    rw [succ]
    rcases Nat.lt_succ_iff_lt_or_eq.mp hn with h | h
    · have hi : (concatBlocks f n).size + j < (concatBlocks f k).size := by
        have := size_le (f := f) (show n + 1 ≤ k by omega)
        rw [succ, Array.size_append] at this; omega
      rw [← ih h, getElem!_pos (concatBlocks f k ++ f k) _ (by simp; omega), getElem!_pos (concatBlocks f k) _ hi,
        Array.getElem_append_left hi]
    · subst h
      rw [getElem!_pos (concatBlocks f n ++ f n) _ (by simp; omega), getElem!_pos (f n) _ hj,
        Array.getElem_append_right (Nat.le_add_right _ _)]
      simp

end concatBlocks

theorem Array.getElem!_push_lt' {α : Type} [Inhabited α] (a : Array α) (x : α) (i : Nat) (h : i < a.size) :
    (a.push x)[i]! = a[i]! := by
  rw [getElem!_pos (a.push x) i (by simp; omega), getElem!_pos a i h, Array.getElem_push_lt]

theorem Array.getElem!_push_eq' {α : Type} [Inhabited α] (a : Array α) (x : α) :
    (a.push x)[a.size]! = x := by
  rw [getElem!_pos (a.push x) a.size (by simp), Array.getElem_push_eq]

/-- The `neRotAdj` block of node `n`. -/
def rotBlock (ne : Nat) (types : Array NodeType) (nvB neB : Array Nat) (ev : Nat → List Nat)
    (mr : Nat → Nat → Array (Option Nat)) (n : Nat) : Array (Option Nat) :=
  layoutRot types[n]! (nvB[n + 1]! - nvB[n]!) neB[n]! neB[n + 1]! (ev n) (mr n) (2 * ne)

/-- Invariant of `planarRelabel` between calls. -/
structure RotInv (g : Graph) (s : PlanarRelabelState) : Prop where
  g_eq : s.base.g = g
  ne_size : s.base.neBounds.size = s.base.types.size + 1
  nv_size : s.base.nvBounds.size = s.base.types.size + 1
  ne_zero : s.base.neBounds[0]! = 0
  ne_last : s.base.neBounds[s.base.types.size]! = s.base.nodeEdges.size
  nv_last : s.base.nvBounds[s.base.types.size]! = s.base.nodeVerts.size
  blocks : ∃ (ev : Nat → List Nat) (mr : Nat → Nat → Array (Option Nat)),
    (∀ n ve, (mr n ve).size = 4) ∧
    (∀ n, n < s.base.types.size → s.base.types[n]! = .R →
      (ev n).length + 1 = s.base.neBounds[n + 1]! - s.base.neBounds[n]!) ∧
    s.aux.neRotAdj = concatBlocks (rotBlock g.ne s.base.types s.base.nvBounds s.base.neBounds ev mr)
      s.base.types.size

theorem RotInv.frame {g : Graph} {s s' : PlanarRelabelState} (hf : PlanarRelabelM.Frame s s')
    (h : RotInv g s) : RotInv g s' := by
  obtain ⟨hg, h1, h2, h3, h4, h5, ev, mr, hmr, hR, hb⟩ := h
  refine ⟨hf.g.trans hg, ?_, ?_, ?_, ?_, ?_, ev, mr, hmr, ?_, ?_⟩ <;>
    simp only [hf.types, hf.neBounds, hf.nvBounds, hf.nodeEdges, hf.nodeVerts, hf.neRotAdj] <;>
    assumption

theorem RotInv.init (g : Graph) (w : PlanarWalkState) : RotInv g (PlanarRelabelState.init g w) :=
  ⟨rfl, rfl, rfl, rfl, rfl, rfl, fun _ => [], fun _ _ => Array.replicate 4 none,
    fun _ _ => by simp, fun n hn => by simp [PlanarRelabelState.init, RelabelState.init] at hn, rfl⟩

end Spqr
