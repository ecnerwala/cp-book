import Spqr.PlanarInv
import Spqr.PlanarEmbed
import Spqr.PlanarShape
import Spqr.PieceSep

/-!
# The gluing invariant of `planarEmbed`

`planarEmbed` runs `embedItem` over the items in reverse preorder. `GluedUpTo i s` is the state
invariant after the items `≥ i` have been processed; the per-item steps are stated separately
(`embedItem_step_*`, under `SpqrTree.WF` and `SpqrTree.ChildShape` of the tree — they are false
for an arbitrary `PlanarSpqrTree`), the fold is proved (`forM_reverse_range_inv`), and the steps
are assembled in `PlanarEmbedFold.lean`.
-/

namespace Spqr

namespace PlanarSpqrTree

variable (t : PlanarSpqrTree)

/-- Original edges of the `Q` items in the subtree of item `i` (items `i ..< subtreeEnd i`). -/
def edgesBelow (i : Nat) : List Nat :=
  (List.range' i (t.subtreeEnd[i]! - i)).filterMap fun j =>
    if t.types[j]! == .Q then t.origId[j]! else none

/-- The piece of `g` below item `i` (its terminals are irrelevant here). -/
def pieceBelow (g : Graph) (i : Nat) : Piece :=
  ⟨t.edgesBelow i, fun e => g.edges[e]!, g.nv, 0, 0⟩

/-- Item `j` is processed (`i ≤ j`) and its parent is not. -/
def Maximal (i j : Nat) : Prop :=
  i ≤ j ∧ j < t.size ∧ ∀ p, t.par[j]? = some (some p) → p < i

/-- The exposed quarter-edges of item `j` in state `s`. -/
def EmbedState.exposedAt (s : EmbedState) (j q : Nat) : Prop :=
  ∃ k : Nat, s.outerE[j]?.bind (fun a : Array (Option Nat) => a[k]?) = some (some q)

/-- The gluing invariant after the items `≥ i` have been processed: for every maximal processed
item `j`, the glued rotation restricted to the real edges below `j` agrees with a planar
embedding `ρ` of that piece; exactly the exposed ends `outerE[j]` are still unset, and each
exposed pair `outerE[j][2k], outerE[j][2k+1]` is a facing pair of `ρ` (so that a `link` is a
transposition conjugation of `ρ`). Quarter-edges of edges not below any processed item are unset,
and unprocessed items have no exposed ends yet. -/
structure GluedPieces (g : Graph) (i : Nat) (s : EmbedState) : Prop where
  rot_size : s.rotAdj.size = 4 * t.ne
  outer_size : s.outerE.size = t.size
  outer_unprocessed : ∀ j, j < i → ∀ q, ¬ s.exposedAt j q
  unset : ∀ q, (∀ j, i ≤ j → j < t.size → QE.edge q ∉ t.edgesBelow j) →
    s.rotAdj[q]? = none ∨ s.rotAdj[q]? = some none
  piece : ∀ j, t.Maximal i j → ∃ ρ : RotationSystem,
    IsPlanarEmbedding (t.pieceBelow g j).es g.nv ρ ∧
    (∀ q r, (t.pieceBelow g j).Mem q → s.rotAdj[q]? = some (some r) →
      ∃ lq lr, (t.pieceBelow g j).loc q = some lq ∧ (t.pieceBelow g j).loc r = some lr ∧
        ρ.get lq = some lr) ∧
    (∀ q, (t.pieceBelow g j).Mem q → (s.rotAdj[q]? = some none ↔ s.exposedAt j q)) ∧
    (∀ k, k < 2 →
      (∀ a, s.outerE[j]?.bind (fun o => o[2 * k]?) = some (some a) →
        ∃ b, s.outerE[j]?.bind (fun o => o[2 * k + 1]?) = some (some b) ∧
          ∃ la lb, (t.pieceBelow g j).loc a = some la ∧ (t.pieceBelow g j).loc b = some lb ∧
            ρ.get la = some lb) ∧
      (∀ b, s.outerE[j]?.bind (fun o => o[2 * k + 1]?) = some (some b) →
        ∃ a, s.outerE[j]?.bind (fun o => o[2 * k]?) = some (some a)))

structure GluedSlots (g : Graph) (i : Nat) (s : EmbedState) : Prop extends t.GluedPieces g i s where
  outer_row_size : ∀ j, j < t.size → s.outerE[j]!.size = 4
  outer_slots : ∀ j k q, s.outerE[j]?.bind (fun o => o[k]?) = some (some q) →
    k < 4 ∧ t.toSpqrTree.type j ≠ .F ∧
      ∀ p, t.toSpqrTree.parent j = some p →
        (t.toSpqrTree.type p = .F ∨ t.toSpqrTree.type p = .V) → k < 2

structure GluedAttachments (g : Graph) (i : Nat) (s : EmbedState) : Prop extends t.GluedSlots g i s where
  outer_at_vertex : ∀ j p v q, t.toSpqrTree.parent j = some p → t.toSpqrTree.type p = .V →
    t.origId[p]! = some v → s.exposedAt j q → QE.vert g.edges.toList q = some v

structure GluedOriented (g : Graph) (i : Nat) (s : EmbedState) : Prop
    extends t.GluedAttachments g i s where
  outer_dir : ∀ (j k q : Nat), s.outerE[j]?.bind (fun o => o[k]?) = some (some q) → q % 2 = k % 2
  outer_present : ∀ j, i ≤ j → j < t.size → ∀ p,
    t.toSpqrTree.parent j = some p → t.toSpqrTree.type p = .V →
    t.edgesBelow j ≠ [] → ∃ q, s.exposedAt j q

structure GluedUpTo (g : Graph) (i : Nat) (s : EmbedState) : Prop
    extends t.GluedOriented g i s where
  outer_vertex : ∀ j, j < t.size → t.toSpqrTree.type j = .V → ∀ v, t.origId[j]! = some v →
    (∀ q, s.exposedAt j q → QE.vert g.edges.toList q = some v) ∧
    (∀ k q, s.outerE[j]?.bind (fun o => o[k]?) = some (some q) → k < 2) ∧
    (i ≤ j → t.edgesBelow j ≠ [] → ∃ q, s.exposedAt j q)

/-- The initial state of `planarEmbed`. -/
def initState : EmbedState :=
  ⟨Array.replicate (4 * t.ne) none, Array.replicate t.size (Array.replicate 4 none)⟩

theorem gluedUpTo_init (g : Graph) : t.GluedUpTo g t.size t.initState where
  rot_size := by simp [initState]
  outer_size := by simp [initState]
  outer_unprocessed := by
    intro j _ q ⟨k, hk⟩
    simp only [initState, Array.getElem?_replicate] at hk
    split at hk <;> simp only [Option.bind_some, Option.bind_none, Array.getElem?_replicate] at hk <;>
      split at hk <;> simp at hk
  unset := by
    intro q _
    simp only [initState, Array.getElem?_replicate]
    split <;> simp
  piece := by
    intro j hj
    obtain ⟨h1, h2, _⟩ := hj
    omega
  outer_row_size := by
    intro j hj
    simp [initState, getElem!_pos, hj]
  outer_slots := by
    intro j k q hk
    by_cases hj : j < t.size <;> by_cases hk' : k < 4 <;>
      simp [initState, Array.getElem?_replicate, hj, hk'] at hk
  outer_at_vertex := by
    intro j p v q _ _ _ ⟨k, hk⟩
    by_cases hj : j < t.size <;> by_cases hk' : k < 4 <;>
      simp [initState, Array.getElem?_replicate, hj, hk'] at hk
  outer_dir := by
    intro j k q hk
    by_cases hj : j < t.size <;> by_cases hk' : k < 4 <;>
      simp [initState, hj, hk'] at hk
  outer_present := by
    intro j hj hjs
    omega
  outer_vertex := by
    intro j hj _ v _
    refine ⟨?_, ?_, ?_⟩
    · intro q ⟨k, hk⟩
      by_cases hk' : k < 4 <;> simp [initState, hj, hk'] at hk
    · intro k q hk
      by_cases hk' : k < 4 <;> simp [initState, hj, hk'] at hk
    · intro hge; omega

/-- `Q` step: the real edge `origId i` is added with its four quarter-edges; the two `I`/`O`
children (the loops / blocks at its endpoints) are 1-summed at the endpoints. Admitted. -/
theorem embedItem_step_Q (g : Graph) (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape)
    (hrep : t.toSpqrTree.Represents g) (hsep : t.toSpqrTree.PieceSep g)
    (i : Nat) (hi : i < t.size) (hty : t.types[i]! = .Q)
    (s : EmbedState) (h : t.GluedUpTo g (i + 1) s) :
    t.GluedUpTo g i ((t.embedItem i).run s).2 := by
  sorry

/-- `S`/`P`/`R` step: the node's local rotation (`nodePlanar_sound`) is 2-summed with each
child's piece through the twin virtual edge (`twoSum_planar`), and the cap's quarter-edges become
the exposed ends. Admitted. -/
theorem embedItem_step_node (g : Graph) (hwf : t.toSpqrTree.WF) (hsh : t.toSpqrTree.ChildShape)
    (hrep : t.toSpqrTree.Represents g) (hsep : t.toSpqrTree.PieceSep g)
    (i : Nat) (hi : i < t.size)
    (hty : t.types[i]! = .S ∨ t.types[i]! = .P ∨ t.types[i]! = .R)
    (hall : t.nodePlanar.all id = true)
    (s : EmbedState) (h : t.GluedUpTo g (i + 1) s) :
    t.GluedUpTo g i ((t.embedItem i).run s).2 := by
  sorry

/-- The reverse-preorder fold of `embedItem` preserves any step-invariant. -/
theorem forM_reverse_range_inv (P : Nat → EmbedState → Prop) (n : Nat)
    (step : ∀ i s, i < n → P (i + 1) s → P i ((t.embedItem i).run s).2)
    (s₀ : EmbedState) (h : P n s₀) :
    P 0 (((List.range n).reverse.forM t.embedItem).run s₀).2 := by
  induction n generalizing s₀ with
  | zero =>
    simp only [List.range_zero, List.reverse_nil, List.forM_eq_forM, List.forM_nil, StateT.run_pure]
    exact h
  | succ n ih =>
    rw [List.range_succ, List.reverse_append, List.reverse_singleton, List.singleton_append]
    simp only [List.forM_eq_forM, List.forM_cons, StateT.run_bind]
    exact ih (fun i s hi => step i s (by omega)) _ (step n s₀ (by omega) h)

end PlanarSpqrTree

end Spqr
