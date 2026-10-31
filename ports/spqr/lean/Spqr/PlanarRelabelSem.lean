import Spqr.PlanarRotInv
import Spqr.Proofs.PlanarWalkFacts
import Spqr.RelabelLoop

/-!
# `planarRelabel`: semantics of the planar steps

`setupNode` is `capLinked`, `applyFlips` is `flipped` (both from `Proofs/PlanarWalkFacts.lean`), and
the R-node initialization of `rotEdgeNe` sends each child virtual edge to `neSt + 1 + idxOf` and
the cap to `neSt`.
-/

namespace Spqr
open PlanarRelabelM

theorem forIn_range4_list {m : Type → Type} [Monad m] {β : Type} (init : β) (f : Nat → β → m (ForInStep β)) :
    forIn [0 : 4] init f = forIn [0, 1, 2, 3] init f := by
  rw [Std.Legacy.Range.forIn_eq_forIn_range']
  have : List.range' (Std.Legacy.Range.start [0:4]) (Std.Legacy.Range.size [0:4]) (Std.Legacy.Range.step [0:4]) = [0, 1, 2, 3] := by
    decide
  rw [this]

theorem capLinked_eq (g : Graph) (q : Qem) (m : Array Nat) :
    capLinked g q m =
      (fun q => (q.set! (8 * g.ne + 3) (some m[3]!)).set! m[3]! (some (8 * g.ne + 3)))
      ((fun q => (q.set! (8 * g.ne + 2) (some m[2]!)).set! m[2]! (some (8 * g.ne + 2)))
      ((fun q => (q.set! (8 * g.ne + 1) (some m[1]!)).set! m[1]! (some (8 * g.ne + 1)))
      ((fun q => (q.set! (8 * g.ne + 0) (some m[0]!)).set! m[0]! (some (8 * g.ne + 0))) q))) := by
  simp only [capLinked, capVe, List.range_succ, List.range_zero, List.nil_append, List.cons_append,
    List.foldl_cons, List.foldl_nil]
  have : ∀ s, 4 * (2 * g.ne) + s = 8 * g.ne + s := by intro s; omega
  simp only [this]

theorem forIn4_run {σ : Type} (f : Nat → σ → σ) (a : σ) :
    (forIn (m := StateM σ) [0, 1, 2, 3] PUnit.unit
      (fun s _ => StateT.bind (modify (f s)) (fun _ => StateT.pure (ForInStep.yield PUnit.unit)))) a =
      (PUnit.unit, f 3 (f 2 (f 1 (f 0 a)))) := by
  rfl

/-- `setupNode` at an S/P/R item recorded planar: links the cap (`capLinked`), returns `true`. -/
theorem setupNode_planar (g : Graph) (type : NodeType) (cur : ItemId) (s : PlanarRelabelState)
    (hty : type = .S ∨ type = .P ∨ type = .R) (m : Array Nat)
    (hm : s.aux.nodePlanarity[cur - (1 + g.nv + g.ne)]! = .planar m) :
    setupNode g type cur s = (true, { s with aux := { s.aux with qem := capLinked g s.aux.qem m } }) := by
  unfold setupNode
  rw [liftAux_eq]
  have key : ∀ (x : StateM PlanarRelabelAux Bool), x s.aux = (true, { s.aux with qem := capLinked g s.aux.qem m }) →
      ((x s.aux).1, { s with aux := (x s.aux).2 }) = (true, { s with aux := { s.aux with qem := capLinked g s.aux.qem m } }) := by
    intro x hx; rw [hx]
  apply key
  rcases hty with rfl | rfl | rfl <;>
  · simp only [bind, StateT.bind, get, getThe, MonadStateOf.get, StateT.get, pure, StateT.pure, hm]
    rw [forIn_range4_list, capLinked_eq,
      forIn4_run (fun s (a : PlanarRelabelAux) => { a with qem := (a.qem.set! (8 * g.ne + s) (some m[s]!)).set! m[s]! (some (8 * g.ne + s)) })]

/-- `setupNode` at an S/P/R item not recorded planar: no change, returns `false`. -/
theorem setupNode_nonplanar (g : Graph) (type : NodeType) (cur : ItemId) (s : PlanarRelabelState)
    (hty : type = .S ∨ type = .P ∨ type = .R)
    (hm : ∀ m, s.aux.nodePlanarity[cur - (1 + g.nv + g.ne)]! ≠ .planar m) :
    setupNode g type cur s = (false, s) := by
  unfold setupNode
  rw [liftAux_eq]
  have key : ∀ (x : StateM PlanarRelabelAux Bool), x s.aux = (false, s.aux) →
      ((x s.aux).1, { s with aux := (x s.aux).2 }) = (false, s) := by
    intro x hx; rw [hx]
  apply key
  rcases hty with rfl | rfl | rfl <;>
  · simp only [bind, StateT.bind, get, getThe, MonadStateOf.get, StateT.get, pure, StateT.pure]

theorem setupNode_other (g : Graph) (type : NodeType) (cur : ItemId) (s : PlanarRelabelState)
    (hty : type ≠ .S ∧ type ≠ .P ∧ type ≠ .R) :
    setupNode g type cur s = (true, s) := by
  unfold setupNode
  rw [liftAux_eq]
  cases type <;> first | exact absurd rfl hty.1 | exact absurd rfl hty.2.1 | exact absurd rfl hty.2.2 | rfl

/-- `rotEdgeNe` after the R-node initialization: node edge `neSt + 1 + idxOf ve` for each child virtual
edge, `neSt` for the cap. -/
theorem foldl_set!_idxOf (vl : List Nat) (n : Nat) (a : Array Nat) (hnd : vl.Nodup)
    (hlt : ∀ v ∈ vl, v < a.size) {v : Nat} (hv : v ∈ vl) :
    ((vl.zipIdx n).foldl (fun a x => a.set! x.1 x.2) a)[v]! = n + vl.idxOf v := by
  obtain ⟨h1, h2⟩ := Ghost.foldl_set!_posOK vl n a hlt v hv
  have hlt' : ((vl.zipIdx n).foldl (fun a x => a.set! x.1 x.2) a)[v]! - n < vl.length := by
    rcases Nat.lt_or_ge (((vl.zipIdx n).foldl (fun a x => a.set! x.1 x.2) a)[v]! - n) vl.length with h | h
    · exact h
    · rw [List.getElem?_eq_none h] at h2; cases h2
  have h3 : vl[((vl.zipIdx n).foldl (fun a x => a.set! x.1 x.2) a)[v]! - n] = v := by
    rw [List.getElem?_eq_getElem hlt'] at h2; exact Option.some.inj h2
  have := hnd.idxOf_getElem _ hlt'
  rw [h3] at this
  omega

theorem rotEdgeNe_init (edgeVes : List Nat) (neSt : Nat) (r : Array Nat) (ne : Nat) (hnd : edgeVes.Nodup)
    (hlt : ∀ v ∈ edgeVes, v < r.size) (hcap : 2 * ne < r.size) (hne : 2 * ne ∉ edgeVes) :
    let r' := ((edgeVes.zipIdx (neSt + 1)).foldl (init := r) (fun r (ve, ne) => r.set! ve ne)).set! (2 * ne) neSt
    r'[2 * ne]! = neSt ∧ ∀ v ∈ edgeVes, r'[v]! = neSt + 1 + edgeVes.idxOf v := by
  intro r'
  refine ⟨?_, ?_⟩
  · show (((edgeVes.zipIdx (neSt + 1)).foldl (fun r x => r.set! x.1 x.2) r).set! (2 * ne) neSt)[2 * ne]! = _
    apply Ghost.set!_get!_self
    rw [Ghost.foldl_set!_size]; exact hcap
  · intro v hv
    have hne' : 2 * ne ≠ v := fun h => hne (h ▸ hv)
    show (((edgeVes.zipIdx (neSt + 1)).foldl (fun r x => r.set! x.1 x.2) r).set! (2 * ne) neSt)[v]! = _
    rw [Ghost.set!_get!_ne _ _ hne']
    exact foldl_set!_idxOf edgeVes (neSt + 1) r hnd hlt hv

theorem applyFlips_eq (g : Graph) (it : Item) (flips : List Bool) (s : PlanarRelabelState) :
    applyFlips g it flips s = ((), { s with aux := { s.aux with qem := flipped g it.ch flips s.aux.qem } }) := by
  unfold applyFlips flipped
  rw [liftAux_eq]
  have h : ∀ (l : List (ItemId × Bool)) (a : PlanarRelabelAux),
      (forIn l () fun (x : ItemId × Bool) (_ : Unit) => (do
        if x.1 ≥ 1 + g.nv && x.2 then
          modify fun (a : PlanarRelabelAux) => { a with qem := (a.qem.swapIfInBounds (4 * (x.1 - (1 + g.nv)) + 0) (4 * (x.1 - (1 + g.nv)) + 1)).swapIfInBounds (4 * (x.1 - (1 + g.nv)) + 2) (4 * (x.1 - (1 + g.nv)) + 3) }
        pure (ForInStep.yield ()) : StateM PlanarRelabelAux (ForInStep Unit))) a =
      ((), { a with qem := l.foldl (fun q (x : ItemId × Bool) => if x.1 ≥ 1 + g.nv && x.2 then
        (q.swapIfInBounds (4 * (x.1 - (1 + g.nv)) + 0) (4 * (x.1 - (1 + g.nv)) + 1)).swapIfInBounds
          (4 * (x.1 - (1 + g.nv)) + 2) (4 * (x.1 - (1 + g.nv)) + 3) else q) a.qem }) := by
    intro l
    induction l with
    | nil => intro a; rfl
    | cons x l ih =>
      intro a
      rw [List.forIn_cons, List.foldl_cons]
      by_cases hc : (decide (x.1 ≥ 1 + g.nv) && x.2) = true
      · simp only [hc, ite_true]
        exact ih _
      · simp only [hc]
        exact ih _
  exact (congrArg (fun p : Unit × PlanarRelabelAux => (p.1, { s with aux := p.2 })) (h (it.ch.zip flips) s.aux))

end Spqr
