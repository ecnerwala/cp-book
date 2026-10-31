import Spqr.PlanarEmbedVLoop

/-!
# Unfolding the `Q` branch of `embedItem`

The three child shapes of a `Q` item (`[]`, `[O]`, `c :: w :: _`) as explicit state functions.
-/

namespace Spqr.PlanarSpqrTree

open EmbedM

variable (t : PlanarSpqrTree)

/-- The four quarter-edges of the `Q` item `i`, in `qes` order. -/
def qes (i : Nat) (k : Nat) : Option Nat :=
  let e := (t.origId[i]!).getD 0
  let flip := t.edgeFlipped[e]!
  #[some (QE.mk e flip false), some (QE.mk e flip true),
    some (QE.mk e (!flip) false), some (QE.mk e (!flip) true)][k]!

/-- The `Q` step's gluing onto the capped child `c`: `q1`, `q2` are linked to `c`'s slots `0`
and `3`, and `c`'s slots `1`, `2` become the new ends (an `I` child is skipped). -/
def qUpper (c : Nat) (q1 q2 : Option Nat) (s : EmbedState) :
    (Option Nat × Option Nat) × EmbedState :=
  if t.types[c]! != .I then
    ((s.outerE[c]![1]!, s.outerE[c]![2]!),
      ((link q2 s.outerE[c]![3]!).run ((link q1 s.outerE[c]![0]!).run s).2).2)
  else ((q1, q2), s)

/-- The `Q` step's gluing onto the lower `V` item `w`: the pair `(a, b)` is closed onto `w`'s
exposed pair, or directly if `w` has none. -/
def qLower (w : Nat) (a b : Option Nat) (s : EmbedState) : EmbedState :=
  if s.outerE[w]![0]!.isSome then
    ((link b s.outerE[w]![1]!).run ((link a s.outerE[w]![0]!).run s).2).2
  else ((link a b).run s).2

/-- Expose all four quarter-edges of the leaf `Q` item `i`. -/
def setOuter4 (i : Nat) (s : EmbedState) : EmbedState :=
  ((setOuter i 3 (t.qes i 3)).run ((setOuter i 2 (t.qes i 2)).run
    ((setOuter i 1 (t.qes i 1)).run ((setOuter i 0 (t.qes i 0)).run s).2).2).2).2

theorem embedItem_Q_nil (i : Nat) (ht : t.types[i]! = .Q) (hc : t.children i = [])
    (s : EmbedState) : ((t.embedItem i).run s).2 = t.setOuter4 i s := by
  unfold embedItem
  rw [ht, hc]
  simp only [Std.Legacy.Range.forIn_eq_forIn_range', Std.Legacy.Range.size, Nat.sub_zero, Nat.add_sub_cancel,
    Nat.div_one, List.range', List.forIn_cons, List.forIn_nil, StateT.run_bind, setOuter4]
  rfl

theorem embedItem_Q_O (i c : Nat) (rest : List Nat) (ht : t.types[i]! = .Q)
    (hc : t.children i = c :: rest) (hO : t.types[c]! = .O) (s : EmbedState) :
    ((t.embedItem i).run s).2 =
      setOuterPair i (t.qes i 0) (t.qes i 3) ((link (t.qes i 1) (t.qes i 2)).run s).2 := by
  unfold embedItem
  rw [ht, hc]
  simp only [hO, beq_self_eq_true, ↓reduceIte, StateT.run_bind, setOuterPair]
  rfl

theorem embedItem_Q_cons (i c : Nat) (rest : List Nat) (ht : t.types[i]! = .Q)
    (hc : t.children i = c :: rest) (hO : t.types[c]! ≠ .O) (s : EmbedState) :
    ((t.embedItem i).run s).2 =
      let r := t.qUpper c (t.qes i 1) (t.qes i 2) s
      setOuterPair i (t.qes i 0) r.1.1 (qLower rest.head! (t.qes i 3) r.1.2 r.2) := by
  have idbind {α β : Type} (a : Id α) (f : α → Id β) : (a >>= f) = f a := rfl
  unfold embedItem
  rw [ht, hc]
  have hO' : (t.types[c]! == .O) = false := by simpa using hO
  simp only [hO', Bool.false_eq_true, ↓reduceIte, StateT.run_bind, outer_run, idbind]
  by_cases hI : t.types[c]! = .I
  · simp only [hI, bne_self_eq_false, Bool.false_eq_true, ↓reduceIte, qUpper, StateT.run_bind,
      outer_run, setOuterPair, qLower, link_outerE, idbind]
    rcases hw : s.outerE[rest.head!]![0]! with _ | w0 <;>
      simp only [hw, Option.isSome_none, Option.isSome_some, Bool.false_eq_true, ↓reduceIte,
        idbind, StateT.run_bind, outer_run, link_outerE] <;> rfl
  · have hI' : (t.types[c]! != .I) = true := by simpa using hI
    simp only [hI', ↓reduceIte, qUpper, StateT.run_bind, outer_run, setOuterPair, qLower,
      link_outerE, idbind]
    rcases hw : s.outerE[rest.head!]![0]! with _ | w0 <;>
      simp only [hw, Option.isSome_none, Option.isSome_some, Bool.false_eq_true, ↓reduceIte,
        idbind, StateT.run_bind, outer_run, link_outerE] <;> rfl

end Spqr.PlanarSpqrTree
