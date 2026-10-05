import Spqr.PlanarInv
import Spqr.PlanarWalkProj

/-!
# Invariant P through the planar walk

`Preserves P m` is the Hoare triple `{P} m {P}` on `PlanarWalkM`. The walk is decomposed into the
steps of `planarFinishEdge`: steps that touch neither the tstack nor the planarity stack are
discharged by `Preserves.frame`; the paired pops by `Preserves.popPair`; the steps that create,
merge, close, prune, flip, fold or unwrap planarity data are the named lemmas `*_inv` (admitted:
each is one case of PROOF.md §8.2–8.4). The plumbing — `planarFinishEdge`, `planarWalkOut`,
`planarWalkTree`, `planarWalkOuts` — is proved from them.
-/

namespace Spqr

/-- Invariant P on every entry, with the two stacks aligned. -/
def WalkInv (s : PlanarWalkState) : Prop :=
  StackInv s ∧ s.base.tstack.length = s.aux.plStack.length

def Preserves (P : PlanarWalkState → Prop) (m : PlanarWalkM α) : Prop := ∀ s, P s → P (m s).2

namespace Preserves

variable {P : PlanarWalkState → Prop}

theorem pure (a : α) : Preserves P (Pure.pure a : PlanarWalkM α) := fun _ h => h

theorem bind {m : PlanarWalkM α} {k : α → PlanarWalkM β} (hm : Preserves P m)
    (hk : ∀ a, Preserves P (k a)) : Preserves P (m >>= k) := fun s h => hk _ _ (hm s h)

theorem seqPair {m₁ : PlanarWalkM α} {m₂ : PlanarWalkM β} {k : β → PlanarWalkM γ}
    (h : Preserves P (m₁ >>= fun _ => m₂)) (hk : ∀ b, Preserves P (k b)) :
    Preserves P (m₁ >>= fun _ => m₂ >>= k) := fun s hs => hk _ _ (h s hs)

theorem fmap {m : PlanarWalkM α} (h : Preserves P m) (g : α → β) : Preserves P (g <$> m) := h

theorem ite (c : Prop) [Decidable c] {t e : PlanarWalkM α} (ht : c → Preserves P t)
    (he : ¬c → Preserves P e) : Preserves P (if c then t else e) := by
  split
  · exact ht ‹_›
  · exact he ‹_›

theorem loop (fuel : Nat) {c : PlanarWalkM Bool} {b : PlanarWalkM Unit} (hc : Preserves P c)
    (hb : Preserves P b) : Preserves P (PlanarWalkM.loop fuel c b) := by
  induction fuel with
  | zero => exact pure _
  | succ n ih =>
    unfold PlanarWalkM.loop
    exact bind hc fun _ => ite _ (fun _ => bind hb fun _ => ih) (fun _ => pure _)

theorem loopFirst (fuel : Nat) {c : PlanarWalkM Bool} {b : Bool → PlanarWalkM Unit}
    (hc : Preserves P c) (hb : ∀ x, Preserves P (b x)) (first : Bool) :
    Preserves P (PlanarWalkM.loopFirst fuel c b first) := by
  induction fuel generalizing first with
  | zero => exact pure _
  | succ n ih =>
    unfold PlanarWalkM.loopFirst
    exact bind hc fun _ => ite _ (fun _ => bind (hb _) fun _ => ih _) (fun _ => pure _)

/-- Steps that leave the tstack, the planarity stack, the matches and the graph alone. -/
theorem frame (m : PlanarWalkM α)
    (h : ∀ s, (m s).2.base.tstack = s.base.tstack ∧ (m s).2.aux.plStack = s.aux.plStack ∧
      (m s).2.aux.qem = s.aux.qem ∧ (m s).2.base.g = s.base.g) :
    Preserves WalkInv m := by
  intro s ⟨hinv, hlen⟩
  obtain ⟨h1, h2, h3, h4⟩ := h s
  refine ⟨?_, ?_⟩
  · unfold StackInv at hinv ⊢
    rw [h1, h2, h3, h4]
    exact hinv
  · rw [h1, h2]; exact hlen

theorem mem_zip_tail {l₁ : List α} {l₂ : List β} {x : α × β} (h : x ∈ l₁.tail.zip l₂.tail) :
    x ∈ l₁.zip l₂ := by
  cases l₁ <;> cases l₂ <;> simp_all

/-- A base pop immediately followed by its planarity pop. -/
theorem popPair {k : TEntry → PlEntry → PlanarWalkM α} (hk : ∀ t tp, Preserves WalkInv (k t tp)) :
    Preserves WalkInv (PlanarWalkM.liftW WalkM.popTstack >>= fun t => PlanarWalkM.popPl >>= fun tp => k t tp) := by
  intro s ⟨hinv, hlen⟩
  refine hk _ _ _ ⟨?_, ?_⟩
  · intro x hx p hp hne
    exact hinv x (mem_zip_tail hx) p hp hne
  · show s.base.tstack.tail.length = s.aux.plStack.tail.length
    simp only [List.length_tail, hlen]

end Preserves

namespace PlanarWalkM

theorem default_ends : (default : Planarity).ends = [] := rfl

/-- A fresh vertex entry has an empty piece and no exposed ends, so it is exempt from
Invariant P. -/
theorem pushVertTstack_inv (v d : Nat) : Preserves WalkInv (pushVertTstack v d) := by
  intro s ⟨hinv, hlen⟩
  refine ⟨?_, ?_⟩
  · intro x hx p hp hne
    rcases List.mem_cons.1 hx with rfl | hx
    · cases hp
      exact absurd default_ends hne
    · exact hinv x hx p hp hne
  · show (_ :: s.base.tstack).length = (_ :: s.aux.plStack).length
    simp only [List.length_cons, hlen]

/-- A fresh single-edge ear (`makeEdgePlanarity`) satisfies Invariant P: the piece is one
virtual edge, its embedding is `bondRot 1`, and the four quarter-edges are the exposed ends
(tree edge: two on each side; back edge: all on side 0, the far end as the `top` return).
Admitted. -/
theorem pushEdgeTstack_inv (v d e : Nat) (isTree : Bool) : Preserves WalkInv (pushEdgeTstack v d e isTree) := by
  sorry

/-- `mergeTstackTops` glues the top piece onto the one below along the shared terminal
(`mergeSide`: 2-sum of the two embeddings at the bottom/top ends); the nesting test is the only
case that yields `none`. Admitted (PROOF.md §8.3). -/
theorem mergeTstackTops_inv : Preserves WalkInv mergeTstackTops := by
  sorry

/-- Reopening a finished node (`unwrapPlanarity`) puts its recorded cap matches back as the
exposed ends of the entry below. Admitted. -/
theorem maybeUnwrapNxt_inv (type : NodeType) (isTree : Bool) : Preserves WalkInv (maybeUnwrapNxt type isTree) := by
  sorry

/-- Finishing the top entry records its exposed ends as the node's cap matches
(`finishMatches`) and replaces the piece by the fresh cap ear. Admitted. -/
theorem finishTstackTop_inv (item : ItemId) (isTree : Bool) : Preserves WalkInv (finishTstackTop item isTree) := by
  sorry

/-- `closeSide` links the inner bottom end to the inner return: the back edges of a closed ear
all return to the current depth, so the piece stays embedded with the same outer face. Admitted. -/
theorem closeBackedges_inv : Preserves WalkInv closeBackedges := by
  sorry

/-- `flipEntry` mirrors the embedding (swaps the two sides); Invariant P is symmetric except for
`minimal`, which the flip conditions re-establish. Admitted. -/
theorem flipBeforeMerge_inv (d fo : Nat) (b : Bool) : Preserves WalkInv (flipBeforeMerge d fo b) := by
  sorry

/-- `pruneSide` closes the returns to depth `d` from the inner end of each side and reopens the
next return. Admitted. -/
theorem pruneBackedges_inv (d : Nat) : Preserves WalkInv (pruneBackedges d) := by
  sorry

/-- Mirrors the entries whose side 0 returns to `lowval`. Admitted (as `flipBeforeMerge_inv`). -/
theorem flipForLowval_inv (lowval origTstack : Nat) : Preserves WalkInv (flipForLowval lowval origTstack) := by
  sorry

/-- Leaving a child: the spans are folded onto one side and `foldPlanarity` glues side 1 around
onto side 0 (the two boundary walks become one, split at the new bottom). Admitted. -/
theorem modifyCur_foldSides_inv (f : TEntry → TEntry) (edgeDir : Bool) (lowval : Nat) :
    Preserves WalkInv (liftW (WalkM.modifyCur f) >>= fun _ => foldSides edgeDir lowval) := by
  sorry

end PlanarWalkM

syntax "pres" : tactic
macro_rules
  | `(tactic| pres) => `(tactic| first
    | with_reducible exact Preserves.pure _
    | with_reducible first
      | exact PlanarWalkM.pushVertTstack_inv _ _ | exact PlanarWalkM.pushEdgeTstack_inv _ _ _ _
      | exact PlanarWalkM.mergeTstackTops_inv | exact PlanarWalkM.maybeUnwrapNxt_inv _ _
      | exact PlanarWalkM.finishTstackTop_inv _ _ | exact PlanarWalkM.closeBackedges_inv
      | exact PlanarWalkM.flipBeforeMerge_inv _ _ _ | exact PlanarWalkM.pruneBackedges_inv _
      | exact PlanarWalkM.flipForLowval_inv _ _
    | with_reducible solve_by_elim only [*]
    | (with_reducible refine Preserves.popPair fun _ _ => ?_; pres)
    | (with_reducible refine Preserves.seqPair (PlanarWalkM.modifyCur_foldSides_inv _ _ _) fun _ => ?_; pres)
    | (with_reducible refine Preserves.bind ?_ fun _ => ?_ <;> pres)
    | (with_reducible refine Preserves.ite _ (fun _ => ?_) (fun _ => ?_) <;> pres)
    | (with_reducible refine Preserves.fmap ?_ _; pres)
    | (with_reducible refine Preserves.loop _ ?_ ?_ <;> pres)
    | (with_reducible refine Preserves.loopFirst _ ?_ (fun _ => ?_) _ <;> pres)
    | exact Preserves.frame _ fun _ => ⟨rfl, rfl, rfl, rfl⟩
    | (dsimp only; pres)
    | (split <;> pres)
    | (unfold letFun; pres)
    | skip)

namespace PlanarWalkM

set_option maxRecDepth 4000 in
theorem planarFinishEdge_inv (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool) :
    Preserves WalkInv (planarFinishEdge curV d o origTstack hasVert) := by
  unfold planarFinishEdge; pres

end PlanarWalkM

mutual
theorem planarWalkTree_inv (t : DfsTree) (d : Nat) : Preserves WalkInv (planarWalkTree t d) := by
  unfold planarWalkTree
  match t with
  | .node v outs =>
    with_reducible refine Preserves.bind (Preserves.frame _ fun _ => ⟨rfl, rfl, rfl, rfl⟩) fun _ => ?_
    with_reducible refine Preserves.bind (planarWalkOuts_inv v d outs false) fun _ => ?_
    pres

theorem planarWalkOuts_inv (v d : Nat) (outs : List DfsOut) (hasVert : Bool) :
    Preserves WalkInv (planarWalkOuts v d outs hasVert) := by
  unfold planarWalkOuts
  match outs with
  | [] => pres
  | o :: rest =>
    with_reducible refine Preserves.bind (planarWalkOut_inv v d o hasVert) fun _ => ?_
    exact planarWalkOuts_inv v d rest _

theorem planarWalkOut_inv (v d : Nat) (o : DfsOut) (hasVert : Bool) :
    Preserves WalkInv (planarWalkOut v d o hasVert) := by
  unfold planarWalkOut
  match o with
  | .tree a b child =>
    pres
    all_goals first
      | exact planarWalkTree_inv _ _
      | exact PlanarWalkM.planarFinishEdge_inv _ _ _ _ _
  | .back a b c =>
    pres
    all_goals exact PlanarWalkM.planarFinishEdge_inv _ _ _ _ _
end

/-- The planar walk preserves Invariant P on every entry (modulo the per-step lemmas `*_inv`). -/
theorem planarWalkOut_stackInv (v d : Nat) (o : DfsOut) (hasVert : Bool) (s : PlanarWalkState)
    (h : StackInv s) (hlen : s.base.tstack.length = s.aux.plStack.length) :
    StackInv ((planarWalkOut v d o hasVert).run s).2 :=
  (planarWalkOut_inv v d o hasVert s ⟨h, hlen⟩).1

end Spqr
