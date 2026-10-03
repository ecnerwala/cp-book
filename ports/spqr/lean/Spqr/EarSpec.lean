import Spqr.Ear
import Spqr.Build
import Spqr.Spec
import Spqr.ItemSpec
import Spqr.RelabelRep
import Spqr.WalkWF
import Spqr.RelabelWF
import Spqr.Frame
import Spqr.Proofs.Dfs
import Mathlib.Data.List.TakeDrop

/-!
# Ear-level statements for phase 2

The walk proof (`PROOF.md`, §4) is stated on the ear-structured walk `walkEarTree`, which is
transported to `walkTree` by `walkEarTree_eq_walkTree`.
-/

namespace Spqr

theorem walkOut_eq_walkOutWith (v d : Nat) (o : DfsOut) (hv : Bool) :
    walkOut v d o hv = walkOutWith walkTree v d o hv := by
  rw [walkOut]; rfl

theorem earOut_succ (fuel v d : Nat) (o : DfsOut) (hv : Bool) :
    earOut (fuel + 1) v d o hv = walkOutWith (walkEarTree fuel) v d o hv := by
  rw [earOut]; rfl

theorem walkOuts_eq_walkOutsWith (v d : Nat) (outs : List DfsOut) (hv : Bool) :
    walkOuts v d outs hv = walkOutsWith walkTree v d outs hv := by
  induction outs generalizing hv with
  | nil => rw [walkOuts]; rfl
  | cons o rest ih => rw [walkOuts]; simp only [walkOutsWith, walkOut_eq_walkOutWith, ih]

theorem earOuts_succ (fuel v d : Nat) (outs : List DfsOut) (hv : Bool) :
    earOuts (fuel + 1) v d outs hv = walkOutsWith (walkEarTree fuel) v d outs hv := by
  induction outs generalizing hv with
  | nil => rw [earOuts]; rfl
  | cons o rest ih => rw [earOuts]; simp only [walkOutsWith, earOut_succ, ih]

theorem walkOutWith_congr {w w' : DfsTree → Nat → WalkM Unit} (v d : Nat) (o : DfsOut) (hv : Bool)
    (h : ∀ e cls c, o = .tree e cls c → w c (d + 1) = w' c (d + 1)) :
    walkOutWith w v d o hv = walkOutWith w' v d o hv := by
  cases o with
  | back => rfl
  | tree e cls c => simp only [walkOutWith, h e cls c rfl]

theorem walkOutsWith_congr {w w' : DfsTree → Nat → WalkM Unit} (v d : Nat) (outs : List DfsOut)
    (hv : Bool) (h : ∀ o ∈ outs, ∀ e cls c, o = .tree e cls c → w c (d + 1) = w' c (d + 1)) :
    walkOutsWith w v d outs hv = walkOutsWith w' v d outs hv := by
  induction outs generalizing hv with
  | nil => rfl
  | cons o rest ih =>
    simp only [walkOutsWith, walkOutWith_congr v d o hv (h o (List.mem_cons_self ..))]
    congr 1; funext hv'
    exact ih hv' fun o ho => h o (List.mem_cons_of_mem _ ho)

theorem walkOutsWith_append (w : DfsTree → Nat → WalkM Unit) (v d : Nat) (l₁ l₂ : List DfsOut)
    (hv : Bool) :
    walkOutsWith w v d (l₁ ++ l₂) hv = walkOutsWith w v d l₁ hv >>= walkOutsWith w v d l₂ := by
  induction l₁ generalizing hv with
  | nil => simp [walkOutsWith]
  | cons o rest ih => simp [walkOutsWith, ih]

theorem ascend_nil (fuel : Nat) : ascend fuel [] = pure () := by
  cases fuel with
  | zero => exact ascend.eq_1 _
  | succ n => exact ascend.eq_2 _ (by omega)

theorem DfsOut.height_le_heightList {outs : List DfsOut} {e : Nat} {cls : OutClass} {c : DfsTree}
    (h : DfsOut.tree e cls c ∈ outs) : c.height ≤ DfsOut.heightList outs := by
  induction outs with
  | nil => simp at h
  | cons o rest ih =>
    rcases List.mem_cons.mp h with rfl | h'
    · exact Nat.le_max_left _ _
    · cases o
      · exact ih h'
      · exact Nat.le_trans (ih h') (Nat.le_max_right _ _)

theorem descend_ascend (n : Nat) : ∀ (t : DfsTree) (F G d : Nat) (acc : List Frame),
    t.height ≤ n → t.height ≤ F → t.height ≤ G →
    (descend F t d acc >>= fun fs => ascend G fs) = (walkTree t d >>= fun _ => ascend G acc) := by
  induction n with
  | zero => intro ⟨v, outs⟩ _ _ _ _ h; simp [DfsTree.height] at h
  | succ n ih =>
    intro ⟨v, outs⟩ F G d acc hn hF hG
    have hw : ∀ (c : DfsTree) (fuel d : Nat), c.height ≤ n → c.height ≤ fuel →
        walkEarTree fuel c d = walkTree c d := by
      intro c fuel d hc hfuel
      rw [walkEarTree.eq_1, ih c fuel fuel d [] hc hfuel hfuel]
      simp [ascend_nil]
    have hch : ∀ o ∈ outs, ∀ e cls c, o = .tree e cls c → c.height ≤ n := by
      intro o ho e cls c rfl
      have := DfsOut.height_le_heightList ho
      simp [DfsTree.height] at hn; omega
    obtain ⟨F, rfl⟩ : ∃ F', F = F' + 1 := ⟨F - 1, by simp [DfsTree.height] at hF; omega⟩
    obtain ⟨G, rfl⟩ : ∃ G', G = G' + 1 := ⟨G - 1, by simp [DfsTree.height] at hG; omega⟩
    have hF' : ∀ o ∈ outs, ∀ e cls c, o = .tree e cls c → walkEarTree F c (d + 1) = walkTree c (d + 1) := by
      intro o ho e cls c h
      have := DfsOut.height_le_heightList (h ▸ ho)
      simp [DfsTree.height] at hF
      exact hw c F (d + 1) (hch o ho e cls c h) (by omega)
    have hG' : ∀ o ∈ outs, ∀ e cls c, o = .tree e cls c → walkEarTree G c (d + 1) = walkTree c (d + 1) := by
      intro o ho e cls c h
      have := DfsOut.height_le_heightList (h ▸ ho)
      simp [DfsTree.height] at hG
      exact hw c G (d + 1) (hch o ho e cls c h) (by omega)
    obtain ⟨boundary, rets, hsp⟩ : ∃ b r, outs.span (fun o => decide (o.cls.lowval d ≥ d)) = (b, r) :=
      ⟨_, _, rfl⟩
    have hbr : boundary ++ rets = outs := by
      rw [List.span_eq_takeWhile_dropWhile] at hsp
      cases hsp; exact List.takeWhile_append_dropWhile
    subst hbr
    rw [descend.eq_2, hsp, walkTree.eq_1]
    simp only [earOuts_succ, walkOuts_eq_walkOutsWith, walkOutsWith_append,
      walkOutsWith_congr v d boundary false (fun o ho => hF' o (List.mem_append_left _ ho)),
      finishVert, bind_assoc]
    congr 1; funext _; congr 1; funext hv
    have hrets := walkOutsWith_congr v d rets hv (fun o ho => hF' o (List.mem_append_right _ ho))
    cases rets with
    | nil => simp only [hrets, bind_assoc, pure_bind]
    | cons o rest =>
      cases o with
      | back e dest cls => simp only [hrets, bind_assoc, pure_bind]
      | tree e cls child =>
        by_cases hc : (decide (OutClass.lowval d cls < d) && !cls.isType1) = true
        · simp only [hc, ite_true, bind_assoc]
          have hmem : DfsOut.tree e cls child ∈ boundary ++ DfsOut.tree e cls child :: rest :=
            List.mem_append_right _ (List.mem_cons_self ..)
          have hh := DfsOut.height_le_heightList hmem
          simp only [DfsTree.height] at hF hG
          have hih : ∀ acc', (descend F child (d + 1) acc' >>= fun fs => ascend (G + 1) fs) =
              (walkTree child (d + 1) >>= fun _ => ascend (G + 1) acc') :=
            fun acc' => ih child F (G + 1) (d + 1) acc' (hch _ hmem e cls child rfl) (by omega) (by omega)
          simp only [hih, ascend.eq_3]
          simp only [earOuts_succ, walkOutsWith, walkOutWith, DfsOut.cls, finishVert, bind_assoc,
            walkOutsWith_congr v d rest _ (fun o ho => hG' o (List.mem_append_right _ (List.mem_cons_of_mem _ ho)))]
          obtain ⟨hlt, ht1⟩ : OutClass.lowval d cls < d ∧ cls.isType1 = false := by
            simpa using hc
          have hnle : ¬ (OutClass.lowval d cls ≥ d) := by omega
          simp only [hnle, ite_false, ht1, Bool.and_false, Bool.false_eq_true, bind_assoc, pure_bind]
        · simp only [hc, Bool.false_eq_true, ite_false, hrets, bind_assoc, pure_bind]

/-- `walkEarTree` is the walk, given enough fuel. -/
theorem walkEarTree_eq_walkTree (t : DfsTree) (d fuel : Nat) (h : t.height ≤ fuel) :
    walkEarTree fuel t d = walkTree t d := by
  rw [walkEarTree.eq_1, descend_ascend _ t fuel fuel d [] (Nat.le_refl _) h h]
  simp [ascend_nil]

theorem foldl_inv {β γ : Type _} (P : β → Prop) (f : β → γ → β) (l : List γ) (init : β)
    (h0 : P init) (hf : ∀ b x, P b → P (f b x)) : P (l.foldl f init) := by
  induction l generalizing init with
  | nil => exact h0
  | cons x l ih => exact ih _ (hf _ _ h0)

theorem DfsOut.heightList_le {outs : List DfsOut} {n : Nat}
    (h : ∀ e cls c, DfsOut.tree e cls c ∈ outs → c.height ≤ n) : DfsOut.heightList outs ≤ n := by
  induction outs with
  | nil => simp [DfsOut.heightList]
  | cons o rest ih =>
    have ih := ih fun e cls c hc => h e cls c (List.mem_cons_of_mem _ hc)
    cases o with
    | back => simpa [DfsOut.heightList] using ih
    | tree e cls c => exact Nat.max_le.mpr ⟨h e cls c (List.mem_cons_self ..), ih⟩

theorem dfsVisit_height (adj : Array (List (Nat × Nat))) :
    ∀ (fuel v d : Nat) (prvE : Option Nat) (depth : Array (Option Nat)),
      (dfsVisit adj fuel v d prvE depth).1.height ≤ fuel + 1 := by
  intro fuel
  induction fuel with
  | zero => intros; simp [dfsVisit, DfsTree.height, DfsOut.heightList]
  | succ fuel ih =>
    intro v d prvE depth
    rw [dfsVisit]
    simp only [DfsTree.height]
    refine Nat.add_le_add_right (DfsOut.heightList_le fun e cls c hc => ?_) 1
    rw [List.mem_mergeSort, List.mem_reverse] at hc
    revert hc
    refine foldl_inv (fun p : List DfsOut × Lowvals × Array (Option Nat) =>
        ∀ e cls c, DfsOut.tree e cls c ∈ p.1 → c.height ≤ fuel + 1) _ _ _ (by simp) ?_ e cls c
    rintro ⟨outs, lv, dp⟩ ⟨nxt, e'⟩ hp e cls c hc
    simp only at hc
    split at hc
    · exact hp e cls c hc
    · split at hc
      · rcases List.mem_cons.mp hc with h | h
        · cases h; exact ih ..
        · exact hp e cls c h
      · rcases List.mem_cons.mp hc with h | h
        · cases h
        · exact hp e cls c h

mutual
theorem DfsTree.height_le_verts_length : ∀ t : DfsTree, t.height ≤ t.verts.length
  | .node _ outs => by
    have := DfsOut.heightList_le_vertsList_length outs
    simp only [DfsTree.height, DfsTree.verts, List.length_cons]; omega
theorem DfsOut.heightList_le_vertsList_length :
    ∀ outs : List DfsOut, DfsOut.heightList outs ≤ (DfsOut.vertsList outs).length
  | [] => by simp [DfsOut.heightList, DfsOut.vertsList]
  | .back .. :: rest => by
    simpa [DfsOut.heightList, DfsOut.vertsList] using DfsOut.heightList_le_vertsList_length rest
  | .tree _ _ c :: rest => by
    have h1 := DfsTree.height_le_verts_length c
    have h2 := DfsOut.heightList_le_vertsList_length rest
    simp only [DfsOut.heightList, DfsOut.vertsList, List.length_append]; omega
end

theorem dfsForest_height (g : Graph) (vo eo : List Nat) :
    ∀ t ∈ g.dfsForest vo eo, t.height ≤ g.nv + 1 := by
  intro t ht
  unfold Graph.dfsForest at ht
  simp only [List.mem_reverse] at ht
  revert ht
  refine foldl_inv (fun p : List DfsTree × Array (Option Nat) => ∀ t ∈ p.1, t.height ≤ g.nv + 1)
    _ _ _ (by simp) ?_ t
  rintro ⟨roots, dp⟩ rt hp t ht
  simp only at ht
  split at ht
  · exact hp t ht
  · rcases List.mem_cons.mp ht with rfl | h
    · exact dfsVisit_height ..
    · exact hp t h

theorem walkEarForest_eq_walkForest (fuel : Nat) (forest : List DfsTree)
    (h : ∀ t ∈ forest, t.height ≤ fuel) : walkEarForest fuel forest = walkForest forest := by
  unfold walkEarForest walkForest
  induction forest with
  | nil => rfl
  | cons t rest ih =>
    rw [List.forM, List.forM, walkEarTree_eq_walkTree t 0 fuel (h t (List.mem_cons_self ..)),
      ih fun t ht => h t (List.mem_cons_of_mem _ ht)]

theorem walkEar_eq_walk (g : Graph) (tern : Bool) (vo eo : List Nat) :
    g.walkEar tern (g.dfsForest vo eo) = g.walk tern (g.dfsForest vo eo) := by
  unfold Graph.walkEar
  rw [walkEarForest_eq_walkForest _ _ (dfsForest_height g vo eo)]; rfl

/-- `m` never looks at the tstack it was started on: it computes what it would on an empty tstack,
leaving the original tstack underneath. -/
def TstackLocal (m : WalkM α) (s : WalkState) : Prop :=
  let (a, s') := m.run s
  let (a', s'') := m.run { s with tstack := [] }
  a = a' ∧ s' = { s'' with tstack := s''.tstack ++ s.tstack }

/-- Frame rule: the walk of a subtree is local to the entries it creates. Entries already on the
tstack that return above the subtree's depth and were not started at its vertices are neither
inspected nor modified. The guards of the walk from the empty tstack are an input (`GuardsTree`;
at a root they are `walkTree_guards`, `EarShape.lean`). -/
theorem walkTree_local (t : DfsTree) (anc : List Nat) (s : WalkState)
    (hg : GuardsTree t anc.length { s with tstack := [] })
    (hS : ∀ e ∈ s.tstack, e.topDepth < anc.length ∧ e.vStart ∉ t.verts) :
    TstackLocal (walkTree t anc.length) s := by
  obtain ⟨a, h, -⟩ := Sim.walkTree (bot := s.tstack) t anc.length hS { s with tstack := [] } hg
  unfold TstackLocal
  rw [show lift s.tstack { s with tstack := [] } = s from rfl] at h
  rw [h]
  rcases (walkTree t anc.length).run { s with tstack := [] } with ⟨a', s''⟩
  exact ⟨Subsingleton.elim _ _, rfl⟩


/-! PROOF.md Lemma 4.4 ("a finished frame's vertex owns exactly one entry", formerly
`ascend_frame_one_entry`) is false for chain frames: a frame `(v, d)` is a type-2 chain edge, so the
ear continues above `d` and its entries returning above `d` stay open. Cycle `0-1-2-3-4-5-0` plus
chord `5-1`: at the frame `(4, 4)` (`origTstack = 0`) the tstack after `finishEdge`, the remaining
out-edges and `finishVert` is `[(4,4), (5,4), (5,1), (5,0), (5,5)]` — five entries. The one-entry
collapse holds only at the ear's top, i.e. at the type-1 edge closing the chain, which is a
non-first out-edge handled by `earOut`: that is `earOut_one_entry`. -/

end Spqr
