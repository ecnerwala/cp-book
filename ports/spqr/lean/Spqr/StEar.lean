import Spqr.StRefEt
import Spqr.Ear

/-!
# Lowvals along an ear

In a well-formed DFS tree (`DfsTree.WF`, established by `dfsVisit_spec`), the lowval of an
out-edge is the minimum of the depths its piece returns to (`DfsOut.lowval_eq_lmin`), and the
first returning out-edge of a child classified `ret l _` returns to exactly `l`
(`first_ret_lowval`). Hence along the first-child chain of an ear every vertex sets `stackDir` to
the same value `!stackDir[l]` (`chain_stackDir_step`): the "semantic half" of the ear uniform-side
fact of `PROOF.md` §7.
-/

namespace Spqr

/-- Two returning out-edges in rank order have non-decreasing lowvals. -/
theorem lowval_le_of_rank_le {d : Nat} {o o' : DfsOut} {l k l' k'} (ho : o.cls = .ret l k)
    (ho' : o'.cls = .ret l' k') (hr : o.cls.rank ≤ o'.cls.rank) :
    o.cls.lowval d ≤ o'.cls.lowval d := by
  rw [ho, ho'] at hr ⊢
  simp only [OutClass.rank, OutClass.lowval] at hr ⊢
  have : k.rank ≤ 2 := by cases k <;> decide
  have : k'.rank ≤ 2 := by cases k' <;> decide
  omega

/-- The first returning out-edge of a child classified `ret l _` returns to exactly `l`
(`classify_eq_ret_iff_tree` gives `l = lmin` of the child's return depths; sortedness by
`OutClass.rank` from `DfsTree.WF` puts the minimal one first). -/
theorem first_ret_lowval (anc : List Nat) (v w l : Nat) (k : RetKind) (outs : List DfsOut)
    (hwf : (DfsTree.node w outs).WF (anc ++ [v]))
    (hcls : classify anc.length true
      (low2 (anc.length + 1) ((DfsTree.node w outs).retDepths (anc.length + 1))) = .ret l k)
    (pre : List DfsOut) (o : DfsOut) (rest : List DfsOut) (houts : outs = pre ++ o :: rest)
    (hpre : ∀ o' ∈ pre, anc.length + 1 ≤ o'.cls.lowval (anc.length + 1))
    (ho : o.cls.lowval (anc.length + 1) < anc.length + 1) :
    o.cls.lowval (anc.length + 1) = l := by
  set d := anc.length with hd
  have hlen : (anc ++ [v]).length = d + 1 := by simp [hd]
  rw [DfsTree.WF] at hwf
  obtain ⟨hsorted, hwfo⟩ := hwf
  obtain ⟨hl, hld, -⟩ := classify_eq_ret_iff_tree.1 hcls
  rw [low2_fst] at hl
  simp only [DfsTree.retDepths, DfsOut.retDepthsList_eq] at hl
  have hwfo' : ∀ o' ∈ outs, o'.cls.lowval (d + 1) = lmin (d + 2) (o'.retDepths (d + 1)) := by
    intro o' ho'; have := DfsOut.lowval_eq_lmin (hwfo o' ho'); rwa [hlen] at this
  have homem : o ∈ outs := by rw [houts]; simp
  -- `l ≤ lowval o`
  have h1 : l ≤ o.cls.lowval (d + 1) := by
    rw [← hl]
    apply lmin_le_of_mem
    rw [List.mem_flatMap]
    refine ⟨o, homem, ?_⟩
    have := DfsOut.lowval_mem_retDepths (hwfo o homem) (by rwa [hlen])
    rwa [hlen] at this
  -- `lowval o ≤ l`
  have h2 : o.cls.lowval (d + 1) ≤ l := by
    have hlmem : l ∈ outs.flatMap (DfsOut.retDepths (d + 1)) := by
      rcases lmin_eq_or_mem (d + 1) (outs.flatMap (DfsOut.retDepths (d + 1))) with h | h
      · omega
      · rwa [hl] at h
    obtain ⟨o'', ho'', hlo''⟩ := List.mem_flatMap.1 hlmem
    have hle : o''.cls.lowval (d + 1) ≤ l := by
      have := DfsOut.lowval_le_retDepths (hwfo o'' ho'') l (by rwa [hlen]); rwa [hlen] at this
    have hret'' : o''.cls.lowval (d + 1) < d + 1 := by omega
    rw [houts, List.mem_append, List.mem_cons] at ho''
    rcases ho'' with hp | rfl | hr
    · have := hpre o'' hp; omega
    · exact hle
    · obtain ⟨lo, ko, hco⟩ := DfsOut.ret_of_lowval_lt (hwfo o homem) (by rwa [hlen])
      obtain ⟨lo', ko', hco'⟩ :=
        DfsOut.ret_of_lowval_lt (hwfo o'' (by rw [houts]; simp [hr])) (by rwa [hlen])
      have hrank : o.cls.rank ≤ o''.cls.rank := by
        rw [houts] at hsorted
        exact (List.pairwise_cons.1 (List.pairwise_append.1 hsorted).2.1).1 o'' hr
      exact Nat.le_trans (lowval_le_of_rank_le hco hco' hrank) hle
  omega

/-- `walkOut`/`earOut` at depth `d + 1` on an out-edge returning to `l` sets `stackDir[d+1]` to
`!stackDir[l]`, which equals the parent's `stackDir[d]` when that was set the same way. -/
theorem chain_stackDir_step (s : WalkState) (d l : Nat) (hl : l < d)
    (hsz : d + 1 < s.stackDir.size) (hd : s.stackDir[d]! = !s.stackDir[l]!) (b : Bool)
    (hb : b = if l ≥ d + 1 then false else !s.stackDir[l]!) :
    (s.stackDir.set! (d + 1) b)[d + 1]! = (s.stackDir.set! (d + 1) b)[d]! := by
  have hne : d + 1 ≠ d := by omega
  have : ¬ l ≥ d + 1 := by omega
  simp only [this, ↓reduceIte] at hb
  rw [Array.getElem!_set!_self _ _ _ hsz, Array.getElem!_set!_ne _ _ _ _ hne, hb, hd]

end Spqr
