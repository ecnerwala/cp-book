import Spqr.StBoundary
import Spqr.StSimLemmas
import Spqr.StFrame
import Spqr.EarSides
import Spqr.EarWalk

/-!
# The st-simulation by the walk induction

`StPre` is the simulation state at an out-edge list (the open stack segment reads the pieces of the
out-edges already walked, the base stack is the segments of the enclosing frames), `StHyps` the
hypotheses threaded through `walkTree.mutual_induct`. `stOut_step` is the one-edge step: a boundary
edge by `finishBoundary_st`, a returning edge by `finishEdge_st` (+ `finishRet_frame_st`).
-/

namespace Spqr

open WalkM WalkState StRefEt

/-- `m` keeps `stackDir.size`. -/
def KeepsSize (m : WalkM α) : Prop :=
  ∀ s, wp m (fun _ s' => s'.stackDir.size = s.stackDir.size) s

namespace KeepsSize

theorem of_all {m : WalkM α} (h : KeepsAll m) : KeepsSize m := fun s => congrArg Array.size (h s)

theorem bind {m : WalkM β} {f : β → WalkM α} (hm : KeepsSize m) (hf : ∀ b, KeepsSize (f b)) :
    KeepsSize (m >>= f) := by
  intro s
  show (((f (m.run s).1).run (m.run s).2)).2.stackDir.size = s.stackDir.size
  rw [hf _ _, hm s]

theorem ite (c : Prop) [Decidable c] {a b : WalkM α} (ha : c → KeepsSize a) (hb : ¬ c → KeepsSize b) :
    KeepsSize (if c then a else b) := by
  split
  · exact ha ‹_›
  · exact hb ‹_›

theorem setStackDir (i : Nat) (b : Bool) : KeepsSize (setStackDir i b) := fun _ => Array.size_set! _ _ _

theorem loop (n : Nat) {cond : WalkM Bool} {body : WalkM Unit} (hc : KeepsSize cond)
    (hb : KeepsSize body) : KeepsSize (loop n cond body) := by
  induction n with
  | zero => exact of_all (KeepsAll.pure ())
  | succ n ih =>
    simp only [WalkM.loop]
    exact bind hc fun b => ite _ (fun _ => bind hb fun _ => ih) fun _ => of_all (KeepsAll.pure ())

theorem loop1Type (d : Nat) (edgeDir : Bool) : KeepsSize (loop1Type d edgeDir) := by
  intro s
  unfold Spqr.loop1Type; wp_simp
  repeat' split
  all_goals first | rfl | simp

theorem loop1Body (d : Nat) (edgeDir : Bool) : KeepsSize (loop1Body d edgeDir) := by
  unfold Spqr.loop1Body
  refine bind (loop1Type d edgeDir) fun ty => bind (of_all (KeepsAll.maybeUnwrapNxt ty)) fun i => ?_
  exact bind (of_all KeepsAll.mergeTstackTops) fun _ => of_all (KeepsAll.finishTstackTop i)

theorem closeEars (nxtV d e : Nat) (edgeDir : Bool) : KeepsSize (closeEars nxtV d e edgeDir) := by
  unfold Spqr.closeEars
  refine bind (of_all (KeepsAll.pushEdgeTstack _ _ _)) fun _ =>
    bind (of_all KeepsAll.tstackSize) fun n => ?_
  exact loop n (of_all (KeepsAll.loop1Cond d)) (loop1Body d edgeDir)

theorem finishTree (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert edgeDir : Bool) :
    KeepsSize (finishTree curV d o origTstack hasVert edgeDir) := by
  unfold Spqr.finishTree
  refine bind (closeEars _ _ _ _) fun _ => bind (of_all (KeepsAll.mergeLate d)) fun isSingle => ?_
  refine ite _ (fun _ => ?_) fun _ => of_all (KeepsAll.finishRest _ _ _ _ _ _)
  exact of_all (KeepsAll.closeVert _ _ _ _ _ fun b => KeepsAll.finishRest _ _ _ _ _ _)

theorem finishEdge (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool) :
    KeepsSize (finishEdge curV d o origTstack hasVert) := by
  rw [finishEdge_eq]
  unfold finishEdge'
  refine bind (of_all KeepsAll.get) fun s => bind (of_all (KeepsAll.stackDir d)) fun edgeDir => ?_
  refine ite _ (fun _ => of_all (KeepsAll.finishBoundary _ _ _ _ _)) fun _ => ?_
  refine bind (of_all (KeepsAll.makeVs _ _)) fun vs => bind (of_all (KeepsAll.modifyItem _ _)) fun _ => ?_
  exact ite _ (fun _ => finishTree _ _ _ _ _ _) fun _ => of_all (KeepsAll.finishBack _ _ _ _)

end KeepsSize

/-- Simulation state after the out-edges `done` of `v` at depth `d` (`d = fs.length`), base stack
`segsStack segs`. -/
structure StPre (g : Graph) (prev : List DfsTree) (fs : List PathFrame)
    (segs : List (List TEntry × List StPiece)) (v d : Nat) (done : List DfsOut) (hasVert : Bool)
    (s : WalkState) : Prop where
  read : ∃ new, s.tstack = new ++ segsStack segs ∧
    StRead s.items new (refOuts g v d (DirsOf s d) done false).1 ∧ (hasVert = false → new = [])
  hv : (refOuts g v d (DirsOf s d) done false).2.2 = hasVert
  segRead : SegRead s.items segs
  vStart : ∀ t ∈ s.tstack,
    t.vStart ∈ v :: DfsOut.vertsList done ∨ ∃ t₀ ∈ segsStack segs, t.vStart = t₀.vStart
  items : StItems g s
    (simBlocks g prev fs (DirsOf s d) ++ (refOuts g v d (DirsOf s d) done false).2.1)

/-- The hypotheses of the out-edge steps (`Full` for the freshness of the new `V`/`Q` items, the
DFS bounds and distinctness, the array size for the subtree, the base entries starting outside). -/
structure StHyps (g : Graph) (prev : List DfsTree) (fs : List PathFrame)
    (segs : List (List TEntry × List StPiece)) (v d : Nat) (done outs : List DfsOut) (hasVert : Bool)
    (P X : ItemId → Prop) (s : WalkState) : Prop where
  full : s.Full g P X
  v_lt : v < g.nv
  w_lt : ∀ w ∈ DfsOut.vertsList outs, w < g.nv
  e_lt : ∀ e ∈ DfsOut.edgesList outs, e < g.ne
  vnodup : (v :: DfsOut.vertsList (done ++ outs)).Nodup
  enodup : (DfsOut.edgesList outs).Nodup
  Pv : ∀ w ∈ DfsOut.vertsList outs, ¬ P (vertItem w)
  Pe : ∀ e ∈ DfsOut.edgesList outs, ¬ P (edgeItem g e)
  Pcur : hasVert = false → ¬ P (vertItem v)
  Pcur' : hasVert = true → P (vertItem v)
  height : d + DfsOut.heightList outs < s.stackDir.size
  base : ∀ t ∈ segsStack segs, t.vStart ∉ v :: DfsOut.vertsList (done ++ outs)
  pre : StPre g prev fs segs v d done hasVert s

abbrev StTreeP (g : Graph) (t : DfsTree) (d : Nat) : Prop :=
  ∀ (s : WalkState) (prev : List DfsTree) (fs : List PathFrame)
    (segs : List (List TEntry × List StPiece)) (P X : ItemId → Prop),
    d = fs.length →
    (∀ v outs, t = .node v outs → ({ s with stackVerts := s.stackVerts.set! d v } : WalkState).Inv' d) →
    Shape s → GuardsTree t d s → BookTree t d s →
    s.Full g P X → (∀ v ∈ t.verts, v < g.nv) → (∀ e ∈ t.edges, e < g.ne) →
    t.verts.Nodup → t.edges.Nodup →
    (∀ v ∈ t.verts, ¬ P (vertItem v)) → (∀ e ∈ t.edges, ¬ P (edgeItem g e)) →
    d + t.height ≤ s.stackDir.size →
    s.tstack = segsStack segs → SegRead s.items segs →
    (∀ t' ∈ segsStack segs, t'.vStart ∉ t.verts) →
    StItems g s (simBlocks g prev fs (DirsOf s d)) →
    wp (walkTree t d) (fun _ s' => s'.g = s.g ∧ s'.stackDir.size = s.stackDir.size ∧
      (∀ t' ∈ s'.tstack, t'.vStart ∈ t.verts ∨ ∃ t₀ ∈ segsStack segs, t'.vStart = t₀.vStart) ∧
      SegRead s'.items segs ∧ StSim g prev fs t (segsStack segs) s') s

abbrev StOutsP (g : Graph) (v d : Nat) (outs : List DfsOut) (hasVert : Bool) : Prop :=
  ∀ (s : WalkState) (prev : List DfsTree) (fs : List PathFrame)
    (segs : List (List TEntry × List StPiece)) (done : List DfsOut) (P X : ItemId → Prop),
    d = fs.length → s.Inv' d → Shape s → GuardsOuts v d outs hasVert s → BookOuts v d outs hasVert s →
    StHyps g prev fs segs v d done outs hasVert P X s →
    wp (walkOuts v d outs hasVert) (fun hv' s' => s'.g = s.g ∧ s'.stackDir.size = s.stackDir.size ∧
      StPre g prev fs segs v d (done ++ outs) hv' s') s

abbrev StOutP (g : Graph) (v d : Nat) (o : DfsOut) (hasVert : Bool) : Prop :=
  ∀ (s : WalkState) (prev : List DfsTree) (fs : List PathFrame)
    (segs : List (List TEntry × List StPiece)) (done : List DfsOut) (P X : ItemId → Prop),
    d = fs.length → s.Inv' d → Shape s → GuardsOut v d o hasVert s → BookOut v d o hasVert s →
    StHyps g prev fs segs v d done [o] hasVert P X s →
    wp (walkOut v d o hasVert) (fun hv' s' => s'.g = s.g ∧ s'.stackDir.size = s.stackDir.size ∧
      StPre g prev fs segs v d (done ++ [o]) hv' s') s

/-! ## `walkOutPre` -/

theorem walkOutPre_run (s : WalkState) (v d : Nat) (o : DfsOut) (hasVert : Bool) :
    (walkOutPre v d o hasVert).run s =
      if (!hasVert && decide (o.cls.lowval d < d) && o.cls.isType1) = true then
        (true, { s with
          stackDir := s.stackDir.set! d (if d ≤ o.cls.lowval d then false else !s.stackDir[o.cls.lowval d]!),
          tstack := ⟨v, d, s.nxtEdgeIdx,
            setSides (s.stackDir.set! d
              (if d ≤ o.cls.lowval d then false else !s.stackDir[o.cls.lowval d]!))[d]!
              [vertItem v] []⟩ :: s.tstack })
      else (hasVert, { s with
        stackDir := s.stackDir.set! d
          (if d ≤ o.cls.lowval d then false else !s.stackDir[o.cls.lowval d]!) }) := by
  by_cases hp : (!hasVert && decide (o.cls.lowval d < d) && o.cls.isType1) = true
  · rw [ite_eq_left hp]; unfold walkOutPre; simp only [hp, ↓reduceIte]; rfl
  · rw [ite_eq_right hp]; unfold walkOutPre; simp only [hp]; rfl

/-- What `walkOutPre` gives (`walkOutPre_st` plus the stack frame). -/
structure PreSpec (g : Graph) (s : WalkState) (v d : Nat) (o : DfsOut) (hasVert : Bool)
    (base : List TEntry) (ps : List StPiece) (blocks : List StBlock) (hv₁ : Bool) (s₁ : WalkState) :
    Prop where
  hv : hv₁ = ((!hasVert && decide (o.cls.lowval d < d) && o.cls.isType1) || hasVert)
  items : s₁.items = s.items
  geq : s₁.g = s.g
  sd : s₁.stackDir = s.stackDir.set! d (if d ≤ o.cls.lowval d then false else !s.stackDir[o.cls.lowval d]!)
  dirs : DirsOf s₁ d = DirsOf s d
  read : ∃ new', s₁.tstack = new' ++ base ∧
    StRead s₁.items new' (ps ++ if (!hasVert && decide (o.cls.lowval d < d) && o.cls.isType1) = true then
      [⟨if d ≤ o.cls.lowval d then false else !s.stackDir[o.cls.lowval d]!, [vertItem v]⟩] else []) ∧
    StItems g s₁ blocks
  tstack : ∀ t ∈ s₁.tstack, t.vStart = v ∨ t ∈ s.tstack
  nopush : (!hasVert && decide (o.cls.lowval d < d) && o.cls.isType1) = false → s₁.tstack = s.tstack

theorem walkOutPre_pre {g : Graph} {blocks : List StBlock} {base new : List TEntry}
    {ps : List StPiece} (s : WalkState) (v d : Nat) (o : DfsOut) (hasVert : Bool)
    (hd : d < s.stackDir.size) (hts : s.tstack = new ++ base) (hR : StRead s.items new ps)
    (hI : StItems g s blocks)
    (hvert : hasVert = false →
      vertItem v < s.items.size ∧ Items.type s.items (vertItem v) = .V ∧
      (∀ p, ¬ Items.IsParent s.items p (vertItem v)) ∧ vertItem v ∉ readStack s.tstack) :
    wp (walkOutPre v d o hasVert) (PreSpec g s v d o hasVert base ps blocks) s := by
  obtain ⟨h1, h2, h3, -, h5, h6, h7⟩ := walkOutPre_st s v d o hasVert hd hts hR hI hvert
  refine ⟨h1, h2, h3, h5, h6, h7, ?_, ?_⟩
  · intro t ht
    rw [walkOutPre_run] at ht
    split at ht
    · rcases List.mem_cons.1 ht with rfl | ht
      · exact Or.inl rfl
      · exact Or.inr ht
    · exact Or.inr ht
  · intro hp
    show ((walkOutPre v d o hasVert).run s).2.tstack = s.tstack
    rw [walkOutPre_run, ite_eq_right (by simp [hp])]

/-! ## Reference pieces of a returning edge, in the implementation's form -/

theorem refOut_ret_back' {g : Graph} {v d : Nat} {dirs : List Bool} {e dest : Nat} {cls : OutClass}
    {hv : Bool} (h : cls.lowval d < d) :
    refOut g v d dirs (.back e dest cls) hv =
      ((if (!hv && cls.isType1) = true then [⟨!dirs.getD (cls.lowval d) false, [vertItem v]⟩] else []) ++
        [⟨dirs.getD (cls.lowval d) false, [edgeItem g e]⟩] ++
        (if (hv || cls.isType1) = true then [] else [⟨!dirs.getD (cls.lowval d) false, [vertItem v]⟩]),
       [], true) := by
  cases hv <;> cases hT : cls.isType1 <;> simp [refOut, DfsOut.cls, Nat.not_le.mpr h, hT]

theorem refOut_ret_tree' {g : Graph} {v d : Nat} {dirs : List Bool} {e : Nat} {cls : OutClass}
    {child : DfsTree} {hv : Bool} (h : cls.lowval d < d) :
    refOut g v d dirs (.tree e cls child) hv =
      ((if (!hv && cls.isType1) = true then [⟨!dirs.getD (cls.lowval d) false, [vertItem v]⟩] else []) ++
        (if (hv || cls.isType1) = true then
          [⟨dirs.getD (cls.lowval d) false, stNest
            ((refTree g child (d + 1) (dirs ++ [!dirs.getD (cls.lowval d) false])).1 ++
              [⟨!dirs.getD (cls.lowval d) false, [edgeItem g e]⟩])⟩]
        else (refTree g child (d + 1) (dirs ++ [!dirs.getD (cls.lowval d) false])).1 ++
              [⟨!dirs.getD (cls.lowval d) false, [edgeItem g e]⟩]) ++
        (if (hv || cls.isType1) = true then [] else [⟨!dirs.getD (cls.lowval d) false, [vertItem v]⟩]),
       (refTree g child (d + 1) (dirs ++ [!dirs.getD (cls.lowval d) false])).2, true) := by
  cases hv <;> cases hT : cls.isType1 <;> simp [refOut, DfsOut.cls, Nat.not_le.mpr h, hT]

/-! ## `finishEdge` of a returning edge at the simulation state -/

/-- The stack above `segsStack segs` reads `finishEdge_st`'s pieces afterwards; the base segments,
the blocks and the directions `≤ d` are kept. -/
theorem stRet_finish {g : Graph} {P X : ItemId → Prop} {D d : Nat} {o : DfsOut} {s : WalkState} {v : Nat}
    {hv₁ : Bool} {sub new₁ : List TEntry} {segs : List (List TEntry × List StPiece)}
    {ps qs : List StPiece} {blocks : List StBlock}
    (hE : s.EarFinish v d o hv₁ sub (new₁ ++ segsStack segs)) (hi : s.Inv' D) (hs : Shape s)
    (hD : D = if o.cls.isTree then d + 1 else d)
    (hg : FinishGuards d o (new₁ ++ segsStack segs).length hv₁ s)
    (hb : FinishBook v d o (new₁ ++ segsStack segs).length hv₁ s)
    (hlt : o.cls.lowval d < d)
    (hB : ∀ t ∈ segsStack segs, t.vStart ≠ v)
    (hpre : hv₁ = false → new₁ = [] ∧ qs = [])
    (hR : StRead s.items sub ps) (hRq : StRead s.items new₁ qs) (hI : StItems g s blocks)
    (hseg : SegRead s.items segs) (hfull : s.Full g P X) :
    let r := (finishEdge v d o (new₁ ++ segsStack segs).length hv₁).run s
    r.1 = true ∧ r.2.g = s.g ∧ DirsOf r.2 d = DirsOf s d ∧
    (∃ new', r.2.tstack = new' ++ segsStack segs ∧
      StRead r.2.items new'
        (qs ++
          (if o.cls.isTree then
            (if hv₁ then
              [⟨!s.stackDir[d]!, stNest (ps ++ [⟨s.stackDir[d]!, [edgeItem s.g o.e]⟩])⟩]
            else ps ++ [⟨s.stackDir[d]!, [edgeItem s.g o.e]⟩])
          else [⟨!s.stackDir[d]!, [edgeItem s.g o.e]⟩]) ++
          (if hv₁ then [] else [⟨s.stackDir[d]!, [vertItem v]⟩]))) ∧
    StItems g r.2 blocks ∧ SegRead r.2.items segs ∧
    ∀ t ∈ r.2.tstack, t.vStart = v ∨ (o.cls.isTree = true ∧ o.cls.lowval d < d ∧ t.vStart = o.dest) ∨
      ∃ t₀ ∈ s.tstack, t.vStart = t₀.vStart := by
  intro r
  obtain ⟨lv, kind, ho, hlow⟩ := WalkState.ret_of_lowval_lt hlt
  have hlv : o.cls.lowval d = lv := by rw [ho]; rfl
  have hok := finishOk_of_guards ho hlow hg hE rfl hi hs hD hb.v_lt hb.e_lt hb.q
    (hb.ends lv kind ho) hb.vert
  have hfr := finishRet_frame_st hE hfull hI hb.v_lt ho hlow hB
  have hsd : s.stackDir[d]! = !s.stackDir[lv]! := hlv ▸ hE.dir_d hlt
  obtain ⟨h1, h2, h3, new', h4, h5, h6⟩ := finishEdge_st hE hi hs hD hok ho hlow hb.v_lt hb.e_lt hb.q
    (hb.ends lv kind ho) (hs.edge _ hb.e_lt) hsd hB hpre hfr.1 hR hRq hI
  exact ⟨h1, h2, DirsOf_congr fun k hk => h3 k (Nat.le_of_lt hk), ⟨new', h4, h5⟩, h6,
    hseg.congr hfr.2, finishEdge_vStart v d o _ hv₁ s⟩

theorem simBlocks_frame {g : Graph} {prev : List DfsTree} {fs : List PathFrame} {s : WalkState}
    {d : Nat} (hd : d = fs.length) (f : PathFrame) (b : Bool) :
    simBlocks g prev (fs ++ [f]) (DirsOf s d ++ [b]) =
      simBlocks g prev fs (DirsOf s d) ++ (refOuts g f.v d (DirsOf s d) f.done false).2.1 := by
  subst hd
  rw [simBlocks_snoc, List.take_left' (DirsOf_length s _),
    simBlocks_congr (dirs' := DirsOf s fs.length) fun j hj =>
      List.take_append_of_le_length (by rw [DirsOf_length]; omega)]

theorem DirsOf_push (s : WalkState) (d : Nat) (b : Bool) (t : TEntry) :
    DirsOf { s with stackDir := s.stackDir.set! d b, tstack := t :: s.tstack } d = DirsOf s d :=
  DirsOf_congr fun k hk => Array.getElem!_set!_ne _ _ _ _ (by omega)

theorem mem_cons_append_left {a : Nat} {v : Nat} {l l' : List Nat} (h : a ∈ v :: l) : a ∈ v :: (l ++ l') := by
  simp only [List.mem_cons, List.mem_append] at h ⊢; tauto

/-! ## The one-edge step -/

theorem stOut_step (g : Graph) (v d : Nat) (o : DfsOut) (hasVert : Bool)
    (ih : match o with | .tree _ _ child => StTreeP g child (d + 1) | .back .. => True) :
    StOutP g v d o hasVert := by
  intro s prev fs segs done P X hd hi hs hg hb hh
  obtain ⟨new, hts, hR, hnew⟩ := hh.pre.read
  have hgeq : s.g = g := hh.full.place.g_eq
  have hvlt : vertItem v < 1 + g.nv + g.ne := by show 1 + v < _; have := hh.v_lt; omega
  have hvert : hasVert = false → vertItem v < s.items.size ∧ Items.type s.items (vertItem v) = .V ∧
      (∀ p, ¬ Items.IsParent s.items p (vertItem v)) ∧ vertItem v ∉ readStack s.tstack := by
    intro hhv
    obtain ⟨hroot, hoff⟩ := Spqr.Place.fresh hh.full.place (by show 0 < 1 + v; omega) hvlt (hh.Pcur hhv)
    have hsz := hs.size
    rw [hgeq] at hsz
    have hty := hs.vert v (hgeq ▸ hh.v_lt)
    exact ⟨Nat.lt_of_lt_of_le hvlt hsz, hty, hroot, hoff⟩
  have hd₀ : d < s.stackDir.size := by have := hh.height; omega
  have hpreW := walkOutPre_pre s v d o hasVert hd₀ hts hR hh.pre.items hvert
  rw [walkOut_eq, wp_bind]
  unfold GuardsOut at hg; unfold BookOut at hb
  refine wp_imp (wp_imp (wp_imp (wp_imp (wp_imp (wp_of_forall
    fun hv₁ s₁ hp₁ ⟨hi₁, hs₁⟩ hg₁ hb₁ ⟨hfull₁, hmono₁⟩ => ?_) hpreW)
    (walkOutPre_inv hi hs hb.1)) hg) hb.2) (walkOutPre_full hh.full hh.v_lt hh.Pcur hh.Pcur')
  obtain ⟨new₁, hts₁, hR₁, hI₁⟩ := hp₁.read
  have hhv₁ := hp₁.hv
  have hsd₁ := hp₁.sd
  have hnopush := hp₁.nopush
  have hdirs₁ : DirsOf s₁ d = DirsOf s d := hp₁.dirs
  have hvs₁ : ∀ t ∈ s₁.tstack,
      t.vStart ∈ v :: DfsOut.vertsList done ∨ ∃ t₀ ∈ segsStack segs, t.vStart = t₀.vStart :=
    fun t ht => (hp₁.tstack t ht).elim (fun h => Or.inl (h ▸ List.mem_cons_self)) (hh.pre.vStart t)
  have hseg₁ : SegRead s₁.items segs := by rw [hp₁.items]; exact hh.pre.segRead
  have hgeq₁ : s₁.g = s.g := hp₁.geq
  have hsgeq : s₁.g = g := hp₁.geq.trans hgeq
  clear hp₁
  generalize hx : (if d ≤ o.cls.lowval d then false else !s.stackDir[o.cls.lowval d]!) = x at hsd₁ hR₁
  generalize hpush : (!hasVert && decide (o.cls.lowval d < d) && o.cls.isType1) = push
    at hhv₁ hR₁ hnopush
  have hxd : s₁.stackDir[d]! = x := by rw [hsd₁]; exact Array.getElem!_set!_self _ _ _ hd₀
  have hsize₁ : s₁.stackDir.size = s.stackDir.size := by rw [hsd₁]; simp
  have hbelow₁ : ∀ k, k < d → s₁.stackDir[k]! = s.stackDir[k]! := fun k hk => by
    rw [hsd₁]; exact Array.getElem!_set!_ne _ _ _ _ (Nat.ne_of_gt hk)
  have hB : ∀ t ∈ segsStack segs, t.vStart ≠ v := fun t ht h => hh.base t ht (h ▸ List.mem_cons_self)
  have hnew₁ : push = false → hasVert = false → new₁ = [] := fun hpf h =>
    (List.append_cancel_right (hts₁.symm.trans ((hnopush hpf).trans hts))).trans (hnew h)
  unfold walkOutRest
  rw [bind_tstackSize]
  cases o with
  | back e dest cls =>
    dsimp only at hg₁ hb₁ ⊢
    obtain ⟨sub, base, hlen, hE⟩ := hb₁.ear
    have hnt : (DfsOut.back e dest cls).cls.isTree = false := by
      cases h : (DfsOut.back e dest cls).cls.isTree
      · rfl
      · obtain ⟨_, _, _, h'⟩ := hb₁.tree.1 h; cases h'
    have hsub : sub = [] := hE.back_nil hnt
    subst hsub
    have hbase : base = new₁ ++ segsStack segs := by simpa [hts₁] using hE.tstack.symm
    subst hbase
    have hD : d = if (DfsOut.back e dest cls).cls.isTree then d + 1 else d := by simp [hnt]
    rw [hts₁] at hg₁ hb₁ ⊢
    have hvl : DfsOut.vertsList (done ++ [DfsOut.back e dest cls]) = DfsOut.vertsList done := by
      simp [DfsOut.vertsList_append, DfsOut.vertsList]
    show _ ∧ _ ∧ StPre g prev fs segs v d (done ++ [.back e dest cls]) _ _
    by_cases hge : d ≤ (DfsOut.back e dest cls).cls.lowval d
    · have hge' : d ≤ cls.lowval d := hge
      have hok := ear_boundary hge hg₁ hi₁ hs₁ hb₁ hD
      obtain ⟨hr1, hrg, hrsd, hrts, hrI, hrfr⟩ :=
        finishBoundary_st hE hi₁ hs₁ hD hok hb₁ hge hsgeq StRead.nil hI₁
      have hpf : push = false := by rw [← hpush]; simp [Nat.not_lt.mpr hge]
      have hdr : DirsOf ((finishEdge v d (.back e dest cls) (new₁ ++ segsStack segs).length hv₁).run s₁).2 d
          = DirsOf s d := by
        rw [DirsOf_congr (s := s₁) (fun k _ => by rw [hrsd]), hdirs₁]
      rw [hpf] at hR₁
      simp only [Bool.false_eq_true, ↓reduceIte, List.append_nil] at hR₁
      rw [hnt] at hrI
      simp only [Bool.false_eq_true, ↓reduceIte, List.append_nil] at hrI
      refine ⟨hrg.trans hgeq₁, (KeepsSize.finishEdge _ _ _ _ _ s₁).trans hsize₁,
        ⟨new₁, hrts, ?_, fun h => hnew₁ hpf ?_⟩, ?_, ?_, ?_, ?_⟩
      · rw [hdr, refOuts_snoc, hh.pre.hv, refOut_boundary_back hge']
        simp only [List.append_nil]
        exact StRead.congr (fun x hx => hrfr x (mem_readStack_append.2 (Or.inl hx))) hR₁
      · rw [hr1, hhv₁, hpf, Bool.false_or] at h; exact h
      · rw [hdr, refOuts_snoc, hh.pre.hv, refOut_boundary_back hge', hr1, hhv₁, hpf, Bool.false_or]
      · exact hseg₁.congr fun x hx => hrfr x (mem_readStack_append.2 (Or.inr hx))
      · rw [hrts, ← hts₁, hvl]; exact hvs₁
      · rw [hdr, refOuts_snoc, hh.pre.hv, refOut_boundary_back hge']
        simpa only [List.append_nil] using hrI
    · have hlt : (DfsOut.back e dest cls).cls.lowval d < d := Nat.lt_of_not_le hge
      have hlt' : cls.lowval d < d := hlt
      have hge' : ¬ d ≤ cls.lowval d := hge
      have hpr : push = (!hasVert && cls.isType1) := by rw [← hpush]; simp [DfsOut.cls, hlt']
      have hhv : hv₁ = (hasVert || cls.isType1) := by
        rw [hhv₁, hpr]; cases hasVert <;> cases cls.isType1 <;> rfl
      have hxv : x = !s.stackDir[cls.lowval d]! := by rw [← hx]; simp [DfsOut.cls, hge']
      have hpre' : hv₁ = false → new₁ = [] ∧
          ((refOuts g v d (DirsOf s d) done false).1 ++
            if push = true then [⟨x, [vertItem v]⟩] else []) = [] := by
        intro h
        rw [hhv, Bool.or_eq_false_iff] at h
        have hpf : push = false := by simp [hpr, h.1, h.2]
        refine ⟨hnew₁ hpf h.1, ?_⟩
        rw [hpf, (refOuts_hv_false _ _ (hh.pre.hv.trans h.1)).2]; rfl
      obtain ⟨hr1, hrg, hdr', ⟨new', hrts, hR'⟩, hrI, hseg', hvs'⟩ :=
        stRet_finish (D := d) hE hi₁ hs₁ hD hg₁ hb₁ hlt hB hpre' StRead.nil hR₁ hI₁ hseg₁ hfull₁
      rw [hdirs₁] at hdr'
      rw [hxd, hxv, Bool.not_not, hsgeq, hpr, hnt] at hR'
      simp only [Bool.false_eq_true, ↓reduceIte] at hR'
      subst hhv
      refine ⟨hrg.trans hgeq₁, (KeepsSize.finishEdge _ _ _ _ _ s₁).trans hsize₁,
        ⟨new', hrts, ?_, fun h => by rw [hr1] at h; cases h⟩, ?_, hseg', ?_, ?_⟩
      · rw [hdr', refOuts_snoc, hh.pre.hv, refOut_ret_back' hlt', DirsOf_getD s hlt']
        simpa only [List.append_assoc, DfsOut.e] using hR'
      · rw [hdr', refOuts_snoc, hh.pre.hv, refOut_ret_back' hlt', hr1]
      · intro t ht
        rw [hvl]
        rcases hvs' t ht with h | ⟨hT, -, -⟩ | ⟨t₀, ht₀, h⟩
        · exact Or.inl (h ▸ List.mem_cons_self)
        · exact Bool.noConfusion ((show cls.isTree = true from hT).symm.trans (show cls.isTree = false from hnt))
        · exact (hvs₁ t₀ ht₀).imp (h ▸ id) fun ⟨t₁, ht₁, h₁⟩ => ⟨t₁, ht₁, h.trans h₁⟩
      · rw [hdr', refOuts_snoc, hh.pre.hv, refOut_ret_back' hlt']
        simpa only [List.append_nil] using hrI
  | tree e cls child =>
    dsimp only at hg₁ hb₁ ih ⊢
    simp only [wp_bind, wp_modify] at hg₁ hb₁ ⊢
    set s₂ : WalkState := { s₁ with firstOccurrence := s₁.firstOccurrence.set! d s₁.g.ne } with hs₂
    have hvl : DfsOut.vertsList (done ++ [DfsOut.tree e cls child]) = DfsOut.vertsList done ++ child.verts := by
      simp [DfsOut.vertsList_append, DfsOut.vertsList]
    have hvn := hh.vnodup
    rw [hvl] at hvn
    obtain ⟨hvnot, hvn'⟩ := List.nodup_cons.1 hvn
    have hel := hh.e_lt
    have hen := hh.enodup
    have hPe := hh.Pe
    simp only [DfsOut.edgesList, List.append_nil, List.mem_cons, List.nodup_cons, forall_eq_or_imp] at hel hen hPe
    have hwl := hh.w_lt
    have hPv := hh.Pv
    simp only [DfsOut.vertsList, List.append_nil] at hwl hPv
    have hhe := hh.height
    simp only [DfsOut.heightList, Nat.max_zero] at hhe
    set qs₁ := (refOuts g v d (DirsOf s d) done false).1 ++
      (if push = true then [⟨x, [vertItem v]⟩] else []) with hqs₁
    have pre : ∀ w outs, child = .node w outs →
        ({ s₂ with stackVerts := s₂.stackVerts.set! (d + 1) w } : WalkState).Inv' (d + 1) :=
      fun w outs _ => (hi₁.frame' (s' := s₂)).setSv w
    have hkb : wp (walkTree child (d + 1))
        (fun _ s₃ => ∀ k, k < d + 1 → s₃.stackDir[k]! = s₂.stackDir[k]!) s₂ :=
      fun k hk => walk_keepsBelow.1 child (d + 1) s₂ k hk
    have hdirs₂ : DirsOf s₂ (d + 1) = DirsOf s d ++ [x] := by
      rw [show DirsOf s₂ (d + 1) = DirsOf s₁ (d + 1) from rfl, DirsOf_succ, hxd, hdirs₁]
    have hih := ih s₂ prev (fs ++ [⟨v, done, .tree e cls child⟩]) ((new₁, qs₁) :: segs)
      (fun i => P i ∨ (hv₁ = true ∧ i = vertItem v)) X (by simp [hd]) pre hs₁.frame' hg₁.1 hb₁.1
      (hfull₁.of_eq rfl rfl rfl) hwl hel.2 hvn'.of_append_right hen.2
      (fun w hw h => h.elim (hPv w hw) fun h =>
        hvnot (List.mem_append_right _ (WalkM.vertItem_inj h.2 ▸ hw)))
      (fun e' he' h => h.elim (hPe.2 e' he') fun h => WalkM.edgeItem_ne_vertItem hh.v_lt e' h.2)
      (by rw [show s₂.stackDir = s₁.stackDir from rfl, hsize₁]; omega)
      (by rw [segsStack_cons]; exact hts₁) (SegRead.cons hR₁ hseg₁)
      (fun t' ht' hmem => by
        rw [segsStack_cons, ← hts₁] at ht'
        rcases hvs₁ t' ht' with h | ⟨t₀, ht₀, h⟩
        · exact List.disjoint_of_nodup_append (l₁ := v :: DfsOut.vertsList done)
            (by simpa using hvn) h hmem
        · refine hh.base t₀ ht₀ ?_
          rw [← h, hvl]
          exact List.mem_cons_of_mem _ (List.mem_append_right _ hmem))
      (by rw [hdirs₂, simBlocks_frame hd]; exact hI₁.congr rfl rfl)
    have hfullT := (walk_full_aux g).1 child (d + 1) (fun i => P i ∨ (hv₁ = true ∧ i = vertItem v)) X s₂
      (hfull₁.of_eq rfl rfl rfl) hwl hel.2 hvn'.of_append_right hen.2
      (fun w hw h => h.elim (hPv w hw) fun h =>
        hvnot (List.mem_append_right _ (WalkM.vertItem_inj h.2 ▸ hw)))
      (fun e' he' h => h.elim (hPe.2 e' he') fun h => WalkM.edgeItem_ne_vertItem hh.v_lt e' h.2)
      (sdTree child (d + 1) s₂ pre hs₁.frame' hg₁.1 hb₁.1)
    refine wp_imp (wp_imp (wp_imp (wp_imp (wp_imp (wp_of_forall
      fun _ s₃ hfull₃ ⟨hg₃, hb₃⟩ ⟨hi₃, hs₃⟩ hkb₃ ⟨hg₃eq, hsz₃, hvs₃, hseg₃, hsim₃⟩ => ?_) hfullT)
      (wp_and hg₁.2 hb₁.2)) (invTree child (d + 1) _ pre hs₁.frame' hg₁.1 hb₁.1)) hkb) hih
    have hsgeq₃ : s₃.g = g := hg₃eq.trans hsgeq
    have hgeq₃ : s₃.g = s.g := hg₃eq.trans hgeq₁
    have hsz₃' : s₃.stackDir.size = s.stackDir.size := hsz₃.trans hsize₁
    have hdirs₃ : DirsOf s₃ (d + 1) = DirsOf s d ++ [x] := by
      rw [DirsOf_congr (s := s₂) hkb₃, hdirs₂]
    have hdirs₃' : DirsOf s₃ d = DirsOf s d := by
      rw [DirsOf_congr (s := s₂) (fun k hk => hkb₃ k (Nat.lt_succ_of_lt hk))]; exact hdirs₁
    have hxd₃ : s₃.stackDir[d]! = x := by rw [hkb₃ d (Nat.lt_succ_self d)]; exact hxd
    obtain ⟨new₃, hts₃, hR₃⟩ := hsim₃.read
    have hI₃ := hsim₃.items
    simp only [List.length_append, List.length_singleton, ← hd, hdirs₃, segsStack_cons] at hts₃ hR₃ hI₃
    rw [simBlocks_frame hd] at hI₃
    have hRq₃ : StRead s₃.items new₁ qs₁ := hseg₃ _ List.mem_cons_self
    have hseg₃' : SegRead s₃.items segs := fun sg hsg => hseg₃ sg (List.mem_cons_of_mem _ hsg)
    rw [hts₁] at hg₃ hb₃ ⊢
    have hit : (DfsOut.tree e cls child).cls.isTree = true := hb₃.tree.2 ⟨_, _, _, rfl⟩
    obtain ⟨sub, base, hlen, hE⟩ := hb₃.ear
    obtain ⟨hsub, hbase⟩ := List.append_inj' (hE.tstack.symm.trans hts₃) hlen
    subst hsub hbase
    have hD : d + 1 = if (DfsOut.tree e cls child).cls.isTree then d + 1 else d := by simp [hit]
    have hvs₃' : ∀ t ∈ s₃.tstack,
        t.vStart ∈ v :: DfsOut.vertsList (done ++ [DfsOut.tree e cls child]) ∨
          ∃ t₀ ∈ segsStack segs, t.vStart = t₀.vStart := by
      intro t ht
      rw [hvl]
      rcases hvs₃ t ht with h | ⟨t₀, ht₀, h⟩
      · exact Or.inl (List.mem_cons_of_mem _ (List.mem_append_right _ h))
      · rw [segsStack_cons, ← hts₁] at ht₀
        exact (hvs₁ t₀ ht₀).imp (fun h' => h ▸ mem_cons_append_left h')
          fun ⟨t₁, ht₁, h₁⟩ => ⟨t₁, ht₁, h.trans h₁⟩
    show _ ∧ _ ∧ StPre g prev fs segs v d (done ++ [.tree e cls child]) _ _
    by_cases hge : d ≤ (DfsOut.tree e cls child).cls.lowval d
    · have hge' : d ≤ cls.lowval d := hge
      have hxv : x = false := by rw [← hx]; simp [DfsOut.cls, hge']
      have hok := ear_boundary hge hg₃ hi₃ hs₃ hb₃ hD
      subst hxv
      obtain ⟨hr1, hrg, hrsd, hrts, hrI, hrfr⟩ :=
        finishBoundary_st hE hi₃ hs₃ hD hok hb₃ hge hsgeq₃ hR₃ hI₃
      have hpf : push = false := by rw [← hpush]; simp [DfsOut.cls, Nat.not_lt.mpr hge']
      have hdr : DirsOf ((finishEdge v d (.tree e cls child) (new₁ ++ segsStack segs).length hv₁).run s₃).2 d
          = DirsOf s d := by
        rw [DirsOf_congr (s := s₃) (fun k _ => by rw [hrsd]), hdirs₃']
      rw [hqs₁, hpf] at hRq₃
      simp only [Bool.false_eq_true, ↓reduceIte, List.append_nil] at hRq₃
      rw [hit] at hrI
      simp only [↓reduceIte] at hrI
      refine ⟨hrg.trans hgeq₃, (KeepsSize.finishEdge _ _ _ _ _ s₃).trans hsz₃',
        ⟨new₁, hrts, ?_, fun h => hnew₁ hpf ?_⟩, ?_, ?_, ?_, ?_⟩
      · rw [hdr, refOuts_snoc, hh.pre.hv, refOut_boundary_tree hge']
        simp only [List.append_nil]
        exact StRead.congr (fun x hx => hrfr x (mem_readStack_append.2 (Or.inl hx))) hRq₃
      · rw [hr1, hhv₁, hpf, Bool.false_or] at h; exact h
      · rw [hdr, refOuts_snoc, hh.pre.hv, refOut_boundary_tree hge', hr1, hhv₁, hpf, Bool.false_or]
      · exact hseg₃'.congr fun x hx => hrfr x (mem_readStack_append.2 (Or.inr hx))
      · rw [hrts, ← hts₁]
        intro t ht; rw [hvl]; exact (hvs₁ t ht).imp_left mem_cons_append_left
      · rw [hdr, refOuts_snoc, hh.pre.hv, refOut_boundary_tree hge']
        simpa only [List.append_assoc, DfsOut.dest] using hrI
    · have hlt : (DfsOut.tree e cls child).cls.lowval d < d := Nat.lt_of_not_le hge
      have hlt' : cls.lowval d < d := hlt
      have hge' : ¬ d ≤ cls.lowval d := hge
      have hpr : push = (!hasVert && cls.isType1) := by rw [← hpush]; simp [DfsOut.cls, hlt']
      have hhv : hv₁ = (hasVert || cls.isType1) := by
        rw [hhv₁, hpr]; cases hasVert <;> cases cls.isType1 <;> rfl
      have hxv : x = !s.stackDir[cls.lowval d]! := by rw [← hx]; simp [DfsOut.cls, hge']
      subst hxv
      have hpre' : hv₁ = false → new₁ = [] ∧ qs₁ = [] := by
        intro h
        rw [hhv, Bool.or_eq_false_iff] at h
        have hpf : push = false := by simp [hpr, h.1, h.2]
        refine ⟨hnew₁ hpf h.1, ?_⟩
        rw [hqs₁, hpf, (refOuts_hv_false _ _ (hh.pre.hv.trans h.1)).2]; rfl
      obtain ⟨hr1, hrg, hdr', ⟨new', hrts, hR'⟩, hrI, hseg', hvs'⟩ :=
        stRet_finish (D := d + 1) hE hi₃ hs₃ hD hg₃ hb₃ hlt hB hpre' hR₃ hRq₃ hI₃ hseg₃' hfull₃
      rw [hdirs₃'] at hdr'
      rw [hxd₃, hsgeq₃, hit, hqs₁, hpr] at hR'
      simp only [↓reduceIte, Bool.not_not] at hR'
      subst hhv
      refine ⟨hrg.trans hgeq₃, (KeepsSize.finishEdge _ _ _ _ _ s₃).trans hsz₃',
        ⟨new', hrts, ?_, fun h => by rw [hr1] at h; cases h⟩, ?_, hseg', ?_, ?_⟩
      · rw [hdr', refOuts_snoc, hh.pre.hv, refOut_ret_tree' hlt', DirsOf_getD s hlt']
        simpa only [List.append_assoc, DfsOut.e] using hR'
      · rw [hdr', refOuts_snoc, hh.pre.hv, refOut_ret_tree' hlt', hr1]
      · intro t ht
        rcases hvs' t ht with h | ⟨-, -, h⟩ | ⟨t₀, ht₀, h⟩
        · exact Or.inl (h ▸ List.mem_cons_self)
        · refine Or.inl (List.mem_cons_of_mem _ ?_)
          rw [h, DfsOut.vertsList_append]
          exact List.mem_append_right _ (by cases child; simp [DfsOut.dest, DfsTree.verts, DfsTree.v])
        · exact (hvs₃' t₀ ht₀).imp (h ▸ id) fun ⟨t₁, ht₁, h₁⟩ => ⟨t₁, ht₁, h.trans h₁⟩
      · rw [hdr', refOuts_snoc, hh.pre.hv, refOut_ret_tree' hlt', DirsOf_getD s hlt']
        simpa only [List.append_assoc, List.append_nil] using hrI


/-! ## The out-edge list and the tree node -/

theorem stOuts_nil (g : Graph) (v d : Nat) (hasVert : Bool) : StOutsP g v d [] hasVert := by
  intro s prev fs segs done P X hd hi hs hg hb hh
  simp only [walkOuts, wp_pure, List.append_nil]
  exact ⟨trivial, trivial, hh.pre⟩

theorem stOuts_cons (g : Graph) (v d : Nat) (hasVert : Bool) (o : DfsOut) (rest : List DfsOut)
    (ih₁ : StOutP g v d o hasVert) (ih₂ : ∀ hv, StOutsP g v d rest hv) :
    StOutsP g v d (o :: rest) hasVert := by
  intro s prev fs segs done P X hd hi hs hg hb hh
  simp only [walkOuts, wp_bind]
  have hvl : DfsOut.vertsList (o :: rest) = DfsOut.vertsList [o] ++ DfsOut.vertsList rest := by
    rw [← List.singleton_append, DfsOut.vertsList_append]
  have hel : DfsOut.edgesList (o :: rest) = DfsOut.edgesList [o] ++ DfsOut.edgesList rest := by
    rw [← List.singleton_append, DfsOut.edgesList_append]
  have hwl := hh.w_lt; have helt := hh.e_lt; have hvn := hh.vnodup; have hen := hh.enodup
  have hPv := hh.Pv; have hPe := hh.Pe
  rw [hvl] at hwl hPv; rw [hel] at helt hPe hen
  rw [DfsOut.vertsList_append, hvl] at hvn
  obtain ⟨hvnot, hvn'⟩ := List.nodup_cons.1 hvn
  have hvo : v ∉ DfsOut.vertsList [o] := fun h => hvnot (List.mem_append_right _ (List.mem_append_left _ h))
  have hfull := (walk_full_aux g).2.2 v d o hasVert P X s hh.full hh.v_lt
    (fun w hw => hwl w (List.mem_append_left _ hw)) (fun e he => helt e (List.mem_append_left _ he))
    (hvn'.of_append_right.of_append_left) hen.of_append_left hvo
    (fun w hw => hPv w (List.mem_append_left _ hw)) (fun e he => hPe e (List.mem_append_left _ he))
    hh.Pcur hh.Pcur' (sdOut v d o hasVert s hi hs hg.1 hb.1)
  have hh₁ : StHyps g prev fs segs v d done [o] hasVert P X s := {
    full := hh.full
    v_lt := hh.v_lt
    w_lt := fun w hw => hwl w (List.mem_append_left _ hw)
    e_lt := fun e he => helt e (List.mem_append_left _ he)
    vnodup := by
      rw [DfsOut.vertsList_append]
      exact hvn.sublist (List.cons_sublist_cons.2 ((List.sublist_append_left _ _).append_left _))
    enodup := hen.of_append_left
    Pv := fun w hw => hPv w (List.mem_append_left _ hw)
    Pe := fun e he => hPe e (List.mem_append_left _ he)
    Pcur := hh.Pcur
    Pcur' := hh.Pcur'
    height := by
      refine Nat.lt_of_le_of_lt ?_ hh.height
      cases o <;> simp [DfsOut.heightList]
    base := fun t ht h => hh.base t ht (by
      simp only [List.mem_cons, DfsOut.vertsList_append, hvl, List.mem_append] at h ⊢; tauto)
    pre := hh.pre }
  refine wp_imp (wp_imp (wp_imp (wp_imp (wp_imp (wp_of_forall
    fun hv₁ s₁ ⟨hfull₁, hhv₁⟩ hb₁ hg₁ ⟨hi₁, hs₁⟩ ⟨hgeq₁, hsize₁, hp₁⟩ => ?_) hfull) hb.2) hg.2)
    (invOut v d o hasVert s hi hs hg.1 hb.1))
    (ih₁ s prev fs segs done P X hd hi hs hg.1 hb.1 hh₁)
  have hdis' : ∀ w ∈ DfsOut.vertsList rest, ∀ w' ∈ DfsOut.vertsList [o], w ≠ w' :=
    fun w hw w' hw' h => (List.disjoint_of_nodup_append hvn'.of_append_right) hw' (h ▸ hw)
  have hh' : StHyps g prev fs segs v d (done ++ [o]) rest hv₁
      (fun i => Pushed g P (DfsOut.vertsList [o]) (DfsOut.edgesList [o]) i ∨ (hv₁ = true ∧ i = vertItem v))
      X s₁ := {
    full := hfull₁.monoP fun i => by
      simp only [Pushed]
      constructor
      · rintro ((h | h) | h | h)
        · exact Or.inl (Or.inl h)
        · exact Or.inr h
        · exact Or.inl (Or.inr (Or.inl h))
        · exact Or.inl (Or.inr (Or.inr h))
      · rintro ((h | h | h) | h)
        · exact Or.inl (Or.inl h)
        · exact Or.inr (Or.inl h)
        · exact Or.inr (Or.inr h)
        · exact Or.inl (Or.inr h)
    v_lt := hh.v_lt
    w_lt := fun w hw => hwl w (List.mem_append_right _ hw)
    e_lt := fun e he => helt e (List.mem_append_right _ he)
    vnodup := by rw [List.append_assoc, List.singleton_append]; exact hh.vnodup
    enodup := hen.of_append_right
    Pv := fun w hw h => by
      rcases h with (h | ⟨w', hw', h⟩ | ⟨e', _, h⟩) | ⟨_, h⟩
      · exact hPv w (List.mem_append_right _ hw) h
      · exact hdis' w hw w' hw' (WalkM.vertItem_inj h)
      · exact WalkM.vertItem_ne_edgeItem (hwl w (List.mem_append_right _ hw)) e' h
      · exact hvnot (List.mem_append_right _ (List.mem_append_right _ (WalkM.vertItem_inj h ▸ hw)))
    Pe := fun e he h => by
      rcases h with (h | ⟨w', hw', h⟩ | ⟨e', he', h⟩) | ⟨_, h⟩
      · exact hPe e (List.mem_append_right _ he) h
      · exact WalkM.edgeItem_ne_vertItem (hwl w' (List.mem_append_left _ hw')) e h
      · exact (List.disjoint_of_nodup_append hen) he' (WalkM.edgeItem_inj h ▸ he)
      · exact WalkM.edgeItem_ne_vertItem hh.v_lt e h
    Pcur := fun h0 h => by
      have hhv : hasVert = false := by
        cases hasVert
        · rfl
        · exact absurd (hhv₁ rfl) (by rw [h0]; decide)
      rcases h with (h | ⟨w', hw', h⟩ | ⟨e', _, h⟩) | ⟨h, _⟩
      · exact hh.Pcur hhv h
      · exact hvo (WalkM.vertItem_inj h ▸ hw')
      · exact WalkM.vertItem_ne_edgeItem hh.v_lt e' h
      · rw [h0] at h; cases h
    Pcur' := fun h => Or.inr ⟨h, rfl⟩
    height := by
      rw [hsize₁]
      refine Nat.lt_of_le_of_lt ?_ hh.height
      cases o <;> simp [DfsOut.heightList]
    base := fun t ht => by rw [List.append_assoc, List.singleton_append]; exact hh.base t ht
    pre := hp₁ }
  refine wp_imp (wp_of_forall fun hv' s' ⟨hg', hsz', hp'⟩ => ?_)
    (ih₂ hv₁ s₁ prev fs segs (done ++ [o]) _ X hd hi₁ hs₁ hg₁ hb₁ hh')
  rw [List.append_assoc, List.singleton_append] at hp'
  exact ⟨hg'.trans hgeq₁, hsz'.trans hsize₁, hp'⟩

theorem stTree_node (g : Graph) (d v : Nat) (outs : List DfsOut) (ih : StOutsP g v d outs false) :
    StTreeP g (.node v outs) d := by
  intro s prev fs segs P X hd hi hs hg hb hfull hvlt helt hvn hen hPv hPe hht hts hseg hbase hI
  simp only [walkTree, wp_bind, wp_modify]
  unfold GuardsTree at hg; unfold BookTree at hb
  simp only [DfsTree.verts, DfsTree.edges, DfsTree.height] at hvlt helt hvn hen hPv hPe hht hbase
  set s₁ : WalkState := { s with stackVerts := s.stackVerts.set! d v } with hs₁
  have hi₁ : s₁.Inv' d := hi v outs rfl
  have hs₁' : Shape s₁ := hs.frame' rfl rfl rfl
  obtain ⟨hvnot, hvn'⟩ := List.nodup_cons.1 hvn
  have hv : v < g.nv := hvlt v List.mem_cons_self
  have hgeq : s.g = g := hfull.place.g_eq
  have hh : StHyps g prev fs segs v d [] outs false P X s₁ := {
    full := hfull.of_eq rfl rfl rfl
    v_lt := hv
    w_lt := fun w hw => hvlt w (List.mem_cons_of_mem _ hw)
    e_lt := helt
    vnodup := by simpa using hvn
    enodup := hen
    Pv := fun w hw => hPv w (List.mem_cons_of_mem _ hw)
    Pe := hPe
    Pcur := fun _ => hPv v List.mem_cons_self
    Pcur' := fun h => by cases h
    height := by show d + _ < s.stackDir.size; omega
    base := fun t ht => by simpa using hbase t ht
    pre := {
      read := ⟨[], by simpa using hts, by rw [refOuts_nil]; exact StRead.nil, fun _ => rfl⟩
      hv := by rw [refOuts_nil]
      segRead := hseg
      vStart := fun t ht => Or.inr ⟨t, by rw [← hts]; exact ht, rfl⟩
      items := by rw [refOuts_nil, List.append_nil]; exact hI.congr rfl rfl } }
  have hfull₁ := (walk_full_aux g).2.1 v d outs false P X s₁ (hfull.of_eq rfl rfl rfl) hv
    (fun w hw => hvlt w (List.mem_cons_of_mem _ hw)) helt hvn' hen hvnot
    (fun w hw => hPv w (List.mem_cons_of_mem _ hw)) hPe (fun _ => hPv v List.mem_cons_self)
    (fun h => by cases h) (sdOuts v d outs false s₁ hi₁ hs₁' hg hb)
  refine wp_imp (wp_imp (wp_imp (wp_of_forall
    fun hv' s₂ ⟨hfull₂, _⟩ ⟨hi₂, hs₂, _⟩ ⟨hgeq₂, hsize₂, hp₂⟩ => ?_) hfull₁)
    (invOuts v d outs false s₁ hi₁ hs₁' hg hb))
    (ih s₁ prev fs segs [] P X hd hi₁ hs₁' hg hb hh)
  have hgeq₂' : s₂.g = s.g := hgeq₂
  have hsgeq₂ : s₂.g = g := hgeq₂'.trans hgeq
  have hsize₂' : s₂.stackDir.size = s.stackDir.size := hsize₂
  simp only [List.nil_append] at hp₂
  obtain ⟨new, hts₂, hR₂, hnew₂⟩ := hp₂.read
  have hhv₂ := hp₂.hv
  have hI₂ := hp₂.items
  have hvs₂ := hp₂.vStart
  have hvs : ∀ t ∈ s₂.tstack, t.vStart ∈ (DfsTree.node v outs).verts ∨
      ∃ t₀ ∈ segsStack segs, t.vStart = t₀.vStart := hvs₂
  cases hv' with
  | true =>
    simp only [↓reduceIte, wp_pure]
    refine ⟨hgeq₂', hsize₂', hvs, hp₂.segRead, ⟨new, hts₂, ?_⟩, ?_⟩
    · rw [refTree_node, ← hd, hhv₂]; exact hR₂
    · rw [refTree_node, ← hd]; exact hI₂
  | false =>
    simp only [Bool.false_eq_true, ↓reduceIte, wp_bind, wp_setStackDir, wp_pushVertTstack]
    have hdlt : d < s₂.stackDir.size := by rw [hsize₂']; omega
    rw [Array.getElem!_set!_self _ _ _ hdlt]
    have hnP : ¬ Pushed g (fun i => P i ∨ (false = true ∧ i = vertItem v))
        (DfsOut.vertsList outs) (DfsOut.edgesList outs) (vertItem v) := by
      rintro ((h | ⟨h, _⟩) | ⟨w, hw, hwe⟩ | ⟨e', _, hee⟩)
      · exact hPv v List.mem_cons_self h
      · cases h
      · exact hvnot (WalkM.vertItem_inj hwe ▸ hw)
      · exact WalkM.vertItem_ne_edgeItem hv e' hee
    obtain ⟨hroot, hoff⟩ := Spqr.Place.fresh hfull₂.place (by show 0 < 1 + v; omega)
      (by show 1 + v < 1 + g.nv + g.ne; omega) hnP
    have hty : Items.type s₂.items (vertItem v) = .V := hs₂.vert v (hsgeq₂ ▸ hv)
    have hsz : vertItem v < s₂.items.size := by
      have := hs₂.size; rw [hsgeq₂] at this; show 1 + v < _; omega
    have hdirs : DirsOf { s₂ with stackDir := s₂.stackDir.set! d true } d = DirsOf s₂ d :=
      DirsOf_congr fun k hk => Array.getElem!_set!_ne _ _ _ _ (by omega)
    refine ⟨hgeq₂', by simp [hsize₂'], ?_, hp₂.segRead,
      ⟨⟨v, d, s₂.nxtEdgeIdx, setSides true [vertItem v] []⟩ :: new, ?_, ?_⟩, ?_⟩
    · intro t ht
      rcases List.mem_cons.1 ht with rfl | ht
      · exact Or.inl List.mem_cons_self
      · exact hvs t ht
    · show _ :: s₂.tstack = _
      rw [hts₂, List.cons_append]
    · rw [refTree_node, ← hd, DirsOf_push, hhv₂]
      simp only [Bool.false_eq_true, ↓reduceIte]
      exact StRead.pushEntry v d s₂.nxtEdgeIdx true (vertItem v) (Or.inl hty) hR₂
    · rw [refTree_node, ← hd, DirsOf_push]
      exact (StItems.pushEntry v d s₂.nxtEdgeIdx true (vertItem v) hsz hroot hoff hI₂).congr rfl rfl


/-- The simulation of the subtree walk, out-edge list and out-edge (modulo `finishBoundary_st`,
`finishRet_frame_st`, `finishEdge_vStart`). -/
theorem stWalk (g : Graph) :
    (∀ (t : DfsTree) (d : Nat), StTreeP g t d) ∧
    (∀ (v d : Nat) (outs : List DfsOut) (hasVert : Bool), StOutsP g v d outs hasVert) ∧
    (∀ (v d : Nat) (o : DfsOut) (hasVert : Bool), StOutP g v d o hasVert) :=
  walkTree.mutual_induct _ _ _ (stTree_node g) (fun o v d hasVert ih => stOut_step g v d o hasVert ih)
    (stOuts_nil g) (fun v d hasVert o rest ih₁ ih₂ => stOuts_cons g v d hasVert o rest ih₁ ih₂)


/-! ## The forest -/

theorem refBlocks_snoc (g : Graph) (pre : List DfsTree) (t : DfsTree) :
    refBlocks g (pre ++ [t]) =
      refBlocks g pre ++ ((refTree g t 0 []).2 ++ [⟨none, stNest (refTree g t 0 []).1⟩]) := by
  simp [refBlocks, List.flatMap_append]

theorem DirsOf_zero (s : WalkState) : DirsOf s 0 = [] := by simp [DirsOf]

/-- The forest walk root by root: after the roots `pre`, every finished S / P / R item is in a block
of `refBlocks g pre` (modulo `rootPop_st` and the admissions of `stWalk`). -/
theorem stForest (g : Graph) : ∀ (forest pre : List DfsTree) (s : WalkState) (P X : ItemId → Prop),
    RootState g pre s → s.Full g P X →
    (∀ v ∈ forest.flatMap DfsTree.verts, ¬ P (vertItem v)) →
    (∀ e ∈ forest.flatMap DfsTree.edges, ¬ P (edgeItem g e)) →
    ForestOK g (pre ++ forest) → (∀ t ∈ forest, t.WF []) → (∀ t ∈ forest, t.Ends g) →
    (∀ t ∈ forest, t.height ≤ g.nv) →
    StItems g s (refBlocks g pre) →
    wp (walkForest forest) (fun _ s' => StItems g s' (refBlocks g (pre ++ forest))) s
  | [], pre, s, P, X, _, _, _, _, _, _, _, _, hI => by
    show StItems g s (refBlocks g (pre ++ []))
    simpa using hI
  | t :: rest, pre, s, P, X, h, hfull, hPv, hPe, hf, hwf, hends, hht, hI => by
    rw [List.flatMap_cons] at hPv hPe
    obtain ⟨hPv₁, hPv₂⟩ := List.forall_mem_append.1 hPv
    obtain ⟨hPe₁, hPe₂⟩ := List.forall_mem_append.1 hPe
    have hb := h.book hf (hwf t (by simp)) (hends t (by simp))
    have hg := gbTree t 0 s hb
    have hi' : ∀ v outs, t = .node v outs →
        ({ s with stackVerts := s.stackVerts.set! 0 v } : WalkState).Inv' 0 :=
      fun _ _ _ => h.inv.stackVerts_of_nil h.tstack _
    have hnv : 0 < g.nv := by
      obtain ⟨v, outs⟩ := t
      exact Nat.lt_of_le_of_lt (Nat.zero_le _) (RootState.hvlt hf v (by simp [DfsTree.verts]))
    have hstep := h.step hf (hwf t (by simp)) (hends t (by simp))
    have hfull₁ := (walk_full_aux g).1 t 0 P X s hfull (RootState.hvlt hf) (RootState.helt hf)
      (RootState.hvn hf).1 (RootState.hen hf).1 hPv₁ hPe₁ (sdTree t 0 s hi' h.shape hg hb)
    have hinv := invTree t 0 s hi' h.shape hg hb
    have hrk : wp (walkTree t 0) (fun _ s' => RootOK s') s :=
      walkTree_rootOK t s hi' h.shape hg hb h.tstack (by rw [h.sd]; exact hnv)
    have hst := (stWalk g).1 t 0 s pre [] [] P X rfl hi' h.shape hg hb hfull (RootState.hvlt hf)
      (RootState.helt hf) (RootState.hvn hf).1 (RootState.hen hf).1 hPv₁ hPe₁
      (by rw [h.sd, Nat.zero_add]; exact hht t (by simp)) (by rw [h.tstack]; rfl) (fun _ h => by cases h)
      (fun _ h => by simp [segsStack] at h) (by simpa [simBlocks, frameBlocks] using hI)
    show wp ((walkTree t 0 >>= fun _ => popTstack >>= fun top =>
      modifyItem rootItem fun it => { it with ch := it.ch ++ top.spans.2 }) >>= fun _ => walkForest rest) _ s
    rw [wp_bind, wp_bind]
    refine wp_mono _ (wp_and hstep (wp_and hfull₁ (wp_and hinv (wp_and hrk hst))))
      fun _ s₁ ⟨hstep₁, hfull₁, ⟨hi₁, hs₁⟩, hrk₁, hg₁, hsd₁, _, _, ⟨new, hts₁, hR₁⟩, hI₁⟩ => ?_
    obtain ⟨tt, htt, hside, hroot⟩ := hrk₁
    simp only [htt, segsStack, List.map_nil, List.flatten_nil, List.append_nil] at hts₁
    subst hts₁
    simp only [wp_bind, wp_popTstack, wp_modifyItem] at hstep₁ ⊢
    simp only [htt, List.head!_cons, List.tail_cons] at hstep₁ ⊢
    have hroot' : ∀ p, ¬ Items.IsParent s₁.items p rootItem := noParent_of_cnt_eq_zero hfull₁.place.root
    have hgeq₁ : s₁.g = g := hg₁.trans h.g_eq
    simp only [List.length_nil, DirsOf_zero] at hR₁ hI₁
    simp only [simBlocks, frameBlocks, List.append_nil] at hI₁
    have hI₂ := rootPop_st htt hgeq₁ hi₁ hs₁ ⟨tt, htt, hside, hroot⟩ hroot' hR₁ hI₁
    have hfull₂ := hfull₁.root_append (by simp [htt, hside]) (by
      simp only [htt, List.head!_cons]
      intro c hc
      obtain ⟨v, hv, rfl⟩ := hroot c hc
      exact ⟨v, by rwa [hgeq₁] at hv, rfl⟩)
    simp only [htt, List.head!_cons, List.tail_cons] at hfull₂
    rw [List.append_cons]
    refine stForest g rest (pre ++ [t]) _ _ X hstep₁ hfull₂ ?_ ?_ (by simpa using hf)
      (fun t' ht' => hwf t' (by simp [ht'])) (fun t' ht' => hends t' (by simp [ht']))
      (fun t' ht' => hht t' (by simp [ht'])) ?_
    · have hdv : ∀ w ∈ t.verts, w ∉ rest.flatMap DfsTree.verts := by
        have h := hf.verts_nodup
        rw [List.flatMap_append, List.flatMap_cons, List.nodup_append, List.nodup_append] at h
        exact fun w hw hw' => h.2.1.2.2 w hw w hw' rfl
      rintro v hv (hh | ⟨w, hw, hvw⟩ | ⟨e, _, hve⟩)
      · exact hPv₂ v hv hh
      · exact hdv w hw (WalkM.vertItem_inj hvw ▸ hv)
      · exact vertItem_ne_edgeItem (hf.verts_lt v (by simp [hv])) e hve
    · have hde : ∀ e ∈ t.edges, e ∉ rest.flatMap DfsTree.edges := by
        have h := hf.edges_nodup
        rw [List.flatMap_append, List.flatMap_cons, List.nodup_append, List.nodup_append] at h
        exact fun e he he' => h.2.1.2.2 e he e he' rfl
      rintro e he (hh | ⟨w, hw, hew⟩ | ⟨e', he', hee⟩)
      · exact hPe₂ e he hh
      · exact vertItem_ne_edgeItem (hf.verts_lt w (by simp [hw])) e hew.symm
      · exact hde e' he' (edgeItem_inj hee ▸ he)
    · rw [refBlocks_snoc]
      rw [← List.append_assoc]; exact hI₂.congr rfl rfl

end Spqr
