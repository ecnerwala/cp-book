import Spqr.WalkInv
import Spqr.EarShape
import Spqr.EarFrontier
import Spqr.RangesWalk
import Spqr.RangesCloseTree
import Spqr.RangesSites
import Spqr.RangesCoverTree
import Spqr.Proofs.RInvWalk
import Spqr.Proofs.RSkelWalk
import Spqr.Proofs.RSide
import Spqr.StInduct
import Spqr.StOpenBlock
import Spqr.StBdPop

/-! # The backbone invariant (PROOF.md §4.7)

One conjunction `WalkInv` at a `walkTree` entry, carrying every layer's state (ear context,
`Inv'`/`Shape`, ranges, close, root coverage, the R frames, the st-simulation hypotheses), and one
theorem `walkTree_inv` producing the conjunction `WalkInvEnd` at the exit.

`WalkInv.sites` is the cross-layer plumbing: the per-site hypotheses the component inductions were
stated against (`BookTree`, `GuardsTree`, `FrontiersTree`, `CsTree`, `RgTree`, `CoverTree`,
`CbTree`, `RSideTree`) are all *produced* from the conjunction, so no layer needs an input that
another induction carries. In this first stage the step bodies are the existing per-layer tree
theorems; `checks/WalkInvCheck.lean` evaluates every field at every site (0 violations, seeds
0..3000 × both modes). -/

namespace Spqr

theorem DfsOut.vertsList_single (o : DfsOut) : DfsOut.vertsList [o] = o.verts := by
  cases o <;> simp [DfsOut.vertsList, DfsOut.verts]

theorem DfsOut.edgesList_single (o : DfsOut) : DfsOut.edgesList [o] = o.edges := by
  cases o <;> simp [DfsOut.edgesList, DfsOut.edges]

theorem DfsOut.heightList_le_of_mem {o : DfsOut} :
    ∀ {outs : List DfsOut}, o ∈ outs → DfsOut.heightList [o] ≤ DfsOut.heightList outs
  | [], h => by simp at h
  | o' :: rest, h => by
    rcases List.mem_cons.1 h with rfl | h
    · cases o <;> simp [DfsOut.heightList]
    · have := DfsOut.heightList_le_of_mem h
      cases o' <;> simp [DfsOut.heightList] <;> omega

theorem DfsOut.verts_disjoint_of_nodup : ∀ {outs : List DfsOut}, (DfsOut.vertsList outs).Nodup →
    ∀ {o o' : DfsOut}, o ∈ outs → o' ∈ outs → o ≠ o' → ∀ {x : Nat}, x ∈ o.verts → x ∈ o'.verts → False
  | [], _, _, _, h, _, _, _, _, _ => by simp at h
  | a :: outs, hn, o, o', ho, ho', hne, x, hx, hx' => by
    rw [DfsOut.vertsList_cons, DfsOut.vertsList_single] at hn
    rcases List.mem_cons.1 ho with rfl | ho₁ <;> rcases List.mem_cons.1 ho' with rfl | ho₁'
    · exact hne rfl
    · exact List.disjoint_of_nodup_append hn hx (WalkState.mem_vertsList_of_verts ho₁' hx')
    · exact List.disjoint_of_nodup_append hn hx' (WalkState.mem_vertsList_of_verts ho₁ hx)
    · exact DfsOut.verts_disjoint_of_nodup hn.of_append_right ho₁ ho₁' hne hx hx'

end Spqr

namespace Spqr
open WalkM
/-- A `V` entry pushed on a live segment keeps it live: everything strictly below a `V` item hangs
under it. -/
theorem StLive.cons_vert {g : Graph} {items : Items} {new : List TEntry} {b : StBlock}
    (h : StLive g items new b) {v d idx : Nat} {dir : Bool} {q : ItemId}
    (hq : Items.type items q = .V) :
    StLive g items (⟨v, d, idx, setSides dir [q] []⟩ :: new) b := by
  intro x hx i hbi hty hh
  rcases (mem_readStack_push v d idx dir q new x).1 hx with rfl | hx
  · exfalso
    by_cases hiq : i = x
    · subst hiq
      rcases hty with h | h | h <;> rw [hq] at h <;> cases h
    · exact hh ⟨x, Relation.ReflTransGen.refl, Ne.symm hiq, Or.inl hq, hbi⟩
  · exact h x hx i hbi hty hh

namespace WalkState

/-- Ghost data of the backbone at a `walkTree` entry: the ear frame (`base`/`bE`/`sv`/`sd`, the
pending edges `pe` of the ancestors), the schedule `σ`/`n`, the DFS data and ancestor frames of the
R layer, the coverage prefix data `sts`/`origs`/`P`/`X`, and the st-simulation path `prev`/`fs`/`segs`. -/
structure TreeGhost where
  g : Graph
  anc : List Nat
  base : List TEntry
  bE : List (Nat → Prop)
  sv : List Nat
  sd : List Bool
  pe : Nat → Prop
  σ : List Nat
  n : Nat
  dfs : DfsData
  F : List RFrame
  sts : List Nat
  origs : List Nat
  P : ItemId → Prop
  X : ItemId → Prop
  prev : List DfsTree
  fs : List PathFrame
  segs : List (List TEntry × List StPiece)

/-- The R-layer context at a non-root entry `(t.v, d)`, `d = dp + 1` (only meaningful on a block). -/
structure RCtx (dfs : DfsData) (F : List RFrame) (t : DfsTree) (dp : Nat) (s : WalkState) : Prop where
  wf : s.g.WF
  spec : dfs.Spec s.g
  rooted : dfs.Rooted s.g
  inv_par : s.Inv' dp
  parent : ∀ v outs, t = .node v outs → dfs.IsParent s.stackVerts[dp]! v
  outs : ∀ t' : DfsTree, t'.Sub t → dfs.outs t'.v = t'.outs
  chain : ∀ k, k ≤ dp → dfs.Anc s.stackVerts[k]! s.stackVerts[dp]! ∧ dfs.depth s.stackVerts[k]! = k
  top : s.RInvTop dfs s.stackVerts[dp]! dp
  frames : ∀ f ∈ F, f.2.2 ≤ s.tstack.length ∧ s.RInvG dfs f.1 f.2.1 f.2.2
  skel : Items.RSkelInv s.g s.items

/-- The backbone invariant at the entry of `walkTree t d` from `s` (before `stackVerts[d] := t.v`). -/
structure WalkInv (G : TreeGhost) (t : DfsTree) (d : Nat) (s : WalkState) : Prop where
  -- frame / DFS facts
  full : s.Full G.g G.P G.X
  d_anc : d = G.anc.length
  wf : t.WF G.anc
  ends : t.Ends G.g
  nodup : (G.anc ++ t.verts).Nodup
  anc_lt : ∀ a ∈ G.anc, a < G.g.nv
  v_lt : ∀ v ∈ t.verts, v < G.g.nv
  e_lt : ∀ e ∈ t.edges, e < G.g.ne
  enodup : t.edges.Nodup
  comp : ∀ e', e' < G.g.ne → ∀ x, G.g.Inc e' x → x ∈ t.verts → e' ∈ t.edges ∨ G.pe e'
  pe_anc : ∀ e', G.pe e' → ∀ x, G.g.Inc e' x → x ∈ G.anc ++ [t.v]
  sv_size : s.stackVerts.size = G.g.nv
  sd_size : s.stackDir.size = G.g.nv
  anc_sv : ∀ k, k < d → G.anc[k]? = some s.stackVerts[k]!
  -- ear
  ear : EarCtx t.v d [] t.outs false G.base G.bE G.sv G.sd { s with stackVerts := s.stackVerts.set! d t.v }
  inv : ({ s with stackVerts := s.stackVerts.set! d t.v } : WalkState).Inv' d
  shape : Shape s
  -- ranges / close
  ranges : ({ s with stackVerts := s.stackVerts.set! d t.v } : WalkState).RangesInv G.σ G.n d
  σ_nodup : G.σ.Nodup
  σ_lt : ∀ e ∈ G.σ, e < s.g.ne
  post : PostAt G.σ G.n t.edgePostorder
  path : AncPath G.g G.σ (G.n + t.edgePostorder.length) (G.anc ++ [t.v])
  close : s.CloseInv
  -- root coverage
  P_past : ∀ e, e < G.g.ne → G.P (edgeItem G.g e) → G.σ.idxOf e < G.n
  owned : OwnedD G.σ (G.sts ++ [G.n]) (G.origs ++ [s.tstack.length]) G.P d G.n
    { s with stackVerts := s.stackVerts.set! d t.v }
  sts_len : G.sts.length = d
  origs_len : G.origs.length = d
  sts_le : ∀ k, k < d → G.sts[k]! ≤ G.n
  origs_le : ∀ k, k < d → G.origs[k]! ≤ s.tstack.length
  Pv : ∀ v ∈ t.verts, ¬ G.P (vertItem v)
  Pe : ∀ e ∈ t.edges, ¬ G.P (edgeItem G.g e)
  -- R (block, non-root entry)
  r : s.g.TwoConnected → ∀ dp, d = dp + 1 → RCtx G.dfs G.F t dp s
  -- st
  d_fs : d = G.fs.length
  height : d + t.height ≤ s.stackDir.size
  tstack : s.tstack = segsStack G.segs
  segRead : SegRead s.items G.segs
  base_out : ∀ t' ∈ segsStack G.segs, t'.vStart ∉ t.verts
  stItems : StItems G.g s (simBlocks G.g G.prev G.fs (DirsOf s d))
  segs_len : G.segs.length = G.fs.length
  live_lower : ∀ j (hj : j < G.segs.length),
    StLive G.g s.items G.segs[j].1
      (openBlock G.g (G.fs.take (G.fs.length - 1 - j)) (DirsOf s (G.fs.length - 1 - j)) G.segs[j].2)

/-- The backbone conclusion at the exit of `walkTree t d` started from `s`. -/
structure WalkInvEnd (G : TreeGhost) (t : DfsTree) (d : Nat) (s s' : WalkState) : Prop where
  treeEnd : TreeEnd t.v d t.outs G.base G.bE G.sv G.sd s'
  inv : s'.Inv' d
  shape : Shape s'
  ranges : RgS G.σ (G.n + t.edgePostorder.length) d s'
  close : s'.CloseInv
  place : s'.Place G.g (WalkM.Pushed G.g G.P t.verts t.edges) G.X
  owned : OwnedD G.σ (G.sts ++ [G.n]) (G.origs ++ [s.tstack.length]) (WalkM.Pushed G.g G.P t.verts t.edges) d
    (G.n + t.edgePostorder.length) s'
  sv_v : s'.stackVerts[d]! = t.v
  vertCover : VertCover t.v s'
  r : s.g.TwoConnected → ∀ dp, d = dp + 1 →
    RWalk G.dfs G.F t.v d s' ∧ BotKeep s.tstack.length s s' ∧
      (∀ k, k < d → s'.stackVerts[k]! = s.stackVerts[k]!) ∧ Items.RSkelInv s'.g s'.items
  g_eq : s'.g = s.g
  sd_size : s'.stackDir.size = s.stackDir.size
  vStart : ∀ t' ∈ s'.tstack, t'.vStart ∈ t.verts ∨ ∃ t₀ ∈ segsStack G.segs, t'.vStart = t₀.vStart
  segRead : SegRead s'.items G.segs
  st : StSim G.g G.prev G.fs t (segsStack G.segs) s'
  live : ∀ new, s'.tstack = new ++ segsStack G.segs →
    StLive G.g s'.items new (openBlock G.g G.fs (DirsOf s' d) (refTree G.g t d (DirsOf s' d)).1)
  live_lower : ∀ j (hj : j < G.segs.length),
    StLive G.g s'.items G.segs[j].1
      (openBlock G.g (G.fs.take (G.fs.length - 1 - j)) (DirsOf s' (G.fs.length - 1 - j)) G.segs[j].2)

/-- The per-site hypotheses of the component inductions, all produced by the conjunction. -/
structure Sites (G : TreeGhost) (t : DfsTree) (d : Nat) (s : WalkState) : Prop where
  ear : EarTree t d s
  book : BookTree t d s
  guards : GuardsTree t d s
  frontiers : FrontiersTree t d s
  cs : CsTree G.σ G.n t d s
  cover : CoverTree G.σ G.n t d s
  rg : RgTree G.σ G.n t d s
  cb : CbTree G.σ G.n t d s
  rside : s.g.TwoConnected → ∀ dp, d = dp + 1 → RSideTree G.dfs t d s

variable {G : TreeGhost} {t : DfsTree} {d : Nat} {s : WalkState}

namespace WalkInv

theorem types (h : WalkInv G t d s) : Types G.g s := h.full.place.types

theorem g_eq (h : WalkInv G t d s) : s.g = G.g := h.full.place.g_eq

theorem inv' (h : WalkInv G t d s) : ∀ v outs, t = .node v outs →
    ({ s with stackVerts := s.stackVerts.set! d v } : WalkState).Inv' d :=
  fun _ _ ht => by subst ht; exact h.inv

theorem ranges' (h : WalkInv G t d s) : ∀ v outs, t = .node v outs →
    ({ s with stackVerts := s.stackVerts.set! d v } : WalkState).RangesInv G.σ G.n d :=
  fun _ _ ht => by subst ht; exact h.ranges

theorem qch (h : WalkInv G t d s) : ∀ e ∈ t.edges, Items.ch s.items (edgeItem G.g e) = [] :=
  fun e he => h.full.qch e (h.e_lt e he) (h.Pe e he)

theorem earTree (h : WalkInv G t d s) :
    EarTree t d s ∧ wp (walkTree t d) (fun _ s' => TreeEnd t.v d t.outs G.base G.bE G.sv G.sd s') s :=
  cTree t d s G.g G.anc G.base G.bE G.sv G.sd G.pe h.types h.d_anc h.wf h.ends h.nodup h.anc_lt
    h.v_lt h.e_lt h.enodup h.comp h.pe_anc h.sv_size h.sd_size h.anc_sv h.ear

theorem coverTree (h : WalkInv G t d s) (hg : GuardsTree t d s) (hb : BookTree t d s)
    (hf : FrontiersTree t d s) (hc : CsTree G.σ G.n t d s) :
    CoverTree G.σ G.n t d s ∧
    wp (walkTree t d) (fun _ s' =>
      s'.Place G.g (WalkM.Pushed G.g G.P t.verts t.edges) G.X ∧
      OwnedD G.σ (G.sts ++ [G.n]) (G.origs ++ [s.tstack.length]) (WalkM.Pushed G.g G.P t.verts t.edges) d
        (G.n + t.edgePostorder.length) s' ∧
      ∀ v outs, t = .node v outs → s'.stackVerts[d]! = v ∧ VertCover v s') s :=
  cvTree G.g G.σ G.n t d G.anc G.sts G.origs G.P G.X s h.ranges' h.shape h.σ_nodup h.σ_lt hg hb hf hc
    h.full.place h.P_past (fun _ _ ht => by subst ht; exact h.owned) h.sts_len h.origs_len h.sts_le
    h.origs_le h.d_anc h.nodup h.anc_lt h.v_lt h.e_lt h.enodup h.Pv h.Pe h.post h.sv_size h.anc_sv

/-- The cross-layer plumbing: every component's site hypotheses from the conjunction. -/
theorem sites (h : WalkInv G t d s) : Sites G t d s := by
  have hE := h.earTree
  have hb : BookTree t d s :=
    bTree t d s G.g G.anc h.types h.d_anc h.wf h.ends h.nodup h.anc_lt h.v_lt h.e_lt h.enodup
      h.sv_size h.anc_sv h.qch hE.1
  have hg : GuardsTree t d s := gbTree t d s hb
  have hf : FrontiersTree t d s := frTree t d s h.inv' h.shape hg hb
  have hc : CsTree G.σ G.n t d s :=
    dsTree G.σ G.n t d s G.g G.anc G.pe h.types h.d_anc h.wf h.ends h.comp h.pe_anc h.nodup h.anc_lt
      h.v_lt h.e_lt h.enodup h.sv_size h.anc_sv h.σ_nodup h.post
      (fun _ _ ht => by subst ht; exact h.path)
  have hcv := h.coverTree hg hb hf hc
  have hr : RgTree G.σ G.n t d s :=
    scheduleTree G.σ G.n t d s h.ranges' h.shape h.σ_nodup h.σ_lt hg hb hf hcv.1 h.post
  refine ⟨hE.1, hb, hg, hf, hc, hcv.1, hr,
    cbTree G.σ G.n t d s h.ranges' h.shape h.σ_nodup h.σ_lt hg hb hf hr hc h.close,
    fun h2 dp hd => ?_⟩
  have R := h.r h2 dp hd
  exact rsTree t d s dp hd R.inv_par h.shape hg hb h2 R.wf R.spec R.rooted
    (by rw [h.g_eq]; exact h.sv_size) R.parent R.outs R.chain R.top

end WalkInv

/-! ## The between-edge level

`WalkInvOut` is the conjunction at a between-edge site of `walkOuts v d` (after `done`, before
`rest`, with the current schedule position `n` and pushed set `P`); `walkOut_inv` steps it over one
`walkOut`, `walkOuts_inv` over the remaining list. -/

/-- The R-layer context at a between-edge site (meaningful on a block, `d = dp + 1`). -/
structure ROutCtx (dfs : DfsData) (F : List RFrame) (B v d : Nat) (outs₀ : List DfsOut) (hasVert : Bool)
    (s : WalkState) : Prop where
  wf : s.g.WF
  spec : dfs.Spec s.g
  rooted : dfs.Rooted s.g
  outs_v : dfs.outs v = outs₀
  sub : ∀ e cls child, DfsOut.tree e cls child ∈ outs₀ → ∀ t' : DfsTree, t'.Sub child → dfs.outs t'.v = t'.outs
  chain : AncChain dfs v d s
  rwalk : RWalk dfs F v d s
  fB : ∀ f ∈ F, f.2.2 ≤ B
  B_le : B ≤ s.tstack.length
  hvB : hasVert = true → B + 1 ≤ s.tstack.length
  skel : Items.RSkelInv s.g s.items

/-- The backbone invariant at a between-edge site. `G.n`/`G.P` are the values at the `walkTree`
entry; `n`/`P` the current ones. -/
structure WalkInvOut (G : TreeGhost) (B v d : Nat) (outs₀ : List DfsOut) (done : List (DfsOut × Bool))
    (rest : List DfsOut) (hasVert : Bool) (n : Nat) (P : ItemId → Prop) (s : WalkState) : Prop where
  -- frame / DFS facts
  full : s.Full G.g P G.X
  d_anc : d = G.anc.length
  split : done.map (·.1) ++ rest = outs₀
  sorted : outs₀.Pairwise (fun a b => a.cls.rank ≤ b.cls.rank)
  wf : ∀ o ∈ outs₀, o.WF G.anc v
  ends : ∀ o ∈ outs₀, DfsOut.Ends G.g v o
  nodup : (G.anc ++ v :: DfsOut.vertsList outs₀).Nodup
  anc_lt : ∀ a ∈ G.anc, a < G.g.nv
  v_lt : v < G.g.nv
  w_lt : ∀ w ∈ DfsOut.vertsList outs₀, w < G.g.nv
  e_lt : ∀ e ∈ DfsOut.edgesList outs₀, e < G.g.ne
  enodup : (DfsOut.edgesList outs₀).Nodup
  comp : ∀ e', e' < G.g.ne → ∀ x, G.g.Inc e' x → x ∈ v :: DfsOut.vertsList outs₀ →
    e' ∈ DfsOut.edgesList outs₀ ∨ G.pe e'
  pe_anc : ∀ e', G.pe e' → ∀ x, G.g.Inc e' x → x ∈ G.anc ++ [v]
  sv_size : s.stackVerts.size = G.g.nv
  sd_size : s.stackDir.size = G.g.nv
  anc_sv : ∀ k, k < d → G.anc[k]? = some s.stackVerts[k]!
  sv_d : s.stackVerts[d]! = v
  -- ear
  ear : EarCtx v d done rest hasVert G.base G.bE G.sv G.sd s
  inv : s.Inv' d
  shape : Shape s
  -- ranges / close
  ranges : s.RangesInv G.σ n d
  σ_nodup : G.σ.Nodup
  σ_lt : ∀ e ∈ G.σ, e < s.g.ne
  post : PostAt G.σ n (DfsOut.edgePostorderList rest)
  path : AncPath G.g G.σ (n + (DfsOut.edgePostorderList rest).length) (G.anc ++ [v])
  close : s.CloseInv
  -- coverage
  P_past : ∀ e, e < G.g.ne → P (edgeItem G.g e) → G.σ.idxOf e < n
  owned : OwnedD G.σ G.sts G.origs P d n s
  sts_len : G.sts.length = d + 1
  origs_len : G.origs.length = d + 1
  sts_le : ∀ k, k ≤ d → G.sts[k]! ≤ G.sts[d]!
  sts_n : G.sts[d]! ≤ n
  origs_le : ∀ k, k ≤ d → G.origs[k]! ≤ G.origs[d]!
  Pv : ∀ w ∈ DfsOut.vertsList rest, ¬ P (vertItem w)
  Pe : ∀ e ∈ DfsOut.edgesList rest, ¬ P (edgeItem G.g e)
  Pcur : hasVert = false → ¬ P (vertItem v)
  Pcur' : hasVert = true → P (vertItem v)
  vcover : hasVert = true → VertCover v s
  -- R
  r : s.g.TwoConnected → ∀ dp, d = dp + 1 → ROutCtx G.dfs G.F B v d outs₀ hasVert s
  -- st
  d_fs : d = G.fs.length
  height : d + DfsOut.heightList outs₀ < s.stackDir.size
  base_out : ∀ t ∈ segsStack G.segs, t.vStart ∉ v :: DfsOut.vertsList outs₀
  stPre : StPre G.g G.prev G.fs G.segs v d (done.map (·.1)) hasVert s
  segs_len : G.segs.length = G.fs.length
  /-- The §7 pairing: the current segment is live in the open block of the current frame … -/
  live_cur : ∀ new, s.tstack = new ++ segsStack G.segs →
    StLive G.g s.items new (openBlock G.g G.fs (DirsOf s d)
      ((refOuts G.g v d (DirsOf s d) (done.map (·.1)) false).1 ++
        if hasVert then [] else [⟨true, [vertItem v]⟩]))
  /-- … and the segment of every lower frame `k` (`segs[j]`, `k = fs.length - 1 - j`) in its own. -/
  live_lower : ∀ j (hj : j < G.segs.length),
    StLive G.g s.items G.segs[j].1
      (openBlock G.g (G.fs.take (G.fs.length - 1 - j)) (DirsOf s (G.fs.length - 1 - j)) G.segs[j].2)

/-- The per-site hypotheses of the component inductions at an out-edge. -/
structure OutSites (G : TreeGhost) (v d : Nat) (o : DfsOut) (hasVert : Bool) (n : Nat) (s : WalkState) :
    Prop where
  ear : EarOut v d o hasVert s
  book : BookOut v d o hasVert s
  guards : GuardsOut v d o hasVert s
  frontiers : FrontiersOut v d o hasVert s
  cs : CsOut G.σ n v d o hasVert s
  cover : CoverOut G.σ n v d o hasVert s
  rg : RgOut G.σ n v d o hasVert s
  cb : CbOut G.σ n v d o hasVert s
  rside : s.g.TwoConnected → ∀ dp, d = dp + 1 → RSideOut G.dfs v d o hasVert s

variable {B v n : Nat} {outs₀ rest : List DfsOut} {done : List (DfsOut × Bool)} {hasVert : Bool}
  {o : DfsOut} {P : ItemId → Prop}

namespace WalkInvOut

theorem types (h : WalkInvOut G B v d outs₀ done rest hasVert n P s) : Types G.g s := h.full.place.types
theorem g_eq (h : WalkInvOut G B v d outs₀ done rest hasVert n P s) : s.g = G.g := h.full.place.g_eq
theorem rgs (h : WalkInvOut G B v d outs₀ done rest hasVert n P s) : RgS G.σ n d s :=
  ⟨h.ranges, h.shape, h.σ_lt⟩

theorem congr {n' : Nat} {P' : ItemId → Prop} (h : WalkInvOut G B v d outs₀ done rest hasVert n P s)
    (hn : n = n') (hP : ∀ i, P i ↔ P' i) : WalkInvOut G B v d outs₀ done rest hasVert n' P' s := by
  subst hn
  exact { h with
    full := h.full.monoP hP
    P_past := fun e he hp => h.P_past e he ((hP _).2 hp)
    owned := h.owned.mono fun i => (hP i).1
    Pv := fun w hw hp => h.Pv w hw ((hP _).2 hp)
    Pe := fun e he hp => h.Pe e he ((hP _).2 hp)
    Pcur := fun h0 hp => h.Pcur h0 ((hP _).2 hp)
    Pcur' := fun h1 => (hP _).1 (h.Pcur' h1) }

section head
variable (h : WalkInvOut G B v d outs₀ done (o :: rest) hasVert n P s)
include h

theorem mem : o ∈ outs₀ := by
  rw [← h.split]; exact List.mem_append_right _ List.mem_cons_self

theorem verts_eq : DfsOut.vertsList outs₀ =
    DfsOut.vertsList (done.map (·.1)) ++ (o.verts ++ DfsOut.vertsList rest) := by
  rw [← h.split, DfsOut.vertsList_append, DfsOut.vertsList_cons, DfsOut.vertsList_single]

theorem edges_eq : DfsOut.edgesList outs₀ =
    DfsOut.edgesList (done.map (·.1)) ++ (o.edges ++ DfsOut.edgesList rest) := by
  rw [← h.split, DfsOut.edgesList_append, DfsOut.edgesList_cons, DfsOut.edgesList_single]

theorem split_facts : (G.anc ++ v :: o.verts).Nodup ∧ (∀ w ∈ o.verts, w ∉ DfsOut.vertsList rest) ∧
    v ∉ o.verts ∧ v ∉ DfsOut.vertsList rest ∧ o.edges.Nodup ∧
    ∀ e ∈ o.edges, e ∉ DfsOut.edgesList rest := by
  have hn := h.nodup; rw [h.verts_eq] at hn
  have he := h.enodup; rw [h.edges_eq] at he
  have hv := (List.nodup_cons.1 hn.of_append_right)
  have hOR := hv.2.of_append_right
  refine ⟨hn.sublist (List.Sublist.append_left (List.Sublist.cons_cons v
      ((List.sublist_append_left _ _).trans (List.sublist_append_right _ _))) _),
    fun w hw hw' => List.disjoint_of_nodup_append hOR hw hw',
    fun hv' => hv.1 (List.mem_append_right _ (List.mem_append_left _ hv')),
    fun hv' => hv.1 (List.mem_append_right _ (List.mem_append_right _ hv')),
    he.of_append_right.of_append_left,
    fun e he' he'' => List.disjoint_of_nodup_append he.of_append_right he' he''⟩

theorem w_lt_o : ∀ w ∈ o.verts, w < G.g.nv :=
  fun w hw => h.w_lt w (by rw [h.verts_eq]; exact List.mem_append_right _ (List.mem_append_left _ hw))
theorem e_lt_o : ∀ e ∈ o.edges, e < G.g.ne :=
  fun e he => h.e_lt e (by rw [h.edges_eq]; exact List.mem_append_right _ (List.mem_append_left _ he))
theorem w_lt_rest : ∀ w ∈ DfsOut.vertsList rest, w < G.g.nv :=
  fun w hw => h.w_lt w (by rw [h.verts_eq]; exact List.mem_append_right _ (List.mem_append_right _ hw))
theorem e_lt_rest : ∀ e ∈ DfsOut.edgesList rest, e < G.g.ne :=
  fun e he => h.e_lt e (by rw [h.edges_eq]; exact List.mem_append_right _ (List.mem_append_right _ he))
theorem Pv_o : ∀ w ∈ o.verts, ¬ P (vertItem w) :=
  fun w hw => h.Pv w (by rw [DfsOut.vertsList_cons, DfsOut.vertsList_single]; exact List.mem_append_left _ hw)
theorem Pe_o : ∀ e ∈ o.edges, ¬ P (edgeItem G.g e) :=
  fun e he => h.Pe e (by rw [DfsOut.edgesList_cons, DfsOut.edgesList_single]; exact List.mem_append_left _ he)
theorem Pv_rest : ∀ w ∈ DfsOut.vertsList rest, ¬ P (vertItem w) :=
  fun w hw => h.Pv w (by rw [DfsOut.vertsList_cons]; exact List.mem_append_right _ hw)
theorem Pe_rest : ∀ e ∈ DfsOut.edgesList rest, ¬ P (edgeItem G.g e) :=
  fun e he => h.Pe e (by rw [DfsOut.edgesList_cons]; exact List.mem_append_right _ he)
theorem qch_o : ∀ e ∈ o.edges, Items.ch s.items (edgeItem G.g e) = [] :=
  fun e he => h.full.qch e (h.e_lt_o e he) (h.Pe_o e he)
theorem post_o : PostAt G.σ n o.block := by
  have := h.post; rw [DfsOut.edgePostorderList_cons] at this; exact this.left
theorem path_o : AncPath G.g G.σ (n + o.block.length) (G.anc ++ [v]) :=
  h.path.mono (by rw [DfsOut.edgePostorderList_cons, List.length_append]; omega)

theorem comp_out : ∀ e', e' < G.g.ne → ∀ x, G.g.Inc e' x → x ∈ o.verts → subEdges o e' ∨ G.pe e' := by
  intro e' he' x hx hxo
  have hxm : x ∈ DfsOut.vertsList outs₀ := mem_vertsList_of_verts h.mem hxo
  rcases h.comp e' he' x hx (List.mem_cons_of_mem _ hxm) with hm | hm
  · obtain ⟨o', ho', hs⟩ := mem_subEdges_edgesList.1 hm
    by_cases heq : o' = o
    · subst heq; exact .inl hs
    · exfalso
      rcases endsOut_wf G.g G.anc v o' (h.wf o' ho') (h.ends o' ho') e' hs x hx with h1 | hxv | h1
      · exact List.disjoint_of_nodup_append h.nodup h1 (List.mem_cons_of_mem _ hxm)
      · subst hxv; exact (List.nodup_cons.1 h.nodup.of_append_right).1 hxm
      · exact DfsOut.verts_disjoint_of_nodup (List.nodup_cons.1 h.nodup.of_append_right).2 ho' h.mem heq h1 hxo
  · exact .inr hm

theorem earOut : EarOut v d o hasVert s ∧ wp (walkOut v d o hasVert)
    (fun hv' s' => ∃ hvF, EarCtx v d (done ++ [(o, hvF)]) rest hv' G.base G.bE G.sv G.sd s') s :=
  cOut v d o hasVert s G.g G.anc outs₀ rest done G.base G.bE G.sv G.sd G.pe h.types h.d_anc h.split
    h.sorted h.wf h.ends h.nodup h.anc_lt h.v_lt h.w_lt h.e_lt h.enodup h.comp h.pe_anc h.sv_size
    h.sd_size h.anc_sv h.ear

theorem coverOut (hg : GuardsOut v d o hasVert s) (hb : BookOut v d o hasVert s)
    (hf : FrontiersOut v d o hasVert s) (hc : CsOut G.σ n v d o hasVert s) :
    CoverOut G.σ n v d o hasVert s ∧
    wp (walkOut v d o hasVert)
      (CvPost G.g G.σ (n + o.block.length) v d G.sts G.origs P G.X hasVert o.verts o.edges) s :=
  cvOut G.g G.σ n v d o hasVert G.anc G.sts G.origs P G.X s h.rgs h.σ_nodup hg hb hf hc h.full.place
    h.P_past h.owned h.sts_len h.origs_len h.sts_le h.sts_n h.origs_le h.d_anc h.split_facts.1
    h.anc_lt h.v_lt h.w_lt_o h.e_lt_o h.split_facts.2.2.2.2.1 h.Pv_o h.Pe_o h.Pcur h.post_o h.sv_size
    h.anc_sv h.sv_d h.vcover

theorem stHyps : StHyps G.g G.prev G.fs G.segs v d (done.map (·.1)) [o] hasVert P G.X s where
  full := h.full
  v_lt := h.v_lt
  w_lt := by rw [DfsOut.vertsList_single]; exact h.w_lt_o
  e_lt := by rw [DfsOut.edgesList_single]; exact h.e_lt_o
  vnodup := by
    have hn := h.nodup.of_append_right
    rw [h.verts_eq] at hn
    rw [DfsOut.vertsList_append, DfsOut.vertsList_single]
    exact hn.sublist (List.cons_sublist_cons.2 ((List.sublist_append_left _ _).append_left _))
  enodup := by rw [DfsOut.edgesList_single]; exact h.split_facts.2.2.2.2.1
  Pv := by rw [DfsOut.vertsList_single]; exact h.Pv_o
  Pe := by rw [DfsOut.edgesList_single]; exact h.Pe_o
  Pcur := h.Pcur
  Pcur' := h.Pcur'
  height := Nat.lt_of_le_of_lt (Nat.add_le_add_left (DfsOut.heightList_le_of_mem h.mem) d) h.height
  base := fun t ht hm => h.base_out t ht (by
    rw [h.verts_eq]
    rw [DfsOut.vertsList_append, DfsOut.vertsList_single] at hm
    rcases List.mem_cons.1 hm with rfl | hm
    · exact List.mem_cons_self
    · exact List.mem_cons_of_mem _ (by
        rcases List.mem_append.1 hm with hm | hm
        · exact List.mem_append_left _ hm
        · exact List.mem_append_right _ (List.mem_append_left _ hm)))
  pre := h.stPre

/-- The cross-layer plumbing at an out-edge. -/
theorem sites : OutSites G v d o hasVert n s := by
  obtain ⟨hndO, -, -, -, hendO, -⟩ := h.split_facts
  have hE := h.earOut
  have hb : BookOut v d o hasVert s :=
    bOut v d o hasVert s G.g G.anc h.types h.d_anc (h.wf o h.mem) (h.ends o h.mem) hndO h.anc_lt h.v_lt
      h.w_lt_o h.e_lt_o hendO h.sv_size h.anc_sv h.sv_d h.qch_o hE.1
  have hg := gbOut v d o hasVert s hb
  have hf : FrontiersOut v d o hasVert s := frOut v d o hasVert s h.inv h.shape hg hb
  have hc : CsOut G.σ n v d o hasVert s :=
    dsOut G.σ n v d o hasVert s G.g G.anc G.pe h.types h.d_anc (h.wf o h.mem) (h.ends o h.mem) h.comp_out
      h.pe_anc hndO h.anc_lt h.v_lt h.w_lt_o h.e_lt_o hendO h.sv_size h.anc_sv h.sv_d h.σ_nodup h.post_o
      h.path_o
  have hcv := h.coverOut hg hb hf hc
  have hr : RgOut G.σ n v d o hasVert s :=
    scheduleOut G.σ n v d o hasVert s h.rgs h.σ_nodup hg hb hf hcv.1 h.post_o
  refine ⟨hE.1, hb, hg, hf, hc, hcv.1, hr,
    cbOut G.σ n v d o hasVert s h.rgs h.σ_nodup hg hb hf hr hc h.close, fun h2 dp hd => ?_⟩
  have R := h.r h2 dp hd
  exact rsOut v d o hasVert s G.F B dp hd h.inv h.shape hg hb h2 R.wf R.spec R.rooted
    (by rw [h.g_eq]; exact h.sv_size) (by rw [R.outs_v]; exact h.mem)
    (fun e cls child ho => R.sub e cls child (ho ▸ h.mem)) R.chain R.rwalk R.fB R.B_le R.hvB

end head
end WalkInvOut

/-- **Named admission** (PROOF.md §4.7, "StLive step"; checker fields `st.live_cur`/`st.live_lower`,
0 violations): one `walkOut` keeps the per-segment `StLive ↔ openBlock` pairing — the current
segment in the open block of the current frame (now with `o` done), every lower segment in its own. -/
theorem walkOut_stLive (h : WalkInvOut G B v d outs₀ done (o :: rest) hasVert n P s) :
    wp (walkOut v d o hasVert) (fun hv' s' =>
      (∀ new, s'.tstack = new ++ segsStack G.segs →
        StLive G.g s'.items new (openBlock G.g G.fs (DirsOf s' d)
          ((refOuts G.g v d (DirsOf s' d) (done.map (·.1) ++ [o]) false).1 ++
            if hv' then [] else [⟨true, [vertItem v]⟩]))) ∧
      ∀ j (hj : j < G.segs.length),
        StLive G.g s'.items G.segs[j].1
          (openBlock G.g (G.fs.take (G.fs.length - 1 - j)) (DirsOf s' (G.fs.length - 1 - j)) G.segs[j].2)) s := by
  sorry

theorem stOutIh (g : Graph) (o : DfsOut) (d : Nat) :
    match o with | .tree _ _ child => StTreeP g child (d + 1) | .back .. => True := by
  cases o with
  | tree e cls child => exact (stWalk g).1 child (d + 1)
  | back e dest cls => exact True.intro

/-- The R part of one `walkOut` step. -/
theorem walkOut_inv_r (h : WalkInvOut G B v d outs₀ done (o :: rest) hasVert n P s) (S : OutSites G v d o hasVert n s) :
    wp (walkOut v d o hasVert) (fun hv' s' => s.g.TwoConnected → ∀ dp, d = dp + 1 →
      ROutCtx G.dfs G.F B v d outs₀ hv' s') s := by
  by_cases h2 : s.g.TwoConnected
  · rcases d with _ | dp
    · exact wp_of_forall fun _ _ _ dp hd => by omega
    · have R := h.r h2 dp rfl
      have hside := S.rside h2 dp rfl
      have hW := rrOut v (dp + 1) o hasVert s G.F B h.inv h.shape S.guards S.book hside h2 R.spec R.rooted
        R.chain R.rwalk R.fB R.B_le R.hvB
      have hK := rkOut v (dp + 1) o hasVert s G.F B h.inv h.shape S.guards S.book hside h2 R.spec R.rooted
        R.chain R.rwalk R.fB R.B_le R.hvB R.skel
      refine wp_mono _ (wp_and hW hK) fun hv' s' ⟨⟨hw, hbk, hhvB, hg, hsv⟩, hskel⟩ _ dp' hd => ?_
      cases hd
      exact {
        wf := by rw [hg]; exact R.wf
        spec := by rw [hg]; exact R.spec
        rooted := by rw [hg]; exact R.rooted
        outs_v := R.outs_v
        sub := R.sub
        chain := ⟨by rw [hsv (dp + 1) (Nat.le_refl _)]; exact R.chain.1,
          fun k hk => by rw [hsv k hk]; exact R.chain.2 k hk⟩
        rwalk := hw
        fB := R.fB
        B_le := hbk.1
        hvB := hhvB
        skel := by rw [hg]; exact hskel }
  · exact wp_of_forall fun _ _ h2' => absurd h2' h2

/-- **One `walkOut` step of the backbone.** -/
theorem walkOut_inv (h : WalkInvOut G B v d outs₀ done (o :: rest) hasVert n P s) :
    wp (walkOut v d o hasVert) (fun hv' s' => (hasVert = true → hv' = true) ∧ ∃ hvF,
      WalkInvOut G B v d outs₀ (done ++ [(o, hvF)]) rest hv' (n + o.block.length)
        (Pushed G.g (fun i => P i ∨ (hv' = true ∧ i = vertItem v)) o.verts o.edges) s') s := by
  obtain ⟨hndO, hdisjV, hvo, hvr, hendO, hdisjE⟩ := h.split_facts
  have S := h.sites
  have hE := h.earOut
  have hK := kOut v d o hasVert s G.g (d + 1) 0 s h.types (Nat.le_refl _) (by omega) (vertItem_ne_zero v)
    (fun w _ => vertItem_ne_zero w) (fun e _ => edgeItem_ne_zero G.g e) Keep.refl
  have hI := invOut v d o hasVert s h.inv h.shape S.guards S.book
  have hRg := rgOut G.σ n v d o hasVert s h.rgs h.σ_nodup S.guards S.book S.rg
  have hC := ccOut G.σ n v d o hasVert s h.rgs h.σ_nodup S.guards S.book S.frontiers S.rg S.cs h.close
  have hV := (h.coverOut S.guards S.book S.frontiers S.cs).2
  have hF := (walk_full_aux G.g).2.2 v d o hasVert P G.X s h.full h.v_lt
    (by rw [DfsOut.vertsList_single]; exact h.w_lt_o) (by rw [DfsOut.edgesList_single]; exact h.e_lt_o)
    (by rw [DfsOut.vertsList_single]; exact (List.nodup_cons.1 hndO.of_append_right).2)
    (by rw [DfsOut.edgesList_single]; exact hendO) (by rw [DfsOut.vertsList_single]; exact hvo)
    (by rw [DfsOut.vertsList_single]; exact h.Pv_o) (by rw [DfsOut.edgesList_single]; exact h.Pe_o)
    h.Pcur h.Pcur' (sdOut v d o hasVert s h.inv h.shape S.guards S.book)
  have hS := stOut_step G.g v d o hasVert (stOutIh G.g o d) s G.prev G.fs G.segs (done.map (·.1)) P G.X h.d_fs h.inv
    h.shape S.guards S.book h.stHyps
  have hRR := walkOut_inv_r h S
  have hL := walkOut_stLive h
  refine wp_mono _ (wp_and hE.2 (wp_and hK (wp_and hI (wp_and hRg (wp_and hC (wp_and hV (wp_and hF
    (wp_and hS (wp_and hRR hL))))))))) fun hv' s' ⟨⟨hvF, hC'⟩, hK', ⟨hi', hs'⟩, ⟨hrg', _, hσ'⟩, hcl',
      ⟨hpl', how', hsv', hvc', hhv'⟩, ⟨hfull', _⟩, ⟨hg', hsd', hp'⟩, hrr', hlc, hll⟩ => ?_
  refine ⟨hhv', hvF, ?_⟩
  have hPiff : ∀ i, Pushed G.g (fun i => P i ∨ (hv' = true ∧ i = vertItem v))
      (DfsOut.vertsList [o]) (DfsOut.edgesList [o]) i ↔
      Pushed G.g (fun i => P i ∨ (hv' = true ∧ i = vertItem v)) o.verts o.edges i := by
    intro i; rw [DfsOut.vertsList_single, DfsOut.edgesList_single]
  have hhf : hv' = false → hasVert = false := fun h0 => by
    cases hasVert
    · rfl
    · exact absurd (hhv' rfl) (by rw [h0]; decide)
  exact {
    full := hfull'.monoP hPiff
    d_anc := h.d_anc
    split := by rw [List.map_append, List.append_assoc]; exact h.split
    sorted := h.sorted
    wf := h.wf
    ends := h.ends
    nodup := h.nodup
    anc_lt := h.anc_lt
    v_lt := h.v_lt
    w_lt := h.w_lt
    e_lt := h.e_lt
    enodup := h.enodup
    comp := h.comp
    pe_anc := h.pe_anc
    sv_size := hK'.sv.trans h.sv_size
    sd_size := hK'.sd.trans h.sd_size
    anc_sv := fun k hk => by rw [hK'.svlo k (by omega)]; exact h.anc_sv k hk
    sv_d := hsv'
    ear := hC'
    inv := hi'
    shape := hs'
    ranges := hrg'
    σ_nodup := h.σ_nodup
    σ_lt := hσ'
    post := by have := h.post; rw [DfsOut.edgePostorderList_cons] at this; exact this.right
    path := by
      have := h.path; rw [DfsOut.edgePostorderList_cons, List.length_append, ← Nat.add_assoc] at this
      exact this
    close := hcl'
    P_past := pushed_past (fun e he hp => by
        rcases hp with hp | ⟨_, hp⟩
        · exact h.P_past e he hp
        · exact absurd hp.symm (vertItem_ne_edgeItem h.v_lt e))
      h.w_lt_o (DfsOut.edges_subset_block o) h.post_o h.σ_nodup
    owned := how'
    sts_len := h.sts_len
    origs_len := h.origs_len
    sts_le := h.sts_le
    sts_n := Nat.le_trans h.sts_n (Nat.le_add_right _ _)
    origs_le := h.origs_le
    Pv := by
      rintro w hw ((hp | ⟨_, hp⟩) | ⟨w', hw', hp⟩ | ⟨e, _, hp⟩)
      · exact h.Pv_rest w hw hp
      · exact hvr (vertItem_inj hp ▸ hw)
      · exact hdisjV w' hw' (vertItem_inj hp ▸ hw)
      · exact vertItem_ne_edgeItem (h.w_lt_rest w hw) e hp
    Pe := by
      rintro e he ((hp | ⟨_, hp⟩) | ⟨w', hw', hp⟩ | ⟨e', he', hp⟩)
      · exact h.Pe_rest e he hp
      · exact vertItem_ne_edgeItem h.v_lt e hp.symm
      · exact vertItem_ne_edgeItem (h.w_lt_o w' hw') e hp.symm
      · exact hdisjE e' he' (edgeItem_inj hp ▸ he)
    Pcur := by
      rintro h0 ((hp | ⟨h1, _⟩) | ⟨w', hw', hp⟩ | ⟨e, _, hp⟩)
      · exact h.Pcur (hhf h0) hp
      · rw [h0] at h1; cases h1
      · exact hvo (vertItem_inj hp ▸ hw')
      · exact vertItem_ne_edgeItem h.v_lt e hp
    Pcur' := fun h1 => Or.inl (Or.inr ⟨h1, rfl⟩)
    vcover := hvc'
    r := fun h2 dp hd => hrr' (by rw [hK'.g] at h2; exact h2) dp hd
    d_fs := h.d_fs
    height := by rw [hsd']; exact h.height
    base_out := h.base_out
    stPre := by rw [List.map_append]; exact hp'
    segs_len := h.segs_len
    live_cur := by rw [List.map_append]; exact hlc
    live_lower := hll }

/-- **The remaining out-edges**: the backbone through `walkOuts v d rest hasVert`. -/
theorem walkOuts_inv : ∀ (rest : List DfsOut) (done : List (DfsOut × Bool)) (hasVert : Bool) (n : Nat)
    (P : ItemId → Prop) (s : WalkState), WalkInvOut G B v d outs₀ done rest hasVert n P s →
    wp (walkOuts v d rest hasVert) (fun hv' s' => (hasVert = true → hv' = true) ∧ ∃ done',
      WalkInvOut G B v d outs₀ done' [] hv' (n + (DfsOut.edgePostorderList rest).length)
        (Pushed G.g (fun i => P i ∨ (hv' = true ∧ i = vertItem v))
          (DfsOut.vertsList rest) (DfsOut.edgesList rest)) s') s
  | [], done, hasVert, n, P, s, h => by
    unfold walkOuts; simp only [wp_pure]
    refine ⟨fun h => h, done, h.congr (by simp [DfsOut.edgePostorderList]) fun i => ?_⟩
    constructor
    · exact fun hp => Or.inl (Or.inl hp)
    · rintro ((hp | ⟨h1, rfl⟩) | ⟨w, hw, _⟩ | ⟨e, he, _⟩)
      · exact hp
      · exact h.Pcur' h1
      · simp [DfsOut.vertsList] at hw
      · simp [DfsOut.edgesList] at he
  | o :: rest, done, hasVert, n, P, s, h => by
    unfold walkOuts; simp only [wp_bind]
    refine wp_mono _ (walkOut_inv h) fun hv₁ s₁ ⟨hhv₁, hvF, h₁⟩ => ?_
    refine wp_mono _ (walkOuts_inv rest _ hv₁ _ _ s₁ h₁) fun hv' s' ⟨hhv', done', h'⟩ =>
      ⟨fun h0 => hhv' (hhv₁ h0), done', h'.congr ?_ fun i => ?_⟩
    · rw [DfsOut.edgePostorderList_cons, List.length_append, Nat.add_assoc]
    · rw [DfsOut.vertsList_cons, DfsOut.vertsList_single, DfsOut.edgesList_cons, DfsOut.edgesList_single]
      constructor
      · exact Pushed.append hhv' i
      · rintro ((hp | hp) | ⟨w, hw, rfl⟩ | ⟨e, he, rfl⟩)
        · exact Or.inl (Or.inl (Or.inl (Or.inl hp)))
        · exact Or.inl (Or.inr hp)
        · rcases List.mem_append.1 hw with hw | hw
          · exact Or.inl (Or.inl (Or.inr (Or.inl ⟨w, hw, rfl⟩)))
          · exact Or.inr (Or.inl ⟨w, hw, rfl⟩)
        · rcases List.mem_append.1 he with he | he
          · exact Or.inl (Or.inl (Or.inr (Or.inr ⟨e, he, rfl⟩)))
          · exact Or.inr (Or.inr ⟨e, he, rfl⟩)


/-! ## Tree-level glue

`WalkInv.toOut` carries the `walkTree` conjunction to the first between-edge site;
`WalkInvOut.exit_true`/`exit_false` carry the last one through the end of `walkTree`
(no vertex entry pushed / `setStackDir d true; pushVertTstack v d`); `walkTree_inv` is the
composition through `walkOuts_inv`. -/

namespace WalkInv

theorem d_lt (h : WalkInv G t d s) : d < G.g.nv := by
  have h1 := h.height; have h2 := h.sd_size
  cases t; simp only [DfsTree.height] at h1; omega

/-- Entry glue: the `walkTree` conjunction at `(.node v outs, d)` gives the between-edge
conjunction before the first out-edge, after `stackVerts.set! d v`. -/
theorem toOut {v : Nat} {outs : List DfsOut} (h : WalkInv G (.node v outs) d s) :
    WalkInvOut { G with sts := G.sts ++ [G.n], origs := G.origs ++ [s.tstack.length] }
      s.tstack.length v d outs [] outs false G.n G.P
      { s with stackVerts := s.stackVerts.set! d v } := by
  have hdlt : d < G.g.nv := h.d_lt
  have hstd : (G.sts ++ [G.n])[d]! = G.n := by
    rw [← h.sts_len]; exact getElem!_concat_length' _ _
  have hord : (G.origs ++ [s.tstack.length])[d]! = s.tstack.length := by
    rw [← h.origs_len]; exact getElem!_concat_length' _ _
  have hwf := h.wf
  simp only [DfsTree.WF] at hwf
  have hends := h.ends
  simp only [DfsTree.Ends] at hends
  have hnd := h.nodup; have hvlt := h.v_lt; have helt := h.e_lt; have hen := h.enodup
  have hcomp := h.comp; have hpe := h.pe_anc; have hpost := h.post; have hpath := h.path
  have hPv := h.Pv; have hPe := h.Pe; have hbase := h.base_out
  simp only [DfsTree.verts, DfsTree.edges] at hnd hvlt helt hen hcomp hpe hPv hPe hbase
  simp only [DfsTree.edgePostorder] at hpost hpath
  have hv : v < G.g.nv := hvlt v (List.mem_cons_self ..)
  have hnew : ∀ new, ({ s with stackVerts := s.stackVerts.set! d v } : WalkState).tstack =
      new ++ segsStack G.segs → new = [] := fun new hts =>
    List.append_cancel_right (hts.symm.trans (by rw [List.nil_append]; exact h.tstack))
  exact {
    full := h.full.of_eq rfl rfl rfl
    d_anc := h.d_anc
    split := by simp
    sorted := hwf.1
    wf := hwf.2
    ends := hends
    nodup := hnd
    anc_lt := h.anc_lt
    v_lt := hv
    w_lt := fun w hw => hvlt w (List.mem_cons_of_mem _ hw)
    e_lt := helt
    enodup := hen
    comp := hcomp
    pe_anc := hpe
    sv_size := by simp [h.sv_size]
    sd_size := h.sd_size
    anc_sv := fun k hk => by
      show G.anc[k]? = some (s.stackVerts.set! d v)[k]!
      rw [getElem!_set!_ne' _ _ _ _ (Nat.ne_of_lt hk)]; exact h.anc_sv k hk
    sv_d := getElem!_set!_self' _ _ _ (by rw [h.sv_size]; exact hdlt)
    ear := h.ear
    inv := h.inv
    shape := h.shape.frame'
    ranges := h.ranges
    σ_nodup := h.σ_nodup
    σ_lt := h.σ_lt
    post := hpost
    path := hpath
    close := h.close.frame rfl rfl (fun _ h => h)
    P_past := h.P_past
    owned := h.owned
    sts_len := by simp [h.sts_len]
    origs_len := by simp [h.origs_len]
    sts_le := fun k hk => by
      rw [hstd]
      rcases Nat.lt_or_eq_of_le hk with hk | rfl
      · rw [getElem!_append_left' (h.sts_len ▸ hk)]; exact h.sts_le k hk
      · rw [hstd]
    sts_n := le_of_eq hstd
    origs_le := fun k hk => by
      rw [hord]
      rcases Nat.lt_or_eq_of_le hk with hk | rfl
      · rw [getElem!_append_left' (h.origs_len ▸ hk)]; exact h.origs_le k hk
      · rw [hord]
    Pv := fun w hw => hPv w (List.mem_cons_of_mem _ hw)
    Pe := hPe
    Pcur := fun _ => hPv v (List.mem_cons_self ..)
    Pcur' := fun h => by cases h
    vcover := fun h => by cases h
    r := fun h2 dp hd => by
      have R := h.r h2 dp hd
      have hr := h.sites.rside h2 dp hd
      subst hd
      unfold RSideTree at hr
      obtain ⟨hanc, hstab, hvs, -⟩ := hr
      simp only [Nat.add_sub_cancel] at hvs
      exact {
        wf := R.wf, spec := R.spec, rooted := R.rooted
        outs_v := R.outs _ (DfsTree.Sub.refl _)
        sub := fun e cls child hm t' hs => R.outs t' (hs.step hm)
        chain := hanc
        rwalk := ⟨fun f hf => ⟨(R.frames f hf).1, ⟨fun t ht hd hne => hstab t (List.mem_of_mem_drop ht)
              ((R.frames f hf).2.entries t ht hd hne), (R.frames f hf).2.disj⟩⟩,
          ⟨fun t ht hd hne => hstab t ht (R.top.entries t ht (by omega) fun h => by
            have := hvs t ht h; omega), R.top.disj⟩⟩
        fB := fun f hf => (R.frames f hf).1
        B_le := le_refl _
        hvB := fun h => nomatch h
        skel := R.skel }
    d_fs := h.d_fs
    height := by
      have := h.height; simp only [DfsTree.height] at this
      show d + _ < s.stackDir.size; omega
    base_out := hbase
    stPre := {
      read := ⟨[], by simpa using h.tstack, by rw [List.map_nil, StRefEt.refOuts_nil]; exact StRead.nil, fun _ => rfl⟩
      hv := by rw [List.map_nil, StRefEt.refOuts_nil]
      segRead := h.segRead
      vStart := fun t ht => Or.inr ⟨t, by rw [← h.tstack]; exact ht, rfl⟩
      items := by rw [List.map_nil, StRefEt.refOuts_nil, List.append_nil]; exact h.stItems.congr rfl rfl }
    segs_len := h.segs_len
    live_cur := fun new hts => by
      rw [hnew new hts]; intro x hx; simp [readStack, readL, readR] at hx
    live_lower := h.live_lower }

end WalkInv

namespace WalkInvOut

variable {outs : List DfsOut} {done' : List (DfsOut × Bool)} {s₂ : WalkState}

/-- Exit glue, `hasVert = true`: the end-of-`walkOuts` conjunction is the `walkTree` conclusion. -/
theorem exit_true (hW : WalkInv G (.node v outs) d s)
    (h : WalkInvOut { G with sts := G.sts ++ [G.n], origs := G.origs ++ [s.tstack.length] }
      s.tstack.length v d outs done' [] true (G.n + (DfsOut.edgePostorderList outs).length)
      (Pushed G.g (fun i => G.P i ∨ (true = true ∧ i = vertItem v))
        (DfsOut.vertsList outs) (DfsOut.edgesList outs)) s₂) :
    WalkInvEnd G (.node v outs) d s s₂ := by
  have hmap : done'.map (·.1) = outs := by simpa using h.split
  obtain ⟨new, hts₂, hR₂, -⟩ := h.stPre.read
  have hhv := h.stPre.hv
  have hI := h.stPre.items
  have hvs := h.stPre.vStart
  have hlc := h.live_cur new hts₂
  rw [hmap] at hhv hR₂ hI hvs hlc
  exact {
    treeEnd := ⟨done', true, s₂, false, [], false, s₂.stackDir, rfl, h.ear, hmap, by simp, rfl,
      fun _ _ => rfl, fun h => absurd h Bool.false_ne_true⟩
    inv := h.inv
    shape := h.shape
    ranges := h.rgs
    close := h.close
    place := h.full.place.mono (fun i hi => Pushed.vert i hi) (fun _ h => h)
    owned := h.owned.mono fun i hi => Pushed.vert i hi
    sv_v := h.sv_d
    vertCover := h.vcover rfl
    r := fun h2 dp hd => by
      have R := h.r (by rw [h.g_eq, ← hW.g_eq]; exact h2) dp hd
      exact ⟨R.rwalk,
        botKeep_of_base (A := []) (A' := new) (by simpa using hW.tstack) hts₂
          (le_of_eq (congrArg List.length hW.tstack)),
        fun k hk => Option.some.inj ((h.anc_sv k hk).symm.trans (hW.anc_sv k hk)), R.skel⟩
    g_eq := h.g_eq.trans hW.g_eq.symm
    sd_size := h.sd_size.trans hW.sd_size.symm
    vStart := hvs
    segRead := h.stPre.segRead
    st := ⟨⟨new, hts₂, by rw [StRefEt.refTree_node, ← h.d_fs, hhv]; exact hR₂⟩,
      by rw [StRefEt.refTree_node, ← h.d_fs]; exact hI⟩
    live := fun new' hts' => by
      obtain rfl : new' = new := List.append_cancel_right (hts'.symm.trans hts₂)
      rw [StRefEt.refTree_node, hhv]; simpa using hlc
    live_lower := h.live_lower }

/-- Exit glue, `hasVert = false`: the end-of-`walkOuts` conjunction gives the `walkTree`
conclusion after `setStackDir d true; pushVertTstack v d`. -/
theorem exit_false (hW : WalkInv G (.node v outs) d s)
    (h : WalkInvOut { G with sts := G.sts ++ [G.n], origs := G.origs ++ [s.tstack.length] }
      s.tstack.length v d outs done' [] false (G.n + (DfsOut.edgePostorderList outs).length)
      (Pushed G.g (fun i => G.P i ∨ (false = true ∧ i = vertItem v))
        (DfsOut.vertsList outs) (DfsOut.edgesList outs)) s₂) :
    WalkInvEnd G (.node v outs) d s
      { s₂ with
        stackDir := s₂.stackDir.set! d true,
        tstack := ⟨v, d, s₂.nxtEdgeIdx,
          setSides (getElem! (s₂.stackDir.set! d true) d) [vertItem v] []⟩ :: s₂.tstack } := by
  have hmap : done'.map (·.1) = outs := by simpa using h.split
  obtain ⟨new, hts₂, hR₂, -⟩ := h.stPre.read
  have hhv := h.stPre.hv
  have hI := h.stPre.items
  have hvs := h.stPre.vStart
  have hlc := h.live_cur new hts₂
  rw [hmap] at hhv hR₂ hI hvs hlc
  have hv : v < s₂.g.nv := by rw [h.g_eq]; exact h.v_lt
  have hdlt : d < s₂.stackDir.size := by have := h.height; omega
  have hDd : (s₂.stackDir.set! d true)[d]! = true := Array.getElem!_set!_self _ _ _ hdlt
  obtain ⟨hc, ha⟩ := h.ear.vert_book rfl
  have hi₂ : ({ s₂ with stackDir := s₂.stackDir.set! d true } : WalkState).Inv' d := h.inv.frame'
  have hs₂ : Shape { s₂ with stackDir := s₂.stackDir.set! d true } := h.shape.frame'
  have hnP : ¬ Pushed G.g (fun i => G.P i ∨ (false = true ∧ i = vertItem v))
      (DfsOut.vertsList outs) (DfsOut.edgesList outs) (vertItem v) := h.Pcur rfl
  have hv0 : 0 < vertItem v := by show 0 < 1 + v; omega
  have hvG : v < G.g.nv := h.v_lt
  have hvlt : vertItem v < 1 + G.g.nv + G.g.ne := by show 1 + v < _; omega
  have hpush : ∀ i, (Pushed G.g (fun i => G.P i ∨ (false = true ∧ i = vertItem v))
      (DfsOut.vertsList outs) (DfsOut.edgesList outs) i ∨ i = vertItem v) →
      Pushed G.g G.P (v :: DfsOut.vertsList outs) (DfsOut.edgesList outs) i := fun i hi => by
    rcases hi with hi | rfl
    · exact Pushed.vert i hi
    · exact Or.inr (Or.inl ⟨v, List.mem_cons_self .., rfl⟩)
  have hSt := Step.pushVert hi₂ hs₂ d hv hc ha
  have hRg := RgStep.pushVert (h.ranges.frame' (s' := { s₂ with stackDir := s₂.stackDir.set! d true }))
    hs₂ d hv hc ha ((h.full.place.set_stackDir _).pushVertR h.v_lt h.P_past)
  obtain ⟨hroot, hoff⟩ := Place.fresh h.full.place hv0 hvlt hnP
  have hty : Items.type s₂.items (vertItem v) = .V := h.shape.vert v hv
  have hsz : vertItem v < s₂.items.size := by
    have := h.shape.size; show 1 + v < _; omega
  have hdirs := DirsOf_push s₂ d true
    ⟨v, d, s₂.nxtEdgeIdx, setSides (getElem! (s₂.stackDir.set! d true) d) [vertItem v] []⟩
  have hd : d = G.fs.length := h.d_fs
  exact {
    treeEnd := ⟨done', false, s₂, (s₂.stackDir.set! d true)[d]!, _, true, s₂.stackDir.set! d true, rfl,
      h.ear, hmap, by simp, rfl, fun k hk => getElem!_set!_ne' s₂.stackDir d k true (Nat.ne_of_lt hk),
      fun _ hd => getElem!_set!_self' _ _ _ (by simpa using hd)⟩
    inv := hSt.inv
    shape := hSt.shape
    ranges := ⟨hRg.ranges, hSt.shape, h.σ_lt⟩
    close := (h.close.frame (s' := { s₂ with stackDir := s₂.stackDir.set! d true }) rfl rfl
      (fun _ h => h)).pushVert v d hv (by have := h.shape.size; show s₂.g.nv < s₂.items.size; omega)
    place := ((h.full.place.set_stackDir (s₂.stackDir.set! d true)).cons_fixed hv0 hvlt hnP v d
      s₂.nxtEdgeIdx (s₂.stackDir.set! d true)[d]!).mono hpush (fun _ h => h)
    owned := (h.owned.pushVert (s' := { s₂ with
          stackDir := s₂.stackDir.set! d true,
          tstack := ⟨v, d, s₂.nxtEdgeIdx, setSides (getElem! (s₂.stackDir.set! d true) d) [vertItem v] []⟩ :: s₂.tstack })
      h.sv_d h.sts_le rfl rfl rfl rfl).mono fun i hi => Pushed.vert i hi
    sv_v := h.sv_d
    vertCover := vertCover_pushVert rfl
    r := fun h2 dp hd => by
      have h2' : s₂.g.TwoConnected := by rw [h.g_eq, ← hW.g_eq]; exact h2
      have R := h.r h2' dp hd
      have hfree₀ : VertFree v s₂ :=
        rSide_vertFree_site h.inv h.shape (fun _ => ⟨hv, hc, ha⟩) h2' R.spec R.rooted R.chain R.rwalk
      have hfree : VertFree v { s₂ with stackDir := s₂.stackDir.set! d true } := hfree₀
      exact ⟨(RWalk.of_eq (s := s₂) (s' := { s₂ with stackDir := s₂.stackDir.set! d true }) rfl rfl rfl rfl
          R.rwalk).pushVert (k := d) hfree,
        (botKeep_of_base (A := []) (A' := new) (by simpa using hW.tstack) hts₂
          (le_of_eq (congrArg List.length hW.tstack))).trans
          (botKeep_cons (s := { s₂ with stackDir := s₂.stackDir.set! d true }) rfl R.B_le),
        fun k hk => Option.some.inj ((h.anc_sv k hk).symm.trans (hW.anc_sv k hk)), R.skel⟩
    g_eq := h.g_eq.trans hW.g_eq.symm
    sd_size := by simp; exact h.sd_size.trans hW.sd_size.symm
    vStart := by
      intro t ht
      rcases List.mem_cons.1 ht with rfl | ht
      · exact Or.inl (List.mem_cons_self ..)
      · exact hvs t ht
    segRead := h.stPre.segRead
    st := ⟨⟨⟨v, d, s₂.nxtEdgeIdx, setSides (getElem! (s₂.stackDir.set! d true) d) [vertItem v] []⟩ :: new,
        by show _ :: s₂.tstack = _; rw [hts₂, List.cons_append],
        by
          rw [StRefEt.refTree_node, ← h.d_fs, hdirs, hhv]
          simp only [Bool.false_eq_true, ↓reduceIte]
          rw [hDd]
          exact StRead.pushEntry v d s₂.nxtEdgeIdx true (vertItem v) (Or.inl hty) hR₂⟩,
      by
        rw [StRefEt.refTree_node, ← h.d_fs, hdirs]
        exact (StItems.pushEntry v d s₂.nxtEdgeIdx _ (vertItem v) hsz hroot hoff hI).congr rfl rfl⟩
    live := fun new' hts' => by
      obtain rfl : new' = _ :: new :=
        List.append_cancel_right (hts'.symm.trans (by show _ :: s₂.tstack = _; rw [hts₂, List.cons_append]))
      rw [StRefEt.refTree_node, hdirs, hhv]
      simp only [Bool.false_eq_true, ↓reduceIte]
      exact StLive.cons_vert (by simpa using hlc) hty
    live_lower := fun j hj => by
      have := h.live_lower j hj
      rwa [DirsOf_congr fun k hk => Array.getElem!_set!_ne _ _ _ _ (by omega)] }

end WalkInvOut

/-- **The backbone theorem**: the conjunction at the entry gives the conjunction at the exit,
through `walkOuts_inv` and the exit glue. -/
theorem walkTree_inv (h : WalkInv G t d s) :
    wp (walkTree t d) (fun _ s' => WalkInvEnd G t d s s') s := by
  obtain ⟨v, outs⟩ := t
  unfold walkTree
  simp only [wp_bind, wp_modify]
  refine wp_mono _ (walkOuts_inv outs [] false G.n G.P _ h.toOut) fun hv' s₂ ⟨_, done', h₂⟩ => ?_
  cases hv'
  · simp only [Bool.false_eq_true, ↓reduceIte, wp_bind, wp_setStackDir, wp_pushVertTstack]
    exact h₂.exit_false h
  · simp only [↓reduceIte, wp_pure]
    exact h₂.exit_true h

end WalkState
end Spqr
