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
open WalkM
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

/-- **The backbone theorem**: the conjunction at the entry gives the conjunction at the exit. -/
theorem walkTree_inv (h : WalkInv G t d s) :
    wp (walkTree t d) (fun _ s' => WalkInvEnd G t d s s') s := by
  have S := h.sites
  have hE := h.earTree
  have hI := invTree t d s h.inv' h.shape S.guards S.book
  have hR := rgTree G.σ G.n t d s h.ranges' h.shape h.σ_nodup h.σ_lt S.guards S.book S.rg
  have hC := ccTree G.σ G.n t d s h.ranges' h.shape h.σ_nodup h.σ_lt S.guards S.book S.frontiers
    S.rg S.cs h.close
  have hV := (h.coverTree S.guards S.book S.frontiers S.cs).2
  have hS := (stWalk G.g).1 t d s G.prev G.fs G.segs G.P G.X h.d_fs h.inv' h.shape S.guards S.book
    h.full h.v_lt h.e_lt h.nodup.of_append_right h.enodup h.Pv h.Pe h.height h.tstack h.segRead
    h.base_out h.stItems
  have hRR : wp (walkTree t d) (fun _ s' => s.g.TwoConnected → ∀ dp, d = dp + 1 →
      RWalk G.dfs G.F t.v d s' ∧ BotKeep s.tstack.length s s' ∧
        (∀ k, k < d → s'.stackVerts[k]! = s.stackVerts[k]!) ∧ Items.RSkelInv s'.g s'.items) s := by
    by_cases h2 : s.g.TwoConnected
    · rcases d with _ | dp
      · exact wp_of_forall fun _ _ _ dp hd => by omega
      · have R := h.r h2 dp rfl
        have hside := S.rside h2 dp rfl
        have hW := rrTree t (dp + 1) s G.F dp rfl h.inv' h.shape S.guards S.book hside h2 R.spec
          R.rooted R.frames R.top
        have hK := rkTree t (dp + 1) s G.F dp rfl h.inv' h.shape S.guards S.book hside h2 R.spec
          R.rooted R.frames R.top R.skel
        refine wp_mono _ (wp_and hW hK) fun _ s' ⟨⟨hw, hk, hg, hsv⟩, hskel⟩ _ dp' hd => ?_
        cases hd
        refine ⟨?_, hk, hsv, by rw [hg]; exact hskel⟩
        cases t with
        | node v outs => exact hw v outs rfl
    · exact wp_of_forall fun _ _ h2' => absurd h2' h2
  refine wp_mono _ (wp_and hE.2 (wp_and hI (wp_and hR (wp_and hC (wp_and hV (wp_and hS hRR))))))
    fun _ s' ⟨hte, ⟨hi, hs⟩, hrs, hcl, ⟨hpl, how, hvc⟩, ⟨hg, hsd, hvs, hsr, hst⟩, hrr⟩ => ?_
  have hvc' : s'.stackVerts[d]! = t.v ∧ VertCover t.v s' := by
    cases t with
    | node v outs => exact hvc v outs rfl
  exact ⟨hte, hi, hs, hrs, hcl, hpl, how, hvc'.1, hvc'.2, hrr, hg, hsd, hvs, hsr, hst⟩

end WalkState
end Spqr
