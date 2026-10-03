import Spqr.Proofs.RSkelWalk
import Spqr.EarWalk
import Spqr.Proofs.Postorder
import Spqr.Proofs.ForestSpec

/-!
# `Items.RSkelInv` at the roots and across `walkForest`

Depth-0 glue for the walk-level invariant: a root out-edge ends in the boundary branch of
`finishEdge` (`keepsR_finishBoundary`), a root tree edge enters `rkTree` at depth 1 with an empty
frame list, and `walkForest`'s pop/append onto `rootItem` touches only the parentless `F` item.
-/

namespace Spqr
open WalkState WalkM
namespace WalkState

theorem init_rSkelInv (g : Graph) (tern : Bool) : Items.RSkelInv g (WalkState.init g tern).items := by
  intro i _ hty
  exfalso
  have h : Items.type (WalkState.init g tern).items i = Items.type (initialItems g) i := rfl
  rw [h, Items.initialItems_type] at hty
  split_ifs at hty

theorem walkOutPre_root (v : Nat) (o : DfsOut) (s : WalkState) (Q : Bool → WalkState → Prop) :
    wp (walkOutPre v 0 o false) Q s = Q false { s with stackDir := s.stackDir.set! 0 false } := by
  unfold walkOutPre
  (try simp only [wp_bind, wp_stackDir, wp_setStackDir, wp_ite, wp_pure])
  (try rfl)

/-- `rootOuts`'s per-edge content: a root out-edge returns `false` and leaves the stack empty. -/
theorem rootOut (v : Nat) (o : DfsOut) (s : WalkState) (hb : BookOut v 0 o false s) (hts : s.tstack = []) :
    wp (walkOut v 0 o false) (fun hv s' => hv = false ∧ s'.tstack = []) s := by
  unfold BookOut at hb
  have hb₁ := hb.2
  rw [walkOut_eq, wp_bind]
  rw [walkOutPre_root] at hb₁ ⊢
  unfold walkOutRest
  rw [wp_bind, wp_tstackSize]
  simp only [hts, List.length_nil] at hb₁ ⊢
  cases o with
  | back e cls dest =>
    exact finishEdge_root hb₁ fun _ => by simp
  | tree e cls child =>
    simp only [wp_bind, wp_modify] at hb₁ ⊢
    exact wp_imp (wp_of_forall fun _ s₃ hb₃ => finishEdge_root hb₃ fun h => by
      rw [hb₃.tree.2 ⟨_, _, _, rfl⟩] at h; cases h) hb₁.2

theorem rkRootOut (dfs : DfsData) (v : Nat) (o : DfsOut) (s : WalkState)
    (hi : s.Inv' 0) (hs : Shape s) (hb : BookOut v 0 o false s) (hts : s.tstack = [])
    (h2 : s.g.TwoConnected) (hsp : dfs.Spec s.g) (hrt : dfs.Rooted s.g)
    (hp : ∀ e cls c couts, o = .tree e cls (.node c couts) →
      dfs.IsParent s.stackVerts[0]! c ∧ couts = dfs.outs c)
    (hinv : Items.RSkelInv s.g s.items) :
    wp (walkOut v 0 o false) (fun _ s' => Items.RSkelInv s.g s'.items) s := by
  have hg := gbOut v 0 o false s hb
  unfold GuardsOut at hg
  unfold BookOut at hb
  have hb₁ := hb.2
  rw [walkOut_eq, wp_bind]
  rw [walkOutPre_root] at hb₁ hg ⊢
  unfold walkOutRest
  rw [wp_bind, wp_tstackSize]
  have hlen : ({ s with stackDir := s.stackDir.set! 0 false } : WalkState).tstack.length = 0 := by
    simp [hts]
  rw [hlen] at hb₁ hg ⊢
  cases o with
  | back e cls dest =>
    try simp only at hg hb₁
    simp only
    have hnt : (DfsOut.back e cls dest).cls.isTree = false := by
      cases h : (DfsOut.back e cls dest).cls.isTree
      · rfl
      · obtain ⟨_, _, _, h'⟩ := hb₁.tree.1 h; cases h'
    have hD : (0 : Nat) = if (DfsOut.back e cls dest).cls.isTree then 0 + 1 else 0 := by simp [hnt]
    have hok := ear_boundary (Nat.zero_le _) hg hi.frame' hs.frame' hb₁ hD
    exact keepsR_finishBoundary v 0 _ 0 false (Nat.zero_le _) hs.frame' hb₁.v_lt hb₁.e_lt hok hinv
  | tree e cls child =>
    obtain ⟨c, couts⟩ := child
    try simp only [wp_modify] at hg hb₁
    simp only [wp_bind, wp_modify]
    set s₂ : WalkState := { s with stackDir := s.stackDir.set! 0 false, firstOccurrence := s.firstOccurrence.set! 0 s.g.ne } with hs₂
    have hi₂ : s₂.Inv' 0 := hi.frame'
    have hs₂' : Shape s₂ := hs.frame'
    have hpre₂ : ∀ w outs, DfsTree.node c couts = .node w outs →
        ({ s₂ with stackVerts := s₂.stackVerts.set! 1 w } : WalkState).Inv' 1 :=
      fun w _ _ => hi₂.setSv w
    have hpar : s₂.RInvTop dfs s₂.stackVerts[0]! 0 :=
      ⟨fun t ht => by simp [hs₂, hts] at ht, by simp [hs₂, hts]⟩
    have hpc := hp e cls c couts rfl
    have hfr := walkTree_frontiers (.node c couts) 1 s₂ hpre₂ hs₂' hg.1 hb₁.1
    have hr := (walkTree_rSide s₂ 0 c couts hi₂ hs₂' hg.1 hfr h2 hsp hrt hpc.1 hpc.2 hpar).2
    have hT : Types s₂.g s₂ := ⟨rfl, hs₂'.size, hs₂'.root, hs₂'.vert, hs₂'.edge⟩
    refine wp_imp (wp_imp (wp_imp (wp_imp (wp_of_forall fun _ s₃ ⟨hg₃, hb₃⟩ ⟨hi₃, hs₃⟩ hinv₃
        ⟨hg₃', _, _, _, _⟩ => ?_)
        (wp_and hg.2 hb₁.2))
        (invTree (.node c couts) 1 s₂ hpre₂ hs₂' hg.1 hb₁.1))
        (rkTree (.node c couts) 1 s₂ [] 0 rfl hpre₂ hs₂' hg.1 hb₁.1 hr h2 hsp hrt
          (fun f hf => by simp at hf) hpar hinv))
        (walkTree_frame (.node c couts) 1 s₂ hT)
    have hD : (1 : Nat) = if (DfsOut.tree e cls (.node c couts)).cls.isTree then 0 + 1 else 0 := by
      simp [hb₃.tree.2 ⟨_, _, _, rfl⟩]
    have hok := ear_boundary (Nat.zero_le _) hg₃ hi₃ hs₃ hb₃ hD
    have hinv₃' : Items.RSkelInv s₃.g s₃.items := by rw [hg₃']; exact hinv₃
    have := keepsR_finishBoundary v 0 _ 0 false (Nat.zero_le _) hs₃ hb₃.v_lt hb₃.e_lt hok hinv₃'
    rw [hg₃'] at this
    exact this

theorem rkRootOuts (dfs : DfsData) : ∀ (v : Nat) (outs : List DfsOut) (s : WalkState),
    s.Inv' 0 → Shape s → BookOuts v 0 outs false s → s.tstack = [] →
    s.g.TwoConnected → dfs.Spec s.g → dfs.Rooted s.g →
    (∀ o ∈ outs, ∀ e cls c couts, o = .tree e cls (.node c couts) →
      dfs.IsParent s.stackVerts[0]! c ∧ couts = dfs.outs c) →
    Items.RSkelInv s.g s.items →
    wp (walkOuts v 0 outs false) (fun _ s' => Items.RSkelInv s.g s'.items) s
  | v, [], s, _, _, _, _, _, _, _, _, hinv => by unfold walkOuts; simpa [wp_pure] using hinv
  | v, o :: rest, s, hi, hs, hb, hts, h2, hsp, hrt, hp, hinv => by
    unfold BookOuts at hb
    have hg := gbOut v 0 o false s hb.1
    have hT : Types s.g s := ⟨rfl, hs.size, hs.root, hs.vert, hs.edge⟩
    have hk := kOut v 0 o false s s.g 1 rootItem s hT (Nat.le_refl _) (by show 0 < 1 + _ + _; omega)
      (by show 1 + v ≠ 0; omega) (fun w _ => by show 1 + w ≠ 0; omega)
      (fun e _ => by show 1 + s.g.nv + e ≠ 0; omega) Keep.refl
    have hout := rkRootOut dfs v o s hi hs hb.1 hts h2 hsp hrt (hp o (List.mem_cons_self ..)) hinv
    unfold walkOuts
    simp only [wp_bind]
    refine wp_imp (wp_imp (wp_imp (wp_imp (wp_imp (wp_of_forall fun hv' s' ⟨hi', hs'⟩ ⟨hv0, hts'⟩ hk'
        hinv' hb' => ?_)
      (invOut v 0 o false s hi hs hg hb.1)) (rootOut v o s hb.1 hts)) hk) hout) hb.2
    subst hv0
    have hg' := hk'.g
    refine wp_mono _ (rkRootOuts dfs v rest s' hi' hs' hb' hts' (by rw [hg']; exact h2)
      (by rw [hg']; exact hsp) (by rw [hg']; exact hrt) ?_ (by rw [hg']; exact hinv'))
      fun _ s'' h => by rw [hg'] at h; exact h
    intro o' ho' e cls c couts heq
    rw [hk'.svlo 0 (by omega)]
    exact hp o' (List.mem_cons_of_mem _ ho') e cls c couts heq

theorem rkRootTree (dfs : DfsData) (t : DfsTree) (s : WalkState) (hi : s.Inv' 0) (hs : Shape s)
    (hb : BookTree t 0 s) (hts : s.tstack = []) (hsz : 0 < s.stackVerts.size)
    (h2 : s.g.TwoConnected) (hsp : dfs.Spec s.g) (hrt : dfs.Rooted s.g)
    (hp : ∀ v outs, t = .node v outs → ∀ o ∈ outs, ∀ e cls c couts,
      o = .tree e cls (.node c couts) → dfs.IsParent v c ∧ couts = dfs.outs c)
    (hinv : Items.RSkelInv s.g s.items) :
    wp (walkTree t 0) (fun _ s' => Items.RSkelInv s.g s'.items) s := by
  obtain ⟨v, outs⟩ := t
  unfold walkTree
  simp only [wp_bind, wp_modify]
  unfold BookTree at hb
  have h₁ := rkRootOuts dfs v outs { s with stackVerts := s.stackVerts.set! 0 v }
    (hi.stackVerts_of_nil hts _) hs.frame' hb hts h2 hsp hrt (fun o ho e cls c couts heq => by
      rw [show ({ s with stackVerts := s.stackVerts.set! 0 v } : WalkState).stackVerts[0]! = v from
        Array.getElem!_set!_self _ _ _ hsz]
      exact hp v outs rfl o ho e cls c couts heq) hinv
  have h₂ := rootOuts v outs _ hb hts
  refine wp_mono _ (wp_and h₁ h₂) fun hv s' ⟨hinv', hv0, hts'⟩ => ?_
  subst hv0
  simp only [Bool.false_eq_true, ↓reduceIte, wp_bind, wp_setStackDir, wp_pushVertTstack]
  exact hinv'

theorem rkForest (g : Graph) (dfs : DfsData) : ∀ (forest pre : List DfsTree) (s : WalkState),
    RootState g pre s → ForestOK g (pre ++ forest) → (∀ t ∈ forest, t.WF []) →
    (∀ t ∈ forest, t.Ends g) → g.TwoConnected → dfs.Spec g → dfs.Rooted g →
    (∀ t ∈ forest, ∀ v outs, t = .node v outs → ∀ o ∈ outs, ∀ e cls c couts,
      o = .tree e cls (.node c couts) → dfs.IsParent v c ∧ couts = dfs.outs c) →
    Items.RSkelInv g s.items → wp (walkForest forest) (fun _ s' => Items.RSkelInv g s'.items) s
  | [], _, _, _, _, _, _, _, _, _, _, hinv => hinv
  | t :: rest, pre, s, h, hf, hwf, hends, h2, hsp, hrt, hp, hinv => by
    have hb := h.book hf (hwf t (by simp)) (hends t (by simp))
    have hg := gbTree t 0 s hb
    have hi' : ∀ v outs, t = .node v outs →
        ({ s with stackVerts := s.stackVerts.set! 0 v } : WalkState).Inv' 0 :=
      fun _ _ _ => h.inv.stackVerts_of_nil h.tstack _
    have hnv : 0 < g.nv := by
      obtain ⟨v, outs⟩ := t
      exact Nat.lt_of_le_of_lt (Nat.zero_le _) (RootState.hvlt hf v (by simp [DfsTree.verts]))
    have hstep := h.step hf (hwf t (by simp)) (hends t (by simp))
    have hplace := (walk_place_aux g).1 t 0 _ _ s h.place (RootState.hvlt hf) (RootState.helt hf)
      (RootState.hvn hf).1 (RootState.hen hf).1 (RootState.hPv hf) (RootState.hPe hf)
    have hinv₁ := rkRootTree dfs t s h.inv h.shape hb h.tstack (by rw [h.sv]; exact hnv)
      (by rw [h.g_eq]; exact h2) (by rw [h.g_eq]; exact hsp) (by rw [h.g_eq]; exact hrt)
      (hp t (by simp)) (by rw [h.g_eq]; exact hinv)
    have hinvT := invTree t 0 s hi' h.shape hg hb
    have hrk : wp (walkTree t 0) (fun _ s' => RootOK s') s :=
      walkTree_rootOK t s hi' h.shape hg hb h.tstack (by rw [h.sd]; exact hnv)
    show wp ((walkTree t 0 >>= fun _ => popTstack >>= fun top =>
      modifyItem rootItem fun it => { it with ch := it.ch ++ top.spans.2 }) >>= fun _ => walkForest rest) _ s
    rw [wp_bind, wp_bind]
    refine wp_mono _ (wp_and hstep (wp_and hplace (wp_and hinv₁ (wp_and hinvT hrk))))
      fun _ s₁ ⟨hst₁, hp₁, hinv₁, ⟨_, hs₁⟩, hrk₁⟩ => ?_
    obtain ⟨tt, hts₁, -, -⟩ := hrk₁
    simp only [wp_bind, wp_popTstack, wp_modifyItem] at hst₁ ⊢
    have hroot : ∀ p, ¬ Items.IsParent s₁.items p rootItem := noParent_of_cnt_eq_zero hp₁.root
    have hrootT : Items.type s₁.items rootItem ≠ .R := by rw [hs₁.root]; decide
    rw [h.g_eq] at hinv₁
    have hinv₂ : Items.RSkelInv g (s₁.items.modify rootItem fun it =>
        { it with ch := it.ch ++ s₁.tstack.head!.spans.2 }) :=
      hinv₁.modify_root_of_ne hroot _ (fun _ => rfl) hrootT
    exact rkForest g dfs rest (pre ++ [t]) _ hst₁ (by simpa using hf)
      (fun t' ht' => hwf t' (by simp [ht'])) (fun t' ht' => hends t' (by simp [ht'])) h2 hsp hrt
      (fun t' ht' => hp t' (by simp [ht'])) hinv₂

end WalkState

/-- **Named admission.** Exact obligation: the DFS data of `g.dfsForest vo eo` is rooted for a
2-connected `g` — `(DfsData.ofForest (g.dfsForest vo eo)).Rooted g`, i.e. every endpoint of
every edge of `g` is a descendant (`Anc`) of `dfs.root`, the first root of the forest. For a
2-connected (hence connected) graph the forest is a single tree whose vertex set is all of
`g.nv` (`dfsForest_spanning'`), so `Anc root x` follows from the tree structure; the connectivity
→ single-tree argument is DFS-layer reasoning outside the R proof files. -/
theorem dfsForest_rooted (g : Graph) (hg : g.WF) (vo eo : List Nat)
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) (h2 : g.TwoConnected) :
    (DfsData.ofForest (g.dfsForest vo eo)).Rooted g := by
  sorry

/-- Every R item of the walk's items is a 3-connected skeleton (`Items.RSkelInv`), threading the
invariant from `WalkState.init` through every root of `g.dfsForest vo eo`. -/
theorem walk_rSkelInv (g : Graph) (hg : g.WF) (tern : Bool) (vo eo : List Nat)
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) (h2 : g.TwoConnected) :
    Items.RSkelInv g (g.walk tern (g.dfsForest vo eo)).items := by
  obtain ⟨hvp, hep⟩ := dfsForest_spanning' hg hvo heo
  have hf : ForestOK g (g.dfsForest vo eo) := ForestOK.of_perm hvp hep
  have hnd : ((g.dfsForest vo eo).flatMap DfsTree.verts).Nodup := hvp.nodup_iff.2 List.nodup_range
  have hspec := dfsForestSpec_of_dfsForest hg hvo heo
  have hrt := dfsForest_rooted g hg vo eo hvo heo h2
  have hp : ∀ t ∈ g.dfsForest vo eo, ∀ v outs, t = .node v outs → ∀ o ∈ outs, ∀ e cls c couts,
      o = .tree e cls (.node c couts) →
      (DfsData.ofForest (g.dfsForest vo eo)).IsParent v c ∧
        couts = (DfsData.ofForest (g.dfsForest vo eo)).outs c := by
    intro t ht v outs htv o ho e cls c couts hoe
    subst htv; subst hoe
    have h1 : (DfsData.ofForest (g.dfsForest vo eo)).outs v = outs :=
      DfsData.ofForest_outs hnd ht (DfsTree.Sub.refl _)
    have h2 : (DfsData.ofForest (g.dfsForest vo eo)).outs c = couts :=
      DfsData.ofForest_outs hnd ht (DfsTree.Sub.step (DfsTree.Sub.refl _) ho)
    exact ⟨⟨_, by rw [h1]; exact ho, rfl, rfl⟩, h2.symm⟩
  exact rkForest g _ (g.dfsForest vo eo) [] (WalkState.init g tern) (rootState_init g tern) hf
    (dfsForest_wf hg hvo heo) (dfsForest_ends g hg hvo heo) h2 hspec.toSpec hrt hp
    (init_rSkelInv g tern)

end Spqr
