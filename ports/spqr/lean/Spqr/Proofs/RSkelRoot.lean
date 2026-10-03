import Spqr.Proofs.RSkelWalk
import Spqr.Proofs.RSide
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
    (hwf : s.g.WF) (hsv : s.stackVerts.size = s.g.nv)
    (hp : ∀ e cls c couts, o = .tree e cls (.node c couts) →
      dfs.IsParent s.stackVerts[0]! c ∧ ∀ t : DfsTree, t.Sub (.node c couts) → dfs.outs t.v = t.outs)
    (hd0 : dfs.depth s.stackVerts[0]! = 0)
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
    have hchain : ∀ k, k ≤ 0 → dfs.Anc s₂.stackVerts[k]! s₂.stackVerts[0]! ∧
        dfs.depth s₂.stackVerts[k]! = k := fun k hk => by
      obtain rfl : k = 0 := Nat.le_zero.1 hk
      exact ⟨Relation.ReflTransGen.refl, hd0⟩
    have hr := walkTree_rSide s₂ 0 c couts hi₂ hs₂' hg.1 hb₁.1 h2 hwf hsp hrt hsv hpc.1 hpc.2 hchain hpar
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
    s.g.TwoConnected → dfs.Spec s.g → dfs.Rooted s.g → s.g.WF → s.stackVerts.size = s.g.nv →
    (∀ o ∈ outs, ∀ e cls c couts, o = .tree e cls (.node c couts) →
      dfs.IsParent s.stackVerts[0]! c ∧ ∀ t : DfsTree, t.Sub (.node c couts) → dfs.outs t.v = t.outs) →
    dfs.depth s.stackVerts[0]! = 0 →
    Items.RSkelInv s.g s.items →
    wp (walkOuts v 0 outs false) (fun _ s' => Items.RSkelInv s.g s'.items) s
  | v, [], s, _, _, _, _, _, _, _, _, _, _, _, hinv => by unfold walkOuts; simpa [wp_pure] using hinv
  | v, o :: rest, s, hi, hs, hb, hts, h2, hsp, hrt, hwf, hsv, hp, hd0, hinv => by
    unfold BookOuts at hb
    have hg := gbOut v 0 o false s hb.1
    have hT : Types s.g s := ⟨rfl, hs.size, hs.root, hs.vert, hs.edge⟩
    have hk := kOut v 0 o false s s.g 1 rootItem s hT (Nat.le_refl _) (by show 0 < 1 + _ + _; omega)
      (by show 1 + v ≠ 0; omega) (fun w _ => by show 1 + w ≠ 0; omega)
      (fun e _ => by show 1 + s.g.nv + e ≠ 0; omega) Keep.refl
    have hout := rkRootOut dfs v o s hi hs hb.1 hts h2 hsp hrt hwf hsv (hp o (List.mem_cons_self ..)) hd0 hinv
    unfold walkOuts
    simp only [wp_bind]
    refine wp_imp (wp_imp (wp_imp (wp_imp (wp_imp (wp_of_forall fun hv' s' ⟨hi', hs'⟩ ⟨hv0, hts'⟩ hk'
        hinv' hb' => ?_)
      (invOut v 0 o false s hi hs hg hb.1)) (rootOut v o s hb.1 hts)) hk) hout) hb.2
    subst hv0
    have hg' := hk'.g
    refine wp_mono _ (rkRootOuts dfs v rest s' hi' hs' hb' hts' (by rw [hg']; exact h2)
      (by rw [hg']; exact hsp) (by rw [hg']; exact hrt) (by rw [hg']; exact hwf)
      (by rw [hk'.sv, hg']; exact hsv) ?_ (by rw [hk'.svlo 0 (by omega)]; exact hd0)
      (by rw [hg']; exact hinv'))
      fun _ s'' h => by rw [hg'] at h; exact h
    intro o' ho' e cls c couts heq
    rw [hk'.svlo 0 (by omega)]
    exact hp o' (List.mem_cons_of_mem _ ho') e cls c couts heq

theorem rkRootTree (dfs : DfsData) (t : DfsTree) (s : WalkState) (hi : s.Inv' 0) (hs : Shape s)
    (hb : BookTree t 0 s) (hts : s.tstack = []) (hsz : 0 < s.stackVerts.size)
    (h2 : s.g.TwoConnected) (hsp : dfs.Spec s.g) (hrt : dfs.Rooted s.g)
    (hwf : s.g.WF) (hsv : s.stackVerts.size = s.g.nv)
    (hp : ∀ v outs, t = .node v outs → ∀ o ∈ outs, ∀ e cls c couts,
      o = .tree e cls (.node c couts) → dfs.IsParent v c ∧ ∀ t : DfsTree, t.Sub (.node c couts) → dfs.outs t.v = t.outs)
    (hd0 : ∀ v outs, t = .node v outs → dfs.depth v = 0)
    (hinv : Items.RSkelInv s.g s.items) :
    wp (walkTree t 0) (fun _ s' => Items.RSkelInv s.g s'.items) s := by
  obtain ⟨v, outs⟩ := t
  unfold walkTree
  simp only [wp_bind, wp_modify]
  unfold BookTree at hb
  have h₁ := rkRootOuts dfs v outs { s with stackVerts := s.stackVerts.set! 0 v }
    (hi.stackVerts_of_nil hts _) hs.frame' hb hts h2 hsp hrt hwf
    (by show (s.stackVerts.set! 0 v).size = s.g.nv; rw [Array.size_set!]; exact hsv)
    (fun o ho e cls c couts heq => by
      rw [show ({ s with stackVerts := s.stackVerts.set! 0 v } : WalkState).stackVerts[0]! = v from
        Array.getElem!_set!_self _ _ _ hsz]
      exact hp v outs rfl o ho e cls c couts heq)
    (by
      rw [show ({ s with stackVerts := s.stackVerts.set! 0 v } : WalkState).stackVerts[0]! = v from
        Array.getElem!_set!_self _ _ _ hsz]
      exact hd0 v outs rfl) hinv
  have h₂ := rootOuts v outs _ hb hts
  refine wp_mono _ (wp_and h₁ h₂) fun hv s' ⟨hinv', hv0, hts'⟩ => ?_
  subst hv0
  simp only [Bool.false_eq_true, ↓reduceIte, wp_bind, wp_setStackDir, wp_pushVertTstack]
  exact hinv'

theorem rkForest (g : Graph) (dfs : DfsData) : ∀ (forest pre : List DfsTree) (s : WalkState),
    RootState g pre s → ForestOK g (pre ++ forest) → (∀ t ∈ forest, t.WF []) →
    (∀ t ∈ forest, t.Ends g) → g.TwoConnected → g.WF → dfs.Spec g → dfs.Rooted g →
    (∀ t ∈ forest, ∀ v outs, t = .node v outs → ∀ o ∈ outs, ∀ e cls c couts,
      o = .tree e cls (.node c couts) → dfs.IsParent v c ∧ ∀ t : DfsTree, t.Sub (.node c couts) → dfs.outs t.v = t.outs) →
    (∀ t ∈ forest, dfs.depth t.v = 0) →
    Items.RSkelInv g s.items → wp (walkForest forest) (fun _ s' => Items.RSkelInv g s'.items) s
  | [], _, _, _, _, _, _, _, _, _, _, _, _, hinv => hinv
  | t :: rest, pre, s, h, hf, hwf, hends, h2, hgw, hsp, hrt, hp, hd0, hinv => by
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
      (by rw [h.g_eq]; exact hgw) (by rw [h.sv, h.g_eq])
      (hp t (by simp)) (fun v outs htv => by subst htv; exact hd0 (.node v outs) (List.mem_cons_self ..))
      (by rw [h.g_eq]; exact hinv)
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
      (fun t' ht' => hwf t' (by simp [ht'])) (fun t' ht' => hends t' (by simp [ht'])) h2 hgw hsp hrt
      (fun t' ht' => hp t' (by simp [ht'])) (fun t' ht' => hd0 t' (by simp [ht'])) hinv₂

end WalkState

/-- The DFS data of `g.dfsForest vo eo`, re-rooted at a depth-0 vertex `r`, is rooted for a
2-connected `g`: every endpoint of every edge is a descendant (`Anc`) of `r`. `r` is the root of
the tree holding the edges — `ofForest`'s own `root` (the *first* tree's root) does not work, since
an isolated vertex ordered first becomes `root` with no out-edges (`checks/RRootedCounter.lean`,
kernel-checked). Proof: an out-edge joins `v` to a child or to an ancestor (`Spec.back_anc`), so
the endpoints of every edge are `Anc`-comparable, and a depth-0 ancestor of one endpoint is an
ancestor of the other (`Anc.comparable`); 2-connectivity joins any edge to edge `0` by a walk, and
the ancestor relation propagates along it. -/
theorem dfsForest_rooted (g : Graph) (hg : g.WF) (vo eo : List Nat)
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) (h2 : g.TwoConnected) :
    ∃ r, (DfsData.ofForest (g.dfsForest vo eo)).depth r = 0 ∧
      ({ DfsData.ofForest (g.dfsForest vo eo) with root := r } : DfsData).Rooted g := by
  obtain ⟨hvp, -⟩ := dfsForest_spanning' hg hvo heo
  have hnd : ((g.dfsForest vo eo).flatMap DfsTree.verts).Nodup := hvp.nodup_iff.2 List.nodup_range
  have hsp : (DfsData.ofForest (g.dfsForest vo eo)).Spec g :=
    (dfsForestSpec_of_dfsForest hg hvo heo).toSpec
  generalize hdfs : DfsData.ofForest (g.dfsForest vo eo) = dfs at hsp ⊢
  have hbound : ∀ {e a b}, g.Joins e a b → a < g.nv ∧ b < g.nv := by
    intro e a b h
    rcases h with h | h
    · obtain ⟨_, h⟩ := Array.getElem?_eq_some_iff.1 h
      exact hg _ (h ▸ Array.getElem_mem _)
    · obtain ⟨_, h⟩ := Array.getElem?_eq_some_iff.1 h
      exact (hg _ (h ▸ Array.getElem_mem _)).symm
  have hpair : ∀ {e a b x y}, g.Joins e a b → g.Joins e x y → (x = a ∧ y = b) ∨ (x = b ∧ y = a) := by
    intro e a b x y hab hxy
    rcases hab with h | h <;> rcases hxy with h' | h' <;>
      have := h'.symm.trans h <;> simp only [Option.some.injEq, Prod.mk.injEq] at this
    · exact .inl this
    · exact .inr ⟨this.2, this.1⟩
    · exact .inr this
    · exact .inl ⟨this.2, this.1⟩
  have hends : ∀ {e a b x}, g.Joins e a b → g.IsEnd e x → x = a ∨ x = b :=
    fun hab ⟨_, hx⟩ => (hpair hab hx).elim (fun h => .inl h.1) (fun h => .inr h.1)
  have hstep : ∀ r x y, dfs.depth r = 0 → dfs.Anc r x → g.Adj x y → dfs.Anc r y := by
    intro r x y hr hx ⟨e, hxy⟩
    obtain ⟨w, o, ho, rfl⟩ := hsp.edge_out e hxy.lt
    have hj := hsp.joins w o ho
    have hcomp : dfs.Anc w o.dest ∨ dfs.Anc o.dest w := by
      cases hdt : o.isTree
      · exact .inr (hsp.back_anc w o ho hdt)
      · exact .inl (DfsData.IsParent.anc ⟨o, ho, hdt, rfl⟩)
    have key : ∀ a b, dfs.Anc r a → (dfs.Anc a b ∨ dfs.Anc b a) → dfs.Anc r b := by
      intro a b hra hab
      rcases hab with hab | hba
      · exact hra.trans hab
      · rcases hra.comparable hsp hba with h | h
        · exact h
        · obtain rfl : b = r := by
            by_contra hne
            have := h.depth_lt hsp hne; omega
          exact DfsData.Anc.refl _
    rcases hpair hj hxy with ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩
    · exact key _ _ hx hcomp
    · exact key _ _ hx hcomp.symm
  by_cases hne : g.ne = 0
  · exact ⟨dfs.root, hsp.depth_root, fun e x hx => absurd hx.lt (by omega)⟩
  obtain ⟨w, o, ho, hoe⟩ := hsp.edge_out 0 (Nat.pos_of_ne_zero hne)
  have hj : g.Joins 0 w o.dest := hoe ▸ hsp.joins w o ho
  obtain ⟨t, ht, hwt⟩ := List.mem_flatMap.1
    (hvp.mem_iff.2 (List.mem_range.2 (hbound hj).1))
  refine ⟨t.v, hdfs ▸ DfsData.ofForest_depth_root hnd ht, fun e x hx => ?_⟩
  have hr0 : dfs.depth t.v = 0 := hdfs ▸ DfsData.ofForest_depth_root hnd ht
  have hrw : dfs.Anc t.v w :=
    hdfs ▸ (DfsData.ofForest_anc_iff hnd ht (DfsTree.Sub.refl _)).2 hwt
  have hend0 : ∀ z, g.IsEnd 0 z → dfs.Anc t.v z := by
    intro z hz
    rcases hends hj hz with rfl | rfl
    · exact hrw
    · exact hstep _ _ _ hr0 hrw ⟨0, hj⟩
  show dfs.Anc t.v x
  rcases h2 g.nv 0 e (Nat.pos_of_ne_zero hne) hx.lt with rfl | ⟨x₀, y, hx₀, hy, hreach⟩
  · exact hend0 x hx
  · have hy' : dfs.Anc t.v y := by
      clear hy
      induction hreach with
      | refl _ => exact hend0 _ hx₀
      | tail _ hadj _ ih => exact hstep _ _ _ hr0 ih hadj
    obtain ⟨y', hyy'⟩ := hy
    rcases hends hyy' hx with rfl | rfl
    · exact hy'
    · exact hstep _ _ _ hr0 hy' ⟨e, hyy'⟩

/-- Every R item of the walk's items is a 3-connected skeleton (`Items.RSkelInv`), threading the
invariant from `WalkState.init` through every root of `g.dfsForest vo eo`. -/
theorem walk_rSkelInv (g : Graph) (hg : g.WF) (tern : Bool) (vo eo : List Nat)
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) (h2 : g.TwoConnected) :
    Items.RSkelInv g (g.walk tern (g.dfsForest vo eo)).items := by
  obtain ⟨hvp, hep⟩ := dfsForest_spanning' hg hvo heo
  have hf : ForestOK g (g.dfsForest vo eo) := ForestOK.of_perm hvp hep
  have hnd : ((g.dfsForest vo eo).flatMap DfsTree.verts).Nodup := hvp.nodup_iff.2 List.nodup_range
  have hspec := dfsForestSpec_of_dfsForest hg hvo heo
  obtain ⟨r, hr0, hrt⟩ := dfsForest_rooted g hg vo eo hvo heo h2
  have hsp' : ({ DfsData.ofForest (g.dfsForest vo eo) with root := r } : DfsData).Spec g :=
    { hspec.toSpec with depth_root := hr0 }
  have hp : ∀ t ∈ g.dfsForest vo eo, ∀ v outs, t = .node v outs → ∀ o ∈ outs, ∀ e cls c couts,
      o = .tree e cls (.node c couts) →
      (DfsData.ofForest (g.dfsForest vo eo)).IsParent v c ∧
        ∀ t : DfsTree, t.Sub (.node c couts) →
          (DfsData.ofForest (g.dfsForest vo eo)).outs t.v = t.outs := by
    intro t ht v outs htv o ho e cls c couts hoe
    subst htv; subst hoe
    have h1 : (DfsData.ofForest (g.dfsForest vo eo)).outs v = outs :=
      DfsData.ofForest_outs hnd ht (DfsTree.Sub.refl _)
    exact ⟨⟨_, by rw [h1]; exact ho, rfl, rfl⟩,
      fun t' hsub => DfsData.ofForest_outs hnd ht (DfsTree.Sub.step hsub ho)⟩
  exact rkForest g _ (g.dfsForest vo eo) [] (WalkState.init g tern) (rootState_init g tern) hf
    (dfsForest_wf hg hvo heo) (dfsForest_ends g hg hvo heo) h2 hg hsp' hrt hp
    (fun t ht => DfsData.ofForest_depth_root hnd ht) (init_rSkelInv g tern)

end Spqr
