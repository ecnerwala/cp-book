import Spqr.RangesCoverTree

/-!
# Range-side site exports

Public per-site conclusions of the Ranges layer, for consumers that cannot run the walk inductions
themselves: `CloseBase` (the range/ear/frontier/R/close exports and the DFS site facts) holds at the
pre-state of every `finishEdge` call of the real walk (`CbTree`/`CbOuts`/`CbOut` on a tree,
`CbForest` on the DFS forest, `walk_closeBase` from the initial state), and from `CloseBase` the range
invariant is carried into `finishEdge`: `rangesInv_l1Iter` (every loop-1 iterate of `closeEars`) and
`rangesInv_feS₂` (the state after `mergeLate`, where a type-1 vertex close starts).
-/

namespace Spqr
open WalkM
namespace WalkState

variable {σ : List Nat} {n D : Nat} {s : WalkState}

/-! ### `CloseBase` at every `finishEdge` pre-state -/

mutual
/-- `CloseBase` at every `finishEdge` pre-state of `walkTree t d` (the structure of `CsTree`). -/
def CbTree (σ : List Nat) (n : Nat) (t : DfsTree) (d : Nat) (s : WalkState) : Prop :=
  match t with
  | .node v outs => CbOuts σ n v d outs false { s with stackVerts := s.stackVerts.set! d v }

def CbOuts (σ : List Nat) (n v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  match outs with
  | [] => True
  | o :: rest => CbOut σ n v d o hasVert s ∧
      wp (walkOut v d o hasVert) (fun hasVert' s' => CbOuts σ (n + o.block.length) v d rest hasVert' s') s

/-- `finishEdge curV d o origTstack hasVert'` is called at `s₃` (tree) / `s₁` (back) with
`origTstack = s₁.tstack.length` (`walkOutRest`). -/
def CbOut (σ : List Nat) (n v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  wp (walkOutPre v d o hasVert) (fun hasVert' s₁ =>
    match o with
    | .tree _ _ child =>
      wp (modify fun s => { s with firstOccurrence := s.firstOccurrence.set! d s.g.ne }) (fun _ s₂ =>
        CbTree σ n child (d + 1) s₂ ∧
        wp (walkTree child (d + 1)) (fun _ s₃ =>
          CloseBase σ (n + child.edgePostorder.length) (d + 1) v d o s₁.tstack.length hasVert' s₃) s₂) s₁
    | .back .. => CloseBase σ n d v d o s₁.tstack.length hasVert' s₁) s
end

abbrev CbInvTree (σ : List Nat) (n : Nat) (t : DfsTree) (d : Nat) (s : WalkState) : Prop :=
  (∀ v outs, t = .node v outs → ({ s with stackVerts := s.stackVerts.set! d v } : WalkState).RangesInv σ n d) →
  Shape s → σ.Nodup → (∀ e ∈ σ, e < s.g.ne) → GuardsTree t d s → BookTree t d s → FrontiersTree t d s →
  RgTree σ n t d s → CsTree σ n t d s → s.CloseInv → CbTree σ n t d s

abbrev CbInvOuts (σ : List Nat) (n v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  RgS σ n d s → σ.Nodup → GuardsOuts v d outs hasVert s → BookOuts v d outs hasVert s →
  FrontiersOuts v d outs hasVert s → RgOuts σ n v d outs hasVert s → CsOuts σ n v d outs hasVert s →
  s.CloseInv → CbOuts σ n v d outs hasVert s

abbrev CbInvOut (σ : List Nat) (n v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  RgS σ n d s → σ.Nodup → GuardsOut v d o hasVert s → BookOut v d o hasVert s →
  FrontiersOut v d o hasVert s → RgOut σ n v d o hasVert s → CsOut σ n v d o hasVert s →
  s.CloseInv → CbOut σ n v d o hasVert s

mutual
theorem cbTree : ∀ (σ : List Nat) (n : Nat) (t : DfsTree) (d : Nat) (s : WalkState), CbInvTree σ n t d s
  | σ, n, .node v outs, d, s => fun hi hs hnd hσ hg hb hf hr hc hcl => by
    unfold CbTree
    unfold GuardsTree at hg; unfold BookTree at hb; unfold FrontiersTree at hf; unfold RgTree at hr
    unfold CsTree at hc
    exact cbOuts σ n v d outs false _ ⟨hi v outs rfl, hs.frame', hσ⟩ hnd hg hb hf hr hc
      (hcl.frame rfl rfl (fun _ h => h))

theorem cbOuts : ∀ (σ : List Nat) (n v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState),
    CbInvOuts σ n v d outs hasVert s
  | σ, n, v, d, [], hasVert, s => fun _ _ _ _ _ _ _ _ => trivial
  | σ, n, v, d, o :: rest, hasVert, s => fun hrs hnd hg hb hf hr hc hcl => by
    unfold GuardsOuts at hg; unfold BookOuts at hb; unfold FrontiersOuts at hf; unfold RgOuts at hr
    unfold CsOuts at hc
    unfold CbOuts
    refine ⟨cbOut σ n v d o hasVert s hrs hnd hg.1 hb.1 hf.1 hr.1 hc.1 hcl, ?_⟩
    exact wp_mono _ (wp_and (rgOut σ n v d o hasVert s hrs hnd hg.1 hb.1 hr.1) (wp_and hg.2 (wp_and hb.2
        (wp_and hf.2 (wp_and hr.2 (wp_and hc.2 (ccOut σ n v d o hasVert s hrs hnd hg.1 hb.1 hf.1 hr.1 hc.1 hcl)))))))
      fun hv' s' ⟨hrs', hg', hb', hf', hr', hc', hcl'⟩ =>
        cbOuts σ (n + o.block.length) v d rest hv' s' hrs' hnd hg' hb' hf' hr' hc' hcl'

theorem cbOut : ∀ (σ : List Nat) (n v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState), CbInvOut σ n v d o hasVert s
  | σ, n, v, d, o, hasVert, s => fun ⟨hi, hs, hσ⟩ hnd hg hb hf hr hc hcl => by
    unfold CbOut
    unfold GuardsOut at hg; unfold BookOut at hb; unfold FrontiersOut at hf; unfold RgOut at hr
    unfold CsOut at hc
    refine wp_mono _ (wp_and (walkOutPre_ranges hi hs hσ hb.1 hr.1) (wp_and hg (wp_and hb.2 (wp_and hf
        (wp_and hr.2 (wp_and hc (walkOutPre_closeInv' hcl hs d o hb.1)))))))
      fun hv' s₁ ⟨⟨hi₁, hs₁, hσ₁⟩, hg₁, hb₁, hf₁, hr₁, hc₁, hcl₁⟩ => ?_
    cases o with
    | back e cls dest =>
      exact ⟨by rw [if_neg fun h => by obtain ⟨_, _, _, h⟩ := hb₁.tree.1 h; cases h],
        hnd, ⟨hi₁, hs₁, hσ₁⟩, hg₁, hb₁, hf₁, hr₁, hcl₁, hc₁⟩
    | tree e cls child =>
      simp only [wp_modify] at hg₁ hb₁ hf₁ hr₁ hc₁ ⊢
      have pre : ∀ w outs, child = .node w outs →
          ({ { s₁ with firstOccurrence := s₁.firstOccurrence.set! d s₁.g.ne } with
            stackVerts := s₁.stackVerts.set! (d + 1) w } : WalkState).RangesInv σ n (d + 1) :=
        fun w outs _ =>
          (hi₁.frame' (s' := { s₁ with firstOccurrence := s₁.firstOccurrence.set! d s₁.g.ne })).setSv w
      refine ⟨cbTree σ n child (d + 1) _ pre hs₁.frame' hnd hσ₁ hg₁.1 hb₁.1 hf₁.1 hr₁.1 hc₁.1
        (hcl₁.frame rfl rfl (fun _ h => h)), ?_⟩
      refine wp_mono _ (wp_and (rgTree σ n child (d + 1) _ pre hs₁.frame' hnd hσ₁ hg₁.1 hb₁.1 hr₁.1)
          (wp_and hg₁.2 (wp_and hb₁.2 (wp_and hf₁.2 (wp_and hr₁.2 (wp_and hc₁.2
          (ccTree σ n child (d + 1) _ pre hs₁.frame' hnd hσ₁ hg₁.1 hb₁.1 hf₁.1 hr₁.1 hc₁.1
            (hcl₁.frame rfl rfl (fun _ h => h)))))))))
        fun _ s₃ ⟨⟨hi₃, hs₃, hσ₃⟩, hg₃, hb₃, hf₃, hr₃, hc₃, hcl₃⟩ => ?_
      exact ⟨by rw [if_pos (hb₃.tree.2 ⟨_, _, _, rfl⟩)], hnd, ⟨hi₃, hs₃, hσ₃⟩, hg₃, hb₃, hf₃, hr₃, hcl₃, hc₃⟩
end

/-- `CloseBase` at every `finishEdge` pre-state of the forest walk. -/
def CbForest (σ : List Nat) (n : Nat) : List DfsTree → WalkState → Prop
  | [], _ => True
  | t :: rest, s => CbTree σ n t 0 s ∧
      wp (walkTree t 0) (fun _ s₁ =>
        wp (popTstack >>= fun top => modifyItem rootItem fun it => { it with ch := it.ch ++ top.spans.2 })
          (fun _ s₂ => CbForest σ (n + t.edgePostorder.length) rest s₂) s₁) s

theorem forest_closeBase {g : Graph} {σ : List Nat} (hnd : σ.Nodup) (hσ : ∀ e ∈ σ, e < g.ne) :
    ∀ forest pre n s, RootState g pre s → s.RangesInv σ n 0 → ForestOK g (pre ++ forest) →
      (∀ t ∈ forest, t.WF []) → (∀ t ∈ forest, t.Ends g) →
      (∀ t ∈ forest, ∀ e, e < g.ne → ∀ x, g.Inc e x → x ∈ t.verts → e ∈ t.edges) →
      RootsCover σ n forest s →
      PostAt σ n (edgePostorderForest forest) → s.CloseInv →
      CbForest σ n forest s
  | [], _, _, _, _, _, _, _, _, _, _, _, _ => trivial
  | t :: rest, pre, n, s, h, hr, hf, hwf, hends, hcomp, hc, hat, hcl => by
    have hb := h.book hf (hwf t (by simp)) (hends t (by simp)) (hcomp t (by simp))
    have hg := gbTree t 0 s hb
    have hi : ∀ v outs, t = .node v outs →
        ({ s with stackVerts := s.stackVerts.set! 0 v } : WalkState).RangesInv σ n 0 :=
      fun v _ _ => hr.stackVerts_of_nil h.tstack _
    have hσ' : ∀ e ∈ σ, e < s.g.ne := by rwa [h.g_eq]
    have hfront := walkTree_frontiers t 0 s (fun v outs ht => (hi v outs ht).inv) h.shape hg hb
    have hat' : PostAt σ n (t.edgePostorder ++ edgePostorderForest rest) := hat
    have hsched := scheduleTree σ n t 0 s hi h.shape hnd hσ' hg hb hfront hc.1 hat'.left
    have hrg := rgTree σ n t 0 s hi h.shape hnd hσ' hg hb hsched
    have hcs : CsTree σ n t 0 s :=
      csTree_of_dfs (anc := []) t 0 s h.place.types rfl (hwf t (by simp)) (hends t (by simp))
        (hcomp t (by simp)) (by simpa using (RootState.hvn hf).1) (by simp) (RootState.hvlt hf) (RootState.helt hf)
        (RootState.hen hf).1 h.sv (fun k hk => absurd hk (Nat.not_lt_zero _)) hnd hat'.left
        (fun v outs _ k hk => absurd hk (by simp))
    have hcc := ccTree σ n t 0 s hi h.shape hnd hσ' hg hb hfront hsched hcs hcl
    have hcb := cbTree σ n t 0 s hi h.shape hnd hσ' hg hb hfront hsched hcs hcl
    have hnv : 0 < s.g.nv := by
      obtain ⟨v, outs⟩ := t
      have hv := RootState.hvlt hf v (by simp [DfsTree.verts])
      rw [h.g_eq]; omega
    have hk : wp (walkTree t 0) (fun _ s' => RootOK s') s :=
      walkTree_rootOK t s (fun v outs ht => (hi v outs ht).inv) h.shape hg hb h.tstack
      (by rw [h.sd, ← h.g_eq]; exact hnv)
    have hp := (walk_place_aux g).1 t 0 _ _ s h.place (RootState.hvlt hf) (RootState.helt hf)
      (RootState.hvn hf).1 (RootState.hen hf).1 (RootState.hPv hf) (RootState.hPe hf)
    have hst := h.step hf (hwf t (by simp)) (hends t (by simp)) (hcomp t (by simp))
    refine ⟨hcb, ?_⟩
    refine wp_mono _ (wp_and hrg (wp_and hk (wp_and hp (wp_and hst (wp_and hc.2 hcc)))))
      fun _ s₁ ⟨hrs, hk, hp, hst, hc, hcl⟩ => ?_
    have hrpop := hrs.1.root_append hrs.2.1 hk (noParent_of_cnt_eq_zero hp.root)
    have hp' : s₁.Place s₁.g (Pushed g (Pushed g (fun _ => False) (pre.flatMap DfsTree.verts)
        (pre.flatMap DfsTree.edges)) t.verts t.edges) (fun _ => False) := by
      rw [hp.g_eq]; exact hp
    have hclpop := rootAppend_closeInv hcl hp' hk
    exact wp_mono (popTstack >>= fun top => modifyItem rootItem fun it => { it with ch := it.ch ++ top.spans.2 })
      (wp_and hrpop (wp_and hst (wp_and hc hclpop))) fun _ s₂ ⟨hr₂, hs₂, hc₂, hcl₂⟩ =>
      forest_closeBase hnd hσ rest (pre ++ [t]) (n + t.edgePostorder.length) s₂ hs₂ hr₂
        (by simpa using hf) (fun t' ht' => hwf t' (by simp [ht']))
        (fun t' ht' => hends t' (by simp [ht'])) (fun t' ht' => hcomp t' (by simp [ht']))
        hc₂ hat'.right hcl₂

/-- `CloseBase` at every `finishEdge` pre-state of the walk of a DFS forest from the initial state. -/
theorem walk_closeBase (g : Graph) (tern : Bool) (forest : List DfsTree)
    (hf : ForestOK g forest) (hwf : ∀ t ∈ forest, t.WF []) (hends : ∀ t ∈ forest, t.Ends g)
    (hecov : ∀ e, e < g.ne → e ∈ forest.flatMap DfsTree.edges) :
    CbForest (edgePostorderForest forest) 0 forest (WalkState.init g tern) := by
  have hnd := DfsData.edgePostorderForest_perm.nodup_iff.2 hf.edges_nodup
  have hσ : ∀ e ∈ edgePostorderForest forest, e < g.ne :=
    fun e he => hf.edges_lt e (DfsData.edgePostorderForest_perm.subset he)
  exact forest_closeBase hnd hσ forest [] 0 (WalkState.init g tern) (rootState_init g tern)
    (init_rangesInv g tern hnd) (by simpa using hf) hwf hends
    (fun t ht => comp_of_forest hf hwf hends hecov ht)
    (walk_rootsCover' g tern forest hf hwf hends hecov)
    ⟨[], [], rfl, by simp⟩ (init_closeInv g tern)

end WalkState

/-- `CloseBase` at every `finishEdge` pre-state of `g.walk tern (g.dfsForest vo eo)`. -/
theorem walk_closeBase (g : Graph) (tern : Bool) (vo eo : List Nat) (hg : g.WF)
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) :
    WalkState.CbForest (edgePostorderForest (g.dfsForest vo eo)) 0 (g.dfsForest vo eo)
      (WalkState.init g tern) := by
  obtain ⟨hvp, hep⟩ := dfsForest_spanning' hg hvo heo
  exact WalkState.walk_closeBase g tern _ (ForestOK.of_perm hvp hep) (dfsForest_wf hg hvo heo)
    (dfsForest_ends g hg hvo heo) (fun e he => hep.mem_iff.2 (List.mem_range.2 he))

end Spqr
