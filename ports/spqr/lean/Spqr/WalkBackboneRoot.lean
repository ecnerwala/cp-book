import Spqr.WalkBackboneStep
import Spqr.WalkTernFrame
import Spqr.StInduct
import Spqr.Proofs.RSkelRoot

/-!
# The backbone at the forest level

`RootInv` is the conjunction between the roots of `walkForest`. `RootInv.step` runs one root
through `walkTree_inv` and then the root pop / append (`popTstack`, `modifyItem rootItem`) from the
exit conjunction `WalkInvEnd` alone: the stack shape from `TreeEnd`, the frame from `KeepF`, the
separation of the new root child from the old ones from `OwnedD`/`Full`. `walk_rootInv` is the
final state; `walk_ternarize`, `walk_canonInv`, `walk_canonical`, `walk_pieceInv`,
`walk_rSkelInv'`, `walk_stItems`, `walk_rangesInv'`, `walk_closeInv''` are its projections
(PROOF.md §4.7 stage 3a).
-/

namespace Spqr
namespace WalkState
open WalkM

theorem Place.parent_unique {g : Graph} {P X : ItemId → Prop} {s : WalkState} (h : s.Place g P X)
    {p p' c : ItemId} (hp : Items.IsParent s.items p c) (hp' : Items.IsParent s.items p' c) :
    p = p' := by
  by_contra hne
  have hle := h.le c
  have h2 : ∑ x ∈ ({p, p'} : Finset ItemId), (Items.ch s.items x).count c ≤ chCount s.items c := by
    unfold chCount
    refine Finset.sum_le_sum_of_subset fun x hx => ?_
    simp only [Finset.mem_insert, Finset.mem_singleton] at hx
    rcases hx with rfl | rfl
    · exact Finset.mem_range.2 (isParent_lt_size hp)
    · exact Finset.mem_range.2 (isParent_lt_size hp')
  rw [Finset.sum_pair hne] at h2
  have hc := List.count_pos_iff.2 hp
  have hc' := List.count_pos_iff.2 hp'
  simp only [cnt] at hle
  omega

theorem Place.cnt_pos_of_parent {g : Graph} {P X : ItemId → Prop} {s : WalkState} (_ : s.Place g P X)
    {p c : ItemId} (hp : Items.IsParent s.items p c) : 0 < s.cnt c := by
  have h1 := List.count_pos_iff.2 hp
  have h2 : (Items.ch s.items p).count c ≤ chCount s.items c :=
    Finset.single_le_sum (f := fun j => (Items.ch s.items j).count c) (fun _ _ => Nat.zero_le _)
      (Finset.mem_range.2 (isParent_lt_size hp))
  simp only [cnt]; omega

theorem below_linear {items : Items}
    (hu : ∀ p p' c, items.IsParent p c → items.IsParent p' c → p = p')
    {a b x : ItemId} (hax : items.Below a x) (hbx : items.Below b x) :
    items.Below a b ∨ items.Below b a := by
  induction hax with
  | refl => exact Or.inr hbx
  | tail hax' hx'x ih =>
    rcases hbx.cases_tail with rfl | ⟨c, hbc, hcx⟩
    · exact Or.inl (Relation.ReflTransGen.tail hax' hx'x)
    · obtain rfl := hu _ _ _ hcx hx'x
      exact ih hbc

theorem Below.eq_of_ch_nil {items : Items} {a x : ItemId} (h : items.ch a = [])
    (hb : items.Below a x) : x = a := by
  rcases hb.cases_head with rfl | ⟨c, hc, -⟩
  · rfl
  · exact absurd hc (by simp [Items.IsParent, h])

theorem idxOf_mid {A B C : List Nat} (hnd : (A ++ B ++ C).Nodup) {e : Nat} (he : e ∈ B) :
    A.length ≤ (A ++ B ++ C).idxOf e ∧ (A ++ B ++ C).idxOf e < A.length + B.length := by
  have hnA : e ∉ A := fun h' =>
    (List.nodup_append.1 (List.nodup_append.1 hnd).1).2.2 e h' e he rfl
  rw [List.idxOf_append_of_mem (List.mem_append_right _ he), List.idxOf_append_of_notMem hnA]
  have := List.idxOf_lt_length_iff.2 he
  omega

theorem mem_mid {A B C : List Nat} (hnd : (A ++ B ++ C).Nodup) {e : Nat}
    (hlo : A.length ≤ (A ++ B ++ C).idxOf e)
    (hhi : (A ++ B ++ C).idxOf e < A.length + B.length) : e ∈ B := by
  have hmem : e ∈ A ++ B ++ C :=
    List.idxOf_lt_length_iff.1 (by simp only [List.length_append]; omega)
  rcases List.mem_append.1 hmem with hm | hm
  · rcases List.mem_append.1 hm with hm | hm
    · rw [List.idxOf_append_of_mem (List.mem_append_left _ hm), List.idxOf_append_of_mem hm] at hlo
      have := List.idxOf_lt_length_iff.2 hm; omega
    · exact hm
  · have hnot : e ∉ A ++ B := fun h' => (List.nodup_append.1 hnd).2.2 e h' e hm rfl
    rw [List.idxOf_append_of_notMem hnot] at hhi
    simp only [List.length_append] at hhi; omega

/-- `TreeEnd` of a root (`d = 0`, empty base): the stack is exactly the root's vertex entry. -/
theorem TreeEnd.root {y : Nat} {outs : List DfsOut} {sv : List Nat} {sd : List Bool}
    {s' : WalkState} (h : TreeEnd y 0 outs [] [] sv sd s') (hsd : 0 < s'.stackDir.size) :
    ∃ idx, s'.tstack = [⟨y, 0, idx, ([], [vertItem y])⟩] := by
  obtain ⟨done', hv', sE, dir', L', push', D₃, rfl, hC, -, hpush, hL, -, hdir⟩ := h
  have hhv : hv' = false := by
    cases hv' with
    | false => rfl
    | true =>
      obtain ⟨o, -, ho⟩ := hC.hv_ret rfl
      exact absurd ho (Nat.not_lt_zero _)
  have hp : push' = true := hpush.2 hhv
  obtain ⟨top, htop, hT⟩ := hC.top
  have htop0 : top = [] := by
    by_contra hne
    obtain ⟨o, -, ho⟩ := hT.ret hne
    exact absurd ho (Nat.not_lt_zero _)
  have hdir' := hdir hp hsd
  subst hL htop0
  refine ⟨sE.nxtEdgeIdx, ?_⟩
  simp [pushEnd, htop, hp, hdir', setSides]

theorem height_le_of_forestOK {g : Graph} {forest : List DfsTree} (hf : ForestOK g forest) :
    ∀ t ∈ forest, t.height ≤ g.nv := by
  intro t ht
  have hnd : t.verts.Nodup := by
    obtain ⟨pre, post, rfl⟩ := List.append_of_mem ht
    have h := hf.verts_nodup
    rw [List.flatMap_append, List.flatMap_cons, List.nodup_append, List.nodup_append] at h
    exact h.2.1.1
  have hsub : t.verts ⊆ List.range g.nv := fun v hv =>
    List.mem_range.2 (hf.verts_lt v (List.mem_flatMap.2 ⟨t, ht, hv⟩))
  have := (List.subperm_of_subset hnd hsub).length_le
  simp at this
  exact (DfsTree.height_le_verts_length t).trans this

/-- The conjunction between roots: the roots `pre` of `forest` are done. -/
structure RootInv (g : Graph) (tern : Bool) (forest pre : List DfsTree) (s : WalkState) : Prop where
  root : RootState g pre s
  full : s.Full g (Pushed g (fun _ => False) (pre.flatMap DfsTree.verts) (pre.flatMap DfsTree.edges))
    (fun _ => False)
  ranges : s.RangesInv (edgePostorderForest forest) (edgePostorderForest pre).length 0
  close : s.CloseInv
  canon : s.CanonInv
  tern : s.ternarize = tern
  piece : PieceFacts s
  rootch : ∀ a, Items.IsParent s.items rootItem a → ∃ w ∈ pre.flatMap DfsTree.verts, a = vertItem w
  skel : g.TwoConnected → Items.RSkelInv g s.items
  stItems : StItems g s (refBlocks g pre)

theorem rootInv_init (g : Graph) (tern : Bool) (forest : List DfsTree)
    (hnd : (edgePostorderForest forest).Nodup) : RootInv g tern forest [] (init g tern) where
  root := rootState_init g tern
  full := (init_full g tern).monoP fun _ => ⟨False.elim, fun h => by simp [Pushed] at h⟩
  ranges := init_rangesInv g tern hnd
  close := init_closeInv g tern
  canon := fun _ p c h => absurd h (by simp [Items.IsParent, init, Items.initialItems_ch])
  tern := rfl
  piece := (PieceInv.init g tern).toPieceFacts
  rootch := fun a h => absurd h (by simp [Items.IsParent, init, Items.initialItems_ch])
  skel := fun _ => init_rSkelInv g tern
  stItems := stItems_init g tern

variable {g : Graph} {tern : Bool} {forest pre rest : List DfsTree} {t : DfsTree} {s : WalkState}

/-- One root: `walkTree_inv`, then the pop / append from `WalkInvEnd`. -/
theorem RootInv.step (h : RootInv g tern forest pre s) (hforest : forest = pre ++ t :: rest)
    (hf : ForestOK g forest) (hwfA : ∀ t ∈ forest, t.WF []) (hendsA : ∀ t ∈ forest, t.Ends g)
    (hcomp : ∀ e, e < g.ne → ∀ x, g.Inc e x → x ∈ t.verts → e ∈ t.edges) (hht : t.height ≤ g.nv)
    (hgw : g.WF) (dfs : DfsData) (hsr : g.TwoConnected → dfs.Spec g ∧ dfs.Rooted g)
    (houts : ∀ t' : DfsTree, t'.Sub t → dfs.outs t'.v = t'.outs) (hd0 : dfs.depth t.v = 0) :
    wp (walkTree t 0 >>= fun _ => popTstack >>= fun top =>
        modifyItem rootItem fun it => { it with ch := it.ch ++ top.spans.2 })
      (fun _ s' => RootInv g tern forest (pre ++ [t]) s') s := by
  subst hforest
  have hwf := hwfA t (by simp)
  have hends := hendsA t (by simp)
  have hvlt := RootState.hvlt hf
  have helt := RootState.helt hf
  have hvn := RootState.hvn hf
  have hen := RootState.hen hf
  have hPv := RootState.hPv hf
  have hPe := RootState.hPe hf
  have hge := h.root.g_eq
  have htv : t.v ∈ t.verts := by obtain ⟨v, outs⟩ := t; exact List.mem_cons_self ..
  have hvt0 : t.v < g.nv := hvlt _ htv
  have hnv : 0 < g.nv := Nat.lt_of_le_of_lt (Nat.zero_le _) hvt0
  have hσ : (edgePostorderForest (pre ++ t :: rest)).Nodup :=
    DfsData.edgePostorderForest_perm.nodup_iff.2 hf.edges_nodup
  have hσsplit : edgePostorderForest (pre ++ t :: rest) =
      edgePostorderForest pre ++ t.edgePostorder ++ edgePostorderForest rest := by
    simp [edgePostorderForest, List.flatMap_append, List.flatMap_cons, List.append_assoc]
  have hσ' : (edgePostorderForest pre ++ t.edgePostorder ++ edgePostorderForest rest).Nodup :=
    hσsplit ▸ hσ
  have hvfresh : ∀ v ∈ t.verts, Items.ch s.items (vertItem v) = [] ∧
      ∀ p, ¬ Items.IsParent s.items p (vertItem v) := fun v hv =>
    ⟨h.root.fresh _ (by show 0 < 1 + v; omega)
      (by have := hvlt v hv; show 1 + v < _; omega)
      (fun w hw hvw => hvn.2 v hv (vertItem_inj hvw ▸ hw))
      (fun e _ hve => vertItem_ne_edgeItem (hvlt v hv) e hve),
     noParent_of_cnt_eq_zero (h.root.place.cnt_eq_zero (by show 0 < 1 + v; omega)
      (by have := hvlt v hv; show 1 + v < _; omega) (hPv v hv))⟩
  have hefresh : ∀ e ∈ t.edges, Items.ch s.items (edgeItem s.g e) = [] ∧
      ∀ p, ¬ Items.IsParent s.items p (edgeItem s.g e) := fun e he => by
    rw [hge]
    exact ⟨h.root.fresh _ (by show 0 < 1 + g.nv + e; omega)
      (by have := helt e he; show 1 + g.nv + e < _; omega)
      (fun w hw hew => vertItem_ne_edgeItem (hf.verts_lt w (by simp [hw])) e hew.symm)
      (fun e' he' hee => hen.2 e he (edgeItem_inj hee ▸ he')),
     noParent_of_cnt_eq_zero (h.root.place.cnt_eq_zero (by show 0 < 1 + g.nv + e; omega)
      (by have := helt e he; show 1 + g.nv + e < _; omega) (hPe e he))⟩
  obtain ⟨sv, sd, hE⟩ := ctx_init_root t s hwf (by rw [hge]; exact hends) (by rw [hge]; exact hvlt)
    (by rw [hge]; exact helt) hvn.1 hen.1 (by rw [hge]; exact h.root.sv) (by rw [hge]; exact h.root.sd)
    (by rw [hge]; exact h.root.fo) h.root.tstack h.root.inv
    h.root.shape hvfresh hefresh
  have hnb : ∀ a, Items.IsParent s.items rootItem a → ¬ Items.Below s.items a (vertItem t.v) := by
    intro a ha hb
    obtain ⟨w, hw, rfl⟩ := h.rootch a ha
    have := (Below_of_root (hvfresh t.v htv).2).1 hb
    exact hvn.2 t.v htv (vertItem_inj this ▸ hw)
  have hsv0 : (s.stackVerts.set! 0 t.v)[0]! = t.v :=
    getElem!_set!_self' s.stackVerts 0 t.v (by rw [h.root.sv]; exact hnv)
  have hnoE : ∀ e, ¬ Items.EdgeBelow s.g s.items (vertItem t.v) e := fun e hb =>
    vertItem_ne_edgeItem (by rw [hge]; exact hvt0) e (Below.eq_of_ch_nil (hvfresh t.v htv).1 hb).symm
  have hepf : edgePostorderForest (pre ++ t :: rest) =
      edgePostorderForest pre ++ edgePostorderForest (t :: rest) := by
    simp [edgePostorderForest, List.flatMap_append]
  have hW : WalkInv
      { g := g, anc := [], base := [], bE := [], sv := sv, sd := sd, pe := (fun _ => False),
        σ := edgePostorderForest (pre ++ t :: rest), n := (edgePostorderForest pre).length,
        dfs := dfs, F := [], sts := [], origs := [],
        P := Pushed g (fun _ => False) (pre.flatMap DfsTree.verts) (pre.flatMap DfsTree.edges),
        X := (fun _ => False), prev := pre, fs := [], segs := [], tern := tern } t 0 s := {
    full := h.full
    d_anc := rfl
    wf := hwf
    ends := hends
    nodup := by simpa using hvn.1
    anc_lt := fun a ha => by simp at ha
    v_lt := hvlt
    e_lt := helt
    enodup := hen.1
    comp := fun e he x hx hxt => Or.inl (hcomp e he x hx hxt)
    pe_anc := fun _ h' => h'.elim
    sv_size := h.root.sv
    sd_size := h.root.sd
    anc_sv := fun k hk => absurd hk (Nat.not_lt_zero _)
    ear := hE
    inv := h.root.inv.stackVerts_of_nil h.root.tstack _
    shape := h.root.shape
    ranges := h.ranges.stackVerts_of_nil h.root.tstack _
    σ_nodup := hσ
    σ_lt := fun e he => by
      rw [hge]; exact hf.edges_lt e (DfsData.edgePostorderForest_perm.subset he)
    post := ⟨edgePostorderForest pre, edgePostorderForest rest, rfl, hσsplit⟩
    path := fun k hk => by simp at hk
    close := h.close
    canon := h.canon
    tern := h.tern
    piece := h.piece.root t.v (by rw [h.root.sv]; exact hnv) hnb
    P_past := fun e he hP => by
      rcases hP with hP | ⟨w, hw, hew⟩ | ⟨e', he', hee⟩
      · exact hP.elim
      · exact absurd hew.symm (vertItem_ne_edgeItem (hf.verts_lt w (by simp [hw])) e)
      · obtain rfl := edgeItem_inj hee
        have hmem : e ∈ edgePostorderForest pre := DfsData.edgePostorderForest_perm.mem_iff.2 he'
        show (edgePostorderForest (pre ++ t :: rest)).idxOf e < (edgePostorderForest pre).length
        rw [hepf, List.idxOf_append_of_mem hmem]
        exact List.idxOf_lt_length_iff.2 hmem
    owned := {
      len := fun k hk => by simp at hk; subst hk; simp
      cover := fun b hb hb' => by simp at hb hb'; omega
      anc := fun k hk => absurd hk (Nat.not_lt_zero _)
      hi := fun e he hb => (hnoE e (by rw [← hsv0]; exact hb)).elim
      lo := fun k hk e he hb => by
        simp at hk; subst hk
        exact (hnoE e (by rw [← hsv0]; exact hb)).elim
      new := fun k hk t' ht' => by simp at hk; subst hk; simp [h.root.tstack] at ht'
      old := fun k hk t' ht' => by simp at hk; subst hk; simp [h.root.tstack] at ht'
      vis := fun t' ht' => by simp [h.root.tstack] at ht'
      fresh := fun w hw hP _ => h.root.fresh _ (by show 0 < 1 + w; omega)
        (by rw [hge] at hw; show 1 + w < _; omega)
        (fun v hv hvw => hP (Or.inr (Or.inl ⟨v, hv, hvw⟩)))
        (fun e _ hve => vertItem_ne_edgeItem (hge ▸ hw) e hve) }
    sts_len := rfl
    origs_len := rfl
    sts_le := fun k hk => absurd hk (Nat.not_lt_zero _)
    origs_le := fun k hk => absurd hk (Nat.not_lt_zero _)
    Pv := hPv
    Pe := hPe
    r := fun h2 => by
      rw [hge] at h2
      exact { wf := by rw [hge]; exact hgw
              spec := by rw [hge]; exact (hsr h2).1
              rooted := by rw [hge]; exact (hsr h2).2
              outs := houts
              par := fun dp hdp => absurd hdp (by omega)
              root := fun _ => ⟨fun v outs htv => by subst htv; exact hd0, h.root.tstack, rfl⟩
              frames := fun f hf => by simp at hf
              skel := by rw [hge]; exact h.skel h2 }
    d_fs := rfl
    height := by rw [h.root.sd]; simpa using hht
    tstack := by rw [h.root.tstack]; rfl
    segRead := fun sg hsg => by simp at hsg
    base_out := fun t' ht' => absurd ht' (by simp [segsStack])
    stItems := by simpa [simBlocks, frameBlocks] using h.stItems
    segs_len := rfl
    live_lower := fun j hj => absurd hj (Nat.not_lt_zero _) }
  rw [wp_bind]
  refine wp_mono _ (walkTree_inv hW) fun _ s₁ hE₁ => ?_
  have hK0 := hE₁.keep 0 (FreshItem.zero _ _ _)
  have hg₁ : s₁.g = g := hK0.g.trans hge
  have hsd₁ : 0 < s₁.stackDir.size := by rw [hK0.sd, h.root.sd]; exact hnv
  obtain ⟨idx, hts₁⟩ := TreeEnd.root hE₁.treeEnd hsd₁
  have hv₁ : t.v < s₁.g.nv := by rw [hg₁]; exact hvt0
  have hroot₁ : ∀ p, ¬ Items.IsParent s₁.items p rootItem :=
    noParent_of_cnt_eq_zero hE₁.full.place.root
  have hrk : RootOK s₁ := ⟨_, hts₁, rfl, fun c hc => ⟨t.v, hv₁, by simpa using hc⟩⟩
  have hs₁ := hE₁.shape
  have hrootS : Items.type s₁.items rootItem ≠ .S := by rw [hs₁.root]; decide
  have hrootP : Items.type s₁.items rootItem ≠ .P := by rw [hs₁.root]; decide
  have hrootR : Items.type s₁.items rootItem ≠ .R := by rw [hs₁.root]; decide
  have hrsz : rootItem < s₁.items.size := by have := hs₁.size; show (0 : Nat) < _; omega
  have hchroot : Items.ch s₁.items rootItem = Items.ch s.items rootItem := hK0.ch
  have hnpv : ∀ p, ¬ Items.IsParent s₁.items p (vertItem t.v) :=
    hE₁.full.place.noParent_of_spans (by simp [spansCount, hts₁])
  have hu : ∀ p p' c, Items.IsParent s₁.items p c → Items.IsParent s₁.items p' c → p = p' :=
    fun _ _ _ hp hp' => hE₁.full.place.parent_unique hp hp'
  have hown := hE₁.owned
  dsimp only at hown
  rw [hσsplit] at hown
  have hbelow_t : ∀ e, e < g.ne → Items.EdgeBelow g s₁.items (vertItem t.v) e → e ∈ t.edges := by
    intro e he hb
    have he₁ : e < s₁.g.ne := by rw [hg₁]; exact he
    have hb' : Items.EdgeBelow s₁.g s₁.items (vertItem s₁.stackVerts[0]!) e := by
      rw [hg₁, hE₁.sv_v]; exact hb
    have hhi := hown.hi e he₁ hb'
    have hlo := hown.lo 0 (le_refl _) e he₁ hb'
    simp only [List.nil_append, List.getElem!_cons_zero] at hlo
    exact t.edgePostorder_perm_edges.subset (mem_mid hσ' hlo hhi)
  have hin_t : ∀ e ∈ t.edges, Items.EdgeBelow g s₁.items (vertItem t.v) e := by
    intro e het
    have hmemB := t.edgePostorder_perm_edges.mem_iff.2 het
    obtain ⟨hlo, hhi⟩ := idxOf_mid hσ' hmemB
    have hmem : e ∈ edgePostorderForest pre ++ t.edgePostorder ++ edgePostorderForest rest := by
      simp [hmemB]
    have hlt := List.idxOf_lt_length_iff.2 hmem
    have hget : (edgePostorderForest pre ++ t.edgePostorder ++ edgePostorderForest rest)[
        (edgePostorderForest pre ++ t.edgePostorder ++ edgePostorderForest rest).idxOf e]! = e := by
      rw [getElem!_of_getElem? (List.getElem?_eq_getElem hlt)]
      exact List.getElem_idxOf hlt
    have hc := hown.cover _ (by simpa using hlo) hhi
    rw [hget] at hc
    rcases hc with ⟨t', ht', hte⟩ | ⟨k, hk, hb⟩
    · rw [hts₁, List.mem_singleton] at ht'; subst ht'
      obtain ⟨i, hi, hb⟩ := hte
      simp at hi; subst hi
      rw [hg₁] at hb; exact hb
    · have : k = 0 := by omega
      subst this
      rw [hg₁, hE₁.sv_v] at hb; exact hb
  have hsep : ∀ a, (Items.IsParent s₁.items rootItem a ∨ a ∈ [vertItem t.v]) →
      ∀ b ∈ [vertItem t.v], a ≠ b →
      ∀ v e e', e < s₁.g.ne → e' < s₁.g.ne → s₁.g.Inc e v → s₁.g.Inc e' v →
      Items.EdgeBelow s₁.g s₁.items a e → Items.EdgeBelow s₁.g s₁.items b e' → False := by
    intro a ha b hb hab v e e' he he' hev he'v hae hbe
    rw [List.mem_singleton] at hb; subst hb
    rw [hg₁] at he he' hev he'v hae hbe
    rcases ha with ha | ha
    · obtain ⟨w, hw, rfl⟩ := h.rootch a (by show a ∈ _; rw [← hchroot]; exact ha)
      have he't : e' ∈ t.edges := hbelow_t e' he' hbe
      have hvt : v ∈ t.verts := tree_inc_verts hwf hends he't he'v
      by_cases het : e ∈ t.edges
      · have hbt := hin_t e het
        rcases below_linear hu hae hbt with hwb | hbw
        · exact hab ((Below_of_root hnpv).1 hwb)
        · rcases hbw.cases_tail with heq | ⟨p, hbp, hp⟩
          · exact hab heq
          · obtain rfl := hu _ _ _ hp ha
            exact absurd ((Below_of_root hroot₁).1 hbp) (vertItem_ne_zero _)
      · have hpar : ∃ p, Items.IsParent s₁.items p (edgeItem g e) := by
          rcases hae.cases_tail with heq | ⟨p, -, hp⟩
          · exact absurd heq.symm (vertItem_ne_edgeItem (hf.verts_lt w (by simp [hw])) e)
          · exact ⟨p, hp⟩
        obtain ⟨p, hp⟩ := hpar
        have hP := hE₁.full.place.fixed _ (by show 0 < 1 + g.nv + e; omega)
          (by show 1 + g.nv + e < 1 + g.nv + g.ne; omega) (hE₁.full.place.cnt_pos_of_parent hp)
        rcases hP with hP | ⟨w', hw', hw''⟩ | ⟨e₀, he₀, hee⟩
        · rcases hP with hP | ⟨w', hw', hw''⟩ | ⟨e₀, he₀, hee⟩
          · exact hP
          · exact vertItem_ne_edgeItem (hf.verts_lt w' (by simp [hw'])) e hw''.symm
          · obtain rfl := edgeItem_inj hee
            obtain ⟨t₀, ht₀, he₀'⟩ := List.mem_flatMap.1 he₀
            have hvt₀ := tree_inc_verts (hwfA t₀ (by simp [ht₀])) (hendsA t₀ (by simp [ht₀])) he₀' hev
            exact hvn.2 v hvt (List.mem_flatMap.2 ⟨t₀, ht₀, hvt₀⟩)
        · exact vertItem_ne_edgeItem (hvlt w' hw') e hw''.symm
        · exact het (edgeItem_inj hee ▸ he₀)
    · rw [List.mem_singleton] at ha; exact hab ha
  have hpop : ({ s₁ with tstack := [] } : WalkState).Inv' 0 :=
    ⟨fun above t below h' => (List.append_ne_nil_of_right_ne_nil above (List.cons_ne_nil t below) h'.symm).elim,
     fun i h1 h2 => ⟨(hE₁.inv.nodes i h1 h2).conn, (hE₁.inv.nodes i h1 h2).attached⟩⟩
  have hinv₂ := hpop.modifyCh rootItem (fun it => { it with ch := it.ch ++ [vertItem t.v] })
    (Nat.lt_of_lt_of_le (Nat.lt_of_lt_of_le Nat.zero_lt_one (Nat.le_add_right 1 _))
      (Nat.le_add_right _ _)) hroot₁ (by simp)
  have hshape₂ : Shape { s₁ with
      items := s₁.items.modify rootItem fun it => { it with ch := it.ch ++ [vertItem t.v] },
      tstack := [] } := by
    refine (hs₁.tstack (l := []) (fun e he => by simp at he)).modify rootItem _
      (fun _ => rfl) fun hj c hc => ?_
    simp only [List.mem_append, List.mem_singleton] at hc
    rcases hc with hc | rfl
    · exact hs₁.ch_lt rootItem c (by simpa [Items.ch, Items.IsParent, hj] using hc)
    · exact hs₁.span ⟨t.v, 0, idx, ([], [vertItem t.v])⟩ (by simp [hts₁]) (vertItem t.v) (by simp)
  have hF₂ := hE₁.full.root_append (by simp [hts₁]) (by
    simp only [hts₁, List.head!_cons]
    intro c hc
    exact ⟨t.v, hvt0, by simpa using hc⟩)
  simp only [hts₁, List.head!_cons, List.tail_cons] at hF₂
  have hFull₂ : Full g (Pushed g (fun _ => False) ((pre ++ [t]).flatMap DfsTree.verts)
      ((pre ++ [t]).flatMap DfsTree.edges)) (fun _ => False) _ :=
    hF₂.monoP fun i => by rw [pushed_append_iff']; simp [List.flatMap_append]
  have hr₂ := hE₁.ranges.1.root_append hs₁ hrk hroot₁
  have hp₁ := hE₁.full.place
  rw [← hg₁] at hp₁
  have hc₂ := rootAppend_closeInv hE₁.close hp₁ hrk
  simp only [wp_bind, wp_popTstack, wp_modifyItem, hts₁, List.head!_cons, List.tail_cons] at hr₂ hc₂
  have hlen : (edgePostorderForest (pre ++ [t])).length =
      (edgePostorderForest pre).length + t.edgePostorder.length := by
    simp [edgePostorderForest, List.flatMap_append]
  obtain ⟨new, hts₁', hR₁⟩ := hE₁.st.read
  simp only [segsStack, List.map_nil, List.flatten_nil, List.append_nil, List.length_nil,
    DirsOf_zero] at hts₁' hR₁
  rw [hts₁] at hts₁'
  subst hts₁'
  have hI₁ := hE₁.st.items
  simp only [simBlocks, frameBlocks, List.append_nil, List.length_nil, DirsOf_zero] at hI₁
  have hI₂ := rootPop_st hts₁ hg₁ hE₁.inv hs₁ hrk hroot₁ hR₁ hI₁
  have hcheq : Items.ch s₁.items rootItem = s₁.items[rootItem].ch := by
    simp [Items.ch, Array.getElem?_eq_getElem hrsz]
  simp only [wp_bind, wp_popTstack, wp_modifyItem, hts₁, List.head!_cons, List.tail_cons]
  exact {
    root := {
      place := hFull₂.place
      g_eq := hg₁
      sv := hK0.sv.trans h.root.sv
      sd := hK0.sd.trans h.root.sd
      fo := hK0.fo.trans h.root.fo
      shape := hshape₂
      tstack := rfl
      inv := hinv₂
      fresh := fun i hi0 hi hv he => by
        simp only [List.flatMap_append, List.flatMap_cons, List.flatMap_nil, List.append_nil] at hv he
        show Items.ch (s₁.items.modify rootItem _) i = []
        rw [Items.ch_modify_of_ne rootItem _ (Nat.ne_of_gt hi0)]
        rw [(hE₁.keep i ⟨hi, fun v hv' => hv v (List.mem_append_right _ hv'),
          fun e he' => he e (List.mem_append_right _ he')⟩).ch]
        exact h.root.fresh i hi0 hi (fun v hv' => hv v (List.mem_append_left _ hv'))
          fun e he' => he e (List.mem_append_left _ he') }
    full := hFull₂
    ranges := by rw [hlen]; exact hr₂
    close := hc₂
    canon := fun ht => (hE₁.canon ht).modify_ch_nonSP rootItem _ hrootS hrootP
    tern := hE₁.tern
    piece := (hE₁.piece.toPieceFacts.rootAppend hs₁ [vertItem t.v] hroot₁
      (fun c hc => by rw [List.mem_singleton] at hc; subst hc; exact hs₁.vert _ hv₁) hsep).frame rfl rfl
    rootch := fun a ha => by
      have ha' : a ∈ Items.ch s₁.items rootItem ++ [vertItem t.v] := by
        change a ∈ Items.ch (s₁.items.modify rootItem _) rootItem at ha
        simp only [Items.ch_modify_at rootItem _ hrsz] at ha
        rw [hcheq]
        exact ha
      simp only [List.mem_append, List.mem_singleton] at ha'
      rcases ha' with ha' | rfl
      · obtain ⟨w, hw, rfl⟩ := h.rootch a (by show a ∈ _; rw [← hchroot]; exact ha')
        exact ⟨w, by simp [hw], rfl⟩
      · exact ⟨t.v, by simp [htv], rfl⟩
    skel := fun h2 => by
      have := (hE₁.skel (by rw [hge]; exact h2)).modify_root_of_ne hroot₁
        (fun it => { it with ch := it.ch ++ [vertItem t.v] }) (fun _ => rfl) hrootR
      rw [hg₁] at this; exact this
    stItems := by rw [refBlocks_snoc, ← List.append_assoc]; exact hI₂ }

theorem forest_inv (g : Graph) (tern : Bool) (forest : List DfsTree) (hf : ForestOK g forest)
    (hwf : ∀ t ∈ forest, t.WF []) (hends : ∀ t ∈ forest, t.Ends g)
    (hcomp : ∀ t ∈ forest, ∀ e, e < g.ne → ∀ x, g.Inc e x → x ∈ t.verts → e ∈ t.edges)
    (hht : ∀ t ∈ forest, t.height ≤ g.nv) (hgw : g.WF) (dfs : DfsData)
    (hsr : g.TwoConnected → dfs.Spec g ∧ dfs.Rooted g)
    (houts : ∀ t ∈ forest, ∀ t' : DfsTree, t'.Sub t → dfs.outs t'.v = t'.outs)
    (hd0 : ∀ t ∈ forest, dfs.depth t.v = 0) :
    ∀ (rest pre : List DfsTree) (s : WalkState), forest = pre ++ rest →
      RootInv g tern forest pre s →
      wp (walkForest rest) (fun _ s' => RootInv g tern forest (pre ++ rest) s') s
  | [], pre, s, _, h => by
    show RootInv g tern forest (pre ++ []) s
    simpa using h
  | t :: rest, pre, s, hforest, h => by
    show wp ((walkTree t 0 >>= fun _ => popTstack >>= fun top =>
      modifyItem rootItem fun it => { it with ch := it.ch ++ top.spans.2 }) >>= fun _ =>
      walkForest rest) _ s
    rw [wp_bind]
    have ht : t ∈ forest := by rw [hforest]; simp
    refine wp_mono _ (h.step hforest hf hwf hends (hcomp t ht) (hht t ht) hgw dfs hsr (houts t ht)
      (hd0 t ht)) fun _ s₁ h₁ => ?_
    rw [List.append_cons]
    exact forest_inv g tern forest hf hwf hends hcomp hht hgw dfs hsr houts hd0 rest (pre ++ [t]) s₁
      (by rw [hforest, List.append_cons]) h₁

/-- The backbone at the final state of the walk. -/
theorem walk_rootInv (g : Graph) (hg : g.WF) (tern : Bool) (vo eo : List Nat)
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) :
    RootInv g tern (g.dfsForest vo eo) (g.dfsForest vo eo) (g.walk tern (g.dfsForest vo eo)) := by
  obtain ⟨hvp, hep⟩ := dfsForest_spanning' hg hvo heo
  have hf : ForestOK g (g.dfsForest vo eo) := ForestOK.of_perm hvp hep
  have hwf := dfsForest_wf hg hvo heo
  have hends := dfsForest_ends g hg hvo heo
  have hecov : ∀ e, e < g.ne → e ∈ (g.dfsForest vo eo).flatMap DfsTree.edges :=
    fun e he => hep.mem_iff.2 (List.mem_range.2 he)
  have hnd : ((g.dfsForest vo eo).flatMap DfsTree.verts).Nodup := hvp.nodup_iff.2 List.nodup_range
  have hσ : (edgePostorderForest (g.dfsForest vo eo)).Nodup :=
    DfsData.edgePostorderForest_perm.nodup_iff.2 hf.edges_nodup
  have hspec := dfsForestSpec_of_dfsForest hg hvo heo
  have key : ∀ dfs : DfsData, (g.TwoConnected → dfs.Spec g ∧ dfs.Rooted g) →
      (∀ t ∈ g.dfsForest vo eo, ∀ t' : DfsTree, t'.Sub t → dfs.outs t'.v = t'.outs) →
      (∀ t ∈ g.dfsForest vo eo, dfs.depth t.v = 0) →
      RootInv g tern (g.dfsForest vo eo) (g.dfsForest vo eo) (g.walk tern (g.dfsForest vo eo)) :=
    fun dfs hsr houts hd0 =>
      forest_inv g tern _ hf hwf hends (fun t ht => comp_of_forest hf hwf hends hecov ht)
        (height_le_of_forestOK hf) hg dfs hsr houts hd0 _ [] _ rfl (rootInv_init g tern _ hσ)
  by_cases h2 : g.TwoConnected
  · obtain ⟨r, hr0, hrt⟩ := dfsForest_rooted g hg vo eo hvo heo h2
    exact key { DfsData.ofForest (g.dfsForest vo eo) with root := r }
      (fun _ => ⟨{ hspec.toSpec with depth_root := hr0 }, hrt⟩)
      (fun t ht t' hsub => DfsData.ofForest_outs hnd ht hsub)
      (fun t ht => DfsData.ofForest_depth_root hnd ht)
  · exact key (DfsData.ofForest _) (fun h => absurd h h2)
      (fun t ht t' hsub => DfsData.ofForest_outs hnd ht hsub)
      (fun t ht => DfsData.ofForest_depth_root hnd ht)

end WalkState

open WalkState in
/-- The walk never writes `ternarize`. -/
theorem walk_ternarize (g : Graph) (tern : Bool) (forest : List DfsTree) :
    (g.walk tern forest).ternarize = tern :=
  tern_walkForest forest (WalkState.init g tern)

theorem walk_canonInv (g : Graph) (hg : g.WF) (tern : Bool) (vo eo : List Nat)
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) :
    (g.walk tern (g.dfsForest vo eo)).CanonInv :=
  (WalkState.walk_rootInv g hg tern vo eo hvo heo).canon

/-- Canonicity of the unternarized walk (`PROOF.md` §4.6; for `tern = true` P under P / S under S
do occur, `check_ranges` seeds 0 and 386). -/
theorem walk_canonical (g : Graph) (hg : g.WF) (vo eo : List Nat) (hvo : OrderOK g.nv vo)
    (heo : OrderOK g.ne eo) :
    Items.Canonical (g.walk false (g.dfsForest vo eo)).items :=
  walk_canonInv g hg false vo eo hvo heo (walk_ternarize g false _)

theorem walk_pieceFacts (g : Graph) (hg : g.WF) (tern : Bool) (vo eo : List Nat)
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) :
    (g.walk tern (g.dfsForest vo eo)).PieceFacts :=
  (WalkState.walk_rootInv g hg tern vo eo hvo heo).piece

theorem walk_rSkelInv' (g : Graph) (hg : g.WF) (tern : Bool) (vo eo : List Nat)
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) (h2 : g.TwoConnected) :
    Items.RSkelInv g (g.walk tern (g.dfsForest vo eo)).items :=
  (WalkState.walk_rootInv g hg tern vo eo hvo heo).skel h2

theorem walk_stItems (g : Graph) (hg : g.WF) (tern : Bool) (vo eo : List Nat)
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) :
    StItems g (g.walk tern (g.dfsForest vo eo)) (refBlocks g (g.dfsForest vo eo)) :=
  (WalkState.walk_rootInv g hg tern vo eo hvo heo).stItems

theorem walk_rangesInv' (g : Graph) (hg : g.WF) (tern : Bool) (vo eo : List Nat)
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) :
    (g.walk tern (g.dfsForest vo eo)).RangesInv (edgePostorderForest (g.dfsForest vo eo))
      (edgePostorderForest (g.dfsForest vo eo)).length 0 :=
  (WalkState.walk_rootInv g hg tern vo eo hvo heo).ranges

theorem walk_closeInv'' (g : Graph) (hg : g.WF) (tern : Bool) (vo eo : List Nat)
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) :
    (g.walk tern (g.dfsForest vo eo)).CloseInv :=
  (WalkState.walk_rootInv g hg tern vo eo hvo heo).close

end Spqr
