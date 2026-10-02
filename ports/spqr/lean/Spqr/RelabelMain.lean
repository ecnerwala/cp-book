import Spqr.RelabelProof

/-!
# `relabel_node_spec`

Glue from the recursive characterization `relabel_spec` (run from `rootItem` on the initial state)
to the top-level per-node statement of `RelabelSpec.lean`.
-/

namespace Spqr.Ghost

variable {g : Graph} {items : Items}

theorem getElem!_replicate_none (n v : Nat) : (Array.replicate n (none : Option Nat))[v]! = none := by
  rw [Array.getElem!_eq_getD_getElem?, Array.getElem?_replicate]
  split <;> rfl

theorem init_order : (RelabelState.init g items).order = #[] := rfl

theorem init_consistent : Consistent g items (RelabelState.init g items) where
  g_eq := rfl
  items_eq := rfl
  order_size := rfl
  order_nodup := List.nodup_nil
  par_size := rfl
  subtreeEnd_size := rfl
  origId_size := rfl
  vertParNv_size := rfl
  chBounds_size := rfl
  nvBounds_size := rfl
  neBounds_size := rfl
  ch_last := rfl
  nv_last := rfl
  ne_last := rfl
  adjBounds_size := rfl
  adjDat_size := rfl
  vertIndex_size := Array.size_replicate
  edgeIndex_size := Array.size_replicate
  edgeFlipped_size := Array.size_replicate
  vertPos_size := Array.size_replicate
  vert_index_mem := fun v h => absurd (getElem!_replicate_none g.nv v) h
  edge_index_mem := fun e h => absurd (getElem!_replicate_none g.ne e) h

theorem relabelRun_post (hwf : items.WF g) :
    CallPost g items rootItem none none none (RelabelState.init g items) (relabelRun g items) :=
  relabel_spec hwf items.size rootItem none none none _
    ⟨init_consistent, by have := hwf.tree.size; show (0 : Nat) < items.size; omega,
      fun _ _ h => by rw [init_order] at h; simp at h⟩
    Items.card_desc_le

theorem relabel_node_spec_proved (hwf : items.WF g) :
    ∃ idx, idx rootItem = 0 ∧ RelabelIdx g items (relabelTree g items) idx ∧
      ∀ i, i < items.size → RelabelNode g items (relabelTree g items) idx i := by
  have hp := relabelRun_post hwf
  have hc := hp.cons
  have ht : relabelTree g items = (relabelRun g items).tree g := relabelTree_eq g items
  have hmem : ∀ i, i ∈ (relabelRun g items).order.toList ↔ i < items.size := fun i => by
    rw [hp.mem, init_order]
    constructor
    · rintro (h | h)
      · simp at h
      · exact (Items.mem_desc.1 h).1
    · intro h; exact Or.inr (Items.mem_desc.2 ⟨h, hwf.tree.reach i h⟩)
  have hsz : (relabelRun g items).types.size = items.size := by
    rw [← hc.order_size, ← Array.length_toList, ← List.toFinset_card_of_nodup hc.order_nodup]
    have : (relabelRun g items).order.toList.toFinset = Finset.range items.size := by
      ext i; rw [List.mem_toFinset, Finset.mem_range, hmem]
    rw [this, Finset.card_range]
  have hidx0 : (relabelRun g items).idx rootItem = 0 := hp.idx_cur
  have hnode : ∀ i, i < items.size → NodeS g items (relabelRun g items) i := fun i hi =>
    (hp.nodes i (Items.mem_desc.2 ⟨hi, hwf.tree.reach i hi⟩)).1
  have hidxlt : ∀ i, i < items.size → (relabelRun g items).idx i < items.size := fun i hi =>
    hsz ▸ hc.idx_lt ((hmem i).2 hi)
  refine ⟨(relabelRun g items).idx, hidx0, ?_, fun i hi => by rw [ht]; exact (hnode i hi).node⟩
  rw [ht]
  exact {
    size := hsz
    root := hidx0
    root_par := by rw [tree_parent_eq]; exact hp.par
    root_par_nv := hp.par_nv
    lt := hidxlt
    inj := fun i j hi _ h => (List.idxOf_inj ((hmem i).2 hi)).1 h
    vert_index := fun v hv => by
      have hn := (hnode (vertItem v) (by have := hwf.tree.size; show (1 + v : Nat) < items.size; omega)).node.vert_index
        (hwf.tree.vert v hv)
      rw [show vertItem v - 1 = v from Nat.add_sub_cancel_left 1 v] at hn
      exact hn
    edge_index := fun e he => by
      have hn := (hnode (edgeItem g e) (by have := hwf.tree.size; show (1 + g.nv + e : Nat) < items.size; omega)).node.edge_index
        (hwf.tree.edge e he)
      rw [show edgeItem g e - 1 - g.nv = e by show (1 + g.nv + e : Nat) - 1 - g.nv = e; omega] at hn
      exact hn.1
    nv := rfl
    ne := rfl
    sizes := {
      par := hc.par_size
      subtreeEnd := hc.subtreeEnd_size
      origId := hc.origId_size
      chBounds := hc.chBounds_size
      nvBounds := hc.nvBounds_size
      neBounds := hc.neBounds_size
      adjBounds := by
        show (relabelRun g items).adjBounds.size = 2 * ((relabelRun g items).nodeVerts.map _).size + 1
        rw [Array.size_map]; exact hc.adjBounds_size
      vertParNv := hc.vertParNv_size
      vertIndex := hc.vertIndex_size
      edgeIndex := hc.edgeIndex_size
      edgeFlipped := hc.edgeFlipped_size }
    ch_zero := by
      show (relabelRun g items).chBounds[0]! = 0
      rw [hp.agree.chBounds.get! (Nat.zero_le _) Nat.zero_lt_one]; rfl
    nv_zero := by
      show (relabelRun g items).nvBounds[0]! = 0
      rw [hp.agree.nvBounds.get! (Nat.zero_le _) Nat.zero_lt_one]; rfl
    ne_zero := by
      show (relabelRun g items).neBounds[0]! = 0
      rw [hp.agree.neBounds.get! (Nat.zero_le _) Nat.zero_lt_one]; rfl
    adj_zero := by
      show (relabelRun g items).adjBounds[0]! = 0
      rw [hp.agree.adjBounds.get! (Nat.zero_le _) Nat.zero_lt_one]; rfl
    ch_last := hc.ch_last
    nv_last := by
      show (relabelRun g items).nvBounds[(relabelRun g items).types.size]! =
        ((relabelRun g items).nodeVerts.map _).size
      rw [Array.size_map]; exact hc.nv_last
    ne_last := hc.ne_last
    adj_dat_size := hc.adjDat_size }

end Spqr.Ghost
