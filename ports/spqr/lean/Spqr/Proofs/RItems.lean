import Spqr.Proofs.RInvFrame
import Spqr.Proofs.ThreeConnected
import Spqr.Correctness

/-!
# R items are 3-connected (PROOF.md §4.5, item level)

`Items.RSkel3 g items i`: the skeleton of the R item `i` (its non-`V` children as pieces, the
complement of its edge set as the parent piece at `i`'s terminals, contracted) is 3-connected.
`RBranch.rSkel3` establishes it for the item Loop 1's R branch closes (`rCloseItems`: `allocItem .R`,
then `finishTstackTop` of the merged entry), from `RBranch.threeConnected`. `items_r_three_connected`
(admitted) is the statement for every R item of the walk's output on a block.
-/

namespace Spqr

theorem findIdx?_congr_mem {α : Type} {l : List α} {p q : α → Bool} (h : ∀ x ∈ l, p x = q x) :
    l.findIdx? p = l.findIdx? q := by
  induction l with
  | nil => rfl
  | cons a l ih =>
    rw [List.findIdx?_cons, List.findIdx?_cons, h a (.head _), ih fun x hx => h x (.tail _ hx)]

theorem Items.type_modify_of_ne {items : Items} (j : ItemId) (f : Item → Item) {p : ItemId}
    (h : p ≠ j) : Items.type (items.modify j f) p = items.type p := by
  simp [Items.type, Array.getElem?_modify, h.symm]

theorem Pieces.ofItems_piece_none (g : Graph) (items : Items) (L : List ItemId) {e : Nat}
    (he : ¬ e < g.ne) : (Pieces.ofItems g items L).piece e = none := by
  simp [Pieces.ofItems, he]

open Classical in
theorem Pieces.ofItems_addParent_edges {g : Graph} {items : Items} {L : List ItemId}
    {U : Nat → Prop} (s t : Nat)
    (hcover : ∀ e, e < g.ne → U e → ∃ i ∈ L, items.EdgeBelow g i e) :
    (((Pieces.ofItems g items L).addParent g U s t).contract g).edges.toList =
      (L.map fun i => ((items.vs i).1.getD 0, (items.vs i).2.getD 0)) ++ [(s, t)] := by
  let P := (Pieces.ofItems g items L).addParent g U s t
  have hnone : ∀ e, e < g.ne → P.piece e ≠ none := by
    intro e he
    dsimp [P, Pieces.addParent]
    by_cases hu : U e
    · simp only [hu, ↓reduceIte, Pieces.ofItems, he]
      intro hn
      obtain ⟨i, hi, hie⟩ := hcover e he hu
      have := List.findIdx?_eq_none_iff.1 hn i hi
      simp [hie] at this
    · simp [hu, he]
  have hempty : ((List.range g.ne).filter fun e => (P.piece e).isNone) = [] := by
    apply List.filter_eq_nil_iff.2
    intro e he
    simpa using hnone e (List.mem_range.1 he)
  change ((P.origins g).map (P.originEnds g)).toArray.toList = _
  rw [List.toList_toArray, Pieces.origins, hempty]
  simp only [List.map_nil, List.nil_append, List.map_map]
  change (List.range (L.length + 1)).map (fun i => (P.x i, P.y i)) = _
  rw [List.range_succ, List.map_append, List.map_singleton]
  congr 1
  · apply List.ext_getElem
    · simp
    · intro i hi hi'
      have hil : i < L.length := by simpa using hi'
      simp [P, Pieces.addParent, Pieces.ofItems, List.getElem?_eq_getElem hil, Nat.ne_of_lt hil]
  · simp [P, Pieces.addParent, Pieces.ofItems]

open Classical in
theorem Items.rSkeleton_perm_contract {g : Graph} {items : Items} {i : ItemId} {s t : Nat}
    (hvs : items.vs i = (some s, some t))
    (hcover : ∀ e, e < g.ne → items.EdgeBelow g i e →
      ∃ c ∈ items.ch i, items.type c ≠ .V ∧ items.EdgeBelow g c e) :
    ((((Pieces.ofItems g items ((items.ch i).filter fun c => decide (items.type c ≠ .V))).addParent g
      (items.EdgeBelow g i) s t).contract g).edges.toList.map
        (fun p => ((items.nvList g i).idxOf p.1, (items.nvList g i).idxOf p.2))).Perm
      (items.rSkeleton g i) := by
  rw [Pieces.ofItems_addParent_edges s t (fun e he hE => by
    obtain ⟨c, hc, ht, hce⟩ := hcover e he hE
    exact ⟨c, List.mem_filter.2 ⟨hc, by simpa using ht⟩, hce⟩)]
  unfold Items.rSkeleton
  rw [hvs]
  exact (List.perm_append_singleton (s, t) (items.virtualEdges i)).map _

open Classical in
theorem Pieces.ofItems_congr {g : Graph} {items items' : Items} {L : List ItemId}
    (hE : ∀ i ∈ L, ∀ e, items'.EdgeBelow g i e ↔ items.EdgeBelow g i e)
    (hvs : ∀ i : Nat, Items.vs items' L[i]! = Items.vs items L[i]!) :
    Pieces.ofItems g items' L = Pieces.ofItems g items L := by
  have hp : (fun e => if e < g.ne then L.findIdx? fun i => decide (items'.EdgeBelow g i e) else none) =
      fun e => if e < g.ne then L.findIdx? fun i => decide (items.EdgeBelow g i e) else none := by
    funext e
    split_ifs with he
    · exact findIdx?_congr_mem fun i hi => by rw [decide_eq_decide]; exact hE i hi e
    · rfl
  have hx : (fun i => (Items.vs items' L[i]!).1.getD 0) = fun i => (Items.vs items L[i]!).1.getD 0 :=
    funext fun i => by rw [hvs i]
  have hy : (fun i => (Items.vs items' L[i]!).2.getD 0) = fun i => (Items.vs items L[i]!).2.getD 0 :=
    funext fun i => by rw [hvs i]
  unfold Pieces.ofItems
  rw [hp, hx, hy]

open Classical in
theorem Pieces.addParent_congr (P : Pieces) (g : Graph) {U U' : Nat → Prop} (s t : Nat)
    (hU : ∀ e, e < g.ne → (U e ↔ U' e)) (hP : ∀ e, ¬ e < g.ne → P.piece e = none) :
    P.addParent g U s t = P.addParent g U' s t := by
  have hpiece : (fun e => if U e then P.piece e else if e < g.ne then some P.k else none) =
      fun e => if U' e then P.piece e else if e < g.ne then some P.k else none := by
    funext e
    by_cases he : e < g.ne
    · rw [hU e he]
    · simp [he, hP e he]
  unfold Pieces.addParent
  rw [hpiece]

/-- The skeleton of the R item `i` in the block `g` is 3-connected: its non-`V` children as
pieces, the complement of its edge set as the parent piece at `i`'s terminals, contracted. -/
def Items.RSkel3 (g : Graph) (items : Items) (i : ItemId) : Prop :=
  ∃ s t, (items.vs i = (some s, some t) ∨ items.vs i = (some t, some s)) ∧
    (((Pieces.ofItems g items ((items.ch i).filter fun c => decide (items.type c ≠ .V))).addParent g
      (items.EdgeBelow g i) s t).contract g).ThreeConnected

theorem Items.RSkel3.rThreeConnected {g : Graph} {items : Items} {i : ItemId}
    (h : Items.RSkel3 g items i) (hwf : items.WF g) (hi : i < items.size)
    (hR : items.type i = .R)
    (hcover : ∀ e, e < g.ne → items.EdgeBelow g i e →
      ∃ c ∈ items.ch i, items.type c ≠ .V ∧ items.EdgeBelow g c e) :
    SpqrTree.ThreeConnected (items.nvList g i).length (items.rSkeleton g i) := by
  classical
  obtain ⟨s, t, hvs, h3⟩ := h
  let P := (Pieces.ofItems g items ((items.ch i).filter fun c => decide (items.type c ≠ .V))).addParent g
    (items.EdgeBelow g i) s t
  obtain ⟨idx, _, hidx, hnode⟩ := relabel_node_spec g items hwf
  let hr : RelabelOK g items (relabelTree g items) idx := ⟨hwf, hidx, hnode⟩
  have hedges : (P.contract g).edges.toList = items.virtualEdges i ++ [(s, t)] :=
    Pieces.ofItems_addParent_edges s t (fun e he hE => by
      obtain ⟨c, hc, ht, hce⟩ := hcover e he hE
      exact ⟨c, List.mem_filter.2 ⟨hc, by simpa using ht⟩, hce⟩)
  have hcap : ∀ v, (items.vs i).1 = some v ∨ (items.vs i).2 = some v → v = s ∨ v = t := by
    intro v hv
    rcases hvs with hvs | hvs <;> simp only [hvs, Option.some.injEq] at hv <;> tauto
  have hcapmem : s ∈ items.nvList g i ∧ t ∈ items.nvList g i := by
    rcases hvs with hvs | hvs <;> simp [RelabelOK.nvList_eq, hvs]
  have hvirt : ∀ q ∈ items.virtualEdges i, q.1 ∈ items.nvList g i ∧ q.2 ∈ items.nvList g i := by
    intro q hq
    obtain ⟨c, hc, rfl⟩ := List.mem_map.1 hq
    obtain ⟨hc, hcV⟩ := List.mem_filter.1 hc
    obtain ⟨u, v, hcv⟩ := (hwf.shapes.r_shape i hi hR).2.2.2.2.2 c hc (by simpa using hcV)
    simp only [hcv, Option.getD_some]
    exact ⟨hr.mem_nvList_of_child (by simp [hR, NodeType.isNode]) hc (by simpa using hcV)
      (.inl (by rw [hcv])),
      hr.mem_nvList_of_child (by simp [hR, NodeType.isNode]) hc (by simpa using hcV)
      (.inr (by rw [hcv]))⟩
  have hends : ∀ e u v, (P.contract g).Joins e u v →
      u ∈ items.nvList g i ∧ v ∈ items.nvList g i := by
    intro e u v he
    rcases he with he | he
    all_goals
      have hm := Array.mem_of_getElem? he
      rw [← Array.mem_toList_iff, hedges] at hm
      rcases List.mem_append.1 hm with hm | hm
      · first | exact hvirt _ hm | exact (hvirt _ hm).symm
      · simp only [List.mem_singleton, Prod.mk.injEq] at hm
        rcases hm with ⟨rfl, rfl⟩
        first | exact hcapmem | exact hcapmem.symm
  have hactive : ∀ v ∈ items.nvList g i, ∃ e, (P.contract g).IsEnd e v := by
    intro v hv
    have hin : ∃ q ∈ (P.contract g).edges.toList, q.1 = v ∨ q.2 = v := by
      rw [RelabelOK.nvList_eq, hr.filter_lt_eq hi] at hv
      rcases List.mem_append.1 hv with hv | hv
      · rcases List.mem_append.1 hv with hv | hv
        · exact ⟨(s, t), by simp [hedges],
            (hcap v (.inl (by simpa using hv))).imp Eq.symm Eq.symm⟩
        · obtain ⟨c, hc, rfl⟩ := List.mem_map.1 hv
          obtain ⟨hc, hcV⟩ := List.mem_filter.1 hc
          have hclt := (hr.type_V_iff hi hc).1 (by simpa using hcV)
          obtain ⟨heq, _⟩ := hr.ch_lt_vert hc hclt
          have hvlt : c - 1 < g.nv := by
            have hc0 : c ≠ 0 := hr.ch_ne_root hc
            have hcLe : c ≤ g.nv := Nat.lt_succ_iff.1 (by simpa [Nat.add_comm] using hclt)
            omega
          have hvpar : items.IsParent i (vertItem (c - 1)) := by rwa [← heq]
          obtain ⟨⟨e, he, hev⟩, hall, hnot⟩ :=
            (hwf.endpoints.interior i (c - 1) hi hvlt (by simp [hR])).1 hvpar
          obtain ⟨j, hj, hjV, hje⟩ := hcover e he (hall e he hev)
          have hother := hnot j hj
          push Not at hother
          obtain ⟨f, hf, hfv, hjf⟩ := hother
          have hterm := hwf.endpoints.separation j (hwf.tree.ch_lt i j hj)
            (by simp [hr.child_type_ne_F hi hj, hjV])
            (c - 1) e f hvlt he hf hev hfv hje hjf
          obtain ⟨x, y, hjvs⟩ := (hwf.shapes.r_shape i hi hR).2.2.2.2.2 j hj hjV
          refine ⟨(x, y), ?_, ?_⟩
          · rw [hedges]
            apply List.mem_append_left
            exact List.mem_map.2 ⟨j, List.mem_filter.2 ⟨hj, by simpa using hjV⟩, by simp [hjvs]⟩
          · simpa [hjvs] using hterm
      · exact ⟨(s, t), by simp [hedges],
          (hcap v (.inr (by simpa using hv))).imp Eq.symm Eq.symm⟩
    obtain ⟨⟨u, w⟩, hq, hqv⟩ := hin
    obtain ⟨e, he⟩ := Graph.exists_joins_of_mem hq
    rcases hqv with rfl | rfl
    · exact ⟨e, he.isEnd⟩
    · exact ⟨e, he.symm.isEnd⟩
  have hne : 3 ≤ (P.contract g).ne := by
    have hlen := congrArg List.length hedges
    have := (hwf.shapes.r_shape i hi hR).2.1
    simp only [List.length_append, List.length_cons, List.length_nil, Array.length_toList] at hlen
    change 3 ≤ (P.contract g).edges.size
    omega
  have hcut := h3.relabel (h3.twoConnected hne) (hr.nv_nodup hi)
    (hr.nvList_R hi hR).choose_spec.choose_spec.2 hactive hends
  apply (SpqrTree.ThreeConnected_congr_undirected (es' := items.rSkeleton g i) ?_).1 hcut
  intro u v
  rw [hedges]
  rcases hvs with hvs | hvs <;>
    simp only [List.map_append, Items.rSkeleton, hvs, Option.getD_some,
      List.map_cons, List.mem_append, List.mem_cons, Prod.mk.injEq] <;> tauto

theorem Pieces.contract_congr {g : Graph} {P Q : Pieces} (hk : P.k = Q.k)
    (hp : ∀ e, e < g.ne → P.piece e = Q.piece e)
    (hxy : ∀ i, i < P.k → (P.x i, P.y i) = (Q.x i, Q.y i)) :
    P.contract g = Q.contract g := by
  have ho : P.origins g = Q.origins g := by
    unfold Pieces.origins
    rw [hk]
    congr 2
    exact List.filter_congr fun e he => by rw [hp e (List.mem_range.1 he)]
  unfold Pieces.contract
  congr 2
  rw [← ho]
  apply List.map_congr_left
  intro o ho
  cases o with
  | inl e => rfl
  | inr i =>
    apply hxy
    simpa [Pieces.origins] using ho

theorem Items.RSkel3.congr {g : Graph} {items items' : Items} {i : ItemId}
    (h : Items.RSkel3 g items i) (hch : items'.ch i = items.ch i)
    (hty : ∀ c ∈ items.ch i, items'.type c = items.type c)
    (hvs : items'.vs i = items.vs i)
    (hcv : ∀ c ∈ items.ch i, items.type c ≠ .V → items'.vs c = items.vs c)
    (hE : ∀ c ∈ items.ch i, items.type c ≠ .V →
      ∀ e, e < g.ne → (items'.EdgeBelow g c e ↔ items.EdgeBelow g c e))
    (hU : ∀ e, e < g.ne → (items'.EdgeBelow g i e ↔ items.EdgeBelow g i e)) :
    Items.RSkel3 g items' i := by
  classical
  obtain ⟨s, t, hv, hc⟩ := h
  refine ⟨s, t, by rwa [hvs], ?_⟩
  have hL : (items'.ch i).filter (fun c => decide (items'.type c ≠ .V)) =
      (items.ch i).filter (fun c => decide (items.type c ≠ .V)) := by
    rw [hch]
    exact List.filter_congr fun c hc => by rw [hty c hc]
  rw [hL]
  convert hc using 1
  refine Pieces.contract_congr ?_ ?_ ?_
  · rfl
  · intro e he
    have hf : ((items.ch i).filter fun c => decide (items.type c ≠ .V)).findIdx?
        (fun c => decide (items'.EdgeBelow g c e)) =
        ((items.ch i).filter fun c => decide (items.type c ≠ .V)).findIdx?
        (fun c => decide (items.EdgeBelow g c e)) := by
      apply findIdx?_congr_mem
      intro c hc
      rw [decide_eq_decide]
      apply hE c
      · exact (List.mem_filter.1 hc).1
      · exact of_decide_eq_true (List.mem_filter.1 hc).2
      · exact he
    dsimp only [Pieces.addParent, Pieces.ofItems]
    simp only [hU e he, hf]
  · intro k hk
    dsimp only [Pieces.addParent, Pieces.ofItems] at hk ⊢
    by_cases hki : k = ((items.ch i).filter fun c => decide (items.type c ≠ .V)).length
    · simp only [hki, ↓reduceIte]
    · have hkl : k < ((items.ch i).filter fun c => decide (items.type c ≠ .V)).length := by omega
      have hm := List.getElem_mem (l := (items.ch i).filter fun c => decide (items.type c ≠ .V)) hkl
      obtain ⟨hm, ht⟩ := List.mem_filter.1 hm
      rw [decide_eq_true_eq] at ht
      simp only [hki, ↓reduceIte, getElem!_pos ((items.ch i).filter fun c => decide (items.type c ≠ .V)) k hkl,
        hcv _ hm ht]

theorem Items.RSkel3.modify_of_not_below {g : Graph} {items : Items} {i j : ItemId}
    (h : Items.RSkel3 g items i) (hj : ¬items.Below i j) (f : Item → Item) :
    Items.RSkel3 g (items.modify j f) i := by
  have hij : i ≠ j := fun he => hj (he ▸ Relation.ReflTransGen.refl)
  have hcj : ∀ c ∈ items.ch i, c ≠ j := fun c hc he =>
    hj (Relation.ReflTransGen.single (he ▸ hc))
  refine h.congr (Items.ch_modify_of_ne j f hij)
    (fun c hc => Items.type_modify_of_ne j f (hcj c hc))
    (Items.vs_modify_of_ne j f hij)
    (fun c hc _ => Items.vs_modify_of_ne j f (hcj c hc)) ?_ ?_
  · intro c hc _ e _
    exact Items.Below_modify_of_not_below j f fun hb =>
      hj (Relation.ReflTransGen.head hc hb)
  · intro e _
    exact Items.Below_modify_of_not_below j f hj

theorem Items.RSkel3.push_nil {g : Graph} {items : Items} {i : ItemId}
    (h : Items.RSkel3 g items i) (hi : i < items.size)
    (hc : ∀ c ∈ items.ch i, c < items.size) (x : Item) (hx : x.ch = []) :
    Items.RSkel3 g (items.push x) i := by
  refine h.congr (Items.ch_push_nil x hx i)
    (fun c hm => Items.type_push_of_ne x (Nat.ne_of_lt (hc c hm)))
    (Items.vs_push_of_ne x (Nat.ne_of_lt hi))
    (fun c hm _ => Items.vs_push_of_ne x (Nat.ne_of_lt (hc c hm))) ?_ ?_
  · intro c _ _ e _
    exact Items.Below_push_nil x hx
  · intro e _
    exact Items.Below_push_nil x hx

namespace WalkState

variable {s : WalkState} {d : Nat} {cur nxt : TEntry} {rest : List TEntry} {dfs : DfsData}

/-- The items after Loop 1's R close of `cur`, `nxt` into the fresh item `s.items.size`
(`allocItem .R`, `mergeTstackTops`, then `finishTstackTop` with the merged spans on side `dir`). -/
def rCloseItems (s : WalkState) (d : Nat) (cur nxt : TEntry) (dir : Bool) : Items :=
  (s.items.push ⟨.R, (none, none), []⟩).modify s.items.size fun it =>
    { it with vs := setSides dir (some s.stackVerts[d]!) (some nxt.vStart),
              ch := getSide (TEntry.mergeInto cur nxt).spans dir }

theorem rU_iff_rItems {e : Nat} :
    s.rU cur nxt e ↔ ∃ i ∈ rItems cur nxt, Items.EdgeBelow s.g s.items i e := by
  constructor
  · rintro (⟨i, hi, hb⟩ | ⟨i, hi, hb⟩)
    · exact ⟨i, by simp only [rItems, List.mem_append] at hi ⊢; tauto, hb⟩
    · exact ⟨i, by simp only [rItems, List.mem_append] at hi ⊢; tauto, hb⟩
  · rintro ⟨i, hi, hb⟩
    exact rU_of_mem_rItems hi hb

/-- The item closed at Loop 1's R branch has a 3-connected skeleton. -/
theorem RBranch.rSkel3 (hsh : Shape s) (hb : s.RBranch d cur nxt rest) (hR : s.RTop dfs cur nxt)
    (h : s.Inv' (d + 1)) (h2 : s.g.TwoConnected) (hsp : dfs.Spec s.g) (hrt : dfs.Rooted s.g)
    (dir : Bool) (hside : getSide (TEntry.mergeInto cur nxt).spans (!dir) = []) :
    Items.RSkel3 s.g (s.rCloseItems d cur nxt dir) s.items.size := by
  unfold rCloseItems
  set n := s.items.size with hn
  set x0 : Item := ⟨.R, (none, none), []⟩ with hx0
  set f : Item → Item := fun it =>
    { it with vs := setSides dir (some s.stackVerts[d]!) (some nxt.vStart),
              ch := getSide (TEntry.mergeInto cur nxt).spans dir } with hf
  set items' := (s.items.push x0).modify n f with hitems'
  have hL : getSide (TEntry.mergeInto cur nxt).spans dir = rItems cur nxt := by
    cases dir
    · simp only [Bool.not_false, getSide, TEntry.mergeInto, Bool.false_eq_true, ↓reduceIte] at hside ⊢
      rw [rItems, hside, List.append_nil]
    · simp only [Bool.not_true, getSide, TEntry.mergeInto, Bool.false_eq_true, ↓reduceIte] at hside ⊢
      rw [rItems, hside, List.nil_append]
  have hmem : ∀ i ∈ rItems cur nxt, i < n := by
    intro i hi
    simp only [rItems, List.mem_append] at hi
    have hc := hsh.span cur (by simp [hb.tstack])
    have hn' := hsh.span nxt (by simp [hb.tstack])
    simp only [List.mem_append] at hc hn'
    rcases hi with (hi | hi) | (hi | hi)
    · exact hc i (.inl hi)
    · exact hn' i (.inl hi)
    · exact hn' i (.inr hi)
    · exact hc i (.inr hi)
  have hroot : ∀ p, ¬ Items.IsParent s.items p n := fun p hp => Nat.lt_irrefl _ (hsh.ch_lt p n hp)
  have hnlt : n < (s.items.push x0).size := by simp [hn]
  have hty : ∀ i, i < n → Items.type items' i = Items.type s.items i := fun i hi => by
    rw [hitems', Items.type_modify_of_ne n f (Nat.ne_of_lt hi),
      Items.type_push_of_ne x0 (Nat.ne_of_lt hi)]
  have hvs : ∀ i, i ≠ n → Items.vs items' i = Items.vs s.items i := fun i hi => by
    rw [hitems', Items.vs_modify_of_ne n f hi, Items.vs_push_of_ne x0 hi]
  have hE : ∀ i, i < n → ∀ e, Items.EdgeBelow s.g items' i e ↔ Items.EdgeBelow s.g s.items i e := by
    intro i hi e
    show Items.Below ((s.items.push x0).modify n f) i (edgeItem s.g e) ↔
      Items.Below s.items i (edgeItem s.g e)
    rw [Items.Below_modify_of_not_below n f, Items.Below_push_nil x0 rfl]
    intro hb
    rw [Items.Below_push_nil x0 rfl] at hb
    exact Nat.ne_of_lt hi (Items.Below.eq_of_no_parent hroot hb)
  have hvs_n : Items.vs items' n = setSides dir (some s.stackVerts[d]!) (some nxt.vStart) := by
    rw [hitems', Items.vs_modify_at n f hnlt]
  have hch_n : Items.ch items' n = rItems cur nxt := by
    rw [hitems', Items.ch_modify_at n f hnlt, ← hL]
  have hEn : ∀ e, e < s.g.ne → (Items.EdgeBelow s.g items' n e ↔ s.rU cur nxt e) := by
    intro e he
    rw [rU_iff_rItems]
    constructor
    · intro hb
      rcases Items.Below.head_cases hb with hb | ⟨c, hc, hb⟩
      · exfalso
        have h1 := hsh.size
        have h2 : (n : Nat) = 1 + s.g.nv + e := hb
        omega
      · rw [Items.IsParent, hch_n] at hc
        exact ⟨c, hc, (hE c (hmem c hc) e).1 hb⟩
    · rintro ⟨c, hc, hb⟩
      exact Relation.ReflTransGen.head (by rw [Items.IsParent, hch_n]; exact hc) ((hE c (hmem c hc) e).2 hb)
  have hlist : (Items.ch items' n).filter (fun c => decide (Items.type items' c ≠ .V)) =
      s.rPieceItems cur nxt := by
    rw [hch_n]
    exact List.filter_congr fun c hc => by rw [hty c (hmem c hc)]
  have hvs' : ∀ i : Nat, Items.vs items' (s.rPieceItems cur nxt)[i]! =
      Items.vs s.items (s.rPieceItems cur nxt)[i]! := by
    intro i
    refine hvs _ ?_
    by_cases hi : i < (s.rPieceItems cur nxt).length
    · rw [getElem!_pos (s.rPieceItems cur nxt) i hi]
      exact Nat.ne_of_lt (hmem _ (List.mem_of_mem_filter (List.getElem_mem hi)))
    · rw [getElem!_neg (s.rPieceItems cur nxt) i hi]
      have h1 := hsh.size
      show (0 : Nat) ≠ (n : Nat)
      omega
  refine ⟨nxt.vStart, s.stackVerts[d]!, ?_, ?_⟩
  · rw [hvs_n]; cases dir <;> simp [setSides]
  · rw [hlist, Pieces.ofItems_congr (L := s.rPieceItems cur nxt)
      (fun i hi e => hE i (hmem i (List.mem_of_mem_filter hi)) e) hvs',
      Pieces.addParent_congr _ _ _ _ hEn (fun e he => Pieces.ofItems_piece_none _ _ _ he)]
    exact hb.threeConnected h h2 hsp hrt hR

end WalkState

/-- Admitted (PROOF.md §4.5, item level): on a block, every R item of the walk's output has a
3-connected skeleton. Every R item is closed either at Loop 1's R branch (`RBranch.rSkel3` under
`loop1_rBranch`) or as the type-1 R close of `finishEdge` (the `isSingle = false` case, the same
argument with `cur` the type-1 entry); the skeleton is then untouched by later steps (items are
only modified when allocated or closed). -/
theorem items_r_three_connected (g : Graph) (hg : g.WF) (tern : Bool) (vo eo : List Nat)
    (hvo : OrderOK g.nv vo) (heo : OrderOK g.ne eo) (h2 : g.TwoConnected) :
    let items := (g.walk tern (g.dfsForest vo eo)).items
    ∀ i, i < items.size → Items.type items i = NodeType.R → Items.RSkel3 g items i := by
  sorry

end Spqr
