import Spqr.StLoop1

/-!
# `closeEars` keeps the st-reading

`feS₁ d o s` (set the tree edge's `vs`, push its `Q` item, run loop 1) for a tree edge with
`lowval < d`: the child's entries `sub` reading as `ps` become entries reading as
`ps ++ [⟨stackDir[d], [Q e]⟩]` above the unchanged `base`.
-/

namespace Spqr
open WalkM WalkState

theorem ExpandsList.leaves {items : Items} (L : List ItemId)
    (h : ∀ x ∈ L, Items.type items x = .V ∨ Items.type items x = .Q) : ExpandsList items L L := by
  induction L with
  | nil => exact .nil
  | cons x xs ih => exact .leaf (h x List.mem_cons_self) (ih fun y hy => h y (List.mem_cons_of_mem _ hy))

theorem mem_readStack_exists {ts : List TEntry} {x : ItemId} (h : x ∈ readStack ts) :
    ∃ u ∈ ts, x ∈ u.spans.1 ++ u.spans.2 := by
  induction ts with
  | nil => simp [readStack, readL, readR] at h
  | cons t ts ih =>
    simp only [readStack, readL, readR, List.mem_append] at h
    rcases h with (h | h) | (h | h)
    · exact ⟨t, List.mem_cons_self, List.mem_append_left _ h⟩
    · obtain ⟨u, hu, hx⟩ := ih (List.mem_append_left _ h); exact ⟨u, List.mem_cons_of_mem _ hu, hx⟩
    · obtain ⟨u, hu, hx⟩ := ih (List.mem_append_right _ h); exact ⟨u, List.mem_cons_of_mem _ hu, hx⟩
    · exact ⟨t, List.mem_cons_self, List.mem_append_right _ h⟩

theorem mem_readStack_push (v d idx : Nat) (dir : Bool) (q : ItemId) (ts : List TEntry) (x : ItemId) :
    x ∈ readStack (⟨v, d, idx, setSides dir [q] []⟩ :: ts) ↔ x = q ∨ x ∈ readStack ts := by
  rw [readStack_cons_setSides]
  cases dir <;> simp [stNestL, stNestR] <;> tauto

theorem nodup_readStack_push (v d idx : Nat) (dir : Bool) (q : ItemId) (ts : List TEntry)
    (hq : q ∉ readStack ts) (hnd : (readStack ts).Nodup) :
    (readStack (⟨v, d, idx, setSides dir [q] []⟩ :: ts)).Nodup := by
  rw [readStack_cons_setSides]
  cases dir <;> simp [stNestL, stNestR, List.nodup_append, hq, hnd] <;>
    exact fun a ha e => hq (e ▸ ha)

theorem StRead.pushEntry {items : Items} (v d idx : Nat) (dir : Bool) (q : ItemId) {new : List TEntry}
    {ps : List StPiece} (hq : Items.type items q = .V ∨ Items.type items q = .Q)
    (h : StRead items new ps) :
    StRead items (⟨v, d, idx, setSides dir [q] []⟩ :: new) (ps ++ [⟨dir, [q]⟩]) := by
  unfold StRead at h ⊢
  rw [readL_cons_setSides, readR_cons_setSides, stNestL_append, stNestR_append]
  have hl : ExpandsList items (stNestL [⟨dir, [q]⟩]) (stNestL [⟨dir, [q]⟩]) := by
    cases dir
    · exact .leaf hq .nil
    · exact .nil
  have hr : ExpandsList items (stNestR [⟨dir, [q]⟩]) (stNestR [⟨dir, [q]⟩]) := by
    cases dir
    · exact .nil
    · exact .leaf hq .nil
  exact ⟨hl.append h.1, h.2.append hr⟩

theorem StItems.modify_root {g : Graph} {s : WalkState} {blocks : List StBlock} (j : ItemId)
    (f : Item → Item) (hf : ∀ it, (f it).type = it.type) (hfc : ∀ it, (f it).ch = it.ch)
    (hr : ∀ p, ¬ Items.IsParent s.items p j)
    (hj : ¬ (Items.type s.items j = .S ∨ Items.type s.items j = .P ∨ Items.type s.items j = .R))
    (h : StItems g s blocks) : StItems g { s with items := s.items.modify j f } blocks := by
  obtain ⟨roots, nodup, bounded, chLt, chNodup, closed, finished⟩ := h
  have hch : ∀ p, Items.ch (s.items.modify j f) p = Items.ch s.items p :=
    Items.ch_modify_ch_eq j f hfc
  have hpar : ∀ p c, Items.IsParent (s.items.modify j f) p c ↔ Items.IsParent s.items p c :=
    fun p c => by unfold Items.IsParent; rw [hch]
  have hbel : ∀ a i, Items.Below (s.items.modify j f) a i ↔ Items.Below s.items a i :=
    fun a i => Items.Below_modify_ch_eq j f hfc
  have hty : ∀ i, Items.type (s.items.modify j f) i = Items.type s.items i :=
    fun i => Items.type_modify s.items j i f hf
  have hsz : (s.items.modify j f).size = s.items.size := Array.size_modify
  refine ⟨fun x hx p hp => roots x hx p ((hpar p x).1 hp), nodup,
    fun x hx y hy => hsz ▸ bounded x hx y ((hbel x y).1 hy),
    fun p c hp => hsz ▸ chLt p c ((hpar p c).1 hp), fun p => hch p ▸ chNodup p, fun i hi hty' => ?_,
    fun x i hVQ hb hsp => ?_⟩
  · rw [hsz] at hi
    rw [hty] at hty'
    rcases closed i hi hty' with ⟨x, hx, hb⟩ | ⟨b, hb, hB⟩
    · exact Or.inl ⟨x, hx, (hbel x i).2 hb⟩
    · exact Or.inr ⟨b, hb, hB.modify_root f hr fun e => hj (e ▸ hty')⟩
  · rw [hty] at hVQ hsp
    obtain ⟨b, hb', hB⟩ := finished x i hVQ ((hbel x i).1 hb) hsp
    exact ⟨b, hb', hB.modify_root f hr fun e => hj (e ▸ hsp)⟩

theorem StItems.pushEntry {g : Graph} {s : WalkState} {blocks : List StBlock} (v d idx : Nat)
    (dir : Bool) (q : ItemId) (hq : q < s.items.size) (hr : ∀ p, ¬ Items.IsParent s.items p q)
    (hfree : q ∉ readStack s.tstack) (h : StItems g s blocks) :
    StItems g { s with tstack := ⟨v, d, idx, setSides dir [q] []⟩ :: s.tstack } blocks := by
  obtain ⟨roots, nodup, bounded, chLt, chNodup, closed, finished⟩ := h
  have hmem := mem_readStack_push v d idx dir q s.tstack
  refine ⟨fun x hx => ?_, nodup_readStack_push v d idx dir q s.tstack hfree nodup,
    fun x hx y hy => ?_, chLt, chNodup, fun i hi hty => ?_, finished⟩
  · rcases (hmem x).1 hx with rfl | hx
    · exact hr
    · exact roots x hx
  · rcases (hmem x).1 hx with rfl | hx
    · rcases hy.cases_tail with rfl | ⟨c, -, hc⟩
      · exact hq
      · exact chLt c y hc
    · exact bounded x hx y hy
  · rcases closed i hi hty with ⟨x, hx, hb⟩ | hB
    · exact Or.inl ⟨x, (hmem x).2 (Or.inr hx), hb⟩
    · exact Or.inr hB

/-- `closeEars` for a tree edge with `lowval < d` keeps the st-reading: the child's entries `sub`
(reading as `ps`) become entries reading as `ps ++ [⟨stackDir[d], [Q e]⟩]` above `base`. -/
theorem closeEars_st {D d : Nat} {o : DfsOut} {s : WalkState} {curV v₀ : Nat} {hasVert : Bool}
    {sub base : List TEntry} {g : Graph} {ps : List StPiece} {blocks : List StBlock}
    (hE : s.EarFinish curV d o hasVert sub base) (hi : s.Inv' D) (hs : Shape s)
    (hD : D = d + 1) (he : o.e < s.g.ne) (hq : Items.ch s.items (edgeItem s.g o.e) = [])
    (hends : Items.PairEq (o.dest, s.stackVerts[d]!) s.g.edges[o.e]!) (hv : v₀ < s.g.nv)
    (ht : o.cls.isTree = true) (hlow : o.cls.lowval d < d)
    (hqty : Items.type s.items (edgeItem s.g o.e) = .Q) {pre B : List TEntry} (hB : base = pre ++ B)
    (hR : StRead s.items (sub ++ pre) ps) (hI : StItems g s blocks) :
    L1StInv g s d (ps ++ [⟨s.stackDir[d]!, [edgeItem s.g o.e]⟩]) blocks B (feS₁ d o s) := by
  obtain ⟨hi', lo, hrange, hc⟩ := L1Ctx.ofEar hE hs hD ht hlow
  have h0 := l1_init hE hi hs hD he hq hends hrange
  show L1StInv g s d _ blocks B
    ((loop (ceS₁ o.dest d o.e (feS₀ d o s)).tstack.length (loop1Cond d) (loop1Body d s.stackDir[d]!)).run
      (ceS₁ o.dest d o.e (feS₀ d o s))).2
  refine l1St_loop hc h0 hv hB ?_ _
  set f : Item → Item :=
    fun it => { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) } with hf
  have hst : ceS₁ o.dest d o.e (feS₀ d o s) =
      { s with items := s.items.modify (edgeItem s.g o.e) f,
               tstack := ⟨o.dest, d, s.nxtEdgeIdx, setSides s.stackDir[d]! [edgeItem s.g o.e] []⟩ ::
                 s.tstack } := rfl
  rw [hst]
  have hqS : ¬ (Items.type s.items (edgeItem s.g o.e) = .S ∨ Items.type s.items (edgeItem s.g o.e) = .P ∨
      Items.type s.items (edgeItem s.g o.e) = .R) := by rw [hqty]; decide
  have hI₀ : StItems g { s with items := s.items.modify (edgeItem s.g o.e) f } blocks :=
    StItems.modify_root _ f (fun _ => rfl) (fun _ => rfl) hE.q_root hqS hI
  have hty₀ : Items.type (s.items.modify (edgeItem s.g o.e) f) (edgeItem s.g o.e) = .Q := by
    rw [Items.type_modify s.items _ _ f (fun _ => rfl)]; exact hqty
  have hroot₀ : ∀ p, ¬ Items.IsParent (s.items.modify (edgeItem s.g o.e) f) p (edgeItem s.g o.e) := by
    intro p hp
    unfold Items.IsParent at hp
    rw [Items.ch_modify_ch_eq _ f (fun _ => rfl)] at hp
    exact hE.q_root p hp
  have hfree : edgeItem s.g o.e ∉ readStack s.tstack := fun hm => by
    obtain ⟨u, hu, hmu⟩ := mem_readStack_exists hm
    exact hE.q_free u hu hmu
  have hqlt : edgeItem s.g o.e < (s.items.modify (edgeItem s.g o.e) f).size := by
    rw [Array.size_modify]; have := hs.size; show 1 + s.g.nv + o.e < _; omega
  have hR₀ : StRead (s.items.modify (edgeItem s.g o.e) f) (sub ++ pre) ps :=
    hR.congr fun x _ y _ =>
      ⟨Items.type_modify s.items _ _ f (fun _ => rfl), Items.ch_modify_ch_eq _ f (fun _ => rfl) y⟩
  refine ⟨⟨⟨o.dest, d, s.nxtEdgeIdx, setSides s.stackDir[d]! [edgeItem s.g o.e] []⟩ :: (sub ++ pre), ?_,
    StRead.pushEntry _ _ _ _ _ (Or.inr hty₀) hR₀⟩, ?_, fun k _ => rfl⟩
  · dsimp only; simp only [hE.tstack, hB, List.cons_append, List.append_assoc]
  · exact StItems.pushEntry (s := { s with items := s.items.modify (edgeItem s.g o.e) f }) o.dest d
      s.nxtEdgeIdx s.stackDir[d]! (edgeItem s.g o.e) hqlt hroot₀ hfree hI₀

end Spqr
