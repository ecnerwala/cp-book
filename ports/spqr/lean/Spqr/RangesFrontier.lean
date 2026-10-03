import Spqr.RangesTree
import Spqr.EarFrontier
import Spqr.Proofs.ForestSpec

namespace Spqr.WalkState
open WalkM

variable {σ : List Nat} {n D : Nat} {s : WalkState}

theorem infix_interval (hnd : σ.Nodup) {xs : List Nat} (hx : xs <:+: σ) :
    ∀ a b c, a ≤ b → b ≤ c → c < σ.length → σ[a]! ∈ xs → σ[c]! ∈ xs → σ[b]! ∈ xs := by
  obtain ⟨pre, post, rfl⟩ := hx
  intro a b c hab hbc hc ha hcc
  have hi : ∀ e ∈ xs, (pre ++ xs ++ post).idxOf e = pre.length + xs.idxOf e := by
    intro e he
    have hn : e ∉ pre := by
      have hd := (List.nodup_append.1 (List.nodup_append.1 hnd).1).2.2
      exact fun hp => hd e hp e he rfl
    rw [List.append_assoc, List.idxOf_append_of_notMem hn, List.idxOf_append_of_mem he]
  have hal := List.idxOf_lt_length_iff.2 ha
  have hcl := List.idxOf_lt_length_iff.2 hcc
  have hia := hi _ ha
  have hic := hi _ hcc
  rw [idxOf_getElem! hnd (by omega : a < (pre ++ xs ++ post).length)] at hia
  rw [idxOf_getElem! hnd hc] at hic
  have hb : b < (pre ++ xs ++ post).length := by omega
  have hpre : pre.length ≤ b := by omega
  have hxs : b < (pre ++ xs).length := by simp only [List.length_append]; omega
  rw [getElem!_pos (pre ++ xs ++ post) b hb, List.getElem_append_left hxs, List.getElem_append_right hpre]
  exact List.getElem_mem _

theorem subEdges_iff_mem_block (o : DfsOut) (e : Nat) : subEdges o e ↔ e ∈ o.block := by
  cases o with
  | back => simp [subEdges, DfsOut.block, DfsOut.e]
  | tree e' cls child =>
    simp only [subEdges, DfsOut.e, DfsOut.block, List.mem_append, List.mem_singleton]
    rw [child.edgePostorder_perm_edges.mem_iff]
    exact or_comm

theorem subEdges_interval (hnd : σ.Nodup) {o : DfsOut} (ho : o.block <:+: σ) :
    ∀ a b c, a ≤ b → b ≤ c → c < σ.length →
      subEdges o σ[a]! → subEdges o σ[c]! → subEdges o σ[b]! := by
  simpa only [subEdges_iff_mem_block] using infix_interval hnd ho

/-- The ear's boundary split owns every gap before its pending tree edge. -/
theorem boundaryAdj_of_book {curV d orig : Nat} {o : DfsOut} {hasVert : Bool}
    (hge : d ≤ o.cls.lowval d) (hb : FinishBook curV d o orig hasVert s)
    (hnd : σ.Nodup) (hσ : ∀ e ∈ σ, e < s.g.ne) (hpos : σ[n]? = some o.e)
    (ho : o.block <:+: σ) : BoundaryAdj σ n d o s := by
  obtain ⟨sub, base, hl, he⟩ := hb.ear
  obtain ⟨hn, hnval⟩ := List.getElem?_eq_some_iff.1 hpos
  have hnval : σ[n]! = o.e := by rw [getElem!_pos σ n hn]; exact hnval
  have hfill : ∀ a b, a ≤ b → b < n →
      (∃ t ∈ sub, t.piece s.g s.items σ[a]!) → ∃ t ∈ sub, t.edges s.g s.items σ[b]! := by
    intro a b hab hbn ha
    obtain ⟨t, ht, ha⟩ := ha
    have hsub := he.sub_edges t ht _ (getElem!_lt hσ (by omega : a < σ.length)) ha.edges
    have hnext : subEdges o σ[n]! := by rw [hnval]; exact Or.inl rfl
    have hmid := subEdges_interval hnd ho a b n hab (Nat.le_of_lt hbn) hn hsub hnext
    apply he.sub_cover _ (getElem!_lt hσ (by omega : b < σ.length)) _ hmid
    intro heq
    have hid := idxOf_getElem! hnd (by omega : b < σ.length)
    rw [heq, ← hnval, idxOf_getElem! hnd hn] at hid
    omega
  refine ⟨hpos, ?_, ?_⟩
  · intro ht hlv t rest hts
    obtain ⟨u, hsub, -, -⟩ := he.bd_bridge ht hlv
    have hu := he.tstack
    rw [hsub, List.singleton_append, hts, List.cons.injEq] at hu
    obtain ⟨rfl, rfl⟩ := hu
    have hs := he.bd_side ht hge
    simp only [hlv, ↓reduceIte, hts, List.head?_cons, Option.mem_def, Option.some.injEq] at hs
    refine ⟨hs _ rfl, ?_⟩
    intro a b hab hbn ha
    obtain ⟨u, hu, huE⟩ := hfill a b hab hbn ⟨t, by simp [hsub], ha⟩
    have hu : u = t := by simpa [hsub] using hu
    exact hu ▸ huE
  · intro ht hlv t₁ t₂ rest hts
    obtain ⟨u₁, u₂, hsub, -, -, -, -⟩ := he.bd_comp ht hge hlv
    have hu := he.tstack
    simp only [hsub, List.cons_append, List.nil_append, hts, List.cons.injEq] at hu
    obtain ⟨rfl, rfl, rfl⟩ := hu
    have hs := he.bd_side ht hge
    simp only [hlv, ↓reduceIte, hts, List.head?_cons, List.tail_cons, Option.mem_def, Option.some.injEq] at hs
    refine ⟨hs.1 _ rfl, hs.2 _ rfl, ?_⟩
    intro a b hab hbn ha
    have hx : ∃ t ∈ sub, t.piece s.g s.items σ[a]! := by simpa [hsub] using ha
    simpa [hsub] using hfill a b hab hbn hx

/-- A covered gap between the top two pieces cannot belong to a lower entry. -/
theorem RangesInv.mergeAdj_of_cover (h : s.RangesInv σ n D) (hnd : σ.Nodup)
    (hσ : ∀ e ∈ σ, e < s.g.ne)
    (hcover : ∀ cur nxt rest, s.tstack = cur :: nxt :: rest →
      ∀ a b c, a ≤ b → b ≤ c → c < σ.length →
        nxt.piece s.g s.items σ[a]! → cur.piece s.g s.items σ[c]! →
        ∃ t ∈ s.tstack, t.edges s.g s.items σ[b]!) : MergeAdj σ s := by
  intro cur nxt rest hts a b c hab hbc hc ha hcc
  obtain ⟨t, ht, he⟩ := hcover cur nxt rest hts a b c hab hbc hc ha hcc
  rw [hts] at ht
  rcases List.mem_cons.1 ht with rfl | ht
  · exact Or.inl he
  rcases List.mem_cons.1 ht with rfl | ht
  · exact Or.inr he
  have hlt := h.ordered [cur] nxt rest hts t ht _ _
    (getElem!_lt hσ (by omega : a < σ.length)) (getElem!_lt hσ (by omega : b < σ.length)) ha he
  rw [idxOf_getElem! hnd (by omega : a < σ.length),
    idxOf_getElem! hnd (by omega : b < σ.length)] at hlt
  omega

/-- Exact frontier ownership of an interval supplies the missing merge adjacency. -/
theorem RangesInv.mergeAdj_of_frontier (h : s.RangesInv σ n D) (hnd : σ.Nodup)
    (hσ : ∀ e ∈ σ, e < s.g.ne) {orig : Nat} {base : List TEntry} {E : Nat → Prop}
    (hf : FrontierOwns orig base E s) (hlen : orig + 2 ≤ s.tstack.length)
    (hinterval : ∀ a b c, a ≤ b → b ≤ c → c < σ.length →
      E σ[a]! → E σ[c]! → E σ[b]!) : MergeAdj σ s := by
  apply h.mergeAdj_of_cover hnd hσ
  intro cur nxt rest hts a b c hab hbc hc ha hcc
  have hlen' : 2 ≤ s.tstack.length - orig := by omega
  have hfront : s.tstack.take (s.tstack.length - orig) =
      cur :: nxt :: rest.take (s.tstack.length - orig - 2) := by
    generalize s.tstack.length - orig = k at hlen' ⊢
    obtain ⟨k, rfl⟩ := Nat.exists_eq_add_of_le hlen'
    rw [hts, show 2 + k = k + 1 + 1 by omega]
    simp
  have he : E σ[b]! := hinterval a b c hab hbc hc
    ((hf.2.2 _ (getElem!_lt hσ (by omega : a < σ.length))).1 ⟨nxt, by simp [hfront], ha.edges⟩)
    ((hf.2.2 _ (getElem!_lt hσ hc)).1 ⟨cur, by simp [hfront], hcc.edges⟩)
  obtain ⟨t, ht, he⟩ := (hf.2.2 _ (getElem!_lt hσ (by omega : b < σ.length))).2 he
  exact ⟨t, List.mem_of_mem_take ht, he⟩

theorem RangesInv.iter_merge_ranges (h : s.RangesInv σ n D) (hnd : σ.Nodup)
    (hσ : ∀ e ∈ σ, e < s.g.ne) {cond : WalkM Bool} {orig : Nat} {base : List TEntry} {E : Nat → Prop}
    (hinterval : ∀ a b c, a ≤ b → b ≤ c → c < σ.length → E σ[a]! → E σ[c]! → E σ[b]!)
    (hok : ∀ k, (∀ j, j ≤ k → result cond (iter mergeTstackTops j s) = true) →
      MergeTopOk D (iter mergeTstackTops k s))
    (hf : ∀ k, (∀ j, j ≤ k → result cond (iter mergeTstackTops j s) = true) →
      FrontierOwns orig base E (iter mergeTstackTops k s) ∧ orig + 2 ≤ (iter mergeTstackTops k s).tstack.length) :
    ∀ k, (∀ j, j < k → result cond (iter mergeTstackTops j s) = true) →
      (iter mergeTstackTops k s).RangesInv σ n D ∧ (iter mergeTstackTops k s).g = s.g := by
  intro k
  induction k with
  | zero => exact fun _ => ⟨h, rfl⟩
  | succ k ih =>
    intro hk
    obtain ⟨hr, hg⟩ := ih fun j hj => hk j (by omega)
    have hk' : ∀ j, j ≤ k → result cond (iter mergeTstackTops j s) = true := fun j hj => hk j (by omega)
    have hσ' : ∀ e ∈ σ, e < (iter mergeTstackTops k s).g.ne := by rwa [hg]
    have ha := hr.mergeAdj_of_frontier hnd hσ' (hf k hk').1 (hf k hk').2 hinterval
    rw [iter_succ']
    exact ⟨hr.mergeTop' hnd hσ' (hok k hk') ha, by rw [run_mergeTstackTops]; exact hg⟩

theorem RangesInv.iter_mergeAdj (h : s.RangesInv σ n D) (hnd : σ.Nodup)
    (hσ : ∀ e ∈ σ, e < s.g.ne) {cond : WalkM Bool} {orig : Nat} {base : List TEntry} {E : Nat → Prop}
    (hinterval : ∀ a b c, a ≤ b → b ≤ c → c < σ.length → E σ[a]! → E σ[c]! → E σ[b]!)
    (hok : ∀ k, (∀ j, j ≤ k → result cond (iter mergeTstackTops j s) = true) →
      MergeTopOk D (iter mergeTstackTops k s))
    (hf : ∀ k, (∀ j, j ≤ k → result cond (iter mergeTstackTops j s) = true) →
      FrontierOwns orig base E (iter mergeTstackTops k s) ∧ orig + 2 ≤ (iter mergeTstackTops k s).tstack.length)
    (k : Nat) (hk : ∀ j, j ≤ k → result cond (iter mergeTstackTops j s) = true) :
    MergeAdj σ (iter mergeTstackTops k s) := by
  obtain ⟨hr, hg⟩ := h.iter_merge_ranges hnd hσ hinterval hok hf k (fun j hj => hk j (by omega))
  exact hr.mergeAdj_of_frontier hnd (by rwa [hg]) (hf k hk).1 (hf k hk).2 hinterval

theorem mergeLateAdj_of_frontier {o : DfsOut} {d orig : Nat}
    (hf : Frontier (o := o) d orig s) (ht : o.cls.isTree = true) (hlow : o.cls.lowval d < d)
    (h : (feS₁ d o s).RangesInv σ n D) (hnd : σ.Nodup)
    (hσ : ∀ e ∈ σ, e < (feS₁ d o s).g.ne) (ho : o.block <:+: σ)
    (hok : MergeLateOk D d (feS₁ d o s)) : MergeLateAdj σ d (feS₁ d o s) := by
  refine ⟨fun hc => h.iter_mergeAdj hnd hσ (orig := orig)
    (base := s.tstack.drop (s.tstack.length - orig)) (subEdges_interval hnd ho) (hok.body hc) ?_⟩
  intro k hk
  obtain ⟨hown, hlen⟩ := hf.loop2 ht hlow k (fun j hj => hk j (by omega))
  exact ⟨hown, hlen (hk k (Nat.le_refl _))⟩

theorem loop3Adj_of_frontier {o : DfsOut} {d orig : Nat}
    (hf : Frontier (o := o) d orig s) (ht : o.cls.isTree = true) (hlow : o.cls.lowval d < d)
    (h : (feS₂ d o s).RangesInv σ n D) (hnd : σ.Nodup)
    (hσ : ∀ e ∈ σ, e < (feS₂ d o s).g.ne) (ho : o.block <:+: σ)
    (hok : CloseVertOk D s.stackVerts[d]! s.stackDir[d]! o.cls.isType1 orig (feSingle d o s) (feS₂ d o s))
    (hty : o.cls.isType1 = false) :
    ∀ k, (∀ j, j ≤ k → result (loop3Cond orig) (iter mergeTstackTops j (feS₂ d o s)) = true) →
      MergeAdj σ (iter mergeTstackTops k (feS₂ d o s)) := by
  apply h.iter_mergeAdj hnd hσ (orig := orig)
    (base := s.tstack.drop (s.tstack.length - orig)) (subEdges_interval hnd ho) (hok.loop3 hty)
  intro k hk
  obtain ⟨hown, hlen⟩ := hf.loop3 ht hlow hty k (fun j hj => hk j (by omega))
  exact ⟨hown, hlen (hk k (Nat.le_refl _))⟩

end Spqr.WalkState
