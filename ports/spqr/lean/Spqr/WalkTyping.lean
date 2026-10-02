import Mathlib.Tactic.SplitIfs
import Spqr.ItemSpec

/-!
# Typing of the walk's items

The part of `Items.WF` that depends only on which items the walk allocates and which `modifyItem`
calls it makes (not on the tstack semantics): item types are fixed at creation, `V` items only get
`Q` children, `I`/`O` items are leaves, and the recorded endpoints `vs` have the per-type shape.
-/

namespace Spqr

/-- `vs_shape` of `Items.Endpoints`, per type. -/
def VsShape (ty : NodeType) (vs : Option Nat × Option Nat) : Prop :=
  match ty with
  | .F | .V => vs = (none, none)
  | .O => ∃ v, vs = (some v, none)
  | .Q => (∃ v, vs = (some v, none)) ∨ (∃ u v, vs = (some u, some v))
  | _ => ∃ u v, vs = (some u, some v)

/-- The weakest-precondition of a `WalkM` program. -/
def WalkM.wp (m : WalkM α) (Q : α → WalkState → Prop) (s : WalkState) : Prop :=
  Q (m.run s).1 (m.run s).2

namespace WalkM

variable {s : WalkState}

theorem wp_pure {Q : α → WalkState → Prop} (a : α) : wp (pure a : WalkM α) Q s = Q a s := rfl
theorem wp_bind {Q : α → WalkState → Prop} (m : WalkM β) (f : β → WalkM α) :
    wp (m >>= f) Q s = wp m (fun b s' => wp (f b) Q s') s := rfl
theorem wp_map {Q : α → WalkState → Prop} (f : β → α) (m : WalkM β) : wp (f <$> m) Q s = wp m (fun b s' => Q (f b) s') s := rfl
theorem wp_get {Q : WalkState → WalkState → Prop} : wp (get : WalkM WalkState) Q s = Q s s := rfl
theorem wp_set {Q : Unit → WalkState → Prop} (s' : WalkState) : wp (set s' : WalkM Unit) Q s = Q () s' := rfl
theorem wp_modify {Q : Unit → WalkState → Prop} (f : WalkState → WalkState) : wp (modify f : WalkM Unit) Q s = Q () (f s) := rfl
theorem wp_ite {Q : α → WalkState → Prop} (c : Prop) [Decidable c] (m₁ m₂ : WalkM α) :
    wp (if c then m₁ else m₂) Q s = if c then wp m₁ Q s else wp m₂ Q s := by split <;> rfl
theorem wp_dite {Q : α → WalkState → Prop} (c : Prop) [Decidable c] (m₁ : c → WalkM α) (m₂ : ¬ c → WalkM α) :
    wp (if h : c then m₁ h else m₂ h) Q s = if h : c then wp (m₁ h) Q s else wp (m₂ h) Q s := by
  split <;> rfl

theorem wp_modifyItem {Q : Unit → WalkState → Prop} (i : Nat) (f : Item → Item) :
    wp (modifyItem i f) Q s = Q () { s with items := s.items.modify i f } := rfl
theorem wp_getItem {Q : Item → WalkState → Prop} (i : ItemId) : wp (getItem i) Q s = Q s.items[i]! s := rfl
theorem wp_allocItem {Q : ItemId → WalkState → Prop} (ty : NodeType) :
    wp (allocItem ty) Q s = Q s.items.size { s with items := s.items.push ⟨ty, (none, none), []⟩ } := rfl
theorem wp_stackDir {Q : Bool → WalkState → Prop} (d : Nat) : wp (stackDir d) Q s = Q s.stackDir[d]! s := rfl
theorem wp_setStackDir {Q : Unit → WalkState → Prop} (d : Nat) (b : Bool) :
    wp (setStackDir d b) Q s = Q () { s with stackDir := s.stackDir.set! d b } := rfl
theorem wp_makeVs {Q : Option Nat × Option Nat → WalkState → Prop} (a b : Nat) :
    wp (makeVs a b) Q s = Q (setSides s.stackDir[b]! (some s.stackVerts[b]!) (some a)) s := rfl
theorem wp_cur {Q : TEntry → WalkState → Prop} : wp cur Q s = Q s.tstack.head! s := rfl
theorem wp_nxt {Q : TEntry → WalkState → Prop} : wp nxt Q s = Q s.tstack.tail.head! s := rfl
theorem wp_tstackSize {Q : Nat → WalkState → Prop} : wp tstackSize Q s = Q s.tstack.length s := rfl
theorem wp_modifyCur {Q : Unit → WalkState → Prop} (f : TEntry → TEntry) :
    wp (modifyCur f) Q s =
      Q () { s with tstack := match s.tstack with | a :: rest => f a :: rest | [] => [] } := rfl
theorem wp_modifyNxt {Q : Unit → WalkState → Prop} (f : TEntry → TEntry) :
    wp (modifyNxt f) Q s =
      Q () { s with tstack := match s.tstack with | a :: b :: rest => a :: f b :: rest | l => l } := rfl
theorem wp_popTstack {Q : TEntry → WalkState → Prop} : wp popTstack Q s = Q s.tstack.head! { s with tstack := s.tstack.tail } := rfl
theorem wp_pushTstack {Q : Unit → WalkState → Prop} (vStart topDepth : Nat) (item : ItemId) :
    wp (pushTstack vStart topDepth item) Q s =
      Q () { s with tstack :=
        ⟨vStart, topDepth, s.nxtEdgeIdx, setSides s.stackDir[topDepth]! [item] []⟩ :: s.tstack } := rfl
theorem wp_pushVertTstack {Q : Unit → WalkState → Prop} (v topDepth : Nat) :
    wp (pushVertTstack v topDepth) Q s =
      Q () { s with tstack :=
        ⟨v, topDepth, s.nxtEdgeIdx, setSides s.stackDir[topDepth]! [vertItem v] []⟩ :: s.tstack } := rfl
theorem wp_pushEdgeTstack {Q : Unit → WalkState → Prop} (vStart topDepth e : Nat) :
    wp (pushEdgeTstack vStart topDepth e) Q s =
      Q () { s with tstack :=
        ⟨vStart, topDepth, s.nxtEdgeIdx, setSides s.stackDir[topDepth]! [edgeItem s.g e] []⟩ :: s.tstack } := rfl

/-- The tstack after `mergeTstackTops`. -/
def mergeTop : List TEntry → List TEntry
  | b :: a :: rest =>
    { a with topDepth := min a.topDepth b.topDepth,
             spans := (b.spans.1 ++ a.spans.1, a.spans.2 ++ b.spans.2) } :: rest
  | _ => []

theorem wp_mergeTstackTops {Q : Unit → WalkState → Prop} : wp mergeTstackTops Q s = Q () { s with tstack := mergeTop s.tstack } := by
  rcases s with ⟨g, tern, items, sv, sd, ni, fo, _ | ⟨b, _ | ⟨a, rest⟩⟩, tb, tl⟩ <;> rfl

theorem wp_mono {Q Q' : α → WalkState → Prop} (m : WalkM α) (h : wp m Q s)
    (hQ : ∀ a s', Q a s' → Q' a s') : wp m Q' s := hQ _ _ h

/-- Invariant rule for `loop`. -/
theorem wp_loop {Q : Unit → WalkState → Prop} (I : WalkState → Prop) (n : Nat) (cond : WalkM Bool) (body : WalkM Unit)
    (hc : ∀ s, I s → wp cond (fun _ s' => I s') s) (hb : ∀ s, I s → wp body (fun _ s' => I s') s)
    (hs : I s) (hQ : ∀ s', I s' → Q () s') : wp (loop n cond body) Q s := by
  induction n generalizing s with
  | zero => exact hQ _ hs
  | succ n ih =>
    simp only [loop, wp_bind, wp_ite]
    refine wp_mono _ (hc s hs) fun b s₁ h₁ => ?_
    split
    · exact wp_mono _ (hb s₁ h₁) fun _ s₂ h₂ => ih h₂
    · exact hQ _ h₁

end WalkM

namespace Items

variable (items : Items)

theorem type_modify (i j : Nat) (f : Item → Item) (hf : ∀ it, (f it).type = it.type) :
    Items.type (items.modify i f) j = items.type j := by
  simp only [type, Array.getElem?_modify]
  split
  · subst i; cases items[j]? <;> simp [hf]
  · rfl

theorem ch_modify_ne (i j : Nat) (f : Item → Item) (h : i ≠ j) :
    Items.ch (items.modify i f) j = items.ch j := by
  simp [ch, Array.getElem?_modify, h]
theorem ch_modify_self (i : Nat) (f : Item → Item) (h : i < items.size) :
    Items.ch (items.modify i f) i = (f items[i]).ch := by
  simp [ch, Array.getElem?_modify, Array.getElem?_eq_getElem h]
theorem ch_modify_of_ch (i j : Nat) (f : Item → Item) (hf : ∀ it, (f it).ch = it.ch) :
    Items.ch (items.modify i f) j = items.ch j := by
  simp only [ch, Array.getElem?_modify]
  split
  · subst i; cases items[j]? <;> simp [hf]
  · rfl

theorem vs_modify_ne (i j : Nat) (f : Item → Item) (h : i ≠ j) :
    Items.vs (items.modify i f) j = items.vs j := by
  simp [vs, Array.getElem?_modify, h]
theorem vs_modify_self (i : Nat) (f : Item → Item) (h : i < items.size) :
    Items.vs (items.modify i f) i = (f items[i]).vs := by
  simp [vs, Array.getElem?_modify, Array.getElem?_eq_getElem h]
theorem vs_modify_of_vs (i j : Nat) (f : Item → Item) (hf : ∀ it, (f it).vs = it.vs) :
    Items.vs (items.modify i f) j = items.vs j := by
  simp only [vs, Array.getElem?_modify]
  split
  · subst i; cases items[j]? <;> simp [hf]
  · rfl

theorem type_push (x : Item) (j : Nat) :
    Items.type (items.push x) j = if j = items.size then x.type else items.type j := by
  simp only [type, Array.getElem?_push]; split <;> rfl
theorem ch_push (x : Item) (j : Nat) :
    Items.ch (items.push x) j = if j = items.size then x.ch else items.ch j := by
  simp only [ch, Array.getElem?_push]; split <;> rfl
theorem vs_push (x : Item) (j : Nat) :
    Items.vs (items.push x) j = if j = items.size then x.vs else items.vs j := by
  simp only [vs, Array.getElem?_push]; split <;> rfl

theorem type_of_le (j : Nat) (h : items.size ≤ j) : items.type j = .F := by
  simp [type, Array.getElem?_eq_none h]
theorem ch_of_le (j : Nat) (h : items.size ≤ j) : items.ch j = [] := by
  simp [ch, Array.getElem?_eq_none h]
theorem vs_of_le (j : Nat) (h : items.size ≤ j) : items.vs j = (none, none) := by
  simp [vs, Array.getElem?_eq_none h]

theorem type_getElem! (j : Nat) (h : j < items.size) : items[j]!.type = items.type j := by
  simp [type, getElem!_pos, h]
theorem ch_getElem! (j : Nat) (h : j < items.size) : items[j]!.ch = items.ch j := by
  simp [ch, getElem!_pos, h]

/-- The typing facts about the items, with `ex` the items whose `vs` may still be unset. -/
structure Typing (g : Graph) (ex : ItemId → Prop) : Prop where
  size : 1 + g.nv + g.ne ≤ items.size
  root : items.type rootItem = .F
  vert : ∀ v, v < g.nv → items.type (vertItem v) = .V
  edge : ∀ e, e < g.ne → items.type (edgeItem g e) = .Q
  node : ∀ i, 1 + g.nv + g.ne ≤ i → i < items.size → items.type i ∉ [NodeType.F, .V, .Q]
  v_children : ∀ v c, v < g.nv → items.IsParent (vertItem v) c → items.type c = .Q
  i_o_leaf : ∀ i, i < items.size → items.type i = .I ∨ items.type i = .O → items.ch i = []
  vs_shape : ∀ i, i < items.size → ¬ ex i → VsShape (items.type i) (items.vs i)
  vs_lt : ∀ i u, (items.vs i).1 = some u ∨ (items.vs i).2 = some u → u < g.nv

variable {items}

theorem Typing.mono {g : Graph} {ex ex' : ItemId → Prop} (h : items.Typing g ex)
    (hex : ∀ i, ex i → ex' i) : items.Typing g ex' :=
  { h with vs_shape := fun i hi hi' => h.vs_shape i hi fun h' => hi' (hex i h') }


theorem type_getElem!' (j : Nat) : items[j]!.type = items.type j := by
  by_cases h : j < items.size
  · simp [type, getElem!_pos, h]
  · simp [type, getElem!_neg, h]; rfl

/-- Allocating a node (`I`, `O`, `S`, `P`, `R`) with unset `vs`. -/
theorem Typing.push {g : Graph} {ex : ItemId → Prop} (h : items.Typing g ex) (ty : NodeType)
    (hty : ty ∉ [NodeType.F, .V, .Q]) :
    Items.Typing (items.push ⟨ty, (none, none), []⟩) g (fun i => ex i ∨ i = items.size) where
  size := by have := h.size; simp; omega
  root := by
    rw [type_push]; split
    · next hj => have := h.size; have hj : (0 : Nat) = items.size := hj; omega
    · exact h.root
  vert v hv := by
    rw [type_push]; split
    · next hj => have := h.size; simp [vertItem] at hj; omega
    · exact h.vert v hv
  edge e he := by
    rw [type_push]; split
    · next hj => have := h.size; simp [edgeItem] at hj; omega
    · exact h.edge e he
  node i hi hi' := by
    rw [type_push]; split
    · exact hty
    · exact h.node i hi (by simp at hi'; omega)
  v_children v c hv hc := by
    rw [type_push]
    simp only [IsParent, ch_push] at hc
    split at hc
    · next hj => have := h.size; simp [vertItem] at hj; omega
    · split
      · exact absurd (h.v_children v c hv hc) (by subst c; simp [type_of_le])
      · exact h.v_children v c hv hc
  i_o_leaf i hi ht := by
    rw [type_push] at ht; rw [ch_push]; split
    · rfl
    · exact h.i_o_leaf i (by simp at hi; omega) (by simpa [*] using ht)
  vs_shape i hi hex := by
    rw [type_push, vs_push]; simp only [not_or] at hex; simp only [hex, ite_false]
    have hi' : i < items.size := by simp at hi; exact Nat.lt_of_le_of_ne (Nat.le_of_lt_succ hi) hex.2
    exact h.vs_shape i hi' hex.1
  vs_lt i u hu := by
    rw [vs_push] at hu; split at hu
    · simp at hu
    · exact h.vs_lt i u hu

/-- Setting `vs` of an item to a value of the right shape. -/
theorem Typing.modify_vs {g : Graph} {ex : ItemId → Prop} (h : items.Typing g ex) (i : Nat)
    (vs : Option Nat × Option Nat) (hshape : VsShape (items.type i) vs)
    (hlt : ∀ u, vs.1 = some u ∨ vs.2 = some u → u < g.nv) :
    Items.Typing (items.modify i fun it => { it with vs := vs }) g (fun j => ex j ∧ j ≠ i) where
  size := by simpa using h.size
  root := by rw [type_modify _ _ _ _ (by intro _; rfl)]; exact h.root
  vert v hv := by rw [type_modify _ _ _ _ (by intro _; rfl)]; exact h.vert v hv
  edge e he := by rw [type_modify _ _ _ _ (by intro _; rfl)]; exact h.edge e he
  node j hj hj' := by rw [type_modify _ _ _ _ (by intro _; rfl)]; exact h.node j hj (by simpa using hj')
  v_children v c hv hc := by
    rw [type_modify _ _ _ _ (by intro _; rfl)]
    rw [IsParent, ch_modify_of_ch _ _ _ _ (by intro _; rfl)] at hc; exact h.v_children v c hv hc
  i_o_leaf j hj ht := by
    rw [type_modify _ _ _ _ (by intro _; rfl)] at ht; rw [ch_modify_of_ch _ _ _ _ (by intro _; rfl)]
    exact h.i_o_leaf j (by simpa using hj) ht
  vs_shape j hj hex := by
    rw [type_modify _ _ _ _ (by intro _; rfl)]
    by_cases hij : i = j
    · subst hij; rw [vs_modify_self _ _ _ (by simpa using hj)]; exact hshape
    · rw [vs_modify_ne _ _ _ _ hij]; exact h.vs_shape j (by simpa using hj) fun h' => hex ⟨h', Ne.symm hij⟩
  vs_lt j u hu := by
    by_cases hij : i = j
    · subst hij
      by_cases hi : i < items.size
      · rw [vs_modify_self _ _ _ hi] at hu; exact hlt u hu
      · rw [vs_of_le _ _ (by simpa using Nat.le_of_not_lt hi)] at hu; simp at hu
    · rw [vs_modify_ne _ _ _ _ hij] at hu; exact h.vs_lt j u hu

/-- Closing an `S`/`P`/`R` node: its `vs` and `ch` are set. -/
theorem Typing.modify_vs_ch {g : Graph} {ex : ItemId → Prop} (h : items.Typing g ex) (i : Nat)
    (hty : items.type i ∈ [NodeType.S, .P, .R]) (vs : Option Nat × Option Nat) (ch : List ItemId)
    (hshape : ∃ u v, vs = (some u, some v))
    (hlt : ∀ u, vs.1 = some u ∨ vs.2 = some u → u < g.nv) :
    Items.Typing (items.modify i fun it => { it with vs := vs, ch := ch }) g (fun j => ex j ∧ j ≠ i) where
  size := by simpa using h.size
  root := by rw [type_modify _ _ _ _ (by intro _; rfl)]; exact h.root
  vert v hv := by rw [type_modify _ _ _ _ (by intro _; rfl)]; exact h.vert v hv
  edge e he := by rw [type_modify _ _ _ _ (by intro _; rfl)]; exact h.edge e he
  node j hj hj' := by rw [type_modify _ _ _ _ (by intro _; rfl)]; exact h.node j hj (by simpa using hj')
  v_children v c hv hc := by
    rw [type_modify _ _ _ _ (by intro _; rfl)]
    have hne : i ≠ vertItem v := fun h' => by rw [h', h.vert v hv] at hty; simp at hty
    rw [IsParent, ch_modify_ne _ _ _ _ hne] at hc; exact h.v_children v c hv hc
  i_o_leaf j hj ht := by
    rw [type_modify _ _ _ _ (by intro _; rfl)] at ht
    have hne : i ≠ j := fun h' => by subst h'; rcases ht with ht | ht <;> rw [ht] at hty <;> simp at hty
    rw [ch_modify_ne _ _ _ _ hne]
    exact h.i_o_leaf j (by simpa using hj) ht
  vs_shape j hj hex := by
    rw [type_modify _ _ _ _ (by intro _; rfl)]
    by_cases hij : i = j
    · subst hij; rw [vs_modify_self _ _ _ (by simpa using hj)]
      simp only [List.mem_cons, List.not_mem_nil, or_false] at hty
      rcases hty with hty | hty | hty <;> rw [hty] <;> exact hshape
    · rw [vs_modify_ne _ _ _ _ hij]; exact h.vs_shape j (by simpa using hj) fun h' => hex ⟨h', Ne.symm hij⟩
  vs_lt j u hu := by
    by_cases hij : i = j
    · subst hij
      by_cases hi : i < items.size
      · rw [vs_modify_self _ _ _ hi] at hu; exact hlt u hu
      · rw [vs_of_le _ _ (by simpa using Nat.le_of_not_lt hi)] at hu; simp at hu
    · rw [vs_modify_ne _ _ _ _ hij] at hu; exact h.vs_lt j u hu

/-- Writing the children of a `Q` item. -/
theorem Typing.modify_q_ch {g : Graph} {ex : ItemId → Prop} (h : items.Typing g ex) (e : Nat)
    (he : e < g.ne) (ch : List ItemId) :
    Items.Typing (items.modify (edgeItem g e) fun it => { it with ch := ch }) g ex where
  size := by simpa using h.size
  root := by rw [type_modify _ _ _ _ (by intro _; rfl)]; exact h.root
  vert v hv := by rw [type_modify _ _ _ _ (by intro _; rfl)]; exact h.vert v hv
  edge e he := by rw [type_modify _ _ _ _ (by intro _; rfl)]; exact h.edge e he
  node j hj hj' := by rw [type_modify _ _ _ _ (by intro _; rfl)]; exact h.node j hj (by simpa using hj')
  v_children v c hv hc := by
    rw [type_modify _ _ _ _ (by intro _; rfl)]
    have hne : edgeItem g e ≠ vertItem v := fun h' => by have : (1 + g.nv + e : Nat) = 1 + v := h'; omega
    rw [IsParent, ch_modify_ne _ _ _ _ hne] at hc; exact h.v_children v c hv hc
  i_o_leaf j hj ht := by
    rw [type_modify _ _ _ _ (by intro _; rfl)] at ht
    have hne : edgeItem g e ≠ j := fun h' => by
      subst h'; rcases ht with ht | ht <;> rw [h.edge e he] at ht <;> simp at ht
    rw [ch_modify_ne _ _ _ _ hne]
    exact h.i_o_leaf j (by simpa using hj) ht
  vs_shape j hj hex := by
    rw [type_modify _ _ _ _ (by intro _; rfl)]; rw [vs_modify_of_vs _ _ _ _ (by intro _; rfl)]
    exact h.vs_shape j (by simpa using hj) hex
  vs_lt j u hu := by rw [vs_modify_of_vs _ _ _ _ (by intro _; rfl)] at hu; exact h.vs_lt j u hu

/-- Appending a `Q` child to a `V` item. -/
theorem Typing.modify_v_ch {g : Graph} {ex : ItemId → Prop} (h : items.Typing g ex) (v e : Nat)
    (hv : v < g.nv) (he : e < g.ne) :
    Items.Typing (items.modify (vertItem v) fun it => { it with ch := it.ch ++ [edgeItem g e] }) g ex where
  size := by simpa using h.size
  root := by rw [type_modify _ _ _ _ (by intro _; rfl)]; exact h.root
  vert v hv := by rw [type_modify _ _ _ _ (by intro _; rfl)]; exact h.vert v hv
  edge e he := by rw [type_modify _ _ _ _ (by intro _; rfl)]; exact h.edge e he
  node j hj hj' := by rw [type_modify _ _ _ _ (by intro _; rfl)]; exact h.node j hj (by simpa using hj')
  v_children w c hw hc := by
    rw [type_modify _ _ _ _ (by intro _; rfl)]
    by_cases hvw : v = w
    · subst hvw
      have hlt : vertItem v < items.size := by have := h.size; show (1 + v : Nat) < _; omega
      simp only [IsParent, ch_modify_self _ _ _ hlt, List.mem_append, List.mem_singleton] at hc
      rcases hc with hc | hc
      · exact h.v_children v c hv (by simpa [IsParent, ch, Array.getElem?_eq_getElem hlt] using hc)
      · subst hc; exact h.edge e he
    · have hne : vertItem v ≠ vertItem w := fun h' => hvw (by have : (1 + v : Nat) = 1 + w := h'; omega)
      rw [IsParent, ch_modify_ne _ _ _ _ hne] at hc; exact h.v_children w c hw hc
  i_o_leaf j hj ht := by
    rw [type_modify _ _ _ _ (by intro _; rfl)] at ht
    have hne : vertItem v ≠ j := fun h' => by
      subst h'; rcases ht with ht | ht <;> rw [h.vert v hv] at ht <;> simp at ht
    rw [ch_modify_ne _ _ _ _ hne]
    exact h.i_o_leaf j (by simpa using hj) ht
  vs_shape j hj hex := by
    rw [type_modify _ _ _ _ (by intro _; rfl)]; rw [vs_modify_of_vs _ _ _ _ (by intro _; rfl)]
    exact h.vs_shape j (by simpa using hj) hex
  vs_lt j u hu := by rw [vs_modify_of_vs _ _ _ _ (by intro _; rfl)] at hu; exact h.vs_lt j u hu

/-- Appending to the root's children. -/
theorem Typing.modify_root_ch {g : Graph} {ex : ItemId → Prop} (h : items.Typing g ex) (l : List ItemId) :
    Items.Typing (items.modify rootItem fun it => { it with ch := it.ch ++ l }) g ex where
  size := by simpa using h.size
  root := by rw [type_modify _ _ _ _ (by intro _; rfl)]; exact h.root
  vert v hv := by rw [type_modify _ _ _ _ (by intro _; rfl)]; exact h.vert v hv
  edge e he := by rw [type_modify _ _ _ _ (by intro _; rfl)]; exact h.edge e he
  node j hj hj' := by rw [type_modify _ _ _ _ (by intro _; rfl)]; exact h.node j hj (by simpa using hj')
  v_children v c hv hc := by
    rw [type_modify _ _ _ _ (by intro _; rfl)]
    have hne : rootItem ≠ vertItem v := fun h' => by have : (0 : Nat) = 1 + v := h'; omega
    rw [IsParent, ch_modify_ne _ _ _ _ hne] at hc; exact h.v_children v c hv hc
  i_o_leaf j hj ht := by
    rw [type_modify _ _ _ _ (by intro _; rfl)] at ht
    have hne : rootItem ≠ j := fun h' => by
      subst h'; rcases ht with ht | ht <;> rw [h.root] at ht <;> simp at ht
    rw [ch_modify_ne _ _ _ _ hne]
    exact h.i_o_leaf j (by simpa using hj) ht
  vs_shape j hj hex := by
    rw [type_modify _ _ _ _ (by intro _; rfl)]; rw [vs_modify_of_vs _ _ _ _ (by intro _; rfl)]
    exact h.vs_shape j (by simpa using hj) hex
  vs_lt j u hu := by rw [vs_modify_of_vs _ _ _ _ (by intro _; rfl)] at hu; exact h.vs_lt j u hu

theorem initialItems_getElem? (g : Graph) (i : Nat) :
    (initialItems g)[i]? =
      if i = 0 then some ⟨.F, (none, none), []⟩
      else if i < 1 + g.nv then some ⟨.V, (none, none), []⟩
      else if i < 1 + g.nv + g.ne then some ⟨.Q, (none, none), []⟩ else none := by
  simp only [initialItems, Array.getElem?_append, Array.size_append, Array.size_replicate,
    Array.getElem?_replicate, List.size_toArray, List.length_cons, List.length_nil]
  split_ifs <;> first | rfl | omega | simp_all

theorem initialItems_size (g : Graph) : (initialItems g).size = 1 + g.nv + g.ne := by
  simp [initialItems]; omega

theorem initialItems_type (g : Graph) (i : Nat) :
    Items.type (initialItems g) i =
      if i = 0 then .F else if i < 1 + g.nv then .V else if i < 1 + g.nv + g.ne then .Q else .F := by
  simp only [Items.type, initialItems_getElem?]; split_ifs <;> rfl

theorem initialItems_ch (g : Graph) (i : Nat) : Items.ch (initialItems g) i = [] := by
  simp only [Items.ch, initialItems_getElem?]; split_ifs <;> rfl

theorem initialItems_vs (g : Graph) (i : Nat) : Items.vs (initialItems g) i = (none, none) := by
  simp only [Items.vs, initialItems_getElem?]; split_ifs <;> rfl

end Items

theorem getElem!_set!_lt {a : Array Nat} {nv i v : Nat} (hv : v < nv) (h : ∀ j : Nat, a[j]! < nv) (j : Nat) :
    (a.set! i v : Array Nat)[j]! < nv := by
  have hj := h j
  rw [Array.set!_eq_setIfInBounds]
  simp only [getElem!_def] at hj ⊢
  rw [Array.getElem?_setIfInBounds]
  by_cases hij : i = j
  · subst hij; simp only [ite_true]
    by_cases hi : i < a.size
    · simp only [hi, ite_true]; exact hv
    · simp only [hi, ite_false]; exact Nat.lt_of_le_of_lt (Nat.zero_le _) hv
  · simp only [hij, ite_false]; exact hj

theorem head!_vStart_lt {l : List TEntry} {nv : Nat} (hnv : 0 < nv) (hl : ∀ t ∈ l, t.vStart < nv) :
    l.head!.vStart < nv := by
  cases l with
  | nil => exact hnv
  | cons t rest => exact hl t (List.mem_cons_self ..)

theorem mergeTop_vStart {P : Nat → Prop} {l : List TEntry} (hl : ∀ t ∈ l, P t.vStart) :
    ∀ t ∈ WalkM.mergeTop l, P t.vStart := by
  match l with
  | b :: a :: rest =>
    intro t ht
    simp only [WalkM.mergeTop, List.mem_cons] at ht
    rcases ht with rfl | ht
    · exact hl a (by simp)
    · exact hl t (by simp [ht])
  | [] | [_] => simp [WalkM.mergeTop]

theorem modifyCur_vStart {P : Nat → Prop} {l : List TEntry} (f : TEntry → TEntry)
    (hf : ∀ t, P t.vStart → P (f t).vStart) (hl : ∀ t ∈ l, P t.vStart) :
    ∀ t ∈ (match l with | a :: rest => f a :: rest | [] => []), P t.vStart := by
  match l with
  | [] => simp
  | a :: rest =>
    intro t ht
    simp only [List.mem_cons] at ht
    rcases ht with rfl | ht
    · exact hf a (hl a (by simp))
    · exact hl t (by simp [ht])

theorem modifyNxt_vStart {P : Nat → Prop} {l : List TEntry} (f : TEntry → TEntry)
    (hf : ∀ t, P t.vStart → P (f t).vStart) (hl : ∀ t ∈ l, P t.vStart) :
    ∀ t ∈ (match l with | a :: b :: rest => a :: f b :: rest | l => l), P t.vStart := by
  match l with
  | [] | [_] => exact hl
  | a :: b :: rest =>
    intro t ht
    simp only [List.mem_cons] at ht
    rcases ht with rfl | rfl | ht
    · exact hl t (by simp)
    · exact hf b (hl b (by simp))
    · exact hl t (by simp [ht])

namespace WalkState

/-- The typing invariant of the walk state; `ex` are the items whose `vs` may still be unset. -/
structure Typing (g : Graph) (ex : ItemId → Prop) (s : WalkState) : Prop where
  g_eq : s.g = g
  items : Items.Typing s.items g ex
  sv_lt : ∀ d : Nat, s.stackVerts[d]! < g.nv
  ts_vStart : ∀ t ∈ s.tstack, t.vStart < g.nv

namespace Typing

variable {g : Graph} {ex : ItemId → Prop} {s : WalkState}

theorem nv_pos (h : s.Typing g ex) : 0 < g.nv := Nat.lt_of_le_of_lt (Nat.zero_le _) (h.sv_lt 0)

theorem mono (h : s.Typing g ex) {ex' : ItemId → Prop} (hex : ∀ i, ex i → ex' i) : s.Typing g ex' :=
  ⟨h.g_eq, h.items.mono hex, h.sv_lt, h.ts_vStart⟩

theorem cur_vStart (h : s.Typing g ex) : s.tstack.head!.vStart < g.nv :=
  head!_vStart_lt h.nv_pos h.ts_vStart
theorem nxt_vStart (h : s.Typing g ex) : s.tstack.tail.head!.vStart < g.nv :=
  head!_vStart_lt h.nv_pos fun t ht => h.ts_vStart t (List.mem_of_mem_tail ht)

theorem makeVs_lt (h : s.Typing g ex) {a b : Nat} (ha : a < g.nv) (u : Nat)
    (hu : (setSides s.stackDir[b]! (some s.stackVerts[b]!) (some a)).1 = some u ∨
      (setSides s.stackDir[b]! (some s.stackVerts[b]!) (some a)).2 = some u) : u < g.nv := by
  have := h.sv_lt b
  unfold setSides at hu
  split at hu <;> simp at hu <;> rcases hu with rfl | rfl <;> assumption

theorem makeVs_shape (s : WalkState) (a b : Nat) :
    ∃ u v, setSides s.stackDir[b]! (some s.stackVerts[b]!) (some a) = (some u, some v) := by
  unfold setSides; split <;> exact ⟨_, _, rfl⟩

theorem of_eq (h : s.Typing g ex) {s' : WalkState} (hg : s'.g = s.g) (hi : s'.items = s.items)
    (hsv : s'.stackVerts = s.stackVerts) (hts : s'.tstack = s.tstack) : s'.Typing g ex where
  g_eq := by rw [hg]; exact h.g_eq
  items := by rw [hi]; exact h.items
  sv_lt := by rw [hsv]; exact h.sv_lt
  ts_vStart := by rw [hts]; exact h.ts_vStart

theorem tstack (h : s.Typing g ex) {l : List TEntry} (hl : ∀ t ∈ l, t.vStart < g.nv) :
    ({ s with tstack := l }).Typing g ex := ⟨h.g_eq, h.items, h.sv_lt, hl⟩

theorem mergeTop (h : s.Typing g ex) : ({ s with tstack := WalkM.mergeTop s.tstack }).Typing g ex :=
  h.tstack (mergeTop_vStart (P := (· < g.nv)) h.ts_vStart)

theorem tail (h : s.Typing g ex) : ({ s with tstack := s.tstack.tail }).Typing g ex :=
  h.tstack fun t ht => h.ts_vStart t (List.mem_of_mem_tail ht)

theorem cons (h : s.Typing g ex) {t : TEntry} (ht : t.vStart < g.nv) :
    ({ s with tstack := t :: s.tstack }).Typing g ex :=
  h.tstack fun t' ht' => by
    rcases List.mem_cons.mp ht' with rfl | ht'
    · exact ht
    · exact h.ts_vStart t' ht'

theorem modifyCur (h : s.Typing g ex) (f : TEntry → TEntry) (hf : ∀ t, t.vStart < g.nv → (f t).vStart < g.nv) :
    ({ s with tstack := match s.tstack with | a :: rest => f a :: rest | [] => [] }).Typing g ex :=
  h.tstack (modifyCur_vStart (P := (· < g.nv)) f hf h.ts_vStart)

theorem modifyNxt (h : s.Typing g ex) (f : TEntry → TEntry) (hf : ∀ t, t.vStart < g.nv → (f t).vStart < g.nv) :
    ({ s with tstack := match s.tstack with | a :: b :: rest => a :: f b :: rest | l => l }).Typing g ex :=
  h.tstack (modifyNxt_vStart (P := (· < g.nv)) f hf h.ts_vStart)

theorem items_eq (h : s.Typing g ex) {items' : Array Item} {ex' : ItemId → Prop} (hi : Items.Typing items' g ex') :
    ({ s with items := items' }).Typing g ex' := ⟨h.g_eq, hi, h.sv_lt, h.ts_vStart⟩

theorem set_stackDir (h : s.Typing g ex) (a : Array Bool) : ({ s with stackDir := a }).Typing g ex :=
  ⟨h.g_eq, h.items, h.sv_lt, h.ts_vStart⟩

theorem set_stackVerts (h : s.Typing g ex) {d v : Nat} (hv : v < g.nv) :
    ({ s with stackVerts := s.stackVerts.set! d v }).Typing g ex :=
  ⟨h.g_eq, h.items, getElem!_set!_lt hv h.sv_lt, h.ts_vStart⟩

end Typing

end WalkState

namespace WalkM

open WalkState

variable {g : Graph} {ex : ItemId → Prop} {s : WalkState}

theorem bind_spec {P : β → WalkState → Prop} {Q : α → WalkState → Prop} {m : WalkM β} {f : β → WalkM α}
    (h₁ : wp m P s) (h₂ : ∀ b s', P b s' → wp (f b) Q s') : wp (m >>= f) Q s := h₂ _ _ h₁

theorem bind_pure {Q : α → WalkState → Prop} {a : β} {f : β → WalkM α} :
    wp (pure a >>= f) Q s = wp (f a) Q s := rfl
theorem bind_get {Q : α → WalkState → Prop} {f : WalkState → WalkM α} :
    wp (get >>= f) Q s = wp (f s) Q s := rfl
theorem bind_stackDir {Q : α → WalkState → Prop} {d : Nat} {f : Bool → WalkM α} :
    wp (stackDir d >>= f) Q s = wp (f s.stackDir[d]!) Q s := rfl
theorem bind_tstackSize {Q : α → WalkState → Prop} {f : Nat → WalkM α} :
    wp (tstackSize >>= f) Q s = wp (f s.tstack.length) Q s := rfl
theorem bind_cur {Q : α → WalkState → Prop} {f : TEntry → WalkM α} :
    wp (cur >>= f) Q s = wp (f s.tstack.head!) Q s := rfl
theorem bind_nxt {Q : α → WalkState → Prop} {f : TEntry → WalkM α} :
    wp (nxt >>= f) Q s = wp (f s.tstack.tail.head!) Q s := rfl

theorem modifyItem_spec (h : s.Typing g ex) {i : Nat} {f : Item → Item} {ex' : ItemId → Prop}
    (hi : Items.Typing (s.items.modify i f) g ex') :
    wp (modifyItem i f) (fun _ s' => s'.Typing g ex') s := h.items_eq hi
theorem modify_spec {f : WalkState → WalkState} {ex' : ItemId → Prop} (hf : (f s).Typing g ex') :
    wp (modify f) (fun _ s' => s'.Typing g ex') s := hf
theorem allocItem_spec (h : s.Typing g ex) (ty : NodeType) (hty : ty ∉ [NodeType.F, .V, .Q]) :
    wp (allocItem ty) (fun item s' => s'.Typing g (fun i => ex i ∨ i = item) ∧ Items.type s'.items item = ty) s :=
  ⟨h.items_eq (h.items.push ty hty), by
    show Items.type (s.items.push _) s.items.size = ty
    rw [Items.type_push]; simp⟩
theorem makeVs_spec (h : s.Typing g ex) {a : Nat} (ha : a < g.nv) (b : Nat) :
    wp (makeVs a b) (fun vs s' => s' = s ∧ (∃ u v, vs = (some u, some v)) ∧
      ∀ u, vs.1 = some u ∨ vs.2 = some u → u < g.nv) s :=
  ⟨rfl, Typing.makeVs_shape s a b, h.makeVs_lt ha⟩
theorem popTstack_spec : wp popTstack (fun _ s' => s' = { s with tstack := s.tstack.tail }) s := rfl
theorem mergeTstackTops_spec : wp mergeTstackTops (fun _ s' => s' = { s with tstack := mergeTop s.tstack }) s := by
  rw [wp_mergeTstackTops]
theorem setStackDir_spec (d : Nat) (b : Bool) :
    wp (setStackDir d b) (fun _ s' => s' = { s with stackDir := s.stackDir.set! d b }) s := rfl
theorem pushVertTstack_spec (v topDepth : Nat) :
    wp (pushVertTstack v topDepth) (fun _ s' => s' = { s with tstack :=
      ⟨v, topDepth, s.nxtEdgeIdx, setSides s.stackDir[topDepth]! [vertItem v] []⟩ :: s.tstack }) s := rfl
theorem pushEdgeTstack_spec (vStart topDepth e : Nat) :
    wp (pushEdgeTstack vStart topDepth e) (fun _ s' => s' = { s with tstack :=
      ⟨vStart, topDepth, s.nxtEdgeIdx, setSides s.stackDir[topDepth]! [edgeItem s.g e] []⟩ :: s.tstack }) s := rfl
theorem modifyCur_spec (f : TEntry → TEntry) :
    wp (modifyCur f) (fun _ s' => s' = { s with tstack := match s.tstack with | a :: rest => f a :: rest | [] => [] }) s :=
  rfl
theorem map_spec {P : β → WalkState → Prop} {m : WalkM β} (f : β → α) (h : wp m P s) :
    wp (f <$> m) (fun b s' => ∃ a, b = f a ∧ P a s') s := ⟨_, rfl, h⟩

theorem loop_typing {n : Nat} {cond : WalkM Bool} {body : WalkM Unit}
    (hc : ∀ s, s.Typing g ex → wp cond (fun _ s' => s'.Typing g ex) s)
    (hb : ∀ s, s.Typing g ex → wp body (fun _ s' => s'.Typing g ex) s) (h : s.Typing g ex) :
    wp (loop n cond body) (fun _ s' => s'.Typing g ex) s :=
  wp_loop (fun s => s.Typing g ex) n cond body hc hb h fun _ h => h

theorem not_fvq_of_spr {ty : NodeType} (h : ty ∈ [NodeType.S, .P, .R]) : ty ∉ [NodeType.F, .V, .Q] := by
  simp only [List.mem_cons, List.not_mem_nil, or_false] at h
  rcases h with rfl | rfl | rfl <;> simp

theorem maybeUnwrapNxt_spec (h : s.Typing g ex) (ty : NodeType) (hty : ty ∈ [NodeType.S, .P, .R]) :
    wp (maybeUnwrapNxt ty)
      (fun item s' => s'.Typing g (fun i => ex i ∨ i = item) ∧ Items.type s'.items item = ty) s := by
  unfold maybeUnwrapNxt
  simp only [wp_bind, wp_get, wp_ite, wp_allocItem, wp_pure, wp_nxt, wp_stackDir, wp_getItem, wp_modifyNxt]
  split
  · exact ⟨h.items_eq (h.items.push ty (not_fvq_of_spr hty)), by simp [Items.type_push]⟩
  · split
    · next heq =>
      refine ⟨(h.modifyNxt _ ?_).mono fun i hi => Or.inl hi, ?_⟩
      · intro t ht; exact ht
      · rw [← Items.type_getElem!']; simpa using heq
    · exact ⟨h.items_eq (h.items.push ty (not_fvq_of_spr hty)), by simp [Items.type_push]⟩

theorem finishTstackTop_spec (h : s.Typing g ex) (item : Nat)
    (hty : Items.type s.items item ∈ [NodeType.S, .P, .R]) :
    wp (finishTstackTop item) (fun _ s' => s'.Typing g (fun i => ex i ∧ i ≠ item)) s := by
  unfold finishTstackTop
  simp only [wp_bind, wp_cur, wp_stackDir, wp_makeVs, wp_modifyItem, wp_modifyCur]
  exact (h.items_eq (h.items.modify_vs_ch item hty _ _ (Typing.makeVs_shape s _ _)
    (h.makeVs_lt h.cur_vStart))).modifyCur _ fun t ht => ht

theorem finishEdge_spec (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool)
    (h : s.Typing g ex) (hv : curV < g.nv) (hd : o.dest < g.nv) (he : o.e < g.ne) :
    wp (finishEdge curV d o origTstack hasVert)
      (fun _ s' => s'.Typing g (fun i => ex i ∧ i ≠ edgeItem g o.e)) s := by
  unfold finishEdge
  rw [bind_get]
  extract_lets g' nxtV e lowval isTree isType1 qItem jpB isT jp4 jpN jpS isF
  rw [bind_stackDir]
  extract_lets jp2 jpT
  have hq : qItem = edgeItem g o.e := congrArg (edgeItem · o.e) h.g_eq
  have hcur : ∀ u, ((some curV, none) : Option Nat × Option Nat).1 = some u ∨
      ((some curV, none) : Option Nat × Option Nat).2 = some u → u < g.nv := by
    rintro u (hu | hu) <;> simp at hu; exact hu ▸ hv
  have hspr : ∀ b : Bool, (if b = true then NodeType.S else NodeType.R) ∈ [NodeType.S, .P, .R] := by
    intro b; cases b <;> simp
  have hjpB : ∀ s', s'.Typing g (fun i => ex i ∧ i ≠ edgeItem g o.e) →
      wp (jpB ()) (fun _ s' => s'.Typing g (fun i => ex i ∧ i ≠ edgeItem g o.e)) s' := by
    intro s' h'
    dsimp -zeta only [jpB]
    rw [hq]
    exact bind_spec (modifyItem_spec h' (h'.items.modify_v_ch curV o.e hv he)) fun _ _ h'' => h''
  have hjpN : ∀ r b s', s'.Typing g (fun i => ex i ∧ i ≠ edgeItem g o.e) →
      wp (jpN r b) (fun _ s' => s'.Typing g (fun i => ex i ∧ i ≠ edgeItem g o.e)) s' := by
    intro r b s' h'
    dsimp -zeta only [jpN]
    rw [bind_tstackSize, bind_nxt, bind_nxt]
    extract_lets jp5
    have hjp5 : ∀ s'', s''.Typing g (fun i => ex i ∧ i ≠ edgeItem g o.e) →
        wp (jp5 ()) (fun _ s' => s'.Typing g (fun i => ex i ∧ i ≠ edgeItem g o.e)) s'' := by
      intro s'' h''
      dsimp -zeta only [jp5]
      split
      · refine bind_spec (pushVertTstack_spec curV d) ?_
        rintro _ _ rfl
        split
        · refine bind_spec mergeTstackTops_spec ?_
          rintro _ _ rfl
          exact (h''.cons hv).mergeTop
        · exact (h''.cons hv)
      · exact h''
    split
    · refine bind_spec (maybeUnwrapNxt_spec h' .P (by simp)) ?_
      rintro item s₁ ⟨h₁, hty⟩
      refine bind_spec mergeTstackTops_spec ?_
      rintro _ _ rfl
      refine bind_spec (finishTstackTop_spec h₁.mergeTop item (by simp [hty])) fun _ s₂ h₂ => ?_
      exact hjp5 _ (h₂.mono fun i hi => hi.1.resolve_right hi.2)
    · exact hjp5 _ h'
  have hjpS : ∀ ty, ty ∈ [NodeType.S, .P, .R] → ∀ s', s'.Typing g (fun i => ex i ∧ i ≠ edgeItem g o.e) →
      wp (jpS ty) (fun _ s' => s'.Typing g (fun i => ex i ∧ i ≠ edgeItem g o.e)) s' := by
    intro ty hty s' h'
    dsimp -zeta only [jpS]
    refine bind_spec (maybeUnwrapNxt_spec h' ty hty) ?_
    rintro item s₁ ⟨h₁, hty₁⟩
    refine bind_spec mergeTstackTops_spec ?_
    rintro _ _ rfl
    exact wp_mono _ (finishTstackTop_spec h₁.mergeTop item (by simp [hty₁, hty])) fun _ s₂ h₂ =>
      h₂.mono fun i hi => hi.1.resolve_right hi.2
  have hjp2 : ∀ r b s', s'.Typing g (fun i => ex i ∧ i ≠ edgeItem g o.e) →
      wp (jp2 r b) (fun _ s' => s'.Typing g (fun i => ex i ∧ i ≠ edgeItem g o.e)) s' := by
    intro r b s' h'
    dsimp -zeta only [jp2]
    extract_lets jp3
    have hjp3 : ∀ item s'', (match item with
        | some i => s''.Typing g (fun j => (ex j ∧ j ≠ edgeItem g o.e) ∨ j = i) ∧
            Items.type s''.items i ∈ [NodeType.S, .P, .R]
        | none => s''.Typing g (fun i => ex i ∧ i ≠ edgeItem g o.e)) →
        wp (jp3 item) (fun _ s' => s'.Typing g (fun i => ex i ∧ i ≠ edgeItem g o.e)) s'' := by
      intro item s'' hi
      dsimp -zeta only [jp3]
      refine bind_spec mergeTstackTops_spec ?_
      rintro _ _ rfl
      refine bind_spec mergeTstackTops_spec ?_
      rintro _ _ rfl
      refine bind_spec (modifyCur_spec _) ?_
      rintro _ _ rfl
      cases item with
      | some i =>
        obtain ⟨hi, hty⟩ := hi
        refine bind_spec (finishTstackTop_spec (hi.mergeTop.mergeTop.modifyCur _ fun _ _ => hv) i
          (by simp [hty])) fun _ s₃ h₃ => ?_
        exact hjpN () _ _ (h₃.mono fun j hj => hj.1.resolve_right hj.2)
      | none => exact hjpN () _ _ (hi.mergeTop.mergeTop.modifyCur _ fun _ _ => hv)
    split
    · refine bind_spec (map_spec some (maybeUnwrapNxt_spec h' _ (hspr b))) ?_
      rintro _ s₁ ⟨i, rfl, h₁, hty⟩
      exact hjp3 (some i) s₁ ⟨h₁, by rw [hty]; exact hspr b⟩
    · exact hjp3 none s' h'
  have hjpT : ∀ r b s', s'.Typing g (fun i => ex i ∧ i ≠ edgeItem g o.e) →
      wp (jpT r b) (fun _ s' => s'.Typing g (fun i => ex i ∧ i ≠ edgeItem g o.e)) s' := by
    intro r b s' h'
    dsimp -zeta only [jpT]
    split
    · split
      · rw [bind_tstackSize]
        refine bind_spec (loop_typing (fun _ h => h)
          (fun s h => by rw [wp_mergeTstackTops]; exact h.mergeTop) h') fun _ s₁ h₁ => ?_
        exact hjp2 () _ _ h₁
      · exact hjp2 () _ _ h'
    · exact hjpN () _ _ h'
  split
  · rw [hq]
    refine bind_spec (modifyItem_spec h (h.items.modify_vs _ (some curV, none)
      (by rw [h.items.edge o.e he]; exact Or.inl ⟨_, rfl⟩) hcur)) fun _ s₁ h₁ => ?_
    refine bind_spec (modify_spec (h₁.of_eq rfl rfl rfl rfl)) fun _ s₂ h₂ => ?_
    split
    · split
      · refine bind_spec (allocItem_spec h₂ .I (by simp)) ?_
        rintro item s₃ ⟨h₃, hty⟩
        refine bind_spec (makeVs_spec h₃ hd d) ?_
        rintro vs _ ⟨rfl, hsh, hlt⟩
        refine bind_spec (modifyItem_spec h₃ (h₃.items.modify_vs item vs (by rw [hty]; exact hsh) hlt))
          fun _ s₄ h₄ => ?_
        refine bind_spec popTstack_spec ?_
        rintro t _ rfl
        refine bind_spec (modifyItem_spec h₄.tail (h₄.tail.items.modify_q_ch o.e he _)) fun _ s₅ h₅ => ?_
        exact hjpB _ (h₅.mono fun i hi => hi.1.resolve_right hi.2)
      · refine bind_spec popTstack_spec ?_
        rintro backedge _ rfl
        refine bind_spec popTstack_spec ?_
        rintro t _ rfl
        refine bind_spec (modifyItem_spec h₂.tail.tail (h₂.tail.tail.items.modify_q_ch o.e he _))
          fun _ s₃ h₃ => ?_
        exact hjpB _ h₃
    · refine bind_spec (modify_spec (h₂.of_eq rfl rfl rfl rfl)) fun _ s₃ h₃ => ?_
      refine bind_spec (allocItem_spec h₃ .O (by simp)) ?_
      rintro item s₄ ⟨h₄, hty⟩
      refine bind_spec (modifyItem_spec h₄ (h₄.items.modify_vs item (some curV, none)
        (by rw [hty]; exact ⟨_, rfl⟩) hcur)) fun _ s₅ h₅ => ?_
      refine bind_spec (modifyItem_spec h₅ (h₅.items.modify_q_ch o.e he _)) fun _ s₆ h₆ => ?_
      exact hjpB _ (h₆.mono fun i hi => hi.1.resolve_right hi.2)
  · refine bind_spec (makeVs_spec h hd d) ?_
    rintro vs _ ⟨rfl, hsh, hlt⟩
    rw [hq]
    refine bind_spec (modifyItem_spec h (h.items.modify_vs _ vs
      (by rw [h.items.edge o.e he]; exact Or.inr hsh) hlt)) fun _ s₁ h₁ => ?_
    split
    · refine bind_spec (pushEdgeTstack_spec nxtV d e) ?_
      rintro _ _ rfl
      rw [bind_tstackSize]
      refine bind_spec (loop_typing (fun _ h => h) ?_ (h₁.cons hd)) fun _ s₂ h₂ => ?_
      · intro s₂ h₂
        rw [bind_nxt]
        split
        · rw [bind_nxt]
          refine bind_spec (setStackDir_spec _ _) ?_
          rintro _ _ rfl
          refine bind_spec mergeTstackTops_spec ?_
          rintro _ _ rfl
          exact hjpS .S (by simp) _ (h₂.set_stackDir _).mergeTop
        · rw [bind_nxt, bind_cur]
          split
          · exact hjpS .P (by simp) _ h₂
          · exact hjpS .R (by simp) _ h₂
      rw [bind_get]
      extract_lets fo
      rw [bind_cur]
      split
      · rw [bind_tstackSize]
        refine bind_spec (loop_typing (fun _ h => h)
          (fun s h => by rw [wp_mergeTstackTops]; exact h.mergeTop) h₂) fun _ s₃ h₃ => ?_
        exact hjpT () _ _ h₃
      · exact hjpT () _ _ h₂
    · refine bind_spec (pushEdgeTstack_spec curV lowval e) ?_
      rintro _ _ rfl
      refine bind_spec (modify_spec ((h₁.cons hv).of_eq rfl rfl rfl rfl)) fun _ s₂ h₂ => ?_
      exact hjpN () _ _ h₂


end WalkM

mutual
/-- All vertices are `< nv` and all edges `< ne`. -/
def DfsTree.Bounded (nv ne : Nat) : DfsTree → Prop
  | .node v outs => v < nv ∧ DfsOut.BoundedList nv ne outs
def DfsOut.Bounded (nv ne : Nat) : DfsOut → Prop
  | .back e dest _ => e < ne ∧ dest < nv
  | .tree e _ child => e < ne ∧ child.Bounded nv ne
def DfsOut.BoundedList (nv ne : Nat) : List DfsOut → Prop
  | [] => True
  | o :: rest => o.Bounded nv ne ∧ DfsOut.BoundedList nv ne rest
end

theorem DfsOut.edgesList_cons (o : DfsOut) (rest : List DfsOut) :
    DfsOut.edgesList (o :: rest) = DfsOut.edgesList [o] ++ DfsOut.edgesList rest := by
  cases o <;> simp [DfsOut.edgesList]

namespace WalkM

open WalkState

variable {g : Graph} {ex : ItemId → Prop} {s : WalkState}

theorem walk_typing_aux (g : Graph) :
    (∀ (t : DfsTree) (d : Nat) (ex : ItemId → Prop) (s : WalkState), s.Typing g ex →
      t.Bounded g.nv g.ne →
      wp (walkTree t d) (fun _ s' => s'.Typing g fun i => ex i ∧ ∀ e ∈ t.edges, i ≠ edgeItem g e) s) ∧
    (∀ (v d : Nat) (outs : List DfsOut) (hasVert : Bool) (ex : ItemId → Prop) (s : WalkState),
      s.Typing g ex → v < g.nv → DfsOut.BoundedList g.nv g.ne outs →
      wp (walkOuts v d outs hasVert)
        (fun _ s' => s'.Typing g fun i => ex i ∧ ∀ e ∈ DfsOut.edgesList outs, i ≠ edgeItem g e) s) ∧
    (∀ (v d : Nat) (o : DfsOut) (hasVert : Bool) (ex : ItemId → Prop) (s : WalkState),
      s.Typing g ex → v < g.nv → o.Bounded g.nv g.ne →
      wp (walkOut v d o hasVert)
        (fun _ s' => s'.Typing g fun i => ex i ∧ ∀ e ∈ DfsOut.edgesList [o], i ≠ edgeItem g e) s) := by
  refine walkTree.mutual_induct _ _ _ ?_ ?_ ?_ ?_
  · intro d v outs ih ex s h hb
    obtain ⟨hv, hb⟩ := hb
    simp only [walkTree, wp_bind, wp_modify]
    refine wp_mono _ (ih ex _ (h.set_stackVerts hv) hv hb) fun hasVert s₁ h₁ => ?_
    split
    · exact h₁
    · refine bind_spec (setStackDir_spec _ _) ?_
      rintro _ _ rfl
      exact wp_mono _ (pushVertTstack_spec v d) (by rintro _ _ rfl; exact (h₁.set_stackDir _).cons hv)
  · intro o v d hasVert ih ex s h hv hb
    cases o with
    | back e dest cls =>
      obtain ⟨he, hd⟩ := hb
      have hfin : ∀ (hasVert₁ : Bool) (origTstack : Nat) (s₁ : WalkState), s₁.Typing g ex →
          wp (finishEdge v d (.back e dest cls) origTstack hasVert₁)
            (fun _ s' => s'.Typing g fun i => ex i ∧
              ∀ e' ∈ DfsOut.edgesList [.back e dest cls], i ≠ edgeItem g e') s₁ :=
        fun hasVert₁ origTstack s₁ h₁ =>
          wp_mono _ (finishEdge_spec v d (.back e dest cls) origTstack hasVert₁ h₁ hv hd he)
            fun _ s₂ h₂ => h₂.mono fun i hi => ⟨hi.1, by simpa [DfsOut.edgesList, DfsOut.e] using hi.2⟩
      unfold walkOut
      dsimp only
      rw [bind_stackDir]
      refine bind_spec (setStackDir_spec _ _) ?_
      rintro _ _ rfl
      split
      · refine bind_spec (pushVertTstack_spec v d) ?_
        rintro _ _ rfl
        rw [bind_pure, bind_tstackSize]
        exact hfin _ _ _ ((h.set_stackDir _).cons hv)
      · rw [bind_pure, bind_tstackSize]
        exact hfin _ _ _ (h.set_stackDir _)
    | tree e cls child =>
      obtain ⟨he, hb⟩ := hb
      obtain ⟨cv, couts⟩ := child
      have hcv : cv < g.nv := hb.1
      dsimp only at ih
      have hfin : ∀ (hasVert₁ : Bool) (origTstack : Nat) (s₁ : WalkState), s₁.Typing g ex →
          wp ((modify fun s => { s with firstOccurrence := s.firstOccurrence.set! d s.g.ne }) >>= fun _ =>
              walkTree (.node cv couts) (d + 1) >>= fun _ =>
              finishEdge v d (.tree e cls (.node cv couts)) origTstack hasVert₁)
            (fun _ s' => s'.Typing g fun i => ex i ∧
              ∀ e' ∈ DfsOut.edgesList [.tree e cls (.node cv couts)], i ≠ edgeItem g e') s₁ := by
        intro hasVert₁ origTstack s₁ h₁
        refine bind_spec (modify_spec (h₁.of_eq rfl rfl rfl rfl)) fun _ s₂ h₂ => ?_
        refine bind_spec (ih ex _ h₂ hb) fun _ s₃ h₃ => ?_
        exact wp_mono _ (finishEdge_spec v d (.tree e cls (.node cv couts)) _ _ h₃ hv hcv he)
          fun _ s₄ h₄ => h₄.mono fun i hi => ⟨hi.1.1, by
            simp only [DfsOut.edgesList, List.append_nil, List.forall_mem_cons]
            exact ⟨hi.2, hi.1.2⟩⟩
      unfold walkOut
      dsimp only
      rw [bind_stackDir]
      refine bind_spec (setStackDir_spec _ _) ?_
      rintro _ _ rfl
      split
      · refine bind_spec (pushVertTstack_spec v d) ?_
        rintro _ _ rfl
        rw [bind_pure, bind_tstackSize]
        exact hfin _ _ _ ((h.set_stackDir _).cons hv)
      · rw [bind_pure, bind_tstackSize]
        exact hfin _ _ _ (h.set_stackDir _)
  · intro v d hasVert ex s h hv hb
    exact h.mono fun i hi => ⟨hi, by simp [DfsOut.edgesList]⟩
  · intro v d hasVert o rest ih₁ ih₂ ex s h hv hb
    simp only [walkOuts]
    refine bind_spec (ih₁ ex s h hv hb.1) fun hasVert₁ s₁ h₁ => ?_
    exact wp_mono _ (ih₂ hasVert₁ _ s₁ h₁ hv hb.2) fun _ s₂ h₂ => h₂.mono fun i hi => by
      rw [DfsOut.edgesList_cons]
      simp only [List.forall_mem_append]
      exact ⟨hi.1.1, hi.1.2, hi.2⟩

theorem walkTree_typing (t : DfsTree) (d : Nat) (h : s.Typing g ex) (hb : t.Bounded g.nv g.ne) :
    wp (walkTree t d) (fun _ s' => s'.Typing g fun i => ex i ∧ ∀ e ∈ t.edges, i ≠ edgeItem g e) s :=
  (walk_typing_aux g).1 t d ex s h hb

theorem walkForest_typing (forest : List DfsTree) (h : s.Typing g ex)
    (hb : ∀ t ∈ forest, t.Bounded g.nv g.ne) :
    wp (walkForest forest)
      (fun _ s' => s'.Typing g fun i => ex i ∧ ∀ t ∈ forest, ∀ e ∈ t.edges, i ≠ edgeItem g e) s := by
  induction forest generalizing ex s with
  | nil => exact h.mono fun i hi => ⟨hi, by simp⟩
  | cons t rest ih =>
    show wp ((walkTree t 0 >>= fun _ => popTstack >>= fun top =>
      modifyItem rootItem fun it => { it with ch := it.ch ++ top.spans.2 }) >>= fun _ => walkForest rest) _ s
    refine bind_spec (P := fun _ s' => s'.Typing g fun i => ex i ∧ ∀ e ∈ t.edges, i ≠ edgeItem g e) ?_
      fun _ s₁ h₁ => ?_
    · refine bind_spec (walkTree_typing t 0 h (hb t (by simp))) fun _ s₁ h₁ => ?_
      refine bind_spec popTstack_spec ?_
      rintro top _ rfl
      exact h₁.tail.items_eq (h₁.tail.items.modify_root_ch _)
    · exact wp_mono _ (ih h₁ fun t ht => hb t (by simp [ht])) fun _ s₂ h₂ => h₂.mono fun i hi => by
        simp only [List.forall_mem_cons]
        exact ⟨hi.1.1, hi.1.2, hi.2⟩

end WalkM

theorem Items.initial_typing (g : Graph) :
    Items.Typing (initialItems g) g fun i => ∃ e, e < g.ne ∧ i = edgeItem g e where
  size := Nat.le_of_eq (Items.initialItems_size g).symm
  root := by rw [Items.initialItems_type]; rfl
  vert v hv := by
    have h1 : (1 + v : Nat) < 1 + g.nv := by omega
    rw [Items.initialItems_type]; simp [vertItem, h1]
  edge e he := by
    have h1 : ¬ (1 + g.nv + e : Nat) < 1 + g.nv := by omega
    have h2 : (1 + g.nv + e : Nat) < 1 + g.nv + g.ne := by omega
    rw [Items.initialItems_type]; simp [edgeItem, h1, h2]
  node i hi hi' := by rw [Items.initialItems_size] at hi'; omega
  v_children v c hv hc := by simp [Items.IsParent, Items.initialItems_ch] at hc
  i_o_leaf i _ _ := Items.initialItems_ch g i
  vs_shape i hi hex := by
    rw [Items.initialItems_vs, Items.initialItems_type]
    split_ifs with h0 h1 h2
    · rfl
    · rfl
    · refine absurd ⟨i - (1 + g.nv), ?_, ?_⟩ hex
      · show (i - (1 + g.nv) : Nat) < g.ne; omega
      · show (i : Nat) = 1 + g.nv + (i - (1 + g.nv)); omega
    · rfl
  vs_lt i u hu := by rw [Items.initialItems_vs] at hu; simp at hu

/-- The initial state is typed, with every `Q` item still unset. -/
theorem WalkState.init_typing (g : Graph) (ternarize : Bool) (hnv : 0 < g.nv) :
    (WalkState.init g ternarize).Typing g fun i => ∃ e, e < g.ne ∧ i = edgeItem g e where
  g_eq := rfl
  items := Items.initial_typing g
  sv_lt d := by
    show (Array.replicate g.nv 0)[d]! < g.nv
    by_cases hd : d < g.nv
    · simpa [getElem!_pos, hd] using hnv
    · simpa [getElem!_neg, hd] using hnv
  ts_vStart _ ht := nomatch ht

/-- The typing / allocation part of `Items.WF` (`Items.Tree.size/root/vert/edge/node/v_children`,
`Items.Shapes.i_o_leaf`, `Items.Endpoints.vs_shape/vs_lt`). -/
structure WalkTyping (g : Graph) (items : Items) : Prop where
  size : 1 + g.nv + g.ne ≤ items.size
  root : items.type rootItem = .F
  vert : ∀ v, v < g.nv → items.type (vertItem v) = .V
  edge : ∀ e, e < g.ne → items.type (edgeItem g e) = .Q
  node : ∀ i, 1 + g.nv + g.ne ≤ i → i < items.size → items.type i ∉ [NodeType.F, .V, .Q]
  v_children : ∀ v c, v < g.nv → items.IsParent (vertItem v) c → items.type c = .Q
  i_o_leaf : ∀ i, i < items.size → items.type i = .I ∨ items.type i = .O → items.ch i = []
  vs_shape : ∀ i, i < items.size →
    match items.type i with
    | .F | .V => items.vs i = (none, none)
    | .O => ∃ v, items.vs i = (some v, none)
    | .Q => (∃ v, items.vs i = (some v, none)) ∨ (∃ u v, items.vs i = (some u, some v))
    | _ => ∃ u v, items.vs i = (some u, some v)
  vs_lt : ∀ i u, (items.vs i).1 = some u ∨ (items.vs i).2 = some u → u < g.nv

theorem walk_typing (g : Graph) (ternarize : Bool) (forest : List DfsTree) (hnv : 0 < g.nv)
    (hb : ∀ t ∈ forest, t.Bounded g.nv g.ne) (hcov : ∀ e, e < g.ne → ∃ t ∈ forest, e ∈ t.edges) :
    WalkTyping g (g.walk ternarize forest).items := by
  have h : (g.walk ternarize forest).Typing g fun i =>
      (∃ e, e < g.ne ∧ i = edgeItem g e) ∧ ∀ t ∈ forest, ∀ e ∈ t.edges, i ≠ edgeItem g e :=
    WalkM.walkForest_typing forest (WalkState.init_typing g ternarize hnv) hb
  have hex : ∀ i, ¬ ((∃ e, e < g.ne ∧ i = edgeItem g e) ∧
      ∀ t ∈ forest, ∀ e ∈ t.edges, i ≠ edgeItem g e) := by
    rintro i ⟨⟨e, he, rfl⟩, hall⟩
    obtain ⟨t, ht, het⟩ := hcov e he
    exact hall t ht e het rfl
  exact
    { h.items with
      vs_shape := fun i hi => by
        have := h.items.vs_shape i hi (hex i)
        unfold VsShape at this
        generalize Items.type (g.walk ternarize forest).items i = ty at this ⊢
        cases ty <;> exact this }

/-- `Items.Shapes.q_children`. A `Q` item's `ch` is written only in the block branch of
`finishEdge`, to `item :: t.spans.2`, `backedge.spans.1 ++ t.spans.2`, or `[item]`; that these have
the shape `[c]` / `[c, vertItem v]` needs the tstack span invariants (PROOF.md §4), not proved here. -/
theorem walk_q_children (g : Graph) (ternarize : Bool) (forest : List DfsTree) (hnv : 0 < g.nv)
    (hb : ∀ t ∈ forest, t.Bounded g.nv g.ne) (hcov : ∀ e, e < g.ne → ∃ t ∈ forest, e ∈ t.edges) :
    ∀ e, e < g.ne → Items.ch (g.walk ternarize forest).items (edgeItem g e) = [] ∨
      ∃ c, Items.type (g.walk ternarize forest).items c ∉ [NodeType.F, .V, .Q] ∧
        (Items.ch (g.walk ternarize forest).items (edgeItem g e) = [c] ∨
          ∃ v, Items.ch (g.walk ternarize forest).items (edgeItem g e) = [c, vertItem v]) := by
  sorry

end Spqr
