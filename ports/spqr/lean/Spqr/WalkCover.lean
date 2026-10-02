import Spqr.Frame
import Spqr.WalkPlace
import Spqr.StWalk

/-!
# Nothing is dropped

`finishTstackTop`, `maybeUnwrapNxt`, the block branch of `finishEdge` and `walkForest` keep one side
of the ear(s) they close and discard the other; `mergeTstackTops` on fewer than two entries drops
the top. `Sides*` below mirrors the walk (in the style of `GuardsTree`) and asserts, at every such
site, that the discarded side is empty (`TEntry.OnSide`) resp. that the entries exist. Under it the
count invariant `Place` is exact (`Full`): every allocated node not loose, and every pushed
`vertItem`/`edgeItem`, occurs exactly once — which gives the "nothing dropped" half of the tree
shape of `Items.WF`.
-/

namespace Spqr
open WalkM WalkState

/-! ### Site conditions -/

/-- `mergeTstackTops` has two entries to merge. -/
def MergeOK (s : WalkState) : Prop := 2 ≤ s.tstack.length

/-- `finishTstackTop`: the top entry exists and the side it discards is empty. -/
def CloseOK (s : WalkState) : Prop :=
  ∃ t rest, s.tstack = t :: rest ∧ t.OnSide s.stackDir[t.topDepth]!

/-- `maybeUnwrapNxt ty` (when it may reuse): if the head item of the entry below the top has type
`ty`, that entry is this single item, on its own side. -/
def UnwrapOK (s : WalkState) (ty : NodeType) : Prop :=
  (ty == .R || s.ternarize) = false →
    ∀ a t rest, s.tstack = a :: t :: rest →
      ∀ x xs, getSide t.spans s.stackDir[t.topDepth]! = x :: xs → Items.type s.items x = ty →
        xs = [] ∧ t.OnSide s.stackDir[t.topDepth]!

/-- The block branch of `finishEdge`: the popped entries' discarded sides are empty. -/
def BoundaryOK (d : Nat) (o : DfsOut) (s : WalkState) : Prop :=
  o.cls.isTree = true →
    if o.cls.lowval d == d + 1 then ∀ t ∈ s.tstack.head?, t.spans.1 = []
    else (∀ b ∈ s.tstack.head?, b.spans.2 = []) ∧ ∀ t ∈ s.tstack.tail.head?, t.spans.1 = []

/-- `walkForest`, after `walkTree t 0`: one entry is left, with everything on its second side. -/
def RootOK (s : WalkState) : Prop := s.tstack.tail = [] ∧ ∀ t ∈ s.tstack.head?, t.spans.1 = []

/-! ### The sites, following the blocks of `finishEdge'` -/

/-- Site conditions `B` at every iteration of `loop n cond body`. -/
def LoopSides (B : WalkState → Prop) (cond : WalkM Bool) (body : WalkM Unit) :
    Nat → WalkState → Prop
  | 0, _ => True
  | n + 1, s => wp cond (fun b s₁ => b = true →
      B s₁ ∧ wp body (fun _ s₂ => LoopSides B cond body n s₂) s₁) s

/-- `maybeUnwrapNxt ty; mergeTstackTops; finishTstackTop`. -/
def UnwrapMergeClose (ty : NodeType) (s : WalkState) : Prop :=
  UnwrapOK s ty ∧ wp (maybeUnwrapNxt ty) (fun _ s₁ =>
    MergeOK s₁ ∧ wp mergeTstackTops (fun _ s₂ => CloseOK s₂) s₁) s

def Loop1Sides (d : Nat) (edgeDir : Bool) (s : WalkState) : Prop :=
  (d < s.tstack.tail.head!.topDepth → MergeOK s) ∧
  wp (loop1Type d edgeDir) (fun ty s₁ => UnwrapMergeClose ty s₁) s

def FinishPSides (curV lowval : Nat) (isType1 : Bool) (s : WalkState) : Prop :=
  wp (condP curV lowval isType1) (fun b s₁ => b = true → UnwrapMergeClose .P s₁) s

def FinishTailSides (curV d : Nat) (hasVert isSingle : Bool) (s : WalkState) : Prop :=
  hasVert = false → isSingle = false → wp (pushVertTstack curV d) (fun _ s₁ => MergeOK s₁) s

def FinishRestSides (curV d lowval : Nat) (isType1 hasVert isSingle : Bool) (s : WalkState) : Prop :=
  FinishPSides curV lowval isType1 s ∧
  wp (finishP curV lowval isType1) (fun _ s₁ => FinishTailSides curV d hasVert isSingle s₁) s

def CloseVertTailSides (curV : Nat) (edgeDir : Bool) (item : Option ItemId) (K : WalkState → Prop)
    (s : WalkState) : Prop :=
  MergeOK s ∧ wp mergeTstackTops (fun _ s₁ => MergeOK s₁ ∧ wp mergeTstackTops (fun _ s₂ =>
    wp (modifyCur fun t => { t with vStart := curV, spans := setSides (!edgeDir) (t.spans.1 ++ t.spans.2) [] })
      (fun _ s₃ => match item with
        | some item => CloseOK s₃ ∧ wp (finishTstackTop item) (fun _ s₄ => K s₄) s₃
        | none => K s₃) s₂) s₁) s

def CloseVertSides (curV : Nat) (edgeDir isType1 : Bool) (origTstack : Nat) (isSingle : Bool)
    (K : Bool → WalkState → Prop) (s : WalkState) : Prop :=
  (isType1 = false →
    LoopSides MergeOK (loop3Cond origTstack) mergeTstackTops s.tstack.length s ∧
    wp (loop s.tstack.length (loop3Cond origTstack) mergeTstackTops)
      (fun _ s₁ => CloseVertTailSides curV edgeDir none (K false) s₁) s) ∧
  (isType1 = true →
    UnwrapOK s (if isSingle then .S else .R) ∧
    wp (maybeUnwrapNxt (if isSingle then .S else .R))
      (fun item s₁ => CloseVertTailSides curV edgeDir (some item) (K true) s₁) s)

def MergeLateSides (d : Nat) (s : WalkState) : Prop :=
  s.firstOccurrence[d]! < s.tstack.head!.firstIdx →
    LoopSides MergeOK (loop2Cond s.firstOccurrence[d]!) mergeTstackTops s.tstack.length s

def FinishTreeSides (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert edgeDir : Bool)
    (s : WalkState) : Prop :=
  wp (pushEdgeTstack o.dest d o.e) (fun _ s₁ =>
    LoopSides (Loop1Sides d edgeDir) (loop1Cond d) (loop1Body d edgeDir) s₁.tstack.length s₁) s ∧
  wp (closeEars o.dest d o.e edgeDir) (fun _ s₂ =>
    MergeLateSides d s₂ ∧
    wp (mergeLate d) (fun isSingle s₃ =>
      (hasVert = true → CloseVertSides curV edgeDir o.cls.isType1 origTstack isSingle
        (fun b => FinishRestSides curV d (o.cls.lowval d) o.cls.isType1 hasVert b) s₃) ∧
      (hasVert = false →
        FinishRestSides curV d (o.cls.lowval d) o.cls.isType1 hasVert isSingle s₃)) s₂) s

def FinishBackSides (curV d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  wp (pushEdgeTstack curV (o.cls.lowval d) o.e) (fun _ s₁ =>
    wp (modify fun s =>
      { s with firstOccurrence := s.firstOccurrence.modify (o.cls.lowval d) (min · s.nxtEdgeIdx),
               nxtEdgeIdx := s.nxtEdgeIdx + 1 })
      (fun _ s₂ => FinishRestSides curV d (o.cls.lowval d) o.cls.isType1 hasVert true s₂) s₁) s

/-- The side conditions of one `finishEdge curV d o origTstack hasVert` call, started at `s`. -/
def FinishSides (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool) (s : WalkState) :
    Prop :=
  (d ≤ o.cls.lowval d → BoundaryOK d o s) ∧
  (o.cls.lowval d < d → wp (finishSetup d o) (fun edgeDir s₁ =>
    (o.cls.isTree = true → FinishTreeSides curV d o origTstack hasVert edgeDir s₁) ∧
    (o.cls.isTree = false → FinishBackSides curV d o hasVert s₁)) s)

mutual
/-- `FinishSides` at every `finishEdge` call of the walk of `t` from `s` (mirrors `GuardsTree`). -/
def SidesTree (t : DfsTree) (d : Nat) (s : WalkState) : Prop :=
  match t with
  | .node v outs => SidesOuts v d outs false { s with stackVerts := s.stackVerts.set! d v }

def SidesOuts (v d : Nat) (outs : List DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  match outs with
  | [] => True
  | o :: rest => SidesOut v d o hasVert s ∧
      wp (walkOut v d o hasVert) (fun hasVert' s' => SidesOuts v d rest hasVert' s') s

def SidesOut (v d : Nat) (o : DfsOut) (hasVert : Bool) (s : WalkState) : Prop :=
  wp (walkOutPre v d o hasVert) (fun hasVert' s₁ =>
    match o with
    | .tree _ _ child =>
      wp (modify fun s => { s with firstOccurrence := s.firstOccurrence.set! d s.g.ne }) (fun _ s₂ =>
        SidesTree child (d + 1) s₂ ∧
        wp (walkTree child (d + 1)) (fun _ s₃ => FinishSides v d o s₁.tstack.length hasVert' s₃) s₂) s₁
    | .back .. => FinishSides v d o s₁.tstack.length hasVert' s₁) s
end

/-- `SidesTree` for every tree of the forest, and `RootOK` after each. -/
def SidesForest : List DfsTree → WalkState → Prop
  | [], _ => True
  | t :: rest, s => SidesTree t 0 s ∧
      wp (walkTree t 0) (fun _ s₁ => RootOK s₁ ∧
        wp (popTstack >>= fun top => modifyItem rootItem fun it => { it with ch := it.ch ++ top.spans.2 })
          (fun _ s₂ => SidesForest rest s₂) s₁) s

/-! ### Exact placement -/

namespace WalkState

/-- `Place`, exactly: pushed fixed items and non-loose allocated nodes are placed (once), and an
unfinished edge's Q item has no children yet. -/
structure Full (g : Graph) (P X : ItemId → Prop) (s : WalkState) : Prop where
  place : s.Place g P X
  fixedP : ∀ i, P i → 0 < i ∧ i < 1 + g.nv + g.ne
  pushed : ∀ i, P i → 0 < s.cnt i
  placed : ∀ i, 1 + g.nv + g.ne ≤ i → i < s.items.size → ¬ X i → 0 < s.cnt i
  qch : ∀ e, e < g.ne → ¬ P (edgeItem g e) → Items.ch s.items (edgeItem g e) = []

namespace Full

variable {g : Graph} {P P' X X' : ItemId → Prop} {s s' : WalkState}

theorem mono (h : s.Full g P X) (hP : ∀ i, P i ↔ P' i) (hX : ∀ i, X i ↔ X' i) : s.Full g P' X' where
  place := h.place.mono (fun i => (hP i).1) fun i => (hX i).2
  fixedP i hi := h.fixedP i ((hP i).2 hi)
  pushed i hi := h.pushed i ((hP i).2 hi)
  placed i hn hi hx := h.placed i hn hi fun h' => hx ((hX i).1 h')
  qch e he hP' := h.qch e he fun h' => hP' ((hP _).1 h')

theorem monoP (h : s.Full g P X) (hP : ∀ i, P i ↔ P' i) : s.Full g P' X := h.mono hP fun _ => Iff.rfl
theorem monoX (h : s.Full g P X) (hX : ∀ i, X i ↔ X' i) : s.Full g P X' := h.mono (fun _ => Iff.rfl) hX

theorem edgeItem_g (h : s.Full g P X) (e : Nat) : edgeItem s.g e = edgeItem g e := h.place.edgeItem_g e

/-- A step that places nothing new and changes no child list of an unfinished edge. -/
theorem of_cnt (h : s.Full g P X) (hp : s'.Place g P X) (hsz : s'.items.size = s.items.size)
    (hq : ∀ e, e < g.ne → ¬ P (edgeItem g e) → Items.ch s'.items (edgeItem g e) = Items.ch s.items (edgeItem g e))
    (hc : ∀ i, s'.cnt i = s.cnt i) : s'.Full g P X where
  place := hp
  fixedP := h.fixedP
  pushed i hi := hc i ▸ h.pushed i hi
  placed i hn hi hx := hc i ▸ h.placed i hn (hsz ▸ hi) hx
  qch e he hP := (hq e he hP).trans (h.qch e he hP)

theorem of_eq (h : s.Full g P X) (hg : s'.g = s.g) (hi : s'.items = s.items) (ht : s'.tstack = s.tstack) :
    s'.Full g P X :=
  h.of_cnt (h.place.of_eq hg hi ht) (by rw [hi]) (fun _ _ _ => by rw [hi]) fun i => by simp [cnt, hi, ht]

theorem tstack (h : s.Full g P X) {ts : List TEntry} (hts : ∀ i, spansCount ts i = spansCount s.tstack i) :
    ({ s with tstack := ts }).Full g P X :=
  h.of_cnt (h.place.tstack_le fun i => Nat.le_of_eq (hts i)) rfl (fun _ _ _ => rfl)
    fun i => by simp only [cnt, hts]

/-- Placing the node `x` (loose before). -/
theorem step_node (h : s.Full g P X) {x : Nat}
    (hp : s'.Place g P (fun i => X i ∧ i ≠ x)) (hsz : s'.items.size = s.items.size)
    (hq : ∀ e, e < g.ne → ¬ P (edgeItem g e) → Items.ch s'.items (edgeItem g e) = Items.ch s.items (edgeItem g e))
    (hc : ∀ i, s'.cnt i = s.cnt i + if i = x then 1 else 0) : s'.Full g P (fun i => X i ∧ i ≠ x) where
  place := hp
  fixedP := h.fixedP
  pushed i hi := by have := h.pushed i hi; rw [hc]; omega
  placed i hn hi hx' := by
    rw [hc]
    by_cases hix : i = x
    · simp [hix]
    · have := h.placed i hn (hsz ▸ hi) fun h' => hx' ⟨h', hix⟩; omega
  qch e he hP := (hq e he hP).trans (h.qch e he hP)

/-- Placing the fixed item `x` (not pushed before). -/
theorem step_fixed (h : s.Full g P X) {x : Nat} (hx0 : 0 < x) (hxlt : x < 1 + g.nv + g.ne)
    (hp : s'.Place g (fun i => P i ∨ i = x) X) (hsz : s'.items.size = s.items.size)
    (hq : ∀ e, e < g.ne → ¬ (P (edgeItem g e) ∨ edgeItem g e = x) →
      Items.ch s'.items (edgeItem g e) = Items.ch s.items (edgeItem g e))
    (hc : ∀ i, s'.cnt i = s.cnt i + if i = x then 1 else 0) : s'.Full g (fun i => P i ∨ i = x) X where
  place := hp
  fixedP i hi := hi.elim (h.fixedP i) fun h' => h' ▸ ⟨hx0, hxlt⟩
  pushed i hi := by
    rw [hc]
    rcases hi with hi | rfl
    · have := h.pushed i hi; omega
    · simp
  placed i hn hi hx' := by have := h.placed i hn (hsz ▸ hi) hx'; rw [hc]; omega
  qch e he hP := (hq e he hP).trans (h.qch e he fun h' => hP (Or.inl h'))

theorem set_stackDir (h : s.Full g P X) (a : Array Bool) : ({ s with stackDir := a }).Full g P X :=
  h.of_eq rfl rfl rfl
theorem set_stackVerts (h : s.Full g P X) (a : Array Nat) : ({ s with stackVerts := a }).Full g P X :=
  h.of_eq rfl rfl rfl

theorem push (h : s.Full g P X) (ty : NodeType) :
    ({ s with items := s.items.push ⟨ty, (none, none), []⟩ }).Full g P
      (fun i => X i ∨ i = s.items.size) where
  place := h.place.push ty
  fixedP := h.fixedP
  pushed i hi := by
    have := h.pushed i hi
    simp only [cnt] at this ⊢; rwa [Items.chCount_push _ _ rfl]
  placed i hn hi hx := by
    rw [Array.size_push] at hi
    rw [not_or] at hx
    have := h.placed i hn (Nat.lt_of_le_of_ne (Nat.le_of_lt_succ hi) hx.2) hx.1
    simp only [cnt] at this ⊢; rwa [Items.chCount_push _ _ rfl]
  qch e he hP := by
    show Items.ch (s.items.push _) _ = []
    rw [Items.ch_push_nil (hx := rfl)]; exact h.qch e he hP

theorem modify_vs (h : s.Full g P X) (x : ItemId) (vs : Option Nat × Option Nat) :
    ({ s with items := s.items.modify x fun it => { it with vs := vs } }).Full g P X :=
  h.of_cnt (h.place.modify_vs x vs) (by simp) (fun _ _ _ => Items.ch_modify_of_ch _ _ _ _ fun _ => rfl)
    fun i => by
      simp only [cnt]
      rw [Items.chCount_modify_of_ch _ x (fun it => { it with vs := vs }) (fun _ => rfl)]

theorem mergeTop (h : s.Full g P X) (hm : MergeOK s) :
    ({ s with tstack := WalkM.mergeTop s.tstack }).Full g P X := by
  refine h.tstack fun i => ?_
  rcases hts : s.tstack with _ | ⟨b, _ | ⟨a, rest⟩⟩
  · simp [MergeOK, hts] at hm
  · simp [MergeOK, hts] at hm
  · simp only [WalkM.mergeTop, spansCount_cons, List.count_append]; omega

theorem modifyCur_reorder (h : s.Full g P X) (f : TEntry → TEntry)
    (hf : ∀ t i, ((f t).spans.1 ++ (f t).spans.2).count i = (t.spans.1 ++ t.spans.2).count i) :
    ({ s with tstack := match s.tstack with | a :: rest => f a :: rest | [] => [] }).Full g P X := by
  refine h.tstack fun i => ?_
  cases s.tstack with
  | nil => rfl
  | cons a rest => simp only [spansCount_cons, hf]

theorem fold (h : s.Full g P X) (curV : Nat) (edgeDir : Bool) :
    ({ s with tstack := match s.tstack with
      | a :: rest => { a with vStart := curV, spans := setSides (!edgeDir) (a.spans.1 ++ a.spans.2) [] } :: rest
      | [] => [] }).Full g P X :=
  h.modifyCur_reorder _ fun t i => by simp only [count_setSides, List.append_nil]

theorem cons (h : s.Full g P X) {x : Nat} (hx0 : 0 < x) (hxlt : x < 1 + g.nv + g.ne) (hP : ¬ P x)
    (vS tD fI : Nat) (dir : Bool) :
    ({ s with tstack := ⟨vS, tD, fI, setSides dir [x] []⟩ :: s.tstack }).Full g (fun i => P i ∨ i = x) X :=
  h.step_fixed hx0 hxlt (h.place.cons_fixed hx0 hxlt hP vS tD fI dir) rfl (fun _ _ _ => rfl) fun i => by
    simp only [cnt, spansCount_cons, count_setSides, List.append_nil]
    by_cases hix : i = x
    · subst hix; simp only [List.count_singleton, beq_self_eq_true, ite_true]; omega
    · simp [Ne.symm hix, hix]

theorem count_getSide_of_onSide (p : List ItemId × List ItemId) (dir : Bool) (h : getSide p (!dir) = [])
    (i : ItemId) : (getSide p dir).count i = (p.1 ++ p.2).count i := by
  cases dir <;> simp_all [getSide]

theorem dropCh_cnt (s : WalkState) {x : Nat} (hx : x < s.items.size) (i : ItemId) :
    (s.dropCh x).cnt i + (Items.ch s.items x).count i = s.cnt i := by
  have := Items.chCount_modify s.items hx (fun it => { it with ch := [] }) i
  simp only [List.count_nil, Nat.add_zero] at this
  simp only [cnt, WalkState.dropCh]; omega

theorem dropCh_fresh (h : s.Full g P X) {x : Nat} (hlt : x < s.items.size) (hx : Items.ch s.items x = []) :
    (s.dropCh x).Full g P X := by
  refine h.of_cnt (h.place.dropCh x) (by simp [WalkState.dropCh]) (fun e _ _ => ?_) fun i => ?_
  · simp only [WalkState.dropCh]
    by_cases hxe : x = edgeItem g e
    · subst hxe; rw [Items.ch_modify_self _ _ _ hlt, hx]
    · rw [Items.ch_modify_ne _ _ _ _ hxe]
  · have := dropCh_cnt s hlt i; rw [hx] at this; simpa using this

end Full
end WalkState
end Spqr
