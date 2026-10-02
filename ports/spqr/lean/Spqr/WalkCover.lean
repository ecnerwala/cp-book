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

theorem head!_spans1_nil {ts : List TEntry} (hb : ∀ t ∈ ts.head?, t.spans.1 = []) : ts.head!.spans.1 = [] := by
  cases ts with
  | nil => rfl
  | cons t ts => exact hb t rfl
theorem head!_spans2_nil {ts : List TEntry} (hb : ∀ t ∈ ts.head?, t.spans.2 = []) : ts.head!.spans.2 = [] := by
  cases ts with
  | nil => rfl
  | cons t ts => exact hb t rfl

theorem node_ne_edgeItem {x : Nat} (hn : 1 + g.nv + g.ne ≤ x) {e : Nat} (he : e < g.ne) : x ≠ edgeItem g e := by
  show x ≠ 1 + g.nv + e; omega

/-- `finishTstackTop x` with the discarded side of the top entry empty. -/
theorem finishTop {x : Nat} (h : (s.dropCh x).Full g P X) (hx : X x) {a : TEntry} {rest : List TEntry}
    (hts : s.tstack = a :: rest) {dir : Bool} (hside : a.OnSide dir) (vs : Option Nat × Option Nat) :
    ({ s with
        items := s.items.modify x fun it => { it with vs := vs, ch := getSide a.spans dir },
        tstack := { a with spans := setSides dir [x] [] } :: rest }).Full g P (fun i => X i ∧ i ≠ x) := by
  have hn := h.place.loose_node x hx
  have hlt : x < s.items.size := by simpa [WalkState.dropCh] using hn.2
  have hp := h.place.finishTop hx dir vs
  rw [hts] at hp
  refine h.step_node hp (by simp [WalkState.dropCh]) (fun e he _ => ?_) fun i => ?_
  · simp only [WalkState.dropCh]
    rw [Items.ch_modify_ne _ _ _ _ (node_ne_edgeItem hn.1 he), Items.ch_modify_ne _ _ _ _ (node_ne_edgeItem hn.1 he)]
  · simp only [cnt, WalkState.dropCh, hts, spansCount_cons, count_setSides, List.append_nil]
    have h1 := Items.chCount_modify s.items hlt (fun it => { it with vs := vs, ch := getSide a.spans dir }) i
    have h2 := Items.chCount_modify s.items hlt (fun it => { it with ch := [] }) i
    simp only [List.count_nil, Nat.add_zero] at h1 h2
    have h3 := count_getSide_of_onSide a.spans dir hside i
    by_cases hix : i = x
    · subst hix; simp only [List.count_singleton, beq_self_eq_true, ite_true]; omega
    · simp only [List.count_singleton, hix, Ne.symm hix, beq_iff_eq, ite_false]; omega

/-- `maybeUnwrapNxt` reusing `x`, the single item of `nxt`. -/
theorem unwrap (h : s.Full g P X) {a t : TEntry} {rest : List TEntry} (hts : s.tstack = a :: t :: rest)
    {x : Nat} {dir : Bool} (hx : getSide t.spans dir = [x]) (hside : t.OnSide dir)
    (hn : 1 + g.nv + g.ne ≤ x) (hlt : x < s.items.size) :
    (({ s with tstack := a :: { t with spans := setSides dir (Items.ch s.items x) [] } :: rest }).dropCh x).Full
      g P (fun i => X i ∨ i = x) := by
  have key : ∀ i, i ≠ x →
      (({ s with tstack := a :: { t with spans := setSides dir (Items.ch s.items x) [] } :: rest }).dropCh x).cnt i
        = s.cnt i := by
    intro i hix
    simp only [cnt, WalkState.dropCh, hts, spansCount_cons, count_setSides, List.append_nil]
    have h2 := Items.chCount_modify s.items hlt (fun it => { it with ch := [] }) i
    simp only [List.count_nil, Nat.add_zero] at h2
    have h3 := count_getSide_of_onSide t.spans dir hside i
    rw [hx] at h3
    simp only [List.count_singleton, beq_iff_eq, Ne.symm hix, ite_false] at h3
    omega
  exact {
    place := h.place.unwrap hts (hx ▸ List.mem_singleton_self x) hn hlt
    fixedP := h.fixedP
    pushed := fun i hi => by
      rw [key i (Nat.ne_of_lt (Nat.lt_of_lt_of_le (h.fixedP i hi).2 hn))]; exact h.pushed i hi
    placed := fun i hn' hi hx' => by
      rw [not_or] at hx'
      rw [key i hx'.2]; exact h.placed i hn' (by simpa [WalkState.dropCh] using hi) hx'.1
    qch := fun e he hP => by
      simp only [WalkState.dropCh]
      rw [Items.ch_modify_ne _ _ _ _ (node_ne_edgeItem hn he)]; exact h.qch e he hP }

/-- Placing node `x` and fixed item `q` at once. -/
theorem step_node_fixed (h : s.Full g P X) {x q : Nat} (hq0 : 0 < q) (hqlt : q < 1 + g.nv + g.ne)
    (hp : s'.Place g (fun i => P i ∨ i = q) (fun i => X i ∧ i ≠ x)) (hsz : s'.items.size = s.items.size)
    (hqc : ∀ e, e < g.ne → ¬ P (edgeItem g e) → edgeItem g e ≠ q →
      Items.ch s'.items (edgeItem g e) = Items.ch s.items (edgeItem g e))
    (hc : ∀ i, s'.cnt i = s.cnt i + (if i = x then 1 else 0) + if i = q then 1 else 0) :
    s'.Full g (fun i => P i ∨ i = q) (fun i => X i ∧ i ≠ x) where
  place := hp
  fixedP i hi := hi.elim (h.fixedP i) fun h' => h' ▸ ⟨hq0, hqlt⟩
  pushed i hi := by
    rw [hc]
    rcases hi with hi | rfl
    · have := h.pushed i hi; omega
    · simp
  placed i hn hi hx' := by
    rw [hc]
    by_cases hix : i = x
    · simp [hix]
    · have := h.placed i hn (hsz ▸ hi) fun h' => hx' ⟨h', hix⟩; omega
  qch e he hP := by
    rw [not_or] at hP
    exact (hqc e he hP.1 hP.2).trans (h.qch e he hP.1)

/-- The bridge case of the block branch: `q` gets `x :: top.spans.2`, then goes under `vertItem curV`. -/
theorem q_cons (h : s.Full g P X) {x q j : Nat} (hx : X x) (hq0 : 0 < q) (hqlt : q < 1 + g.nv + g.ne)
    (hP : ¬ P q) (hqch : Items.ch s.items q = []) (hj : j < s.items.size)
    (hjne : ∀ e, e < g.ne → j ≠ edgeItem g e) (hb : ∀ t ∈ s.tstack.head?, t.spans.1 = []) :
    ({ s with
        items := (s.items.modify q fun it => { it with ch := x :: s.tstack.head!.spans.2 }).modify j
          fun it => { it with ch := it.ch ++ [q] },
        tstack := s.tstack.tail }).Full g (fun i => P i ∨ i = q) (fun i => X i ∧ i ≠ x) := by
  have hq : q < s.items.size := Nat.lt_of_lt_of_le hqlt h.place.size
  have hn := h.place.loose_node x hx
  have hxq : x ≠ q := by omega
  refine h.step_node_fixed hq0 hqlt
    ((h.place.q_cons hx hq).append_fixed (by simpa using hj) hq0 hqlt hP) (by simp) (fun e he _ hne => ?_) fun i => ?_
  · rw [Items.ch_modify_ne _ _ _ _ (hjne e he), Items.ch_modify_ne _ _ _ _ hne.symm]
  · have hb' := head!_spans1_nil hb
    simp only [cnt]
    have h1 := Items.chCount_modify s.items hq (fun it => { it with ch := x :: s.tstack.head!.spans.2 }) i
    have h2 := Items.chCount_modify (s.items.modify q fun it => { it with ch := x :: s.tstack.head!.spans.2 })
      (by simpa using hj) (fun it => { it with ch := it.ch ++ [q] }) i
    have h3 := spansCount_head_tail s.tstack i
    rw [← Place.ch_eq_getElem _ (by simpa using hj)] at h2
    simp only [hqch, hb', List.count_nil, List.count_cons, List.count_append, beq_iff_eq,
      List.nil_append] at h1 h2 h3
    rcases eq_or_ne i x with rfl | hix <;> rcases eq_or_ne i q with rfl | hiq
    · exact absurd rfl hxq
    · simp only [ite_true, ite_eq_right hiq, ite_eq_right (Ne.symm hiq)] at h1 h2 ⊢; omega
    · simp only [ite_true, ite_eq_right hix, ite_eq_right (Ne.symm hix)] at h1 h2 ⊢; omega
    · simp only [ite_eq_right hix, ite_eq_right (Ne.symm hix), ite_eq_right hiq, ite_eq_right (Ne.symm hiq)] at h1 h2 ⊢; omega


/-- The block branch without a bridge: `q` gets `backedge.spans.1 ++ t.spans.2`, then goes under `j`. -/
theorem q_merge (h : s.Full g P X) {q j : Nat} (hq0 : 0 < q) (hqlt : q < 1 + g.nv + g.ne) (hP : ¬ P q)
    (hqch : Items.ch s.items q = []) (hj : j < s.items.size) (hjne : ∀ e, e < g.ne → j ≠ edgeItem g e)
    (hb : ∀ b ∈ s.tstack.head?, b.spans.2 = []) (ht : ∀ t ∈ s.tstack.tail.head?, t.spans.1 = []) :
    ({ s with
        items := (s.items.modify q fun it =>
            { it with ch := s.tstack.head!.spans.1 ++ s.tstack.tail.head!.spans.2 }).modify j
          fun it => { it with ch := it.ch ++ [q] },
        tstack := s.tstack.tail.tail }).Full g (fun i => P i ∨ i = q) X := by
  have hq : q < s.items.size := Nat.lt_of_lt_of_le hqlt h.place.size
  refine h.step_fixed hq0 hqlt ((h.place.q_merge hq).append_fixed (by simpa using hj) hq0 hqlt hP) (by simp)
    (fun e he hP' => ?_) fun i => ?_
  · rw [not_or] at hP'
    rw [Items.ch_modify_ne _ _ _ _ (hjne e he), Items.ch_modify_ne _ _ _ _ (Ne.symm hP'.2)]
  · have hb' := head!_spans2_nil hb
    have ht' := head!_spans1_nil ht
    simp only [cnt]
    have h1 := Items.chCount_modify s.items hq
      (fun it => { it with ch := s.tstack.head!.spans.1 ++ s.tstack.tail.head!.spans.2 }) i
    have h2 := Items.chCount_modify (s.items.modify q fun it =>
      { it with ch := s.tstack.head!.spans.1 ++ s.tstack.tail.head!.spans.2 })
      (by simpa using hj) (fun it => { it with ch := it.ch ++ [q] }) i
    have h3 := spansCount_head_tail s.tstack i
    have h4 := spansCount_head_tail s.tstack.tail i
    rw [← Place.ch_eq_getElem _ (by simpa using hj)] at h2
    simp only [hb', ht', hqch, List.count_nil, List.count_append, List.count_singleton, beq_iff_eq,
      List.append_nil, List.nil_append] at h1 h2 h3 h4
    rcases eq_or_ne i q with rfl | hiq
    · simp only [ite_true] at h2 ⊢; omega
    · simp only [ite_eq_right hiq, ite_eq_right (Ne.symm hiq)] at h2 ⊢; omega


/-- The self-loop case: `q` gets `[x]`, then goes under `j`. -/
theorem q_single (h : s.Full g P X) {x q j : Nat} (hx : X x) (hq0 : 0 < q) (hqlt : q < 1 + g.nv + g.ne)
    (hP : ¬ P q) (hqch : Items.ch s.items q = []) (hj : j < s.items.size)
    (hjne : ∀ e, e < g.ne → j ≠ edgeItem g e) :
    ({ s with
        items := (s.items.modify q fun it => { it with ch := [x] }).modify j
          fun it => { it with ch := it.ch ++ [q] } }).Full g (fun i => P i ∨ i = q) (fun i => X i ∧ i ≠ x) := by
  have hq : q < s.items.size := Nat.lt_of_lt_of_le hqlt h.place.size
  have hn := h.place.loose_node x hx
  have hxq : x ≠ q := by omega
  refine h.step_node_fixed hq0 hqlt
    ((h.place.q_single hx hq).append_fixed (by simpa using hj) hq0 hqlt hP) (by simp) (fun e he _ hne => ?_) fun i => ?_
  · rw [Items.ch_modify_ne _ _ _ _ (hjne e he), Items.ch_modify_ne _ _ _ _ hne.symm]
  · simp only [cnt]
    have h1 := Items.chCount_modify s.items hq (fun it => { it with ch := [x] }) i
    have h2 := Items.chCount_modify (s.items.modify q fun it => { it with ch := [x] })
      (by simpa using hj) (fun it => { it with ch := it.ch ++ [q] }) i
    rw [← Place.ch_eq_getElem _ (by simpa using hj)] at h2
    simp only [hqch, List.count_nil, List.count_append, List.count_singleton, beq_iff_eq] at h1 h2
    rcases eq_or_ne i x with rfl | hix <;> rcases eq_or_ne i q with rfl | hiq
    · exact absurd rfl hxq
    · simp only [ite_true, ite_eq_right hiq, ite_eq_right (Ne.symm hiq)] at h1 h2 ⊢; omega
    · simp only [ite_true, ite_eq_right hix, ite_eq_right (Ne.symm hix)] at h1 h2 ⊢; omega
    · simp only [ite_eq_right hix, ite_eq_right (Ne.symm hix), ite_eq_right hiq, ite_eq_right (Ne.symm hiq)]
        at h1 h2 ⊢; omega

/-- `walkForest`: the last entry's second side goes under the root. -/
theorem root_append (h : s.Full g P X) (hb : ∀ t ∈ s.tstack.head?, t.spans.1 = []) :
    ({ s with
        items := s.items.modify rootItem fun it => { it with ch := it.ch ++ s.tstack.head!.spans.2 },
        tstack := s.tstack.tail }).Full g P X := by
  have h0 : rootItem < s.items.size := Nat.lt_of_lt_of_le (show (0 : Nat) < 1 + g.nv + g.ne by omega) h.place.size
  refine h.of_cnt h.place.root_append (by simp) (fun e _ _ => ?_) fun i => ?_
  · exact Items.ch_modify_ne _ _ _ _ (by show (0 : Nat) ≠ 1 + g.nv + e; omega)
  · have hb' := head!_spans1_nil hb
    simp only [cnt]
    have h1 := Items.chCount_modify s.items h0 (fun it => { it with ch := it.ch ++ s.tstack.head!.spans.2 }) i
    have h2 := spansCount_head_tail s.tstack i
    simp only [List.count_append, hb', List.nil_append, ← Place.ch_eq_getElem s.items h0] at h1 h2
    omega

theorem dropCh_mergeTop {x : Nat} (h : (s.dropCh x).Full g P X) (hm : MergeOK s) :
    (({ s with tstack := WalkM.mergeTop s.tstack }).dropCh x).Full g P X := h.mergeTop hm
theorem dropCh_fold {x : Nat} (h : (s.dropCh x).Full g P X) (curV : Nat) (edgeDir : Bool) :
    (({ s with tstack := match s.tstack with
      | a :: rest => { a with vStart := curV, spans := setSides (!edgeDir) (a.spans.1 ++ a.spans.2) [] } :: rest
      | [] => [] }).dropCh x).Full g P X := h.fold curV edgeDir

end Full
end WalkState

/-! ### The primitives, as `wp` facts -/

open WalkState

variable {g : Graph} {P X : ItemId → Prop} {s : WalkState}

theorem wp_and {A B : α → WalkState → Prop} {m : WalkM α} (ha : wp m A s) (hb : wp m B s) :
    wp m (fun a s' => A a s' ∧ B a s') s := ⟨ha, hb⟩

theorem modifyItem_spec' (i : Nat) (f : Item → Item) :
    wp (modifyItem i f) (fun _ s' => s' = { s with items := s.items.modify i f }) s := rfl
theorem modify_spec' (f : WalkState → WalkState) : wp (modify f : WalkM Unit) (fun _ s' => s' = f s) s := rfl
theorem allocItem_spec' (ty : NodeType) :
    wp (allocItem ty) (fun item s' => item = s.items.size ∧
      s' = { s with items := s.items.push ⟨ty, (none, none), []⟩ }) s := ⟨rfl, rfl⟩

theorem mergeTstackTops_full (h : s.Full g P X) (hm : MergeOK s) :
    wp mergeTstackTops (fun _ s' => s'.Full g P X) s := by
  rw [wp_mergeTstackTops]; exact h.mergeTop hm

theorem loop_full {I B : WalkState → Prop} {cond : WalkM Bool} {body : WalkM Unit}
    (hc : ∀ s, wp cond (fun _ s' => s' = s) s)
    (hb : ∀ s, I s → B s → wp body (fun _ s' => I s') s) :
    ∀ (n : Nat) (s : WalkState), I s → LoopSides B cond body n s → wp (loop n cond body) (fun _ s' => I s') s
  | 0, _, h, _ => h
  | n + 1, s, h, hl => by
    simp only [loop, wp_bind, wp_ite]
    refine wp_mono _ (wp_and (hc s) hl) ?_
    rintro b s' ⟨rfl, hl'⟩
    split
    · next hb' =>
      obtain ⟨hB, hrest⟩ := hl' hb'
      exact wp_mono _ (wp_and (hb s' h hB) hrest) fun _ s₂ ⟨h₂, hl₂⟩ => loop_full hc hb n s₂ h₂ hl₂
    · exact h

theorem maybeUnwrapNxt_full (h : s.Full g P X) (ty : NodeType) (hty : ty ∈ [NodeType.S, .P, .R])
    (hu : UnwrapOK s ty) :
    wp (maybeUnwrapNxt ty)
      (fun item s' => (s'.dropCh item).Full g P (fun i => X i ∨ i = item) ∧ ¬ X item) s := by
  have halloc : wp (allocItem ty)
      (fun item s' => (s'.dropCh item).Full g P (fun i => X i ∨ i = item) ∧ ¬ X item) s :=
    ⟨(h.push ty).dropCh_fresh (x := s.items.size) (by simp) (Items.ch_push_size _ rfl),
      fun hx => Nat.lt_irrefl _ (h.place.loose_node _ hx).2⟩
  have hroot : Items.type s.items 0 ≠ ty := by
    have hr : Items.type s.items 0 = .F := h.place.root_type
    rw [hr]; rintro rfl; simp at hty
  unfold maybeUnwrapNxt
  simp only [wp_bind, wp_get, wp_ite, wp_pure, wp_nxt, wp_stackDir, wp_getItem, wp_modifyNxt]
  split
  · exact halloc
  · next hcond =>
    split
    · next heq =>
      rw [beq_iff_eq, Items.type_getElem!'] at heq
      rcases hts : s.tstack with _ | ⟨a, _ | ⟨t, rest⟩⟩
      · rw [hts] at heq
        exact absurd heq (by cases s.stackDir[([] : List TEntry).tail.head!.topDepth]! <;> exact hroot)
      · rw [hts] at heq
        exact absurd heq (by cases s.stackDir[[a].tail.head!.topDepth]! <;> exact hroot)
      · rw [hts] at heq
        simp only [List.tail_cons, head!_cons] at heq ⊢
        cases hl : getSide t.spans s.stackDir[t.topDepth]! with
        | nil => rw [hl] at heq; exact absurd heq hroot
        | cons x xs =>
          rw [hl, List.head!_cons] at heq
          obtain ⟨rfl, hside⟩ := hu ((Bool.not_eq_true _).mp hcond) a t rest hts x xs hl heq
          simp only [List.head!_cons]
          have hlt : x < s.items.size := by
            by_contra hge
            rw [Items.type_of_le _ _ (Nat.le_of_not_lt hge)] at heq
            rw [← heq] at hty; simp at hty
          have hn := h.place.node_ge (x := x) (heq ▸ hty)
          rw [Items.ch_getElem! _ _ hlt]
          refine ⟨h.unwrap hts hl hside hn hlt, fun hx => ?_⟩
          have := h.place.loose _ hx
          have h2 := count_getSide_le t.spans s.stackDir[t.topDepth]! x
          have h3 : spansCount s.tstack x ≤ s.cnt x := Nat.le_add_right _ _
          rw [hts] at h3
          simp only [spansCount_cons, hl, List.count_singleton, beq_self_eq_true, ite_true] at h2 h3
          omega
    · exact halloc

theorem finishTstackTop_full {x : Nat} (h : (s.dropCh x).Full g P X) (hx : X x) (hc : CloseOK s) :
    wp (finishTstackTop x) (fun _ s' => s'.Full g P (fun i => X i ∧ i ≠ x)) s := by
  obtain ⟨a, rest, hts, hside⟩ := hc
  unfold finishTstackTop
  simp only [wp_bind, wp_cur, wp_stackDir, wp_makeVs, wp_modifyItem, wp_modifyCur, hts, head!_cons]
  exact h.finishTop hx hts hside _

end Spqr
