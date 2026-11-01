import Spqr.RelabelPieceSep
import Spqr.RangesCanon
import Spqr.RangesOwned

/-! # The PieceSep item facts as a walk-state invariant (`PieceFacts`, `PieceInv`, `ClosePiece`)

`spqrTree_pieceSep` needs five facts about the final items (`Items.QUpper`/`RootSep`/`RootV`/
`QChildVs`/`PChildVs`, `RelabelPieceSep.lean`). `PieceFacts s` states them for the items of a walk
state; `PieceInv d s` adds the frame clause `root_path` (no vertex item of the open DFS path lies
below a root child), which is what keeps `RootSep` through the boundary writes. They are kept by
every walk primitive given the per-site record `ClosePiece` (the orientation of the closed node at
a boundary and `FinishPiece` at the three `finishTstackTop` sites: a P's new children carry its
`vs`); `finishEdge_piece` is the step lemma the walk induction (`WalkBackbone.lean`) uses, and the
executable checker (`checks/WalkInvCheck/Ranges.lean`, `checkPieceInv`/`checkPiece`) evaluates
every clause (`piece_*`). -/

namespace Spqr
namespace WalkState
open WalkM

/-- The five final-items facts of `spqrTree_pieceSep` at a walk state. -/
structure PieceFacts (s : WalkState) : Prop where
  q_upper : Items.QUpper s.g s.items
  root_sep : Items.RootSep s.g s.items
  root_v : Items.RootV s.items
  q_child_vs : Items.QChildVs s.g s.items
  p_child_vs : Items.PChildVs s.items

/-- `PieceFacts` plus the frame clause: no vertex item of the open path `stackVerts[0..d]` lies below
a root child (root children are the finished DFS trees). -/
structure PieceInv (d : Nat) (s : WalkState) : Prop extends PieceFacts s where
  root_path : ∀ a, Items.IsParent s.items rootItem a → ∀ k, k ≤ d →
    ¬ Items.Below s.items a (vertItem s.stackVerts[k]!)

/-- `finishTstackTop x` from `r` keeps `PChildVs`: if `x` is a P, every non-V item on the closing
side of the top entry `t` already carries `makeVs t.vStart t.topDepth`, the `vs` the close writes
to `x`. -/
def FinishPiece (x : ItemId) (r : WalkState) : Prop :=
  Items.type r.items x = .P →
    ∀ c ∈ getSide (curE r).spans r.stackDir[(curE r).topDepth]!, Items.type r.items c ≠ .V →
      Items.vs r.items c =
        setSides r.stackDir[(curE r).topDepth]! (some r.stackVerts[(curE r).topDepth]!) (some (curE r).vStart)

/-- The PieceSep site facts of one `finishEdge curV d o origTstack hasVert` call from `s`: at a
boundary the new I node (`bd_bridge`) or the closed block node (`bd_node`) is oriented
`(curV, dest)`; `FinishPiece` at the P close (`finishP`, whenever `condP` holds at `feRest`), the
type-1 vertex close (`cvS₅`) and every loop-1 iteration whose first `k` conditions held. -/
structure ClosePiece (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool)
    (s : WalkState) : Prop where
  bd_bridge : o.cls.isTree = true → d ≤ o.cls.lowval d → o.cls.lowval d = d + 1 →
    setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) = (some curV, some o.dest)
  bd_node : o.cls.isTree = true → d ≤ o.cls.lowval d → o.cls.lowval d ≠ d + 1 →
    ∀ b ∈ s.tstack.head?, ∀ c, b.spans.1 = [c] → Items.vs s.items c = (some curV, some o.dest)
  p_site : o.cls.lowval d < d → o.cls.isType1 = true →
    result (condP curV (o.cls.lowval d) true) (feRest curV d o origTstack hasVert s) = true →
    FinishPiece (pItem (feRest curV d o origTstack hasVert s)) (pPre (feRest curV d o origTstack hasVert s))
  v_site : o.cls.isTree = true → o.cls.lowval d < d → hasVert = true → o.cls.isType1 = true →
    FinishPiece ((maybeUnwrapNxt (if feSingle d o s then NodeType.S else .R)).run (feS₂ d o s)).1
      (cvS₅ curV s.stackDir[d]! true origTstack (feSingle d o s) (feS₂ d o s))
  l1_site : o.cls.isTree = true → o.cls.lowval d < d → ∀ k,
    (∀ j, j ≤ k → result (loop1Cond d) (l1Iter d o s j) = true) →
    FinishPiece (l1Node d o s k) (l1Pre d o s k)

end WalkState

/-! ## Preservation through the item primitives -/

namespace Items
variable {items : Items}

theorem IsParent_modify_iff (j : ItemId) (f : Item → Item) (p c : ItemId) :
    Items.IsParent (items.modify j f) p c ↔
      if p = j then ∃ hj : j < items.size, c ∈ (f items[j]).ch else items.IsParent p c := by
  by_cases hp : p = j
  · subst hp
    simp only [IsParent, ite_true]
    by_cases hj : p < items.size
    · rw [ch_modify_at p f hj]; exact ⟨fun h => ⟨hj, h⟩, fun ⟨_, h⟩ => h⟩
    · have hn : items[p]? = none := Array.getElem?_eq_none_iff.mpr (Nat.le_of_not_lt hj)
      have : Items.ch (items.modify p f) p = [] := by simp [ch, Array.getElem?_modify, hn]
      rw [this]; simp [hj]
  · rw [if_neg hp]; simp only [IsParent, ch_modify_of_ne j f hp]

/-- Modifying an item with no parent changes no `Below` relation from another item. -/
theorem Below_modify_of_no_parent {x : ItemId} (f : Item → Item) (hx : ∀ p, ¬ items.IsParent p x)
    {a i : ItemId} (ha : a ≠ x) : Items.Below (items.modify x f) a i ↔ items.Below a i := by
  have key : ∀ i, Items.Below (items.modify x f) a i → items.Below a i ∧ i ≠ x := by
    intro i h
    induction h with
    | refl => exact ⟨.refl, ha⟩
    | tail _ hpc ih =>
      have hc := (IsParent_modify_iff x f _ _).1 hpc
      rw [if_neg ih.2] at hc
      exact ⟨ih.1.tail hc, fun h => hx _ (h ▸ hc)⟩
  have key' : ∀ i, items.Below a i → Items.Below (items.modify x f) a i ∧ i ≠ x := by
    intro i h
    induction h with
    | refl => exact ⟨.refl, ha⟩
    | @tail b c _ hpc ih =>
      have hc : Items.IsParent (items.modify x f) b c := by
        rw [IsParent_modify_iff, if_neg ih.2]; exact hpc
      exact ⟨ih.1.tail hc, fun h => hx _ (h ▸ hpc)⟩
  exact ⟨fun h => (key i h).1, fun h => (key' i h).1⟩

theorem ch_eq_nil_of_le {p : ItemId} (hp : items.size ≤ p) : items.ch p = [] := by
  simp [ch, Array.getElem?_eq_none_iff.mpr hp]

end Items

namespace WalkState
open WalkM
variable {s s' : WalkState} {d : Nat}

theorem PieceFacts.congr (h : PieceFacts s) (hg : s'.g = s.g) (hsz : s'.items.size = s.items.size)
    (hty : ∀ p, Items.type s'.items p = Items.type s.items p)
    (hch : ∀ p, Items.ch s'.items p = Items.ch s.items p)
    (hvs : ∀ p, Items.vs s'.items p = Items.vs s.items p) : PieceFacts s' where
  q_upper := by
    intro v c hv hpc
    rw [hg] at hv; rw [Items.IsParent_congr hch] at hpc; rw [hvs]
    exact h.q_upper v c hv hpc
  root_sep := by
    intro a b ha hb hab v e e' he he' hi hi' hea heb
    rw [hg] at he he' hi hi'
    rw [Items.IsParent_congr hch] at ha hb
    unfold Items.EdgeBelow at hea heb
    rw [hg, Items.Below_congr hch] at hea heb
    exact h.root_sep a b ha hb hab v e e' he he' hi hi' hea heb
  root_v := by
    intro c hc; rw [Items.IsParent_congr hch] at hc; rw [hty]; exact h.root_v c hc
  q_child_vs := by
    intro e he; rw [hg] at he ⊢; simp only [hch, hvs]; exact h.q_child_vs e he
  p_child_vs := by
    intro i c hi hti hpc htc
    rw [hsz] at hi; rw [hty] at hti htc; rw [Items.IsParent_congr hch] at hpc; rw [hvs, hvs]
    exact h.p_child_vs i c hi hti hpc htc

theorem PieceFacts.frame (h : PieceFacts s) (hg : s'.g = s.g) (hi : s'.items = s.items) : PieceFacts s' :=
  h.congr hg (by rw [hi]) (fun _ => by rw [hi]) (fun _ => by rw [hi]) (fun _ => by rw [hi])

theorem PieceInv.congr (h : PieceInv d s) (hg : s'.g = s.g) (hsz : s'.items.size = s.items.size)
    (hty : ∀ p, Items.type s'.items p = Items.type s.items p)
    (hch : ∀ p, Items.ch s'.items p = Items.ch s.items p)
    (hvs : ∀ p, Items.vs s'.items p = Items.vs s.items p)
    (hsv : ∀ k, k ≤ d → s'.stackVerts[k]! = s.stackVerts[k]!) : PieceInv d s' where
  toPieceFacts := h.toPieceFacts.congr hg hsz hty hch hvs
  root_path := by
    intro a ha k hk
    rw [Items.IsParent_congr hch] at ha; rw [hsv k hk, Items.Below_congr hch]
    exact h.root_path a ha k hk

theorem PieceInv.frame (h : PieceInv d s) (hg : s'.g = s.g) (hi : s'.items = s.items)
    (hsv : s'.stackVerts = s.stackVerts) : PieceInv d s' :=
  h.congr hg (by rw [hi]) (fun _ => by rw [hi]) (fun _ => by rw [hi]) (fun _ => by rw [hi])
    (fun _ _ => by rw [hsv])

theorem PieceInv.exit (h : PieceInv (d + 1) s) : PieceInv d s :=
  ⟨h.toPieceFacts, fun a ha k hk => h.root_path a ha k (by omega)⟩

/-- Entering the child `v` at depth `d + 1`. -/
theorem PieceInv.setVert (h : PieceInv d s) (v : Nat) (hd : d + 1 < s.stackVerts.size)
    (hnb : ∀ a, Items.IsParent s.items rootItem a → ¬ Items.Below s.items a (vertItem v)) :
    PieceInv (d + 1) { s with stackVerts := s.stackVerts.set! (d + 1) v } where
  toPieceFacts := h.toPieceFacts.frame rfl rfl
  root_path := by
    intro a ha k hk
    show ¬ Items.Below s.items a (vertItem (s.stackVerts.set! (d + 1) v)[k]!)
    by_cases hkd : k = d + 1
    · subst hkd; rw [getElem!_set!_self _ _ _ hd]; exact hnb a ha
    · rw [getElem!_set!_ne _ _ _ _ hkd]; exact h.root_path a ha k (by omega)

/-- Starting a new DFS tree at the root `v` (depth 0): `PieceFacts` between trees becomes `PieceInv 0`
once `vertItem v` (unvisited: no root child lies above it) is placed at `stackVerts[0]`. -/
theorem PieceFacts.root (h : PieceFacts s) (v : Nat) (hd : 0 < s.stackVerts.size)
    (hnb : ∀ a, Items.IsParent s.items rootItem a → ¬ Items.Below s.items a (vertItem v)) :
    PieceInv 0 { s with stackVerts := s.stackVerts.set! 0 v } where
  toPieceFacts := h.frame rfl rfl
  root_path := by
    intro a ha k hk
    have hk0 : k = 0 := by omega
    subst hk0
    show ¬ Items.Below s.items a (vertItem (s.stackVerts.set! 0 v)[0]!)
    rw [getElem!_set!_self _ _ _ hd]; exact hnb a ha

theorem PieceInv.init (g : Graph) (tern : Bool) : PieceInv 0 (init g tern) := by
  have hch : ∀ p, Items.ch (WalkState.init g tern).items p = [] := fun p => Items.initialItems_ch g p
  have hnp : ∀ p c, ¬ Items.IsParent (WalkState.init g tern).items p c := fun p c h => by
    simp [Items.IsParent, hch] at h
  refine ⟨⟨?_, ?_, ?_, ?_, ?_⟩, ?_⟩
  · intro v c _ h; exact absurd h (hnp _ _)
  · intro a _ ha; exact absurd ha (hnp _ _)
  · intro c hc; exact absurd hc (hnp _ _)
  · intro e _; exact ⟨fun c w h => by simp [hch] at h, fun c h => by simp [hch] at h⟩
  · intro i c _ _ h; exact absurd h (hnp _ _)
  · intro a ha; exact absurd ha (hnp _ _)

/-- Pushing a fresh childless item. -/
theorem PieceInv.alloc (h : PieceInv d s) (hs : Shape s) (ty : NodeType) :
    PieceInv d { s with items := s.items.push ⟨ty, (none, none), []⟩ } := by
  have hch : ∀ p, Items.ch (s.items.push ⟨ty, (none, none), []⟩) p = Items.ch s.items p :=
    Items.ch_push_nil _ rfl
  have hpar : ∀ p c, Items.IsParent (s.items.push ⟨ty, (none, none), []⟩) p c ↔ Items.IsParent s.items p c :=
    fun _ _ => Items.IsParent_congr hch
  have hclt : ∀ p c, Items.IsParent s.items p c → c ≠ s.items.size := fun p c hpc =>
    Nat.ne_of_lt (hs.ch_lt p c hpc)
  have hplt : ∀ p c, Items.IsParent s.items p c → p ≠ s.items.size := fun p c hpc =>
    Nat.ne_of_lt (Items.parent_lt hpc)
  refine ⟨⟨?_, ?_, ?_, ?_, ?_⟩, ?_⟩
  · intro v c hv hpc
    rw [hpar] at hpc
    show (Items.vs (s.items.push _) c).1 = _
    rw [Items.vs_push_of_ne _ (hclt _ _ hpc)]; exact h.q_upper v c hv hpc
  · intro a b ha hb hab v e e' he he' hi hi' hea heb
    rw [hpar] at ha hb
    unfold Items.EdgeBelow at hea heb
    rw [Items.Below_push_nil _ rfl] at hea heb
    exact h.root_sep a b ha hb hab v e e' he he' hi hi' hea heb
  · intro c hc
    rw [hpar] at hc
    show Items.type (s.items.push _) c = _
    rw [Items.type_push_of_ne _ (hclt _ _ hc)]; exact h.root_v c hc
  · intro e he
    change e < s.g.ne at he
    have hesz : edgeItem s.g e ≠ s.items.size := by
      have := hs.size; show 1 + s.g.nv + e ≠ _; omega
    show (∀ c w, Items.ch (s.items.push _) (edgeItem s.g e) = [c, vertItem w] →
        Items.vs (s.items.push _) c = ((Items.vs (s.items.push _) (edgeItem s.g e)).1, some w)) ∧
      ∀ c, Items.ch (s.items.push _) (edgeItem s.g e) = [c] →
        Items.vs (s.items.push _) c = ((Items.vs (s.items.push _) (edgeItem s.g e)).1, none)
    rw [hch, Items.vs_push_of_ne _ hesz]
    constructor
    · intro c w hcw
      have hpc : Items.IsParent s.items (edgeItem s.g e) c := by simp [Items.IsParent, hcw]
      rw [Items.vs_push_of_ne _ (hclt _ _ hpc)]; exact (h.q_child_vs e he).1 c w hcw
    · intro c hcw
      have hpc : Items.IsParent s.items (edgeItem s.g e) c := by simp [Items.IsParent, hcw]
      rw [Items.vs_push_of_ne _ (hclt _ _ hpc)]; exact (h.q_child_vs e he).2 c hcw
  · intro i c hi hti hpc htc
    rw [hpar] at hpc
    have hi' := Items.parent_lt hpc
    show Items.vs (s.items.push _) c = Items.vs (s.items.push _) i
    rw [Items.vs_push_of_ne _ (hclt _ _ hpc), Items.vs_push_of_ne _ (hplt _ _ hpc)]
    rw [show Items.type (s.items.push ⟨ty, (none, none), []⟩) i = Items.type s.items i from
      Items.type_push_of_ne _ (hplt _ _ hpc)] at hti
    rw [show Items.type (s.items.push ⟨ty, (none, none), []⟩) c = Items.type s.items c from
      Items.type_push_of_ne _ (hclt _ _ hpc)] at htc
    exact h.p_child_vs i c hi' hti hpc htc
  · intro a ha k hk
    rw [hpar] at ha
    show ¬ Items.Below (s.items.push _) a (vertItem s.stackVerts[k]!)
    rw [Items.Below_push_nil _ rfl]; exact h.root_path a ha k hk

/-- Rewriting an item `x` with no parent that is neither the root nor a vertex item, given the
`QChildVs`/`PChildVs` clauses at `x` for the new record. -/
theorem PieceInv.modify (h : PieceInv d s) (hs : Shape s) (x : ItemId) (f : Item → Item)
    (hf : ∀ it, (f it).type = it.type)
    (hx : ∀ p, ¬ Items.IsParent s.items p x)
    (hF : Items.type s.items x ≠ .F) (hV : Items.type s.items x ≠ .V)
    (hQ : ∀ e, e < s.g.ne → x = edgeItem s.g e →
      (∀ c w, Items.ch (s.items.modify x f) x = [c, vertItem w] →
        Items.vs (s.items.modify x f) c = ((Items.vs (s.items.modify x f) x).1, some w)) ∧
      ∀ c, Items.ch (s.items.modify x f) x = [c] →
        Items.vs (s.items.modify x f) c = ((Items.vs (s.items.modify x f) x).1, none))
    (hP : Items.type s.items x = .P → ∀ c ∈ Items.ch (s.items.modify x f) x,
      Items.type s.items c ≠ .V → Items.vs (s.items.modify x f) c = Items.vs (s.items.modify x f) x) :
    PieceInv d { s with items := s.items.modify x f } := by
  have hty : ∀ p, Items.type (s.items.modify x f) p = Items.type s.items p :=
    Items.type_modify_type_eq x f hf
  have hpar : ∀ p c, Items.IsParent (s.items.modify x f) p c → p ≠ x → Items.IsParent s.items p c := by
    intro p c hpc hpx
    have := (Items.IsParent_modify_iff (items := s.items) x f p c).1 hpc
    rwa [if_neg hpx] at this
  have hne : ∀ p c, Items.IsParent s.items p c → c ≠ x := fun p c hpc hcx => hx p (hcx ▸ hpc)
  have hvs : ∀ c, c ≠ x → Items.vs (s.items.modify x f) c = Items.vs s.items c := fun c hc =>
    Items.vs_modify_of_ne x f hc
  have hch : ∀ p, p ≠ x → Items.ch (s.items.modify x f) p = Items.ch s.items p := fun p hp =>
    Items.ch_modify_of_ne x f hp
  have hroot : rootItem ≠ x := fun h => hF (h ▸ hs.root)
  have hvert : ∀ v, v < s.g.nv → vertItem v ≠ x := fun v hv h => hV (h ▸ hs.vert v hv)
  refine ⟨⟨?_, ?_, ?_, ?_, ?_⟩, ?_⟩
  · intro v c hv hpc
    have hpc' := hpar _ _ hpc (hvert v hv)
    show (Items.vs (s.items.modify x f) c).1 = _
    rw [hvs c (hne _ _ hpc')]; exact h.q_upper v c hv hpc'
  · intro a b ha hb hab v e e' he he' hi hi' hea heb
    have ha' := hpar _ _ ha hroot
    have hb' := hpar _ _ hb hroot
    unfold Items.EdgeBelow at hea heb
    rw [Items.Below_modify_of_no_parent f hx (hne _ _ ha')] at hea
    rw [Items.Below_modify_of_no_parent f hx (hne _ _ hb')] at heb
    exact h.root_sep a b ha' hb' hab v e e' he he' hi hi' hea heb
  · intro c hc
    show Items.type (s.items.modify x f) c = _
    rw [hty]; exact h.root_v c (hpar _ _ hc hroot)
  · intro e he
    change e < s.g.ne at he
    by_cases hxe : x = edgeItem s.g e
    · subst hxe; exact hQ e he rfl
    · show (∀ c w, Items.ch (s.items.modify x f) (edgeItem s.g e) = [c, vertItem w] →
          Items.vs (s.items.modify x f) c = ((Items.vs (s.items.modify x f) (edgeItem s.g e)).1, some w)) ∧
        ∀ c, Items.ch (s.items.modify x f) (edgeItem s.g e) = [c] →
          Items.vs (s.items.modify x f) c = ((Items.vs (s.items.modify x f) (edgeItem s.g e)).1, none)
      rw [hch _ (Ne.symm hxe), hvs _ (Ne.symm hxe)]
      constructor
      · intro c w hcw
        have hpc : Items.IsParent s.items (edgeItem s.g e) c := by simp [Items.IsParent, hcw]
        rw [hvs c (hne _ _ hpc)]; exact (h.q_child_vs e he).1 c w hcw
      · intro c hcw
        have hpc : Items.IsParent s.items (edgeItem s.g e) c := by simp [Items.IsParent, hcw]
        rw [hvs c (hne _ _ hpc)]; exact (h.q_child_vs e he).2 c hcw
  · intro i c hi hti hpc htc
    have hi' : i < s.items.size := by simpa [Array.size_modify] using hi
    rw [show Items.type (s.items.modify x f) i = Items.type s.items i from hty i] at hti
    rw [show Items.type (s.items.modify x f) c = Items.type s.items c from hty c] at htc
    by_cases hix : i = x
    · subst hix; exact hP hti c hpc htc
    · have hpc' := hpar _ _ hpc hix
      show Items.vs (s.items.modify x f) c = Items.vs (s.items.modify x f) i
      rw [hvs c (hne _ _ hpc'), hvs i hix]
      exact h.p_child_vs i c hi' hti hpc' htc
  · intro a ha k hk
    have ha' := hpar _ _ ha hroot
    show ¬ Items.Below (s.items.modify x f) a (vertItem s.stackVerts[k]!)
    rw [Items.Below_modify_of_no_parent f hx (hne _ _ ha')]; exact h.root_path a ha' k hk

/-- Writing the `vs` of a childless item with no parent (not the root, not a vertex item). -/
theorem PieceInv.modifyVs (h : PieceInv d s) (hs : Shape s) (x : ItemId)
    (vsv : Option Nat × Option Nat) (hx : ∀ p, ¬ Items.IsParent s.items p x)
    (hF : Items.type s.items x ≠ .F) (hV : Items.type s.items x ≠ .V) (hch : Items.ch s.items x = []) :
    PieceInv d { s with items := s.items.modify x fun it => { it with vs := vsv } } := by
  have hch' : Items.ch (s.items.modify x fun it => { it with vs := vsv }) x = [] := by
    rw [Items.ch_modify_ch_eq x (fun it => { it with vs := vsv }) (fun _ => rfl)]; exact hch
  refine h.modify hs x _ (fun _ => rfl) hx hF hV ?_ ?_
  · intro e _ _; rw [hch']; exact ⟨fun _ _ h => by simp at h, fun _ h => by simp at h⟩
  · intro _ c hc; rw [hch'] at hc; simp at hc

/-- Appending the block root `q` (oriented `(curV, _)`) under `vertItem curV`, where no root child
reaches `vertItem curV`. -/
theorem PieceInv.vAppend (h : PieceInv d s) (hs : Shape s) {curV : Nat} (hv : curV < s.g.nv) (q : ItemId)
    (hq : ∀ p, ¬ Items.IsParent s.items p q) (hqt : Items.type s.items q = .Q)
    (hqvs : (Items.vs s.items q).1 = some curV)
    (hnb : ∀ a, Items.IsParent s.items rootItem a → ¬ Items.Below s.items a (vertItem curV)) :
    PieceInv d { s with items := s.items.modify (vertItem curV) fun it => { it with ch := it.ch ++ [q] } } := by
  set M := s.items.modify (vertItem curV) fun it => { it with ch := it.ch ++ [q] } with hM
  have hvsz : vertItem curV < s.items.size := by have := hs.size; show 1 + curV < _; omega
  have hty : ∀ p, Items.type M p = Items.type s.items p :=
    Items.type_modify_type_eq _ _ fun _ => rfl
  have hvs : ∀ p, Items.vs M p = Items.vs s.items p := fun p =>
    Items.vs_modify_of_vs _ _ p _ fun _ => rfl
  have hch : ∀ p, p ≠ vertItem curV → Items.ch M p = Items.ch s.items p := fun p hp =>
    Items.ch_modify_of_ne _ _ hp
  have hpar : ∀ p c, Items.IsParent M p c → p ≠ vertItem curV → Items.IsParent s.items p c := by
    intro p c hpc hpx
    have := (Items.IsParent_modify_iff (items := s.items) _ _ p c).1 hpc
    rwa [if_neg hpx] at this
  have hparV : ∀ c, Items.IsParent M (vertItem curV) c → Items.IsParent s.items (vertItem curV) c ∨ c = q := by
    intro c hpc
    have := (Items.IsParent_modify_iff (items := s.items) _ _ _ c).1 hpc
    rw [if_pos rfl] at this
    obtain ⟨_, hc⟩ := this
    simp only [List.mem_append, List.mem_singleton] at hc
    rcases hc with hc | hc
    · left; show c ∈ Items.ch s.items (vertItem curV); rw [Items.ch_eq_getElem hvsz]; exact hc
    · right; exact hc
  have hbel : ∀ a i, Items.Below M a i ↔ Items.Below s.items a i ∨
      (Items.Below s.items a (vertItem curV) ∧ ∃ q' ∈ [q], Items.Below s.items q' i) := fun a i =>
    Below_modify_ch_append _ _ hvsz
  have hbel' : ∀ a, Items.IsParent s.items rootItem a → ∀ i, Items.Below M a i ↔ Items.Below s.items a i := by
    intro a ha i; rw [hbel]
    exact ⟨fun h => h.elim id fun h => absurd h.1 (hnb a ha), Or.inl⟩
  have hroot : rootItem ≠ vertItem curV := by show 0 ≠ 1 + curV; omega
  refine ⟨⟨?_, ?_, ?_, ?_, ?_⟩, ?_⟩
  · intro v c hv' hpc
    show (Items.vs M c).1 = _
    rw [hvs]
    by_cases hvv : vertItem v = vertItem curV
    · have hvc : v = curV := by simp [vertItem] at hvv; omega
      subst hvc
      rcases hparV c hpc with hpc' | rfl
      · exact h.q_upper v c hv hpc'
      · exact hqvs
    · exact h.q_upper v c hv' (hpar _ _ hpc hvv)
  · intro a b ha hb hab v e e' he he' hi hi' hea heb
    have ha' := hpar _ _ ha hroot
    have hb' := hpar _ _ hb hroot
    unfold Items.EdgeBelow at hea heb
    rw [hbel' a ha'] at hea; rw [hbel' b hb'] at heb
    exact h.root_sep a b ha' hb' hab v e e' he he' hi hi' hea heb
  · intro c hc
    show Items.type M c = _
    rw [hty]; exact h.root_v c (hpar _ _ hc hroot)
  · intro e he
    change e < s.g.ne at he
    have hev : edgeItem s.g e ≠ vertItem curV := by show 1 + s.g.nv + e ≠ 1 + curV; omega
    show (∀ c w, Items.ch M (edgeItem s.g e) = [c, vertItem w] →
        Items.vs M c = ((Items.vs M (edgeItem s.g e)).1, some w)) ∧
      ∀ c, Items.ch M (edgeItem s.g e) = [c] → Items.vs M c = ((Items.vs M (edgeItem s.g e)).1, none)
    rw [hch _ hev]; simp only [hvs]; exact h.q_child_vs e he
  · intro i c hi hti hpc htc
    have hi' : i < s.items.size := by simpa [hM, Array.size_modify] using hi
    rw [show Items.type M i = Items.type s.items i from hty i] at hti
    rw [show Items.type M c = Items.type s.items c from hty c] at htc
    have hiv : i ≠ vertItem curV := fun h => by rw [h, hs.vert curV hv] at hti; cases hti
    show Items.vs M c = Items.vs M i
    rw [hvs, hvs]; exact h.p_child_vs i c hi' hti (hpar _ _ hpc hiv) htc
  · intro a ha k hk
    have ha' := hpar _ _ ha hroot
    show ¬ Items.Below M a (vertItem s.stackVerts[k]!)
    rw [hbel' a ha']; exact h.root_path a ha' k hk

/-- Appending the finished trees' vertex items `L` under the root (`walkForest`). -/
theorem PieceFacts.rootAppend (h : PieceFacts s) (hs : Shape s) (L : List ItemId)
    (hroot : ∀ p, ¬ Items.IsParent s.items p rootItem)
    (hLv : ∀ c ∈ L, Items.type s.items c = .V)
    (hsep : ∀ a, (Items.IsParent s.items rootItem a ∨ a ∈ L) → ∀ b ∈ L, a ≠ b →
      ∀ v e e', e < s.g.ne → e' < s.g.ne → s.g.Inc e v → s.g.Inc e' v →
        Items.EdgeBelow s.g s.items a e → Items.EdgeBelow s.g s.items b e' → False) :
    PieceFacts { s with items := s.items.modify rootItem fun it => { it with ch := it.ch ++ L } } := by
  set M := s.items.modify rootItem fun it => { it with ch := it.ch ++ L } with hM
  have hrsz : rootItem < s.items.size := by have := hs.size; show (0 : Nat) < s.items.size; omega
  have hty : ∀ p, Items.type M p = Items.type s.items p :=
    Items.type_modify_type_eq _ _ fun _ => rfl
  have hvs : ∀ p, Items.vs M p = Items.vs s.items p := fun p =>
    Items.vs_modify_of_vs _ _ p _ fun _ => rfl
  have hch : ∀ p, p ≠ rootItem → Items.ch M p = Items.ch s.items p := fun p hp =>
    Items.ch_modify_of_ne _ _ hp
  have hpar : ∀ p c, Items.IsParent M p c → p ≠ rootItem → Items.IsParent s.items p c := by
    intro p c hpc hpx
    have := (Items.IsParent_modify_iff (items := s.items) _ _ p c).1 hpc
    rwa [if_neg hpx] at this
  have hparR : ∀ c, Items.IsParent M rootItem c → Items.IsParent s.items rootItem c ∨ c ∈ L := by
    intro c hpc
    have := (Items.IsParent_modify_iff (items := s.items) _ _ _ c).1 hpc
    rw [if_pos rfl] at this
    obtain ⟨_, hc⟩ := this
    simp only [List.mem_append] at hc
    rcases hc with hc | hc
    · left; show c ∈ Items.ch s.items rootItem; rw [Items.ch_eq_getElem hrsz]; exact hc
    · right; exact hc
  have hner : ∀ a, Items.IsParent s.items rootItem a ∨ a ∈ L → a ≠ rootItem := by
    intro a ha h; subst h
    rcases ha with ha | ha
    · exact hroot _ ha
    · have := hLv _ ha; rw [hs.root] at this; cases this
  have hbel : ∀ a, Items.IsParent s.items rootItem a ∨ a ∈ L → ∀ i, Items.Below M a i ↔ Items.Below s.items a i := by
    intro a ha i
    rw [Below_modify_ch_append _ _ hrsz]
    refine ⟨fun h => h.elim id fun h => absurd ((Below_of_root hroot).1 h.1) (hner a ha), Or.inl⟩
  refine ⟨?_, ?_, ?_, ?_, ?_⟩
  · intro v c hv hpc
    have hvr : vertItem v ≠ rootItem := by show 1 + v ≠ 0; omega
    show (Items.vs M c).1 = _
    rw [hvs]; exact h.q_upper v c hv (hpar _ _ hpc hvr)
  · intro a b ha hb hab v e e' he he' hi hi' hea heb
    have ha' := hparR a ha
    have hb' := hparR b hb
    unfold Items.EdgeBelow at hea heb
    rw [hbel a ha'] at hea; rw [hbel b hb'] at heb
    rcases hb' with hb' | hb'
    · rcases ha' with ha' | ha'
      · exact h.root_sep a b ha' hb' hab v e e' he he' hi hi' hea heb
      · exact hsep b (Or.inl hb') a ha' (Ne.symm hab) v e' e he' he hi' hi heb hea
    · exact hsep a ha' b hb' hab v e e' he he' hi hi' hea heb
  · intro c hc
    show Items.type M c = _
    rw [hty]
    rcases hparR c hc with hc' | hc'
    · exact h.root_v c hc'
    · exact hLv c hc'
  · intro e he
    change e < s.g.ne at he
    have her : edgeItem s.g e ≠ rootItem := by show 1 + s.g.nv + e ≠ 0; omega
    show (∀ c w, Items.ch M (edgeItem s.g e) = [c, vertItem w] →
        Items.vs M c = ((Items.vs M (edgeItem s.g e)).1, some w)) ∧
      ∀ c, Items.ch M (edgeItem s.g e) = [c] → Items.vs M c = ((Items.vs M (edgeItem s.g e)).1, none)
    rw [hch _ her]; simp only [hvs]; exact h.q_child_vs e he
  · intro i c hi hti hpc htc
    have hi' : i < s.items.size := by simpa [hM, Array.size_modify] using hi
    rw [show Items.type M i = Items.type s.items i from hty i] at hti
    rw [show Items.type M c = Items.type s.items c from hty c] at htc
    have hir : i ≠ rootItem := fun h => by rw [h, hs.root] at hti; cases hti
    show Items.vs M c = Items.vs M i
    rw [hvs, hvs]; exact h.p_child_vs i c hi' hti (hpar _ _ hpc hir) htc

/-! ## Preservation through the walk primitives -/

theorem PieceInv.mergeTop (h : PieceInv d s) : PieceInv d (after mergeTstackTops s) := by
  rw [after_mergeTstackTops]; exact h.frame rfl rfl rfl

theorem PieceInv.loop (cond : WalkM Bool) (body : WalkM Unit) (fuel : Nat)
    (hcond : ∀ s, (cond.run s).2 = s)
    (hbody : ∀ k, (∀ j, j ≤ k → (cond.run (iter body j s)).1 = true) →
      PieceInv d (iter body k s) → PieceInv d (after body (iter body k s)))
    (h : PieceInv d s) : PieceInv d (after (WalkM.loop fuel cond body) s) := by
  induction fuel generalizing s with
  | zero => exact h
  | succ fuel ih =>
    show PieceInv d ((WalkM.loop (fuel + 1) cond body).run s).2
    rw [loop_succ_run fuel cond body s (hcond s)]
    by_cases hc : (cond.run s).1 = true
    · simp only [hc, ↓reduceIte]
      have h₁ := hbody 0 (fun j hj => by rw [Nat.le_zero.1 hj]; exact hc) h
      exact ih (s := (body.run s).2) (fun k hk => hbody (k + 1) fun j hj => by
        cases j with
        | zero => exact hc
        | succ j => exact hk j (Nat.le_of_succ_le_succ hj)) h₁
    · simp only [hc]; exact h

theorem PieceInv.mergeLoop (h : PieceInv d s) (cond : WalkM Bool) (hcond : ∀ s, (cond.run s).2 = s)
    (fuel : Nat) : PieceInv d (after (WalkM.loop fuel cond mergeTstackTops) s) :=
  PieceInv.loop cond mergeTstackTops fuel hcond (fun _ _ hk => hk.mergeTop) h

theorem PieceInv.mergeLate (h : PieceInv d s) (d' : Nat) : PieceInv d (after (Spqr.mergeLate d') s) := by
  show PieceInv d ((Spqr.mergeLate d').run s).2
  rw [mergeLate_run]
  by_cases hc : (curE s).firstIdx > s.firstOccurrence[d']!
  · simp only [hc, ↓reduceIte]
    exact h.mergeLoop _ (fun _ => rfl) _
  · simp only [hc, ↓reduceIte]; exact h

theorem PieceInv.vertPre (h : PieceInv d s) (isType1 : Bool) (origTstack : Nat) (isSingle : Bool) :
    PieceInv d (cvS₁ isType1 origTstack isSingle s) := by
  cases isType1
  · show PieceInv d ((WalkState.vertPre false origTstack isSingle).run s).2
    simp only [WalkState.vertPre, Bool.not_false, ↓reduceIte, WalkM.run_bind, run_tstackSize,
      WalkM.pure_run]
    exact h.mergeLoop _ (fun _ => rfl) _
  · exact h

theorem PieceInv.retarget (h : PieceInv d s) (curV : Nat) (dir : Bool) :
    PieceInv d (after (WalkState.retarget curV dir) s) := by
  show PieceInv d ((WalkState.retarget curV dir).run s).2
  rw [WalkState.retarget, run_modifyCur]; exact h.frame rfl rfl rfl

theorem PieceInv.pushVert (h : PieceInv d s) (v d' : Nat) : PieceInv d (after (pushVertTstack v d') s) := by
  show PieceInv d ((pushVertTstack v d').run s).2
  rw [pushVertTstack, run_pushTstack]; exact h.frame rfl rfl rfl

theorem PieceInv.finishTail (h : PieceInv d s) (curV d' : Nat) (hasVert isSingle : Bool) :
    PieceInv d (after (Spqr.finishTail curV d' hasVert isSingle) s) := by
  cases hasVert
  · have h₁ := h.pushVert curV d'
    cases isSingle
    · show PieceInv d (after mergeTstackTops (after (pushVertTstack curV d') s))
      exact h₁.mergeTop
    · exact h₁
  · exact h

theorem PieceInv.l1Type (h : PieceInv d s) (d' : Nat) (dir : Bool) : PieceInv d (l1S₁ d' dir s) := by
  show PieceInv d ((loop1Type d' dir).run s).2
  rw [loop1Type_run]
  split
  · show PieceInv d (after mergeTstackTops { s with stackDir := s.stackDir.set! (nxtE s).topDepth dir })
    exact (h.frame (s' := { s with stackDir := s.stackDir.set! (nxtE s).topDepth dir }) rfl rfl rfl).mergeTop
  · split <;> exact h

theorem lt_of_type_ne_F {items : Items} {x : ItemId} (h : Items.type items x ≠ .F) : x < items.size :=
  Nat.lt_of_not_le fun hle => h (Items.type_of_le items x hle)

theorem Shape.mergeTop (hs : Shape s) : Shape (after mergeTstackTops s) := by
  rw [after_mergeTstackTops]
  refine ⟨hs.size, hs.root, hs.vert, hs.edge, hs.ch_lt, fun t ht i hi => ?_⟩
  rcases hts : s.tstack with _ | ⟨b, _ | ⟨a, rest⟩⟩
  · rw [hts] at ht; simp [WalkM.mergeTop] at ht
  · rw [hts] at ht; simp [WalkM.mergeTop] at ht
  · rw [hts] at ht
    simp only [WalkM.mergeTop, List.mem_cons] at ht
    rcases ht with rfl | ht
    · simp only [List.mem_append] at hi
      rcases hi with (hi | hi) | (hi | hi)
      · exact hs.span b (by rw [hts]; simp) i (List.mem_append_left _ hi)
      · exact hs.span a (by rw [hts]; simp) i (List.mem_append_left _ hi)
      · exact hs.span a (by rw [hts]; simp) i (List.mem_append_right _ hi)
      · exact hs.span b (by rw [hts]; simp) i (List.mem_append_right _ hi)
    · exact hs.span t (by rw [hts]; simp [ht]) i hi

theorem Shape.maybeUnwrap (hs : Shape s) (ty : NodeType) {a b : TEntry} {rest : List TEntry}
    (hts : s.tstack = a :: b :: rest) : Shape (after (maybeUnwrapNxt ty) s) := by
  show Shape ((maybeUnwrapNxt ty).run s).2
  rw [maybeUnwrapNxt_run_eq ty s a b rest hts _ rfl _ rfl]
  split
  · rw [run_allocItem]; exact hs.push ty
  · split
    · refine ⟨hs.size, hs.root, hs.vert, hs.edge, hs.ch_lt, fun t ht i hi => ?_⟩
      simp only [List.mem_cons] at ht
      rcases ht with rfl | rfl | ht
      · exact hs.span t (by rw [hts]; simp) i hi
      · have key : ∀ (dir : Bool) (L : List ItemId), i ∈ (setSides dir L []).1 ++ (setSides dir L []).2 →
            i ∈ L := by
          intro dir L h; cases dir <;> simpa [setSides] using h
        have hi' : i ∈ s.items[(getSide b.spans s.stackDir[b.topDepth]!).head!]!.ch := key _ _ hi
        rw [← Items.ch_eq_getElem!] at hi'
        exact hs.ch_lt _ _ hi'
      · exact hs.span t (by rw [hts]; simp [ht]) i hi
    · rw [run_allocItem]; exact hs.push ty

/-- `maybeUnwrapNxt ty`: the invariant, and the returned item is a parentless item of type `ty`. -/
theorem PieceInv.maybeUnwrap (h : PieceInv d s) (hs : Shape s) (ty : NodeType) (hok : UnwrapOk ty s)
    {a b : TEntry} {rest : List TEntry} (hts : s.tstack = a :: b :: rest) :
    PieceInv d (after (maybeUnwrapNxt ty) s) ∧
    (∀ p, ¬ Items.IsParent (after (maybeUnwrapNxt ty) s).items p (result (maybeUnwrapNxt ty) s)) ∧
    Items.type (after (maybeUnwrapNxt ty) s).items (result (maybeUnwrapNxt ty) s) = ty := by
  have hn : nxtE s = b := by rw [nxtE, hts]; rfl
  have hd : nxtDir s = s.stackDir[b.topDepth]! := by rw [nxtDir, hn]
  have hh : nxtHead s = (getSide b.spans s.stackDir[b.topDepth]!).head! := by rw [nxtHead, hn, hd]
  show PieceInv d ((maybeUnwrapNxt ty).run s).2 ∧
    (∀ p, ¬ Items.IsParent ((maybeUnwrapNxt ty).run s).2.items p ((maybeUnwrapNxt ty).run s).1) ∧
    Items.type ((maybeUnwrapNxt ty).run s).2.items ((maybeUnwrapNxt ty).run s).1 = ty
  have alloc : PieceInv d ((allocItem ty).run s).2 ∧
      (∀ p, ¬ Items.IsParent ((allocItem ty).run s).2.items p ((allocItem ty).run s).1) ∧
      Items.type ((allocItem ty).run s).2.items ((allocItem ty).run s).1 = ty := by
    rw [run_allocItem]
    refine ⟨h.alloc hs ty, fun p hp => ?_, Items.type_push_size _⟩
    have := hs.ch_lt p _ ((Items.IsParent_congr (Items.ch_push_nil _ rfl)).1 hp)
    exact Nat.lt_irrefl _ this
  rw [maybeUnwrapNxt_run_eq ty s a b rest hts _ rfl _ rfl]
  split
  · exact alloc
  · split
    · rename_i h1 h2
      have hu : UnwrapAt s := hok.unwrap h1 (by rw [hh]; exact h2)
      refine ⟨h.frame rfl rfl rfl, fun p hp => hu.root p (by rw [hh]; exact hp), ?_⟩
      show Items.type s.items _ = ty
      rw [Items.type_eq_getElem!]; exact h2
    · exact alloc

/-- `finishTstackTop x` for a parentless S/P/R item `x`, given `FinishPiece`. -/
theorem PieceInv.finishTop (h : PieceInv d s) (hs : Shape s) (x : ItemId)
    (hx : ∀ p, ¬ Items.IsParent s.items p x) (hxt : Items.type s.items x ∉ [NodeType.F, .V, .Q])
    (hp : FinishPiece x s) {t : TEntry} {rest : List TEntry} (hts : s.tstack = t :: rest) :
    PieceInv d (after (finishTstackTop x) s) := by
  have hF : Items.type s.items x ≠ .F := fun h => hxt (by rw [h]; simp)
  have hV : Items.type s.items x ≠ .V := fun h => hxt (by rw [h]; simp)
  have hQ : Items.type s.items x ≠ .Q := fun h => hxt (by rw [h]; simp)
  have hxsz : x < s.items.size := lt_of_type_ne_F hF
  have hc : curE s = t := by rw [curE, hts, List.head!_cons']
  show PieceInv d ((finishTstackTop x).run s).2
  rw [finishTstackTop_run_eq s x t rest hts]
  refine PieceInv.frame (h.modify hs x (fun it => { it with
      vs := setSides s.stackDir[t.topDepth]! (some s.stackVerts[t.topDepth]!) (some t.vStart),
      ch := getSide t.spans s.stackDir[t.topDepth]! }) (fun _ => rfl) hx hF hV ?_ ?_) rfl rfl rfl
  · intro e he hxe; exfalso; exact hQ (hxe ▸ hs.edge e he)
  · intro hP c hcm htc
    rw [Items.ch_modify_at x _ hxsz] at hcm
    rw [Items.vs_modify_at x _ hxsz]
    by_cases hcx : c = x
    · subst hcx; rw [Items.vs_modify_at c _ hxsz]
    · rw [Items.vs_modify_of_ne x _ hcx]
      have := hp hP c (by rw [hc]; exact hcm) htc
      rw [hc] at this; exact this

theorem PieceInv.vertUnwrap (h : PieceInv d s) (hs : Shape s) {isType1 isSingle : Bool}
    (hok : isType1 = true → UnwrapOk (if isSingle then .S else .R) s) :
    PieceInv d (after (vertUnwrap isType1 isSingle) s) ∧
    (isType1 = true →
      (∀ p, ¬ Items.IsParent (after (vertUnwrap isType1 isSingle) s).items p
        ((maybeUnwrapNxt (if isSingle then NodeType.S else .R)).run s).1) ∧
      Items.type (after (vertUnwrap isType1 isSingle) s).items
        ((maybeUnwrapNxt (if isSingle then NodeType.S else .R)).run s).1 = if isSingle then .S else .R) := by
  cases isType1
  · exact ⟨h, fun h => by cases h⟩
  · obtain ⟨a, b, rest, hts⟩ := two_entries_of_le (hok rfl).two
    have hu := h.maybeUnwrap hs (if isSingle then .S else .R) (hok rfl) hts
    have heq : after (WalkState.vertUnwrap true isSingle) s =
        after (maybeUnwrapNxt (if isSingle then .S else .R)) s := by
      show ((some <$> maybeUnwrapNxt _).run s).2 = _
      rw [WalkM.map_run]; rfl
    rw [heq]
    exact ⟨hu.1, fun _ => ⟨hu.2.1, hu.2.2⟩⟩

theorem PieceInv.pSite (h : PieceInv d s) (hs : Shape s) {curV lv : Nat} {b : Bool}
    (hok : result (condP curV lv b) s = true → UnwrapOk .P s)
    (hpc : result (condP curV lv b) s = true → FinishPiece (pItem s) (pPre s)) :
    PieceInv d (after (finishP curV lv b) s) := by
  show PieceInv d ((finishP curV lv b).run s).2
  simp only [finishP, WalkM.run_bind, run_condP]
  by_cases hc : result (condP curV lv b) s = true
  · have hc' : (b && decide (s.tstack.length ≥ 2) && (s.tstack.tail.head!.vStart == curV) &&
        (s.tstack.tail.head!.topDepth == lv)) = true := hc
    simp only [hc', ↓reduceIte, WalkM.run_bind]
    have h2 : 2 ≤ s.tstack.length :=
      of_decide_eq_true (Bool.and_eq_true_iff.1 (Bool.and_eq_true_iff.1 (Bool.and_eq_true_iff.1 hc').1).1).2
    obtain ⟨a, b', rest, hts⟩ := two_entries_of_le h2
    obtain ⟨h₁, hroot₁, hty₁⟩ := h.maybeUnwrap hs .P (hok hc) hts
    have h₂ : PieceInv d (pPre s) := h₁.mergeTop
    obtain ⟨y', items, heq, -, -, -⟩ := maybeUnwrapNxt_run .P s hts
    have hts₁ : (after (maybeUnwrapNxt .P) s).tstack = a :: y' :: rest := by
      show ((maybeUnwrapNxt .P).run s).2.tstack = _; rw [heq]
    obtain ⟨y'', heq₂, -, -, -⟩ := mergeTstackTops_run _ hts₁
    have hts₂ : (pPre s).tstack = y'' :: rest := by
      show (mergeTstackTops.run _).2.tstack = _; rw [heq₂]
    have hs₂ : Shape (pPre s) := (Shape.maybeUnwrap hs .P hts).mergeTop
    refine h₂.finishTop hs₂ (pItem s) ?_ ?_ (hpc hc) hts₂
    · show ∀ p, ¬ Items.IsParent (after mergeTstackTops (after (maybeUnwrapNxt .P) s)).items p _
      rw [after_mergeTstackTops]; exact hroot₁
    · show Items.type (after mergeTstackTops (after (maybeUnwrapNxt .P) s)).items _ ∉ _
      rw [after_mergeTstackTops]
      show Items.type (after (maybeUnwrapNxt .P) s).items (result (maybeUnwrapNxt .P) s) ∉ _
      rw [hty₁]; simp
  · have hc' : (b && decide (s.tstack.length ≥ 2) && (s.tstack.tail.head!.vStart == curV) &&
        (s.tstack.tail.head!.topDepth == lv)) = false := Bool.eq_false_iff.2 hc
    simp only [hc', Bool.false_eq_true, ↓reduceIte]
    exact h

theorem Shape.tail (hs : Shape s) : Shape { s with tstack := s.tstack.tail } :=
  ⟨hs.size, hs.root, hs.vert, hs.edge, hs.ch_lt, fun t ht => hs.span t (List.mem_of_mem_tail ht)⟩

theorem Shape.modify_vs (hs : Shape s) (j : ItemId) (vsv : Option Nat × Option Nat) :
    Shape { s with items := s.items.modify j fun it => { it with vs := vsv } } :=
  hs.modify j _ (fun _ => rfl) fun hj c hc => hs.ch_lt j c (by
    show c ∈ Items.ch s.items j
    simpa [Items.ch, Array.getElem?_eq_getElem hj] using hc)

theorem Shape.modify_ch (hs : Shape s) (j : ItemId) (L : List ItemId) (hL : ∀ c ∈ L, c < s.items.size) :
    Shape { s with items := s.items.modify j fun it => { it with ch := L } } :=
  hs.modify j _ (fun _ => rfl) fun _ c hc => hL c hc

theorem Shape.retarget (hs : Shape s) (v : Nat) (dir : Bool) :
    Shape (after (WalkState.retarget v dir) s) := by
  show Shape ((WalkState.retarget v dir).run s).2
  rw [WalkState.retarget, run_modifyCur]
  refine ⟨hs.size, hs.root, hs.vert, hs.edge, hs.ch_lt, fun t ht i hi => ?_⟩
  rcases hts : s.tstack with _ | ⟨a, rest⟩
  · rw [hts] at ht; simp at ht
  · rw [hts] at ht; simp only [List.mem_cons] at ht
    rcases ht with rfl | ht
    · have : i ∈ a.spans.1 ++ a.spans.2 := by cases dir <;> simpa [setSides] using hi
      exact hs.span a (by rw [hts]; simp) i this
    · exact hs.span t (by rw [hts]; simp [ht]) i hi

theorem Shape.vertUnwrap (hs : Shape s) {isType1 isSingle : Bool}
    (hok : isType1 = true → UnwrapOk (if isSingle then .S else .R) s) :
    Shape (after (WalkState.vertUnwrap isType1 isSingle) s) := by
  cases isType1
  · exact hs
  · obtain ⟨a, b, rest, hts⟩ := two_entries_of_le (hok rfl).two
    show Shape ((some <$> maybeUnwrapNxt _).run s).2
    rw [WalkM.map_run]
    exact hs.maybeUnwrap _ hts

theorem items_retarget (v : Nat) (dir : Bool) : (after (WalkState.retarget v dir) s).items = s.items := by
  show ((WalkState.retarget v dir).run s).2.items = _
  rw [WalkState.retarget, run_modifyCur]

theorem items_mergeTop : (after mergeTstackTops s).items = s.items := by
  rw [after_mergeTstackTops]

/-- Closing a block at the Q item of edge `e`: write its children `L` (oriented `(curV, _)`) and
append it under `vertItem curV`. -/
theorem PieceInv.qClose {r : WalkState} (h : PieceInv d r) (hs : Shape r) {curV e : Nat}
    (he : e < r.g.ne) (hv : curV < r.g.nv)
    (hqroot : ∀ p, ¬ Items.IsParent r.items p (edgeItem r.g e))
    (hqvs : (Items.vs r.items (edgeItem r.g e)).1 = some curV)
    (hsvd : r.stackVerts[d]! = curV) (L : List ItemId) (hL : ∀ c ∈ L, c < r.items.size)
    (hQL : edgeItem r.g e ∉ L)
    (hL1 : ∀ c w, L = [c, vertItem w] → Items.vs r.items c = (some curV, some w))
    (hL2 : ∀ c, L = [c] → Items.vs r.items c = (some curV, none)) :
    PieceInv d { r with
      items := Array.modify (Array.modify r.items (edgeItem r.g e) (fun it => { it with ch := L }))
                 (vertItem curV) (fun it => { it with ch := it.ch ++ [edgeItem r.g e] }) } := by
  have hqt : Items.type r.items (edgeItem r.g e) = .Q := hs.edge e he
  have hqsz : edgeItem r.g e < r.items.size := by
    show 1 + r.g.nv + e < _; have := hs.size; omega
  have hvs₁ : ∀ c, Items.vs (r.items.modify (edgeItem r.g e) fun it => { it with ch := L }) c =
      Items.vs r.items c := fun c => Items.vs_modify_of_vs _ _ c _ (fun _ => rfl)
  have hty₁ : ∀ c, Items.type (r.items.modify (edgeItem r.g e) fun it => { it with ch := L }) c =
      Items.type r.items c := Items.type_modify_type_eq _ _ (fun _ => rfl)
  have h₁ : PieceInv d { r with items := r.items.modify (edgeItem r.g e) fun it => { it with ch := L } } := by
    refine h.modify hs _ _ (fun _ => rfl) hqroot (by rw [hqt]; decide) (by rw [hqt]; decide) ?_ ?_
    · intro e' _ hee
      have hee' : e' = e := by
        have : 1 + r.g.nv + e = 1 + r.g.nv + e' := hee; omega
      subst hee'
      rw [Items.ch_modify_at _ _ hqsz]
      simp only [hvs₁, hqvs]
      exact ⟨hL1, hL2⟩
    · intro hP; rw [hqt] at hP; cases hP
  have hs₁ : Shape { r with items := r.items.modify (edgeItem r.g e) fun it => { it with ch := L } } :=
    hs.modify_ch _ L hL
  refine h₁.vAppend hs₁ hv (edgeItem r.g e) ?_ (by rw [hty₁]; exact hqt) (by rw [hvs₁]; exact hqvs) ?_
  · intro p hp
    by_cases hpq : p = edgeItem r.g e
    · subst hpq
      apply hQL
      have := hp
      simp only [Items.IsParent, Items.ch_modify_at _ _ hqsz] at this
      exact this
    · have : Items.IsParent r.items p (edgeItem r.g e) := by
        have := hp
        simp only [Items.IsParent, Items.ch_modify_of_ne _ _ hpq] at this
        exact this
      exact hqroot p this
  · intro a ha
    have := h₁.root_path a ha d le_rfl
    change ¬ Items.Below _ a (vertItem r.stackVerts[d]!) at this
    rwa [hsvd] at this

/-- The boundary branch (`d ≤ lowval`): the Q item gets `vs := (some curV, none)`, its children (the
block entry's vertex item plus a bridge/loop/closed node, oriented `(curV, dest)`), and is appended
under `vertItem curV`. -/
theorem PieceInv.boundary (h : PieceInv d s) (hs : Shape s) {curV : Nat} {o : DfsOut} {hasVert : Bool}
    (he : o.e < s.g.ne) (hv : curV < s.g.nv) (hdest : o.dest < s.g.nv)
    (hq : Items.ch s.items (edgeItem s.g o.e) = [])
    (hqroot : ∀ p, ¬ Items.IsParent s.items p (edgeItem s.g o.e))
    (hsvd : s.stackVerts[d]! = curV)
    (hbridge : o.cls.isTree = true → o.cls.lowval d = d + 1 →
      s.tstack.head!.spans.2 = [vertItem o.dest] ∧
        setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) = (some curV, some o.dest))
    (hcomp : o.cls.isTree = true → o.cls.lowval d ≠ d + 1 →
      ∃ c, s.tstack.head!.spans.1 = [c] ∧ s.tstack.tail.head!.spans.2 = [vertItem o.dest] ∧
        c < s.items.size ∧ c ≠ edgeItem s.g o.e ∧ Items.vs s.items c = (some curV, some o.dest)) :
    PieceInv d (after (finishBoundary curV d o (edgeItem s.g o.e) hasVert) s) := by
  have hqt : Items.type s.items (edgeItem s.g o.e) = .Q := hs.edge o.e he
  have hqsz : edgeItem s.g o.e < s.items.size := by
    show 1 + s.g.nv + o.e < _; have := hs.size; omega
  have hqne_v : edgeItem s.g o.e ≠ vertItem o.dest := by
    show 1 + s.g.nv + o.e ≠ 1 + o.dest; omega
  have hvsz : vertItem o.dest < s.items.size := by
    show 1 + o.dest < _; have := hs.size; omega
  have h₀ : PieceInv d { s with
      items := s.items.modify (edgeItem s.g o.e) fun it => { it with vs := (some curV, none) },
      totBlocks := s.totBlocks + 1 } :=
    (h.modifyVs hs _ _ hqroot (by rw [hqt]; decide) (by rw [hqt]; decide) hq).frame rfl rfl rfl
  have hs₀ : Shape { s with
      items := s.items.modify (edgeItem s.g o.e) fun it => { it with vs := (some curV, none) },
      totBlocks := s.totBlocks + 1 } :=
    (hs.modify_vs (edgeItem s.g o.e) (some curV, none)).frame (by rfl) (by rfl) (by rfl)
  have hch₀ : ∀ p, Items.ch (s.items.modify (edgeItem s.g o.e) fun it => { it with vs := (some curV, none) }) p =
      Items.ch s.items p := Items.ch_modify_vs _ _
  have hvs₀ : Items.vs (s.items.modify (edgeItem s.g o.e) fun it => { it with vs := (some curV, none) })
      (edgeItem s.g o.e) = (some curV, none) := by rw [Items.vs_modify_at _ _ hqsz]
  have hsz₀ : (s.items.modify (edgeItem s.g o.e) fun it => { it with vs := (some curV, none) }).size =
      s.items.size := Array.size_modify ..
  show wp (finishBoundary curV d o (edgeItem s.g o.e) hasVert) (fun _ s' => PieceInv d s') s
  unfold finishBoundary
  simp only [wp_bind, wp_modifyItem, wp_modify, wp_ite, wp_allocItem, wp_makeVs, wp_popTstack, wp_pure]
  split
  · split
    · rename_i ht hl
      have hl' : o.cls.lowval d = d + 1 := by simpa using hl
      obtain ⟨ht2, hmk⟩ := hbridge ht hl'
      have hch₁ := Items.ch_push_nil
        (items := s.items.modify (edgeItem s.g o.e) fun it => { it with vs := (some curV, none) })
        ⟨.I, (none, none), []⟩ rfl
      have hroot₁ : ∀ p, ¬ Items.IsParent
          ((s.items.modify (edgeItem s.g o.e) fun it => { it with vs := (some curV, none) }).push
            ⟨.I, (none, none), []⟩) p
          (s.items.modify (edgeItem s.g o.e) fun it => { it with vs := (some curV, none) }).size := by
        intro p hp
        have := hs₀.ch_lt p _ ((Items.IsParent_congr hch₁).1 hp)
        exact Nat.lt_irrefl _ this
      have h₁ := (h₀.alloc hs₀ .I).modifyVs ((hs₀.push .I)) _
        (setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest)) hroot₁
        (by rw [Items.type_push_size]; decide) (by rw [Items.type_push_size]; decide)
        (by rw [hch₁]; exact Items.ch_of_le _ _ le_rfl)
      have hs₁ := ((hs₀.push .I).modify_vs (s.items.modify (edgeItem s.g o.e) fun it => { it with vs := (some curV, none) }).size
        (setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest))).tail
      have hqne_i : edgeItem s.g o.e ≠
          (s.items.modify (edgeItem s.g o.e) fun it => { it with vs := (some curV, none) }).size := by
        rw [hsz₀]; exact Nat.ne_of_lt hqsz
      have hch₂ : ∀ p, Items.ch (((s.items.modify (edgeItem s.g o.e) fun it => { it with vs := (some curV, none) }).push
            ⟨.I, (none, none), []⟩).modify
            (s.items.modify (edgeItem s.g o.e) fun it => { it with vs := (some curV, none) }).size
            fun it => { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) }) p =
          Items.ch s.items p := fun p => by rw [Items.ch_modify_vs, hch₁, hch₀]
      refine PieceInv.qClose (r := { s with
          items := ((s.items.modify (edgeItem s.g o.e) fun it => { it with vs := (some curV, none) }).push ⟨.I, (none, none), []⟩).modify (s.items.modify (edgeItem s.g o.e) fun it => { it with vs := (some curV, none) }).size fun it =>
            { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) },
          totBlocks := s.totBlocks + 1, tstack := s.tstack.tail })
        (h₁.frame rfl rfl rfl) hs₁ he hv ?_ ?_ hsvd _ ?_ ?_ ?_ ?_
      · intro p hp; exact hqroot p ((Items.IsParent_congr hch₂).1 hp)
      · show (Items.vs _ (edgeItem s.g o.e)).1 = some curV
        rw [Items.vs_modify_of_ne _ _ hqne_i, Items.vs_push_of_ne _ hqne_i, hvs₀]
      · show ∀ c ∈ _ :: s.tstack.head!.spans.2, c < _
        rw [ht2]
        intro c hc
        simp only [List.mem_cons, List.not_mem_nil, or_false] at hc
        rcases hc with rfl | rfl
        · dsimp only; simp only [Array.size_modify, Array.size_push]; exact Nat.lt_succ_self _
        · dsimp only; simp only [Array.size_modify, Array.size_push]; exact Nat.lt_succ_of_lt hvsz
      · show edgeItem s.g o.e ∉ _ :: s.tstack.head!.spans.2
        rw [ht2]
        simp only [List.mem_cons, List.not_mem_nil, or_false, not_or]
        exact ⟨hqne_i, hqne_v⟩
      · show ∀ c w, _ :: s.tstack.head!.spans.2 = [c, vertItem w] → _
        rw [ht2]
        intro c w hcw
        simp only [List.cons.injEq, and_true] at hcw
        obtain ⟨rfl, hw⟩ := hcw
        have hw' : w = o.dest := by
          have : 1 + o.dest = 1 + w := hw; omega
        subst hw'
        show Items.vs (Array.modify _ _ _) _ = _
        rw [Items.vs_modify_at _ _ (by simp [Array.size_push])]
        exact hmk
      · show ∀ c, _ :: s.tstack.head!.spans.2 = [c] → _
        rw [ht2]; intro c hc; simp at hc
    · rename_i ht hl
      have hl' : o.cls.lowval d ≠ d + 1 := by simpa using hl
      obtain ⟨c, hc1, ht2, hclt, hcq, hcvs⟩ := hcomp ht hl'
      refine PieceInv.qClose (r := { s with
          items := (s.items.modify (edgeItem s.g o.e) fun it => { it with vs := (some curV, none) }),
          totBlocks := s.totBlocks + 1, tstack := s.tstack.tail.tail })
        (h₀.frame rfl rfl rfl) hs₀.tail.tail he hv ?_ ?_ hsvd _ ?_ ?_ ?_ ?_
      · intro p hp; exact hqroot p ((Items.IsParent_congr hch₀).1 hp)
      · show (Items.vs _ (edgeItem s.g o.e)).1 = some curV
        rw [hvs₀]
      · show ∀ c' ∈ s.tstack.head!.spans.1 ++ s.tstack.tail.head!.spans.2, c' < _
        rw [hc1, ht2]
        intro c' hc'
        simp only [List.cons_append, List.nil_append, List.mem_cons, List.mem_singleton,
          List.not_mem_nil, or_false] at hc'
        rw [hsz₀]
        rcases hc' with rfl | rfl
        · exact hclt
        · exact hvsz
      · show edgeItem s.g o.e ∉ s.tstack.head!.spans.1 ++ s.tstack.tail.head!.spans.2
        rw [hc1, ht2]
        simp only [List.cons_append, List.nil_append, List.mem_cons, List.mem_singleton,
          List.not_mem_nil, or_false, not_or]
        exact ⟨Ne.symm hcq, hqne_v⟩
      · show ∀ c' w, s.tstack.head!.spans.1 ++ s.tstack.tail.head!.spans.2 = [c', vertItem w] → _
        rw [hc1, ht2]
        intro c' w hcw
        simp only [List.cons_append, List.nil_append, List.cons.injEq, and_true] at hcw
        obtain ⟨rfl, hw⟩ := hcw
        have hw' : w = o.dest := by
          have : 1 + o.dest = 1 + w := hw; omega
        subst hw'
        show Items.vs (Array.modify _ _ _) _ = _
        rw [Items.vs_modify_of_ne _ _ hcq]; exact hcvs
      · show ∀ c', s.tstack.head!.spans.1 ++ s.tstack.tail.head!.spans.2 = [c'] → _
        rw [hc1, ht2]; intro c' hc'; simp at hc'
  · rename_i ht
    have hch₁ := Items.ch_push_nil
      (items := s.items.modify (edgeItem s.g o.e) fun it => { it with vs := (some curV, none) })
      ⟨.O, (none, none), []⟩ rfl
    have hroot₁ : ∀ p, ¬ Items.IsParent
        ((s.items.modify (edgeItem s.g o.e) fun it => { it with vs := (some curV, none) }).push
          ⟨.O, (none, none), []⟩) p
        (s.items.modify (edgeItem s.g o.e) fun it => { it with vs := (some curV, none) }).size := by
      intro p hp
      have := hs₀.ch_lt p _ ((Items.IsParent_congr hch₁).1 hp)
      exact Nat.lt_irrefl _ this
    have h₀' : PieceInv d { s with
        items := s.items.modify (edgeItem s.g o.e) fun it => { it with vs := (some curV, none) },
        totBlocks := s.totBlocks + 1, totSelfLoops := s.totSelfLoops + 1 } := h₀.frame rfl rfl rfl
    have hs₀' : Shape { s with
        items := s.items.modify (edgeItem s.g o.e) fun it => { it with vs := (some curV, none) },
        totBlocks := s.totBlocks + 1, totSelfLoops := s.totSelfLoops + 1 } := hs₀.frame (by rfl) (by rfl) (by rfl)
    have h₁ := (h₀'.alloc hs₀' .O).modifyVs ((hs₀'.push .O)) _ (some curV, none) hroot₁
      (by rw [Items.type_push_size]; decide) (by rw [Items.type_push_size]; decide)
      (by rw [hch₁]; exact Items.ch_of_le _ _ le_rfl)
    have hs₁ := (hs₀'.push .O).modify_vs (s.items.modify (edgeItem s.g o.e) fun it => { it with vs := (some curV, none) }).size (some curV, none)
    have hqne_i : edgeItem s.g o.e ≠
        (s.items.modify (edgeItem s.g o.e) fun it => { it with vs := (some curV, none) }).size := by
      rw [hsz₀]; exact Nat.ne_of_lt hqsz
    have hch₂ : ∀ p, Items.ch (((s.items.modify (edgeItem s.g o.e) fun it => { it with vs := (some curV, none) }).push
          ⟨.O, (none, none), []⟩).modify
          (s.items.modify (edgeItem s.g o.e) fun it => { it with vs := (some curV, none) }).size
          fun it => { it with vs := (some curV, none) }) p =
        Items.ch s.items p := fun p => by rw [Items.ch_modify_vs, hch₁, hch₀]
    refine PieceInv.qClose (r := { s with
        items := ((s.items.modify (edgeItem s.g o.e) fun it => { it with vs := (some curV, none) }).push ⟨.O, (none, none), []⟩).modify (s.items.modify (edgeItem s.g o.e) fun it => { it with vs := (some curV, none) }).size fun it => { it with vs := (some curV, none) },
        totBlocks := s.totBlocks + 1, totSelfLoops := s.totSelfLoops + 1 })
      h₁ hs₁ he hv ?_ ?_ hsvd _ ?_ ?_ ?_ ?_
    · intro p hp; exact hqroot p ((Items.IsParent_congr hch₂).1 hp)
    · show (Items.vs _ (edgeItem s.g o.e)).1 = some curV
      rw [Items.vs_modify_of_ne _ _ hqne_i, Items.vs_push_of_ne _ hqne_i, hvs₀]
    · intro c hc
      simp only [List.mem_singleton] at hc
      subst hc
      dsimp only; simp only [Array.size_modify, Array.size_push]; exact Nat.lt_succ_self _
    · simp only [List.mem_singleton]; exact hqne_i
    · intro c w hcw; simp at hcw
    · intro c hc
      simp only [List.cons.injEq, and_true] at hc
      subst hc
      show Items.vs (Array.modify _ _ _) _ = _
      rw [Items.vs_modify_at _ _ (by simp [Array.size_push])]

theorem l1Ty_mem (d : Nat) (dir : Bool) (s : WalkState) : l1Ty d dir s ∉ [NodeType.F, .V, .Q] := by
  show ((loop1Type d dir).run s).1 ∉ _
  rw [loop1Type_run]
  split
  · simp
  · split <;> simp

/-- One loop-1 iteration from its merged pre-close state. -/
theorem PieceInv.l1Body (h : PieceInv d s) (d' : Nat) (dir : Bool) (hs : Shape (l1S₁ d' dir s))
    (hok : UnwrapOk (l1Ty d' dir s) (l1S₁ d' dir s))
    (hp : FinishPiece (result (maybeUnwrapNxt (l1Ty d' dir s)) (l1S₁ d' dir s))
      (after mergeTstackTops (l1S₂ d' dir s))) :
    PieceInv d (after (loop1Body d' dir) s) := by
  obtain ⟨a, b, rest, hts⟩ := two_entries_of_le hok.two
  obtain ⟨h₁, hroot, hty⟩ := (h.l1Type d' dir).maybeUnwrap hs (l1Ty d' dir s) hok hts
  have h₂ : PieceInv d (after mergeTstackTops (l1S₂ d' dir s)) := h₁.mergeTop
  obtain ⟨y', items, heq, -, -, -⟩ := maybeUnwrapNxt_run (l1Ty d' dir s) _ hts
  have hts₁ : (l1S₂ d' dir s).tstack = a :: y' :: rest := by
    show ((maybeUnwrapNxt _).run _).2.tstack = _; rw [heq]
  obtain ⟨y'', heq₂, -, -, -⟩ := mergeTstackTops_run _ hts₁
  have hts₂ : (after mergeTstackTops (l1S₂ d' dir s)).tstack = y'' :: rest := by
    show (mergeTstackTops.run _).2.tstack = _; rw [heq₂]
  have hs₂ : Shape (after mergeTstackTops (l1S₂ d' dir s)) := (hs.maybeUnwrap _ hts).mergeTop
  refine h₂.finishTop hs₂ _ ?_ ?_ hp hts₂
  · show ∀ p, ¬ Items.IsParent (after mergeTstackTops (l1S₂ d' dir s)).items p _
    rw [items_mergeTop]; exact hroot
  · show Items.type (after mergeTstackTops (l1S₂ d' dir s)).items _ ∉ _
    rw [items_mergeTop]
    show Items.type (after (maybeUnwrapNxt (l1Ty d' dir s)) (l1S₁ d' dir s)).items
      (result (maybeUnwrapNxt (l1Ty d' dir s)) (l1S₁ d' dir s)) ∉ _
    rw [hty]; exact l1Ty_mem d' dir s

theorem PieceInv.rest (h : PieceInv d s) (hs : Shape s) {curV d' lv : Nat}
    {isType1 hasVert isSingle : Bool}
    (hok : result (condP curV lv isType1) s = true → UnwrapOk .P s)
    (hp : result (condP curV lv isType1) s = true → FinishPiece (pItem s) (pPre s)) :
    PieceInv d (after (finishRest curV d' lv isType1 hasVert isSingle) s) := by
  show PieceInv d ((finishRest curV d' lv isType1 hasVert isSingle).run s).2
  simp only [finishRest, WalkM.run_bind]
  exact (h.pSite hs hok hp).finishTail curV d' hasVert isSingle

/-- `finishEdge` preserves `PieceInv`, from `CloseBase` + `CloseContent` (block entry contents) +
`ClosePiece` (orientations at the boundary / P / V / loop-1 sites). -/
theorem finishEdge_piece (h : CloseBase σ n D curV d o origTstack hasVert s)
    (hc : CloseContent curV d o origTstack hasVert s) (hp : ClosePiece curV d o origTstack hasVert s)
    (hi : PieceInv d s) :
    PieceInv d (after (finishEdge curV d o origTstack hasVert) s) := by
  have hs := h.rgs.2.1
  obtain ⟨sub, base, hlen, hE⟩ := h.book.ear
  by_cases hge : d ≤ o.cls.lowval d
  · have hge' : o.cls.lowval d ≥ d := hge
    show PieceInv d ((finishEdge curV d o origTstack hasVert).run s).2
    rw [finishEdge_eq]
    simp only [finishEdge', WalkM.run_bind, WalkM.get_run, run_stackDir, hge', ↓reduceIte]
    refine hi.boundary hs h.book.e_lt h.book.v_lt h.site.dest_lt h.book.q hE.q_root hE.sv_d ?_ ?_
    · intro ht hl
      obtain ⟨t, hsub, -, -⟩ := hE.bd_bridge ht hl
      have hts : s.tstack = t :: base := by rw [hE.tstack, hsub]; rfl
      refine ⟨?_, hp.bd_bridge ht hge hl⟩
      have := hc.bd_vert ht hge t (by rw [if_pos hl, hts]; exact Option.mem_def.2 rfl)
      rw [hts]; exact this
    · intro ht hl
      obtain ⟨t₁, t₂, hsub, -, -, -, -⟩ := hE.bd_comp ht hge hl
      have hts : s.tstack = t₁ :: t₂ :: base := by rw [hE.tstack, hsub]; rfl
      obtain ⟨c, hc1, -, -, -⟩ := hc.bd_node ht hge hl t₁ (by rw [hts]; exact Option.mem_def.2 rfl)
      have hmem : c ∈ t₁.spans.1 ++ t₁.spans.2 := List.mem_append_left _ (by rw [hc1]; simp)
      refine ⟨c, ?_, ?_, ?_, ?_, ?_⟩
      · rw [hts]; exact hc1
      · have := hc.bd_vert ht hge t₂ (by rw [if_neg hl, hts]; exact Option.mem_def.2 rfl)
        rw [hts]; exact this
      · exact hs.span t₁ (by rw [hts]; simp) c hmem
      · intro hcq; apply hE.q_free t₁ (by rw [hts]; simp); rw [← hcq]; exact hmem
      · exact hp.bd_node ht hge hl t₁ (by rw [hts]; exact Option.mem_def.2 rfl) c hc1
  have hlow : o.cls.lowval d < d := Nat.lt_of_not_le hge
  obtain ⟨lv, kind, ho, hl⟩ := ret_of_lowval_lt hlow
  have hok := h.finishOk ho hl
  have hlv : o.cls.lowval d = lv := by rw [ho]; rfl
  have hge' : ¬ (lv ≥ d) := by omega
  have hv := h.book.v_lt
  have hj : edgeItem s.g o.e < 1 + s.g.nv + s.g.ne := by
    show 1 + s.g.nv + o.e < _; have := hok.e_lt; omega
  have hqt : Items.type s.items (edgeItem s.g o.e) = .Q := hs.edge o.e hok.e_lt
  have st₀ : Step D curV s (feS₀ d o s) := Step.modifyVs h.rgs.1.inv hs (edgeItem s.g o.e) _ hj
  have hv₀ : curV < (feS₀ d o s).g.nv := by rw [st₀.g]; exact hv
  have h₀ : PieceInv d (feS₀ d o s) :=
    hi.modifyVs hs _ _ hE.q_root (by rw [hqt]; decide) (by rw [hqt]; decide) h.book.q
  have pSite : ∀ r, feRest curV d o origTstack hasVert s = r →
      result (condP curV lv o.cls.isType1) r = true → FinishPiece (pItem r) (pPre r) := by
    intro r hr hcond
    have hc' : (o.cls.isType1 && decide (r.tstack.length ≥ 2) && (r.tstack.tail.head!.vStart == curV) &&
        (r.tstack.tail.head!.topDepth == lv)) = true := hcond
    have hb : o.cls.isType1 = true :=
      (Bool.and_eq_true_iff.1 (Bool.and_eq_true_iff.1 (Bool.and_eq_true_iff.1 hc').1).1).1
    rw [hb] at hcond
    have := hp.p_site hlow hb (by rw [hr, hlv]; exact hcond)
    rwa [hr] at this
  show PieceInv d ((finishEdge curV d o origTstack hasVert).run s).2
  rw [finishEdge_eq]
  simp only [finishEdge', finishTree, finishBack, hlv, hge', ↓reduceIte, WalkM.run_bind, WalkM.get_run,
    run_stackDir, run_makeVs, run_modifyItem]
  by_cases ht : o.cls.isTree = true
  · simp only [ht, ↓reduceIte, closeVert_eq, WalkM.run_bind]
    have st₁ : Step D curV _ (feS₁ d o s) := Step.closeEars st₀.inv st₀.shape hv₀ (hok.ears ht)
    have hv₁ : curV < (feS₁ d o s).g.nv := by rw [st₁.g]; exact hv₀
    have st₂ : Step D curV _ (feS₂ d o s) := Step.mergeLate st₁.inv st₁.shape hv₁ (hok.late ht)
    have hv₂ : curV < (feS₂ d o s).g.nv := by rw [st₂.g]; exact hv₁
    have h₁ : PieceInv d (feS₁ d o s) := by
      show PieceInv d (after (loop _ (loop1Cond d) (loop1Body d s.stackDir[d]!))
        (ceS₁ o.dest d o.e (feS₀ d o s)))
      refine PieceInv.loop _ _ _ (fun _ => rfl) (fun k hk hk' => ?_) ?_
      · have hok₁ : Loop1BodyOk D d s.stackDir[d]! (l1Iter d o s k) := (hok.ears ht).body k hk
        have hadj : Loop1BodyAdj σ d s.stackDir[d]! (l1Iter d o s k) :=
          ((h.finishR.1 lv kind ho hl).ears ht).body k hk
        have st := rangesInv_l1Iter h ht hlow k (fun j hj => hk j hj.le)
        have st₁ : RgStep σ (n + 1) D curV _ (l1S₁ d s.stackDir[d]! (l1Iter d o s k)) :=
          RgStep.loop1Type st.ranges st.step.shape h.nodup (st.hσ h.rgs.2.2) hok₁.mergeS hadj.mergeS
        exact hk'.l1Body d _ st₁.step.shape hok₁.unwrap (hp.l1_site ht hlow k hk)
      · show PieceInv d ((pushEdgeTstack o.dest d o.e).run (feS₀ d o s)).2
        rw [run_pushEdgeTstack]; exact h₀.frame rfl rfl rfl
    have h₂ : PieceInv d (feS₂ d o s) := h₁.mergeLate d
    cases hasVert
    · simp only [Bool.false_eq_true, ↓reduceIte]
      have hr : feRest curV d o origTstack false s = feS₂ d o s := by simp [feRest, ht]
      exact h₂.rest (st₀.trans (st₁.trans st₂)).shape
        (fun hc' => ((hok.rest_tree ht rfl).p.ok hc').1) (pSite _ hr)
    · simp only [↓reduceIte, WalkM.run_bind]
      have hcv := hok.vert ht rfl
      have st₃ : Step D curV _ (cvS₁ o.cls.isType1 origTstack (feSingle d o s) (feS₂ d o s)) :=
        Step.vertPre st₂.inv st₂.shape hv₂ hcv.loop3
      have h₃ := h₂.vertPre o.cls.isType1 origTstack (feSingle d o s)
      have hok' : o.cls.isType1 = true →
          UnwrapOk (if cvB₁ o.cls.isType1 origTstack (feSingle d o s) (feS₂ d o s) then .S else .R)
            (cvS₁ o.cls.isType1 origTstack (feSingle d o s) (feS₂ d o s)) :=
        fun h => by rw [h]; exact hcv.unwrap h
      have h₄' := h₃.vertUnwrap st₃.shape hok'
      have h₄ : PieceInv d (cvS₂ o.cls.isType1 origTstack (feSingle d o s) (feS₂ d o s)) := h₄'.1
      have hs₄ : Shape (cvS₂ o.cls.isType1 origTstack (feSingle d o s) (feS₂ d o s)) :=
        st₃.shape.vertUnwrap hok'
      have h₅ : PieceInv d (cvS₃ o.cls.isType1 origTstack (feSingle d o s) (feS₂ d o s)) := by
        show PieceInv d (after mergeTstackTops _); exact h₄.mergeTop
      have hs₅ : Shape (cvS₃ o.cls.isType1 origTstack (feSingle d o s) (feS₂ d o s)) := by
        show Shape (after mergeTstackTops _); exact hs₄.mergeTop
      have h₆ : PieceInv d (cvS₄ o.cls.isType1 origTstack (feSingle d o s) (feS₂ d o s)) := by
        show PieceInv d (after mergeTstackTops _); exact h₅.mergeTop
      have hs₆ : Shape (cvS₄ o.cls.isType1 origTstack (feSingle d o s) (feS₂ d o s)) := by
        show Shape (after mergeTstackTops _); exact hs₅.mergeTop
      have h₇ := h₆.retarget curV s.stackDir[d]!
      have hs₇ := hs₆.retarget curV s.stackDir[d]!
      have h₈ : PieceInv d (feS₃ curV d o origTstack s) := by
        cases h1 : o.cls.isType1
        · have hS₃ : feS₃ curV d o origTstack s =
              cvS₅ curV s.stackDir[d]! false origTstack (feSingle d o s) (feS₂ d o s) := by
            simp only [feS₃, h1, closeVert', cvS₅, cvS₄, cvS₃, cvS₂, cvS₁, after, vertPre,
              vertUnwrap, vertFinish, Bool.not_false, Bool.false_eq_true, ↓reduceIte, WalkM.run_bind,
              WalkM.pure_run, run_tstackSize]
          rw [hS₃]; rw [h1] at h₇; exact h₇
        · have hS₃ : feS₃ curV d o origTstack s =
              after (finishTstackTop
                ((maybeUnwrapNxt (if feSingle d o s then .S else .R)).run (feS₂ d o s)).1)
                (cvS₅ curV s.stackDir[d]! true origTstack (feSingle d o s) (feS₂ d o s)) := by
            simp only [feS₃, h1, closeVert', cvS₅, cvS₄, cvS₃, cvS₂, cvS₁, cvB₁, after, result, vertPre,
              vertUnwrap, vertFinish, Bool.not_true, Bool.false_eq_true, ↓reduceIte, WalkM.run_bind,
              WalkM.pure_run, WalkM.map_run]
            rfl
          have hne := (hcv.finish h1).nonempty
          rw [h1] at h₇ hne hs₇ h₄'
          obtain ⟨t, rest, hts⟩ := List.exists_cons_of_ne_nil hne
          obtain ⟨hroot, hty⟩ := h₄'.2 rfl
          have hi₅ : (cvS₅ curV s.stackDir[d]! true origTstack (feSingle d o s) (feS₂ d o s)).items =
              (cvS₂ true origTstack (feSingle d o s) (feS₂ d o s)).items := by
            show (after (WalkState.retarget _ _) (after mergeTstackTops (after mergeTstackTops _))).items = _
            rw [items_retarget, items_mergeTop, items_mergeTop]
          rw [hS₃]
          refine h₇.finishTop hs₇ _ ?_ ?_ (hp.v_site ht hlow rfl h1) hts
          · show ∀ p, ¬ Items.IsParent
              (cvS₅ curV s.stackDir[d]! true origTstack (feSingle d o s) (feS₂ d o s)).items p _
            rw [hi₅]; exact hroot
          · show Items.type
              (cvS₅ curV s.stackDir[d]! true origTstack (feSingle d o s) (feS₂ d o s)).items _ ∉ _
            rw [hi₅]
            show Items.type (after
                (WalkState.vertUnwrap true (cvB₁ true origTstack (feSingle d o s) (feS₂ d o s)))
                (cvS₁ true origTstack (feSingle d o s) (feS₂ d o s))).items
              ((maybeUnwrapNxt (if cvB₁ true origTstack (feSingle d o s) (feS₂ d o s) then NodeType.S
                else .R)).run (cvS₁ true origTstack (feSingle d o s) (feS₂ d o s))).1 ∉ _
            rw [hty]; split <;> simp
      have st₈ : Step D curV _ (feS₃ curV d o origTstack s) :=
        Step.closeVert' st₂.inv st₂.shape hv₂ hcv
      have hr : feRest curV d o origTstack true s = feS₃ curV d o origTstack s := by simp [feRest, ht]
      exact h₈.rest (st₀.trans (st₁.trans (st₂.trans st₈))).shape
        (fun hc' => ((hok.rest_vert ht rfl).p.ok hc').1) (pSite _ hr)
  · have ht' : o.cls.isTree = false := Bool.eq_false_iff.2 ht
    simp only [ht', Bool.false_eq_true, ↓reduceIte, WalkM.run_bind]
    have hq : Items.ch (feS₀ d o s).items (edgeItem (feS₀ d o s).g o.e) = [] := by
      show Items.ch (s.items.modify _ _) (edgeItem s.g o.e) = []
      rw [Items.ch_modify_ch_eq (edgeItem s.g o.e)
        (fun it => { it with vs := setSides s.stackDir[d]! (some s.stackVerts[d]!) (some o.dest) })
        (fun _ => rfl)]
      exact hok.q ht'
    have st₁ : Step D curV _ (after (pushEdgeTstack curV lv o.e) (feS₀ d o s)) :=
      Step.pushEdge st₀.inv st₀.shape curV lv o.e hok.e_lt hq (hok.ends ht') (hok.lv_le ht')
    have st₂ : Step D curV _ (feBack curV lv d o s) :=
      Step.frame st₁.inv st₁.shape rfl rfl rfl rfl
    have hB : PieceInv d (feBack curV lv d o s) := by
      have h₁ : PieceInv d (after (pushEdgeTstack curV lv o.e) (feS₀ d o s)) := by
        show PieceInv d ((pushEdgeTstack curV lv o.e).run (feS₀ d o s)).2
        rw [run_pushEdgeTstack]; exact h₀.frame rfl rfl rfl
      exact h₁.frame rfl rfl rfl
    have hr : feRest curV d o origTstack hasVert s = feBack curV lv d o s := by simp [feRest, ht', hlv]
    exact hB.rest (st₀.trans (st₁.trans st₂)).shape
      (fun hc' => ((hok.rest_back ht').p.ok hc').1) (pSite _ hr)

end WalkState
end Spqr
