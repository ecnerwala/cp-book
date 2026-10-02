import Lean
import Spqr.Walk
import Spqr.CatList
import Spqr.Refine
import Mathlib.Data.List.Induction

/-!
# Phase 2, linear time

`Spqr.Walk` with the ear spans as `Span`s (`O(1)` append, flattened once when they become a
node's children), the tstack as an `Array`, and a step counter `ticks`. Everything else is the
same program, and `WalkState.toSlow` projects a fast state onto the reference one;
`walkFast_toSlow` shows the two walks agree under that projection.
-/

namespace Spqr.Fast

structure Item where
  type : NodeType
  vs : Option Nat × Option Nat := (none, none)
  ch : Span ItemId := .nil
deriving Repr, Inhabited

def Item.toSlow (it : Item) : Spqr.Item := ⟨it.type, it.vs, it.ch.toList⟩

def initialItems (g : Graph) : Array Item :=
  #[⟨.F, (none, none), .nil⟩]
    ++ Array.replicate g.nv ⟨.V, (none, none), .nil⟩
    ++ Array.replicate g.ne ⟨.Q, (none, none), .nil⟩

structure TEntry where
  vStart : Nat
  topDepth : Nat
  firstIdx : Nat
  spans : Span ItemId × Span ItemId
deriving Repr, Inhabited

def TEntry.toSlow (t : TEntry) : Spqr.TEntry :=
  ⟨t.vStart, t.topDepth, t.firstIdx, Prod.map Span.toList Span.toList t.spans⟩

structure WalkState where
  g : Graph
  ternarize : Bool
  items : Array Item
  stackVerts : Array Nat
  stackDir : Array Bool
  nxtEdgeIdx : Nat := 0
  firstOccurrence : Array Nat
  /-- The ear stack; the last entry is the top (`cur`), the one before it is `nxt`. -/
  tstack : Array TEntry := #[]
  totBlocks : Nat := 0
  totSelfLoops : Nat := 0
  /-- Steps: loop iterations, tstack pushes / pops, item allocations. Only read by the
  complexity theorems. -/
  ticks : Nat := 0

namespace WalkState

def init (g : Graph) (ternarize : Bool) : WalkState where
  g := g
  ternarize := ternarize
  items := initialItems g
  stackVerts := Array.replicate g.nv 0
  stackDir := Array.replicate g.nv false
  firstOccurrence := Array.replicate g.nv 0

end WalkState

/-- The reference tstack (top first) of an array tstack (top last). -/
def tstackList (a : Array TEntry) : List Spqr.TEntry := (a.toList.map TEntry.toSlow).reverse

def WalkState.toSlow (s : WalkState) : Spqr.WalkState :=
  { g := s.g, ternarize := s.ternarize, items := s.items.map Item.toSlow, stackVerts := s.stackVerts,
    stackDir := s.stackDir, nxtEdgeIdx := s.nxtEdgeIdx, firstOccurrence := s.firstOccurrence,
    tstack := tstackList s.tstack, totBlocks := s.totBlocks, totSelfLoops := s.totSelfLoops }

abbrev WalkM := StateM WalkState

namespace WalkM

def tick : WalkM Unit := modify fun s => { s with ticks := s.ticks + 1 }

def modifyItem (item : ItemId) (f : Item → Item) : WalkM Unit :=
  modify fun s => { s with items := s.items.modify item f }

def getItem (item : ItemId) : WalkM Item := do return (← get).items[item]!

def allocItem (type : NodeType) : WalkM ItemId :=
  modifyGet fun s =>
    (s.items.size, { s with items := s.items.push ⟨type, (none, none), .nil⟩, ticks := s.ticks + 1 })

def stackDir (d : Nat) : WalkM Bool := do return (← get).stackDir[d]!
def setStackDir (d : Nat) (b : Bool) : WalkM Unit :=
  modify fun s => { s with stackDir := s.stackDir.set! d b }

def makeVs (vStart topDepth : Nat) : WalkM (Option Nat × Option Nat) := do
  let s ← get
  return setSides s.stackDir[topDepth]! (some s.stackVerts[topDepth]!) (some vStart)

def cur : WalkM TEntry := do return (← get).tstack.back!
def nxt : WalkM TEntry := do
  let s ← get
  return if 2 ≤ s.tstack.size then s.tstack[s.tstack.size - 2]! else default
def tstackSize : WalkM Nat := do return (← get).tstack.size
def modifyCur (f : TEntry → TEntry) : WalkM Unit :=
  modify fun s => { s with tstack := s.tstack.modify (s.tstack.size - 1) f }
def modifyNxt (f : TEntry → TEntry) : WalkM Unit :=
  modify fun s => { s with
    tstack := if 2 ≤ s.tstack.size then s.tstack.modify (s.tstack.size - 2) f else s.tstack }
def popTstack : WalkM TEntry :=
  modifyGet fun s => (s.tstack.back!, { s with tstack := s.tstack.pop, ticks := s.ticks + 1 })

def pushTstack (vStart topDepth : Nat) (item : ItemId) : WalkM Unit :=
  modify fun s =>
    { s with tstack := s.tstack.push ⟨vStart, topDepth, s.nxtEdgeIdx, setSides s.stackDir[topDepth]! (.single item) .nil⟩,
             ticks := s.ticks + 1 }
def pushVertTstack (v topDepth : Nat) : WalkM Unit := pushTstack v topDepth (vertItem v)
def pushEdgeTstack (vStart topDepth e : Nat) : WalkM Unit := do
  pushTstack vStart topDepth (edgeItem (← get).g e)

/-- Merge the top ear into the one below it. -/
def mergeTstackTops : WalkM Unit := do
  let b ← popTstack
  modifyCur fun a =>
    { a with topDepth := min a.topDepth b.topDepth,
             spans := (b.spans.1 ++ a.spans.1, a.spans.2 ++ b.spans.2) }

/-- Allocate a node of the given type for the ear below the top, unless (when not ternarizing)
that ear is a single node of the same type, which is then reopened and reused. -/
def maybeUnwrapNxt (type : NodeType) : WalkM ItemId := do
  if type == .R || (← get).ternarize then return ← allocItem type
  let t ← nxt
  let topDir ← stackDir t.topDepth
  let item := (getSide t.spans topDir).head!
  let it ← getItem item
  if it.type == type then
    modifyNxt fun t => { t with spans := setSides topDir it.ch .nil }
    return item
  else
    allocItem type

/-- Close the top ear into `item`: it becomes the ear's single child. -/
def finishTstackTop (item : ItemId) : WalkM Unit := do
  let t ← cur
  let topDir ← stackDir t.topDepth
  let vs ← makeVs t.vStart t.topDepth
  modifyItem item fun it => { it with vs := vs, ch := getSide t.spans topDir }
  modifyCur fun t => { t with spans := setSides topDir (.single item) .nil }

/-- `while cond do body`, with `fuel` bounding the number of iterations. -/
def loop : Nat → WalkM Bool → WalkM Unit → WalkM Unit
  | 0, _, _ => pure ()
  | fuel + 1, cond, body => do
    if ← cond then
      body
      tick
      loop fuel cond body

/-- Loop 1 of `finishEdge`: while the ear below the top stays at or below depth `d`, close it. -/
def closeCond (d : Nat) : WalkM Bool := do return (← tstackSize) ≥ 2 && (← nxt).topDepth ≥ d
def closeBody (d : Nat) (edgeDir : Bool) : WalkM Unit := do
  let type ← do
    if (← nxt).topDepth > d then
      setStackDir (← nxt).topDepth edgeDir
      mergeTstackTops
      pure NodeType.S
    else if (← nxt).vStart == (← cur).vStart then pure .P
    else pure .R
  let item ← maybeUnwrapNxt type
  mergeTstackTops
  finishTstackTop item

/-- Loop 2 of `finishEdge`: while the top ear's first back edge is after `fo`, merge it down. -/
def firstIdxCond (fo : Nat) : WalkM Bool := do return (← cur).firstIdx > fo
/-- Loop 3 of `finishEdge`: while more than `n` ears are stacked, merge the top two. -/
def sizeCond (n : Nat) : WalkM Bool := do return (← tstackSize) > n

end WalkM

open Spqr.Fast.WalkM in
/-- `Spqr.finishEdge`, with the three loops named. -/
def finishEdge (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool) : WalkM Bool := do
  tick
  let g := (← get).g
  let nxtV := o.dest
  let e := o.e
  let lowval := o.cls.lowval d
  let isTree := o.cls.isTree
  let isType1 := o.cls.isType1
  let edgeDir ← stackDir d
  let qItem := edgeItem g e

  if lowval ≥ d then
    modifyItem qItem fun it => { it with vs := (some curV, none) }
    modify fun s => { s with totBlocks := s.totBlocks + 1 }
    if isTree then
      if lowval == d + 1 then
        let item ← allocItem .I
        let vs ← makeVs nxtV d
        modifyItem item fun it => { it with vs := vs }
        let t ← popTstack
        modifyItem qItem fun it => { it with ch := .single item ++ t.spans.2 }
      else
        let backedge ← popTstack
        let t ← popTstack
        modifyItem qItem fun it => { it with ch := backedge.spans.1 ++ t.spans.2 }
    else
      modify fun s => { s with totSelfLoops := s.totSelfLoops + 1 }
      let item ← allocItem .O
      modifyItem item fun it => { it with vs := (some curV, none) }
      modifyItem qItem fun it => { it with ch := .single item }
    modifyItem (vertItem curV) fun it => { it with ch := it.ch ++ .single qItem }
    return hasVert

  let vs ← makeVs nxtV d
  modifyItem qItem fun it => { it with vs := vs }

  let mut isSingle := true
  if isTree then
    pushEdgeTstack nxtV d e
    loop (← tstackSize) (closeCond d) (closeBody d edgeDir)
    let fo := (← get).firstOccurrence[d]!
    if (← cur).firstIdx > fo then
      loop (← tstackSize) (firstIdxCond fo) mergeTstackTops
      isSingle := false
    if hasVert then
      if !isType1 then
        loop (← tstackSize) (sizeCond (origTstack + 3)) mergeTstackTops
        isSingle := false
      let item ← if isType1 then some <$> maybeUnwrapNxt (if isSingle then .S else .R) else pure none
      mergeTstackTops
      mergeTstackTops
      modifyCur fun t => { t with vStart := curV, spans := setSides (!edgeDir) (t.spans.1 ++ t.spans.2) .nil }
      if let some item := item then
        finishTstackTop item
        isSingle := true
  else
    pushEdgeTstack curV lowval e
    modify fun s =>
      { s with firstOccurrence := s.firstOccurrence.modify lowval (min · s.nxtEdgeIdx),
               nxtEdgeIdx := s.nxtEdgeIdx + 1 }

  if isType1 && (← tstackSize) ≥ 2 && (← nxt).vStart == curV && (← nxt).topDepth == lowval then
    let item ← maybeUnwrapNxt .P
    mergeTstackTops
    finishTstackTop item

  if !hasVert then
    pushVertTstack curV d
    if !isSingle then mergeTstackTops
    return true
  return hasVert

open Spqr.Fast.WalkM in
mutual
def walkTree (t : DfsTree) (d : Nat) : WalkM Unit := do
  match t with
  | .node v outs =>
    tick
    modify fun s => { s with stackVerts := s.stackVerts.set! d v }
    let hasVert ← walkOuts v d outs false
    unless hasVert do
      setStackDir d true
      pushVertTstack v d

def walkOuts (v d : Nat) (outs : List DfsOut) (hasVert : Bool) : WalkM Bool := do
  match outs with
  | [] => return hasVert
  | o :: rest =>
    let hasVert ← walkOut v d o hasVert
    walkOuts v d rest hasVert

def walkOut (v d : Nat) (o : DfsOut) (hasVert : Bool) : WalkM Bool := do
  let lowval := o.cls.lowval d
  let lowDir ← stackDir lowval
  setStackDir d (if lowval ≥ d then false else !lowDir)
  let hasVert ← do
    if !hasVert && lowval < d && o.cls.isType1 then
      pushVertTstack v d
      pure true
    else pure hasVert
  let origTstack ← tstackSize
  match o with
  | .tree _ _ child =>
    modify fun s => { s with firstOccurrence := s.firstOccurrence.set! d s.g.ne }
    walkTree child (d + 1)
  | .back .. => pure ()
  finishEdge v d o origTstack hasVert
end

open Spqr.Fast.WalkM in
def walkForest (forest : List DfsTree) : WalkM Unit :=
  forest.forM fun t => do
    walkTree t 0
    let top ← popTstack
    modifyItem rootItem fun it => { it with ch := it.ch ++ top.spans.2 }

end Fast

/-- Phase 2 entry point, linear time. -/
def Graph.walkFast (g : Graph) (ternarize : Bool) (forest : List DfsTree) : Fast.WalkState :=
  (Fast.walkForest forest).run (Fast.WalkState.init g ternarize) |>.2

/-! ## `walkFast` simulates `walk` -/

namespace Fast
open StateRun

attribute [sim_simp] Item.toSlow TEntry.toSlow WalkState.toSlow Prod.map_fst Prod.map_snd Span.head!_toList

@[sim_simp] theorem getSide_map {α β : Type} (f : α → β) (p : α × α) (d : Bool) :
    getSide (Prod.map f f p) d = f (getSide p d) := by
  cases d <;> rfl

@[sim_simp] theorem Prod.map_setSides {α β : Type} (f : α → β) (d : Bool) (a b : α) :
    Prod.map f f (setSides d a b) = setSides d (f a) (f b) := by
  cases d <;> rfl

section Lemmas

@[sim_simp] theorem WalkState.toSlow_g (s : WalkState) : s.toSlow.g = s.g := rfl
@[sim_simp] theorem WalkState.toSlow_ternarize (s : WalkState) : s.toSlow.ternarize = s.ternarize := rfl
@[sim_simp] theorem WalkState.toSlow_items (s : WalkState) : s.toSlow.items = s.items.map Item.toSlow := rfl
@[sim_simp] theorem WalkState.toSlow_stackVerts (s : WalkState) : s.toSlow.stackVerts = s.stackVerts := rfl
@[sim_simp] theorem WalkState.toSlow_stackDir (s : WalkState) : s.toSlow.stackDir = s.stackDir := rfl
@[sim_simp] theorem WalkState.toSlow_nxtEdgeIdx (s : WalkState) : s.toSlow.nxtEdgeIdx = s.nxtEdgeIdx := rfl
@[sim_simp] theorem WalkState.toSlow_firstOccurrence (s : WalkState) :
    s.toSlow.firstOccurrence = s.firstOccurrence := rfl
@[sim_simp] theorem WalkState.toSlow_tstack (s : WalkState) : s.toSlow.tstack = tstackList s.tstack := rfl
@[sim_simp] theorem WalkState.toSlow_totBlocks (s : WalkState) : s.toSlow.totBlocks = s.totBlocks := rfl
@[sim_simp] theorem WalkState.toSlow_totSelfLoops (s : WalkState) : s.toSlow.totSelfLoops = s.totSelfLoops := rfl

theorem Item.toSlow_default : (default : Item).toSlow = default := rfl
theorem TEntry.toSlow_default : (default : TEntry).toSlow = default := rfl

theorem initialItems_toSlow (g : Graph) : (initialItems g).map Item.toSlow = Spqr.initialItems g := by
  simp [initialItems, Spqr.initialItems, Item.toSlow]

theorem Array.getElem!_map' {α β : Type} [Inhabited α] [Inhabited β] (f : α → β) (hf : f default = default)
    (a : Array α) (i : Nat) : (a.map f)[i]! = f a[i]! := by
  by_cases h : i < a.size
  · simp [getElem!_pos, h]
  · simp [getElem!_neg, h, hf]

theorem Array.map_modify {α β : Type} (f : α → β) (a : Array α) (i : Nat) (gf : α → α) (g : β → β)
    (h : ∀ x, f (gf x) = g (f x)) : (a.modify i gf).map f = (a.map f).modify i g := by
  apply Array.ext
  · simp
  · intro j h1 h2
    simp only [Array.getElem_map, Array.getElem_modify]
    split <;> simp [h]

@[simp] theorem tstackList_push (a : Array TEntry) (t : TEntry) :
    tstackList (a.push t) = t.toSlow :: tstackList a := by
  simp [tstackList]

@[simp] theorem tstackList_length (a : Array TEntry) : (tstackList a).length = a.size := by
  simp [tstackList]

@[simp] theorem tstackList_empty : tstackList #[] = [] := rfl

theorem tstackList_pop (a : Array TEntry) : tstackList a.pop = (tstackList a).tail := by
  simp [tstackList, List.tail_reverse, List.map_dropLast]

theorem Array.eq_empty_or_push {α : Type} (a : Array α) :
    a = #[] ∨ ∃ (b : Array α) (x : α), a = b.push x := by
  rcases a with ⟨l⟩
  induction l using List.reverseRecOn with
  | nil => left; rfl
  | append_singleton l x _ => right; exact ⟨⟨l⟩, x, Array.ext' (by simp)⟩

theorem tstackList_back! (a : Array TEntry) : a.back!.toSlow = (tstackList a).head! := by
  rcases Array.eq_empty_or_push a with rfl | ⟨b, x, rfl⟩
  · rfl
  · rw [Array.back!_push, tstackList_push]; rfl

theorem tstackList_nxt (a : Array TEntry) :
    (if 2 ≤ a.size then a[a.size - 2]! else default).toSlow = (tstackList a).tail.head! := by
  rcases Array.eq_empty_or_push a with rfl | ⟨b, x, rfl⟩
  · rfl
  rcases Array.eq_empty_or_push b with rfl | ⟨c, y, rfl⟩
  · rfl
  have h2 : 2 ≤ c.size + 1 + 1 := by omega
  simp only [Array.size_push, tstackList_push, List.tail_cons, h2, ↓reduceIte]
  have h3 : c.size < ((c.push y).push x).size := by simp only [Array.size_push]; omega
  have h4 : c.size < (c.push y).size := by simp only [Array.size_push]; omega
  rw [show c.size + 1 + 1 - 2 = c.size by omega, getElem!_pos ((c.push y).push x) c.size h3,
    Array.getElem_push_lt h4, Array.getElem_push_eq]
  rfl

theorem Array.modify_push_last {α : Type} (a : Array α) (x : α) (f : α → α) :
    (a.push x).modify a.size f = a.push (f x) := by
  apply Array.ext
  · simp
  · intro j h1 h2
    simp only [Array.getElem_modify, Array.getElem_push]
    split <;> split <;> first | rfl | omega | (simp_all; try omega)

theorem Array.modify_push_push {α : Type} (a : Array α) (x y : α) (f : α → α) :
    ((a.push y).push x).modify a.size f = (a.push (f y)).push x := by
  apply Array.ext
  · simp
  · intro j h1 h2
    simp only [Array.getElem_modify, Array.getElem_push, Array.size_push]
    split <;> split <;> (try split) <;> first | rfl | omega

theorem tstackList_modifyCur (a : Array TEntry) (ff : TEntry → TEntry) (f : Spqr.TEntry → Spqr.TEntry)
    (h : ∀ t, (ff t).toSlow = f t.toSlow) :
    tstackList (a.modify (a.size - 1) ff) =
      match tstackList a with | t :: rest => f t :: rest | [] => [] := by
  rcases Array.eq_empty_or_push a with rfl | ⟨b, x, rfl⟩
  · rfl
  simp only [Array.size_push, Nat.add_sub_cancel, Array.modify_push_last, tstackList_push, h]

theorem tstackList_modifyNxt (a : Array TEntry) (ff : TEntry → TEntry) (f : Spqr.TEntry → Spqr.TEntry)
    (h : ∀ t, (ff t).toSlow = f t.toSlow) :
    tstackList (if 2 ≤ a.size then a.modify (a.size - 2) ff else a) =
      match tstackList a with | t :: u :: rest => t :: f u :: rest | l => l := by
  rcases Array.eq_empty_or_push a with rfl | ⟨b, x, rfl⟩
  · rfl
  rcases Array.eq_empty_or_push b with rfl | ⟨c, y, rfl⟩
  · rfl
  have h2 : 2 ≤ c.size + 1 + 1 := by omega
  simp only [Array.size_push, h2, ↓reduceIte]
  rw [show c.size + 1 + 1 - 2 = c.size by omega, Array.modify_push_push]
  simp [h]

end Lemmas

namespace WalkM

theorem sim_modifyItem (i : ItemId) (ff : Item → Item) (f : Spqr.Item → Spqr.Item)
    (h : ∀ it, (ff it).toSlow = f it.toSlow) :
    Refine.Sim WalkState.toSlow Eq (modifyItem i ff) (Spqr.WalkM.modifyItem i f) := fun s =>
  ⟨by simp only [run_simp, modifyItem, Spqr.WalkM.modifyItem, WalkState.toSlow, Array.map_modify _ _ _ _ _ h],
   rfl⟩

theorem sim_getItem (i : ItemId) :
    Refine.Sim WalkState.toSlow (fun a b => b = a.toSlow) (getItem i) (Spqr.WalkM.getItem i) := fun s =>
  ⟨rfl, Array.getElem!_map' _ Item.toSlow_default s.items i⟩

theorem sim_allocItem (t : NodeType) :
    Refine.Sim WalkState.toSlow Eq (allocItem t) (Spqr.WalkM.allocItem t) := fun s =>
  ⟨by simp [run_simp, allocItem, Spqr.WalkM.allocItem, WalkState.toSlow, Item.toSlow] <;> rfl,
   by show s.items.size = (s.items.map Item.toSlow).size; simp⟩

theorem sim_stackDir (d : Nat) : Refine.Sim WalkState.toSlow Eq (stackDir d) (Spqr.WalkM.stackDir d) :=
  fun _ => ⟨rfl, rfl⟩
theorem sim_setStackDir (d : Nat) (b : Bool) :
    Refine.Sim WalkState.toSlow Eq (setStackDir d b) (Spqr.WalkM.setStackDir d b) := fun _ => ⟨rfl, rfl⟩
theorem sim_makeVs (v d : Nat) : Refine.Sim WalkState.toSlow Eq (makeVs v d) (Spqr.WalkM.makeVs v d) :=
  fun _ => ⟨rfl, rfl⟩
theorem sim_tstackSize : Refine.Sim WalkState.toSlow Eq tstackSize Spqr.WalkM.tstackSize := fun s =>
  ⟨rfl, (tstackList_length s.tstack).symm⟩

theorem sim_cur : Refine.Sim WalkState.toSlow (fun a b => b = a.toSlow) cur Spqr.WalkM.cur := fun s =>
  ⟨rfl, (tstackList_back! s.tstack).symm⟩

theorem sim_nxt : Refine.Sim WalkState.toSlow (fun a b => b = a.toSlow) nxt Spqr.WalkM.nxt := fun s =>
  ⟨rfl, (tstackList_nxt s.tstack).symm⟩

theorem sim_modifyCur (ff : TEntry → TEntry) (f : Spqr.TEntry → Spqr.TEntry)
    (h : ∀ t, (ff t).toSlow = f t.toSlow) :
    Refine.Sim WalkState.toSlow Eq (modifyCur ff) (Spqr.WalkM.modifyCur f) := fun s =>
  ⟨by simp only [run_simp, modifyCur, Spqr.WalkM.modifyCur, WalkState.toSlow, tstackList_modifyCur _ _ _ h]
      <;> rfl,
   rfl⟩

theorem sim_modifyNxt (ff : TEntry → TEntry) (f : Spqr.TEntry → Spqr.TEntry)
    (h : ∀ t, (ff t).toSlow = f t.toSlow) :
    Refine.Sim WalkState.toSlow Eq (modifyNxt ff) (Spqr.WalkM.modifyNxt f) := fun s =>
  ⟨by simp only [run_simp, modifyNxt, Spqr.WalkM.modifyNxt, WalkState.toSlow, tstackList_modifyNxt _ _ _ h]
      <;> rfl,
   rfl⟩

theorem sim_popTstack :
    Refine.Sim WalkState.toSlow (fun a b => b = a.toSlow) popTstack Spqr.WalkM.popTstack := fun s =>
  ⟨by simp only [run_simp, popTstack, Spqr.WalkM.popTstack, WalkState.toSlow, tstackList_pop] <;> rfl,
   (tstackList_back! s.tstack).symm⟩

theorem sim_pushTstack (v d i : Nat) :
    Refine.Sim WalkState.toSlow Eq (pushTstack v d i) (Spqr.WalkM.pushTstack v d i) := fun s =>
  ⟨by simp only [run_simp, pushTstack, Spqr.WalkM.pushTstack, WalkState.toSlow, tstackList_push]
      simp [TEntry.toSlow, setSides]
      split <;> simp_all,
   rfl⟩

theorem sim_tick (s : WalkState) : (tick.run s).2.toSlow = s.toSlow := rfl

end WalkM

sim_silent Spqr.Fast.WalkM.tick

macro_rules
  | `(tactic| sim_leaf) => `(tactic| first
      | apply WalkM.sim_modifyItem | apply WalkM.sim_getItem | apply WalkM.sim_allocItem
      | apply WalkM.sim_stackDir | apply WalkM.sim_setStackDir | apply WalkM.sim_makeVs
      | apply WalkM.sim_tstackSize | apply WalkM.sim_cur | apply WalkM.sim_nxt
      | apply WalkM.sim_modifyCur | apply WalkM.sim_modifyNxt | apply WalkM.sim_popTstack
      | apply WalkM.sim_pushTstack)

namespace WalkM

theorem sim_pushVertTstack (v d : Nat) :
    Refine.Sim WalkState.toSlow Eq (pushVertTstack v d) (Spqr.WalkM.pushVertTstack v d) := by
  unfold pushVertTstack Spqr.WalkM.pushVertTstack; sim_auto
theorem sim_pushEdgeTstack (v d e : Nat) :
    Refine.Sim WalkState.toSlow Eq (pushEdgeTstack v d e) (Spqr.WalkM.pushEdgeTstack v d e) := by
  unfold pushEdgeTstack Spqr.WalkM.pushEdgeTstack; sim_auto
theorem sim_mergeTstackTops :
    Refine.Sim WalkState.toSlow Eq mergeTstackTops Spqr.WalkM.mergeTstackTops := by
  unfold mergeTstackTops Spqr.WalkM.mergeTstackTops; sim_auto
theorem sim_maybeUnwrapNxt (t : NodeType) :
    Refine.Sim WalkState.toSlow Eq (maybeUnwrapNxt t) (Spqr.WalkM.maybeUnwrapNxt t) := by
  unfold maybeUnwrapNxt Spqr.WalkM.maybeUnwrapNxt; sim_auto
theorem sim_finishTstackTop (i : ItemId) :
    Refine.Sim WalkState.toSlow Eq (finishTstackTop i) (Spqr.WalkM.finishTstackTop i) := by
  unfold finishTstackTop Spqr.WalkM.finishTstackTop; sim_auto

theorem sim_loop {cf : WalkM Bool} {c : Spqr.WalkM Bool} {bf : WalkM Unit} {b : Spqr.WalkM Unit}
    (hc : Refine.Sim WalkState.toSlow Eq cf c) (hb : Refine.Sim WalkState.toSlow Eq bf b) :
    ∀ n, Refine.Sim WalkState.toSlow Eq (loop n cf bf) (Spqr.WalkM.loop n c b)
  | 0 => Refine.Sim.pure _ _ rfl
  | n + 1 => by
    have ih := sim_loop hc hb n
    unfold loop Spqr.WalkM.loop
    sim_auto

end WalkM

macro_rules
  | `(tactic| sim_leaf) => `(tactic| first
      | apply WalkM.sim_loop | apply WalkM.sim_pushVertTstack | apply WalkM.sim_pushEdgeTstack
      | apply WalkM.sim_mergeTstackTops | apply WalkM.sim_maybeUnwrapNxt
      | apply WalkM.sim_finishTstackTop)

namespace WalkM

set_option maxHeartbeats 4000000 in
set_option sim_auto.trace true in
theorem sim_finishEdge (curV d : Nat) (o : DfsOut) (origTstack : Nat) (hasVert : Bool) :
    Refine.Sim WalkState.toSlow Eq (finishEdge curV d o origTstack hasVert)
      (Spqr.finishEdge curV d o origTstack hasVert) := by
  unfold finishEdge Spqr.finishEdge
  sim_auto

end WalkM

macro_rules
  | `(tactic| sim_leaf) => `(tactic| apply WalkM.sim_finishEdge)

mutual
theorem sim_walkTree (t : DfsTree) (d : Nat) :
    Refine.Sim WalkState.toSlow Eq (walkTree t d) (Spqr.walkTree t d) := by
  match t with
  | .node v outs =>
    have ih := sim_walkOuts v d outs false
    unfold walkTree Spqr.walkTree
    sim_auto
theorem sim_walkOuts (v d : Nat) (outs : List DfsOut) (hasVert : Bool) :
    Refine.Sim WalkState.toSlow Eq (walkOuts v d outs hasVert) (Spqr.walkOuts v d outs hasVert) := by
  match outs with
  | [] => unfold walkOuts Spqr.walkOuts; sim_auto
  | o :: rest =>
    have ih1 := sim_walkOut v d o hasVert
    have ih2 := fun hv => sim_walkOuts v d rest hv
    unfold walkOuts Spqr.walkOuts
    sim_auto
theorem sim_walkOut (v d : Nat) (o : DfsOut) (hasVert : Bool) :
    Refine.Sim WalkState.toSlow Eq (walkOut v d o hasVert) (Spqr.walkOut v d o hasVert) := by
  match o with
  | .tree a b child =>
    have ih := sim_walkTree child (d + 1)
    unfold walkOut Spqr.walkOut
    sim_auto
  | .back .. =>
    unfold walkOut Spqr.walkOut
    sim_auto
end

theorem sim_walkForest (forest : List DfsTree) :
    Refine.Sim WalkState.toSlow Eq (walkForest forest) (Spqr.walkForest forest) := by
  have := fun t d => sim_walkTree t d
  unfold walkForest Spqr.walkForest
  sim_auto

theorem WalkState.init_toSlow (g : Graph) (tern : Bool) :
    (WalkState.init g tern).toSlow = Spqr.WalkState.init g tern := by
  simp [WalkState.toSlow, WalkState.init, Spqr.WalkState.init, initialItems_toSlow, tstackList]

end Fast

theorem Graph.walkFast_toSlow (g : Graph) (tern : Bool) (forest : List DfsTree) :
    (g.walkFast tern forest).toSlow = g.walk tern forest := by
  have h := (Fast.sim_walkForest forest (Fast.WalkState.init g tern)).1
  rw [Fast.WalkState.init_toSlow] at h
  exact h

/-- The fast walk produces the same items as `Graph.walk`, up to flattening the children spans. -/
theorem Graph.walkFast_items (g : Graph) (tern : Bool) (forest : List DfsTree) :
    (g.walkFast tern forest).items.map Fast.Item.toSlow = (g.walk tern forest).items := by
  rw [← Graph.walkFast_toSlow]; rfl

end Spqr
