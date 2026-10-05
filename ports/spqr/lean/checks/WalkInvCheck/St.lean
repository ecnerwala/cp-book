import WalkInvCheck.Common
/-! Verbatim copy of the executable mirrors of `CheckStSim.lean` (`StRead`, `StItems` m1–m7, the
truncated-reference fields m8/m10/m11, `StLiveCtx` m16/m17 (`truncCtxB`) and m13/m15 (`truncSiteB`);
m9/m12/m18 stay off: they are false, see PROOF.md §7). Keep in sync with that file. -/
namespace WalkInvCheck.St
open Spqr WalkM

/-! Differential test of the simulation relation `StSim` / `StSimOuts` (`Spqr/StSim.lean`): the walk
is re-run step by step (mirroring `walkTree` / `walkOuts` / `walkOut`, calling the real
`finishEdge`) and the relation is checked at the start of every out-edge, before every `finishEdge`
and at the end of every `walkTree`. -/

instance (p q : Nat × Nat) : Decidable (Items.PairEq p q) := by unfold Items.PairEq; infer_instance

structure Stats where
  checks : Nat := 0
  bad : Nat := 0

instance : BEq TEntry :=
  ⟨fun a b => a.vStart == b.vStart && a.topDepth == b.topDepth && a.firstIdx == b.firstIdx && a.spans == b.spans⟩

def isSPR (items : Items) (i : ItemId) : Bool :=
  match Items.type items i with
  | .S | .P | .R => true
  | _ => false

def isInfixB (l m : List Nat) : Bool :=
  (List.range (m.length + 1)).any fun k => (m.drop k).take l.length == l

def descOf (items : Items) : Nat → ItemId → List ItemId
  | 0, i => [i]
  | fuel + 1, i => i :: (Items.ch items i).flatMap (descOf items fuel)

def leavesB (items : Items) (i : ItemId) : List ItemId := Items.leaves items items.size i

def stReadB (items : Items) (new : List TEntry) (ps : List StPiece) : Bool :=
  ((readL new).flatMap (leavesB items) == stNestL ps && (readR new).flatMap (leavesB items) == stNestR ps)

def vsOrientedAtB (g : Graph) (items : Items) (b : StBlock) (i : ItemId) : Bool :=
  let lv := leavesB items i
  let ds := descOf items items.size i
  lv.all (· ∈ b.items) &&
  ((List.range g.ne).all fun e =>
    let x := edgeItem g e
    !(x ∈ b.items || b.root.any fun r => decide (Items.PairEq g.edges[e]! r)) || !(x ∈ ds) || x ∈ lv) &&
  decide (Oriented (b.seq g) (Items.vs items i)) &&
  (match Items.vs items i with
    | (some s, some t) => (Items.ch items i).all fun c => Items.type items c ≠ .V ||
        (decide (Precedes (b.seq g) s (c - 1)) && decide (Precedes (b.seq g) (c - 1) t))
    | _ => false) &&
  ((Items.ch items i).all fun c => Items.type items c = .V || decide (Oriented (b.seq g) (Items.vs items c))) &&
  ((Items.ch items i).all fun c => Items.type items c = .V ||
    match Items.vs items c with
    | (some u, some v) =>
      let dc := descOf items items.size c
      (List.range g.ne).all fun e =>
        !(edgeItem g e ∈ b.items || b.root.any fun r => decide (Items.PairEq g.edges[e]! r)) ||
        !(edgeItem g e ∈ dc) ||
        [(g.edges[e]!).1, (g.edges[e]!).2].all fun y =>
          (u == y || decide (Precedes (b.seq g) u y)) && (y == v || decide (Precedes (b.seq g) y v))
    | _ => false)

def inBlockB (g : Graph) (items : Items) (b : StBlock) (i : ItemId) : Bool :=
  isInfixB (leavesB items i) b.items && vsOrientedAtB g items b i

/-- Violations of `StItems`. -/
def stItemsB (g : Graph) (s : WalkState) (blocks : List StBlock) : List String :=
  let rs := readStack s.tstack
  let items := s.items
  let m1 := if rs.all (fun x => (List.range items.size).all fun p => !(x ∈ Items.ch items p)) then []
    else ["roots: a stack item has a parent"]
  let m2 := if decide rs.Nodup then [] else [s!"nodup: {rs}"]
  let m4 := if rs.all (fun x => (descOf items items.size x).all (· < items.size)) then []
    else ["bounded: a stack item's subtree is out of range"]
  let m3 := (List.range items.size).filterMap fun i =>
    if !isSPR items i then none
    else if rs.any (fun x => i ∈ descOf items items.size x) then none
    else if blocks.any (inBlockB g items · i) then none
    else some s!"item {i} ({repr (Items.type items i)}) vs {repr (Items.vs items i)} ch {Items.ch items i} leaves {leavesB items i} not live and in no block {blocks.map (·.items)}"
  let m5 := if (List.range items.size).all (fun p => (Items.ch items p).all (· < items.size)) then []
    else ["chLt: a child is out of range"]
  let m6 := if (List.range items.size).all (fun p => decide (Items.ch items p).Nodup) then []
    else ["chNodup: a child list repeats"]
  let m7 := (List.range items.size).filterMap fun x =>
    if !(Items.type items x == .V || Items.type items x == .Q) then none
    else match (descOf items items.size x).find? fun i =>
        isSPR items i && !blocks.any (inBlockB g items · i) with
      | none => none
      | some i => some s!"finished: item {i} below the leaf-type item {x} is in no block"
  m1 ++ m2 ++ m4 ++ m5 ++ m6 ++ m3 ++ m7

def report (st : IO.Ref Stats) (ok : Bool) (msg : String) : IO Unit := do
  st.modify fun x => { x with checks := x.checks + 1, bad := if ok then x.bad else x.bad + 1 }
  unless ok do IO.println s!"BAD {msg}"

def above (s : WalkState) (base : List TEntry) : List TEntry := s.tstack.take (s.tstack.length - base.length)
/-- The per-entry terminal checks m9 / m12 (`(stackVerts[t.topDepth], t.vStart)` oriented, V span items
strictly between) are FALSE as stated (seeds 6, 9, 12, 17, 554: buried entries of finished sibling
subtrees keep a `topDepth` whose `stackVerts` is stale; two-sided merged entries and the fold retarget
`vStart := curV`); kept here as the record of what was tried, off by default. -/
def truncEntryChecks : Bool := false

/-- Candidate invariant (sim-4): against the blocks of the *truncated* forest `prev ++ [truncTree fs t]`
(the open blocks along the path read with the pieces walked so far), every S / P / R item is
`InBlock` (m8), and every live non-V item with two endpoints is oriented in its block's st-order
(m10). -/
def truncItemsB (g : Graph) (s : WalkState) (tblocks : List StBlock) : List String :=
  let items := s.items
  let live := (readStack s.tstack).flatMap (descOf items items.size)
  let m8 := (List.range items.size).filterMap fun i =>
    if !isSPR items i then none
    else if tblocks.any (inBlockB g items · i) then none
    else some s!"trunc: S/P/R item {i} is in no truncated block"
  let m10 := live.filterMap fun y =>
    if Items.type items y = .V then none
    else match Items.vs items y with
      | (some a, some b) =>
        if tblocks.any (fun b' => decide (Oriented (b'.seq g) (some a, some b))) then none
        else some s!"trunc: live item {y} vs ({a}, {b}) oriented in no truncated block"
      | _ => none
  let m9 := s.tstack.filterMap fun t =>
    let a := s.stackVerts[t.topDepth]!
    let dir := s.stackDir[t.topDepth]!
    let σ := if getSide t.spans dir != [] then dir else !dir
    if a == t.vStart || vertItem t.vStart ∈ t.spans.1 ++ t.spans.2 then none
    else if tblocks.any (fun b' => decide (Oriented (b'.seq g) (setSides σ (some a) (some t.vStart)))) then none
    else some s!"trunc: entry (vStart {t.vStart}, topDepth {t.topDepth}, sv {a}, dir {σ}) spans {t.spans} unoriented"
  let m11 := live.flatMap fun y =>
    if !isSPR items y then []
    else match Items.vs items y with
      | (some u, some v) =>
        let dy := descOf items items.size y
        let b? := tblocks.find? (inBlockB g items · y)
        match b? with
        | none => []
        | some b' =>
          (List.range g.ne).filterMap fun e =>
            if !(edgeItem g e ∈ dy) then none
            else if !(edgeItem g e ∈ b'.items || b'.root.any fun r => decide (Items.PairEq g.edges[e]! r)) then none
            else if [(g.edges[e]!).1, (g.edges[e]!).2].all (fun z =>
                (u == z || decide (Precedes (b'.seq g) u z)) && (z == v || decide (Precedes (b'.seq g) z v))) then none
            else some s!"trunc: live S/P/R item {y} vs ({u}, {v}) edge {e} {g.edges[e]!} outside; seq {b'.seq g} leaves {leavesB items y} stack {s.tstack.map fun t => (t.vStart, t.topDepth, t.spans)}"
      | _ => []
  let m12 := s.tstack.flatMap fun t =>
    let a := s.stackVerts[t.topDepth]!
    let dir := s.stackDir[t.topDepth]!
    let σ := if getSide t.spans dir != [] then dir else !dir
    let (a', b') := setSides σ a t.vStart
    if a == t.vStart || vertItem t.vStart ∈ t.spans.1 ++ t.spans.2 || (t.spans.1 != [] && t.spans.2 != []) then [] else
    (t.spans.1 ++ t.spans.2).filterMap fun c =>
      if Items.type items c != .V then none
      else if tblocks.any (fun b'' => decide (Precedes (b''.seq g) a' (c - 1)) && decide (Precedes (b''.seq g) (c - 1) b')) then none
      else some s!"trunc: V item {c} under entry (vStart {t.vStart}, topDepth {t.topDepth}, sv {a}, dir {σ}) (span item) not strictly between; seqs {tblocks.map (·.seq g)} stack {s.tstack.map fun t => (t.vStart, t.topDepth, t.spans)} sv {s.stackVerts.toList.take 6}"
  m8 ++ m10 ++ m11 ++ (if truncEntryChecks then m9 ++ m12 else [])

def betweenAny (xs : List Nat) (x y z : Nat) : Bool :=
  (decide (Precedes xs x z) && decide (Precedes xs z y)) ||
  (decide (Precedes xs y z) && decide (Precedes xs z x))

/-- `some false`: `(a, b)` oriented in a block; `some true`: `(b, a)`; `none`: neither. -/
def orientedDir (tblocks : List StBlock) (g : Graph) (a b : Nat) : Option Bool :=
  if tblocks.any (fun b' => decide (Oriented (b'.seq g) (some a, some b))) then some false
  else if tblocks.any (fun b' => decide (Oriented (b'.seq g) (some b, some a))) then some true
  else none

def vSpans (items : Items) (t : TEntry) : List ItemId :=
  (t.spans.1 ++ t.spans.2).filter (Items.type items · == .V)

def showStack (ts : List TEntry) : String := s!"{ts.map fun t => (t.vStart, t.topDepth, t.spans)}"

/-- Entry-level candidate (sim-5), only for the entries `(v, l)`, `l < d`, started at the current
vertex `v` (`EarCtx.above`, whose `stackVerts[l]` is current): `(stackVerts[l], v)` is oriented by
`stackDir[l]` (m16) and the V span items are strictly between (m17). -/
def truncCtxB (g : Graph) (s : WalkState) (tblocks : List StBlock) (v d : Nat) : List String :=
  s.tstack.flatMap fun t =>
    if t.vStart != v || d ≤ t.topDepth then [] else
    let a := s.stackVerts[t.topDepth]!
    let dir := s.stackDir[t.topDepth]!
    let m16 := match orientedDir tblocks g a v with
      | none => [s!"ctx: entry (v {v}, l {t.topDepth}, sv {a}) unoriented; stack {showStack s.tstack}"]
      | some flip => if flip != dir then
          [s!"ctx: entry (v {v}, l {t.topDepth}, sv {a}) oriented against stackDir {dir}; stack {showStack s.tstack}"]
        else []
    let m17 := (vSpans s.items t).filterMap fun c =>
      if tblocks.any (fun b' => betweenAny (b'.seq g) a v (c - 1)) then none
      else some s!"ctx: V item {c} of entry (v {v}, l {t.topDepth}, sv {a}) not strictly between; stack {showStack s.tstack}"
    m16 ++ m17

/-- Site-level candidate (sim-5), before `finishEdge o` of `v` at depth `d` with the child's leftover
`sub`: for a returning edge the ear's terminals `(stackVerts[lowval], v)` are oriented by
`stackDir[lowval]` and every V span item of `sub` is strictly between (m13/m14); for every loop-1
entry `t'` of depth `d` the close `(v, t'.vStart)` is oriented by `stackDir[d]` and the V span items
of the range down to `t'` are strictly between (m15). (m18, `topDepth ≤ d + 1` for loop-1 entries,
is FALSE: seeds 154, 195, 442 — buried chain entries keep a deeper `topDepth`; off.) -/
def truncSiteB (g : Graph) (s : WalkState) (tblocks : List StBlock) (v d : Nat) (o : DfsOut)
    (sub : List TEntry) : List String :=
  let lowval := o.cls.lowval d
  let m13 := if d ≤ lowval then [] else
    let a := s.stackVerts[lowval]!
    let dir := s.stackDir[lowval]!
    (match orientedDir tblocks g a v with
      | none => [s!"ret: ({a}, {v}) lowval {lowval} unoriented; sub {showStack sub}"]
      | some flip => if flip != dir then
          [s!"ret: ({a}, {v}) lowval {lowval} oriented against stackDir {dir}; sub {showStack sub}"] else []) ++
    ((sub.flatMap (vSpans s.items)).filter (· != vertItem v)).filterMap fun c =>
      if tblocks.any (fun b' => betweenAny (b'.seq g) a v (c - 1)) then none
      else some s!"ret: V item {c} not strictly between ({a}, {v}) lowval {lowval}; sub {showStack sub}"
  let hi := sub.takeWhile (fun t => d ≤ t.topDepth)
  let m15 := (List.range hi.length).flatMap fun j =>
    let t' := hi[j]!
    if t'.topDepth != d then [] else
    let bottom := t'.vStart
    let dir := s.stackDir[d]!
    (match orientedDir tblocks g v bottom with
      | none => [s!"l1: ({v}, {bottom}) unoriented; hi {showStack hi}"]
      | some flip => if flip != dir then
          [s!"l1: ({v}, {bottom}) oriented against stackDir {dir}; hi {showStack hi}"] else []) ++
    (((hi.take (j + 1)).flatMap (vSpans s.items)).filter (· != vertItem v)).filterMap fun c =>
      if tblocks.any (fun b' => betweenAny (b'.seq g) v bottom (c - 1)) then none
      else some s!"l1: V item {c} not strictly between ({v}, {bottom}); hi {showStack hi}"
  let m18 := hi.filterMap fun t =>
    if d + 1 < t.topDepth then some s!"l1: entry topDepth {t.topDepth} > d+1; hi {showStack hi}" else none
  m13 ++ m15 ++ (if truncEntryChecks then m18 else [])

def baseOk (s : WalkState) (base : List TEntry) : Bool :=
  base.length ≤ s.tstack.length && s.tstack.drop (s.tstack.length - base.length) == base

end WalkInvCheck.St
