import WalkInvCheck.R
/-! Executable mirrors of the `RCloseContent`/`REntryContent` fields (`Spqr/Proofs/RSiteContent.lean`),
one reported field each (`rc.<field>`). -/
open Spqr WalkM

namespace WalkInvCheck.RC
open WalkInvCheck.R

/-- The HT content of `RCloseShape` (`single`, `maximal`, `bond`, `type1`, `type2`) for the pieces
`Pi` with union `Et` closed at `(a, b)`; `entryR` for a union of entries. -/
def unionR (D : Dfs) (s : WalkState) (Pi : List ItemId) (Et : List Nat) (a b : Nat) :
    List (String × String) := Id.run do
  let g := s.g
  let mut bad := []
  let P := Pi.map fun i => (i, belowList s i)
  if !twoAttached g Et.contains a b && !oneClass D a b Et then bad := ("single", "") :: bad
  for (i, Ei) in P do
    match vsOf s i with
    | none => bad := ("maximal", s!"vs i={i}") :: bad
    | some (x, y) =>
      if !oneClass D x y ((List.range g.ne).filter fun e => !Ei.contains e) then
        bad := ("maximal", s!"i={i}") :: bad
  for e₁ in Et do
    for e₂ in Et do
      if e₁ != e₂ && (g.edges[e₁]! == g.edges[e₂]! || g.edges[e₁]! == (g.edges[e₂]!.2, g.edges[e₂]!.1)) then
        if !P.any (fun (_, Ei) => Ei.contains e₁ && Ei.contains e₂) then
          bad := ("bond", s!"e={e₁},{e₂}") :: bad
  let sp (a' b' : Nat) : Bool :=
    P.all (fun (i, Ei) => match vsOf s i with
      | some (x, y) => (!(Ei.any fun e => inc g e b') || b' == x || b' == y) &&
          !((a' == x && b' == y) || (a' == y && b' == x))
      | none => true) &&
    !((a' == a && b' == b) || (a' == b && b' == a))
  for (a', b', K) in D.t1 do
    if sp a' b' && !laminar s P Et K then bad := ("type1", s!"a={a'} b={b'}") :: bad
  for (a', b', K) in D.t2 do
    if sp a' b' && !laminar s P Et K then bad := ("type2", s!"a={a'} b={b'}") :: bad
  return bad

/-- `vert_free`: `VertFree v s` (no open entry owns an edge below `vertItem v`). -/
def vertFree (v : Nat) (s : WalkState) : List (String × String) :=
  let vE := belowList s (vertItem v)
  s.tstack.flatMap fun t =>
    if (t.spans.1 ++ t.spans.2).any fun i => (belowList s i).any vE.contains then
      [("vert_free", s!"v={v} {showT t}")] else []

/-- Loop 1 iterate by iterate from `ceS₁`; at every `.R` iterate (state `rl1Iter d o s k`, stack
`cur :: nxt :: _`) the fields `l1_mid` (`RBranch.mid`), `l1_fields` (`RBranchFields`) and `l1_top`
(`RTop`: `EntryR cur`, `EntryR nxt`, disjoint). -/
partial def l1Emu (D : Dfs) (d : Nat) (edgeDir : Bool) (s : WalkState) (fuel : Nat) :
    List (String × String) :=
  if fuel = 0 then [] else
  if WalkState.result (loop1Cond d) s then
    let bad := if WalkState.l1Ty d edgeDir s == .R then match s.tstack with
      | cur :: nxt :: _ =>
        ((rBranch d s cur nxt).filter fun b => !b.startsWith "interior").map (fun b =>
          (if b.startsWith "mid" then "l1_mid" else "l1_fields", b)) ++
        (entryR D s cur).map (("l1_top", "cur " ++ ·)) ++
        (entryR D s nxt).map (("l1_top", "nxt " ++ ·)) ++
        (let Ec := (cur.spans.1 ++ cur.spans.2).flatMap (belowList s)
         let En := (nxt.spans.1 ++ nxt.spans.2).flatMap (belowList s)
         if Ec.any En.contains then [("l1_top", s!"disj {showT cur} {showT nxt}")] else [])
      | _ => [("l1_mid", "short stack")]
      else []
    bad ++ l1Emu D d edgeDir ((loop1Body d edgeDir).run s).2 (fuel - 1)
  else []

/-- The `RCloseContent` fields at the pre-state `s` of `finishEdge v d o orig hv` for a returning
tree edge (`lowval < d`): `settled` (the `feS₁` frontier entries above `orig` are `EntryR` unless
exempt), `s2_top` (`hv = false`: the `feS₂` top with `topDepth = d`, `vStart ≠ v` is `EntryR`),
`vert_own` (`hv = false`: `VertFree v` after the P-check), the loop-1 fields (`l1Emu`), and the
type-1 vertex close with `feSingle = false` (`v_entry`: `EntryR` of `c`, `py` at `feS₂`;
`v_single`/`v_maximal`/`v_bond`/`v_type1`/`v_type2`: the HT content of `vU c py vy` at
`(v, stackVerts[lv])`). -/
def closeContent (D : Dfs) (v d : Nat) (o : DfsOut) (orig : Nat) (hv : Bool) (s : WalkState) :
    List (String × String) := Id.run do
  let mut bad := []
  let s₁ := WalkState.feS₁ d o s
  for t in s₁.tstack.tail.take (s₁.tstack.length - 1 - orig) do
    for b in settledEntry D v d s₁ t "" do bad := ("settled", b) :: bad
  let s₂ := WalkState.feS₂ d o s
  if !hv then
    match s₂.tstack with
    | c :: _ =>
      if c.topDepth == d && c.vStart != v then
        for b in entryR D s₂ c do bad := ("s2_top", b) :: bad
    | [] => pure ()
  if !hv then
    let sP := ((Spqr.finishP v (o.cls.lowval d) o.cls.isType1).run s₂).2
    for (_, b) in vertFree v sP do bad := ("vert_own", b) :: bad
  let sp := WalkState.ceS₁ o.dest d o.e (WalkState.feS₀ d o s)
  bad := l1Emu D d s.stackDir[d]! sp sp.tstack.length ++ bad
  if hv && o.cls.isType1 && !WalkState.feSingle d o s then
    match s₂.tstack with
    | c :: py :: vy :: _ =>
      for b in entryR D s₂ c do bad := ("v_entry", s!"c: {b}") :: bad
      for b in entryR D s₂ py do bad := ("v_entry", s!"py: {b}") :: bad
      let vItems := ((c.spans.1 ++ py.spans.1) ++ vy.spans.1) ++ (vy.spans.2 ++ (py.spans.2 ++ c.spans.2))
      let Pi := vItems.filter fun i => Items.type s₂.items i != .V
      let Et := vItems.flatMap (belowList s₂)
      let lv := o.cls.lowval d
      for (f, b) in unionR D s₂ Pi Et v s₂.stackVerts[lv]! do
        bad := (s!"v_{f}", s!"lv={lv} c={showT c} py={showT py} vy={showT vy} {b}") :: bad
    | _ => bad := ("v_entry", "short stack at feS₂") :: bad
  return bad.reverse

/-- The `REntryContent` fields at the entry of `walkTree c dep` (`s` before
`stackVerts.set! dep c`, `s'` after; parent `p = stackVerts[dep - 1]`): `entry_parent_top` (no
entry starts at `p` with `topDepth > dep - 1`), `entry_stab` (an `EntryR` entry with
`topDepth = dep` stays `EntryR` under the write). -/
def entryContent (D : Dfs) (dep c : Nat) (s s' : WalkState) : List (String × String) :=
  if dep = 0 then [] else
  let p := s.stackVerts[dep - 1]!
  s.tstack.flatMap fun t =>
    (if t.vStart == p && t.topDepth > dep - 1 then [("entry_parent_top", s!"p={p} {showT t}")] else []) ++
    (if t.topDepth == dep && (entryR D s t).isEmpty && !(entryR D s' t).isEmpty then
      [("entry_stab", s!"c={c} {showT t} loses {entryR D s' t}")] else [])

end WalkInvCheck.RC
