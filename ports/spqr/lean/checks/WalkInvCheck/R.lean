import WalkInvCheck.Common
/-! Verbatim copy of the executable mirrors of `checks/RFinishEdgeCheck.lean` (`EntryR` = `RInvTop`/
`RInvFront` contract B, `RBranch`, loop-1 emulation, `FinishRShape`, `RInvG` base frame, `vertOwn`).
Keep in sync with that file. -/
open Spqr WalkM

namespace WalkInvCheck.R

instance : Inhabited DfsOut := ⟨.back 0 0 .selfLoop⟩

def inc (g : Graph) (e v : Nat) : Bool := g.edges[e]!.1 == v || g.edges[e]!.2 == v

/-- Component labels of the vertices satisfying `ok` in the subgraph of the edges `E` with both
endpoints `ok` (min-label propagation, `nv + 1` rounds). -/
def comps (g : Graph) (E : Nat → Bool) (ok : Nat → Bool) : Array Nat := Id.run do
  let mut lab := (List.range g.nv).toArray
  let es := (List.range g.ne).filter fun e => E e && ok g.edges[e]!.1 && ok g.edges[e]!.2
  for _ in List.range (g.nv + 1) do
    for e in es do
      let (x, y) := g.edges[e]!
      let m := min lab[x]! lab[y]!
      lab := (lab.set! x m).set! y m
  return lab

/-- The DFS data of a block, with the separation-class tables precomputed. -/
structure Dfs where
  g : Graph
  parent : Array (Option Nat)
  depth : Array Nat
  outs : Nat → List DfsOut
  /-- `sepLab[a * nv + b]`: component labels of `g − {a, b}`. -/
  sepLab : Array (Array Nat)
  /-- `(a, b, K)` for each type-1 child edge `o` of `b` to `depth a`: `K = EndIn o.dest`. -/
  t1 : List (Nat × Nat × Array Bool)
  /-- `(a, b, K)` for each `Type2Pair a b` and tree edge `o` of `a` towards `b`:
  `K = SepClass a b o.e`. -/
  t2 : List (Nat × Nat × Array Bool)

partial def fillTree (p : Option Nat) (d : Nat) : DfsTree → Array (Option Nat) × Array Nat →
    Array (Option Nat) × Array Nat
  | .node v outs, (par, dep) =>
    outs.foldl (fun acc o => match o with
      | .tree _ _ c => fillTree (some v) (d + 1) c acc
      | .back .. => acc) (par.set! v p, dep.set! v d)

partial def ancP (parent : Array (Option Nat)) (a b : Nat) : Bool :=
  if a == b then true else match parent[b]! with
    | some p => ancP parent a p
    | none => false

/-- The edge class of `e` for `SepClass a b` under the labels `lab`: the label of an endpoint
outside `{a, b}` (`none`: `e` is a class of its own). -/
def rep (g : Graph) (a b : Nat) (lab : Array Nat) (e : Nat) : Option Nat :=
  let (x, y) := g.edges[e]!
  if x != a && x != b then some lab[x]! else if y != a && y != b then some lab[y]! else none

def sepClass (g : Graph) (a b : Nat) (lab : Array Nat) (e e' : Nat) : Bool :=
  e == e' || match rep g a b lab e, rep g a b lab e' with
    | some r, some r' => r == r'
    | _, _ => false

def mkDfs (g : Graph) (forest : List DfsTree) : Dfs := Id.run do
  let (par, dep) := forest.foldl (fun acc t => fillTree none 0 t acc)
    (Array.replicate g.nv none, Array.replicate g.nv 0)
  let dd := DfsData.ofForest forest
  let anc := ancP par
  let backD := (List.range g.nv).toArray.map fun v => (dd.outs v).filterMap fun o =>
    if o.isTree then none else some dep[o.dest]!
  let returns := fun c l => (List.range g.nv).any fun u => anc c u && backD[u]!.contains l
  let noBothSides := fun a b => ((dd.outs b).filter (·.isTree)).all fun o =>
    !((List.range dep[a]!).any (returns o.dest)) ||
      (List.range dep[b]!).all fun l => !(dep[a]! < l) || !(returns o.dest l)
  let betweenStays := fun a a' b => (List.range g.nv).all fun u => !(anc a' u) || anc b u ||
    (dd.outs u).all fun o => o.isTree || dep[a]! ≤ dep[o.dest]!
  let sepLab := (List.range (g.nv * g.nv)).toArray.map fun k =>
    let a := k / g.nv; let b := k % g.nv
    comps g (fun _ => true) fun v => v != a && v != b
  let mut t1 := []
  for b in List.range g.nv do
    for o in dd.outs b do
      match o.cls with
      | .ret l .type1Child =>
        for a in List.range g.nv do
          if dep[a]! == l && anc a b then
            t1 := (a, b, (List.range g.ne).toArray.map fun e =>
              [g.edges[e]!.1, g.edges[e]!.2].any fun x => anc o.dest x) :: t1
      | _ => pure ()
  let mut t2 := []
  for a in List.range g.nv do
    if par[a]!.isSome then
      for o in dd.outs a do
        if o.isTree then
          for b in List.range g.nv do
            if anc o.dest b && o.dest != b && noBothSides a b && betweenStays a o.dest b then
              let lab := sepLab[a * g.nv + b]!
              t2 := (a, b, (List.range g.ne).toArray.map fun e => sepClass g a b lab o.e e) :: t2
  return ⟨g, par, dep, dd.outs, sepLab, t1, t2⟩

/-- `∀ e e' ∈ C, g.SepClass a b e e'`. -/
def oneClass (D : Dfs) (a b : Nat) (C : List Nat) : Bool :=
  match C with
  | [] => true
  | e :: rest =>
    let lab := D.sepLab[a * D.g.nv + b]!
    rest.all fun e' => sepClass D.g a b lab e e'

/-- `g.ConnEdges E` for a nonempty `E`. -/
def connEdges (g : Graph) (E : List Nat) : Bool :=
  match E with
  | [] => true
  | e₀ :: _ =>
    let lab := comps g E.contains fun _ => true
    E.all fun e => lab[g.edges[e]!.1]! == lab[g.edges[e₀]!.1]!

def twoAttached (g : Graph) (E : Nat → Bool) (a b : Nat) : Bool :=
  (List.range g.nv).all fun v => v == a || v == b ||
    !((List.range g.ne).any fun e => E e && inc g e v) ||
    !((List.range g.ne).any fun e => !E e && inc g e v)

/-- The edges `e < ne` below item `i`. -/
partial def belowList (s : WalkState) (i : ItemId) : List Nat :=
  (if 1 + s.g.nv ≤ i ∧ i < 1 + s.g.nv + s.g.ne then [i - 1 - s.g.nv] else []) ++
    (Items.ch s.items i).flatMap (belowList s)

def pieceItems (s : WalkState) (t : TEntry) : List ItemId :=
  (t.spans.1 ++ t.spans.2).filter fun i => Items.type s.items i != .V

def vsOf (s : WalkState) (i : ItemId) : Option (Nat × Nat) :=
  match Items.vs s.items i with
  | (some x, some y) => some (x, y)
  | _ => none

/-- `s.EntrySkelPair t a b`. -/
def skelPair (s : WalkState) (t : TEntry) (P : List (ItemId × List Nat)) (a b : Nat) : Bool :=
  P.all (fun (i, Ei) => match vsOf s i with
    | some (x, y) => (!(Ei.any fun e => inc s.g e b) || b == x || b == y) &&
        !((a == x && b == y) || (a == y && b == x))
    | none => true) &&
  !((a == t.vStart && b == s.stackVerts[t.topDepth]!) ||
    (a == s.stackVerts[t.topDepth]! && b == t.vStart))

/-- `s.EntryLaminar t K` (`K` over `e < ne`). -/
def laminar (s : WalkState) (P : List (ItemId × List Nat)) (Et : List Nat) (K : Array Bool) : Bool :=
  let Ks := (List.range s.g.ne).filter (K[·]!)
  P.any (fun (_, Ei) => Ks.all Ei.contains) || Ks.all (fun e => !Et.contains e) || Et.all (K[·]!)

def showT (t : TEntry) : String := s!"({t.vStart},{t.topDepth},{t.firstIdx},{t.spans})"

/-- The failing fields of `s.EntryR dfs t`. -/
def entryR (D : Dfs) (s : WalkState) (t : TEntry) : List String := Id.run do
  let g := s.g
  let mut bad := []
  let Pi := pieceItems s t
  let P := Pi.map fun i => (i, belowList s i)
  let Et := (t.spans.1 ++ t.spans.2).flatMap (belowList s)
  let a := t.vStart
  let b := s.stackVerts[t.topDepth]!
  if !Pi.Nodup then bad := "pieces.nodup" :: bad
  for (i, Ei) in P do
    match vsOf s i with
    | none => bad := s!"pieces.vs i={i}" :: bad
    | some (x, y) =>
      if Ei.isEmpty then bad := s!"pieces.ne i={i}" :: bad
      if !connEdges g Ei then bad := s!"pieces.conn i={i}" :: bad
      if !twoAttached g Ei.contains x y then bad := s!"pieces.attached i={i}" :: bad
      if !oneClass D x y ((List.range g.ne).filter fun e => !Ei.contains e) then
        bad := s!"maximal i={i}" :: bad
    for (j, Ej) in P do
      if i != j && Ei.any Ej.contains then bad := s!"pieces.disj i={i} j={j}" :: bad
  for e₁ in Et do
    for e₂ in Et do
      if e₁ != e₂ && (g.edges[e₁]! == g.edges[e₂]! || g.edges[e₁]! == (g.edges[e₂]!.2, g.edges[e₂]!.1)) then
        if !P.any (fun (_, Ei) => Ei.contains e₁ && Ei.contains e₂) then
          bad := s!"bond e={e₁},{e₂}" :: bad
  if !twoAttached g Et.contains a b && !oneClass D a b Et then bad := "single" :: bad
  for (a', b', K) in D.t1 do
    if skelPair s t P a' b' && !laminar s P Et K then bad := s!"type1 a={a'} b={b'}" :: bad
  for (a', b', K) in D.t2 do
    if skelPair s t P a' b' && !laminar s P Et K then bad := s!"type2 a={a'} b={b'}" :: bad
  return bad.map (· ++ s!" t={showT t}")

def disjoint (s : WalkState) : List String := Id.run do
  let mut bad := []
  let E := s.tstack.map fun t => (t.spans.1 ++ t.spans.2).flatMap (belowList s)
  for i in List.range s.tstack.length do
    for j in List.range i do
      if E[i]!.any E[j]!.contains then
        bad := s!"disj {showT s.tstack[i]!} {showT s.tstack[j]!}" :: bad
  return bad

/-- `EntryR` of the entries of `ts` kept by `keep`, tagged with the entry's position relative to `d`. -/
def checkEntries (D : Dfs) (d : Nat) (ts : List TEntry) (s : WalkState) (keep : TEntry → Bool) : List String :=
  ts.flatMap fun t => if keep t then (entryR D s t).map
    (· ++ s!" [td={t.topDepth} d={d} depV={D.depth[t.vStart]!}]") else []

def tagged (tag : String) (l : List String) : List String :=
  l.map fun b => if b.startsWith "stat:" then b else tag ++ b

/-- `RBranch`'s content fields at an R branch of loop 1 (`Spqr/RInv.lean`): `mid` (the child is the
bottom or interior to `cur`; the bottom-is-the-child form `cur_c` fails after an S merge,
`checks/RBranchCounter.lean`), `interior`, `cur_piece`, `cur_vs`, `nxt_ne`, `proper`,
`nxt_touch_top`, `nxt_touch_bot`, `nxt_no_cu`, read on the state after `loop1Type`. -/
def rBranch (d : Nat) (s : WalkState) (cur nxt : TEntry) : List String := Id.run do
  let g := s.g
  let mut bad := []
  let top := s.stackVerts[d]!
  let Ecur := (cur.spans.1 ++ cur.spans.2).flatMap (belowList s)
  let Enxt := (nxt.spans.1 ++ nxt.spans.2).flatMap (belowList s)
  let c := s.stackVerts[d + 1]!
  if cur.vStart != c && !(List.range g.ne).all (fun e => !inc g e c || Ecur.contains e) then
    bad := "mid" :: bad
  if !(List.range g.ne).all (fun e => !inc g e cur.vStart || Ecur.contains e || Enxt.contains e) then
    bad := "interior" :: bad
  let Pi := pieceItems s cur
  if Pi.isEmpty then bad := "cur_piece" :: bad
  for i in Pi do
    if vsOf s i != some (cur.vStart, top) && vsOf s i != some (top, cur.vStart) then
      bad := s!"cur_vs i={i}" :: bad
  if Enxt.isEmpty then bad := "nxt_ne" :: bad
  if (List.range g.ne).all (fun e => Ecur.contains e || Enxt.contains e) then bad := "proper" :: bad
  if !Enxt.any (inc g · top) then bad := "nxt_touch_top" :: bad
  if !Enxt.any (inc g · nxt.vStart) then bad := "nxt_touch_bot" :: bad
  if Enxt.any (fun e => inc g e cur.vStart && inc g e top && cur.vStart != top) then
    bad := "nxt_no_cu" :: bad
  return bad

/-- Loop 1 step by step from the state after `pushEdgeTstack`; at every R branch, `RTop` and the
`RBranch` content fields (`rBranch`) of `cur`, `nxt` (on the state after `loop1Type`). -/
partial def loop1Emu (D : Dfs) (d : Nat) (edgeDir : Bool) (s : WalkState) (fuel : Nat) :
    WalkState × List String :=
  if fuel = 0 then (s, []) else
  if WalkState.result (loop1Cond d) s then
    let (ty, s₁) := (loop1Type d edgeDir).run s
    let bad := if ty == .R then match s₁.tstack with
      | cur :: nxt :: _ => "stat:rclose" :: (entryR D s₁ cur).map ("rclose-cur " ++ ·) ++
          (entryR D s₁ nxt).map ("rclose-nxt " ++ ·) ++ (rBranch d s₁ cur nxt).map ("rbranch " ++ ·)
      | _ => ["rclose-short"]
      else []
    let (item, s₂) := (maybeUnwrapNxt ty).run s₁
    let s₃ := ((finishTstackTop item).run (mergeTstackTops.run s₂).2).2
    let (s', bad') := loop1Emu D d edgeDir s₃ (fuel - 1)
    (s', bad ++ bad')
  else (s, [])

/-- `Exempt v d t`, decidably. -/
def exempt (v d : Nat) (t : TEntry) : Bool := t.topDepth < d || t.vStart == v

/-- `d ≤ t.topDepth → t.vStart ≠ v → EntryR t`, as failure lines. -/
def settledEntry (D : Dfs) (v d : Nat) (s : WalkState) (t : TEntry) (tag : String) : List String :=
  if exempt v d t then [] else (entryR D s t).map fun b => s!"{tag} {showT t}: {b}"

/-- `FinishRShape dfs v d o orig hv s` (`Proofs/RInvTree.lean`) plus the admitted
`feS₂_top_entryR` (`s2top`: the `feS₂` top with `topDepth = d`, `vStart ≠ v` is `EntryR`, first-edge
case) and `finishEdge_tree_top_settled_first` (`ptop`: the top after the P-check), computed on the
library's `feS₁`, `feS₂`, `feP`. -/
def shapeCheck (D : Dfs) (v d : Nat) (o : DfsOut) (orig : Nat) (hv : Bool) (s : WalkState) :
    List String := Id.run do
  let mut bad := []
  for t in s.tstack do
    if (t.spans.1 ++ t.spans.2).any fun i => (belowList s i).contains o.e then
      bad := s!"pend owned by {showT t}" :: bad
  let s₁ := WalkState.feS₁ d o s
  for t in s₁.tstack.tail.take (s₁.tstack.length - 1 - orig) do
    bad := settledEntry D v d s₁ t "settled" ++ bad
  let s₂ := WalkState.feS₂ d o s
  if hv && o.cls.isType1 then
    match s₂.tstack with
    | _ :: b :: _ => if !exempt v d b then bad := s!"unwrap nonexempt nxt={showT b}" :: bad
    | _ => pure ()
  if !hv then
    match s₂.tstack with
    | c :: _ =>
      if c.topDepth == d && c.vStart != v then
        bad := "stat:s2top-nonexempt" :: (entryR D s₂ c).map (fun b => s!"s2top {showT c}: {b}") ++ bad
    | [] => pure ()
    let sP := WalkState.feP v d o s
    let vE := belowList sP (vertItem v)
    for t in sP.tstack do
      if (t.spans.1 ++ t.spans.2).any fun i => (belowList sP i).any vE.contains then
        bad := s!"vert_own {showT t}" :: bad
    match sP.tstack with
    | t :: _ =>
      bad := settledEntry D v d sP t "ptop" ++ bad
      if !exempt v d t then bad := "stat:ptop-nonexempt" :: bad
    | [] => pure ()
  return bad.reverse

def isBlock (g : Graph) : Bool := (List.range g.nv).all fun v => Id.run do
  let mut seen := if g.ne == 0 then [] else [0]
  for _ in List.range g.ne do
    seen := (List.range g.ne).filter fun e => seen.contains e || seen.any fun f =>
      [g.edges[e]!.1, g.edges[e]!.2].any fun u =>
        u != v && (g.edges[f]!.1 == u || g.edges[f]!.2 == u)
  return seen.length == g.ne

/-- No stack entry owns an edge below `vertItem v` (the vertex pushes of `walkOutPre` and the
`walkTree` tail). -/
def vertOwn (v : Nat) (s : WalkState) (tag : String) : List String := Id.run do
  let vE := belowList s (vertItem v)
  let mut bad := ["stat:vertown-site"]
  for t in s.tstack do
    if (t.spans.1 ++ t.spans.2).any fun i => (belowList s i).any vE.contains then
      bad := s!"{tag}: vert_own {showT t}" :: bad
  return bad

end WalkInvCheck.R
