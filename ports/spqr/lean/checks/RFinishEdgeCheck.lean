import Spqr.Proofs.RInvTree

/-!
# Empirical check of the provisional-frontier `finishEdge` R contract (`finishEdge_rInvTop`)

Run with `lake env lean --run checks/RFinishEdgeCheck.lean < cases` from `ports/spqr/lean` (input:
concatenated `gen.py` cases prefixed by their number; not part of the library; 6000 extra random
multigraphs are generated inside). Re-implements `walkTree` around the library's `finishEdge` (as
`checks/RInvReturnCheck.lean`) and, on every block, at every
`finishEdge curV d o origTstack hasVert` with `o.cls = .ret lv _`, `lv < d`, evaluates two
contracts, each as `pre` (entering `finishEdge`, i.e. the child-return state `RReturn` for tree
edges) and `post` (after `finishEdge`):
* contract A (old, refuted): every base entry (`pre`) / every entry (`post`) with `vStart ≠ curV`
  is `EntryR` (all six fields: `pieces`, `maximal`, `bond`, `single`, `type1`, `type2`);
* contract B (`RInvFront`/`RInvTop`, current): the same restricted to entries with
  `topDepth ≥ d`;
* `D`: the whole stack is edge-disjoint.
A-failures are reported as `expected-old-contract` and do not count: entries topping out above `d`
are provisional until the walk returns to their top (`checks/RFinishEdgeCounter.lean`,
`Deep`/`Base`). B/D failures, Loop-1 emulation mismatches, the R-branch `RTop` checks and `FinishRShape`
(`Proofs/RInvTree.lean`, the call-site hypotheses of the tree-edge branch: the pending tree edge
is unowned, after loop 1 every frontier entry below the top is settled, the type-1 `closeVert`
unwraps an exempt entry, the first-edge vertex entry takes unowned edges; plus the admitted
`finishEdge_tree_top_settled` (a): the top after the P-check is settled; `shape` lines) count.
* `base` lines (`finishEdge_rInvG_base`, `Proofs/RInvBase.lean`): for every ancestor frame
  `(p, dp, n₀)` of the site (`p = stackVerts[dp]`, `n₀` the stack size when `p`'s child subtree
  was entered), the bottom `n₀` entries are unchanged by `finishEdge`, their non-exempt entries
  (`dp ≤ topDepth`, `vStart ≠ p`) are `EntryR` before and after, and the positional side
  conditions hold: `n₀ ≤ origTstack`, `hasVert → n₀ + 1 ≤ origTstack`, and for a tree edge
  `origTstack + 3 ≤ length` after loop 2 (`hclose`).

`dfs` is `DfsData.ofForest forest`; `Anc` is read off the parent map, `Type2Pair` is evaluated
through `NoBothSides`/`BetweenStays`, `SepClass a b` through the components of `g − {a, b}`.
Everything "edge ∈ item" is restricted to `e < g.ne`.
-/

open Spqr WalkM

namespace RFinishEdgeCheck

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

/-- Loop 1 step by step from the state after `pushEdgeTstack`; at every R branch, `RTop` of
`cur`, `nxt` (on the state after `loop1Type`). -/
partial def loop1Emu (D : Dfs) (d : Nat) (edgeDir : Bool) (s : WalkState) (fuel : Nat) :
    WalkState × List String :=
  if fuel = 0 then (s, []) else
  if WalkState.result (loop1Cond d) s then
    let (ty, s₁) := (loop1Type d edgeDir).run s
    let bad := if ty == .R then match s₁.tstack with
      | cur :: nxt :: _ => "stat:rclose" :: (entryR D s₁ cur).map ("rclose-cur " ++ ·) ++
          (entryR D s₁ nxt).map ("rclose-nxt " ++ ·)
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
`finishEdge_tree_top_settled` (a) (the top after the P-check, first-edge case), computed on the
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
    let sP := WalkState.feP v d o s
    let vE := belowList sP (vertItem v)
    for t in sP.tstack do
      if (t.spans.1 ++ t.spans.2).any fun i => (belowList sP i).any vE.contains then
        bad := s!"vert_own {showT t}" :: bad
    match sP.tstack with
    | t :: _ => bad := settledEntry D v d sP t "ptop" ++ bad
    | [] => pure ()
  return bad.reverse

def isBlock (g : Graph) : Bool := (List.range g.nv).all fun v => Id.run do
  let mut seen := if g.ne == 0 then [] else [0]
  for _ in List.range g.ne do
    seen := (List.range g.ne).filter fun e => seen.contains e || seen.any fun f =>
      [g.edges[e]!.1, g.edges[e]!.2].any fun u =>
        u != v && (g.edges[f]!.1 == u || g.edges[f]!.2 == u)
  return seen.length == g.ne

mutual
partial def tree (D : Dfs) (t : DfsTree) (d : Nat) (s : WalkState) (fr : List (Nat × Nat × Nat)) :
    WalkState × List String × Nat :=
  match t with
  | .node v os =>
    let s := { s with stackVerts := s.stackVerts.set! d v }
    let (hv, s, bad, n) := outs D v d os false s fr
    let s := if hv then s else (setStackDir d true *> pushVertTstack v d).run s |>.2
    (s, bad, n)

partial def outs (D : Dfs) (v d : Nat) (os : List DfsOut) (hv : Bool) (s : WalkState)
    (fr : List (Nat × Nat × Nat)) : Bool × WalkState × List String × Nat :=
  match os with
  | [] => (hv, s, [], 0)
  | o :: rest =>
    let (hv, s) := (walkOutPre v d o hv).run s
    let orig := s.tstack.length
    let (s, bad, n) := match o with
      | .back .. => (s, [], 0)
      | .tree _ _ child => tree D child (d + 1) { s with firstOccurrence := s.firstOccurrence.set! d s.g.ne }
        ((v, d, orig) :: fr)
    let site : Bool := match o.cls with | .ret lv _ => decide (lv < d) | _ => false
    let tag := s!"v={v} d={d} e={o.e} {if o.isTree then "tree" else "back"} orig={orig} hv={hv} t1={o.cls.isType1}: "
    let base := s.tstack.drop (s.tstack.length - orig)
    let pre := if site then
      tagged s!"preA {tag}" (checkEntries D d base s (fun t => t.vStart != v)) ++
      tagged s!"preB {tag}" ("stat:preB-site" :: checkEntries D d base s (fun t => t.vStart != v && t.topDepth ≥ d)) ++
      tagged s!"preD {tag}" (disjoint s) else []
    let rcl := if site && o.isTree then
      let sp := WalkState.after (pushEdgeTstack o.dest d o.e) (WalkState.feS₀ d o s)
      let (sEnd, rbad) := loop1Emu D d s.stackDir[d]! sp sp.tstack.length
      (if sEnd.tstack.map showT != (WalkState.feS₁ d o s).tstack.map showT then [s!"emu-mismatch {tag}"] else []) ++
        tagged "" (rbad.map fun b => if b.startsWith "stat:" then b else s!"{b} {tag}") else []
    let shp := if site && o.isTree then
      tagged s!"shape {tag}" ("stat:shape-site" :: shapeCheck D v d o orig hv s) else []
    let sPre := s
    let basePre := if site then fr.flatMap fun (p, dp, n0) =>
      let tagB := s!"base p={p} dp={dp} n0={n0} {tag}"
      let bot := s.tstack.drop (s.tstack.length - n0)
      (if n0 > orig then [s!"{tagB}n0 > orig"] else []) ++
      (if hv && n0 + 1 > orig then [s!"{tagB}hasVert but n0 + 1 > orig"] else []) ++
      (if o.isTree && orig + 3 > (WalkState.feS₂ d o s).tstack.length then
        [s!"{tagB}close: orig + 3 > len(feS₂)={(WalkState.feS₂ d o s).tstack.length}"] else []) ++
      "stat:base-site" :: bot.flatMap (fun t => settledEntry D p dp s t s!"{tagB}pre") else []
    let (hv, s) := (finishEdge v d o orig hv).run s
    let basePost := if site then fr.flatMap fun (p, dp, n0) =>
      let tagB := s!"base p={p} dp={dp} n0={n0} {tag}"
      let bot := sPre.tstack.drop (sPre.tstack.length - n0)
      let bot' := s.tstack.drop (s.tstack.length - n0)
      (if bot'.map showT != bot.map showT then [s!"{tagB}bottom {n0} changed: {bot.map showT} -> {bot'.map showT}"] else []) ++
      bot'.flatMap (fun t => settledEntry D p dp s t s!"{tagB}post") else []
    let nB := (s.tstack.filter fun t => t.vStart != v && t.topDepth ≥ d).length
    let post := if site then
      tagged s!"postA {tag}" (checkEntries D d s.tstack s (fun t => t.vStart != v)) ++
      tagged s!"postB {tag}" ((List.replicate nB "stat:postB-entry") ++ checkEntries D d s.tstack s (fun t => t.vStart != v && t.topDepth ≥ d)) ++
      tagged s!"postD {tag}" (disjoint s) else []
    let (hv, s, bad', n') := outs D v d rest hv s fr
    (hv, s, bad ++ pre ++ rcl ++ shp ++ basePre ++ post ++ basePost ++ bad', n + n' + (if site then 1 else 0))
end

def runGraph (g : Graph) (vo eo : List Nat) (tern : Bool) : List String × Nat := Id.run do
  let forest := g.dfsForest vo eo
  let D := mkDfs g forest
  let mut s := WalkState.init g tern
  let mut bad := []
  let mut n := 0
  for t in forest do
    let (s', bs, n') := tree D t 0 s []
    bad := bad ++ bs
    n := n + n'
    s := (do
      let top ← popTstack
      modifyItem rootItem fun it => { it with ch := it.ch ++ top.spans.2 }).run s' |>.2
  let ref := g.walk tern forest
  if s.items.toList.map (fun i => (i.type, i.vs, i.ch)) !=
      ref.items.toList.map (fun i => (i.type, i.vs, i.ch)) then
    bad := "instrumented walk differs" :: bad
  return (bad, n)

end RFinishEdgeCheck

/-- Input is concatenated `gen.py` cases, prefixed by the number of cases. -/
def main : IO UInt32 := do
  let input ← (← IO.getStdin).readToEnd
  let toks := ((input.splitOn " ").flatMap (·.splitOn "\n") |>.filter (· != "") |>.map String.toNat!).toArray
  let mut p := 1
  let mut fails := 0
  let mut blocks := 0
  let mut sites := 0
  let mut expectedA := 0
  let mut stats : Std.HashMap String Nat := {}
  let mut cases : List (Nat × Graph × List Nat × List Nat) := []
  for seed in [0:toks[0]!] do
    let nv := toks[p]!; let ne := toks[p+1]!
    p := p + 3
    let edges := (List.range ne).map fun i => (toks[p+2*i]!, toks[p+2*i+1]!)
    p := p + 2*ne
    let k := toks[p]!; p := p + 1
    let vo := (List.range k).map fun i => toks[p+i]!
    p := p + k
    let k := toks[p]!; p := p + 1
    let eo := (List.range k).map fun i => toks[p+i]!
    p := p + k
    cases := (seed, ⟨nv, edges.toArray⟩, vo, eo) :: cases
  -- the triangle with a tripled edge (`checks/RInvReturnCheck.lean`), as a fixed regression
  cases := (1000, ⟨3, #[(0, 1), (1, 2), (1, 2), (1, 2), (0, 2)]⟩, [], []) :: cases
  -- extra corpus: 6000 small random multigraphs (the second half denser) (LCG), kept below when they are blocks
  let mut x : Nat := 12345
  let step := fun (x : Nat) => (x * 1103515245 + 12345) % 2147483648
  for i in [0:6000] do
    x := step x
    let nv := 3 + (x / 65536) % 6
    x := step x
    let ne := if i < 3000 then nv + (x / 65536) % (nv + 5) else 2 * nv + (x / 65536) % (nv + 2)
    let mut es : List (Nat × Nat) := []
    for _ in [0:ne] do
      x := step x
      let u := (x / 65536) % nv
      x := step x
      let w := (x / 65536) % nv
      es := (u, w) :: es
    let vo := if i % 2 == 1 then (List.range nv).reverse else []
    let eo := if i % 3 == 1 then (List.range ne).reverse else if i % 3 == 2 then (List.range ne).rotateLeft (ne / 2) else []
    cases := (2000 + i, ⟨nv, es.toArray⟩, vo, eo) :: cases
  for (seed, g, vo, eo) in cases.reverse do
    if RFinishEdgeCheck.isBlock g && g.ne ≥ 2 then
      blocks := blocks + 1
      for tern in [false, true] do
        let (bad, n) := RFinishEdgeCheck.runGraph g vo eo tern
        sites := sites + n
        for b in bad do
          if b.startsWith "stat:" then stats := stats.insert b (stats.getD b 0 + 1)
        let bad := bad.filter fun b => !(b.startsWith "stat:")
        let isA := fun (b : String) => b.startsWith "preA" || b.startsWith "postA"
        let badA := bad.filter isA
        let bad := bad.filter fun b => !isA b
        if !badA.isEmpty then
          expectedA := expectedA + 1
          IO.println s!"expected-old-contract seed={seed} tern={tern}: nv={g.nv} edges={g.edges} vo={vo} eo={eo}"
          for b in badA do IO.println s!"  {b}"
        if !bad.isEmpty then
          fails := fails + 1
          IO.println s!"seed={seed} tern={tern}: nv={g.nv} edges={g.edges} vo={vo} eo={eo}"
          for b in bad do IO.println s!"  {b}"
  IO.println s!"stats={stats.toList}"
  IO.println s!"cases={toks[0]!} + fixed + 6000 extra; block cases={blocks}; finishEdge sites checked (both ternarize)={sites}; old-contract-A runs failing (expected)={expectedA}; failures (contract B / disjointness / loop-1 emulation / R-branch RTop / FinishRShape / RInvG base frame)={fails}"
  return if fails == 0 then 0 else 1
