import Spqr.Ear
import Spqr.ItemSpec
/-!
# Empirical check of the walk's attachment invariant (`Spqr/EarInv.lean`)

Run with `lake env lean checks/InvCheck.lean` from `ports/spqr/lean` (not part of the library).

`iTree`/`iOuts`/`iOut` re-implement `walkTree`/`walkOuts`/`walkOut` around the library's
`finishEdge` and record, at every `walkTree` start/end and after every `walkOut`, the open entries
whose attachment vertices (boundary of `Items.EdgeBelow` over the spans) violate
* `A(d)`   : `TEntry.Term d`     (`WalkSpec.EntryInv d`, `d` the current depth);
* `A(d+1)` : `TEntry.Term (d+1)` (`WalkSpec.EntryInv (d+1)`, the index `walkTree_inv'` uses under a
             tree edge);
* `B(d)`   : `TEntry.Term' d`    (`EarInv.EntryInv' d`: `Term d` plus `vStart`s of the entries above);
* `conn`   : connectivity of the entry's edge set.
`cexBuried` is the hand counterexample of `EarInv.lean`; `summarize` is the random sweep.
-/
open Spqr WalkM
instance : Inhabited Spqr.DfsTree := ⟨.node 0 []⟩
instance : Inhabited Spqr.DfsOut := ⟨.back 0 0 .selfLoop⟩
instance : Inhabited Spqr.Frame := ⟨⟨0, 0, .back 0 0 .selfLoop, 0, false, []⟩⟩

partial def edgesBelow (s : WalkState) (i : ItemId) : List Nat :=
  (if 1 + s.g.nv ≤ i ∧ i < 1 + s.g.nv + s.g.ne then [i - 1 - s.g.nv] else []) ++
    (s.items[i]!.ch).flatMap (edgesBelow s)
def entryEdges (s : WalkState) (t : TEntry) : List Nat := (t.spans.1 ++ t.spans.2).flatMap (edgesBelow s)
def inc (s : WalkState) (e v : Nat) : Bool := s.g.edges[e]!.1 == v || s.g.edges[e]!.2 == v
def boundary (s : WalkState) (E : List Nat) : List Nat :=
  (List.range s.g.nv).filter fun v => E.any (fun e => inc s e v) && (List.range s.g.ne).any fun e' => !E.contains e' && inc s e' v
def termA (D : Nat) (s : WalkState) (t : TEntry) (v : Nat) : Bool :=
  v == t.vStart || (List.range (D+1)).any fun k => t.topDepth ≤ k && s.stackVerts[k]! == v
def termB (D : Nat) (s : WalkState) (t : TEntry) (above : List TEntry) (v : Nat) : Bool :=
  termA D s t v || above.any fun t' => t'.vStart == v
partial def connected (s : WalkState) (E : List Nat) : Bool :=
  match E with
  | [] => true
  | e :: _ =>
    let rec go (seen : List Nat) (todo : List Nat) (fuel : Nat) : List Nat :=
      match fuel, todo with
      | 0, _ | _, [] => seen
      | f + 1, x :: rest =>
        let nb := E.filter fun e' => !seen.contains e' && (inc s e' x) 
        let seen := seen ++ nb
        let vs := nb.flatMap fun e' => [s.g.edges[e']!.1, s.g.edges[e']!.2]
        go seen (rest ++ vs) f
    let seen := go [e] [s.g.edges[e]!.1, s.g.edges[e]!.2] (4 * E.length + 4)
    E.all seen.contains

structure V where
  site : String
  d : Nat
  idx : Nat
  entry : Nat × Nat × (List Nat × List Nat)
  edges : List Nat
  bad : List Nat
  kind : String
deriving Repr

def check (site : String) (d : Nat) (s : WalkState) : List V := Id.run do
  let mut out : List V := []
  let mut above : List TEntry := []
  let mut i := 0
  for t in s.tstack do
    let E := entryEdges s t
    let B := boundary s E
    let mk (k : String) (bad : List Nat) : V := ⟨site, d, i, (t.vStart, t.topDepth, t.spans), E, bad, k⟩
    if !connected s E then out := mk "conn" [] :: out
    let badA := B.filter fun v => !termA d s t v
    if badA ≠ [] then out := mk "A(d)" badA :: out
    let badA1 := B.filter fun v => !termA (d+1) s t v
    if badA1 ≠ [] then out := mk "A(d+1)" badA1 :: out
    let badB := B.filter fun v => !termB d s t above v
    if badB ≠ [] then out := mk "B(d)" badB :: out
    above := above ++ [t]
    i := i + 1
  return out

mutual
partial def iTree (t : DfsTree) (d : Nat) (s : WalkState) : WalkState × List V :=
  match t with
  | .node v outs =>
    let s := { s with stackVerts := s.stackVerts.set! d v }
    let vs0 := check "treeStart" d s
    let (hv, s, vs) := iOuts v d outs false s
    let s := if hv then s else ((setStackDir d true *> pushVertTstack v d).run s).2
    (s, vs0 ++ vs ++ check "treeEnd" d s)
partial def iOuts (v d : Nat) (outs : List DfsOut) (hv : Bool) (s : WalkState) : Bool × WalkState × List V :=
  match outs with
  | [] => (hv, s, [])
  | o :: rest =>
    let (hv, s, vs) := iOut v d o hv s
    let (hv', s', vs') := iOuts v d rest hv s
    (hv', s', vs ++ vs')
partial def iOut (v d : Nat) (o : DfsOut) (hv : Bool) (s : WalkState) : Bool × WalkState × List V :=
  let lowval := o.cls.lowval d
  let s := ((do let lowDir ← stackDir lowval; setStackDir d (if lowval ≥ d then false else !lowDir) : WalkM Unit).run s).2
  let (hv, s) := if !hv && lowval < d && o.cls.isType1 then (true, ((pushVertTstack v d).run s).2) else (hv, s)
  let orig := s.tstack.length
  let (s, vs) := match o with
    | .tree _ _ child => iTree child (d+1) { s with firstOccurrence := s.firstOccurrence.set! d s.g.ne }
    | .back .. => (s, [])
  let (hv', s) := (finishEdge v d o orig hv).run s
  (hv', s, vs ++ check "afterOut" d s)
end

def iForest (forest : List DfsTree) (s : WalkState) : WalkState × List V :=
  forest.foldl (fun (s, vs) t =>
    let (s, vs') := iTree t 0 s
    let s := ((popTstack >>= fun top => modifyItem rootItem fun it => { it with ch := it.ch ++ top.spans.2 }).run s).2
    (s, vs ++ vs')) (s, [])

/-! ### The hand counterexample -/

namespace cexBuried
def g : Graph := ⟨8, #[(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(6,0),(6,1),(5,2),(4,7),(7,3)]⟩
def f := g.dfsForest [] []
/-- The subtree at vertex `2` (depth 2): the type-2 chain `2→3→4→5→6`. -/
def t2 : DfsTree :=
  match f.head! with
  | .node _ [.tree _ _ (.node _ [.tree _ _ c])] => c
  | _ => default
def s0 : WalkState :=
  let s := WalkState.init g false
  { s with stackVerts := (s.stackVerts.set! 0 0).set! 1 1,
           firstOccurrence := (s.firstOccurrence.set! 0 g.ne).set! 1 g.ne }
def fr : List Frame × WalkState := (descend 20 t2 2 []).run s0
/-- After the frame `(5, 5)`. -/
def s1 : WalkState := ((ascend 20 (fr.1.take 1)).run fr.2).2
def f4 : Frame := fr.1[1]!
/-- After `finishEdge` of the frame `(4, 4)` (type-2 edge `4→5`). -/
def s2 : WalkState := ((finishEdge f4.v f4.d f4.o f4.origTstack f4.hasVert).run s1).2
/-- The start of `walkTree (node 7) 5` for the sibling `4→7`. -/
def s3 : WalkState := { s2 with stackVerts := s2.stackVerts.set! 5 7 }
def showT (s : WalkState) : List (Nat × Nat × Nat × (List Nat × List Nat)) :=
  s.tstack.map fun t => (t.vStart, t.topDepth, t.firstIdx, t.spans)
#eval (fr.1.map fun f => (f.v, f.d, f.origTstack, f.hasVert), showT s2, s2.stackVerts)
-- `(6, 5, [], [Q(5,6)])` violates `Term 4` after the frame, and `Term D` for every `D` at `s3`:
#eval (check "s2" 4 s2).map fun v => (v.kind, v.entry, v.bad)
#eval ((check "s3" 5 s3).map fun v => (v.kind, v.entry, v.bad), (check "s3" 6 s3).map fun v => (v.kind, v.entry, v.bad))
end cexBuried

/-! ### Random sweep -/

def lcg (x : Nat) : Nat := (x * 6364136223846793005 + 1442695040888963407) % 2^64
def randGraph (seed : Nat) : Graph := Id.run do
  let mut x := lcg (seed + 12345)
  let nv := 2 + (x >>> 33) % 11
  x := lcg x
  let ne := 1 + (x >>> 33) % 24
  let mut es : Array (Nat × Nat) := #[]
  for _ in List.range ne do
    x := lcg x
    let a := (x >>> 33) % nv
    x := lcg x
    let b := (x >>> 33) % nv
    es := es.push (a, b)
  return ⟨nv, es⟩

def runSeed (seed : Nat) : Bool × Nat × List V :=
  let g := randGraph seed
  let f := g.dfsForest [] []
  let (s, vs) := iForest f (WalkState.init g false)
  let ref := g.walk false f
  (s.items.toList.map (fun it => (it.vs, it.ch)) == ref.items.toList.map (fun it => (it.vs, it.ch)) && s.tstack.length == ref.tstack.length, s.tstack.length, vs)

def summarize (lo hi : Nat) : IO Unit := do
  let mut counts : List (String × Nat) := []
  let mut okAll := true
  let mut shown : List String := []
  for seed in List.range' lo (hi - lo) do
    let (ok, _, vs) := runSeed seed
    if !ok then okAll := false
    for v in vs do
      counts := match counts.find? (·.1 == v.kind) with
        | some _ => counts.map fun (k, n) => if k == v.kind then (k, n+1) else (k, n)
        | none => counts ++ [(v.kind, 1)]
      if (v.kind == "B(d)" || v.kind == "conn") && shown.length < 4 then
        shown := shown ++ [s!"seed {seed} {repr v} g={repr (randGraph seed).edges}"]
  IO.println s!"instrumentation matches library walk: {okAll}"
  IO.println s!"violations by kind: {counts}"
  for l in shown do IO.println l

#eval summarize 0 6000
