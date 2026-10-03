import Spqr.Walk
import Spqr.SepPair
import Spqr.Ranges
/-!
# Empirical check of `WalkState.RangesInv` (`Spqr/RangesInv.lean`)

Run with `lake env lean checks/RangesInvCheck.lean` from `ports/spqr/lean` (not part of the
library). Re-implements `walkTree` around the library's `finishEdge` (as `checks/EarCheck.lean`)
and, right before every `finishEdge curV d o origTstack hasVert`, evaluates each field of
`RangesInv σ n (d+1) s` with `σ = edgePostorderForest forest`, `n = pos o.e` (the edges processed
so far are exactly `σ.take n`), and `Saturated k s` with `k = tstack.length - origTstack` (the
`Frontier` split).
-/
open Spqr WalkM
instance : Inhabited Spqr.DfsTree := ⟨.node 0 []⟩
instance : Inhabited Spqr.DfsOut := ⟨.back 0 0 .selfLoop⟩

partial def edgesBelow (s : WalkState) (i : ItemId) : List Nat :=
  (if 1 + s.g.nv ≤ i ∧ i < 1 + s.g.nv + s.g.ne then [i - 1 - s.g.nv] else []) ++
    (s.items[i]!.ch).flatMap (edgesBelow s)
partial def pieceBelow (s : WalkState) (i : ItemId) : List Nat :=
  (if 1 + s.g.nv ≤ i ∧ i < 1 + s.g.nv + s.g.ne then [i - 1 - s.g.nv] else []) ++
    ((s.items[i]!.ch).filter fun c => s.items[c]!.type != .V).flatMap (pieceBelow s)
def spanItems (t : TEntry) : List ItemId := t.spans.1 ++ t.spans.2
def entryEdges (s : WalkState) (t : TEntry) : List Nat := (spanItems t).flatMap (edgesBelow s)
def entryPiece (s : WalkState) (t : TEntry) : List Nat :=
  ((spanItems t).filter fun i => s.items[i]!.type != .V).flatMap (pieceBelow s)
def inc (s : WalkState) (e v : Nat) : Bool := s.g.edges[e]!.1 == v || s.g.edges[e]!.2 == v
def touches (s : WalkState) (E : List Nat) (v : Nat) : Bool := E.any (inc s · v)
/-- Attachment vertices of `E` (touching an edge in `E` and one outside). -/
def atts (s : WalkState) (E : List Nat) : List Nat :=
  (List.range s.g.nv).filter fun v =>
    touches s E v && (List.range s.g.ne).any fun e => inc s e v && !E.contains e
def showT (t : TEntry) : String := s!"({t.vStart},{t.topDepth},{t.firstIdx},{t.spans})"

structure V where
  seed : Nat
  tern : Bool
  curV : Nat
  d : Nat
  o : String
  kind : String
  info : String
deriving Repr

def check (seed : Nat) (σ : List Nat) (curV d : Nat) (o : DfsOut) (orig : Nat) (s : WalkState) : List V := Id.run do
  let pos := fun e => σ.idxOf e
  let n := pos o.e
  let D := d + 1
  let oStr := s!"{if o.cls.isTree then "tree" else "back"} e={o.e} dest={o.dest} lv={o.cls.lowval d}"
  let mut out : List V := []
  let bad (k : String) (info : String) : V :=
    ⟨seed, s.ternarize, curV, d, oStr, k, s!"{info} | orig={orig} n={n} sv={s.stackVerts.toList.take (d+2)} σ={σ} stack={s.tstack.map showT}"⟩
  let E := fun t => entryEdges s t
  let P := fun t => entryPiece s t
  let stk := s.tstack
  -- processed
  for t in stk do
    for e in E t do
      if !(pos e < n) then out := bad "processed" s!"{showT t} e={e}" :: out
  -- ordered (top-first: pieces of higher entries are later)
  let rec ord : List TEntry → List V
    | [] => []
    | t :: rest =>
      (rest.flatMap fun t' => (P t).flatMap fun e => (P t').filterMap fun e' =>
        if pos e' < pos e then none else some (bad "ordered" s!"{showT t} e={e} below {showT t'} e'={e'}")) ++ ord rest
  out := ord stk ++ out
  -- open convexity (holes are nobody's piece); candidate: holes are below the entry
  for t in stk do
    let ps := (P t).map pos
    match ps.min?, ps.max? with
    | some lo, some hi =>
      for b in List.range' lo (hi + 1 - lo) do
        let e := σ[b]!
        if !(E t).contains e && stk.any (fun t' => (P t').contains e) then
          out := bad "cand_piece_consecutive" s!"{showT t} hole e={e}" :: out
        if !(E t).contains e then out := bad "convex" s!"{showT t} hole e={e}" :: out
    | _, _ => pure ()
  -- cover: processed edges are in some entry or below the root
  for b in List.range n do
    let e := σ[b]!
    if !stk.any (fun t => (E t).contains e) && !(edgesBelow s rootItem).contains e &&
        !(List.range (d+1)).any (fun k => (edgesBelow s (vertItem s.stackVerts[k]!)).contains e) then
      out := bad "cand_cover" s!"e={e}" :: out
  -- terms: ≤ 2 attachments, in Term D (vStart or path[topDepth..D]); candidate variants
  for t in stk do
    let A := atts s (E t)
    let term := fun v => v == t.vStart || (List.range (D+1)).any fun k => t.topDepth ≤ k && v == s.stackVerts[k]!
    let path := fun v => v == t.vStart || (List.range (D+1)).any fun k => v == s.stackVerts[k]!
    if A.length > 2 then out := bad "cand_terms_two" s!"{showT t} atts={A}" :: out
    if !A.all term then out := bad "cand_terms_term" s!"{showT t} atts={A}" :: out
    if !A.all path then out := bad "cand_terms_path" s!"{showT t} atts={A}" :: out
  -- closed items
  for i in List.range s.items.size do
    let it := s.items[i]!
    if 1 + s.g.nv + s.g.ne ≤ i then
      match it.vs with
      | (some u, some v) =>
        let A := atts s (edgesBelow s i)
        if !A.all (fun x => x == u || x == v) then out := bad "closed_att" s!"i={i} vs={it.vs} atts={A}" :: out
      | _ => pure ()
    if it.type != .F && it.type != .V then
      let ps := (pieceBelow s i).map pos
      let below := edgesBelow s i
      match ps.min?, ps.max? with
      | some lo, some hi =>
        for b in List.range' lo (hi + 1 - lo) do
          if !below.contains σ[b]! then out := bad "closed_convex" s!"i={i} {repr it.type} hole e={σ[b]!}" :: out
      | _, _ => pure ()
  -- saturation below the frontier, over runs (and the pairwise variant, for comparison)
  let k := stk.length - orig
  let base := stk.drop k
  let rec runs : List TEntry → List (List TEntry)
    | [] => []
    | t :: rest => ((List.range rest.length).map fun j => t :: rest.take (j+1)) ++ runs rest
  for r in runs base do
    let U := r.flatMap E
    let A := atts s U
    if U ≠ [] && A.length ≤ 2 then
      if r.length == 2 then out := bad "cand_sat_pair" s!"run={r.map showT} atts={A}" :: out
      else out := bad "cand_sat_run" s!"run={r.map showT} atts={A}" :: out
  return out

mutual
partial def iTree (seed : Nat) (σ : List Nat) (t : DfsTree) (d : Nat) (s : WalkState) : WalkState × List V :=
  match t with
  | .node v outs =>
    let s := { s with stackVerts := s.stackVerts.set! d v }
    let (hv, s, vs) := iOuts seed σ v d outs false s
    let s := if hv then s else ((setStackDir d true *> pushVertTstack v d).run s).2
    (s, vs)
partial def iOuts (seed : Nat) (σ : List Nat) (v d : Nat) (outs : List DfsOut) (hv : Bool) (s : WalkState) : Bool × WalkState × List V :=
  match outs with
  | [] => (hv, s, [])
  | o :: rest =>
    let (hv, s, vs) := iOut seed σ v d o hv s
    let (hv', s', vs') := iOuts seed σ v d rest hv s
    (hv', s', vs ++ vs')
partial def iOut (seed : Nat) (σ : List Nat) (v d : Nat) (o : DfsOut) (hv : Bool) (s : WalkState) : Bool × WalkState × List V :=
  let lowval := o.cls.lowval d
  let s := ((do let lowDir ← stackDir lowval; setStackDir d (if lowval ≥ d then false else !lowDir) : WalkM Unit).run s).2
  let (hv, s) := if !hv && lowval < d && o.cls.isType1 then (true, ((pushVertTstack v d).run s).2) else (hv, s)
  let orig := s.tstack.length
  let (s, vs) := match o with
    | .tree _ _ child => iTree seed σ child (d+1) { s with firstOccurrence := s.firstOccurrence.set! d s.g.ne }
    | .back .. => (s, [])
  let vs := vs ++ check seed σ v d o orig s
  let (hv', s) := (finishEdge v d o orig hv).run s
  (hv', s, vs)
end

def iForest (seed : Nat) (σ : List Nat) (forest : List DfsTree) (s : WalkState) : WalkState × List V :=
  forest.foldl (fun (s, vs) t =>
    let (s, vs') := iTree seed σ t 0 s
    let s := ((popTstack >>= fun top => modifyItem rootItem fun it => { it with ch := it.ch ++ top.spans.2 }).run s).2
    (s, vs ++ vs')) (s, [])

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

def runGraph (seed : Nat) (g : Graph) (tern : Bool) : Bool × List V :=
  let f := g.dfsForest [] []
  let σ := edgePostorderForest f
  let (s, vs) := iForest seed σ f (WalkState.init g tern)
  let ref := g.walk tern f
  (s.items.toList.map (fun it => (it.vs, it.ch)) == ref.items.toList.map (fun it => (it.vs, it.ch)) && s.tstack.length == ref.tstack.length, vs)

def tiny : List Graph :=
  [⟨2, #[(0,1)]⟩, ⟨1, #[(0,0)]⟩, ⟨3, #[]⟩, ⟨4, #[(0,1),(0,2),(0,3)]⟩, ⟨5, #[(0,1),(1,2),(2,0),(2,3),(3,4),(4,2)]⟩,
   ⟨3, #[(0,1),(1,2),(2,0)]⟩, ⟨2, #[(0,1),(1,0)]⟩, ⟨4, #[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]⟩]

def summarize (lo hi : Nat) : IO Unit := do
  let mut counts : List (String × Nat) := []
  let mut okAll := true
  let mut shown : List (String × String) := []
  let gs := (tiny.map fun g => (1000000, g)) ++ (List.range' lo (hi - lo)).map fun seed => (seed, randGraph seed)
  for (seed, g) in gs do
    for tern in [false, true] do
      let (ok, vs) := runGraph seed g tern
      if !ok then okAll := false
      for v in vs do
        counts := match counts.find? (·.1 == v.kind) with
          | some _ => counts.map fun (k, n) => if k == v.kind then (k, n+1) else (k, n)
          | none => counts ++ [(v.kind, 1)]
        if (shown.filter (·.1 == v.kind)).length < 2 then
          shown := shown ++ [(v.kind, s!"{repr v} g={repr g.edges}")]
  IO.println s!"instrumentation matches library walk: {okAll}"
  let req := counts.filter fun (k, _) => !k.startsWith "cand_"
  IO.println s!"required failures: {req.foldl (fun a (_, n) => a + n) 0} {req}"
  IO.println s!"candidate observations (not fields of RangesInv, see PROOF.md §4.6): {counts.filter fun (k, _) => k.startsWith "cand_"}"
  for (_, l) in shown do IO.println l
#eval summarize 0 400
