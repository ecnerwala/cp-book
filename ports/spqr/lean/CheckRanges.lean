import Spqr
import Spqr.Ranges

open Spqr

/-! Empirical check of every field of `Items.Ranges` (plus some candidate strengthenings, reported
as `cand:` and not counted) on the walk's output; input format of `gen.py`.  `--dump` prints the
raw items with their edge intervals. -/

def pairEq (p q : Nat × Nat) : Bool := p == q || p == (q.2, q.1)

partial def pieceEdges (items : Items) (nv ne : Nat) (i : ItemId) : List Nat :=
  (if nv + 1 ≤ i && i < nv + 1 + ne then [i - nv - 1] else []) ++
    ((Items.ch items i).filter fun c => Items.type items c != .V).flatMap (pieceEdges items nv ne)

partial def edgesBelow (items : Items) (nv ne : Nat) (i : ItemId) : List Nat :=
  (if nv + 1 ≤ i && i < nv + 1 + ne then [i - nv - 1] else []) ++ (Items.ch items i).flatMap (edgesBelow items nv ne)

def main (args : List String) : IO Unit := do
  let input ← (← IO.getStdin).readToEnd
  let toks := (input.splitOn " ").flatMap (·.splitOn "\n") |>.filter (· ≠ "") |>.map String.toNat!
  let toks := toks.toArray
  let mut p := 0
  let next := fun (p : Nat) => (toks[p]!, p + 1)
  let (nv, p1) := next p; let (ne, p2) := next p1; let (tern, p3) := next p2
  p := p3
  let mut edges : Array (Nat × Nat) := #[]
  for _ in [0:ne] do
    edges := edges.push (toks[p]!, toks[p+1]!)
    p := p + 2
  let k := toks[p]!; p := p + 1
  let vertOrder := (List.range k).map fun i => toks[p + i]!
  p := p + k
  let l := toks[p]!; p := p + 1
  let edgeOrder := (List.range l).map fun i => toks[p + i]!
  let g : Graph := ⟨nv, edges⟩
  let forest := g.dfsForest vertOrder edgeOrder
  let items := (g.walk (tern != 0) forest).items
  let σ := edgePostorderForest forest
  let n := items.size
  let ty := fun i => Items.type items i
  let ch := fun i => Items.ch items i
  let vs := fun i => Items.vs items i
  let isNode := fun i => ty i != .F && ty i != .V
  let isSPR := fun i => [NodeType.S, .P, .R].contains (ty i)
  -- edges below, as a membership table
  let mut belowTbl : Array (Array Bool) := #[]
  for i in [0:n] do
    let mut row := Array.replicate ne false
    for e in edgesBelow items nv ne i do row := row.set! e true
    belowTbl := belowTbl.push row
  let below := fun i e => belowTbl[i]![e]!
  let mut pieceTbl : Array (Array Bool) := #[]
  for i in [0:n] do
    let mut row := Array.replicate ne false
    for e in pieceEdges items nv ne i do row := row.set! e true
    pieceTbl := pieceTbl.push row
  let piece := fun i e => pieceTbl[i]![e]!
  let inc := fun e v => edges[e]!.1 == v || edges[e]!.2 == v
  let isVs := fun i v => (vs i).1 == some v || (vs i).2 == some v
  let att := fun i v => (List.range ne).any (fun e => inc e v && below i e) &&
    (List.range ne).any (fun e => inc e v && !below i e)
  let inner := fun i v => (List.range ne).any (fun e => inc e v) &&
    (List.range ne).all (fun e => !inc e v || below i e)
  let allBelow := fun i v => (List.range ne).all (fun e => !inc e v || below i e)
  let pos : Array Nat := Id.run do
    let mut a := Array.replicate ne ne
    for j in [0:σ.length] do a := a.set! σ[j]! j
    return a
  let interval := fun i => ((List.range ne).filter (piece i ·)).map (pos[·]!)
  let fullInterval := fun i => ((List.range ne).filter (below i ·)).map (pos[·]!)
  let isConvex := fun (ps : List Nat) => ps.isEmpty || ps.max?.getD 0 + 1 - ps.min?.getD 0 == ps.length
  let convexAt := fun i =>
    let ps := interval i
    ps.isEmpty || (List.range ne).all fun b =>
      !(ps.min?.getD 0 ≤ pos[b]! && pos[b]! ≤ ps.max?.getD 0) || below i b
  if args.contains "--dump" then
    IO.println s!"sigma {σ}"
    for i in [0:n] do
      IO.println s!"{i}: {repr (ty i)} vs={vs i} ch={ch i} E={(List.range ne).filter (below i ·)} piecePos={interval i}"
  let mut bad := 0
  let mut cand := 0
  let fail := fun (msg : String) => IO.println s!"FAIL {msg}"
  -- sanity: σ is a permutation of the edges
  if σ.length != ne || (List.range ne).any (fun e => pos[e]! == ne) then
    bad := bad + 1; fail s!"sigma not a permutation: {σ}"
  for i in [0:n] do
    if isNode i then
      if !convexAt i then bad := bad + 1; fail s!"convex: i={i} {repr (ty i)} pos={interval i}"
      if !isConvex (interval i) then cand := cand + 1; IO.println s!"cand: strict piece convexity: i={i} {repr (ty i)}"
      for v in [0:nv] do
        if att i v && !isVs i v then bad := bad + 1; fail s!"att_vs: i={i} {repr (ty i)} v={v} vs={vs i}"
    else
      pure ()
    if isSPR i then
      for v in [0:nv] do
        if isVs i v && !att i v then bad := bad + 1; fail s!"vs_att: i={i} {repr (ty i)} v={v}"
        let rhs := inner i v && (ch i).all fun c => !allBelow c v
        if (ch i).contains (vertItem v) != rhs then
          bad := bad + 1; fail s!"interior: i={i} {repr (ty i)} v={v} rhs={rhs}"
      match vs i with
      | (some u, some v) => if u == v then bad := bad + 1; fail s!"vs_ne: i={i} vs={vs i}"
      | _ => bad := bad + 1; fail s!"vs_ne shape: i={i} vs={vs i}"
      for c in ch i do
        if ty c != .V then
          match vs c with
          | (some a, some b) => if a == b then bad := bad + 1; fail s!"child_two: p={i} c={c} vs={vs c}"
          | _ => bad := bad + 1; fail s!"child_two: p={i} c={c} {repr (ty c)} vs={vs c}"
    else if ty i == .Q && ch i == [] then
      for v in [0:nv] do
        if isVs i v && !att i v then bad := bad + 1; fail s!"vs_att at Q leaf: i={i} v={v}"
    for c in ch i do
      if (ty c == .I || ty c == .O) && ty i != .Q then bad := bad + 1; fail s!"io_parent: p={i} c={c}"
      if ty i == .V && ch c == [] then bad := bad + 1; fail s!"q_under_v: v={i} c={c}"
      if ty c == .S && ty i == .S then cand := cand + 1; IO.println s!"cand: canonical S/S: p={i} c={c} tern={tern}"
      if ty c == .P && ty i == .P then cand := cand + 1; IO.println s!"cand: canonical P/P: p={i} c={c} tern={tern}"
    if ty i == .P then
      let ve := Items.virtualEdges items i
      if ve.length < 2 then bad := bad + 1; fail s!"p_shape len: i={i} ch={ch i}"
      match vs i with
      | (some u, some v) =>
        for q in ve do
          if !pairEq q (u, v) then bad := bad + 1; fail s!"p_shape pair: i={i} q={q} vs={vs i}"
      | _ => bad := bad + 1; fail s!"p_shape vs: i={i}"
      if (ch i).any (fun c => ty c == .V) then bad := bad + 1; fail s!"p_shape: P with V child {i}"
    if ty i == .S then
      let xs := ((ch i).filter fun c => ty c == .V).map (· - 1)
      match vs i with
      | (some u, some v) =>
        if xs.length < 1 || Items.virtualEdges items i != List.zip (u :: xs) (xs ++ [v]) then
          bad := bad + 1; fail s!"s_order: i={i} vs={vs i} xs={xs} virt={Items.virtualEdges items i}"
      | _ => bad := bad + 1; fail s!"s_order vs: i={i} vs={vs i}"
    if ty i == .R then
      let ve := Items.virtualEdges items i
      let xs := (ch i).filter fun c => ty c == .V
      let keys := ve.map fun q => min q.1 q.2 + nv * max q.1 q.2
      if xs.length < 2 || ve.length < 5 || keys.Nodup = false then
        bad := bad + 1; fail s!"r_shape: i={i} xs={xs.length} ve={ve}"
      match vs i with
      | (some u, some v) => for q in ve do if pairEq q (u, v) then bad := bad + 1; fail s!"r_shape par: i={i} q={q}"
      | _ => bad := bad + 1; fail s!"r_shape vs: i={i}"
  for e in [0:ne] do
    let q := edgeItem g e
    let (x, y) := edges[e]!
    match ch q with
    | [] =>
      match vs q with
      | (some a, some b) => if a == b || !pairEq (a, b) (x, y) then bad := bad + 1; fail s!"q_leaf: e={e} vs={vs q} edge={(x,y)}"
      | _ => bad := bad + 1; fail s!"q_leaf vs: e={e} vs={vs q}"
    | cs =>
      match vs q with
      | (some u, none) =>
        if !inc e u then bad := bad + 1; fail s!"q_root inc: e={e} u={u}"
        match cs with
        | [c] =>
          if x != y || !isNode c || (ty c == .Q && ch c != []) || vs c != (some u, none) then
            bad := bad + 1; fail s!"q_root loop: e={e} c={c} {repr (ty c)} vs c={vs c} edge={(x,y)}"
          if ty c != .O then cand := cand + 1; IO.println s!"cand: loop child not O: e={e}"
        | [c, w] =>
          let ok := x != y && isNode c && (ty c != .Q || ch c == []) && w ≥ 1 && pairEq (u, w - 1) (x, y) &&
            (match vs c with | (some a, some b) => pairEq (a, b) (u, w - 1) | _ => false)
          if !ok then bad := bad + 1; fail s!"q_root pair: e={e} c={c} {repr (ty c)} vs c={vs c} w={w} vs={vs q} edge={(x,y)}"
          if ty w != .V then bad := bad + 1; fail s!"q_root w not V: e={e}"
        | _ => bad := bad + 1; fail s!"q_root ch: e={e} ch={cs}"
      | _ => bad := bad + 1; fail s!"q_root vs: e={e} vs={vs q}"
  IO.println s!"items {n} bad {bad} cand {cand}"
