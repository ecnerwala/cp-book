import Spqr.Walk
import Spqr.SepPair
import Spqr.Ranges
import Spqr.RangesTree
/-!
# Empirical check of `WalkState.RangesInv` (`Spqr/RangesInv.lean`)

Run with `lake env lean checks/RangesInvCheck.lean` from `ports/spqr/lean` (not part of the
library). Re-implements `walkTree` around the library's `finishEdge` (as `checks/EarCheck.lean`)
and, right before every `finishEdge curV d o origTstack hasVert`, evaluates each field of
`RangesInv σ n (d+1) s` with `σ = edgePostorderForest forest`, `n = pos o.e` (the edges processed
so far are exactly `σ.take n`), and `Saturated k s` with `k = tstack.length - origTstack` (the
`Frontier` split).
-/
open Spqr WalkM WalkState
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

def mergeSites (cond : WalkM Bool) (body : WalkM Unit) (s : WalkState) : List WalkState :=
  go (s.tstack.length + 1) s
where
  go : Nat → WalkState → List WalkState
    | 0, _ => []
    | k + 1, s => if result cond s then s :: go k (after body s) else []

def adjacencySites (curV d : Nat) (o : DfsOut) (orig : Nat) (hv : Bool) (s : WalkState) :
    List (String × WalkState) := Id.run do
  let lv := o.cls.lowval d
  if d ≤ lv then return []
  let mut sites := []
  let mut rest := feBack curV lv d o s
  let mut single := true
  if o.cls.isTree then
    let edgeDir := s.stackDir[d]!
    for st in mergeSites (loop1Cond d) (loop1Body d edgeDir) (ceS₁ o.dest d o.e (feS₀ d o s)) do
      if (nxtE st).topDepth > d then sites := sites ++ [("loop1.S", st)]
      sites := sites ++ [("loop1.close", l1S₂ d edgeDir st)]
    let st := feS₁ d o s
    if (curE st).firstIdx > st.firstOccurrence[d]! then
      sites := sites ++ (mergeSites (loop2Cond st.firstOccurrence[d]!) mergeTstackTops st).map ("loop2", ·)
    let st := feS₂ d o s
    let b := feSingle d o s
    if hv then
      if !o.cls.isType1 then
        sites := sites ++ (mergeSites (loop3Cond orig) mergeTstackTops st).map ("loop3", ·)
      sites := sites ++ [("vertex.1", cvS₂ o.cls.isType1 orig b st), ("vertex.2", cvS₃ o.cls.isType1 orig b st)]
      rest := feS₃ curV d o orig s
      single := feB₃ curV d o orig s
    else
      rest := st
      single := b
  if result (condP curV lv o.cls.isType1) rest then
    sites := sites ++ [("P", after (maybeUnwrapNxt .P) rest)]
  if !hv && !single then
    sites := sites ++ [("tail", after (pushVertTstack curV d) (after (finishP curV lv o.cls.isType1) rest))]
  return sites

def checkAdj (seed : Nat) (σ : List Nat) (curV d : Nat) (o : DfsOut) (orig : Nat) (hv : Bool)
    (s : WalkState) : List V := Id.run do
  let mut out := []
  for (site, st) in adjacencySites curV d o orig hv s do
    let cur := st.tstack.head!
    let nxt := st.tstack.tail.head!
    if site == "P" then
      let hi := σ.idxOf o.e + 1
      let lo := hi - o.block.length
      for b in List.range' lo (hi - lo) do
        if !st.tstack.any (fun t => (entryEdges st t).contains σ[b]!) then
          out := ⟨seed, s.ternarize, curV, d, s!"e={o.e}", "p_child_cover", s!"lo={lo} hi={hi} b={b}"⟩ :: out
      for a in (entryPiece st nxt).map σ.idxOf do
        for b in List.range' a (lo - a) do
          if !st.tstack.any (fun t => (entryEdges st t).contains σ[b]!) then
            out := ⟨seed, s.ternarize, curV, d, s!"e={o.e}", "p_base_cover", s!"lo={lo} hi={hi} a={a} b={b}"⟩ :: out
    let low := ((entryPiece st nxt).map σ.idxOf).min?
    let high := ((entryPiece st cur).map σ.idxOf).max?
    match low, high with
    | some lo, some hi =>
      for b in List.range' lo (hi + 1 - lo) do
        if !(entryEdges st cur).contains σ[b]! && !(entryEdges st nxt).contains σ[b]! then
          out := ⟨seed, s.ternarize, curV, d, s!"e={o.e} lv={o.cls.lowval d}", "merge_adj",
            s!"site={site} gap={σ[b]!} orig={orig} σ={σ} stack={st.tstack.map showT}"⟩ :: out
    | _, _ => pure ()
  return out

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
  if d ≤ o.cls.lowval d && o.cls.isTree then
    let front := if o.cls.lowval d == d + 1 then stk.take 1 else stk.take 2
    let sideOk := if o.cls.lowval d == d + 1 then
      stk.head!.spans.1.isEmpty
    else stk.head!.spans.2.isEmpty && stk.tail.head!.spans.1.isEmpty
    if !sideOk then out := bad "boundary_side" s!"front={front.map showT}" :: out
    for t in front do
      for e in P t do
        for b in List.range' (pos e) (n - pos e) do
          if !front.any (fun u => (E u).contains σ[b]!) then
            out := bad "boundary_fill" s!"front={front.map showT} e={e} gap={σ[b]!}" :: out
  -- processed
  for t in stk do
    for e in E t do
      if !(pos e < n) then out := bad "processed" s!"{showT t} e={e}" :: out
  -- ordered (top-first: pieces of higher entries are later)
  let rec ord : List TEntry → List V
    | [] => []
    | t :: rest =>
      (rest.flatMap fun t' => (P t).flatMap fun e => (E t').filterMap fun e' =>
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
  -- candidate: adjacent entries are σ-adjacent (the gap between their pieces is their own edges),
  -- and the variant allowing edges below the unpushed vertex items of the DFS path
  let unpushed := fun (w : Nat) => !stk.any fun t => (spanItems t).contains (vertItem w)
  let rec gaps (above : List TEntry) : List TEntry → List V
    | t :: t' :: rest =>
      (match ((P t').map pos).max?, ((P t).map pos).min? with
      | some hi', some lo =>
        (List.range' (hi' + 1) (lo - hi' - 1)).flatMap fun b =>
          let e := σ[b]!
          if (E t).contains e || (E t').contains e then []
          else
            let inAbove := above.any fun u => (E u).contains e
            let inUnpushed := (List.range (D+1)).any fun k =>
              unpushed s.stackVerts[k]! && (edgesBelow s (vertItem s.stackVerts[k]!)).contains e
            [bad "cand_gapfree" s!"{showT t} / {showT t'} gap e={e} above={inAbove} unpushed={inUnpushed}"] ++
            (if inAbove || inUnpushed then [] else [bad "cand_gapfree_gen" s!"{showT t} / {showT t'} gap e={e}"])
      | _, _ => []) ++ gaps (above ++ [t]) (t' :: rest)
    | _ => []
  out := gaps [] stk ++ out
  -- candidate: the top entry's edges reach the edge about to be pushed (positions from its first
  -- piece edge up to n are its own edges)
  if o.cls.lowval d < d then
    match stk with
    | t :: _ =>
      match ((P t).map pos).min? with
      | some lo =>
        for b in List.range' lo (n - lo) do
          if !(E t).contains σ[b]! then out := bad "cand_top_reaches" s!"{showT t} e={σ[b]!}" :: out
      | none => pure ()
    | [] => pure ()
    -- generalized: from any entry's first piece edge up to n, every position is owned by an entry
    -- at or above it
    for j in List.range stk.length do
      let t' := stk[j]!
      let run := stk.take (j+1)
      match ((P t').map pos).min? with
      | some lo =>
        for b in List.range' lo (n - lo) do
          if !run.any (fun u => (E u).contains σ[b]!) then
            out := bad "cand_reach_all" s!"{showT t'} j={j} e={σ[b]!}" :: out
      | none => pure ()
  -- candidate: runs of consecutive entries are σ-convex (holes are edges of the run)
  for j in List.range stk.length do
    for j' in List.range' (j+1) (stk.length - j - 1) do
      let t := stk[j]!
      let t' := stk[j']!
      let run := (stk.drop j).take (j' - j + 1)
      match ((P t').map pos).max?, ((P t).map pos).min? with
      | some hi', some lo =>
        for b in List.range' (hi' + 1) (lo - hi' - 1) do
          if !run.any (fun u => (E u).contains σ[b]!) then
            out := bad "cand_runs" s!"{showT t} j={j} / {showT t'} j'={j'} e={σ[b]!}" :: out
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

def checkVertPast (seed : Nat) (σ : List Nat) (v n : Nat) (s : WalkState) : List V :=
  if (edgesBelow s (vertItem v)).all (fun e => σ.idxOf e < n) then []
  else [⟨seed, s.ternarize, v, 0, "vertex", "vertex_past", s!"n={n} edges={edgesBelow s (vertItem v)} σ={σ}"⟩]

def checkClose (seed : Nat) (s : WalkState) : List V := Id.run do
  let ty := fun i => s.items[i]!.type
  let ch := fun i => s.items[i]!.ch
  let vs := fun i => s.items[i]!.vs
  let live := fun i => s.tstack.any (fun t => (spanItems t).contains i) || s.items.any (fun it => it.ch.contains i)
  let spr := fun i => [NodeType.S, .P, .R].contains (ty i)
  let pair := fun (p q : Nat × Nat) => p == q || p == (q.2, q.1)
  let two := fun i => match vs i with | (some a, some b) => a != b | _ => false
  let isVs := fun i v => (vs i).1 == some v || (vs i).2 == some v
  let below := ((List.range s.items.size).map (edgesBelow s)).toArray
  let attachments := below.map (atts s)
  let allBelow := fun i v => (List.range s.g.ne).all fun e => !inc s e v || below[i]!.contains e
  let inner := fun i v => touches s (List.range s.g.ne) v && allBelow i v
  let att := fun i v => attachments[i]!.contains v
  let mut out := []
  for i in List.range s.items.size do
    if live i then
      let bad := fun k => (⟨seed, s.ternarize, 0, 0, "close", "close_" ++ k,
        s!"i={i} type={repr (ty i)} vs={vs i} ch={ch i} stack={s.tstack.map showT}"⟩ : V)
      if spr i then
        if !two i then out := bad "vs_ne" :: out
        for v in List.range s.g.nv do
          if isVs i v && !att i v then out := bad "vs_att" :: out
          if (ch i).contains (vertItem v) != (inner i v && (ch i).all (fun c => !allBelow c v)) then
            out := bad "interior" :: out
        for c in ch i do
          if ty c != .V && !two c then out := bad "child_two" :: out
      for c in ch i do
        if (ty c == .I || ty c == .O) && ty i != .Q then out := bad "io_parent" :: out
        if ty i == .V && (ch c).isEmpty then out := bad "q_under_v" :: out
      let ve := Items.virtualEdges s.items i
      let xs := ((ch i).filter (fun c => ty c == .V)).map (· - 1)
      if ty i == .P then
        let ok := ve.length ≥ 2 && xs.isEmpty &&
          (match vs i with | (some u, some v) => ve.all (fun q => pair q (u, v)) | _ => false)
        if !ok then out := bad "p_shape" :: out
      if ty i == .S then
        let ok := xs.length ≥ 1 &&
          (match vs i with | (some u, some v) => ve == List.zip (u :: xs) (xs ++ [v]) | _ => false)
        if !ok then out := bad "s_order" :: out
      if ty i == .R then
        let keys := ve.map fun q => min q.1 q.2 + s.g.nv * max q.1 q.2
        let ok := xs.length ≥ 2 && ve.length ≥ 5 && keys.Nodup &&
          (match vs i with | (some u, some v) => ve.all (fun q => !pair q (u, v)) | _ => false)
        if !ok then out := bad "r_shape" :: out
      if ty i == .Q then
        let e := i - 1 - s.g.nv
        let edge := s.g.edges[e]!
        for v in List.range s.g.nv do
          if att i v && !isVs i v then out := bad "q_att_vs" :: out
          if (ch i).isEmpty && isVs i v && !att i v then out := bad "vs_att_q" :: out
        if (ch i).isEmpty then
          let ok := two i && (match vs i with | (some a, some b) => pair (a, b) edge | _ => false)
          if !ok then out := bad "q_leaf" :: out
        else
          let ok := match vs i, ch i with
            | (some u, none), [c] => edge.1 == edge.2 && inc s e u && ty c != .F && ty c != .V &&
                (ty c != .Q || (ch c).isEmpty) && vs c == (some u, none)
            | (some u, none), [c, w] => edge.1 != edge.2 && 1 ≤ w && w ≤ s.g.nv && inc s e u &&
                ty c != .F && ty c != .V && (ty c != .Q || (ch c).isEmpty) && pair (u, w - 1) edge &&
                (match vs c with | (some a, some b) => pair (a, b) (u, w - 1) | _ => false)
            | _, _ => false
          if !ok then out := bad "q_root" :: out
  return out

mutual
partial def iTree (seed : Nat) (σ : List Nat) (t : DfsTree) (d : Nat) (s : WalkState) : WalkState × List V :=
  match t with
  | .node v outs =>
    let n := σ.idxOf t.edgePostorder.head!
    let s := { s with stackVerts := s.stackVerts.set! d v }
    let (hv, s, vs) := iOuts seed σ v d outs false s
    let vs := vs ++ if hv then [] else checkVertPast seed σ v (if t.edgePostorder.isEmpty then σ.length else n + t.edgePostorder.length) s
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
  let vp := if hv then [] else checkVertPast seed σ v (σ.idxOf o.block.head!) s
  let lowval := o.cls.lowval d
  let s := ((do let lowDir ← stackDir lowval; setStackDir d (if lowval ≥ d then false else !lowDir) : WalkM Unit).run s).2
  let (hv, s) := if !hv && lowval < d && o.cls.isType1 then (true, ((pushVertTstack v d).run s).2) else (hv, s)
  let orig := s.tstack.length
  let (s, vs) := match o with
    | .tree _ _ child => iTree seed σ child (d+1) { s with firstOccurrence := s.firstOccurrence.set! d s.g.ne }
    | .back .. => (s, [])
  let vs := vp ++ vs ++ check seed σ v d o orig s ++ checkAdj seed σ v d o orig hv s ++ checkClose seed s
  let vs := vs ++ if hv then [] else checkVertPast seed σ v (σ.idxOf o.e) s
  let (hv', s) := (finishEdge v d o orig hv).run s
  (hv', s, vs ++ checkClose seed s)
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
   ⟨3, #[(0,1),(1,2),(2,0)]⟩, ⟨2, #[(0,1),(1,0)]⟩, ⟨4, #[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]⟩,
   ⟨4, #[(0,1),(1,2),(1,0),(2,3),(2,0)]⟩]

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
#eval summarize 0 401
