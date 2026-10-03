import WalkInvCheck.Common
/-! Verbatim copy of the executable mirrors of `checks/RangesInvCheck.lean` (`RangesInv`, `FinishAdj`,
`CloseInv`/`CloseAt`, `CloseCtx`, the P / vertex / loop-1 / `RangesInv`-iterate sites, `OwnedD`,
`VertCover`). Keep in sync with that file. -/
namespace WalkInvCheck.Ranges
open Spqr WalkM WalkState

theorem closeFacts_needs_graph_wf :
    let g : Graph := ⟨0, #[(0, 0)]⟩
    ¬ Items.CloseFacts g (g.walk false (g.dfsForest [] [])).items := by
  intro g h
  obtain ⟨a, b, hvs, -, -⟩ := h.q_leaf 0 (by decide) (by decide)
  have hn : Items.vs (g.walk false (g.dfsForest [] [])).items (edgeItem g 0) = (none, none) := by decide
  rw [hn] at hvs
  cases hvs

#print axioms closeFacts_needs_graph_wf


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

/-- The block boundaries of one `finishEdge` at which `RangesCloseSites.lean` asserts `CloseInv`,
labelled by the theorem (`closeEars_closeAt`, `loop1Body_closeAt`, `closeVertTail_closeAt`,
`finishP_closeAt`, `finishBack_closeAt`, `finishBoundary_closeAt`) or frame block responsible. -/
def closeSites (curV d : Nat) (o : DfsOut) (orig : Nat) (hv : Bool) (s : WalkState) :
    List (String × WalkState) := Id.run do
  let lv := o.cls.lowval d
  let fin := after (finishEdge curV d o orig hv) s
  if d ≤ lv then return [("finishBoundary_closeAt", fin)]
  let mut sites := []
  let mut rest := feBack curV lv d o s
  if o.cls.isTree then
    let edgeDir := s.stackDir[d]!
    let s₁ := ceS₁ o.dest d o.e (feS₀ d o s)
    sites := [("feS₀", feS₀ d o s), ("closeEars_closeAt", s₁)]
    for st in mergeSites (loop1Cond d) (loop1Body d edgeDir) s₁ do
      sites := sites ++ [("loop1Body_closeAt", after (loop1Body d edgeDir) st)]
    sites := sites ++ [("mergeLate", feS₂ d o s)]
    if hv then
      sites := sites ++
        [("closeVert'.pre", cvS₅ curV edgeDir o.cls.isType1 orig (feSingle d o s) (feS₂ d o s)),
         ("closeVertTail_closeAt", feS₃ curV d o orig s)]
      rest := feS₃ curV d o orig s
    else
      rest := feS₂ d o s
  else
    sites := [("finishBack_closeAt", rest)]
  return sites ++ [("finishP_closeAt", after (finishP curV lv o.cls.isType1) rest), ("finishTail", fin)]

/-- The boundary fields of `CloseCtx` that `finishBoundary_closeAt` consumes, evaluated at the pre-state of
every block-boundary `finishEdge`: `dest_lt`/`bd_loop` (proved, `dsOut`) and the admissions
`closeCtx_bd_vert`/`closeCtx_bd_node`. -/
def checkCtx (seed : Nat) (curV d : Nat) (o : DfsOut) (s : WalkState) : List V := Id.run do
  let lv := o.cls.lowval d
  if !(d ≤ lv) then return []
  let ty : Nat → NodeType := fun i => s.items[i]!.type
  let ch : Nat → List ItemId := fun i => s.items[i]!.ch
  let vs : Nat → Option Nat × Option Nat := fun i => s.items[i]!.vs
  let bad := fun th k => (⟨seed, s.ternarize, curV, d, th, "ctx_" ++ k,
    s!"o.e={o.e} dest={o.dest} lv={lv} stack={s.tstack.map showT}"⟩ : V)
  let mut out := []
  if !(o.dest < s.g.nv) then out := bad "dsOut" "dest_lt" :: out
  if !o.cls.isTree && o.dest != curV then out := bad "dsOut" "bd_loop" :: out
  if o.cls.isTree then
    match (if lv == d + 1 then s.tstack.head? else s.tstack.tail.head?) with
    | some t => if t.spans.2 != [vertItem o.dest] then out := bad "closeCtx_bd_vert" "bd_vert" :: out
    | none => out := bad "closeCtx_bd_vert" "bd_vert_none" :: out
    if lv != d + 1 then
      match s.tstack.head? with
      | some b =>
        let ok := match b.spans.1 with
          | [c] => ty c != .F && ty c != .V && (ty c != .Q || (ch c).isEmpty) &&
              (match vs c with
                | (some a, some b') => (a, b') == (curV, o.dest) || (a, b') == (o.dest, curV)
                | _ => false)
          | _ => false
        if !ok then out := bad "closeCtx_bd_node" "bd_node" :: out
      | none => out := bad "closeCtx_bd_node" "bd_node_none" :: out
  return out

/-- The `PSite` fields of `CloseCtx.p_site`, evaluated at the state `finishP` runs from whenever
`condP` holds there. -/
def checkP (seed : Nat) (curV d : Nat) (o : DfsOut) (orig : Nat) (hv : Bool) (s : WalkState) :
    List V := Id.run do
  let lv := o.cls.lowval d
  if !(lv < d) then return []
  let r := if o.cls.isTree then (if hv then feS₃ curV d o orig s else feS₂ d o s)
    else feBack curV lv d o s
  if !(result (condP curV lv o.cls.isType1) r) then return []
  let ty : Nat → NodeType := fun i => r.items[i]!.type
  let ch : Nat → List ItemId := fun i => r.items[i]!.ch
  let vs : Nat → Option Nat × Option Nat := fun i => r.items[i]!.vs
  let u := r.stackVerts[lv]!
  let dir := r.stackDir[lv]!
  let top := r.tstack.take 2
  let bad := fun k => (⟨seed, r.ternarize, curV, d, "closeCtx_p_site", "psite_" ++ k,
    s!"o.e={o.e} lv={lv} u={u} stack={r.tstack.map showT}"⟩ : V)
  let mut out := []
  match r.tstack with
  | cur :: _ => if cur.vStart != curV || cur.topDepth != lv then out := bad "stack" :: out
  | [] => out := bad "stack" :: out
  if curV == u then out := bad "ne" :: out
  let pair := fun (p q : Nat × Nat) => p == q || p == (q.2, q.1)
  let kind := fun j => [NodeType.S, .P, .R, .Q].contains (ty j) && (ty j != .Q || (ch j).isEmpty)
  let allE := List.range r.g.ne
  for t in top do
    let E := entryEdges r t
    let interior := fun w => allE.all fun e => !inc r e w || E.contains e
    match getSide t.spans dir, getSide t.spans (!dir) with
    | [j], [] =>
      match vs j with
      | (some a, some b) => if !pair (a, b) (u, curV) then out := bad "single_vs" :: out
      | _ => out := bad "single_vs" :: out
    | _, _ => out := bad "single" :: out
    for w in List.range r.g.nv do
      if touches r E w && w != curV && w != u && !interior w then out := bad "att" :: out
    if !(touches r E curV && touches r E u) then out := bad "touch" :: out
    for j in spanItems t do
      for j' in j :: (if ty j == .P then ch j else []) do
        if !kind j' then out := bad "kinds" :: out
      let cnt := (r.tstack.map fun t' => (spanItems t').count j).sum
      if cnt != 1 || r.items.any (fun it => it.ch.contains j) then out := bad "once" :: out
  let Eall := top.flatMap (entryEdges r)
  for w in [curV, u] do
    if !(allE.any fun e => inc r e w && !Eall.contains e) then out := bad "pend" :: out
  return out

def checkVSite (bad : String → V) (curV x : Nat) (r : WalkState) : List V := Id.run do
  let ty : Nat → NodeType := fun i => r.items[i]!.type
  let vs : Nat → Option Nat × Option Nat := fun i => r.items[i]!.vs
  match r.tstack with
  | [] => return [bad "stack"]
  | t :: _ =>
    let mut out := []
    let dir := r.stackDir[t.topDepth]!
    let u := r.stackVerts[t.topDepth]!
    let cs := getSide t.spans dir
    let E := entryEdges r t
    let allE := List.range r.g.ne
    let interior := fun (F : List Nat) (w : Nat) => allE.all fun e => !inc r e w || F.contains e
    if t.vStart != curV then out := bad "vstart" :: out
    if !(getSide t.spans (!dir)).isEmpty then out := bad "side" :: out
    if !(x < r.items.size && 1 + r.g.nv + r.g.ne ≤ x && r.items.all (fun it => !it.ch.contains x)
        && r.tstack.all (fun t' => !(spanItems t').contains x)) then out := bad "free" :: out
    if curV == u then out := bad "ne" :: out
    for c in cs do
      if ![NodeType.S, .P, .R, .Q, .V].contains (ty c) then out := bad "kinds" :: out
      if ty c != .V then
        match vs c with
        | (some a, some b) => if a == b then out := bad "two" :: out
        | _ => out := bad "two" :: out
    for w in List.range r.g.nv do
      if touches r E w && w != curV && w != u && !interior E w then out := bad "att" :: out
      let inner := touches r E w && interior E w && cs.all fun c => !interior (edgesBelow r c) w
      if cs.contains (vertItem w) != inner then out := bad "inner" :: out
    if !(touches r E curV && touches r E u) then out := bad "touch" :: out
    for w in [curV, u] do
      if !(allE.any fun e => inc r e w && !E.contains e) then out := bad "pend" :: out
    let terms := setSides dir u curV
    let pair := fun (p q : Nat × Nat) => p == q || p == (q.2, q.1)
    let xs := (cs.filter (fun c => ty c == .V)).map (· - 1)
    let ve := (cs.filter (fun c => ty c != .V)).map fun c => ((vs c).1.getD 0, (vs c).2.getD 0)
    if ty x == .S then
      if !(xs.length ≥ 1 && ve == List.zip (terms.1 :: xs) (xs ++ [terms.2])) then
        out := bad "s_order" :: out
    if ty x == .R then
      let keys := ve.map fun q => min q.1 q.2 + r.g.nv * max q.1 q.2
      if !(xs.length ≥ 2 && ve.length ≥ 5 && keys.Nodup && ve.all fun q => !pair q terms) then
        out := bad "r_shape" :: out
    if ty x == .P then
      if !(ve.length ≥ 2 && xs.isEmpty && ve.all fun q => pair q terms) then
        out := bad "p_shape" :: out
    return out

def checkV (seed : Nat) (curV d : Nat) (o : DfsOut) (orig : Nat) (hv : Bool) (s : WalkState) :
    List V :=
  let lv := o.cls.lowval d
  if !(o.cls.isTree && lv < d && hv && o.cls.isType1) then [] else
  let x := ((maybeUnwrapNxt (if feSingle d o s then NodeType.S else .R)).run (feS₂ d o s)).1
  let r := cvS₅ curV s.stackDir[d]! true orig (feSingle d o s) (feS₂ d o s)
  checkVSite (fun k => ⟨seed, r.ternarize, curV, d, "closeCtx_v_site", "vsite_" ++ k,
    s!"o.e={o.e} x={x} stack={r.tstack.map showT}"⟩) curV x r

def checkL1 (seed : Nat) (curV d : Nat) (o : DfsOut) (s : WalkState) : List V := Id.run do
  if !(o.cls.isTree && o.cls.lowval d < d) then return []
  let dir := s.stackDir[d]!
  let mut out := []
  for st in mergeSites (loop1Cond d) (loop1Body d dir) (ceS₁ o.dest d o.e (feS₀ d o s)) do
    let s₁ := l1S₁ d dir st
    let x := result (maybeUnwrapNxt (l1Ty d dir st)) s₁
    let r := after mergeTstackTops (l1S₂ d dir st)
    let bad := fun k => (⟨seed, r.ternarize, curV, d, "closeCtx_l1_site", "l1site_" ++ k,
      s!"o.e={o.e} x={x} stack={r.tstack.map showT}"⟩ : V)
    if s₁.tstack.length < 2 then out := bad "two" :: out
    match r.tstack with
    | t :: _ => out := checkVSite bad t.vStart x r ++ out
    | [] => out := bad "stack" :: out
  return out

/-- The `PContent` fields of `CloseContent.p_site` (`RangesCloseContent.lean`). -/
def checkPContent (bad : String → V) (curV lv : Nat) (r : WalkState) : List V := Id.run do
  let ty : Nat → NodeType := fun i => r.items[i]!.type
  let ch : Nat → List ItemId := fun i => r.items[i]!.ch
  let vs : Nat → Option Nat × Option Nat := fun i => r.items[i]!.vs
  let u := r.stackVerts[lv]!
  let dir := r.stackDir[lv]!
  let top := r.tstack.take 2
  let mut out := []
  match r.tstack with
  | cur :: _ => if cur.vStart != curV || cur.topDepth != lv then out := bad "p_stack" :: out
  | [] => out := bad "p_stack" :: out
  let pair := fun (p q : Nat × Nat) => p == q || p == (q.2, q.1)
  let kind := fun j => [NodeType.S, .P, .R, .Q].contains (ty j) && (ty j != .Q || (ch j).isEmpty)
  let allE := List.range r.g.ne
  for t in top do
    let E := entryEdges r t
    let interior := fun w => allE.all fun e => !inc r e w || E.contains e
    match getSide t.spans dir, getSide t.spans (!dir) with
    | [j], [] =>
      match vs j with
      | (some a, some b) => if !pair (a, b) (u, curV) then out := bad "p_single" :: out
      | _ => out := bad "p_single" :: out
    | _, _ => out := bad "p_single" :: out
    for w in List.range r.g.nv do
      if touches r E w && w != curV && w != u && !interior w then out := bad "p_att" :: out
    if !(touches r E curV && touches r E u) then out := bad "p_touch" :: out
    for j in spanItems t do
      for j' in j :: (if ty j == .P then ch j else []) do
        if !kind j' then out := bad "p_kinds" :: out
      let cnt := (r.tstack.map fun t' => (spanItems t').count j).sum
      if cnt != 1 || r.items.any (fun it => it.ch.contains j) then out := bad "p_once" :: out
  return out

/-- The `VContent` fields of `CloseContent.v_site`/`l1_site` for the top entry `t` of `r`
(`pre` prefixes the kind: `v_` or `l1_`). -/
def checkVContent (bad : String → V) (pre : String) (curV x : Nat) (r : WalkState) : List V := Id.run do
  let ty : Nat → NodeType := fun i => r.items[i]!.type
  let vs : Nat → Option Nat × Option Nat := fun i => r.items[i]!.vs
  match r.tstack with
  | [] => return [bad (pre ++ "stack")]
  | t :: _ =>
    let mut out := []
    let dir := r.stackDir[t.topDepth]!
    let u := r.stackVerts[t.topDepth]!
    let cs := getSide t.spans dir
    let E := entryEdges r t
    let allE := List.range r.g.ne
    let interior := fun (F : List Nat) (w : Nat) => allE.all fun e => !inc r e w || F.contains e
    if t.vStart != curV then out := bad (pre ++ "vstart") :: out
    if !(getSide t.spans (!dir)).isEmpty then out := bad (pre ++ "side") :: out
    if !(x < r.items.size && 1 + r.g.nv + r.g.ne ≤ x && r.items.all (fun it => !it.ch.contains x)
        && r.tstack.all (fun t' => !(spanItems t').contains x)) then out := bad (pre ++ "free") :: out
    for c in cs do
      if ![NodeType.S, .P, .R, .Q, .V].contains (ty c) then out := bad (pre ++ "kinds") :: out
      if ty c != .V then
        match vs c with
        | (some a, some b) => if a == b then out := bad (pre ++ "two") :: out
        | _ => out := bad (pre ++ "two") :: out
    if !(touches r E curV && touches r E u) then out := bad (pre ++ "touch") :: out
    for w in List.range r.g.nv do
      let inner := touches r E w && interior E w && cs.all fun c => !interior (edgesBelow r c) w
      if cs.contains (vertItem w) != inner then out := bad (pre ++ "inner") :: out
    let terms := setSides dir u curV
    let pair := fun (p q : Nat × Nat) => p == q || p == (q.2, q.1)
    let xs := (cs.filter (fun c => ty c == .V)).map (· - 1)
    let ve := (cs.filter (fun c => ty c != .V)).map fun c => ((vs c).1.getD 0, (vs c).2.getD 0)
    if ty x == .S then
      if !(xs.length ≥ 1 && ve == List.zip (terms.1 :: xs) (xs ++ [terms.2])) then
        out := bad (pre ++ "s_order") :: out
    if ty x == .R then
      let keys := ve.map fun q => min q.1 q.2 + r.g.nv * max q.1 q.2
      if !(xs.length ≥ 2 && ve.length ≥ 5 && keys.Nodup && ve.all fun q => !pair q terms) then
        out := bad (pre ++ "r_shape") :: out
    if ty x == .P then
      if !(ve.length ≥ 2 && xs.isEmpty && ve.all fun q => pair q terms) then
        out := bad (pre ++ "p_shape") :: out
    return out

/-- Every field of `CloseContent curV d o orig hv s` (`RangesCloseContent.lean`), evaluated at its
site state and reported separately (`content_*`). -/
def checkContent (seed : Nat) (curV d : Nat) (o : DfsOut) (orig : Nat) (hv : Bool) (s : WalkState) :
    List V := Id.run do
  let lv := o.cls.lowval d
  let bad := fun (r : WalkState) k => (⟨seed, s.ternarize, curV, d, "CloseContent", "content_" ++ k,
    s!"o.e={o.e} dest={o.dest} lv={lv} stack={r.tstack.map showT}"⟩ : V)
  let mut out := []
  if o.cls.isTree && d ≤ lv then
    let ty : Nat → NodeType := fun i => s.items[i]!.type
    let ch : Nat → List ItemId := fun i => s.items[i]!.ch
    let vs : Nat → Option Nat × Option Nat := fun i => s.items[i]!.vs
    match (if lv == d + 1 then s.tstack.head? else s.tstack.tail.head?) with
    | some t => if t.spans.2 != [vertItem o.dest] then out := bad s "bd_vert" :: out
    | none => out := bad s "bd_vert_none" :: out
    if lv != d + 1 then
      match s.tstack.head? with
      | some b =>
        let ok := match b.spans.1 with
          | [c] => ty c != .F && ty c != .V && (ty c != .Q || (ch c).isEmpty) &&
              (match vs c with
                | (some a, some b') => (a, b') == (curV, o.dest) || (a, b') == (o.dest, curV)
                | _ => false)
          | _ => false
        if !ok then out := bad s "bd_node" :: out
      | none => out := bad s "bd_node_none" :: out
  if lv < d && o.cls.isType1 then
    let r := feRest curV d o orig hv s
    if result (condP curV lv true) r then out := checkPContent (bad r) curV lv r ++ out
  if o.cls.isTree && lv < d && hv && o.cls.isType1 then
    let x := ((maybeUnwrapNxt (if feSingle d o s then NodeType.S else .R)).run (feS₂ d o s)).1
    let r := cvS₅ curV s.stackDir[d]! true orig (feSingle d o s) (feS₂ d o s)
    out := checkVContent (bad r) "v_" curV x r ++ out
  if o.cls.isTree && lv < d then
    let dir := s.stackDir[d]!
    for st in mergeSites (loop1Cond d) (loop1Body d dir) (ceS₁ o.dest d o.e (feS₀ d o s)) do
      let s₁ := l1S₁ d dir st
      let x := result (maybeUnwrapNxt (l1Ty d dir st)) s₁
      let r := after mergeTstackTops (l1S₂ d dir st)
      match r.tstack with
      | t :: _ =>
        let u := r.stackVerts[t.topDepth]!
        if t.vStart == u then out := bad r "l1_ne" :: out
        let E := entryEdges r t
        for w in [t.vStart, u] do
          if !((List.range r.g.ne).any fun e => inc r e w && !E.contains e) then out := bad r "l1_pend" :: out
        out := checkVContent (bad r) "l1_" t.vStart x r ++ out
      | [] => out := bad r "l1_stack" :: out
  return out

/-- `rangesInv_l1Iter`/`rangesInv_feS₂`: the `RangesInv` clauses (`o.e` processed, `D = d + 1`) at every
loop-1 iterate of `closeEars` and at `feS₂` of a returning tree edge. -/
def checkRI (seed : Nat) (σ : List Nat) (curV d : Nat) (o : DfsOut) (s : WalkState) : List V := Id.run do
  if !(o.cls.isTree && o.cls.lowval d < d) then return []
  let pos := fun e => σ.idxOf e
  let n := pos o.e + 1
  let dir := s.stackDir[d]!
  let iters := mergeSites (loop1Cond d) (loop1Body d dir) (ceS₁ o.dest d o.e (feS₀ d o s))
  let sites := iters.map (fun st => ("rangesInv_l1Iter", st)) ++ [("rangesInv_feS₂", feS₂ d o s)]
  let mut out := []
  for (site, st) in sites do
    let bad := fun k info => (⟨seed, st.ternarize, curV, d, site, k,
      s!"o.e={o.e} {info} stack={st.tstack.map showT}"⟩ : V)
    let E := fun t => entryEdges st t
    let P := fun t => entryPiece st t
    let stk := st.tstack
    for t in stk do
      for e in E t do
        if !(pos e < n) then out := bad "processed" s!"{showT t} e={e}" :: out
    let rec ord : List TEntry → List V
      | [] => []
      | t :: rest =>
        (rest.flatMap fun t' => (P t).flatMap fun e => (E t').filterMap fun e' =>
          if pos e' < pos e then none else some (bad "ordered" s!"{showT t} e={e} below {showT t'} e'={e'}")) ++
        ord rest
    out := ord stk ++ out
    for t in stk do
      let ps := (P t).map pos
      match ps.min?, ps.max? with
      | some lo, some hi =>
        for b in List.range' lo (hi + 1 - lo) do
          if !(E t).contains σ[b]! then out := bad "convex" s!"{showT t} hole e={σ[b]!}" :: out
      | _, _ => pure ()
    for i in List.range st.items.size do
      let it := st.items[i]!
      if it.type != .F && it.type != .V then
        let ps := (pieceBelow st i).map pos
        let below := edgesBelow st i
        match ps.min?, ps.max? with
        | some lo, some hi =>
          for b in List.range' lo (hi + 1 - lo) do
            if !below.contains σ[b]! then out := bad "closed_convex" s!"i={i} hole e={σ[b]!}" :: out
        | _, _ => pure ()
  return out

def checkClose (seed : Nat) (site : String) (s : WalkState) : List V := Id.run do
  let ty := fun i => s.items[i]!.type
  let ch := fun i => s.items[i]!.ch
  let vs := fun i => s.items[i]!.vs
  let live := fun i => i == rootItem || (1 ≤ i && i ≤ s.g.nv) ||
    s.tstack.any (fun t => (spanItems t).contains i) || s.items.any (fun it => it.ch.contains i)
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
      let bad := fun k => (⟨seed, s.ternarize, 0, 0, site, "close_" ++ k,
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

/-- Depth-indexed ownership invariant `OwnedD` (`RangesOwned.lean`) at the path vertex of depth `d`
(`c.sts[k]`/`c.origs[k]` = schedule start / stack length at the entry of the path vertex of depth
`k`, `c.visited` = placed vertices); `n` = number of processed edges. Checked at every vertex entry
(`entry`), `finishEdge` pre- and post-state (`pre`, `post` with `n + 1`), P site (`P`) and
vertex end (`end`); `finishEdge_ownedD` is `pre` → `post`. -/
structure OwnCtx where
  sts : List Nat
  origs : List Nat
  visited : List Nat

def checkOwned (seed : Nat) (σ : List Nat) (c : OwnCtx) (n v d : Nat) (site : String)
    (s : WalkState) : List V := Id.run do
  let mut out := []
  let fail := fun (kind msg : String) => (⟨seed, s.ternarize, v, 0, site, kind, msg⟩ : V)
  let owned := s.tstack.flatMap (entryEdges s)
  let len := s.tstack.length
  if s.stackVerts[d]! != v then out := fail "own_sv" s!"d={d}" :: out
  for k in List.range (d + 1) do
    if len < c.origs[k]! then out := fail "own_len" s!"k={k} len={len} orig={c.origs[k]!}" :: out
  for b in List.range' c.sts[0]! (n - c.sts[0]!) do
    let e := σ[b]!
    if !owned.contains e &&
        !(List.range (d + 1)).any (fun k => (edgesBelow s (vertItem s.stackVerts[k]!)).contains e) then
      out := fail "own_cover" s!"b={b} e={e} n={n} stack={s.tstack.map showT}" :: out
  for k in List.range d do
    for e in edgesBelow s (vertItem s.stackVerts[k]!) do
      if !(σ.idxOf e < c.sts[k + 1]!) then out := fail "own_anc" s!"k={k} e={e}" :: out
  for e in edgesBelow s (vertItem v) do
    if !(σ.idxOf e < n) then out := fail "own_hi" s!"e={e}" :: out
  for k in List.range (d + 1) do
    for e in edgesBelow s (vertItem s.stackVerts[k]!) do
      if σ.idxOf e < c.sts[k]! then out := fail "own_lo" s!"k={k} e={e}" :: out
    for t in s.tstack.take (len - c.origs[k]!) do
      for e in entryEdges s t do
        if σ.idxOf e < c.sts[k]! then out := fail "own_new" s!"k={k} e={e} t={showT t}" :: out
    for t in s.tstack.drop (len - c.origs[k]!) do
      if t.vStart == s.stackVerts[k]! then out := fail "own_old" s!"k={k} t={showT t}" :: out
  for t in s.tstack do
    if !c.visited.contains t.vStart then out := fail "own_vis" s!"t={showT t}" :: out
  for w in List.range s.g.nv do
    if !c.visited.contains w && s.items[vertItem w]!.ch != [] then out := fail "own_fresh" s!"w={w}" :: out
  return out

def checkVCover (seed : Nat) (v : Nat) (site : String) (s : WalkState) : List V := Id.run do
  let mut out := []
  let owned := s.tstack.flatMap (entryEdges s)
  for e in edgesBelow s (vertItem v) do
    if !owned.contains e then
      out := (⟨seed, s.ternarize, v, 0, site, "own_vcover", s!"e={e} stack={s.tstack.map showT}"⟩ : V) :: out
  return out

end WalkInvCheck.Ranges
