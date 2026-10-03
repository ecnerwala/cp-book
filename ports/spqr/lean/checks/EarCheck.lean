import Spqr.Ear
import Spqr.ItemSpec
import Spqr.EarInv
/-!
# Empirical check of `WalkState.EarFinish` (`Spqr/EarInv.lean`)

Run with `lake env lean checks/EarCheck.lean` from `ports/spqr/lean` (not part of the library).

`iTree`/`iOuts`/`iOut` re-implement `walkTree`/`walkOuts`/`walkOut` around the library's
`finishEdge` (as in `checks/InvCheck.lean`) and, right before every `finishEdge curV d o
origTstack hasVert`, evaluate each field of `EarFinish curV d o hasVert sub base` with
`sub = tstack.take (len - origTstack)`, `base = tstack.drop (len - origTstack)` (the only split
`EarAt` allows), plus the ear-bottom anchor (`eb_*`: `EarFinish.bottom`) and the shape after loops 1–2 (`t_*`/`t1_*`/`t2_*`:
`EarFinish.loops`), recorded per chain bottom in `EB`.
-/
open Spqr WalkM
instance : Inhabited Spqr.DfsTree := ⟨.node 0 []⟩
instance : Inhabited Spqr.DfsOut := ⟨.back 0 0 .selfLoop⟩

partial def edgesBelow (s : WalkState) (i : ItemId) : List Nat :=
  (if 1 + s.g.nv ≤ i ∧ i < 1 + s.g.nv + s.g.ne then [i - 1 - s.g.nv] else []) ++
    (s.items[i]!.ch).flatMap (edgesBelow s)
def entryEdges (s : WalkState) (t : TEntry) : List Nat := (t.spans.1 ++ t.spans.2).flatMap (edgesBelow s)
def inc (s : WalkState) (e v : Nat) : Bool := s.g.edges[e]!.1 == v || s.g.edges[e]!.2 == v
def touches (s : WalkState) (E : List Nat) (v : Nat) : Bool := E.any (inc s · v)
def isParent (s : WalkState) (p c : ItemId) : Bool := (s.items[p]!.ch).contains c
def hasParent (s : WalkState) (c : ItemId) : Bool :=
  (List.range s.items.size).any fun p => isParent s p c
def subEdgesL (o : DfsOut) : List Nat :=
  match o with
  | .tree e _ child => e :: child.edges
  | .back e _ _ => [e]
def spanItems (t : TEntry) : List ItemId := t.spans.1 ++ t.spans.2
def onSide (t : TEntry) (dir : Bool) : Bool := getSide t.spans (!dir) == []
/-- `Graph.ConnEdges E` literally: every `.1` endpoint of `E` reaches every other through `E`. -/
def connB (s : WalkState) (E : List Nat) : Bool := Id.run do
  match E with
  | [] => return true
  | e₀ :: _ =>
    let mut seen : List Nat := [s.g.edges[e₀]!.1]
    let mut changed := true
    while changed do
      changed := false
      for e in E do
        let (a, b) := s.g.edges[e]!
        if seen.contains a && !seen.contains b then seen := b :: seen; changed := true
        if seen.contains b && !seen.contains a then seen := a :: seen; changed := true
    return E.all fun e => seen.contains s.g.edges[e]!.1
/-- `Graph.TwoAttached E a b` literally. -/
def twoAttB (s : WalkState) (E : List Nat) (a b : Nat) : Bool :=
  (List.range s.g.nv).all fun v =>
    !(touches s E v && (List.range s.g.ne).any (fun e => !E.contains e && inc s e v)) || v == a || v == b
def showT (t : TEntry) : String := s!"({t.vStart},{t.topDepth},{t.firstIdx},{t.spans})"
def sameT (a b : TEntry) : Bool := showT a == showT b

structure V where
  seed : Nat
  curV : Nat
  d : Nat
  o : String
  hasVert : Bool
  kind : String
  info : String
deriving Repr

def traceOn : Bool := false
def earCheck (seed curV d : Nat) (o : DfsOut) (orig : Nat) (hv : Bool) (s : WalkState) : List V := Id.run do
  let n := s.tstack.length
  let sub := s.tstack.take (n - orig)
  let base := s.tstack.drop (n - orig)
  let dir := s.stackDir[d]!
  let lowval := o.cls.lowval d
  let oStr := s!"{if o.cls.isTree then "tree" else "back"} e={o.e} dest={o.dest} lv={lowval} t1={o.cls.isType1}"
  let mut out : List V := []
  if traceOn then out := [⟨seed, curV, d, oStr, hv, "TRACE", s!"orig={orig} dir={s.stackDir.toList.take (d+2)} sv={s.stackVerts.toList.take (d+2)} stack={s.tstack.map showT}"⟩]
  let bad (k : String) (info : String) : V :=
    ⟨seed, curV, d, oStr, hv, k, s!"{info} | orig={orig} dir={s.stackDir.toList.take (d+2)} sv={s.stackVerts.toList.take (d+2)} nv={s.g.nv} stack={s.tstack.map showT}"⟩
  let E := fun t => entryEdges s t
  let SE := subEdgesL o
  if !o.cls.isTree && sub ≠ [] then out := bad "back_nil" (toString (sub.map showT)) :: out
  for t in base do
    if o.cls.isTree && t.vStart == o.dest then out := bad "base_bot" (showT t) :: out
    if (E t).any (SE.contains ·) then out := bad "base_disj" (showT t) :: out
  for t in sub do
    if (E t).any (fun e => !SE.contains e) then out := bad "sub_edges" (showT t) :: out
    if t.vStart == curV then out := bad "sub_bot" (showT t) :: out
  for e in SE do
    if e ≠ o.e && !sub.any (fun t => (E t).contains e) then out := bad "sub_cover" s!"e={e}" :: out
  for k in List.range (d+1) do
    for k' in List.range (d+1) do
      if k < k' && s.stackVerts[k]! == s.stackVerts[k']! then out := bad "path" s!"{k} {k'}" :: out
  let rec pw : List TEntry → Bool
    | [] => true
    | t :: rest => rest.all (fun t' => !(E t).any ((E t').contains ·)) && pw rest
  if !pw s.tstack then out := bad "disj" "" :: out
  let rec pws : List TEntry → Bool
    | [] => true
    | t :: rest => rest.all (fun t' => !(spanItems t).any ((spanItems t').contains ·)) && pws rest
  if !pws s.tstack then out := bad "span_disj" "" :: out
  let hi := sub.takeWhile fun t => d ≤ t.topDepth
  let ret := lowval < d
  if ret then
    for t in hi do
      if !onSide t dir then out := bad "loop1_side" (showT t) :: out
      if (E t) ≠ [] && t.topDepth ≤ d + 1 && !touches s (E t) s.stackVerts[t.topDepth]! then
        out := bad "loop1_touch" (showT t) :: out
  for t in s.tstack do
    if (E t) ≠ [] && !touches s (E t) t.vStart then out := bad "touch_bot" (showT t) :: out
    if (spanItems t).contains (vertItem curV) then
      if !hv then out := bad "vert_free.hv" (showT t) :: out
      if !base.any (sameT t) then out := bad "vert_free.base" (showT t) :: out
    if (spanItems t).contains (edgeItem s.g o.e) then out := bad "q_free" (showT t) :: out
  if hv && !base.any (fun t => t.topDepth ≤ d && (spanItems t).contains (vertItem curV)) then
    out := bad "vert" "" :: out
  if hasParent s (edgeItem s.g o.e) then out := bad "q_root" "" :: out
  if hasParent s (vertItem curV) then out := bad "v_root" "" :: out
  if ret && o.cls.isType1 then
    for t in base do
      if t.vStart == curV && t.topDepth == lowval then
        if !onSide t s.stackDir[lowval]! then out := bad "p_entry.side" (showT t) :: out
        match spanItems t with
        | [i] => if hasParent s i then out := bad "p_entry.root" (showT t) :: out
        | _ => out := bad "p_entry.single" (showT t) :: out
        for v in List.range s.g.nv do
          if touches s (E t) v && v ≠ curV && v ≠ s.stackVerts[lowval]! &&
              (List.range s.g.ne).any (fun e => inc s e v && !(E t).contains e) then
            out := bad "p_entry.att" s!"{showT t} v={v}" :: out
  if d ≤ lowval then
    for t in sub do
      for u in base do
        for v in List.range s.g.nv do
          if touches s (E t) v && touches s (E u) v && v ≠ curV then
            out := bad "boundary" s!"{showT t} {showT u} v={v}" :: out
  if hv && d ≤ lowval then out := bad "bd_noVert" "" :: out
  if o.cls.isTree then
    for e in List.range s.g.ne do
      if inc s e o.dest && !SE.contains e then out := bad "dest_edges" s!"e={e}" :: out
    for u in base do
      if touches s (E u) o.dest then out := bad "base_touch" (showT u) :: out
    if d ≤ lowval then
      if lowval == d + 1 then
        match sub with
        | [t] => if t.vStart ≠ o.dest || t.topDepth ≠ d + 1 then out := bad "bd_bridge" (showT t) :: out
        | _ => out := bad "bd_bridge" (toString (sub.map showT)) :: out
      else
        match sub with
        | [t₁, t₂] =>
          if t₁.vStart ≠ o.dest || t₁.topDepth ≠ lowval || t₂.vStart ≠ o.dest || t₂.topDepth ≠ d + 1 then
            out := bad "bd_comp" (toString (sub.map showT)) :: out
        | _ => out := bad "bd_comp" (toString (sub.map showT)) :: out
      for u in s.tstack.tail do
        if touches s (E u) curV && u.vStart ≠ curV && d < u.topDepth then out := bad "bd_term" (showT u) :: out
      if lowval == d + 1 then
        match s.tstack.head? with
        | some t => if t.spans.1 ≠ [] then out := bad "bd_side.bridge" (showT t) :: out
        | none => pure ()
      else
        match s.tstack with
        | b :: t :: _ =>
          if b.spans.2 ≠ [] then out := bad "bd_side.back" (showT b) :: out
          if t.spans.1 ≠ [] then out := bad "bd_side.vert" (showT t) :: out
        | _ => pure ()
    let fe := WalkState.after (finishEdge curV d o orig hv) s
    if fe.g.ne ≠ s.g.ne || fe.g.edges ≠ s.g.edges || fe.stackVerts ≠ s.stackVerts then out := bad "lower_frame" "" :: out
    for i in List.range fe.tstack.length do
      let t := fe.tstack[i]!
      let above := fe.tstack.take i
      let Et := entryEdges fe t
      if touches fe Et o.dest && !(List.range s.g.ne).all (fun e => !inc fe e o.dest || Et.contains e) &&
          t.vStart ≠ o.dest && !above.any (fun t' => t'.vStart == o.dest) then
        out := bad "lower" s!"{showT t} final={fe.tstack.map showT}" :: out
  return out

/-! ### Loop 1: `Loop1BodyOk` (D = d+1) at every iterate -/

def interiorB (s : WalkState) (Es : List Nat) (v : Nat) : Bool :=
  (List.range s.g.ne).all fun e => !inc s e v || Es.contains e
/-! ### `EarFinish.loop1` (`Loop1Spec`) literally -/
def l1CloseB (s : WalkState) (d v : Nat) (E : List Nat) (t : TEntry) (withMid : Bool) : List String :=
  let Et := entryEdges s t
  let U := E ++ Et
  (if Et == [] || (List.range s.g.nv).any (fun x => touches s E x && touches s Et x) then [] else ["share"]) ++
  (if v == t.vStart || (List.range (d+2)).any (fun k => d ≤ k && v == s.stackVerts[k]!) || interiorB s U v then [] else ["bottom"]) ++
  (if !withMid then [] else
    let x := s.stackVerts[d+1]!
    if x == t.vStart || interiorB s U x || !touches s U x then [] else ["mid"])
def unwrapB (s : WalkState) (ty : NodeType) (t : TEntry) : List String :=
  [false, true].flatMap fun dir =>
    let h := (getSide t.spans dir).headD 0
    if Items.type s.items h ≠ ty then [] else
    (if getSide t.spans dir ≠ [h] then ["single"] else []) ++
    (if hasParent s h then ["root"] else []) ++
    (if (Items.ch s.items h).any (fun c => s.tstack.any fun u => (spanItems u).contains c) then ["child"] else [])
/-- Walk the reached splits of `hi`, checking `Loop1Spec`'s per-split facts. -/
def specWalk (s : WalkState) (d : Nat) (o : DfsOut) (bad : String → String → V) :
    Nat → List TEntry → List TEntry → List V
  | 0, _, _ => []
  | _ + 1, _, [] => []
  | fuel + 1, done, t :: rest' =>
    let v := (done.getLast?.map TEntry.vStart).getD o.dest
    let E := o.e :: done.flatMap (entryEdges s)
    if t.topDepth == d then
      (l1CloseB s d v E t true).map (fun m => bad s!"close_{m}" (showT t)) ++
      (if t.vStart == v then (unwrapB s .P t).map (fun m => bad s!"unwrapP_{m}" (showT t)) else []) ++
      specWalk s d o bad fuel (done ++ [t]) rest'
    else
      (l1CloseB s d v E t false).map (fun m => bad s!"merge_{m}" (showT t)) ++
      match rest' with
      | t' :: rest'' =>
        (l1CloseB s d t.vStart (E ++ entryEdges s t) t' true).map (fun m => bad s!"close2_{m}" s!"{showT t} {showT t'}") ++
        (unwrapB s .S t').map (fun m => bad s!"unwrapS_{m}" s!"{showT t} {showT t'}") ++
        specWalk s d o bad fuel (done ++ [t, t']) rest''
      | [] => [bad "series_last" (showT t)]
def specCheck (seed curV d : Nat) (o : DfsOut) (orig : Nat) (hv : Bool) (s : WalkState) : List V :=
  let lowval := o.cls.lowval d
  if !(o.cls.isTree && lowval < d) then [] else
  let oStr := s!"tree e={o.e} dest={o.dest} lv={lowval} t1={o.cls.isType1}"
  let n := s.tstack.length
  let sub := s.tstack.take (n - orig)
  let hi := sub.takeWhile fun t => d ≤ t.topDepth
  let bad (kind : String) (info : String) : V :=
    ⟨seed, curV, d, oStr, hv, s!"spec_{kind}", s!"{info} | orig={orig} sv={s.stackVerts.toList.take (d+2)} stack={s.tstack.map showT}"⟩
  specWalk s d o bad (hi.length + 1) [] hi

/-! ### `EarFinish.late` / `EarFinish.close` (`EarLate`, `EarClose`) and the path fields -/
def termB (s : WalkState) (D : Nat) (t : TEntry) (v : Nat) : Bool :=
  v == t.vStart || (List.range (D+1)).any fun k => t.topDepth ≤ k && v == s.stackVerts[k]!
/-- `MergeOk D s cur nxt` literally. -/
def mergeOkB (s : WalkState) (D : Nat) (cur nxt : TEntry) : List String :=
  let Ec := entryEdges s cur
  let En := entryEdges s nxt
  (if Ec == [] || En == [] || (List.range s.g.nv).any (fun x => touches s Ec x && touches s En x) then [] else ["share"]) ++
  (if termB s D (TEntry.mergeInto cur nxt) cur.vStart || interiorB s (Ec ++ En) cur.vStart then [] else ["bottom"])
def pairwiseDisjB (s : WalkState) : List TEntry → Bool
  | [] => true
  | t :: rest => rest.all (fun u => (entryEdges s t).all fun e => !(entryEdges s u).contains e) && pairwiseDisjB s rest
def pairwiseSpanDisjB : List TEntry → Bool
  | [] => true
  | t :: rest => rest.all (fun u => (spanItems t).all fun i => !(spanItems u).contains i) && pairwiseSpanDisjB rest
def sameEdges (a b : List Nat) : Bool := a.all b.contains && b.all a.contains
/-- `FoldSpec D s cur R` literally (`l2Cur` at every split). -/
def foldWalk (s : WalkState) (D : Nat) (pre : String) (bad : String → String → V) : TEntry → List TEntry → List V
  | _, [] => []
  | cur, t :: rest =>
    (mergeOkB s D cur t).map (fun m => bad s!"{pre}_{m}" s!"cur={showT cur} t={showT t}") ++
    foldWalk s D pre bad (TEntry.mergeInto cur t) rest
/-- `EarLate` at `feS₁`: the splits reached while the loop condition held. -/
def lateWalk (s : WalkState) (D fo : Nat) (bad : String → String → V) : TEntry → List TEntry → List V
  | _, [] => []
  | cur, t :: rest =>
    (mergeOkB s D cur t).map (fun m => bad s!"late_{m}" s!"cur={showT cur} t={showT t}") ++
    (if fo < t.firstIdx then lateWalk s D fo bad (TEntry.mergeInto cur t) rest else [])
def closeCheck (seed curV d : Nat) (o : DfsOut) (orig : Nat) (hv : Bool) (s : WalkState) : List V := Id.run do
  let lowval := o.cls.lowval d
  let oStr := s!"{if o.cls.isTree then "tree" else "back"} e={o.e} dest={o.dest} lv={lowval} t1={o.cls.isType1}"
  let n := s.tstack.length
  let base := s.tstack.drop (n - orig)
  let mut out : List V := []
  let bad (k : String) (info : String) : V :=
    ⟨seed, curV, d, oStr, hv, k, s!"{info} | orig={orig} dir={s.stackDir.toList.take (d+2)} sv={s.stackVerts.toList.take (d+2)} nv={s.g.nv} stack={s.tstack.map showT}"⟩
  -- the path fields
  if s.stackVerts[d]! ≠ curV then out := bad "sv_d" "" :: out
  if o.cls.isTree then
    if s.stackVerts[d+1]! ≠ o.dest then out := bad "sv_child" "" :: out
    if (List.range (d+1)).any (fun k => s.stackVerts[k]! == o.dest) then out := bad "path_child" "" :: out
  if lowval < d && s.stackDir[d]! ≠ !s.stackDir[lowval]! then out := bad "dir_d" "" :: out
  if !(o.cls.isTree && lowval < d) then return out
  -- `EarLate` at `feS₁`
  let s₁ := WalkState.feS₁ d o s
  let bad₁ (k : String) (info : String) : V :=
    ⟨seed, curV, d, oStr, hv, k, s!"{info} | orig={orig} fo={s₁.firstOccurrence[d]!} sv={s.stackVerts.toList.take (d+2)} nv={s.g.nv} stack₁={s₁.tstack.map showT} stack={s.tstack.map showT}"⟩
  match s₁.tstack with
  | [] => out := bad₁ "late_empty" "" :: out
  | c₀ :: R =>
    if !pairwiseDisjB s₁ s₁.tstack then out := bad₁ "late_disj" "" :: out
    let fo := s₁.firstOccurrence[d]!
    if fo < c₀.firstIdx then out := lateWalk s₁ (d+1) fo bad₁ c₀ R ++ out
    match s₁.tstack.getLast? with
    | some e => if fo < e.firstIdx then out := bad₁ "late_fo" (showT e) :: out
    | none => pure ()
  -- `EarClose` at `feS₂`
  let s₂ := WalkState.feS₂ d o s
  let bad₂ (k : String) (info : String) : V :=
    ⟨seed, curV, d, oStr, hv, k, s!"{info} | orig={orig} sv={s.stackVerts.toList.take (d+2)} nv={s.g.nv} stack₂={s₂.tstack.map showT} stack={s.tstack.map showT}"⟩
  let n₂ := s₂.tstack.length
  let sub₂ := s₂.tstack.take (n₂ - orig)
  let base₂ := s₂.tstack.drop (n₂ - orig)
  let E₂ := fun t => entryEdges s₂ t
  if !(base₂.map showT == base.map showT) then out := bad₂ "close_base" "" :: out
  if s₂.g.edges ≠ s.g.edges || s₂.g.nv ≠ s.g.nv then out := bad₂ "close_g" "" :: out
  if s₂.stackVerts ≠ s.stackVerts then out := bad₂ "close_sv" "" :: out
  if s₂.stackDir[lowval]! ≠ s.stackDir[lowval]! then out := bad₂ "close_dir_l" "" :: out
  for t in base do
    if !sameEdges (E₂ t) (entryEdges s t) then out := bad₂ "close_base_edges" (showT t) :: out
    for i in spanItems t do
      if !hasParent s i && hasParent s₂ i then out := bad₂ "close_base_root" (showT t) :: out
  if !pairwiseDisjB s₂ s₂.tstack then out := bad₂ "close_disj" "" :: out
  if !pairwiseSpanDisjB s₂.tstack then out := bad₂ "close_span_disj" "" :: out
  let subE := subEdgesL o
  for t in sub₂ do
    if !(E₂ t).all subE.contains then out := bad₂ "close_sub_edges" (showT t) :: out
  for e in subE do
    if !sub₂.any (fun t => (E₂ t).contains e) then out := bad₂ "close_sub_cover" s!"e={e}" :: out
  match sub₂ with
  | c :: rest =>
    match rest.reverse with
    | vy :: py :: midr =>
      let mid := midr.reverse
      if c.topDepth < lowval || d < c.topDepth then out := bad₂ "close_c_top" (showT c) :: out
      if vy.topDepth ≤ d then out := bad₂ "close_vy_top" (showT vy) :: out
      if !(E₂ c).contains o.e then out := bad₂ "close_c_edge" (showT c) :: out
      if (List.range s.g.ne).any (fun e => inc s e vy.vStart && !subE.contains e) then out := bad₂ "close_y_edges" (showT vy) :: out
      match spanItems py with
      | [i] => if hasParent s₂ i then out := bad₂ "close_py_root" (showT py) :: out
      | _ => out := bad₂ "close_py_single" (showT py) :: out
      if !onSide py s.stackDir[lowval]! then out := bad₂ "close_py_side" (showT py) :: out
      if py.vStart ≠ vy.vStart || py.topDepth ≠ lowval then out := bad₂ "close_py" (showT py) :: out
      out := foldWalk s₂ (d+1) "fold" bad₂ c (mid ++ [py, vy]) ++ out
      if o.cls.isType1 then
        if mid ≠ [] then out := bad₂ "close_t1_mid" "" :: out
        let U := E₂ c ++ E₂ py ++ E₂ vy
        for v in List.range s.g.nv do
          if touches s₂ U v && v ≠ curV && v ≠ s.stackVerts[lowval]! && !interiorB s₂ U v then
            out := bad₂ "close_t1_touch" s!"v={v}" :: out
    | _ => out := bad₂ "close_short" "" :: out
  | [] => out := bad₂ "close_short" "" :: out
  if !hv then
    let EV := edgesBelow s₂ (vertItem curV)
    for t in s₂.tstack do
      if (E₂ t).any EV.contains then out := bad₂ "close_vert_disj" (showT t) :: out
    if EV ≠ [] && !touches s₂ EV curV then out := bad₂ "close_vert_touch" "" :: out
  return out

/-- Ear bottoms: `(y, (y,l) piece, V y)` recorded when the walk of `y` ends, if `y`'s own entries
begin (bottom-up) with its vertex entry. -/
abbrev EB := List (Nat × TEntry × TEntry)

/-- Checks of the ear-bottom claim at a `finishEdge`: for a returning tree edge the child's entries
end with `[(y, lowval) piece, V y]` as recorded for `y`; every recorded pair still on the stack is
untouched (and adjacent). -/
def ebCheck (seed curV d : Nat) (o : DfsOut) (orig : Nat) (hv : Bool) (eb : EB) (s : WalkState) : List V := Id.run do
  let n := s.tstack.length
  let sub := s.tstack.take (n - orig)
  let lowval := o.cls.lowval d
  let oStr := s!"{if o.cls.isTree then "tree" else "back"} e={o.e} dest={o.dest} lv={lowval} t1={o.cls.isType1}"
  let bad (k : String) (info : String) : V :=
    ⟨seed, curV, d, oStr, hv, k, s!"{info} | orig={orig} sv={s.stackVerts.toList.take (d+2)} nv={s.g.nv} stack={s.tstack.map showT}"⟩
  let mut out : List V := []
  if o.cls.isTree && lowval < d then
    match sub.reverse with
    | vy :: py :: _ =>
      match eb.find? (·.1 == vy.vStart) with
      | some (_, p, v) =>
        if !(sameT vy v && sameT py p) then out := bad "eb_sub_bottom" s!"rec=({showT p},{showT v})" :: out
        if py.topDepth ≠ lowval then out := bad "eb_lowval" (showT py) :: out
        if !onSide py s.stackDir[lowval]! then out := bad "eb_py_side" (showT py) :: out
        if (spanItems py).any (hasParent s) then out := bad "eb_p_root" (showT py) :: out
        if d ≥ vy.topDepth then out := bad "eb_vy_top" (showT vy) :: out
        if !touches s (entryEdges s py) vy.vStart || !touches s (entryEdges s py) s.stackVerts[lowval]! then out := bad "eb_py_touch" (showT py) :: out
        if py.vStart ≠ vy.vStart then out := bad "eb_vstart" (showT py) :: out
        if spanItems vy ≠ [vertItem vy.vStart] then out := bad "eb_vspans" (showT vy) :: out
        match spanItems py with
        | [_] => pure ()
        | _ => out := bad "eb_psingle" (showT py) :: out
      | none => out := bad "eb_unrecorded" (showT vy ++ showT py) :: out
    | _ => out := bad "eb_short" (toString (sub.map showT)) :: out
  -- type-1 three-entry shape after loop 1 / loop 2 (`feS₂`), and the type-2 shape
  if o.cls.isTree && lowval < d then
    let s₂ := WalkState.feS₂ d o s
    let n₂ := s₂.tstack.length
    let sub₂ := s₂.tstack.take (n₂ - orig)
    let E₂ := fun t => entryEdges s₂ t
    match sub₂.reverse, sub.reverse with
    | vy₂ :: py₂ :: _, vy :: py :: _ =>
      if !(sameT vy₂ vy && sameT py₂ py) then out := bad "t_bottom_same" s!"{showT py₂} {showT vy₂}" :: out
    | _, _ => pure ()
    if o.cls.isType1 then
      match sub₂ with
      | [c, py, vy] =>
        if py.topDepth ≠ lowval || py.vStart ≠ vy.vStart then out := bad "t1_py" (showT py) :: out
        if spanItems vy ≠ [vertItem vy.vStart] then out := bad "t1_vy" (showT vy) :: out
        match spanItems py with
        | [_] => pure ()
        | _ => out := bad "t1_py_single" (showT py) :: out
        for v in List.range s₂.g.nv do
          if touches s₂ (E₂ vy) v && v ≠ vy.vStart && !((List.range s₂.g.ne).all fun e => !inc s₂ e v || (E₂ vy).contains e) then
            out := bad "t1_vy_bd" s!"v={v} {showT vy}" :: out
        if c.topDepth < lowval || d < c.topDepth then out := bad "t1_cur_range" (showT c) :: out
        -- cur ∪ py ∪ vy = subEdges
        for e in subEdgesL o do
          if !((E₂ c).contains e || (E₂ py).contains e || (E₂ vy).contains e) then out := bad "t1_cover" s!"e={e}" :: out
        -- py's edges all end at stackVerts[lowval] or are interior; cur touches curV and o.dest
        if !touches s₂ (E₂ c) curV then out := bad "t1_cur_touch_cur" (showT c) :: out
        if !touches s₂ (E₂ c) o.dest then out := bad "t1_cur_touch_dest" (showT c) :: out
        if !touches s₂ (E₂ py) s₂.stackVerts[lowval]! then out := bad "t1_py_touch_top" (showT py) :: out
        if !touches s₂ (E₂ py) py.vStart then out := bad "t1_py_touch_bot" (showT py) :: out
        if !onSide py s.stackDir[lowval]! then out := bad "t1_py_side" (showT py) :: out
      | _ => out := bad "t1_shape" s!"hv={hv} {sub₂.map showT}" :: out
    else
      match sub₂.reverse with
      | vy :: py :: _ =>
        if py.topDepth ≠ lowval || py.vStart ≠ vy.vStart then out := bad "t2_py" (showT py) :: out
        if spanItems vy ≠ [vertItem vy.vStart] then out := bad "t2_vy" (showT vy) :: out
        match sub₂ with
        | c :: _ => if c.topDepth < lowval || d < c.topDepth then out := bad "t2_cur_range" (showT c) :: out
        | _ => pure ()
      | _ => out := bad "t2_shape" s!"hv={hv} {sub₂.map showT}" :: out
  -- untouched: every recorded `V y` still present as its own entry has the recorded piece above it
  let rec scan : List TEntry → List V
    | [] => []
    | [_] => []
    | a :: b :: rest =>
      (match eb.find? (fun (y, _, _) => spanItems b == [vertItem y]) with
       | some (_, p, v) => if sameT b v && sameT a p then [] else [bad "eb_touched" s!"{showT a} {showT b} rec=({showT p},{showT v})"]
       | none => []) ++ scan (b :: rest)
  out := scan s.tstack ++ out
  for (y, _, v) in eb do
    match s.tstack.getLast? with
    | some b => if spanItems b == [vertItem y] && !sameT b v then out := bad "eb_touched_last" (showT b) :: out
    | none => pure ()
  return out

/-! ### `EarCtx` (PROOF.md §4.2b): the between-edges invariant of `walkOuts v d`, checked at the
start of `walkOuts` and after every `walkOut` return (`done` = finished outs, `rest` = the others;
`base₀`/`bE₀`/`sv₀` = the tstack, its entries' edge sets and `stackVerts[0..d]` at entry to `v`). -/
def ctxCheck (seed v d : Nat) (done : List (DfsOut × Bool)) (rest : List DfsOut) (hv : Bool) (base₀ : List TEntry)
    (bE₀ : List (List Nat)) (sv₀ : List Nat) (sd₀ : List Bool) (s : WalkState) (endPush : Bool := false) :
    List V := Id.run do
  let n := s.tstack.length
  let top := s.tstack.take (n - base₀.length)
  let base := s.tstack.drop (n - base₀.length)
  let oStr := s!"ctx done={done.length} rest={rest.length}"
  let mut out : List V := []
  let bad (k : String) (info : String) : V :=
    ⟨seed, v, d, oStr, hv, s!"ctx_{k}", s!"{info} | sv={s.stackVerts.toList.take (d+2)} nv={s.g.nv} top={top.map showT} base={base.map showT}"⟩
  let E := fun t => entryEdges s t
  -- (i) frame: `base` untouched, the path
  if base.length ≠ base₀.length || !(List.zip base base₀).all (fun (a, b) => sameT a b) then out := bad "base" "" :: out
  if !(List.zip base bE₀).all (fun (a, e₀) => sameEdges (E a) e₀) then out := bad "base_edges" "" :: out
  if s.stackVerts.toList.take (d+1) ≠ sv₀ then out := bad "sv" "" :: out
  if s.stackDir.toList.take d ≠ sd₀ then out := bad "sd" "" :: out
  if s.stackVerts[d]! ≠ v then out := bad "sv_d" "" :: out
  for k in List.range (d+1) do
    for k' in List.range (d+1) do
      if k < k' && s.stackVerts[k]! == s.stackVerts[k']! then out := bad "path" s!"{k} {k'}" :: out
  -- (ii) `top = above ++ vt :: below` once `hv`, `vt` the (unique) entry holding `V v`,
  -- `topDepth ≤ d` (after a type-2 close with the vertex entry it is merged: `vStart ≠ v`);
  -- `above` are the entries of the outs finished after the vertex push (`doneV`)
  let isV (t : TEntry) : Bool := (spanItems t).contains (vertItem v)
  let above := if hv then top.takeWhile (fun t => !isV t) else []
  let vt := if hv then (top.dropWhile (fun t => !isV t)).head? else none
  let below := if hv then (top.dropWhile (fun t => !isV t)).drop 1 else top
  if hv then
    match vt with
    | some t => if d < t.topDepth then out := bad "vert_le" (showT t) :: out
    | none => out := bad "vert_missing" "" :: out
  if hasParent s (vertItem v) then out := bad "v_root" "" :: out
  for t in s.tstack do
    if isV t && !(hv && vt.any (sameT t)) then out := bad "vert_free" (showT t) :: out
  for t in below do
    if t.vStart == v then out := bad "open_bot" (showT t) :: out
  if top ≠ [] && !endPush && !done.any (fun o => o.1.cls.lowval d < d) then out := bad "top_ret" "" :: out
  for t in base do
    if t.vStart == v then out := bad "base_bot" (showT t) :: out
  -- (ii') the entry holding `V v`: if still started at `v` it is the untouched vertex entry
  match vt with
  | some t =>
    if t.vStart == v && t.topDepth ≠ d then out := bad "vt_vstart" (showT t) :: out
    if t.vStart == v && spanItems t ≠ [vertItem v] then out := bad "vt_single" (showT t) :: out
  | none => pure ()
  -- candidate clauses for `tree_comp_shape`
  let retDone := done.filter (fun o => o.1.cls.lowval d < d)
  let allT1 := retDone.all (fun o => o.1.cls.isType1)
  let oneLow := match retDone with
    | [] => true
    | o :: _ => retDone.all (fun o' => o'.1.cls.lowval d == o.1.cls.lowval d)
  if retDone ≠ [] && !hv then out := bad "ret_hv" "" :: out
  for o in retDone do
    if o.1.cls.isType1 && !o.2 then out := bad "t1_flag" "" :: out
  for i in List.range above.length do
    for j in List.range above.length do
      if allT1 && i < j && above[i]!.topDepth ≤ above[j]!.topDepth then
        out := bad "t1_above_distinct" s!"{showT above[i]!} {showT above[j]!}" :: out
  if hv && allT1 && below ≠ [] then out := bad "t1_below" (toString (below.map showT)) :: out
  if hv && retDone ≠ [] && allT1 && oneLow && !(above.length == 1 && below.isEmpty && vt.any (·.vStart == v)) then
    out := bad "t1_top" (toString (top.map showT)) :: out
  match vt with
  | some t =>
    if hv && allT1 && t.vStart ≠ v then out := bad "t1_vt" (showT t) :: out
    if t.vStart == v && t.spans ≠ ([], [vertItem v]) &&
        !retDone.any (fun o => t.spans == setSides (!s.stackDir[o.1.cls.lowval d]!) [vertItem v] []) then
      out := bad "vt_side" (showT t) :: out
    if !endPush && t.vStart == v &&
        !retDone.any (fun o => t.spans == setSides (!s.stackDir[o.1.cls.lowval d]!) [vertItem v] []) then
      out := bad "vt_side_strict" (showT t) :: out
  | none => pure ()
  -- every open entry touches its bottom; every span item is a root; `hasVert` only after a return
  for t in s.tstack do
    if (E t) ≠ [] && !touches s (E t) t.vStart then out := bad "touch_bot" (showT t) :: out
    for i in spanItems t do
      if i ≥ s.items.size then out := bad "span_lt" (showT t) :: out
    for i in spanItems t do
      if hasParent s i then out := bad "span_root" s!"{showT t} i={i}" :: out
  for j in List.range s.items.size do
    for c in s.items[j]!.ch do
      if c ≥ s.items.size then out := bad "ch_lt" s!"{j}" :: out
  if hv && !endPush && !done.any (fun o => o.1.cls.lowval d < d) then out := bad "hv_ret" "" :: out
  for t in above do
    if t.firstIdx ≥ s.nxtEdgeIdx then out := bad "above_firstIdx" (showT t) :: out
  for t in s.tstack do
    if t.vStart == v && !isV t && t.firstIdx ≥ s.nxtEdgeIdx then out := bad "vfirst" (showT t) :: out
  if !hv && done.any (·.2) then out := bad "noVert_after" "" :: out
  if !hv && !(connB s (edgesBelow s (vertItem v)) && twoAttB s (edgesBelow s (vertItem v)) v v) then
    out := bad "vert_book" "" :: out
  if !hv then
    for t in s.tstack do
      if (E t).any (edgesBelow s (vertItem v)).contains then out := bad "vert_disj" (showT t) :: out
  -- `top` was made by this vertex's outs: bottoms off the path above, `V`/`Q` items of the finished
  -- subtrees only, edges exactly the finished outs' sub-ear edges (returning outs covered by `top`)
  let doneVerts := done.flatMap fun o => match o.1 with | .tree _ _ c => c.verts | .back .. => []
  let doneE := done.flatMap (subEdgesL ·.1)
  for t in top do
    if (List.range d).any (fun k => s.stackVerts[k]! == t.vStart) then out := bad "top_bot" (showT t) :: out
    for i in spanItems t do
      if 1 ≤ i && i < 1 + s.g.nv && i - 1 ≠ v && !doneVerts.contains (i - 1) then out := bad "top_vitems" s!"{showT t} i={i}" :: out
      if 1 + s.g.nv ≤ i && i < 1 + s.g.nv + s.g.ne && !doneE.contains (i - 1 - s.g.nv) then out := bad "top_qitems" s!"{showT t} i={i}" :: out
    if (E t).any (fun e => !doneE.contains e) then out := bad "top_edges" (showT t) :: out
  for o in done do
    if o.1.cls.lowval d < d then
      for e in subEdgesL o.1 do
        if !top.any (fun t => (E t).contains e) then out := bad "done_cover" s!"e={e}" :: out
  let bdDone := (done.map (·.1)).filter fun o => d ≤ o.cls.lowval d
  let afterV := (done.filter (·.2)).map (·.1)
  let lvs := afterV.map (·.cls.lowval d)
  let allT1 (l : Nat) : Bool := afterV.all fun o => o.cls.lowval d ≠ l || o.cls.isType1
  for o in afterV do
    if !(o.cls.lowval d < d) then out := bad "afterV_ret" s!"lv={o.cls.lowval d}" :: out
  let tops := above
  for t in tops do
    if t.vStart ≠ v then out := bad "top_vstart" (showT t) :: out
    if d ≤ t.topDepth then out := bad "top_depth" (showT t) :: out
    if !lvs.contains t.topDepth then out := bad "top_lowval" (showT t) :: out
    if (E t) == [] then out := bad "top_nonempty" (showT t) :: out
    if !touches s (E t) v then out := bad "top_touch_bot" (showT t) :: out
    if !touches s (E t) s.stackVerts[t.topDepth]! then out := bad "top_touch_top" (showT t) :: out
    if !onSide t s.stackDir[t.topDepth]! then out := bad "top_side" (showT t) :: out
    if allT1 t.topDepth then
      match spanItems t with
      | [i] => if hasParent s i then out := bad "top_root" (showT t) :: out
      | _ => out := bad "top_single" (showT t) :: out
    for x in List.range s.g.nv do
      if touches s (E t) x && x ≠ v && !interiorB s (E t) x then
        if !(List.range (d+1)).any (fun k => t.topDepth ≤ k && s.stackVerts[k]! == x) then
          out := bad "top_att" s!"{showT t} x={x}" :: out
        if allT1 t.topDepth && x ≠ s.stackVerts[t.topDepth]! then out := bad "top_att1" s!"{showT t} x={x}" :: out
  for l in lvs do
    if !tops.any (fun t => t.topDepth == l) && vt.all (fun t => t.topDepth ≠ l) then out := bad "top_present" s!"lv={l}" :: out
  let rec noninc : List TEntry → Bool
    | a :: b :: r => a.topDepth ≥ b.topDepth && noninc (b :: r)
    | _ => true
  let rec fdec : List TEntry → Bool
    | a :: b :: r => a.firstIdx > b.firstIdx && fdec (b :: r)
    | _ => true
  if !noninc tops then out := bad "top_noninc" "" :: out
  if !fdec tops then out := bad "top_first_dec" "" :: out
  -- ownership: `tops` hold exactly the returning outs' sub-ear edges, `V v` the boundary blocks
  let retE := afterV.flatMap subEdgesL
  let bdE := bdDone.flatMap subEdgesL
  let EVt := vt.map (fun t => E t) |>.getD []
  let EV := edgesBelow s (vertItem v)
  for t in tops do
    if (E t).any (fun e => !retE.contains e) then out := bad "top_sub_edges" (showT t) :: out
  for e in retE do
    if !tops.any (fun t => (E t).contains e) && !EVt.contains e then out := bad "top_cover" s!"e={e}" :: out
  for e in bdE do
    if !EV.contains e then out := bad "bd_under_vert" s!"e={e}" :: out
  if EV.any (fun e => !bdE.contains e) then out := bad "vert_edges" "" :: out
  if !pairwiseDisjB s s.tstack then out := bad "disj" "" :: out
  if !pairwiseSpanDisjB s.tstack then out := bad "span_disj" "" :: out
  -- freshness of the remaining outs' items and vertices
  for o in rest do
    for e in subEdgesL o do
      if hasParent s (edgeItem s.g e) then out := bad "q_root" s!"e={e}" :: out
      if (Items.ch s.items (edgeItem s.g e)) ≠ [] then out := bad "q_ch" s!"e={e}" :: out
      if s.tstack.any (fun t => (spanItems t).contains (edgeItem s.g e)) then out := bad "q_free" s!"e={e}" :: out
    match o with
    | .tree _ _ child =>
      for y in child.verts do
        if hasParent s (vertItem y) then out := bad "vy_root" s!"y={y}" :: out
        if (Items.ch s.items (vertItem y)) ≠ [] then out := bad "vy_ch" s!"y={y}" :: out
        if s.tstack.any (fun t => (spanItems t).contains (vertItem y)) then out := bad "vy_free" s!"y={y}" :: out
        if (List.range (d+1)).any (fun k => s.stackVerts[k]! == y) then out := bad "path_rest" s!"y={y}" :: out
        if s.tstack.any (fun t => t.vStart == y) then out := bad "vy_bot" s!"y={y}" :: out
        if s.tstack.any (fun t => touches s (E t) y) then out := bad "vy_touch" s!"y={y}" :: out
    | .back .. => pure ()
  return out

/-- Two-state facts of a finished subtree walk (`s₀` at entry to `walkTree y D`, `s` at its end):
a fixed item that was a parentless non-span item and is not a `V`/`Q` item of the subtree is still
parentless. -/
def leftCheck (seed y D : Nat) (t : DfsTree) (s₀ s : WalkState) : List V := Id.run do
  let mut out : List V := []
  let bad (k : String) (info : String) : V :=
    ⟨seed, y, D, "left", false, s!"left_{k}", s!"{info} | stack={s.tstack.map showT}"⟩
  for i in List.range (1 + s₀.g.nv + s₀.g.ne) do
    if 0 < i && !hasParent s₀ i &&
        !t.verts.any (fun w => vertItem w == i) && !t.edges.any (fun e => edgeItem s₀.g e == i) &&
        hasParent s i then
      out := bad "fresh_root" s!"i={i}" :: out
  return out

/-- Frame facts of the child's walk used by the assembly: `edgeItem e` and `vertItem v` stay roots,
and the edges below `vertItem v` are unchanged. -/
def keptCheck (seed v d : Nat) (o : DfsOut) (hv : Bool) (s₂ s : WalkState) : List V :=
  match o with
  | .back .. => []
  | .tree e cls child => Id.run do
    let mut out : List V := []
    let bad (k : String) (info : String) : V := ⟨seed, v, d, s!"tree e={e}", hv, s!"kept_{k}", info⟩
    if d ≤ cls.lowval d then
      for e' in e :: child.edges do
        let (a, b) := s.g.edges[e']!
        if !(v :: child.verts).contains a || !(v :: child.verts).contains b then
          out := bad "ends" s!"e'={e'}" :: out
    if hasParent s (edgeItem s.g e) then out := bad "q_root" "" :: out
    if hasParent s (vertItem v) then out := bad "v_root" "" :: out
    if edgesBelow s (vertItem v) ≠ edgesBelow s₂ (vertItem v) then out := bad "v_below" "" :: out
    let childItem (j : Nat) : Bool :=
      child.verts.any (fun w => vertItem w == j) || child.edges.any (fun e' => edgeItem s.g e' == j)
    let parents (s : WalkState) (j : Nat) : List Nat := (List.range s.items.size).filter (isParent s · j)
    for j in List.range s₂.items.size do
      if !childItem j then
        if s.items[j]!.ch ≠ s₂.items[j]!.ch then out := bad "items_ch" s!"j={j}" :: out
        if parents s j ≠ parents s₂ j then out := bad "items_par" s!"j={j}" :: out
    return out

/-- Freshness of the top entries (`EarCtx` over `base₀`, at every `walkOuts` step): every item in a
top entry's spans is a `V`/`Q` item of the enclosing tree (`v`, its outs' vertices and edges) or was
allocated after the enclosing `walkTree` entry (`sz₀` = items then). -/
def freshCheck (seed v d : Nat) (outs₀ : List DfsOut) (sz₀ : Nat) (base₀ : List TEntry) (s : WalkState) :
    List V := Id.run do
  let mut out : List V := []
  let top := s.tstack.take (s.tstack.length - base₀.length)
  let treeItem (j : Nat) : Bool :=
    (v :: DfsOut.vertsList outs₀).any (fun w => vertItem w == j) ||
    (DfsOut.edgesList outs₀).any (fun e => edgeItem s.g e == j)
  for t in top do
    for i in spanItems t do
      if !treeItem i && i < sz₀ then
        out := ⟨seed, v, d, "", false, "fresh_span", s!"i={i} t={showT t}"⟩ :: out
  return out

/-- Per-site frame of `walkOut v d o` (from its entry state `sIn` to the state after `finishEdge`):
every item allocated at entry other than `V v`, `Q o.e`, the out's subtree items and the items in
the spans of the top entries (above `base₀`) keeps its children and its parents. -/
def siteKeptCheck (seed v d : Nat) (o : DfsOut) (hv : Bool) (base₀ : List TEntry) (sIn s : WalkState) :
    List V := Id.run do
  let mut out : List V := []
  let bad (k : String) (info : String) : V := ⟨seed, v, d, s!"e={o.e}", hv, s!"site_{k}", info⟩
  let top := sIn.tstack.take (sIn.tstack.length - base₀.length)
  let overts := match o with | .tree _ _ child => child.verts | .back .. => []
  let excl (j : Nat) : Bool :=
    j == vertItem v || j == edgeItem sIn.g o.e || overts.any (fun w => vertItem w == j) ||
    (subEdgesL o).any (fun e' => edgeItem sIn.g e' == j) || top.any (fun t => (spanItems t).contains j)
  let parents (s : WalkState) (j : Nat) : List Nat := (List.range s.items.size).filter (isParent s · j)
  for j in List.range sIn.items.size do
    if !excl j then
      if s.items[j]!.ch ≠ sIn.items[j]!.ch then out := bad "ch" s!"j={j}" :: out
      if parents s j ≠ parents sIn j then out := bad "par" s!"j={j}" :: out
  if (s.tstack.drop (s.tstack.length - base₀.length)).map showT ≠ base₀.map showT then
    out := bad "base" "" :: out
  let top' := s.tstack.take (s.tstack.length - base₀.length)
  for t in top' do
    for i in spanItems t do
      if !excl i && i < sIn.items.size then out := bad "span" s!"i={i} t={showT t}" :: out
  return out

/-- The returning tree-edge admissions (`TreeSite.tree_ret_vert`/`tree_ret_noVert`): relative to
the pre-child state `s₂` (after the first-edge push), the state `sX` after loops 1–2 (and, with a
vertex entry, after `closeVert'`) is `R ++ s₂.tstack` with `s₂.tstack` untouched; the items of `R`
are the child's, `Q e` or fresh, roots, and own exactly the out's edges; the out's edges touch only
`v`, the child's vertices and the path at depths `[lowval, d]`. With a vertex entry `R` is the single
`(v, lowval)` entry (one item if type 1); without (type 2) its entries start at child vertices. -/
def retCheck (seed v d : Nat) (o : DfsOut) (hv : Bool) (s₂ s : WalkState) : List V :=
  match o with
  | .back .. => []
  | .tree e cls child => Id.run do
    let lowval := cls.lowval d
    if d ≤ lowval then return []
    let mut out : List V := []
    let bad (k : String) (info : String) : V :=
      ⟨seed, v, d, s!"tree e={e} lv={lowval} t1={cls.isType1}", hv, s!"ret_{k}", info⟩
    let orig := s₂.tstack.length
    let sX := if hv then
        WalkState.after (WalkState.closeVert' v s.stackDir[d]! cls.isType1 orig (WalkState.feSingle d o s))
          (WalkState.feS₂ d o s)
      else WalkState.feS₂ d o s
    let n := sX.tstack.length
    if n < orig then return [bad "short" ""]
    let R := sX.tstack.take (n - orig)
    let baseX := sX.tstack.drop (n - orig)
    let info := s!"R={R.map showT} base={s₂.tstack.map showT}"
    if !(baseX.length == orig && (List.zip baseX s₂.tstack).all fun (a, b) => sameT a b) then
      out := bad "base" info :: out
    if sX.g.edges ≠ s₂.g.edges || sX.g.nv ≠ s₂.g.nv then out := bad "g" "" :: out
    if (List.range (d+1)).any (fun k => sX.stackVerts[k]! ≠ s₂.stackVerts[k]!) then out := bad "sv" "" :: out
    if (List.range d).any (fun k => sX.stackDir[k]! ≠ s₂.stackDir[k]!) then out := bad "sd" "" :: out
    if sX.stackDir[d]! ≠ !s₂.stackDir[lowval]! then out := bad "sdd" "" :: out
    if sX.items.size < s₂.items.size then out := bad "size" "" :: out
    if sX.nxtEdgeIdx < s₂.nxtEdgeIdx then out := bad "nxt" "" :: out
    let childItem (j : Nat) : Bool :=
      child.verts.any (fun w => vertItem w == j) || child.edges.any (fun e' => edgeItem s₂.g e' == j) ||
        edgeItem s₂.g e == j
    let parents (s : WalkState) (j : Nat) : List Nat := (List.range s.items.size).filter (isParent s · j)
    for j in List.range s₂.items.size do
      if !childItem j then
        if sX.items[j]!.ch ≠ s₂.items[j]!.ch then out := bad "keep_ch" s!"j={j}" :: out
        if parents sX j ≠ parents s₂ j then out := bad "keep_par" s!"j={j}" :: out
    for j in List.range sX.items.size do
      for c in sX.items[j]!.ch do
        if c ≥ sX.items.size then out := bad "ch_lt" s!"{j}" :: out
    let oldSpans := s₂.tstack.flatMap spanItems
    let subE := subEdgesL o
    for t in R do
      for i in spanItems t do
        if sX.items.size ≤ i then out := bad "span_lt" s!"i={i}" :: out
        if hasParent sX i then out := bad "span_root" s!"i={i}" :: out
        if oldSpans.contains i then out := bad "span_old" s!"i={i}" :: out
        if i < s₂.items.size && !childItem i then out := bad "span_owned" s!"i={i}" :: out
      let E := entryEdges sX t
      if E ≠ [] && !touches sX E t.vStart then out := bad "touch_bot" (showT t) :: out
    let ER := R.flatMap (entryEdges sX)
    if !sameEdges ER subE then out := bad "edges" s!"ER={ER} sub={subE} {info}" :: out
    if !pairwiseDisjB sX R then out := bad "disj" info :: out
    if !pairwiseSpanDisjB R then out := bad "span_disj" info :: out
    for t in R do
      if t.vStart ≠ v && !child.verts.contains t.vStart then out := bad "vstart" (showT t) :: out
    for x in List.range s₂.g.nv do
      if touches sX ER x && x ≠ v && !child.verts.contains x &&
          !(List.range (d+1)).any (fun k => lowval ≤ k && s₂.stackVerts[k]! == x) then
        out := bad "touch" s!"x={x}" :: out
    if hv then
      match R with
      | [x] =>
        if x.vStart ≠ v || x.topDepth ≠ lowval then out := bad "x" (showT x) :: out
        if !onSide x s₂.stackDir[lowval]! then out := bad "x_side" (showT x) :: out
        if x.firstIdx < s₂.nxtEdgeIdx || sX.nxtEdgeIdx ≤ x.firstIdx then out := bad "x_first" (showT x) :: out
        if cls.isType1 && (spanItems x).length ≠ 1 then out := bad "x_single" (showT x) :: out
        if !touches sX ER sX.stackVerts[lowval]! then out := bad "x_touch_top" (showT x) :: out
        let S₂ := WalkState.feS₂ d o s
        let R₂ := S₂.tstack.take (S₂.tstack.length - orig)
        match R₂.getLast? with
        | some vy =>
          if vy.firstIdx < s₂.nxtEdgeIdx || sX.nxtEdgeIdx ≤ vy.firstIdx then out := bad "x_first_pre" (showT vy) :: out
        | none => out := bad "x_pre_empty" info :: out
        if !touches S₂ (R₂.flatMap (entryEdges S₂)) s₂.stackVerts[lowval]! then out := bad "x_top_pre" info :: out
        if R₂.any (fun t => t.topDepth < lowval) then out := bad "x_mid_top" info :: out
        for w in List.range s₂.g.nv do
          if touches sX ER w && w ≠ v && !interiorB sX ER w &&
              !(List.range (d+1)).any (fun k => lowval ≤ k && sX.stackVerts[k]! == w) then
            out := bad "x_att" s!"w={w} {showT x}" :: out
          if cls.isType1 && touches sX ER w && w ≠ v && w ≠ sX.stackVerts[lowval]! && !interiorB sX ER w then
            out := bad "x_single_att" s!"w={w} {showT x}" :: out
      | _ => out := bad "vert_shape" info :: out
    else
      if cls.isType1 then out := bad "noVert_t1" "" :: out
      if R.isEmpty then out := bad "noVert_empty" info :: out
      match R with
      | c :: _ => if !(entryEdges sX c).contains e then out := bad "noVert_head_e" info :: out
      | [] => pure ()
      for t in R do
        if t.vStart == v then out := bad "noVert_v" (showT t) :: out
        if (List.range (d+1)).any (fun k => s₂.stackVerts[k]! == t.vStart) then
          out := bad "noVert_path" (showT t) :: out
    return out

/-- The child's end-of-outs stack at a component edge (`TreeSite.comp_shape`): `hasVert`, and exactly
`[(y, d-1), V y]` above the parent's stack, the `(y, d-1)` entry on side 1, `V y` on side 2. -/
def compEndCheck (seed v d : Nat) (hv : Bool) (top : List TEntry) : List V := Id.run do
  let mut out : List V := []
  let bad (k : String) (info : String) : V := ⟨seed, v, d, "compEnd", hv, s!"compEnd_{k}", info⟩
  if !hv then out := bad "hv" "" :: out
  match top with
  | [t₁, t₂] =>
    if t₁.vStart ≠ v || t₁.topDepth ≠ d - 1 || t₁.spans.2 ≠ [] then out := bad "t1" (showT t₁) :: out
    if t₂.vStart ≠ v || t₂.topDepth ≠ d || t₂.spans ≠ ([], [vertItem v]) then out := bad "t2" (showT t₂) :: out
  | _ => out := bad "shape" (toString (top.map showT)) :: out
  return out

mutual
partial def iTree (seed : Nat) (t : DfsTree) (d : Nat) (eb : EB) (s : WalkState) (pcls : Option OutClass) : WalkState × EB × List V :=
  match t with
  | .node v outs =>
    let orig := s.tstack.length
    let s₀ := s
    let s := { s with stackVerts := s.stackVerts.set! d v }
    let base₀ := s.tstack
    let bE₀ := s.tstack.map (entryEdges s)
    let sv₀ := s.stackVerts.toList.take (d+1)
    let sd₀ := s.stackDir.toList.take d
    let (hv, s, eb, vs, done) := iOuts seed v d outs false eb s [] base₀ bE₀ sv₀ sd₀ s₀.items.size
    let vs := vs ++ (match pcls with
      | some .component => compEndCheck seed v d hv (s.tstack.take (s.tstack.length - orig))
      | _ => [])
    let s := if hv then s else ((setStackDir d true *> pushVertTstack v d).run s).2
    let vs := vs ++ ctxCheck seed v d done [] true base₀ bE₀ sv₀ sd₀ s true ++ leftCheck seed v d t s₀ s
    let subv := s.tstack.take (s.tstack.length - orig)
    let eb := match subv.reverse with
      | vy :: py :: _ => if spanItems vy == [vertItem v] then (v, py, vy) :: eb else eb
      | _ => eb
    (s, eb, vs)
partial def iOuts (seed v d : Nat) (outs : List DfsOut) (hv : Bool) (eb : EB) (s : WalkState)
    (done : List (DfsOut × Bool)) (base₀ : List TEntry) (bE₀ : List (List Nat)) (sv₀ : List Nat) (sd₀ : List Bool)
    (sz₀ : Nat) :
    Bool × WalkState × EB × List V × List (DfsOut × Bool) :=
  let cvs := ctxCheck seed v d done outs hv base₀ bE₀ sv₀ sd₀ s ++
    freshCheck seed v d (done.map (·.1) ++ outs) sz₀ base₀ s
  match outs with
  | [] => (hv, s, eb, cvs, done)
  | o :: rest =>
    let hvF := hv || (o.cls.lowval d < d && o.cls.isType1)
    let (hv, s, eb, vs) := iOut seed v d o hv eb s base₀
    let (hv', s', eb, vs', done') := iOuts seed v d rest hv eb s (done ++ [(o, hvF)]) base₀ bE₀ sv₀ sd₀ sz₀
    (hv', s', eb, cvs ++ vs ++ vs', done')
partial def iOut (seed v d : Nat) (o : DfsOut) (hv : Bool) (eb : EB) (s : WalkState) (base₀ : List TEntry) :
    Bool × WalkState × EB × List V :=
  let sIn := s
  let hvIn := hv
  let lowval := o.cls.lowval d
  let s := ((do let lowDir ← stackDir lowval; setStackDir d (if lowval ≥ d then false else !lowDir) : WalkM Unit).run s).2
  let (hv, s) := if !hv && lowval < d && o.cls.isType1 then (true, ((pushVertTstack v d).run s).2) else (hv, s)
  let orig := s.tstack.length
  let s₂ := s
  let (s, eb, vs) := match o with
    | .tree _ cls child => iTree seed child (d+1) eb { s with firstOccurrence := s.firstOccurrence.set! d s.g.ne } (some cls)
    | .back .. => (s, eb, [])
  let vs := vs ++ keptCheck seed v d o hv s₂ s
  let vs := vs ++ retCheck seed v d o hv s₂ s
  let vs := vs ++ earCheck seed v d o orig hv s ++ ebCheck seed v d o orig hv eb s ++ specCheck seed v d o orig hv s ++ closeCheck seed v d o orig hv s
  let (hv', s) := (finishEdge v d o orig hv).run s
  let vs := vs ++ siteKeptCheck seed v d o hvIn base₀ sIn s
  (hv', s, eb, vs)
end

def iForest (seed : Nat) (forest : List DfsTree) (s : WalkState) : WalkState × List V :=
  forest.foldl (fun (s, vs) t =>
    let (s, _, vs') := iTree seed t 0 [] s none
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

def runSeed (seed : Nat) (tern : Bool := false) : Bool × List V :=
  let g := randGraph seed
  let f := g.dfsForest [] []
  let (s, vs) := iForest seed f (WalkState.init g tern)
  let ref := g.walk tern f
  (s.items.toList.map (fun it => (it.vs, it.ch)) == ref.items.toList.map (fun it => (it.vs, it.ch)) && s.tstack.length == ref.tstack.length, vs)

def summarize (lo hi : Nat) (tern : Bool := false) : IO Unit := do
  let mut counts : List (String × Nat) := []
  let mut okAll := true
  let mut shown : List (String × String) := []
  let mut checks := 0
  for seed in List.range' lo (hi - lo) do
    let (ok, vs) := runSeed seed tern
    if !ok then okAll := false
    checks := checks + 1
    for v in vs do
      counts := match counts.find? (·.1 == v.kind) with
        | some _ => counts.map fun (k, n) => if k == v.kind then (k, n+1) else (k, n)
        | none => counts ++ [(v.kind, 1)]
      if (shown.filter (·.1 == v.kind)).length < 1 then
        shown := shown ++ [(v.kind, s!"{repr v} g={repr (randGraph seed).edges}")]
  IO.println s!"tern={tern} instrumentation matches library walk: {okAll}"
  IO.println s!"violations by kind: {counts}"
  for (_, l) in shown do IO.println l


/-! ### Trace of one seed -/
mutual
partial def showTree (t : DfsTree) (ind : String) : String :=
  match t with
  | .node v outs => s!"{ind}{v}\n" ++ String.join (outs.map (showOut · ind))
partial def showOut (o : DfsOut) (ind : String) : String :=
  match o with
  | .back e dest cls => s!"{ind}  back e{e}->{dest} {repr cls}\n"
  | .tree e cls child => s!"{ind}  tree e{e} {repr cls}\n" ++ showTree child (ind ++ "    ")
end
def traceSeed (seed : Nat) : IO Unit := do
  let g := randGraph seed
  IO.println s!"g = {repr g.edges} nv={g.nv}"
  for t in g.dfsForest [] [] do IO.println (showTree t "")
  let (_, vs) := runSeed seed
  for v in vs do
    if v.kind == "TRACE" then IO.println s!"finishEdge curV={v.curV} d={v.d} {v.o} hv={v.hasVert} {v.info}"
#eval summarize 0 3000
#eval summarize 0 3000 true

