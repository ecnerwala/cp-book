import WalkInvCheck.Ear
import WalkInvCheck.Ranges
import WalkInvCheck.R
import WalkInvCheck.RContent
import WalkInvCheck.St
import WalkInvCheck.Extra
import WalkInvCheck.E3False
import WalkInvCheck.Planar
/-!
# `WalkInvCheck`: the executable conjunction of the backbone invariant (PROOF.md §4.7)

Build with `lake build check_walkinv`; run `./.lake/build/bin/check_walkinv [lo hi]` (default seeds
`0..3000`, i.e. `lo = 0`, `hi = 3001`), BOTH `ternarize` modes, graphs `randGraph seed` of
`checks/EarCheck.lean`. The walk is re-run by ONE instrumented traversal (`iTree`/`iOuts`/`iOut`
mirror `walkTree`/`walkOuts`/`walkOut` around the library's `finishEdge`) and every field of the
proposed `WalkInv` is evaluated as a `Bool` mirror at the sites its layer states it at:

* between-edge sites (`walkOuts v d` with `done` finished, before every `walkOut` and at the end of
  the outs): `inv.*` (`Inv' d`), `ear.ctx_*` (`EarCtx`), `ranges.own_*` (`OwnedD`, at the tree
  entry / end), `st.*_cur` (`StRead`/`StItems`/m8–m11/m16–m17 and the new per-segment
  `StLive ↔ openBlock` pairing `st.live_cur`, `st.live_lower`), `rskel.*`, `e4.*`;
* the child-entry site (`walkTree child (d+1)`): `e1.entry`, `e3.stab_single` (the original E3
  `buried_vacuous` is refuted in `WalkInvCheck/E3False.lean`), `r.stab`, `r.anc`;
* before every `finishEdge` (and on its `feS₀`/loop-1 iterates/`feS₁`/`feS₂`/P/vertex/tail sub-sites
  the Ranges/R/ear contracts are stated at): `ear.*` (`EarFinish`, `Loop1Spec`, `EarLate`/`EarClose`,
  ear bottoms, `Frontier`), `ranges.*` (`RangesInv`, `FinishAdj`, `CloseInv`/`CloseCtx`, the P /
  vertex / loop-1 sites), `r.*` (contract B `EntryR`, disjointness, loop-1 emulation `RTop`/`RBranch`,
  `FinishRShape`, `RInvG` base frame — on block graphs, as `checks/RFinishEdgeCheck.lean`),
  `st.*_pre`, `st.site_*` (m13/m15), `st.live_pre`, `e2.bd_free`, `e5.*`;
* after every `finishEdge`: `ranges.own_*` (post), `r.postB`/`r.base_post`, `st.boundary`;
* the end of every `walkTree`: `ear.compEnd_*`, `ear.left_*`, `st.*_end`, `st.live_end`,
  `ranges.own_vcover`/`vert_past`, `e2.bd_free` (tail).

Violations are counted PER FIELD with the smallest violating seed; the instrumented walk is compared
with `Graph.walk` (items and tstack length) on every run. `check_walkinv minimize <field> <seed>`
edge-minimises a violating graph (and reports every violation on it); `check_walkinv small <field>
<n>` searches `n` small random graphs for the smallest violation.
-/
namespace WalkInvCheck
open Spqr WalkM

structure Ctx where
  seed : Nat
  g : Graph
  σ : List Nat
  D : R.Dfs
  block : Bool

/-- One frame of the open path: the vertex at depth `k`, its finished outs and the open one, the
bookkeeping of the layers (`nEntry`/`baseLen`: Ranges `OwnedD`; `origLen`: R base frames and the
segment split; `qs`: the pieces of the lower segment `new_k`). -/
structure Frame where
  v : Nat
  done : List DfsOut
  o : DfsOut
  nEntry : Nat
  baseLen : Nat
  origLen : Nat
  qs : List StPiece
deriving Inhabited

def pathFrames (fs : List Frame) : List PathFrame := fs.map fun f => ⟨f.v, f.done, f.o⟩

def ofEar (tern : Bool) (v : Ear.V) : Viol :=
  ⟨v.seed, tern, s!"ear.{v.kind}", s!"v={v.curV} d={v.d} {v.o} hv={v.hasVert} {v.info}"⟩
/-- The `cand_*` kinds of `checks/RangesInvCheck.lean` are candidate observations, not fields of
`RangesInv` (PROOF.md §4.6); they are dropped here. -/
def ofRanges (l : List Ranges.V) : List Viol :=
  l.filterMap fun v => if v.kind.startsWith "cand_" then none else
    some ⟨v.seed, v.tern, s!"ranges.{v.kind}", s!"v={v.curV} d={v.d} {v.o} {v.info}"⟩
def ofStrs (c : Ctx) (tern : Bool) (field : String) (site : String) (l : List String) : List Viol :=
  l.filterMap fun b => if b.startsWith "stat:" then none else some ⟨c.seed, tern, field, s!"{site} {b}"⟩

/-- The `(field, info)` pairs of `RC` as `rc.<field>`. -/
def ofRC (c : Ctx) (tern : Bool) (site : String) (l : List (String × String)) : List Viol :=
  l.map fun (f, b) => ⟨c.seed, tern, s!"rc.{f}", s!"{site} {b}"⟩

/-- The R strings of `loop1Emu` carry their field as a prefix. -/
def ofEmu (c : Ctx) (tern : Bool) (site : String) (l : List String) : List Viol :=
  l.filterMap fun b =>
    if b.startsWith "stat:" then none
    else if b.startsWith "rclose" then some ⟨c.seed, tern, "r.rclose", s!"{site} {b}"⟩
    else if b.startsWith "rbranch" then some ⟨c.seed, tern, "r.rbranch", s!"{site} {b}"⟩
    else some ⟨c.seed, tern, "r.emu", s!"{site} {b}"⟩

/-- The St fields at a site with the current segment `new`, the completed blocks `blocks`
(`StItems`) and the truncated reference's blocks `tblocks`. -/
def stSite (c : Ctx) (site : String) (v d : Nat) (new : List TEntry) (ps : List StPiece)
    (blocks tblocks : List StBlock) (s : WalkState) : List Viol :=
  let g := c.g
  (if St.stReadB s.items new ps then [] else
    [⟨c.seed, s.ternarize, s!"st.read_{site}", s!"v={v} d={d} read={(readStack new).flatMap (St.leavesB s.items)} ps={stNest ps} stack={St.showStack s.tstack}"⟩]) ++
  ofStrs c s.ternarize s!"st.items_{site}" s!"v={v} d={d}" (St.stItemsB g s blocks) ++
  ofStrs c s.ternarize s!"st.trunc_{site}" s!"v={v} d={d}" (St.truncItemsB g s tblocks) ++
  ofStrs c s.ternarize s!"st.ctx_{site}" s!"v={v} d={d}" (St.truncCtxB g s tblocks v d)

/-- The per-segment `StLive ↔ openBlock` pairing of the lower segments `k < fs.length`. -/
def liveLower (c : Ctx) (site : String) (fs : List Frame) (s : WalkState) : List Viol :=
  let pfs := pathFrames fs
  (List.range fs.length).flatMap fun k =>
    let f := fs[k]!
    let n := s.tstack.length
    let new := (s.tstack.drop (n - f.origLen)).take (f.origLen - f.baseLen)
    let b := openBlock c.g (pfs.take k) (DirsOf s k) f.qs
    (if f.origLen ≤ n then [] else [⟨c.seed, s.ternarize, "st.seg_lower", s!"{site} k={k} origLen={f.origLen} > len={n}"⟩]) ++
    (if St.stReadB s.items new f.qs then [] else
      [⟨c.seed, s.ternarize, "st.segread_lower", s!"{site} k={k} v={f.v} new={new.map Extra.showT} qs={stNest f.qs}"⟩]) ++
    ofStrs c s.ternarize "st.live_lower" s!"{site} k={k} v={f.v} new={new.map Extra.showT}" (Extra.stLiveB c.g s.items new b)

/-- The fields stated at every site. -/
def everySite (c : Ctx) (site : String) (D : Nat) (s : WalkState) : List Viol :=
  Extra.invCheck c.seed site D s ++ Extra.rskelCheck c.seed site s ++ Extra.e4Check c.seed site s

mutual
partial def iTree (c : Ctx) (prev : List DfsTree) (fs : List Frame) (n : Nat) (t : DfsTree) (d : Nat)
    (eb : Ear.EB) (visited : List Nat) (s : WalkState) : WalkState × Ear.EB × List Nat × List Viol :=
  match t with
  | .node v outs =>
    let tern := s.ternarize
    let g := c.g
    let s₀ := s
    let baseLen := s.tstack.length
    let vs := (Extra.e3Check c.seed d v s).map fun x =>
      { x with info := s!"block={c.block} frames(v,done,e)={fs.map fun f => (f.v, f.done.length, f.o.e)} {x.info}" }
    let s' := { s with stackVerts := s.stackVerts.set! d v }
    let vs := vs ++ (if !c.block then [] else
      s.tstack.flatMap fun t =>
        if (R.entryR c.D s t).isEmpty && !(R.entryR c.D s' t).isEmpty then
          [⟨c.seed, tern, "r.stab", s!"v={v} d={d} {R.showT t} loses EntryR under stackVerts.set! d v"⟩] else [])
    let vs := vs ++ ((List.range (d + 1)).flatMap fun k =>
      let a := s'.stackVerts[k]!
      if c.D.depth[a]! == k && R.ancP c.D.parent a v then [] else
        [⟨c.seed, tern, "r.anc", s!"v={v} d={d} k={k} sv[k]={a} depth={c.D.depth[a]!}"⟩])
    let vs := vs ++ (if !c.block || d = 0 then [] else ofRC c tern s!"entry v={v} d={d}" (RC.entryContent c.D d v s s'))
    let s := s'
    let base₀ := s.tstack
    let bE₀ := s.tstack.map (Ear.entryEdges s)
    let sv₀ := s.stackVerts.toList.take (d+1)
    let sd₀ := s.stackDir.toList.take d
    let own : Ranges.OwnCtx := ⟨fs.map (·.nEntry) ++ [n], fs.map (·.baseLen) ++ [baseLen], visited ++ [v]⟩
    let vs := vs ++ ofRanges (Ranges.checkOwned c.seed c.σ own n v d "entry" s ++ Ranges.checkPieceInv c.seed "entry" d s)
    let (hv, s, eb, visited, done, vs') := iOuts c prev fs n v d outs [] false base₀ bE₀ sv₀ sd₀ eb own s
    let vs := vs ++ vs'
    let own := { own with visited := visited }
    let vs := vs ++ (match fs.getLast? with
      | some f => match f.o.cls with
        | .component => (Ear.compEndCheck c.seed v d hv (s.tstack.take (s.tstack.length - baseLen))).map (ofEar tern)
        | _ => []
      | none => [])
    let nEnd := n + t.edgePostorder.length
    let vs := vs ++ ofRanges (Ranges.checkOwned c.seed c.σ own nEnd v d "end" s ++ Ranges.checkPieceInv c.seed "end" d s)
    let vs := vs ++ (if hv then ofRanges (Ranges.checkVCover c.seed v "end" s)
      else ofRanges (Ranges.checkVertPast c.seed c.σ v nEnd s) ++ Extra.e2Check c.seed "tail" v s ++
        (if !c.block then [] else ofRC c tern s!"tail v={v} d={d}" (RC.vertFree v s)))
    let s := if hv then s else ((setStackDir d true *> pushVertTstack v d).run s).2
    let vs := vs ++ ofRanges (Ranges.checkPieceInv c.seed "tail" d s)
    let vs := vs ++ (Ear.ctxCheck c.seed v d done [] true base₀ bE₀ sv₀ sd₀ s true ++ Ear.leftCheck c.seed v d t s₀ s).map (ofEar tern)
    let pfs := pathFrames fs
    let dirs := DirsOf s d
    let (ps, bl) := refTree g t d dirs
    let new := St.above s base₀
    let vs := vs ++ (if St.baseOk s base₀ then [] else [⟨c.seed, tern, "st.base_end", s!"v={v} d={d}"⟩])
    let vs := vs ++ stSite c "end" v d new ps (simBlocks g prev pfs dirs ++ bl) (refBlocks g (prev ++ [truncTree pfs t])) s
    let vs := vs ++ ofStrs c tern "st.live_end" s!"v={v} d={d} new={new.map Extra.showT}" (Extra.stLiveB g s.items new (openBlock g pfs dirs ps))
    let vs := vs ++ liveLower c s!"end v={v} d={d}" fs s
    let vs := vs ++ everySite c s!"end v={v}" d s
    let subv := s.tstack.take (s.tstack.length - baseLen)
    let eb := match subv.reverse with
      | vy :: py :: _ => if Ear.spanItems vy == [vertItem v] then (v, py, vy) :: eb else eb
      | _ => eb
    (s, eb, visited, vs)

partial def iOuts (c : Ctx) (prev : List DfsTree) (fs : List Frame) (n v d : Nat) (outs : List DfsOut)
    (done : List (DfsOut × Bool)) (hv : Bool) (base₀ : List TEntry) (bE₀ : List (List Nat)) (sv₀ : List Nat)
    (sd₀ : List Bool) (eb : Ear.EB) (own : Ranges.OwnCtx) (s : WalkState) :
    Bool × WalkState × Ear.EB × List Nat × List (DfsOut × Bool) × List Viol :=
  let tern := s.ternarize
  let g := c.g
  let site := s!"between v={v} d={d} done={done.length}"
  let cvs := (Ear.ctxCheck c.seed v d done outs hv base₀ bE₀ sv₀ sd₀ s).map (ofEar tern)
  let pfs := pathFrames fs
  let dirs := DirsOf s d
  let doneO := done.map (·.1)
  let (ps, bl, hvRef) := refOuts g v d dirs doneO false
  let new := St.above s base₀
  let cvs := cvs ++ (if St.baseOk s base₀ then [] else [⟨c.seed, tern, "st.base_cur", site⟩])
  let cvs := cvs ++ (if hvRef == hv then [] else [⟨c.seed, tern, "st.hv", s!"{site} ref={hvRef} walk={hv}"⟩])
  let cvs := cvs ++ stSite c "cur" v d new ps (simBlocks g prev pfs dirs ++ bl) (refBlocks g (prev ++ [truncTree pfs (.node v doneO)])) s
  let cvs := cvs ++ ofStrs c tern "st.live_cur" s!"{site} new={new.map Extra.showT}"
    ([false, true].flatMap fun b =>
      Extra.stLiveB g s.items new (openBlock g pfs dirs (ps ++ if hv then [] else [⟨b, [vertItem v]⟩])))
  let cvs := cvs ++ liveLower c site fs s
  let cvs := cvs ++ everySite c site d s
  match outs with
  | [] => (hv, s, eb, own.visited, done, cvs)
  | o :: rest =>
    let hvF := hv || (o.cls.lowval d < d && o.cls.isType1)
    let (hv, s, eb, visited, vs) := iOut c prev fs n v d o doneO hv eb own s
    let (hv', s', eb, visited, done', vs') :=
      iOuts c prev fs n v d rest (done ++ [(o, hvF)]) hv base₀ bE₀ sv₀ sd₀ eb { own with visited := visited } s
    (hv', s', eb, visited, done', cvs ++ vs ++ vs')

partial def iOut (c : Ctx) (prev : List DfsTree) (fs : List Frame) (n v d : Nat) (o : DfsOut)
    (doneO : List DfsOut) (hv : Bool) (eb : Ear.EB) (own : Ranges.OwnCtx) (s : WalkState) :
    Bool × WalkState × Ear.EB × List Nat × List Viol :=
  let tern := s.ternarize
  let g := c.g
  let σ := c.σ
  let lowval := o.cls.lowval d
  let pushNow := !hv && lowval < d && o.cls.isType1
  let vs := if hv then [] else ofRanges (Ranges.checkVertPast c.seed σ v (σ.idxOf o.block.head!) s) ++
    (if !c.block then [] else ofRC c tern s!"out-entry v={v} d={d} e={o.e}" (RC.vertFree v s))
  let vs := vs ++ (if d = 0 && !(s.tstack.isEmpty && !hv) then
    [⟨c.seed, tern, "r.root", s!"v={v} d=0 e={o.e} tstack={s.tstack.length} hv={hv}"⟩] else [])
  let vs := vs ++ (if pushNow then Extra.e2Check c.seed s!"pre-push e={o.e}" v s else [])
  let (hv₀, s) := (walkOutPre v d o hv).run s
  let hv := hv₀
  let orig := s.tstack.length
  let s₂ := s
  let pfs := pathFrames fs
  let f : Frame := ⟨v, doneO, o, n, own.origs.getLast!, orig,
    (refOuts g v d (DirsOf s d) doneO false).1 ++ (if pushNow then [⟨s.stackDir[d]!, [vertItem v]⟩] else [])⟩
  let (s, eb, visited, vs') := match o with
    | .tree _ _ child => iTree c prev (fs ++ [f]) (σ.idxOf o.e - child.edgePostorder.length) child (d+1) eb own.visited
        { s with firstOccurrence := s.firstOccurrence.set! d s.g.ne }
    | .back .. => (s, eb, own.visited, [])
  let vs := vs ++ vs'
  let own := { own with visited := visited }
  let oStr := s!"{if o.cls.isTree then "tree" else "back"} e={o.e} dest={o.dest} lv={lowval} t1={o.cls.isType1}"
  let site := s!"pre v={v} d={d} {oStr} hv={hv} orig={orig}"
  -- ear
  let vs := vs ++ (Ear.keptCheck c.seed v d o hv s₂ s ++ Ear.earCheck c.seed v d o orig hv s ++
    Ear.ebCheck c.seed v d o orig hv eb s ++ Ear.specCheck c.seed v d o orig hv s ++
    Ear.closeCheck c.seed v d o orig hv s ++ Ear.frontierCheck c.seed v d o orig hv s).map (ofEar tern)
  -- ranges
  let vs := vs ++ ofRanges (Ranges.checkOwned c.seed σ own (σ.idxOf o.e) v d "pre" s ++
    ((Ranges.adjacencySites v d o orig hv s).filter (·.1 == "P")).flatMap (fun (_, st) =>
      Ranges.checkOwned c.seed σ own (σ.idxOf o.e + 1) v d "P" st) ++
    Ranges.check c.seed σ v d o orig s ++ Ranges.checkAdj c.seed σ v d o orig hv s ++ Ranges.checkClose c.seed "pre" s ++
    Ranges.checkCtx c.seed v d o s ++ Ranges.checkP c.seed v d o orig hv s ++ Ranges.checkV c.seed v d o orig hv s ++
    Ranges.checkL1 c.seed v d o s ++ Ranges.checkRI c.seed σ v d o s ++
    Ranges.checkContent c.seed v d o orig hv s ++
    Ranges.checkCanon c.seed v d o orig hv s ++
    Ranges.checkPieceInv c.seed "pre" d s ++ Ranges.checkPiece c.seed v d o orig hv s ++
    Ranges.checkCanonInv c.seed v d s ++
    (Ranges.closeSites v d o orig hv s).flatMap (fun (st, r) => Ranges.checkClose c.seed st r) ++
    (if hv then [] else Ranges.checkVertPast c.seed σ v (σ.idxOf o.e) s))
  -- R (block graphs), before `finishEdge`
  let rsite := lowval < d
  let base := s.tstack.drop (s.tstack.length - orig)
  let rfr := (List.range fs.length).map fun k => (fs[k]!.v, k, fs[k]!.origLen)
  let vs := vs ++ (if !c.block || !rsite then [] else
    ofStrs c tern "r.preB" site (R.checkEntries c.D d base s (fun t => t.vStart != v && t.topDepth ≥ d)) ++
    ofStrs c tern "r.disj" site (R.disjoint s) ++
    (if !o.isTree then [] else
      let sp := WalkState.after (pushEdgeTstack o.dest d o.e) (WalkState.feS₀ d o s)
      let (sEnd, rbad) := R.loop1Emu c.D d s.stackDir[d]! sp sp.tstack.length
      (if sEnd.tstack.map R.showT != (WalkState.feS₁ d o s).tstack.map R.showT then
        [⟨c.seed, tern, "r.emu", s!"{site} loop-1 emulation mismatch"⟩] else []) ++
      ofEmu c tern site rbad ++
      ofStrs c tern "r.shape" site (R.shapeCheck c.D v d o orig hv s)) ++
    rfr.flatMap (fun (p, dp, n0) =>
      let bot := s.tstack.drop (s.tstack.length - n0)
      (if n0 > orig then [⟨c.seed, tern, "r.base_pre", s!"{site} p={p} dp={dp} n0={n0} > orig"⟩] else []) ++
      (if hv && n0 + 1 > orig then [⟨c.seed, tern, "r.base_pre", s!"{site} p={p} dp={dp} n0={n0} hasVert but n0+1 > orig"⟩] else []) ++
      (if o.isTree && orig + 3 > (WalkState.feS₂ d o s).tstack.length then
        [⟨c.seed, tern, "r.base_pre", s!"{site} p={p} dp={dp} close: orig+3 > len(feS₂)"⟩] else []) ++
      ofStrs c tern "r.base_pre" s!"{site} p={p} dp={dp} n0={n0}" (bot.flatMap fun t => R.settledEntry c.D p dp s t "pre")))
  let vs := vs ++ (if !c.block || !rsite || !o.isTree then [] else ofRC c tern site (RC.closeContent c.D v d o orig hv s))
  -- St, before `finishEdge`
  let dirs := DirsOf s d
  let (psChild, blChild) := match o with
    | .tree _ _ child => refTree g child (d + 1) (DirsOf s (d + 1))
    | .back .. => ([], [])
  let sub := s.tstack.take (s.tstack.length - orig)
  let tblocks := refBlocks g (prev ++ [truncTree pfs (.node v (doneO ++ [o]))])
  let vs := vs ++ (if orig ≤ s.tstack.length then [] else [⟨c.seed, tern, "st.base_pre", site⟩])
  let vs := vs ++ stSite c "pre" v d sub psChild (simBlocks g prev pfs dirs ++ (refOuts g v d dirs doneO false).2.1 ++ blChild) tblocks s
  let vs := vs ++ (St.truncSiteB g s tblocks v d o sub).filterMap fun b =>
    some ⟨c.seed, tern, if b.startsWith "ret" then "st.site_m13" else "st.site_m15", s!"{site} {b}"⟩
  let vs := vs ++ ofStrs c tern "st.live_pre" s!"{site} new={sub.map Extra.showT}"
    (Extra.stLiveB g s.items sub (openBlock g (pfs ++ [⟨v, doneO, o⟩]) (DirsOf s (d + 1)) psChild))
  let vs := vs ++ liveLower c site (fs ++ [f]) s
  -- the extra fields
  let vs := vs ++ everySite c site (if o.cls.isTree then d + 1 else d) s ++ Extra.e5Check c.seed v d o orig s
  let vs := vs ++ (if hv then [] else Extra.e2Check c.seed site v s)
  let sPre := s
  let (hv', s) := (finishEdge v d o orig hv).run s
  -- after `finishEdge`
  let post := s!"post v={v} d={d} {oStr} hv={hv'}"
  let vs := vs ++ ofRanges (Ranges.checkOwned c.seed σ own (σ.idxOf o.e + 1) v d "post" s ++ Ranges.checkPieceInv c.seed "post" d s)
  let vs := vs ++ (if hv' then ofRanges (Ranges.checkVCover c.seed v "post" s) else [])
  let vs := vs ++ (if !c.block || !rsite then [] else
    ofStrs c tern "r.postB" post (R.checkEntries c.D d s.tstack s (fun t => t.vStart != v && t.topDepth ≥ d)) ++
    ofStrs c tern "r.disj_post" post (R.disjoint s) ++
    rfr.flatMap (fun (p, dp, n0) =>
      let bot := sPre.tstack.drop (sPre.tstack.length - n0)
      let bot' := s.tstack.drop (s.tstack.length - n0)
      (if bot'.map R.showT != bot.map R.showT then
        [⟨c.seed, tern, "r.base_post", s!"{post} p={p} dp={dp} bottom {n0} changed"⟩] else []) ++
      ofStrs c tern "r.base_post" s!"{post} p={p} dp={dp} n0={n0}" (bot'.flatMap fun t => R.settledEntry c.D p dp s t "post")))
  let vs := vs ++ (if d ≤ lowval && !(s.tstack == s₂.tstack && hv' == hv && s.stackDir == sPre.stackDir) then
    [⟨c.seed, tern, "st.boundary", s!"{post} tstack/hasVert/stackDir changed"⟩] else [])
  (hv', s, eb, visited, vs)
end

def iForest (c : Ctx) (forest : List DfsTree) (s : WalkState) : WalkState × List Viol :=
  let (s, vs, _, _, _) := forest.foldl (fun (s, vs, visited, n, prev) t =>
    let (s, _, visited, vs') := iTree c prev [] n t 0 [] visited s
    let s := ((popTstack >>= fun top => modifyItem rootItem fun it => { it with ch := it.ch ++ top.spans.2 }).run s).2
    let prev := prev ++ [t]
    let vs' := vs' ++ (if s.tstack == [] then [] else [⟨c.seed, s.ternarize, "st.root_pop", s!"tstack not empty after tree {t.v}"⟩])
    let vs' := vs' ++ ofRanges (Ranges.checkPieceInv c.seed "root" 0 s false)
    let vs' := vs' ++ ofStrs c s.ternarize "st.items_root" s!"after tree {t.v}" (St.stItemsB c.g s (simBlocks c.g prev [] []))
    (s, vs ++ vs', visited, n + t.edgePostorder.length, prev))
    (s, ([] : List Viol), ([] : List Nat), 0, ([] : List DfsTree))
  (s, vs)

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

/-- One run: whether the instrumented walk equals `g.walk`, whether `g` is a block (R fields
evaluated), the violations. -/
def runGraph (seed : Nat) (g : Graph) (tern : Bool) : Bool × Bool × List Viol :=
  let f := g.dfsForest [] []
  let c : Ctx := ⟨seed, g, edgePostorderForest f, R.mkDfs g f, R.isBlock g && g.ne ≥ 2⟩
  let (s, vs) := iForest c f (WalkState.init g tern)
  let vs := vs ++ ofRanges (Ranges.checkFinal seed s)
  let vs := vs ++ Planar.finishB seed g tern f
  let ref := g.walk tern f
  (s.items.toList.map (fun it => (it.type, it.vs, it.ch)) == ref.items.toList.map (fun it => (it.type, it.vs, it.ch)) &&
    s.tstack.length == ref.tstack.length, c.block, vs)
def runSeed (seed : Nat) (tern : Bool) : Bool × Bool × List Viol := runGraph seed (randGraph seed) tern

/-- Small random multigraphs (for minimal counterexamples): `nv ∈ 2..5`, `ne ∈ 1..8`. -/
def smallGraph (seed : Nat) : Graph := Id.run do
  let mut x := lcg (seed + 777)
  let nv := 2 + (x >>> 33) % 4
  x := lcg x
  let ne := 1 + (x >>> 33) % 8
  let mut es : Array (Nat × Nat) := #[]
  for _ in List.range ne do
    x := lcg x
    let a := (x >>> 33) % nv
    x := lcg x
    let b := (x >>> 33) % nv
    es := es.push (a, b)
  return ⟨nv, es⟩

/-- `check_walkinv small <field> <n>`: the smallest (by `ne`, then `nv`) small random multigraph
among `smallGraph 0..n-1` violating `field` (both modes). -/
def smallest (field : String) (n : Nat) : IO UInt32 := do
  let mut best : Option (Nat × Nat × Nat × Bool × Viol) := none
  for seed in List.range n do
    let g := smallGraph seed
    for tern in [false, true] do
      let (_, _, vs) := runGraph seed g tern
      for v in vs do
        if v.field == field then
          match best with
          | some (ne, nv, _, _, _) =>
            if g.ne < ne || (g.ne == ne && g.nv < nv) then best := some (g.ne, g.nv, seed, tern, v)
          | none => best := some (g.ne, g.nv, seed, tern, v)
  match best with
  | some (_, _, seed, tern, v) =>
    IO.println s!"smallest {field}: seed={seed} tern={tern} nv={(smallGraph seed).nv} edges={repr (smallGraph seed).edges}"
    IO.println s!"  {repr v}"
    return 0
  | none => IO.println s!"no violation of {field} in {n} small graphs"; return 1

/-- `check_walkinv minimize <field> <seed>`: greedily delete edges of `randGraph seed` (then drop
unused trailing vertices) while some run (either mode) still violates `field`. -/
def minimize (field : String) (seed : Nat) : IO UInt32 := do
  let viol (g : Graph) : Option Viol := Id.run do
    for tern in [false, true] do
      let (_, _, vs) := runGraph seed g tern
      match vs.find? (·.field.startsWith field) with
      | some v => return some v
      | none => pure ()
    return none
  let mut g := randGraph seed
  if (viol g).isNone then IO.println "no violation"; return 1
  let mut changed := true
  while changed do
    changed := false
    for e in List.range g.edges.size do
      if e < g.edges.size then
        let g' : Graph := ⟨g.nv, g.edges.eraseIdx! e⟩
        if (viol g').isSome then g := g'; changed := true
  let mut nv := g.nv
  while nv > 1 && !(g.edges.any fun (a, b) => a == nv - 1 || b == nv - 1) do nv := nv - 1
  g := ⟨nv, g.edges⟩
  IO.println s!"minimal {field} from seed {seed}: nv={g.nv} edges={repr g.edges} block={R.isBlock g}"
  for tern in [false, true] do
    let (ok, _, vs) := runGraph seed g tern
    IO.println s!"  tern={tern} matches={ok} violations:"
    for v in vs do IO.println s!"    {v.field}: {v.info}"
  return 0

structure FieldStat where
  count : Nat
  minSeed : Nat
  minTern : Bool
  sample : String

def summarize (lo hi : Nat) (all : Bool := false) : IO UInt32 := do
  let mut stats : List (String × FieldStat) := []
  let mut okAll := true
  let mut mismatches : List (Nat × Bool) := []
  let mut runs := 0
  let mut blocks := 0
  for seed in List.range' lo (hi - lo) do
    for tern in [false, true] do
      let (ok, block, vs) := runSeed seed tern
      runs := runs + 1
      if block then blocks := blocks + 1
      if !ok then
        okAll := false
        mismatches := mismatches ++ [(seed, tern)]
      for v in vs do
        if all then IO.println s!"{repr v} g={repr (randGraph seed).edges}"
        match stats.find? (·.1 == v.field) with
        | some (_, st) =>
          let st := { st with count := st.count + 1 }
          let st := if seed < st.minSeed then { st with minSeed := seed, minTern := tern, sample := s!"{repr v} g={repr (randGraph seed).edges}" } else st
          stats := stats.map fun (k, x) => if k == v.field then (k, st) else (k, x)
        | none => stats := stats ++ [(v.field, ⟨1, seed, tern, s!"{repr v} g={repr (randGraph seed).edges}"⟩)]
  IO.println s!"instrumentation matches library walk: {okAll}{if okAll then "" else s!" mismatches={mismatches.take 5}"}"
  let total := stats.foldl (fun a (_, st) => a + st.count) 0
  IO.println s!"seeds {lo}..{hi - 1} both modes: runs={runs} block runs (R fields)={blocks} fields violated={stats.length} violations={total}"
  for (k, st) in stats do
    IO.println s!"FIELD {k}: count={st.count} smallest seed={st.minSeed} tern={st.minTern}"
    IO.println s!"  {st.sample}"
  return if okAll && stats.isEmpty then 0 else 1

end WalkInvCheck

def main (args : List String) : IO UInt32 := do
  if args[0]? == some "minimize" then
    return ← WalkInvCheck.minimize (args.getD 1 "") ((args[2]?.bind String.toNat?).getD 0)
  if args[0]? == some "small" then
    return ← WalkInvCheck.smallest (args.getD 1 "") ((args[2]?.bind String.toNat?).getD 20000)
  let lo := args[0]?.bind String.toNat? |>.getD 0
  let hi := args[1]?.bind String.toNat? |>.getD 3001
  WalkInvCheck.summarize lo hi (args.getD 2 "" == "all")
