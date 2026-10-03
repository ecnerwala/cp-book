import Spqr.Build
import Spqr.Spec
import Spqr.PlanarEmbedSteps

open Spqr

def checkCapFace (g : Graph) (t : PlanarSpqrTree) (j : Nat)
    (s : PlanarSpqrTree.EmbedState) : Bool := Id.run do
  let p := t.pieceBelow g j
  let o := s.outerE[j]!
  let mut adj := #[]
  for l in [0:4 * p.ves.length] do
    let q := 4 * p.ves[l / 4]! + l % 4
    let r := if o[0]! == some q then o[1]!
      else if o[1]! == some q then o[0]!
      else if o[2]! == some q then o[3]!
      else if o[3]! == some q then o[2]!
      else s.rotAdj[q]!
    adj := adj.push (r.bind p.loc)
  match o[0]!.bind p.loc, o[2]!.bind p.loc with
  | some a, some b =>
    let ρ : RotationSystem := ⟨adj⟩
    let mut q := a
    for _ in [0:ρ.size] do
      if q == b then return true
      q := (ρ.faceStep q).getD q
    return false
  | _, _ => return p.ves.isEmpty

def edgeIn (t : SpqrTree) (i e : Nat) : Bool :=
  match t.edgeIndex[e]! with
  | none => false
  | some j => i ≤ j && j < t.subtreeEnd[i]!

def incident (g : Graph) (v e : Nat) : Bool :=
  g.edges[e]!.1 == v || g.edges[e]!.2 == v

def touches (t : SpqrTree) (g : Graph) (i v : Nat) : Bool :=
  (List.range g.ne).any fun e => edgeIn t i e && incident g v e

def checkOuter (g : Graph) (t : PlanarSpqrTree) (i : Nat)
    (s : PlanarSpqrTree.EmbedState) : IO Nat := do
  let mut bad := 0
  for j in [0:t.size] do
    if let some ne := t.toSpqrTree.capNe j then
      if i <= j && ((t.toSpqrTree.parent j).all (· < i)) && !checkCapFace g t j s then
        bad := bad + 1
        IO.println s!"cap_cofacial: i={i} j={j}"
      if let some (u, v) := t.toSpqrTree.neOrig ne then
        for k in [0:4] do
          if let some q := s.outerE[j]![k]! then
            if QE.vert g.edges.toList q != some (if k < 2 then u else v) then
              bad := bad + 1
              IO.println s!"outer_cap_vertex: j={j} k={k} q={q}"
          else if i <= j && !(t.edgesBelow j).isEmpty then
            bad := bad + 1
            IO.println s!"outer_cap_present: i={i} j={j} k={k}"
    if t.toSpqrTree.type j == .V then
      if i <= j && !(t.edgesBelow j).isEmpty && !(s.outerE[j]!.any Option.isSome) then
        bad := bad + 1
        IO.println s!"outer_vertex_present: i={i} j={j}"
      for k in [0:s.outerE[j]!.size] do
        if let some q := s.outerE[j]![k]! then
          if k >= 2 || QE.vert g.edges.toList q != t.origId[j]! then
            bad := bad + 1
            IO.println s!"outer_vertex: j={j} k={k} q={q}"
    if i <= j && ((t.toSpqrTree.parent j).map t.toSpqrTree.type == some NodeType.V) &&
        !(t.edgesBelow j).isEmpty && !(s.outerE[j]!.any Option.isSome) then
      bad := bad + 1
      IO.println s!"outer_present: i={i} j={j}"
    if s.outerE[j]!.size != 4 then
      bad := bad + 1
      IO.println s!"outer_row_size: j={j}"
    for k in [0:s.outerE[j]!.size] do
      if s.outerE[j]![k]!.isSome then
        if (s.outerE[j]![k]!).map (· % 2) != some (k % 2) then
          bad := bad + 1
          IO.println s!"outer_dir: j={j} k={k}"
        let parentType := (t.toSpqrTree.parent j).map t.toSpqrTree.type
        if k >= 4 || t.toSpqrTree.type j == .F ||
            ((parentType == some NodeType.F || parentType == some NodeType.V) && k >= 2) then
          bad := bad + 1
          IO.println s!"outer_slots: j={j} k={k}"
        if let some p := t.toSpqrTree.parent j then
          if t.toSpqrTree.type p == .V then
            if let some v := t.origId[p]! then
              if (s.outerE[j]![k]!).bind (QE.vert g.edges.toList) != some v then
                bad := bad + 1
                IO.println s!"outer_at_vertex: j={j} k={k} parent={p} vertex={v}"
  return bad

def check (g : Graph) (tern : Bool) (vo eo : List Nat) : IO Nat := do
  let t := relabelTree g (g.walk tern (g.dfsForest vo eo)).items
  let mut bad := 0
  for a in t.children 0 do
    for b in t.children 0 do
      for v in [0:g.nv] do
        if a != b && touches t g a v && touches t g b v then
          bad := bad + 1
          IO.println s!"root_disjoint: a={a} b={b} v={v}"
  for i in [0:t.size] do
    if t.type i == .V then
      for a in t.children i do
        if !((List.range g.ne).any fun e => edgeIn t a e) then
          bad := bad + 1
          IO.println s!"v_nonempty: i={i} a={a}"
        for b in t.children i do
          for v in [0:g.nv] do
            if a != b && touches t g a v && touches t g b v && t.origId[i]! != some v then
              bad := bad + 1
              IO.println s!"v_attach: i={i} a={a} b={b} v={v}"
    if let some ne := t.capNe i then
      if (t.neOrig ne).isNone then
        bad := bad + 1
        IO.println s!"cap_orig: i={i} ne={ne}"
    if t.hasCap i && t.type i != .I && t.type i != .O && !((List.range g.ne).any fun e => edgeIn t i e) then
      bad := bad + 1
      IO.println s!"cap_nonempty: i={i} type={repr (t.type i)} children={t.children i} parent={t.parent i}"
    if t.type i == .Q then
      match t.children i, t.origId[i]! with
      | [_, w], some e =>
        if t.type w == .V then
          for v in [0:g.nv] do
            for e' in [0:g.ne] do
              if touches t g i v && incident g v e' && !edgeIn t i e' && !incident g v e then
                bad := bad + 1
                IO.println s!"q_root_attach: i={i} e={e} v={v} outside={e'}"
      | _, _ => pure ()
      let lower (e : Nat) : Nat :=
        if t.edgeFlipped[e]! then (g.edges[e]!).1 else (g.edges[e]!).2
      let orient (e : Nat) : Nat × Nat :=
        if t.edgeFlipped[e]! then ((g.edges[e]!).2, (g.edges[e]!).1) else g.edges[e]!
      let shapeOk := match t.children i with
        | [] => true
        | [c] => t.type c == .O
        | [c, w] => t.type c != .V && t.type c != .O && t.hasCap c && t.type w == .V
        | _ => false
      if !shapeOk then
        bad := bad + 1
        IO.println s!"q_shape: i={i} children={t.children i}"
      if t.children i == [] then
        if let some p := t.parent i then
          if t.type p == .F || t.type p == .V then
            bad := bad + 1
            IO.println s!"q_leaf_parent: i={i} p={p}"
      if let some e := t.origId[i]! then
        let upper := if t.edgeFlipped[e]! then (g.edges[e]!).2 else (g.edges[e]!).1
        let isLoop := (g.edges[e]!).1 == (g.edges[e]!).2
        if isLoop != (t.children i).any (fun c => t.type c == .O) then
          bad := bad + 1
          IO.println s!"q_loop: i={i} e={e}"
        if let some p := t.parent i then
          if t.type p == .V && t.origId[p]! != some upper then
            bad := bad + 1
            IO.println s!"q_upper: i={i} p={p}"
        match t.children i with
        | [c, w] =>
          if t.origId[w]! != some (lower e) then
            bad := bad + 1
            IO.println s!"q_lower: i={i} w={w}"
          for v in [0:g.nv] do
            if (touches t g c v || incident g v e) && touches t g w v && v != lower e then
              bad := bad + 1
              IO.println s!"q_lower_attach: i={i} c={c} w={w} v={v}"
        | _ => pure ()
        for c in t.children i do
          if let some ne := t.capNe c then
            if t.neOrig ne != some (orient e) then
              bad := bad + 1
              IO.println s!"q_cap_orient: i={i} c={c}"
  for i in [0:t.size] do
    let (neSt, neEn) := t.neRange i
    for ne in [neSt:neEn] do
      if let some tw := t.twin ne then
        if t.neOrig ne != t.neOrig tw then
          bad := bad + 1
          IO.println s!"twin_orient: i={i} ne={ne} tw={tw}"
    if t.type i == .S || t.type i == .P || t.type i == .R then
      for ne in [neSt + (if t.hasCap i then 1 else 0):neEn] do
        match t.twin ne with
        | none =>
          bad := bad + 1
          IO.println s!"node_twin: i={i} ne={ne}"
        | some tw =>
          match t.nodeOfNe tw with
          | some c =>
            if !(t.children i).contains c || t.capNe c != some tw then
              bad := bad + 1
              IO.println s!"node_twin_child: i={i} ne={ne} tw={tw} c={c}"
          | none =>
            bad := bad + 1
            IO.println s!"node_twin_node: i={i} ne={ne} tw={tw}"
  for i in [0:t.size] do
    if t.type i == .S || t.type i == .P || t.type i == .R then
      let (neSt, neEn) := t.neRange i
      let (nvSt, nvEn) := t.nvRange i
      for nv in [nvSt:nvEn] do
        for nv' in [nvSt:nv] do
          if t.nvOrig nv == t.nvOrig nv' then
            bad := bad + 1
            IO.println s!"nv_orig_inj: i={i} nv={nv} nv'={nv'}"
      let nvInc (c nv : Nat) : Bool :=
        (t.nodeVerts[nv]!).vert == c ||
          (List.range (neEn - neSt)).any fun k =>
            let ne := neSt + k
            t.twin ne == t.capNe c && t.capNe c != none &&
              ((t.nodeEdges[ne]!).nvs.1 == nv || (t.nodeEdges[ne]!).nvs.2 == nv)
      for c in t.children i do
        if t.type c != .V && (!t.hasCap c || t.type c == .I || t.type c == .O) then
          bad := bad + 1
          IO.println s!"node_child_cap: i={i} c={c} type={repr (t.type c)}"
        for w in [0:g.nv] do
          if touches t g c w then
            for nv in [nvSt:nvEn] do
              if t.nvOrig nv == some w && !nvInc c nv then
                bad := bad + 1
                IO.println s!"node_touch: i={i} c={c} w={w} nv={nv}"
      for a in t.children i do
        for b in t.children i do
          if a != b then
            for w in [0:g.nv] do
              if touches t g a w && touches t g b w then
                if !((List.range (nvEn - nvSt)).any fun k => t.nvOrig (nvSt + k) == some w) then
                  bad := bad + 1
                  IO.println s!"node_attach: i={i} a={a} b={b} w={w}"
  let pt := g.planarSpqrTree tern vo eo
  for i in [0:pt.size] do
    let ty := pt.toSpqrTree.type i
    if (ty == .S || ty == .P || ty == .R) && pt.nodePlanar[i]! then
      let (neSt, neEn) := pt.toSpqrTree.neRange i
      let (nvSt, nvEn) := pt.toSpqrTree.nvRange i
      let cap := pt.toSpqrTree.nodeEdges[neSt]!
      for nv in [nvSt:nvEn] do
        let mut cnt := 0
        let mut cntRev := 0
        for ta in [4 * (neSt + 1):4 * neEn] do
          if ta % 4 == 2 && (pt.toSpqrTree.nodeEdges[ta / 4]!).nvs.2 == nv then
            if let some tb := pt.neRotAdj[ta]! then
              if tb % 4 == 1 then
                if ta < tb then cnt := cnt + 1 else cntRev := cntRev + 1
        let v := (pt.toSpqrTree.nodeVerts[nv]!).vert
        let capEnd := cap.nvs.1 == nv || cap.nvs.2 == nv
        let hasPieces := !(pt.edgesBelow v).isEmpty
        let isChild := pt.toSpqrTree.parent v == some i
        if capEnd && (cnt + cntRev != 0) then
          bad := bad + 1
          IO.println s!"corner_capend: i={i} nv={nv} cnt={cnt} rev={cntRev}"
        if !capEnd && (cnt != 1 || cntRev != 0) then
          bad := bad + 1
          IO.println s!"corner_inner: i={i} ty={repr ty} nv={nv} cnt={cnt} rev={cntRev} pieces={hasPieces}"
        if !capEnd && !isChild then
          bad := bad + 1
          IO.println s!"corner_vitem_not_child: i={i} nv={nv} v={v}"
        if capEnd && isChild then
          bad := bad + 1
          IO.println s!"corner_capend_child: i={i} nv={nv} v={v} pieces={hasPieces}"
      if let some p := pt.toSpqrTree.parent i then
        if pt.toSpqrTree.type p == .F || pt.toSpqrTree.type p == .V then
          bad := bad + 1
          IO.println s!"node_parent: i={i} p={p}"
      for c in pt.toSpqrTree.children i do
        if pt.toSpqrTree.type c == .V then
          let nvs := (List.range (nvEn - nvSt)).filter fun k =>
            (pt.toSpqrTree.nodeVerts[nvSt + k]!).vert == c
          if nvs.length != 1 then
            bad := bad + 1
            IO.println s!"corner_vchild_nv: i={i} c={c} n={nvs.length}"
      for nv in [nvSt:nvEn] do
        for nv' in [nvSt:nv] do
          if (pt.toSpqrTree.nodeVerts[nv]!).vert == (pt.toSpqrTree.nodeVerts[nv']!).vert then
            bad := bad + 1
            IO.println s!"corner_nv_inj: i={i} nv={nv} nv'={nv'}"
  let mut s := pt.initState
  bad := bad + (← checkOuter g pt pt.size s)
  if pt.nodePlanar.all id then
    for i in (List.range pt.size).reverse do
      s := ((pt.embedItem i).run s).2
      bad := bad + (← checkOuter g pt i s)
  return bad

def main : IO UInt32 := do
  let input ← (← IO.getStdin).readToEnd
  let toks := (input.splitOn " ").flatMap (·.splitOn "\n") |>.filter (· ≠ "") |>.map String.toNat!
  let toks := toks.toArray
  let nv := toks[0]!; let ne := toks[1]!; let tern := toks[2]! != 0
  let edges := ((List.range ne).map fun e => (toks[3 + 2 * e]!, toks[4 + 2 * e]!)).toArray
  let p := 3 + 2 * ne
  let k := toks[p]!
  let vo := (List.range k).map fun j => toks[p + 1 + j]!
  let p := p + 1 + k
  let l := toks[p]!
  let eo := (List.range l).map fun j => toks[p + 1 + j]!
  let bad ← check ⟨nv, edges⟩ tern vo eo
  IO.println s!"bad {bad}"
  return if bad == 0 then 0 else 1
