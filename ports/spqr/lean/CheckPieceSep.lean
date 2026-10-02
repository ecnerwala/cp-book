import Spqr.Build
import Spqr.Spec
import Spqr.PlanarEmbedSteps

open Spqr

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
  let pt := g.planarSpqrTree tern vo eo
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
