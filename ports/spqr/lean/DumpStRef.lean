import Spqr

open Spqr

/-! Dump the lowval-sorted DFS forest and the walk's items, for prototyping the st-order
reference (`Spqr/StRef.lean`). -/

def optStr : Option Nat → String
  | none => "-1"
  | some x => toString x

def NodeType.char : NodeType → String
  | .F => "F" | .V => "V" | .Q => "Q" | .I => "I" | .O => "O" | .S => "S" | .P => "P" | .R => "R"

mutual
partial def dumpTree (t : DfsTree) (d : Nat) : IO Unit := do
  match t with
  | .node v outs =>
    IO.println s!"N {v} {d} {outs.length}"
    for o in outs do dumpOut o d
partial def dumpOut (o : DfsOut) (d : Nat) : IO Unit := do
  IO.println s!"O {o.e} {o.cls.lowval d} {if o.cls.isTree then 1 else 0} {if o.cls.isType1 then 1 else 0} {o.dest}"
  match o with
  | .tree _ _ child => dumpTree child (d + 1)
  | .back .. => pure ()
end

def main : IO Unit := do
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
  IO.println s!"G {nv} {ne} {tern}"
  for (u, v) in edges do IO.println s!"E {u} {v}"
  for t in forest do
    IO.println "T"
    dumpTree t 0
  let items := (g.walk (tern != 0) forest).items
  for i in [0:items.size] do
    let it := items[i]!
    let ch := String.intercalate " " (it.ch.map toString)
    IO.println s!"I {i} {NodeType.char it.type} {optStr it.vs.1} {optStr it.vs.2} {ch}"
