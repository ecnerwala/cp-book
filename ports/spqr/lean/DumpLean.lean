import Spqr

open Spqr

def optStr : Option Nat → String
  | none => "-1"
  | some x => toString x

def line (name : String) (xs : List String) : String :=
  name ++ ":" ++ String.join (xs.map (" " ++ ·)) ++ "\n"

def NodeType.char : NodeType → String
  | .F => "F" | .V => "V" | .Q => "Q" | .I => "I" | .O => "O" | .S => "S" | .P => "P" | .R => "R"

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
  let ear := (← IO.getEnv "SPQR_EAR").isSome
  let t := if ear then relabelTree g (g.walkEar (tern != 0) (g.dfsForest vertOrder edgeOrder)).items
    else g.spqrTree (tern != 0) vertOrder edgeOrder
  let nats (xs : Array Nat) := xs.toList.map toString
  let opts (xs : Array (Option Nat)) := xs.toList.map optStr
  let out := String.join [
    line "vert_index" (opts t.vertIndex),
    line "edge_index" (opts t.edgeIndex),
    line "edge_flipped" (t.edgeFlipped.toList.map fun b => if b then "1" else "0"),
    line "par" (opts t.par),
    line "subtree_end" (nats t.subtreeEnd),
    line "types" (t.types.toList.map NodeType.char),
    line "orig_id" (opts t.origId),
    line "ch.bounds" (nats t.chBounds),
    line "ch.dat" (nats t.chDat),
    line "node_nvs.bounds" (nats t.nvBounds),
    line "node_verts" (t.nodeVerts.toList.map fun x => s!"{x.node},{x.vert}"),
    line "vert_par_nv" (opts t.vertParNv),
    line "node_nes.bounds" (nats t.neBounds),
    line "node_edges" (t.nodeEdges.toList.map fun x => s!"{x.node},{optStr x.twin},{x.nvs.1},{x.nvs.2}"),
    line "node_adj.bounds" (nats t.adjBounds),
    line "node_adj.dat" (t.adjDat.toList.map fun x => s!"{x.ne},{x.destNv}")]
  (← IO.getStdout).putStr out
