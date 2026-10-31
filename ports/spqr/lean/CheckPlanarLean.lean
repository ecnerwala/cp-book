import Spqr

open Spqr

/-- Reads a `gen.py` instance, builds the planar SPQR tree and checks the planarity output
against the specification: every planar S/P/R node's local rotation system is a planar
embedding of its skeleton, and the glued embedding (if any) is a planar embedding of `g`. -/
def main : IO UInt32 := do
  let input ← (← IO.getStdin).readToEnd
  let toks := (input.splitOn " ").flatMap (·.splitOn "\n") |>.filter (· ≠ "") |>.map String.toNat!
  let toks := toks.toArray
  let mut p := 0
  let (nv, ne, tern) := (toks[0]!, toks[1]!, toks[2]!)
  p := 3
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
  let t := g.planarSpqrTree (tern != 0) vertOrder edgeOrder
  let mut ok := true
  let mut nodes := 0
  for i in [0 : t.size] do
    let ty := t.toSpqrTree.type i
    if ty == .S || ty == .P || ty == .R then
      if t.isPlanar i then
        nodes := nodes + 1
        unless IsPlanarEmbedding (t.localSkeleton i) (t.toSpqrTree.nVerts i) (t.nodeRot i) do
          IO.println s!"FAIL node {i}"
          ok := false
      else
        unless ty == .R && (t.nodeRot i).rotAdj.all (·.isNone) do
          IO.println s!"FAIL nonplanar node {i} is not an R node with all entries unset"
          ok := false
  match t.planarEmbed with
  | some rs =>
    unless IsPlanarEmbedding edges.toList nv rs do
      IO.println "FAIL glued embedding"
      ok := false
    unless t.nodePlanar.all id do
      IO.println "FAIL embedding despite nonplanar node"
      ok := false
    IO.println s!"ok nodes={nodes} embedded"
  | none =>
    if t.nodePlanar.all id then
      IO.println "FAIL no embedding although all nodes are planar"
      ok := false
    IO.println s!"ok nodes={nodes} nonplanar"
  return if ok then 0 else 1
