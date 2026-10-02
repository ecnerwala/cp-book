import Spqr
import Spqr.StRef

open Spqr

/-! Differential test of `walk_st'`, `walk_vsOriented` and `refBlocks_st`: the walk's S / P / R child
lists equal `restrictCh` of `refOrder`, every S / P / R item lies in one reference block with its
`vs` and its non-V children's `vs` oriented along the block's sequence, and every reference block
is st-numbered. -/

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
  let items := (g.walk (tern != 0) forest).items
  let blocks := refBlocks g forest
  let order := refOrder g forest
  let mut bad := 0
  let mut n := 0
  for b in blocks do
    if ¬ b.St g then
      bad := bad + 1
      IO.println s!"BLOCK not st: root {repr b.root} seq {b.seq g} edges {b.edges g}"
  for i in [0:items.size] do
    match Items.type items i with
    | .S | .P | .R =>
      n := n + 1
      let r := restrictCh items items.size order i
      if r ≠ Items.ch items i then
        bad := bad + 1
        IO.println s!"MISMATCH item {i}: walk {Items.ch items i} ref {r}"
      let lv := Items.leaves items items.size i
      let ok := blocks.any fun b => lv.all (· ∈ b.items) && decide (Oriented (b.seq g) (Items.vs items i)) &&
        (Items.ch items i).all fun c => Items.type items c = .V || decide (Oriented (b.seq g) (Items.vs items c))
      if !ok then
        bad := bad + 1
        IO.println s!"ORIENT item {i}: vs {repr (Items.vs items i)} ch {(Items.ch items i).map fun c => (c, repr (Items.vs items c))} leaves {lv}"
    | _ => pure ()
  IO.println s!"items {n} blocks {blocks.length} bad {bad}"
