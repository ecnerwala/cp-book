import Spqr
import Spqr.StRef

open Spqr

/-! Differential test of `walk_st'`: the walk's S / P / R child lists equal `restrictCh` of `refOrder`. -/

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
  let order := refOrder g forest
  let mut bad := 0
  let mut n := 0
  for i in [0:items.size] do
    match Items.type items i with
    | .S | .P | .R =>
      n := n + 1
      let r := restrictCh items items.size order i
      if r ≠ Items.ch items i then
        bad := bad + 1
        IO.println s!"MISMATCH item {i}: walk {Items.ch items i} ref {r}"
    | _ => pure ()
  IO.println s!"items {n} bad {bad}"
