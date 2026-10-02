import Spqr
import Spqr.ItemSpec

open Spqr

/-! Empirical check of the item facts beyond the former `Items.WF` (`q_root`, `o_parent`, `s_order`)
on the walk's output; input format of `gen.py`. -/

def pairEq (p q : Nat × Nat) : Bool := p == q || p == (q.2, q.1)

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
  let mut bad := 0
  for e in [0:ne] do
    let q := edgeItem g e
    match Items.ch items q, (Items.vs items q).1 with
    | [_], _ =>
      if edges[e]!.1 != edges[e]!.2 then
        bad := bad + 1; IO.println s!"q_root loop: e={e} ch={Items.ch items q} edge={edges[e]!}"
    | [_, w], some u =>
      if w < 1 || !pairEq (u, w - 1) edges[e]! then
        bad := bad + 1; IO.println s!"q_root pair: e={e} ch={Items.ch items q} vs={Items.vs items q} edge={edges[e]!}"
    | [_, _], none =>
      bad := bad + 1; IO.println s!"q_root pair (no vs.1): e={e}"
    | _, _ => pure ()
  for i in [0:items.size] do
    for c in Items.ch items i do
      if Items.type items c == .O then
        if Items.type items i != .Q || Items.ch items i != [c] then
          bad := bad + 1; IO.println s!"o_parent: p={i} type={repr (Items.type items i)} ch={Items.ch items i}"
    if Items.type items i == .S then
      let xs := ((Items.ch items i).filter fun c => Items.type items c == .V).map (· - 1)
      match Items.vs items i with
      | (some u, some v) =>
        if Items.virtualEdges items i != List.zip (u :: xs) (xs ++ [v]) then
          bad := bad + 1
          IO.println s!"s_order: i={i} vs={Items.vs items i} xs={xs} virt={Items.virtualEdges items i}"
      | _ => bad := bad + 1; IO.println s!"s_order vs: i={i} vs={Items.vs items i}"
  IO.println s!"items {items.size} bad {bad}"
