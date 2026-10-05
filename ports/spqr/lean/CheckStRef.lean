import Spqr
import Spqr.StRef

open Spqr

/-! Differential test of `walk_st'`, `walk_vsOriented` and `refBlocks_st`: the walk's S / P / R child
lists equal `restrictCh` of `refOrder`, every S / P / R item lies in one reference block with its
`vs` and its non-V children's `vs` oriented along the block's sequence, and every reference block
is st-numbered. -/

instance (p q : Nat × Nat) : Decidable (Items.PairEq p q) := by unfold Items.PairEq; infer_instance

/-- The descendants-or-self of `i` (`Items.Below i`), by fuel. -/
def descOf (items : Items) : Nat → ItemId → List ItemId
  | 0, i => [i]
  | fuel + 1, i => i :: (Items.ch items i).flatMap (descOf items fuel)

/-- Boolean `Items.StItem`. -/
def stItemB (items : Items) (i : ItemId) : Bool :=
  match Items.vs items i with
  | (some s, some t) =>
    let xs := Items.vertList items i
    let es := (s, t) :: Items.virtualEdges items i
    xs.Nodup && es.all (fun p => p.1 ∈ xs && p.2 ∈ xs && p.1 ≠ p.2) &&
    xs.all (fun x => xs.head? == some x || xs.getLast? == some x ||
      (es.any (fun p => (p.1 == x && xs.idxOf p.2 < xs.idxOf x) || (p.2 == x && xs.idxOf p.1 < xs.idxOf x)) &&
       es.any (fun p => (p.1 == x && xs.idxOf x < xs.idxOf p.2) || (p.2 == x && xs.idxOf x < xs.idxOf p.1)))) &&
    (Items.virtualEdges items i).all fun p => xs.idxOf p.1 < xs.idxOf p.2
  | _ => false

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
  if order.eraseDups.length != order.length then
    bad := bad + 1
    IO.println s!"ORDER not nodup: {order}"
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
      let ds := descOf items items.size i
      let ok := blocks.any fun b => lv.all (· ∈ b.items) &&
        ((List.range g.ne).all fun e =>
          let x := edgeItem g e
          !(x ∈ b.items || b.root.any fun r => decide (Items.PairEq g.edges[e]! r)) || !(x ∈ ds) || x ∈ lv) &&
        decide (Oriented (b.seq g) (Items.vs items i)) &&
        (match Items.vs items i with
          | (some s, some t) => (Items.ch items i).all fun c => Items.type items c ≠ .V ||
              (decide (Precedes (b.seq g) s (c - 1)) && decide (Precedes (b.seq g) (c - 1) t))
          | _ => false) &&
        ((Items.ch items i).all fun c => Items.type items c = .V || decide (Oriented (b.seq g) (Items.vs items c))) &&
        ((Items.ch items i).all fun c => Items.type items c = .V ||
          match Items.vs items c with
          | (some u, some v) =>
            let dc := descOf items items.size c
            (List.range g.ne).all fun e =>
              !(edgeItem g e ∈ b.items || b.root.any fun r => decide (Items.PairEq g.edges[e]! r)) ||
              !(edgeItem g e ∈ dc) ||
              [(g.edges[e]!).1, (g.edges[e]!).2].all fun y =>
                (u == y || decide (Precedes (b.seq g) u y)) && (y == v || decide (Precedes (b.seq g) y v))
          | _ => false)
      if !stItemB items i then
        bad := bad + 1
        IO.println s!"NOT ST item {i}: vs {repr (Items.vs items i)} vertList {Items.vertList items i} ve {Items.virtualEdges items i}"
      if !ok then
        bad := bad + 1
        IO.println s!"ORIENT item {i}: vs {repr (Items.vs items i)} ch {(Items.ch items i).map fun c => (c, repr (Items.vs items c))} leaves {lv}"
    | _ => pure ()
  IO.println s!"items {n} blocks {blocks.length} bad {bad}"
