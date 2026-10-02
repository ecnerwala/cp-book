import Spqr.Build
import Spqr.Spec

open Spqr

def edgeIn (t : SpqrTree) (i e : Nat) : Bool :=
  match t.edgeIndex[e]! with
  | none => false
  | some j => i ≤ j && j < t.subtreeEnd[i]!

def incident (g : Graph) (v e : Nat) : Bool :=
  g.edges[e]!.1 == v || g.edges[e]!.2 == v

def touches (t : SpqrTree) (g : Graph) (i v : Nat) : Bool :=
  (List.range g.ne).any fun e => edgeIn t i e && incident g v e

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
