import WalkInvCheck.Common
/-!
# Open-stack content (`content.*`)

`EntryContent` (`Spqr/WalkContent.lean`) evaluated on every open entry at every site: the span items
are S/P/R/Q/V roots (Q leaves), non-V ones have two distinct terminals among the entry's terminals
(`Term'`: `vStart`, the open path above `topDepth`, the `vStart`s of the entries above) or its V
items, and a vertex item belongs to the entry iff it is a non-terminal vertex all of whose edges lie
in the entry or the entries above it, outside every single span item.
-/
namespace WalkInvCheck.Content
open Spqr WalkM

partial def edgesBelow (s : WalkState) (i : ItemId) : List Nat :=
  (if 1 + s.g.nv ≤ i ∧ i < 1 + s.g.nv + s.g.ne then [i - 1 - s.g.nv] else []) ++
    (s.items[i]!.ch).flatMap (edgesBelow s)
def spanItems (t : TEntry) : List ItemId := t.spans.1 ++ t.spans.2
def entryEdges (s : WalkState) (t : TEntry) : List Nat := (spanItems t).flatMap (edgesBelow s)
def inc (s : WalkState) (e v : Nat) : Bool := s.g.edges[e]!.1 == v || s.g.edges[e]!.2 == v
def touches (s : WalkState) (E : List Nat) (v : Nat) : Bool := E.any (inc s · v)
def interior (s : WalkState) (E : List Nat) (v : Nat) : Bool :=
  (List.range s.g.ne).all fun e => !inc s e v || E.contains e
def showT (t : TEntry) : String := s!"({t.vStart},{t.topDepth},{t.firstIdx},{t.spans})"

/-- `EntryContent` on every open entry of `s`. -/
def check (seed : Nat) (site : String) (D : Nat) (s : WalkState) : List Viol := Id.run do
  let ty : Nat → NodeType := fun i => s.items[i]!.type
  let ch : Nat → List ItemId := fun i => s.items[i]!.ch
  let vs : Nat → Option Nat × Option Nat := fun i => s.items[i]!.vs
  let mut out : List Viol := []
  let bad (k info : String) : Viol :=
    ⟨seed, s.ternarize, s!"content.{k}", s!"{site} {info} | sv={s.stackVerts.toList.take 8} stack={s.tstack.map showT}"⟩
  for (above, t) in (List.range s.tstack.length).map (fun i => (s.tstack.take i, s.tstack[i]!)) do
    let u := s.stackVerts[t.topDepth]!
    let items := spanItems t
    let E := entryEdges s t
    for j in items do
      if ![NodeType.S, .P, .R, .Q, .V].contains (ty j) || (ty j == .Q && !(ch j).isEmpty) then
        out := bad "kinds" s!"{showT t} j={j} ty={repr (ty j)}" :: out
      if ty j != .V then
        match vs j with
        | (some a, some b) => if a == b then out := bad "two" s!"{showT t} j={j}" :: out
        | _ => out := bad "two" s!"{showT t} j={j} vs={repr (vs j)}" :: out
      if s.items.any (fun it => it.ch.contains j) then
        out := bad "root" s!"{showT t} j={j}" :: out
      if 1 ≤ j && j < 1 + s.g.nv then
        let w := j - 1
        let onPath := (List.range (D + 1)).any fun k => t.topDepth ≤ k && s.stackVerts[k]! == w
        let term := w == t.vStart || onPath || above.any (·.vStart == w)
        let Eab := E ++ above.flatMap (entryEdges s)
        if !(term || interior s Eab w) then
          out := bad "vabove" s!"{showT t} w={w} u={u}" :: out
        if !term && !(touches s E w) then
          out := bad "vtouch" s!"{showT t} w={w} u={u}" :: out
      if ty j != .V then
        let path := (List.range (D + 1)).filterMap fun k => if t.topDepth ≤ k then some s.stackVerts[k]! else none
        let ok (a : Nat) := a == t.vStart || path.contains a || above.any (·.vStart == a) ||
          items.contains (vertItem a)
        match vs j with
        | (some a, some b) =>
          if !(ok a && ok b) then out := bad "ends" s!"{showT t} j={j} vs=({a},{b}) u={u}" :: out
          if !(touches s (edgesBelow s j) a && touches s (edgesBelow s j) b) then
            out := bad "ends_touch" s!"{showT t} j={j}" :: out
        | _ => pure ()
    for w in List.range s.g.nv do
      let inner := touches s E w && interior s E w &&
        items.all fun c => !interior s (edgesBelow s c) w
      if inner && !items.contains (vertItem w) then
        out := bad "vinner" s!"{showT t} w={w}" :: out
  return out

end WalkInvCheck.Content
