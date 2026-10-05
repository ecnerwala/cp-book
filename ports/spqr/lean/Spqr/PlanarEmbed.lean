import Spqr.PlanarRelabel

/-!
# Gluing the local embeddings into an embedding of the original graph

Bottom-up over the preorder, every node reports the two exposed quarter-edges of its cap on
each side (`outerE`), and S/P/R nodes link their children's exposed ends according to their local
rotation system, through the twin of every non-cap node-edge. Vertex items splice the loops /
blocks hanging off a vertex into the vertex's rotation.
-/

namespace Spqr

namespace PlanarSpqrTree

variable (t : PlanarSpqrTree)

def size : Nat := t.types.size
def children (i : Nat) : List Nat :=
  (List.range (t.chBounds[i + 1]! - t.chBounds[i]!)).map fun k => t.chDat[t.chBounds[i]! + k]!

structure EmbedState where
  rotAdj : Array (Option Nat)
  /-- Per item, the exposed quarter-edges `[side 0 outer, side 0 inner, side 1 outer, side 1 inner]`. -/
  outerE : Array (Array (Option Nat))

abbrev EmbedM := StateM EmbedState

namespace EmbedM

def link (a b : Option Nat) : EmbedM Unit := do
  match a, b with
  | some a, some b => modify fun s => { s with rotAdj := (s.rotAdj.set! a (some b)).set! b (some a) }
  | _, _ => pure ()

def outer (i k : Nat) : EmbedM (Option Nat) := do return (← get).outerE[i]![k]!
def setOuter (i k : Nat) (q : Option Nat) : EmbedM Unit :=
  modify fun s => { s with outerE := s.outerE.modify i (·.set! k q) }

end EmbedM

open EmbedM in
def embedItem (i : Nat) : EmbedM Unit := do
  match t.types[i]! with
  | .F =>
    for j in t.children i do
      let a ← outer j 0
      if a.isSome then link a (← outer j 1)
  | .V =>
    let mut qes : Option Nat × Option Nat := (none, none)
    for j in t.children i do
      if qes.1.isNone then
        qes := (← outer j 0, ← outer j 1)
      else
        link qes.2 (← outer j 0)
        qes := (qes.1, ← outer j 1)
    setOuter i 0 qes.1
    setOuter i 1 qes.2
  | .Q =>
    let e := (t.origId[i]!).getD 0
    let flip := t.edgeFlipped[e]!
    let mut qes : Array (Option Nat) :=
      #[some (QE.mk e flip false), some (QE.mk e flip true), some (QE.mk e (!flip) false), some (QE.mk e (!flip) true)]
    match t.children i with
    | [] => for k in [0 : 4] do setOuter i k qes[k]!
    | j :: rest =>
      if t.types[j]! == .O then
        link qes[1]! qes[2]!
        setOuter i 0 qes[0]!
        setOuter i 1 qes[3]!
      else
        if t.types[j]! != .I then
          link qes[1]! (← outer j 0)
          qes := qes.set! 1 (← outer j 1)
          link qes[2]! (← outer j 3)
          qes := qes.set! 2 (← outer j 2)
        let k := rest.head!
        if (← outer k 0).isSome then
          link qes[3]! (← outer k 0)
          link qes[2]! (← outer k 1)
        else
          link qes[3]! qes[2]!
        setOuter i 0 qes[0]!
        setOuter i 1 qes[1]!
  | .O | .I => pure ()
  | _ =>
    let neSt := t.neBounds[i]!
    let neEn := t.neBounds[i + 1]!
    let treeQeToQe (q : Nat) : EmbedM (Option Nat) := do
      let twin := (t.nodeEdges[QE.edge q]!.twin).getD 0
      outer (t.nodeEdges[twin]!.node) (2 * QE.side q + QE.dir q)
    for ta in [4 * neSt : 4 * neEn] do
      if let some tb := t.neRotAdj[ta]! then
        if tb < ta then continue
        let qb ← treeQeToQe tb
        if ta < 4 * (neSt + 1) then
          setOuter i (2 * QE.side ta + (1 - QE.dir ta)) qb
        else
          let qa ← treeQeToQe ta
          if ta % 4 == 2 && tb % 4 == 1 then
            let v := t.nodeVerts[(t.nodeEdges[QE.edge ta]!).nvs.2]!.vert
            if (← outer v 0).isSome then
              link qa (← outer v 1)
              link qb (← outer v 0)
            else link qa qb
          else link qa qb

/-- Glue the node embeddings into a rotation system of the original graph, or `none` if some
node is nonplanar. -/
def planarEmbed : Option RotationSystem :=
  if t.nodePlanar.all id then
    let s := ((List.range t.size).reverse.forM t.embedItem).run
      ⟨Array.replicate (4 * t.ne) none, Array.replicate t.size (Array.replicate 4 none)⟩ |>.2
    some ⟨s.rotAdj⟩
  else none

end PlanarSpqrTree

end Spqr
