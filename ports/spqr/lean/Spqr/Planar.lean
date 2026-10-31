import Mathlib.Logic.Relation
import Mathlib.Logic.Function.Iterate
import Spqr.Graph

/-!
# Specification: rotation systems and planar embeddings

An embedding of a graph is given on *quarter-edges*: edge `e` has the four quarter-edges
`4 * e + 2 * side + dir` (`side` picks the endpoint, `dir` the two directions around it).
Three involutions act on them: `QE.flipDir` (`q ^^^ 1`, the other direction at the same
endpoint), `QE.across` (`q ^^^ 3`, the opposite corner of the same edge) and the rotation system's
`rotAdj` (the facing quarter-edge). Walking around a vertex alternates `flipDir` and `rotAdj`;
walking around a face alternates `across` and `rotAdj`. A rotation system is a planar embedding
when Euler's formula holds, counting every vertex/face orbit twice (once per direction).
-/

namespace Spqr

namespace QE

def edge (q : Nat) : Nat := q / 4
def side (q : Nat) : Nat := q / 2 % 2
def dir (q : Nat) : Nat := q % 2
def flipDir (q : Nat) : Nat := q ^^^ 1
def across (q : Nat) : Nat := q ^^^ 3
/-- The quarter-edge of edge `e` at endpoint `side` in direction `dir`. -/
def mk (e : Nat) (side dir : Bool) : Nat := 4 * e + 2 * side.toNat + dir.toNat

/-- The endpoint vertex of `q` in the edge list `es`. -/
def vert (es : List (Nat × Nat)) (q : Nat) : Option Nat :=
  (es[edge q]?).map fun p => if side q = 0 then p.1 else p.2

end QE

/-- A (partial) rotation system: `rotAdj[q]` is the quarter-edge facing `q`, `none` if unset. -/
structure RotationSystem where
  rotAdj : Array (Option Nat)
deriving Repr, Inhabited

/-- The total extension of a partial step map: unset points are fixed. -/
def stepFn (step : Nat → Option Nat) (q : Nat) : Nat := (step q).getD q

/-- `q` is the least point among the first `n` iterates of `step` starting at `q`. -/
def isOrbitMin (step : Nat → Option Nat) (n q : Nat) : Bool :=
  ((List.range n).foldl (fun (p : Nat × Bool) _ => (stepFn step p.1, p.2 && decide (q ≤ p.1))) (q, true)).2

/-- Number of orbits of `step` on `[0, n)`: the number of points that are minimal on their
`n`-step forward orbit. For a permutation of `[0, n)` this is the number of cycles. -/
def numOrbits (step : Nat → Option Nat) (n : Nat) : Nat :=
  ((List.range n).filter (isOrbitMin step n)).length

/-- Vertex `v` has an incident edge in `es`. -/
def nonIsolated (es : List (Nat × Nat)) (v : Nat) : Bool := es.any fun p => p.1 == v || p.2 == v

def numNonIsolated (es : List (Nat × Nat)) (nVerts : Nat) : Nat :=
  ((List.range nVerts).filter (nonIsolated es)).length

/-- One relaxation round: both endpoints of every edge receive the smaller of their labels
(a vertex without a label entry is its own label). -/
def relaxLabels (es : List (Nat × Nat)) (l : List Nat) : List Nat :=
  es.foldl (fun l p => let m := min (l[p.1]?.getD p.1) (l[p.2]?.getD p.2); (l.set p.1 m).set p.2 m) l

/-- Component labels: after `nVerts` rounds of relaxation every vertex carries the least vertex
of its connected component. -/
def compLabels (es : List (Nat × Nat)) (nVerts : Nat) : List Nat :=
  (relaxLabels es)^[nVerts] (List.range nVerts)

/-- Number of connected components of `es` with at least one edge: the non-isolated vertices
that are the least vertex of their component. -/
def numComponents (es : List (Nat × Nat)) (nVerts : Nat) : Nat :=
  let l := compLabels es nVerts
  ((List.range nVerts).filter fun v => nonIsolated es v && l[v]?.getD v == v).length

namespace RotationSystem

variable (rs : RotationSystem)

def size : Nat := rs.rotAdj.size
def get (q : Nat) : Option Nat := (rs.rotAdj[q]?).bind id

/-- One step around a vertex. -/
def vertexStep (q : Nat) : Option Nat := rs.get (QE.flipDir q)
/-- One step around a face. -/
def faceStep (q : Nat) : Option Nat := rs.get (QE.across q)

def Total : Prop := ∀ q, q < rs.size → (rs.get q).isSome
def Involution : Prop := ∀ q, q < rs.size → ∀ r ∈ rs.get q, r < rs.size ∧ rs.get r = some q
def OppositeDir : Prop := ∀ q, q < rs.size → ∀ r ∈ rs.get q, QE.dir r ≠ QE.dir q
def SameVertex (es : List (Nat × Nat)) : Prop :=
  ∀ q, q < rs.size → ∀ r ∈ rs.get q, QE.vert es q = QE.vert es r

/-- `q` and `r` lie on the same vertex orbit. -/
def SameVertexOrbit (q r : Nat) : Prop :=
  Relation.ReflTransGen (fun a b => rs.vertexStep a = some b) q r
/-- `q` and `r` lie on the same face orbit. -/
def SameFaceOrbit (q r : Nat) : Prop :=
  Relation.ReflTransGen (fun a b => rs.faceStep a = some b) q r

def numVertexOrbits : Nat := numOrbits rs.vertexStep rs.size
def numFaceOrbits : Nat := numOrbits rs.faceStep rs.size

instance : Decidable rs.Total := by unfold Total; infer_instance
instance : Decidable rs.Involution := by unfold Involution; infer_instance
instance : Decidable rs.OppositeDir := by unfold OppositeDir; infer_instance
instance (es) : Decidable (rs.SameVertex es) := by unfold SameVertex; infer_instance

end RotationSystem

/-- `rs` is a rotation system of the graph with edges `es` on vertices `[0, nVerts)`: a total
involution on the `4 * |es|` quarter-edges pairing quarter-edges at the same vertex with opposite
direction, whose vertex orbits are exactly two per non-isolated vertex (one per direction), i.e.
a single cyclic order of the edges around each vertex. -/
structure IsEmbedding (es : List (Nat × Nat)) (nVerts : Nat) (rs : RotationSystem) : Prop where
  size : rs.size = 4 * es.length
  verts : ∀ p ∈ es, p.1 < nVerts ∧ p.2 < nVerts
  total : rs.Total
  involution : rs.Involution
  opposite_dir : rs.OppositeDir
  same_vertex : rs.SameVertex es
  vertex_orbits : rs.numVertexOrbits = 2 * numNonIsolated es nVerts

/-- Euler's formula, each orbit counted twice: `F + 2V = 2(2C + E)` with `C` the number of
components with edges (`C = 1` for a connected edge set). -/
def EulerFormula (es : List (Nat × Nat)) (nVerts : Nat) (rs : RotationSystem) : Prop :=
  rs.numFaceOrbits + 2 * numNonIsolated es nVerts = 2 * (2 * numComponents es nVerts + es.length)

structure IsPlanarEmbedding (es : List (Nat × Nat)) (nVerts : Nat) (rs : RotationSystem) : Prop
    extends IsEmbedding es nVerts rs where
  euler : EulerFormula es nVerts rs

instance (es nVerts rs) : Decidable (EulerFormula es nVerts rs) := by unfold EulerFormula; infer_instance

/-- The graph `es` on `[0, nVerts)` is planar. -/
def Planar (es : List (Nat × Nat)) (nVerts : Nat) : Prop := ∃ rs, IsPlanarEmbedding es nVerts rs

instance (es nVerts rs) : Decidable (IsEmbedding es nVerts rs) :=
  decidable_of_iff (rs.size = 4 * es.length ∧ (∀ p ∈ es, p.1 < nVerts ∧ p.2 < nVerts) ∧ rs.Total ∧
      rs.Involution ∧ rs.OppositeDir ∧ rs.SameVertex es ∧ rs.numVertexOrbits = 2 * numNonIsolated es nVerts)
    ⟨fun ⟨a, b, c, d, e, f, g⟩ => ⟨a, b, c, d, e, f, g⟩, fun ⟨a, b, c, d, e, f, g⟩ => ⟨a, b, c, d, e, f, g⟩⟩

instance (es nVerts rs) : Decidable (IsPlanarEmbedding es nVerts rs) :=
  decidable_of_iff (IsEmbedding es nVerts rs ∧ EulerFormula es nVerts rs)
    ⟨fun ⟨a, b⟩ => ⟨a, b⟩, fun ⟨a, b⟩ => ⟨a, b⟩⟩

end Spqr
