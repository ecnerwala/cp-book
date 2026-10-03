import Spqr.Planar
/-!
# Empirical check of the genus inequality (`Spqr/Proofs/PlanarGenus.lean`)

Run with `lake env lean checks/PlanarGenusCheck.lean` from `ports/spqr/lean` (not part of the library).

Random edge lists (loops and parallel edges allowed, isolated vertices allowed) with a random
cyclic order at every vertex give random `IsEmbedding`s; we check
`F + 2 * numNonIsolated ≤ 2 * (2 * numComponents + E)`, and for the deletion of the last edge
(`delLast`, the vertex-orbit shortcut around its four quarter-edges) the bookkeeping the proof
uses: the result is an `IsEmbedding`, `F` changes by exactly `-2`/`+2` when both ends are
non-isolated in the rest according to whether the two sides `4e`/`4e+2` are in distinct/same
face orbits, components/non-isolated counts move as expected.
-/
open Spqr

def lcg (s : Nat) : Nat := (s * 6364136223846793005 + 1442695040888963407) % 2^64
structure Rng where s : Nat
def Rng.next (r : Rng) (k : Nat) : Nat × Rng := let s := lcg r.s; ((s / 2^33) % k, ⟨s⟩)

def shuffle (r : Rng) (l : List Nat) : List Nat × Rng := Id.run do
  let mut a := l.toArray
  let mut r := r
  for i in [0:a.size] do
    let (j, r') := r.next (i + 1)
    r := r'
    a := a.swapIfInBounds i j
  return (a.toList, r)

/-- Random `IsEmbedding` of random edges. -/
def genEmb (r : Rng) (n m : Nat) : List (Nat × Nat) × RotationSystem × Rng := Id.run do
  let mut r := r
  let mut es : List (Nat × Nat) := []
  for _ in [0:m] do
    let (u, r1) := r.next n
    let (v, r2) := r1.next n
    r := r2
    es := es ++ [(u, v)]
  let mut rot : Array (Option Nat) := Array.replicate (4 * m) none
  for x in [0:n] do
    -- dir-0 quarter-edges at x
    let qs := (List.range m).flatMap fun e =>
      (if es[e]!.1 == x then [4 * e] else []) ++ (if es[e]!.2 == x then [4 * e + 2] else [])
    let (cyc, r') := shuffle r qs
    r := r'
    let k := cyc.length
    for i in [0:k] do
      let a := cyc[i]!
      let b := cyc[(i + 1) % k]!
      rot := rot.set! (a ^^^ 1) (some b)
      rot := rot.set! b (some (a ^^^ 1))
  return (es, ⟨rot⟩, r)

def isEmb (es : List (Nat × Nat)) (n : Nat) (rs : RotationSystem) : Bool :=
  rs.size == 4 * es.length && es.all (fun p => p.1 < n && p.2 < n) &&
  decide rs.Total && decide rs.Involution && decide rs.OppositeDir && decide (rs.SameVertex es) &&
  rs.numVertexOrbits == 2 * numNonIsolated es n

def genus (es : List (Nat × Nat)) (n : Nat) (rs : RotationSystem) : Bool :=
  rs.numFaceOrbits + 2 * numNonIsolated es n ≤ 2 * (2 * numComponents es n + es.length)

def sameFace (rs : RotationSystem) (a b : Nat) : Bool := Id.run do
  let mut q := a
  for _ in [0:rs.size + 1] do
    if q == b then return true
    q := (rs.faceStep q).getD q
  return false

/-- Delete the last edge: shortcut the vertex orbits around its four quarter-edges. -/
def delLast (rs : RotationSystem) : RotationSystem := Id.run do
  let m := rs.size / 4 - 1
  let inX (q : Nat) := 4 * m ≤ q
  let exit (q : Nat) : Nat := Id.run do
    let mut q := q
    for _ in [0:3] do
      if inX q then q := (rs.get (q ^^^ 1)).getD q
    return q
  let mut rot : Array (Option Nat) := Array.replicate (4 * m) none
  for a in [0:4 * m] do
    rot := rot.set! a (some (exit ((rs.get a).getD a)))
  return ⟨rot⟩

/-- Delete the first edge: shortcut the vertex orbits around its four quarter-edges and shift. -/
def delFirst (rs : RotationSystem) : RotationSystem := Id.run do
  let m := rs.size / 4 - 1
  let exit (q : Nat) : Nat := Id.run do
    let mut q := q
    for _ in [0:3] do
      if q < 4 then q := (rs.get (q ^^^ 1)).getD q
    return q
  let mut rot : Array (Option Nat) := Array.replicate (4 * m) none
  for a in [0:4 * m] do
    rot := rot.set! a (some (exit ((rs.get (4 + a)).getD (4 + a)) - 4))
  return ⟨rot⟩

/-- `EdgesConn es u v`, by iterated relaxation. -/
def conn (es : List (Nat × Nat)) (n u v : Nat) : Bool := Id.run do
  let mut comp : Array Nat := Array.range n
  for _ in [0:n + 1] do
    for (a, b) in es do
      let c := min comp[a]! comp[b]!
      comp := comp.set! a c
      comp := comp.set! b c
  return comp[u]! == comp[v]!

def main : IO Unit := do
  let mut r : Rng := ⟨12345⟩
  let mut fails := 0
  let mut cnt := 0
  let mut distinct := 0
  let mut same := 0
  let mut planars := 0
  let mut uninserts := 0
  for _ in [0:4000] do
    let (n1, r1) := r.next 7
    let (m1, r2) := r1.next 9
    let n := n1 + 1
    let m := m1
    let (es, rs, r3) := genEmb r2 n m
    r := r3
    if !isEmb es n rs then
      fails := fails + 1
      IO.println s!"not an embedding: {es} {repr rs.rotAdj}"
    else
      cnt := cnt + 1
      if !genus es n rs then
        fails := fails + 1
        IO.println s!"genus: {es} n={n} F={rs.numFaceOrbits} C={numComponents es n}"
      if m > 0 then
        let e := m - 1
        let es' := es.take e
        let rs' := delLast rs
        let (u, v) := es[e]!
        if !isEmb es' n rs' then
          fails := fails + 1
          IO.println s!"delLast not an embedding: {es} {repr rs.rotAdj}"
        else
          let F := rs.numFaceOrbits
          let F' := rs'.numFaceOrbits
          let C := numComponents es n
          let C' := numComponents es' n
          let NI := numNonIsolated es n
          let NI' := numNonIsolated es' n
          let uOld := nonIsolated es' u
          let vOld := nonIsolated es' v
          let sides := sameFace rs (4 * e) (4 * e + 2)
          let planar := rs.numFaceOrbits + 2 * NI == 2 * (2 * C + m)
          if planar then planars := planars + 1
          let mirror := sameFace rs (4 * e + 1) (4 * e + 3)
          if sides != mirror then
            fails := fails + 1
            IO.println s!"mirror: {es} {repr rs.rotAdj}"
          if uOld && vOld then
            if sides then same := same + 1 else distinct := distinct + 1
            if (if sides then F' != F + 2 else F' + 2 != F) then
              fails := fails + 1
              IO.println s!"F change: {es} {repr rs.rotAdj} F={F} F'={F'} sides={sides}"
            if NI' != NI then fails := fails + 1; IO.println "NI both old"
            if !(C' ≤ C + 1 && C ≤ C') then fails := fails + 1; IO.println s!"C both old {C} {C'}"
            if !sides && C' != C then
              fails := fails + 1
              IO.println s!"distinct sides but bridge: {es} {repr rs.rotAdj}"
            if planar && sides && C' != C + 1 then
              fails := fails + 1
              IO.println s!"planar, same sides, not a bridge: {es} {repr rs.rotAdj}"
          else if uOld || vOld then
            if u == v then fails := fails + 1; IO.println "loop one old?"
            if F' != F then fails := fails + 1; IO.println s!"F one old {F} {F'}"
            if NI' + 1 != NI then fails := fails + 1; IO.println "NI one old"
            if C' != C then fails := fails + 1; IO.println s!"C one old {C} {C'}"
          else
            if (if u == v then F' + 4 != F else F' + 2 != F) then
              fails := fails + 1; IO.println s!"F both new {F} {F'}"
            if (if u == v then NI' + 1 != NI else NI' + 2 != NI) then
              fails := fails + 1; IO.println "NI both new"
            if C' + 1 != C then fails := fails + 1; IO.println s!"C both new {C} {C'}"
      -- `IsPlanarEmbedding.uninsert`: first edge, both ends non-isolated in the rest, planar ⇒
      -- sides distinct, the deleted system is planar, pairs `rot 0 ↔ rot 1`, `rot 2 ↔ rot 3`,
      -- and `rot 1 - 4`, `rot 3 - 4` are cofacial.
      if m > 0 then
        let es' := es.drop 1
        let (u, v) := es[0]!
        let rot (q : Nat) : Nat := (rs.get q).getD q
        let A := rot 0; let B := rot 1; let C := rot 2; let D := rot 3
        let C0 := numComponents es n
        let planar := rs.numFaceOrbits + 2 * numNonIsolated es n == 2 * (2 * C0 + m)
        if planar && u != v && conn es' n u v then
          if !(4 ≤ A && 4 ≤ B && 4 ≤ C && 4 ≤ D) then
            fails := fails + 1; IO.println s!"uninsert: neighbour inside edge {es} {repr rs.rotAdj}"
          else
            uninserts := uninserts + 1
            if sameFace rs 0 2 then
              fails := fails + 1; IO.println s!"uninsert: sides cofacial {es} {repr rs.rotAdj}"
            let rs' := delFirst rs
            let C' := numComponents es' n
            if !isEmb es' n rs' then
              fails := fails + 1; IO.println s!"uninsert: not an embedding {es} {repr rs.rotAdj}"
            if rs'.numFaceOrbits + 2 * numNonIsolated es' n != 2 * (2 * C' + (m - 1)) then
              fails := fails + 1; IO.println s!"uninsert: not planar {es} {repr rs.rotAdj}"
            if !(rs'.get (A - 4) == some (B - 4) && rs'.get (B - 4) == some (A - 4) &&
                rs'.get (C - 4) == some (D - 4) && rs'.get (D - 4) == some (C - 4)) then
              fails := fails + 1; IO.println s!"uninsert: pairs {es} {repr rs.rotAdj}"
            if !sameFace rs' (B - 4) (D - 4) then
              fails := fails + 1; IO.println s!"uninsert: B D not cofacial {es} {repr rs.rotAdj}"
            if sameFace rs' (B - 4) (A - 4) then
              fails := fails + 1; IO.println s!"uninsert: B A cofacial {es} {repr rs.rotAdj}"
            for q in [0:rs'.size] do
              if 4 + q != A && 4 + q != B && 4 + q != C && 4 + q != D then
                if rs'.get q != (rs.get (4 + q)).map (· - 4) then
                  fails := fails + 1; IO.println s!"uninsert: frame {es} {repr rs.rotAdj}"
  IO.println s!"embeddings={cnt} planar={planars} distinct={distinct} same={same} uninserts={uninserts} fails={fails}"
#eval main
