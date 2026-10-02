import Spqr

open Spqr

/-- Per-phase wall-clock timing harness: `prof_lean < input`. -/
def main : IO Unit := do
  let t0 ← IO.monoMsNow
  let input ← (← IO.getStdin).readToEnd
  let toks := (input.splitOn " ").flatMap (·.splitOn "\n") |>.filter (· ≠ "") |>.map String.toNat!
  let toks := toks.toArray
  let nv := toks[0]!; let ne := toks[1]!; let tern := toks[2]!
  let mut edges : Array (Nat × Nat) := #[]
  for i in [0:ne] do
    edges := edges.push (toks[3 + 2 * i]!, toks[4 + 2 * i]!)
  let g : Graph := ⟨nv, edges⟩
  let t1 ← IO.monoMsNow
  IO.println s!"parse: {t1 - t0} ms (nv={g.nv} ne={g.ne})"
  let forest := g.dfsForest [] []
  let t2 ← IO.monoMsNow
  IO.println s!"dfs: {t2 - t1} ms (roots={forest.length})"
  let w := g.walk (tern != 0) forest
  let t3 ← IO.monoMsNow
  IO.println s!"walk: {t3 - t2} ms (items={w.items.size} ticks={w.ticks})"
  let t := relabelTree g w.items
  let t4 ← IO.monoMsNow
  IO.println s!"relabel: {t4 - t3} ms (nodes={t.types.size} ticks={(relabelRun g w.items).ticks})"
  let out := String.join ((t.adjDat.toList.map fun x => s!"{x.ne},{x.destNv}").map (" " ++ ·))
  let t5 ← IO.monoMsNow
  IO.println s!"format adjDat: {t5 - t4} ms ({out.length} chars)"
