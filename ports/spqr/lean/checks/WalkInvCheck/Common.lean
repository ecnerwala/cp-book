import Spqr
import Spqr.StSim
import Spqr.StOpenBlock
import Spqr.StBdPop
import Spqr.Proofs.RInvTree
import Spqr.Proofs.RSkelInv
import Spqr.Ranges
import Spqr.RangesTree
import Spqr.RangesClose
/-!
# `WalkInvCheck`: shared declarations

The executable conjunction of the proposed backbone invariant `WalkInv` (PROOF.md §4.7). The
per-layer modules `WalkInvCheck.Ear` / `.Ranges` / `.R` / `.St` are verbatim copies of the
executable mirrors of `checks/EarCheck.lean` (+ `frontierCheck` of `checks/FrontierCheck.lean`),
`checks/RangesInvCheck.lean`, `checks/RFinishEdgeCheck.lean` and `CheckStSim.lean` (each in its own
namespace, so the helper names may repeat); `WalkInvCheck.Extra` adds the fields no checker had
(`Inv'`, `RSkelInv`, E2–E5, the per-segment `StLive ↔ openBlock` pairing); the root module
`WalkInvCheck` is the single instrumented traversal that evaluates all of them.
-/
namespace WalkInvCheck
open Spqr

instance : Inhabited Spqr.DfsTree := ⟨.node 0 []⟩
instance : Inhabited Spqr.DfsOut := ⟨.back 0 0 .selfLoop⟩

/-- One violation: the field of `WalkInv` (prefixed by its layer), the seed / mode, the site. -/
structure Viol where
  seed : Nat
  tern : Bool
  field : String
  info : String
deriving Repr

end WalkInvCheck
