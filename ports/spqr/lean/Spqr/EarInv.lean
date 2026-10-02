import Spqr.WalkSpec

/-!
# Ear invariant (PROOF.md §4.1–4.3)

The corrected attachment invariant `WalkState.Inv'` (`TEntry.Term'`, `WalkState.Stack`; PROOF.md
§4.2b correction) and its primitive steps live in `Spqr/WalkSpec.lean`, where `Step` is stated
through it. This file is reserved for the ear invariant (span ownership + graph facts per entry).
-/
