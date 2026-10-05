import Spqr.SepPair

/-!
# Blocks relative to a lowpoint-sorted DFS forest

Blocks of `g` are the classes of `Graph.SameBlock`: no single vertex separates two edges of a
block. On the DFS forest (`DfsData`) they are read off the `lowval ≥ d` classes (PROOF.md §2): the
*block roots* are the roots and the child ends of `component` / `bridge` tree edges, and every
non-bridge, non-loop edge belongs to the block of the deepest block root above its deeper
endpoint. `Spqr.Proofs.Blocks` proves `SameBlock` equal to this description.
-/

namespace Spqr

namespace Graph

variable (g : Graph)

/-- No vertex separates `e` from `e'`: the two edges lie in a common block. -/
def SameBlock (e e' : Nat) : Prop := ∀ v, g.EdgeConn (· ≠ v) e e'

end Graph

namespace DfsOut

/-- The deeper endpoint of the out-edge `o` of `v`: the child for a tree edge, `v` itself for a
back edge or loop. -/
def deep (v : Nat) : DfsOut → Nat
  | .back .. => v
  | .tree _ _ child => child.v

/-- The shallower endpoint: `v` for a tree edge, the destination for a back edge or loop. -/
def shallow (v : Nat) : DfsOut → Nat
  | .back _ dest _ => dest
  | .tree .. => v

/-- A bridge or a loop is a block by itself. -/
def isSingleton (o : DfsOut) : Prop := o.cls = .bridge ∨ o.cls = .selfLoop

end DfsOut

namespace DfsData

variable (d : DfsData)

def IsRoot (v : Nat) : Prop := ∀ p, ¬d.IsParent p v

/-- `c` tops a block: it is a root, or is entered by a `component` or `bridge` tree edge. -/
def BlockRoot (c : Nat) : Prop :=
  d.IsRoot c ∨ ∃ p, ∃ o ∈ d.outs p, o.isTree = true ∧ o.dest = c ∧
    (o.cls = .component ∨ o.cls = .bridge)

/-- `c` is the deepest block root above (or at) `x`. -/
def BlockTop (c x : Nat) : Prop :=
  d.BlockRoot c ∧ d.Anc c x ∧ ∀ c', d.BlockRoot c' → d.Anc c' x → d.Anc c' c

/-- Edge `e` has its deeper endpoint in `T_c`. -/
def DeepIn (c e : Nat) : Prop := ∃ v, ∃ o ∈ d.outs v, o.e = e ∧ d.Anc c (o.deep v)

/-- Edge `e` lies in the block topped by `c`: it is neither a bridge nor a loop, and `c` is the
block top of its deeper endpoint. -/
def InBlock (c e : Nat) : Prop :=
  ∃ v, ∃ o ∈ d.outs v, o.e = e ∧ ¬o.isSingleton ∧ d.BlockTop c (o.deep v)

end DfsData

end Spqr
