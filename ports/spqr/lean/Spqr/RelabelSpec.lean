import Spqr.Spec
import Spqr.ItemSpec

/-!
# Per-node characterization of `relabel` (the interface of phase 3)

`relabelTree g items` numbers the items in preorder (`idx : ItemId → Nat`) and, for the item `i`
numbered `idx i`, writes exactly what `RelabelNode g items t idx i` says: its type / `origId`,
its node-vert list `Items.nvList`, its children in `Items.ordered` order, and the skeleton of
`layoutNode` for those children (edges, the two adjacency rows of every node-vert) with the
`twin` back-patches. Everything a node writes into its *children's* slots (`par`, `vertParNv`, the
cap `twin`) is stated in the parent's record; the root's slots are in `RelabelIdx`.

The vertex positions used by `orderedChildren` and `edgeChildren` are the scratch array `vertPos`;
for an R node they are the positions in its own node-vert list (`PosOK`), for other nodes they are
whatever an earlier R node left behind, so the record only exposes them existentially
(`RelabelNode.layout`); `layoutNode` ignores them for non-R types anyway.

Not stated here: `adjBounds[2 nvSt] = 2 neSt` (the CSR start of a node's rows). It is the previous
node's last bound, and for an R node the last bound of `layoutNode .R` is `2 neEn` only when every
edge child `(a, c)` is oriented `a < c` along the node-vert order (the counting phase writes to
bounds `2 a + 2` and `2 c + 1`), which `Items.WF` does not promise. See `relabel_adj_spec`.
-/

namespace Spqr

namespace Items

variable (g : Graph) (items : Items)

/-- The vertices of item `i`'s node in node-vert order: first cap endpoint, the V children in
`ch` order, second cap endpoint (`nodeVerts` of `relabel`, before re-indexing through
`vertIndex`). -/
def nvList (i : ItemId) : List Nat :=
  (items.vs i).1.toList ++ ((items.ch i).filter (· < 1 + g.nv)).map (· - 1) ++ (items.vs i).2.toList

/-- Sort key of `RelabelM.orderedChildren` with vertex positions `pos`. -/
def loc (nvSt : Nat) (pos : Nat → Nat) (c : ItemId) : Nat :=
  if c < 1 + g.nv then 2 * (pos (c - 1) - nvSt)
  else (pos ((items.vs c).1.getD 0) - nvSt) + (pos ((items.vs c).2.getD 0) - nvSt)

/-- The children of `i` in output order: `ch`, stably sorted by `loc` for R nodes. -/
def ordered (i : ItemId) (nvSt : Nat) (pos : Nat → Nat) : List ItemId :=
  if items.type i ≠ .R then items.ch i
  else (items.ch i).mergeSort fun a b => items.loc g nvSt pos a ≤ items.loc g nvSt pos b

/-- The non-V items of `children` as pairs of vertex positions (`edgeChildren` of `relabel`). -/
def edgeChildren (pos : Nat → Nat) (children : List ItemId) : List (Nat × Nat) :=
  (children.filter (· ≥ 1 + g.nv)).map fun c => (pos ((items.vs c).1.getD 0), pos ((items.vs c).2.getD 0))

/-- `pos` sends every vertex of `vl` (listed from node-vert `nvSt`) to a node-vert holding it. -/
def PosOK (nvSt : Nat) (vl : List Nat) (pos : Nat → Nat) : Prop :=
  ∀ v ∈ vl, nvSt ≤ pos v ∧ vl[pos v - nvSt]? = some v

/-- The node has a cap edge: it is a node and not a block-root Q. -/
def hasCap (i : ItemId) : Bool :=
  (items.type i).isNode && !(items.type i == .Q && !(items.ch i).isEmpty)

def capCount (i : ItemId) : Nat := if items.hasCap i then 1 else 0

/-- Number of node-edges: one per non-V child, plus the cap. -/
def nEdges (i : ItemId) : Nat :=
  (if (items.type i).isNode then (items.ch i).countP (· ≥ 1 + g.nv) else 0) + items.capCount i

/-- `origId` of item `i`. -/
def origOf (i : ItemId) : Option Nat :=
  match items.type i with
  | .V => some (i - 1)
  | .Q => some (i - 1 - g.nv)
  | _ => none

/-- Every edge child of an R node is oriented along its node-vert order (and the node-vert list
has no repeats). `Items.StNumbered` implies this; `Items.WF` does not. -/
def ROriented : Prop :=
  ∀ i, i < items.size → items.type i = .R →
    (items.nvList g i).Nodup ∧
    ∀ p ∈ items.virtualEdges i, (items.nvList g i).idxOf p.1 < (items.nvList g i).idxOf p.2

end Items

/-- The skeleton `relabel` lays out for item `i` (numbered `idx i`) with vertex positions `pos`. -/
def nodeLayout (g : Graph) (items : Items) (t : SpqrTree) (idx : ItemId → Nat) (i : ItemId)
    (pos : Nat → Nat) : Layout :=
  layoutNode (items.type i) (idx i) (t.nvRange (idx i)).1 (t.nvRange (idx i)).2
    (t.neRange (idx i)).1 (t.neRange (idx i)).2
    (items.edgeChildren g pos (items.ordered g i (t.nvRange (idx i)).1 pos))

/-- The part of item `i`'s record that depends on the vertex positions `pos` seen by `relabel`
when it processed `i`: children order, skeleton, back-patched twins, preorder of the children. -/
structure RelabelLayout (g : Graph) (items : Items) (t : SpqrTree) (idx : ItemId → Nat) (i : ItemId)
    (pos : Nat → Nat) : Prop where
  pos_ok : items.type i = .R → Items.PosOK (t.nvRange (idx i)).1 (items.nvList g i) pos
  children : t.children (idx i) = (items.ordered g i (t.nvRange (idx i)).1 pos).map idx
  /-- Node-edge `k` of the node is edge `k` of `layoutNode` (up to its `twin`). -/
  edge_node : ∀ k, k < items.nEdges g i →
    (t.nodeEdges[(t.neRange (idx i)).1 + k]!).node = ((nodeLayout g items t idx i pos).edges[k]!).node
  edge_nvs : ∀ k, k < items.nEdges g i →
    (t.nodeEdges[(t.neRange (idx i)).1 + k]!).nvs = ((nodeLayout g items t idx i pos).edges[k]!).nvs
  /-- Global adjacency bound `2 nvSt + j`, `j ∈ [1, 2 nVerts]`, is local bound `j`. -/
  adj_bounds : ∀ j, 1 ≤ j → j ≤ 2 * (items.nvList g i).length →
    t.adjBounds[2 * (t.nvRange (idx i)).1 + j]! = (nodeLayout g items t idx i pos).adjBounds[j]!
  /-- Global adjacency slot `2 neSt + j`, `j < 2 nEdges`, is local slot `j`. -/
  adj_dat : ∀ j, j < 2 * items.nEdges g i →
    t.adjDat[2 * (t.neRange (idx i)).1 + j]! = (nodeLayout g items t idx i pos).adjDat[j]!
  /-- The `k`-th V child (in output order) is told its node-vert `nvSt + |vs.1| + k`. -/
  vert_par_nv : ∀ k (hk : k < ((items.ordered g i (t.nvRange (idx i)).1 pos).filter (· < 1 + g.nv)).length),
    t.vertParNv[idx ((items.ordered g i (t.nvRange (idx i)).1 pos).filter (· < 1 + g.nv))[k]]! =
      some ((t.nvRange (idx i)).1 + (items.vs i).1.toList.length + k)
  /-- Node-edge `capCount + k` of a node is twinned with the cap (first node-edge) of its `k`-th
  non-V child, both ways. -/
  twin : (items.type i).isNode →
    ∀ k (hk : k < ((items.ordered g i (t.nvRange (idx i)).1 pos).filter (· ≥ 1 + g.nv)).length),
      t.twin ((t.neRange (idx i)).1 + items.capCount i + k) =
        some (t.neRange (idx ((items.ordered g i (t.nvRange (idx i)).1 pos).filter (· ≥ 1 + g.nv))[k])).1 ∧
      t.twin (t.neRange (idx ((items.ordered g i (t.nvRange (idx i)).1 pos).filter (· ≥ 1 + g.nv))[k])).1 =
        some ((t.neRange (idx i)).1 + items.capCount i + k)
  /-- Preorder: the first child is numbered next, each later child right after its predecessor's
  subtree. -/
  child_idx : ∀ k (hk : k < (items.ordered g i (t.nvRange (idx i)).1 pos).length),
    idx (items.ordered g i (t.nvRange (idx i)).1 pos)[k] =
      if k = 0 then idx i + 1 else t.subtreeEnd[idx (items.ordered g i (t.nvRange (idx i)).1 pos)[k - 1]]!
  subtree_end : t.subtreeEnd[idx i]! =
    match (items.ordered g i (t.nvRange (idx i)).1 pos).getLast? with
    | none => idx i + 1
    | some c => t.subtreeEnd[idx c]!

/-- What `relabel` writes for item `i`, numbered `idx i`, into `t = relabelTree g items`: its own
slots, its ranges, its node-verts (re-indexed through `vertIndex`, as `relabelTree` does), and the
slots it writes for its children (`par`, `vertParNv`, the cap twin of children of non-nodes). -/
structure RelabelNode (g : Graph) (items : Items) (t : SpqrTree) (idx : ItemId → Nat) (i : ItemId) :
    Prop where
  idx_lt : idx i < t.size
  type : t.type (idx i) = items.type i
  orig : t.origId[idx i]! = items.origOf g i
  vert_index : items.type i = .V → t.vertIndex[i - 1]! = some (idx i)
  edge_index : items.type i = .Q → t.edgeIndex[i - 1 - g.nv]! = some (idx i) ∧
    t.edgeFlipped[i - 1 - g.nv]! = ((items.vs i).1 != some (g.edges[i - 1 - g.nv]!).1)
  ch_range : (t.chRange (idx i)).2 = (t.chRange (idx i)).1 + (items.ch i).length
  nv_range : (t.nvRange (idx i)).2 = (t.nvRange (idx i)).1 + (items.nvList g i).length
  ne_range : (t.neRange (idx i)).2 = (t.neRange (idx i)).1 + items.nEdges g i
  /-- Node-vert `nvSt + k` is vertex `nvList[k]`, named by the index of its V item. -/
  node_verts : ∀ k (hk : k < (items.nvList g i).length),
    t.nodeVerts[(t.nvRange (idx i)).1 + k]! = ⟨idx i, (t.vertIndex[(items.nvList g i)[k]]!).getD 0⟩
  child_lt : ∀ c ∈ items.ch i, idx c < t.size
  child_par : ∀ c ∈ items.ch i, t.parent (idx c) = some (idx i)
  child_par_nv_none : ∀ c ∈ items.ch i, 1 + g.nv ≤ c → t.vertParNv[idx c]! = none
  /-- Children of F / V items get no cap twin. -/
  child_cap_twin_none : ¬ (items.type i).isNode → ∀ c ∈ items.ch i, items.hasCap c →
    t.twin (t.neRange (idx c)).1 = none
  layout : ∃ pos, RelabelLayout g items t idx i pos

/-- The global facts: `idx` is a bijection `[0, items.size) → [0, t.size)` with the root at `0`,
`vertIndex` / `edgeIndex` are `idx` on V / Q items, array sizes agree and the CSR arrays start at
`0` and end at the data sizes. -/
structure RelabelIdx (g : Graph) (items : Items) (t : SpqrTree) (idx : ItemId → Nat) : Prop where
  size : t.size = items.size
  root : idx rootItem = 0
  root_par : t.parent 0 = none
  root_par_nv : t.vertParNv[0]! = none
  lt : ∀ i, i < items.size → idx i < items.size
  inj : ∀ i j, i < items.size → j < items.size → idx i = idx j → i = j
  vert_index : ∀ v, v < g.nv → t.vertIndex[v]! = some (idx (vertItem v))
  edge_index : ∀ e, e < g.ne → t.edgeIndex[e]! = some (idx (edgeItem g e))
  nv : t.nv = g.nv
  ne : t.ne = g.ne
  sizes : t.Sizes
  ch_zero : t.chBounds[0]! = 0
  nv_zero : t.nvBounds[0]! = 0
  ne_zero : t.neBounds[0]! = 0
  adj_zero : t.adjBounds[0]! = 0
  ch_last : t.chBounds[t.size]! = t.chDat.size
  nv_last : t.nvBounds[t.size]! = t.nodeVerts.size
  ne_last : t.neBounds[t.size]! = t.nodeEdges.size
  adj_dat_size : t.adjDat.size = 2 * t.nodeEdges.size

/-- Phase 3, per node: `relabelTree` numbers a well-formed item tree in preorder and writes each
node's record. -/
theorem relabel_node_spec (g : Graph) (items : Items) (h : items.WF g) :
    ∃ idx, idx rootItem = 0 ∧ RelabelIdx g items (relabelTree g items) idx ∧
      ∀ i, i < items.size → RelabelNode g items (relabelTree g items) idx i := by
  sorry

/-- The CSR facts about `adjBounds` that need the R children oriented: each node's rows start at
`2 neSt`, and the last bound is the end of `adjDat`. -/
theorem relabel_adj_spec (g : Graph) (items : Items) (h : items.WF g) (hor : items.ROriented g) :
    (∀ n, n < (relabelTree g items).size →
      (relabelTree g items).adjBounds[2 * ((relabelTree g items).nvRange n).1]! =
        2 * ((relabelTree g items).neRange n).1) ∧
    (relabelTree g items).adjBounds[2 * (relabelTree g items).nodeVerts.size]! =
      (relabelTree g items).adjDat.size := by
  sorry

end Spqr
