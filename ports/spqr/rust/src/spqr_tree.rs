//! SPQR tree. Port of `wala::spqr_tree` / `wala::planar_spqr_tree` (cp-book `graph/spqr_tree.hpp`).
//!
//! All ids are `i32` with `-1` as the "none" sentinel, exactly like the C++ original.
//! The C++ `[&]` lambdas are expressed as methods on small state structs (`LowvalDfs`, `Builder`, `Relabel`);
//! `with_planarity` is a const generic, so the non-planar build compiles the planarity code out.
// The arithmetic mirrors the C++ literally (e.g. `2 * x + 0` for quarter-edge ids), so silence the lints about it.
#![allow(clippy::identity_op, clippy::erasing_op, clippy::int_plus_one, clippy::needless_range_loop)]

use std::fmt;
use std::ops::{Deref, DerefMut, Index, IndexMut, Range};

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct CsrIndex {
	pub bounds: Vec<i32>,
}

impl CsrIndex {
	pub fn indices(&self, i: i32) -> Range<i32> {
		self.bounds[i as usize]..self.bounds[i as usize + 1]
	}
	pub fn slice<'a, T>(&self, i: i32, base: &'a [T]) -> &'a [T] {
		&base[self.bounds[i as usize] as usize..self.bounds[i as usize + 1] as usize]
	}
	pub fn num_rows(&self) -> i32 {
		if self.bounds.is_empty() { 0 } else { self.bounds.len() as i32 - 1 }
	}
	pub fn num_entries(&self) -> i32 {
		self.bounds.last().copied().unwrap_or(0)
	}
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Csr<T> {
	pub bounds: Vec<i32>,
	pub dat: Vec<T>,
}

impl<T> Csr<T> {
	pub fn indices(&self, i: i32) -> Range<i32> {
		self.bounds[i as usize]..self.bounds[i as usize + 1]
	}
	pub fn size(&self) -> i32 {
		self.bounds.len() as i32 - 1
	}
}

impl<T> Index<i32> for Csr<T> {
	type Output = [T];
	fn index(&self, i: i32) -> &[T] {
		&self.dat[self.bounds[i as usize] as usize..self.bounds[i as usize + 1] as usize]
	}
}

impl<T> IndexMut<i32> for Csr<T> {
	fn index_mut(&mut self, i: i32) -> &mut [T] {
		&mut self.dat[self.bounds[i as usize] as usize..self.bounds[i as usize + 1] as usize]
	}
}

#[derive(Clone, Debug, Default)]
pub struct CsrBuilder<T> {
	pub bounds: Vec<i32>,
	pub dat: Vec<T>,
}

impl<T: Copy + Default> CsrBuilder<T> {
	pub fn new(n: i32) -> Self {
		Self { bounds: vec![0; n as usize + 1], dat: Vec::new() }
	}

	pub fn count(&mut self, k: i32) {
		self.bounds[k as usize + 1] += 1;
	}
	pub fn allocate(&mut self) {
		let mut l = 0;
		for i in 1..self.bounds.len() {
			let c = self.bounds[i];
			self.bounds[i] = l;
			l += c;
		}
		self.dat.resize(l as usize, T::default());
	}
	#[must_use]
	pub fn push(&mut self, k: i32) -> &mut T {
		let idx = self.bounds[k as usize + 1];
		self.bounds[k as usize + 1] += 1;
		&mut self.dat[idx as usize]
	}
	#[must_use]
	pub fn finalize(self) -> Csr<T> {
		Csr { bounds: self.bounds, dat: self.dat }
	}
}

/// The SPQR tree of a graph is a canonical/"maximal" decomposition of the graph by 2-vertex cuts.
/// The tree consists of nodes which are graphs of virtual edges (vedges), corresponding to nontrivial 2-vertex cuts.
/// Virtual edges are paired, and we can reassemble the graph by gluing nodes at their matching vedges (and removing the vedge).
/// Real edges are represented as special Q nodes which each contain exactly 1 real edge and exactly 1 vedge.
///
/// Traditionally, the SPQR tree is defined for each biconnected component,
/// but we will embed the SPQR decompositions inside the block-cut tree to get a (rooted) decomposition of the entire graph.
///
/// As such, we will have a tree of "items", which consist of SPQR nodes, real vertices, and a special "forest root" item:
///  - Each vertex will be a child of the topmost node which contains it (or the forest root).
///  - Each block will be a subtree of nodes rooted at a Q edge, which is the child of one of its vertices.
///
/// Item types:
///  F - forest - a root node corresponding to the whole forest.
///  V - vertex - not really a node, just there because they're mixed into the tree a la block/cut tree
///  Q - real edge - has exactly 1 vedge and 1 real edge
///  I - bridge - has exactly 1 vedge connecting to a bridge Q node
///  O - self-loop - has exactly 1 vedge connecting to a self-loop Q node
///  S - series - a cycle of >= 3 vedges; note that any 2 vertices of the cycle form a cut
///  P - parallel - a parallel group of >= 3 vedges with the same endpoints
///  R - rigid - a 3-vertex-connected component
///
/// Q nodes occur in 2 places: block roots and block leaves.
/// Block leaf Q's simply have no children.
/// Block root Q's have 2 children: their vedge, and their deeper vertex (unless it's a self-loop).
///
/// Degenerate blocks:
///  - a block consisting of a self-loop is a Q node connected to an O node.
///  - a block consisting of a bridge is a Q node connected to an I node.
///  - a block consisting of exactly 2 parallel edges is represented by 2 glued Q nodes.
///
/// We have several id spaces:
///  - items are in preorder
///  - node_verts (nv's) are each node's vertices, given in node order then s-t order.
///  - node_edges (ne's) are each node's vedges, given in node order then a s-t order.
///  - node_adj is each node_vert's incident vedges, given as 2 lists per nv: left/rightwards edges each in reverse s-t order.
///  - original verts and original edges can be converted to items as vert_item / edge_item
///
/// Children of a node will be sorted in s-t order.
/// Specifically vertices are sorted, and edges are guaranteed to satisfy the strong "dominance" partial order:
/// if a.nvs[0] <= b.nvs[0] and a.nvs[1] <= b.nvs[1], then a <= b. (In practice, we'll sort by midpoint.)
/// Adjacency lists are sorted as "center-is-longest", which helps make laminar/bracket cases clean.
///   (5->4) (5->3) (5->2) (5->1) *vertex 5* (5->9) (5->8) (5->7) (5->6)
/// More specifically, node_adj contains two lists per vertex: 2*nv+0 is leftwards and 2*nv+1 is rightwards.
///
/// All id's are item indices unless clearly nv/ne id's.
///
/// In general, there are 2 ways to use the SPQR tree: the rooted view and the unrooted view.
///  - The rooted view uses par / ch walks, and either treats the tree as 1 top-down big decomposition, or walks in paths up/down the tree with LCA-like queries.
///  - The unrooted view mostly uses nv/ne/nd lists and works locally within a node/sometimes jumps between them.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SpqrTree {
	pub vert_index: Vec<i32>,
	pub edge_index: Vec<i32>,
	/// Whether edge e is stored in its node with its endpoints swapped
	/// (`item_vs[edge_item(e)][0] != edges[e][0]`).
	pub edge_flipped: Vec<bool>,

	pub par: Vec<i32>,
	pub subtree_end: Vec<i32>,
	pub types: Vec<NodeType>,
	pub orig_id: Vec<i32>,

	pub ch: Csr<i32>,
	pub node_verts: Vec<NodeVert>,
	pub node_nvs: CsrIndex,
	/// The nv index of a vertex within its parent node
	pub vert_par_nv: Vec<i32>,
	// TODO: Should we store a vert_nodes CSR?
	pub node_edges: Vec<NodeEdge>,
	pub node_nes: CsrIndex,
	pub node_adj: Csr<NodeAdj>,
}

#[repr(u8)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum NodeType {
	#[default]
	F = b'F',
	V = b'V',
	Q = b'Q',
	I = b'I',
	O = b'O',
	S = b'S',
	P = b'P',
	R = b'R',
}

impl fmt::Display for NodeType {
	fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
		write!(f, "{}", *self as u8 as char)
	}
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct NodeVert {
	pub node: i32,
	pub vert: i32,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct NodeEdge {
	pub node: i32,
	pub twin_ne: i32,
	// TODO: Should we store the twin node, the twin node type, and/or twin node type == Q?
	pub nvs: [i32; 2],
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct NodeAdj {
	pub ne: i32,
	pub dest_nv: i32,
}

impl SpqrTree {
	pub fn size(&self) -> i32 {
		self.par.len() as i32
	}

	/// `vert_order` and `edge_order` are (prefixes of) permutations of vertex / edge ids;
	/// listed ids are visited first in the given order, then the rest in id order.
	/// Roots are the first unvisited vertices, and DFS children are explored in edge order.
	/// Use [`PlanarSpqrTree::build`] to also compute the planar embeddings.
	pub fn build(nv: i32, edges: &[[i32; 2]], ternarize: bool, vert_order: &[i32], edge_order: &[i32]) -> SpqrTree {
		build_impl::<false>(nv, edges, ternarize, vert_order, edge_order).tree
	}
}

/// Quarter-edges are indexed according to 4 * edge + 2 * side + dir, where side is v0 vs v1, and dir is cw vs ccw:
///
///       1     2
///    v0 ---e--- v1
///       0     3
///
/// A planar embedding is 3 involutions on quarter-edges:
/// * qe <-> qe ^ 1 maps quarter edges to their opposite side around the endpoint vertex.
/// * qe <-> qe ^ 3 maps quarter edges to their opposite side along the edge (around the face).
/// * qe <-> rot_adj[qe] maps quarter edges to their facing pair.
///
/// Walking around a vertex is alternating qe ^ 1 and rot_adj[qe], and walking around a face is qe ^ 3 and rot_adj[qe].
///
/// Partial embeddings are represented with -1's in the rot_adj array.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PlanarEmbedding {
	pub rot_adj: Vec<i32>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PlanarSpqrTree {
	pub tree: SpqrTree,
	pub node_planar: Vec<bool>,
	/// Planarity adjacencies: ne_embedding.rot_adj is an involution of facing quarter-edges, indexed according to:
	/// ne_embedding.rot_adj[4 * node_edge + 2 * side + dir]
	/// Nonplanar nodes have all entries -1.
	pub ne_embedding: PlanarEmbedding,
}

impl Deref for PlanarSpqrTree {
	type Target = SpqrTree;
	fn deref(&self) -> &SpqrTree {
		&self.tree
	}
}

impl DerefMut for PlanarSpqrTree {
	fn deref_mut(&mut self) -> &mut SpqrTree {
		&mut self.tree
	}
}

impl PlanarSpqrTree {
	pub fn build(nv: i32, edges: &[[i32; 2]], ternarize: bool, vert_order: &[i32], edge_order: &[i32]) -> PlanarSpqrTree {
		build_impl::<true>(nv, edges, ternarize, vert_order, edge_order)
	}
}

fn min(a: i32, b: i32) -> i32 {
	if a < b { a } else { b }
}
fn setmin(a: &mut i32, b: i32) {
	if b < *a {
		*a = b;
	}
}

/// Calls f(i) for i in order, then for the remaining i in [0, n) in increasing order.
fn for_each_in_order(n: i32, order: &[i32], mut f: impl FnMut(i32)) {
	for &i in order {
		f(i);
	}
	if order.len() as i32 == n {
		return;
	}
	if order.is_empty() {
		for i in 0..n {
			f(i);
		}
	} else if order.len() == 1 {
		for i in 0..n {
			if i != order[0] {
				f(i);
			}
		}
	} else {
		let mut listed = vec![false; n as usize];
		for &i in order {
			listed[i as usize] = true;
		}
		for i in 0..n {
			if !listed[i as usize] {
				f(i);
			}
		}
	}
}

// Helpers for working with [T; 2] - these compile to cmov's better than direct index access.

/// return arr[dir] == a, arr[!dir] == b
fn set_sides<T>(dir: bool, a: T, b: T) -> [T; 2] {
	if dir { [b, a] } else { [a, b] }
}
fn get_side<T: Copy>(a: [T; 2], dir: bool) -> T {
	if dir { a[1] } else { a[0] }
}

// ---------------------------------------------------------------------------------------------
// Phase 1: build a sorted skeleton

#[derive(Clone, Copy, Default)]
struct OutEdge {
	src: i32,
	dest: i32,
	e: i32,
	key: i32,
}

#[derive(Clone, Copy, Default)]
struct AdjEdge {
	dest: i32,
	e: i32,
}

/// Return the 2 lowvals from this subtree
struct LowvalStack {
	cur: i32,
	prv_e: i32,
	lowvals: [i32; 2],
	ch_idx: i32,
	ch_end: i32,
}

struct LowvalDfs<'a> {
	adj: &'a Csr<AdjEdge>,
	depth: Vec<i32>,
	all_outedges: Vec<OutEdge>,
	stk: Vec<LowvalStack>,
}

impl LowvalDfs<'_> {
	fn push_vert(&mut self, cur: i32, prv_e: i32) {
		let d = self.stk.len() as i32;
		self.depth[cur as usize] = d;
		self.stk.push(LowvalStack {
			cur,
			prv_e,
			lowvals: [d, d],
			ch_idx: self.adj.bounds[cur as usize],
			ch_end: self.adj.bounds[cur as usize + 1],
		});
	}

	fn finish_edge(&mut self, is_tree: bool, n_lowvals: [i32; 2]) {
		let d = self.stk.len() as i32 - 1;
		let s = self.stk.last_mut().unwrap();
		let cur = s.cur;
		debug_assert!(s.ch_idx < s.ch_end);
		let AdjEdge { dest: nxt, e } = self.adj.dat[s.ch_idx as usize];
		let lowvals = &mut s.lowvals;
		s.ch_idx += 1;

		{
			// Extra bit is 0 for type-1 children, 1 for backedges, 2 for children with lowval2
			// Bridges have lowval -2 (kind 0), and components loops have lowval -1 (components are kind 0, loops are kind 1)
			// We don't really need to distinguish backedges vs type-1 children, but do it just for fun?
			let mut lowval = n_lowvals[0];
			if lowval >= d {
				lowval = !(lowval - d);
			}
			let kind = 2 * (n_lowvals[1] < d) as i32 + (!is_tree) as i32;
			self.all_outedges.push(OutEdge { src: cur, dest: nxt, e, key: 3 * (lowval + 2) + kind });
		}

		// Keep the 2 distinct mins
		if n_lowvals[0] < lowvals[0] {
			*lowvals = [n_lowvals[0], min(n_lowvals[1], lowvals[0])];
		} else {
			lowvals[1] = min(lowvals[1], if n_lowvals[0] == lowvals[0] { n_lowvals[1] } else { n_lowvals[0] });
		}
	}

	fn start_edge(&mut self) {
		let d = self.stk.len() as i32 - 1;
		let s = self.stk.last_mut().unwrap();
		debug_assert!(s.ch_idx < s.ch_end);
		let AdjEdge { dest: nxt, e } = self.adj.dat[s.ch_idx as usize];

		if e == s.prv_e || self.depth[nxt as usize] > d {
			// skip the edge
			s.ch_idx += 1;
			return;
		}

		let is_tree = self.depth[nxt as usize] == -1;
		if is_tree {
			self.push_vert(nxt, e);
		} else {
			self.finish_edge(false, [self.depth[nxt as usize], d]);
		}
	}

	fn pop_vert(&mut self) -> [i32; 2] {
		self.stk.pop().unwrap().lowvals
	}
}

// ---------------------------------------------------------------------------------------------
// Phase 2: do the big ear-decomposition-like walk

// We're going to build a tree of all SPQR *nodes* + all original *vertices* (collectively *items*).
// Vertices will hang off the first SPQR node containing them, and blocks will be rooted at a topmost Q node for the top edge.

// As we build, we will represent the children of our nodes/vertices as linked lists.
const ROOT_ITEM: i32 = 0;

#[derive(Clone, Copy, Debug)]
struct ItemList {
	// Items are actually 2 * item + planarity_flip (always 0 without planarity)
	v: [i32; 2],
}

impl Default for ItemList {
	fn default() -> Self {
		ItemList { v: [-1, -1] }
	}
}

impl ItemList {
	#[must_use]
	fn empty(&self) -> bool {
		self.v[0] < 0
	}
}

fn unit_list(item: i32) -> ItemList {
	ItemList { v: [item << 1, item << 1] }
}

#[derive(Clone, Copy, Debug, Default)]
struct NonplanarityCertificate;

#[derive(Clone, Copy, Debug)]
struct TstackPlanaritySide {
	// For each side, store pointers to the "linked lists" of the edges inside.
	// v[0] is the outer / longer edges and v[1] is the inner / shorter edges, matching the outside-in sort order.

	// bot_ends are the outer/innermost exposed pieces of the walk down the ear in the tree (they're connected to the bottommost/topmost vertices of the tree path)
	bot_ends: [i32; 2],
	// top_ends are the outer/innermost exposed backedges
	top_ends: [i32; 2],
	// depths should be increasing going inwards
	top_depths: [i32; 2],
}

impl Default for TstackPlanaritySide {
	fn default() -> Self {
		TstackPlanaritySide { bot_ends: [-1, -1], top_ends: [-1, -1], top_depths: [-1, -1] }
	}
}

#[derive(Clone, Copy, Debug, Default)]
struct TstackPlanarity {
	// The convention is that sides[0].top_depths[0] == top_depth, i.e. at least one minimal return lives on side 0
	sides: [TstackPlanaritySide; 2],
}

#[derive(Clone, Copy, Debug, Default)]
struct TstackNonplanarity {
	// TODO: What's the nonplanarity certificate look like?
}

// NB: unlike the C++ `std::conditional_t<with_planarity, std::expected<...>, std::monostate>`,
// this is stored (but never touched) in the non-planar build too.
type TstackMaybePlanarity = Result<TstackPlanarity, TstackNonplanarity>;

#[derive(Clone, Copy, Debug)]
struct Tstack {
	v_start: i32,
	top_depth: i32,
	first_idx: i32,
	spans: [ItemList; 2],
	planarity: TstackMaybePlanarity,
}

#[derive(Clone, Copy)]
struct Key {
	lowval: i32,
	is_tree: bool,
	is_type_1: bool,
}

fn decode_key(cur_depth: i32, key: i32) -> Key {
	let mut lowval = key / 3 - 2;
	if lowval < 0 {
		lowval = cur_depth + !lowval;
	}
	let kind = key % 3;
	let is_tree = kind != 1;
	let is_type_1 = kind <= 1;
	Key { lowval, is_tree, is_type_1 }
}

#[derive(Clone, Copy)]
struct WalkStack {
	has_vert_tstack: bool,
	ch_idx: i32,
	ch_end: i32,
	orig_tstack: i32,
}

fn merge_planarity(quarter_edge_matches: &mut [i32], a: &mut TstackMaybePlanarity, b: &TstackMaybePlanarity) {
	let Ok(ap) = a else { return };
	let Ok(bp) = b else {
		*a = *b;
		return;
	};
	for z in 0..2 {
		let as_ = &mut ap.sides[z];
		let bs = &bp.sides[z];
		// If there's no bottom edges, then we must be an isolated vertex, so we can end early.
		if bs.bot_ends[0] == -1 {
			// Do nothing
		} else if as_.bot_ends[0] == -1 {
			*as_ = *bs;
		} else {
			quarter_edge_matches[as_.bot_ends[1] as usize] = bs.bot_ends[0];
			quarter_edge_matches[bs.bot_ends[0] as usize] = as_.bot_ends[1];
			as_.bot_ends[1] = bs.bot_ends[1];

			if bs.top_ends[0] == -1 {
				// Do nothing
			} else if as_.top_ends[0] == -1 {
				as_.top_ends = bs.top_ends;
				as_.top_depths = bs.top_depths;
			} else if as_.top_depths[1] > bs.top_depths[0] {
				// TODO: Certificate
				*a = Err(TstackNonplanarity {});
				return;
			} else {
				quarter_edge_matches[as_.top_ends[1] as usize] = bs.top_ends[0];
				quarter_edge_matches[bs.top_ends[0] as usize] = as_.top_ends[1];
				as_.top_ends[1] = bs.top_ends[1];
				as_.top_depths[1] = bs.top_depths[1];
			}
		}
	}
}

struct Builder<const WP: bool> {
	nv: i32,
	ne: i32,
	ternarize: bool,
	outedges: Csr<OutEdge>,

	ch_nxt: Vec<i32>,
	item_vs: Vec<[i32; 2]>,
	item_ch: Vec<ItemList>,
	item_types: Vec<NodeType>,

	// Quarter edges for planar embedding building.
	// Each vedge has 4 entries by 4 * vedge_id + 2 * source_vert + is_cw (is_cw is arbitrary)
	// vedges are identified with what item they cap, numbered by (item - 1 - NV)
	quarter_edge_matches: Vec<i32>,
	node_planarity: Vec<Result<[i32; 4], NonplanarityCertificate>>,

	tot_blocks: i32,
	tot_self_loops: i32,

	// Declare these here: most of our code will be in terms of v_start / top_depth, so we'll want to read these out
	stack_verts: Vec<i32>,
	stack_dir: Vec<bool>,

	nxt_edge_idx: i32, // Counts backedges only
	first_occurrence: Vec<i32>, // First backedge to this depth

	edge_top_depths: Vec<i32>,

	tstack: Vec<Tstack>,
	stk: Vec<WalkStack>,
}

impl<const WP: bool> Builder<WP> {
	fn vert_item(&self, v: i32) -> i32 {
		1 + v
	}
	fn edge_item(&self, e: i32) -> i32 {
		1 + self.nv + e
	}

	fn concat(&mut self, a: ItemList, b: ItemList) -> ItemList {
		if b.empty() {
			return a;
		}
		if a.empty() {
			return b;
		}
		self.ch_nxt[(a.v[1] >> 1) as usize] = b.v[0] ^ (a.v[1] & 1);
		ItemList { v: [a.v[0], b.v[1]] }
	}

	fn alloc_item(&mut self, ty: NodeType) -> i32 {
		let item = self.item_vs.len() as i32;
		self.item_vs.push([0, 0]);
		self.item_ch.push(ItemList::default());
		self.item_types.push(ty);
		self.ch_nxt.push(-1);
		if WP {
			self.node_planarity.push(Ok([0; 4]));
		}
		item
	}

	fn make_vs(&self, v_start: i32, top_depth: i32) -> [i32; 2] {
		set_sides(self.stack_dir[top_depth as usize], self.stack_verts[top_depth as usize], v_start)
	}

	fn make_edge_planarity(&mut self, item: i32, top_depth: i32, is_tree: bool) -> TstackMaybePlanarity {
		if WP {
			debug_assert!(item >= 1 + self.nv);
			let ve = item - (1 + self.nv);
			let top_dir = self.stack_dir[top_depth as usize] as i32;
			self.edge_top_depths[ve as usize] = top_depth;
			let mut p = TstackPlanarity::default();
			if is_tree {
				p.sides[0].bot_ends = [4 * ve + 2 * (1 - top_dir) + 0, 4 * ve + 2 * top_dir + 1];
				p.sides[1].bot_ends = [4 * ve + 2 * (1 - top_dir) + 1, 4 * ve + 2 * top_dir + 0];
			} else {
				p.sides[0].bot_ends = [4 * ve + 2 * (1 - top_dir) + 0, 4 * ve + 2 * (1 - top_dir) + 1];
				p.sides[0].top_ends = [4 * ve + 2 * top_dir + 1, 4 * ve + 2 * top_dir + 0];
				p.sides[0].top_depths = [top_depth, top_depth];
			}
			Ok(p)
		} else {
			Ok(TstackPlanarity::default())
		}
	}

	fn cur_tstack(&self) -> &Tstack {
		&self.tstack[self.tstack.len() - 1]
	}
	fn nxt_tstack(&self) -> &Tstack {
		&self.tstack[self.tstack.len() - 2]
	}

	fn push_tstack(&mut self, v_start: i32, top_depth: i32, item: i32, planarity: TstackMaybePlanarity) {
		let spans = set_sides(self.stack_dir[top_depth as usize], unit_list(item), ItemList::default());
		self.tstack.push(Tstack { v_start, top_depth, first_idx: self.nxt_edge_idx, spans, planarity });
	}
	fn push_vert_tstack(&mut self, v: i32, top_depth: i32) {
		let item = self.vert_item(v);
		self.push_tstack(v, top_depth, item, Ok(TstackPlanarity::default()));
	}
	fn push_edge_tstack(&mut self, v_start: i32, top_depth: i32, e: i32, is_tree: bool) {
		let item = self.edge_item(e);
		let planarity = self.make_edge_planarity(item, top_depth, is_tree);
		self.push_tstack(v_start, top_depth, item, planarity);
	}
	fn flip_tstack_planarity(a: &mut Tstack) {
		if WP {
			a.spans[0].v[0] ^= 1;
			a.spans[0].v[1] ^= 1;
			a.spans[1].v[0] ^= 1;
			a.spans[1].v[1] ^= 1;
			if let Ok(p) = &mut a.planarity {
				p.sides.swap(0, 1);
			}
		}
	}
	fn merge_tstack_tops(&mut self) {
		let b = self.tstack.pop().unwrap();
		let ai = self.tstack.len() - 1;
		setmin(&mut self.tstack[ai].top_depth, b.top_depth);
		let s0 = self.concat(b.spans[0], self.tstack[ai].spans[0]);
		self.tstack[ai].spans[0] = s0;
		let s1 = self.concat(self.tstack[ai].spans[1], b.spans[1]);
		self.tstack[ai].spans[1] = s1;
		if WP {
			merge_planarity(&mut self.quarter_edge_matches, &mut self.tstack[ai].planarity, &b.planarity);
		}
	}

	fn maybe_unwrap_nxt(&mut self, ty: NodeType, is_tree: bool) -> i32 {
		let ti = self.tstack.len() - 2;

		if ty == NodeType::R {
			return self.alloc_item(ty);
		}

		debug_assert!(ty == NodeType::P || ty == NodeType::S);

		// If we want to ternarize, never reuse.
		if self.ternarize {
			return self.alloc_item(ty);
		}

		let t = self.tstack[ti];
		let top_dir = self.stack_dir[t.top_depth as usize];
		debug_assert!(get_side(t.spans, !top_dir).empty());
		let item = get_side(t.spans, top_dir).v[0] >> 1;
		debug_assert!(item == (get_side(t.spans, top_dir).v[1] >> 1));
		if self.item_types[item as usize] == ty {
			self.tstack[ti].spans = set_sides(top_dir, self.item_ch[item as usize], ItemList::default());
			if WP {
				// Unwrap the planarity data
				// We don't really need to maintain this at all because S/P nodes are known to be trivially planar
				// The current state is just make_edge_planarity(wrapped), which means that it has the right shape, just needs to be relabelled.
				let matches = self.node_planarity[(item - (1 + self.nv + self.ne)) as usize]
					.expect("unwrapped S/P nodes are always planar");
				let p = self.tstack[ti].planarity.as_mut().expect("unwrapped tstack is a single edge");
				let top_dir = top_dir as usize;
				if is_tree {
					p.sides[0].bot_ends[0] = matches[2 * (1 - top_dir) + 1];
					p.sides[0].bot_ends[1] = matches[2 * top_dir + 0];
					p.sides[1].bot_ends[0] = matches[2 * (1 - top_dir) + 0];
					p.sides[1].bot_ends[1] = matches[2 * top_dir + 1];
				} else {
					p.sides[0].bot_ends[0] = matches[2 * (1 - top_dir) + 1];
					p.sides[0].bot_ends[1] = matches[2 * (1 - top_dir) + 0];
					p.sides[0].top_ends[0] = matches[2 * top_dir + 0];
					p.sides[0].top_ends[1] = matches[2 * top_dir + 1];
				}
			}
			item
		} else {
			self.alloc_item(ty)
		}
	}

	fn finish_tstack_top(&mut self, item: i32, is_tree: bool) {
		let ti = self.tstack.len() - 1;
		let t = self.tstack[ti];
		let top_dir = self.stack_dir[t.top_depth as usize];
		debug_assert!(get_side(t.spans, !top_dir).empty());

		if WP {
			let np = (item - (1 + self.nv + self.ne)) as usize;
			if let Ok(p) = &t.planarity {
				let top_dir = top_dir as usize;
				let mut matches = [0; 4];
				if is_tree {
					matches[2 * (1 - top_dir) + 1] = p.sides[0].bot_ends[0];
					matches[2 * top_dir + 0] = p.sides[0].bot_ends[1];
					matches[2 * (1 - top_dir) + 0] = p.sides[1].bot_ends[0];
					matches[2 * top_dir + 1] = p.sides[1].bot_ends[1];
				} else {
					matches[2 * (1 - top_dir) + 1] = p.sides[0].bot_ends[0];
					matches[2 * (1 - top_dir) + 0] = p.sides[0].bot_ends[1];
					matches[2 * top_dir + 0] = p.sides[0].top_ends[0];
					matches[2 * top_dir + 1] = p.sides[0].top_ends[1];
				}
				self.node_planarity[np] = Ok(matches);
			} else {
				debug_assert!(self.item_types[item as usize] == NodeType::R);
				self.node_planarity[np] = Err(NonplanarityCertificate);
			}
		}
		self.item_vs[item as usize] = self.make_vs(t.v_start, t.top_depth);
		self.item_ch[item as usize] = get_side(t.spans, top_dir);

		self.tstack[ti].spans = set_sides(top_dir, unit_list(item), ItemList::default());
		self.tstack[ti].planarity = self.make_edge_planarity(item, t.top_depth, is_tree);
	}

	fn push_vert(&mut self, cur: i32) {
		self.stk.push(WalkStack {
			has_vert_tstack: false,
			ch_idx: self.outedges.bounds[cur as usize],
			ch_end: self.outedges.bounds[cur as usize + 1],
			orig_tstack: -1,
		});
		let cur_depth = self.stk.len() - 1;
		self.stack_verts[cur_depth] = cur;
	}

	/// Some(nxt) means jump to push_vert(nxt), None means jump to finish_edge
	fn start_edge(&mut self) -> Option<i32> {
		let cur_depth = self.stk.len() as i32 - 1;
		let si = self.stk.len() - 1;
		let s = self.stk[si];
		let cur = self.stack_verts[cur_depth as usize];
		debug_assert!(s.ch_idx < s.ch_end);
		let OutEdge { dest: nxt, key, .. } = self.outedges.dat[s.ch_idx as usize];
		let Key { lowval, is_tree, is_type_1 } = decode_key(cur_depth, key);

		// edge_dir convention: false is forwards, true is backwards.
		// That means that cur is on the edge_dir side and nxt is on the !edge_dir side.
		self.stack_dir[cur_depth as usize] = if lowval >= cur_depth { false } else { !self.stack_dir[lowval as usize] };

		if !s.has_vert_tstack && lowval < cur_depth && is_type_1 {
			// Do this with the correct stack_dir set
			self.push_vert_tstack(cur, cur_depth);
			self.stk[si].has_vert_tstack = true;
		}

		self.stk[si].orig_tstack = self.tstack.len() as i32;
		if is_tree {
			self.first_occurrence[cur_depth as usize] = self.ne;
			Some(nxt)
		} else {
			None
		}
	}

	fn finish_edge(&mut self) {
		let cur_depth = self.stk.len() as i32 - 1;
		let si = self.stk.len() - 1;
		let cur = self.stack_verts[cur_depth as usize];
		debug_assert!(self.stk[si].ch_idx < self.stk[si].ch_end);

		let OutEdge { dest: nxt, e, key, .. } = self.outedges.dat[self.stk[si].ch_idx as usize];
		self.stk[si].ch_idx += 1;

		let Key { lowval, is_tree, is_type_1 } = decode_key(cur_depth, key);

		let orig_tstack = self.stk[si].orig_tstack;
		let edge_dir = self.stack_dir[cur_depth as usize];
		let ei = self.edge_item(e) as usize;

		if lowval >= cur_depth {
			// There's no planarity handling for this because it's just a Q node. I/O nodes also don't need any tracking.
			self.item_vs[ei] = [cur, -1];
			self.tot_blocks += 1;
			if is_tree {
				// Bridges and components
				if lowval == cur_depth + 1 {
					// tstack[tstack_size-1] is currently just smuggling out the child vertex, prepend the bridge component
					// This is just a shortcut for allocating a full I-type tstack
					let item = self.alloc_item(NodeType::I);
					self.item_vs[item as usize] = self.make_vs(nxt, cur_depth);
					let t = self.tstack.pop().unwrap();
					self.item_ch[ei] = self.concat(unit_list(item), t.spans[1]);
				} else {
					// tstack[tstack_size-2] is the vertex and tstack[tstack_size-1] is the backedge
					let backedge = self.tstack.pop().unwrap().spans[0];
					let t = self.tstack.pop().unwrap();
					self.item_ch[ei] = self.concat(backedge, t.spans[1]);
				}
			} else {
				// self loops
				debug_assert!(nxt == cur);
				self.tot_self_loops += 1;
				let item = self.alloc_item(NodeType::O);
				// Make sure the nxt is -1 as well
				self.item_vs[item as usize] = [cur, -1];
				self.item_ch[ei] = unit_list(item);
			}
			let vi = self.vert_item(cur) as usize;
			let l = self.concat(self.item_ch[vi], unit_list(ei as i32));
			self.item_ch[vi] = l;
			return;
		}
		debug_assert!(lowval < cur_depth);

		self.item_vs[ei] = self.make_vs(nxt, cur_depth);

		// Whether cur_tstack() is a single edge
		let mut is_single = true;
		if is_tree {
			// The span lives on side edge_dir
			self.push_edge_tstack(nxt, cur_depth, e, true);
			while self.tstack.len() >= 2 && self.nxt_tstack().top_depth >= cur_depth {
				let ty;
				if self.nxt_tstack().top_depth > cur_depth {
					// Just backfill this for maybe_unwrap
					let td = self.nxt_tstack().top_depth;
					self.stack_dir[td as usize] = edge_dir;

					// The tstack currently contains a tree-edge followed by a vertex; merge the vertex first
					self.merge_tstack_tops();

					ty = NodeType::S;
				} else if self.nxt_tstack().v_start == self.cur_tstack().v_start {
					// This will be a P node
					ty = NodeType::P;
				} else {
					ty = NodeType::R;
				}
				let item = self.maybe_unwrap_nxt(ty, ty == NodeType::S);
				self.merge_tstack_tops();
				if WP {
					let ci = self.tstack.len() - 1;
					if let Ok(p) = &mut self.tstack[ci].planarity {
						// Merge all backedges into the component
						for side in &mut p.sides {
							debug_assert!(side.bot_ends[1] != -1);
							if side.top_ends[1] == -1 {
								continue;
							}
							debug_assert!(side.top_depths[0] == cur_depth);
							debug_assert!(side.top_depths[1] == cur_depth);
							self.quarter_edge_matches[side.bot_ends[1] as usize] = side.top_ends[1];
							self.quarter_edge_matches[side.top_ends[1] as usize] = side.bot_ends[1];
							side.bot_ends[1] = side.top_ends[0];
							side.top_depths = [-1, -1];
							side.top_ends = [-1, -1];
						}
					}
				}
				self.finish_tstack_top(item, true);
			}

			if self.cur_tstack().first_idx > self.first_occurrence[cur_depth as usize] {
				while self.cur_tstack().first_idx > self.first_occurrence[cur_depth as usize] {
					if WP {
						let n = self.tstack.len();
						if self.tstack[n - 2].first_idx > self.first_occurrence[cur_depth as usize] {
							// We will put cur_depth on side 1 until the bottom
							if self.tstack[n - 2].top_depth == cur_depth {
								Self::flip_tstack_planarity(&mut self.tstack[n - 2]);
							}
						} else if !is_single {
							debug_assert!(self.tstack[n - 1].top_depth < cur_depth);
							if let Ok(p) = self.tstack[n - 2].planarity {
								if p.sides[0].top_depths[1] == cur_depth {
									// We need to flip cur_tstack and nxt_tstack relative to each other.
									// Flip the one with worse top_depth.
									let which = if self.tstack[n - 1].top_depth < self.tstack[n - 2].top_depth { n - 2 } else { n - 1 };
									Self::flip_tstack_planarity(&mut self.tstack[which]);
								} else {
									debug_assert!(p.sides[1].top_depths[1] == cur_depth);
								}
							}
						}
					}
					self.merge_tstack_tops();
					is_single = false;
				}
				if WP {
					let ci = self.tstack.len() - 1;
					if let Ok(p) = &mut self.tstack[ci].planarity {
						// Prune off finished cur-side things
						for side in &mut p.sides {
							debug_assert!(side.bot_ends[1] != -1);
							while side.top_depths[1] == cur_depth {
								{
									// Link these to bot_ends[1]
									self.quarter_edge_matches[side.bot_ends[1] as usize] = side.top_ends[1];
									self.quarter_edge_matches[side.top_ends[1] as usize] = side.bot_ends[1];
									side.bot_ends[1] = side.top_ends[1] ^ 1;
								}
								side.top_ends[1] = std::mem::replace(&mut self.quarter_edge_matches[side.bot_ends[1] as usize], -1);
								if side.top_ends[1] != -1 {
									self.quarter_edge_matches[side.top_ends[1] as usize] = -1;
									side.top_depths[1] = self.edge_top_depths[(side.top_ends[1] >> 2) as usize];
								} else {
									side.top_depths = [-1, -1];
									side.top_ends = [-1, -1];
								}
							}
						}
					}
				}
			}

			if is_type_1 {
				debug_assert!(self.stk[si].has_vert_tstack);
			}
			if self.stk[si].has_vert_tstack {
				// NB: tstack[orig_size] is the vertex and tstack[orig_size+1] is the backedge; maybe we should reverse them?
				debug_assert!(self.tstack.len() as i32 >= orig_tstack + 3);

				if !is_type_1 {
					if WP {
						// The lowval side should be side 1, everything else goes on side 0.
						// The exception is tstack[orig_tstack + 2], which could be == lowval on one/both sides,
						// but is guaranteed to have *something* > lowval by non-type-1-ness
						let ti = (orig_tstack + 2) as usize;
						if let Ok(p) = self.tstack[ti].planarity {
							debug_assert!(p.sides[0].top_depths[0] == self.tstack[ti].top_depth);
							if p.sides[0].top_depths[1] == lowval {
								Self::flip_tstack_planarity(&mut self.tstack[ti]);
							}
							let p = self.tstack[ti].planarity.as_ref().unwrap();
							debug_assert!(p.sides[0].top_depths[1] != -1);
							debug_assert!(p.sides[0].top_depths[1] > lowval);
						}
						for i in (orig_tstack + 3) as usize..self.tstack.len() {
							if self.tstack[i].top_depth == lowval {
								Self::flip_tstack_planarity(&mut self.tstack[i]);
							}
						}
					}
					while self.tstack.len() as i32 > orig_tstack + 3 {
						self.merge_tstack_tops();
						is_single = false;
					}
					debug_assert!(!is_single);
				}

				debug_assert!(self.tstack.len() as i32 == orig_tstack + 3);
				let item = if is_type_1 {
					self.maybe_unwrap_nxt(if is_single { NodeType::S } else { NodeType::R }, false)
				} else {
					// Just for the type checker
					-1
				};
				// Merge with the backedge
				self.merge_tstack_tops();
				// Merge with the vertex
				self.merge_tstack_tops();

				let ci = self.tstack.len() - 1;
				self.tstack[ci].v_start = cur;
				debug_assert!(self.tstack[ci].top_depth == lowval);

				// Fold everything to the correct side now that we're leaving the child.
				// The entire subtree should go to the !edge_dir side.
				let all = self.concat(self.tstack[ci].spans[0], self.tstack[ci].spans[1]);
				self.tstack[ci].spans = set_sides(!edge_dir, all, ItemList::default());

				if WP {
					if let Ok(p) = &mut self.tstack[ci].planarity {
						// precondition: side 1 should be the lowval only side
						let [s0, s1] = &mut p.sides;
						self.quarter_edge_matches[s0.bot_ends[0] as usize] = s1.bot_ends[0];
						self.quarter_edge_matches[s1.bot_ends[0] as usize] = s0.bot_ends[0];
						s0.bot_ends[0] = s1.bot_ends[1];
						let mut nonplanar = false;
						if s1.top_ends[0] != -1 {
							if s1.top_depths[1] != lowval {
								debug_assert!(!is_type_1);
								nonplanar = true;
							} else {
								debug_assert!(s1.top_depths[0] == lowval);
								self.quarter_edge_matches[s0.top_ends[0] as usize] = s1.top_ends[0];
								self.quarter_edge_matches[s1.top_ends[0] as usize] = s0.top_ends[0];
								s0.top_ends[0] = s1.top_ends[1];
								// Already true since the backedge was on side 0
								debug_assert!(s0.top_depths[0] == lowval);
							}
						}
						if nonplanar {
							self.tstack[ci].planarity = Err(TstackNonplanarity {});
						} else {
							*s1 = TstackPlanaritySide::default();
						}
					}
				}

				if is_type_1 {
					self.finish_tstack_top(item, false);
					is_single = true;
				}
			}
		} else {
			debug_assert!(is_type_1);
			// The span lives on side !edge_dir
			self.push_edge_tstack(cur, lowval, e, false);
			let idx = self.nxt_edge_idx;
			self.nxt_edge_idx += 1;
			setmin(&mut self.first_occurrence[lowval as usize], idx);
		}

		// NB: We can do this check in lots of ways, maybe there's a cleaner check
		if is_type_1 && self.tstack.len() >= 2 && self.nxt_tstack().v_start == cur && self.nxt_tstack().top_depth == lowval {
			// This will be a P node
			let item = self.maybe_unwrap_nxt(NodeType::P, false);
			self.merge_tstack_tops();
			self.finish_tstack_top(item, false);
		}

		if !self.stk[si].has_vert_tstack {
			// Throw cur_vert_node onto the tstack so it'll get interleaved correctly
			self.push_vert_tstack(cur, cur_depth);
			self.stk[si].has_vert_tstack = true;
			debug_assert!(!is_type_1);
			if !is_single {
				// Just eagerly merge the vertex into the R to avoid a later spurious finish_tstack
				self.merge_tstack_tops();
			}
		}
	}

	fn pop_vert(&mut self) {
		let cur_depth = self.stk.len() as i32 - 1;
		let si = self.stk.len() - 1;
		let cur = self.stack_verts[cur_depth as usize];
		debug_assert!(self.stk[si].ch_idx == self.stk[si].ch_end);
		if !self.stk[si].has_vert_tstack {
			// Either our parent is a bridge edge, or we're just a root.
			// We'll just leave it on tstack for future cleanup, it'll just get popped of immediately.
			// edge_dir == !stack_dir[lowval == cur_depth - 1] == true
			self.stack_dir[cur_depth as usize] = true;
			self.push_vert_tstack(cur, cur_depth);
			self.stk[si].has_vert_tstack = true;
		}
		self.stk.pop();
	}
}

// ---------------------------------------------------------------------------------------------
// Phase 3: relabel the full tree in preorder

#[derive(Clone, Copy, Default)]
struct ChBuf {
	loc: i32,
	item_id: i32,
}

#[derive(Clone, Copy)]
struct RelabelStack {
	cur_idx: i32,
	ch_idx: i32,
	ch_end: i32,
	cur_nv: i32,
	cur_ne: i32,
}

struct Relabel<'a, const WP: bool> {
	b: Builder<WP>,
	edges: &'a [[i32; 2]],

	vert_index: Vec<i32>,
	edge_index: Vec<i32>,
	edge_flipped: Vec<bool>,
	par: Vec<i32>,
	subtree_end: Vec<i32>,
	types: Vec<NodeType>,
	orig_id: Vec<i32>,
	ch: Csr<i32>,
	node_verts: Vec<NodeVert>,
	node_nvs: CsrIndex,
	vert_par_nv: Vec<i32>,
	node_edges: Vec<NodeEdge>,
	node_nes: CsrIndex,
	node_adj: Csr<NodeAdj>,
	node_planar: Vec<bool>,
	ne_rot_adj: Vec<i32>,

	vert_pos_buf: Vec<i32>,
	cnts_buf: Vec<i32>,
	ch_buf: Vec<ChBuf>,
	rot_edge_ne: Vec<i32>,

	nxt_unassigned_idx: i32,
	stk: Vec<RelabelStack>,
}

impl<const WP: bool> Relabel<'_, WP> {
	fn set_ne(&mut self, cur_idx: i32, ne: i32, nvs: [i32; 2], nds: [i32; 2], rot_adjs: [i32; 4]) {
		self.node_edges[ne as usize].node = cur_idx;
		self.node_edges[ne as usize].nvs = nvs;
		self.node_adj.dat[nds[0] as usize] = NodeAdj { ne, dest_nv: nvs[1] };
		self.node_adj.dat[nds[1] as usize] = NodeAdj { ne, dest_nv: nvs[0] };
		if WP {
			for z in 0..4 {
				self.ne_rot_adj[(4 * ne) as usize + z] = rot_adjs[z];
			}
		}
	}

	fn map_rot_edge(&self, planar: bool, ve: i32) -> [i32; 4] {
		if !WP {
			return [-1, -1, -1, -1];
		}
		if !planar {
			return [-1, -1, -1, -1];
		}
		let mut res = [0; 4];
		for z in 0..4 {
			let o = self.b.quarter_edge_matches[(4 * ve) as usize + z];
			debug_assert!(o != -1);
			res[z] = (self.rot_edge_ne[(o >> 2) as usize] << 2) + (o & 2) + ((z & 1) == 0) as i32;
		}
		res
	}

	fn push_item(&mut self, cur_item: i32) {
		let nv = self.b.nv;
		let ne = self.b.ne;
		let cur_idx = self.nxt_unassigned_idx;
		self.nxt_unassigned_idx += 1;
		let cur_type = self.b.item_types[cur_item as usize];
		self.types[cur_idx as usize] = cur_type;
		let mut planar = true;
		if cur_type == NodeType::F {
			debug_assert!(cur_item == 0);
		} else if cur_type == NodeType::V {
			debug_assert!(1 <= cur_item && cur_item < 1 + nv);
			let orig_vert = cur_item - 1;
			self.orig_id[cur_idx as usize] = orig_vert;
			self.vert_index[orig_vert as usize] = cur_idx;
		} else if cur_type == NodeType::Q {
			debug_assert!(1 + nv <= cur_item && cur_item < 1 + nv + ne);
			let orig_edge = cur_item - 1 - nv;
			self.orig_id[cur_idx as usize] = orig_edge;
			self.edge_index[orig_edge as usize] = cur_idx;
			debug_assert!(self.b.item_vs[cur_item as usize][0] != -1);
			self.edge_flipped[orig_edge as usize] = self.b.item_vs[cur_item as usize][0] != self.edges[orig_edge as usize][0];
		} else {
			debug_assert!(1 + nv + ne <= cur_item);
			if WP {
				if cur_type == NodeType::O || cur_type == NodeType::I {
					// No planarity data was set up
				} else if cur_type == NodeType::S || cur_type == NodeType::P || cur_type == NodeType::R {
					let p = self.b.node_planarity[(cur_item - (1 + nv + ne)) as usize];
					if let Ok(p) = p {
						// Make sure this runs before our planarity_flip checks
						for s in 0..4 {
							let a = 8 * ne + s as i32;
							let b = p[s];
							self.b.quarter_edge_matches[a as usize] = b;
							self.b.quarter_edge_matches[b as usize] = a;
						}
					} else {
						// TODO: Any certificate stuff
						planar = false;
					}
				} else {
					debug_assert!(false);
				}
			}
		}
		if WP {
			self.node_planar[cur_idx as usize] = planar;
		}

		// HACK: Fill ch and vert_items in with orig items / orig verts for now,
		// because we don't have the final item id's yet.
		let ch_st = self.ch.bounds[cur_idx as usize];
		let mut ch_en = ch_st;
		let nv_st = self.node_nvs.bounds[cur_idx as usize];
		let mut nv_en = nv_st;
		let mut n_edges = 0;
		let cur_item_vs = self.b.item_vs[cur_item as usize];
		if cur_item_vs[0] != -1 {
			self.node_verts[nv_en as usize] = NodeVert { node: cur_idx, vert: cur_item_vs[0] };
			nv_en += 1;
		}
		let cur_item_ch = self.b.item_ch[cur_item as usize];
		if !cur_item_ch.empty() {
			let mut planarity_flip = (cur_item_ch.v[0] & 1) != 0;
			let mut ch_item = cur_item_ch.v[0] >> 1;
			loop {
				self.ch.dat[ch_en as usize] = ch_item;
				ch_en += 1;
				debug_assert!(ch_item >= 1);
				if ch_item < 1 + nv {
					self.node_verts[nv_en as usize] = NodeVert { node: cur_idx, vert: ch_item - 1 };
					nv_en += 1;
				} else {
					if WP {
						if cur_type != NodeType::R {
							debug_assert!(!planarity_flip);
						} else {
							// Fix the planarity direction right here: reverse quarter_edge_matches upfront;
							// this breaks the involution property, but from here on we'll never read the low bits anyways.
							let ve = (ch_item - (1 + nv)) as usize;
							if planarity_flip {
								self.b.quarter_edge_matches.swap(4 * ve + 0, 4 * ve + 1);
								self.b.quarter_edge_matches.swap(4 * ve + 2, 4 * ve + 3);
							}
						}
					}
					n_edges += 1;
				}
				if ch_item == (cur_item_ch.v[1] >> 1) {
					debug_assert!(self.b.ch_nxt[ch_item as usize] == -1);
					break;
				}
				let nxt = self.b.ch_nxt[ch_item as usize];
				planarity_flip ^= (nxt & 1) != 0;
				ch_item = nxt >> 1;
			}
			planarity_flip ^= (cur_item_ch.v[1] & 1) != 0;
			debug_assert!(!planarity_flip);
		}
		if cur_item_vs[1] != -1 {
			self.node_verts[nv_en as usize] = NodeVert { node: cur_idx, vert: cur_item_vs[1] };
			nv_en += 1;
		}
		self.ch.bounds[cur_idx as usize + 1] = ch_en;
		self.node_nvs.bounds[cur_idx as usize + 1] = nv_en;

		let n_verts = nv_en - nv_st;

		let is_node = cur_type != NodeType::F && cur_type != NodeType::V;
		let has_cap = is_node && !(cur_type == NodeType::Q && ch_en - ch_st > 0);

		if !is_node {
			n_edges = 0;
		}
		if has_cap {
			n_edges += 1;
		}

		let ne_st = self.node_nes.bounds[cur_idx as usize];
		let ne_en = ne_st + n_edges;
		self.node_nes.bounds[cur_idx as usize + 1] = ne_en;

		if cur_type == NodeType::F {
			// Just set node_adj bounds and we're good
			for i in 2 * nv_st + 1..=2 * nv_en {
				self.node_adj.bounds[i as usize] = 2 * ne_st;
			}
		} else if cur_type == NodeType::V {
			// Nothing to do
		} else if n_verts == 1 {
			debug_assert!(cur_type == NodeType::Q || cur_type == NodeType::O);
			debug_assert!(n_edges == 1);
			self.node_adj.bounds[(2 * nv_st + 1) as usize] = 2 * ne_st + 1 * n_edges;
			self.node_adj.bounds[(2 * nv_st + 2) as usize] = 2 * ne_st + 2 * n_edges;
			self.set_ne(cur_idx, ne_st, [nv_st, nv_st], [2 * ne_st + 1, 2 * ne_st], [4 * ne_st + 3, 4 * ne_st + 2, 4 * ne_st + 1, 4 * ne_st + 0]);
		} else if cur_type == NodeType::Q || cur_type == NodeType::I {
			debug_assert!(n_verts == 2);
			debug_assert!(n_edges == 1);
			self.node_adj.bounds[(2 * nv_st + 1) as usize] = 2 * ne_st + 0 * n_edges;
			self.node_adj.bounds[(2 * nv_st + 2) as usize] = 2 * ne_st + 1 * n_edges;
			self.node_adj.bounds[(2 * nv_st + 3) as usize] = 2 * ne_st + 2 * n_edges;
			self.node_adj.bounds[(2 * nv_st + 4) as usize] = 2 * ne_st + 2 * n_edges;
			self.set_ne(cur_idx, ne_st, [nv_st, nv_st + 1], [2 * ne_st, 2 * ne_st + 1], [4 * ne_st + 1, 4 * ne_st + 0, 4 * ne_st + 3, 4 * ne_st + 2]);
		} else if cur_type == NodeType::P {
			// Special case: tiebreak the parallel edges so they're reversed
			debug_assert!(n_verts == 2);
			debug_assert!(n_edges >= 3);
			self.node_adj.bounds[(2 * nv_st + 1) as usize] = 2 * ne_st + 0 * n_edges;
			self.node_adj.bounds[(2 * nv_st + 2) as usize] = 2 * ne_st + 1 * n_edges;
			self.node_adj.bounds[(2 * nv_st + 3) as usize] = 2 * ne_st + 2 * n_edges;
			self.node_adj.bounds[(2 * nv_st + 4) as usize] = 2 * ne_st + 2 * n_edges;
			for ne in ne_st..ne_en {
				let ne_prv = (if ne == ne_st { ne_en } else { ne }) - 1;
				let ne_nxt = if ne + 1 == ne_en { ne_st } else { ne + 1 };
				let rot_adjs = [4 * ne_prv + 1, 4 * ne_nxt + 0, 4 * ne_nxt + 3, 4 * ne_prv + 2];
				self.set_ne(cur_idx, ne, [nv_st, nv_st + 1], [2 * ne_st + (ne - ne_st), 2 * ne_en - 1 - (ne - ne_st)], rot_adjs);
			}
		} else if cur_type == NodeType::S {
			debug_assert!(n_verts == n_edges);
			debug_assert!(n_verts >= 3);
			for i in 2 * nv_st + 1..=2 * nv_en {
				self.node_adj.bounds[i as usize] = i + 2 * (ne_st - nv_st);
			}
			// Fix bounds for the cap
			self.node_adj.bounds[(2 * nv_st + 1) as usize] -= 1;
			self.node_adj.bounds[(2 * nv_en - 1) as usize] += 1;
			self.set_ne(cur_idx, ne_st, [nv_st, nv_en - 1], [2 * ne_st, 2 * ne_en - 1], [4 * (ne_st + 1) + 1, 4 * (ne_st + 1) + 0, 4 * (ne_en - 1) + 3, 4 * (ne_en - 1) + 2]);
			for i in 1..n_edges {
				let ne = ne_st + i;
				let mut rot_adjs = [4 * (ne - 1) + 3, 4 * (ne - 1) + 2, 4 * (ne + 1) + 1, 4 * (ne + 1) + 0];
				if ne - 1 == ne_st {
					rot_adjs[0] = 4 * ne_st + 1;
					rot_adjs[1] = 4 * ne_st + 0;
				}
				if ne + 1 == ne_en {
					rot_adjs[2] = 4 * ne_st + 3;
					rot_adjs[3] = 4 * ne_st + 2;
				}
				self.set_ne(cur_idx, ne, [nv_st + i - 1, nv_st + i], [2 * ne - 1, 2 * ne], rot_adjs);
			}
		} else if cur_type == NodeType::R {
			// Bucketsort the children by the midpoint
			for nv_ in nv_st..nv_en {
				self.vert_pos_buf[self.node_verts[nv_ as usize].vert as usize] = nv_;
			}
			self.cnts_buf.clear();
			self.cnts_buf.resize((n_verts * 2 - 1) as usize, 0);
			self.ch_buf.clear();

			debug_assert!(has_cap);

			// Cap node_adj bounds
			self.node_adj.bounds[(2 * nv_st + 2) as usize] += 1;
			self.node_adj.bounds[(2 * nv_en - 1) as usize] += 1;

			for i in ch_st..ch_en {
				let item = self.ch.dat[i as usize];
				debug_assert!(item >= 1);
				let nvs: [i32; 2];
				if item < 1 + nv {
					nvs = [self.vert_pos_buf[(item - 1) as usize], self.vert_pos_buf[(item - 1) as usize]];
				} else {
					let vs = self.b.item_vs[item as usize];
					nvs = [self.vert_pos_buf[vs[0] as usize], self.vert_pos_buf[vs[1] as usize]];
					debug_assert!(nvs[0] < nvs[1]);
					self.node_adj.bounds[(2 * nvs[0] + 2) as usize] += 1;
					self.node_adj.bounds[(2 * nvs[1] + 1) as usize] += 1;
				}
				let loc = (nvs[0] - nv_st) + (nvs[1] - nv_st);
				self.ch_buf.push(ChBuf { loc, item_id: item });
				self.cnts_buf[loc as usize] += 1;
			}
			let mut offset = ch_st;
			for cnt in self.cnts_buf.iter_mut() {
				offset += *cnt;
				*cnt = offset;
			}
			for &ChBuf { loc, item_id: n } in self.ch_buf.iter().rev() {
				self.cnts_buf[loc as usize] -= 1;
				self.ch.dat[self.cnts_buf[loc as usize] as usize] = n;
			}

			if WP {
				// Set up the reverse mapping for ourselves
				let mut nxt_ne = ne_en;
				for i in (ch_st..ch_en).rev() {
					let item = self.ch.dat[i as usize];
					debug_assert!(item >= 1);
					if item < 1 + nv {
						continue;
					}
					nxt_ne -= 1;
					self.rot_edge_ne[(item - (1 + nv)) as usize] = nxt_ne;
				}
				debug_assert!(nxt_ne == ne_st + 1);
				self.rot_edge_ne[(2 * ne) as usize] = ne_st;
			}

			{
				let mut off = 2 * ne_st;
				for i in 2 * nv_st + 1..=2 * nv_en {
					off += std::mem::replace(&mut self.node_adj.bounds[i as usize], off);
				}
				debug_assert!(off == 2 * ne_en);
			}

			// Fill in node_edges and node_adj.
			// Reverse order to get the adj in bracket ordering.
			{
				// Handle cap as special: it's first in the node_edges, which means it's in the wrong place for the left endpoint.
				self.node_adj.bounds[(2 * nv_st + 2) as usize] += 1;

				let mut nxt_ne = ne_en;
				for i in (ch_st..ch_en).rev() {
					let item = self.ch.dat[i as usize];
					debug_assert!(item >= 1);
					if item < 1 + nv {
						continue;
					}
					nxt_ne -= 1;
					let [v0, v1] = self.b.item_vs[item as usize];
					// TODO: Reuse this from the ch pass?
					let nvs = [self.vert_pos_buf[v0 as usize], self.vert_pos_buf[v1 as usize]];
					let nd0 = self.node_adj.bounds[(2 * nvs[0] + 2) as usize];
					self.node_adj.bounds[(2 * nvs[0] + 2) as usize] += 1;
					let nd1 = self.node_adj.bounds[(2 * nvs[1] + 1) as usize];
					self.node_adj.bounds[(2 * nvs[1] + 1) as usize] += 1;
					let rot_adjs = self.map_rot_edge(planar, item - (1 + nv));
					self.set_ne(cur_idx, nxt_ne, nvs, [nd0, nd1], rot_adjs);
				}
				debug_assert!(nxt_ne == ne_st + 1);

				// Insert the cap / bump its bound
				let rot_adjs = self.map_rot_edge(planar, 2 * ne);
				self.set_ne(cur_idx, ne_st, [nv_st, nv_en - 1], [2 * ne_st, 2 * ne_en - 1], rot_adjs);
				self.node_adj.bounds[(2 * nv_en - 1) as usize] += 1;
			}
		} else {
			debug_assert!(false);
		}

		let cur_nv = nv_st + (cur_item_vs[0] != -1) as i32;
		let cur_ne = ne_st + has_cap as i32;
		self.stk.push(RelabelStack { cur_idx, ch_idx: ch_st, ch_end: ch_en, cur_nv, cur_ne });
	}

	fn start_child(&mut self) -> i32 {
		let si = self.stk.len() - 1;
		let RelabelStack { cur_idx, ch_idx, ch_end: ch_en, cur_nv, cur_ne } = self.stk[si];
		debug_assert!(ch_idx < ch_en);
		let nxt_item = self.ch.dat[ch_idx as usize];
		let nxt_idx = self.nxt_unassigned_idx;
		self.ch.dat[ch_idx as usize] = nxt_idx;
		self.par[nxt_idx as usize] = cur_idx;
		let nxt_ne = self.node_nes.bounds[nxt_idx as usize];
		if nxt_item < 1 + self.b.nv {
			self.vert_par_nv[nxt_idx as usize] = cur_nv;
			self.stk[si].cur_nv += 1;
		} else if self.types[cur_idx as usize] != NodeType::F && self.types[cur_idx as usize] != NodeType::V {
			self.node_edges[cur_ne as usize].twin_ne = nxt_ne;
			self.node_edges[nxt_ne as usize].twin_ne = cur_ne;
			self.stk[si].cur_ne += 1;
		}

		self.stk[si].ch_idx += 1;
		nxt_item
	}

	fn pop_item(&mut self) {
		let RelabelStack { cur_idx, ch_idx, ch_end: ch_en, .. } = self.stk.pop().unwrap();
		debug_assert!(ch_idx == ch_en);
		self.subtree_end[cur_idx as usize] = self.nxt_unassigned_idx;
	}
}

fn build_impl<const WP: bool>(nv: i32, edges: &[[i32; 2]], ternarize: bool, vert_order: &[i32], edge_order: &[i32]) -> PlanarSpqrTree {
	let ne = edges.len() as i32;
	assert!(vert_order.len() as i32 <= nv);
	assert!(edge_order.len() as i32 <= ne);

	let mut roots: Vec<i32> = Vec::with_capacity(nv as usize);
	let outedges: Csr<OutEdge>;

	// Phase 1: build a sorted skeleton
	{
		// 1a: build a normal adjacency list for the initial lowval dfs
		let mut adj_builder = CsrBuilder::<AdjEdge>::new(nv);
		for &[u, v] in edges {
			adj_builder.count(u);
			if u != v {
				adj_builder.count(v);
			}
		}
		adj_builder.allocate();
		for_each_in_order(ne, edge_order, |e| {
			let [u, v] = edges[e as usize];
			*adj_builder.push(u) = AdjEdge { dest: v, e };
			if u != v {
				*adj_builder.push(v) = AdjEdge { dest: u, e };
			}
		});
		let adj = adj_builder.finalize();

		let mut dfs = LowvalDfs {
			adj: &adj,
			depth: vec![-1; nv as usize],
			all_outedges: Vec::with_capacity(ne as usize),
			stk: Vec::with_capacity(nv as usize),
		};
		for_each_in_order(nv, vert_order, |rt| {
			if dfs.depth[rt as usize] == -1 {
				roots.push(rt);
				dfs.push_vert(rt, -1);
				loop {
					let s = dfs.stk.last().unwrap();
					if s.ch_idx == s.ch_end {
						let lowvals = dfs.pop_vert();
						if dfs.stk.is_empty() {
							break;
						}
						dfs.finish_edge(true, lowvals);
					} else {
						dfs.start_edge();
					}
				}
			}
		});
		let all_outedges = dfs.all_outedges;

		let mut by_key_builder = CsrBuilder::<OutEdge>::new(3 * nv + 6);
		for edge in &all_outedges {
			by_key_builder.count(edge.key);
		}
		by_key_builder.allocate();
		for &edge in &all_outedges {
			*by_key_builder.push(edge.key) = edge;
		}
		let by_key = by_key_builder.finalize();

		let mut by_src_builder = CsrBuilder::<OutEdge>::new(nv);
		// Hack to reuse memory
		by_src_builder.dat = all_outedges;
		for edge in &by_key.dat {
			by_src_builder.count(edge.src);
		}
		by_src_builder.allocate();
		for &edge in &by_key.dat {
			*by_src_builder.push(edge.src) = edge;
		}
		outedges = by_src_builder.finalize();
	}

	// Phase 2: do the big ear-decomposition-like walk
	let n_items0 = (1 + nv + ne) as usize;
	let mut b = Builder::<WP> {
		nv,
		ne,
		ternarize,
		outedges,
		ch_nxt: {
			let mut v = Vec::with_capacity(n_items0 + ne as usize);
			v.resize(n_items0, -1);
			v
		},
		item_vs: {
			let mut v = Vec::with_capacity(n_items0 + ne as usize);
			v.resize(n_items0, [-1, -1]);
			v
		},
		item_ch: {
			let mut v = Vec::with_capacity(n_items0 + ne as usize);
			v.resize(n_items0, ItemList::default());
			v
		},
		item_types: {
			let mut v = Vec::with_capacity(n_items0 + ne as usize);
			v.resize(1, NodeType::F);
			v.resize(1 + nv as usize, NodeType::V);
			v.resize(n_items0, NodeType::Q);
			v
		},
		quarter_edge_matches: vec![-1; if WP { (8 * ne + 4) as usize } else { 0 }],
		node_planarity: Vec::with_capacity(if WP { ne as usize } else { 0 }),
		tot_blocks: 0,
		tot_self_loops: 0,
		stack_verts: vec![0; nv as usize],
		stack_dir: vec![false; nv as usize],
		nxt_edge_idx: 0,
		first_occurrence: vec![0; nv as usize],
		edge_top_depths: vec![-1; if WP { (2 * ne) as usize } else { 0 }],
		tstack: Vec::with_capacity((nv + ne) as usize),
		stk: Vec::with_capacity(nv as usize),
	};

	for &rt in &roots {
		b.push_vert(rt);
		loop {
			let s = b.stk[b.stk.len() - 1];
			if s.ch_idx == s.ch_end {
				b.pop_vert();
				if b.stk.is_empty() {
					break;
				}
				b.finish_edge();
			} else if let Some(nxt) = b.start_edge() {
				b.push_vert(nxt);
			} else {
				b.finish_edge();
			}
		}
		let t = b.tstack.pop().unwrap();
		let l = b.concat(b.item_ch[ROOT_ITEM as usize], t.spans[1]);
		b.item_ch[ROOT_ITEM as usize] = l;
	}

	// Phase 3: relabel the full tree in preorder
	let tot_items = b.item_types.len() as i32;
	let tot_blocks = b.tot_blocks;
	let tot_self_loops = b.tot_self_loops;

	// Each node is a child, and additionally most non-block node has 2 cap verts; blocks have 1, and O nodes have 1
	let tot_node_verts = nv + (tot_items - 1 - nv) * 2 - tot_blocks - tot_self_loops;
	let tot_node_edges = (tot_items - 1 - nv - tot_blocks) * 2;

	let mut r = Relabel::<WP> {
		b,
		edges,
		vert_index: vec![-1; nv as usize],
		edge_index: vec![-1; ne as usize],
		edge_flipped: vec![false; ne as usize],
		par: vec![-1; tot_items as usize],
		subtree_end: vec![-1; tot_items as usize],
		types: vec![NodeType::F; tot_items as usize],
		orig_id: vec![-1; tot_items as usize],
		ch: Csr { bounds: vec![0; tot_items as usize + 1], dat: vec![0; (tot_items - 1) as usize] },
		node_verts: vec![NodeVert::default(); tot_node_verts as usize],
		node_nvs: CsrIndex { bounds: vec![0; tot_items as usize + 1] },
		vert_par_nv: vec![-1; tot_items as usize],
		node_edges: vec![NodeEdge::default(); tot_node_edges as usize],
		node_nes: CsrIndex { bounds: vec![0; tot_items as usize + 1] },
		node_adj: Csr { bounds: vec![0; (tot_node_verts * 2 + 1) as usize], dat: vec![NodeAdj::default(); (tot_node_edges * 2) as usize] },
		node_planar: vec![false; if WP { tot_items as usize } else { 0 }],
		ne_rot_adj: vec![-1; if WP { (4 * tot_node_edges) as usize } else { 0 }],
		vert_pos_buf: vec![-1; nv as usize],
		cnts_buf: vec![-1; (2 * nv) as usize],
		ch_buf: Vec::with_capacity(tot_items as usize),
		rot_edge_ne: vec![0; if WP { (2 * ne + 1) as usize } else { 0 }],
		nxt_unassigned_idx: 0,
		stk: Vec::with_capacity(tot_items as usize),
	};

	r.par[r.nxt_unassigned_idx as usize] = -1;
	r.push_item(ROOT_ITEM);
	loop {
		let s = r.stk[r.stk.len() - 1];
		if s.ch_idx == s.ch_end {
			r.pop_item();
			if r.stk.is_empty() {
				break;
			}
		} else {
			let nxt = r.start_child();
			r.push_item(nxt);
		}
	}

	debug_assert!(r.nxt_unassigned_idx == tot_items);
	debug_assert!(*r.ch.bounds.last().unwrap() == r.ch.dat.len() as i32);
	debug_assert!(*r.node_nvs.bounds.last().unwrap() == r.node_verts.len() as i32);
	debug_assert!(*r.node_nes.bounds.last().unwrap() == r.node_edges.len() as i32);
	debug_assert!(*r.node_adj.bounds.last().unwrap() == r.node_adj.dat.len() as i32);

	// Rewrite node_vertices to the correct index
	for v in r.node_verts.iter_mut() {
		v.vert = r.vert_index[v.vert as usize];
	}

	PlanarSpqrTree {
		tree: SpqrTree {
			vert_index: r.vert_index,
			edge_index: r.edge_index,
			edge_flipped: r.edge_flipped,
			par: r.par,
			subtree_end: r.subtree_end,
			types: r.types,
			orig_id: r.orig_id,
			ch: r.ch,
			node_verts: r.node_verts,
			node_nvs: r.node_nvs,
			vert_par_nv: r.vert_par_nv,
			node_edges: r.node_edges,
			node_nes: r.node_nes,
			node_adj: r.node_adj,
		},
		node_planar: r.node_planar,
		ne_embedding: PlanarEmbedding { rot_adj: r.ne_rot_adj },
	}
}
