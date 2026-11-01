//! SPQR tree — idiomatic Rust version of cp-book's `graph/spqr_tree.hpp`.
//!
//! Same algorithm and identical output as [`crate::spqr_tree`] (the line-by-line port), restructured the
//! way one would write it in Rust from scratch:
//!  * ids are `usize` locally and [`Idx`] (a `u32` with `u32::MAX` as niche) in stored arrays, so
//!    `Option<Idx>` is 4 bytes and `None` replaces every `-1` sentinel;
//!  * the C++ `[&]` lambdas become methods on small objects with disjoint data ([`Items`], [`TStacks`],
//!    [`Planarity`]) that a per-phase driver composes, instead of one struct with 30 fields;
//!  * linear passes are iterators ([`in_order`], [`Csr::bucket`]), and packed bit tricks become structs
//!    (`FlipItem`, `Span`, `Key`);
//!  * `with_planarity` is a const generic; the planarity data lives in its own object and is empty otherwise.

use std::fmt;
use std::iter;
use std::num::NonZeroU32;
use std::ops::{Deref, DerefMut, Index, IndexMut, Range};

// ---------------------------------------------------------------------------------------------
// Ids

/// A `u32` index that is never `u32::MAX`, so `Option<Idx>` is 4 bytes with `None` taking the slot the
/// C++ uses for `-1`. (Stored inverted in a `NonZeroU32`, which is what gives the niche.)
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Idx(NonZeroU32);

// Stored complemented, so ordering must go through `get`.
impl PartialOrd for Idx {
	fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
		Some(self.cmp(other))
	}
}
impl Ord for Idx {
	fn cmp(&self, other: &Self) -> std::cmp::Ordering {
		self.get().cmp(&other.get())
	}
}

impl Idx {
	pub const fn new(i: u32) -> Option<Idx> {
		match NonZeroU32::new(!i) {
			Some(n) => Some(Idx(n)),
			None => None,
		}
	}
	pub const fn get(self) -> u32 {
		!self.0.get()
	}
	pub const fn usize(self) -> usize {
		self.get() as usize
	}
}

impl From<usize> for Idx {
	fn from(i: usize) -> Idx {
		u32::try_from(i).ok().and_then(Idx::new).expect("index out of Idx range")
	}
}

impl fmt::Debug for Idx {
	fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
		write!(f, "{}", self.get())
	}
}
impl fmt::Display for Idx {
	fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
		write!(f, "{}", self.get())
	}
}

fn idx(i: usize) -> Idx {
	Idx::from(i)
}

fn u32_of(i: usize) -> u32 {
	u32::try_from(i).expect("index fits in u32")
}

// ---------------------------------------------------------------------------------------------
// CSR (jagged array)

/// The row bounds of a jagged array (n + 1 offsets) whose entries live in a separate flat `Vec`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct CsrIndex {
	pub bounds: Vec<usize>,
}

impl CsrIndex {
	pub fn len(&self) -> usize {
		self.bounds.len().saturating_sub(1)
	}
	pub fn is_empty(&self) -> bool {
		self.len() == 0
	}
	pub fn range(&self, i: usize) -> Range<usize> {
		self.bounds[i]..self.bounds[i + 1]
	}
	pub fn slice<'a, T>(&self, i: usize, base: &'a [T]) -> &'a [T] {
		&base[self.range(i)]
	}
}

/// A jagged array stored as `bounds` (n + 1 offsets) into `dat`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Csr<T> {
	pub bounds: Vec<usize>,
	pub dat: Vec<T>,
}

impl<T> Csr<T> {
	pub fn len(&self) -> usize {
		self.bounds.len() - 1
	}
	pub fn is_empty(&self) -> bool {
		self.len() == 0
	}
	pub fn range(&self, i: usize) -> Range<usize> {
		self.bounds[i]..self.bounds[i + 1]
	}
	pub fn rows(&self) -> impl Iterator<Item = &[T]> + '_ {
		self.bounds.windows(2).map(move |w| &self.dat[w[0]..w[1]])
	}
}

impl<T: Copy> Csr<T> {
	/// Counting sort: bucket `items` by `key` (in `0..n`), keeping the iteration order within each bucket.
	/// Bounds of bucket `i` as `u32`, for compact cursors.
	fn range_u32(&self, i: usize) -> Range<u32> {
		u32_of(self.bounds[i])..u32_of(self.bounds[i + 1])
	}

	pub fn bucket<I>(n: usize, items: I, key: impl Fn(&T) -> usize) -> Csr<T>
	where
		I: Iterator<Item = T> + Clone,
	{
		let mut bounds = vec![0; n + 1];
		for it in items.clone() {
			bounds[key(&it) + 1] += 1;
		}
		let total = prefix_sums(&mut bounds);
		// Fill with any item as a placeholder, then scatter with a cursor per bucket.
		let Some(first) = items.clone().next() else {
			return Csr { bounds, dat: Vec::new() };
		};
		let mut dat = vec![first; total];
		let mut cursor = bounds.clone();
		for it in items {
			let k = key(&it);
			dat[cursor[k]] = it;
			cursor[k] += 1;
		}
		Csr { bounds, dat }
	}
}

impl<T> Index<usize> for Csr<T> {
	type Output = [T];
	fn index(&self, i: usize) -> &[T] {
		&self.dat[self.range(i)]
	}
}
impl<T> IndexMut<usize> for Csr<T> {
	fn index_mut(&mut self, i: usize) -> &mut [T] {
		let r = self.range(i);
		&mut self.dat[r]
	}
}

// ---------------------------------------------------------------------------------------------
// Public output types (see `spqr_tree.rs` for the full description of the decomposition)

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum NodeType {
	/// forest root
	#[default]
	F,
	/// original vertex
	V,
	/// real edge: 1 vedge + 1 real edge
	Q,
	/// bridge
	I,
	/// self-loop
	O,
	/// series (cycle of >= 3 vedges)
	S,
	/// parallel (>= 3 vedges)
	P,
	/// rigid (3-connected)
	R,
}

impl NodeType {
	pub fn letter(self) -> char {
		match self {
			NodeType::F => 'F',
			NodeType::V => 'V',
			NodeType::Q => 'Q',
			NodeType::I => 'I',
			NodeType::O => 'O',
			NodeType::S => 'S',
			NodeType::P => 'P',
			NodeType::R => 'R',
		}
	}
	/// SPQR nodes proper (everything except the forest root and vertex items).
	pub fn is_node(self) -> bool {
		!matches!(self, NodeType::F | NodeType::V)
	}
}

impl fmt::Display for NodeType {
	fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
		write!(f, "{}", self.letter())
	}
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct NodeVert {
	pub node: Idx,
	pub vert: Idx,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct NodeEdge {
	pub node: Idx,
	pub twin_ne: Idx,
	pub nvs: [Idx; 2],
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct NodeAdj {
	pub ne: Idx,
	pub dest_nv: Idx,
}

/// The SPQR forest of a graph, embedded in its block-cut tree; items are in preorder.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SpqrTree {
	/// item of each original vertex / edge
	pub vert_index: Vec<Idx>,
	pub edge_index: Vec<Idx>,
	/// whether edge e is stored in its node with its endpoints swapped (`node_verts[..].vert` of its first nv != `edges[e][0]`)
	pub edge_flipped: Vec<bool>,

	pub par: Vec<Option<Idx>>,
	/// preorder end of each item's subtree
	pub subtree_end: Vec<Idx>,
	pub types: Vec<NodeType>,
	/// original vertex / edge id for V / Q items
	pub orig_id: Vec<Option<Idx>>,

	pub ch: Csr<Idx>,
	/// nv's of each node: `node_nvs.slice(i, &node_verts)`
	pub node_verts: Vec<NodeVert>,
	pub node_nvs: CsrIndex,
	/// The nv index of a vertex within its parent node
	pub vert_par_nv: Vec<Option<Idx>>,
	/// ne's of each node: `node_nes.slice(i, &node_edges)`
	pub node_edges: Vec<NodeEdge>,
	pub node_nes: CsrIndex,
	pub node_adj: Csr<NodeAdj>,
}

impl SpqrTree {
	pub fn len(&self) -> usize {
		self.par.len()
	}
	pub fn is_empty(&self) -> bool {
		self.par.is_empty()
	}

	/// `vert_order` and `edge_order` are (prefixes of) permutations of vertex / edge ids;
	/// listed ids are visited first in the given order, then the rest in id order.
	/// Roots are the first unvisited vertices, and DFS children are explored in edge order.
	pub fn build(nv: usize, edges: &[[usize; 2]], ternarize: bool, vert_order: &[usize], edge_order: &[usize]) -> SpqrTree {
		build::<false>(nv, edges, ternarize, vert_order, edge_order).tree
	}
}

/// Quarter-edges are indexed by `4 * edge + 2 * side + dir` (side: v0 vs v1, dir: cw vs ccw); `qe ^ 1` is the
/// other side around the endpoint, `qe ^ 3` the other side along the edge, and `rot_adj[qe]` the facing quarter-edge.
/// `None` entries mark non-embedded edges (a partial embedding).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PlanarEmbedding {
	pub rot_adj: Vec<Option<Idx>>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PlanarSpqrTree {
	pub tree: SpqrTree,
	pub node_planar: Vec<bool>,
	/// Embedding of each node's vedges, indexed by `4 * ne + 2 * side + dir`; `None` for nonplanar nodes.
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
	pub fn build(nv: usize, edges: &[[usize; 2]], ternarize: bool, vert_order: &[usize], edge_order: &[usize]) -> PlanarSpqrTree {
		build::<true>(nv, edges, ternarize, vert_order, edge_order)
	}
}

// ---------------------------------------------------------------------------------------------
// Helpers

/// Which of the remaining ids (those not in an explicit order prefix) still need visiting.
#[derive(Clone)]
enum Rest {
	All,
	Except(usize),
	Unlisted(Vec<bool>),
}

/// `order` first, then the remaining ids in `0..n` in increasing order.
fn in_order(n: usize, order: &[usize]) -> impl Iterator<Item = usize> + Clone + '_ {
	let rest = match order {
		_ if order.len() == n => Rest::Unlisted(vec![true; n]),
		[] => Rest::All,
		&[o] => Rest::Except(o),
		_ => {
			let mut listed = vec![false; n];
			for &i in order {
				listed[i] = true;
			}
			Rest::Unlisted(listed)
		}
	};
	let visits = move |i: &usize| match &rest {
		Rest::All => true,
		Rest::Except(o) => *i != *o,
		Rest::Unlisted(listed) => !listed[*i],
	};
	order.iter().copied().chain((0..n).filter(visits))
}

/// `[T; 2]` indexed by a side: `arr[dir] == a`, `arr[!dir] == b`. (Compiles to cmovs.)
fn by_side<T>(dir: bool, a: T, b: T) -> [T; 2] {
	if dir { [b, a] } else { [a, b] }
}
fn side<T: Copy>(arr: [T; 2], dir: bool) -> T {
	if dir { arr[1] } else { arr[0] }
}

// ---------------------------------------------------------------------------------------------
// Phase 1: lowpoint DFS, producing the out-edges of the DFS tree sorted by (lowpoint class, kind)

#[derive(Clone, Copy)]
struct AdjEdge {
	dest: usize,
	e: usize,
}

/// Sort key of a DFS out-edge as seen from its source at `cur_depth`. Encodes, in this order:
/// first the two "local" classes (bridges / whole components: lowval == cur_depth + 1, then self loops and
/// non-tree parents: lowval == cur_depth), then real lowvals increasing; ties broken by kind.
#[derive(Clone, Copy)]
struct Key(u32);

#[derive(Clone, Copy)]
struct EdgeKind {
	lowval: usize,
	is_tree: bool,
	/// type-1 = the subtree's second lowval is also above cur (single return), or a back edge
	is_type_1: bool,
}

impl Key {
	fn encode(cur_depth: usize, child_lowvals: [usize; 2], is_tree: bool) -> Key {
		let d = cur_depth as i64;
		let mut lowval = child_lowvals[0] as i64;
		if lowval >= d {
			// Bridges have lowval -2, and loops/components -1
			lowval = !(lowval - d);
		}
		let kind = 2 * i64::from(child_lowvals[1] < cur_depth) + i64::from(!is_tree);
		Key(u32::try_from(3 * (lowval + 2) + kind).expect("key fits"))
	}
	fn decode(self, cur_depth: usize) -> EdgeKind {
		let mut lowval = i64::from(self.0 / 3) - 2;
		if lowval < 0 {
			lowval = cur_depth as i64 + !lowval;
		}
		let kind = self.0 % 3;
		EdgeKind { lowval: lowval as usize, is_tree: kind != 1, is_type_1: kind <= 1 }
	}
	fn range(nv: usize) -> usize {
		3 * nv + 6
	}
}

#[derive(Clone, Copy)]
struct OutEdge {
	src: Idx,
	dest: Idx,
	e: Idx,
	key: Key,
}

struct LowvalFrame {
	cur: Idx,
	prv_e: Option<Idx>,
	/// The two smallest distinct depths reachable from the subtree so far
	lowvals: [Idx; 2],
	/// Unvisited part of `adj[cur]`
	edges: Range<u32>,
}

struct LowvalDfs<'a> {
	adj: &'a Csr<AdjEdge>,
	depth: Vec<Option<Idx>>,
	outedges: Vec<OutEdge>,
	stk: Vec<LowvalFrame>,
}

impl<'a> LowvalDfs<'a> {
	fn new(nv: usize, ne: usize, adj: &'a Csr<AdjEdge>) -> Self {
		LowvalDfs { adj, depth: vec![None; nv], outedges: Vec::with_capacity(ne), stk: Vec::with_capacity(nv) }
	}

	fn visited(&self, v: usize) -> bool {
		self.depth[v].is_some()
	}

	fn push_vert(&mut self, cur: usize, prv_e: Option<usize>) {
		let d = self.stk.len();
		self.depth[cur] = Some(idx(d));
		self.stk.push(LowvalFrame { cur: idx(cur), prv_e: prv_e.map(idx), lowvals: [idx(d); 2], edges: self.adj.range_u32(cur) });
	}

	/// Account for the edge at the top frame's cursor leading to a subtree (or back edge) with `child_lowvals`.
	fn finish_edge(&mut self, is_tree: bool, child_lowvals: [usize; 2]) {
		let d = self.stk.len() - 1;
		let s = self.stk.last_mut().expect("in a dfs");
		let AdjEdge { dest, e } = self.adj.dat[s.edges.next().expect("edge pending") as usize];
		self.outedges.push(OutEdge { src: s.cur, dest: idx(dest), e: idx(e), key: Key::encode(d, child_lowvals, is_tree) });

		// Keep the 2 smallest distinct values
		let [l0, l1] = s.lowvals;
		let [c0, c1] = child_lowvals.map(idx);
		s.lowvals = if c0 < l0 {
			[c0, c1.min(l0)]
		} else if c0 == l0 {
			[l0, l1.min(c1)]
		} else {
			[l0, l1.min(c0)]
		};
	}

	/// Look at the next edge of the top frame: skip it, descend (returns `Some(child)`), or fold in a back edge.
	fn start_edge(&mut self) -> Option<usize> {
		let d = self.stk.len() - 1;
		let s = self.stk.last_mut().expect("in a dfs");
		let AdjEdge { dest, e } = self.adj.dat[s.edges.start as usize];
		if Some(idx(e)) == s.prv_e || self.depth[dest].is_some_and(|dd| dd.usize() > d) {
			// parent edge, or already handled from the other side
			s.edges.next();
			return None;
		}
		match self.depth[dest] {
			None => {
				self.push_vert(dest, Some(e));
				Some(dest)
			}
			Some(dd) => {
				self.finish_edge(false, [dd.usize(), d]);
				None
			}
		}
	}

	fn run(&mut self, root: usize) {
		self.push_vert(root, None);
		loop {
			let s = self.stk.last().expect("in a dfs");
			if s.edges.is_empty() {
				let lowvals = self.stk.pop().expect("in a dfs").lowvals.map(Idx::usize);
				if self.stk.is_empty() {
					break;
				}
				self.finish_edge(true, lowvals);
			} else {
				self.start_edge();
			}
		}
	}
}

/// Every DFS out-edge (tree and back edges), grouped by source, each group sorted by [`Key`].
/// Turn per-bucket counts (stored at `counts[k + 1]`) into CSR bounds; returns the total.
fn prefix_sums(counts: &mut [usize]) -> usize {
	let mut acc = 0;
	for b in &mut counts[1..] {
		acc += *b;
		*b = acc;
	}
	acc
}

/// Adjacency lists; `edge_order` decides the order within each list, self-loops appear once.
fn adjacency(nv: usize, edges: &[[usize; 2]], edge_order: &[usize]) -> Csr<AdjEdge> {
	let mut bounds = vec![0; nv + 1];
	for &[u, v] in edges {
		bounds[u + 1] += 1;
		if u != v {
			bounds[v + 1] += 1;
		}
	}
	let total = prefix_sums(&mut bounds);
	let mut dat = vec![AdjEdge { dest: 0, e: 0 }; total];
	let mut cursor = bounds.clone();
	for e in in_order(edges.len(), edge_order) {
		let [u, v] = edges[e];
		dat[cursor[u]] = AdjEdge { dest: v, e };
		cursor[u] += 1;
		if u != v {
			dat[cursor[v]] = AdjEdge { dest: u, e };
			cursor[v] += 1;
		}
	}
	Csr { bounds, dat }
}

fn sorted_skeleton(nv: usize, edges: &[[usize; 2]], vert_order: &[usize], edge_order: &[usize]) -> (Vec<usize>, Csr<OutEdge>) {
	let ne = edges.len();
	let adj = adjacency(nv, edges, edge_order);

	let mut dfs = LowvalDfs::new(nv, ne, &adj);
	let mut roots = Vec::new();
	for rt in in_order(nv, vert_order) {
		if !dfs.visited(rt) {
			roots.push(rt);
			dfs.run(rt);
		}
	}
	debug_assert!(dfs.outedges.len() == ne);
	debug_assert!((0..nv).all(|v| dfs.visited(v)));

	let by_key = Csr::bucket(Key::range(nv), dfs.outedges.iter().copied(), |o| o.key.0 as usize);
	let by_src = Csr::bucket(nv, by_key.dat.iter().copied(), |o| o.src.usize());
	(roots, by_src)
}

// ---------------------------------------------------------------------------------------------
// Phase 2: the ear-decomposition-like walk that discovers the SPQR nodes
//
// We build a tree of all SPQR *nodes* + all original *vertices* (collectively *items*). Vertices hang off the
// first SPQR node containing them, and blocks are rooted at a topmost Q node for the top edge.
// Children are kept as linked lists (`Items::ch_nxt`) until phase 3 lays them out.

const ROOT_ITEM: usize = 0;

/// A link in a child list. `flip` toggles the planar orientation of everything after it (always false
/// without planarity).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct FlipItem {
	item: Idx,
	flip: bool,
}

impl FlipItem {
	fn new(item: usize) -> Self {
		FlipItem { item: idx(item), flip: false }
	}
	fn flipped(self) -> Self {
		FlipItem { flip: !self.flip, ..self }
	}
}

/// A non-empty child list, by its two ends.
#[derive(Clone, Copy, Debug)]
struct Span {
	first: FlipItem,
	last: FlipItem,
}

type ItemList = Option<Span>;

fn unit(item: usize) -> ItemList {
	let f = FlipItem::new(item);
	Some(Span { first: f, last: f })
}

fn flipped(l: ItemList) -> ItemList {
	l.map(|Span { first, last }| Span { first: first.flipped(), last: last.flipped() })
}

/// The item tree under construction.
struct Items {
	nv: usize,
	ne: usize,
	/// The endpoints (original vertices) of each node: two for regular nodes, one for blocks / loops.
	vs: Vec<[Option<Idx>; 2]>,
	ch: Vec<ItemList>,
	types: Vec<NodeType>,
	ch_nxt: Vec<Option<FlipItem>>,
}

impl Items {
	fn new(nv: usize, ne: usize) -> Items {
		let n0 = 1 + nv + ne;
		let cap = n0 + ne;
		let mut types = Vec::with_capacity(cap);
		types.push(NodeType::F);
		types.extend(iter::repeat_n(NodeType::V, nv));
		types.extend(iter::repeat_n(NodeType::Q, ne));
		let mut vs = Vec::with_capacity(cap);
		vs.resize(n0, [None; 2]);
		let mut ch = Vec::with_capacity(cap);
		ch.resize(n0, None);
		let mut ch_nxt = Vec::with_capacity(cap);
		ch_nxt.resize(n0, None);
		Items { nv, ne, vs, ch, types, ch_nxt }
	}

	fn len(&self) -> usize {
		self.types.len()
	}
	fn vert(&self, v: usize) -> usize {
		1 + v
	}
	fn edge(&self, e: usize) -> usize {
		1 + self.nv + e
	}
	fn is_vert(&self, item: usize) -> bool {
		(1..=self.nv).contains(&item)
	}
	/// Vedges are identified with the item they cap: Q items and nodes.
	fn vedge(&self, item: usize) -> usize {
		item - (1 + self.nv)
	}
	/// Index among allocated (non-Q) nodes.
	fn node(&self, item: usize) -> usize {
		item - (1 + self.nv + self.ne)
	}

	fn alloc(&mut self, ty: NodeType) -> usize {
		self.vs.push([None; 2]);
		self.ch.push(None);
		self.ch_nxt.push(None);
		self.types.push(ty);
		self.types.len() - 1
	}

	fn concat(&mut self, a: ItemList, b: ItemList) -> ItemList {
		match (a, b) {
			(a, None) => a,
			(None, b) => b,
			(Some(a), Some(b)) => {
				self.ch_nxt[a.last.item.usize()] = Some(FlipItem { item: b.first.item, flip: b.first.flip ^ a.last.flip });
				Some(Span { first: a.first, last: b.last })
			}
		}
	}

	fn append_children(&mut self, item: usize, l: ItemList) {
		let cur = self.ch[item];
		self.ch[item] = self.concat(cur, l);
	}
}

/// Planar-embedding bookkeeping (empty without planarity).
/// Each vedge `ve` has 4 quarter-edges `4 * ve + 2 * source_vert + is_cw`; `matches` pairs facing ones.
#[derive(Default)]
struct Planarity {
	matches: Vec<Option<Idx>>,
	edge_top_depths: Vec<Option<Idx>>,
	node_planarity: Vec<Result<[Idx; 4], Nonplanar>>,
}

impl Planarity {
	fn new(nv: usize, ne: usize) -> Self {
		let _ = nv;
		Planarity { matches: vec![None; 8 * ne + 4], edge_top_depths: vec![None; 2 * ne], node_planarity: Vec::with_capacity(ne) }
	}

	fn quarter(ve: usize, side: usize, dir: usize) -> Idx {
		idx(4 * ve + 2 * side + dir)
	}
	fn link(&mut self, a: Idx, b: Idx) {
		self.matches[a.usize()] = Some(b);
		self.matches[b.usize()] = Some(a);
	}
}

/// TODO: What's the nonplanarity certificate look like?
#[derive(Clone, Copy, Debug, Default)]
struct Nonplanar;

#[derive(Clone, Copy, Debug, Default)]
struct PlanaritySide {
	// For each side, pointers to the "linked lists" of the edges inside.
	// [0] is the outer / longer edges and [1] the inner / shorter edges, matching the outside-in sort order.

	// exposed pieces of the walk down the ear in the tree (connected to the bottommost/topmost vertices of the tree path)
	bot_ends: [Option<Idx>; 2],
	// exposed backedges
	top_ends: [Option<Idx>; 2],
	// depths should be increasing going inwards
	top_depths: [Option<Idx>; 2],
}

/// Convention: `sides[0].top_depths[0] == top_depth`, i.e. at least one minimal return lives on side 0.
#[derive(Clone, Copy, Debug, Default)]
struct TstackPlanarity {
	sides: [PlanaritySide; 2],
}

type MaybePlanarity = Result<TstackPlanarity, Nonplanar>;

fn merge_planarity(pl: &mut Planarity, a: &mut MaybePlanarity, b: &MaybePlanarity) {
	let Ok(ap) = a else { return };
	let Ok(bp) = b else {
		*a = Err(Nonplanar);
		return;
	};
	for (as_, bs) in ap.sides.iter_mut().zip(&bp.sides) {
		// No bottom edges means an isolated vertex: nothing to merge.
		let (Some(a_bot1), Some(b_bot0)) = (as_.bot_ends[1], bs.bot_ends[0]) else {
			if as_.bot_ends[0].is_none() {
				*as_ = *bs;
			}
			continue;
		};
		pl.link(a_bot1, b_bot0);
		as_.bot_ends[1] = bs.bot_ends[1];

		let (Some(a_top1), Some(b_top0)) = (as_.top_ends[1], bs.top_ends[0]) else {
			if as_.top_ends[0].is_none() {
				as_.top_ends = bs.top_ends;
				as_.top_depths = bs.top_depths;
			}
			continue;
		};
		if as_.top_depths[1] > bs.top_depths[0] {
			// TODO: Certificate
			*a = Err(Nonplanar);
			return;
		}
		pl.link(a_top1, b_top0);
		as_.top_ends[1] = bs.top_ends[1];
		as_.top_depths[1] = bs.top_depths[1];
	}
}

/// An open ear / partial component on the triconnectivity stack.
#[derive(Clone, Copy, Debug)]
struct Tstack {
	v_start: usize,
	top_depth: usize,
	/// Back-edge counter at creation, to tell which tstacks were pushed after a given back edge
	first_idx: usize,
	/// Children on the two sides of the ear
	spans: [ItemList; 2],
}

/// The triconnectivity stack, innermost on top. `planarity` shadows `stack` element-wise when `WP`.
struct TStacks<const WP: bool> {
	stack: Vec<Tstack>,
	planarity: Vec<MaybePlanarity>,
}

impl<const WP: bool> TStacks<WP> {
	fn with_capacity(n: usize) -> Self {
		TStacks { stack: Vec::with_capacity(n), planarity: Vec::with_capacity(if WP { n } else { 0 }) }
	}
	fn len(&self) -> usize {
		self.stack.len()
	}
	fn cur(&self) -> &Tstack {
		self.stack.last().expect("nonempty tstack")
	}
	fn cur_mut(&mut self) -> &mut Tstack {
		self.stack.last_mut().expect("nonempty tstack")
	}
	fn nxt(&self) -> &Tstack {
		&self.stack[self.stack.len() - 2]
	}
	fn push(&mut self, t: Tstack, p: MaybePlanarity) {
		self.stack.push(t);
		if WP {
			self.planarity.push(p);
		}
	}
	fn pop(&mut self) -> (Tstack, MaybePlanarity) {
		let t = self.stack.pop().expect("nonempty tstack");
		let p = if WP { self.planarity.pop().expect("shadow stack in sync") } else { Ok(TstackPlanarity::default()) };
		(t, p)
	}
	fn planarity_mut(&mut self, i: usize) -> &mut MaybePlanarity {
		debug_assert!(WP);
		&mut self.planarity[i]
	}
	fn cur_planarity_mut(&mut self) -> &mut MaybePlanarity {
		let i = self.len() - 1;
		self.planarity_mut(i)
	}

	/// Mirror tstack `i`: swap which side is which.
	fn flip(&mut self, i: usize) {
		if WP {
			let t = &mut self.stack[i];
			t.spans = t.spans.map(flipped);
			if let Ok(p) = &mut self.planarity[i] {
				p.sides.swap(0, 1);
			}
		}
	}

	/// Merge the top (inner) tstack into the one below it.
	#[inline]
	fn merge_tops(&mut self, items: &mut Items, pl: &mut Planarity) {
		let (b, bp) = self.pop();
		let a = self.cur_mut();
		a.top_depth = a.top_depth.min(b.top_depth);
		let [a0, a1] = a.spans;
		let spans = [items.concat(b.spans[0], a0), items.concat(a1, b.spans[1])];
		self.cur_mut().spans = spans;
		if WP {
			merge_planarity(pl, self.cur_planarity_mut(), &bp);
		}
	}
}

struct WalkFrame {
	has_vert_tstack: bool,
	/// Unvisited part of `outedges[cur]`
	edges: Range<u32>,
	/// tstack height when the current edge was started
	orig_tstack: u32,
}

struct Builder<const WP: bool> {
	ternarize: bool,
	outedges: Csr<OutEdge>,

	items: Items,
	tstacks: TStacks<WP>,
	planarity: Planarity,

	tot_blocks: usize,
	tot_self_loops: usize,

	// Per depth of the current DFS path. Most of the code is in terms of v_start / top_depth, and reads these out.
	stack_verts: Vec<usize>,
	stack_dir: Vec<bool>,
	/// First back edge (by `nxt_edge_idx`) to each depth
	first_occurrence: Vec<usize>,
	/// Counts back edges only
	nxt_edge_idx: usize,

	stk: Vec<WalkFrame>,
}

impl<const WP: bool> Builder<WP> {
	fn new(nv: usize, ne: usize, ternarize: bool, outedges: Csr<OutEdge>) -> Self {
		Builder {
			ternarize,
			outedges,
			items: Items::new(nv, ne),
			tstacks: TStacks::with_capacity(nv + ne),
			planarity: if WP { Planarity::new(nv, ne) } else { Planarity::default() },
			tot_blocks: 0,
			tot_self_loops: 0,
			stack_verts: vec![0; nv],
			stack_dir: vec![false; nv],
			first_occurrence: vec![0; nv],
			nxt_edge_idx: 0,
			stk: Vec::with_capacity(nv),
		}
	}

	fn alloc(&mut self, ty: NodeType) -> usize {
		let item = self.items.alloc(ty);
		if WP {
			// Set by finish_tstack_top before it is ever read
			self.planarity.node_planarity.push(Err(Nonplanar));
		}
		item
	}

	fn make_vs(&self, v_start: usize, top_depth: usize) -> [Option<Idx>; 2] {
		by_side(self.stack_dir[top_depth], Some(idx(self.stack_verts[top_depth])), Some(idx(v_start)))
	}

	fn make_edge_planarity(&mut self, item: usize, top_depth: usize, is_tree: bool) -> MaybePlanarity {
		let mut p = TstackPlanarity::default();
		if WP {
			let ve = self.items.vedge(item);
			let top_dir = usize::from(self.stack_dir[top_depth]);
			self.planarity.edge_top_depths[ve] = Some(idx(top_depth));
			let q = |s: usize, d: usize| Some(Planarity::quarter(ve, s, d));
			if is_tree {
				p.sides[0].bot_ends = [q(1 - top_dir, 0), q(top_dir, 1)];
				p.sides[1].bot_ends = [q(1 - top_dir, 1), q(top_dir, 0)];
			} else {
				p.sides[0].bot_ends = [q(1 - top_dir, 0), q(1 - top_dir, 1)];
				p.sides[0].top_ends = [q(top_dir, 1), q(top_dir, 0)];
				p.sides[0].top_depths = [Some(idx(top_depth)); 2];
			}
		}
		Ok(p)
	}

	#[inline]
	fn push_tstack(&mut self, v_start: usize, top_depth: usize, item: usize, planarity: MaybePlanarity) {
		let spans = by_side(self.stack_dir[top_depth], unit(item), None);
		self.tstacks.push(Tstack { v_start, top_depth, first_idx: self.nxt_edge_idx, spans }, planarity);
	}
	fn push_vert_tstack(&mut self, v: usize, top_depth: usize) {
		let item = self.items.vert(v);
		self.push_tstack(v, top_depth, item, Ok(TstackPlanarity::default()));
	}
	fn push_edge_tstack(&mut self, v_start: usize, top_depth: usize, e: usize, is_tree: bool) {
		let item = self.items.edge(e);
		let planarity = self.make_edge_planarity(item, top_depth, is_tree);
		self.push_tstack(v_start, top_depth, item, planarity);
	}
	fn merge_tstack_tops(&mut self) {
		self.tstacks.merge_tops(&mut self.items, &mut self.planarity);
	}

	/// The node that will absorb the second-to-top tstack: a fresh one, or (unless ternarizing) the S/P node
	/// that tstack consists of, re-opened.
	fn maybe_unwrap_nxt(&mut self, ty: NodeType, is_tree: bool) -> usize {
		if ty == NodeType::R || self.ternarize {
			return self.alloc(ty);
		}
		debug_assert!(matches!(ty, NodeType::P | NodeType::S));

		let ti = self.tstacks.len() - 2;
		let t = self.tstacks.stack[ti];
		let top_dir = self.stack_dir[t.top_depth];
		debug_assert!(side(t.spans, !top_dir).is_none());
		let Span { first, last } = side(t.spans, top_dir).expect("tstack has a child");
		let item = first.item.usize();
		debug_assert!(item == last.item.usize());
		if self.items.types[item] != ty {
			return self.alloc(ty);
		}

		self.tstacks.stack[ti].spans = by_side(top_dir, self.items.ch[item], None);
		if WP {
			// Unwrap the planarity data. S/P nodes are trivially planar, so the current state is just
			// make_edge_planarity(wrapped): right shape, just needs relabelling.
			let m = self.planarity.node_planarity[self.items.node(item)].expect("unwrapped S/P nodes are always planar");
			let p = self.tstacks.planarity_mut(ti).as_mut().expect("unwrapped tstack is a single edge");
			let td = usize::from(top_dir);
			let m = |i: usize| Some(m[i]);
			if is_tree {
				p.sides[0].bot_ends = [m(2 * (1 - td) + 1), m(2 * td)];
				p.sides[1].bot_ends = [m(2 * (1 - td)), m(2 * td + 1)];
			} else {
				p.sides[0].bot_ends = [m(2 * (1 - td) + 1), m(2 * (1 - td))];
				p.sides[0].top_ends = [m(2 * td), m(2 * td + 1)];
			}
		}
		item
	}

	/// Close the top tstack into `item`, leaving `item` on the tstack as a single edge.
	fn finish_tstack_top(&mut self, item: usize, is_tree: bool) {
		let ti = self.tstacks.len() - 1;
		let t = self.tstacks.stack[ti];
		let top_dir = self.stack_dir[t.top_depth];
		debug_assert!(side(t.spans, !top_dir).is_none());

		if WP {
			let np = self.items.node(item);
			self.planarity.node_planarity[np] = match &self.tstacks.planarity[ti] {
				Ok(p) => {
					let td = usize::from(top_dir);
					let [s0, s1] = &p.sides;
					let mut m = [None; 4];
					if is_tree {
						m[2 * (1 - td) + 1] = s0.bot_ends[0];
						m[2 * td] = s0.bot_ends[1];
						m[2 * (1 - td)] = s1.bot_ends[0];
						m[2 * td + 1] = s1.bot_ends[1];
					} else {
						m[2 * (1 - td) + 1] = s0.bot_ends[0];
						m[2 * (1 - td)] = s0.bot_ends[1];
						m[2 * td] = s0.top_ends[0];
						m[2 * td + 1] = s0.top_ends[1];
					}
					Ok(m.map(|q| q.expect("all quarter edges of a finished node are exposed")))
				}
				Err(Nonplanar) => {
					debug_assert!(self.items.types[item] == NodeType::R);
					Err(Nonplanar)
				}
			};
		}
		self.items.vs[item] = self.make_vs(t.v_start, t.top_depth);
		self.items.ch[item] = side(t.spans, top_dir);

		self.tstacks.stack[ti].spans = by_side(top_dir, unit(item), None);
		if WP {
			*self.tstacks.planarity_mut(ti) = self.make_edge_planarity(item, t.top_depth, is_tree);
		}
	}

	fn push_vert(&mut self, cur: usize) {
		self.stk.push(WalkFrame { has_vert_tstack: false, edges: self.outedges.range_u32(cur), orig_tstack: 0 });
		let cur_depth = self.stk.len() - 1;
		self.stack_verts[cur_depth] = cur;
	}

	/// Set up for the top frame's next edge; `Some(child)` means descend into it, otherwise finish it directly.
	fn start_edge(&mut self) -> Option<usize> {
		let cur_depth = self.stk.len() - 1;
		let cur = self.stack_verts[cur_depth];
		let s = self.stk.last_mut().expect("in a walk");
		let OutEdge { dest, key, .. } = self.outedges.dat[s.edges.start as usize];
		let dest = dest.usize();
		let EdgeKind { lowval, is_tree, is_type_1 } = key.decode(cur_depth);

		// edge_dir convention: false is forwards, true is backwards,
		// i.e. cur is on the edge_dir side and nxt on the !edge_dir side.
		self.stack_dir[cur_depth] = lowval < cur_depth && !self.stack_dir[lowval];

		if !s.has_vert_tstack && lowval < cur_depth && is_type_1 {
			// Do this with the correct stack_dir set
			s.has_vert_tstack = true;
			self.push_vert_tstack(cur, cur_depth);
		}

		let orig_tstack = self.tstacks.len();
		let s = self.stk.last_mut().expect("in a walk");
		s.orig_tstack = u32_of(orig_tstack);
		if is_tree {
			self.first_occurrence[cur_depth] = self.items.ne;
			Some(dest)
		} else {
			None
		}
	}

	fn finish_edge(&mut self) {
		let cur_depth = self.stk.len() - 1;
		let cur = self.stack_verts[cur_depth];
		let s = self.stk.last_mut().expect("in a walk");
		let OutEdge { dest: nxt, e, key, .. } = self.outedges.dat[s.edges.next().expect("edge pending") as usize];
		let (nxt, e) = (nxt.usize(), e.usize());
		let EdgeKind { lowval, is_tree, is_type_1 } = key.decode(cur_depth);
		let orig_tstack = s.orig_tstack as usize;
		let has_vert_tstack = s.has_vert_tstack;
		let edge_dir = self.stack_dir[cur_depth];
		let ei = self.items.edge(e);

		if lowval >= cur_depth {
			// A whole block: no planarity handling needed, it's just a Q node with an I/O node.
			self.items.vs[ei] = [Some(idx(cur)), None];
			self.tot_blocks += 1;
			let children = if is_tree {
				if lowval == cur_depth + 1 {
					// Bridge: the top tstack is just smuggling out the child vertex; prepend the bridge component.
					// (A shortcut for allocating a full I-type tstack.)
					let item = self.alloc(NodeType::I);
					self.items.vs[item] = self.make_vs(nxt, cur_depth);
					let (t, _) = self.tstacks.pop();
					self.items.concat(unit(item), t.spans[1])
				} else {
					// Component: the top tstack is the backedge and the one below it the vertex
					let (backedge, _) = self.tstacks.pop();
					let (t, _) = self.tstacks.pop();
					self.items.concat(backedge.spans[0], t.spans[1])
				}
			} else {
				debug_assert!(nxt == cur);
				self.tot_self_loops += 1;
				let item = self.alloc(NodeType::O);
				self.items.vs[item] = [Some(idx(cur)), None];
				unit(item)
			};
			self.items.ch[ei] = children;
			let vi = self.items.vert(cur);
			self.items.append_children(vi, unit(ei));
			return;
		}

		self.items.vs[ei] = self.make_vs(nxt, cur_depth);

		// Whether the top tstack is a single edge
		let mut is_single = true;
		if is_tree {
			// The span lives on side edge_dir
			self.push_edge_tstack(nxt, cur_depth, e, true);
			while self.tstacks.len() >= 2 && self.tstacks.nxt().top_depth >= cur_depth {
				let ty = if self.tstacks.nxt().top_depth > cur_depth {
					// Backfill this for maybe_unwrap
					self.stack_dir[self.tstacks.nxt().top_depth] = edge_dir;
					// The tstack currently holds a tree edge followed by a vertex; merge the vertex first
					self.merge_tstack_tops();
					NodeType::S
				} else if self.tstacks.nxt().v_start == self.tstacks.cur().v_start {
					NodeType::P
				} else {
					NodeType::R
				};
				let item = self.maybe_unwrap_nxt(ty, ty == NodeType::S);
				self.merge_tstack_tops();
				if WP {
					if let Ok(p) = self.tstacks.cur_planarity_mut() {
						// Merge all backedges into the component
						for side in &mut p.sides {
							let bot1 = side.bot_ends[1].expect("open ear");
							let Some(top1) = side.top_ends[1] else { continue };
							debug_assert!(side.top_depths == [Some(idx(cur_depth)); 2]);
							self.planarity.link(bot1, top1);
							side.bot_ends[1] = side.top_ends[0];
							side.top_depths = [None; 2];
							side.top_ends = [None; 2];
						}
					}
				}
				self.finish_tstack_top(item, true);
			}

			if self.tstacks.cur().first_idx > self.first_occurrence[cur_depth] {
				while self.tstacks.cur().first_idx > self.first_occurrence[cur_depth] {
					if WP {
						let n = self.tstacks.len();
						let (below, top) = (&self.tstacks.stack[n - 2], &self.tstacks.stack[n - 1]);
						if below.first_idx > self.first_occurrence[cur_depth] {
							// We will put cur_depth on side 1 until the bottom
							if below.top_depth == cur_depth {
								self.tstacks.flip(n - 2);
							}
						} else if !is_single {
							debug_assert!(top.top_depth < cur_depth);
							if let Ok(p) = &self.tstacks.planarity[n - 2] {
								if p.sides[0].top_depths[1] == Some(idx(cur_depth)) {
									// Flip cur and nxt relative to each other: the one with worse top_depth.
									let which = if top.top_depth < below.top_depth { n - 2 } else { n - 1 };
									self.tstacks.flip(which);
								} else {
									debug_assert!(p.sides[1].top_depths[1] == Some(idx(cur_depth)));
								}
							}
						}
					}
					self.merge_tstack_tops();
					is_single = false;
				}
				if WP {
					if let Ok(p) = self.tstacks.cur_planarity_mut() {
						// Prune off finished cur-side things
						for side in &mut p.sides {
							debug_assert!(side.bot_ends[1].is_some());
							while side.top_depths[1] == Some(idx(cur_depth)) {
								let (bot1, top1) = (side.bot_ends[1].expect("open ear"), side.top_ends[1].expect("has a top"));
								self.planarity.link(bot1, top1);
								// The other quarter edge of the same (vedge, side) is now the exposed bottom.
								let new_bot1 = idx(top1.usize() ^ 1);
								side.bot_ends[1] = Some(new_bot1);
								side.top_ends[1] = self.planarity.matches[new_bot1.usize()].take();
								if let Some(t1) = side.top_ends[1] {
									self.planarity.matches[t1.usize()] = None;
									side.top_depths[1] = self.planarity.edge_top_depths[t1.usize() >> 2];
								} else {
									side.top_depths = [None; 2];
									side.top_ends = [None; 2];
								}
							}
						}
					}
				}
			}

			if is_type_1 {
				debug_assert!(has_vert_tstack);
			}
			if has_vert_tstack {
				// NB: tstack[orig_tstack] is the vertex and tstack[orig_tstack + 1] is the backedge.
				debug_assert!(self.tstacks.len() >= orig_tstack + 3);

				if !is_type_1 {
					if WP {
						// The lowval side should be side 1, everything else goes on side 0.
						// The exception is tstack[orig_tstack + 2], which could be == lowval on one/both sides,
						// but is guaranteed to have *something* > lowval by non-type-1-ness.
						let ti = orig_tstack + 2;
						if let Ok(p) = self.tstacks.planarity[ti] {
							debug_assert!(p.sides[0].top_depths[0] == Some(idx(self.tstacks.stack[ti].top_depth)));
							if p.sides[0].top_depths[1] == Some(idx(lowval)) {
								self.tstacks.flip(ti);
							}
							let p = self.tstacks.planarity[ti].expect("unchanged by flip");
							debug_assert!(p.sides[0].top_depths[1].is_some_and(|d| d.usize() > lowval));
						}
						for i in orig_tstack + 3..self.tstacks.len() {
							if self.tstacks.stack[i].top_depth == lowval {
								self.tstacks.flip(i);
							}
						}
					}
					while self.tstacks.len() > orig_tstack + 3 {
						self.merge_tstack_tops();
						is_single = false;
					}
					debug_assert!(!is_single);
				}

				debug_assert!(self.tstacks.len() == orig_tstack + 3);
				let item = is_type_1.then(|| self.maybe_unwrap_nxt(if is_single { NodeType::S } else { NodeType::R }, false));
				// Merge with the backedge, then with the vertex
				self.merge_tstack_tops();
				self.merge_tstack_tops();

				let t = self.tstacks.cur_mut();
				t.v_start = cur;
				debug_assert!(t.top_depth == lowval);

				// Fold everything to the correct side now that we're leaving the child:
				// the entire subtree goes to the !edge_dir side.
				let [s0, s1] = t.spans;
				let all = self.items.concat(s0, s1);
				self.tstacks.cur_mut().spans = by_side(!edge_dir, all, None);

				if WP {
					if let Ok(p) = self.tstacks.cur_planarity_mut() {
						// precondition: side 1 is the lowval-only side
						let [s0, s1] = &mut p.sides;
						self.planarity.link(s0.bot_ends[0].expect("open ear"), s1.bot_ends[0].expect("open ear"));
						s0.bot_ends[0] = s1.bot_ends[1];
						let mut nonplanar = false;
						if let Some(s1_top0) = s1.top_ends[0] {
							if s1.top_depths[1] != Some(idx(lowval)) {
								debug_assert!(!is_type_1);
								nonplanar = true;
							} else {
								debug_assert!(s1.top_depths[0] == Some(idx(lowval)));
								self.planarity.link(s0.top_ends[0].expect("backedge on side 0"), s1_top0);
								s0.top_ends[0] = s1.top_ends[1];
								// Already true since the backedge was on side 0
								debug_assert!(s0.top_depths[0] == Some(idx(lowval)));
							}
						}
						if nonplanar {
							*self.tstacks.cur_planarity_mut() = Err(Nonplanar);
						} else {
							*s1 = PlanaritySide::default();
						}
					}
				}

				if let Some(item) = item {
					self.finish_tstack_top(item, false);
					is_single = true;
				}
			}
		} else {
			debug_assert!(is_type_1);
			// The span lives on side !edge_dir
			self.push_edge_tstack(cur, lowval, e, false);
			let fo = &mut self.first_occurrence[lowval];
			*fo = (*fo).min(self.nxt_edge_idx);
			self.nxt_edge_idx += 1;
		}

		// NB: We can do this check in lots of ways, maybe there's a cleaner check
		if is_type_1 && self.tstacks.len() >= 2 && self.tstacks.nxt().v_start == cur && self.tstacks.nxt().top_depth == lowval {
			let item = self.maybe_unwrap_nxt(NodeType::P, false);
			self.merge_tstack_tops();
			self.finish_tstack_top(item, false);
		}

		if !has_vert_tstack {
			// Throw the vertex onto the tstack so it'll get interleaved correctly
			self.push_vert_tstack(cur, cur_depth);
			self.stk.last_mut().expect("in a walk").has_vert_tstack = true;
			debug_assert!(!is_type_1);
			if !is_single {
				// Eagerly merge the vertex into the R to avoid a later spurious finish_tstack
				self.merge_tstack_tops();
			}
		}
	}

	fn pop_vert(&mut self) {
		let cur_depth = self.stk.len() - 1;
		let cur = self.stack_verts[cur_depth];
		let s = self.stk.pop().expect("in a walk");
		debug_assert!(s.edges.is_empty());
		if !s.has_vert_tstack {
			// Either our parent is a bridge edge, or we're just a root.
			// Leave it on the tstack for future cleanup; it'll get popped off immediately.
			// edge_dir == !stack_dir[lowval == cur_depth - 1] == true
			self.stack_dir[cur_depth] = true;
			self.push_vert_tstack(cur, cur_depth);
		}
	}

	/// Walk the DFS tree rooted at `root`, hanging its blocks off the forest root.
	fn run(&mut self, root: usize) {
		self.push_vert(root);
		loop {
			let s = self.stk.last().expect("in a walk");
			if s.edges.is_empty() {
				self.pop_vert();
				if self.stk.is_empty() {
					break;
				}
				self.finish_edge();
			} else if let Some(nxt) = self.start_edge() {
				self.push_vert(nxt);
			} else {
				self.finish_edge();
			}
		}
		let (t, _) = self.tstacks.pop();
		debug_assert!(self.tstacks.len() == 0);
		self.items.append_children(ROOT_ITEM, t.spans[1]);
	}
}

// ---------------------------------------------------------------------------------------------
// Phase 3: relabel the item tree in preorder and lay out the CSR outputs

const ZERO: Idx = match Idx::new(0) {
	Some(i) => i,
	None => unreachable!(),
};

struct RelabelFrame {
	cur_idx: Idx,
	/// Unassigned children of `cur_idx` in `ch.dat`
	ch: Range<u32>,
	/// Next nv / ne slot of `cur_idx` to hand to a child
	cur_nv: u32,
	cur_ne: u32,
}

/// Preorder numbering + CSR layout. Output arrays are preallocated to their exact final sizes and filled in
/// as items are numbered (`ch.dat` / `node_verts.dat[..].vert` temporarily hold item / original-vertex ids).
struct Relabel<'a, const WP: bool> {
	items: Items,
	planarity: Planarity,
	edges: &'a [[usize; 2]],

	vert_index: Vec<Option<Idx>>,
	edge_index: Vec<Option<Idx>>,
	edge_flipped: Vec<bool>,
	par: Vec<Option<Idx>>,
	subtree_end: Vec<Idx>,
	types: Vec<NodeType>,
	orig_id: Vec<Option<Idx>>,
	ch: Csr<Idx>,
	node_verts: Vec<NodeVert>,
	node_nvs: CsrIndex,
	vert_par_nv: Vec<Option<Idx>>,
	node_edges: Vec<NodeEdge>,
	node_nes: CsrIndex,
	node_adj: Csr<NodeAdj>,
	node_planar: Vec<bool>,
	ne_rot_adj: Vec<Option<Idx>>,

	// Scratch for R nodes
	vert_pos_buf: Vec<usize>,
	cnts_buf: Vec<usize>,
	ch_buf: Vec<(usize, Idx)>,
	/// ne of each vedge of the R node being laid out (`2 * ne` is its cap)
	rot_edge_ne: Vec<usize>,

	nxt_idx: usize,
	stk: Vec<RelabelFrame>,
}

impl<'a, const WP: bool> Relabel<'a, WP> {
	fn new(b: Builder<WP>, edges: &'a [[usize; 2]]) -> Self {
		let Builder { items, planarity, tot_blocks, tot_self_loops, .. } = b;
		let (nv, ne) = (items.nv, items.ne);
		let tot_items = items.len();
		// Each node is a child, and additionally most nodes have 2 cap verts; blocks and O nodes have 1
		let tot_node_verts = nv + (tot_items - 1 - nv) * 2 - tot_blocks - tot_self_loops;
		let tot_node_edges = (tot_items - 1 - nv - tot_blocks) * 2;
		fn csr<T: Clone>(n: usize, m: usize, fill: T) -> Csr<T> {
			Csr { bounds: vec![0; n + 1], dat: vec![fill; m] }
		}
		let wp = |n: usize| if WP { n } else { 0 };
		Relabel {
			items,
			planarity,
			edges,
			vert_index: vec![None; nv],
			edge_index: vec![None; ne],
			edge_flipped: vec![false; ne],
			par: vec![None; tot_items],
			subtree_end: vec![ZERO; tot_items],
			types: vec![NodeType::F; tot_items],
			orig_id: vec![None; tot_items],
			ch: csr(tot_items, tot_items - 1, ZERO),
			node_verts: vec![NodeVert { node: ZERO, vert: ZERO }; tot_node_verts],
			node_nvs: CsrIndex { bounds: vec![0; tot_items + 1] },
			vert_par_nv: vec![None; tot_items],
			node_edges: vec![NodeEdge { node: ZERO, twin_ne: ZERO, nvs: [ZERO; 2] }; tot_node_edges],
			node_nes: CsrIndex { bounds: vec![0; tot_items + 1] },
			node_adj: csr(2 * tot_node_verts, 2 * tot_node_edges, NodeAdj { ne: ZERO, dest_nv: ZERO }),
			node_planar: vec![false; wp(tot_items)],
			ne_rot_adj: vec![None; wp(4 * tot_node_edges)],
			vert_pos_buf: vec![0; nv],
			cnts_buf: Vec::with_capacity(2 * nv),
			ch_buf: Vec::with_capacity(tot_items),
			rot_edge_ne: vec![0; wp(2 * ne + 1)],
			nxt_idx: 0,
			stk: Vec::with_capacity(tot_items),
		}
	}

	fn set_ne(&mut self, cur_idx: usize, ne: usize, nvs: [usize; 2], nds: [usize; 2], rot_adjs: [Option<Idx>; 4]) {
		let nvs = nvs.map(idx);
		self.node_edges[ne].node = idx(cur_idx);
		self.node_edges[ne].nvs = nvs;
		self.node_adj.dat[nds[0]] = NodeAdj { ne: idx(ne), dest_nv: nvs[1] };
		self.node_adj.dat[nds[1]] = NodeAdj { ne: idx(ne), dest_nv: nvs[0] };
		if WP {
			self.ne_rot_adj[4 * ne..4 * ne + 4].copy_from_slice(&rot_adjs);
		}
	}

	/// Quarter-edge neighbours of vedge `ve` of the current R node, in ne_rot_adj numbering.
	fn map_rot_edge(&self, planar: bool, ve: usize) -> [Option<Idx>; 4] {
		if !WP || !planar {
			return [None; 4];
		}
		std::array::from_fn(|z| {
			let o = self.planarity.matches[4 * ve + z].expect("planar node: all quarter edges matched").usize();
			Some(idx((self.rot_edge_ne[o >> 2] << 2) + (o & 2) + usize::from(z & 1 == 0)))
		})
	}

	fn push_item(&mut self, cur_item: usize) {
		let ne = self.items.ne;
		let cur_idx = self.nxt_idx;
		self.nxt_idx += 1;
		let cur_type = self.items.types[cur_item];
		self.types[cur_idx] = cur_type;
		let mut planar = true;
		match cur_type {
			NodeType::F => debug_assert!(cur_item == ROOT_ITEM),
			NodeType::V => {
				let orig_vert = cur_item - 1;
				self.orig_id[cur_idx] = Some(idx(orig_vert));
				self.vert_index[orig_vert] = Some(idx(cur_idx));
			}
			NodeType::Q => {
				let orig_edge = self.items.vedge(cur_item);
				self.orig_id[cur_idx] = Some(idx(orig_edge));
				self.edge_index[orig_edge] = Some(idx(cur_idx));
				let v0 = self.items.vs[cur_item][0].expect("edge endpoint");
				self.edge_flipped[orig_edge] = v0.usize() != self.edges[orig_edge][0];
			}
			NodeType::I | NodeType::O => {
				// No planarity data was set up
			}
			NodeType::S | NodeType::P | NodeType::R => {
				if WP {
					match self.planarity.node_planarity[self.items.node(cur_item)] {
						Ok(p) => {
							// Expose our cap as vedge `2 * ne`. Make sure this runs before the planarity_flip checks.
							for (s, &b) in p.iter().enumerate() {
								self.planarity.link(idx(8 * ne + s), b);
							}
						}
						// TODO: Any certificate stuff
						Err(Nonplanar) => planar = false,
					}
				}
			}
		}
		if WP {
			self.node_planar[cur_idx] = planar;
		}

		// Fill ch and node_verts with items / original verts for now, since the children aren't numbered yet.
		let ch_st = self.ch.bounds[cur_idx];
		let mut ch_en = ch_st;
		let nv_st = self.node_nvs.bounds[cur_idx];
		let mut nv_en = nv_st;
		let mut n_edges = 0;
		let mut push_nv = |node_verts: &mut Vec<NodeVert>, vert: usize| {
			node_verts[nv_en] = NodeVert { node: idx(cur_idx), vert: idx(vert) };
			nv_en += 1;
		};
		let cur_vs = self.items.vs[cur_item];
		if let Some(v) = cur_vs[0] {
			push_nv(&mut self.node_verts, v.usize());
		}
		if let Some(Span { first, last }) = self.items.ch[cur_item] {
			let mut planarity_flip = first.flip;
			let mut ch_item = first.item.usize();
			loop {
				self.ch.dat[ch_en] = idx(ch_item);
				ch_en += 1;
				if self.items.is_vert(ch_item) {
					push_nv(&mut self.node_verts, ch_item - 1);
				} else {
					if WP {
						if cur_type != NodeType::R {
							debug_assert!(!planarity_flip);
						} else if planarity_flip {
							// Fix the planarity direction right here: reverse quarter_edge_matches upfront;
							// this breaks the involution property, but from here on we'll never read the low bits anyways.
							let m = &mut self.planarity.matches[4 * self.items.vedge(ch_item)..][..4];
							m.swap(0, 1);
							m.swap(2, 3);
						}
					}
					n_edges += 1;
				}
				if ch_item == last.item.usize() {
					debug_assert!(self.items.ch_nxt[ch_item].is_none());
					break;
				}
				let nxt = self.items.ch_nxt[ch_item].expect("child list continues to its last item");
				planarity_flip ^= nxt.flip;
				ch_item = nxt.item.usize();
			}
			planarity_flip ^= last.flip;
			debug_assert!(!planarity_flip);
		}
		if let Some(v) = cur_vs[1] {
			push_nv(&mut self.node_verts, v.usize());
		}
		self.ch.bounds[cur_idx + 1] = ch_en;
		self.node_nvs.bounds[cur_idx + 1] = nv_en;

		let n_verts = nv_en - nv_st;
		let is_node = cur_type.is_node();
		// Every node but a leaf Q has a vedge to its parent
		let has_cap = is_node && !(cur_type == NodeType::Q && ch_en > ch_st);
		if !is_node {
			n_edges = 0;
		}
		if has_cap {
			n_edges += 1;
		}

		let ne_st = self.node_nes.bounds[cur_idx];
		let ne_en = ne_st + n_edges;
		self.node_nes.bounds[cur_idx + 1] = ne_en;

		let adj_bounds = |this: &mut Self, bounds: [usize; 4]| {
			this.node_adj.bounds[2 * nv_st + 1..2 * nv_st + 5].copy_from_slice(&bounds);
		};
		match cur_type {
			NodeType::F => {
				// Just set node_adj bounds and we're good
				self.node_adj.bounds[2 * nv_st + 1..=2 * nv_en].fill(2 * ne_st);
			}
			NodeType::V => {}
			_ if n_verts == 1 => {
				debug_assert!(matches!(cur_type, NodeType::Q | NodeType::O));
				debug_assert!(n_edges == 1);
				self.node_adj.bounds[2 * nv_st + 1] = 2 * ne_st + 1;
				self.node_adj.bounds[2 * nv_st + 2] = 2 * ne_st + 2;
				let q = |k: usize| Some(idx(4 * ne_st + k));
				self.set_ne(cur_idx, ne_st, [nv_st, nv_st], [2 * ne_st + 1, 2 * ne_st], [q(3), q(2), q(1), q(0)]);
			}
			NodeType::O => unreachable!("O nodes have a single vertex"),
			NodeType::Q | NodeType::I => {
				debug_assert!(n_verts == 2);
				debug_assert!(n_edges == 1);
				adj_bounds(self, [2 * ne_st, 2 * ne_st + 1, 2 * ne_st + 2, 2 * ne_st + 2]);
				let q = |k: usize| Some(idx(4 * ne_st + k));
				self.set_ne(cur_idx, ne_st, [nv_st, nv_st + 1], [2 * ne_st, 2 * ne_st + 1], [q(1), q(0), q(3), q(2)]);
			}
			NodeType::P => {
				// Special case: tiebreak the parallel edges so they're reversed
				debug_assert!(n_verts == 2);
				debug_assert!(n_edges >= 3);
				adj_bounds(self, [2 * ne_st, 2 * ne_st + n_edges, 2 * ne_en, 2 * ne_en]);
				for (i, ne) in (ne_st..ne_en).enumerate() {
					let ne_prv = if ne == ne_st { ne_en - 1 } else { ne - 1 };
					let ne_nxt = if ne + 1 == ne_en { ne_st } else { ne + 1 };
					let q = |e: usize, k: usize| Some(idx(4 * e + k));
					let rot_adjs = [q(ne_prv, 1), q(ne_nxt, 0), q(ne_nxt, 3), q(ne_prv, 2)];
					self.set_ne(cur_idx, ne, [nv_st, nv_st + 1], [2 * ne_st + i, 2 * ne_en - 1 - i], rot_adjs);
				}
			}
			NodeType::S => {
				debug_assert!(n_verts == n_edges);
				debug_assert!(n_verts >= 3);
				for i in 2 * nv_st + 1..=2 * nv_en {
					self.node_adj.bounds[i] = i - 2 * nv_st + 2 * ne_st;
				}
				// Fix bounds for the cap
				self.node_adj.bounds[2 * nv_st + 1] -= 1;
				self.node_adj.bounds[2 * nv_en - 1] += 1;
				let q = |e: usize, k: usize| Some(idx(4 * e + k));
				self.set_ne(cur_idx, ne_st, [nv_st, nv_en - 1], [2 * ne_st, 2 * ne_en - 1], [q(ne_st + 1, 1), q(ne_st + 1, 0), q(ne_en - 1, 3), q(ne_en - 1, 2)]);
				for i in 1..n_edges {
					let ne = ne_st + i;
					// The cycle closes through the cap, which sits at ne_st
					let (prv, nxt) = (ne - 1, if ne + 1 == ne_en { ne_st } else { ne + 1 });
					let rot_adjs = if prv == ne_st {
						[q(ne_st, 1), q(ne_st, 0), q(nxt, 1), q(nxt, 0)]
					} else if nxt == ne_st {
						[q(prv, 3), q(prv, 2), q(ne_st, 3), q(ne_st, 2)]
					} else {
						[q(prv, 3), q(prv, 2), q(nxt, 1), q(nxt, 0)]
					};
					self.set_ne(cur_idx, ne, [nv_st + i - 1, nv_st + i], [2 * ne - 1, 2 * ne], rot_adjs);
				}
			}
			NodeType::R => {
				debug_assert!(has_cap);
				// Bucketsort the children by the midpoint
				for nv_ in nv_st..nv_en {
					self.vert_pos_buf[self.node_verts[nv_].vert.usize()] = nv_;
				}
				self.cnts_buf.clear();
				self.cnts_buf.resize(n_verts * 2 - 1, 0);
				self.ch_buf.clear();

				// Cap node_adj bounds
				self.node_adj.bounds[2 * nv_st + 2] += 1;
				self.node_adj.bounds[2 * nv_en - 1] += 1;

				for i in ch_st..ch_en {
					let item = self.ch.dat[i];
					let nvs = if self.items.is_vert(item.usize()) {
						[self.vert_pos_buf[item.usize() - 1]; 2]
					} else {
						let nvs = self.items.vs[item.usize()].map(|v| self.vert_pos_buf[v.expect("node endpoint").usize()]);
						debug_assert!(nvs[0] < nvs[1]);
						self.node_adj.bounds[2 * nvs[0] + 2] += 1;
						self.node_adj.bounds[2 * nvs[1] + 1] += 1;
						nvs
					};
					let loc = (nvs[0] - nv_st) + (nvs[1] - nv_st);
					self.ch_buf.push((loc, item));
					self.cnts_buf[loc] += 1;
				}
				let mut offset = ch_st;
				for cnt in &mut self.cnts_buf {
					offset += *cnt;
					*cnt = offset;
				}
				for &(loc, item) in self.ch_buf.iter().rev() {
					self.cnts_buf[loc] -= 1;
					self.ch.dat[self.cnts_buf[loc]] = item;
				}

				if WP {
					// Set up the reverse mapping for ourselves
					let mut nxt_ne = ne_en;
					for i in (ch_st..ch_en).rev() {
						let item = self.ch.dat[i].usize();
						if self.items.is_vert(item) {
							continue;
						}
						nxt_ne -= 1;
						self.rot_edge_ne[self.items.vedge(item)] = nxt_ne;
					}
					debug_assert!(nxt_ne == ne_st + 1);
					self.rot_edge_ne[2 * ne] = ne_st;
				}

				let mut off = 2 * ne_st;
				for b in &mut self.node_adj.bounds[2 * nv_st + 1..=2 * nv_en] {
					off += std::mem::replace(b, off);
				}
				debug_assert!(off == 2 * ne_en);

				// Fill in node_edges and node_adj, in reverse order to get the adj in bracket ordering.
				// The cap is special: it's first in node_edges, which means it's in the wrong place for the left endpoint.
				self.node_adj.bounds[2 * nv_st + 2] += 1;
				let mut nxt_ne = ne_en;
				for i in (ch_st..ch_en).rev() {
					let item = self.ch.dat[i].usize();
					if self.items.is_vert(item) {
						continue;
					}
					nxt_ne -= 1;
					let nvs = self.items.vs[item].map(|v| self.vert_pos_buf[v.expect("node endpoint").usize()]);
					let nd0 = self.node_adj.bounds[2 * nvs[0] + 2];
					self.node_adj.bounds[2 * nvs[0] + 2] += 1;
					let nd1 = self.node_adj.bounds[2 * nvs[1] + 1];
					self.node_adj.bounds[2 * nvs[1] + 1] += 1;
					let rot_adjs = self.map_rot_edge(planar, self.items.vedge(item));
					self.set_ne(cur_idx, nxt_ne, nvs, [nd0, nd1], rot_adjs);
				}
				debug_assert!(nxt_ne == ne_st + 1);

				// Insert the cap / bump its bound
				let rot_adjs = self.map_rot_edge(planar, 2 * ne);
				self.set_ne(cur_idx, ne_st, [nv_st, nv_en - 1], [2 * ne_st, 2 * ne_en - 1], rot_adjs);
				self.node_adj.bounds[2 * nv_en - 1] += 1;
			}
		}

		let cur_nv = nv_st + usize::from(cur_vs[0].is_some());
		let cur_ne = ne_st + usize::from(has_cap);
		self.stk.push(RelabelFrame { cur_idx: idx(cur_idx), ch: u32_of(ch_st)..u32_of(ch_en), cur_nv: u32_of(cur_nv), cur_ne: u32_of(cur_ne) });
	}

	/// Number the next child of the top frame and hook it up to its parent; returns the child item.
	fn start_child(&mut self) -> usize {
		let s = self.stk.last_mut().expect("in a walk");
		let cur_idx = s.cur_idx.usize();
		let ch_idx = s.ch.next().expect("child pending") as usize;
		let nxt_item = self.ch.dat[ch_idx].usize();
		let nxt_idx = self.nxt_idx;
		self.ch.dat[ch_idx] = idx(nxt_idx);
		self.par[nxt_idx] = Some(idx(cur_idx));
		if self.items.is_vert(nxt_item) {
			self.vert_par_nv[nxt_idx] = Some(idx(s.cur_nv as usize));
			s.cur_nv += 1;
		} else if self.types[cur_idx].is_node() {
			let (cur_ne, nxt_ne) = (s.cur_ne as usize, self.node_nes.bounds[nxt_idx]);
			self.node_edges[cur_ne].twin_ne = idx(nxt_ne);
			self.node_edges[nxt_ne].twin_ne = idx(cur_ne);
			s.cur_ne += 1;
		}
		nxt_item
	}

	fn pop_item(&mut self) {
		let s = self.stk.pop().expect("in a walk");
		debug_assert!(s.ch.is_empty());
		self.subtree_end[s.cur_idx.usize()] = idx(self.nxt_idx);
	}

	fn run(mut self) -> PlanarSpqrTree {
		self.push_item(ROOT_ITEM);
		while let Some(s) = self.stk.last() {
			if s.ch.is_empty() {
				self.pop_item();
			} else {
				let nxt = self.start_child();
				self.push_item(nxt);
			}
		}

		debug_assert!(self.nxt_idx == self.items.len());
		debug_assert!(self.ch.bounds.last() == Some(&self.ch.dat.len()));
		debug_assert!(self.node_nvs.bounds.last() == Some(&self.node_verts.len()));
		debug_assert!(self.node_nes.bounds.last() == Some(&self.node_edges.len()));
		debug_assert!(self.node_adj.bounds.last() == Some(&self.node_adj.dat.len()));

		let vert_index: Vec<Idx> = self.vert_index.into_iter().map(|i| i.expect("every vertex is an item")).collect();
		let edge_index = self.edge_index.into_iter().map(|i| i.expect("every edge is an item")).collect();
		// Rewrite node_verts to item indices
		for v in &mut self.node_verts {
			v.vert = vert_index[v.vert.usize()];
		}

		PlanarSpqrTree {
			tree: SpqrTree {
				vert_index,
				edge_index,
				edge_flipped: self.edge_flipped,
				par: self.par,
				subtree_end: self.subtree_end,
				types: self.types,
				orig_id: self.orig_id,
				ch: self.ch,
				node_verts: self.node_verts,
				node_nvs: self.node_nvs,
				vert_par_nv: self.vert_par_nv,
				node_edges: self.node_edges,
				node_nes: self.node_nes,
				node_adj: self.node_adj,
			},
			node_planar: self.node_planar,
			ne_embedding: PlanarEmbedding { rot_adj: self.ne_rot_adj },
		}
	}
}

fn build<const WP: bool>(nv: usize, edges: &[[usize; 2]], ternarize: bool, vert_order: &[usize], edge_order: &[usize]) -> PlanarSpqrTree {
	assert!(vert_order.len() <= nv);
	assert!(edge_order.len() <= edges.len());

	let (roots, outedges) = sorted_skeleton(nv, edges, vert_order, edge_order);

	let mut b = Builder::<WP>::new(nv, edges.len(), ternarize, outedges);
	for rt in roots {
		b.run(rt);
	}

	Relabel::new(b, edges).run()
}
