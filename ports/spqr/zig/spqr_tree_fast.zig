//! SPQR tree (performance variant: see `Stack` / `fill` below; otherwise identical to spqr_tree.zig). Port of `wala::spqr_tree` / `wala::planar_spqr_tree` (cp-book `graph/spqr_tree.hpp`).
//!
//! All ids are `i32` with `-1` as the "none" sentinel, exactly like the C++ original.
//! The C++ `[&]` lambdas are expressed as methods on small state structs (`LowvalDfs`, `Builder`, `Relabel`);
//! `with_planarity` is a comptime bool, so the non-planar build compiles the planarity code out (like `if constexpr`).
//! `std::expected<T, Nonplanarity>` is expressed as `?T` (null == nonplanar); the certificate is empty in the original too.
//!
//! All output arrays are owned slices allocated from the `gpa` passed to `build`; free them with `deinit`.

const std = @import("std");
const Allocator = std.mem.Allocator;
const ArrayList = std.ArrayList;
const assert = std.debug.assert;

inline fn ix(i: i32) usize {
    return @intCast(i);
}
inline fn bi(b: bool) i32 {
    return @intFromBool(b);
}

/// Performance variant: fixed-capacity stack over a preallocated slice (a plain `[]T` plus a length).
/// Every stack in this file has a known exact upper bound, so this replaces `ArrayList`: pushes never
/// allocate (and thus never fail), and everything inlines. Bounds are asserted in debug builds.
fn Stack(comptime T: type) type {
    return struct {
        buf: []T,
        len: usize = 0,

        const Self = @This();

        fn init(gpa: Allocator, cap: usize) Allocator.Error!Self {
            return .{ .buf = try gpa.alloc(T, cap) };
        }
        fn deinit(self: *Self, gpa: Allocator) void {
            gpa.free(self.buf);
        }
        inline fn push(self: *Self, v: T) void {
            assert(self.len < self.buf.len);
            self.buf[self.len] = v;
            self.len += 1;
        }
        fn pushNTimes(self: *Self, v: T, n: usize) void {
            assert(self.len + n <= self.buf.len);
            fill(self.buf[self.len .. self.len + n], v);
            self.len += n;
        }
        inline fn pop(self: *Self) T {
            assert(self.len > 0);
            self.len -= 1;
            return self.buf[self.len];
        }
        inline fn top(self: *Self) *T {
            assert(self.len > 0);
            return &self.buf[self.len - 1];
        }
        inline fn at(self: *Self, i: usize) *T {
            assert(i < self.len);
            return &self.buf[i];
        }
        inline fn items(self: *Self) []T {
            return self.buf[0..self.len];
        }
    };
}

/// `@memset` lowers to compiler_rt's memset, which is several times slower than glibc's for
/// non-zero patterns; for 4-byte elements do the fill with explicit vector stores instead.
fn fill(s: anytype, val: anytype) void {
    const T = @TypeOf(s[0]);
    if (@sizeOf(T) == 4 and @bitSizeOf(T) == 32) {
        const V = @Vector(8, T);
        const vv: V = @splat(val);
        var i: usize = 0;
        while (i + 8 <= s.len) : (i += 8) {
            s[i..][0..8].* = vv;
            // Keep LLVM's loop-idiom pass from turning this back into a memset call.
            asm volatile ("" ::: .{ .memory = true });
        }
        while (i < s.len) : (i += 1) s[i] = val;
    } else {
        @memset(s, val);
    }
}

/// The row bounds of a jagged array whose entries live in a separate flat slice.
pub const CsrIndex = struct {
    bounds: []i32,

    pub fn size(self: CsrIndex) i32 {
        return @intCast(self.bounds.len - 1);
    }
    /// Equivalent of `csr_index[i]`: the half-open range of row `i`.
    pub fn range(self: CsrIndex, i: i32) [2]i32 {
        return .{ self.bounds[ix(i)], self.bounds[ix(i) + 1] };
    }
    /// Row `i` of `base`.
    pub fn slice(self: CsrIndex, i: i32, base: anytype) @TypeOf(base) {
        return base[ix(self.bounds[ix(i)])..ix(self.bounds[ix(i) + 1])];
    }
    pub fn deinit(self: *CsrIndex, gpa: Allocator) void {
        gpa.free(self.bounds);
        self.* = undefined;
    }
};

pub fn Csr(comptime T: type) type {
    return struct {
        const Self = @This();

        bounds: []i32,
        dat: []T,

        pub fn size(self: Self) i32 {
            return @intCast(self.bounds.len - 1);
        }
        /// Equivalent of `csr[i]`.
        pub fn slice(self: Self, i: i32) []T {
            return self.dat[ix(self.bounds[ix(i)])..ix(self.bounds[ix(i) + 1])];
        }
        pub fn deinit(self: *Self, gpa: Allocator) void {
            gpa.free(self.bounds);
            gpa.free(self.dat);
            self.* = undefined;
        }
    };
}

pub fn CsrBuilder(comptime T: type) type {
    return struct {
        const Self = @This();

        bounds: []i32,
        dat: []T = &.{},

        pub fn init(gpa: Allocator, n: i32) Allocator.Error!Self {
            const bounds = try gpa.alloc(i32, ix(n) + 1);
            fill(bounds, 0);
            return .{ .bounds = bounds };
        }

        pub fn count(self: *Self, k: i32) void {
            self.bounds[ix(k) + 1] += 1;
        }
        pub fn allocate(self: *Self, gpa: Allocator) Allocator.Error!void {
            var l: i32 = 0;
            for (self.bounds[1..]) |*b| {
                const c = b.*;
                b.* = l;
                l += c;
            }
            self.dat = try gpa.realloc(self.dat, ix(l));
        }
        pub fn push(self: *Self, k: i32) *T {
            const idx = self.bounds[ix(k) + 1];
            self.bounds[ix(k) + 1] += 1;
            return &self.dat[ix(idx)];
        }
        pub fn finalize(self: Self) Csr(T) {
            return .{ .bounds = self.bounds, .dat = self.dat };
        }
    };
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
pub const SpqrTree = struct {
    vert_index: []i32,
    edge_index: []i32,
    /// Whether edge e is stored in its Q node with its endpoints swapped: `node_verts[node_nvs.bounds[edge_index[e]]].vert != edges[e][0]`
    edge_flipped: []bool,

    par: []i32,
    subtree_end: []i32,
    types: []NodeType,
    orig_id: []i32,

    ch: Csr(i32),
    /// nv's of each node: `node_nvs.slice(i, node_verts)`
    node_verts: []NodeVert,
    node_nvs: CsrIndex,
    /// The nv index of a vertex within its parent node
    vert_par_nv: []i32,
    // TODO: Should we store a vert_nodes CSR?
    /// ne's of each node: `node_nes.slice(i, node_edges)`
    node_edges: []NodeEdge,
    node_nes: CsrIndex,
    node_adj: Csr(NodeAdj),

    pub fn size(self: SpqrTree) i32 {
        return @intCast(self.par.len);
    }

    /// `vert_order` and `edge_order` are (prefixes of) permutations of vertex / edge ids;
    /// listed ids are visited first in the given order, then the rest in id order.
    /// Roots are the first unvisited vertices, and DFS children are explored in edge order.
    /// Use `PlanarSpqrTree.build` to also compute the planar embeddings.
    pub fn build(gpa: Allocator, NV: i32, edges: []const [2]i32, ternarize: bool, vert_order: []const i32, edge_order: []const i32) Allocator.Error!SpqrTree {
        const t = try buildImpl(false, gpa, NV, edges, ternarize, vert_order, edge_order);
        gpa.free(t.node_planar);
        gpa.free(t.ne_embedding.rot_adj);
        return t.tree;
    }

    pub fn deinit(self: *SpqrTree, gpa: Allocator) void {
        gpa.free(self.vert_index);
        gpa.free(self.edge_index);
        gpa.free(self.edge_flipped);
        gpa.free(self.par);
        gpa.free(self.subtree_end);
        gpa.free(self.types);
        gpa.free(self.orig_id);
        self.ch.deinit(gpa);
        gpa.free(self.node_verts);
        self.node_nvs.deinit(gpa);
        gpa.free(self.vert_par_nv);
        gpa.free(self.node_edges);
        self.node_nes.deinit(gpa);
        self.node_adj.deinit(gpa);
        self.* = undefined;
    }
};

pub const NodeType = enum(u8) {
    F = 'F',
    V = 'V',
    Q = 'Q',
    I = 'I',
    O = 'O',
    S = 'S',
    P = 'P',
    R = 'R',

    pub fn char(self: NodeType) u8 {
        return @intFromEnum(self);
    }
};

pub const NodeVert = struct {
    node: i32,
    vert: i32,
};

pub const NodeEdge = struct {
    node: i32,
    twin_ne: i32,
    // TODO: Should we store the twin node, the twin node type, and/or twin node type == Q?
    nvs: [2]i32,
};

pub const NodeAdj = struct {
    ne: i32,
    dest_nv: i32,
};

/// A planar embedding of a graph, as the rotation system of its quarter-edges.
/// Quarter-edges are indexed by 4 * edge + 2 * side + dir (side: v0 vs v1, dir: cw vs ccw):
/// qe ^ 1 is the other side around the endpoint, qe ^ 3 the other side along the edge, and rot_adj[qe] the facing quarter-edge.
/// Partial embeddings are represented with -1's in rot_adj.
pub const PlanarEmbedding = struct {
    rot_adj: []i32,
};

pub const PlanarSpqrTree = struct {
    tree: SpqrTree,
    node_planar: []bool,
    /// Planarity adjacencies of each node's vedges, indexed according to:
    /// ne_embedding.rot_adj[4 * node_edge + 2 * side + dir]
    /// Nonplanar nodes have all entries -1.
    ne_embedding: PlanarEmbedding,

    pub fn build(gpa: Allocator, NV: i32, edges: []const [2]i32, ternarize: bool, vert_order: []const i32, edge_order: []const i32) Allocator.Error!PlanarSpqrTree {
        return buildImpl(true, gpa, NV, edges, ternarize, vert_order, edge_order);
    }

    pub fn deinit(self: *PlanarSpqrTree, gpa: Allocator) void {
        self.tree.deinit(gpa);
        gpa.free(self.node_planar);
        gpa.free(self.ne_embedding.rot_adj);
        self.* = undefined;
    }
};

fn setmin(a: *i32, b: i32) void {
    if (b < a.*) a.* = b;
}

/// Yields the ids in `order`, then the remaining ids in [0, n) in increasing order.
/// (The C++ `for_each_in_order` helper, as an iterator since Zig has no closures.)
const OrderIter = struct {
    n: i32,
    order: []const i32,
    listed: []bool,
    pos: usize = 0,
    rest: i32 = 0,

    fn init(gpa: Allocator, n: i32, order: []const i32) Allocator.Error!OrderIter {
        var listed: []bool = &.{};
        if (order.len >= 2 and order.len != ix(n)) {
            listed = try gpa.alloc(bool, ix(n));
            fill(listed, false);
            for (order) |i| listed[ix(i)] = true;
        }
        return .{ .n = n, .order = order, .listed = listed };
    }
    fn deinit(self: *OrderIter, gpa: Allocator) void {
        gpa.free(self.listed);
    }
    fn next(self: *OrderIter) ?i32 {
        if (self.pos < self.order.len) {
            const i = self.order[self.pos];
            self.pos += 1;
            return i;
        }
        if (self.order.len == ix(self.n)) return null;
        while (self.rest < self.n) {
            const i = self.rest;
            self.rest += 1;
            const skip = if (self.order.len == 0) false else if (self.order.len == 1) i == self.order[0] else self.listed[ix(i)];
            if (!skip) return i;
        }
        return null;
    }
};

// Helpers for working with [2]T - these compile to cmov's better than direct index access.

/// return arr[dir] == a, arr[!dir] == b
fn setSides(dir: bool, a: anytype, b: @TypeOf(a)) [2]@TypeOf(a) {
    const r: [2]@TypeOf(a) = if (dir) .{ b, a } else .{ a, b };
    return r;
}
fn getSide(a: anytype, dir: bool) @TypeOf(a[0]) {
    return if (dir) a[1] else a[0];
}

// ---------------------------------------------------------------------------------------------
// Phase 1: build a sorted skeleton

const OutEdge = struct {
    src: i32,
    dest: i32,
    e: i32,
    key: i32,
};

const AdjEdge = struct {
    dest: i32,
    e: i32,
};

/// Return the 2 lowvals from this subtree
const LowvalStack = struct {
    cur: i32,
    prv_e: i32,
    lowvals: [2]i32,
    ch_idx: i32,
    ch_end: i32,
};

const LowvalDfs = struct {
    gpa: Allocator,
    adj: *const Csr(AdjEdge),
    depth: []i32,
    all_outedges: Stack(OutEdge),
    stk: Stack(LowvalStack),

    fn pushVert(self: *LowvalDfs, cur: i32, prv_e: i32) void {
        const d: i32 = @intCast(self.stk.len);
        self.depth[ix(cur)] = d;
        self.stk.push(.{
            .cur = cur,
            .prv_e = prv_e,
            .lowvals = .{ d, d },
            .ch_idx = self.adj.bounds[ix(cur)],
            .ch_end = self.adj.bounds[ix(cur) + 1],
        });
    }

    fn finishEdge(self: *LowvalDfs, is_tree: bool, n_lowvals: [2]i32) void {
        const d: i32 = @intCast(self.stk.len - 1);
        const s = &self.stk.buf[self.stk.len - 1];
        const cur = s.cur;
        assert(s.ch_idx < s.ch_end);
        const ae = self.adj.dat[ix(s.ch_idx)];
        const nxt = ae.dest;
        const e = ae.e;
        const lowvals = &s.lowvals;
        s.ch_idx += 1;

        {
            // Extra bit is 0 for type-1 children, 1 for backedges, 2 for children with lowval2
            // Bridges have lowval -2 (kind 0), and components loops have lowval -1 (components are kind 0, loops are kind 1)
            // We don't really need to distinguish backedges vs type-1 children, but do it just for fun?
            var lowval = n_lowvals[0];
            if (lowval >= d) {
                lowval = ~(lowval - d);
            }
            const kind = 2 * bi(n_lowvals[1] < d) + bi(!is_tree);
            self.all_outedges.push(.{ .src = cur, .dest = nxt, .e = e, .key = 3 * (lowval + 2) + kind });
        }

        // Keep the 2 distinct mins
        if (n_lowvals[0] < lowvals[0]) {
            // NB: two statements, not an array literal: Zig's result-location semantics would alias lowvals[0].
            lowvals[1] = @min(n_lowvals[1], lowvals[0]);
            lowvals[0] = n_lowvals[0];
        } else {
            lowvals[1] = @min(lowvals[1], if (n_lowvals[0] == lowvals[0]) n_lowvals[1] else n_lowvals[0]);
        }
    }

    fn startEdge(self: *LowvalDfs) void {
        const d: i32 = @intCast(self.stk.len - 1);
        const s = &self.stk.buf[self.stk.len - 1];
        assert(s.ch_idx < s.ch_end);
        const ae = self.adj.dat[ix(s.ch_idx)];
        const nxt = ae.dest;
        const e = ae.e;

        if (e == s.prv_e or self.depth[ix(nxt)] > d) {
            // skip the edge
            s.ch_idx += 1;
            return;
        }

        const is_tree = self.depth[ix(nxt)] == -1;
        if (is_tree) {
            self.pushVert(nxt, e);
        } else {
            self.finishEdge(false, .{ self.depth[ix(nxt)], d });
        }
    }

    fn popVert(self: *LowvalDfs) [2]i32 {
        return self.stk.pop().lowvals;
    }
};

// ---------------------------------------------------------------------------------------------
// Phase 2: do the big ear-decomposition-like walk

// We're going to build a tree of all SPQR *nodes* + all original *vertices* (collectively *items*).
// Vertices will hang off the first SPQR node containing them, and blocks will be rooted at a topmost Q node for the top edge.

// As we build, we will represent the children of our nodes/vertices as linked lists.
const ROOT_ITEM: i32 = 0;

const ItemList = struct {
    // Items are actually 2 * item + planarity_flip (always 0 without planarity)
    v: [2]i32 = .{ -1, -1 },

    const empty: ItemList = .{};

    fn isEmpty(self: ItemList) bool {
        return self.v[0] < 0;
    }
};

fn unitList(item: i32) ItemList {
    return .{ .v = .{ item << 1, item << 1 } };
}

const TstackPlanaritySide = struct {
    // For each side, store pointers to the "linked lists" of the edges inside.
    // v[0] is the outer / longer edges and v[1] is the inner / shorter edges, matching the outside-in sort order.

    // bot_ends are the outer/innermost exposed pieces of the walk down the ear in the tree (they're connected to the bottommost/topmost vertices of the tree path)
    bot_ends: [2]i32 = .{ -1, -1 },
    // top_ends are the outer/innermost exposed backedges
    top_ends: [2]i32 = .{ -1, -1 },
    // depths should be increasing going inwards
    top_depths: [2]i32 = .{ -1, -1 },
};

const TstackPlanarity = struct {
    // The convention is that sides[0].top_depths[0] == top_depth, i.e. at least one minimal return lives on side 0
    sides: [2]TstackPlanaritySide = .{ .{}, .{} },
};

// TODO: What's the nonplanarity certificate look like?
// `null` means nonplanar.
const TstackMaybePlanarity = ?TstackPlanarity;

/// Like the C++ `std::conditional_t<with_planarity, std::expected<...>, std::monostate>` member,
/// the planarity data is zero-sized in the non-planar build, keeping this at 28 bytes there.
fn Tstack(comptime WP: bool) type {
    return struct {
        v_start: i32,
        top_depth: i32,
        first_idx: i32,
        spans: [2]ItemList,
        planarity: if (WP) TstackMaybePlanarity else void,
    };
}

const Key = struct {
    lowval: i32,
    is_tree: bool,
    is_type_1: bool,
};

fn decodeKey(cur_depth: i32, key: i32) Key {
    var lowval = @divTrunc(key, 3) - 2;
    if (lowval < 0) {
        lowval = cur_depth + ~lowval;
    }
    const kind = @rem(key, 3);
    return .{ .lowval = lowval, .is_tree = kind != 1, .is_type_1 = kind <= 1 };
}

const WalkStack = struct {
    has_vert_tstack: bool,
    ch_idx: i32,
    ch_end: i32,
    orig_tstack: i32,
};

fn mergePlanarity(quarter_edge_matches: []i32, a: *TstackMaybePlanarity, b: *const TstackMaybePlanarity) void {
    const ap: *TstackPlanarity = if (a.*) |*ap| ap else return;
    const bp: *const TstackPlanarity = if (b.*) |*bp| bp else {
        a.* = null;
        return;
    };
    for (0..2) |z| {
        const as_ = &ap.sides[z];
        const bs = &bp.sides[z];
        // If there's no bottom edges, then we must be an isolated vertex, so we can end early.
        if (bs.bot_ends[0] == -1) {
            // Do nothing
        } else if (as_.bot_ends[0] == -1) {
            as_.* = bs.*;
        } else {
            quarter_edge_matches[ix(as_.bot_ends[1])] = bs.bot_ends[0];
            quarter_edge_matches[ix(bs.bot_ends[0])] = as_.bot_ends[1];
            as_.bot_ends[1] = bs.bot_ends[1];

            if (bs.top_ends[0] == -1) {
                // Do nothing
            } else if (as_.top_ends[0] == -1) {
                as_.top_ends = bs.top_ends;
                as_.top_depths = bs.top_depths;
            } else if (as_.top_depths[1] > bs.top_depths[0]) {
                // TODO: Certificate
                a.* = null;
                return;
            } else {
                quarter_edge_matches[ix(as_.top_ends[1])] = bs.top_ends[0];
                quarter_edge_matches[ix(bs.top_ends[0])] = as_.top_ends[1];
                as_.top_ends[1] = bs.top_ends[1];
                as_.top_depths[1] = bs.top_depths[1];
            }
        }
    }
}

fn Builder(comptime WP: bool) type {
    return struct {
        const TstackT = Tstack(WP);
        const Self = @This();

        gpa: Allocator,
        NV: i32,
        NE: i32,
        ternarize: bool,
        outedges: Csr(OutEdge),

        ch_nxt: Stack(i32),
        item_vs: Stack([2]i32),
        item_ch: Stack(ItemList),
        item_types: Stack(NodeType),

        // Quarter edges for planar embedding building.
        // Each vedge has 4 entries by 4 * vedge_id + 2 * source_vert + is_cw (is_cw is arbitrary)
        // vedges are identified with what item they cap, numbered by (item - 1 - NV)
        quarter_edge_matches: []i32,
        // null == nonplanar (the C++ NonplanarityCertificate is empty)
        node_planarity: Stack(?[4]i32),

        tot_blocks: i32 = 0,
        tot_self_loops: i32 = 0,

        // Declare these here: most of our code will be in terms of v_start / top_depth, so we'll want to read these out
        stack_verts: []i32,
        stack_dir: []bool,

        nxt_edge_idx: i32 = 0, // Counts backedges only
        first_occurrence: []i32, // First backedge to this depth

        edge_top_depths: []i32,

        tstack: Stack(TstackT),
        stk: Stack(WalkStack),

        fn vertItem(self: *const Self, v: i32) i32 {
            _ = self;
            return 1 + v;
        }
        fn edgeItem(self: *const Self, e: i32) i32 {
            return 1 + self.NV + e;
        }

        fn concat(self: *Self, a: ItemList, b: ItemList) ItemList {
            if (b.isEmpty()) return a;
            if (a.isEmpty()) return b;
            self.ch_nxt.buf[ix(a.v[1] >> 1)] = b.v[0] ^ (a.v[1] & 1);
            // Build the result in a local so the result location can alias `a` / `b` at the call site.
            const r: ItemList = .{ .v = .{ a.v[0], b.v[1] } };
            return r;
        }

        fn allocItem(self: *Self, ty: NodeType) i32 {
            const item: i32 = @intCast(self.item_vs.len);
            self.item_vs.push(.{ 0, 0 });
            self.item_ch.push(.{});
            self.item_types.push(ty);
            self.ch_nxt.push(-1);
            if (WP) {
                self.node_planarity.push(.{ 0, 0, 0, 0 });
            }
            return item;
        }

        fn makeVs(self: *const Self, v_start: i32, top_depth: i32) [2]i32 {
            return setSides(self.stack_dir[ix(top_depth)], self.stack_verts[ix(top_depth)], v_start);
        }

        fn makeEdgePlanarity(self: *Self, item: i32, top_depth: i32, is_tree: bool) TstackMaybePlanarity {
            if (WP) {
                assert(item >= 1 + self.NV);
                const ve = item - (1 + self.NV);
                const top_dir = bi(self.stack_dir[ix(top_depth)]);
                self.edge_top_depths[ix(ve)] = top_depth;
                var p: TstackPlanarity = .{};
                if (is_tree) {
                    p.sides[0].bot_ends = .{ 4 * ve + 2 * (1 - top_dir) + 0, 4 * ve + 2 * top_dir + 1 };
                    p.sides[1].bot_ends = .{ 4 * ve + 2 * (1 - top_dir) + 1, 4 * ve + 2 * top_dir + 0 };
                } else {
                    p.sides[0].bot_ends = .{ 4 * ve + 2 * (1 - top_dir) + 0, 4 * ve + 2 * (1 - top_dir) + 1 };
                    p.sides[0].top_ends = .{ 4 * ve + 2 * top_dir + 1, 4 * ve + 2 * top_dir + 0 };
                    p.sides[0].top_depths = .{ top_depth, top_depth };
                }
                return p;
            } else {
                return TstackPlanarity{};
            }
        }

        fn curTstack(self: *Self) *TstackT {
            return &self.tstack.buf[self.tstack.len - 1];
        }
        fn nxtTstack(self: *Self) *TstackT {
            return &self.tstack.buf[self.tstack.len - 2];
        }

        fn pushTstack(self: *Self, v_start: i32, top_depth: i32, item: i32, planarity: TstackMaybePlanarity) void {
            const spans = setSides(self.stack_dir[ix(top_depth)], unitList(item), ItemList.empty);
            self.tstack.push(.{ .v_start = v_start, .top_depth = top_depth, .first_idx = self.nxt_edge_idx, .spans = spans, .planarity = if (WP) planarity else {} });
        }
        fn pushVertTstack(self: *Self, v: i32, top_depth: i32) void {
            const item = self.vertItem(v);
            self.pushTstack(v, top_depth, item, TstackPlanarity{});
        }
        fn pushEdgeTstack(self: *Self, v_start: i32, top_depth: i32, e: i32, is_tree: bool) void {
            const item = self.edgeItem(e);
            const planarity = self.makeEdgePlanarity(item, top_depth, is_tree);
            self.pushTstack(v_start, top_depth, item, planarity);
        }
        fn flipTstackPlanarity(a: *TstackT) void {
            if (WP) {
                a.spans[0].v[0] ^= 1;
                a.spans[0].v[1] ^= 1;
                a.spans[1].v[0] ^= 1;
                a.spans[1].v[1] ^= 1;
                if (a.planarity) |*p| {
                    std.mem.swap(TstackPlanaritySide, &p.sides[0], &p.sides[1]);
                }
            }
        }
        fn mergeTstackTops(self: *Self) void {
            const b = self.tstack.pop();
            const a = self.curTstack();
            setmin(&a.top_depth, b.top_depth);
            a.spans[0] = self.concat(b.spans[0], a.spans[0]);
            a.spans[1] = self.concat(a.spans[1], b.spans[1]);
            if (WP) {
                mergePlanarity(self.quarter_edge_matches, &a.planarity, &b.planarity);
            }
        }

        fn maybeUnwrapNxt(self: *Self, ty: NodeType, is_tree: bool) i32 {
            if (ty == .R) {
                return self.allocItem(ty);
            }

            assert(ty == .P or ty == .S);

            // If we want to ternarize, never reuse.
            if (self.ternarize) {
                return self.allocItem(ty);
            }

            const t = self.nxtTstack();
            const top_dir = self.stack_dir[ix(t.top_depth)];
            assert(getSide(t.spans, !top_dir).isEmpty());
            const item = getSide(t.spans, top_dir).v[0] >> 1;
            assert(item == (getSide(t.spans, top_dir).v[1] >> 1));
            if (self.item_types.buf[ix(item)] == ty) {
                t.spans = setSides(top_dir, self.item_ch.buf[ix(item)], ItemList.empty);
                if (WP) {
                    // Unwrap the planarity data
                    // We don't really need to maintain this at all because S/P nodes are known to be trivially planar
                    // The current state is just makeEdgePlanarity(wrapped), which means that it has the right shape, just needs to be relabelled.
                    const matches = self.node_planarity.buf[ix(item - (1 + self.NV + self.NE))].?; // unwrapped S/P nodes are always planar
                    const p: *TstackPlanarity = if (t.planarity) |*p| p else unreachable; // unwrapped tstack is a single edge
                    const td: usize = @intFromBool(top_dir);
                    if (is_tree) {
                        p.sides[0].bot_ends[0] = matches[2 * (1 - td) + 1];
                        p.sides[0].bot_ends[1] = matches[2 * td + 0];
                        p.sides[1].bot_ends[0] = matches[2 * (1 - td) + 0];
                        p.sides[1].bot_ends[1] = matches[2 * td + 1];
                    } else {
                        p.sides[0].bot_ends[0] = matches[2 * (1 - td) + 1];
                        p.sides[0].bot_ends[1] = matches[2 * (1 - td) + 0];
                        p.sides[0].top_ends[0] = matches[2 * td + 0];
                        p.sides[0].top_ends[1] = matches[2 * td + 1];
                    }
                }
                return item;
            } else {
                return self.allocItem(ty);
            }
        }

        fn finishTstackTop(self: *Self, item: i32, is_tree: bool) void {
            const t = self.curTstack();
            const top_dir = self.stack_dir[ix(t.top_depth)];
            assert(getSide(t.spans, !top_dir).isEmpty());

            if (WP) {
                const np = ix(item - (1 + self.NV + self.NE));
                if (t.planarity) |p| {
                    const td: usize = @intFromBool(top_dir);
                    var matches: [4]i32 = .{ 0, 0, 0, 0 };
                    if (is_tree) {
                        matches[2 * (1 - td) + 1] = p.sides[0].bot_ends[0];
                        matches[2 * td + 0] = p.sides[0].bot_ends[1];
                        matches[2 * (1 - td) + 0] = p.sides[1].bot_ends[0];
                        matches[2 * td + 1] = p.sides[1].bot_ends[1];
                    } else {
                        matches[2 * (1 - td) + 1] = p.sides[0].bot_ends[0];
                        matches[2 * (1 - td) + 0] = p.sides[0].bot_ends[1];
                        matches[2 * td + 0] = p.sides[0].top_ends[0];
                        matches[2 * td + 1] = p.sides[0].top_ends[1];
                    }
                    self.node_planarity.buf[np] = matches;
                } else {
                    assert(self.item_types.buf[ix(item)] == .R);
                    self.node_planarity.buf[np] = null;
                }
            }
            self.item_vs.buf[ix(item)] = self.makeVs(t.v_start, t.top_depth);
            self.item_ch.buf[ix(item)] = getSide(t.spans, top_dir);

            t.spans = setSides(top_dir, unitList(item), ItemList.empty);
            if (WP) t.planarity = self.makeEdgePlanarity(item, t.top_depth, is_tree);
        }

        fn pushVert(self: *Self, cur: i32) void {
            self.stk.push(.{
                .has_vert_tstack = false,
                .ch_idx = self.outedges.bounds[ix(cur)],
                .ch_end = self.outedges.bounds[ix(cur) + 1],
                .orig_tstack = -1,
            });
            const cur_depth = self.stk.len - 1;
            self.stack_verts[cur_depth] = cur;
        }

        /// Returns nxt to jump to pushVert(nxt), or null to jump to finishEdge
        fn startEdge(self: *Self) ?i32 {
            const cur_depth: i32 = @intCast(self.stk.len - 1);
            const s = &self.stk.buf[self.stk.len - 1];
            const cur = self.stack_verts[ix(cur_depth)];
            assert(s.ch_idx < s.ch_end);
            const oe = self.outedges.dat[ix(s.ch_idx)];
            const nxt = oe.dest;
            const k = decodeKey(cur_depth, oe.key);

            // edge_dir convention: false is forwards, true is backwards.
            // That means that cur is on the edge_dir side and nxt is on the !edge_dir side.
            self.stack_dir[ix(cur_depth)] = if (k.lowval >= cur_depth) false else !self.stack_dir[ix(k.lowval)];

            if (!s.has_vert_tstack and k.lowval < cur_depth and k.is_type_1) {
                // Do this with the correct stack_dir set
                self.pushVertTstack(cur, cur_depth);
                s.has_vert_tstack = true;
            }

            s.orig_tstack = @intCast(self.tstack.len);
            if (k.is_tree) {
                self.first_occurrence[ix(cur_depth)] = self.NE;
                return nxt;
            } else {
                return null;
            }
        }

        fn finishEdge(self: *Self) void {
            const cur_depth: i32 = @intCast(self.stk.len - 1);
            const s = &self.stk.buf[self.stk.len - 1];
            const cur = self.stack_verts[ix(cur_depth)];
            assert(s.ch_idx < s.ch_end);

            const oe = self.outedges.dat[ix(s.ch_idx)];
            const nxt = oe.dest;
            const e = oe.e;
            s.ch_idx += 1;

            const k = decodeKey(cur_depth, oe.key);
            const lowval = k.lowval;
            const is_tree = k.is_tree;
            const is_type_1 = k.is_type_1;

            const orig_tstack = s.orig_tstack;
            const edge_dir = self.stack_dir[ix(cur_depth)];
            const ei = ix(self.edgeItem(e));

            if (lowval >= cur_depth) {
                // There's no planarity handling for this because it's just a Q node. I/O nodes also don't need any tracking.
                self.item_vs.buf[ei] = .{ cur, -1 };
                self.tot_blocks += 1;
                if (is_tree) {
                    // Bridges and components
                    if (lowval == cur_depth + 1) {
                        // tstack[tstack_size-1] is currently just smuggling out the child vertex, prepend the bridge component
                        // This is just a shortcut for allocating a full I-type tstack
                        const item = self.allocItem(.I);
                        self.item_vs.buf[ix(item)] = self.makeVs(nxt, cur_depth);
                        const t = self.tstack.pop();
                        self.item_ch.buf[ei] = self.concat(unitList(item), t.spans[1]);
                    } else {
                        // tstack[tstack_size-2] is the vertex and tstack[tstack_size-1] is the backedge
                        const backedge = self.tstack.pop().spans[0];
                        const t = self.tstack.pop();
                        self.item_ch.buf[ei] = self.concat(backedge, t.spans[1]);
                    }
                } else {
                    // self loops
                    assert(nxt == cur);
                    self.tot_self_loops += 1;
                    const item = self.allocItem(.O);
                    // Make sure the nxt is -1 as well
                    self.item_vs.buf[ix(item)] = .{ cur, -1 };
                    self.item_ch.buf[ei] = unitList(item);
                }
                const vi = ix(self.vertItem(cur));
                self.item_ch.buf[vi] = self.concat(self.item_ch.buf[vi], unitList(@intCast(ei)));
                return;
            }
            assert(lowval < cur_depth);

            self.item_vs.buf[ei] = self.makeVs(nxt, cur_depth);

            // Whether curTstack() is a single edge
            var is_single = true;
            if (is_tree) {
                // The span lives on side edge_dir
                self.pushEdgeTstack(nxt, cur_depth, e, true);
                while (self.tstack.len >= 2 and self.nxtTstack().top_depth >= cur_depth) {
                    const ty: NodeType = if (self.nxtTstack().top_depth > cur_depth) blk: {
                        // Just backfill this for maybeUnwrap
                        self.stack_dir[ix(self.nxtTstack().top_depth)] = edge_dir;

                        // The tstack currently contains a tree-edge followed by a vertex; merge the vertex first
                        self.mergeTstackTops();

                        break :blk .S;
                    } else if (self.nxtTstack().v_start == self.curTstack().v_start)
                        // This will be a P node
                        .P
                    else
                        .R;
                    const item = self.maybeUnwrapNxt(ty, ty == .S);
                    self.mergeTstackTops();
                    if (WP) {
                        if (self.curTstack().planarity) |*p| {
                            // Merge all backedges into the component
                            for (&p.sides) |*side| {
                                assert(side.bot_ends[1] != -1);
                                if (side.top_ends[1] == -1) continue;
                                assert(side.top_depths[0] == cur_depth);
                                assert(side.top_depths[1] == cur_depth);
                                self.quarter_edge_matches[ix(side.bot_ends[1])] = side.top_ends[1];
                                self.quarter_edge_matches[ix(side.top_ends[1])] = side.bot_ends[1];
                                side.bot_ends[1] = side.top_ends[0];
                                side.top_depths = .{ -1, -1 };
                                side.top_ends = .{ -1, -1 };
                            }
                        }
                    }
                    self.finishTstackTop(item, true);
                }

                if (self.curTstack().first_idx > self.first_occurrence[ix(cur_depth)]) {
                    while (self.curTstack().first_idx > self.first_occurrence[ix(cur_depth)]) {
                        if (WP) {
                            const n = self.tstack.len;
                            const ts = self.tstack.buf;
                            if (ts[n - 2].first_idx > self.first_occurrence[ix(cur_depth)]) {
                                // We will put cur_depth on side 1 until the bottom
                                if (ts[n - 2].top_depth == cur_depth) {
                                    flipTstackPlanarity(&ts[n - 2]);
                                }
                            } else if (!is_single) {
                                assert(ts[n - 1].top_depth < cur_depth);
                                if (ts[n - 2].planarity) |p| {
                                    if (p.sides[0].top_depths[1] == cur_depth) {
                                        // We need to flip curTstack and nxtTstack relative to each other.
                                        // Flip the one with worse top_depth.
                                        const which: usize = if (ts[n - 1].top_depth < ts[n - 2].top_depth) n - 2 else n - 1;
                                        flipTstackPlanarity(&ts[which]);
                                    } else {
                                        assert(p.sides[1].top_depths[1] == cur_depth);
                                    }
                                }
                            }
                        }
                        self.mergeTstackTops();
                        is_single = false;
                    }
                    if (WP) {
                        if (self.curTstack().planarity) |*p| {
                            // Prune off finished cur-side things
                            for (&p.sides) |*side| {
                                assert(side.bot_ends[1] != -1);
                                while (side.top_depths[1] == cur_depth) {
                                    {
                                        // Link these to bot_ends[1]
                                        self.quarter_edge_matches[ix(side.bot_ends[1])] = side.top_ends[1];
                                        self.quarter_edge_matches[ix(side.top_ends[1])] = side.bot_ends[1];
                                        side.bot_ends[1] = side.top_ends[1] ^ 1;
                                    }
                                    side.top_ends[1] = self.quarter_edge_matches[ix(side.bot_ends[1])];
                                    self.quarter_edge_matches[ix(side.bot_ends[1])] = -1;
                                    if (side.top_ends[1] != -1) {
                                        self.quarter_edge_matches[ix(side.top_ends[1])] = -1;
                                        side.top_depths[1] = self.edge_top_depths[ix(side.top_ends[1] >> 2)];
                                    } else {
                                        side.top_depths = .{ -1, -1 };
                                        side.top_ends = .{ -1, -1 };
                                    }
                                }
                            }
                        }
                    }
                }

                if (is_type_1) {
                    assert(s.has_vert_tstack);
                }
                if (s.has_vert_tstack) {
                    // NB: tstack[orig_size] is the vertex and tstack[orig_size+1] is the backedge; maybe we should reverse them?
                    assert(self.tstack.len >= ix(orig_tstack + 3));

                    if (!is_type_1) {
                        if (WP) {
                            // The lowval side should be side 1, everything else goes on side 0.
                            // The exception is tstack[orig_tstack + 2], which could be == lowval on one/both sides,
                            // but is guaranteed to have *something* > lowval by non-type-1-ness
                            const t = &self.tstack.buf[ix(orig_tstack + 2)];
                            if (t.planarity) |p| {
                                assert(p.sides[0].top_depths[0] == t.top_depth);
                                if (p.sides[0].top_depths[1] == lowval) {
                                    flipTstackPlanarity(t);
                                }
                                assert(t.planarity.?.sides[0].top_depths[1] != -1);
                                assert(t.planarity.?.sides[0].top_depths[1] > lowval);
                            }
                            for (self.tstack.items()[ix(orig_tstack + 3)..]) |*ti| {
                                if (ti.top_depth == lowval) {
                                    flipTstackPlanarity(ti);
                                }
                            }
                        }
                        while (self.tstack.len > ix(orig_tstack + 3)) {
                            self.mergeTstackTops();
                            is_single = false;
                        }
                        assert(!is_single);
                    }

                    assert(self.tstack.len == ix(orig_tstack + 3));
                    const item: i32 = if (is_type_1)
                        self.maybeUnwrapNxt(if (is_single) .S else .R, false)
                    else
                        // Just for the type checker
                        -1;
                    // Merge with the backedge
                    self.mergeTstackTops();
                    // Merge with the vertex
                    self.mergeTstackTops();

                    const c = self.curTstack();
                    c.v_start = cur;
                    assert(c.top_depth == lowval);

                    // Fold everything to the correct side now that we're leaving the child.
                    // The entire subtree should go to the !edge_dir side.
                    const all = self.concat(c.spans[0], c.spans[1]);
                    c.spans = setSides(!edge_dir, all, ItemList.empty);

                    if (WP) {
                        if (c.planarity) |*p| {
                            // precondition: side 1 should be the lowval only side
                            const s0 = &p.sides[0];
                            const s1 = &p.sides[1];
                            self.quarter_edge_matches[ix(s0.bot_ends[0])] = s1.bot_ends[0];
                            self.quarter_edge_matches[ix(s1.bot_ends[0])] = s0.bot_ends[0];
                            s0.bot_ends[0] = s1.bot_ends[1];
                            var nonplanar = false;
                            if (s1.top_ends[0] != -1) {
                                if (s1.top_depths[1] != lowval) {
                                    assert(!is_type_1);
                                    nonplanar = true;
                                } else {
                                    assert(s1.top_depths[0] == lowval);
                                    self.quarter_edge_matches[ix(s0.top_ends[0])] = s1.top_ends[0];
                                    self.quarter_edge_matches[ix(s1.top_ends[0])] = s0.top_ends[0];
                                    s0.top_ends[0] = s1.top_ends[1];
                                    // Already true since the backedge was on side 0
                                    assert(s0.top_depths[0] == lowval);
                                }
                            }
                            if (nonplanar) {
                                c.planarity = null;
                            } else {
                                s1.* = .{};
                            }
                        }
                    }

                    if (is_type_1) {
                        self.finishTstackTop(item, false);
                        is_single = true;
                    }
                }
            } else {
                assert(is_type_1);
                // The span lives on side !edge_dir
                self.pushEdgeTstack(cur, lowval, e, false);
                const idx = self.nxt_edge_idx;
                self.nxt_edge_idx += 1;
                setmin(&self.first_occurrence[ix(lowval)], idx);
            }

            // NB: We can do this check in lots of ways, maybe there's a cleaner check
            if (is_type_1 and self.tstack.len >= 2 and self.nxtTstack().v_start == cur and self.nxtTstack().top_depth == lowval) {
                // This will be a P node
                const item = self.maybeUnwrapNxt(.P, false);
                self.mergeTstackTops();
                self.finishTstackTop(item, false);
            }

            if (!s.has_vert_tstack) {
                // Throw cur_vert_node onto the tstack so it'll get interleaved correctly
                self.pushVertTstack(cur, cur_depth);
                s.has_vert_tstack = true;
                assert(!is_type_1);
                if (!is_single) {
                    // Just eagerly merge the vertex into the R to avoid a later spurious finishTstack
                    self.mergeTstackTops();
                }
            }
        }

        fn popVert(self: *Self) void {
            const cur_depth: i32 = @intCast(self.stk.len - 1);
            const s = &self.stk.buf[self.stk.len - 1];
            const cur = self.stack_verts[ix(cur_depth)];
            assert(s.ch_idx == s.ch_end);
            if (!s.has_vert_tstack) {
                // Either our parent is a bridge edge, or we're just a root.
                // We'll just leave it on tstack for future cleanup, it'll just get popped of immediately.
                // edge_dir == !stack_dir[lowval == cur_depth - 1] == true
                self.stack_dir[ix(cur_depth)] = true;
                self.pushVertTstack(cur, cur_depth);
                s.has_vert_tstack = true;
            }
            _ = self.stk.pop();
        }

        fn deinit(self: *Self) void {
            const gpa = self.gpa;
            self.outedges.deinit(gpa);
            self.ch_nxt.deinit(gpa);
            self.item_vs.deinit(gpa);
            self.item_ch.deinit(gpa);
            self.item_types.deinit(gpa);
            gpa.free(self.quarter_edge_matches);
            self.node_planarity.deinit(gpa);
            gpa.free(self.stack_verts);
            gpa.free(self.stack_dir);
            gpa.free(self.first_occurrence);
            gpa.free(self.edge_top_depths);
            self.tstack.deinit(gpa);
            self.stk.deinit(gpa);
        }
    };
}

// ---------------------------------------------------------------------------------------------
// Phase 3: relabel the full tree in preorder

const ChBuf = struct {
    loc: i32,
    item_id: i32,
};

const RelabelStack = struct {
    cur_idx: i32,
    ch_idx: i32,
    ch_end: i32,
    cur_nv: i32,
    cur_ne: i32,
};

fn Relabel(comptime WP: bool) type {
    return struct {
        const Self = @This();

        gpa: Allocator,
        b: *Builder(WP),
        edges: []const [2]i32,

        vert_index: []i32,
        edge_index: []i32,
        edge_flipped: []bool,
        par: []i32,
        subtree_end: []i32,
        types: []NodeType,
        orig_id: []i32,
        ch: Csr(i32),
        node_verts: []NodeVert,
        node_nvs: CsrIndex,
        vert_par_nv: []i32,
        node_edges: []NodeEdge,
        node_nes: CsrIndex,
        node_adj: Csr(NodeAdj),
        node_planar: []bool,
        ne_rot_adj: []i32,

        vert_pos_buf: []i32,
        cnts_buf: Stack(i32),
        ch_buf: Stack(ChBuf),
        rot_edge_ne: []i32,

        nxt_unassigned_idx: i32 = 0,
        stk: Stack(RelabelStack),

        fn setNe(self: *Self, cur_idx: i32, ne: i32, nvs: [2]i32, nds: [2]i32, rot_adjs: [4]i32) void {
            self.node_edges[ix(ne)].node = cur_idx;
            self.node_edges[ix(ne)].nvs = nvs;
            self.node_adj.dat[ix(nds[0])] = .{ .ne = ne, .dest_nv = nvs[1] };
            self.node_adj.dat[ix(nds[1])] = .{ .ne = ne, .dest_nv = nvs[0] };
            if (WP) {
                for (0..4) |z| {
                    self.ne_rot_adj[ix(4 * ne) + z] = rot_adjs[z];
                }
            }
        }

        fn mapRotEdge(self: *const Self, planar: bool, ve: i32) [4]i32 {
            if (!WP) return .{ -1, -1, -1, -1 };
            if (!planar) return .{ -1, -1, -1, -1 };
            var res: [4]i32 = undefined;
            for (0..4) |z| {
                const o = self.b.quarter_edge_matches[ix(4 * ve) + z];
                assert(o != -1);
                res[z] = (self.rot_edge_ne[ix(o >> 2)] << 2) + (o & 2) + bi((z & 1) == 0);
            }
            return res;
        }

        fn pushItem(self: *Self, cur_item: i32) void {
            const NV = self.b.NV;
            const NE = self.b.NE;
            const cur_idx = self.nxt_unassigned_idx;
            self.nxt_unassigned_idx += 1;
            const cur_type = self.b.item_types.buf[ix(cur_item)];
            self.types[ix(cur_idx)] = cur_type;
            var planar = true;
            if (cur_type == .F) {
                assert(cur_item == 0);
            } else if (cur_type == .V) {
                assert(1 <= cur_item and cur_item < 1 + NV);
                const orig_vert = cur_item - 1;
                self.orig_id[ix(cur_idx)] = orig_vert;
                self.vert_index[ix(orig_vert)] = cur_idx;
            } else if (cur_type == .Q) {
                assert(1 + NV <= cur_item and cur_item < 1 + NV + NE);
                const orig_edge = cur_item - 1 - NV;
                self.orig_id[ix(cur_idx)] = orig_edge;
                self.edge_index[ix(orig_edge)] = cur_idx;
                assert(self.b.item_vs.buf[ix(cur_item)][0] != -1);
                self.edge_flipped[ix(orig_edge)] = self.b.item_vs.buf[ix(cur_item)][0] != self.edges[ix(orig_edge)][0];
            } else {
                assert(1 + NV + NE <= cur_item);
                if (WP) {
                    if (cur_type == .O or cur_type == .I) {
                        // No planarity data was set up
                    } else if (cur_type == .S or cur_type == .P or cur_type == .R) {
                        if (self.b.node_planarity.buf[ix(cur_item - (1 + NV + NE))]) |p| {
                            // Make sure this runs before our planarity_flip checks
                            for (0..4) |s| {
                                const a: i32 = 8 * NE + @as(i32, @intCast(s));
                                const bb = p[s];
                                self.b.quarter_edge_matches[ix(a)] = bb;
                                self.b.quarter_edge_matches[ix(bb)] = a;
                            }
                        } else {
                            // TODO: Any certificate stuff
                            planar = false;
                        }
                    } else {
                        unreachable;
                    }
                }
            }
            if (WP) {
                self.node_planar[ix(cur_idx)] = planar;
            }

            // HACK: Fill ch and vert_items in with orig items / orig verts for now,
            // because we don't have the final item id's yet.
            const ch_st = self.ch.bounds[ix(cur_idx)];
            var ch_en = ch_st;
            const nv_st = self.node_nvs.bounds[ix(cur_idx)];
            var nv_en = nv_st;
            var n_edges: i32 = 0;
            const cur_item_vs = self.b.item_vs.buf[ix(cur_item)];
            if (cur_item_vs[0] != -1) {
                self.node_verts[ix(nv_en)] = .{ .node = cur_idx, .vert = cur_item_vs[0] };
                nv_en += 1;
            }
            const cur_item_ch = self.b.item_ch.buf[ix(cur_item)];
            if (!cur_item_ch.isEmpty()) {
                var planarity_flip = (cur_item_ch.v[0] & 1) != 0;
                var ch_item = cur_item_ch.v[0] >> 1;
                while (true) {
                    self.ch.dat[ix(ch_en)] = ch_item;
                    ch_en += 1;
                    assert(ch_item >= 1);
                    if (ch_item < 1 + NV) {
                        self.node_verts[ix(nv_en)] = .{ .node = cur_idx, .vert = ch_item - 1 };
                        nv_en += 1;
                    } else {
                        if (WP) {
                            if (cur_type != .R) {
                                assert(!planarity_flip);
                            } else {
                                // Fix the planarity direction right here: reverse quarter_edge_matches upfront;
                                // this breaks the involution property, but from here on we'll never read the low bits anyways.
                                const ve = ix(ch_item - (1 + NV));
                                if (planarity_flip) {
                                    const qem = self.b.quarter_edge_matches;
                                    std.mem.swap(i32, &qem[4 * ve + 0], &qem[4 * ve + 1]);
                                    std.mem.swap(i32, &qem[4 * ve + 2], &qem[4 * ve + 3]);
                                }
                            }
                        }
                        n_edges += 1;
                    }
                    if (ch_item == (cur_item_ch.v[1] >> 1)) {
                        assert(self.b.ch_nxt.buf[ix(ch_item)] == -1);
                        break;
                    }
                    const nxt = self.b.ch_nxt.buf[ix(ch_item)];
                    planarity_flip = planarity_flip != ((nxt & 1) != 0);
                    ch_item = nxt >> 1;
                }
                planarity_flip = planarity_flip != ((cur_item_ch.v[1] & 1) != 0);
                assert(!planarity_flip);
            }
            if (cur_item_vs[1] != -1) {
                self.node_verts[ix(nv_en)] = .{ .node = cur_idx, .vert = cur_item_vs[1] };
                nv_en += 1;
            }
            self.ch.bounds[ix(cur_idx) + 1] = ch_en;
            self.node_nvs.bounds[ix(cur_idx) + 1] = nv_en;

            const n_verts = nv_en - nv_st;

            const is_node = cur_type != .F and cur_type != .V;
            const has_cap = is_node and !(cur_type == .Q and ch_en - ch_st > 0);

            if (!is_node) {
                n_edges = 0;
            }
            if (has_cap) {
                n_edges += 1;
            }

            const ne_st = self.node_nes.bounds[ix(cur_idx)];
            const ne_en = ne_st + n_edges;
            self.node_nes.bounds[ix(cur_idx) + 1] = ne_en;

            const nab = self.node_adj.bounds;
            if (cur_type == .F) {
                // Just set node_adj bounds and we're good
                var i = 2 * nv_st + 1;
                while (i <= 2 * nv_en) : (i += 1) {
                    nab[ix(i)] = 2 * ne_st;
                }
            } else if (cur_type == .V) {
                // Nothing to do
            } else if (n_verts == 1) {
                assert(cur_type == .Q or cur_type == .O);
                assert(n_edges == 1);
                nab[ix(2 * nv_st + 1)] = 2 * ne_st + 1 * n_edges;
                nab[ix(2 * nv_st + 2)] = 2 * ne_st + 2 * n_edges;
                self.setNe(cur_idx, ne_st, .{ nv_st, nv_st }, .{ 2 * ne_st + 1, 2 * ne_st }, .{ 4 * ne_st + 3, 4 * ne_st + 2, 4 * ne_st + 1, 4 * ne_st + 0 });
            } else if (cur_type == .Q or cur_type == .I) {
                assert(n_verts == 2);
                assert(n_edges == 1);
                nab[ix(2 * nv_st + 1)] = 2 * ne_st + 0 * n_edges;
                nab[ix(2 * nv_st + 2)] = 2 * ne_st + 1 * n_edges;
                nab[ix(2 * nv_st + 3)] = 2 * ne_st + 2 * n_edges;
                nab[ix(2 * nv_st + 4)] = 2 * ne_st + 2 * n_edges;
                self.setNe(cur_idx, ne_st, .{ nv_st, nv_st + 1 }, .{ 2 * ne_st, 2 * ne_st + 1 }, .{ 4 * ne_st + 1, 4 * ne_st + 0, 4 * ne_st + 3, 4 * ne_st + 2 });
            } else if (cur_type == .P) {
                // Special case: tiebreak the parallel edges so they're reversed
                assert(n_verts == 2);
                assert(n_edges >= 3);
                nab[ix(2 * nv_st + 1)] = 2 * ne_st + 0 * n_edges;
                nab[ix(2 * nv_st + 2)] = 2 * ne_st + 1 * n_edges;
                nab[ix(2 * nv_st + 3)] = 2 * ne_st + 2 * n_edges;
                nab[ix(2 * nv_st + 4)] = 2 * ne_st + 2 * n_edges;
                var ne = ne_st;
                while (ne < ne_en) : (ne += 1) {
                    const ne_prv = (if (ne == ne_st) ne_en else ne) - 1;
                    const ne_nxt = if (ne + 1 == ne_en) ne_st else ne + 1;
                    const rot_adjs: [4]i32 = .{ 4 * ne_prv + 1, 4 * ne_nxt + 0, 4 * ne_nxt + 3, 4 * ne_prv + 2 };
                    self.setNe(cur_idx, ne, .{ nv_st, nv_st + 1 }, .{ 2 * ne_st + (ne - ne_st), 2 * ne_en - 1 - (ne - ne_st) }, rot_adjs);
                }
            } else if (cur_type == .S) {
                assert(n_verts == n_edges);
                assert(n_verts >= 3);
                {
                    var i = 2 * nv_st + 1;
                    while (i <= 2 * nv_en) : (i += 1) {
                        nab[ix(i)] = i + 2 * (ne_st - nv_st);
                    }
                }
                // Fix bounds for the cap
                nab[ix(2 * nv_st + 1)] -= 1;
                nab[ix(2 * nv_en - 1)] += 1;
                self.setNe(cur_idx, ne_st, .{ nv_st, nv_en - 1 }, .{ 2 * ne_st, 2 * ne_en - 1 }, .{ 4 * (ne_st + 1) + 1, 4 * (ne_st + 1) + 0, 4 * (ne_en - 1) + 3, 4 * (ne_en - 1) + 2 });
                var i: i32 = 1;
                while (i < n_edges) : (i += 1) {
                    const ne = ne_st + i;
                    var rot_adjs: [4]i32 = .{ 4 * (ne - 1) + 3, 4 * (ne - 1) + 2, 4 * (ne + 1) + 1, 4 * (ne + 1) + 0 };
                    if (ne - 1 == ne_st) {
                        rot_adjs[0] = 4 * ne_st + 1;
                        rot_adjs[1] = 4 * ne_st + 0;
                    }
                    if (ne + 1 == ne_en) {
                        rot_adjs[2] = 4 * ne_st + 3;
                        rot_adjs[3] = 4 * ne_st + 2;
                    }
                    self.setNe(cur_idx, ne, .{ nv_st + i - 1, nv_st + i }, .{ 2 * ne - 1, 2 * ne }, rot_adjs);
                }
            } else if (cur_type == .R) {
                // Bucketsort the children by the midpoint
                {
                    var nv_ = nv_st;
                    while (nv_ < nv_en) : (nv_ += 1) {
                        self.vert_pos_buf[ix(self.node_verts[ix(nv_)].vert)] = nv_;
                    }
                }
                self.cnts_buf.len = 0;
                self.cnts_buf.pushNTimes(0, ix(n_verts * 2 - 1));
                self.ch_buf.len = 0;

                assert(has_cap);

                // Cap node_adj bounds
                nab[ix(2 * nv_st + 2)] += 1;
                nab[ix(2 * nv_en - 1)] += 1;

                {
                    var i = ch_st;
                    while (i < ch_en) : (i += 1) {
                        const item = self.ch.dat[ix(i)];
                        assert(item >= 1);
                        var nvs: [2]i32 = undefined;
                        if (item < 1 + NV) {
                            nvs = .{ self.vert_pos_buf[ix(item - 1)], self.vert_pos_buf[ix(item - 1)] };
                        } else {
                            const vs = self.b.item_vs.buf[ix(item)];
                            nvs = .{ self.vert_pos_buf[ix(vs[0])], self.vert_pos_buf[ix(vs[1])] };
                            assert(nvs[0] < nvs[1]);
                            nab[ix(2 * nvs[0] + 2)] += 1;
                            nab[ix(2 * nvs[1] + 1)] += 1;
                        }
                        const loc = (nvs[0] - nv_st) + (nvs[1] - nv_st);
                        self.ch_buf.push(.{ .loc = loc, .item_id = item });
                        self.cnts_buf.buf[ix(loc)] += 1;
                    }
                }
                {
                    var offset = ch_st;
                    for (self.cnts_buf.items()) |*cnt| {
                        offset += cnt.*;
                        cnt.* = offset;
                    }
                }
                {
                    var i = self.ch_buf.len;
                    while (i > 0) {
                        i -= 1;
                        const cb = self.ch_buf.buf[i];
                        self.cnts_buf.buf[ix(cb.loc)] -= 1;
                        self.ch.dat[ix(self.cnts_buf.buf[ix(cb.loc)])] = cb.item_id;
                    }
                }

                if (WP) {
                    // Set up the reverse mapping for ourselves
                    var nxt_ne = ne_en;
                    var i = ch_en;
                    while (i > ch_st) {
                        i -= 1;
                        const item = self.ch.dat[ix(i)];
                        assert(item >= 1);
                        if (item < 1 + NV) continue;
                        nxt_ne -= 1;
                        self.rot_edge_ne[ix(item - (1 + NV))] = nxt_ne;
                    }
                    assert(nxt_ne == ne_st + 1);
                    self.rot_edge_ne[ix(2 * NE)] = ne_st;
                }

                {
                    var off = 2 * ne_st;
                    var i = 2 * nv_st + 1;
                    while (i <= 2 * nv_en) : (i += 1) {
                        const old = nab[ix(i)];
                        nab[ix(i)] = off;
                        off += old;
                    }
                    assert(off == 2 * ne_en);
                }

                // Fill in node_edges and node_adj.
                // Reverse order to get the adj in bracket ordering.
                {
                    // Handle cap as special: it's first in the node_edges, which means it's in the wrong place for the left endpoint.
                    nab[ix(2 * nv_st + 2)] += 1;

                    var nxt_ne = ne_en;
                    var i = ch_en;
                    while (i > ch_st) {
                        i -= 1;
                        const item = self.ch.dat[ix(i)];
                        assert(item >= 1);
                        if (item < 1 + NV) continue;
                        nxt_ne -= 1;
                        const vs = self.b.item_vs.buf[ix(item)];
                        // TODO: Reuse this from the ch pass?
                        const nvs: [2]i32 = .{ self.vert_pos_buf[ix(vs[0])], self.vert_pos_buf[ix(vs[1])] };
                        const nd0 = nab[ix(2 * nvs[0] + 2)];
                        nab[ix(2 * nvs[0] + 2)] += 1;
                        const nd1 = nab[ix(2 * nvs[1] + 1)];
                        nab[ix(2 * nvs[1] + 1)] += 1;
                        const rot_adjs = self.mapRotEdge(planar, item - (1 + NV));
                        self.setNe(cur_idx, nxt_ne, nvs, .{ nd0, nd1 }, rot_adjs);
                    }
                    assert(nxt_ne == ne_st + 1);

                    // Insert the cap / bump its bound
                    const rot_adjs = self.mapRotEdge(planar, 2 * NE);
                    self.setNe(cur_idx, ne_st, .{ nv_st, nv_en - 1 }, .{ 2 * ne_st, 2 * ne_en - 1 }, rot_adjs);
                    nab[ix(2 * nv_en - 1)] += 1;
                }
            } else {
                unreachable;
            }

            const cur_nv = nv_st + bi(cur_item_vs[0] != -1);
            const cur_ne = ne_st + bi(has_cap);
            self.stk.push(.{ .cur_idx = cur_idx, .ch_idx = ch_st, .ch_end = ch_en, .cur_nv = cur_nv, .cur_ne = cur_ne });
        }

        fn startChild(self: *Self) i32 {
            const s = &self.stk.buf[self.stk.len - 1];
            const cur_idx = s.cur_idx;
            assert(s.ch_idx < s.ch_end);
            const nxt_item = self.ch.dat[ix(s.ch_idx)];
            const nxt_idx = self.nxt_unassigned_idx;
            self.ch.dat[ix(s.ch_idx)] = nxt_idx;
            self.par[ix(nxt_idx)] = cur_idx;
            const nxt_ne = self.node_nes.bounds[ix(nxt_idx)];
            if (nxt_item < 1 + self.b.NV) {
                self.vert_par_nv[ix(nxt_idx)] = s.cur_nv;
                s.cur_nv += 1;
            } else if (self.types[ix(cur_idx)] != .F and self.types[ix(cur_idx)] != .V) {
                self.node_edges[ix(s.cur_ne)].twin_ne = nxt_ne;
                self.node_edges[ix(nxt_ne)].twin_ne = s.cur_ne;
                s.cur_ne += 1;
            }

            s.ch_idx += 1;
            return nxt_item;
        }

        fn popItem(self: *Self) void {
            const s = self.stk.pop();
            assert(s.ch_idx == s.ch_end);
            self.subtree_end[ix(s.cur_idx)] = self.nxt_unassigned_idx;
        }
    };
}

fn allocFilled(gpa: Allocator, comptime T: type, n: usize, val: T) Allocator.Error![]T {
    const s = try gpa.alloc(T, n);
    fill(s, val);
    return s;
}

fn buildImpl(comptime WP: bool, gpa: Allocator, NV: i32, edges: []const [2]i32, ternarize: bool, vert_order: []const i32, edge_order: []const i32) Allocator.Error!PlanarSpqrTree {
    const NE: i32 = @intCast(edges.len);
    assert(vert_order.len <= ix(NV));
    assert(edge_order.len <= ix(NE));

    var roots: Stack(i32) = try .init(gpa, ix(NV));
    defer roots.deinit(gpa);
    var outedges: Csr(OutEdge) = undefined;

    // Phase 1: build a sorted skeleton
    {
        // 1a: build a normal adjacency list for the initial lowval dfs
        var adj_builder = try CsrBuilder(AdjEdge).init(gpa, NV);
        for (edges) |uv| {
            adj_builder.count(uv[0]);
            if (uv[0] != uv[1]) adj_builder.count(uv[1]);
        }
        try adj_builder.allocate(gpa);
        {
            var it = try OrderIter.init(gpa, NE, edge_order);
            defer it.deinit(gpa);
            while (it.next()) |e| {
                const u = edges[ix(e)][0];
                const v = edges[ix(e)][1];
                adj_builder.push(u).* = .{ .dest = v, .e = e };
                if (u != v) adj_builder.push(v).* = .{ .dest = u, .e = e };
            }
        }
        var adj = adj_builder.finalize();
        defer adj.deinit(gpa);

        var dfs: LowvalDfs = .{
            .gpa = gpa,
            .adj = &adj,
            .depth = try allocFilled(gpa, i32, ix(NV), -1),
            .all_outedges = try .init(gpa, ix(NE)),
            .stk = try .init(gpa, ix(NV)),
        };
        defer gpa.free(dfs.depth);
        defer dfs.stk.deinit(gpa);
        {
            var it = try OrderIter.init(gpa, NV, vert_order);
            defer it.deinit(gpa);
            while (it.next()) |rt| {
                if (dfs.depth[ix(rt)] == -1) {
                    roots.push(rt);
                    dfs.pushVert(rt, -1);
                    while (true) {
                        const s = dfs.stk.top().*;
                        if (s.ch_idx == s.ch_end) {
                            const lowvals = dfs.popVert();
                            if (dfs.stk.len == 0) break;
                            dfs.finishEdge(true, lowvals);
                        } else {
                            dfs.startEdge();
                        }
                    }
                }
            }
        }
        // Every edge produced exactly one outedge, so the stack is full; hand its buffer over.
        assert(dfs.all_outedges.len == dfs.all_outedges.buf.len);
        const all_outedges = dfs.all_outedges.buf;

        var by_key_builder = try CsrBuilder(OutEdge).init(gpa, 3 * NV + 6);
        for (all_outedges) |edge| by_key_builder.count(edge.key);
        try by_key_builder.allocate(gpa);
        for (all_outedges) |edge| by_key_builder.push(edge.key).* = edge;
        var by_key = by_key_builder.finalize();
        defer by_key.deinit(gpa);

        var by_src_builder = try CsrBuilder(OutEdge).init(gpa, NV);
        // Hack to reuse memory
        by_src_builder.dat = all_outedges;
        for (by_key.dat) |edge| by_src_builder.count(edge.src);
        try by_src_builder.allocate(gpa);
        for (by_key.dat) |edge| by_src_builder.push(edge.src).* = edge;
        outedges = by_src_builder.finalize();
    }

    // Phase 2: do the big ear-decomposition-like walk
    const n_items0 = ix(1 + NV + NE);
    var b: Builder(WP) = .{
        .gpa = gpa,
        .NV = NV,
        .NE = NE,
        .ternarize = ternarize,
        .outedges = outedges,
        .ch_nxt = try .init(gpa, n_items0 + ix(NE)),
        .item_vs = try .init(gpa, n_items0 + ix(NE)),
        .item_ch = try .init(gpa, n_items0 + ix(NE)),
        .item_types = try .init(gpa, n_items0 + ix(NE)),
        .quarter_edge_matches = try allocFilled(gpa, i32, if (WP) ix(8 * NE + 4) else 0, -1),
        .node_planarity = try .init(gpa, if (WP) ix(NE) else 0),
        .stack_verts = try allocFilled(gpa, i32, ix(NV), 0),
        .stack_dir = try allocFilled(gpa, bool, ix(NV), false),
        .first_occurrence = try allocFilled(gpa, i32, ix(NV), 0),
        .edge_top_depths = try allocFilled(gpa, i32, if (WP) ix(2 * NE) else 0, -1),
        .tstack = try .init(gpa, ix(NV + NE)),
        .stk = try .init(gpa, ix(NV)),
    };
    defer b.deinit();
    b.ch_nxt.pushNTimes(-1, n_items0);
    b.item_vs.pushNTimes(.{ -1, -1 }, n_items0);
    b.item_ch.pushNTimes(.{}, n_items0);
    b.item_types.pushNTimes(.F, 1);
    b.item_types.pushNTimes(.V, ix(NV));
    b.item_types.pushNTimes(.Q, ix(NE));

    for (roots.items()) |rt| {
        b.pushVert(rt);
        while (true) {
            const s = b.stk.top().*;
            if (s.ch_idx == s.ch_end) {
                b.popVert();
                if (b.stk.len == 0) break;
                b.finishEdge();
            } else if (b.startEdge()) |nxt| {
                b.pushVert(nxt);
            } else {
                b.finishEdge();
            }
        }
        const t = b.tstack.pop();
        b.item_ch.buf[ix(ROOT_ITEM)] = b.concat(b.item_ch.buf[ix(ROOT_ITEM)], t.spans[1]);
    }

    // Phase 3: relabel the full tree in preorder
    const tot_items: i32 = @intCast(b.item_types.len);
    const tot_blocks = b.tot_blocks;
    const tot_self_loops = b.tot_self_loops;

    // Each node is a child, and additionally most non-block node has 2 cap verts; blocks have 1, and O nodes have 1
    const tot_node_verts = NV + (tot_items - 1 - NV) * 2 - tot_blocks - tot_self_loops;
    const tot_node_edges = (tot_items - 1 - NV - tot_blocks) * 2;

    var r: Relabel(WP) = .{
        .gpa = gpa,
        .b = &b,
        .edges = edges,
        .vert_index = try allocFilled(gpa, i32, ix(NV), -1),
        .edge_index = try allocFilled(gpa, i32, ix(NE), -1),
        .edge_flipped = try allocFilled(gpa, bool, ix(NE), false),
        .par = try allocFilled(gpa, i32, ix(tot_items), -1),
        .subtree_end = try allocFilled(gpa, i32, ix(tot_items), -1),
        .types = try allocFilled(gpa, NodeType, ix(tot_items), .F),
        .orig_id = try allocFilled(gpa, i32, ix(tot_items), -1),
        .ch = .{ .bounds = try allocFilled(gpa, i32, ix(tot_items) + 1, 0), .dat = try allocFilled(gpa, i32, ix(tot_items - 1), 0) },
        .node_verts = try allocFilled(gpa, NodeVert, ix(tot_node_verts), .{ .node = 0, .vert = 0 }),
        .node_nvs = .{ .bounds = try allocFilled(gpa, i32, ix(tot_items) + 1, 0) },
        .vert_par_nv = try allocFilled(gpa, i32, ix(tot_items), -1),
        .node_edges = try allocFilled(gpa, NodeEdge, ix(tot_node_edges), .{ .node = 0, .twin_ne = 0, .nvs = .{ 0, 0 } }),
        .node_nes = .{ .bounds = try allocFilled(gpa, i32, ix(tot_items) + 1, 0) },
        .node_adj = .{ .bounds = try allocFilled(gpa, i32, ix(tot_node_verts * 2 + 1), 0), .dat = try allocFilled(gpa, NodeAdj, ix(tot_node_edges * 2), .{ .ne = 0, .dest_nv = 0 }) },
        .node_planar = try allocFilled(gpa, bool, if (WP) ix(tot_items) else 0, false),
        .ne_rot_adj = try allocFilled(gpa, i32, if (WP) ix(4 * tot_node_edges) else 0, -1),
        .vert_pos_buf = try allocFilled(gpa, i32, ix(NV), -1),
        .cnts_buf = try .init(gpa, ix(2 * NV)),
        .ch_buf = try .init(gpa, ix(tot_items)),
        .rot_edge_ne = try allocFilled(gpa, i32, if (WP) ix(2 * NE + 1) else 0, 0),
        .stk = try .init(gpa, ix(tot_items)),
    };
    defer gpa.free(r.vert_pos_buf);
    defer r.cnts_buf.deinit(gpa);
    defer r.ch_buf.deinit(gpa);
    defer gpa.free(r.rot_edge_ne);
    defer r.stk.deinit(gpa);

    r.par[ix(r.nxt_unassigned_idx)] = -1;
    r.pushItem(ROOT_ITEM);
    while (true) {
        const s = r.stk.top().*;
        if (s.ch_idx == s.ch_end) {
            r.popItem();
            if (r.stk.len == 0) break;
        } else {
            const nxt = r.startChild();
            r.pushItem(nxt);
        }
    }

    assert(r.nxt_unassigned_idx == tot_items);
    assert(r.ch.bounds[r.ch.bounds.len - 1] == @as(i32, @intCast(r.ch.dat.len)));
    assert(r.node_nvs.bounds[r.node_nvs.bounds.len - 1] == @as(i32, @intCast(r.node_verts.len)));
    assert(r.node_nes.bounds[r.node_nes.bounds.len - 1] == @as(i32, @intCast(r.node_edges.len)));
    assert(r.node_adj.bounds[r.node_adj.bounds.len - 1] == @as(i32, @intCast(r.node_adj.dat.len)));

    // Rewrite node_vertices to the correct index
    for (r.node_verts) |*v| {
        v.vert = r.vert_index[ix(v.vert)];
    }

    return .{
        .tree = .{
            .vert_index = r.vert_index,
            .edge_index = r.edge_index,
            .edge_flipped = r.edge_flipped,
            .par = r.par,
            .subtree_end = r.subtree_end,
            .types = r.types,
            .orig_id = r.orig_id,
            .ch = r.ch,
            .node_verts = r.node_verts,
            .node_nvs = r.node_nvs,
            .vert_par_nv = r.vert_par_nv,
            .node_edges = r.node_edges,
            .node_nes = r.node_nes,
            .node_adj = r.node_adj,
        },
        .node_planar = r.node_planar,
        .ne_embedding = .{ .rot_adj = r.ne_rot_adj },
    };
}
