const std = @import("std");
const spqr = @import("spqr_tree_fast.zig");

fn dumpInts(w: *std.Io.Writer, name: []const u8, v: []const i32) !void {
    try w.print("{s}:", .{name});
    for (v) |x| try w.print(" {d}", .{x});
    try w.writeAll("\n");
}

fn bytesEql(a: anytype, b: @TypeOf(a)) bool {
    return std.mem.eql(u8, std.mem.sliceAsBytes(a), std.mem.sliceAsBytes(b));
}

fn treesEqual(a: spqr.SpqrTree, b: spqr.SpqrTree) bool {
    return bytesEql(a.vert_index, b.vert_index) and bytesEql(a.edge_index, b.edge_index) and
        bytesEql(a.par, b.par) and bytesEql(a.subtree_end, b.subtree_end) and
        bytesEql(a.types, b.types) and bytesEql(a.orig_id, b.orig_id) and
        bytesEql(a.ch.bounds, b.ch.bounds) and bytesEql(a.ch.dat, b.ch.dat) and
        bytesEql(a.node_verts.bounds, b.node_verts.bounds) and bytesEql(a.node_verts.dat, b.node_verts.dat) and
        bytesEql(a.vert_par_nv, b.vert_par_nv) and
        bytesEql(a.node_edges.bounds, b.node_edges.bounds) and bytesEql(a.node_edges.dat, b.node_edges.dat) and
        bytesEql(a.node_adj.bounds, b.node_adj.bounds) and bytesEql(a.node_adj.dat, b.node_adj.dat);
}

pub fn main(init: std.process.Init) !void {
    const io = init.io;
    const gpa = init.gpa;
    var rbuf: [1 << 16]u8 = undefined;
    var stdin = std.Io.File.stdin().reader(io, &rbuf);
    const input = try stdin.interface.allocRemaining(gpa, .unlimited);
    defer gpa.free(input);

    var toks = std.mem.tokenizeAny(u8, input, " \t\r\n");
    const Next = struct {
        fn int(t: *std.mem.TokenIterator(u8, .any)) !i32 {
            return std.fmt.parseInt(i32, t.next().?, 10);
        }
    };
    const NV = try Next.int(&toks);
    const NE = try Next.int(&toks);
    const ternarize = (try Next.int(&toks)) != 0;
    const edges = try gpa.alloc([2]i32, @intCast(NE));
    defer gpa.free(edges);
    for (edges) |*e| e.* = .{ try Next.int(&toks), try Next.int(&toks) };
    const k = try Next.int(&toks);
    const vert_order = try gpa.alloc(i32, @intCast(k));
    defer gpa.free(vert_order);
    for (vert_order) |*x| x.* = try Next.int(&toks);
    const l = try Next.int(&toks);
    const edge_order = try gpa.alloc(i32, @intCast(l));
    defer gpa.free(edge_order);
    for (edge_order) |*x| x.* = try Next.int(&toks);

    var t = try spqr.PlanarSpqrTree.build(gpa, NV, edges, ternarize, vert_order, edge_order);
    defer t.deinit(gpa);

    var wbuf: [1 << 16]u8 = undefined;
    var stdout = std.Io.File.stdout().writer(io, &wbuf);
    const w = &stdout.interface;
    try dumpInts(w, "vert_index", t.tree.vert_index);
    try dumpInts(w, "edge_index", t.tree.edge_index);
    try dumpInts(w, "par", t.tree.par);
    try dumpInts(w, "subtree_end", t.tree.subtree_end);
    try w.writeAll("types:");
    for (t.tree.types) |x| try w.print(" {c}", .{x.char()});
    try w.writeAll("\n");
    try dumpInts(w, "orig_id", t.tree.orig_id);
    try dumpInts(w, "ch.bounds", t.tree.ch.bounds);
    try dumpInts(w, "ch.dat", t.tree.ch.dat);
    try dumpInts(w, "node_verts.bounds", t.tree.node_verts.bounds);
    try w.writeAll("node_verts.dat:");
    for (t.tree.node_verts.dat) |x| try w.print(" {d},{d}", .{ x.node, x.vert });
    try w.writeAll("\n");
    try dumpInts(w, "vert_par_nv", t.tree.vert_par_nv);
    try dumpInts(w, "node_edges.bounds", t.tree.node_edges.bounds);
    try w.writeAll("node_edges.dat:");
    for (t.tree.node_edges.dat) |x| try w.print(" {d},{d},{d},{d}", .{ x.node, x.twin_ne, x.nvs[0], x.nvs[1] });
    try w.writeAll("\n");
    try dumpInts(w, "node_adj.bounds", t.tree.node_adj.bounds);
    try w.writeAll("node_adj.dat:");
    for (t.tree.node_adj.dat) |x| try w.print(" {d},{d}", .{ x.ne, x.dest_nv });
    try w.writeAll("\n");
    try w.writeAll("node_planar:");
    for (t.node_planar) |x| try w.print(" {d}", .{@intFromBool(x)});
    try w.writeAll("\n");
    try dumpInts(w, "ne_rot_adj", t.ne_rot_adj);

    var s = try spqr.SpqrTree.build(gpa, NV, edges, ternarize, vert_order, edge_order);
    defer s.deinit(gpa);
    try w.print("nonplanar_build_same: {d}\n", .{@intFromBool(treesEqual(s, t.tree))});
    try w.flush();
}
