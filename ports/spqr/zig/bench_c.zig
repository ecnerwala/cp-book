const std = @import("std");
const spqr = @import("spqr_tree.zig");

pub fn main(init: std.process.Init) !void {
    const io = init.io;
    const gpa = std.heap.c_allocator;
    var args = try std.process.Args.Iterator.initAllocator(init.minimal.args, gpa);
    defer args.deinit();
    _ = args.next();
    const R: usize = if (args.next()) |a| try std.fmt.parseInt(usize, a, 10) else 5;

    var rbuf: [1 << 16]u8 = undefined;
    var stdin = std.Io.File.stdin().reader(io, &rbuf);
    const input = try stdin.interface.allocRemaining(gpa, .unlimited);
    defer gpa.free(input);
    var toks = std.mem.tokenizeAny(u8, input, " \t\r\n");
    const NV = try std.fmt.parseInt(i32, toks.next().?, 10);
    const NE = try std.fmt.parseInt(i32, toks.next().?, 10);
    const edges = try gpa.alloc([2]i32, @intCast(NE));
    defer gpa.free(edges);
    for (edges) |*e| e.* = .{ try std.fmt.parseInt(i32, toks.next().?, 10), try std.fmt.parseInt(i32, toks.next().?, 10) };

    var wbuf: [1 << 12]u8 = undefined;
    var stdout = std.Io.File.stdout().writer(io, &wbuf);
    const w = &stdout.interface;

    inline for (.{ false, true }) |planar| {
        var best: f64 = std.math.floatMax(f64);
        var sink: usize = 0;
        for (0..R) |_| {
            const t0 = std.Io.Timestamp.now(io, .awake);
            if (planar) {
                var t = try spqr.PlanarSpqrTree.build(gpa, NV, edges, false, &.{}, &.{});
                const ns: i96 = t0.durationTo(std.Io.Timestamp.now(io, .awake)).nanoseconds;
                sink += t.tree.par.len;
                t.deinit(gpa);
                best = @min(best, @as(f64, @floatFromInt(ns)) / 1e6);
            } else {
                var t = try spqr.SpqrTree.build(gpa, NV, edges, false, &.{}, &.{});
                const ns: i96 = t0.durationTo(std.Io.Timestamp.now(io, .awake)).nanoseconds;
                sink += t.par.len;
                t.deinit(gpa);
                best = @min(best, @as(f64, @floatFromInt(ns)) / 1e6);
            }
        }
        try w.print("zig  {s:<10} {d:8.2} ms  (items={d})\n", .{ if (planar) "planar" else "spqr", best, sink / R });
    }
    try w.flush();
}
