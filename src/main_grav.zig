// Simulates the case where the particles aren't affected by gravity

const std = @import("std");
const sim = @import("sim.zig");
const math = std.math;
const Allocator = std.mem.Allocator;
const alloc = std.heap.page_allocator;

const default_n = 1000;
const default_t = 1000.0;

pub const Col = sim.Col;
pub const ColType = sim.ColType;

pub fn main() !void {
    var timer = try std.time.Timer.start();
    const t0 = timer.read();

    const params = try getArgs();
    const n = params.n;
    const max_time = params.t;

    var pct_elapsed: u8 = 0;
    var progress = std.Progress{};
    const root_node = progress.start("Simulating", 100);
    defer root_node.end();

    std.debug.print("Particles: {d}\tWorldtime: {d}s\n", .{ n, max_time });

    var xs = try alloc.alloc(f64, n);
    var vs = try alloc.alloc(f64, n);
    var cols = try alloc.alloc(Col, n);
    defer {
        alloc.free(xs);
        alloc.free(vs);
        alloc.free(cols);
    }

    sim.initArrays(0.001, xs, vs, cols);

    var t: f64 = 0;
    var ct: usize = 0; // collision count

    while (t < max_time) {
        const dt = stepForward(cols, xs, vs);
        updateProgress(root_node, &pct_elapsed, t, max_time);
        t += dt;
        ct += 1;
    }

    const time_taken = @as(f64, @floatFromInt(timer.read() - t0));
    std.debug.print("\x1b[2K\x1b[0G", .{}); // clear the line of updates
    std.debug.print("Collisions: {}\n", .{ct});
    std.debug.print("Time taken: {d:.2}s\n", .{time_taken / 1000000000.0});
}

fn getArgs() !struct { n: u32, t: f64 } {
    var arg_it = std.process.args();
    _ = arg_it.skip();
    var arg = arg_it.next() orelse return .{ .n = default_n, .t = default_t };
    const n = try std.fmt.parseInt(u32, arg, 10);
    arg = arg_it.next() orelse return .{ .n = n, .t = default_t };
    const t = try std.fmt.parseFloat(f64, arg);
    return .{ .n = n, .t = t };
}

fn updateProgress(root_node: *std.Progress.Node, pct_elapsed: *u8, t: f64, max_time: f64) void {
    const new_pct_elapsed = 100 * t / max_time;
    const int_pct_elapsed = @as(u8, @intFromFloat(new_pct_elapsed));
    if (int_pct_elapsed > pct_elapsed.*) {
        pct_elapsed.* = int_pct_elapsed;
        root_node.completeOne();
    }
}

fn stepForward(cols: []Col, xs: []f64, vs: []f64) f64 {
    const next_col = sim.nextCollision(cols);
    const dt = next_col.dt;
    const j = next_col.j;
    const col_type = next_col.type;

    sim.advanceXsPsWithGravity(next_col.dt, xs, vs);

    if (col_type == ColType.ground) {
        sim.computeGroundColWithGravity(j, dt, xs, vs, cols);
    } else {
        sim.computePistColWithGravity(j, dt, xs, vs, cols);
    }
    return dt;
}
