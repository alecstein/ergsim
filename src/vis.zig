// Visualization with zig-webui

const std = @import("std");
const sim = @import("sim.zig");

const webui = @import("webui");
const html = @embedFile("index.html");

const math = std.math;
const Allocator = std.mem.Allocator;
const alloc = std.heap.page_allocator;

const default_n = 1000;

pub const Col = sim.Col;
pub const ColType = sim.ColType;

pub fn main() !void {
    var nwin = webui.newWindow;
    _ = nwin.show(html);

    var timer = try std.time.Timer.start();
    const t0 = timer.read();

    const params = try getArgs();
    const n = params.n;

    var pct_elapsed: u8 = 0;
    _ = pct_elapsed;
    var progress = std.Progress{};
    const root_node = progress.start("Simulating", 100);
    defer root_node.end();

    std.debug.print("Particles: {d}\n", .{n});

    var xs = try alloc.alloc(f64, n);
    var vs = try alloc.alloc(f64, n);
    var cols = try alloc.alloc(Col, n);
    defer {
        alloc.free(xs);
        alloc.free(vs);
        alloc.free(cols);
    }

    sim.initArrays(xs, vs, cols);

    var t: f64 = 0;
    var ct: usize = 0; // collision count

    while (t < max_time) {
        const dt = stepForward(cols, xs, vs);
        t += dt;
        ct += 1;

        const xsBytes = std.mem.sliceAsBytes(xs);
        const vsBytes = std.mem.sliceAsBytes(vs);

        if (ct % 500 == 0) {
            nwin.sendRaw(
                "updateGasDensityHistogram",
                xsBytes,
            );

            nwin.sendRaw(
                "updateMomentumHistogram",
                vsBytes,
            );
        }

        // webui.wait();

        webui.clean();
    }

    const time_taken = @as(f64, @floatFromInt(timer.read() - t0));
    std.debug.print("\x1b[2K\x1b[0G", .{}); // clear the line of updates
    std.debug.print("Collisions: {}\n", .{ct});
    std.debug.print("Time taken: {d:.2}s\n", .{time_taken / 1000000000.0});
}

fn getArgs() !struct { n: u32 } {
    var arg_it = std.process.args();
    _ = arg_it.skip();
    var arg = arg_it.next() orelse return .{ .n = default_n };
    const n = try std.fmt.parseInt(u32, arg, 10);
    return .{ .n = n };
}

fn stepForward(cols: []Col, xs: []f64, vs: []f64) f64 {
    const next_col = sim.nextCollision(cols);
    const dt = next_col.dt;
    const j = next_col.j;
    const col_type = next_col.type;

    sim.advanceXsVs(next_col.dt, xs, vs);

    if (col_type == ColType.ground) {
        sim.computeGroundCol(j, dt, xs, vs, cols);
    } else {
        sim.computePistCol(j, dt, xs, vs, cols);
    }
    return dt;
}
