// Visualization with zig-webui

const std = @import("std");
const sim = @import("sim.zig");

const webui = @import("webui");
const html = @embedFile("index.html");
var win: webui = undefined;

const math = std.math;
const Allocator = std.mem.Allocator;
const alloc = std.heap.page_allocator;

const default_n = 1000;
const x_bins = 100;
const v_bins = 200;

pub const Col = sim.Col;
pub const ColType = sim.ColType;

const refresh_ns = 20000000;

pub fn main() !void {
    win = webui.newWindow();
    _ = win.show(html);

    var timer = try std.time.Timer.start();

    const params = try getArgs();
    const n = params.n;

    const p_max = sim.mu * 0.1;
    _ = p_max;

    var buf = try std.fmt.allocPrint(alloc, "setConstants({}, {});", .{
        n,
        sim.mu,
    });
    win.run(buf);

    var xs = try alloc.alloc(f64, n);
    var ps = try alloc.alloc(f64, n);
    var cols = try alloc.alloc(Col, n);
    defer {
        alloc.free(xs);
        alloc.free(ps);
        alloc.free(cols);
    }

    sim.initArrays(xs, ps, cols);

    while (true) {
        stepForward(cols, xs, ps);

        var xHist = buildHistogram(xs, 0, 1.5, x_bins);
        var xHistBytes = std.mem.asBytes(&xHist);

        var pHist = buildHistogram(ps, -math.sqrt(sim.mu) * 0.1, math.sqrt(sim.mu) * 0.1, v_bins);
        var pHistBytes = std.mem.asBytes(&pHist);

        if (timer.read() > refresh_ns) {
            timer.reset();
            win.sendRaw(
                "updateDensityHist",
                xHistBytes,
            );

            win.sendRaw(
                "updateMomentumHist",
                pHistBytes,
            );
        }
    }

    // webui.wait();

    webui.clean();
}

fn getArgs() !struct { n: u32 } {
    var arg_it = std.process.args();
    _ = arg_it.skip();
    var arg = arg_it.next() orelse return .{ .n = default_n };
    const n = try std.fmt.parseInt(u32, arg, 10);
    return .{ .n = n };
}

fn stepForward(cols: []Col, xs: []f64, ps: []f64) void {
    const next_col = sim.nextCollision(cols);
    const dt = next_col.dt;
    const j = next_col.j;
    const col_type = next_col.type;

    sim.advanceXsPs(next_col.dt, xs, ps);

    if (col_type == ColType.ground) {
        sim.computeGroundCol(j, dt, xs, ps, cols);
    } else {
        sim.computePistCol(j, dt, xs, ps, cols);
    }
}

fn buildHistogram(arr: []f64, lower: f64, upper: f64, comptime n: usize) [n]u32 {
    const step = (upper - lower) / @as(f64, n);
    var hist: [n]u32 = undefined;
    @memset(&hist, 0);

    var sum: usize = 0;

    for (arr) |x| {
        if (x < lower or x > upper) {
            continue;
        }
        // std.debug.print("x: {}", .{x});

        const i = @as(usize, @intFromFloat((x - lower) / step));
        if (i >= 0 and i < n) {
            hist[i] += 1;
            sum += 1;
        }
    }
    return hist;
}
