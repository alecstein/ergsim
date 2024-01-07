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
const v_bins = 100;

pub const Col = sim.Col;
pub const ColType = sim.ColType;

pub fn main() !void {
    win = webui.newWindow();
    _ = win.show(html);

    var timer = std.time.Timer();

    const params = try getArgs();
    const n = params.n;

    var buf = try std.fmt.allocPrint(alloc, "setConstants({}, {}, {}, {}, {});", .{ n, sim.m_particle, sim.m_piston, sim.piston_x, 1 });
    win.run(buf);

    std.debug.print("Particles: {d}\n", .{n});

    var xs = try alloc.alloc(f64, n);
    var ps = try alloc.alloc(f64, n);
    var cols = try alloc.alloc(Col, n);
    defer {
        alloc.free(xs);
        alloc.free(ps);
        alloc.free(cols);
    }

    sim.initArrays(xs, ps, cols);

    var t: f64 = 0;
    var ct: usize = 0; // collision count
    var t_ct = timer.read();

    while (true) {
        const dt = stepForward(cols, xs, ps);
        t += dt;
        ct += 1;

        const xsBytes = std.mem.sliceAsBytes(xs);
        _ = xsBytes;
        const psBytes = std.mem.sliceAsBytes(ps);
        _ = psBytes;

        var xHist = buildHistogram(xs, 0, 1, x_bins);
        // std.debug.print("xHist: {any}\n", .{xHist});
        var xHistBytes = std.mem.asBytes(&xHist);
        // std.debug.print("xHistBytes: {any}\n", .{xHistBytes});

        var pHist = buildHistogram(ps, -sim.init_v * 4, sim.init_v * 4, v_bins);
        var pHistBytes = std.mem.asBytes(&pHist);

        if (timer.read() > t_ct + 50000000) {
            t_ct = timer.read();
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

fn stepForward(cols: []Col, xs: []f64, ps: []f64) f64 {
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
    return dt;
}

fn buildHistogram(arr: []f64, lower: f64, upper: f64, comptime n: usize) [n]u32 {
    const step = (upper - lower) / @as(f64, n);
    var hist: [n]u32 = undefined;
    @memset(&hist, 0);

    for (arr) |x| {
        if (x < lower or x > upper) {
            continue;
        }
        // std.debug.print("x: {}", .{x});

        const i = @as(usize, @intFromFloat((x - lower) / step));
        if (i >= 0 and i < n) {
            hist[i] += 1;
        }
    }
    return hist;
}
