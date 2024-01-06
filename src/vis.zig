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

pub const Col = sim.Col;
pub const ColType = sim.ColType;

pub fn main() !void {
    win = webui.newWindow();
    _ = win.show(html);

    const params = try getArgs();
    const n = params.n;
    const x_bins = 100;
    const v_bins = 100;

    var xs = try alloc.alloc(f64, n);
    var vs = try alloc.alloc(f64, n);
    var cols = try alloc.alloc(Col, n);
    defer {
        alloc.free(xs);
        alloc.free(vs);
        alloc.free(cols);
    }

    sim.initArrays(xs, vs, cols);

    var init_ke_particles: f64 = 0;
    for (vs) |v| {
        init_ke_particles += sim.m_particle * math.pow(f64, v, 2) / 2;
    }
    const init_ke_piston = sim.m_piston * math.pow(f64, sim.piston_v, 2) / 2;
    const init_pe_particles: f64 = 0;
    const init_pe_piston = sim.g * sim.m_piston * sim.piston_x;
    const energy = init_ke_particles + init_ke_piston + init_pe_particles + init_pe_piston;

    var t: f64 = 0;
    var ct: usize = 0; // collision count

    // get real system time
    var timer = try std.time.Timer.start();
    var t_ct = timer.read();

    var buf = try std.fmt.allocPrint(alloc, "setConstants({}, {}, {}, {}, {}, {}, {});", .{ n, sim.m_particle, sim.m_piston, sim.piston_x, x_bins, v_bins, energy });
    win.run(buf);

    while (true) {
        const dt = stepForward(cols, xs, vs);
        t += dt;
        ct += 1;

        const xsBytes = std.mem.sliceAsBytes(xs);
        _ = xsBytes;
        const vsBytes = std.mem.sliceAsBytes(vs);
        _ = vsBytes;

        var xHist = buildHistogram(xs, 0, 2, x_bins);
        // std.debug.print("xHist: {any}\n", .{xHist});
        var xHistBytes = std.mem.asBytes(&xHist);
        // std.debug.print("xHistBytes: {any}\n", .{xHistBytes});

        var vHist = buildHistogram(vs, -5, 5, v_bins);
        var vHistBytes = std.mem.asBytes(&vHist);

        if (timer.read() > t_ct + 50000000) {
            t_ct = timer.read();
            win.sendRaw(
                "updateDensityHist",
                xHistBytes,
            );

            win.sendRaw(
                "updateMomentumHist",
                vHistBytes,
            );
        }
    }

    webui.wait();

    webui.clean();
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
