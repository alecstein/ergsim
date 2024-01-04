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

    var buf = try std.fmt.allocPrint(alloc, "setConstants({}, {}, {}, {}, {});", .{ n, sim.m_particle, sim.m_piston, sim.piston_x, 1 });
    win.run(buf);

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

    while (true) {
        const dt = stepForward(cols, xs, vs);
        t += dt;
        ct += 1;

        const xsBytes = std.mem.sliceAsBytes(xs);
        const vsBytes = std.mem.sliceAsBytes(vs);

        if (ct % (n / 50) == 0) {
            win.sendRaw(
                "updateGasDensityHistogram",
                xsBytes,
            );

            win.sendRaw(
                "updateMomentumHistogram",
                vsBytes,
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
