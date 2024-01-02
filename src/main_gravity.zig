// Simulates the case where the particles are operating under the
// same gravitational field as the piston

const std = @import("std");
const math = std.math;
const Allocator = std.mem.Allocator;
const alloc = std.heap.page_allocator;

const g = 9.81;
const m_particle = 0.001;
const m_piston = 1.0;
const default_n = 1000;
const default_t = 1000.0;
var piston_x: f64 = 1;
var piston_v: f64 = 0;

const ColType = enum {
    ground,
    piston,
};

const Col = struct {
    time: f64,
    type: ColType,
};

pub fn main() !void {
    var timer = try std.time.Timer.start();
    const t0 = timer.read();

    const params = try getArgs();
    const n = params.n;
    const max_time = params.t;

    std.debug.print("Particles: {d}\tWorldtime: {d}s\n", .{ n, max_time });

    var xs = try alloc.alloc(f64, n);
    var vs = try alloc.alloc(f64, n);
    var cols = try alloc.alloc(Col, n);
    defer {
        alloc.free(xs);
        alloc.free(vs);
        alloc.free(cols);
    }

    initArrays(xs, vs, cols);

    const est_total_items = @as(usize, @intFromFloat(@divFloor(max_time, 10)));
    var pct_done: u8 = 0;
    var progress = std.Progress{};
    const root_node = progress.start("Simulating", est_total_items);
    defer root_node.end();

    var t: f64 = 0;
    var ct: usize = 0; // collision count

    while (t < max_time) {
        const next_col = nextCollision(cols);
        const dt = next_col.dt;
        const j = next_col.j;
        const col_type = next_col.type;

        advanceXsVs(next_col.dt, xs, vs);

        if (col_type == ColType.ground) {
            computeGroundCol(j, dt, xs, vs, cols);
        } else {
            computePistCol(j, dt, xs, vs, cols);
        }

        updateProgress(root_node, &pct_done, t, max_time, est_total_items);

        t += dt;
        ct += 1;
    }

    const time_taken = @as(f64, @floatFromInt(timer.read() - t0));
    std.debug.print("\nCollisions: {}\n", .{ct});
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

fn initVel() f64 {
    const v = 1.0;
    return v;
}

fn initArrays(xs: []f64, vs: []f64, cols: []Col) void {
    // PRNG always seeded with 0 for reproducible results
    var prng = std.rand.DefaultPrng.init(0);
    const rand = prng.random();

    for (xs, vs, cols) |*x, *v, *col| {
        x.* = piston_x * rand.float(f32);
        v.* = initVel();

        const t_p = getTimeToPiston(x.*, v.*);
        const t_g = getTimeToGround(x.*, v.*);

        const c = Col{
            .time = @min(t_p, t_g),
            .type = if (t_g < t_p) ColType.ground else ColType.piston,
        };
        col.* = c;
    }
}

fn getTimeToGround(x: f64, v: f64) f64 {
    const a = -0.5 * g;
    const b = v;
    const c = x;
    const disc = b * b - 4 * a * c;
    const t = (-b - math.sqrt(disc)) / (2 * a);
    return t;
}

fn getTimeToPiston(x: f64, v: f64) f64 {
    if (piston_v - v >= 0) return math.inf(f64);
    const t = (x - piston_x) / (piston_v - v);
    std.debug.assert(t >= 0);
    return t;
}

fn elasticCol(m1: f64, v1: f64, m2: f64, v2: f64) struct { v1_prime: f64, v2_prime: f64 } {
    const v1_prime = (m1 - m2) * v1 / (m1 + m2) + 2 * m2 * v2 / (m1 + m2);
    const v2_prime = 2 * m1 * v1 / (m1 + m2) - (m1 - m2) * v2 / (m1 + m2);
    return .{ .v1_prime = v1_prime, .v2_prime = v2_prime };
}

fn nextCollision(cols: []Col) struct { dt: f64, j: usize, type: ColType } {
    // "j" is the index of the colliding particle

    var dt = math.inf(f64);
    var j: usize = undefined;
    var col_type: ColType = undefined;

    for (cols, 0..) |c, k| {
        const trial_t = c.time;
        if (trial_t < dt) {
            dt = trial_t;
            j = k;
            col_type = c.type;
        }
    }

    return .{ .dt = dt, .j = j, .type = col_type };
}

fn advanceXsVs(dt: f64, xs: []f64, vs: []f64) void {
    piston_x += piston_v * dt - g * dt * dt / 2;
    piston_v -= g * dt;

    for (xs, vs) |*x, *v| {
        x.* += v.* * dt - g * dt * dt / 2;
        v.* -= g * dt;
    }
}

fn computeGroundCol(j: usize, dt: f64, xs: []f64, vs: []f64, cols: []Col) void {
    vs[j] = -vs[j];

    // Every collision gets advanced dt in time
    for (cols) |*c| {
        c.*.time -= dt;
    }

    const t_g = getTimeToGround(xs[j], vs[j]);
    const t_p = getTimeToPiston(xs[j], vs[j]);
    cols[j].time = @min(t_g, t_p);
    cols[j].type = if (t_g < t_p) ColType.ground else ColType.piston;
}

fn computePistCol(j: usize, dt: f64, xs: []f64, vs: []f64, cols: []Col) void {
    const new_vs = elasticCol(m_particle, vs[j], m_piston, piston_v);
    vs[j] = new_vs.v1_prime;
    piston_v = new_vs.v2_prime;

    // A small optimization is used here.
    // We know that if a particle is going to collide with a ground,
    // some *other* particle hitting the piston won't change that.
    // Another particle hitting the piston will only slow it down. So
    // any particles that would hit the ground before this collision will
    // still hit the ground after it.

    for (cols, xs, vs) |*c, x, v| {
        if (c.*.type == ColType.ground) {
            c.*.time -= dt;
        } else {
            const t_g = getTimeToGround(x, v);
            const t_p = getTimeToPiston(x, v);
            c.*.time = @min(t_g, t_p);
            c.*.type = if (t_g < t_p) ColType.ground else ColType.piston;
        }
    }
}

fn updateProgress(root_node: *std.Progress.Node, pct_done: *u8, t: f64, max_time: f64, est_total_items: usize) void {
    const frac_time = (t * @as(f64, @floatFromInt(est_total_items))) / max_time;
    const int_frac_time = @as(u8, @intFromFloat(frac_time));

    const new_pct_done = @as(u8, int_frac_time);
    if (new_pct_done > pct_done.*) {
        pct_done.* = new_pct_done;
        root_node.completeOne();
    }
}
