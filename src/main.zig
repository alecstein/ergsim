const std = @import("std");
const math = std.math;
const Allocator = std.mem.Allocator;
const alloc = std.heap.page_allocator;

// const tracer = @import("tracer");
// pub const tracer_impl = tracer.spall;

const g = 9.81;
const m_particle = 0.001;
const m_piston = 1.0;
const default_n = 1000;
const default_t = 1000.0;
const gravity = true;
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
    // try tracer.init();
    // defer tracer.deinit();

    // try tracer.init_thread(std.fs.cwd());
    // defer tracer.deinit_thread();

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

    // var filename: []const u8 = try std.fmt.allocPrint(
    //     alloc,
    //     "N_{d}_Time_{d}_x.txt",
    //     .{ n, max_time },
    // );

    // const file = try std.fs.cwd().createFile(
    //     filename,
    //     .{ .read = true },
    // );
    // var writer = file.writer();
    // defer file.close();

    while (t < max_time) {
        ct += 1;

        const int = getNextInteraction(cols);

        // advance the particles forward
        advanceSystem(int.dt, xs, vs);

        if (int.col_type == ColType.ground) {
            computeGroundCol(int.j, int.dt, xs, vs, cols);
        } else {
            computePistCol(int.j, int.dt, xs, vs, cols);
        }

        updateProgress(root_node, &pct_done, t, max_time, est_total_items);

        t += int.dt;

        // if (ct % 10000 == 0) {
        //     for (vs) |v| {
        //         try std.fmt.format(writer, "{}\n", .{v});
        //     }
        // }
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
    // const vAvg: f64 = 0.5 * math.sqrt(2 * m_piston * g * piston_x / m_particle / @as(f64, @floatFromInt(n)));
    const v = 1.0;
    return v;
}

fn initArrays(xs: []f64, vs: []f64, cols: []Col) void {
    // without a seed, this will always produce the same sequence
    // this is good for debugging, but a more robust
    // solution might be better long-term.
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
    // const t_ = tracer.trace(@src(), "getTimeToGround", .{});
    // defer t_.end();

    @setFloatMode(.Optimized);

    if (!gravity) {
        if (v >= 0) return math.inf(f64);
        return -x / v;
    }
    const a = -0.5 * g;
    const b = v;
    const c = x;
    const disc = b * b - 4 * a * c;
    const t = (-b - math.sqrt(disc)) / (2 * a);
    return t;
}

fn getTimeToPiston(x: f64, v: f64) f64 {
    // const t_ = tracer.trace(@src(), "getTimeToPiston", .{});
    // defer t_.end();

    @setFloatMode(.Optimized);

    if (!gravity) {
        const a = -0.5 * g;
        const b = piston_v - v;
        const c = piston_x - x;
        const disc = b * b - 4 * a * c;
        return (-b - math.sqrt(disc)) / (2 * a);
    }

    if (piston_v - v >= 0) return math.inf(f64);
    const t = (x - piston_x) / (piston_v - v);
    std.debug.assert(t >= 0);
    return t;
}

fn elasticCol(m1: f64, v1: f64, m2: f64, v2: f64) struct { v1_prime: f64, v2_prime: f64 } {
    // const t_ = tracer.trace(@src(), "elasticCol", .{});
    // defer t_.end();

    const v1_prime = (m1 - m2) * v1 / (m1 + m2) + 2 * m2 * v2 / (m1 + m2);
    const v2_prime = 2 * m1 * v1 / (m1 + m2) - (m1 - m2) * v2 / (m1 + m2);
    return .{ .v1_prime = v1_prime, .v2_prime = v2_prime };
}

fn getNextInteraction(cols: []Col) struct { dt: f64, j: usize, col_type: ColType } {
    // const t_ = tracer.trace(@src(), "getNext", .{});
    // defer t_.end();

    // this is always O(N) time, so there's no prettier way to do it
    var dt = math.inf(f64);
    var j: usize = undefined;
    var col_type: ColType = undefined;

    for (cols, 0..) |c, i| {
        const trial_t = c.time;
        if (trial_t < dt) {
            dt = trial_t;
            j = i;
            col_type = c.type;
        }
    }

    return .{ .dt = dt, .j = j, .col_type = col_type };
}

fn advanceSystem(dt: f64, xs: []f64, vs: []f64) void {
    // const t_ = tracer.trace(@src(), "advanceSys", .{});
    // defer t_.end();
    piston_x += piston_v * dt - g * dt * dt / 2;
    piston_v -= g * dt;

    for (xs, vs) |*x, *v| {
        if (!gravity) {
            x.* += v.* * dt;
        } else {
            x.* += v.* * dt - g * dt * dt / 2;
            v.* -= g * dt;
        }
    }
}

fn computeGroundCol(j: usize, dt: f64, xs: []f64, vs: []f64, cols: []Col) void {
    // const t_ = tracer.trace(@src(), "computeGroundCol", .{});
    // defer t_.end();

    // reverse the velocity of the colliding particle
    vs[j] = -vs[j];

    // every collision gets dt closer
    for (cols) |*c| {
        c.*.time -= dt;
    }

    // update the next interaction for the particle that just collided with the ground
    if (!gravity) {
        // without gravity, the next collision is _always_ with the piston
        cols[j].time = math.inf(f64);
        cols[j].type = ColType.piston;
    } else {
        // under gravity it could be either
        const t_g = getTimeToGround(xs[j], vs[j]);
        const t_p = getTimeToPiston(xs[j], vs[j]);
        cols[j].time = @min(t_g, t_p);
        cols[j].type = if (t_g < t_p) ColType.ground else ColType.piston;
    }
}

fn computePistCol(j: usize, dt: f64, xs: []f64, vs: []f64, cols: []Col) void {
    // const t_ = tracer.trace(@src(), "computePistCol", .{});
    // defer t_.end();

    // if a particle is going to hit the ground, there's no need to recompute its
    // collision time with the piston. A piston colliding with any particle is not
    // going to make the piston move downward faster.

    const new_vs = elasticCol(m_particle, vs[j], m_piston, piston_v);
    vs[j] = new_vs.v1_prime;
    piston_v = new_vs.v2_prime;

    for (cols, 0..) |*c, i| {
        if (c.*.type == ColType.ground) {
            c.*.time -= dt;
        } else {
            const t_g = getTimeToGround(xs[i], vs[i]);
            const t_p = getTimeToPiston(xs[i], vs[i]);
            c.*.time = @min(t_g, t_p);
            c.*.type = if (t_g < t_p) ColType.ground else ColType.piston;
        }
    }
}

fn updateProgress(root_node: *std.Progress.Node, pct_done: *u8, t: f64, max_time: f64, est_total_items: usize) void {
    // const t_ = tracer.trace(@src(), "updateProgress", .{});
    // defer t_.end();

    const frac_time = (t * @as(f64, @floatFromInt(est_total_items))) / max_time;
    const int_frac_time = @as(u8, @intFromFloat(frac_time));

    const new_pct_done = @as(u8, int_frac_time);
    if (new_pct_done > pct_done.*) {
        pct_done.* = new_pct_done;
        root_node.completeOne();
    }
}
