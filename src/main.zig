const std = @import("std");
const math = std.math;
const Allocator = std.mem.Allocator;
const global_allocator = std.heap.page_allocator;

const g = 9.81;
const m_particle = 0.001;
const m_piston = 1.0;
const default_n = 1000;
const default_t = 1000.0;
const gravity = false;
var piston_x: f64 = 1;
var piston_v: f64 = 0;

pub fn main() !void {
    var timer = try std.time.Timer.start();
    const t0 = timer.read();

    const params = try getArgs();
    const n = params.n;
    const max_time = params.t;

    std.debug.print("Particles: {d}\tWorldtime: {d}s\n", .{ n, max_time });

    var xs = try global_allocator.alloc(f64, n);
    var vs = try global_allocator.alloc(f64, n);
    var t_grounds = try global_allocator.alloc(f64, n);
    var t_pistons = try global_allocator.alloc(f64, n);
    defer {
        global_allocator.free(xs);
        global_allocator.free(vs);
        global_allocator.free(t_grounds);
        global_allocator.free(t_pistons);
    }

    initializeArrays(&xs, &vs, &t_grounds, &t_pistons, n);

    const est_total_items = @as(usize, @intFromFloat(@divFloor(max_time, 10)));
    var pct_done: u8 = 0;
    var progress = std.Progress{};
    const root_node = progress.start("Simulating", est_total_items);
    defer root_node.end();

    var int = getNextInteraction(t_grounds, t_pistons, n);
    var dt = int.dt;
    var j = int.j;
    var is_ground_col = int.is_ground_col;

    var t: f64 = 0;
    var ct: usize = 0; // collision count

    while (t < max_time) {
        t += dt;
        ct += 1;

        try stepForward(&xs, &vs, &t_grounds, n, dt);

        if (is_ground_col) {
            try handleGroundCollision(&xs, &vs, &t_grounds, &t_pistons, n, j, dt);
        } else {
            try handlePistonCollision(&xs, &vs, &t_grounds, &t_pistons, n, j);
        }

        int = getNextInteraction(t_grounds, t_pistons, n);
        dt = int.dt;
        j = int.j;
        is_ground_col = int.is_ground_col;

        updateProgress(root_node, &pct_done, t, max_time, est_total_items);
    }

    const t1 = timer.read();
    const time_taken = @as(f64, @floatFromInt(t1 - t0));

    std.debug.print("\nCollisions: {}\n", .{ct});
    std.debug.print("Time taken: {d:.2}s\n", .{time_taken / 1000000000.0});
}

fn initializeArrays(xs: *[]f64, vs: *[]f64, t_grounds: *[]f64, t_pistons: *[]f64, n: usize) void {
    // without a seed, this will always produce the same sequence
    // this is good for debugging, but a more robust
    // solution might be better long-term.
    var prng = std.rand.DefaultPrng.init(0);
    const rand = prng.random();

    for (0..n) |i| {
        xs.*[i] = piston_x * rand.float(f32);
        vs.*[i] = initVel();
        t_grounds.*[i] = getTimeToGround(xs.*[i], vs.*[i]);
        t_pistons.*[i] = getTimeToPiston(xs.*[i], vs.*[i]);
    }
}

// get time to hit ground
inline fn getTimeToGround(x: f64, v: f64) f64 {
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

// get time to collide with piston
fn getTimeToPiston(x: f64, v: f64) f64 {
    if (!gravity) {
        const a = -0.5 * g;
        const b = piston_v - v;
        const c = piston_x - x;
        const disc = b * b - 4 * a * c;
        return (-b - math.sqrt(disc)) / (2 * a);
    }
    const t = (x - piston_x) / (piston_v - v);
    if (piston_v - v >= 0) return math.inf(f64);
    if (t == 0) return math.inf(f64);
    return t;
}

// compute velocities after elastic collision
fn elasticCol(m1: f64, v1: f64, m2: f64, v2: f64) struct { v1_prime: f64, v2_prime: f64 } {
    const v1_prime = (m1 - m2) * v1 / (m1 + m2) + 2 * m2 * v2 / (m1 + m2);
    const v2_prime = 2 * m1 * v1 / (m1 + m2) - (m1 - m2) * v2 / (m1 + m2);
    return .{ .v1_prime = v1_prime, .v2_prime = v2_prime };
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

fn getNextInteraction(t_grounds: []f64, t_pistons: []f64, n: usize) struct { dt: f64, j: usize, is_ground_col: bool } {
    // this is always O(N) time, so there's no prettier way to do it
    var dt = math.inf(f64);
    var j: usize = undefined;
    var is_ground_col: bool = true;

    for (0..n) |i| {
        const trial_t = @min(t_grounds[i], t_pistons[i]);
        if (trial_t < dt) {
            dt = trial_t;
            j = i;
            is_ground_col = t_grounds[i] < t_pistons[i];
        }
    }

    return .{ .dt = dt, .j = j, .is_ground_col = is_ground_col };
}

fn initVel() f64 {
    //     const vAvg: f64 = 0.5 * math.sqrt(2 * m_piston * g * piston_x / m_particle / @as(f64, @floatFromInt(n)));
    const v = 1.0;
    return v;
}

fn stepForward(xs: *[]f64, vs: *[]f64, t_grounds: *[]f64, n: usize, dt: f64) !void {
    for (0..n) |i| {
        if (!gravity) {
            xs.*[i] += vs.*[i] * dt;
        } else {
            xs.*[i] += vs.*[i] * dt - g * dt * dt / 2;
            vs.*[i] -= g * dt;
        }
        t_grounds.*[i] -= dt;
    }
    piston_x += piston_v * dt - g * dt * dt / 2;
    piston_v -= g * dt;
}

fn handleGroundCollision(xs: *[]f64, vs: *[]f64, t_grounds: *[]f64, t_pistons: *[]f64, n: usize, j: usize, dt: f64) !void {
    vs.*[j] = -vs.*[j];
    if (!gravity) {
        t_grounds.*[j] = math.inf(f64);
    } else {
        t_grounds.*[j] = getTimeToGround(xs.*[j], vs.*[j]);
    }

    for (0..n) |i| {
        t_pistons.*[i] -= dt;
    }

    t_pistons.*[j] = getTimeToPiston(xs.*[j], vs.*[j]);
}

fn handlePistonCollision(xs: *[]f64, vs: *[]f64, t_grounds: *[]f64, t_pistons: *[]f64, n: usize, j: usize) !void {
    const new_vs = elasticCol(m_particle, vs.*[j], m_piston, piston_v);
    vs.*[j] = new_vs.v1_prime;
    piston_v = new_vs.v2_prime;
    t_grounds.*[j] = getTimeToGround(xs.*[j], vs.*[j]);
    for (0..n) |i| {
        t_pistons.*[i] = getTimeToPiston(xs.*[i], vs.*[i]);
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
