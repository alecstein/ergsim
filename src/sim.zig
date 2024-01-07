const std = @import("std");
const math = std.math;

const Allocator = std.mem.Allocator;
const alloc = std.heap.page_allocator;

const mu = 0.001;
const T = 0.1;
pub var piston_p: f64 = 0;
pub var piston_x: f64 = 1;
pub var piston_v: f64 = 0;

pub const g = 0;
pub const m_particle = 0.0001;
pub const m_piston = 1.0;

pub const ColType = enum {
    ground,
    piston,
};

pub const Col = struct {
    time: f64,
    type: ColType,
};

pub fn initVel() f64 {
    const v = -0.00001;
    return v;
}

pub fn initArrays(xs: []f64, ps: []f64, cols: []Col) void {
    // PRNG always seeded with 0 for reproducible results

    var prng = std.rand.DefaultPrng.init(0);
    const rand = prng.random();

    for (xs, ps, cols, 0..) |*x, *p, *col, i| {
        x.* = rand.float(f32);
        p.* = initVel();

        const t_p = getTimeToPiston(x.*, p.*);
        const t_g = getTimeToGround(x.*, p.*);

        const c = Col{
            .time = @min(t_p, t_g),
            .type = if (t_g < t_p) ColType.ground else ColType.piston,
        };
        col.* = c;

        std.debug.print("{}: {}, {}, {}\n", .{ i, c, ps[i], xs[i] });
    }
}

pub fn getTimeToGround(x: f64, p: f64) f64 {
    if (p >= 0) return math.inf(f64);
    return -mu * T * x / p;
}

pub fn getTimeToPiston(x: f64, p: f64) f64 {
    const dp = piston_p - p / mu;
    const dx = piston_x - x;
    return 2 * T * dp + 2 * T * math.sqrt(dp * dp + dx);
}

// Returns the momentum of the particle after collision
pub fn elasticCol(cur_piston_p: f64, p: f64) struct { new_piston_p: f64, new_p: f64 } {
    const new_piston_p = 2 * math.sqrt(mu) * p / (1 + mu * mu) + (1 - mu * mu) * cur_piston_p / (1 + mu * mu);
    const new_p = 2 * mu * math.sqrt(mu) * cur_piston_p / (1 + mu * mu) - (1 - mu * mu) * p / (1 + mu * mu);
    return .{ .new_piston_p = new_piston_p, .new_p = new_p };
}

pub fn nextCollision(cols: []Col) struct { dt: f64, j: usize, type: ColType } {
    // "j" is the index of the colliding particle

    var dt = math.inf(f64);
    var j: usize = undefined;
    var col_type: ColType = undefined;

    for (cols, 0..) |c, k| {
        const trial_t = c.time;

        std.debug.print("{d} col: {any}\n", .{ k, c });

        std.debug.assert(dt != 0);

        if (trial_t < dt) {
            dt = trial_t;
            j = k;
            col_type = c.type;
        }
    }

    std.debug.print("dt: {}\n", .{dt});
    std.debug.print("j: {}\n", .{j});
    std.debug.print("col_type: {}\n", .{col_type});

    // std.debug.print("minimum j: {}\n", .{j});

    return .{ .dt = dt, .j = j, .type = col_type };
}

pub fn advanceXsPs(dt: f64, xs: []f64, ps: []f64) void {
    piston_x += piston_p * dt / T - dt * dt / (4 * T * T);
    piston_p += -dt / (2 * T);

    for (xs, ps) |*x, *p| {
        x.* += p.* * dt / (mu * T);
    }
}

pub fn computeGroundCol(j: usize, dt: f64, xs: []f64, ps: []f64, cols: []Col) void {
    // std.debug.print("j: {}\n", .{j});
    ps[j] = -ps[j];

    // Ground collisions don't affect the collision times of other particles
    // Every collision gets advanced forward by dt
    for (cols) |*c| {
        c.*.time -= dt;
    }

    cols[j].time = getTimeToPiston(xs[j], ps[j]);
    cols[j].type = ColType.piston;
}

pub fn computePistCol(j: usize, dt: f64, xs: []f64, ps: []f64, cols: []Col) void {
    std.debug.print("j: {}\n", .{j});
    // print out old and new ps to check for errors
    std.debug.print("old ps[j]: {}\n", .{ps[j]});
    const new_ps = elasticCol(piston_p, ps[j]);
    ps[j] = new_ps.new_p;
    piston_p = new_ps.new_piston_p;
    std.debug.print("new ps[j]: {}\n", .{ps[j]});

    // A small optimization is used here.
    // We know that if a particle is going to collide with a ground,
    // some *other* particle hitting the piston won't change that.
    // Another particle hitting the piston will only slow it down. So
    // any particles that would hit the ground before this collision will
    // still hit the ground after it.

    for (cols, xs, ps) |*c, x, v| {
        if (c.*.type == ColType.ground) {
            c.*.time -= dt;
        } else if (v > 0) {
            // If the particles are traveling upwards (and not subject to gravity)
            // they'll never hit the ground
            c.*.time = getTimeToPiston(x, v);
            c.*.type = ColType.piston;
        } else {
            // This is for those cases where the particles are traveling downward,
            // but still are scheduled to get hit by the piston. These are the
            // ambiguous cases we need to solve
            const t_g = getTimeToGround(x, v);
            const t_p = getTimeToPiston(x, v);
            c.*.time = @min(t_g, t_p);
            c.*.type = if (t_g < t_p) ColType.ground else ColType.piston;
        }
    }
}

pub fn getTimeToGroundWithGravity(x: f64, v: f64) f64 {
    const a = -0.5 * g;
    const b = v;
    const c = x;
    const disc = b * b - 4 * a * c;
    const t = (-b - math.sqrt(disc)) / (2 * a);
    return t;
}

pub fn getTimeToPistonWithGravity(x: f64, v: f64) f64 {
    if (piston_v - v >= 0) return math.inf(f64);
    const t = (x - piston_x) / (piston_v - v);
    std.debug.assert(t >= 0);
    return t;
}

pub fn advanceXsPsWithGravity(dt: f64, xs: []f64, ps: []f64) void {
    piston_x += piston_v * dt - g * dt * dt / 2;
    piston_v -= g * dt;

    for (xs, ps) |*x, *v| {
        x.* += v.* * dt - g * dt * dt / 2;
        v.* -= g * dt;
    }
}

pub fn computeGroundColWithGravity(j: usize, dt: f64, xs: []f64, ps: []f64, cols: []Col) void {
    ps[j] = -ps[j];

    // Ground collisions don't affect the collision times of other particles
    // Every collision gets advanced forward by dt
    for (cols) |*c| {
        c.*.time -= dt;
    }

    const t_g = getTimeToGroundWithGravity(xs[j], ps[j]);
    const t_p = getTimeToPistonWithGravity(xs[j], ps[j]);
    cols[j].time = @min(t_g, t_p);
    cols[j].type = if (t_g < t_p) ColType.ground else ColType.piston;
}

pub fn computePistColWithGravity(j: usize, dt: f64, xs: []f64, ps: []f64, cols: []Col) void {
    const new_ps = elasticCol(piston_p, ps[j]);
    ps[j] = new_ps.new_p;
    piston_p = new_ps.new_piston_p;

    // A small optimization is used here.
    // We know that if a particle is going to collide with a ground,
    // some *other* particle hitting the piston won't change that.
    // Another particle hitting the piston will only slow it down. So
    // any particles that would hit the ground before this collision will
    // still hit the ground after it.

    for (cols, xs, ps) |*c, x, v| {
        if (c.*.type == ColType.ground) {
            c.*.time -= dt;
        } else if (piston_v - v >= 0) {
            // This particle will never hit the piston
            c.*.time = getTimeToGroundWithGravity(x, v);
            c.*.type = ColType.ground;
        } else {
            const t_g = getTimeToGroundWithGravity(x, v);
            const t_p = getTimeToPistonWithGravity(x, v);
            c.*.time = @min(t_g, t_p);
            c.*.type = if (t_g < t_p) ColType.ground else ColType.piston;
        }
    }
}
