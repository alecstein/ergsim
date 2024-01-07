const std = @import("std");
const math = std.math;

const Allocator = std.mem.Allocator;
const alloc = std.heap.page_allocator;

pub const g = 9.81;
pub const m_particle = 0.001;
pub const m_piston = 1.0;
pub var piston_x: f64 = 1;
pub var piston_v: f64 = 0;

pub const ColType = enum {
    ground,
    piston,
};

pub const Col = struct {
    time: f64,
    type: ColType,
};

pub fn initVel() f64 {
    const v = 1.0;
    return v;
}

pub fn initArrays(xs: []f64, vs: []f64, cols: []Col) void {
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

pub fn getTimeToGround(x: f64, v: f64) f64 {
    if (v >= 0) return math.inf(f64);
    return -x / v;
}

pub fn getTimeToPiston(x: f64, v: f64) f64 {
    const a = -0.5 * g;
    const b = piston_v - v;
    const c = piston_x - x;
    const disc = b * b - 4 * a * c;
    return (-b - math.sqrt(disc)) / (2 * a);
}

pub fn elasticCol(m1: f64, v1: f64, m2: f64, v2: f64) struct { v1_prime: f64, v2_prime: f64 } {
    const v1_prime = (m1 - m2) * v1 / (m1 + m2) + 2 * m2 * v2 / (m1 + m2);
    const v2_prime = 2 * m1 * v1 / (m1 + m2) - (m1 - m2) * v2 / (m1 + m2);
    return .{ .v1_prime = v1_prime, .v2_prime = v2_prime };
}

pub fn nextCollision(cols: []Col) struct { dt: f64, j: usize, type: ColType } {
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

pub fn advanceXsVs(dt: f64, xs: []f64, vs: []f64) void {
    piston_x += piston_v * dt - g * dt * dt / 2;
    piston_v -= g * dt;

    for (xs, vs) |*x, *v| {
        x.* += v.* * dt;
    }
}

pub fn computeGroundCol(j: usize, dt: f64, xs: []f64, vs: []f64, cols: []Col) void {
    vs[j] = -vs[j];

    // Ground collisions don't affect the collision times of other particles
    // Every collision gets advanced forward by dt
    for (cols) |*c| {
        c.*.time -= dt;
    }

    cols[j].time = getTimeToPiston(xs[j], vs[j]);
    cols[j].type = ColType.piston;
}

pub fn computePistCol(j: usize, dt: f64, xs: []f64, vs: []f64, cols: []Col) void {
    const new_vs = elasticCol(m_particle, vs[j], m_piston, piston_v);
    vs[j] = new_vs.v1_prime;
    piston_v = new_vs.v2_prime;

    // A small optimization is used here.
    // We know that if a particle is going to collide with a ground,
    // some *other* particle hitting the piston won't change that.
    // Another particle hitting the piston will only slow it down. So
    // any particles that would hit the ground before this collision will
    // still hit the ground after it.
    std.debug.assert(new_vs.v2_prime >= piston_v);

    for (cols, xs, vs) |*c, x, v| {
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

pub fn advanceXsVsWithGravity(dt: f64, xs: []f64, vs: []f64) void {
    piston_x += piston_v * dt - g * dt * dt / 2;
    piston_v -= g * dt;

    for (xs, vs) |*x, *v| {
        x.* += v.* * dt - g * dt * dt / 2;
        v.* -= g * dt;
    }
}

pub fn computeGroundColWithGravity(j: usize, dt: f64, xs: []f64, vs: []f64, cols: []Col) void {
    vs[j] = -vs[j];

    // Ground collisions don't affect the collision times of other particles
    // Every collision gets advanced forward by dt
    for (cols) |*c| {
        c.*.time -= dt;
    }

    const t_g = getTimeToGroundWithGravity(xs[j], vs[j]);
    const t_p = getTimeToPistonWithGravity(xs[j], vs[j]);
    cols[j].time = @min(t_g, t_p);
    cols[j].type = if (t_g < t_p) ColType.ground else ColType.piston;
}

pub fn computePistColWithGravity(j: usize, dt: f64, xs: []f64, vs: []f64, cols: []Col) void {
    const new_vs = elasticCol(m_particle, vs[j], m_piston, piston_v);
    vs[j] = new_vs.v1_prime;
    piston_v = new_vs.v2_prime;

    // A small optimization is used here.
    // We know that if a particle is going to collide with a ground,
    // some *other* particle hitting the piston won't change that.
    // Another particle hitting the piston will only slow it down. So
    // any particles that would hit the ground before this collision will
    // still hit the ground after it.
    std.debug.assert(new_vs.v2_prime >= piston_v);

    for (cols, xs, vs) |*c, x, v| {
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
