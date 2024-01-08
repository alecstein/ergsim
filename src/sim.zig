const std = @import("std");
const math = std.math;

const Allocator = std.mem.Allocator;
const alloc = std.heap.page_allocator;

pub const mu = 0.001;

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
    const v = 1;
    return v;
}

pub fn initArrays(xs: []f64, ps: []f64, cols: []Col) void {
    // PRNG always seeded with 0 for reproducible results

    var prng = std.rand.DefaultPrng.init(0);
    const rand = prng.random();

    var energy = piston_x;

    for (ps) |*p| {
        p.* = rand.float(f64) * 1;
        energy += p.* * p.* / mu;
    }

    for (ps) |*p| {
        p.* /= math.sqrt(energy);
    }

    piston_x /= math.sqrt(energy);

    for (xs, ps, cols) |*x, *p, *col| {
        x.* = rand.float(f32) * piston_x;
        // p.* = rand.float(f32) * 1;

        const t_p = getTimeToPiston(x.*, p.*);
        const t_g = getTimeToGround(x.*, p.*);

        energy += p.* * p.* / mu;

        const c = Col{
            .time = @min(t_p, t_g),
            .type = if (t_g < t_p) ColType.ground else ColType.piston,
        };
        col.* = c;
    }
}

pub fn getTimeToGround(x: f64, p: f64) f64 {
    if (p >= 0) return math.inf(f64);
    return -mu * x / p;
}

pub fn getTimeToPiston(x: f64, p: f64) f64 {
    const dp = piston_p - p / mu;
    const dx = piston_x - x;
    return 2 * dp + 2 * math.sqrt(dp * dp + dx);
}

// Returns the momentum of the particle after collision
pub fn elasticCol(cur_piston_p: f64, p: f64) struct { new_piston_p: f64, new_p: f64 } {
    const new_piston_p = 2 * p / (mu + 1) - (mu - 1) * cur_piston_p / (mu + 1);
    const new_p = (mu - 1) * p / (mu + 1) + 2 * mu * cur_piston_p / (mu + 1);
    return .{ .new_piston_p = new_piston_p, .new_p = new_p };
}

pub fn nextCollision(cols: []Col) struct { dt: f64, j: usize, type: ColType } {
    // "j" is the index of the colliding particle

    var dt = math.inf(f64);
    var j: usize = undefined;
    var col_type: ColType = undefined;

    for (cols, 0..) |c, k| {
        const trial_t = c.time;

        std.debug.assert(dt != 0);

        if (trial_t < dt) {
            dt = trial_t;
            j = k;
            col_type = c.type;
        }
    }

    return .{ .dt = dt, .j = j, .type = col_type };
}

pub fn advanceXsPs(dt: f64, xs: []f64, ps: []f64) void {
    piston_x += piston_p * dt - dt * dt / 4;
    piston_p += -dt / 2;

    for (xs, ps) |*x, *p| {
        x.* += p.* * dt / (mu);
    }
}

pub fn computeGroundCol(j: usize, dt: f64, xs: []f64, ps: []f64, cols: []Col) void {
    ps[j] = -ps[j];

    for (cols) |*c| {
        c.*.time -= dt;
    }

    cols[j].time = getTimeToPiston(xs[j], ps[j]);
    cols[j].type = ColType.piston;
}

pub fn computePistCol(j: usize, dt: f64, xs: []f64, ps: []f64, cols: []Col) void {
    const new_ps = elasticCol(piston_p, ps[j]);
    ps[j] = new_ps.new_p;
    piston_p = new_ps.new_piston_p;

    // A small optimization is used here.
    // We know that if a particle is going to collide with a ground,
    // some *other* particle hitting the piston won't change that.
    // Another particle hitting the piston will only slow it down. So
    // any particles that would hit the ground before this collision will
    // still hit the ground after it.

    for (cols, xs, ps) |*c, x, p| {
        if (c.*.type == ColType.ground) {
            c.*.time -= dt;
        } else if (p > 0) {
            // If the particles are traveling upwards (and not subject to gravity)
            // they'll never hit the ground
            c.*.time = getTimeToPiston(x, p);
            c.*.type = ColType.piston;
        } else {
            // This is for those cases where the particles are traveling downward,
            // but still are scheduled to get hit by the piston. These are the
            // ambiguous cases we need to solve
            const t_g = getTimeToGround(x, p);
            const t_p = getTimeToPiston(x, p);
            std.debug.assert(t_p != 0);
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
