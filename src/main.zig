const std = @import("std");
const math = std.math;

const g = 9.81;
const m_particle = 0.001;
const m_piston = 1.0;

var piston_x: f64 = 1;
var piston_v: f64 = 0;

const Particle = struct {
    x: f64,
    v: f64,
    t_ground: f64,
    t_piston: f64,
};

const ParticleList = std.MultiArrayList(Particle);

// get time to hit ground
fn tGround(x: f64, v: f64) f64 {
    if (v >= 0) return math.inf(f64);
    return -x / v;
}

// get time to collide with piston
fn tPiston(x: f64, v: f64) f64 {
    const a = 0.5 * g;
    const b = v - piston_v;
    const c = x - piston_x;
    const disc = b * b - 4 * a * c;
    return (-b + math.sqrt(disc)) / (2 * a);
}

// compute velocities after elastic collision
fn elasticCol(m1: f64, v1: f64, m2: f64, v2: f64) struct { v1_prime: f64, v2_prime: f64 } {
    const v1_prime = (m1 - m2) * v1 / (m1 + m2) + 2 * m2 * v2 / (m1 + m2);
    const v2_prime = 2 * m1 * v1 / (m1 + m2) - (m1 - m2) * v2 / (m1 + m2);
    return .{ .v1_prime = v1_prime, .v2_prime = v2_prime };
}

pub fn main() !void {
    const allocator = std.heap.page_allocator;
    const args = try std.process.argsAlloc(allocator);

    var timer = try std.time.Timer.start();
    const t0 = timer.read();

    if (args.len < 3) return error.ExpectedArgument;
    const n = try std.fmt.parseInt(u32, args[1], 0);
    const max_time = try std.fmt.parseFloat(f64, args[2]);

    var particles = ParticleList{};
    defer particles.deinit(allocator);

    const vAvg: f64 = 0.5 * math.sqrt(2 * m_piston * g * piston_x / m_particle / @as(f64, @floatFromInt(n)));
    _ = vAvg;

    var prng = std.rand.DefaultPrng.init(0);
    const r = prng.random();

    for (0..n) |i| {
        _ = i;
        // TODO write this as init method
        var p = Particle{
            .x = piston_x * r.float(f64),
            .v = 1,
            .t_ground = undefined,
            .t_piston = undefined,
        };
        p.t_ground = tGround(p.x, p.v);
        p.t_piston = tPiston(p.x, p.v);
        try particles.append(allocator, p);
    }

    var ps = particles.slice();
    var xs = ps.items(.x);
    var vs = ps.items(.v);
    var t_gs = ps.items(.t_ground);
    var t_ps = ps.items(.t_piston);

    // progress indicator
    var pct_done: u8 = 0;
    var progress = std.Progress{};
    const root_node = progress.start("Simulating", 100);
    defer root_node.end();

    var t: f64 = 0;
    var ct: usize = 0; // collision count

    // std.debug.print("Particles: {d}\tWorldtime: {d}s\n", .{ n, max_time });

    while (t < max_time) {
        ct += 1;

        // update progress indicator
        const new_pct_done = @as(u8, @intFromFloat((t * 100) / max_time));
        if (new_pct_done > pct_done) {
            pct_done = new_pct_done;
            root_node.completeOne();
            progress.refresh();
        }

        // get the next interaction
        // this is always O(N) time, so there's no prettier way to do it
        var dt = math.inf(f64);
        var j: usize = undefined;
        var is_ground_col: bool = true;

        for (t_gs, t_ps, 0..n) |t_g, t_p, i| {
            const trial_t = @min(t_g, t_p);
            if (trial_t < dt) {
                dt = trial_t;
                j = i;
                is_ground_col = t_g < t_p;
            }
        }

        // step forward in time
        piston_x += piston_v * dt - g * dt * dt / 2;
        piston_v -= g * dt;
        t += dt;

        for (xs, vs, t_gs) |*x, v, *t_g| {
            x.* += v * dt;
            t_g.* -= dt;
        }

        var pj = ps.get(j);
        if (is_ground_col) {
            // handle ground collision
            pj.v = -pj.v;
            pj.t_ground = math.inf(f64);
            for (t_ps) |*t_p| {
                t_p.* -= dt;
            }
        } else {

            // handle piston collision
            const new_vs = elasticCol(m_particle, pj.v, m_piston, piston_v);
            piston_v = new_vs.v2_prime;

            pj.v = new_vs.v1_prime;
            pj.t_ground = tGround(pj.x, pj.v);
            pj.t_piston = tPiston(pj.x, pj.v);
            for (t_ps, xs, vs) |*t_p, x, v| {
                t_p.* = tPiston(x, v);
            }
        }
        particles.set(j, pj);
    }

    const t1 = timer.read();
    const dt = @as(f64, @floatFromInt(t1 - t0));

    std.debug.print("\nCollisions: {}\n", .{ct});
    std.debug.print("Time taken: {d:.2}s\n", .{dt / 1000000000.0});
}
