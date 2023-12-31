const std = @import("std");
const math = std.math;

const g = 9.81;
const m_particle = 0.001;
const m_piston = 1.0;

var piston_x: f64 = 1;
var piston_v: f64 = 0;

const gravity: bool = false;

// get time to hit ground
fn getTimeToGround(x: f64, v: f64) f64 {
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

pub fn main() !void {
    const allocator = std.heap.page_allocator;
    const args = try std.process.argsAlloc(allocator);

    var timer = try std.time.Timer.start();
    const t0 = timer.read();

    if (args.len < 3) return error.ExpectedArgument;
    const n = try std.fmt.parseInt(u32, args[1], 0);
    const max_time = try std.fmt.parseFloat(f64, args[2]);

    var xs = try allocator.alloc(f64, n);
    defer allocator.free(xs);
    var vs = try allocator.alloc(f64, n);
    defer allocator.free(vs);
    var t_grounds = try allocator.alloc(f64, n);
    defer allocator.free(t_grounds);
    var t_pistons = try allocator.alloc(f64, n);
    defer allocator.free(t_pistons);

    const vAvg: f64 = 0.5 * math.sqrt(2 * m_piston * g * piston_x / m_particle / @as(f64, @floatFromInt(n)));
    _ = vAvg;

    var prng = std.rand.DefaultPrng.init(0);
    const rand = prng.random();

    for (0..n) |i| {
        xs[i] = piston_x * rand.float(f32);
        vs[i] = 1;
        t_grounds[i] = getTimeToGround(xs[i], vs[i]);
        t_pistons[i] = getTimeToPiston(xs[i], vs[i]);
    }

    // progress indicator
    var pct_done: u8 = 0;
    var progress = std.Progress{};
    const root_node = progress.start("Simulating", 100);
    defer root_node.end();

    var t: f64 = 0;
    var ct: usize = 0; // collision count

    std.debug.print("Particles: {d}\tWorldtime: {d}s\n", .{ n, max_time });

    var buf: [100]u8 = undefined;
    const path = try std.fmt.bufPrint(&buf, "N_{d}_Time_{d}_out.txt", .{ n, max_time });

    const file = try std.fs.cwd().createFile(
        path,
        .{ .read = true },
    );
    defer file.close();
    var writer = file.writer();

    while (t < max_time) {
        ct += 1;

        // update progress indicator
        const frac_time = (t * 100) / max_time;
        const int_frac_time = @as(u8, @intFromFloat(frac_time));

        const new_pct_done = @as(u8, int_frac_time);
        if (new_pct_done > pct_done) {
            pct_done = new_pct_done;
            root_node.completeOne();
            progress.refresh();

            for (vs) |v| {
                try std.fmt.format(writer, "{}\n", .{v});
            }
        }

        // get the next interaction
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

        // step forward in time
        for (0..n) |i| {
            if (!gravity) {
                xs[i] += vs[i] * dt;
            } else {
                xs[i] += vs[i] * dt - g * dt * dt / 2;
                vs[i] -= g * dt;
            }
            t_grounds[i] -= dt;
        }
        piston_x += piston_v * dt - g * dt * dt / 2;
        piston_v -= g * dt;
        t += dt;

        if (is_ground_col) {
            vs[j] = -vs[j];
            t_grounds[j] = math.inf(f64);

            for (0..n) |i| {
                t_pistons[i] -= dt;
            }

            t_pistons[j] = getTimeToPiston(xs[j], vs[j]);
        } else {
            const new_vs = elasticCol(m_particle, vs[j], m_piston, piston_v);
            vs[j] = new_vs.v1_prime;
            piston_v = new_vs.v2_prime;
            t_grounds[j] = getTimeToGround(xs[j], vs[j]);
            for (0..n) |i| {
                t_pistons[i] = getTimeToPiston(xs[i], vs[i]);
            }
        }
    }

    const t1 = timer.read();
    const dt = @as(f64, @floatFromInt(t1 - t0));

    std.debug.print("\nCollisions: {}\n", .{ct});
    std.debug.print("Time taken: {d:.2}s\n", .{dt / 1000000000.0});
}
