//! Off-targets of random guides on the reference and on a haplotype, and the
//! sites the haplotype gains or loses.
//!
//! usage: kerf-demo REF.fa HAP.vcf GUIDES K SEED
//!
//! HAP.vcf lists the variants of one haplotype (CHROM POS ID REF ALT, sorted,
//! biallelic). Guides are random 20-mers of the reference followed by NGG.
//! K is the largest edit cost (mismatches and bulges count 1 each).

const std = @import("std");
const kerf = @import("kerf");
const nucleo = @import("nucleo");

fn now(io: std.Io) i96 {
    return std.Io.Timestamp.now(io, .awake).nanoseconds;
}

const Cluster = struct { guide: u32, strand: kerf.Strand, pos: u32, cost: i64 };

/// One entry per site: hits of the same guide and strand within 3 bases merge.
fn cluster(gpa: std.mem.Allocator, sites: []kerf.Site) ![]Cluster {
    std.mem.sort(kerf.Site, sites, {}, struct {
        fn f(_: void, a: kerf.Site, b: kerf.Site) bool {
            if (a.guide != b.guide) return a.guide < b.guide;
            if (a.strand != b.strand) return @intFromEnum(a.strand) < @intFromEnum(b.strand);
            return a.ref_pos < b.ref_pos;
        }
    }.f);
    var out: std.ArrayList(Cluster) = .empty;
    for (sites) |s| {
        if (out.items.len > 0) {
            const l = &out.items[out.items.len - 1];
            if (l.guide == s.guide and l.strand == s.strand and s.ref_pos -| l.pos <= 3) {
                l.pos = s.ref_pos;
                l.cost = @min(l.cost, s.cost);
                continue;
            }
        }
        try out.append(gpa, .{ .guide = s.guide, .strand = s.strand, .pos = s.ref_pos, .cost = s.cost });
    }
    return out.toOwnedSlice(gpa);
}

fn near(list: []const Cluster, c: Cluster) bool {
    for (list) |o| {
        if (o.guide == c.guide and o.strand == c.strand and (if (o.pos > c.pos) o.pos - c.pos else c.pos - o.pos) <= 8) return true;
    }
    return false;
}

pub fn main(init: std.process.Init) !void {
    const io = init.io;
    const gpa = init.arena.allocator();
    const args = try init.minimal.args.toSlice(gpa);
    if (args.len < 6) return error.Usage;
    const fa = try std.Io.Dir.cwd().readFileAlloc(io, args[1], gpa, .unlimited);
    var fp = nucleo.fastx.FastaParser{ .buf = fa };
    const rec = (try fp.next()) orelse return error.EmptyFasta;
    const ref = try gpa.dupe(u8, rec.seq);
    for (ref) |*c| c.* = std.ascii.toUpper(c.*);
    const contig = std.mem.sliceTo(rec.name, ' ');

    const vt = try std.Io.Dir.cwd().readFileAlloc(io, args[2], gpa, .unlimited);
    var vp = nucleo.vcf.Parser{ .buf = vt };
    var vars: std.ArrayList(kerf.Variant) = .empty;
    while (try vp.next()) |r| {
        if (!std.mem.eql(u8, r.chrom, contig) or std.mem.indexOfScalar(u8, r.alt, ',') != null) continue;
        var k: usize = 0;
        while (k < r.ref.len and k < r.alt.len and r.ref[k] == r.alt[k]) k += 1;
        if (k == r.ref.len and k == r.alt.len) continue;
        const alt = try gpa.dupe(u8, r.alt[k..]);
        if (std.mem.indexOfNone(u8, alt, "ACGT") != null) continue;
        try vars.append(gpa, .{ .pos = r.pos - 1 + k, .ref_len = r.ref.len - k, .alt = alt });
    }
    const nguides = try std.fmt.parseInt(usize, args[3], 10);
    const k = try std.fmt.parseInt(i64, args[4], 10);
    var prng = std.Random.DefaultPrng.init(try std.fmt.parseInt(u64, args[5], 10));
    const rnd = prng.random();

    const guides = try gpa.alloc([]const u8, nguides);
    for (guides) |*g| while (true) {
        const at = rnd.uintLessThan(usize, ref.len - 30);
        const s = ref[at .. at + 23];
        if (std.mem.indexOfNone(u8, s, "ACGT") == null and s[21] == 'G' and s[22] == 'G') {
            g.* = s[0..20];
            break;
        }
    };

    const ref_hap = try kerf.haplotype(gpa, ref, &.{});
    const alt_hap = try kerf.haplotype(gpa, ref, vars.items);
    std.debug.print("{s}: {d} bp, {d} variants applied, {d} guides, K={d}\n", .{ contig, ref.len, vars.items.len, nguides, k });

    var ref_sites: std.ArrayList(kerf.Site) = .empty;
    var alt_sites: std.ArrayList(kerf.Site) = .empty;
    const t0 = now(io);
    try kerf.scan(gpa, ref_hap, guides, k, &ref_sites);
    const t1 = now(io);
    try kerf.scan(gpa, alt_hap, guides, k, &alt_sites);
    const t2 = now(io);
    const rs = try cluster(gpa, ref_sites.items);
    const as = try cluster(gpa, alt_sites.items);
    var gained: usize = 0;
    var lost: usize = 0;
    var gained_lt: usize = 0;
    for (as) |c| if (!near(rs, c)) {
        gained += 1;
    };
    for (rs) |c| if (!near(as, c)) {
        lost += 1;
    };
    for (as) |c| if (!near(rs, c) and c.cost < k) {
        gained_lt += 1;
    };
    const cells = @as(f64, @floatFromInt(ref.len)) * 2 * 20 * @as(f64, @floatFromInt(nguides));
    std.debug.print("reference: {d} sites in {d:.2}s ({d:.1} GCUPS)\n", .{ rs.len, @as(f64, @floatFromInt(t1 - t0)) * 1e-9, cells / (@as(f64, @floatFromInt(t1 - t0)) * 1e-9) * 1e-9 });
    std.debug.print("haplotype: {d} sites in {d:.2}s\n", .{ as.len, @as(f64, @floatFromInt(t2 - t1)) * 1e-9 });
    std.debug.print("gained by the haplotype: {d} ({d} below the cost limit), lost: {d}\n", .{ gained, gained_lt, lost });
}
