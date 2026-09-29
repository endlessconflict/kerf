//! kerf: CRISPR off-target search over haplotypes. A haplotype is the
//! reference with a set of variants applied; guides are scanned along it
//! (both strands) in search mode with unit edit costs, so mismatches and
//! bulges both count, and a hit needs an NGG PAM right after the protospacer.
//! The scan is bitdp's packed kernel, several guides per lane word.

const std = @import("std");
const bitdp = @import("bitdp");

/// Replace `ref_len` reference bases at 0-based `pos` by `alt` (possibly empty).
pub const Variant = struct { pos: u64, ref_len: u64, alt: []const u8 };

pub const Haplotype = struct {
    seq: []u8,
    /// Reference position of each haplotype base; inserted bases take the
    /// position of the variant that inserted them.
    map: []u32,

    pub fn deinit(h: *Haplotype, gpa: std.mem.Allocator) void {
        gpa.free(h.seq);
        gpa.free(h.map);
    }
};

/// The reference with `variants` applied. They must be sorted by position;
/// a variant that overlaps the previous one is skipped.
pub fn haplotype(gpa: std.mem.Allocator, ref: []const u8, variants: []const Variant) !Haplotype {
    var seq: std.ArrayList(u8) = .empty;
    errdefer seq.deinit(gpa);
    var map: std.ArrayList(u32) = .empty;
    errdefer map.deinit(gpa);
    var at: u64 = 0; // next unread reference base
    for (variants) |v| {
        if (v.pos < at or v.pos + v.ref_len > ref.len) continue;
        for (at..v.pos) |i| {
            try seq.append(gpa, ref[i]);
            try map.append(gpa, @intCast(i));
        }
        for (v.alt) |c| {
            try seq.append(gpa, c);
            try map.append(gpa, @intCast(v.pos));
        }
        at = v.pos + v.ref_len;
    }
    for (at..ref.len) |i| {
        try seq.append(gpa, ref[i]);
        try map.append(gpa, @intCast(i));
    }
    return .{ .seq = try seq.toOwnedSlice(gpa), .map = try map.toOwnedSlice(gpa) };
}

pub const Strand = enum { plus, minus };

/// An off-target site: guide `guide` aligns with cost `cost` to the protospacer
/// ending at `ref_pos` (the reference position of its 3' base on the strand
/// searched), followed by an NGG PAM.
pub const Site = struct { guide: u32, strand: Strand, ref_pos: u32, cost: i64 };

const Search = bitdp.Kernel(.{ .mode = .search });

fn isPam(text: []const u8, end: usize) bool {
    return end + 3 <= text.len and text[end + 1] == 'G' and text[end + 2] == 'G';
}

fn revcomp(gpa: std.mem.Allocator, s: []const u8) ![]u8 {
    const out = try gpa.alloc(u8, s.len);
    for (s, 0..) |c, i| out[s.len - 1 - i] = switch (c) {
        'A' => 'T',
        'C' => 'G',
        'G' => 'C',
        'T' => 'A',
        else => 'N',
    };
    return out;
}

/// All sites of `guides` on `hap` with cost at most `max_cost`, appended to `out`.
/// Guides are 1..62 bases over ACGT. Several hits can belong to one site
/// (neighbouring end positions); each is reported.
pub fn scan(gpa: std.mem.Allocator, hap: Haplotype, guides: []const []const u8, max_cost: i64, out: *std.ArrayList(Site)) !void {
    const rc = try revcomp(gpa, hap.seq);
    defer gpa.free(rc);
    var hits: std.ArrayList(Search.Hit) = .empty;
    defer hits.deinit(gpa);
    inline for ([_]Strand{ .plus, .minus }) |strand| {
        const text: []const u8 = if (strand == .plus) hap.seq else rc;
        var start: usize = 0;
        while (start < guides.len) {
            const pk = Search.Packed.init(guides[start..]);
            hits.clearRetainingCapacity();
            try pk.scan(gpa, text, max_cost, &hits);
            for (hits.items) |h| {
                if (!isPam(text, h.end)) continue;
                const idx = if (strand == .plus) h.end - 1 else text.len - h.end;
                try out.append(gpa, .{ .guide = @intCast(start + h.lane), .strand = strand, .ref_pos = hap.map[idx], .cost = h.cost });
            }
            start += pk.count;
        }
    }
}

// ---------------------------------------------------------------- tests

/// Brute force: DP of every guide over the text, minimum cost ending at each position.
fn brute(gpa: std.mem.Allocator, hap: Haplotype, guides: []const []const u8, max_cost: i64, out: *std.ArrayList(Site)) !void {
    const rc = try revcomp(gpa, hap.seq);
    defer gpa.free(rc);
    for ([_]Strand{ .plus, .minus }) |strand| {
        const text: []const u8 = if (strand == .plus) hap.seq else rc;
        for (guides, 0..) |g, gi| {
            const col = try gpa.alloc(i64, g.len + 1);
            defer gpa.free(col);
            for (col, 0..) |*x, i| x.* = @intCast(i);
            for (text, 1..) |c, end| {
                var diag = col[0];
                col[0] = 0;
                for (g, 1..) |gc, i| {
                    const cell = @min(diag + @intFromBool(gc != c), @min(col[i] + 1, col[i - 1] + 1));
                    diag = col[i];
                    col[i] = cell;
                }
                if (col[g.len] <= max_cost and isPam(text, end)) {
                    const idx = if (strand == .plus) end - 1 else text.len - end;
                    try out.append(gpa, .{ .guide = @intCast(gi), .strand = strand, .ref_pos = hap.map[idx], .cost = col[g.len] });
                }
            }
        }
    }
}

fn lessSite(_: void, a: Site, b: Site) bool {
    if (a.guide != b.guide) return a.guide < b.guide;
    if (a.strand != b.strand) return @intFromEnum(a.strand) < @intFromEnum(b.strand);
    if (a.ref_pos != b.ref_pos) return a.ref_pos < b.ref_pos;
    return a.cost < b.cost;
}

test "packed scan equals a brute-force DP over every haplotype" {
    const gpa = std.testing.allocator;
    var arena_state = std.heap.ArenaAllocator.init(gpa);
    defer arena_state.deinit();
    const arena = arena_state.allocator();
    var prng = std.Random.DefaultPrng.init(21);
    const rnd = prng.random();
    for (0..40) |_| {
        _ = arena_state.reset(.retain_capacity);
        const ref = try arena.alloc(u8, 400);
        for (ref) |*c| c.* = "ACGT"[rnd.int(u2)];
        // Variants: SNVs, small insertions and deletions, sorted and spaced.
        var vs: std.ArrayList(Variant) = .empty;
        var pos: u64 = 5;
        while (pos < 380) : (pos += rnd.intRangeAtMost(u64, 8, 30)) {
            const kind = rnd.uintLessThan(u8, 3);
            const alt = try arena.alloc(u8, if (kind == 1) 0 else if (kind == 2) rnd.intRangeAtMost(usize, 2, 3) else 1);
            for (alt) |*c| c.* = "ACGT"[rnd.int(u2)];
            try vs.append(arena, .{ .pos = pos, .ref_len = if (kind == 1) 2 else 1, .alt = alt });
        }
        // Guides drawn from the reference, some mutated, some random; the packed
        // kernel is exercised with more guides than fit one word.
        var guides: std.ArrayList([]const u8) = .empty;
        for (0..9) |_| {
            const m = rnd.intRangeAtMost(usize, 10, 22);
            const g = try arena.alloc(u8, m);
            const at = rnd.uintLessThan(usize, ref.len - m - 3);
            @memcpy(g, ref[at..][0..m]);
            g[rnd.uintLessThan(usize, m)] = "ACGT"[rnd.int(u2)];
            @memcpy(ref[at + m ..][0..3], "AGG"); // a PAM after the source of the guide
            try guides.append(arena, g);
        }
        for ([_][]const Variant{ &.{}, vs.items }) |set| {
            var hap = try haplotype(gpa, ref, set);
            defer hap.deinit(gpa);
            var want: std.ArrayList(Site) = .empty;
            var got: std.ArrayList(Site) = .empty;
            try brute(arena, hap, guides.items, 2, &want);
            try scan(arena, hap, guides.items, 2, &got);
            try std.testing.expect(want.items.len > 0);
            std.mem.sort(Site, want.items, {}, lessSite);
            std.mem.sort(Site, got.items, {}, lessSite);
            try std.testing.expectEqual(want.items.len, got.items.len);
            for (want.items, got.items) |a, b| try std.testing.expectEqual(a, b);
        }
    }
}

test "a variant can create an off-target that the reference lacks" {
    const gpa = std.testing.allocator;
    // Reference: guide with three mismatches and a PAM; an SNV fixes one of them.
    const guide = "ACGTACGTACGTACGTACGT";
    var ref: [80]u8 = @splat('T');
    @memcpy(ref[20..][0..20], guide);
    ref[25] = 'T'; // two mismatches against the guide, at 5 and 15
    ref[35] = 'A';
    ref[40] = 'A';
    ref[41] = 'G';
    ref[42] = 'G';
    var none = try haplotype(gpa, &ref, &.{});
    defer none.deinit(gpa);
    var out: std.ArrayList(Site) = .empty;
    defer out.deinit(gpa);
    try scan(gpa, none, &.{guide}, 1, &out);
    try std.testing.expectEqual(@as(usize, 0), out.items.len);
    // An SNV restores the guide base at 35, leaving one mismatch.
    var alt = try haplotype(gpa, &ref, &.{.{ .pos = 35, .ref_len = 1, .alt = "T" }});
    defer alt.deinit(gpa);
    try scan(gpa, alt, &.{guide}, 1, &out);
    try std.testing.expect(out.items.len > 0);
}
