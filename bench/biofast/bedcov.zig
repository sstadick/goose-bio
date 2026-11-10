//! Zig port of biofast's `bedcov_c1_cgr.c` and the parts of `cgranges.c` it
//! uses: load one BED file into an implicit interval tree, then print, for
//! every interval of a second BED file, how many loaded intervals overlap it
//! and how many of its bases they cover.
//!
//! The port follows the C line by line where that reads naturally. Lines are
//! slices of the reader's buffer, so a line must fit in `buffer_size` bytes.

const std = @import("std");
const Io = std.Io;
const Allocator = std.mem.Allocator;

const buffer_size = 0x10000;

// ---------------------------------------------------------------------------
// Radix sort
// ---------------------------------------------------------------------------

const min_size = 64;
const radix_bits = 8;
const bucket_count = 1 << radix_bits;

const Bucket = struct {
    begin: [*]Interval,
    end: [*]Interval,
};

inline fn digit(interval: Interval, shift: u6) usize {
    return @intCast((interval.x >> shift) & (bucket_count - 1));
}

fn insertionSort(begin: [*]Interval, end: [*]Interval) void {
    var current = begin + 1;
    while (@intFromPtr(current) < @intFromPtr(end)) : (current += 1) {
        if (current[0].x < (current - 1)[0].x) {
            const temporary = current[0];
            var destination = current;
            while (destination != begin and temporary.x < (destination - 1)[0].x) : (destination -= 1) {
                destination[0] = (destination - 1)[0];
            }
            destination[0] = temporary;
        }
    }
}

/// `rs_sort_cr_intv`.
fn radixSortInner(begin: [*]Interval, end: [*]Interval, shift: u6) void {
    var buckets: [bucket_count]Bucket = undefined;

    for (&buckets) |*bucket| bucket.* = .{ .begin = begin, .end = begin };
    var item = begin;
    while (item != end) : (item += 1) buckets[digit(item[0], shift)].end += 1;

    for (1..bucket_count) |bucket| {
        const previous_end = buckets[bucket - 1].end;
        buckets[bucket].end += previous_end - begin;
        buckets[bucket].begin = previous_end;
    }

    var bucket: usize = 0;
    while (bucket < bucket_count) {
        const current = &buckets[bucket];
        if (current.begin != current.end) {
            var target = digit(current.begin[0], shift);
            if (target != bucket) {
                var temporary = current.begin[0];
                while (true) {
                    const destination = &buckets[target];
                    const displaced = destination.begin[0];
                    destination.begin[0] = temporary;
                    destination.begin += 1;
                    temporary = displaced;
                    target = digit(temporary, shift);
                    if (target == bucket) break;
                }
                current.begin[0] = temporary;
            }
            current.begin += 1;
        } else {
            bucket += 1;
        }
    }

    buckets[0].begin = begin;
    for (1..bucket_count) |k| buckets[k].begin = buckets[k - 1].end;

    if (shift != 0) {
        const next_shift = shift -| radix_bits;
        for (buckets) |current| {
            const length = current.end - current.begin;
            if (length > min_size) {
                radixSortInner(current.begin, current.end, next_shift);
            } else if (length > 1) {
                insertionSort(current.begin, current.end);
            }
        }
    }
}

fn radixSort(intervals: []Interval) void {
    const begin = intervals.ptr;
    const end = begin + intervals.len;
    if (intervals.len <= min_size) {
        insertionSort(begin, end);
    } else {
        radixSortInner(begin, end, 64 - radix_bits);
    }
}

// ---------------------------------------------------------------------------
// cgranges
// ---------------------------------------------------------------------------

/// `cr_intv_t`. Before indexing, `x` is `contig << 32 | start`; afterwards it
/// is `start << 32 | end` and `y` is the largest end in the interval's subtree.
const Interval = extern struct {
    x: u64,
    y: i32,
    /// Carried, as in C, so a record is 16 bytes.
    label: i32,

    inline fn start(interval: Interval) i32 {
        return @bitCast(@as(u32, @truncate(interval.x >> 32)));
    }

    inline fn end(interval: Interval) i32 {
        return @bitCast(@as(u32, @truncate(interval.x)));
    }
};

/// `cr_ctg_t` without the name, which the map owns, and the unused length.
const Contig = struct {
    /// Height of the contig's tree.
    root_k: u6 = 0,
    /// The contig's intervals are `intervals[offset..][0..count]`.
    offset: usize = 0,
    count: usize = 0,
};

/// One node of the traversal in `overlap`: the subtree rooted at index `x`
/// with height `k`, and whether its left child has been handled.
const StackEntry = struct {
    x: usize,
    k: u6,
    left_done: bool,
};

const Cgranges = struct {
    intervals: std.ArrayList(Interval) = .empty,
    contigs: std.ArrayList(Contig) = .empty,
    contig_ids: std.StringHashMapUnmanaged(u32) = .empty,

    fn deinit(ranges: *Cgranges, gpa: Allocator) void {
        var names = ranges.contig_ids.keyIterator();
        while (names.next()) |name| gpa.free(name.*);
        ranges.contig_ids.deinit(gpa);
        ranges.contigs.deinit(gpa);
        ranges.intervals.deinit(gpa);
    }

    /// `cr_add_ctg`.
    fn addContig(ranges: *Cgranges, gpa: Allocator, name: []const u8) !u32 {
        const entry = try ranges.contig_ids.getOrPut(gpa, name);
        if (!entry.found_existing) {
            errdefer ranges.contig_ids.removeByPtr(entry.key_ptr);
            entry.key_ptr.* = try gpa.dupe(u8, name);
            entry.value_ptr.* = @intCast(ranges.contigs.items.len);
            try ranges.contigs.append(gpa, .{});
        }
        return entry.value_ptr.*;
    }

    /// `cr_add`.
    fn add(ranges: *Cgranges, gpa: Allocator, contig: []const u8, start: i32, end: i32, label: i32) !void {
        if (start > end) return;
        const id = try ranges.addContig(gpa, contig);
        try ranges.intervals.append(gpa, .{
            .x = @as(u64, id) << 32 | @as(u32, @bitCast(start)),
            .y = end,
            .label = label,
        });
    }

    /// `cr_index`: sort by contig and start, then build each contig's tree.
    fn index(ranges: *Cgranges) void {
        // cr_index_prepare
        const intervals = ranges.intervals.items;
        if (!std.sort.isSorted(Interval, intervals, {}, lessThanX)) radixSort(intervals);
        var first: usize = 0;
        for (1..intervals.len + 1) |i| {
            if (i == intervals.len or intervals[i].x >> 32 != intervals[first].x >> 32) {
                const contig = &ranges.contigs.items[intervals[first].x >> 32];
                contig.offset = first;
                contig.count = i - first;
                first = i;
            }
        }
        for (intervals) |*interval| {
            interval.x = interval.x << 32 | @as(u32, @bitCast(interval.y));
            interval.y = 0;
        }

        for (ranges.contigs.items) |*contig| {
            contig.root_k = indexContig(intervals[contig.offset..][0..contig.count]);
        }
    }

    /// `cr_overlap`: replace `hits` with the indices in `intervals`, in start
    /// order, of the intervals of `contig` that overlap `start..end`.
    fn overlap(
        ranges: *const Cgranges,
        gpa: Allocator,
        contig_name: []const u8,
        start: i32,
        end: i32,
        hits: *std.ArrayList(usize),
    ) !void {
        hits.clearRetainingCapacity();
        const contig = ranges.contigs.items[ranges.contig_ids.get(contig_name) orelse return];
        const tree = ranges.intervals.items[contig.offset..][0..contig.count];
        if (tree.len == 0) return;

        var stack: [64]StackEntry = undefined;
        var depth: usize = 1;
        stack[0] = .{ .x = (@as(usize, 1) << contig.root_k) - 1, .k = contig.root_k, .left_done = false };
        while (depth > 0) {
            depth -= 1;
            const node = stack[depth];
            if (node.k <= 3) {
                // The subtree has at most (1 << (k + 1)) - 1 intervals: scan it.
                const first = node.x >> node.k << node.k;
                const last = @min(first + (@as(usize, 1) << (node.k + 1)) - 1, tree.len);
                for (tree[first..last], first..) |interval, i| {
                    if (interval.start() >= end) break;
                    if (start < interval.end()) try hits.append(gpa, contig.offset + i);
                }
            } else if (!node.left_done) {
                const left = node.x - (@as(usize, 1) << (node.k - 1));
                stack[depth] = .{ .x = node.x, .k = node.k, .left_done = true };
                depth += 1;
                if (left >= tree.len or tree[left].y > start) {
                    stack[depth] = .{ .x = left, .k = node.k - 1, .left_done = false };
                    depth += 1;
                }
            } else if (node.x < tree.len and tree[node.x].start() < end) {
                if (start < tree[node.x].end()) try hits.append(gpa, contig.offset + node.x);
                const right = node.x + (@as(usize, 1) << (node.k - 1));
                stack[depth] = .{ .x = right, .k = node.k - 1, .left_done = false };
                depth += 1;
            }
        }
    }
};

fn lessThanX(_: void, a: Interval, b: Interval) bool {
    return a.x < b.x;
}

/// `cr_index1`: store in every `y` the largest end in that interval's subtree
/// and return the tree's height.
fn indexContig(tree: []Interval) u6 {
    if (tree.len == 0) return 0;
    var last_i: usize = 0;
    var last: i32 = 0;
    var i: usize = 0;
    while (i < tree.len) : (i += 2) {
        last_i = i;
        last = tree[i].end();
        tree[i].y = last;
    }
    var k: u6 = 1;
    while (@as(usize, 1) << k <= tree.len) : (k += 1) {
        const x = @as(usize, 1) << (k - 1);
        i = (x << 1) - 1;
        while (i < tree.len) : (i += x << 2) {
            const left = tree[i - x].y;
            const right = if (i + x < tree.len) tree[i + x].y else last;
            tree[i].y = @max(tree[i].end(), left, right);
        }
        last_i = if (last_i >> k & 1 != 0) last_i - x else last_i + x;
        if (last_i < tree.len and tree[last_i].y > last) last = tree[last_i].y;
    }
    return k - 1;
}

// ---------------------------------------------------------------------------
// bedcov
// ---------------------------------------------------------------------------

/// C's `atol` narrowed to 32 bits: optional blanks and sign, then digits.
fn atoi(text: []const u8) i32 {
    var i: usize = 0;
    while (i < text.len and std.ascii.isWhitespace(text[i])) i += 1;
    var negative = false;
    if (i < text.len and (text[i] == '+' or text[i] == '-')) {
        negative = text[i] == '-';
        i += 1;
    }
    var value: i64 = 0;
    while (i < text.len and std.ascii.isDigit(text[i])) : (i += 1) {
        value = value *% 10 +% (text[i] - '0');
    }
    return @truncate(if (negative) -%value else value);
}

const Bed = struct {
    contig: []const u8,
    start: i32,
    end: i32,
};

/// `parse_bed`: the contig, start, and end of a line with at least three
/// tab-separated fields.
fn parseBed(line: []const u8) ?Bed {
    var fields = std.mem.splitScalar(u8, line, '\t');
    const contig = fields.next() orelse return null;
    const start = atoi(fields.next() orelse return null);
    const end = atoi(fields.next() orelse return null);
    return .{ .contig = contig, .start = start, .end = end };
}

/// The next line without its terminator, or null at the end of the input.
fn takeLine(reader: *Io.Reader) !?[]const u8 {
    const line = try reader.takeDelimiter('\n') orelse return null;
    if (line.len > 1 and line[line.len - 1] == '\r') return line[0 .. line.len - 1];
    return line;
}

fn readBed(io: Io, gpa: Allocator, path: []const u8, buffer: []u8) !Cgranges {
    const file = try Io.Dir.cwd().openFile(io, path, .{});
    defer file.close(io);
    var file_reader = file.reader(io, buffer);
    const reader = &file_reader.interface;

    var ranges: Cgranges = .{};
    errdefer ranges.deinit(gpa);
    var label: i32 = 0;
    while (try takeLine(reader)) |line| {
        const bed = parseBed(line) orelse continue;
        try ranges.add(gpa, bed.contig, bed.start, bed.end, label);
        label += 1;
    }
    return ranges;
}

pub fn main(init: std.process.Init) !void {
    const io = init.io;
    const gpa = init.gpa;
    var stdout_buffer: [buffer_size]u8 = undefined;
    var stdout_writer: Io.File.Writer = .init(.stdout(), io, &stdout_buffer);
    const out = &stdout_writer.interface;

    const arguments = try init.minimal.args.toSlice(init.arena.allocator());
    if (arguments.len != 3) {
        try out.writeAll("Usage: bedcov <loaded.bed> <streamed.bed>\n");
        return out.flush();
    }

    var read_buffer: [buffer_size]u8 = undefined;
    var ranges = try readBed(io, gpa, arguments[1], &read_buffer);
    defer ranges.deinit(gpa);
    ranges.index();

    const file = try Io.Dir.cwd().openFile(io, arguments[2], .{});
    defer file.close(io);
    var file_reader = file.reader(io, &read_buffer);
    const reader = &file_reader.interface;

    var hits: std.ArrayList(usize) = .empty;
    defer hits.deinit(gpa);
    while (try takeLine(reader)) |line| {
        const bed = parseBed(line) orelse continue;
        try ranges.overlap(gpa, bed.contig, bed.start, bed.end, &hits);
        var coverage: i64 = 0;
        var covered_start: i64 = 0;
        var covered_end: i64 = 0;
        for (hits.items) |hit| {
            const interval = ranges.intervals.items[hit];
            const start0: i64 = @max(interval.start(), bed.start);
            const end0: i64 = @min(interval.end(), bed.end);
            if (start0 > covered_end) {
                coverage += covered_end - covered_start;
                covered_start = start0;
                covered_end = end0;
            } else {
                covered_end = @max(covered_end, end0);
            }
        }
        coverage += covered_end - covered_start;
        try out.print("{s}\t{d}\t{d}\t{d}\t{d}\n", .{ bed.contig, bed.start, bed.end, hits.items.len, coverage });
    }
    try out.flush();
}
