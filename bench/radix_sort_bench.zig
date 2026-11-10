const std = @import("std");
const Io = std.Io;

const min_size = 64;
const radix_bits = 8;
const bucket_count = 1 << radix_bits;

const Record = extern struct {
    key: u64,
    payload: u64,
};

const Bucket = struct {
    begin: [*]Record,
    end: [*]Record,
};

extern fn c_setup(len: i64) void;
extern fn c_reset() void;
extern fn c_radix_sort() void;
extern fn c_validate() bool;
extern fn c_checksum() u64;
extern fn c_teardown() void;

inline fn digit(record: Record, shift: u6) usize {
    return @intCast((record.key >> shift) & (bucket_count - 1));
}

fn insertionSort(begin: [*]Record, end: [*]Record) void {
    var current = begin + 1;
    while (@intFromPtr(current) < @intFromPtr(end)) : (current += 1) {
        if (current[0].key < (current - 1)[0].key) {
            const temporary = current[0];
            var destination = current;
            while (destination != begin and temporary.key < (destination - 1)[0].key) : (destination -= 1) {
                destination[0] = (destination - 1)[0];
            }
            destination[0] = temporary;
        }
    }
}

// A direct port of rs_sort_record in radix_sort_reference.c: every recursive
// call keeps its own 4 KiB table of begin/end record pointers on the stack.
fn radixSortInner(begin: [*]Record, end: [*]Record, shift: u6) void {
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

fn radixSort(records: []Record) void {
    const begin = records.ptr;
    const end = begin + records.len;
    if (records.len <= min_size) {
        insertionSort(begin, end);
    } else {
        radixSortInner(begin, end, 64 - radix_bits);
    }
}

fn nextRandom(state: *u64) u64 {
    var x = state.*;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    state.* = x;
    return x;
}

fn zigSample(io: Io, work: []Record) i96 {
    const start = Io.Clock.awake.now(io);
    radixSort(work);
    return start.durationTo(Io.Clock.awake.now(io)).toNanoseconds();
}

fn cSample(io: Io) i96 {
    const start = Io.Clock.awake.now(io);
    c_radix_sort();
    return start.durationTo(Io.Clock.awake.now(io)).toNanoseconds();
}

fn checksum(work: []const Record) u64 {
    return work[0].key ^ work[work.len / 2].payload ^ work[work.len - 1].key;
}

fn runCase(io: Io, gpa: std.mem.Allocator, out: *Io.Writer, n: usize, repetitions: usize) !void {
    const payload_mask = 0x9E3779B97F4A7C15;
    const source = try gpa.alloc(Record, n);
    defer gpa.free(source);
    const work = try gpa.alloc(Record, n);
    defer gpa.free(work);
    var state: u64 = 0x4D595DF4D0F33173;
    for (source) |*record| {
        const key = nextRandom(&state);
        record.* = .{ .key = key, .payload = key ^ payload_mask };
    }

    // Warm and validate both implementations before collecting samples.
    @memcpy(work, source);
    radixSort(work);
    for (work, 0..) |record, i| {
        if (record.payload != record.key ^ payload_mask) return error.PayloadDetached;
        if (i != 0 and work[i - 1].key > record.key) return error.NotSorted;
    }
    c_setup(@intCast(n));
    c_reset();
    c_radix_sort();
    if (!c_validate()) return error.ReferenceNotSorted;

    var zig_checksum: u64 = 0;
    var reference_checksum: u64 = 0;
    for (0..repetitions) |repetition| {
        // Alternate order to balance short-term frequency and temperature changes.
        for (0..2) |turn| {
            if ((repetition + turn) % 2 == 0) {
                @memcpy(work, source);
                const elapsed = zigSample(io, work);
                zig_checksum +%= checksum(work);
                try out.print("zig,{d},{d}\n", .{ n, elapsed });
            } else {
                c_reset();
                const elapsed = cSample(io);
                reference_checksum +%= c_checksum();
                try out.print("c,{d},{d}\n", .{ n, elapsed });
            }
        }
    }
    if (zig_checksum != reference_checksum) return error.ChecksumMismatch;
    c_teardown();
}

pub fn main(init: std.process.Init) !void {
    var stdout_buffer: [1 << 16]u8 = undefined;
    var stdout_writer: Io.File.Writer = .init(.stdout(), init.io, &stdout_buffer);
    const out = &stdout_writer.interface;

    try runCase(init.io, init.gpa, out, 1024, 1000);
    try runCase(init.io, init.gpa, out, 16384, 200);
    try runCase(init.io, init.gpa, out, 262144, 20);
    try runCase(init.io, init.gpa, out, 1048576, 6);
    try out.flush();
}
