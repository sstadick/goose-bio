//! Zig port of biofast's `fqcnt_c1_kseq.c` and the FASTA/FASTQ parser of
//! `kseq.h`: count the records, bases, and quality values of a file.
//!
//! `std.Io.Reader` plays the part of `kstream_t`; the parser follows
//! `kseq_read`, including multi-line records.

const std = @import("std");
const Io = std.Io;
const Allocator = std.mem.Allocator;

const buffer_size = 16384;

/// `kseq_t`: the most recently read record and the stream it came from.
const Kseq = struct {
    name: std.ArrayList(u8) = .empty,
    comment: std.ArrayList(u8) = .empty,
    seq: std.ArrayList(u8) = .empty,
    qual: std.ArrayList(u8) = .empty,
    /// The `>` or `@` of the next record, when reading a sequence consumed it.
    last_char: ?u8 = null,
    reader: *Io.Reader,
    gpa: Allocator,

    const Error = Allocator.Error || error{ ReadFailed, TruncatedQuality };

    fn deinit(kseq: *Kseq) void {
        kseq.name.deinit(kseq.gpa);
        kseq.comment.deinit(kseq.gpa);
        kseq.seq.deinit(kseq.gpa);
        kseq.qual.deinit(kseq.gpa);
    }

    /// `ks_getc`: the next byte, or null at the end of the input.
    fn getc(kseq: *Kseq) error{ReadFailed}!?u8 {
        return kseq.reader.takeByte() catch |err| switch (err) {
            error.EndOfStream => null,
            error.ReadFailed => |e| e,
        };
    }

    /// The bytes read ahead, after reading more if there are none. Empty at
    /// the end of the input.
    fn buffered(kseq: *Kseq) error{ReadFailed}![]const u8 {
        while (kseq.reader.bufferedLen() == 0) {
            kseq.reader.fillMore() catch |err| switch (err) {
                error.EndOfStream => break,
                error.ReadFailed => |e| return e,
            };
        }
        return kseq.reader.buffered();
    }

    /// `ks_getuntil2` with `KS_SEP_SPACE`: append the bytes before the next
    /// whitespace to `text` and return that whitespace, or null at the end of
    /// the input.
    fn readWord(kseq: *Kseq, text: *std.ArrayList(u8)) Error!?u8 {
        while (true) {
            const bytes = try kseq.buffered();
            if (bytes.len == 0) return null;
            const i = std.mem.findAny(u8, bytes, &std.ascii.whitespace) orelse bytes.len;
            try text.appendSlice(kseq.gpa, bytes[0..i]);
            if (i < bytes.len) {
                const delimiter = bytes[i];
                kseq.reader.toss(i + 1);
                return delimiter;
            }
            kseq.reader.toss(i);
        }
    }

    /// `ks_getuntil2` with `KS_SEP_LINE`: append the rest of the line to
    /// `text`, without its terminator. Returns whether the input had anything
    /// left to read.
    fn readLine(kseq: *Kseq, text: *std.ArrayList(u8)) Error!bool {
        var got_any = false;
        while (true) {
            const bytes = try kseq.buffered();
            if (bytes.len == 0) break;
            got_any = true;
            const i = std.mem.findScalar(u8, bytes, '\n') orelse bytes.len;
            try text.appendSlice(kseq.gpa, bytes[0..i]);
            if (i < bytes.len) {
                kseq.reader.toss(i + 1);
                break;
            }
            kseq.reader.toss(i);
        }
        if (got_any and text.items.len > 1 and text.items[text.items.len - 1] == '\r') text.items.len -= 1;
        return got_any;
    }

    /// `kseq_read`: read the next record into `kseq`. Returns false at the
    /// end of the input and `error.TruncatedQuality` for a FASTQ record whose
    /// quality string is missing or not as long as its sequence.
    fn read(kseq: *Kseq) Error!bool {
        if (kseq.last_char == null) {
            // Jump to the next header line.
            while (true) {
                const byte = try kseq.getc() orelse return false;
                if (byte == '>' or byte == '@') {
                    kseq.last_char = byte;
                    break;
                }
            }
        }
        kseq.name.clearRetainingCapacity();
        kseq.comment.clearRetainingCapacity();
        kseq.seq.clearRetainingCapacity();
        kseq.qual.clearRetainingCapacity();

        const delimiter = try kseq.readWord(&kseq.name);
        if (delimiter == null and kseq.name.items.len == 0) return false;
        if (delimiter != '\n') _ = try kseq.readLine(&kseq.comment);

        const next = while (true) {
            const byte = try kseq.getc() orelse break null;
            switch (byte) {
                '>', '+', '@' => break byte,
                '\n' => {}, // Skip empty lines.
                else => {
                    try kseq.seq.append(kseq.gpa, byte);
                    _ = try kseq.readLine(&kseq.seq);
                },
            }
        };
        if (next == '>' or next == '@') kseq.last_char = next; // The next record's header has begun.
        if (next != '+') return true; // FASTA

        // Skip the rest of the '+' line.
        while (true) {
            const byte = try kseq.getc() orelse return error.TruncatedQuality;
            if (byte == '\n') break;
        }
        while (try kseq.readLine(&kseq.qual) and kseq.qual.items.len < kseq.seq.items.len) {}
        kseq.last_char = null; // The next header line is still ahead.
        if (kseq.seq.items.len != kseq.qual.items.len) return error.TruncatedQuality;
        return true;
    }
};

pub fn main(init: std.process.Init) !void {
    const io = init.io;
    var stdout_buffer: [256]u8 = undefined;
    var stdout_writer: Io.File.Writer = .init(.stdout(), io, &stdout_buffer);
    const out = &stdout_writer.interface;

    const arguments = try init.minimal.args.toSlice(init.arena.allocator());
    if (arguments.len != 2) {
        try out.writeAll("Usage: fqcnt <in.fq>");
        return out.flush();
    }

    const file = try Io.Dir.cwd().openFile(io, arguments[1], .{});
    defer file.close(io);
    var read_buffer: [buffer_size]u8 = undefined;
    var file_reader = file.reader(io, &read_buffer);
    var record: Kseq = .{ .reader = &file_reader.interface, .gpa = init.gpa };
    defer record.deinit();

    var count: u64 = 0;
    var seq_length: usize = 0;
    var qual_length: usize = 0;
    const malformed = while (true) {
        const more = record.read() catch |err| switch (err) {
            error.TruncatedQuality => break true,
            else => |e| return e,
        };
        if (!more) break false;
        count += 1;
        seq_length += record.seq.items.len;
        qual_length += record.qual.items.len;
    };
    try out.print("{d}\t{d}\t{d}\n", .{ count, seq_length, qual_length });
    try out.flush();
    if (malformed) std.debug.print("ERROR: malformated FASTX\n", .{});
}
