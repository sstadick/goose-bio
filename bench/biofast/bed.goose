// BED reading for the Goose `bedcov` program: `parse_bed` and the line loop
// of `bedcov_c1_cgr.c`, with goose-io's buffered reader standing in for
// `kstream_t`.
import std;
import os;

import io;

namespace bed;

const BUFFER_SIZE = 0x10000;

// `isspace` in the C locale; std's `is_space` leaves out '\v' and '\f'.
inline fn is_c_space(byte: u8) -> bool { byte == ' ' || (byte >= '\t' && byte <= '\r') }

// C's `atol` narrowed to 32 bits: optional blanks and sign, then digits.
fn atoi(text: const u8[:]) -> i32 {
    var digits = text;
    while digits.len > 0 && is_c_space(digits[0]) { digits = digits[1..]; }
    var negative = false;
    if digits.len > 0 && (digits[0] == '-' || digits[0] == '+') {
        negative = digits[0] == '-';
        digits = digits[1..];
    }
    var value = 0;
    for digit in digits {
        guard digit >= '0' && digit <= '9' else { break; }
        value = value * 10 + (digit - '0');
    }
    (if negative { -value } else { value }) as! i32
}

// `parse_bed`: the contig, start, and end of a line, and whether the line has
// the three tab-separated fields that make it a BED record.
fn parse_bed(line: const u8[:]) -> const u8[:], i32, i32, bool {
    var contig = line[..0];
    var start: i32 = -1;
    var end:   i32 = -1;
    var fields = 0;
    each_split(line, '\t') { field =>
        if fields == 0 {
            contig = field;
        } else if fields == 1 {
            start = atoi(field);
        } else if fields == 2 {
            end = atoi(field);
        }
        fields += 1;
    };
    return contig, start, end, fields >= 3;
}

// Call F with the contig, start, and end of every BED record of the file at
// path. The contig is only valid during the call. Returns whether the file
// could be read to its end.
fn each_bed_record<F: fn(const u8[:], i32, i32)>(path: const u8[:]) -> bool {
    let file, opened = io::open_read(path);
    guard opened else { return false; }
    let reader = io::buffered_reader(file, BUFFER_SIZE);
    defer io::close(reader);

    loop {
        let line, status = io::read_line(reader);
        guard status == io::LineStatus.Line else { return status == io::LineStatus.Eof; }
        let contig, start, end, is_record = parse_bed(line);
        if is_record { F(contig, start, end); }
    }
}

// Report that the file at path could not be read and exit.
fn fail(path: const u8[:]) {
    write_stderr(str("bedcov: cannot read ", path, "\n"));
    exit(1);
}
