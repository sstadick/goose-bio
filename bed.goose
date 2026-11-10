import std;
import dictionary;

import io from "github.com/sstadick/goose-io@v0.2.2";
import bstr from "github.com/sstadick/goose-bstr@v0.3.1";
import simd_bytes from "github.com/sstadick/goose-simd@v0.2.0";

namespace bed;

root struct ChrLookup {
    map:      dictionary<const u8[:], u16> = dictionary<const u8[:], u16> {},
    chrs:     u8[varint][>..]              = [],
    chrs_idx: u8[varint]&[>..]             = [],
}

inline fn get_id(map: ChrLookup&, chr: const u8[:]) -> u16 {
    var next_id = map.chrs.len;
    if !map.map.contains(chr) {
        map.chrs_idx.push(map.chrs.push(chr));
        map.map.insert(map.chrs_idx[map.chrs_idx.len - 1], (map.chrs.len - 1) as u16);
    }
    let idx = map.map.get(&chr);
    assert(idx);
    idx
}

struct BedRecord {
    chrom: u16,
    start: i32,
    stop:  i32,
}

fn read_bedfile<F: fn(BedRecord&, const u8[:][:])>(bedfile: u8[:], chr_map: ChrLookup&) -> bool, const u8[:] {
    let fh, opened = io::open_read(bedfile);
    guard opened else { return false, "Failed to open bedfile"; }
    defer io::close(fh);

    let reader = io::buffered_reader(fh);

    loop {
        let line, status = io::read_line(reader);
        if status == io::LineStatus.Eof { break; }
        guard status == io::LineStatus.Line else { return false, "Invalid LineStatus"; }

        var fields: const u8[:][>..] = [];
        bstr::split_str(line, "\t") { field => fields.push(field) };
        guard fields.len >= 3 else { return false, "Invalid bed record."; }

        let chr = fields[0];
        let start = parse_pos(fields[1]);
        let stop = parse_pos(fields[2]);
        let chr_idx = chr_map.get_id(chr);

        let record = BedRecord {
            chrom: chr_idx,
            start: start,
            stop:  stop,
        };
        F(&record, fields[3..]);
    }

    return true, "";
}

inline fn parse_pos(bytes: const u8[:]) -> i32 {
    var num: i32 = 0;
    for b in bytes {
        num = (num * 10) + (b - '0');
    }
    num
}
