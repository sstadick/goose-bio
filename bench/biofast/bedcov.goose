// Goose port of biofast's `bedcov_c1_cgr.c`: load one BED file into an
// implicit interval tree, then print, for every interval of a second BED
// file, how many loaded intervals overlap it and how many of its bases they
// cover.
//
// The interval tree, a port of `cgranges.c`, is this repository's `iitree`
// module, and `bed.goose` reads the files. The overlaps arrive through a
// block instead of cgranges' array of indices.
import std;
import os;

import io from "github.com/sstadick/goose-io@v0.2.2";
import bstr from "github.com/sstadick/goose-bstr@v0.3.1";
import simd_bytes from "github.com/sstadick/goose-simd@v0.2.0";
import bio from "../..";
import bio.iitree;
import bed from ".";

fn main() {
    let arguments = args();
    if arguments.len != 3 {
        print("Usage: bedcov <loaded.bed> <streamed.bed>");
        return;
    }

    let tree = iitree::IITree<i32> { .. };
    var label: i32 = 0;
    let loaded = bed::each_bed_record(arguments[1]) { contig, start, end =>
        tree.add(contig, start, end, label);
        label += 1;
    };
    if !loaded { bed::fail(arguments[1]); }
    tree.index();

    var out: u8[>..] = [];
    let streamed = bed::each_bed_record(arguments[2]) { contig, start1, end1 =>
        var count = 0;
        var coverage = 0;
        var covered_start = 0;
        var covered_end = 0;
        tree.overlap(contig, start1, end1) { interval =>
            // Overlaps arrive by start, so each one extends the covered run
            // or starts a new one.
            let start0 = max(interval.start, start1);
            let end0 = min(interval.stop, end1);
            if start0 > covered_end {
                coverage += covered_end - covered_start;
                covered_start = start0;
                covered_end = end0;
            } else if end0 > covered_end {
                covered_end = end0;
            }
            count += 1;
        };
        coverage += covered_end - covered_start;

        format(out, contig, "\t", start1, "\t", end1, "\t", count, "\t", coverage, "\n");
        if out.len >= bed::BUFFER_SIZE {
            write_stdout(out);
            out.clear();
        }
    };
    write_stdout(out);
    if !streamed { bed::fail(arguments[2]); }
}
