import std;
import os;
import bio from "..";
import bio.bed;
import bio.lapper;

import opts from "github.com/sstadick/goose-opts@v0.1.0";
import io from "github.com/sstadick/goose-io@v0.2.2";
import bstr from "github.com/sstadick/goose-bstr@v0.3.1";
import simd_bytes from "github.com/sstadick/goose-simd@v0.2.0";

struct Options {
    side_a: opts::Field<const u8[:]> = opts::Field<const u8[:]> {
        short: "a", required: true, help: "Input file",
    },
    side_b: opts::Field<const u8[:]> = opts::Field<const u8[:]> {
        short: "b", required: true, help: "Input file",
    },
    datastructure: opts::Field<DataStructure> = opts::Field<DataStructure> {
        short: "d", required: true, help: "The datastructure to use",
        value: DataStructure.Lapper,
    },
}

enum DataStructure {
    Lapper,
    IITree,
}

// Basically bedcov, read a into mem, stream b and ouput bedcov-like oputput
// Follow https://github.com/lh3/biofast/tree/master/bedcov
fn main() {
    let argv = args();
    let options, help_requested, error = opts::parse<Options>(argv[1..]);
    if error.len > 0 {
        write_stderr(str(error, "\n"));
        exit(2);
    }
    if help_requested {
        write_stdout(opts::help<Options>("demo"));
        return;
    }

    match options.datastructure.value {
        Lapper => lapper(options),
        IITree => print("IITree not implemented yet"),
    }
}

fn lapper(options: Options&) {
    let chr_map = bed::ChrLookup {};

    // TODO: enforce that data lives in a global data array?
    let intervals: lapper::Interval<i32, i32>[>..] = [];
    bed::read_bedfile(options.side_a.value, chr_map) { record, rest =>
        intervals.push(lapper::Interval {
            record.start, record.stop, 0 as i32,
        });
    };
    var index: i32[>..] = [];
    let lapper = lapper::new_lapper(intervals, index);

    bed::read_bedfile(options.side_b.value, chr_map) { record, rest =>
        // TODO: don't really need to do the mapping here
        var count = 0;
        lapper::find(lapper, record.start, record.stop) { iv =>
            count += 1;
        };
        print(record.chrom, "\t", record.start, "\t", record.stop, "\t", count);
    };
}
