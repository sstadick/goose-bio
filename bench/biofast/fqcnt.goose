// Goose port of biofast's `fqcnt_c1_kseq.c`: count the records, bases, and
// quality values of a FASTA or FASTQ file.
//
// The parser, a port of `kseq.h`, is this repository's `fastx` module.
import std;
import os;

import io from "github.com/sstadick/goose-io@v0.2.2";
import bstr from "github.com/sstadick/goose-bstr@v0.3.1";
import simd_bytes from "github.com/sstadick/goose-simd@v0.2.0";
import bio from "../..";
import bio.fastx;

fn main() {
    let arguments = args();
    if arguments.len < 2 {
        write_stdout("Usage: fqcnt <in.fq>");
        return;
    }
    let file, opened = io::open_read(arguments[1]);
    if !opened {
        write_stderr(str("fqcnt: cannot read ", arguments[1], "\n"));
        exit(1);
    }
    let record = fastx::Reader { stream: fastx::Stream { file: file, .. } };
    defer io::close(record.stream.file);

    var count = 0;
    var seq_length = 0;
    var qual_length = 0;
    var status = record.read();
    while status == fastx::Status.Record {
        count += 1;
        seq_length += record.seq.len;
        qual_length += record.qual.len;
        status = record.read();
    }
    print(count, "\t", seq_length, "\t", qual_length);
    if status != fastx::Status.End { write_stderr("ERROR: malformated FASTX\n"); }
}
