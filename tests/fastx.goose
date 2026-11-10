import std;
import os;

import io from "github.com/sstadick/goose-io@v0.2.2";
import bstr from "github.com/sstadick/goose-bstr@v0.3.1";
import simd_bytes from "github.com/sstadick/goose-simd@v0.2.0";
import bio from "..";
import bio.fastx;

let PATH = "fastx-test.tmp";

// Parse text as a file. Appends `name|comment|seq|qual\n` to out for every
// record and returns the status that ended the file.
fn parse(text: const u8[:], out: u8[>..]&) -> fastx::Status {
    assert(write_file(PATH, text));
    defer delete_file(PATH);
    let file, opened = io::open_read(PATH);
    assert(opened);
    let reader = fastx::Reader { stream: fastx::Stream { file: file, .. } };
    defer io::close(reader.stream.file);

    var status = reader.read();
    while status == fastx::Status.Record {
        format(out, reader.name, "|", reader.comment, "|", reader.seq, "|", reader.qual, "\n");
        status = reader.read();
    }
    status
}

fn expect(text: const u8[:], records: const u8[:], status: fastx::Status) {
    var out: u8[>..] = [];
    let found = parse(text, out);
    assert(out == records);
    assert(found == status);
}

fn test_empty_input_has_no_records() {
    expect("", "", fastx::Status.End);
    expect("\n\n", "", fastx::Status.End);
    expect("no header here\n", "", fastx::Status.End);
}

fn test_fastq_records() {
    expect(
        "@r1 first read\nACGT\n+\nIIII\n@r2\nGG\n+r2\nHH\n",
        "r1|first read|ACGT|IIII\nr2||GG|HH\n",
        fastx::Status.End,
    );
    // A tab also ends the name, and the comment keeps its inner blanks.
    expect("@r1\tx  y\nA\n+\nI\n", "r1|x  y|A|I\n", fastx::Status.End);
    // Quality that begins with a header character is still quality.
    expect("@r1\nACGT\n+\n@>+I\n@r2\nA\n+\n>\n", "r1||ACGT|@>+I\nr2||A|>\n", fastx::Status.End);
}

fn test_fasta_and_mixed_records() {
    expect(
        ">chr1 one\nACGT\nAC\n\nGT\n>chr2\nTT\n",
        "chr1|one|ACGTACGT|\nchr2||TT|\n",
        fastx::Status.End,
    );
    expect(
        ">a\nAC\n@b\nGT\n+\nII\n>c\nT\n",
        "a||AC|\nb||GT|II\nc||T|\n",
        fastx::Status.End,
    );
    // A header without a sequence is an empty record.
    expect(">only\n", "only|||\n", fastx::Status.End);
}

fn test_wrapped_fastq_and_line_endings() {
    expect("@r1\nACGT\nACGT\n+\nIIII\nHHHH\n", "r1||ACGTACGT|IIIIHHHH\n", fastx::Status.End);
    expect("@r1 c\r\nACGT\r\n+\r\nIIII\r\n", "r1|c|ACGT|IIII\n", fastx::Status.End);
    // The last line needs no terminator.
    expect("@r1\nACGT\n+\nIIII", "r1||ACGT|IIII\n", fastx::Status.End);
    expect(">a\nACGT", "a||ACGT|\n", fastx::Status.End);
}

fn test_truncated_fastq() {
    expect("@r1\nACGT\n+\nII\n", "", fastx::Status.Truncated);
    expect("@r1\nACGT\n+\n", "", fastx::Status.Truncated);
    expect("@r1\nACGT\n+", "", fastx::Status.Truncated);
    expect("@r1\nA\n+\nI\n@r2\nACGT\n+\nIIIII\n", "r1||A|I\n", fastx::Status.Truncated);
}

// Records longer than the stream's 16 KiB buffer, and headers and line ends
// that land on every side of its boundary.
fn test_records_across_the_buffer() {
    for padding in 16370..16400 {
        var text:     u8[>..] = [];
        var expected: u8[>..] = [];
        text.append("@pad\n");
        for _ in 0..padding { text.push('A'); }
        text.append("\n+\n");
        for _ in 0..padding { text.push('I'); }
        text.append("\n@next read two\nACGT\nTG\n+\nIIII\nHH\n>tail\n");
        for line in 0..700 {
            for i in 0..60 { text.push("ACGT"[(line + i) % 4]); }
            text.push('\n');
        }

        expected.append("pad||");
        for _ in 0..padding { expected.push('A'); }
        expected.push('|');
        for _ in 0..padding { expected.push('I'); }
        expected.append("\nnext|read two|ACGTTG|IIIIHH\ntail||");
        for line in 0..700 {
            for i in 0..60 { expected.push("ACGT"[(line + i) % 4]); }
        }
        expected.append("|\n");
        expect(text, expected, fastx::Status.End);
    }
}

fn main() {
    test_empty_input_has_no_records();
    test_fastq_records();
    test_fasta_and_mixed_records();
    test_wrapped_fastq_and_line_endings();
    test_truncated_fastq();
    test_records_across_the_buffer();
    print("fastx tests passed");
}
