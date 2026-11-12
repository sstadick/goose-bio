import std;
import os;
import io from "github.com/sstadick/goose-io@c3554678191cd3e1ff04d2c5207ddd526f3535fe";
import bstr from "github.com/sstadick/goose-bstr@v0.3.1";
import simd_bytes from "github.com/sstadick/goose-simd@v0.2.0";
import bio from "..";
import bio.faidx;

let PATH = "faidx-test.tmp";
let INDEX = "faidx-index.tmp";

fn test_queries() {
    let reader = faidx::Reader { .. };
    defer reader.close();
    var out: u8[>..] = [];
    assert(reader.open("tests/example.fa") == faidx::Status.Ok);
    assert(reader.references.len == 2 && reader.reference_name(1) == "chr2");
    let expected = "ACGTacgtnRYkm";
    // Every short interval, including empty, complete, line ends and final base.
    for begin in 0..expected.len + 1 {
        for end in begin..expected.len + 1 {
            assert(reader.fetch("chr1", begin, end, out) == faidx::Status.Ok);
            assert(out == expected[begin..end]);
        }
    }
    assert(reader.fetch(1, 1, 5, out) == faidx::Status.Ok && out == "TGCA");
    assert(reader.fetch(0, 0, 2, out) == faidx::Status.Ok && out == "AC");
    var saved: u8[>..] = [];
    assert(reader.fetch(0, 8, 13, saved) == faidx::Status.Ok);
    for range in [[-1, 1], [2, 1], [0, 14], [14, 14]] {
        assert(reader.fetch(0, range[0], range[1], out) == faidx::Status.OutOfRange);
        assert(out.len == 0);
    }
    assert(reader.fetch("nope", 0, 1, out) == faidx::Status.UnknownReference);
    assert(reader.fetch(-1, 0, 1, out) == faidx::Status.UnknownReference);
    assert(reader.fetch(0, 0, 1, out) == faidx::Status.Ok && out == "A");
    assert(reader.close());
    assert(saved == "nRYkm");
    assert(reader.reference_id("chr1") == -1);
    assert(reader.fetch(0, 0, 1, out) == faidx::Status.Closed && out.len == 0);
    assert(reader.close());
}

fn test_wrapping() {
    defer delete_file(PATH);
    defer delete_file(INDEX);
    // Different line lengths and terminators, including a newline split across
    // the read cache and a base run longer than that cache. No final newline.
    for width in [1, 7, 60, 65535, 65536, 65537] {
        for newline in ["\n", "\r\n", "\r"] {
            var sequence: u8[>..] = [];
            var text: u8[>..] = [];
            text.append(">wrapped description\n");
            let offset = text.len;
            for i in 0..131143 {
                let base = "ACgTNrykmBDHV"[(i * 7 + i / 19) % 13];
                sequence.push(base);
                if i > 0 && i % width == 0 { text.append(newline); }
                text.push(base);
            }
            assert(write_file(PATH, text));
            var index: u8[>..] = [];
            format(index, "wrapped\t", sequence.len, "\t", offset, "\t", width, "\t", width + newline.len, "\r\n");
            assert(write_file(INDEX, index));
            let reader = faidx::Reader { .. };
            defer reader.close();
            assert(reader.open(PATH, INDEX) == faidx::Status.Ok);
            var out: u8[>..] = [];
            assert(reader.fetch(0, 0, sequence.len, out) == faidx::Status.Ok);
            assert(out == sequence);
            // Alternating backward and forward queries must not reuse stale bytes.
            for begin in [131130, 0, 65534, 3, 65535, 65536, 65537, 17] {
                assert(reader.fetch(0, begin, begin + 13, out) == faidx::Status.Ok);
                assert(out == sequence[begin..begin + 13]);
            }
        }
    }
}

fn test_index_validation() {
    defer delete_file(PATH);
    defer delete_file(INDEX);
    assert(write_file(PATH, ">a\nAC\nGT\n"));
    let reader = faidx::Reader { .. };
    defer reader.close();
    var out: u8[>..] = [];
    for index in ["", "\n", "a\t4\t3\t2", "a\t4\t3\t2\t3\t9",
                  "a\t4\t3\t2\t3\t", "\t4\t3\t2\t3", "a b\t4\t3\t2\t3",
                  "a\t-4\t3\t2\t3", "a\t4x\t3\t2\t3", "a\t\t3\t2\t3",
                  "a\t4\t3\t0\t1", "a\t4\t3\t2\t1", "a\t4\t3\t2\t5",
                  "a\t4\t3\t5\t6", "a\t4\t3\t2\t2", "a\t0\t3\t1\t1",
                  "a\t4\t9223372036854775808\t2\t3", "a\t4\t9223372036854775805\t2\t3",
                  "a\t9223372036854775807\t3\t1\t2", "a\t4\t3\t2\t3\na\t4\t3\t2\t3\n"] {
        let candidate = faidx::Reader { .. };
        defer candidate.close();
        assert(candidate.open("tests/example.fa") == faidx::Status.Ok);
        assert(write_file(INDEX, index));
        assert(candidate.open(PATH, INDEX) == faidx::Status.InvalidIndex);
        assert(candidate.references.len == 0 && candidate.reference_id("chr1") == -1);
        var empty: u8[>..] = [];
        assert(candidate.fetch(0, 0, 1, empty) == faidx::Status.Closed);
    }
    // Unterminated FAI line and single-line sequence with no terminator.
    assert(write_file(PATH, ">a\nAcgT"));
    assert(write_file(INDEX, "a\t4\t3\t4\t4"));
    assert(reader.open(PATH, INDEX) == faidx::Status.Ok);
    assert(reader.fetch(0, 0, 4, out) == faidx::Status.Ok && out == "AcgT");
    reader.close();
    assert(write_file(INDEX, "a\t0\t3\t0\t0\n"));
    assert(reader.open(PATH, INDEX) == faidx::Status.Ok);
    assert(reader.fetch(0, 0, 0, out) == faidx::Status.Ok && out.len == 0);
    assert(reader.fetch(0, 0, 1, out) == faidx::Status.OutOfRange);
    assert(reader.open("tests/example.fa", "no-such-index") == faidx::Status.IoError);
    assert(reader.open("no-such-fasta") == faidx::Status.IoError);
}

fn test_bad_data_and_reopen() {
    defer delete_file(PATH);
    defer delete_file(INDEX);
    let reader = faidx::Reader { .. };
    defer reader.close();
    var out: u8[>..] = [];
    assert(write_file(INDEX, "a\t4\t3\t2\t3\n"));
    assert(write_file(PATH, ">a\nAC!GT\n"));
    assert(reader.open(PATH, INDEX) == faidx::Status.Ok);
    assert(reader.fetch(0, 0, 4, out) == faidx::Status.InvalidData && out.len == 0);
    assert(reader.fetch(0, 0, 2, out) == faidx::Status.Ok && out == "AC");
    reader.close();
    assert(write_file(PATH, ">a\nAC\nG"));
    assert(reader.open(PATH, INDEX) == faidx::Status.Ok);
    assert(reader.fetch(0, 0, 4, out) == faidx::Status.Truncated && out.len == 0);
    reader.close();
    assert(write_file(PATH, "\x1f\x8brest"));
    assert(reader.open(PATH, INDEX) == faidx::Status.UnsupportedFormat);
    assert(write_file(PATH, "@fastq\n"));
    assert(reader.open(PATH, INDEX) == faidx::Status.InvalidData);
    assert(reader.open("tests/example.fa") == faidx::Status.Ok);
    assert(reader.fetch(1, 0, 5, out) == faidx::Status.Ok && out == "TTGCA");
}

fn test_many_references() {
    defer delete_file(PATH);
    defer delete_file(INDEX);
    var fasta: u8[>..] = [];
    var index: u8[>..] = [];
    for i in 0..1000 {
        format(fasta, ">reference_", i, " description\n");
        format(index, "reference_", i, "\t1\t", fasta.len, "\t1\t2\n");
        fasta.push("ACGT"[i % 4]);
        fasta.push('\n');
    }
    assert(write_file(PATH, fasta));
    assert(write_file(INDEX, index));
    let reader = faidx::Reader { .. };
    defer reader.close();
    assert(reader.open(PATH, INDEX) == faidx::Status.Ok);
    var name: u8[>..] = [];
    var out: u8[>..] = [];
    for j in 0..1000 {
        let i = 999 - j;
        name.clear();
        format(name, "reference_", i);
        assert(reader.reference_id(name) == i);
        assert(reader.reference_name(i) == name);
        assert(reader.fetch(name, 0, 1, out) == faidx::Status.Ok);
        assert(out.len == 1 && out[0] == "ACGT"[i % 4]);
    }
}

fn main() {
    test_queries();
    test_wrapping();
    test_index_validation();
    test_bad_data_and_reopen();
    test_many_references();
    // The Python runner supplies a sparse fixture beyond 4 GiB. This tests
    // the native seek offset without allocating or reading several GiB.
    let args = args();
    if args.len == 2 {
        let reader = faidx::Reader { .. };
        defer reader.close();
        assert(reader.open(args[1]) == faidx::Status.Ok);
        assert(reader.references[0].offset > 4294967296);
        var out: u8[>..] = [];
        assert(reader.fetch(0, 1, 8, out) == faidx::Status.Ok && out == "cGTnRYk");
        assert(reader.fetch(0, 0, 1, out) == faidx::Status.Ok && out == "A");
    }
    print("faidx tests PASSED");
}
