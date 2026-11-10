// FASTA and FASTQ reader.
//
// A `Reader` parses a file one record at a time into buffers it reuses:
// `read` fills name, comment, seq, and qual and says what it found. Records
// may wrap their sequence and quality over several lines, and FASTA and FASTQ
// records may be mixed in one file.
//
// The implementation is a port of the parser in Heng Li's kseq.h:
// https://github.com/lh3/biofast/blob/4177ef282f4e1addc5b20beb7ee2f6435806f844/lib/kseq.h
import std;

import io from "github.com/sstadick/goose-io@v0.2.2";
import bstr from "github.com/sstadick/goose-bstr@v0.3.1";

namespace fastx;

const BUFFER_SIZE = 16384;

// A file read through a buffer whose unread bytes are `buffer[begin..end]`;
// kseq.h's `kstream_t`. Construct one with `Stream { file: file, .. }`. It
// does not close the file.
//
// ```goose
// let file, opened = io::open_read("tests/example.fq");
// assert(opened);
// let reader = fastx::Reader { stream: fastx::Stream { file: file, .. } };
// defer io::close(reader.stream.file);
// assert(reader.read() == fastx::Status.Record);
// assert(!reader.stream.failed);
// ```
struct Stream {
    file:   io::FileHandle,
    begin:  i64  = 0,
    end:    i64  = 0,
    is_eof: bool = false,
    failed: bool = false,
    buffer: u8[BUFFER_SIZE],
}

// Make sure the buffer has unread bytes. Returns false once the input has
// ended or failed.
inline fn fill(stream: Stream&) -> bool {
    stream.begin < stream.end || refill(stream)
}

fn refill(stream: Stream&) -> bool {
    guard !stream.is_eof else { return false; }
    let count = io::read(stream.file, stream.buffer);
    if count <= 0 {
        stream.is_eof = true;
        stream.failed = count < 0;
        return false;
    }
    stream.begin = 0;
    stream.end = count;
    true
}

// `ks_getc`: the next byte, or -1 at the end of the input.
inline fn getc(stream: Stream&) -> i32 {
    guard stream.fill() else { return -1; }
    let byte = stream.buffer[stream.begin];
    stream.begin += 1;
    byte
}

// `isspace` in the C locale; std's `is_space` leaves out '\v' and '\f'.
inline fn is_c_space(byte: u8) -> bool { byte == ' ' || (byte >= '\t' && byte <= '\r') }

// `ks_getuntil2` with `KS_SEP_SPACE`: append the bytes before the next
// whitespace to text. Returns that whitespace, or -1 if the input ended first.
fn read_word(stream: Stream&, text: u8[>..]&) -> i32 {
    while stream.fill() {
        let unread = stream.buffer[stream.begin..stream.end];
        var length = unread.len;
        for byte, i in unread {
            if is_c_space(byte) { length = i; break; }
        }
        text.append(unread[..length]);
        stream.begin += length + 1;
        if length < unread.len { return unread[length]; }
    }
    -1
}

// `ks_getuntil2` with `KS_SEP_LINE`: append the rest of the line to text,
// without its terminator. Returns whether the input had anything left.
fn read_line(stream: Stream&, text: u8[>..]&) -> bool {
    var got_any = false;
    while stream.fill() {
        got_any = true;
        let unread = stream.buffer[stream.begin..stream.end];
        let newline = bstr::find_byte(unread, '\n');
        let length = if newline >= 0 { newline } else { unread.len };
        text.append(unread[..length]);
        stream.begin += length + 1;
        if newline >= 0 { break; }
    }
    if got_any && text.len > 1 && text[text.len - 1] == '\r' { text.pop(); }
    got_any
}

// What `read` found.
//
// ```goose
// let file, opened = io::open_read("tests/example.fq");
// assert(opened);
// let reader = fastx::Reader { stream: fastx::Stream { file: file, .. } };
// defer io::close(reader.stream.file);
// var records = 0;
// var status = reader.read();
// while status == fastx::Status.Record {
//     records += 1;
//     status = reader.read();
// }
// assert(records == 2 && status == fastx::Status.End);
// ```
enum Status {
    Record,
    // The input has no more records.
    End,
    // A FASTQ record's quality is missing or not as long as its sequence.
    Truncated,
    // The file could not be read.
    Failed,
}

// The most recently read record and the stream it came from; kseq.h's
// `kseq_t`. name is the header up to its first whitespace and comment the
// rest of the header line; qual is empty for a FASTA record. `read` replaces
// all four.
//
// ```goose
// let file, opened = io::open_read("tests/example.fq");
// assert(opened);
// let reader = fastx::Reader { stream: fastx::Stream { file: file, .. } };
// defer io::close(reader.stream.file);
// assert(reader.read() == fastx::Status.Record);
// assert(reader.name == "read1" && reader.comment == "first read");
// assert(reader.seq == "ACGTACGT" && reader.qual == "IIIIHHHH");
// ```
root struct Reader {
    name:    u8[>..] = [],
    comment: u8[>..] = [],
    seq:     u8[>..] = [],
    qual:    u8[>..] = [],
    // The `>` or `@` of the next record when reading a sequence consumed it,
    // otherwise 0.
    last_char: i32 = 0,
    stream:    Stream,
}

// Read the next record into reader, as `kseq_read` does, and return what was
// found. After `Truncated` the fields hold what the record had; after `End`
// or `Failed` they are not meaningful.
//
// ```goose
// let file, opened = io::open_read("tests/example.fq");
// assert(opened);
// let reader = fastx::Reader { stream: fastx::Stream { file: file, .. } };
// defer io::close(reader.stream.file);
// var bases = 0;
// while reader.read() == fastx::Status.Record { bases += reader.seq.len; }
// assert(bases == 12);
// ```
fn read(reader: Reader&) -> Status {
    let stream .= reader.stream;
    if reader.last_char == 0 {
        // Jump to the next header line.
        var byte = stream.getc();
        while byte >= 0 && byte != '>' && byte != '@' { byte = stream.getc(); }
        guard byte >= 0 else { return stream.end_status(); }
        reader.last_char = byte;
    }
    reader.name.clear();
    reader.comment.clear();
    reader.seq.clear();
    reader.qual.clear();

    let delimiter = stream.read_word(reader.name);
    guard delimiter >= 0 || reader.name.len > 0 else { return stream.end_status(); }
    if delimiter != '\n' { stream.read_line(reader.comment); }

    var next = stream.getc();
    while next >= 0 && next != '>' && next != '+' && next != '@' {
        if next != '\n' { // Skip empty lines.
            reader.seq.push(next as! u8);
            stream.read_line(reader.seq);
        }
        next = stream.getc();
    }
    if next == '>' || next == '@' { reader.last_char = next; } // The next record's header has begun.
    guard !stream.failed else { return Status.Failed; }
    guard next == '+' else { return Status.Record; } // FASTA

    // Skip the rest of the '+' line.
    var skipped = stream.getc();
    while skipped >= 0 && skipped != '\n' { skipped = stream.getc(); }
    guard skipped >= 0 else { return Status.Truncated; }
    while stream.read_line(reader.qual) && reader.qual.len < reader.seq.len {}
    guard !stream.failed else { return Status.Failed; }
    reader.last_char = 0; // The next header line is still ahead.
    if reader.seq.len == reader.qual.len { Status.Record } else { Status.Truncated }
}

fn end_status(stream: Stream&) -> Status {
    if stream.failed { Status.Failed } else { Status.End }
}
