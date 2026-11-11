// Random access to uncompressed FASTA through a standard five-column .fai.
namespace faidx;
import std;
import dictionary;
import io from "github.com/sstadick/goose-io@v0.2.2";
import bstr from "github.com/sstadick/goose-bstr@v0.3.1";

extern fn seek_raw(fd: i32, offset: i64) -> i64 = "goose_bio_faidx_seek";

// Results of opening an index or fetching bases. Queries do not clamp ranges.
// InvalidIndex includes malformed numbers, duplicate names and offset overflow;
// InvalidData means a non-FASTA header or a queried line separator mismatch.
// UnsupportedFormat means gzip/BGZF input; use an uncompressed FASTA instead.
// Errors are not sticky: another query may succeed. A failed fetch clears output.
//
// ```goose
// let reader = faidx::Reader { .. };
// var bases: u8[>..] = [];
// assert(reader.fetch("chr1", 0, 1, bases) == faidx::Status.Closed);
// assert(bases.len == 0);
// ```
enum Status { Ok, Closed, IoError, InvalidIndex, UnsupportedFormat, UnknownReference, OutOfRange, Truncated, InvalidData }

// One FAI entry, in index order. length counts bases; offset and line_width
// count file bytes. line_bases excludes the line terminator. name/name_len
// address the reader's name slab; use reference_name to borrow a name.
// Treat all metadata as read-only; next belongs to the name lookup table.
//
// ```goose
// let reader = faidx::Reader { .. };
// defer reader.close();
// assert(reader.open("tests/example.fa") == faidx::Status.Ok);
// assert(reader.reference_name(0) == "chr1");
// assert(reader.references[0].length == 13);
// ```
struct Reference { name: i64, name_len: i64, length: i64, offset: i64, line_bases: i64, line_width: i64, next: i64 }

// Owns the FASTA file, parsed index and a 64 KiB read cache. Construct with
// Reader { .. } and close explicitly, normally with defer. No whole sequence
// is loaded on open. Use a separate reader per worker; do not share or copy an
// open reader. Query outputs belong to the caller and survive subsequent queries.
// Implementation fields other than references are not a public mutation API.
//
// ```goose
// let reader = faidx::Reader { .. };
// defer reader.close();
// assert(reader.open("tests/example.fa") == faidx::Status.Ok);
// var bases: u8[>..] = [];
// assert(reader.fetch("chr1", 3, 10, bases) == faidx::Status.Ok);
// assert(bases == "TacgtnR");
// ```
root struct Reader {
    file: io::FileHandle = io::FileHandle { fd: -1 },
    references: Reference[>..] = [],
    names: u8[>..] = [],
    ids: dictionary<u64, i64> = dictionary<u64, i64> {},
    buffer: u8[65536],
    buffer_offset: i64 = 0,
    buffered: i64 = 0,
}

// Close the file and clear the index/cache. Repeated close is safe. Return
// false only if closing the OS handle failed; it is invalidated either way.
//
// ```goose
// let reader = faidx::Reader { .. };
// assert(reader.open("tests/example.fa") == faidx::Status.Ok);
// assert(reader.close());
// assert(reader.close());
// assert(reader.references.len == 0);
// ```
fn close(reader: Reader&) -> bool {
    reader.ids.clear();
    reader.references.clear();
    reader.names.clear();
    reader.buffered = 0;
    io::close(reader.file)
}

// Open FASTA and an existing index; an omitted index_path means path + ".fai".
// No index is generated or modified. Failure closes the file and clears any
// previous index. Accept LF/CRLF index lines, exactly five tab-separated fields,
// unique nonempty names and nonnegative signed-64-bit values. FASTQ indexes
// and compressed FASTA are unsupported. The index must match the unchanged
// FASTA: FAI has no identity checksum, and open does not scan all sequence data.
//
// ```goose
// let reader = faidx::Reader { .. };
// defer reader.close();
// assert(reader.open("tests/example.fa", "tests/example.fa.fai") == faidx::Status.Ok);
// assert(reader.open("") == faidx::Status.IoError);
// assert(reader.references.len == 0);
// ```
fn open(reader: Reader&, path: const u8[:], index_path: const u8[:] = "") -> Status {
    reader.close();
    let file, ok = io::open_read(path);
    if !ok { return Status.IoError; }
    reader.file = file;
    file.fd = -1;
    let status = reader.open_index(path, index_path);
    if status != Status.Ok { reader.close(); }
    status
}

fn number(bytes: const u8[:]) -> i64, bool {
    if bytes.len == 0 { return 0, false; }
    var n: i64 = 0;
    for b in bytes {
        if b < '0' || b > '9' { return 0, false; }
        let digit = (b - '0') as i64;
        if n > (9223372036854775807 - digit) / 10 { return 0, false; }
        n = n * 10 + digit;
    }
    return n, true;
}

fn open_index(reader: Reader&, path: const u8[:], index_path: const u8[:]) -> Status {
    // Probe without changing the cache; the first query seeks to its FAI offset.
    var magic: u8[2] = [0; 2];
    var count: i64 = 0;
    while count < 2 {
        let n = io::read(reader.file, magic[count..]);
        if n < 0 { return Status.IoError; }
        if n == 0 { break; }
        count += n;
    }
    if count == 2 && magic[0] == 31 && magic[1] == 139 { return Status.UnsupportedFormat; }
    if count == 0 || magic[0] != '>' { return Status.InvalidData; }
    var automatic: u8[>..] = [];
    if index_path.len == 0 { automatic.append(path); automatic.append(".fai"); }
    let file, ok = io::open_read(if index_path.len == 0 { automatic } else { index_path });
    if !ok { return Status.IoError; }
    let input = io::buffered_reader(file);
    defer io::close(input);
    loop {
        let line, status = io::read_line(input);
        if status == io::LineStatus.Eof { break; }
        if status != io::LineStatus.Line { return Status.IoError; }
        var fields: const u8[:][..5] = [];
        var too_many = false;
        bstr::split_str(line, "\t") { field =>
            if fields.len == 5 { too_many = true; } else { fields.push(field); }
        };
        if too_many || fields.len != 5 || fields[0].len == 0 { return Status.InvalidIndex; }
        for b in fields[0] { if b <= 32 || b == 127 { return Status.InvalidIndex; } }
        if reader.reference_id(fields[0]) >= 0 { return Status.InvalidIndex; }
        var values: i64[4] = [0; 4];
        for i in 4 {
            let n, valid = number(fields[i + 1]);
            if !valid { return Status.InvalidIndex; }
            values[i] = n;
        }
        let length = values[0]; let offset = values[1];
        let bases = values[2]; let width = values[3];
        if length == 0 {
            if bases != 0 || width != 0 { return Status.InvalidIndex; }
        } else {
            if bases == 0 || bases > length || width < bases || width - bases > 2 ||
               (width == bases && length > bases) { return Status.InvalidIndex; }
            let column = (length - 1) % bases;
            if offset > 9223372036854775807 - column ||
               (length - 1) / bases > (9223372036854775807 - offset - column) / width {
                return Status.InvalidIndex;
            }
        }
        // Store offsets, not self-borrowed dictionary keys, so reopening can
        // clear the names safely. Compare names within each hash collision chain.
        let key = hash(fields[0]);
        let head = reader.ids.get(key);
        let next = if head { head } else { -1 };
        reader.ids.insert(key, reader.references.len);
        reader.references.push(Reference { name: reader.names.len, name_len: fields[0].len,
                                          length: length, offset: offset,
                                          line_bases: bases, line_width: width, next: next });
        reader.names.append(fields[0]);
    }
    if reader.references.len == 0 { Status.InvalidIndex } else { Status.Ok }
}

// Look up a reference by its exact FAI name. Return -1 when absent or closed.
// IDs are zero-based index order and remain valid only until close/open.
//
// ```goose
// let reader = faidx::Reader { .. };
// defer reader.close();
// assert(reader.open("tests/example.fa") == faidx::Status.Ok);
// assert(reader.reference_id("chr1") == 0);
// assert(reader.reference_id("missing") == -1);
// ```
fn reference_id(reader: Reader&, name: const u8[:]) -> i64 {
    let head = reader.ids.get(hash(name));
    var id = if head { head } else { -1 };
    while id >= 0 {
        if reader.reference_name(id) == name { return id; }
        id = reader.references[id].next;
    }
    -1
}

// Borrow the name of an existing reference ID until close/open. id must be
// in [0, references.len). Names are the exact first-column FAI bytes.
//
// ```goose
// let reader = faidx::Reader { .. };
// defer reader.close();
// assert(reader.open("tests/example.fa") == faidx::Status.Ok);
// assert(reader.reference_name(1) == "chr2");
// ```
fn reference_name(reader: const Reader&, id: i64) -> const u8[:] {
    let ref = reader.references[id];
    reader.names[ref.name..ref.name + ref.name_len]
}

// Replace output with bases in zero-based [begin, end), preserving case and
// ambiguity codes and removing line endings. Empty intervals, including at
// sequence end, succeed. Negative, reversed or past-end ranges fail rather
// than clamp. Output is cleared on any failure. Its capacity is reused; no
// borrowed sequence views escape. The FASTA must not change while open.
// LF, CRLF and CR sequence lines are supported, as is a final unterminated line.
//
// ```goose
// let reader = faidx::Reader { .. };
// defer reader.close();
// assert(reader.open("tests/example.fa") == faidx::Status.Ok);
// var bases: u8[>..] = [];
// assert(reader.fetch("chr1", 3, 10, bases) == faidx::Status.Ok);
// assert(bases == "TacgtnR");
// assert(reader.fetch(0, 13, 13, bases) == faidx::Status.Ok && bases.len == 0);
// assert(reader.fetch(0, 0, 14, bases) == faidx::Status.OutOfRange);
// ```
fn fetch(reader: Reader&, name: const u8[:], begin: i64, end: i64, output: u8[>..]&) -> Status {
    reader.fetch(reader.reference_id(name), begin, end, output)
}

// ID overload of fetch; cache reference_id outside a query loop when convenient.
//
// ```goose
// let reader = faidx::Reader { .. };
// defer reader.close();
// assert(reader.open("tests/example.fa") == faidx::Status.Ok);
// let id = reader.reference_id("chr2");
// var bases: u8[>..] = [];
// assert(reader.fetch(id, 1, 4, bases) == faidx::Status.Ok && bases == "TGC");
// ```
fn fetch(reader: Reader&, id: i64, begin: i64, end: i64, output: u8[>..]&) -> Status {
    output.clear();
    if reader.file.fd < 0 { return Status.Closed; }
    if id < 0 || id >= reader.references.len { return Status.UnknownReference; }
    let ref = reader.references[id];
    if begin < 0 || end < begin || end > ref.length { return Status.OutOfRange; }
    if begin == end { return Status.Ok; }
    output.resize(end - begin, 0);
    let status = reader.read_bases(ref, begin, output);
    if status != Status.Ok { output.clear(); }
    status
}

fn fill(reader: Reader&, offset: i64) -> Status {
    if offset >= reader.buffer_offset && offset - reader.buffer_offset < reader.buffered { return Status.Ok; }
    reader.buffered = 0;
    if seek_raw(reader.file.fd, offset) != offset { return Status.IoError; }
    let n = io::read(reader.file, reader.buffer);
    if n < 0 { return Status.IoError; }
    if n == 0 { return Status.Truncated; }
    reader.buffer_offset = offset;
    reader.buffered = n;
    Status.Ok
}

fn read_bases(reader: Reader&, ref: Reference, begin: i64, output: u8[:]) -> Status {
    var column = begin % ref.line_bases;
    var offset = ref.offset + (begin / ref.line_bases) * ref.line_width + column;
    var written: i64 = 0;
    while written < output.len {
        let status = reader.fill(offset);
        if status != Status.Ok { return status; }
        let at = offset - reader.buffer_offset;
        let n = min(output.len - written, min(ref.line_bases - column, reader.buffered - at));
        output[written..written + n].move_from(reader.buffer[at..at + n]);
        written += n;
        if written == output.len { break; }
        offset += n;
        column += n;
        if column == ref.line_bases {
            let terminator = ref.line_width - ref.line_bases;
            for i in terminator {
                let status = reader.fill(offset);
                if status != Status.Ok { return status; }
                let b = reader.buffer[offset - reader.buffer_offset];
                if terminator == 2 {
                    if b != (if i == 0 { '\r' } else { '\n' }) { return Status.InvalidData; }
                } else if b != '\n' && b != '\r' { return Status.InvalidData; }
                offset += 1;
            }
            column = 0;
        }
    }
    Status.Ok
}
