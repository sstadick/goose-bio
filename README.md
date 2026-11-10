# goose-bio

Bioinformatics data structures and parsers for
[Goose](https://github.com/aardappel/goose): an implicit interval tree, a
block-indexed interval lookup, a FASTA and FASTQ reader, and the radix sort
both interval structures are built with.

**[API reference and examples](API.md)** · [Benchmarks](bench/BIOFAST.md)

## Get started

Declare the package once in your program's root file, then import the modules
you use. Each module has its own namespace, and Goose downloads the package:

```goose
import bio from "github.com/sstadick/goose-bio@main";
import bio.iitree;

fn main() {
    let tree = iitree::IITree<u8[..8]> { .. };
    tree.add("chr1", 10, 20, "geneA");
    tree.add("chr1", 15, 30, "geneB");
    tree.add("chr2", 10, 20, "geneC");
    tree.index();

    var names: u8[>..] = [];
    tree.overlap("chr1", 18, 25) { interval => names.append(interval.value); };
    assert(names == "geneAgeneB");
}
```

```sh
goose run main.goose
goose build main.goose -O2 -o demo
```

| Module | Namespace | Needs |
| --- | --- | --- |
| `bio.iitree` | `iitree` | nothing else |
| `bio.lapper` | `lapper` | nothing else |
| `bio.radix_sort` | `radix_sort` | nothing else |
| `bio.fastx` | `fastx` | goose-io, goose-bstr, goose-simd |

The examples follow `main`. For reproducible builds, pin the import to a
[release](https://github.com/sstadick/goose-bio/releases) tag, as in
`"github.com/sstadick/goose-bio@v0.2.0"`.

## Implicit interval tree

`iitree.goose` is a port of Heng Li's
[cgranges](https://github.com/lh3/cgranges): an interval tree over named
contigs that stores nothing beyond its sorted intervals. An `IITree<V>` holds
`Interval<V> { start, stop, value }`. Add intervals, index once, then pass a
block to `overlap`, which calls it with each overlapping interval in order of
start, as in the program above.

## Lapper

`lapper.goose` is a sorted-interval lookup after
[rust-lapper](https://github.com/sstadick/rust-lapper). A `Lapper<T, V>`
covers one sequence: `new_lapper` sorts a slice of
`Interval<T, V> { start, stop, val }` in place, for any integer coordinate
type `T`, and `find` calls a block with each overlapping interval, ordered by
start, then stop.

```goose
import bio from "github.com/sstadick/goose-bio@main";
import bio.lapper;

fn main() {
    let exons: lapper::Interval<i32, u8>[>..] = [
        lapper::Interval<i32, u8> { 20, 30, 'c' },
        lapper::Interval<i32, u8> { 0, 10, 'a' },
        lapper::Interval<i32, u8> { 8, 12, 'b' },
    ];
    var index: i32[>..] = [];
    let chr1 = lapper::new_lapper(exons, index);

    var found: u8[>..] = [];
    chr1.find(9, 25) { interval => found.push(interval.val); };
    assert(found == "abc");
}
```

The index groups the sorted intervals in blocks of 32 and keeps a copy of
every start and stop, each block's smallest and largest stop, and a binary
max tree over the blocks. A query binary-searches for the first block that
can overlap it, tests each block it visits with vector compares, and uses the
tree to jump over blocks that end before it, so one long interval among short
ones does not slow queries down.

`new_lapper` appends that index, two coordinates per interval plus five per
block, to an array the program owns. A `Lapper` only holds slices of the
intervals and of that array, so one index array can serve a `Lapper` per
chromosome, kept in an ordinary array. Queries are faster than the interval
tree's; the tree stores nothing beyond its intervals and looks contigs up by
name.

## FASTA and FASTQ

`fastx.goose` is a port of the parser in Heng Li's
[kseq.h](https://github.com/lh3/biofast/blob/4177ef282f4e1addc5b20beb7ee2f6435806f844/lib/kseq.h).
A `Reader` reads one record at a time into buffers it reuses, including
records wrapped over several lines. It needs
[goose-io](https://github.com/sstadick/goose-io) and
[goose-bstr](https://github.com/sstadick/goose-bstr), which the program
declares along with goose-simd for goose-bstr:

```goose
import std;
import os;

import io from "github.com/sstadick/goose-io@v0.2.2";
import bstr from "github.com/sstadick/goose-bstr@v0.3.1";
import simd_bytes from "github.com/sstadick/goose-simd@v0.2.0";
import bio from "github.com/sstadick/goose-bio@main";
import bio.fastx;

fn main() {
    let file, opened = io::open_read("reads.fq");
    if !opened { exit(1); }
    let reader = fastx::Reader { stream: fastx::Stream { file: file, .. } };
    defer io::close(reader.stream.file);

    var bases = 0;
    var status = reader.read();
    while status == fastx::Status.Record {
        bases += reader.seq.len;
        status = reader.read();
    }
    if status != fastx::Status.End { exit(1); }
    print(bases);
}
```

The [biofast benchmarks](bench/BIOFAST.md) compare the interval tree and
this reader with C, Rust, Zig, and Mojo ports.

## Radix sort

`radix_sort.goose` provides in-place insertion sort and MSD radix sort for all
signed and unsigned integer widths. Radix sort uses eight-bit digits, falls
back to stable insertion sort for inputs and buckets of at most 64 elements,
and is itself unstable. A keyed overload moves complete values while sorting
by an integer key.

```goose
import bio from "github.com/sstadick/goose-bio@main";
import bio.radix_sort;

fn main() {
    let values: i16[] = [256, -4, 1, 0, -20];
    radix_sort::radix_sort(values);
    let expected: i16[] = [-20, -4, 0, 1, 256];
    assert(values == expected);
}
```

See the [API reference](API.md) for both direct and keyed sorting, and the
[benchmark against C, Rust, Zig, and Mojo ports](bench/RADIX_SORT.md) for
methodology and current performance.

## Unfinished module

`bed.goose` (a BED reader) is an earlier experiment. It is not part of the
documented API and may change or go away.

## Development

Run the implementation tests and executable documentation from the repository
root:

```sh
python3 tests/run.py --cc gcc clang
scripts/docs.sh --check
scripts/doctest.sh
```

`tests/run.py` runs every test with TinyCC and with each named C compiler.
Regenerate `API.md` with `scripts/docs.sh`. Both documentation scripts use
[goose-doc](https://github.com/sstadick/goose-doc); set `GOOSE_DOC` when its
`main.goose` is not at `../goose-doc/main.goose`. The benchmarks are run with
[pixi](https://pixi.sh) from `bench/`, as their reports describe.

CI runs the tests, the documentation check, and the examples on Linux
x86-64/ARM64 (GCC and Clang), macOS ARM64 (Clang), and Windows x86-64 (MSVC
and clang-cl). Pushing a `vMAJOR.MINOR.PATCH` tag publishes a source release
only after the entire matrix passes, attaching an API reference with
commit-pinned links.

## License

MIT; see [LICENSE](LICENSE), which also carries the notices of cgranges and
klib.
