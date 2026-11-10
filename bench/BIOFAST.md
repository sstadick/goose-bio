# biofast benchmarks

These are ports of the two [biofast](https://github.com/lh3/biofast)
benchmarks, measured against biofast's own C:

- `bedcov` loads one BED file into an implicit interval tree, then prints,
  for every interval of a second BED file, how many loaded intervals overlap
  it and how many of its bases they cover. The C is `bedcov_c1_cgr.c` over
  `cgranges.c`.
- `fqcnt` parses a FASTA or FASTQ file and prints its record, base, and
  quality counts. The C is `fqcnt_c1_kseq.c` over `kseq.h`.

`biofast/c/` holds the upstream C unchanged. Each port is two command-line
programs that take the same arguments as the C and must print the same bytes:

| Port | Sources | Interval tree | Input | Output |
|:---|:---|:---|:---|:---|
| `c` | `biofast/c/` | `cgranges.c`; overlaps fill an index array | `kstream_t` over zlib | `printf` |
| `goose` | `biofast/bedcov.goose`, `bed.goose`, `fqcnt.goose`, `../iitree.goose`, `../fastx.goose`, `../radix_sort.goose` | the repository's `iitree` module; overlaps arrive through a block | goose-io `read_line` for `bedcov`; the repository's `fastx` module, with its ported `kstream_t`, for `fqcnt` | `format` into a buffer |
| `rust` | `biofast/*.rs` | in `bedcov.rs`; overlaps fill a `Vec<usize>` | `BufReader` | `write!` to a `BufWriter` |
| `zig` | `biofast/*.zig` | in `bedcov.zig`; overlaps fill an `ArrayList` | `std.Io.Reader` | `print` to a buffered `std.Io.Writer` |
| `mojo` | `biofast/*.mojo` | in `bedcov.mojo`; overlaps fill a `List` | a ported `kstream_t` | integers formatted into a `List` buffer |

Every port follows the C's algorithms: the cgranges radix sort, the implicit
tree and its explicit stack, and the `kseq_read` state machine with
multi-line records. Around them each port uses its own language's
containers, hash map, and buffered I/O where the standard library has them.
`rust` is safe Rust on the standard library, `goose` and `mojo` keep their
bounds checks, and `zig` builds as `ReleaseFast`, which removes them. Mojo
1.2.0 has no buffered reader, so `mojo` ports `kstream_t` for both programs.

`goose` differs from the C and the other ports in one call. It passes
`overlap` a block that updates the coverage as each overlap is found, where
they fill an array of indices, as `cr_overlap` does, and then read the
intervals back through it. [Why Goose is ahead](#why-goose-is-ahead-of-c-on-bedcov)
measures what that is worth.

The ports leave out what the benchmark does not use: gzip input,
`bedcov`'s `-c` and `-C` options, and with them `cr_contain`,
`cr_min_start`, and the contig lengths cgranges tracks. Every port parses
numbers with a hand-written `atoi` that matches C `atol` up to 64 bits.

On two kinds of input the ports agree with one another and not with the C.
A number of more than 19 digits wraps in the ports, where `atol` stops at
`LONG_MAX`. A NUL byte ends the line for the C, whose lines are C strings,
and is an ordinary byte to the ports. biofast's data has neither.

## Running

The data is biofast's release archive, 567 MB:

```sh
mkdir -p bench/data && cd bench/data
curl -LO https://github.com/lh3/biofast/releases/download/biofast-data-v1/biofast-data-v1.tar.gz
tar xzf biofast-data-v1.tar.gz
```

`bench/pixi.toml` supplies Zig and Mojo; Goose, `rustc`, and the C compiler
come from `PATH`:

```sh
cd bench && pixi run biofast --rounds 3 --cpu 4
```

Three rounds of all five ports take about three minutes and `--rounds 1`
about 80 seconds. `--ports goose,zig` selects ports (C always
runs), `--cases g2r,fqcnt` selects cases, `--cc clang` changes the C
compiler, `--data` names another data directory, and `--cpu` requires Linux
`taskset`.

The script builds every port, then runs each case once per port and compares
the output's MD5 with C's from the same run; any difference stops the run.
It then times the cases in rounds, one port after another, with standard
output sent to `/dev/null`. A sample is the process's user plus system CPU
seconds and its peak resident set size, read by `measure.c` from `wait4`.
Python does not start the programs itself because Linux carries the peak
resident set size across `exec`, so a child of the Python process never
reports less than Python's own 44 MB.

Every port is built for the native CPU at its highest optimization level:
`cc -O3 -march=native`; Goose `-O2` with the same C flags for its generated
C; `rustc -Copt-level=3 -Ctarget-cpu=native -Ccodegen-units=1`; `zig
build-exe -O ReleaseFast -mcpu native`; and `mojo build -O3 --march native`.

| Case | Command | Work |
|:---|:---|:---|
| `g2r` | `bedcov ex-rna.bed ex-anno.bed` | index 8.9M RNA-seq intervals, query 1.2M annotations |
| `r2g` | `bedcov ex-anno.bed ex-rna.bed` | index 1.2M annotations, query 8.9M RNA-seq intervals |
| `fqcnt` | `fqcnt M_abscessus_HiSeq.fq` | parse 5.7M FASTQ records, 1.3 GB |

## Ryzen 9 3950X

Linux x86-64, Goose `v2026.3-269-g74c641c5`, rustc 1.97.1, Zig 0.17.0, Mojo
1.2.0, three rounds pinned to one CPU. Cells are the median CPU seconds, then
the ratio to C; below 1 is faster than C.

GCC 13.3.0 compiles the C and Goose's generated C:

| Port  |          g2r |          r2g |        fqcnt | Geomean port/C |
|:------|-------------:|-------------:|-------------:|---------------:|
| c     | 2.96 (1.00x) | 5.77 (1.00x) | 0.52 (1.00x) |                |
| goose | 2.01 (0.68x) | 3.40 (0.59x) | 0.56 (1.08x) |         0.756x |
| rust  | 2.79 (0.94x) | 5.33 (0.92x) | 0.56 (1.09x) |         0.983x |
| zig   | 2.49 (0.84x) | 3.89 (0.67x) | 0.53 (1.02x) |         0.834x |
| mojo  | 3.62 (1.23x) | 5.66 (0.98x) | 0.63 (1.22x) |         1.137x |

Ubuntu Clang 18.1.3 compiles them (`--cc clang`):

| Port  |          g2r |          r2g |        fqcnt | Geomean port/C |
|:------|-------------:|-------------:|-------------:|---------------:|
| c     | 2.96 (1.00x) | 5.65 (1.00x) | 0.53 (1.00x) |                |
| goose | 1.96 (0.66x) | 3.28 (0.58x) | 0.48 (0.92x) |         0.707x |
| rust  | 2.81 (0.95x) | 5.29 (0.94x) | 0.56 (1.06x) |         0.982x |
| zig   | 2.51 (0.85x) | 3.89 (0.69x) | 0.53 (1.00x) |         0.836x |
| mojo  | 3.63 (1.23x) | 5.66 (1.00x) | 0.63 (1.20x) |         1.138x |

Peak resident set size in MB, the same in both runs to within 0.5 MB:

| Port  |   g2r |  r2g | fqcnt |
|:------|------:|-----:|------:|
| c     | 139.2 | 20.2 |   1.5 |
| goose | 138.0 | 19.8 |   1.5 |
| rust  | 139.5 | 20.2 |   1.8 |
| zig   | 139.5 | 20.2 |   1.8 |
| mojo  | 300.5 | 56.8 |   9.0 |

## Why Goose is ahead of C on `bedcov`

Goose compiles to C and the same compiler builds it and the reference, so
the lead is not code generation. Four programs that are not in this
repository separate what it is: two patched copies of `bedcov_c1_cgr.c`,
`goose` with every bounds check kept, and `goose-cgr`. `goose-cgr` is `goose`
with an `overlap` that fills an `i64` array with the overlaps' indices, as
`cr_overlap` does, and a coverage loop that reads the intervals back through
them, as the C's does. `iitree` had that overload when these were measured
and has since dropped it, and its `Interval` then packed start and end into
one `u64`; the `goose` row here predates both changes. CPU seconds, the
median of three pinned runs:

| Program | g2r, GCC | r2g, GCC | g2r, Clang | r2g, Clang |
|:---|---:|---:|---:|---:|
| `c` | 2.91 | 5.64 | 2.93 | 5.59 |
| C with `atol` replaced by a digit loop | 2.63 | 5.24 | | |
| … and `printf` by integers formatted by hand into a 64 KiB buffer | 2.51 | 3.97 | 2.38 | 3.87 |
| `goose-cgr` | 2.50 | 3.93 | 2.22 | 3.59 |
| `goose` | 2.08 | 3.59 | 1.89 | 3.35 |
| `goose` built with `--no-bce`, which keeps every bounds check | 2.06 | 3.63 | | |

libc is the larger part: 0.4 s of `g2r` and 1.7 s of `r2g` through GCC.
`printf` interprets its format string for every line and `strtol` handles
bases, overflow, and the locale, and `r2g` prints 8.9M lines. Goose's
`format` takes typed arguments, so there is no format string to interpret,
and every port's `atoi` is a digit loop.

With libc out of the way the C and `goose-cgr` are within 1% through GCC, on
the same tree walk and the same radix sort. `goose-cgr` gets there with its
bounds checks: Goose reports that it keeps 9 of the 12 checks in the tree
walk, and keeping every check in the program with `--no-bce` changes nothing
measurable.

Through Clang `goose-cgr` is 6–7% ahead of the patched C, and that is libc
again. zlib is a small part: the C reads plain files through `gzread`, and
calling `read` instead takes 1–3% off it. The rest is what the C does with
each line. `ks_getuntil` copies the line into a `kstring_t`, `parse_bed`
ends its fields with NULs, and khash and the output then call `strcmp` and
`strlen` on them. Goose's `read_line` returns a view of the reader's buffer
and a view knows its length, so none of that happens. With `gzread` already
replaced, libc is 14% of the C's samples while it loads `ex-rna.bed` and 6%
of a whole `r2g` run; loading alone takes the C 0.93 s and `goose-cgr`
0.79 s, which is the whole `g2r` difference. Through GCC the same C is 1–2%
ahead of `goose-cgr`; in that build Goose's copy and compare helpers are
calls, 2% of the samples, and in the Clang build they are inlined.

The block is the rest: 0.2–0.4 s in both cases. Both cases report the same
305M overlaps, 256 per query on `g2r` and 34 on `r2g`, so the index array
costs about 1 ns per overlap to store, bounds-check, and load back. Goose
specializes `overlap` for the caller's block and inlines the coverage
update into the traversal, so `goose` never stores an index. C can do the
same with a macro or an inlined callback; cgranges' interface does not.

`zig` fills an array as the C does and matches `goose-cgr` through GCC, 2.49
and 3.89 s against 2.50 and 3.93.

## The other ports

`rust` is 5–8% ahead of C on `bedcov` and takes 36–47% longer on `r2g` than
`goose-cgr`, which fills the same array, did. A sampling profile gives Rust's
tree walk and coverage loop about as many samples as `goose-cgr`'s. About two
thirds of the gap is in two standard library defaults: `core::fmt` and the
`BufWriter` adapter take 15% of Rust's `r2g` samples, three times what
Goose's formatter takes, and `HashMap`'s SipHash takes 4%, about ten times
what Goose spends hashing. The port keeps both because it is meant to be the
Rust a reader would write first.

`fqcnt` is the closer comparison, because the C there is one tight
`kseq_read` loop over `memchr`. `zig` is within 2% of C, `rust` is 6–9%
behind, and `goose` is 8% behind through GCC and 8% ahead through Clang.

`mojo` is level with C on `r2g` and 20–23% behind on `g2r` and `fqcnt`:

- Mojo 1.2.0 stores the arguments of a bounds-check failure message before
  every checked index, as [the radix sort
  notes](RADIX_SORT.md#matching-c-from-indices) describe. Scanning the read
  buffer with `for byte in buffer[begin:end]` and a slice search in place of
  `buffer[i]` took `fqcnt` from 1.00 to 0.47 user seconds.
- Building an earlier `bedcov.mojo` with `-D ASSERT=none`, which removes
  every bounds check, measured 3.16 s on `g2r` and 5.34 s on `r2g` in single
  runs. Without any checks `mojo` was still behind C on `g2r`.
- `mojo` peaks at twice everyone's memory on `g2r`. The likely cause, which
  these notes did not confirm, is `List` growth: growing allocates a new
  buffer and copies, so the old and new buffers of the 8.9M-interval list
  are resident together, where the other ports peak at the list's final
  size.

`zig` is the only port without bounds checks, so its lead over `rust` and
`mojo` includes whatever those checks cost. `goose-cgr` matches it with its
checks on.

A few percent in these tables is noise: [Placement
noise](RADIX_SORT.md#placement-noise) in the radix sort notes describes how
much the same code moves from one executable to the next.

## A bug in the C

Comparing each port with the C on random BED files found that
`bedcov_c1_cgr.c` corrupts its heap when a start coordinate is negative.
`cr_add` packs `(uint64_t)contig << 32 | start` with a signed 32-bit start,
which sign-extends and overwrites the contig number, and indexing then
writes through the wrong contig. The ports keep the low 32 bits of the
start and agree with one another on such input. biofast's data has no
negative coordinates.
