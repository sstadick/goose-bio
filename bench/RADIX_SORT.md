# Radix sort benchmark

This benchmark compares ports of the cgranges `KRADIX_SORT_INIT` radix sort
with the original C. Every port generates the same random 16-byte records by
the same PRNG and sorts them by a `u64` key. Each port's executable links the
C reference, which owns its own storage; no records cross the FFI boundary.
The port and the reference run in one process, use the same monotonic clock,
alternate execution order, and exclude resetting the input from each timed
sample. Every warmup is checked for sorted order and both implementations must
produce matching checksums.

| Port | Source | Record cursors | Bucket tables | Tables live in |
|:---|:---|:---|:---|:---|
| `goose` | `radix_sort_bench.goose`, calling `radix_sort.goose` | bounds-checked `i64` indices | separate `begins` and `ends` | each stack frame |
| `rust-safe` | `radix_sort_bench.rs` | bounds-checked `usize` indices | separate `begins` and `ends` | each stack frame |
| `rust-unsafe` | `radix_sort_bench.rs` | raw addresses | begin/end pairs | each stack frame |
| `zig` | `radix_sort_bench.zig` | raw pointers | begin/end pairs | each stack frame |
| `mojo-safe` | `radix_sort_bench.mojo` | bounds-checked `Int` indices into a `Span` | separate `begins` and `ends` | each stack frame |
| `mojo` | `radix_sort_bench.mojo` | `Int` indices into unchecked pointer offsets | separate `begins` and `ends` | each stack frame |
| `mojo-opt` | `radix_sort_bench_optimized.mojo` | raw addresses | begin/end pairs | one heap workspace per sort |

`zig` follows `rs_sort_record` in `radix_sort_reference.c` step by step: every
bucket starts as two pointers to the first record, and counting advances the
end pointer. Every other port zeroes its counts with one `memset` and counts
into them. The prefix sum then leaves slice indices in `goose`, `rust-safe`,
`mojo-safe`, and `mojo`, and record addresses in `rust-unsafe` and `mojo-opt`.

The four index ports also drop two things the C does: the 16-byte begin/end
pair and the loop that restores every `begin` after the
permutation. [Matching C from indices](#matching-c-from-indices) measures
both. `mojo-opt` is the earlier tuned Mojo port, kept unchanged; [its
section](#tuned-mojo-port) describes it.

`rust-unsafe` and `zig` each use the bucket form that measured fastest for
them. With both forms in one Rust executable, zeroed counts are 7% faster at
1K and 17% faster at 16K than pointer buckets, and equal on the two large
cases. In Zig, against an identical reference, zeroed counts are 1–5% slower
than pointer buckets, and zeroed indices with unchecked `data[index]` access
are up to 10% slower.

Run every port whose compiler is on `PATH` from the repository root:

```sh
bench/run_radix_sort_bench.py --rounds 5 --cpu 0
```

`bench/pixi.toml` supplies Zig and Mojo; Goose, `rustc`, and the C compiler
come from `PATH`:

```sh
cd bench && pixi run radix --rounds 5 --cpu 0
```

The script builds the ports in parallel, runs the executables in turn for each
round, and prints one table. Five rounds of all seven ports take about 35
seconds and `--rounds 1` about 15. `--ports goose,rust-safe` selects ports,
`--cc clang` changes the C compiler, and `--cpu` requires Linux `taskset`.

Every port is built for the native CPU at its highest optimization level:
Goose `-O2` (its specialization and inlining pass) with `-O3 -march=native`
appended when the C compiler builds the generated C, the runtime, and the
reference; `rustc -Copt-level=3 -Ctarget-cpu=native -Ccodegen-units=1`;
`zig build-exe -O ReleaseFast -mcpu native`; and `mojo build -O3 --march
native`.

## Ryzen 9 3950X

Linux x86-64, Goose `v2026.3-269-g74c641c5`, rustc 1.97.1, Zig 0.17.0, Mojo
1.2.0, five rounds pinned to one CPU. Port cells are the median ns/item, then
the ratio to the C reference's median in the same executable. The `c` row is
the median over every executable, then the range of the per-executable
medians. Each cell holds 5,000, 1,000, 100, and 30 samples from left to right.

GCC 13.3.0 compiles the reference and Goose's generated C:

| Port        |          1,024 |            16,384 |           262,144 |         1,048,576 | Geomean port/C |
|:------------|---------------:|------------------:|------------------:|------------------:|---------------:|
| c           | 6.78 (6.6-7.0) | 18.83 (18.3-20.8) | 17.73 (17.6-17.8) | 21.02 (20.8-21.6) |                |
| goose       |   6.93 (1.01x) |     18.95 (1.00x) |     17.64 (0.99x) |     21.53 (1.03x) |         1.007x |
| rust-safe   |   5.49 (0.81x) |     17.42 (0.93x) |     15.24 (0.86x) |     19.12 (0.91x) |         0.877x |
| rust-unsafe |   5.35 (0.79x) |     16.82 (0.91x) |     15.32 (0.86x) |     18.46 (0.88x) |         0.861x |
| zig         |   5.63 (0.80x) |     19.45 (0.94x) |     15.18 (0.85x) |     19.01 (0.88x) |         0.866x |
| mojo-safe   |   6.31 (0.93x) |     17.48 (0.92x) |     16.13 (0.91x) |     19.26 (0.92x) |         0.922x |
| mojo        |   5.49 (0.83x) |     15.76 (0.83x) |     15.05 (0.85x) |     18.55 (0.89x) |         0.851x |
| mojo-opt    |   5.78 (0.86x) |     20.69 (1.13x) |     15.16 (0.86x) |     17.70 (0.84x) |         0.913x |

Ubuntu Clang 18.1.3 compiles them (`--cc clang`):

| Port        |          1,024 |            16,384 |           262,144 |         1,048,576 | Geomean port/C |
|:------------|---------------:|------------------:|------------------:|------------------:|---------------:|
| c           | 5.61 (5.5-5.9) | 18.24 (17.9-19.7) | 15.34 (15.2-15.5) | 17.97 (17.9-18.8) |                |
| goose       |   5.25 (0.95x) |     16.26 (0.91x) |     15.23 (1.00x) |     18.71 (1.05x) |         0.974x |
| rust-safe   |   5.47 (0.98x) |     16.85 (0.93x) |     15.18 (0.99x) |     18.52 (1.03x) |         0.983x |
| rust-unsafe |   5.29 (0.94x) |     16.84 (0.93x) |     15.24 (0.99x) |     18.38 (1.03x) |         0.972x |
| zig         |   5.57 (0.94x) |     19.25 (0.98x) |     15.15 (0.99x) |     19.00 (1.01x) |         0.977x |
| mojo-safe   |   6.32 (1.14x) |     17.58 (0.96x) |     16.27 (1.06x) |     20.14 (1.12x) |         1.069x |
| mojo        |   5.28 (0.95x) |     16.21 (0.89x) |     14.94 (0.98x) |     18.29 (1.02x) |         0.960x |
| mojo-opt    |   5.84 (1.02x) |     21.36 (1.17x) |     15.08 (0.97x) |     18.23 (1.01x) |         1.042x |

Against Clang's C, five ports are 2–4% ahead by geomean and within noise of
each other: `mojo`, `rust-unsafe`, `goose`, `zig`, and `rust-safe`. Three of
them keep integer cursors, `goose` and `rust-safe` keep their bounds checks,
and `rust-safe` contains no `unsafe`. Each is 7–11% ahead at 16K except
`zig`, which ties there, and none is more than 5% behind at any size.

`mojo-safe` is 7% behind by geomean, 14% at 1K and 12% at 1M, and ahead only
at 16K. Its gap to `mojo` is mostly how Mojo 1.2.0 compiles a `Span` bounds check,
not the check itself; [Matching C from indices](#matching-c-from-indices)
shows the measurement.

`mojo-opt` is 11% behind `mojo` at 1K and 32% behind at 16K, and within 1% of
it on the two large cases; beside GCC's reference it is 5% ahead at 1M. It
keeps the scan loop and the restore loop that `mojo` dropped, and its heap
workspace does not make up for them.

The C compiler still decides where Goose stands. Through Clang, Goose is 3%
ahead of Clang's C by geomean. Through GCC, Goose is within 3% of GCC's C at
every size, and slower than every LLVM port there except `zig` and
`mojo-opt` at 16K.
GCC's C is itself 16–21% slower than Clang's C at every size but 16K, so the
Rust, Zig, and Mojo ratios below 1 in the first table measure GCC against
LLVM more than they measure those languages.

### Placement noise

Small differences in these tables are not language differences. The `c` row's
range is the same C source, compiled with the same flags, linked into
different executables: at 16K its median varies by up to 14%. Rebuilding an
earlier Rust executable in another directory moved every function by 16
bytes and changed the reference's 16K median in that executable from about
20.0 to 18.5 ns/item, which turned that `rust-safe` from 0.98x to about
1.05x without any change to either sort.

The ports are subject to the same effect. The same `mojo-opt` source and
flags take 17.7 ns/item at 1M beside GCC's reference and 18.2 beside Clang's.
Replacing a signed division in a Mojo port's recursive dispatch with a shift
removed ten instructions per bucket and left the inlined insertion sort's
instructions identical, yet a 4,096-record case, which is one pass and then
256 insertion sorts, slowed from 9.5 to 10.7 ns/item. Whole runs also drift
by about 5% from one invocation to the next.

The 1K and 16K cases are the most sensitive; [the 16K cutoff
cliff](#the-16k-cutoff-cliff) describes what the 16K case stresses. Treat
differences under about 10% at those sizes and under about 5% at the two
large sizes as ties, and compare absolute ns/item across ports before
trusting a ratio.

## Matching C from indices

An index costs a shift that a pointer does not: a record is 16 bytes, x86
addressing scales by at most 8, so `data[index]` needs a shift before the
load. That shift sits on the permutation's dependency chain, where each
swap's load address depends on the record the previous swap loaded. Three
changes, none of which needs a pointer, more than pay for it. The rows below
are successive Mojo variants in one interleaved run beside Clang's reference,
which took about 5.7, 18.3, 15.4, and 18.0 ns/item:

| Mojo variant | 1,024 | 16,384 | 262,144 | 1,048,576 |
|:---|---:|---:|---:|---:|
| Raw addresses, begin/end pairs: the previous `mojo` | 6.02 | 21.36 | 15.74 | 18.12 |
| … with the scan as `for` over `while` | 5.59 | 18.05 | 15.54 | 19.00 |
| … and no restore loop | 5.23 | 16.27 | 15.17 | 18.34 |
| Indices, begin/end pairs, no restore loop | 5.87 | 16.83 | 16.39 | 19.93 |
| Indices, separate `begins` and `ends`: `mojo` | 5.37 | 16.97 | 15.12 | 18.26 |
| The same over a bounds-checked `Span`: `mojo-safe` | 6.28 | 17.59 | 16.24 | 20.11 |

The scan loop. cgranges writes the permutation as one `while` that either
places a record or, when the bucket is full, advances to the next bucket.
Mojo 1.2.0 compiles the advance as `sete` and `add`, so the bucket number
depends on the comparison of two values just loaded from the table, and each
step of the scan over 256 buckets waits about nine cycles for the previous
one. Rust and Clang compile the same source shape to a branch, which the
processor predicts. A call that sorts 70 records still scans all 256
buckets, so this dominated the 16K case, whose top level makes 115 such
calls: timing each phase put Mojo's depth-1 permutation at 2.5 times Rust's,
its dispatch to the next level at 1.25 times, and every other phase within
about 10%. A `for` over buckets with a `while` inside compiles to branches
in Mojo too.

That row also replaces a signed division in the dispatch with a shift and
loads the insertion sort's previous record once, which together were worth
about 0.5 ns/item at 16K. Its 1M column rises although the swap loop's
instructions are identical in both variants; [Placement
noise](#placement-noise) describes a like case.

The restore loop. The permutation advances every bucket's `begin` to its
`end`, and the C then rewrites all 256 `begin`s before it recurses. A bucket
begins where the previous one ends, so the recursion can carry one running
`begin` instead.

Separate tables. With 8-byte entries, `begins[digit]` is one scaled address,
so the shift that the 16-byte pair needed on the same dependency chain goes
away and the index's shift replaces it. This is what recovers the 8–9% that
indices otherwise lose on the two large cases. The counts also need only a
2 KiB `memset`, and the swap loop touches half as many table cache lines.

The same three changes took `rust-safe` from 5.87, 19.14, 16.39, and 19.74
ns/item to 5.50, 16.97, 15.20, and 19.19 in one interleaved run; dropping
only the restore loop accounted for 0.5 ns/item of that at 16K. In
`radix_sort.goose` they took `goose` from 5.62, 17.90, 15.97, and 19.34
ns/item to 5.27, 16.37, 15.25, and 18.71 beside Clang's reference, also in
one interleaved run. Across full table runs that moved `goose` from 4% behind
Clang's C by geomean to 3% ahead, and from 8% behind GCC's C to 1% behind.

`mojo-safe` shows what Mojo's bounds checks cost today. `Array` checks on
the tables disappear, because the digit and the bucket number are provably
below 256. Each `Span` access in the swap loop keeps its check, and Mojo
1.2.0 stores the six words of the failure message's arguments to the stack
before every comparison instead of on the failure path. Writing the
permutation's checks by hand, as an unsigned comparison and `abort()` before
each unchecked access, measured 5.45, 17.10, 15.21, and 19.29 ns/item in a
run where `mojo` took 5.37, 16.99, 15.16, and 18.33 and `mojo-safe` 6.29,
17.60, 16.28, and 20.10. The comparison there costs under 2% at three sizes
and 5% at 1M, which points at how the checks are compiled, not at the
comparisons, for most of `mojo-safe`'s gap.

## Tuned Mojo port

`radix_sort_bench_optimized.mojo`, the `mojo-opt` port, is the earlier tuned
Mojo port, kept unchanged. It was the result of the tuning this section
describes, and `mojo` has since passed it at 1K and 16K; [Matching C from
indices](#matching-c-from-indices) explains why. The tables in this section
were measured with Mojo 1.1.0, before the other ports existed.

The timed wrapper allocates one eight-level heap workspace per sort, and each
recursive depth clears and reuses its own 256-entry bucket table. Buckets
start as zero integer counts; after the prefix sum their cursors become raw
record addresses so the hot permutation advances pointers without bounds
checks or repeated base/index scaling. The implementation also uses
compile-time radix constants, forced inlining for digit extraction and the
public wrapper, and an out-of-line specialized recursive kernel. Input reset is
an unchecked `unsafe_memcpy`. There are no collection bounds checks in the sort.

Build and run it with Mojo 1.1 or newer:

```sh
bench/run_radix_sort_bench.py --ports mojo-opt --rounds 5 --cpu 0
```

The runner compiles Mojo at `-O3 --march native`, compiles the unchanged C
reference at GCC `-O3 -march=native`, links them into one executable, and uses
the same alternating-order, validation, checksum, and excluded-reset protocol
as every other port. A command with arguments can be supplied when Mojo is
managed by Pixi, for example:

```sh
bench/run_radix_sort_bench.py --ports mojo-opt --rounds 5 --cpu 0 \
  --mojo 'pixi run mojo'
```

On the same Ryzen 9 3950X, using Mojo 1.1.0 and GCC 13.3.0:

| Items | C median | Mojo median | Mojo/C | C ns/item | Mojo ns/item | Samples each |
|---:|---:|---:|---:|---:|---:|---:|
| 1,024 | 0.007 ms | 0.006 ms | 0.863x | 6.63 | 5.72 | 5,000 |
| 16,384 | 0.291 ms | 0.334 ms | 1.150x | 17.74 | 20.39 | 1,000 |
| 262,144 | 4.561 ms | 3.915 ms | 0.858x | 17.40 | 14.94 | 100 |
| 1,048,576 | 21.516 ms | 18.178 ms | 0.845x | 20.52 | 17.34 | 30 |

For these inputs Mojo beats the cgranges/GCC reference by about 14–16% on the
two large cases, but loses by about 15% at 16K. This compares whole toolchains,
not just languages: Mojo uses its LLVM pipeline while the reference uses GCC,
and Mojo specializes the fixed radix width directly into its recursive kernel.
The measured recursive bodies were 621 bytes for Mojo and 720 bytes for GCC's
constant-propagated C clone on this build.

To compare frontends on LLVM-family pipelines, the same five-round run was also
made with Ubuntu Clang 18.1.3:

```sh
bench/run_radix_sort_bench.py --ports mojo-opt --rounds 5 --cpu 0 --cc clang
```

| Items | Clang C median | Mojo median | Mojo/C | C ns/item | Mojo ns/item | Samples each |
|---:|---:|---:|---:|---:|---:|---:|
| 1,024 | 0.006 ms | 0.006 ms | 0.989x | 6.17 | 6.10 | 5,000 |
| 16,384 | 0.336 ms | 0.365 ms | 1.085x | 20.53 | 22.28 | 1,000 |
| 262,144 | 4.409 ms | 4.141 ms | 0.939x | 16.82 | 15.80 | 100 |
| 1,048,576 | 20.888 ms | 20.840 ms | 0.998x | 19.92 | 19.87 | 30 |

Against Clang C, Mojo is at parity or faster in three cases and is about 9%
slower at the 16K cutoff cliff. Thus the large-input win above is specifically a win over
this GCC build, not evidence that Mojo beats optimized C independently of the C
compiler. Clang and Mojo may still use different LLVM revisions and pipelines,
so this is closer to a frontend comparison rather than a perfectly controlled
one.

### The 16K cutoff cliff

The 16K case distributes exactly 64 records per top-byte bucket on average. In
the deterministic input used here, 115 of 256 buckets contain more than the
64-record cutoff and therefore clear and process another depth-specific 4 KiB
table; the other 141 use insertion sort directly. A threshold sweep makes the source
of the anomaly visible:

| Cutoff | Recursive top-byte buckets | Clang C ns/item | Mojo ns/item | Mojo/C |
|---:|---:|---:|---:|---:|
| 32 | 256 | 22.74 | 26.67 | 1.173x |
| 64 | 115 | 19.38 | 20.94 | 1.081x |
| 128 | 0 | 14.50 | 14.46 | 0.998x |

The generated standalone insertion loops are nearly identical: 95 bytes for
Clang C and 111 bytes for Mojo, with the same 128-bit record-move loop. Moving
the per-call 4 KiB local `Array` into one caller-owned heap workspace reduced
the matched-cutoff gap from about 21% to 8%. Raising the cutoff to 128 removes
the remaining tiny recursive calls and reaches parity, but the benchmark keeps
64 to match cgranges exactly.

The gap that remained at the matched cutoff was not the cutoff. It was the
scan over the 256 buckets in every one of those small recursive calls, which
Mojo compiled without a branch; [Matching C from
indices](#matching-c-from-indices) describes it. With that loop rewritten,
`mojo` keeps the local `Array` and is 11% ahead of Clang's C at 16K.

Experiments with integer cursors throughout, a second indexed kernel for small
recursive calls, fill batches from 8 through 64, explicit pointer buckets, and
an iterative DFS worklist all failed to improve the matched-cutoff result.
The retained design combines zeroed integer buckets, a heap workspace indexed
by recursion depth, and raw-address cursors in the hot permutation.

## Goose and Mojo for performance work

This comparison predates the all-port tables: Mojo here is `mojo-opt` on Mojo
1.1.0, and Goose's C was compiled with `-O3` but without `-march=native`.

The same algorithm does not produce one universal winner:

| Items | Goose ns/item | Goose/C | Mojo ns/item | Mojo/C | Faster port |
|---:|---:|---:|---:|---:|---:|
| 1,024 | 5.42 | 0.930x | 6.10 | 0.989x | Goose |
| 16,384 | 16.80 | 0.866x | 22.28 | 1.085x | Goose |
| 262,144 | 15.65 | 1.035x | 15.80 | 0.939x | Mojo |
| 1,048,576 | 18.77 | 1.022x | 19.87 | 0.998x | Mojo |

The ratios use each port's separately linked, paired Clang C run, which is more
reliable than comparing absolute times across binaries. Goose obtains excellent
small-input code from safe integer cursors and Clang's mature C optimization,
but cannot express unchecked pointer induction when its bounds analysis stops.
Mojo permits zero-cost unsafe access, raw addresses, SIMD, and representation
experiments, reaching C parity on large inputs, but those controls do not make
optimization automatic: several plausible low-level rewrites regressed this
kernel.

For this workload Goose offers a strong performance floor with substantially
stronger safety and less tuning. Mojo offers the higher controllable ceiling and
the tools needed when the compiler cannot prove an invariant, at the cost of
more code, more benchmarking, and an unsafe correctness burden. These results
should not be generalized to all workloads, especially SIMD, accelerator, or
FFI-heavy kernels where Mojo exposes capabilities Goose intentionally does not.

## Assembly investigation

This section predates the separate `begins` and `ends` tables and the
removal of the restore loop in `radix_sort.goose`; the assembly has not been
compared again since.

The explicit zeroing loop over `ends` is intentional. Declaring a fixed array
without an initializer and filling it in place gives GCC one array. Using a
repeated array value made Goose construct a temporary array before assigning
it to the local, adding avoidable stack initialization.

The remaining assembly is not identical:

- The C reference stores direct begin/end pointers. Goose stores safe `i64`
  cursors, so its cycle loop loads an index, checks it against the slice length,
  scales it by the 16-byte record size, and adds the slice base.
- Goose retains checks for cursor-derived element accesses and recursive slice
  bounds. C has no corresponding checks.
- Goose's recursive body is consequently larger and has different register and
  control flow from the pointer-based C body.

Bounds checks are therefore part of the difference, especially for small
inputs, but the assembly also differs in cursor representation and surrounding
control flow. The benchmark does not isolate one cause as the entire gap.
