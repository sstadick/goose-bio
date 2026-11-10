// Sorted-interval lookup for one sequence, after rust-lapper:
// https://github.com/sstadick/rust-lapper
//
// `new_lapper` sorts a slice of intervals by start, then stop, and indexes it
// in blocks of 32; `find` calls a block with every interval that overlaps a
// query, in sorted order. A Lapper is a view: it holds slices of the
// intervals it was given and of the index storage, so it is a fixed-size
// value that can sit in an array, one per chromosome, and both arrays must
// outlive it.
//
// The index keeps, for every block, its first start, its smallest and largest
// stop, and the largest stop of all blocks up to it, plus a binary max tree
// over the blocks' largest stops and a copy of every start and stop. A query
// binary-searches the running maximum for the first block that can overlap
// it, walks blocks while they start before the query's stop, and uses the
// tree to jump over runs of blocks whose intervals all stop before it.
// Because of the tree, one long interval among short ones does not slow
// queries down.
import std;
import radix_sort;

namespace lapper;

// The number of intervals in a block: one bit each in a u32 mask.
const BLOCK = 32;

// A half-open interval `start..stop` carrying a value. T is an integer type.
//
// ```goose
// let exon = lapper::Interval<i32, u8> { 10, 20, 'a' };
// assert(exon.start == 10 && exon.stop == 20 && exon.val == 'a');
// ```
struct Interval<T, V> {
    // 0-based inclusive start
    start: T,
    // 0-based exclusive end
    stop: T,
    // A stored value.
    val: V,
}

// Whether a and b share a position. Intervals that only touch do not overlap.
//
// ```goose
// let a = lapper::Interval<i32, u8> { 10, 20, 'a' };
// let b = lapper::Interval<i32, u8> { 19, 30, 'b' };
// let c = lapper::Interval<i32, u8> { 20, 30, 'c' };
// assert(lapper::overlap(a, b));
// assert(!lapper::overlap(a, c));
// ```
inline fn overlap<T, V>(a: Interval<T, V>&, b: Interval<T, V>&) -> bool {
    a.start < b.stop && a.stop > b.start
}

// Whether a shares a position with `b_start..b_stop`.
//
// ```goose
// let a = lapper::Interval<i32, u8> { 10, 20, 'a' };
// assert(lapper::overlap(a, 19, 30));
// assert(!lapper::overlap(a, 0, 10));
// ```
inline fn overlap<T, V>(a: Interval<T, V>&, b_start: T, b_stop: T) -> bool {
    a.start < b_stop && a.stop > b_start
}

// One sequence's intervals, sorted and indexed for `find`. Make one with
// `new_lapper`. Its fields are slices of the two arrays passed to
// `new_lapper`, so a Lapper is a plain value: keep one per chromosome in an
// array, and keep both arrays alive and their intervals unchanged while it
// is in use. Everything but intervals is the index.
//
// ```goose
// let chr1: lapper::Interval<i32, u8>[>..] = [
//     lapper::Interval<i32, u8> { 0, 10, 'a' },
//     lapper::Interval<i32, u8> { 5, 15, 'b' },
// ];
// let chr2: lapper::Interval<i32, u8>[>..] = [lapper::Interval<i32, u8> { 100, 110, 'c' }];
// var index: i32[>..] = [];
// var lappers: lapper::Lapper<i32, u8>[>..] = [];
// lappers.push(lapper::new_lapper(chr1, index));
// lappers.push(lapper::new_lapper(chr2, index));
//
// var found: u8[>..] = [];
// lappers[0].find(7, 9) { interval => found.push(interval.val); };
// lappers[1].find(105, 106) { interval => found.push(interval.val); };
// assert(found == "abc");
// ```
struct Lapper<T, V> {
    // The intervals by start, then stop. Block b is intervals[b * BLOCK..(b + 1) * BLOCK].
    intervals:   Interval<T, V>[:],
    // The start and stop of every interval, padded to whole blocks.
    starts:      T[:],
    stops:       T[:],
    // Each block's first start, its smallest stop, and the largest stop of
    // the blocks up to and including it.
    block_start: T[:],
    block_min:   T[:],
    prefix_max:  T[:],
    // A binary tree of the blocks' largest stops with leaves leaves: node i
    // has the children 2 * i and 2 * i + 1, and block b is node leaves + b.
    tree:        T[:],
    leaves:      i64,
}

checkfree fn start_of<T, V>(interval: Interval<T, V>) -> T { interval.start }

// Whether intervals is ordered by start, then stop. Counting instead of stopping at the
// first pair out of order keeps the loop free of branches.
fn in_order<T, V>(intervals: Interval<T, V>[:]) -> bool {
    guard intervals.len > 0 else { return true; }
    var disorder = 0;
    var start = intervals[0].start;
    var stop = intervals[0].stop;
    for &interval in intervals {
        if interval.start < start || (interval.start == start && interval.stop < stop) { disorder += 1; }
        start = interval.start;
        stop = interval.stop;
    }
    disorder == 0
}

// Sort intervals in place by start, then stop, and index them. The index is
// appended to index, 2 T per interval plus 5 T per block of 32; the result
// holds slices of both arguments. One index array can serve any number of
// Lappers. Intervals that are already in order are not sorted again.
//
// ```goose
// let intervals: lapper::Interval<i32, u8>[>..] = [
//     lapper::Interval<i32, u8> { 20, 30, 'c' },
//     lapper::Interval<i32, u8> { 0, 10, 'b' },
//     lapper::Interval<i32, u8> { 0, 5, 'a' },
// ];
// var index: i32[>..] = [];
// let sorted = lapper::new_lapper(intervals, index);
// assert(sorted.intervals.len == 3 && index.len > 0);
// assert(intervals[0].val == 'a' && intervals[1].val == 'b' && intervals[2].val == 'c');
// ```
fn new_lapper<T, V>(intervals: Interval<T, V>[:], index: T[>..]&) -> Lapper<T, V> {
    if !in_order(intervals) {
        radix_sort::radix_sort(intervals, start_of);
        // The radix sort orders by start alone: order equal starts by stop.
        var first = 0;
        for i in 1..intervals.len + 1 {
            if i == intervals.len || intervals[i].start != intervals[first].start {
                if i - first > 1 { sort(intervals[first..i]) { a, b => a.stop < b.stop }; }
                first = i;
            }
        }
    }
    let blocks = (intervals.len + BLOCK - 1) / BLOCK;
    var leaves = 1;
    while leaves < blocks { leaves *= 2; }

    // The padding and the tree's unused leaves hold default values: `find`
    // masks the padding off, and a tree walk that lands on an unused leaf has
    // passed every block.
    let padding = blocks * BLOCK - intervals.len;
    let base = index.len;
    for &interval in intervals { index.push(interval.start); }
    push_n(index, default<T>(), padding);
    for &interval in intervals { index.push(interval.stop); }
    push_n(index, default<T>(), padding + 3 * blocks + 2 * leaves);

    let lanes = blocks * BLOCK;
    let starts: T[:] = index[base..base + lanes];
    let stops: T[:] = index[base + lanes..base + 2 * lanes];
    let rest: T[:] = index[base + 2 * lanes..];
    let block_start: T[:] = rest[..blocks];
    let block_min: T[:] = rest[blocks..2 * blocks];
    let prefix_max: T[:] = rest[2 * blocks..3 * blocks];
    let tree: T[:] = rest[3 * blocks..];

    for b in 0..blocks {
        let run: T[:] = stops[b * BLOCK..min(intervals.len, (b + 1) * BLOCK)];
        var largest = run[0];
        var smallest = run[0];
        for stop in run {
            if stop > largest { largest = stop; }
            if stop < smallest { smallest = stop; }
        }
        block_start[b] = starts[b * BLOCK];
        block_min[b] = smallest;
        tree[leaves + b] = largest;
        prefix_max[b] = if b > 0 && prefix_max[b - 1] > largest { prefix_max[b - 1] } else { largest };
    }
    var i = leaves - 1;
    while i >= 1 {
        tree[i] = if tree[2 * i] > tree[2 * i + 1] { tree[2 * i] } else { tree[2 * i + 1] };
        i -= 1;
    }
    Lapper { intervals, starts, stops, block_start, block_min, prefix_max, tree, leaves }
}

// The first block after b that holds a stop above start, or a number past
// the last block if there is none.
inline fn next_above<T>(tree: T[:], leaves: i64, b: i64, start: T) -> i64 {
    var i = leaves + b;
    while true {
        // Leave every subtree that ends at i for the next one to its right.
        while i & 1 == 1 { i >>= 1; }
        if i == 0 { return leaves; }
        i += 1;
        if tree[i] > start { break; }
    }
    while i < leaves {
        i *= 2;
        if !(tree[i] > start) { i += 1; }
    }
    i - leaves
}

// The number of leading values at or below key in a non-decreasing slice.
inline fn count_up_to<T>(values: T[:], key: T) -> i64 {
    guard values.len > 0 else { return 0; }
    var low = 0;
    var size = values.len;
    while size > 1 {
        let half = size / 2;
        low = if values[low + half - 1] <= key { low + half } else { low };
        size -= half;
    }
    if values[low] <= key { low + 1 } else { low }
}

// Bit i is set where starts[i]..stops[i] overlaps start..stop. Both slices
// are one block; the loop has no branches so that it compiles to vector
// compares.
inline fn overlap_mask<T>(starts: T[:], stops: T[:], start: T, stop: T) -> u32 {
    var mask: u32 = 0;
    for i in 0..BLOCK {
        let below: u32 = if starts[i] < stop { 1 } else { 0 };
        let above: u32 = if stops[i] > start { 1 } else { 0 };
        mask |= (below & above) << i;
    }
    mask
}

// Bit i is set where stops[i] is above start. stops is one block.
inline fn stop_mask<T>(stops: T[:], start: T) -> u32 {
    var mask: u32 = 0;
    for i in 0..BLOCK {
        let bit: u32 = if stops[i] > start { 1 } else { 0 };
        mask |= bit << i;
    }
    mask
}

// Call F with the intervals of run whose bit in mask is set, in order: by
// testing every bit when most are set, else by jumping from set bit to set
// bit.
inline fn each_set<T, V, F: fn(Interval<T, V>&)>(run: Interval<T, V>[:], mask: u32) {
    var bits = mask;
    if popcount(bits as u64) > BLOCK / 2 {
        for &interval in run {
            if bits & 1 != 0 { F(interval); }
            bits >>= 1;
        }
    } else {
        while bits != 0 {
            F(run[ctz(bits as u64)]);
            bits &= bits - 1;
        }
    }
}

// Call F with every interval that overlaps `start..stop`, ordered by start,
// then stop.
//
// ```goose
// let intervals: lapper::Interval<i32, u8>[>..] = [
//     lapper::Interval<i32, u8> { 0, 10, 'a' },
//     lapper::Interval<i32, u8> { 8, 12, 'b' },
//     lapper::Interval<i32, u8> { 20, 30, 'c' },
//     lapper::Interval<i32, u8> { 30, 40, 'd' },
// ];
// var index: i32[>..] = [];
// let exons = lapper::new_lapper(intervals, index);
//
// var found: u8[>..] = [];
// exons.find(9, 25) { interval => found.push(interval.val); };
// assert(found == "abc");
//
// found.clear();
// exons.find(12, 20) { interval => found.push(interval.val); };
// assert(found.len == 0);
// ```
fn find<T, V, F: fn(Interval<T, V>&)>(lapper: Lapper<T, V>&, start: T, stop: T) {
    let intervals = lapper.intervals;
    let firsts = lapper.block_start;
    let blocks = firsts.len;
    // Blocks before b hold no stop above start; block b holds one.
    var b = count_up_to(lapper.prefix_max, start);
    while b < blocks && firsts[b] < stop {
        let first = b * BLOCK;
        if b + 1 < blocks && firsts[b + 1] < stop {
            // Every start in the block is below stop.
            if lapper.block_min[b] > start {
                for &interval in intervals[first..first + BLOCK] { F(interval); }
            } else {
                each_set(intervals[first..first + BLOCK], stop_mask(lapper.stops[first..first + BLOCK], start), F);
            }
        } else {
            // The last block to look at: later blocks start at or after stop.
            var mask = overlap_mask(lapper.starts[first..first + BLOCK], lapper.stops[first..first + BLOCK], start, stop);
            if intervals.len - first < BLOCK { mask &= ((1 as u32) << (intervals.len - first)) - 1; }
            each_set(intervals[first..min(intervals.len, first + BLOCK)], mask, F);
            return;
        }
        b += 1;
        if !(lapper.tree[lapper.leaves + b] > start) {
            b = next_above(lapper.tree, lapper.leaves, b, start);
        }
    }
}
