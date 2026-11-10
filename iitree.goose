// Implicit interval tree over named contigs.
//
// An `IITree` collects intervals with `add`, sorts and indexes them once with
// `index`, and then calls a block with the intervals that overlap a query in
// `overlap`. The tree is implicit: the sorted array of intervals is its own
// balanced binary tree, so indexing only stores, in every interval, the
// largest stop in its subtree.
//
// The implementation is a port of Heng Li's cgranges:
// https://github.com/lh3/biofast/blob/4177ef282f4e1addc5b20beb7ee2f6435806f844/lib/cgranges.c
import std;
import dictionary;
import radix_sort;

namespace iitree;

// A half-open interval `start..stop` on one contig, carrying a value.
//
// max belongs to the tree: until `index` it is the number of the interval's
// contig, and afterwards the largest stop in the interval's subtree. start
// and max are adjacent so that they load as one sort key.
//
// ```goose
// let tree = iitree::IITree<u8> { .. };
// tree.add("chr1", 10, 20, 'a');
// tree.index();
// tree.overlap("chr1", 15, 16) { interval =>
//     assert(interval.start == 10 && interval.stop == 20 && interval.value == 'a');
// };
// ```
struct Interval<V> {
    start: i32,
    max:   i32,
    stop:  i32,
    value: V,
}

// One contig's part of the tree: `intervals[offset..offset + count]`, a tree
// of height root_k. root_k is -1 for a contig without intervals.
//
// ```goose
// let tree = iitree::IITree<i32> { .. };
// tree.add("chr1", 10, 20, 0);
// tree.add("chr2", 5, 8, 1);
// tree.add("chr2", 6, 9, 2);
// tree.index();
// assert(tree.contigs[1].offset == 1 && tree.contigs[1].count == 2);
// ```
struct Contig {
    root_k: i32,
    offset: i64,
    count:  i64,
}

// Intervals with values of type V, grouped by contig. Construct one with
// `IITree<V> { .. }`, `add` every interval, call `index` once, and then query
// with `overlap`.
//
// Contigs are numbered in the order `add` first sees their names; contig_ids
// maps a name, kept in names, to that number.
//
// ```goose
// let tree = iitree::IITree<i32> { .. };
// tree.add("chr1", 10, 20, 0);
// tree.add("chr1", 30, 40, 1);
// tree.index();
// var hits = 0;
// tree.overlap("chr1", 19, 31) { interval => hits += 1; };
// assert(hits == 2);
// ```
root struct IITree<V> {
    intervals:  Interval<V>[>..]             = [],
    contigs:    Contig[>..]                  = [],
    names:      u8[varint][>..]              = [],
    contig_ids: dictionary<const u8[:], i32> = dictionary<const u8[:], i32> {},
}

// Return the number of the contig called name, adding the contig if it is new.
//
// ```goose
// let tree = iitree::IITree<i32> { .. };
// assert(tree.add_contig("chr1") == 0);
// assert(tree.add_contig("chr2") == 1);
// assert(tree.add_contig("chr1") == 0);
// ```
fn add_contig<V>(tree: IITree<V>&, name: const u8[:]) -> i32 {
    let known = tree.contig_ids.get(name);
    if known { return known; }
    let id = tree.contigs.len as i32;
    tree.contigs.push(Contig { root_k: -1, offset: 0, count: 0 });
    // The dictionary keys are views of the names this tree keeps.
    tree.contig_ids.insert(tree.names.push(name), id);
    id
}

// Add the interval `start..stop` of contig, with its value, to a tree that is
// not yet indexed. An interval whose start is after its stop is ignored.
//
// ```goose
// let tree = iitree::IITree<f64> { .. };
// tree.add("chr1", 10, 20, 0.5);
// tree.add("chr1", 20, 10, 1.5); // Ignored.
// assert(tree.intervals.len == 1 && tree.intervals[0].value == 0.5);
// ```
fn add<V>(tree: IITree<V>&, contig: const u8[:], start: i32, stop: i32, value: V) {
    guard start <= stop;
    let id = tree.add_contig(contig);
    tree.intervals.push(Interval<V> { start: start, max: id, stop: stop, value: value });
}

// Contig, then start: max still holds the contig's number.
checkfree fn sort_key<V>(interval: Interval<V>) -> u64 {
    (interval.max as! u64) << 32 | (interval.start as! u32)
}

// Sort the intervals by contig and start, then build each contig's tree.
// Call it once, after the last `add` and before the first `overlap`.
//
// ```goose
// let tree = iitree::IITree<i32> { .. };
// tree.add("chr1", 30, 40, 1);
// tree.add("chr1", 10, 20, 0);
// tree.index();
// assert(tree.intervals[0].value == 0 && tree.intervals[1].value == 1);
// ```
fn index<V>(tree: IITree<V>&) {
    let intervals: Interval<V>[:] = tree.intervals;
    var sorted = true;
    for i in 1..intervals.len {
        if sort_key(intervals[i - 1]) > sort_key(intervals[i]) { sorted = false; break; }
    }
    if !sorted { radix_sort::radix_sort(intervals, sort_key); }

    // Find each contig's run of intervals.
    var first = 0;
    for i in 1..intervals.len + 1 {
        if i == intervals.len || intervals[i].max != intervals[first].max {
            let contig .= tree.contigs[intervals[first].max];
            contig.offset = first;
            contig.count = i - first;
            first = i;
        }
    }

    // The contig is now implied by the run, which frees max for the tree.
    for &contig in tree.contigs {
        contig.root_k = index_contig(intervals[contig.offset..contig.offset + contig.count]);
    }
}

// Store in every max the largest stop in that interval's subtree and return
// the tree's height, or -1 for an empty tree.
fn index_contig<V>(tree: Interval<V>[:]) -> i32 {
    let count = tree.len;
    if count == 0 { return -1; }
    // last is the largest stop in the, possibly incomplete, rightmost subtree
    // of the current level, whose root would sit at last_i.
    var last_i = 0;
    var last: i32 = 0;
    var i = 0;
    while i < count {
        last_i = i;
        last = tree[i].stop;
        tree[i].max = last;
        i += 2;
    }
    var k: i32 = 1;
    while (1 << k) <= count {
        let x = 1 << (k - 1);
        i = (x << 1) - 1;
        while i < count {
            let left = tree[i - x].max;
            let right = if i + x < count { tree[i + x].max } else { last };
            var largest = tree[i].stop;
            if left > largest { largest = left; }
            if right > largest { largest = right; }
            tree[i].max = largest;
            i += x << 2;
        }
        last_i = if (last_i >> k) & 1 != 0 { last_i - x } else { last_i + x };
        if last_i < count && tree[last_i].max > last { last = tree[last_i].max; }
        k += 1;
    }
    k - 1
}

// One step of the traversal in `overlap`: the subtree rooted at index x with
// height k, and whether its left child has been handled.
struct StackEntry {
    x:         i64,
    k:         i32,
    left_done: bool,
}

// Call F with every interval of contig that overlaps `start..stop`, in order
// of interval start. An unknown contig has no overlaps. The tree must be
// indexed.
//
// ```goose
// let tree = iitree::IITree<i32> { .. };
// tree.add("chr1", 10, 20, 0);
// tree.add("chr1", 15, 25, 1);
// tree.add("chr1", 30, 40, 2);
// tree.add("chr2", 10, 20, 3);
// tree.index();
// var values: i32[>..] = [];
// tree.overlap("chr1", 19, 30) { interval => values.push(interval.value); };
// let expected: i32[] = [0, 1]; // Half-open: 30..40 only touches the query.
// assert(values == expected);
// ```
fn overlap<V, F: fn(Interval<V>&)>(tree: IITree<V>&, contig: const u8[:], start: i32, stop: i32) {
    let id = tree.contig_ids.get(contig);
    guard id;
    let found = tree.contigs[id];
    guard found.root_k >= 0;
    let intervals: Interval<V>[:] = tree.intervals[found.offset..found.offset + found.count];
    let count = intervals.len;

    var stack: StackEntry[64];
    var depth = 1;
    stack[0] = StackEntry { x: (1 << found.root_k) - 1, k: found.root_k, left_done: false };
    while depth > 0 {
        depth -= 1;
        let node = stack[depth];
        if node.k <= 3 {
            // The subtree has at most (1 << (k + 1)) - 1 intervals: scan it.
            let first = node.x >> node.k << node.k;
            var last = first + (1 << (node.k + 1)) - 1;
            if last > count { last = count; }
            for &interval in intervals[first..last] {
                if interval.start >= stop { break; }
                if start < interval.stop { F(interval); }
            }
        } else if !node.left_done {
            // Come back to this node after its left subtree.
            let left = node.x - (1 << (node.k - 1));
            stack[depth] = StackEntry { x: node.x, k: node.k, left_done: true };
            depth += 1;
            if left >= count || intervals[left].max > start {
                stack[depth] = StackEntry { x: left, k: node.k - 1, left_done: false };
                depth += 1;
            }
        } else if node.x < count && intervals[node.x].start < stop {
            if start < intervals[node.x].stop { F(intervals[node.x]); }
            stack[depth] = StackEntry { x: node.x + (1 << (node.k - 1)), k: node.k - 1, left_done: false };
            depth += 1;
        }
    }
}
