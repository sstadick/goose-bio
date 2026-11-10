import std;
import bio from "..";
import bio.iitree;

struct Query {
    contig: i64,
    start:  i32,
    end:    i32,
}

fn next_random(state: u64&) -> u64 {
    var x = state;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    state = x;
    return x;
}

fn test_empty_tree() {
    let tree = iitree::IITree<i32> { .. };
    tree.index();
    var hits = 0;
    tree.overlap("chr1", 0, 100) { interval => hits += 1; };
    assert(hits == 0);
}

fn test_add_numbers_contigs_and_skips_reversed_intervals() {
    let tree = iitree::IITree<i32> { .. };
    assert(tree.add_contig("chr2") == 0);
    tree.add("chr1", 10, 20, 0);
    tree.add("chr2", 5, 5, 1);   // An empty interval is kept.
    tree.add("chr1", 20, 10, 2); // A reversed one is not.
    assert(tree.add_contig("chr1") == 1);
    assert(tree.contigs.len == 2);
    assert(tree.intervals.len == 2);
}

fn test_index_sorts_and_groups_by_contig() {
    let tree = iitree::IITree<i32> { .. };
    tree.add("chr2", 30, 40, 0);
    tree.add("chr1", 50, 60, 1);
    tree.add("chr2", 10, 20, 2);
    tree.add("chr1", 5, 7, 3);
    tree.add_contig("chr3"); // A contig without intervals.
    tree.index();

    // chr2 was seen first, so its intervals come first, each run by start.
    let values: i32[4] = [2, 0, 3, 1];
    for interval, i in tree.intervals { assert(interval.value == values[i]); }
    assert(tree.intervals[0].start == 10 && tree.intervals[0].stop == 20);
    assert(tree.contigs[0].offset == 0 && tree.contigs[0].count == 2);
    assert(tree.contigs[1].offset == 2 && tree.contigs[1].count == 2);
    assert(tree.contigs[2].root_k == -1);

    var hits = 0;
    tree.overlap("chr3", 0, 100) { interval => hits += 1; };
    tree.overlap("chrUn", 0, 100) { interval => hits += 1; };
    assert(hits == 0);
}

fn test_overlap_is_half_open_and_in_start_order() {
    let tree = iitree::IITree<i32> { .. };
    tree.add("chr1", 30, 40, 2);
    tree.add("chr1", 10, 20, 0);
    tree.add("chr1", 15, 25, 1);
    tree.add("chr2", 10, 20, 3);
    tree.index();

    var values: i32[>..] = [];
    tree.overlap("chr1", 19, 30) { interval => values.push(interval.value); };
    let expected: i32[] = [0, 1];
    assert(values == expected);

    // Touching endpoints do not overlap.
    var hits = 0;
    tree.overlap("chr1", 25, 30) { interval => hits += 1; };
    tree.overlap("chr1", 0, 10) { interval => hits += 1; };
    assert(hits == 0);

    tree.overlap("chr2", 0, 100) { interval => hits += 1; };
    assert(hits == 1);
}

struct Gene {
    name:   u8[..8],
    strand: u8,
}

// A value wider than the key keeps its interval through the sort.
fn test_values_of_another_type_follow_their_intervals() {
    let tree = iitree::IITree<Gene> { .. };
    tree.add("chr1", 300, 400, Gene { name: "third", strand: '+' });
    tree.add("chr1", 100, 200, Gene { name: "first", strand: '-' });
    tree.add("chr1", 150, 350, Gene { name: "second", strand: '+' });
    tree.index();

    var names: u8[>..] = [];
    tree.overlap("chr1", 180, 320) { interval =>
        names.append(interval.value.name);
        names.push(interval.value.strand);
    };
    assert(names == "first-second+third+");
}

// Compare every query against a scan of all intervals. The sizes cover an
// empty contig, trees shallow enough to be scanned outright (height <= 3),
// complete trees, and trees with an incomplete right side.
fn test_overlap_matches_a_linear_scan() {
    let sizes: i64[] = [0, 1, 2, 15, 16, 17, 31, 100, 1000, 1023, 1024, 1025];
    var state: u64 = 0x4D595DF4D0F33173;

    for size in sizes {
        let tree = iitree::IITree<i32> { .. };
        let names: u8[..4][2] = ["chrA", "chrB"];
        // starts[c] and ends[c] hold contig c's intervals; the value of an
        // interval is its index there.
        var starts_a: i32[>..] = [];
        var ends_a:   i32[>..] = [];
        var starts_b: i32[>..] = [];
        var ends_b:   i32[>..] = [];
        for _ in 0..size {
            let start = (next_random(state) % 10000) as i32;
            let end = start + (next_random(state) % 300) as i32;
            if next_random(state) % 4 == 0 {
                tree.add(names[1], start, end, starts_b.len as i32);
                starts_b.push(start);
                ends_b.push(end);
            } else {
                tree.add(names[0], start, end, starts_a.len as i32);
                starts_a.push(start);
                ends_a.push(end);
            }
        }
        tree.index();

        for _ in 0..200 {
            let on_b = next_random(state) % 4 == 0;
            let start = (next_random(state) % 10400) as i32 - 100;
            let end = start + (next_random(state) % 500) as i32;
            let starts: i32[:] = if on_b { starts_b } else { starts_a };
            let ends:   i32[:] = if on_b { ends_b } else { ends_a };

            var expected = 0;
            for i in 0..starts.len {
                if starts[i] < end && start < ends[i] { expected += 1; }
            }

            var hits = 0;
            var previous_start: i32 = -1;
            tree.overlap(names[if on_b { 1 } else { 0 }], start, end) { interval =>
                assert(interval.start == starts[interval.value]);
                assert(interval.stop == ends[interval.value]);
                assert(interval.start < end && start < interval.stop);
                assert(interval.start >= previous_start);
                previous_start = interval.start;
                hits += 1;
            };
            assert(hits == expected);
        }
    }
}

fn main() {
    test_empty_tree();
    test_add_numbers_contigs_and_skips_reversed_intervals();
    test_index_sorts_and_groups_by_contig();
    test_overlap_is_half_open_and_in_start_order();
    test_values_of_another_type_follow_their_intervals();
    test_overlap_matches_a_linear_scan();
    print("iitree tests passed");
}
