import std;
import bio from "..";
import bio.lapper;

fn test_overlap() {
    let a = lapper::Interval { 10, 20, 0 };
    let overlaps_end = lapper::Interval { 15, 25, 0 };
    let touches_start = lapper::Interval { 0, 10, 0 };

    assert(lapper::overlap(a, overlaps_end));
    assert(lapper::overlap(a, 5, 15));
    assert(lapper::overlap(a, 12, 18));
    assert(lapper::overlap(a, 5, 25));
    assert(lapper::overlap(a, 15, 25));

    // Intervals are half-open: touching endpoints do not overlap.
    assert(!lapper::overlap(a, touches_start));
    assert(!lapper::overlap(a, 20, 30));
}

fn test_new_lapper_sorts_by_start_then_stop() {
    let intervals: lapper::Interval<i64, i64>[>..] = [
        lapper::Interval { 5, 7, 3 },
        lapper::Interval { 0, 5, 0 },
        lapper::Interval { 4, 6, 2 },
        lapper::Interval { 0, 3, 1 },
    ];
    var storage: i64[>..] = [];
    let index = lapper::new_lapper(intervals, storage);

    assert(index.intervals.len == 4);
    assert(index.intervals[0].start == 0 && index.intervals[0].stop == 3);
    assert(index.intervals[1].start == 0 && index.intervals[1].stop == 5);
    assert(index.intervals[2].start == 4 && index.intervals[2].stop == 6);
    assert(index.intervals[3].start == 5 && index.intervals[3].stop == 7);
}

fn test_empty_lapper() {
    let intervals: lapper::Interval<i64, i64>[>..] = [];
    var storage: i64[>..] = [];
    let index = lapper::new_lapper(intervals, storage);

    var hits = 0;
    lapper::find(index, 0, 100) { iv => hits += 1; };
    assert(hits == 0);
}

// One index array serves several Lappers, which are plain values in an array.
fn test_lappers_share_index_storage() {
    let first: lapper::Interval<i64, i64>[>..] = [lapper::Interval { 0, 10, 0 }, lapper::Interval { 5, 15, 1 }];
    let second: lapper::Interval<i64, i64>[>..] = [lapper::Interval { 100, 110, 2 }];
    var storage: i64[>..] = [];
    var lappers: lapper::Lapper<i64, i64>[>..] = [];
    lappers.push(lapper::new_lapper(first, storage));
    lappers.push(lapper::new_lapper(second, storage));

    var found: i64[>..] = [];
    lappers[0].find(7, 105) { iv => found.push(iv.val); };
    assert(found == [0, 1]);

    found.clear();
    lappers[1].find(7, 105) { iv => found.push(iv.val); };
    assert(found == [2]);
}

fn test_nonoverlapping_boundaries() {
    let intervals: lapper::Interval<i64, i64>[>..] = [
        lapper::Interval { 0, 10, 0 },
        lapper::Interval { 20, 30, 1 },
        lapper::Interval { 40, 50, 2 },
        lapper::Interval { 60, 70, 3 },
        lapper::Interval { 80, 90, 4 },
    ];
    var storage: i64[>..] = [];
    let index = lapper::new_lapper(intervals, storage);

    var found: i64[>..] = [];
    lapper::find(index, 15, 20) { iv => found.push(iv.val); };
    assert(found.len == 0); // query stop equals interval start

    found.clear();
    lapper::find(index, 30, 35) { iv => found.push(iv.val); };
    assert(found.len == 0); // query start equals interval stop

    found.clear();
    lapper::find(index, 15, 25) { iv => found.push(iv.val); };
    assert(found == [1]); // overlaps the interval start

    found.clear();
    lapper::find(index, 25, 35) { iv => found.push(iv.val); };
    assert(found == [1]); // overlaps the interval stop
}

fn test_containment_queries() {
    let intervals: lapper::Interval<i64, i64>[>..] = [
        lapper::Interval { 0, 10, 0 },
        lapper::Interval { 20, 30, 1 },
        lapper::Interval { 40, 50, 2 },
    ];
    var storage: i64[>..] = [];
    let index = lapper::new_lapper(intervals, storage);

    var found: i64[>..] = [];
    lapper::find(index, 22, 27) { iv => found.push(iv.val); };
    assert(found == [1]); // interval envelops query

    found.clear();
    lapper::find(index, 15, 35) { iv => found.push(iv.val); };
    assert(found == [1]); // query envelops interval

    found.clear();
    lapper::find(index, -5, 55) { iv => found.push(iv.val); };
    assert(found == [0, 1, 2]);
}

fn test_overlapping_intervals() {
    let intervals: lapper::Interval<i64, i64>[>..] = [
        lapper::Interval { 0, 15, 0 },
        lapper::Interval { 10, 25, 1 },
        lapper::Interval { 20, 35, 2 },
        lapper::Interval { 30, 45, 3 },
    ];
    var storage: i64[>..] = [];
    let index = lapper::new_lapper(intervals, storage);

    var found: i64[>..] = [];
    lapper::find(index, 8, 20) { iv => found.push(iv.val); };
    assert(found == [0, 1]);

    found.clear();
    lapper::find(index, 14, 31) { iv => found.push(iv.val); };
    assert(found == [0, 1, 2, 3]);
}

// This shape inspired rust-lapper's "bad lapper": one long interval pushes
// the lower-bound search far behind a query, while shorter intervals must
// still be returned correctly.
fn test_find_behind_first_match() {
    let intervals: lapper::Interval<i64, i64>[>..] = [
        lapper::Interval { 70, 120, 9 },
        lapper::Interval { 10, 15, 0 },
        lapper::Interval { 10, 15, 1 },
        lapper::Interval { 12, 15, 2 },
        lapper::Interval { 14, 16, 3 },
        lapper::Interval { 40, 45, 4 },
        lapper::Interval { 50, 55, 5 },
        lapper::Interval { 60, 65, 6 },
        lapper::Interval { 68, 71, 7 },
        lapper::Interval { 70, 75, 8 },
    ];
    var storage: i64[>..] = [];
    let index = lapper::new_lapper(intervals, storage);

    var found: i64[>..] = [];
    lapper::find(index, 50, 55) { iv => found.push(iv.val); };
    assert(found == [5]);

    found.clear();
    lapper::find(index, 70, 72) { iv => found.push(iv.val); };
    assert(found == [7, 8, 9]);
}

fn test_find_overlaps_in_large_intervals() {
    let intervals: lapper::Interval<i64, i64>[>..] = [
        lapper::Interval { 0, 8, 0 },
        lapper::Interval { 1, 10, 1 },
        lapper::Interval { 2, 5, 2 },
        lapper::Interval { 3, 8, 3 },
        lapper::Interval { 4, 7, 4 },
        lapper::Interval { 5, 8, 5 },
        lapper::Interval { 8, 8, 6 },
        lapper::Interval { 9, 11, 7 },
        lapper::Interval { 10, 13, 8 },
        lapper::Interval { 100, 200, 9 },
        lapper::Interval { 110, 120, 10 },
        lapper::Interval { 110, 124, 11 },
        lapper::Interval { 111, 160, 12 },
        lapper::Interval { 150, 200, 13 },
    ];
    var storage: i64[>..] = [];
    let index = lapper::new_lapper(intervals, storage);

    var found: i64[>..] = [];
    lapper::find(index, 8, 11) { iv => found.push(iv.val); };
    assert(found == [1, 7, 8]);

    found.clear();
    lapper::find(index, 145, 151) { iv => found.push(iv.val); };
    assert(found == [9, 12, 13]);
}

fn test_long_interval_spanning_short_intervals() {
    let intervals: lapper::Interval<i64, i64>[>..] = [
        lapper::Interval { 25264912, 25264986, 0 },
        lapper::Interval { 27273024, 27273065, 1 },
        lapper::Interval { 27440273, 27440318, 2 },
        lapper::Interval { 27488033, 27488125, 3 },
        lapper::Interval { 27938410, 27938470, 4 },
        lapper::Interval { 27959118, 27959171, 5 },
        lapper::Interval { 28866309, 33141404, 6 },
    ];
    var storage: i64[>..] = [];
    let index = lapper::new_lapper(intervals, storage);

    var found: i64[>..] = [];
    lapper::find(index, 28974798, 33141355) { iv => found.push(iv.val); };
    assert(found == [6]);
}

fn test_single_interval() {
    let intervals: lapper::Interval<i64, i64>[>..] = [
        lapper::Interval { 10, 35, 7 },
    ];
    var storage: i64[>..] = [];
    let index = lapper::new_lapper(intervals, storage);

    var found: i64[>..] = [];
    lapper::find(index, 10, 35) { iv => found.push(iv.val); };
    assert(found == [7]);

    found.clear();
    lapper::find(index, 0, 10) { iv => found.push(iv.val); };
    assert(found.len == 0);

    found.clear();
    lapper::find(index, 35, 40) { iv => found.push(iv.val); };
    assert(found.len == 0);
}

fn next_random(state: u64&) -> u64 {
    var x = state;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    state = x;
    return x;
}

fn below(state: u64&, n: i64) -> i64 { (next_random(state) % (n as u64)) as i64 }

// Build count intervals of one shape, shifted by offset, and compare `find`
// with a scan of every interval for many queries.
fn check_against_scan<T>(shape: i64, count: i64, offset: i64, state: u64&) {
    let span = 20 * count + 100;
    var intervals: lapper::Interval<T, i32>[>..] = [];
    for i in 0..count {
        var start = below(state, span);
        var stop = start + 1 + below(state, 50);
        if shape == 1 {
            // Now and then an interval over nearly everything.
            if below(state, 100) == 0 { start = below(state, 10); stop = span - below(state, 10); }
        } else if shape == 2 {
            // Nested: every interval holds all later ones.
            start = i;
            stop = 2 * count - i;
        } else if shape == 3 {
            // One start, many stops.
            start = 100;
            stop = 101 + below(state, span);
        } else if shape == 4 {
            // Copies of one interval.
            start = 100;
            stop = 200;
        } else if shape == 5 {
            // Any length up to everything.
            stop = start + below(state, span);
        } else if shape == 6 {
            // Already sorted by start, with ties and empty intervals.
            start = 3 * (i / 2);
            stop = start + below(state, 40);
        } else if shape == 7 {
            // Clusters separated by wide gaps.
            start = (i / 40) * 5000 + below(state, 60);
            stop = start + 1 + below(state, 30);
        }
        intervals.push(lapper::Interval { start: (start + offset) as T, stop: (stop + offset) as T, val: i as i32 });
    }

    var storage: T[>..] = [];
    storage.push(default<T>()); // Something ahead of this Lapper's index.
    let index = lapper::new_lapper(intervals, storage);
    for i in 1..intervals.len {
        let a, b = intervals[i - 1], intervals[i];
        assert(a.start < b.start || (a.start == b.start && a.stop <= b.stop));
    }

    var found: i32[>..] = [];
    var expected: i32[>..] = [];
    for q in 0..60 {
        var start = below(state, span + 40) - 20;
        var stop = start + below(state, 60);
        let kind = q % 5;
        if kind == 1 { stop = start + below(state, span); }
        if kind == 2 { start = -20; stop = span + 20; }
        if kind == 3 { stop = start; }
        if offset == 0 && start < 0 { stop -= start; start = 0; }
        let s = (start + offset) as T;
        let e = (stop + offset) as T;
        found.clear();
        expected.clear();
        index.find(s, e) { iv => found.push(iv.val); };
        for &iv in intervals {
            if iv.start < e && iv.stop > s { expected.push(iv.val); }
        }
        assert(found == expected);
    }
}

// Interval counts on both sides of a block of 32 and of the powers of two
// that the tree of blocks grows at.
fn test_find_matches_scan() {
    var state: u64 = 0x9E3779B97F4A7C15;
    let counts: i64[] = [1, 2, 31, 32, 33, 63, 64, 65, 96, 97, 128, 129, 1000, 1024, 1025, 2049, 4100];
    for shape in 0..8 {
        for count in counts {
            check_against_scan<i32>(shape, count, -777, state);
            check_against_scan<u32>(shape, count, 0, state);
            check_against_scan<i64>(shape, count, -5000000000, state);
        }
    }
}

fn main() {
    test_overlap();
    test_new_lapper_sorts_by_start_then_stop();
    test_empty_lapper();
    test_lappers_share_index_storage();
    test_nonoverlapping_boundaries();
    test_containment_queries();
    test_overlapping_intervals();
    test_find_behind_first_match();
    test_find_overlaps_in_large_intervals();
    test_long_interval_spanning_short_intervals();
    test_single_interval();
    test_find_matches_scan();
    print("lapper tests PASSED");
}
