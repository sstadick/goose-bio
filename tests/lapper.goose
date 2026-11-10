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

fn test_new_lapper_sorts_and_tracks_longest() {
    let intervals: lapper::Interval<i64, i64>[>..] = [
        lapper::Interval { 5, 7, 3 },
        lapper::Interval { 0, 5, 0 },
        lapper::Interval { 4, 6, 2 },
        lapper::Interval { 0, 3, 1 },
    ];
    let index = lapper::new_lapper(intervals);

    assert(index.longest == 5);
    assert(index.intervals.len == 4);
    assert(index.intervals[0].start == 0 && index.intervals[0].stop == 3);
    assert(index.intervals[1].start == 0 && index.intervals[1].stop == 5);
    assert(index.intervals[2].start == 4 && index.intervals[2].stop == 6);
    assert(index.intervals[3].start == 5 && index.intervals[3].stop == 7);
}

fn test_lower_bound() {
    let intervals: lapper::Interval<i64, i64>[>..] = [
        lapper::Interval { 0, 3, 0 },
        lapper::Interval { 2, 3, 1 },
        lapper::Interval { 4, 6, 2 },
        lapper::Interval { 5, 7, 3 },
        lapper::Interval { 9, 11, 4 },
    ];

    assert(lapper::lower_bound(-5, intervals) == 0);
    assert(lapper::lower_bound(0, intervals) == 0);
    assert(lapper::lower_bound(1, intervals) == 1);
    assert(lapper::lower_bound(2, intervals) == 1);
    assert(lapper::lower_bound(4, intervals) == 2);
    assert(lapper::lower_bound(6, intervals) == 4);
    assert(lapper::lower_bound(9, intervals) == 4);
    assert(lapper::lower_bound(10, intervals) == 5);

    let empty: lapper::Interval<i64, i64>[>..] = [];
    assert(lapper::lower_bound(0, empty) == -1);
}

fn test_nonoverlapping_boundaries() {
    let intervals: lapper::Interval<i64, i64>[>..] = [
        lapper::Interval { 0, 10, 0 },
        lapper::Interval { 20, 30, 1 },
        lapper::Interval { 40, 50, 2 },
        lapper::Interval { 60, 70, 3 },
        lapper::Interval { 80, 90, 4 },
    ];
    let index = lapper::new_lapper(intervals);

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
    let index = lapper::new_lapper(intervals);

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
    let index = lapper::new_lapper(intervals);

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
    let index = lapper::new_lapper(intervals);
    assert(index.longest == 50);

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
    let index = lapper::new_lapper(intervals);

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
    let index = lapper::new_lapper(intervals);

    var found: i64[>..] = [];
    lapper::find(index, 28974798, 33141355) { iv => found.push(iv.val); };
    assert(found == [6]);
}

fn test_single_interval() {
    let intervals: lapper::Interval<i64, i64>[>..] = [
        lapper::Interval { 10, 35, 7 },
    ];
    let index = lapper::new_lapper(intervals);

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

fn main() {
    test_overlap();
    test_new_lapper_sorts_and_tracks_longest();
    test_lower_bound();
    test_nonoverlapping_boundaries();
    test_containment_queries();
    test_overlapping_intervals();
    test_find_behind_first_match();
    test_find_overlaps_in_large_intervals();
    test_long_interval_spanning_short_intervals();
    test_single_interval();
    print("lapper tests PASSED");
}
