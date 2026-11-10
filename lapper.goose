namespace lapper;

struct Interval<T, V> {
    // 0-based inclusive start
    start: T,
    // 0-based exclusive end
    stop: T,
    // A stored value.
    val: V,
}

inline fn overlap<T, V>(a: Interval<T, V>&, b: Interval<T, V>&) -> bool {
    a.start < b.stop && a.stop > b.start
}
inline fn overlap<T, V>(a: Interval<T, V>&, b_start: T, b_stop: T) -> bool {
    a.start < b_stop && a.stop > b_start
}

inline fn less<T, V>(a: Interval<T, V>&, b: Interval<T, V>&) -> bool {
    if a.start < b.start {
        true
    } else if a.start == b.start && a.stop < b.stop {
        true
    } else {
        false
    }
}

struct Lapper<T, V> {
    longest:   i64,
    intervals: Interval<T, V>[:],
}

fn new_lapper<T, V>(intervals: Interval<T, V>[:]) -> Lapper<T, V> {
    sort(intervals, less);
    var longest = 0;
    for iv in intervals {
        let len = iv.stop - iv.start;
        longest = if len > longest { len } else { longest };
    }
    Lapper { longest, intervals }
}

fn find<T, V, F: fn(Interval<T, V>&)>(lapper: Lapper<T, V>&, start: T, stop: T) {
    var low = lower_bound((start - lapper.longest) as T, lapper.intervals);
    for iv in lapper.intervals[low..] {
        if iv.overlap(start, stop) {
            F(iv)
        } else if iv.start > stop {
            break
        }
    }
}

// Assumes longest inteval has been subracted form start
fn lower_bound<T, V>(start: T, ivs: Interval<T, V>[:]) -> i64 {
    // find the i that is to the left of the passed in stop
    if ivs.len == 0 {
        return -1;
    }

    var size = ivs.len;
    var low = 0;

    while size > 0 {
        let half = size / 2;
        let other_half = size - half;
        let probe = low + half;
        let other_probe = low + other_half;
        let v .= ivs[probe];
        size = half;
        low = if v.start < start { other_probe } else { low };
    }
    low
}

inline fn bool_as<T>(cond: bool) -> T {
    if cond {
        1
    } else {
        0
    }
}
