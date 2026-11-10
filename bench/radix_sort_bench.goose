import os;
import bio from "..";
import bio.radix_sort;

struct Record {
    key:     u64,
    payload: u64,
}

extern fn c_setup(len: i64);
extern fn c_reset();
extern fn c_radix_sort();
extern fn c_validate() -> bool;
extern fn c_checksum() -> u64;
extern fn c_teardown();

checkfree fn record_key(record: Record) -> u64 { record.key }

fn next_random(state: u64&) -> u64 {
    var x = state;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    state = x;
    return x;
}

fn reset(source: Record[:], work: Record[:]) {
    for i in 0..source.len { work[i] = source[i]; }
}

fn goose_sample(work: Record[:]) -> i64 {
    let start = clock_ns();
    radix_sort::radix_sort(work, record_key);
    return clock_ns() - start;
}

fn c_sample() -> i64 {
    let start = clock_ns();
    c_radix_sort();
    return clock_ns() - start;
}

fn run_case(n: i64, repetitions: i64) {
    var source: Record[>..] = [];
    var work:   Record[>..] = [];
    var state:  u64         = 0x4D595DF4D0F33173;
    for _ in 0..n {
        let key = next_random(state);
        source.push(Record { key: key, payload: key ^ 0x9E3779B97F4A7C15 });
        work.push(Record { .. });
    }

    // Warm and validate both implementations before collecting samples.
    reset(source, work);
    radix_sort::radix_sort(work, record_key);
    for i in 1..n { assert(work[i - 1].key <= work[i].key); }
    c_setup(n);
    c_reset();
    c_radix_sort();
    assert(c_validate());

    var goose_checksum: u64 = 0;
    var c_checksum:     u64 = 0;
    for repetition in 0..repetitions {
        // Alternate order to balance short-term frequency and temperature changes.
        if repetition % 2 == 0 {
            reset(source, work);
            let elapsed = goose_sample(work);
            goose_checksum += work[0].key ^ work[n / 2].payload ^ work[n - 1].key;
            print("goose,", n, ",", elapsed);

            c_reset();
            let c_elapsed = c_sample();
            c_checksum += ::c_checksum();
            print("c,", n, ",", c_elapsed);
        } else {
            c_reset();
            let c_elapsed = c_sample();
            c_checksum += ::c_checksum();
            print("c,", n, ",", c_elapsed);

            reset(source, work);
            let elapsed = goose_sample(work);
            goose_checksum += work[0].key ^ work[n / 2].payload ^ work[n - 1].key;
            print("goose,", n, ",", elapsed);
        }
    }
    assert(goose_checksum == c_checksum);
    c_teardown();
    write_stderr(str("checksum ", n, ": ", goose_checksum, "\n"));
}

fn main() {
    run_case(1024, 1000);
    run_case(16384, 200);
    run_case(262144, 20);
    run_case(1048576, 6);
}
