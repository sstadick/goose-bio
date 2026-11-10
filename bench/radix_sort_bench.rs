//! Two Rust ports of the cgranges radix sort, timed against the C reference.
//!
//! `safe` keeps integer cursors into a bounds-checked slice, like the Goose
//! and Mojo ports. `unsafe` keeps raw record addresses, like the Zig port.

use std::env;
use std::io::{self, BufWriter, Write};
use std::process;
use std::time::Instant;

const MIN_SIZE: usize = 64;
const RADIX_BITS: u32 = 8;
const BUCKET_COUNT: usize = 1 << RADIX_BITS;

#[derive(Clone, Copy)]
#[repr(C)]
struct Record {
    key: u64,
    payload: u64,
}

unsafe extern "C" {
    fn c_setup(len: i64);
    fn c_reset();
    fn c_radix_sort();
    fn c_validate() -> bool;
    fn c_checksum() -> u64;
    fn c_teardown();
}

#[inline(always)]
fn digit(record: Record, shift: u32) -> usize {
    ((record.key >> shift) & (BUCKET_COUNT as u64 - 1)) as usize
}

mod safe {
    use super::{BUCKET_COUNT, MIN_SIZE, RADIX_BITS, Record, digit};
    use std::mem;

    fn insertion_sort(slice: &mut [Record]) {
        for i in 1..slice.len() {
            if slice[i].key < slice[i - 1].key {
                let temporary = slice[i];
                let mut j = i;
                while j > 0 && temporary.key < slice[j - 1].key {
                    slice[j] = slice[j - 1];
                    j -= 1;
                }
                slice[j] = temporary;
            }
        }
    }

    /// Unlike `rs_sort_record` in `radix_sort_reference.c`, this keeps the
    /// bucket cursors and ends in two tables of 8-byte indices, so a table
    /// lookup needs no shift, which pays for the shift that turns an index
    /// into an address. A bucket begins where the previous one ends, so
    /// nothing restores the `begins` that the permutation advanced.
    fn radix_sort_inner(slice: &mut [Record], shift: u32) {
        let mut begins = [0usize; BUCKET_COUNT];
        let mut ends = [0usize; BUCKET_COUNT];
        for record in slice.iter() {
            ends[digit(*record, shift)] += 1;
        }
        let mut previous_end = 0;
        for bucket in 0..BUCKET_COUNT {
            begins[bucket] = previous_end;
            previous_end += ends[bucket];
            ends[bucket] = previous_end;
        }

        for bucket in 0..BUCKET_COUNT {
            while begins[bucket] != ends[bucket] {
                let origin = begins[bucket];
                let mut temporary = slice[origin];
                let mut target = digit(temporary, shift);
                while target != bucket {
                    let destination = &mut begins[target];
                    temporary = mem::replace(&mut slice[*destination], temporary);
                    *destination += 1;
                    target = digit(temporary, shift);
                }
                slice[origin] = temporary;
                begins[bucket] = origin + 1;
            }
        }

        if shift != 0 {
            let next_shift = shift.saturating_sub(RADIX_BITS);
            let mut begin = 0;
            for end in ends {
                let part = &mut slice[begin..end];
                if part.len() > MIN_SIZE {
                    radix_sort_inner(part, next_shift);
                } else if part.len() > 1 {
                    insertion_sort(part);
                }
                begin = end;
            }
        }
    }

    pub fn radix_sort(slice: &mut [Record]) {
        if slice.len() <= MIN_SIZE {
            insertion_sort(slice);
        } else {
            radix_sort_inner(slice, u64::BITS - RADIX_BITS);
        }
    }
}

mod raw {
    use super::{BUCKET_COUNT, MIN_SIZE, RADIX_BITS, Record, digit};

    /// A record count while counting, then the addresses of a bucket's next
    /// unplaced record and of its end.
    #[derive(Clone, Copy)]
    struct Bucket {
        begin: usize,
        end: usize,
    }

    const RECORD_SIZE: usize = size_of::<Record>();

    /// # Safety
    ///
    /// `begin..end` must be one live, exclusively owned range of at least two
    /// records.
    unsafe fn insertion_sort(begin: *mut Record, end: *mut Record) {
        // SAFETY: `current` stays in `begin + 1..end` and `destination` in
        // `begin..=current`, so every access is inside the caller's range.
        unsafe {
            let mut current = begin.add(1);
            while current < end {
                if (*current).key < (*current.sub(1)).key {
                    let temporary = *current;
                    let mut destination = current;
                    while destination > begin && temporary.key < (*destination.sub(1)).key {
                        *destination = *destination.sub(1);
                        destination = destination.sub(1);
                    }
                    *destination = temporary;
                }
                current = current.add(1);
            }
        }
    }

    /// Every recursive call keeps its own 4 KiB bucket table on the stack, as
    /// `rs_sort_record` in `radix_sort_reference.c` does. Unlike the C, the
    /// table starts as zeroed counts, so one `memset` initializes it, and the
    /// prefix sum turns the counts into record addresses for the permutation.
    ///
    /// # Safety
    ///
    /// `begin..end` must be one live, exclusively owned range of records.
    unsafe fn radix_sort_inner(begin: *mut Record, end: *mut Record, shift: u32) {
        let mut buckets = [Bucket { begin: 0, end: 0 }; BUCKET_COUNT];
        // SAFETY: the prefix sum makes the bucket ranges partition
        // `begin..end`. Each cursor then only advances within its own bucket,
        // so every address is inside `begin..end` and takes its provenance
        // from `begin`.
        unsafe {
            let mut item = begin;
            while item != end {
                buckets[digit(*item, shift)].end += 1;
                item = item.add(1);
            }
            let mut previous_end = begin;
            for bucket in &mut buckets {
                bucket.begin = previous_end.addr();
                previous_end = previous_end.add(bucket.end);
                bucket.end = previous_end.addr();
            }

            let mut bucket = 0;
            while bucket < BUCKET_COUNT {
                let cursor = buckets[bucket].begin;
                if cursor != buckets[bucket].end {
                    let origin = begin.with_addr(cursor);
                    let mut temporary = *origin;
                    let mut target = digit(temporary, shift);
                    if target != bucket {
                        loop {
                            let destination = &mut buckets[target].begin;
                            temporary = begin.with_addr(*destination).replace(temporary);
                            *destination += RECORD_SIZE;
                            target = digit(temporary, shift);
                            if target == bucket {
                                break;
                            }
                        }
                        *origin = temporary;
                    }
                    buckets[bucket].begin = cursor + RECORD_SIZE;
                } else {
                    bucket += 1;
                }
            }

            buckets[0].begin = begin.addr();
            for k in 1..BUCKET_COUNT {
                buckets[k].begin = buckets[k - 1].end;
            }

            if shift != 0 {
                let next_shift = shift.saturating_sub(RADIX_BITS);
                for bucket in &buckets {
                    let length = (bucket.end - bucket.begin) / RECORD_SIZE;
                    let part_begin = begin.with_addr(bucket.begin);
                    let part_end = begin.with_addr(bucket.end);
                    if length > MIN_SIZE {
                        radix_sort_inner(part_begin, part_end, next_shift);
                    } else if length > 1 {
                        insertion_sort(part_begin, part_end);
                    }
                }
            }
        }
    }

    pub fn radix_sort(slice: &mut [Record]) {
        if slice.len() < 2 {
            return;
        }
        let range = slice.as_mut_ptr_range();
        // SAFETY: the range is the exclusively borrowed slice, which holds at
        // least two records.
        unsafe {
            if slice.len() <= MIN_SIZE {
                insertion_sort(range.start, range.end);
            } else {
                radix_sort_inner(range.start, range.end, u64::BITS - RADIX_BITS);
            }
        }
    }
}

fn next_random(state: &mut u64) -> u64 {
    let mut x = *state;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    *state = x;
    x
}

fn checksum(work: &[Record]) -> u64 {
    work[0].key ^ work[work.len() / 2].payload ^ work[work.len() - 1].key
}

fn run_case(
    label: &str,
    sort: fn(&mut [Record]),
    n: usize,
    repetitions: usize,
    out: &mut impl Write,
) -> io::Result<()> {
    const PAYLOAD_MASK: u64 = 0x9E3779B97F4A7C15;
    let mut state = 0x4D595DF4D0F33173;
    let source: Vec<Record> = (0..n)
        .map(|_| {
            let key = next_random(&mut state);
            Record {
                key,
                payload: key ^ PAYLOAD_MASK,
            }
        })
        .collect();
    let mut work = source.clone();

    // Warm and validate both implementations before collecting samples.
    sort(&mut work);
    assert!(
        work.iter()
            .all(|record| record.payload == record.key ^ PAYLOAD_MASK)
    );
    assert!(work.is_sorted_by_key(|record| record.key));
    // SAFETY: the reference only touches the buffers that `c_setup` allocates
    // for `n` records, and nothing else calls it.
    unsafe {
        c_setup(n as i64);
        c_reset();
        c_radix_sort();
        assert!(c_validate());
    }

    let mut rust_checksum: u64 = 0;
    let mut reference_checksum: u64 = 0;
    for repetition in 0..repetitions {
        // Alternate order to balance short-term frequency and temperature changes.
        for turn in 0..2 {
            if (repetition + turn) % 2 == 0 {
                work.copy_from_slice(&source);
                let start = Instant::now();
                sort(&mut work);
                let elapsed = start.elapsed().as_nanos();
                rust_checksum = rust_checksum.wrapping_add(checksum(&work));
                writeln!(out, "{label},{n},{elapsed}")?;
            } else {
                // SAFETY: as above; `c_setup` has run for this case.
                let elapsed = unsafe {
                    c_reset();
                    let start = Instant::now();
                    c_radix_sort();
                    let elapsed = start.elapsed().as_nanos();
                    reference_checksum = reference_checksum.wrapping_add(c_checksum());
                    elapsed
                };
                writeln!(out, "c,{n},{elapsed}")?;
            }
        }
    }
    assert_eq!(rust_checksum, reference_checksum);
    // SAFETY: as above; no reference buffer is used after this.
    unsafe { c_teardown() };
    Ok(())
}

fn main() -> io::Result<()> {
    let (label, sort): (&str, fn(&mut [Record])) = match env::args().nth(1).as_deref() {
        Some("safe") => ("rust-safe", safe::radix_sort),
        Some("unsafe") => ("rust-unsafe", raw::radix_sort),
        _ => {
            eprintln!("usage: radix-sort-bench safe|unsafe");
            process::exit(2);
        }
    };

    let mut out = BufWriter::new(io::stdout().lock());
    for (n, repetitions) in [(1024, 1000), (16384, 200), (262144, 20), (1048576, 6)] {
        run_case(label, sort, n, repetitions, &mut out)?;
    }
    out.flush()
}
