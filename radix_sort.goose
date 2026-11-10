// In-place insertion sort and most-significant-digit radix sort.
//
// `radix_sort` sorts signed or unsigned integer keys eight bits at a time. It
// uses insertion sort for inputs and buckets of at most 64 elements, requires
// constant auxiliary space per recursion level, and is not stable. A key
// function moves whole values while ordering them by the returned key.
//
// The radix algorithm is adapted from Heng Li's cgranges implementation:
// https://github.com/lh3/biofast/blob/4177ef282f4e1addc5b20beb7ee2f6435806f844/lib/cgranges.c#L23
namespace radix_sort;

const RS_MIN_SIZE = 64;
const RS_MAX_BITS = 8;

checkfree fn identity<V>(x: V) -> V { x }

// Map each integer representation to u64 while preserving its numeric order.
// Flipping a signed value's sign bit places negatives before nonnegatives.
checkfree fn ordered_key(key: u8) -> u64 { key }
checkfree fn ordered_key(key: u16) -> u64 { key }
checkfree fn ordered_key(key: u32) -> u64 { key }
checkfree fn ordered_key(key: u64) -> u64 { key }
checkfree fn ordered_key(key: i8) -> u64 { (key as! u8) ^ 0x80 }
checkfree fn ordered_key(key: i16) -> u64 { (key as! u16) ^ 0x8000 }
checkfree fn ordered_key(key: i32) -> u64 { (key as! u32) ^ 0x80000000 }
checkfree fn ordered_key(key: i64) -> u64 { (key as! u64) ^ ((1 as u64) << 63) }

// Return the most-significant eight-bit digit's offset for each key width.
checkfree fn top_shift(key: u8) -> u64 { 0 }
checkfree fn top_shift(key: i8) -> u64 { 0 }
checkfree fn top_shift(key: u16) -> u64 { 8 }
checkfree fn top_shift(key: i16) -> u64 { 8 }
checkfree fn top_shift(key: u32) -> u64 { 24 }
checkfree fn top_shift(key: i32) -> u64 { 24 }
checkfree fn top_shift(key: u64) -> u64 { 56 }
checkfree fn top_shift(key: i64) -> u64 { 56 }

// Sort slice in place with insertion sort, comparing its values directly.
// This sort is stable and works well for small or nearly sorted inputs.
//
// ```goose
// let values = [4, 1, 3, 1];
// radix_sort::insert_sort(values);
// assert(values == [1, 1, 3, 4]);
// ```
checkfree fn insert_sort<V>(slice: V[:]) {
    insert_sort(slice, identity);
}

// Sort slice in place with insertion sort, comparing values through Key.
// Key must be callable from a `checkfree` function. Equal-key values retain
// their original order.
//
// Define key functions at top level so they can be `checkfree`:
//
// ```goose ignore
// checkfree fn low_byte(value: u16) -> u16 { value & 0xFF }
// ```
//
// ```goose
// let values: u16[] = [0x102, 0x001, 0x100];
// radix_sort::insert_sort(values, low_byte);
// let expected: u16[] = [0x100, 0x001, 0x102];
// assert(values == expected);
// ```
checkfree fn insert_sort<V, Key: fn(V) -> T>(slice: V[:]) {
    for i in 1..slice.len {
        if Key(slice[i]) < Key(slice[i - 1]) {
            let temp = slice[i];
            var j = i;
            while j > 0 && Key(temp) < Key(slice[j - 1]) {
                slice[j] = slice[j - 1];
                j -= 1;
            }
            slice[j] = temp;
        }
    }
}

// Sort signed or unsigned integers in place with an MSD radix sort.
//
// The sort is not stable. Inputs of at most 64 elements use insertion sort;
// larger inputs process the most-significant byte first and recursively refine
// each bucket. Empty slices are accepted.
//
// ```goose
// let values: u16[] = [500, 1, 256, 0, 255];
// radix_sort::radix_sort(values);
// let expected: u16[] = [0, 1, 255, 256, 500];
// assert(values == expected);
// ```
fn radix_sort<V>(slice: V[:]) {
    radix_sort(slice, identity);
}

// Sort values in place by the integer returned by Key.
//
// Key must return an `i8`–`i64` or `u8`–`u64` and must be callable from a
// `checkfree` function. Signed keys are normalized so negative values sort
// before nonnegative values. The sort is not stable, so values with equal keys
// may be reordered. The complete value is moved during partitioning; payload
// data remains attached to its key.
//
// Define key functions at top level so they can be `checkfree`:
//
// ```goose ignore
// checkfree fn low_byte(value: u16) -> u16 { value & 0xFF }
// ```
//
// ```goose
// let values: u16[] = [0x102, 0x001, 0x100];
// radix_sort::radix_sort(values, low_byte);
// let expected: u16[] = [0x100, 0x001, 0x102];
// assert(values == expected);
// ```
fn radix_sort<V, Key: fn(V) -> T>(slice: V[:]) {
    if slice.len <= RS_MIN_SIZE {
        insert_sort(slice, Key);
        return;
    }

    radix_sort_inner(slice, top_shift(Key(slice[0])), Key);
}

recursive fn radix_sort_inner<V, Key: fn(V) -> T>(slice: V[:], s: u64) {
    let size: i64 = 1 << RS_MAX_BITS;
    let mask: u64 = (1 << RS_MAX_BITS) - 1;
    // Bucket k occupies slice[begins[k]..ends[k]]. Two tables of 8-byte
    // indices instead of one of begin/end pairs let a lookup scale its index
    // in the address, with no shift on the cycle loop's dependency chain.
    let begins: i64[1 << RS_MAX_BITS];
    let ends:   i64[1 << RS_MAX_BITS];

    // Phase 1: count digits into ends. Zero in place to avoid a temporary array.
    for i in 0..size { ends[i] = 0; }
    for i in 0..slice.len {
        ends[(ordered_key(Key(slice[i])) >> s) & mask] += 1;
    }

    // Phase 2: turn the counts into each bucket's indices in slice
    var previous_end: i64 = 0;
    for k in 0..size {
        begins[k] = previous_end;
        previous_end += ends[k];
        ends[k] = previous_end;
    }

    // Phase 3: Fill the buckets, using begins as each bucket's cursor.
    // This partitions the slice by the current digit. Keep the bucket loop
    // separate from the loop that places records: one loop that does either
    // can compile to a bucket number that depends on the table loads.
    for b in 0..size {
        while begins[b] != ends[b] {
            let origin = begins[b];
            var tmp = slice[origin];
            var decoded = ((ordered_key(Key(tmp)) >> s) & mask) as i64;
            while decoded != b { // record belongs elsewhere
                let destination = begins[decoded];
                let displaced = slice[destination];
                slice[destination] = tmp;
                begins[decoded] = destination + 1;
                tmp = displaced;
                decoded = ((ordered_key(Key(tmp)) >> s) & mask) as i64;
            }
            slice[origin] = tmp;
            begins[b] = origin + 1;
        }
    }

    // Phase 4: refine each bucket, whose records now share the current digit.
    // Phase 3 advanced every begin to its end, but a bucket begins where the
    // previous one ends, so the ranges need no restoring.
    if s != 0 {
        let new_s = if s > RS_MAX_BITS { s - RS_MAX_BITS } else { 0 };
        var bucket_begin: i64 = 0;
        for k in 0..size { // Sort within each bucket
            let bucket_end = ends[k];
            if bucket_end - bucket_begin > RS_MIN_SIZE {
                radix_sort_inner(slice[bucket_begin..bucket_end], new_s, Key);
            } else if bucket_end - bucket_begin > 1 {
                insert_sort(slice[bucket_begin..bucket_end], Key);
            }
            bucket_begin = bucket_end;
        }
    }
}
