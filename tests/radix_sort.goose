import bio from "..";
import bio.radix_sort;

struct Item {
    key:     u16,
    payload: i64,
}

struct SignedItem {
    key:     i32,
    payload: i64,
}

checkfree fn item_key(item: Item) -> u16 { item.key }
checkfree fn signed_item_key(item: SignedItem) -> i32 { item.key }

fn test_insert_sort() {
    let empty: i64[] = [];
    radix_sort::insert_sort(empty);
    assert(empty.len == 0);

    let values = [5, 7, 4, 2, 3, 0, 11, 15, 13, 3];
    radix_sort::insert_sort(values);
    assert(values == [0, 2, 3, 3, 4, 5, 7, 11, 13, 15]);

    let items = [
        Item { key: 2, payload: 20 },
        Item { key: 1, payload: 10 },
        Item { key: 2, payload: 21 },
    ];
    radix_sort::insert_sort(items, item_key);
    assert(items[0].payload == 10);
    assert(items[1].payload == 20);
    assert(items[2].payload == 21); // Insertion sort is stable.
}

fn test_radix_sort_small_inputs() {
    let empty: u16[] = [];
    radix_sort::radix_sort(empty);
    assert(empty.len == 0);

    let singleton: u16[] = [42];
    radix_sort::radix_sort(singleton);
    assert(singleton.len == 1 && singleton[0] == 42);

    // The cutoff itself takes the insertion-sort path.
    var values: u16[>..] = [];
    for i in 0..64 { values.push((63 - i) as u16); }
    radix_sort::radix_sort(values);
    for i in 0..values.len { assert(values[i] == i as u16); }
}

fn test_radix_sort_multiple_digits() {
    // More than 64 values exercises radix partitioning and recursive refinement.
    var values: u16[>..] = [];
    for i in 0..300 { values.push((299 - i) as u16); }

    radix_sort::radix_sort(values);

    for i in 0..values.len { assert(values[i] == i as u16); }
}

fn test_radix_sort_signed_values() {
    var values: i16[>..] = [];
    for i in 0..300 { values.push((149 - i) as i16); }

    radix_sort::radix_sort(values);

    for i in 0..values.len { assert(values[i] == (i - 150) as i16); }
}

fn test_radix_sort_integer_widths_and_extrema() {
    var u8s: u8[>..] = [];
    for i in 0..200 { u8s.push((199 - i) as u8); }
    u8s.push(255);
    radix_sort::radix_sort(u8s);
    for i in 0..200 { assert(u8s[i] == i as u8); }
    assert(u8s[200] == 255);

    var i8s: i8[>..] = [];
    for i in 0..200 { i8s.push((99 - i) as i8); }
    i8s.push(-128);
    i8s.push(127);
    radix_sort::radix_sort(i8s);
    assert(i8s[0] == -128 && i8s[201] == 127);
    for i in 1..201 { assert(i8s[i] == (i - 101) as i8); }

    var u32s: u32[>..] = [];
    for i in 0..130 { u32s.push((129 - i) as u32); }
    u32s.push(0xFFFFFFFF);
    radix_sort::radix_sort(u32s);
    for i in 0..130 { assert(u32s[i] == i as u32); }
    assert(u32s[130] == 0xFFFFFFFF);

    var i32s: i32[>..] = [];
    for i in 0..130 { i32s.push((64 - i) as i32); }
    i32s.push(-2147483647 - 1);
    i32s.push(2147483647);
    radix_sort::radix_sort(i32s);
    assert(i32s[0] == -2147483647 - 1 && i32s[131] == 2147483647);
    for i in 1..131 { assert(i32s[i] == (i - 66) as i32); }

    var i64s: i64[>..] = [];
    for i in 0..130 { i64s.push(64 - i); }
    i64s.push(-0x7FFFFFFFFFFFFFFF - 1);
    i64s.push(0x7FFFFFFFFFFFFFFF);
    radix_sort::radix_sort(i64s);
    assert(i64s[0] == -0x7FFFFFFFFFFFFFFF - 1);
    assert(i64s[131] == 0x7FFFFFFFFFFFFFFF);
    for i in 1..131 { assert(i64s[i] == i - 66); }
}

fn test_radix_sort_keyed_payloads() {
    var items: Item[>..] = [];
    for i in 0..300 {
        let key = ((i * 73) % 300) as u16;
        items.push(Item { key: key, payload: key as i64 * 11 + 7 });
    }

    radix_sort::radix_sort(items, item_key);

    for i in 0..items.len {
        assert(items[i].key == i as u16);
        assert(items[i].payload == i * 11 + 7);
    }
}

fn test_radix_sort_signed_keys() {
    var items: SignedItem[>..] = [];
    for i in 0..300 {
        let key = ((i * 73) % 300 - 150) as i32;
        items.push(SignedItem { key: key, payload: key as i64 * 13 - 5 });
    }

    radix_sort::radix_sort(items, signed_item_key);

    for i in 0..items.len {
        let key = i - 150;
        assert(items[i].key == key as i32);
        assert(items[i].payload == key * 13 - 5);
    }
}

fn test_radix_sort_duplicates() {
    var items: Item[>..] = [];
    for i in 0..300 {
        items.push(Item { key: ((i * 37) % 17) as u16, payload: i });
    }

    radix_sort::radix_sort(items, item_key);

    let seen:   bool[300] = [false; 300];
    let counts: i64[17]   = [0; 17];
    for i in 0..items.len {
        if i > 0 { assert(items[i - 1].key <= items[i].key); }
        assert(!seen[items[i].payload]);
        seen[items[i].payload] = true;
        counts[items[i].key] += 1;
    }
    for present in seen { assert(present); }
    var total = 0;
    for count in counts { total += count; }
    assert(total == 300);
}

fn test_radix_sort_high_u64_digit() {
    var values:       u64[>..] = [];
    var expected_xor: u64      = 0;
    for i in 0..130 {
        let value = if i % 2 == 0 { (1 as u64) << 63 | i as u64 } else { i as u64 };
        values.push(value);
        expected_xor ^= value;
    }

    radix_sort::radix_sort(values);

    var actual_xor: u64 = 0;
    for i in 0..values.len {
        if i > 0 { assert(values[i - 1] <= values[i]); }
        actual_xor ^= values[i];
    }
    assert(actual_xor == expected_xor);
}

fn main() {
    test_insert_sort();
    test_radix_sort_small_inputs();
    test_radix_sort_multiple_digits();
    test_radix_sort_signed_values();
    test_radix_sort_integer_widths_and_extrema();
    test_radix_sort_keyed_payloads();
    test_radix_sort_signed_keys();
    test_radix_sort_duplicates();
    test_radix_sort_high_u64_digit();
}
