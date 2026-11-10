from std.ffi import external_call
from std.memory import unsafe_memcpy
from std.memory.alloc import unsafe_alloc
from std.sys import argv, exit, stderr
from std.testing import assert_equal, assert_true
from std.time import monotonic


comptime MIN_SIZE = 64
comptime RADIX_BITS = 8
comptime BUCKET_COUNT = 1 << RADIX_BITS


@fieldwise_init
struct Record(ImplicitlyCopyable):
    var key: UInt64
    var payload: UInt64


comptime RecordPtr = Pointer[mut=True, Record, MutUntrackedOrigin]


@always_inline("nodebug")
def digit(record: Record, shift: Int) -> Int:
    return Int((record.key >> UInt64(shift)) & UInt64(BUCKET_COUNT - 1))


# The `mojo` port: integer cursors into unchecked pointer offsets.


def insertion_sort(data: RecordPtr, length: Int):
    for current in range(1, length):
        var temporary = data[unsafe_offset=current]
        if temporary.key < data[unsafe_offset=current - 1].key:
            var destination = current
            while destination > 0:
                var previous = data[unsafe_offset=destination - 1]
                if temporary.key >= previous.key:
                    break
                data[unsafe_offset=destination] = previous
                destination -= 1
            data[unsafe_offset=destination] = temporary


# Every recursive call keeps its own 4 KiB of bucket tables on the stack, as
# rs_sort_record in radix_sort_reference.c does. Three things differ from the
# C, and RADIX_SORT.md measures each:
#
# - The tables hold record indices, not addresses. Splitting them into
#   `begins` and `ends` makes each entry 8 bytes, so a table lookup needs no
#   shift, which pays for the shift that turns an index into an address.
# - `ends` starts as zeroed counts, which one memset initializes.
# - A bucket begins where the previous one ends, so nothing restores the
#   `begins` that the permutation advanced.
#
# The permutation is a `for` over buckets with a `while` inside. Written as
# one `while` that either places a record or advances the bucket, Mojo
# compiles the advance without a branch, and every step of the scan over the
# 256 buckets then waits for the previous step's table load.
@no_inline
def radix_sort_inner(data: RecordPtr, length: Int, shift: Int):
    var begins = Array[Int, BUCKET_COUNT](uninitialized=True)
    var ends = Array[Int, BUCKET_COUNT](fill=0)

    for i in range(length):
        ends.unsafe_get(digit(data[unsafe_offset=i], shift)) += 1

    var previous_end = 0
    for bucket in range(BUCKET_COUNT):
        begins.unsafe_get(bucket) = previous_end
        previous_end += ends.unsafe_get(bucket)
        ends.unsafe_get(bucket) = previous_end

    for bucket in range(BUCKET_COUNT):
        while begins.unsafe_get(bucket) != ends.unsafe_get(bucket):
            var cursor = begins.unsafe_get(bucket)
            var temporary = data[unsafe_offset=cursor]
            var target = digit(temporary, shift)
            while target != bucket:
                var destination = begins.unsafe_get(target)
                var displaced = data[unsafe_offset=destination]
                data[unsafe_offset=destination] = temporary
                begins.unsafe_get(target) = destination + 1
                temporary = displaced
                target = digit(temporary, shift)
            data[unsafe_offset=cursor] = temporary
            begins.unsafe_get(bucket) = cursor + 1

    if shift != 0:
        var next_shift = max(0, shift - RADIX_BITS)
        var bucket_begin = 0
        for bucket in range(BUCKET_COUNT):
            var bucket_end = ends.unsafe_get(bucket)
            var bucket_length = bucket_end - bucket_begin
            if bucket_length > MIN_SIZE:
                radix_sort_inner(
                    data.unsafe_offset(bucket_begin), bucket_length, next_shift
                )
            elif bucket_length > 1:
                insertion_sort(data.unsafe_offset(bucket_begin), bucket_length)
            bucket_begin = bucket_end


@always_inline
def radix_sort(data: RecordPtr, length: Int):
    if length <= MIN_SIZE:
        insertion_sort(data, length)
    else:
        radix_sort_inner(data, length, 64 - RADIX_BITS)


# The `mojo-safe` port: the same sort over a bounds-checked Span. The digit
# and the bucket number are provably below 256, so only the record accesses
# keep their checks.


def safe_insertion_sort[origin: MutOrigin](data: Span[Record, origin]):
    for current in range(1, len(data)):
        var temporary = data[current]
        if temporary.key < data[current - 1].key:
            var destination = current
            while destination > 0:
                var previous = data[destination - 1]
                if temporary.key >= previous.key:
                    break
                data[destination] = previous
                destination -= 1
            data[destination] = temporary


@no_inline
def safe_radix_sort_inner[
    origin: MutOrigin
](data: Span[Record, origin], shift: Int):
    var begins = Array[Int, BUCKET_COUNT](fill=0)
    var ends = Array[Int, BUCKET_COUNT](fill=0)

    for record in data:
        ends[digit(record, shift)] += 1

    var previous_end = 0
    for bucket in range(BUCKET_COUNT):
        begins[bucket] = previous_end
        previous_end += ends[bucket]
        ends[bucket] = previous_end

    for bucket in range(BUCKET_COUNT):
        while begins[bucket] != ends[bucket]:
            var cursor = begins[bucket]
            var temporary = data[cursor]
            var target = digit(temporary, shift)
            while target != bucket:
                var destination = begins[target]
                var displaced = data[destination]
                data[destination] = temporary
                begins[target] = destination + 1
                temporary = displaced
                target = digit(temporary, shift)
            data[cursor] = temporary
            begins[bucket] = cursor + 1

    if shift != 0:
        var next_shift = max(0, shift - RADIX_BITS)
        var bucket_begin = 0
        for bucket in range(BUCKET_COUNT):
            var bucket_end = ends[bucket]
            var bucket_length = bucket_end - bucket_begin
            if bucket_length > MIN_SIZE:
                safe_radix_sort_inner(data[bucket_begin:bucket_end], next_shift)
            elif bucket_length > 1:
                safe_insertion_sort(data[bucket_begin:bucket_end])
            bucket_begin = bucket_end


@always_inline
def safe_radix_sort[origin: MutOrigin](data: Span[Record, origin]):
    if len(data) <= MIN_SIZE:
        safe_insertion_sort(data)
    else:
        safe_radix_sort_inner(data, 64 - RADIX_BITS)


@always_inline
def sort[safe: Bool](data: RecordPtr, length: Int):
    comptime if safe:
        safe_radix_sort(
            Span[Record, MutUntrackedOrigin](unsafe_ptr=data, length=length)
        )
    else:
        radix_sort(data, length)


@always_inline("nodebug")
def next_random(mut state: UInt64) -> UInt64:
    var value = state
    value ^= value << 13
    value ^= value >> 7
    value ^= value << 17
    state = value
    return value


def reset(source: RecordPtr, work: RecordPtr, length: Int):
    unsafe_memcpy(dest=work, src=source, count=length)


def run_case[safe: Bool](length: Int, repetitions: Int) raises:
    var source = unsafe_alloc[Record](length)
    var work = unsafe_alloc[Record](length)
    var state = UInt64(0x4D595DF4D0F33173)
    var source_sum = UInt64(0)
    var source_xor = UInt64(0)
    for i in range(length):
        var key = next_random(state)
        var record = Record(key, key ^ UInt64(0x9E3779B97F4A7C15))
        source.unsafe_offset(i).unsafe_write(record)
        work.unsafe_offset(i).unsafe_write(record)
        source_sum += key
        source_xor ^= key

    reset(source, work, length)
    sort[safe](work, length)
    var output_sum = UInt64(0)
    var output_xor = UInt64(0)
    for i in range(length):
        var record = work[unsafe_offset=i]
        assert_equal(record.payload, record.key ^ UInt64(0x9E3779B97F4A7C15))
        output_sum += record.key
        output_xor ^= record.key
        if i != 0:
            assert_true(work[unsafe_offset=i - 1].key <= record.key)
    assert_equal(output_sum, source_sum)
    assert_equal(output_xor, source_xor)

    external_call["c_setup", NoneType](length)
    external_call["c_reset", NoneType]()
    external_call["c_radix_sort", NoneType]()
    assert_true(external_call["c_validate", Bool]())

    var mojo_checksum = UInt64(0)
    var c_checksum = UInt64(0)
    for repetition in range(repetitions):
        if repetition % 2 == 0:
            reset(source, work, length)
            var start = monotonic()
            sort[safe](work, length)
            var elapsed = monotonic() - start
            mojo_checksum += (
                work[].key
                ^ work[unsafe_offset=length // 2].payload
                ^ work[unsafe_offset=length - 1].key
            )
            print("mojo,", length, ",", elapsed, sep="")

            external_call["c_reset", NoneType]()
            start = monotonic()
            external_call["c_radix_sort", NoneType]()
            elapsed = monotonic() - start
            c_checksum += external_call["c_checksum", UInt64]()
            print("c,", length, ",", elapsed, sep="")
        else:
            external_call["c_reset", NoneType]()
            var start = monotonic()
            external_call["c_radix_sort", NoneType]()
            var elapsed = monotonic() - start
            c_checksum += external_call["c_checksum", UInt64]()
            print("c,", length, ",", elapsed, sep="")

            reset(source, work, length)
            start = monotonic()
            sort[safe](work, length)
            elapsed = monotonic() - start
            mojo_checksum += (
                work[].key
                ^ work[unsafe_offset=length // 2].payload
                ^ work[unsafe_offset=length - 1].key
            )
            print("mojo,", length, ",", elapsed, sep="")

    assert_equal(mojo_checksum, c_checksum)
    external_call["c_teardown", NoneType]()
    source.unsafe_free()
    work.unsafe_free()


def run_cases[safe: Bool]() raises:
    run_case[safe](1024, 1000)
    run_case[safe](16384, 200)
    run_case[safe](262144, 20)
    run_case[safe](1048576, 6)


def main() raises:
    var arguments = argv()
    if len(arguments) == 2 and arguments[1] == "unsafe":
        run_cases[False]()
    elif len(arguments) == 2 and arguments[1] == "safe":
        run_cases[True]()
    else:
        print("usage: radix-sort-bench safe|unsafe", file=stderr)
        exit(2)
