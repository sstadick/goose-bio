from std.ffi import external_call
from std.memory import unsafe_memcpy, unsafe_memset_zero
from std.memory.alloc import unsafe_alloc
from std.sys import size_of
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


@fieldwise_init
struct Bucket(ImplicitlyCopyable):
    var begin: Int
    var end: Int


comptime BucketPtr = Pointer[mut=True, Bucket, MutUntrackedOrigin]
comptime RADIX_LEVELS = 8


@always_inline("nodebug")
def digit(record: Record, shift: Int) -> Int:
    return Int((record.key >> UInt64(shift)) & UInt64(BUCKET_COUNT - 1))


def insertion_sort(begin: RecordPtr, end: RecordPtr):
    var current = begin.unsafe_offset(1)
    while current < end:
        if current[].key < current.unsafe_offset(-1)[].key:
            var temporary = current[]
            var destination = current
            while (
                destination > begin
                and temporary.key < destination.unsafe_offset(-1)[].key
            ):
                destination[] = destination.unsafe_offset(-1)[]
                destination = destination.unsafe_offset(-1)
            destination[] = temporary
        current = current.unsafe_offset(1)


@no_inline
def radix_sort_inner(
    data: RecordPtr,
    length: Int,
    shift: Int,
    workspace: BucketPtr,
    depth: Int,
):
    # Keep one table per radix depth in a sort-owned heap workspace. This
    # removes the 4 KiB local Array from every recursive stack frame.
    var bucket_ptr = workspace.unsafe_offset(depth * BUCKET_COUNT)
    unsafe_memset_zero(bucket_ptr, BUCKET_COUNT)
    var buckets = Span(unsafe_ptr=bucket_ptr, length=BUCKET_COUNT)

    var item = data
    var item_end = data.unsafe_offset(length)
    while item != item_end:
        buckets.unsafe_get(digit(item[], shift)).end += 1
        item = item.unsafe_offset(1)

    var bucket = 0
    var cumulative = 0
    var previous_end = Int(data)
    while bucket < BUCKET_COUNT:
        cumulative += buckets.unsafe_get(bucket).end
        buckets.unsafe_get(bucket).begin = previous_end
        previous_end = Int(data.unsafe_offset(cumulative))
        buckets.unsafe_get(bucket).end = previous_end
        bucket += 1

    bucket = 0
    while bucket < BUCKET_COUNT:
        var cursor = buckets.unsafe_get(bucket).begin
        if cursor != buckets.unsafe_get(bucket).end:
            var cursor_ptr = RecordPtr(unsafe_from_address=cursor)
            var target = digit(cursor_ptr[], shift)
            if target != bucket:
                var temporary = cursor_ptr[]
                while True:
                    var destination = buckets.unsafe_get(target).begin
                    var destination_ptr = RecordPtr(
                        unsafe_from_address=destination
                    )
                    var displaced = destination_ptr[]
                    destination_ptr[] = temporary
                    buckets.unsafe_get(target).begin = Int(
                        destination_ptr.unsafe_offset(1)
                    )
                    temporary = displaced
                    target = digit(temporary, shift)
                    if target == bucket:
                        break
                cursor_ptr[] = temporary
            buckets.unsafe_get(bucket).begin = Int(cursor_ptr.unsafe_offset(1))
        else:
            bucket += 1

    buckets.unsafe_get(0).begin = Int(data)
    bucket = 1
    while bucket < BUCKET_COUNT:
        buckets.unsafe_get(bucket).begin = buckets.unsafe_get(bucket - 1).end
        bucket += 1

    if shift != 0:
        var next_shift = max(0, shift - RADIX_BITS)
        bucket = 0
        while bucket < BUCKET_COUNT:
            var bucket_begin = buckets.unsafe_get(bucket).begin
            var bucket_end = buckets.unsafe_get(bucket).end
            var bucket_length = (bucket_end - bucket_begin) // size_of[Record]()
            var bucket_data = RecordPtr(unsafe_from_address=bucket_begin)
            if bucket_length > MIN_SIZE:
                radix_sort_inner(
                    bucket_data,
                    bucket_length,
                    next_shift,
                    workspace,
                    depth + 1,
                )
            elif bucket_length > 1:
                insertion_sort(
                    bucket_data, RecordPtr(unsafe_from_address=bucket_end)
                )
            bucket += 1


@always_inline
def radix_sort(data: RecordPtr, length: Int):
    if length <= 1:
        return
    if length <= MIN_SIZE:
        insertion_sort(data, data.unsafe_offset(length))
    else:
        var workspace = unsafe_alloc[Bucket](RADIX_LEVELS * BUCKET_COUNT)
        radix_sort_inner(data, length, 64 - RADIX_BITS, workspace, 0)
        workspace.unsafe_free()


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


def run_case(length: Int, repetitions: Int) raises:
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
    radix_sort(work, length)
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
            radix_sort(work, length)
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
            radix_sort(work, length)
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


def main() raises:
    run_case(1024, 1000)
    run_case(16384, 200)
    run_case(262144, 20)
    run_case(1048576, 6)
