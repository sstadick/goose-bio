"""Mojo port of biofast's `bedcov_c1_cgr.c` and the parts of `cgranges.c` and
`kseq.h` it uses: load one BED file into an implicit interval tree, then
print, for every interval of a second BED file, how many loaded intervals
overlap it and how many of its bases they cover.

The port follows the C line by line where that reads naturally and uses
`List`, `Dict`, and `Span` with their bounds checks. The standard library has
no buffered reader or writer, so the file ports `kstream_t` and buffers its
own output. Only the radix sort indexes through an unchecked pointer.
"""

from std.io.file_descriptor import RawFileDescriptor
from std.sys import argv, stdout


comptime BUFFER_SIZE = 0x10000

comptime TAB = UInt8(ord("\t"))
comptime CARRIAGE_RETURN = UInt8(ord("\r"))
comptime SPACE = UInt8(ord(" "))
comptime PLUS = UInt8(ord("+"))
comptime MINUS = UInt8(ord("-"))
comptime ZERO = UInt8(ord("0"))
comptime NINE = UInt8(ord("9"))

# Radix sort


comptime MIN_SIZE = 64
comptime RADIX_BITS = 8
comptime BUCKET_COUNT = 1 << RADIX_BITS

comptime IntervalPtr = Pointer[mut=True, Interval, MutUntrackedOrigin]


@always_inline("nodebug")
def digit(interval: Interval, shift: Int) -> Int:
    return Int((interval.x >> UInt64(shift)) & UInt64(BUCKET_COUNT - 1))


def insertion_sort(data: IntervalPtr, length: Int):
    for current in range(1, length):
        var temporary = data[unsafe_offset=current]
        if temporary.x < data[unsafe_offset=current - 1].x:
            var destination = current
            while destination > 0:
                var previous = data[unsafe_offset=destination - 1]
                if temporary.x >= previous.x:
                    break
                data[unsafe_offset=destination] = previous
                destination -= 1
            data[unsafe_offset=destination] = temporary


# `rs_sort_cr_intv`, with bucket cursors and ends as two tables of indices and
# the permutation as a `for` over buckets (see `bench/RADIX_SORT.md` for why
# this layout matches C's pointer pairs).
@no_inline
def radix_sort_inner(data: IntervalPtr, length: Int, shift: Int):
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


def radix_sort[origin: MutOrigin](intervals: Span[Interval, origin]):
    var data = intervals.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()
    if len(intervals) <= MIN_SIZE:
        insertion_sort(data, len(intervals))
    else:
        radix_sort_inner(data, len(intervals), 64 - RADIX_BITS)


# cgranges


@fieldwise_init
struct Interval(ImplicitlyCopyable):
    """`cr_intv_t`. Before indexing, `x` is `contig << 32 | start`; afterwards
    it is `start << 32 | end` and `y` is the largest end in the interval's
    subtree."""

    var x: UInt64
    var y: Int32
    var label: Int32
    """Carried, as in C, so a record is 16 bytes."""

    @always_inline
    def start(self) -> Int:
        return Int(Int32(self.x >> 32))

    @always_inline
    def end(self) -> Int:
        return Int(Int32(self.x & 0xFFFFFFFF))


@fieldwise_init
struct Contig(ImplicitlyCopyable):
    """`cr_ctg_t` without the name, which the map owns, and the unused length."""

    var root_k: Int
    """Height of the contig's tree, or -1 for a contig without intervals."""
    var offset: Int
    """The contig's intervals are `intervals[offset : offset + count]`."""
    var count: Int


@fieldwise_init
struct StackEntry(ImplicitlyCopyable):
    """One node of the traversal in `overlap`: the subtree rooted at index `x`
    with height `k`, and whether its left child has been handled."""

    var x: Int
    var k: Int
    var left_done: Bool


struct Cgranges(Movable):
    var intervals: List[Interval]
    var contigs: List[Contig]
    var contig_ids: Dict[String, Int]

    def __init__(out self):
        self.intervals = List[Interval]()
        self.contigs = List[Contig]()
        self.contig_ids = Dict[String, Int]()

    def add_contig(mut self, name: StringSlice) -> Int:
        """`cr_add_ctg`."""
        var key = String(name)
        var found = self.contig_ids.find(key)
        if found:
            return found.value()
        var id = len(self.contigs)
        self.contigs.append(Contig(-1, 0, 0))
        self.contig_ids[key^] = id
        return id

    def add(mut self, contig: StringSlice, start: Int, end: Int, label: Int):
        """`cr_add`."""
        if start > end:
            return
        var id = self.add_contig(contig)
        var x = UInt64(id) << 32 | (UInt64(start) & 0xFFFFFFFF)
        self.intervals.append(Interval(x, Int32(end), Int32(label)))

    def index(mut self):
        """`cr_index`: sort by contig and start, then build each contig's
        tree."""
        # cr_index_prepare
        var count = len(self.intervals)
        var is_sorted = True
        for i in range(1, count):
            if self.intervals[i - 1].x > self.intervals[i].x:
                is_sorted = False
                break
        if not is_sorted:
            radix_sort(self.intervals[0:count])
        var first = 0
        for i in range(1, count + 1):
            if (
                i == count
                or self.intervals[i].x >> 32 != self.intervals[first].x >> 32
            ):
                var id = Int(self.intervals[first].x >> 32)
                self.contigs[id].offset = first
                self.contigs[id].count = i - first
                first = i
        for ref interval in self.intervals:
            interval.x = interval.x << 32 | (UInt64(interval.y) & 0xFFFFFFFF)
            interval.y = 0

        for id in range(len(self.contigs)):
            var contig = self.contigs[id]
            self.contigs[id].root_k = index_contig(
                self.intervals[contig.offset : contig.offset + contig.count]
            )

    def overlap(
        self, contig_name: StringSlice, start: Int, end: Int, mut hits: List[Int]
    ):
        """`cr_overlap`: replace `hits` with the indices in `intervals`, in
        start order, of the intervals of the contig that overlap
        `start..end`."""
        hits.clear()
        var found = self.contig_ids.find(String(contig_name))
        if not found:
            return
        var contig = self.contigs[found.value()]
        if contig.root_k < 0:
            return
        var count = contig.count
        var tree = self.intervals[contig.offset : contig.offset + count]

        var stack = Array[StackEntry, 64](uninitialized=True)
        var depth = 1
        stack[0] = StackEntry((1 << contig.root_k) - 1, contig.root_k, False)
        while depth > 0:
            depth -= 1
            var node = stack[depth]
            if node.k <= 3:
                # The subtree has at most (1 << (k + 1)) - 1 intervals: scan it.
                var first = node.x >> node.k << node.k
                var last = min(first + (1 << (node.k + 1)) - 1, count)
                var index = contig.offset + first
                for interval in tree[first:last]:
                    if interval.start() >= end:
                        break
                    if start < interval.end():
                        hits.append(index)
                    index += 1
            elif not node.left_done:
                var left = node.x - (1 << (node.k - 1))
                stack[depth] = StackEntry(node.x, node.k, True)
                depth += 1
                if left >= count or Int(tree[left].y) > start:
                    stack[depth] = StackEntry(left, node.k - 1, False)
                    depth += 1
            elif node.x < count and tree[node.x].start() < end:
                if start < tree[node.x].end():
                    hits.append(contig.offset + node.x)
                var right = node.x + (1 << (node.k - 1))
                stack[depth] = StackEntry(right, node.k - 1, False)
                depth += 1


def index_contig[origin: MutOrigin](tree: Span[Interval, origin]) -> Int:
    """`cr_index1`: store in every `y` the largest end in that interval's
    subtree and return the tree's height."""
    var count = len(tree)
    if count == 0:
        return -1
    var last_i = 0
    var last = Int32(0)
    for i in range(0, count, 2):
        last_i = i
        last = Int32(tree[i].end())
        tree[i].y = last
    var k = 1
    while 1 << k <= count:
        var x = 1 << (k - 1)
        for i in range((x << 1) - 1, count, x << 2):
            var left = tree[i - x].y
            var right = tree[i + x].y if i + x < count else last
            tree[i].y = max(Int32(tree[i].end()), max(left, right))
        last_i = last_i - x if (last_i >> k) & 1 != 0 else last_i + x
        if last_i < count and tree[last_i].y > last:
            last = tree[last_i].y
        k += 1
    return k - 1


# kstream


struct LineReader(Movable):
    """`kstream_t` and the `KS_SEP_LINE` case of `ks_getuntil`."""

    var file: FileHandle
    var buffer: List[UInt8]
    var begin: Int
    var end: Int
    var is_eof: Bool

    def __init__(out self, path: String) raises:
        self.file = open(path, "r")
        self.buffer = List[UInt8](length=BUFFER_SIZE, fill=0)
        self.begin = 0
        self.end = 0
        self.is_eof = False

    def read_line(mut self, mut line: List[UInt8]) raises -> Bool:
        """Replace `line` with the next line, without its terminator. Returns
        whether the input had anything left to read."""
        line.clear()
        var got_any = False
        while True:
            if self.begin >= self.end:
                if self.is_eof:
                    break
                self.begin = 0
                self.end = self.file.read(self.buffer[0:BUFFER_SIZE])
                if self.end == 0:
                    self.is_eof = True
                    break
            got_any = True
            var rest = self.buffer[self.begin : self.end]
            var newline = StringSlice(unsafe_from_utf8=rest).find("\n")
            if newline < 0:
                line.extend(rest)
                self.begin = self.end
            else:
                line.extend(rest[0:newline])
                self.begin += newline + 1
                break
        if got_any and len(line) > 1 and line[len(line) - 1] == CARRIAGE_RETURN:
            _ = line.pop()
        return got_any


struct BufferedWriter(Movable, Writer):
    """Collects text and writes it to a file `BUFFER_SIZE` bytes at a time."""

    var file: RawFileDescriptor
    var buffer: List[UInt8]

    def __init__(out self, file: RawFileDescriptor):
        self.file = file
        self.buffer = List[UInt8](capacity=BUFFER_SIZE)

    def write_string(mut self, string: StringSlice):
        var bytes = string.as_bytes()
        if len(self.buffer) + len(bytes) > BUFFER_SIZE:
            self.flush()
        self.buffer.extend(bytes)

    def flush(mut self):
        self.file.write_bytes(self.buffer[0 : len(self.buffer)])
        self.buffer.clear()


# bedcov


def atoi(text: Span[Byte, _]) -> Int:
    """C's `atol` narrowed to 32 bits: optional blanks and sign, then
    digits."""
    var value = 0
    var negative = False
    var started = False  # Past the blanks: a sign or a digit has been read.
    for byte in text:
        if ZERO <= byte <= NINE:
            value = value * 10 + Int(byte - ZERO)
            started = True
        elif started:
            break
        elif byte == PLUS or byte == MINUS:
            negative = byte == MINUS
            started = True
        elif not (byte == SPACE or TAB <= byte <= CARRIAGE_RETURN):
            break
    return Int(Int32(-value if negative else value))


@fieldwise_init
struct Bed(ImplicitlyCopyable):
    """The fields of a BED line; its contig is the line's first
    `contig_length` bytes."""

    var contig_length: Int
    var start: Int
    var end: Int


def field_end(line: Span[Byte, _], begin: Int) -> Int:
    """The index of the first tab at or after `begin`, or the line's length."""
    var i = begin
    for byte in line[begin:]:
        if byte == TAB:
            break
        i += 1
    return i


def parse_bed(line: Span[Byte, _]) -> Optional[Bed]:
    """`parse_bed`: the contig, start, and end of a line with at least three
    tab-separated fields."""
    var contig_end = field_end(line, 0)
    if contig_end == len(line):
        return None
    var start_end = field_end(line, contig_end + 1)
    if start_end == len(line):
        return None
    var end_end = field_end(line, start_end + 1)
    return Bed(
        contig_end,
        atoi(line[contig_end + 1 : start_end]),
        atoi(line[start_end + 1 : end_end]),
    )


def read_bed(path: String) raises -> Cgranges:
    var ranges = Cgranges()
    var reader = LineReader(path)
    var line = List[UInt8]()
    var label = 0
    while reader.read_line(line):
        var parsed = parse_bed(line[0 : len(line)])
        if parsed:
            var bed = parsed.value()
            var contig = StringSlice(unsafe_from_utf8=line[0 : bed.contig_length])
            ranges.add(contig, bed.start, bed.end, label)
            label += 1
    return ranges^


def main() raises:
    var arguments = argv()
    if len(arguments) != 3:
        print("Usage: bedcov <loaded.bed> <streamed.bed>")
        return

    var ranges = read_bed(String(arguments[1]))
    ranges.index()

    var out = BufferedWriter(stdout)
    var reader = LineReader(String(arguments[2]))
    var line = List[UInt8]()
    var hits = List[Int]()
    while reader.read_line(line):
        var parsed = parse_bed(line[0 : len(line)])
        if not parsed:
            continue
        var bed = parsed.value()
        var contig = StringSlice(unsafe_from_utf8=line[0 : bed.contig_length])
        ranges.overlap(contig, bed.start, bed.end, hits)
        var coverage = 0
        var covered_start = 0
        var covered_end = 0
        for hit in hits:
            var interval = ranges.intervals[hit]
            var start0 = max(interval.start(), bed.start)
            var end0 = min(interval.end(), bed.end)
            if start0 > covered_end:
                coverage += covered_end - covered_start
                covered_start = start0
                covered_end = end0
            else:
                covered_end = max(covered_end, end0)
        coverage += covered_end - covered_start
        out.write(
            contig, "\t", bed.start, "\t", bed.end, "\t", len(hits), "\t",
            coverage, "\n",
        )
    out.flush()
