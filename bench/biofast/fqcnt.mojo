"""Mojo port of biofast's `fqcnt_c1_kseq.c` and the FASTA/FASTQ parser of
`kseq.h`: count the records, bases, and quality values of a file.

The standard library has no buffered reader, so `KStream` ports `kstream_t`;
`Kseq.read` follows `kseq_read`, including multi-line records. Buffers are
`List`s indexed with bounds checks.
"""

from std.sys import argv, stderr


comptime BUFFER_SIZE = 16384

comptime END_OF_FILE = -1
"""What `KStream.getc` and the delimiter results return at the end of the
input, as in C."""

comptime NEWLINE = ord("\n")
comptime CARRIAGE_RETURN = UInt8(ord("\r"))
comptime SPACE = UInt8(ord(" "))
comptime TAB = UInt8(ord("\t"))
comptime FASTA_HEADER = ord(">")
comptime FASTQ_HEADER = ord("@")
comptime QUALITY_HEADER = ord("+")


struct KStream(Movable):
    """`kstream_t`: a file read through a buffer."""

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

    @always_inline
    def fill(mut self) raises -> Bool:
        """Read more when the buffer is used up. Returns whether it then
        holds unread bytes."""
        if self.begin < self.end:
            return True
        return self.refill()

    @no_inline
    def refill(mut self) raises -> Bool:
        if self.is_eof:
            return False
        self.begin = 0
        self.end = self.file.read(self.buffer[0:BUFFER_SIZE])
        self.is_eof = self.end == 0
        return not self.is_eof

    @always_inline
    def getc(mut self) raises -> Int:
        """`ks_getc`: the next byte, or `END_OF_FILE`."""
        if not self.fill():
            return END_OF_FILE
        var byte = self.buffer[self.begin]
        self.begin += 1
        return Int(byte)

    def read_word(mut self, mut text: List[UInt8]) raises -> Int:
        """`ks_getuntil2` with `KS_SEP_SPACE`: append the bytes before the
        next whitespace to `text` and return that whitespace, or
        `END_OF_FILE`."""
        while self.fill():
            var i = self.begin
            for byte in self.buffer[self.begin : self.end]:
                if byte == SPACE or TAB <= byte <= CARRIAGE_RETURN:
                    break
                i += 1
            text.extend(self.buffer[self.begin : i])
            self.begin = i + 1
            if i < self.end:
                return Int(self.buffer[i])
        return END_OF_FILE

    def read_line(mut self, mut text: List[UInt8]) raises -> Bool:
        """`ks_getuntil2` with `KS_SEP_LINE`: append the rest of the line to
        `text`, without its terminator. Returns whether the input had
        anything left to read."""
        var got_any = False
        while self.fill():
            got_any = True
            var rest = self.buffer[self.begin : self.end]
            var newline = StringSlice(unsafe_from_utf8=rest).find("\n")
            if newline >= 0:
                text.extend(rest[0:newline])
                self.begin += newline + 1
                break
            text.extend(rest)
            self.begin = self.end
        if got_any and len(text) > 1 and text[len(text) - 1] == CARRIAGE_RETURN:
            _ = text.pop()
        return got_any


struct Kseq(Movable):
    """`kseq_t`: the most recently read record and the stream it came from."""

    var name: List[UInt8]
    var comment: List[UInt8]
    var seq: List[UInt8]
    var qual: List[UInt8]
    var last_char: Int
    """The `>` or `@` of the next record when reading a sequence consumed it,
    otherwise 0."""
    var stream: KStream

    def __init__(out self, var stream: KStream):
        self.name = List[UInt8]()
        self.comment = List[UInt8]()
        self.seq = List[UInt8](capacity=256)
        self.qual = List[UInt8]()
        self.last_char = 0
        self.stream = stream^

    def read(mut self) raises -> Bool:
        """`kseq_read`: read the next record into `self`. Returns False at the
        end of the input and raises for a FASTQ record whose quality string
        is missing or not as long as its sequence."""
        if self.last_char == 0:
            # Jump to the next header line.
            while True:
                var byte = self.stream.getc()
                if byte == END_OF_FILE:
                    return False
                if byte == FASTA_HEADER or byte == FASTQ_HEADER:
                    self.last_char = byte
                    break
        self.name.clear()
        self.comment.clear()
        self.seq.clear()
        self.qual.clear()

        var delimiter = self.stream.read_word(self.name)
        if delimiter == END_OF_FILE and len(self.name) == 0:
            return False
        if delimiter != NEWLINE:
            _ = self.stream.read_line(self.comment)

        var next: Int
        while True:
            next = self.stream.getc()
            if (
                next == END_OF_FILE
                or next == FASTA_HEADER
                or next == QUALITY_HEADER
                or next == FASTQ_HEADER
            ):
                break
            if next == NEWLINE:
                continue  # Skip empty lines.
            self.seq.append(UInt8(next))
            _ = self.stream.read_line(self.seq)
        if next == FASTA_HEADER or next == FASTQ_HEADER:
            self.last_char = next  # The next record's header has begun.
        if next != QUALITY_HEADER:
            return True  # FASTA

        # Skip the rest of the '+' line.
        while True:
            var byte = self.stream.getc()
            if byte == END_OF_FILE:
                raise Error("truncated quality string")
            if byte == NEWLINE:
                break
        while self.stream.read_line(self.qual) and len(self.qual) < len(self.seq):
            pass
        self.last_char = 0  # The next header line is still ahead.
        if len(self.seq) != len(self.qual):
            raise Error("truncated quality string")
        return True


def main() raises:
    var arguments = argv()
    if len(arguments) != 2:
        print("Usage: fqcnt <in.fq>", end="")
        return

    var record = Kseq(KStream(String(arguments[1])))
    var count = 0
    var seq_length = 0
    var qual_length = 0
    var malformed = False
    try:
        while record.read():
            count += 1
            seq_length += len(record.seq)
            qual_length += len(record.qual)
    except:
        malformed = True
    print(count, seq_length, qual_length, sep="\t")
    if malformed:
        print("ERROR: malformated FASTX", file=stderr)
