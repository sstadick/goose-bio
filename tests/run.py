#!/usr/bin/env python3
"""Run library tests with TinyCC and native C compilers on the host architecture."""
import argparse
import os
from pathlib import Path
import platform
import struct
import subprocess
import tempfile
import zlib

ROOT = Path(__file__).resolve().parents[1]
# bed is an unfinished module; its test keeps it compiling.
SUITES = ("radix_sort", "iitree", "fastx", "faidx", "lapper", "bed")


def run(*args):
    print("+", " ".join(map(str, args)), flush=True)
    subprocess.run(list(map(str, args)), cwd=ROOT, check=True)


def compressed_fixture(directory):
    sequence = bytes(b"ACgTNrykmBDHV"[(i * 7 + i // 19) % 13] for i in range(131143))
    prefix = b">wrapped description\n"
    data = prefix + b"\r\n".join(sequence[i:i + 61] for i in range(0, len(sequence), 61))
    data += b"\n>other\n"; other_offset = len(data); data += b"nRy\n"
    path = Path(directory) / "wrapped.fa.gz"
    offsets = []; encoded = bytearray(); at = 0
    sizes = (1, 17, 65536, 3, 1024, 45000)
    while at < len(data):
        size = sizes[len(offsets) % len(sizes)]
        chunk = data[at:at + size]
        offsets.append((len(encoded), at))
        compressor = zlib.compressobj(wbits=-15)
        payload = compressor.compress(chunk) + compressor.flush()
        encoded += (b"\x1f\x8b\x08\x04" + bytes(6) + struct.pack("<H", 6) + b"BC" +
                    struct.pack("<HH", 2, len(payload) + 25) + payload +
                    struct.pack("<II", zlib.crc32(chunk), len(chunk)))
        at += len(chunk)
    path.write_bytes(encoded)  # A complete last block without an EOF marker is valid.
    Path(str(path) + ".fai").write_text(f"wrapped\t{len(sequence)}\t{len(prefix)}\t61\t63\nother\t3\t{other_offset}\t3\t4\n")
    Path(str(path) + ".gzi").write_bytes(struct.pack("<Q", len(offsets) - 1) +
                                       b"".join(struct.pack("<QQ", *pair) for pair in offsets[1:]))
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cc", nargs="+", required=True, help="native C compilers to test")
    args = parser.parse_args()
    goose = os.environ.get("GOOSE", "goose")
    print(f"Testing {platform.system()} {platform.machine()}", flush=True)
    with tempfile.TemporaryDirectory(prefix="goose-bio-tests-") as directory:
        packed = compressed_fixture(directory)
        fasta = Path(directory) / "large.fa"
        with fasta.open("wb") as output:
            # Mark the file sparse on Windows before extending it.
            if os.name == "nt":
                import ctypes
                import msvcrt
                returned = ctypes.c_ulong()
                handle = ctypes.c_void_p(msvcrt.get_osfhandle(output.fileno()))
                if not ctypes.windll.kernel32.DeviceIoControl(
                        handle, 0x900C4, None, 0, None, 0, ctypes.byref(returned), None):
                    raise ctypes.WinError()
            output.write(b">prefix\n")
            output.seek((1 << 32) + 17)
            output.write(b">high\n")
            offset = output.tell()
            output.write(b"AcGTn\nRYkm\n")
        Path(str(fasta) + ".fai").write_text(f"high\t9\t{offset}\t5\t6\n")
        for level in (0, 2):
            for suite in SUITES:
                extra = ["--", fasta, packed] if suite == "faidx" else []
                run(goose, "run", f"tests/{suite}.goose", "--debug", f"-O{level}", *extra)
        exe = Path(directory) / ("test.exe" if os.name == "nt" else "test")
        for cc in args.cc:
            for suite in SUITES:
                run(goose, "build", f"tests/{suite}.goose", "--cc", cc,
                    "--debug", "-O2", "-o", exe)
                run(exe, *([fasta, packed] if suite == "faidx" else []))


if __name__ == "__main__":
    main()
