#!/usr/bin/env python3
"""Run library tests with TinyCC and native C compilers on the host architecture."""
import argparse
import os
from pathlib import Path
import platform
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
# bed is an unfinished module; its test keeps it compiling.
SUITES = ("radix_sort", "iitree", "fastx", "faidx", "lapper", "bed")


def run(*args):
    print("+", " ".join(map(str, args)), flush=True)
    subprocess.run(list(map(str, args)), cwd=ROOT, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cc", nargs="+", required=True, help="native C compilers to test")
    args = parser.parse_args()
    goose = os.environ.get("GOOSE", "goose")
    print(f"Testing {platform.system()} {platform.machine()}", flush=True)
    with tempfile.TemporaryDirectory(prefix="goose-bio-tests-") as directory:
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
                extra = ["--", fasta] if suite == "faidx" else []
                run(goose, "run", f"tests/{suite}.goose", "--debug", f"-O{level}", *extra)
        exe = Path(directory) / ("test.exe" if os.name == "nt" else "test")
        for cc in args.cc:
            for suite in SUITES:
                run(goose, "build", f"tests/{suite}.goose", "--cc", cc,
                    "--debug", "-O2", "-o", exe)
                run(exe, *([fasta] if suite == "faidx" else []))


if __name__ == "__main__":
    main()
