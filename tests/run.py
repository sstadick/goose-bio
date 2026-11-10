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
SUITES = ("radix_sort", "iitree", "fastx", "lapper", "bed")


def run(*args):
    print("+", " ".join(map(str, args)), flush=True)
    subprocess.run(list(map(str, args)), cwd=ROOT, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cc", nargs="+", required=True, help="native C compilers to test")
    args = parser.parse_args()
    goose = os.environ.get("GOOSE", "goose")
    print(f"Testing {platform.system()} {platform.machine()}", flush=True)
    for level in (0, 2):
        for suite in SUITES:
            run(goose, "run", f"tests/{suite}.goose", "--debug", f"-O{level}")
    with tempfile.TemporaryDirectory(prefix="goose-bio-tests-") as directory:
        exe = Path(directory) / ("test.exe" if os.name == "nt" else "test")
        for cc in args.cc:
            for suite in SUITES:
                run(goose, "build", f"tests/{suite}.goose", "--cc", cc,
                    "--debug", "-O2", "-o", exe)
                run(exe)


if __name__ == "__main__":
    main()
