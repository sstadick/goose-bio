#!/usr/bin/env python3
"""Build the biofast ports and compare each with biofast's C implementation.

A port is a pair of command-line programs, `bedcov` and `fqcnt`, that read
the same files and must print exactly what the C programs print. A case is
one program with its arguments. Each port's output is checked once against
C's, then every case is timed with its output discarded.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
import hashlib
import math
import os
from pathlib import Path
import shlex
import shutil
import statistics
import subprocess
import sys
import tempfile

BENCH = Path(__file__).resolve().parent
ROOT = BENCH.parent
SOURCES = BENCH / "biofast"
DATA_URL = (
    "https://github.com/lh3/biofast/releases/download/"
    "biofast-data-v1/biofast-data-v1.tar.gz"
)

# Compiler option that builds each port. C is the reference and always runs.
PORTS = {
    "c": "cc",
    "goose": "goose",
    "rust": "rustc",
    "zig": "zig",
    "mojo": "mojo",
}
PROGRAMS = ("bedcov", "fqcnt")


@dataclass(frozen=True)
class Case:
    program: str
    files: tuple[str, ...]
    description: str


CASES = {
    "g2r": Case(
        "bedcov",
        ("ex-rna.bed", "ex-anno.bed"),
        "index 8.9M RNA-seq intervals, query 1.2M annotations",
    ),
    "r2g": Case(
        "bedcov",
        ("ex-anno.bed", "ex-rna.bed"),
        "index 1.2M annotations, query 8.9M RNA-seq intervals",
    ),
    "fqcnt": Case(
        "fqcnt", ("M_abscessus_HiSeq.fq",), "parse 5.7M FASTQ records, 1.3 GB"
    ),
}

VERSION_ARGUMENTS = {"zig": ["version"]}


@dataclass(frozen=True)
class Sample:
    cpu: float  # User plus system seconds.
    peak: float  # Peak resident set size in MB.


def build_command(
    port: str, program: str, tool: list[str], cc_wrapper: Path, executable: Path
) -> list[str]:
    """Return the command that builds one program at its highest optimization level."""
    match port:
        case "c":
            name = "bedcov_c1_cgr.c" if program == "bedcov" else "fqcnt_c1_kseq.c"
            sources = [name, "cgranges.c"] if program == "bedcov" else [name]
            return [
                str(cc_wrapper),
                *(str(SOURCES / "c" / source) for source in sources),
                "-lz",
                "-o",
                str(executable),
            ]
        case "goose":
            return [
                *tool,
                "build",
                str(SOURCES / f"{program}.goose"),
                "--cc",
                str(cc_wrapper),
                "-O2",
                "-o",
                str(executable),
            ]
        case "rust":
            return [
                *tool,
                "--edition=2024",
                "-Copt-level=3",
                "-Ctarget-cpu=native",
                "-Ccodegen-units=1",
                str(SOURCES / f"{program}.rs"),
                "-o",
                str(executable),
            ]
        case "zig":
            return [
                *tool,
                "build-exe",
                str(SOURCES / f"{program}.zig"),
                "-O",
                "ReleaseFast",
                "-mcpu",
                "native",
                "-lc",
                f"-femit-bin={executable}",
            ]
        case "mojo":
            return [
                *tool,
                "build",
                "-O3",
                "--march",
                "native",
                str(SOURCES / f"{program}.mojo"),
                "-o",
                str(executable),
            ]
    raise ValueError(port)


def run(command: list[str], cwd: Path) -> str:
    result = subprocess.run(
        command, cwd=cwd, text=True, capture_output=True, check=False
    )
    if result.returncode != 0:
        sys.exit(f"{shlex.join(command)} failed:\n{result.stdout}{result.stderr}")
    return result.stdout


def version(name: str, tool: list[str]) -> str:
    output = run([*tool, *VERSION_ARGUMENTS.get(name, ["--version"])], ROOT)
    return output.splitlines()[0].strip()


def table(header: list[str], rows: list[list[str]]) -> str:
    """Format a Markdown table whose columns also line up as plain text."""
    widths = [max(len(row[i]) for row in [header, *rows]) for i in range(len(header))]
    lines = [
        "| "
        + " | ".join(
            cell.ljust(width) if i == 0 else cell.rjust(width)
            for i, (cell, width) in enumerate(zip(row, widths))
        )
        + " |"
        for row in [header, *rows]
    ]
    rule = (
        "|"
        + "|".join(
            ":" + "-" * (width + 1) if i == 0 else "-" * (width + 1) + ":"
            for i, width in enumerate(widths)
        )
        + "|"
    )
    return "\n".join([lines[0], rule, *lines[1:]])


def parse_arguments() -> tuple[
    argparse.Namespace, dict[str, list[str]], list[str], list[str]
]:
    """Return the options, each compiler's command, the ports, and the cases."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--ports",
        help=f"comma-separated ports from {', '.join(PORTS)} "
        "(default: every port whose compiler is found; c always runs)",
    )
    parser.add_argument(
        "--cases",
        default=",".join(CASES),
        help=f"comma-separated cases from {', '.join(CASES)} (default: all)",
    )
    parser.add_argument(
        "--rounds", type=int, default=3, help="timed runs of each case (default: 3)"
    )
    parser.add_argument("--cpu", type=int, help="pin each timed run with taskset")
    parser.add_argument(
        "--data",
        type=Path,
        default=BENCH / "data" / "biofast-data-v1",
        help="directory of the unpacked biofast-data-v1 files",
    )
    parser.add_argument(
        "--cc",
        default=os.environ.get("CC", "gcc"),
        help="C compiler for the C programs and for Goose's C (default: gcc)",
    )
    for tool in dict.fromkeys(PORTS.values()):
        if tool != "cc":
            parser.add_argument(
                f"--{tool}",
                default=os.environ.get(tool.upper(), tool),
                help=f"{tool} command; shell-style arguments are allowed",
            )
    args = parser.parse_args()
    if args.rounds < 1:
        parser.error("--rounds must be positive")
    if args.cpu is not None and shutil.which("taskset") is None:
        parser.error("--cpu requires taskset")

    cases = args.cases.split(",")
    for name in cases:
        if name not in CASES:
            parser.error(f"unknown case {name}; choose from {', '.join(CASES)}")
    missing = sorted(
        {
            file
            for name in cases
            for file in CASES[name].files
            if not (args.data / file).is_file()
        }
    )
    if missing:
        parser.error(
            f"{args.data} lacks {', '.join(missing)}; unpack {DATA_URL} "
            f"in {args.data.parent}, or pass --data"
        )

    tools = {tool: shlex.split(getattr(args, tool)) for tool in PORTS.values()}
    for name, tool in tools.items():
        if not tool:
            parser.error(f"--{name} must not be empty")
    found = {name for name, tool in tools.items() if shutil.which(tool[0])}
    if "cc" not in found:
        parser.error(f"C compiler {tools['cc'][0]} not found")
    if args.ports is None:
        ports = [name for name, tool in PORTS.items() if tool in found]
        for name in sorted(tools.keys() - found):
            print(f"skipping {name}: {tools[name][0]} not found", file=sys.stderr)
    else:
        ports = list(dict.fromkeys(["c", *args.ports.split(",")]))
        for name in ports:
            if name not in PORTS:
                parser.error(f"unknown port {name}; choose from {', '.join(PORTS)}")
            if PORTS[name] not in found:
                parser.error(
                    f"port {name} needs {tools[PORTS[name]][0]}, which was not found"
                )
    return args, tools, ports, cases


def build(
    tools: dict[str, list[str]], ports: list[str], programs: list[str], directory: Path
) -> None:
    """Build `measure` and each port's `<program>-<port>` executables in directory."""
    # Goose runs the C compiler itself, so the flags go in a wrapper.
    cc_wrapper = directory / "cc"
    cc_wrapper.write_text(
        f'#!/bin/sh\nexec {shlex.join(tools["cc"])} "$@" -O3 -march=native\n'
    )
    cc_wrapper.chmod(0o755)
    measure = directory / "measure"
    run([str(cc_wrapper), str(BENCH / "measure.c"), "-o", str(measure)], ROOT)

    def build_with(tool: str) -> None:
        # One compiler's programs build in turn: they share package and module
        # caches.
        for port in ports:
            for program in programs:
                if PORTS[port] == tool:
                    executable = directory / f"{program}-{port}"
                    run(
                        build_command(
                            port, program, tools[tool], cc_wrapper, executable
                        ),
                        ROOT,
                    )

    print(f"building {', '.join(ports)}", file=sys.stderr)
    with ThreadPoolExecutor() as pool:
        list(pool.map(build_with, dict.fromkeys(PORTS[port] for port in ports)))


def command(directory: Path, port: str, case: Case, data: Path) -> list[str]:
    return [
        str(directory / f"{case.program}-{port}"),
        *(str(data / file) for file in case.files),
    ]


def output_digest(command: list[str]) -> str:
    """Run command and return a digest of what it wrote to standard output."""
    digest = hashlib.md5()
    with subprocess.Popen(command, stdout=subprocess.PIPE) as process:
        while chunk := process.stdout.read(1 << 20):
            digest.update(chunk)
    if process.returncode != 0:
        sys.exit(f"{shlex.join(command)} failed with status {process.returncode}")
    return digest.hexdigest()


def check(directory: Path, ports: list[str], cases: list[str], data: Path) -> None:
    """Exit unless every port prints what C prints for every case."""
    print("checking every port's output against c", file=sys.stderr)
    jobs = [(port, name) for name in cases for port in ports]
    with ThreadPoolExecutor(max_workers=os.cpu_count()) as pool:
        digests = dict(
            zip(
                jobs,
                pool.map(
                    lambda job: output_digest(
                        command(directory, job[0], CASES[job[1]], data)
                    ),
                    jobs,
                ),
            )
        )
    wrong = [
        f"{port} {name}"
        for port, name in jobs
        if digests[port, name] != digests["c", name]
    ]
    if wrong:
        sys.exit(f"output differs from c: {', '.join(wrong)}")


def measure(directory: Path, prefix: list[str], command: list[str]) -> Sample:
    """Run command with its output discarded and return what it cost."""
    # measure.c explains why the child's resource usage is not read from here.
    report = directory / "report"
    command = [*prefix, str(directory / "measure"), str(report), *command]
    status = subprocess.run(command, stdout=subprocess.DEVNULL, check=False).returncode
    if status != 0:
        sys.exit(f"{shlex.join(command)} failed with status {status}")
    user, system, peak_kb = report.read_text().split()
    return Sample(float(user) + float(system), int(peak_kb) / 1024)


def report(
    ports: list[str], cases: list[str], samples: dict[tuple[str, str], list[Sample]]
) -> None:
    def median(port: str, case: str, field: str) -> float:
        return statistics.median(
            getattr(sample, field) for sample in samples[port, case]
        )

    def cpu(port: str, case: str) -> str:
        seconds = median(port, case, "cpu")
        return f"{seconds:.2f} ({seconds / median('c', case, 'cpu'):.2f}x)"

    rows = []
    for port in ports:
        ratios = [
            median(port, case, "cpu") / median("c", case, "cpu") for case in cases
        ]
        geomean = (
            "" if port == "c" else f"{math.prod(ratios) ** (1 / len(ratios)):.3f}x"
        )
        rows.append([port, *(cpu(port, case) for case in cases), geomean])
    print(table(["Port", *cases, "Geomean port/C"], rows))
    print()
    print("Median CPU seconds (user plus system), then the ratio to c; below 1 is")
    print(f"faster than C. Runs per cell: {len(samples[ports[0], cases[0]])}.")
    print()

    def peak(port: str, case: str) -> str:
        return f"{median(port, case, 'peak'):.1f}"

    rows = [[port, *(peak(port, case) for case in cases)] for port in ports]
    print(table(["Peak memory, MB", *cases], rows))
    print()

    for name in cases:
        case = CASES[name]
        print(f"- {name}: `{case.program} {' '.join(case.files)}`, {case.description}")


def main() -> None:
    args, tools, ports, cases = parse_arguments()
    programs = [
        program
        for program in PROGRAMS
        if any(CASES[name].program == program for name in cases)
    ]
    samples: dict[tuple[str, str], list[Sample]] = {
        (port, name): [] for port in ports for name in cases
    }

    with tempfile.TemporaryDirectory(prefix="goose-bio-biofast-bench-") as tmp:
        directory = Path(tmp)
        build(tools, ports, programs, directory)
        check(directory, ports, cases, args.data)
        prefix = [] if args.cpu is None else ["taskset", "-c", str(args.cpu)]
        # Interleave the ports so slow drift affects every port alike.
        for round_number in range(1, args.rounds + 1):
            print(f"round {round_number}/{args.rounds}", file=sys.stderr)
            for name in cases:
                for port in ports:
                    samples[port, name].append(
                        measure(
                            directory,
                            prefix,
                            command(directory, port, CASES[name], args.data),
                        )
                    )

    report(ports, cases, samples)
    print()
    for name in dict.fromkeys(PORTS[port] for port in ports):
        print(f"- {name}: {version(name, tools[name])}")


if __name__ == "__main__":
    main()
