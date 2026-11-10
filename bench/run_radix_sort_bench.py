#!/usr/bin/env python3
"""Build the radix sort ports and compare each with the cgranges C reference.

Every port's executable links the same C reference, alternates the two sorts
in one process, and prints one `implementation,items,nanoseconds` row per
sample. Ratios therefore compare a port with the C samples taken beside it.
"""

import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
import csv
from dataclasses import dataclass
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


@dataclass(frozen=True)
class Port:
    tool: str  # Compiler option that builds it.
    target: str  # Executable it runs; ports of one language may share one.
    arguments: tuple[str, ...] = ()


PORTS = {
    "goose": Port("goose", "goose"),
    "rust-safe": Port("rustc", "rust", ("safe",)),
    "rust-unsafe": Port("rustc", "rust", ("unsafe",)),
    "zig": Port("zig", "zig"),
    "mojo-safe": Port("mojo", "mojo", ("safe",)),
    "mojo": Port("mojo", "mojo", ("unsafe",)),
    "mojo-opt": Port("mojo", "mojo-opt"),
}

VERSION_ARGUMENTS = {"zig": ["version"]}

# Elapsed nanoseconds of every sample, by port name and then by item count.
Samples = dict[str, dict[int, list[int]]]


def build_command(
    target: str, tool: list[str], reference: Path, cc_wrapper: Path, executable: Path
) -> list[str]:
    """Return the command that builds one target at its highest optimization level."""
    match target:
        case "goose":
            # The .build file beside the source compiles the C reference.
            source = BENCH / "radix_sort_bench.goose"
            return [
                *tool,
                "build",
                str(source),
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
                f"-Clink-arg={reference}",
                str(BENCH / "radix_sort_bench.rs"),
                "-o",
                str(executable),
            ]
        case "zig":
            source = BENCH / "radix_sort_bench.zig"
            return [
                *tool,
                "build-exe",
                str(source),
                str(reference),
                "-O",
                "ReleaseFast",
                "-mcpu",
                "native",
                "-lc",
                f"-femit-bin={executable}",
            ]
        case "mojo" | "mojo-opt":
            name = (
                "radix_sort_bench.mojo"
                if target == "mojo"
                else "radix_sort_bench_optimized.mojo"
            )
            return [
                *tool,
                "build",
                "-O3",
                "--march",
                "native",
                str(BENCH / name),
                "-Xlinker",
                str(reference),
                "-o",
                str(executable),
            ]
    raise ValueError(target)


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


def parse_arguments() -> tuple[argparse.Namespace, dict[str, list[str]], list[str]]:
    """Return the options, each compiler's command, and the ports to compare."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--ports",
        help=f"comma-separated ports from {', '.join(PORTS)} "
        "(default: every port whose compiler is found)",
    )
    parser.add_argument(
        "--rounds", type=int, default=5, help="runs of each executable (default: 5)"
    )
    parser.add_argument("--cpu", type=int, help="pin each run with taskset")
    parser.add_argument(
        "--cc",
        default=os.environ.get("CC", "gcc"),
        help="C compiler for the reference and for Goose's C (default: gcc)",
    )
    for tool in sorted({port.tool for port in PORTS.values()}):
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

    tools = {
        name: shlex.split(getattr(args, name))
        for name in ("cc", *(port.tool for port in PORTS.values()))
    }
    for name, tool in tools.items():
        if not tool:
            parser.error(f"--{name} must not be empty")
    found = {name for name, tool in tools.items() if shutil.which(tool[0])}
    if "cc" not in found:
        parser.error(f"C compiler {tools['cc'][0]} not found")
    if args.ports is None:
        ports = [name for name, port in PORTS.items() if port.tool in found]
        for name in sorted(tools.keys() - found):
            print(f"skipping {name} ports: {tools[name][0]} not found", file=sys.stderr)
        if not ports:
            parser.error("no port's compiler was found")
    else:
        ports = args.ports.split(",")
        for name in ports:
            if name not in PORTS:
                parser.error(f"unknown port {name}; choose from {', '.join(PORTS)}")
            if PORTS[name].tool not in found:
                tool = tools[PORTS[name].tool][0]
                parser.error(f"port {name} needs {tool}, which was not found")
    return args, tools, ports


def build(tools: dict[str, list[str]], ports: list[str], directory: Path) -> None:
    """Build each port's executable in directory, named after its target."""
    c_flags = ["-O3", "-march=native"]
    source = BENCH / "radix_sort_reference.c"
    reference = directory / "radix-sort-reference.o"
    run([*tools["cc"], *c_flags, "-c", str(source), "-o", str(reference)], ROOT)
    # Goose runs the C compiler itself, so the flags go in a wrapper.
    cc_wrapper = directory / "cc"
    cc_wrapper.write_text(
        f'#!/bin/sh\nexec {shlex.join(tools["cc"])} "$@" {shlex.join(c_flags)}\n'
    )
    cc_wrapper.chmod(0o755)

    commands = {
        port.target: build_command(
            port.target,
            tools[port.tool],
            reference,
            cc_wrapper,
            directory / port.target,
        )
        for port in map(PORTS.get, ports)
    }
    print(f"building {', '.join(commands)}", file=sys.stderr)
    with ThreadPoolExecutor() as pool:
        list(pool.map(lambda command: run(command, ROOT), commands.values()))


def report(ports: list[str], port_samples: Samples, c_samples: Samples) -> None:
    counts = sorted(port_samples[ports[0]])
    shape = {count: len(port_samples[ports[0]][count]) for count in counts}
    for name in ports:
        for samples in (port_samples[name], c_samples[name]):
            if {count: len(values) for count, values in samples.items()} != shape:
                sys.exit(f"port {name} does not take the same samples as {ports[0]}")

    c_cells = []
    for count in counts:
        pooled = statistics.median(
            [sample for name in ports for sample in c_samples[name][count]]
        )
        medians = [statistics.median(c_samples[name][count]) for name in ports]
        low, high = min(medians) / count, max(medians) / count
        c_cells.append(f"{pooled / count:.2f} ({low:.1f}-{high:.1f})")
    rows = [["c", *c_cells, ""]]
    for name in ports:
        cells = []
        ratios = []
        for count in counts:
            median = statistics.median(port_samples[name][count])
            ratios.append(median / statistics.median(c_samples[name][count]))
            cells.append(f"{median / count:.2f} ({ratios[-1]:.2f}x)")
        rows.append([name, *cells, f"{math.prod(ratios) ** (1 / len(ratios)):.3f}x"])

    print(table(["Port", *(f"{count:,}" for count in counts), "Geomean port/C"], rows))
    print()
    print(
        "Port cells are median ns/item, then the ratio to the C reference's median in the"
    )
    print(
        "same executable; below 1 is faster than C. The c row is the median over every"
    )
    print(
        "executable, then the range of the per-executable medians: the difference that"
    )
    print("code and data placement alone make to the same C code.")
    print(
        f"Samples per cell: {', '.join(f'{samples:,}' for samples in shape.values())}."
    )


def main() -> None:
    args, tools, ports = parse_arguments()
    port_samples: Samples = {name: defaultdict(list) for name in ports}
    c_samples: Samples = {name: defaultdict(list) for name in ports}

    with tempfile.TemporaryDirectory(prefix="goose-bio-radix-bench-") as tmp:
        build(tools, ports, Path(tmp))
        prefix = [] if args.cpu is None else ["taskset", "-c", str(args.cpu)]
        # Interleave the executables so slow drift affects every port alike.
        for round_number in range(1, args.rounds + 1):
            print(f"round {round_number}/{args.rounds}", file=sys.stderr)
            for name in ports:
                port = PORTS[name]
                output = run(
                    [*prefix, str(Path(tmp) / port.target), *port.arguments], ROOT
                )
                for implementation, count, elapsed in csv.reader(output.splitlines()):
                    samples = c_samples if implementation == "c" else port_samples
                    samples[name][int(count)].append(int(elapsed))

    report(ports, port_samples, c_samples)
    print()
    for name in dict.fromkeys(["cc", *(PORTS[name].tool for name in ports)]):
        print(f"- {name}: {version(name, tools[name])}")


if __name__ == "__main__":
    main()
