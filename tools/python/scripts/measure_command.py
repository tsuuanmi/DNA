"""Measure one Linux command for release-validation evidence."""

from __future__ import annotations

import argparse
import json
import os
import platform
import resource
import subprocess
import sys
import time
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

SCHEMA_VERSION = "dna.performance/v1"


def parser() -> argparse.ArgumentParser:
    built = argparse.ArgumentParser(description=__doc__)
    built.add_argument(
        "--output",
        type=Path,
        required=True,
        help="JSON evidence path",
    )
    built.add_argument(
        "command",
        nargs=argparse.REMAINDER,
        help="command to measure; prefix with -- to separate wrapper arguments",
    )
    return built


def normalized_command(values: list[str]) -> list[str]:
    command = values[1:] if values and values[0] == "--" else values
    if not command:
        raise ValueError("a command is required after --")
    return command


def measure(command: list[str]) -> tuple[int, dict[str, Any]]:
    if sys.platform != "linux":
        raise RuntimeError("performance evidence measurement currently supports Linux only")

    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    started_at = datetime.now(UTC)
    started = time.perf_counter()
    completed = subprocess.run(command, check=False)
    elapsed = time.perf_counter() - started
    after = resource.getrusage(resource.RUSAGE_CHILDREN)

    evidence: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "measured_at": started_at.isoformat(),
        "platform": {
            "system": platform.system(),
            "release": platform.release(),
            "machine": platform.machine(),
        },
        "command": command,
        "exit_code": completed.returncode,
        "wall_seconds": round(elapsed, 6),
        "user_cpu_seconds": round(after.ru_utime - before.ru_utime, 6),
        "system_cpu_seconds": round(after.ru_stime - before.ru_stime, 6),
        "peak_rss_kib": after.ru_maxrss,
    }
    return completed.returncode, evidence


def write_evidence(path: Path, evidence: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(
        json.dumps(evidence, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        command = normalized_command(args.command)
        status, evidence = measure(command)
        write_evidence(args.output, evidence)
    except (OSError, RuntimeError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    print(args.output)
    return status


if __name__ == "__main__":
    raise SystemExit(main())
