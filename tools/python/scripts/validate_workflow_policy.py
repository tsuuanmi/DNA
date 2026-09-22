"""Validate GitHub Actions supply-chain and aggregate-gate policy."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
WORKFLOW_ROOT = ROOT / ".github" / "workflows"
USES = re.compile(r"^\s*(?:-\s*)?uses:\s*([^@\s]+)@([^\s#]+)")
FULL_SHA = re.compile(r"[0-9a-f]{40}")
JOB_ID = re.compile(r"^  ([A-Za-z_][A-Za-z0-9_-]*):\s*$")
NEED = re.compile(r"^      - ([A-Za-z_][A-Za-z0-9_-]*)\s*$")
FORBIDDEN = (
    "pull_request_target:",
    "workflow_run:",
    "permissions: write-all",
)
AGGREGATE_JOB = "ci-success"


def action_violations(root: Path) -> list[str]:
    found: list[str] = []
    for path in sorted(root.glob("*.y*ml")):
        for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            stripped = line.strip()
            for token in FORBIDDEN:
                if token in stripped:
                    found.append(
                        f"{path}:{number}: forbidden workflow construct {token!r}"
                    )
            match = USES.match(line)
            if match and not FULL_SHA.fullmatch(match.group(2)):
                found.append(
                    f"{path}:{number}: action {match.group(1)!r} must use a full "
                    "40-character commit SHA"
                )
    return found


def aggregate_violations(path: Path) -> list[str]:
    if not path.is_file():
        return [f"{path}: CI workflow is missing"]

    lines = path.read_text(encoding="utf-8").splitlines()
    if "jobs:" not in lines:
        return []

    jobs_index = lines.index("jobs:")
    job_lines = lines[jobs_index + 1 :]
    jobs = {match.group(1) for line in job_lines if (match := JOB_ID.match(line))}
    if AGGREGATE_JOB not in jobs:
        return [f"{path}: missing aggregate job {AGGREGATE_JOB!r}"]

    aggregate_index = lines.index(f"  {AGGREGATE_JOB}:")
    needs_index = next(
        (
            index
            for index in range(aggregate_index + 1, len(lines))
            if lines[index] == "    needs:"
        ),
        None,
    )
    if needs_index is None:
        return [f"{path}: aggregate job {AGGREGATE_JOB!r} has no needs list"]

    needs: set[str] = set()
    for line in lines[needs_index + 1 :]:
        match = NEED.match(line)
        if match:
            needs.add(match.group(1))
            continue
        if not line.strip():
            continue
        break

    mandatory = jobs - {AGGREGATE_JOB}
    missing = sorted(mandatory - needs)
    unexpected = sorted(needs - mandatory)
    errors: list[str] = []
    if missing:
        errors.append(
            f"{path}: aggregate job is missing mandatory dependencies: "
            + ", ".join(missing)
        )
    if unexpected:
        errors.append(
            f"{path}: aggregate job references unknown dependencies: "
            + ", ".join(unexpected)
        )
    return errors


def violations(root: Path) -> list[str]:
    return action_violations(root) + aggregate_violations(root / "ci.yml")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=WORKFLOW_ROOT)
    args = parser.parse_args(argv)

    errors = violations(args.root)
    for error in errors:
        print(f"FAIL: {error}", file=sys.stderr)
    if errors:
        return 1
    print(
        "OK: GitHub Actions use immutable pins, safe triggers, and complete CI aggregation"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
