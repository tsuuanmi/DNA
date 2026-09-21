"""Validate GitHub Actions supply-chain policy."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
WORKFLOW_ROOT = ROOT / ".github" / "workflows"
USES = re.compile(r"^\s*(?:-\s*)?uses:\s*([^@\s]+)@([^\s#]+)")
FULL_SHA = re.compile(r"[0-9a-f]{40}")
FORBIDDEN = (
    "pull_request_target:",
    "workflow_run:",
    "permissions: write-all",
)


def violations(root: Path) -> list[str]:
    found: list[str] = []
    for path in sorted(root.glob("*.y*ml")):
        for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            stripped = line.strip()
            for token in FORBIDDEN:
                if token in stripped:
                    found.append(\n                        f"{path}:{number}: forbidden workflow construct {token!r}"\n                    )
            match = USES.match(line)
            if match and not FULL_SHA.fullmatch(match.group(2)):
                found.append(
                    f"{path}:{number}: action {match.group(1)!r} must use a full 40-character commit SHA"
                )
    return found


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=WORKFLOW_ROOT)
    args = parser.parse_args(argv)

    errors = violations(args.root)
    for error in errors:
        print(f"FAIL: {error}", file=sys.stderr)
    if errors:
        return 1
    print("OK: GitHub Actions use immutable pins and safe trigger policy")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
