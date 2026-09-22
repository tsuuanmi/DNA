"""Reject explicit legacy/compatibility scaffolding in first-party Rust source."""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SOURCE_ROOT = ROOT / "src"

DEPRECATED_ATTRIBUTE = re.compile(r"#\s*\[\s*deprecated(?:\s*[=(]|\s*\])")
LINT_SUPPRESSION = re.compile(r"\b(?:allow|expect)\s*\(([^)]*)\)")
DECLARATION = re.compile(
    r"^\s*(?:pub(?:\([^)]*\))?\s+)?(?:async\s+)?"
    r"(?:fn|struct|enum|trait|type|mod|const|static)\s+"
    r"([A-Za-z_][A-Za-z0-9_]*)",
    re.MULTILINE,
)
USE_ALIAS = re.compile(
    r"^\s*(?:pub(?:\([^)]*\))?\s+)?use\s+[^;\n]*\bas\s+"
    r"([A-Za-z_][A-Za-z0-9_]*)",
    re.MULTILINE,
)
FEATURE_GATE = re.compile(
    r"feature\s*=\s*[\"']([^\"']*(?:legacy|compat|backward)[^\"']*)[\"']",
    re.IGNORECASE,
)

FORBIDDEN_SUPPRESSED_LINTS = {
    "dead_code",
    "deprecated",
    "unreachable_code",
    "unused",
    "unused_attributes",
    "unused_imports",
    "unused_mut",
    "unused_variables",
}

LEGACY_IDENTIFIER_PARTS = {"legacy", "compat", "compatibility"}


@dataclass(frozen=True)
class Violation:
    path: Path
    line: int
    rule: str
    detail: str

    def render(self) -> str:
        relative = (
            self.path.relative_to(ROOT) if self.path.is_relative_to(ROOT) else self.path
        )
        return f"{relative}:{self.line}: {self.rule}: {self.detail}"


def line_number(text: str, offset: int) -> int:
    return text.count("\n", 0, offset) + 1


def identifier_is_legacy(identifier: str) -> bool:
    lowered = identifier.lower()
    parts = lowered.split("_")
    return (
        any(part in LEGACY_IDENTIFIER_PARTS for part in parts)
        or "backward_compat" in lowered
        or "backwards_compat" in lowered
    )


def lint_names(match: re.Match[str]) -> set[str]:
    found: set[str] = set()
    for raw in match.group(1).split(","):
        lint = raw.strip().split("::")[-1]
        if lint in FORBIDDEN_SUPPRESSED_LINTS or lint.startswith("unused_"):
            found.add(lint)
    return found


def inspect_source(path: Path) -> list[Violation]:
    text = path.read_text(encoding="utf-8")
    violations: list[Violation] = []

    for match in DEPRECATED_ATTRIBUTE.finditer(text):
        violations.append(
            Violation(
                path,
                line_number(text, match.start()),
                "deprecated-api",
                "first-party production Rust must not declare #[deprecated] compatibility APIs",
            )
        )

    for match in LINT_SUPPRESSION.finditer(text):
        for lint in sorted(lint_names(match)):
            violations.append(
                Violation(
                    path,
                    line_number(text, match.start()),
                    "lint-suppression",
                    f"must not suppress {lint}; remove the stale code or fix the warning",
                )
            )

    for match in FEATURE_GATE.finditer(text):
        violations.append(
            Violation(
                path,
                line_number(text, match.start()),
                "compat-feature",
                f"compatibility feature gate {match.group(1)!r} is not permitted",
            )
        )

    for pattern in (DECLARATION, USE_ALIAS):
        for match in pattern.finditer(text):
            identifier = match.group(1)
            if identifier_is_legacy(identifier):
                violations.append(
                    Violation(
                        path,
                        line_number(text, match.start()),
                        "legacy-identifier",
                        f"explicit legacy/compatibility identifier {identifier!r} is not permitted",
                    )
                )

    return violations


def validate(root: Path = SOURCE_ROOT) -> list[Violation]:
    violations: list[Violation] = []
    for path in sorted(root.rglob("*.rs")):
        violations.extend(inspect_source(path))
    return violations


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=SOURCE_ROOT,
        help="Rust source root to validate (default: repository src/)",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    violations = validate(args.root)
    for violation in violations:
        print(f"FAIL: {violation.render()}", file=sys.stderr)
    if violations:
        print(
            f"{len(violations)} Rust source-policy violation(s) found",
            file=sys.stderr,
        )
        return 1
    print(
        "OK: first-party Rust source contains no explicit legacy compatibility scaffolding"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
