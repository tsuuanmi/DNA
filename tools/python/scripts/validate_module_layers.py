"""Enforce the plugin-family crate map of first-party Rust source (ADR-0069).

Every module belongs to one crate of the plugin-first workspace:

- ``kernel``: shared contracts;
- ``core``: the core caller;
- ``sanger``: the Sanger modality;
- ``post``: post-calling plugins;
- ``dna``: the facade that composes them.

``model`` children and ``input`` children are assigned individually. A module
may depend only on its own crate or on a crate its crate is allowed to depend
on: the plugin crates depend only on the kernel, and only the facade composes
them. The module graph must also be acyclic.
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SOURCE_ROOT = ROOT / "src"

# Crates each crate may depend on, besides itself.
ALLOWED: dict[str, frozenset[str]] = {
    "kernel": frozenset(),
    "core": frozenset({"kernel"}),
    "sanger": frozenset({"kernel"}),
    "post": frozenset({"kernel"}),
    "dna": frozenset({"kernel", "core", "sanger", "post"}),
}

# Crate of every module path under the source root; the longest matching
# prefix wins. Keep in sync with PROP-0002 phase 5.
MODULES: dict[str, str] = {
    "error": "kernel",
    "checksum": "kernel",
    "bounds": "kernel",
    "read_evidence": "kernel",
    "plugin": "kernel",
    "variant": "kernel",
    "reference": "kernel",
    "profile": "kernel",
    "model::nucleotide": "kernel",
    "model::reference": "kernel",
    "alignment": "core",
    "variant_calling": "core",
    "sample": "core",
    "read_call": "core",
    "model::alignment": "core",
    "model::variant": "core",
    "model::coordinate": "core",
    "model::called_read": "core",
    "model::reference_call": "core",
    "model::sample_evidence": "core",
    "locus": "sanger",
    "basecalling": "sanger",
    "signal_processing": "sanger",
    "callability": "sanger",
    "quality_control": "sanger",
    "read_processing": "sanger",
    "input::sanger::abif": "sanger",
    "model::basecalls": "sanger",
    "model::callability": "sanger",
    "model::locus_evidence": "sanger",
    "model::quality": "sanger",
    "model::sanger": "sanger",
    "model::signal": "sanger",
    "model::attachment": "sanger",
    "variant_representation": "post",
    "variant_normalization": "post",
    "variant_nomenclature": "post",
    "conformance": "post",
    "cli": "dna",
    "pipeline": "dna",
    "report": "dna",
    "operation_log": "dna",
    "config": "dna",
    "input": "dna",
    "variant_analysis": "dna",
    "model::read_observation": "dna",
    "model::result": "dna",
    "model::sample_result": "dna",
    "model::basecall_result": "dna",
    "model::variants_result": "dna",
    "model::notation_result": "dna",
}

# Files that only declare child modules or compose the crate.
CONTAINERS = {"lib.rs", "main.rs", "model/mod.rs"}

TEST_MODULE = re.compile(
    r"#\[cfg\(test\)\]\s*(?:pub(?:\([^)]*\))?\s+)?mod\s+\w+\s*\{", re.MULTILINE
)
LINE_COMMENT = re.compile(r"//.*$", re.MULTILINE)
CRATE_USE = re.compile(r"\bcrate::")
TOKEN = re.compile(r"\s*([A-Za-z_][A-Za-z0-9_]*|::|\{|\}|,|\*)")


@dataclass(frozen=True)
class Edge:
    source: str
    target: str
    path: Path
    line: int


def module_path(root: Path, path: Path) -> str:
    parts = list(path.relative_to(root).with_suffix("").parts)
    if parts[-1] == "mod":
        parts.pop()
    return "::".join(parts)


def unit(path: str) -> str | None:
    """The longest mapped module prefix of a ``::`` path, or None."""
    segments = path.split("::")
    for length in range(len(segments), 0, -1):
        prefix = "::".join(segments[:length])
        if prefix in MODULES:
            return prefix
    return None


def production_text(path: Path) -> str:
    """Source without inline test modules (which run to the end of the file) or comments."""
    text = path.read_text(encoding="utf-8")
    test_module = TEST_MODULE.search(text)
    if test_module:
        text = text[: test_module.start()]
    return LINE_COMMENT.sub("", text)


def path_tokens(text: str, start: int) -> list[str]:
    """Tokens of the path that begins at ``start``, through any ``{...}`` group."""
    tokens: list[str] = []
    position = start
    depth = 0
    expect_segment = True
    while match := TOKEN.match(text, position):
        token = match.group(1)
        if token == "{":
            depth += 1
        elif token == "}":
            if depth == 0:
                break
            depth -= 1
        elif token == "::":
            expect_segment = True
        elif token == ",":
            if depth == 0:
                break
            expect_segment = True
        elif not expect_segment and depth == 0:
            break
        else:
            expect_segment = False
        tokens.append(token)
        position = match.end()
        if depth == 0 and token == "}":
            break
    return tokens


def expand(tokens: list[str]) -> list[str]:
    """Full module paths of a path with nested ``{...}`` groups."""
    paths: list[str] = []

    def walk(index: int, prefix: str) -> int:
        current = prefix
        while index < len(tokens):
            token = tokens[index]
            if token == "::":
                index += 1
            elif token == "{":
                index += 1
                while index < len(tokens) and tokens[index] != "}":
                    index = walk(index, current)
                    if index < len(tokens) and tokens[index] == ",":
                        index += 1
                return index + 1
            elif token in {",", "}"}:
                break
            else:
                if token not in {"self", "*"}:
                    current = f"{current}::{token}" if current else token
                index += 1
        paths.append(current)
        return index

    walk(0, "")
    return [path for path in paths if path]


def edges(root: Path) -> list[Edge]:
    found: list[Edge] = []
    for path in sorted(root.rglob("*.rs")):
        if path.relative_to(root).as_posix() in CONTAINERS or path.name == "tests.rs":
            continue
        source = module_path(root, path)
        text = production_text(path)
        for match in CRATE_USE.finditer(text):
            line = text.count("\n", 0, match.start()) + 1
            for target in expand(path_tokens(text, match.end())):
                found.append(Edge(source, target, path, line))
    return found


def find_cycle(graph: dict[str, set[str]]) -> list[str] | None:
    """Return one dependency cycle as a module path, or None when acyclic."""
    state: dict[str, int] = {}
    stack: list[str] = []

    def visit(module: str) -> list[str] | None:
        state[module] = 1
        stack.append(module)
        for target in sorted(graph.get(module, set())):
            if state.get(target) == 1:
                return [*stack[stack.index(target) :], target]
            if target not in state and (cycle := visit(target)):
                return cycle
        stack.pop()
        state[module] = 2
        return None

    for module in sorted(graph):
        if module not in state and (cycle := visit(module)):
            return cycle
    return None


def relative(path: Path) -> Path:
    return path.relative_to(ROOT) if path.is_relative_to(ROOT) else path


def validate(root: Path = SOURCE_ROOT) -> list[str]:
    failures: list[str] = []
    for path in sorted(root.rglob("*.rs")):
        if path.relative_to(root).as_posix() in CONTAINERS or path.name == "tests.rs":
            continue
        source = module_path(root, path)
        if unit(source) is None:
            failures.append(f"{source}: unmapped module; assign it a crate in MODULES")
    graph: dict[str, set[str]] = {}
    for edge in edges(root):
        source_unit = unit(edge.source)
        if source_unit is None:
            continue
        target_unit = unit(edge.target)
        if target_unit is None:
            failures.append(
                f"{relative(edge.path)}:{edge.line}: unmapped dependency crate::{edge.target}"
            )
            continue
        if target_unit != source_unit:
            graph.setdefault(source_unit, set()).add(target_unit)
        source_crate = MODULES[source_unit]
        target_crate = MODULES[target_unit]
        if target_crate != source_crate and target_crate not in ALLOWED[source_crate]:
            failures.append(
                f"{relative(edge.path)}:{edge.line}: {source_crate} module {source_unit} "
                f"must not depend on {target_crate} module {target_unit}"
            )
    if cycle := find_cycle(graph):
        failures.append(f"module dependency cycle: {' -> '.join(cycle)}")
    return failures


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
    failures = validate(args.root)
    for failure in failures:
        print(f"FAIL: {failure}", file=sys.stderr)
    if failures:
        print(f"{len(failures)} crate-map violation(s) found", file=sys.stderr)
        return 1
    print(
        "OK: every module belongs to a plugin-family crate, depends only on its own "
        "crate or an allowed one, and the module graph is acyclic"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
