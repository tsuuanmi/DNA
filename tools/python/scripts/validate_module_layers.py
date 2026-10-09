"""Enforce the plugin-family crate map of the workspace (ADR-0069, ADR-0070).

The workspace has one facade crate at the root (``src/``, crate ``dna``) and
one member crate per plugin family under ``crates/``:

- ``dna-kernel``: shared contracts;
- ``dna-core``: the core caller;
- ``dna-sanger``: the Sanger modality;
- ``dna-post``: post-calling plugins.

The plugin crates depend only on the kernel, and only the facade composes
them. Both source paths (``dna_<crate>::…``) and manifest dependencies must
follow that graph. Within each crate the module graph must be acyclic.
"""

from __future__ import annotations

import argparse
import re
import sys
import tomllib
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]

# Crates each crate may depend on, besides itself.
ALLOWED: dict[str, frozenset[str]] = {
    "kernel": frozenset(),
    "core": frozenset({"kernel"}),
    "sanger": frozenset({"kernel"}),
    "post": frozenset({"kernel"}),
    "dna": frozenset({"kernel", "core", "sanger", "post"}),
}
# Source path prefix of each member crate.
PATH_PREFIXES = {f"dna_{name}": name for name in ALLOWED if name != "dna"}
# Modules whose children are tracked as separate nodes of the module graph.
SPLIT_MODULES = {"model", "input"}

TEST_MODULE = re.compile(
    r"#\[cfg\(test\)\]\s*(?:pub(?:\([^)]*\))?\s+)?mod\s+\w+\s*\{", re.MULTILINE
)
LINE_COMMENT = re.compile(r"//.*$", re.MULTILINE)
PATH_START = re.compile(r"\b(crate|dna_kernel|dna_core|dna_sanger|dna_post)::")
TOKEN = re.compile(r"\s*([A-Za-z_][A-Za-z0-9_]*|::|\{|\}|,|\*)")


@dataclass(frozen=True)
class Edge:
    source: str
    target: str
    path: Path
    line: int


@dataclass(frozen=True)
class Crate:
    name: str
    source: Path
    manifest: Path


def crates(root: Path) -> tuple[list[Crate], list[str]]:
    """The facade and every member crate, plus unknown member directories."""
    found = [Crate("dna", root / "src", root / "Cargo.toml")]
    unknown: list[str] = []
    for directory in sorted((root / "crates").glob("*")):
        if not directory.is_dir():
            continue
        name = directory.name.removeprefix("dna-")
        if not directory.name.startswith("dna-") or name not in ALLOWED:
            unknown.append(directory.name)
            continue
        found.append(Crate(name, directory / "src", directory / "Cargo.toml"))
    return found, unknown


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


def unit(path: str) -> str:
    """Module-graph node of a crate-relative path."""
    segments = path.split("::")
    if segments[0] in SPLIT_MODULES and len(segments) > 1:
        return "::".join(segments[:2])
    return segments[0]


def module_of(source: Path, path: Path) -> str | None:
    """Crate-relative module of a source file, or None for the crate root."""
    parts = list(path.relative_to(source).with_suffix("").parts)
    if parts[-1] == "mod":
        parts.pop()
    if parts in (["lib"], ["main"]):
        return None
    return unit("::".join(parts)) if parts else None


def edges(crate: Crate) -> list[tuple[str, Edge]]:
    """Every dependency of the crate's production source, with its target crate."""
    found: list[tuple[str, Edge]] = []
    for path in sorted(crate.source.rglob("*.rs")):
        if path.name == "tests.rs":
            continue
        source = module_of(crate.source, path) or "crate root"
        text = production_text(path)
        for match in PATH_START.finditer(text):
            prefix = match.group(1)
            target_crate = crate.name if prefix == "crate" else PATH_PREFIXES[prefix]
            line = text.count("\n", 0, match.start()) + 1
            for target in expand(path_tokens(text, match.end())):
                found.append((target_crate, Edge(source, target, path, line)))
    return found


def manifest_dependencies(manifest: Path) -> list[str]:
    """Workspace crates a manifest depends on, including dev-dependencies."""
    if not manifest.is_file():
        return []
    document = tomllib.loads(manifest.read_text(encoding="utf-8"))
    names: list[str] = []
    for table in ("dependencies", "dev-dependencies"):
        names.extend(
            name.removeprefix("dna-")
            for name in document.get(table, {})
            if name.startswith("dna-")
        )
    return names


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


def relative(root: Path, path: Path) -> Path:
    return path.relative_to(root) if path.is_relative_to(root) else path


def validate(root: Path = ROOT) -> list[str]:
    failures: list[str] = []
    members, unknown = crates(root)
    failures.extend(
        f"crates/{name}: unknown crate; add it to the crate map (ADR-0070)"
        for name in unknown
    )
    for crate in members:
        for dependency in manifest_dependencies(crate.manifest):
            if dependency not in ALLOWED[crate.name]:
                failures.append(
                    f"{relative(root, crate.manifest)}: {crate.name} must not depend on "
                    f"dna-{dependency}"
                )
        graph: dict[str, set[str]] = {}
        for target_crate, edge in edges(crate):
            if target_crate != crate.name:
                if target_crate not in ALLOWED[crate.name]:
                    failures.append(
                        f"{relative(root, edge.path)}:{edge.line}: {crate.name} module "
                        f"{edge.source} must not depend on {target_crate} module "
                        f"{unit(edge.target)}"
                    )
                continue
            target = unit(edge.target)
            if edge.source != "crate root" and target != edge.source:
                graph.setdefault(edge.source, set()).add(target)
        if cycle := find_cycle(graph):
            failures.append(
                f"{crate.name}: module dependency cycle: {' -> '.join(cycle)}"
            )
    return failures


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=ROOT,
        help="workspace root to validate (default: this repository)",
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
        "OK: crates depend only on the kernel or, for the facade, on every plugin "
        "crate, in source and manifests, and each crate's module graph is acyclic"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
