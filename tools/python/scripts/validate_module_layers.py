"""Enforce the crate map of the workspace (ADR-0069, ADR-0070, ADR-0075, ADR-0079).

The workspace has one facade crate at the root (``src/``, crate ``dna``), one
member crate per plugin family or core step under ``crates/``:

- ``dna-kernel``: shared contracts;
- ``dna-align``, ``dna-call``, ``dna-sample``, ``dna-consensus``: the core's
  steps, a chain in which each depends only on the kernel and the steps
  before it (``dna-align`` <- ``dna-call`` <- ``dna-sample``,
  ``dna-consensus``);
- ``dna-sanger``: the Sanger modality;
- ``dna-trace``: trace signal engineering, a library only ``dna-sanger``
  depends on (ADR-0084); the harness also uses its ``simulate`` feature, the
  forward-model renderer of the simulated bench (PROP-0015);
- ``dna-normalize``: post-calling representation, normalization, and
  nomenclature;
- ``dna-empop``: post-calling notation conformance;

and the evaluation harness: its library ``dna-eval-lib`` and the thin
``dna-eval`` binary over it, which records what it was built from.

Inside the facade, the modules that stay modality-neutral (ADR-0069, the
``NEUTRAL_FACADE`` paths) may not name a modality crate (``dna-sanger``,
``dna-trace``): the contracts' records, the readers of neutral documents,
the neutral reports and pipelines, and the review's rules but its Sanger
part. A modality value a contract reports has its own facade record,
projected in ``model::sanger_read``.

The modality and post-calling crates depend only on the kernel, the Sanger
modality on ``dna-trace`` as well, and only the facade composes the plugins. The harness may depend on the facade and the
crates whose records it reads, and no crate but its binary may depend on the
harness. Both source paths (``dna::…``, ``dna_<crate>::…``) and manifest
dependencies must follow that graph. Within each crate the module graph must
be acyclic.

The harness caches what a stage of DNA made, keyed by the stage's code
identity: the git object ids of the sources ``crates/dna-eval/stages.txt``
lists for it. Each list must be closed, or a cache would serve output that
code outside the list changed: every module a listed source uses (through
``crate::``, ``super::``, ``self::``, a child module, a ``use`` binding, or a
re-export), every module file above a listed source, every workspace crate a
listed source or crate depends on, and the build's own inputs, is listed. A
file listed for some items only (``path: a, b``) must define or re-export
them, and what those items' definitions and ``impl`` blocks use (directly,
through the file's ``use`` bindings, or through the file's own top-level
items they name, transitively) is checked like a whole file's uses; what its
other items use is not followed. A listed file may not glob-import its own
crate (``use super::*``), whose names the check cannot follow. A path marked ``!path``
is used by a listed source but is not on the stage's path; it must be used.
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
    "align": frozenset({"kernel"}),
    "call": frozenset({"kernel", "align"}),
    "sample": frozenset({"kernel", "align", "call"}),
    "consensus": frozenset({"kernel", "align", "call"}),
    "trace": frozenset({"kernel"}),
    "sanger": frozenset({"kernel", "trace"}),
    "normalize": frozenset({"kernel"}),
    "empop": frozenset({"kernel"}),
    "dna": frozenset(
        {
            "kernel",
            "align",
            "call",
            "sample",
            "consensus",
            "sanger",
            "normalize",
            "empop",
        }
    ),
    "eval": frozenset({"eval-lib"}),
    "eval-lib": frozenset(
        {
            "dna",
            "kernel",
            "align",
            "sample",
            "consensus",
            "sanger",
            "trace",
            "normalize",
            "empop",
        }
    ),
}
# Modality crates the neutral facade modules may not name.
MODALITY_CRATES = frozenset({"sanger", "trace"})
# Facade source files (globs relative to src/) that stay modality-neutral.
NEUTRAL_FACADE = (
    "model/result.rs",
    "model/*_result.rs",
    "model/run_context.rs",
    "input/consensus.rs",
    "input/run_context.rs",
    "input/sequence.rs",
    "input/variants.rs",
    "pipeline/call.rs",
    "pipeline/notation.rs",
    "pipeline/path.rs",
    "pipeline/represent.rs",
    "pipeline/sample_call.rs",
    "report/aligned_calls.rs",
    "report/atomic.rs",
    "report/call.rs",
    "report/consensus.rs",
    "report/json.rs",
    "report/notation.rs",
    "report/variant.rs",
    "report/review/*.rs",
)
# The review's Sanger part: its run signal input (mod.rs) and second peaks.
NEUTRAL_EXCEPTIONS = frozenset({"report/review/mod.rs", "report/review/second.rs"})
# Source path prefix of each crate.
PATH_PREFIXES = {
    ("dna" if name == "dna" else f"dna_{name.replace('-', '_')}"): name
    for name in ALLOWED
}
# Modules whose children are tracked as separate nodes of the module graph.
SPLIT_MODULES = {"model", "input"}

TEST_MODULE = re.compile(
    r"#\[cfg\(test\)\]\s*(?:pub(?:\([^)]*\))?\s+)?mod\s+\w+\s*\{", re.MULTILINE
)
LINE_COMMENT = re.compile(r"//.*$", re.MULTILINE)
PATH_START = re.compile(
    r"\b(crate|" + "|".join(sorted(PATH_PREFIXES)) + r")::",
)
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
        for name in document.get(table, {}):
            if name == "dna":
                names.append(name)
            elif name.startswith("dna-"):
                names.append(name.removeprefix("dna-"))
    return names


def package(name: str) -> str:
    """Package name of a crate-map name."""
    return name if name == "dna" else f"dna-{name}"


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


# SRS-COMPAT-003
# SRS-CONS-002: the consensus crate may not depend on a modality crate.
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
                    f"{package(dependency)}"
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
    failures.extend(neutrality_violations(root))
    return failures


def neutral_files(source: Path) -> list[Path]:
    """The facade's modality-neutral source files."""
    found = {
        path
        for pattern in NEUTRAL_FACADE
        for path in source.glob(pattern)
        if path.relative_to(source).as_posix() not in NEUTRAL_EXCEPTIONS
    }
    return sorted(found)


# ADR-0069: the facade's neutral modules name no modality crate.
def neutrality_violations(root: Path) -> list[str]:
    source = root / "src"
    failures: list[str] = []
    for path in neutral_files(source):
        text = production_text(path)
        for match in PATH_START.finditer(text):
            prefix = match.group(1)
            if prefix != "crate" and PATH_PREFIXES[prefix] in MODALITY_CRATES:
                line = text.count("\n", 0, match.start()) + 1
                failures.append(
                    f"{relative(root, path)}:{line}: a modality-neutral facade module "
                    f"must not name {package(PATH_PREFIXES[prefix])}"
                )
    return failures


# The sources of each cached stage of the evaluation harness.
STAGES = Path("crates/dna-eval/stages.txt")
# What every build reads besides the sources.
BUILD_INPUTS = ("Cargo.lock", "Cargo.toml", "rust-toolchain.toml")
MOD_DECLARATION = re.compile(r"\bmod\s+([a-z_][a-z0-9_]*)\s*;")
USE = re.compile(r"\buse\s+([^;]*);")
PATH_HEAD = re.compile(r"(?<!::)\b(crate|super|self|[a-z_][a-z0-9_]*)\s*::")
WORD = re.compile(r"[A-Za-z_][A-Za-z0-9_]*|::|\{|\}|,|\*")
# A named item a file defines.
ITEM_DEFINITION = re.compile(
    r"(?<!')\b(?:fn|struct|enum|union|const|static|type|trait)\s+([A-Za-z_][A-Za-z0-9_]*)"
)
WORD_NAME = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")
# A glob import, whose names the closure cannot follow: its path.
GLOB_USE = re.compile(r"\buse\s+([^;]*\*)\s*\}?\s*;")
# How deep re-exports are followed.
MAX_REEXPORTS = 8


@dataclass
class Stage:
    name: str
    line: int
    # Listed path -> the items it is listed for, or None for all of it.
    listed: dict[str, tuple[str, ...] | None]
    excluded: list[str]


def parse_stages(text: str) -> tuple[list[Stage], list[str]]:
    """The stages of a stage file, and its syntax errors."""
    stages: list[Stage] = []
    failures: list[str] = []
    for number, raw in enumerate(text.splitlines(), start=1):
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            stages.append(Stage(line[1:-1].strip(), number, {}, []))
        elif not stages:
            failures.append(f"{STAGES}:{number}: a path before any [stage]")
        elif line.startswith("!"):
            stages[-1].excluded.append(line[1:].strip())
        else:
            path, colon, items = line.partition(":")
            names = tuple(item.strip() for item in items.split(",") if item.strip())
            if colon and not names:
                failures.append(f"{STAGES}:{number}: {path.strip()}: lists no item")
            stages[-1].listed[path.strip()] = names if colon else None
    return stages, failures


@dataclass(frozen=True)
class Tree:
    """The module tree of one crate's source directory."""

    source: Path

    def file(self, module: list[str]) -> Path | None:
        """The file of a module, if it exists."""
        if not module:
            candidates = [self.source / "lib.rs", self.source / "main.rs"]
        else:
            base = self.source.joinpath(*module)
            candidates = [base.with_suffix(".rs"), base / "mod.rs"]
        return next((path for path in candidates if path.is_file()), None)

    def module(self, path: Path) -> list[str]:
        """The module a source file holds."""
        parts = list(path.relative_to(self.source).with_suffix("").parts)
        if parts in (["lib"], ["main"]):
            return []
        if parts[-1] == "mod":
            parts.pop()
        return parts


def use_bindings(text: str) -> list[tuple[str, list[str]]]:
    """Every name a ``use`` binds, with the path it binds."""
    bindings: list[tuple[str, list[str]]] = []
    for statement in USE.finditer(text):
        for path in expand(WORD.findall(statement.group(1))):
            segments = path.split("::")
            if "as" in segments:
                index = segments.index("as")
                if index + 1 < len(segments):
                    bindings.append((segments[index + 1], segments[:index]))
            elif segments:
                bindings.append((segments[-1], segments))
    return bindings


def resolve(
    tree: Tree, module: list[str], rest: list[str], depth: int = 0
) -> Path | None:
    """The file that holds ``rest`` of ``module``, through child modules and
    re-exports."""
    while rest and tree.file([*module, rest[0]]):
        module, rest = [*module, rest[0]], rest[1:]
    path = tree.file(module)
    if path is None or not rest or depth >= MAX_REEXPORTS:
        return path
    text = production_text(path)
    children = set(MOD_DECLARATION.findall(text))
    for name, bound in use_bindings(text):
        if name != rest[0] or not bound:
            continue
        target = relative_path(module, children, bound)
        if target is not None:
            return resolve(tree, target[0], target[1] + rest[1:], depth + 1)
    return path


def relative_path(
    module: list[str], children: set[str], segments: list[str]
) -> tuple[list[str], list[str]] | None:
    """The (module, rest) a crate-relative path starts from, or None for a
    path into another crate."""
    head, rest = segments[0], segments[1:]
    if head == "crate":
        return [], rest
    if head == "self":
        return module, rest
    if head == "super":
        base = module[:-1]
        while rest and rest[0] == "super":
            base, rest = base[:-1], rest[1:]
        return base, rest
    if head in children:
        return module, segments
    return None


def item_spans(text: str, item: str) -> list[str]:
    """The text of every definition, ``impl`` block, and ``use`` of ``item``:
    from its keyword to its matching closing brace, or to its ``;``."""
    spans: list[str] = []
    starts = re.finditer(
        rf"\b(?:fn|struct|enum|mod|const|static|type|trait)\s+{item}\b"
        rf"|\bimpl\b(?:\s*<[^{{;]*?>)?\s+(?:[\w:<>, ]+\s+for\s+)?{item}\b"
        rf"|\buse\b[^;]*\b{item}\b",
        text,
    )
    for start in starts:
        position = start.start()
        brace = text.find("{", start.end())
        end = text.find(";", start.end())
        if start.group(0).startswith("use") or brace < 0 or 0 <= end < brace:
            spans.append(text[position : (end + 1 if end >= 0 else len(text))])
            continue
        depth = 0
        for index in range(brace, len(text)):
            if text[index] == "{":
                depth += 1
            elif text[index] == "}":
                depth -= 1
                if depth == 0:
                    spans.append(text[position : index + 1])
                    break
        else:
            spans.append(text[position:])
    return spans


def top_level_items(text: str) -> set[str]:
    """The names a file defines at its top level, outside every block (a
    method of an ``impl`` is not one)."""
    depths: list[int] = []
    depth = 0
    for character in text:
        depths.append(depth)
        if character == "{":
            depth += 1
        elif character == "}":
            depth -= 1
    return {
        match.group(1)
        for match in ITEM_DEFINITION.finditer(text)
        if depths[match.start()] == 0 and match.group(1) != "fn"
    }


def item_closure(text: str, items: tuple[str, ...]) -> str:
    """The spans of ``items`` and, transitively, of every top-level item of
    the same file their spans name (a private helper a listed function
    calls). A method call (``.name(``), a path segment (``a::name``), a
    child module's path (``name::``), and a ``use`` statement name no item
    of the file."""
    defined = top_level_items(text)
    children = set(MOD_DECLARATION.findall(text))
    names = set(items)
    pending = list(items)
    spans: list[str] = []
    while pending:
        for span in item_spans(text, pending.pop()):
            spans.append(span)
            if span.startswith("use"):
                continue
            for match in WORD_NAME.finditer(span):
                word = match.group(0)
                before = span[: match.start()].rstrip()
                after = span[match.end() :].lstrip()
                if before.endswith((".", "::")) or (
                    word in children and after.startswith("::")
                ):
                    continue
                if word in defined and word not in names:
                    names.add(word)
                    pending.append(word)
    return "\n".join(spans)


def stage_uses(
    tree: Tree, path: Path, items: tuple[str, ...] | None = None
) -> tuple[set[Path], set[str]]:
    """The files of its own crate and the crates a source file uses, or,
    with ``items``, what those items' definitions use (through the file's
    ``use`` bindings too)."""
    text = production_text(path)
    module = tree.module(path)
    children = set(MOD_DECLARATION.findall(text))
    bindings = dict(use_bindings(text))
    scanned = text if items is None else item_closure(text, items)
    files: set[Path] = set()
    crates_used: set[str] = set()
    if items is not None:
        words = set(re.findall(r"[A-Za-z_][A-Za-z0-9_]*", scanned))
        for name, bound in bindings.items():
            if name not in words or not bound:
                continue
            if bound[0] in PATH_PREFIXES:
                crates_used.add(PATH_PREFIXES[bound[0]])
                continue
            start = relative_path(module, children, bound)
            if start is not None and (found := resolve(tree, *start)) is not None:
                files.add(found)
    for match in PATH_HEAD.finditer(scanned):
        head = match.group(1)
        for full in expand(path_tokens(scanned, match.end())):
            segments = [head, *full.split("::")]
            if head in PATH_PREFIXES:
                crates_used.add(PATH_PREFIXES[head])
                continue
            start = relative_path(module, children, segments)
            if start is None and head in bindings:
                start = relative_path(
                    module, children, [*bindings[head], *segments[1:]]
                )
            if start is None:
                continue
            if (found := resolve(tree, *start)) is not None:
                files.add(found)
    files.discard(path)
    return files, crates_used


def crate_directory(name: str) -> str:
    """The directory of a crate-map name."""
    return "src" if name == "dna" else f"crates/{package(name)}"


def validate_stages(root: Path = ROOT) -> list[str]:
    """Every failure of the stage lists' closure."""
    if not (root / STAGES).is_file():
        return []
    stages, failures = parse_stages((root / STAGES).read_text(encoding="utf-8"))
    for stage in stages:
        failures.extend(check_stage(root, stage))
    return failures


def check_stage(root: Path, stage: Stage) -> list[str]:
    failures: list[str] = []
    where = f"{STAGES}:{stage.line}: [{stage.name}]"
    listed = {root / path: items for path, items in stage.listed.items()}
    excluded = [root / path for path in stage.excluded]
    used: set[Path] = set()

    def covered(path: Path) -> bool:
        for exclusion in excluded:
            if path == exclusion or exclusion in path.parents:
                used.add(exclusion)
                return True
        return any(path == entry or entry in path.parents for entry in listed)

    def show(path: Path) -> str:
        return relative(root, path).as_posix()

    for path in [*listed, *excluded]:
        if not path.exists():
            failures.append(f"{where} {show(path)} does not exist")
    failures.extend(
        f"{where} lists no {name}" for name in BUILD_INPUTS if root / name not in listed
    )
    sources: list[tuple[Path, tuple[str, ...] | None]] = []
    for path, items in listed.items():
        if path.parent == root / "crates" and path.is_dir():
            for dependency in manifest_dependencies(path / "Cargo.toml"):
                if dependency in ALLOWED and not covered(
                    root / crate_directory(dependency)
                ):
                    failures.append(
                        f"{where} {show(path)} depends on {package(dependency)}, "
                        "which is not listed"
                    )
        elif path.is_dir():
            sources.extend((source, None) for source in sorted(path.rglob("*.rs")))
        elif path.suffix == ".rs":
            sources.append((path, items))
    for path, items in sources:
        source = next(
            (
                directory
                for directory in [
                    root / "src",
                    *sorted((root / "crates").glob("*/src")),
                ]
                if directory in path.parents
            ),
            None,
        )
        if source is None:
            continue
        tree = Tree(source)
        module = tree.module(path)
        for depth in range(len(module)):
            parent = tree.file(module[:depth])
            if parent is not None and not covered(parent):
                failures.append(
                    f"{where} {show(path)} is in {show(parent)}, which is not listed"
                )
        for glob in GLOB_USE.finditer(production_text(path)):
            # A glob of another crate's module is that crate's, which the
            # crate check covers; a glob of this crate hides what it uses.
            head = WORD_NAME.match(glob.group(1).strip())
            if head is None or head.group(0) not in PATH_PREFIXES:
                failures.append(
                    f"{where} {show(path)} has a glob import (`use …::*`), whose "
                    "names the check cannot follow; import the names"
                )
        if items is not None:
            text = production_text(path)
            failures.extend(
                f"{where} {show(path)} has no item {item}"
                for item in items
                if not item_spans(text, item)
            )
        files, crates_used = stage_uses(tree, path, items)
        failures.extend(
            f"{where} {show(path)} uses {show(file)}, which is not listed"
            for file in sorted(files)
            if not covered(file)
        )
        own = (
            "dna" if source == root / "src" else source.parent.name.removeprefix("dna-")
        )
        failures.extend(
            f"{where} {show(path)} uses {package(name)}, which is not listed"
            for name in sorted(crates_used)
            if name != own and not covered(root / crate_directory(name))
        )
    failures.extend(
        f"{where} !{show(path)} is used by no listed source"
        for path in excluded
        if path not in used
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
    failures = validate(args.root) + validate_stages(args.root)
    for failure in failures:
        print(f"FAIL: {failure}", file=sys.stderr)
    if failures:
        print(f"{len(failures)} crate-map violation(s) found", file=sys.stderr)
        return 1
    print(
        "OK: the core steps form a chain on the kernel, the other plugin crates "
        "depend only on the kernel, the facade on the plugin crates, and nothing "
        "on the evaluation harness but its binary, in source and manifests, each "
        "crate's module graph is acyclic, the facade's neutral modules name no "
        "modality crate, and the harness's stage lists are closed"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
