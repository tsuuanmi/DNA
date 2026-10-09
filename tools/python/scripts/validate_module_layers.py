"""Enforce the crate-ready module layering of first-party Rust source (ADR-0064).

Every top-level module under ``src/`` belongs to one layer. Production code may
depend only on modules in the same or a lower layer, and the module graph must
be acyclic, so each layer can later become a crate without redesign.

Modality-neutral modules (the core caller, the evidence contract, and
post-calling representation, ADR-0069) must not depend on the Sanger modality:
neither on its modules nor on the Sanger children of ``model``.
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SOURCE_ROOT = ROOT / "src"

# Layer of every top-level module: 0 core, 1 target data, 2 science,
# 3 adapters/capabilities, 4 delivery. Keep in sync with ADR-0064.
LAYERS: dict[str, int] = {
    "error": 0,
    "checksum": 0,
    "model": 0,
    "locus": 0,
    "variant": 0,
    "read_evidence": 0,
    "config": 1,
    "reference": 1,
    "profile": 1,
    "basecalling": 2,
    "signal_processing": 2,
    "callability": 2,
    "quality_control": 2,
    "read_processing": 2,
    "alignment": 2,
    "variant_calling": 2,
    "sample": 2,
    "variant_representation": 2,
    "variant_normalization": 2,
    "variant_nomenclature": 2,
    "input": 3,
    "variant_analysis": 3,
    "report": 4,
    "operation_log": 4,
    "pipeline": 4,
    "cli": 4,
}

# Modules that must stay independent of any sequencing modality (ADR-0069).
NEUTRAL: frozenset[str] = frozenset(
    {
        "read_evidence",
        "variant",
        "alignment",
        "variant_calling",
        "variant_representation",
        "variant_normalization",
        "variant_nomenclature",
    }
)
# Modules that implement the Sanger modality.
SANGER: frozenset[str] = frozenset(
    {
        "basecalling",
        "signal_processing",
        "callability",
        "quality_control",
        "read_processing",
        "input",
        "locus",
        "sample",
    }
)
# Children of ``model`` that neutral modules may use; every other child is
# Sanger-specific or a delivery contract.
NEUTRAL_MODEL_CHILDREN: frozenset[str] = frozenset(
    {"alignment", "coordinate", "nucleotide", "reference", "variant"}
)

# Crate roots compose every module and are not part of any layer.
CRATE_ROOTS = {"lib.rs", "main.rs"}

TEST_MODULE = re.compile(
    r"#\[cfg\(test\)\]\s*(?:pub(?:\([^)]*\))?\s+)?mod\s+\w+\s*\{", re.MULTILINE
)
LINE_COMMENT = re.compile(r"//.*$", re.MULTILINE)
CRATE_PATH = re.compile(r"\bcrate::([a-z_][a-z0-9_]*)")
CRATE_GROUP = re.compile(r"\bcrate::\{([^}]*)\}")
GROUP_MEMBER = re.compile(r"(?:^|,)\s*([a-z_][a-z0-9_]*)")
MODEL_PATH = re.compile(r"\bcrate::model::([a-z_][a-z0-9_]*)")
MODEL_GROUP = re.compile(r"\bcrate::model::\{([^}]*)\}")


@dataclass(frozen=True)
class Edge:
    source: str
    target: str
    path: Path
    line: int


def module_of(root: Path, path: Path) -> str | None:
    relative = path.relative_to(root)
    if len(relative.parts) == 1:
        return None if relative.name in CRATE_ROOTS else relative.stem
    return relative.parts[0]


def production_text(path: Path) -> str:
    """Source without inline test modules (which run to the end of the file) or comments."""
    text = path.read_text(encoding="utf-8")
    test_module = TEST_MODULE.search(text)
    if test_module:
        text = text[: test_module.start()]
    return LINE_COMMENT.sub("", text)


def edges(root: Path) -> list[Edge]:
    found: list[Edge] = []
    for path in sorted(root.rglob("*.rs")):
        source = module_of(root, path)
        if source is None or path.name == "tests.rs":
            continue
        text = production_text(path)
        targets: list[tuple[str, int]] = []
        for match in CRATE_GROUP.finditer(text):
            line = text.count("\n", 0, match.start()) + 1
            targets.extend(
                (member.group(1), line)
                for member in GROUP_MEMBER.finditer(match.group(1))
            )
        for match in CRATE_PATH.finditer(text):
            targets.append((match.group(1), text.count("\n", 0, match.start()) + 1))
        for target, line in targets:
            if target != source:
                found.append(Edge(source, target, path, line))
    return found


def model_children(root: Path) -> list[Edge]:
    """References from modules to children of ``model``, as ``model::child`` edges."""
    found: list[Edge] = []
    for path in sorted(root.rglob("*.rs")):
        source = module_of(root, path)
        if source is None or path.name == "tests.rs":
            continue
        text = production_text(path)
        for match in MODEL_GROUP.finditer(text):
            line = text.count("\n", 0, match.start()) + 1
            found.extend(
                Edge(source, member.group(1), path, line)
                for member in GROUP_MEMBER.finditer(match.group(1))
            )
        for match in MODEL_PATH.finditer(text):
            found.append(
                Edge(
                    source, match.group(1), path, text.count("\n", 0, match.start()) + 1
                )
            )
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
    graph: dict[str, set[str]] = {}
    modules = {
        module
        for path in root.rglob("*.rs")
        if (module := module_of(root, path)) is not None
    }
    failures.extend(
        f"{module}: unmapped module; assign it a layer in LAYERS (ADR-0064)"
        for module in sorted(modules - LAYERS.keys())
    )
    for edge in edges(root):
        graph.setdefault(edge.source, set()).add(edge.target)
        source_layer = LAYERS.get(edge.source)
        target_layer = LAYERS.get(edge.target)
        if target_layer is None:
            if edge.target in modules:
                continue
            failures.append(
                f"{relative(edge.path)}:{edge.line}: unknown module crate::{edge.target}"
            )
        elif source_layer is not None and target_layer > source_layer:
            failures.append(
                f"{relative(edge.path)}:{edge.line}: layer {source_layer} module "
                f"{edge.source} must not depend on layer {target_layer} module {edge.target}"
            )
    for edge in edges(root):
        if edge.source in NEUTRAL and edge.target in SANGER:
            failures.append(
                f"{relative(edge.path)}:{edge.line}: modality-neutral module {edge.source} "
                f"must not depend on Sanger module {edge.target} (ADR-0069)"
            )
    for edge in model_children(root):
        if edge.source in NEUTRAL and edge.target not in NEUTRAL_MODEL_CHILDREN:
            failures.append(
                f"{relative(edge.path)}:{edge.line}: modality-neutral module {edge.source} "
                f"must not depend on model::{edge.target} (ADR-0069)"
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
        print(f"{len(failures)} module-layering violation(s) found", file=sys.stderr)
        return 1
    print(
        "OK: Rust modules depend only on the same or lower layers, without cycles, "
        "and modality-neutral modules stay independent of Sanger"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
