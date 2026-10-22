"""Check that the shipped configuration matches its reference table.

`config/dna.toml` is the configuration DNA ships; the table under the
"Schema version N" heading of `docs/reference/configuration.md` documents it
key by key. This check fails when they drift apart: a key in one and not the
other, a different value, or a schema version that differs between the
heading, the table's `schema_version` row, the TOML, and the version the
build accepts (`SCHEMA_VERSION` in `src/config/types.rs`).
"""

from __future__ import annotations

import argparse
import re
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
CONFIG = Path("config/dna.toml")
REFERENCE = Path("docs/reference/configuration.md")
TYPES = Path("src/config/types.rs")
HEADING = re.compile(r"^## Schema version (\d+)\s*$")
SOURCE_VERSION = re.compile(r"^const SCHEMA_VERSION: u32 = (\d+);", re.MULTILINE)
CODE = re.compile(r"^`([^`]*)`$")
ROOT_SECTION = "root"

type Value = object
type Key = tuple[str, str]


def parse_value(text: str) -> Value:
    """A table value as TOML reads it; a bare path stays a string."""
    try:
        return tomllib.loads(f"value = {text}")["value"]
    except tomllib.TOMLDecodeError:
        return text


def same(left: Value, right: Value) -> bool:
    """Equal values of the same TOML type (`1` is not `1.0`)."""
    if isinstance(left, list) and isinstance(right, list):
        return len(left) == len(right) and all(map(same, left, right))
    return type(left) is type(right) and left == right


def reference_table(text: str) -> tuple[int | None, dict[Key, Value], list[str]]:
    """The heading's schema version and the table's keys with their values.

    A row whose section cell is empty belongs to the section above it.
    """
    errors: list[str] = []
    version: int | None = None
    rows: dict[Key, Value] = {}
    lines = text.splitlines()
    start = next(
        (index for index, line in enumerate(lines) if HEADING.match(line)), None
    )
    if start is None:
        return None, rows, ["reference: no '## Schema version N' heading"]
    heading = HEADING.match(lines[start])
    assert heading is not None
    version = int(heading.group(1))
    section = ""
    for line in lines[start + 1 :]:
        if line.startswith("## "):
            break
        if not line.startswith("|"):
            if rows:
                break
            continue
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        if len(cells) < 3 or cells[1] in ("Key", "") or set(cells[1]) <= set("-:"):
            continue
        section_cell, key_cell, value_cell = cells[:3]
        if section_cell:
            section_match = CODE.match(section_cell)
            section = section_match.group(1) if section_match else section_cell
        if not section:
            errors.append(f"reference: row {key_cell} has no section")
            continue
        key_match = CODE.match(key_cell)
        value_match = CODE.match(value_cell)
        if key_match is None or value_match is None:
            errors.append(f"reference: row {key_cell} is not `key` | `value`")
            continue
        key = (section, key_match.group(1))
        if key in rows:
            errors.append(f"reference: {section}.{key[1]} is listed twice")
        rows[key] = parse_value(value_match.group(1))
    return version, rows, errors


def config_keys(config: dict[str, object]) -> dict[Key, Value]:
    """The shipped configuration's keys, root keys under `root`."""
    keys: dict[Key, Value] = {}
    for name, value in config.items():
        if isinstance(value, dict):
            for key, inner in value.items():
                keys[(name, key)] = inner
        else:
            keys[(ROOT_SECTION, name)] = value
    return keys


def violations(root: Path) -> list[str]:
    config = tomllib.loads((root / CONFIG).read_text(encoding="utf-8"))
    version, rows, errors = reference_table(
        (root / REFERENCE).read_text(encoding="utf-8")
    )
    shipped = config_keys(config)
    for key in sorted(shipped.keys() - rows.keys()):
        errors.append(f"{CONFIG}: {key[0]}.{key[1]} is not in {REFERENCE}")
    for key in sorted(rows.keys() - shipped.keys()):
        errors.append(f"{REFERENCE}: {key[0]}.{key[1]} is not in {CONFIG}")
    for key in sorted(rows.keys() & shipped.keys()):
        if not same(rows[key], shipped[key]):
            errors.append(
                f"{REFERENCE}: {key[0]}.{key[1]} is {rows[key]!r}, "
                f"{CONFIG} has {shipped[key]!r}"
            )
    source = SOURCE_VERSION.search((root / TYPES).read_text(encoding="utf-8"))
    versions = {
        f"{CONFIG} schema_version": config.get("schema_version"),
        f"{REFERENCE} heading": version,
        f"{TYPES} SCHEMA_VERSION": int(source.group(1)) if source else None,
    }
    if len(set(versions.values())) != 1:
        stated = ", ".join(f"{name} {value}" for name, value in versions.items())
        errors.append(f"schema versions differ: {stated}")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    arguments = parser.parse_args()
    errors = violations(arguments.root)
    for error in errors:
        print(error, file=sys.stderr)
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
