import re
from pathlib import Path

from docs_graph import load_profile

ROOT = Path(__file__).resolve().parents[3]
DOCS = ROOT / "docs"
PROFILE = load_profile(ROOT)
# Source trees and forbidden names come from the repository profile
# (docs/governance/documentation.md); this script enforces DOC-STR-* and DOC-LNK-001.
SOURCES = [
    path for pattern in PROFILE.source_roots for path in sorted(ROOT.glob(pattern))
]
FORBIDDEN_NAMES = set(PROFILE.forbidden_names)
MARKDOWN_LINK = re.compile(r"!?\[[^\]]*\]\(([^)]+)\)")
# PROP-0011: each crate's source README names the contract its step writes.
OUTPUT_CONTRACT = "Output contract:"
SCHEMA_FILE = re.compile(r"[A-Za-z0-9_.-]+\.schema\.json")


def output_contract_violations(root: Path) -> list[str]:
    """Check the "Output contract:" paragraphs of the source READMEs.

    Every crate's source README has one; it names at least one schema unless it
    states "none"; every schema it names exists under docs/reference/schemas.
    """
    errors: list[str] = []
    schemas = root / "docs" / "reference" / "schemas"
    crate_readmes = sorted(root.glob("crates/*/src/README.md"))
    readmes = sorted(
        {
            *(root / "src").rglob("README.md"),
            *root.glob("crates/*/src/**/README.md"),
        }
    )
    for readme in readmes:
        relative = readme.relative_to(root)
        paragraphs = [
            paragraph.strip()
            for paragraph in re.split(r"\n\s*\n", readme.read_text(encoding="utf-8"))
            if paragraph.lstrip().startswith(OUTPUT_CONTRACT)
        ]
        if not paragraphs and readme in crate_readmes:
            errors.append(f"crate README names no output contract: {relative}")
        for paragraph in paragraphs:
            names = SCHEMA_FILE.findall(paragraph)
            stated = paragraph.removeprefix(OUTPUT_CONTRACT).strip().lower()
            if not names and not stated.startswith("none"):
                errors.append(
                    f"output contract names no schema and is not 'none': {relative}"
                )
            errors.extend(
                f"output contract names a missing schema: {relative} -> {name}"
                for name in names
                if not (schemas / name).is_file()
            )
    return errors


def validate_markdown_links(path: Path, errors: list[str]) -> None:
    text = path.read_text(encoding="utf-8")
    for raw_target in MARKDOWN_LINK.findall(text):
        target = raw_target.strip().split(maxsplit=1)[0].strip("<>")
        if not target or target.startswith(("#", "http://", "https://", "mailto:")):
            continue

        target = target.split("#", 1)[0]
        if not target:
            continue

        resolved = (path.parent / target).resolve()
        try:
            resolved.relative_to(ROOT)
        except ValueError:
            errors.append(
                f"markdown link escapes repository: {path.relative_to(ROOT)} -> {raw_target}"
            )
            continue

        if not resolved.exists():
            errors.append(
                f"broken markdown link: {path.relative_to(ROOT)} -> {raw_target}"
            )


def main() -> None:
    errors: list[str] = []

    if not (ROOT / "README.md").is_file():
        errors.append("README.md is required")
    if not (DOCS / "README.md").is_file():
        errors.append("docs/README.md is required")
    if (DOCS / "src").exists():
        errors.append(
            "docs/src is not allowed; implementation ownership belongs in "
            "source-directory README.md files"
        )

    documentation_paths = sorted(DOCS.rglob("*"))

    for path in documentation_paths:
        if path.is_file() and path.stem.lower() in FORBIDDEN_NAMES:
            errors.append(
                f"legacy/temporary documentation file is not allowed: "
                f"{path.relative_to(ROOT)}"
            )

    for directory in (path for path in documentation_paths if path.is_dir()):
        relative = directory.relative_to(ROOT)

        if directory.name.lower() in FORBIDDEN_NAMES:
            errors.append(f"legacy/temporary docs directory is not allowed: {relative}")

        if not (directory / "README.md").is_file():
            errors.append(
                f"documentation folder is missing README.md index: {relative}"
            )

        companion = directory.parent / f"{directory.name}.md"
        if companion.is_file():
            errors.append(
                "duplicate documentation entry points are not allowed: "
                f"{companion.relative_to(ROOT)} and {relative / 'README.md'}"
            )

    source_directories = [
        directory
        for source in SOURCES
        for directory in [
            source,
            *sorted(path for path in source.rglob("*") if path.is_dir()),
        ]
    ]
    for directory in source_directories:
        if not (directory / "README.md").is_file():
            errors.append(
                "source module/package directory is missing README.md: "
                f"{directory.relative_to(ROOT)}"
            )
    errors.extend(output_contract_violations(ROOT))

    markdown_files = [
        ROOT / "README.md",
        ROOT / "AGENTS.md",
        *sorted(DOCS.rglob("*.md")),
        *sorted((ROOT / "crates").rglob("README.md")),
        *sorted((ROOT / "src").rglob("README.md")),
    ]
    for markdown_file in markdown_files:
        validate_markdown_links(markdown_file, errors)

    if errors:
        raise SystemExit("\n".join(errors))


if __name__ == "__main__":
    main()
