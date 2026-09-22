import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
DOCS = ROOT / "docs"
SOURCE = ROOT / "src"
FORBIDDEN_NAMES = {
    "archive",
    "archives",
    "deprecated",
    "history",
    "legacy",
    "old",
    "temp",
    "temporary",
}
MARKDOWN_LINK = re.compile(r"!?\[[^\]]*\]\(([^)]+)\)")


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
        SOURCE,
        *sorted(path for path in SOURCE.rglob("*") if path.is_dir()),
    ]
    for directory in source_directories:
        if not (directory / "README.md").is_file():
            errors.append(
                "source module/package directory is missing README.md: "
                f"{directory.relative_to(ROOT)}"
            )

    markdown_files = [
        ROOT / "README.md",
        ROOT / "AGENTS.md",
        *sorted(DOCS.rglob("*.md")),
        *sorted(SOURCE.rglob("README.md")),
    ]
    for markdown_file in markdown_files:
        validate_markdown_links(markdown_file, errors)

    if errors:
        raise SystemExit("\n".join(errors))


if __name__ == "__main__":
    main()
