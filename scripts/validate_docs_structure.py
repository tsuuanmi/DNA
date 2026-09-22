from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
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

    if errors:
        raise SystemExit("\n".join(errors))


if __name__ == "__main__":
    main()
