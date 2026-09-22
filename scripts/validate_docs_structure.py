from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
SOURCE_MIRROR = DOCS / "src"
FORBIDDEN_DIRECTORY_NAMES = {
    "archive",
    "archives",
    "deprecated",
    "history",
    "legacy",
    "old",
    "temp",
    "temporary",
}


def is_source_mirror(path: Path) -> bool:
    return path == SOURCE_MIRROR or SOURCE_MIRROR in path.parents


def main() -> None:
    errors: list[str] = []

    if not (DOCS / "README.md").is_file():
        errors.append("docs/README.md is required")

    documentation_paths = sorted(\n        path for path in DOCS.rglob("*") if not is_source_mirror(path)\n    )

    forbidden_files = [
        path
        for path in documentation_paths
        if path.is_file() and path.stem.lower() in FORBIDDEN_DIRECTORY_NAMES
    ]
    for path in forbidden_files:
        errors.append(
            f"legacy/temporary documentation file is not allowed: {path.relative_to(ROOT)}"
        )

    directories = [path for path in documentation_paths if path.is_dir()]
    for directory in directories:
        if is_source_mirror(directory):
            continue

        relative = directory.relative_to(ROOT)

        if directory.name.lower() in FORBIDDEN_DIRECTORY_NAMES:
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

    if errors:
        raise SystemExit("\n".join(errors))


if __name__ == "__main__":
    main()
