#!/usr/bin/env python3
"""Cleanly analyze selected manifest samples into per-sample result folders."""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SAMPLE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
type Workload = dict[str, list[Path]]


def resolved(path: Path) -> Path:
    """Resolve a path relative to the repository root."""
    return path.resolve() if path.is_absolute() else (ROOT / path).resolve()


def regular_file(path: Path, label: str) -> Path:
    """Require an existing regular file."""
    path = resolved(path)
    if not path.is_file():
        raise ValueError(f"{label} is not a regular file: {path}")
    return path


def displayed(path: Path) -> Path:
    """Prefer a repository-relative path for status output."""
    try:
        return path.relative_to(ROOT)
    except ValueError:
        return path


def sample_ids(manifest: Path, limit: int) -> list[str]:
    """Read and validate the first non-empty sample IDs."""
    if limit < 1:
        raise ValueError("--limit must be at least 1")
    samples = [
        line.strip() for line in manifest.read_text(encoding="utf-8").splitlines()
    ]
    samples = [sample for sample in samples if sample]
    selected = samples[:limit]
    if len(selected) < limit:
        raise ValueError(
            f"manifest contains {len(selected)} non-empty sample IDs; {limit} requested"
        )
    invalid = [sample for sample in selected if SAMPLE_ID.fullmatch(sample) is None]
    if invalid:
        raise ValueError(f"invalid sample ID: {invalid[0]!r}")
    if len(set(selected)) != len(selected):
        raise ValueError("the selected sample IDs are not unique")
    return selected


def discover_workload(trace_dir: Path, selected: list[str]) -> Workload:
    """Discover every selected trace and reject ambiguous batch identities."""
    workload: Workload = {}
    owners: dict[Path, str] = {}
    for sample in selected:
        candidates = sorted(trace_dir.glob(f"*_{sample}_*.ab1"))
        symlink = next((trace for trace in candidates if trace.is_symlink()), None)
        if symlink is not None:
            raise ValueError(f"refusing symlinked trace: {symlink}")
        traces = [trace.resolve() for trace in candidates if trace.is_file()]
        if not traces:
            raise ValueError(f"{sample}: no matching AB1 files")
        for trace in traces:
            if owner := owners.get(trace):
                raise ValueError(
                    f"trace {trace.name} matches multiple samples: {owner} and {sample}"
                )
            owners[trace] = sample
        workload[sample] = traces
    return workload


def paths_overlap(left: Path, right: Path) -> bool:
    """Return whether either resolved path contains the other."""
    return left == right or left in right.parents or right in left.parents


def validate_cleanup_roots(
    output_dir: Path, log_dir: Path, protected: tuple[Path, ...]
) -> None:
    """Reject destructive roots that overlap each other or protected inputs."""
    if paths_overlap(output_dir, log_dir):
        raise ValueError(f"output and log roots overlap: {output_dir} and {log_dir}")
    for root, label in ((output_dir, "output"), (log_dir, "log")):
        for protected_path in protected:
            if paths_overlap(root, protected_path):
                raise ValueError(
                    f"{label} root overlaps protected path: {root} and {protected_path}"
                )


def cleanup_targets(
    output_dir: Path,
    log_dir: Path,
    workload: Workload,
    protected: tuple[Path, ...] = (),
) -> tuple[list[Path], list[Path]]:
    """Validate and return selected result directories and log files to remove."""
    validate_cleanup_roots(output_dir, log_dir, protected)
    result_targets: list[Path] = []
    log_targets: set[Path] = set()
    for sample in workload:
        result = output_dir / sample
        if result.is_symlink():
            raise ValueError(f"refusing symlinked result directory: {result}")
        if result.exists():
            if not result.is_dir():
                raise ValueError(f"result target is not a directory: {result}")
            if result.parent != output_dir:
                raise ValueError(f"result target escapes output root: {result}")
            result_targets.append(result)
        sample_log = log_dir / f"{sample}.log"
        if sample_log.exists() or sample_log.is_symlink():
            log_targets.add(sample_log)

    for log in log_targets:
        if log.parent != log_dir:
            raise ValueError(f"log target escapes log root: {log}")
        if log.is_symlink():
            raise ValueError(f"refusing symlinked log target: {log}")
        if log.exists() and not log.is_file():
            raise ValueError(f"log target is not a file: {log}")
    return result_targets, sorted(log_targets)


def clean_previous_results(
    output_dir: Path,
    log_dir: Path,
    workload: Workload,
    protected: tuple[Path, ...] = (),
) -> tuple[int, int]:
    """Remove only preflighted selected-sample result directories and logs."""
    result_targets, log_targets = cleanup_targets(
        output_dir, log_dir, workload, protected
    )
    for result in result_targets:
        shutil.rmtree(result)
        print(f"CLEAN results {displayed(result)}")
    for log in log_targets:
        log.unlink(missing_ok=True)
        print(f"CLEAN log     {displayed(log)}")
    return len(result_targets), len(log_targets)


def parser() -> argparse.ArgumentParser:
    built = argparse.ArgumentParser(description=__doc__)
    built.add_argument(
        "--manifest",
        type=Path,
        default=Path("data/MS_010426_001.txt"),
        help="sample ID manifest (default: data/MS_010426_001.txt)",
    )
    built.add_argument(
        "--trace-dir",
        type=Path,
        default=Path("data/raw/MS_010426_001"),
        help="directory containing AB1 files",
    )
    built.add_argument(
        "--reference",
        type=Path,
        default=Path("references/rCRS.fasta"),
        help="reference FASTA",
    )
    built.add_argument(
        "--config",
        type=Path,
        default=Path("config/dna.toml"),
        help="DNA TOML configuration",
    )
    built.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results"),
        help="result root directory",
    )
    built.add_argument(
        "--log-dir",
        type=Path,
        default=Path("logs"),
        help="Rust log root directory",
    )
    built.add_argument(
        "--binary",
        type=Path,
        default=Path("target/release/dna"),
        help="DNA executable",
    )
    built.add_argument("--limit", type=int, default=89, help="number of sample IDs")
    built.add_argument(
        "--no-build",
        action="store_true",
        help="use the existing binary without running cargo build --release",
    )
    return built


def sync_directory(directory: Path) -> None:
    """Synchronize a directory after an atomic link or rollback."""
    descriptor = os.open(directory, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def create_result_directory(directory: Path) -> list[Path]:
    """Create a result directory and return parents needing synchronization."""
    if directory.is_symlink():
        raise OSError(f"refusing symlinked result directory: {directory}")
    missing: list[Path] = []
    current = directory
    while not current.exists():
        missing.append(current)
        if current.parent == current:
            break
        current = current.parent
    if current.exists() and not current.is_dir():
        raise OSError(f"result directory ancestor is not a directory: {current}")
    directory.mkdir(parents=True, exist_ok=True)
    if directory.is_symlink() or not directory.is_dir():
        raise OSError(f"result directory is not a regular directory: {directory}")
    return [created.parent for created in missing]


def publish_result(generated: Path, destination: Path) -> tuple[bool, str]:
    """Publish one generated result atomically without overwrite."""
    try:
        directory_parents = create_result_directory(destination.parent)
    except OSError as directory_error:
        return False, f"could not prepare {destination.parent}: {directory_error}"
    staged: Path | None = None
    published = False
    error = ""
    try:
        with (
            generated.open("rb") as source,
            tempfile.NamedTemporaryFile(
                mode="wb",
                prefix=f".{destination.stem}.",
                suffix=".tmp",
                dir=destination.parent,
                delete=False,
            ) as target,
        ):
            staged = Path(target.name)
            shutil.copyfileobj(source, target)
            target.flush()
            os.fsync(target.fileno())
        os.link(staged, destination)
        published = True
        staged.unlink()
        staged = None
        sync_directory(destination.parent)
        for parent in dict.fromkeys(directory_parents):
            sync_directory(parent)
    except FileExistsError:
        error = f"result appeared while analysis was running: {destination}"
    except OSError as publication_error:
        error = f"could not publish {destination}: {publication_error}"
        if published:
            try:
                destination.unlink()
                sync_directory(destination.parent)
            except OSError as rollback_error:
                error = f"{error}; rollback failed: {rollback_error}"
    finally:
        if staged is not None:
            try:
                staged.unlink(missing_ok=True)
            except OSError as cleanup_error:
                if not error:
                    error = f"could not remove staged result {staged}: {cleanup_error}"
    return not error, error


def run_analysis(
    binary: Path,
    trace: Path,
    reference: Path,
    config: Path,
    destination: Path,
) -> tuple[bool, str]:
    """Run one analysis in isolation and publish its JSON without overwrite."""
    with tempfile.TemporaryDirectory(prefix="dna-batch-") as temporary:
        work = Path(temporary)
        environment = os.environ.copy()
        environment["DNA_CONFIG"] = str(config)
        environment["DNA_LOG_DIR"] = str(work / "logs")
        completed = subprocess.run(
            [str(binary), "analyze", str(trace), "--reference", str(reference)],
            cwd=work,
            env=environment,
            check=False,
            capture_output=True,
            text=True,
        )
        if completed.returncode != 0:
            detail = completed.stderr.strip() or f"exit status {completed.returncode}"
            return False, detail
        generated = work / "results" / f"{trace.stem}.json"
        if not generated.is_file():
            return False, f"analysis succeeded but did not create {generated}"
        return publish_result(generated, destination)


def run_sample(
    binary: Path,
    sample: str,
    traces: list[Path],
    reference: Path,
    config: Path,
    log_dir: Path,
    destination: Path,
) -> tuple[bool, str]:
    """Aggregate one sample and publish <sample>/<sample>.json without overwrite."""
    with tempfile.TemporaryDirectory(prefix="dna-sample-") as temporary:
        work = Path(temporary)
        environment = os.environ.copy()
        environment["DNA_CONFIG"] = str(config)
        environment["DNA_LOG_DIR"] = str(log_dir)
        completed = subprocess.run(
            [
                str(binary),
                "sample",
                sample,
                *(str(trace) for trace in traces),
                "--reference",
                str(reference),
            ],
            cwd=work,
            env=environment,
            check=False,
            capture_output=True,
            text=True,
        )
        if completed.returncode != 0:
            detail = completed.stderr.strip() or f"exit status {completed.returncode}"
            return False, detail
        generated = work / "results" / f"{sample}.sample.json"
        if not generated.is_file():
            return False, f"sample analysis succeeded but did not create {generated}"
        return publish_result(generated, destination)


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        manifest = regular_file(args.manifest, "manifest")
        reference = regular_file(args.reference, "reference")
        config = regular_file(args.config, "config")
        trace_dir = resolved(args.trace_dir)
        output_dir = resolved(args.output_dir)
        log_dir = resolved(args.log_dir)
        binary = resolved(args.binary)
        if not trace_dir.is_dir():
            raise ValueError(f"trace directory does not exist: {trace_dir}")
        selected = sample_ids(manifest, args.limit)
        workload = discover_workload(trace_dir, selected)
        protected = (manifest, reference, config, trace_dir, binary)
        cleanup_targets(output_dir, log_dir, workload, protected)
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2

    if not args.no_build:
        build = subprocess.run(["cargo", "build", "--release"], cwd=ROOT, check=False)
        if build.returncode != 0:
            return build.returncode
    if not binary.is_file():
        print(f"error: DNA binary is not a regular file: {binary}", file=sys.stderr)
        return 2

    try:
        cleaned_results, cleaned_logs = clean_previous_results(
            output_dir, log_dir, workload, protected
        )
    except (OSError, ValueError) as error:
        print(f"error: cleanup failed: {error}", file=sys.stderr)
        return 2

    completed_count = 0
    failed_count = 0
    sample_completed_count = 0
    sample_failed_count = 0
    trace_count = 0
    for sample, traces in workload.items():
        sample_trace_failed = False
        for trace in traces:
            trace_count += 1
            destination = output_dir / sample / f"{trace.stem}.json"
            succeeded, detail = run_analysis(
                binary, trace, reference, config, destination
            )
            if succeeded:
                print(f"OK   {displayed(destination)}")
                completed_count += 1
            else:
                print(f"FAIL {trace.name}: {detail}", file=sys.stderr)
                failed_count += 1
                sample_trace_failed = True

        if sample_trace_failed:
            print(
                f"SKIP {sample}: sample JSON not generated because a trace failed",
                file=sys.stderr,
            )
            continue

        sample_destination = output_dir / sample / f"{sample}.json"
        succeeded, detail = run_sample(
            binary,
            sample,
            traces,
            reference,
            config,
            log_dir,
            sample_destination,
        )
        if succeeded:
            print(f"OK   {displayed(sample_destination)}")
            sample_completed_count += 1
        else:
            print(f"FAIL {sample}: {detail}", file=sys.stderr)
            sample_failed_count += 1

    print(
        "Summary: "
        f"samples={len(selected)} traces={trace_count} completed={completed_count} "
        f"failed={failed_count} sample_json={sample_completed_count} "
        f"sample_failed={sample_failed_count} cleaned_results={cleaned_results} "
        f"cleaned_logs={cleaned_logs}"
    )
    return 1 if failed_count or sample_failed_count else 0


if __name__ == "__main__":
    raise SystemExit(main())
