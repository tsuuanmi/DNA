"""Validate DNA analysis, basecall, and sample-evidence result contracts."""

from __future__ import annotations

import argparse
import copy
import json
import sys
from pathlib import Path
from typing import Any

from jsonschema import Draft202012Validator, SchemaError, ValidationError

ROOT = Path(__file__).resolve().parents[3]
CONTRACTS = ROOT / "docs" / "reference"
ANALYSIS_SCHEMA = CONTRACTS / "schemas" / "analysis-v7.schema.json"
ANALYSIS_EXAMPLE = CONTRACTS / "examples" / "analysis-v7.example.json"
BASECALL_SCHEMA = CONTRACTS / "schemas" / "basecalls-v2.schema.json"
BASECALL_EXAMPLE = CONTRACTS / "examples" / "basecalls-v2.example.json"
SAMPLE_SCHEMA = CONTRACTS / "schemas" / "sample-evidence-v8.schema.json"
SAMPLE_EXAMPLE = CONTRACTS / "examples" / "sample-evidence-v8.example.json"


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def validator(path: Path, errors: list[str]) -> Draft202012Validator:
    schema = load_json(path)
    built = Draft202012Validator(schema)
    try:
        built.check_schema(schema)
    except SchemaError as exc:
        errors.append(f"schema {path} is not valid Draft 2020-12: {exc.message}")
    return built


def validate_documents(
    built: Draft202012Validator, paths: list[Path], errors: list[str]
) -> None:
    for path in paths:
        try:
            built.validate(load_json(path))
        except ValidationError as exc:
            errors.append(f"{path}: invalid: {exc.message}")


def analysis_call_shapes(example: dict[str, Any]) -> tuple[dict[str, Any], ...]:
    """Return valid SNV, insertion, and deletion analysis documents."""
    snv = copy.deepcopy(example)
    supporting = copy.deepcopy(snv["variants"][0]["calls"][0])
    flanking = copy.deepcopy(supporting)
    flanking["role"] = "flanking"
    flanking["base"] = "A"

    insertion = copy.deepcopy(example)
    insertion["variants"][0].update(
        {
            "reference": "A",
            "alternate": "AG",
            "kind": "INS",
            "calls": [flanking, supporting],
        }
    )
    deletion = copy.deepcopy(example)
    deletion["variants"][0].update(
        {"reference": "AG", "alternate": "A", "kind": "DEL", "calls": [flanking]}
    )
    return snv, insertion, deletion


def rejected_analysis_shapes(
    example: dict[str, Any],
) -> list[tuple[str, dict[str, Any]]]:
    snv, _, deletion = analysis_call_shapes(example)
    supporting = copy.deepcopy(snv["variants"][0]["calls"][0])
    flanking = copy.deepcopy(deletion["variants"][0]["calls"][0])

    def document(kind: str, calls: list[dict[str, Any]]) -> dict[str, Any]:
        built = copy.deepcopy(example)
        built["variants"][0]["kind"] = kind
        built["variants"][0]["calls"] = calls
        return built

    missing_peaks = copy.deepcopy(supporting)
    missing_peaks.pop("peaks")
    missing_quality = copy.deepcopy(supporting)
    missing_quality.pop("quality")
    legacy_pointer = copy.deepcopy(supporting)
    legacy_pointer["ploc"] = 123
    missing_metric = copy.deepcopy(example)
    missing_metric["signal_quality"]["noisy_regions"][0].pop("minimum_primary_snr")
    negative_metric = copy.deepcopy(example)
    negative_metric["signal_quality"]["noisy_regions"][0]["minimum_primary_snr"] = -1
    unknown_field = copy.deepcopy(example)
    unknown_field["signal_quality"]["windows"] = []
    removed_section = copy.deepcopy(example)
    removed_section["sequence"] = {"primary": "ACGT"}
    removed_software_version = copy.deepcopy(example)
    removed_software_version["provenance"]["software_version"] = "0.1.0"
    old_schema = copy.deepcopy(example)
    old_schema["schema_version"] = "dna.analysis/v6"
    missing_integrity = copy.deepcopy(example)
    missing_integrity["signal_quality"].pop("integrity")
    invalid_integrity_ratio = copy.deepcopy(example)
    invalid_integrity_ratio["signal_quality"]["integrity"][
        "maximum_to_median_event_signal_ratio"
    ] = 0.5
    excessive_vendor_mismatches = copy.deepcopy(example)
    excessive_vendor_mismatches["warnings"]["ploc_vendor_length_mismatches"] = 3

    return [
        ("SNV with no calls", document("SNV", [])),
        ("SNV flanking call", document("SNV", [flanking])),
        ("INS with only flanking calls", document("INS", [flanking])),
        ("DEL supporting call", document("DEL", [supporting])),
        ("variant call without peaks", document("SNV", [missing_peaks])),
        ("variant call without quality", document("SNV", [missing_quality])),
        ("variant call with legacy pointer", document("SNV", [legacy_pointer])),
        ("noisy region missing primary SNR", missing_metric),
        ("noisy region with negative primary SNR", negative_metric),
        ("signal quality with removed windows", unknown_field),
        ("document with removed sequence section", removed_section),
        ("analysis provenance with removed software version", removed_software_version),
        ("analysis using old schema version", old_schema),
        ("analysis without trace integrity", missing_integrity),
        ("analysis with invalid event-signal ratio", invalid_integrity_ratio),
        (
            "analysis with more than two vendor length mismatches",
            excessive_vendor_mismatches,
        ),
    ]


def rejected_basecall_shapes(
    example: dict[str, Any],
) -> list[tuple[str, dict[str, Any]]]:
    invalid_primary = copy.deepcopy(example)
    invalid_primary["read"]["primary"] = "ACGTX"
    empty_retained = copy.deepcopy(example)
    empty_retained["read"]["retained"] = ""
    negative_trim = copy.deepcopy(example)
    negative_trim["read"]["trim"]["start"] = -1
    unknown_field = copy.deepcopy(example)
    unknown_field["read"]["relative_quality"] = []
    reference = copy.deepcopy(example)
    reference["provenance"]["reference"] = {"name": "unexpected"}
    software_version = copy.deepcopy(example)
    software_version["provenance"]["software_version"] = "0.1.0"
    old_schema = copy.deepcopy(example)
    old_schema["schema_version"] = "dna.basecalls/v1"
    missing_integrity = copy.deepcopy(example)
    missing_integrity["signal_quality"].pop("integrity")
    invalid_single_ploc_spacing = copy.deepcopy(example)
    integrity = invalid_single_ploc_spacing["signal_quality"]["integrity"]
    integrity["ploc_count"] = 1
    excessive_vendor_mismatches = copy.deepcopy(example)
    excessive_vendor_mismatches["warnings"]["ploc_vendor_length_mismatches"] = 3

    return [
        ("basecall primary with unsupported symbol", invalid_primary),
        ("basecall empty retained sequence", empty_retained),
        ("basecall negative trim start", negative_trim),
        ("basecall read with unknown field", unknown_field),
        ("basecall provenance with reference", reference),
        ("basecall provenance with software version", software_version),
        ("basecall using old schema version", old_schema),
        ("basecall without trace integrity", missing_integrity),
        ("single-PLOC basecall carrying spacing summary", invalid_single_ploc_spacing),
        (
            "basecall with more than two vendor length mismatches",
            excessive_vendor_mismatches,
        ),
    ]


def validate_sample_support_topology_document(
    document: dict[str, Any], label: str, errors: list[str]
) -> None:
    reads = document.get("reads")
    if not isinstance(reads, list):
        return

    orientations: dict[str, str] = {}
    for read in reads:
        if not isinstance(read, dict):
            continue
        name = read.get("name")
        alignment = read.get("alignment")
        orientation = (
            alignment.get("orientation") if isinstance(alignment, dict) else None
        )
        if isinstance(name, str) and orientation in {"forward", "reverse"}:
            if name in orientations:
                errors.append(f"{label}: duplicate read name {name!r}")
                continue
            orientations[name] = orientation

    loci = document.get("locus_differences")
    if isinstance(loci, list):
        for index, locus in enumerate(loci):
            if not isinstance(locus, dict):
                continue
            observations = locus.get("observations")
            topology = locus.get("support_topology")
            if not isinstance(observations, list) or not isinstance(topology, dict):
                continue

            expected = {
                "reads": len(observations),
                "forward_reads": 0,
                "reverse_reads": 0,
                "reference_reads": 0,
                "alternate_reads": 0,
                "unresolved_reads": 0,
                "deletion_reads": 0,
                "profile_reads": 0,
                "profile_forward_reads": 0,
                "profile_reverse_reads": 0,
            }
            seen_reads: set[str] = set()
            valid = True
            for observation in observations:
                if not isinstance(observation, dict):
                    valid = False
                    continue
                read_name = observation.get("read")
                state = observation.get("state")
                if not isinstance(read_name, str) or read_name not in orientations:
                    errors.append(
                        f"{label}: locus {index} observation references unknown read {read_name!r}"
                    )
                    valid = False
                    continue
                if read_name in seen_reads:
                    errors.append(
                        f"{label}: locus {index} repeats read {read_name!r} in observations"
                    )
                    valid = False
                    continue
                seen_reads.add(read_name)
                orientation = orientations[read_name]
                expected[f"{orientation}_reads"] += 1
                if state in {"reference", "alternate", "unresolved", "deletion"}:
                    expected[f"{state}_reads"] += 1
                else:
                    valid = False

                profile = observation.get("profile")
                if isinstance(profile, dict):
                    expected["profile_reads"] += 1
                    expected[f"profile_{orientation}_reads"] += 1
                    weights: list[float] = []
                    for base in ("A", "C", "G", "T"):
                        value = profile.get(base)
                        if not isinstance(value, (int, float)):
                            valid = False
                            break
                        weights.append(float(value))
                    if len(weights) == 4 and abs(sum(weights) - 1.0) > 1e-9:
                        errors.append(
                            f"{label}: locus {index} observation profile does not sum to one"
                        )

            if valid and any(
                topology.get(key) != value for key, value in expected.items()
            ):
                errors.append(
                    f"{label}: locus {index} support_topology does not match observations/read orientation"
                )

    variants = document.get("variants")
    if not isinstance(variants, list):
        return

    for index, variant in enumerate(variants):
        if not isinstance(variant, dict):
            continue
        support = variant.get("support")
        topology = variant.get("support_topology")
        if not isinstance(support, list) or not isinstance(topology, dict):
            continue

        expected = {
            "reads": len(support),
            "eligible_reads": 0,
            "forward_reads": 0,
            "reverse_reads": 0,
            "eligible_forward_reads": 0,
            "eligible_reverse_reads": 0,
        }
        seen_support_reads: set[str] = set()
        valid = True
        for item in support:
            if not isinstance(item, dict):
                valid = False
                continue
            read_name = item.get("read")
            if not isinstance(read_name, str) or read_name not in orientations:
                errors.append(
                    f"{label}: variant {index} support references unknown read {read_name!r}"
                )
                valid = False
                continue
            if read_name in seen_support_reads:
                errors.append(
                    f"{label}: variant {index} repeats read {read_name!r} in support"
                )
                valid = False
                continue
            seen_support_reads.add(read_name)
            orientation = orientations[read_name]
            expected[f"{orientation}_reads"] += 1
            if item.get("eligible") is True:
                expected["eligible_reads"] += 1
                expected[f"eligible_{orientation}_reads"] += 1

        if valid and any(topology.get(key) != value for key, value in expected.items()):
            errors.append(
                f"{label}: variant {index} support_topology does not match support/read orientation"
            )


def validate_sample_support_topology(paths: list[Path], errors: list[str]) -> None:
    for path in paths:
        document = load_json(path)
        if isinstance(document, dict):
            validate_sample_support_topology_document(document, str(path), errors)


def rejected_sample_shapes(
    example: dict[str, Any],
) -> list[tuple[str, dict[str, Any]]]:
    missing_reads = copy.deepcopy(example)
    missing_reads["reads"] = []

    missing_overlaps = copy.deepcopy(example)
    missing_overlaps.pop("overlaps")

    eligible_overlap_with_reason = copy.deepcopy(example)
    eligible_overlap_with_reason["overlaps"][0]["exclusion_reasons"] = [
        "comparable_bases_below_minimum"
    ]

    ineligible_overlap_without_reason = copy.deepcopy(example)
    ineligible_overlap_without_reason["overlaps"][0]["eligible"] = False

    invalid_overlap_agreement = copy.deepcopy(example)
    invalid_overlap_agreement["overlaps"][0]["agreement"] = 1.1

    missing_overlap_agreement = copy.deepcopy(example)
    missing_overlap_agreement["overlaps"][0].pop("agreement")

    zero_comparable_with_agreement = copy.deepcopy(example)
    zero_comparable_with_agreement["overlaps"][0]["comparable_bases"] = 0
    zero_comparable_with_agreement["overlaps"][0]["agreements"] = 0
    zero_comparable_with_agreement["overlaps"][0]["conflicts"] = 0

    old_sample_schema = copy.deepcopy(example)
    old_sample_schema["schema_version"] = "dna.sample_evidence/v7"

    missing_locus_support_topology = copy.deepcopy(example)
    missing_locus_support_topology["locus_differences"][0].pop("support_topology")
    missing_locus_noisy_context = copy.deepcopy(example)
    missing_locus_noisy_context["locus_differences"][0]["observations"][0].pop(
        "in_noisy_region"
    )
    out_of_range_profile = copy.deepcopy(example)
    out_of_range_profile["locus_differences"][0]["observations"][0]["profile"]["A"] = (
        1.1
    )

    missing_support_topology = copy.deepcopy(example)
    missing_support_topology["variants"][0].pop("support_topology")
    zero_support_topology_reads = copy.deepcopy(example)
    zero_support_topology_reads["variants"][0]["support_topology"]["reads"] = 0

    missing_coverage = copy.deepcopy(example)
    missing_coverage.pop("coverage")
    empty_coverage = copy.deepcopy(example)
    empty_coverage["coverage"] = []
    zero_coverage_depth = copy.deepcopy(example)
    zero_coverage_depth["coverage"][0]["read_depth"] = 0

    missing_read_integrity = copy.deepcopy(example)
    missing_read_integrity["reads"][0].pop("integrity")

    invalid_sample_id = copy.deepcopy(example)
    invalid_sample_id["sample_id"] = "../sample"

    all_reference_locus = copy.deepcopy(example)
    for observation in all_reference_locus["locus_differences"][0]["observations"]:
        observation["state"] = "reference"
        observation["base"] = all_reference_locus["locus_differences"][0]["reference"]

    verbose_deletion = copy.deepcopy(example)
    observation = verbose_deletion["locus_differences"][0]["observations"][0]
    observation["state"] = "deletion"

    empty_read = copy.deepcopy(example)
    empty_read["variants"][0]["support"][0]["read"] = ""

    repeated_identity = copy.deepcopy(example)
    repeated_identity["variants"][0]["support"][0]["read_name"] = "legacy.ab1"

    legacy_loci = copy.deepcopy(example)
    legacy_loci["loci"] = []
    legacy_loci.pop("locus_differences")

    unknown_field = copy.deepcopy(example)
    unknown_field["consensus"] = "ACGT"

    empty_support = copy.deepcopy(example)
    empty_support["variants"][0]["support"] = []

    eligible_with_reason = copy.deepcopy(example)
    eligible_with_reason["variants"][0]["support"][0]["exclusion_reasons"] = [
        "peak_below_minimum"
    ]

    ineligible_without_reason = copy.deepcopy(example)
    ineligible_without_reason["variants"][0]["support"][0]["eligible"] = False

    missing_call_peaks = copy.deepcopy(example)
    missing_call_peaks["variants"][0]["support"][0]["calls"][0].pop("peaks")

    empty_variant_calls = copy.deepcopy(example)
    empty_variant_calls["variants"][0]["support"][0]["calls"] = []

    unknown_exclusion_reason = copy.deepcopy(example)
    unknown_exclusion_reason["variants"][0]["support"][0]["eligible"] = False
    unknown_exclusion_reason["variants"][0]["support"][0]["exclusion_reasons"] = [
        "unsupported_reason"
    ]

    return [
        ("sample evidence with no reads", missing_reads),
        ("sample evidence without overlap graph", missing_overlaps),
        (
            "eligible overlap with exclusion reason",
            eligible_overlap_with_reason,
        ),
        (
            "ineligible overlap without exclusion reason",
            ineligible_overlap_without_reason,
        ),
        ("overlap agreement above one", invalid_overlap_agreement),
        ("overlap with comparable bases but no agreement", missing_overlap_agreement),
        ("zero-comparable overlap with agreement", zero_comparable_with_agreement),
        ("sample evidence using old schema version", old_sample_schema),
        ("sample locus without support topology", missing_locus_support_topology),
        ("sample called locus without noisy context", missing_locus_noisy_context),
        ("sample locus profile with out-of-range channel", out_of_range_profile),
        ("sample variant without support topology", missing_support_topology),
        ("sample variant topology with zero reads", zero_support_topology_reads),
        ("sample evidence without coverage topology", missing_coverage),
        ("sample evidence with empty coverage topology", empty_coverage),
        ("sample coverage with zero read depth", zero_coverage_depth),
        ("sample read without trace integrity", missing_read_integrity),
        ("sample evidence with invalid sample id", invalid_sample_id),
        (
            "sparse difference locus with only reference observations",
            all_reference_locus,
        ),
        ("deletion difference carrying called-base fields", verbose_deletion),
        ("sample evidence with empty read reference", empty_read),
        ("sample support with repeated read identity", repeated_identity),
        ("sample evidence using removed loci field", legacy_loci),
        ("sample evidence with consensus field", unknown_field),
        ("sample variant with no supporting reads", empty_support),
        ("eligible sample variant support with exclusion reason", eligible_with_reason),
        (
            "ineligible sample variant support without exclusion reason",
            ineligible_without_reason,
        ),
        ("sample variant call without peaks", missing_call_peaks),
        ("sample variant support without mapped calls", empty_variant_calls),
        (
            "sample variant support with unknown exclusion reason",
            unknown_exclusion_reason,
        ),
    ]


def assert_rejected(
    built: Draft202012Validator,
    rejected: list[tuple[str, dict[str, Any]]],
    errors: list[str],
) -> None:
    for label, document in rejected:
        try:
            built.validate(document)
        except ValidationError:
            continue
        errors.append(f"expected {label} to be rejected, but it validated")


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--analysis",
        action="append",
        type=Path,
        default=[],
        metavar="RESULT",
        help="analysis result to validate; may be repeated",
    )
    parser.add_argument(
        "--basecalls",
        action="append",
        type=Path,
        default=[],
        metavar="RESULT",
        help="basecall result to validate; may be repeated",
    )
    parser.add_argument(
        "--sample-evidence",
        action="append",
        type=Path,
        default=[],
        metavar="RESULT",
        help="sample-evidence result to validate; may be repeated",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    errors: list[str] = []
    analysis_validator = validator(ANALYSIS_SCHEMA, errors)
    basecall_validator = validator(BASECALL_SCHEMA, errors)
    sample_validator = validator(SAMPLE_SCHEMA, errors)
    analysis_paths = args.analysis or [ANALYSIS_EXAMPLE]
    basecall_paths = args.basecalls or [BASECALL_EXAMPLE]
    sample_paths = args.sample_evidence or [SAMPLE_EXAMPLE]
    validate_documents(analysis_validator, analysis_paths, errors)
    validate_documents(basecall_validator, basecall_paths, errors)
    validate_documents(sample_validator, sample_paths, errors)
    validate_sample_support_topology(sample_paths, errors)

    analysis_example = load_json(ANALYSIS_EXAMPLE)
    valid_shapes = analysis_call_shapes(analysis_example)
    for index, document in enumerate(valid_shapes, start=1):
        try:
            analysis_validator.validate(document)
        except ValidationError as exc:
            errors.append(
                f"valid analysis call shape {index} was rejected: {exc.message}"
            )
    rejected_analysis = rejected_analysis_shapes(analysis_example)
    rejected_basecalls = rejected_basecall_shapes(load_json(BASECALL_EXAMPLE))
    sample_example = load_json(SAMPLE_EXAMPLE)
    rejected_samples = rejected_sample_shapes(sample_example)
    invalid_profile_mass = copy.deepcopy(sample_example)
    invalid_profile_mass["locus_differences"][0]["observations"][0]["profile"] = {
        "A": 0.4,
        "C": 0.4,
        "G": 0.4,
        "T": 0.0,
    }
    invalid_profile_errors: list[str] = []
    validate_sample_support_topology_document(
        invalid_profile_mass,
        "synthetic invalid locus profile mass",
        invalid_profile_errors,
    )
    if not invalid_profile_errors:
        errors.append("expected invalid sample locus profile mass to be rejected")

    inconsistent_locus_topology = copy.deepcopy(sample_example)
    inconsistent_locus_topology["locus_differences"][0]["support_topology"][
        "forward_reads"
    ] += 1
    inconsistent_locus_errors: list[str] = []
    validate_sample_support_topology_document(
        inconsistent_locus_topology,
        "synthetic inconsistent locus support topology",
        inconsistent_locus_errors,
    )
    if not inconsistent_locus_errors:
        errors.append(
            "expected inconsistent sample locus support topology to be rejected"
        )

    inconsistent_topology = copy.deepcopy(sample_example)
    inconsistent_topology["variants"][0]["support_topology"]["forward_reads"] += 1
    semantic_rejection_errors: list[str] = []
    validate_sample_support_topology_document(
        inconsistent_topology,
        "synthetic inconsistent support topology",
        semantic_rejection_errors,
    )
    if not semantic_rejection_errors:
        errors.append("expected inconsistent sample support topology to be rejected")

    duplicate_support = copy.deepcopy(sample_example)
    duplicate_support["variants"][0]["support"].append(
        copy.deepcopy(duplicate_support["variants"][0]["support"][0])
    )
    duplicate_support_errors: list[str] = []
    validate_sample_support_topology_document(
        duplicate_support,
        "synthetic duplicate variant support",
        duplicate_support_errors,
    )
    if not duplicate_support_errors:
        errors.append("expected duplicate sample variant support read to be rejected")

    unknown_support = copy.deepcopy(sample_example)
    unknown_support["variants"][0]["support"][0]["read"] = "unknown-read"
    unknown_support_errors: list[str] = []
    validate_sample_support_topology_document(
        unknown_support,
        "synthetic unknown support read",
        unknown_support_errors,
    )
    if not unknown_support_errors:
        errors.append("expected unknown sample variant support read to be rejected")

    duplicate_read_name = copy.deepcopy(sample_example)
    duplicate_read_name["reads"][1]["name"] = duplicate_read_name["reads"][0]["name"]
    duplicate_read_errors: list[str] = []
    validate_sample_support_topology_document(
        duplicate_read_name,
        "synthetic duplicate read name",
        duplicate_read_errors,
    )
    if not duplicate_read_errors:
        errors.append("expected duplicate sample read name to be rejected")

    assert_rejected(analysis_validator, rejected_analysis, errors)
    assert_rejected(basecall_validator, rejected_basecalls, errors)
    assert_rejected(sample_validator, rejected_samples, errors)

    for error in errors:
        print(f"FAIL: {error}", file=sys.stderr)
    if errors:
        print(f"{len(errors)} check(s) failed", file=sys.stderr)
        return 1
    print(
        f"OK: validated {len(analysis_paths)} analysis, {len(basecall_paths)} basecall, "
        f"and {len(sample_paths)} sample-evidence document(s), {len(valid_shapes)} "
        f"analysis call shapes; rejected "
        f"{len(rejected_analysis) + len(rejected_basecalls) + len(rejected_samples)} "
        "invalid shape(s)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
