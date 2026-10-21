"""Validate DNA analysis, basecall, sample-evidence, variants, notation, consensus, review, run signal, run context, read evidence, and evaluation run contracts."""

from __future__ import annotations

import argparse
import copy
import itertools
import json
import re
import sys
import tempfile
from pathlib import Path
from typing import Any

from jsonschema import Draft202012Validator, SchemaError, ValidationError
from referencing import Registry, Resource

ROOT = Path(__file__).resolve().parents[3]
CONTRACTS = ROOT / "docs" / "reference"
ANALYSIS_SCHEMA = CONTRACTS / "schemas" / "analysis-v10.schema.json"
ANALYSIS_EXAMPLE = CONTRACTS / "examples" / "analysis-v10.example.json"
BASECALL_SCHEMA = CONTRACTS / "schemas" / "basecalls-v3.schema.json"
BASECALL_EXAMPLE = CONTRACTS / "examples" / "basecalls-v3.example.json"
SAMPLE_SCHEMA = CONTRACTS / "schemas" / "sample-evidence-v10.schema.json"
SAMPLE_EXAMPLE = CONTRACTS / "examples" / "sample-evidence-v10.example.json"
VARIANTS_SCHEMA = CONTRACTS / "schemas" / "variants-v1.schema.json"
VARIANTS_EXAMPLE = CONTRACTS / "examples" / "variants-v1.example.json"
NOTATION_SCHEMA = CONTRACTS / "schemas" / "notation-v1.schema.json"
NOTATION_EXAMPLE = CONTRACTS / "examples" / "notation-v1.example.json"
CONSENSUS_SCHEMA = CONTRACTS / "schemas" / "consensus-v1.schema.json"
CONSENSUS_EXAMPLE = CONTRACTS / "examples" / "consensus-v1.example.json"
REVIEW_SCHEMA = CONTRACTS / "schemas" / "review-v1.schema.json"
REVIEW_EXAMPLE = CONTRACTS / "examples" / "review-v1.example.json"
# The labelling rules' one source: the `review` plugin's descriptor.
REVIEW_SOURCE = ROOT / "src" / "report" / "review.rs"
REVIEW_REFERENCE = CONTRACTS / "review.md"
REVIEW_PLUGIN = re.compile(
    r'PluginDescriptor = PluginDescriptor \{\s*id: "review",[^}]*?version: (\d+),'
)
RULES_HEADING = re.compile(r"^## Labels \(rules version (\d+)\)$", re.MULTILINE)
RANK_ROW = re.compile(r"^\| (\d) \| (.+) \|$", re.MULTILINE)
SIGNAL_SCHEMA = CONTRACTS / "schemas" / "signal-v1.schema.json"
SIGNAL_EXAMPLE = CONTRACTS / "examples" / "signal-v1.example.json"
RUN_CONTEXT_SCHEMA = CONTRACTS / "schemas" / "run-context-v1.schema.json"
RUN_CONTEXT_EXAMPLE = CONTRACTS / "examples" / "run-context-v1.example.json"
EVIDENCE_SCHEMA = CONTRACTS / "schemas" / "read-evidence-v1.schema.json"
EVIDENCE_EXAMPLE = CONTRACTS / "examples" / "read-evidence-v1.example.json"
REJECTED_EVIDENCE_EXAMPLE = (
    CONTRACTS / "examples" / "read-evidence-v1.rejected.example.json"
)
PLACED_SCHEMA = CONTRACTS / "schemas" / "placed-read-v1.schema.json"
PLACED_EXAMPLE = CONTRACTS / "examples" / "placed-read-v1.example.json"
EVAL_RUN_SCHEMA = CONTRACTS / "schemas" / "eval-run-v1.schema.json"
EVAL_RUN_EXAMPLE = CONTRACTS / "examples" / "eval-run-v1.example.json"


def review_label_table(schema_path: Path) -> dict[str, tuple[int, bool]]:
    """The review labels, in the summary's order, with their rank and whether
    they make an item a check, as the schema's summary lists them (`x-rank`,
    `x-check`): the published copy of DNA's label table (PROP-0004)."""
    labels = load_json(schema_path)["$defs"]["summary"]["properties"]["labels"]
    return {
        name: (
            labels["properties"][name]["x-rank"],
            labels["properties"][name]["x-check"],
        )
        for name in labels["required"]
    }


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def schema_registry() -> Registry[Any]:
    """Every result schema by its `$id`, so a schema can reuse another's definitions."""
    resources = []
    for path in sorted((CONTRACTS / "schemas").glob("*.schema.json")):
        schema = load_json(path)
        resources.append((schema["$id"], Resource.from_contents(schema)))
    return Registry().with_resources(resources)


def validator(path: Path, errors: list[str]) -> Draft202012Validator:
    schema = load_json(path)
    built = Draft202012Validator(schema, registry=schema_registry())
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
    old_schema["schema_version"] = "dna.analysis/v9"
    missing_profile = copy.deepcopy(example)
    missing_profile["provenance"].pop("profile")
    invalid_profile_id = copy.deepcopy(example)
    invalid_profile_id["provenance"]["profile"]["id"] = "Human mtDNA"
    missing_plugins = copy.deepcopy(example)
    missing_plugins["provenance"].pop("plugins")
    unknown_family = copy.deepcopy(example)
    unknown_family["provenance"]["plugins"][0]["family"] = "kernel"
    missing_integrity = copy.deepcopy(example)
    missing_integrity["signal_quality"].pop("integrity")
    invalid_integrity_ratio = copy.deepcopy(example)
    invalid_integrity_ratio["signal_quality"]["integrity"][
        "maximum_to_median_event_signal_ratio"
    ] = 0.5
    excessive_vendor_mismatches = copy.deepcopy(example)
    excessive_vendor_mismatches["warnings"]["ploc_vendor_length_mismatches"] = 3
    missing_callability = copy.deepcopy(example)
    missing_callability["read"].pop("callability")
    unknown_phase_state = copy.deepcopy(example)
    unknown_phase_state["read"]["callability"]["segments"][0]["state"] = "noisy"
    segment_without_attribution = copy.deepcopy(example)
    segment_without_attribution["read"]["callability"]["segments"][0].pop(
        "after_repeat"
    )
    callability_with_mask_array = copy.deepcopy(example)
    callability_with_mask_array["read"]["callability"]["mask"] = [False]
    callability_without_segments = copy.deepcopy(example)
    callability_without_segments["read"]["callability"]["segments"] = []
    alignment_without_masked_bases = copy.deepcopy(example)
    alignment_without_masked_bases["alignment"].pop("masked_bases")
    alignment_without_callable_segments = copy.deepcopy(example)
    alignment_without_callable_segments["alignment"].pop("callable_reference_segments")
    without_calls = copy.deepcopy(example)
    without_calls.pop("calls")
    trusted_withheld = copy.deepcopy(example)
    trusted = next(
        call for call in trusted_withheld["calls"] if call["trust"]["trusted"]
    )
    trusted["trust"]["withheld_by"] = [{"kind": "read_end"}]
    masked_without_reason = copy.deepcopy(example)
    masked_without_reason["calls"][0]["trust"]["withheld_by"] = [{"kind": "masked"}]
    untrimmed_with_reason = copy.deepcopy(example)
    untrimmed_with_reason["calls"][0]["trust"]["withheld_by"] = [
        {"kind": "untrimmed", "reason": "weak_signal"}
    ]
    unknown_withholder = copy.deepcopy(example)
    unknown_withholder["calls"][0]["trust"]["withheld_by"] = [{"kind": "edge"}]
    column_zero = copy.deepcopy(example)
    column_zero["calls"][0]["column"]["position"] = 0
    unknown_scope = copy.deepcopy(example)
    unknown_scope["calls"][0]["trust"]["scope"] = "length"
    scored_with_reason = copy.deepcopy(example)
    scored_with_reason["calls"][0]["trust"]["placement"]["reason"] = "weak_signal"
    segments = example["read"]["callability"]["segments"]
    dephased = next(
        index
        for index, segment in enumerate(segments)
        if segment["state"] == "dephased"
    )
    in_phase = next(
        index
        for index, segment in enumerate(segments)
        if segment["state"] == "in_phase"
    )
    shadow_cases: list[tuple[str, dict[str, Any]]] = []
    for description, index, offsets in [
        ("dephased segment without shadow offsets", dephased, None),
        ("in-phase segment with shadow offsets", in_phase, [-1]),
        ("shadow offset zero", dephased, [0]),
        ("duplicate shadow offsets", dephased, [-1, -1]),
        ("empty shadow offsets", dephased, []),
        ("shadow offsets without a one-call shadow", dephased, [2]),
    ]:
        mutated = copy.deepcopy(example)
        segment = mutated["read"]["callability"]["segments"][index]
        if offsets is None:
            segment.pop("shadow_offsets")
        else:
            segment["shadow_offsets"] = offsets
        shadow_cases.append((f"analysis callability {description}", mutated))

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
        ("analysis without profile identity", missing_profile),
        ("analysis with an invalid profile id", invalid_profile_id),
        ("analysis without plugin provenance", missing_plugins),
        ("analysis with an unknown plugin family", unknown_family),
        ("analysis without trace integrity", missing_integrity),
        ("analysis with invalid event-signal ratio", invalid_integrity_ratio),
        (
            "analysis with more than two vendor length mismatches",
            excessive_vendor_mismatches,
        ),
        ("analysis read without callability", missing_callability),
        ("analysis callability with unknown phase state", unknown_phase_state),
        (
            "analysis callability segment without repeat attribution",
            segment_without_attribution,
        ),
        ("analysis callability with per-position mask", callability_with_mask_array),
        ("analysis callability without segments", callability_without_segments),
        ("analysis alignment without masked bases", alignment_without_masked_bases),
        (
            "analysis alignment without callable reference segments",
            alignment_without_callable_segments,
        ),
        ("analysis without aligned calls", without_calls),
        ("analysis trusted call with a withholder", trusted_withheld),
        ("analysis masked withholder without its reason", masked_without_reason),
        ("analysis untrimmed withholder with a reason", untrimmed_with_reason),
        ("analysis withholder of an unknown kind", unknown_withholder),
        ("analysis call column at position zero", column_zero),
        ("analysis call with an unknown trust scope", unknown_scope),
        ("analysis scored call with a placement reason", scored_with_reason),
        *shadow_cases,
    ]


def validate_analysis_calls(document: Any, label: str, errors: list[str]) -> None:
    """Check that an analysis's aligned calls are distinct, in call order, and calls of the read."""
    if not isinstance(document, dict):
        return
    indices = [call.get("index") for call in document.get("calls", [])]
    if any(
        not isinstance(index, int) or not isinstance(after, int) or index >= after
        for index, after in itertools.pairwise(indices)
    ):
        errors.append(f"{label}: aligned calls are not distinct and in call order")
    call_count = document.get("read", {}).get("call_count", 0)
    if any(isinstance(index, int) and index >= call_count for index in indices):
        errors.append(f"{label}: an aligned call lies beyond the read's calls")


def rejected_variants_shapes(
    example: dict[str, Any],
) -> list[tuple[str, dict[str, Any]]]:
    no_reads = copy.deepcopy(example)
    no_reads["reads"] = []
    missing_plugins = copy.deepcopy(example)
    missing_plugins["provenance"].pop("plugins")
    modality_field = copy.deepcopy(example)
    modality_field["reads"][0]["variants"][0]["peaks"] = {"A": 1}
    invalid_reason = copy.deepcopy(example)
    invalid_reason["reads"][0]["variants"][0]["exclusion_reasons"] = ["Read End"]
    invalid_allele = copy.deepcopy(example)
    invalid_allele["reads"][0]["variants"][0]["alternate"] = "N"
    old_schema = copy.deepcopy(example)
    old_schema["schema_version"] = "dna.variants/v0"
    return [
        ("variants without reads", no_reads),
        ("variants without plugin provenance", missing_plugins),
        ("called variant carrying modality evidence", modality_field),
        ("called variant with an invalid reason label", invalid_reason),
        ("called variant with a non-canonical allele", invalid_allele),
        ("variants using an unknown schema version", old_schema),
    ]


def rejected_notation_shapes(
    example: dict[str, Any],
) -> list[tuple[str, dict[str, Any]]]:
    missing_source = copy.deepcopy(example)
    missing_source["provenance"].pop("source")
    other_source = copy.deepcopy(example)
    other_source["provenance"]["source"]["schema_version"] = "dna.sample_evidence/v10"
    unknown_rule = copy.deepcopy(example)
    unknown_rule["conformance"]["rules"].append("minimize_differences")
    empty_finding = copy.deepcopy(example)
    empty_finding["conformance"]["findings"][0]["calls"] = []
    missing_notation = copy.deepcopy(example)
    missing_notation.pop("notation")
    return [
        ("notation without its source document", missing_source),
        ("notation from a non-variants document", other_source),
        ("notation with an unknown conformance rule", unknown_rule),
        ("notation finding without calls", empty_finding),
        ("notation document without notation", missing_notation),
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
    old_schema["schema_version"] = "dna.basecalls/v2"
    missing_plugins = copy.deepcopy(example)
    missing_plugins["provenance"].pop("plugins")
    missing_integrity = copy.deepcopy(example)
    missing_integrity["signal_quality"].pop("integrity")
    invalid_single_ploc_spacing = copy.deepcopy(example)
    integrity = invalid_single_ploc_spacing["signal_quality"]["integrity"]
    integrity["ploc_count"] = 1
    excessive_vendor_mismatches = copy.deepcopy(example)
    excessive_vendor_mismatches["warnings"]["ploc_vendor_length_mismatches"] = 3
    missing_callability = copy.deepcopy(example)
    missing_callability["read"].pop("callability")
    negative_masked_calls = copy.deepcopy(example)
    negative_masked_calls["read"]["callability"]["masked_calls"] = -1

    return [
        ("basecall primary with unsupported symbol", invalid_primary),
        ("basecall empty retained sequence", empty_retained),
        ("basecall negative trim start", negative_trim),
        ("basecall read with unknown field", unknown_field),
        ("basecall provenance with reference", reference),
        ("basecall provenance with software version", software_version),
        ("basecall using old schema version", old_schema),
        ("basecall without plugin provenance", missing_plugins),
        ("basecall without trace integrity", missing_integrity),
        ("single-PLOC basecall carrying spacing summary", invalid_single_ploc_spacing),
        (
            "basecall with more than two vendor length mismatches",
            excessive_vendor_mismatches,
        ),
        ("basecall read without callability", missing_callability),
        ("basecall callability with negative masked calls", negative_masked_calls),
    ]


def validate_callability_view(
    callability: Any, call_count: Any, label: str, errors: list[str]
) -> None:
    """Check that segments partition the read and the counts agree with them."""
    if not isinstance(callability, dict) or not isinstance(call_count, int):
        return
    segments = callability.get("segments")
    span = callability.get("callable_span")
    if not isinstance(segments, list) or not isinstance(span, dict):
        return
    expected_start = 0
    masked = 0
    unmasked: list[int] = []
    for segment in segments:
        if not isinstance(segment, dict):
            return
        calls = segment.get("calls")
        if not isinstance(calls, dict):
            return
        start, end = calls.get("start"), calls.get("end")
        if not isinstance(start, int) or not isinstance(end, int):
            return
        if start != expected_start or end <= start:
            errors.append(f"{label}: callability segments do not partition the read")
            return
        expected_start = end
        offsets = segment.get("shadow_offsets")
        if isinstance(offsets, list) and offsets != sorted(offsets):
            errors.append(f"{label}: callability shadow offsets are not ascending")
        if segment.get("state") == "in_phase":
            unmasked.extend((start, end))
        else:
            masked += end - start
    if expected_start != call_count:
        errors.append(f"{label}: callability segments do not cover every call")
        return
    if callability.get("masked_calls") != masked:
        errors.append(f"{label}: callability masked_calls disagrees with segments")
    expected_span = (
        {"start": unmasked[0], "end": unmasked[-1]}
        if unmasked
        else {"start": call_count, "end": call_count}
    )
    if span != expected_span:
        errors.append(f"{label}: callability callable_span disagrees with segments")


def validate_callable_segments(alignment: Any, label: str, errors: list[str]) -> None:
    """Check that every callable reference segment lies inside one mapped segment."""
    if not isinstance(alignment, dict):
        return
    mapped = alignment.get("reference_segments")
    callable_segments = alignment.get("callable_reference_segments")
    if not isinstance(mapped, list) or not isinstance(callable_segments, list):
        return
    for segment in callable_segments:
        if not isinstance(segment, dict):
            return
        start, end = segment.get("start"), segment.get("end")
        if not isinstance(start, int) or not isinstance(end, int) or end <= start:
            errors.append(f"{label}: empty or invalid callable reference segment")
            return
        if not any(
            isinstance(outer, dict)
            and isinstance(outer.get("start"), int)
            and isinstance(outer.get("end"), int)
            and outer["start"] <= start
            and end <= outer["end"]
            for outer in mapped
        ):
            errors.append(
                f"{label}: callable reference segment {start}..{end} lies outside the mapped segments"
            )


def validate_read_callability(paths: list[Path], errors: list[str]) -> None:
    for path in paths:
        document = load_json(path)
        if not isinstance(document, dict):
            continue
        read = document.get("read")
        if isinstance(read, dict):
            validate_callability_view(
                read.get("callability"), read.get("call_count"), str(path), errors
            )
        validate_callable_segments(document.get("alignment"), str(path), errors)
        reads = [
            *(document.get("reads") or []),
            *(document.get("rejected_reads") or []),
        ]
        for item in reads:
            if isinstance(item, dict):
                validate_callable_segments(item.get("alignment"), str(path), errors)
                integrity = item.get("integrity")
                call_count = (
                    integrity.get("ploc_count") if isinstance(integrity, dict) else None
                )
                validate_callability_view(
                    item.get("callability"), call_count, str(path), errors
                )


def validate_short_rejection(
    rejection: dict[str, Any], label: str, errors: list[str]
) -> None:
    """Check that a too-short read has fewer calls than its minimum."""
    calls = rejection.get("calls")
    minimum_calls = rejection.get("minimum_calls")
    if (
        isinstance(calls, int)
        and isinstance(minimum_calls, int)
        and calls >= minimum_calls
    ):
        errors.append(f"{label}: too-short read reaches minimum_calls")


def validate_consensus_rejections(
    document: dict[str, Any], label: str, errors: list[str]
) -> None:
    """Check the counts of a consensus document's too-short rejected reads."""
    rejected = document.get("rejected_reads")
    if not isinstance(rejected, list):
        return
    for read in rejected:
        rejection = read.get("rejection") if isinstance(read, dict) else None
        if isinstance(rejection, dict) and rejection.get("reason") == "too_short":
            validate_short_rejection(rejection, label, errors)


def validate_rejection(read: dict[str, Any], label: str, errors: list[str]) -> None:
    """Check that a rejected read's counts agree with its callability view."""
    rejection = read.get("rejection")
    if isinstance(rejection, dict) and rejection.get("reason") == "too_short":
        validate_short_rejection(rejection, label, errors)
        return
    callability = read.get("callability")
    integrity = read.get("integrity")
    if not (
        isinstance(rejection, dict)
        and isinstance(callability, dict)
        and isinstance(integrity, dict)
    ):
        return
    callable_calls = rejection.get("callable_calls")
    minimum = rejection.get("minimum_callable_calls")
    calls = integrity.get("ploc_count")
    masked = callability.get("masked_calls")
    if not (
        isinstance(callable_calls, int)
        and isinstance(minimum, int)
        and isinstance(calls, int)
        and isinstance(masked, int)
    ):
        return
    if callable_calls != calls - masked:
        errors.append(
            f"{label}: rejected read callable_calls disagrees with its callability"
        )
    if callable_calls >= minimum:
        errors.append(f"{label}: rejected read reaches minimum_callable_calls")


def expected_opposition(
    document: dict[str, Any], variant: dict[str, Any]
) -> dict[str, Any] | None:
    """Recompute a variant's opposition from the read registry and its support."""
    reads = document.get("reads")
    support = variant.get("support")
    position = variant.get("position")
    reference = variant.get("reference")
    if not (
        isinstance(reads, list)
        and isinstance(support, list)
        and isinstance(position, int)
        and isinstance(reference, str)
    ):
        return None
    supporting = {item.get("read") for item in support if isinstance(item, dict)}
    start = position - 1
    end = start + len(reference) + (1 if variant.get("kind") == "INS" else 0)
    names: list[str] = []
    counts = {"forward": 0, "reverse": 0}
    for read in reads:
        if not isinstance(read, dict) or read.get("name") in supporting:
            continue
        alignment = read.get("alignment")
        if not isinstance(alignment, dict):
            return None
        segments = alignment.get("callable_reference_segments")
        orientation = alignment.get("orientation")
        if not isinstance(segments, list) or orientation not in counts:
            return None
        if all(
            any(
                isinstance(segment, dict)
                and isinstance(segment.get("start"), int)
                and isinstance(segment.get("end"), int)
                and segment["start"] <= reference_index < segment["end"]
                for segment in segments
            )
            for reference_index in range(start, end)
        ):
            names.append(str(read.get("name")))
            counts[orientation] += 1
    return {
        "reads": names,
        "forward_reads": counts["forward"],
        "reverse_reads": counts["reverse"],
    }


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

    rejected = document.get("rejected_reads")
    if isinstance(rejected, list):
        for read in rejected:
            if not isinstance(read, dict):
                continue
            name = read.get("name")
            if isinstance(name, str) and name in orientations:
                errors.append(
                    f"{label}: rejected read {name!r} is also an admitted read"
                )
            validate_rejection(read, label, errors)

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
                "masked_reads": 0,
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
                if state in {
                    "reference",
                    "alternate",
                    "unresolved",
                    "deletion",
                    "masked",
                }:
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
        opposition = expected_opposition(document, variant)
        if opposition is not None and variant.get("opposition") != opposition:
            errors.append(
                f"{label}: variant {index} opposition does not match callable coverage and support"
            )
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
    old_sample_schema["schema_version"] = "dna.sample_evidence/v9"
    sample_without_plugins = copy.deepcopy(example)
    sample_without_plugins["provenance"].pop("plugins")
    sample_plugin_version_zero = copy.deepcopy(example)
    sample_plugin_version_zero["provenance"]["plugins"][0]["version"] = 0

    unresolved_notation_call = copy.deepcopy(example)
    unresolved_notation_call["notation"]["calls"][0]["call"] = "150N"
    zero_insertion_ordinal = copy.deepcopy(example)
    zero_insertion_ordinal["notation"]["calls"][0]["call"] = "309.0C"
    notation_call_without_reads = copy.deepcopy(example)
    notation_call_without_reads["notation"]["calls"][0]["reads"] = []
    unknown_notation_style = copy.deepcopy(example)
    unknown_notation_style["notation"]["style"] = "hgvs"
    legacy_notation_policy = copy.deepcopy(example)
    legacy_notation_policy["notation"]["policy"] = "rcrs_right_aligned_control_region"
    sample_without_profile = copy.deepcopy(example)
    sample_without_profile["provenance"].pop("profile")
    extra_notation_field = copy.deepcopy(example)
    extra_notation_field["notation"]["consensus"] = []

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
    missing_read_callability = copy.deepcopy(example)
    missing_read_callability["reads"][0].pop("callability")

    invalid_sample_id = copy.deepcopy(example)
    invalid_sample_id["sample_id"] = "../sample"

    all_reference_locus = copy.deepcopy(example)
    for observation in all_reference_locus["locus_differences"][0]["observations"]:
        observation["state"] = "reference"
        observation["base"] = all_reference_locus["locus_differences"][0]["reference"]

    short_index, callable_index = (
        next(
            index
            for index, read in enumerate(example["rejected_reads"])
            if read["rejection"]["reason"] == reason
        )
        for reason in ("too_short", "callable_calls_below_minimum")
    )
    rejected_with_alignment = copy.deepcopy(example)
    rejected_with_alignment["rejected_reads"][callable_index]["alignment"] = (
        copy.deepcopy(example["reads"][0]["alignment"])
    )
    rejected_with_unknown_reason = copy.deepcopy(example)
    rejected_with_unknown_reason["rejected_reads"][callable_index]["rejection"][
        "reason"
    ] = "noisy"
    short_with_alignment = copy.deepcopy(example)
    short_with_alignment["rejected_reads"][short_index]["alignment"] = copy.deepcopy(
        example["reads"][0]["alignment"]
    )
    missing_rejected_reads = copy.deepcopy(example)
    missing_rejected_reads.pop("rejected_reads")
    short_with_integrity = copy.deepcopy(example)
    short_with_integrity["rejected_reads"][short_index]["integrity"] = copy.deepcopy(
        example["rejected_reads"][callable_index]["integrity"]
    )
    short_with_callable_calls = copy.deepcopy(example)
    short_with_callable_calls["rejected_reads"][short_index]["rejection"][
        "callable_calls"
    ] = 0
    callable_without_callability = copy.deepcopy(example)
    callable_without_callability["rejected_reads"][callable_index].pop("callability")

    all_masked_locus = copy.deepcopy(example)
    for observation in all_masked_locus["locus_differences"][0]["observations"]:
        observation["state"] = "masked"
    missing_masked_bases = copy.deepcopy(example)
    missing_masked_bases["reads"][0]["alignment"].pop("masked_bases")
    missing_opposition = copy.deepcopy(example)
    missing_opposition["variants"][0].pop("opposition")
    missing_masked_reads = copy.deepcopy(example)
    missing_masked_reads["locus_differences"][0]["support_topology"].pop("masked_reads")

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
        ("sample evidence without plugin provenance", sample_without_plugins),
        ("sample evidence with plugin version zero", sample_plugin_version_zero),
        ("sample notation with an unresolved call base", unresolved_notation_call),
        ("sample notation with a zero insertion ordinal", zero_insertion_ordinal),
        ("sample notation call without supporting reads", notation_call_without_reads),
        ("sample notation with an unknown style", unknown_notation_style),
        ("sample notation with the removed policy field", legacy_notation_policy),
        ("sample evidence without profile identity", sample_without_profile),
        ("sample notation with an extra field", extra_notation_field),
        ("sample locus without support topology", missing_locus_support_topology),
        ("sample called locus without noisy context", missing_locus_noisy_context),
        ("sample locus profile with out-of-range channel", out_of_range_profile),
        ("sample variant without support topology", missing_support_topology),
        ("sample variant topology with zero reads", zero_support_topology_reads),
        ("sample evidence without coverage topology", missing_coverage),
        ("sample evidence with empty coverage topology", empty_coverage),
        ("sample coverage with zero read depth", zero_coverage_depth),
        ("sample read without trace integrity", missing_read_integrity),
        ("sample read without callability", missing_read_callability),
        ("sample locus retained only by masked observations", all_masked_locus),
        ("sample rejected read with an alignment", rejected_with_alignment),
        ("sample rejected read with an unknown reason", rejected_with_unknown_reason),
        ("sample too-short read with an alignment", short_with_alignment),
        ("sample without rejected reads", missing_rejected_reads),
        ("sample too-short read with integrity", short_with_integrity),
        ("sample too-short rejection with callable calls", short_with_callable_calls),
        (
            "sample callable-calls rejection without callability",
            callable_without_callability,
        ),
        ("sample read alignment without masked bases", missing_masked_bases),
        ("sample locus topology without masked reads", missing_masked_reads),
        ("sample variant without opposition", missing_opposition),
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


def rejected_consensus_shapes(
    example: dict[str, Any],
) -> list[tuple[str, dict[str, Any]]]:
    unknown_state = copy.deepcopy(example)
    unknown_state["sites"][0]["state"] = "adjudicated"
    lowercase_sequence = copy.deepcopy(example)
    lowercase_sequence["segments"][0]["sequence"] = "gtattt"
    iupac_call = copy.deepcopy(example)
    iupac_call["sites"][0]["call"] = "R"
    no_reads = copy.deepcopy(example)
    no_reads["reads"] = []
    missing_summary = copy.deepcopy(example)
    missing_summary.pop("summary")
    unknown_evidence = copy.deepcopy(example)
    unknown_evidence["sites"][0]["runs"] = [
        {"base": "C", "length": 9, "length_evidence": "dominant"}
    ]
    missing_runs = copy.deepcopy(example)
    missing_runs["sites"][0].pop("runs")
    short_without_minimum = copy.deepcopy(example)
    for read in short_without_minimum["rejected_reads"]:
        if read["rejection"]["reason"] == "too_short":
            read["rejection"].pop("minimum_calls")
    return [
        ("consensus site with an unknown state", unknown_state),
        ("consensus segment with lower-case bases", lowercase_sequence),
        ("consensus site with an IUPAC call", iupac_call),
        ("consensus without admitted reads", no_reads),
        ("consensus without a summary", missing_summary),
        ("consensus run with an unknown length evidence", unknown_evidence),
        ("consensus site without runs", missing_runs),
        ("consensus too-short rejection without minimum", short_without_minimum),
    ]


def review_rules_errors() -> list[str]:
    """The rules version the reference and the example state is the `review`
    plugin's version in the source, its one source."""
    plugin = REVIEW_PLUGIN.search(REVIEW_SOURCE.read_text(encoding="utf-8"))
    if plugin is None:
        return [f"{REVIEW_SOURCE}: no review plugin descriptor"]
    version = int(plugin.group(1))
    errors: list[str] = []
    heading = RULES_HEADING.search(REVIEW_REFERENCE.read_text(encoding="utf-8"))
    if heading is None or int(heading.group(1)) != version:
        errors.append(
            f"{REVIEW_REFERENCE}: the labels heading must state rules version {version}"
        )
    if load_json(REVIEW_EXAMPLE)["provenance"]["rules"] != version:
        errors.append(f"{REVIEW_EXAMPLE}: rules must be {version}")
    # The reference's rank table lists every label at the schema's rank.
    listed: dict[str, int] = {}
    for rank, cell in RANK_ROW.findall(REVIEW_REFERENCE.read_text(encoding="utf-8")):
        for name in re.findall(r"^`([a-z_]+)`$", cell.replace(", ", "\n"), re.M):
            listed[name] = int(rank)
    ranks = {
        name: rank for name, (rank, _) in review_label_table(REVIEW_SCHEMA).items()
    }
    if listed != ranks:
        errors.append(f"{REVIEW_REFERENCE}: the rank table disagrees with the schema")
    return errors


def validate_review(document: Any, label: str, errors: list[str]) -> None:
    """Check a review document's item order, kinds, and summary counts."""
    if not isinstance(document, dict):
        return
    items = document.get("items")
    summary = document.get("summary")
    if not isinstance(items, list) or not isinstance(summary, dict):
        return
    provenance = document.get("provenance", {})
    plugins = {
        plugin.get("id"): plugin.get("version")
        for plugin in provenance.get("plugins", [])
    }
    if provenance.get("rules") != plugins.get("review"):
        errors.append(f"{label}: rules is not the review plugin's version")
    table = review_label_table(REVIEW_SCHEMA)
    # An explained possible heteroplasmy is context: the context rank, no check.
    context = max(rank for rank, _ in table.values())
    previous: tuple[int, int, int] | None = None
    known_calls = {entry.get("call") for entry in document.get("calls", [])}
    counts = dict.fromkeys(table, 0)
    checks = 0
    for index, item in enumerate(items):
        labels = [
            (entry.get("label"), entry.get("explained") is True)
            for entry in item.get("labels", [])
            if entry.get("label") in table
        ]
        names = {name for name, _ in labels}
        priority = min(
            (context if explained else table[name][0] for name, explained in labels),
            default=context + 1,
        )
        first = item.get("first", {})
        order = (priority, first.get("position", 0), first.get("ordinal", 0))
        if previous is not None and order < previous:
            errors.append(f"{label}: item {index} is out of priority order")
        previous = order
        # Without notation calls, a difference is the consensus's decision.
        differs = bool(item.get("calls")) or "calls" not in document
        # A difference inside a read's poor signal is a check too (SRS-CALL-024).
        unassured = any(
            table[name][1] and not explained for name, explained in labels
        ) or (differs and "poor_signal" in names)
        # Without a difference, an explained possible heteroplasmy or a
        # tract_length label makes a note.
        noted = any(explained or name == "tract_length" for name, explained in labels)
        if unassured:
            kind = "check"
        elif differs:
            kind = "variant"
        else:
            kind = "note" if noted else None
        if item.get("kind") != kind:
            errors.append(
                f"{label}: item {index} kind disagrees with its calls and labels"
            )
        if not set(item.get("calls", [])) <= known_calls:
            errors.append(
                f"{label}: item {index} holds a call the document does not list"
            )
        if kind == "check":
            checks += 1
            if not names:
                errors.append(f"{label}: check item {index} has no label")
        for name in names:
            if name in counts:
                counts[name] += 1
    if summary.get("items") != len(items) or summary.get("checks") != checks:
        errors.append(f"{label}: summary item counts disagree with the items")
    variants = sum(1 for item in items if item.get("kind") == "variant")
    if summary.get("variants") != variants:
        errors.append(f"{label}: summary variant count disagrees with the items")
    if summary.get("labels") != counts:
        errors.append(f"{label}: summary label counts disagree with the items")
    positions = summary.get("positions", {})
    held = {
        position
        for item in items
        if item.get("kind") == "check"
        for position in range(
            item.get("first", {}).get("position", 0),
            item.get("last", {}).get("position", -1) + 1,
        )
    }
    if positions.get("check") != len(held) or len(held) > positions.get(
        "reportable", 0
    ):
        errors.append(f"{label}: summary positions disagree with the check items")
    clear = document.get("status") == "called" and checks == 0
    if summary.get("clear") != clear:
        errors.append(f"{label}: summary clear disagrees with the status and checks")
    confirmed = document.get("reference_only", [])
    if isinstance(confirmed, list):
        if confirmed != sorted(set(confirmed)):
            errors.append(f"{label}: reference_only positions are not increasing")
        if positions.get("reference_only") != len(confirmed):
            errors.append(
                f"{label}: summary reference_only count disagrees with the positions"
            )
    left_out = document.get("not_analysed", [])
    if isinstance(left_out, list):
        held_out = sum(stretch["last"] - stretch["first"] + 1 for stretch in left_out)
        if positions.get("not_analysed") != held_out:
            errors.append(
                f"{label}: summary not_analysed count disagrees with the stretches"
            )
    validate_zones(document, items, positions, summary, label, errors)
    validate_mixtures(document, label, errors)


def validate_mixtures(document: dict[str, Any], label: str, errors: list[str]) -> None:
    """Check that each read's mixtures lead with the dominant sequence."""
    for read in document.get("reads", []):
        for index, mixture in enumerate(read.get("mixtures", [])):
            components = mixture.get("components", [])
            offsets = [component.get("offset") for component in components]
            if not offsets or offsets[0] != 0 or len(set(offsets)) != len(offsets):
                errors.append(
                    f"{label}: mixture {index} of {read.get('name')} does not lead "
                    "with one dominant sequence"
                )
            total = sum(component.get("fraction", 0) for component in components)
            if total > 1.0 + 1e-6:
                errors.append(
                    f"{label}: mixture {index} of {read.get('name')} has fractions "
                    "above one"
                )
            unexplained = mixture.get("unexplained", [])
            if unexplained != sorted(set(unexplained)):
                errors.append(
                    f"{label}: mixture {index} of {read.get('name')} lists unexplained "
                    "positions out of order"
                )


def validate_zones(
    document: dict[str, Any],
    items: list[Any],
    positions: dict[str, Any],
    summary: dict[str, Any],
    label: str,
    errors: list[str],
) -> None:
    """Check that zones tile the analysed reportable positions and count the checks."""
    zones = document.get("zones", [])
    if not isinstance(zones, list):
        return
    covered = 0
    previous: dict[str, Any] | None = None
    for index, zone in enumerate(zones):
        first, last = zone.get("first", 0), zone.get("last", -1)
        covered += last - first + 1
        if previous is not None:
            if first <= previous.get("last", 0):
                errors.append(f"{label}: zone {index} is out of order or overlaps")
            elif (
                first == previous.get("last", 0) + 1
                and zone.get("kind") == previous.get("kind")
                and zone.get("name") == previous.get("name")
            ):
                errors.append(f"{label}: zone {index} continues the zone before it")
        previous = zone
        touching = sum(
            1
            for item in items
            if item.get("kind") == "check"
            and item.get("first", {}).get("position", 0) <= last
            and first <= item.get("last", {}).get("position", -1)
        )
        if zone.get("checks") != touching:
            errors.append(f"{label}: zone {index} check count disagrees with the items")
    # The stretches the call does not analyse (SRS-CONS-012) belong to no
    # zone, so the zones tile the reportable positions without them.
    analysed = positions.get("reportable", 0) - positions.get("not_analysed", 0)
    if zones and covered != analysed:
        errors.append(f"{label}: zones do not tile the analysed reportable positions")
    counts = summary.get("zones", {})
    flagged = sum(1 for zone in zones if zone.get("checks", 0) > 0)
    if counts.get("total") != len(zones) or counts.get("check") != flagged:
        errors.append(f"{label}: summary zone counts disagree with the zones")


def rejected_review_shapes(
    example: dict[str, Any],
) -> list[tuple[str, dict[str, Any]]]:
    unknown_label = copy.deepcopy(example)
    unknown_label["items"][0]["labels"].append({"label": "watched"})
    failed_without_error = copy.deepcopy(example)
    failed_without_error["status"] = "call_failed"
    for key in (
        "calls",
        "callable",
        "reference_only",
        "not_analysed",
        "empop",
        "zones",
    ):
        failed_without_error.pop(key)
    called_with_error = copy.deepcopy(example)
    called_with_error["error"] = "the consensus has no segment"
    rejected_with_items = copy.deepcopy(example)
    rejected_with_items["status"] = "rejected"
    for key in (
        "calls",
        "callable",
        "reference_only",
        "not_analysed",
        "empop",
        "zones",
    ):
        rejected_with_items.pop(key)
    unlisted_call = copy.deepcopy(example)
    for item in unlisted_call["items"]:
        if item["calls"]:
            item["calls"].append("999X")
            break
    missing_count = copy.deepcopy(example)
    missing_count["summary"]["labels"].pop("read_edge")
    compared_unexplained = copy.deepcopy(example)
    for read in compared_unexplained["reads"]:
        for mixture in read.get("mixtures", []):
            mixture["explained"] = False
    uncompared_explained = copy.deepcopy(example)
    for read in uncompared_explained["reads"]:
        for mixture in read.get("mixtures", []):
            mixture["components"][0].pop("calls")  # partly compared
    no_mixtures = copy.deepcopy(example)
    for read in no_mixtures["reads"]:
        if read.get("admitted"):
            read.pop("mixtures")
            break
    unnamed_repeat = copy.deepcopy(example)
    for zone in unnamed_repeat["zones"]:
        if zone["kind"] == "repeat":
            zone.pop("name")
            break
    named_clean = copy.deepcopy(example)
    for zone in named_clean["zones"]:
        if zone["kind"] == "clean":
            zone["name"] = "HVS-I"
            break
    inserted_evidence = copy.deepcopy(example)
    inserted_evidence["items"][0]["evidence"][0]["positions"][0]["ordinal"] = 1
    no_threshold = copy.deepcopy(example)
    for item in no_threshold["items"]:
        for entry in item["labels"]:
            if entry["label"] == "possible_heteroplasmy":
                entry["reads"][0].pop("threshold")
    observed_length = copy.deepcopy(example)
    for item in observed_length["items"]:
        for entry in item["labels"]:
            if entry["label"] == "estimated_length":
                entry["runs"][0]["length_evidence"] = "in_phase"
    no_rules = copy.deepcopy(example)
    no_rules["provenance"].pop("rules")
    return [
        ("review label outside the rules", unknown_label),
        ("failed review call without its error", failed_without_error),
        ("called review with an error", called_with_error),
        ("rejected review with items", rejected_with_items),
        ("review item call the document does not list", unlisted_call),
        ("review summary without a label count", missing_count),
        ("review repeat zone without its window name", unnamed_repeat),
        ("review clean zone with a window name", named_clean),
        ("admitted review read without mixtures", no_mixtures),
        ("unexplained mixture with compared components", compared_unexplained),
        ("mixture component compared without its calls", uncompared_explained),
        ("review evidence at an inserted column", inserted_evidence),
        ("possible heteroplasmy without its threshold", no_threshold),
        ("estimated length that a read shows", observed_length),
        ("review without its rules version", no_rules),
    ]


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
    parser.add_argument(
        "--variants",
        action="append",
        type=Path,
        default=[],
        metavar="RESULT",
        help="variants result to validate; may be repeated",
    )
    parser.add_argument(
        "--notation",
        action="append",
        type=Path,
        default=[],
        metavar="RESULT",
        help="notation result to validate; may be repeated",
    )
    parser.add_argument(
        "--consensus",
        action="append",
        type=Path,
        default=[],
        metavar="RESULT",
        help="consensus result to validate; may be repeated",
    )
    parser.add_argument(
        "--review",
        action="append",
        type=Path,
        default=[],
        metavar="RESULT",
        help="review result to validate; may be repeated",
    )
    parser.add_argument(
        "--signal",
        action="append",
        type=Path,
        default=[],
        metavar="RESULT",
        help="run signal result to validate; may be repeated",
    )
    parser.add_argument(
        "--run-context",
        action="append",
        type=Path,
        default=[],
        metavar="RESULT",
        help="run context result to validate; may be repeated",
    )
    parser.add_argument(
        "--evidence",
        action="append",
        type=Path,
        default=[],
        metavar="RESULT",
        help="read evidence result to validate; may be repeated",
    )
    parser.add_argument(
        "--eval-run",
        action="append",
        type=Path,
        default=[],
        metavar="DIRECTORY",
        help="dna-eval run output directory to validate; may be repeated",
    )
    return parser.parse_args(argv)


def validate_rejected_evidence(
    document: dict[str, Any], label: str, errors: list[str]
) -> None:
    """Check a rejected read's numbers: below its minimum, and callable and masked calls that add up."""
    rejection = document.get("rejection", {})
    if rejection.get("reason") == "too_short":
        if not rejection.get("calls", 0) < rejection.get("minimum_calls", 0):
            errors.append(f"{label}: a too-short read reaches its minimum calls")
        return
    callable_calls = rejection.get("callable_calls", 0)
    if not callable_calls < rejection.get("minimum_callable_calls", 0):
        errors.append(f"{label}: a read reaches its minimum callable calls")
    callability = document.get("callability", {})
    segments = callability.get("segments", [])
    spans = [
        (segment.get("calls", {}).get("start"), segment.get("calls", {}).get("end"))
        for segment in segments
    ]
    calls = spans[-1][1] if spans else 0
    if [start for start, _ in spans] != [0] + [end for _, end in spans[:-1]]:
        errors.append(f"{label}: callability segments do not partition the calls")
    if callable_calls + callability.get("masked_calls", 0) != calls:
        errors.append(f"{label}: callable and masked calls do not add up to the calls")


def validate_evidence(document: Any, label: str, errors: list[str]) -> None:
    """Check a read evidence's call order, intervals, warnings, and listed reasons, or a rejected read's numbers."""
    if not isinstance(document, dict):
        return
    if "rejection" in document:
        validate_rejected_evidence(document, label, errors)
        return
    read = document.get("read", {})
    calls = document.get("calls", [])
    if read.get("call_count") != len(calls):
        errors.append(f"{label}: read call_count disagrees with the calls")
    if [call.get("index") for call in calls] != list(range(len(calls))):
        errors.append(f"{label}: calls are not in call order")
    for name in ("informative", "trim"):
        interval = read.get(name, {})
        if not 0 <= interval.get("start", 0) <= interval.get("end", 0) <= len(calls):
            errors.append(f"{label}: {name} interval exceeds the calls")
    warnings = document.get("warnings", {})
    integrity = document.get("signal_quality", {}).get("integrity", {})
    unresolved = sum(1 for call in calls if call.get("base") == "N")
    if warnings.get("unresolved_primary_calls") != unresolved:
        errors.append(f"{label}: unresolved_primary_calls disagrees with the calls")
    if warnings.get("clipped_channel_samples") != integrity.get(
        "clipped_channel_samples"
    ):
        errors.append(f"{label}: clipped_channel_samples disagrees with the integrity")
    for index, mixture in enumerate(read.get("mixtures", [])):
        span = mixture.get("calls", {})
        offsets = [
            component.get("offset") for component in mixture.get("components", [])
        ]
        if not 0 <= span.get("start", 0) < span.get("end", 0) <= len(calls):
            errors.append(f"{label}: mixture {index} stretch exceeds the calls")
        if not offsets or offsets[0] != 0 or len(set(offsets)) != len(offsets):
            errors.append(
                f"{label}: mixture {index} does not lead with one dominant sequence"
            )
        if any(
            not span.get("start", 0) <= call < span.get("end", 0)
            for call in mixture.get("unexplained", [])
        ):
            errors.append(
                f"{label}: mixture {index} lists an unexplained call outside it"
            )
    for index, tract in enumerate(read.get("tracts", [])):
        validate_tract(tract, len(calls), f"{label}: tract {index}", errors)
    listed = set(read.get("reasons", []))
    for call in calls:
        held = set(call.get("identity", {}).get("withheld", [])) | set(
            call.get("indel", {}).get("withheld", [])
        )
        if not held <= listed:
            errors.append(
                f"{label}: call {call.get('index')} holds a reason read.reasons does not list"
            )


def validate_marginal(
    marginal: dict[str, Any], longest: int, label: str, errors: list[str]
) -> None:
    """Check a tract marginal as the kernel does (SRS-CALL-025)."""
    lengths = marginal.get("lengths", [])
    shares = marginal.get("shares", [])
    if lengths != sorted(set(lengths)) or (lengths and lengths[-1] > longest):
        errors.append(f"{label}: lengths are not ascending within the candidates")
    if len(shares) != len(lengths) or len(marginal.get("errors", [])) != len(lengths):
        errors.append(f"{label}: shares or errors do not match the lengths")
    if abs(sum(shares) - 1.0) > 1e-6:
        errors.append(f"{label}: shares do not sum to one")
    dominant, following = marginal.get("dominant"), marginal.get("next")
    if dominant not in lengths or (
        following is not None and (following not in lengths or following == dominant)
    ):
        errors.append(f"{label}: dominant or next length is not among the lengths")
    undefined = following is not None and marginal.get("z_adjusted") is None
    if (marginal.get("margin") == "undefined") != undefined:
        errors.append(f"{label}: margin is undefined unless two lengths lack one")


def validate_tract(
    tract: dict[str, Any], count: int, label: str, errors: list[str]
) -> None:
    """Check a tract decode as the kernel does (SRS-CALL-025)."""
    span = tract.get("calls", {})
    start, end = span.get("start", 0), span.get("end", 0)
    if not start < tract.get("onset", 0) < count or end > count:
        errors.append(f"{label}: does not fit the calls")
    if end - start < tract.get("extent", 0):
        errors.append(f"{label}: fits fewer calls than its extent")
    candidates = tract.get("candidates", [])
    marks = [candidate.get("mark") for candidate in candidates]
    order = [None] + [base for base in "ACGT" if base != tract.get("run")]
    if marks != [mark for mark in order if mark in marks] or marks[:1] != [None]:
        errors.append(f"{label}: candidates are not in evaluation order")
    chosen = [c for c in candidates if c.get("mark") == tract.get("mark")]
    if not chosen or chosen[0].get("residual") != tract.get("residual"):
        errors.append(f"{label}: the chosen structure is not a candidate")
    elif any(c.get("residual", 0) < tract.get("residual", 0) for c in candidates):
        errors.append(f"{label}: a candidate fits better than the chosen one")
    total = tract.get("total", {})
    lengths = total.get("lengths", [])
    if lengths != list(range(lengths[0], lengths[0] + len(lengths)) if lengths else []):
        errors.append(f"{label}: total lengths are not consecutive")
    longest = lengths[-1] if lengths else 0
    validate_marginal(total, longest, f"{label} total", errors)
    validate_marginal(tract.get("first", {}), longest, f"{label} first", errors)
    cleaned = tract.get("cleaned")
    if cleaned is not None:
        validate_marginal(cleaned.get("first", {}), longest, f"{label} cleaned", errors)


def rejected_evidence_shapes(
    example: dict[str, Any],
) -> list[tuple[str, dict[str, Any]]]:
    scored_with_reason = copy.deepcopy(example)
    scored_with_reason["calls"][0]["placement"]["reason"] = "weak_signal"
    masked_without_reason = copy.deepcopy(example)
    masked_without_reason["calls"][5]["placement"].pop("reason")
    unknown_cause = copy.deepcopy(example)
    unknown_cause["calls"][4]["signal"]["cause"] = "context"
    wide_share = copy.deepcopy(example)
    wide_share["calls"][0]["second"]["share"] = 0.6
    no_signal = copy.deepcopy(example)
    no_signal["calls"][0].pop("signal")
    no_trim = copy.deepcopy(example)
    no_trim["read"].pop("trim")
    no_signal_quality = copy.deepcopy(example)
    no_signal_quality.pop("signal_quality")
    no_warnings = copy.deepcopy(example)
    no_warnings.pop("warnings")
    three_peaks = copy.deepcopy(example)
    three_peaks["calls"][0]["signal"]["peaks"].pop("T")
    no_channels = copy.deepcopy(example)
    no_channels["calls"][0]["signal"].pop("channels")
    no_name = copy.deepcopy(example)
    no_name["provenance"]["input"].pop("name")
    path_name = copy.deepcopy(example)
    path_name["provenance"]["input"]["name"] = "runs/trace.ab1"
    read_and_rejection = copy.deepcopy(example)
    read_and_rejection["rejection"] = {
        "reason": "too_short",
        "calls": 5,
        "minimum_calls": 10,
    }
    rejected = load_json(REJECTED_EVIDENCE_EXAMPLE)
    short_with_callability = copy.deepcopy(rejected)
    short_with_callability["rejection"] = {
        "reason": "too_short",
        "calls": 5,
        "minimum_calls": 10,
    }
    callable_without_integrity = copy.deepcopy(rejected)
    callable_without_integrity.pop("integrity")
    placement_rejection = copy.deepcopy(rejected)
    placement_rejection["rejection"] = {"reason": "ambiguous_placement"}
    return [
        ("evidence scored call with a placement reason", scored_with_reason),
        ("evidence masked call without its reason", masked_without_reason),
        ("evidence second-peak cause the modality does not give", unknown_cause),
        ("evidence second share above one half", wide_share),
        ("evidence call without its signal", no_signal),
        ("evidence without the trim", no_trim),
        ("evidence without its signal quality", no_signal_quality),
        ("evidence without its warnings", no_warnings),
        ("evidence peaks without one channel", three_peaks),
        ("evidence call without its window peaks", no_channels),
        ("evidence without the trace's name", no_name),
        ("evidence naming its trace by a path", path_name),
        ("evidence of a read with a rejection", read_and_rejection),
        ("too-short evidence with callability", short_with_callability),
        ("rejected evidence without its integrity", callable_without_integrity),
        ("evidence rejected by the core", placement_rejection),
    ]


def inconsistent_tract_shapes(
    example: dict[str, Any],
) -> list[tuple[str, dict[str, Any]]]:
    """Read evidence whose tract the schema admits but the kernel refuses."""
    past_calls = copy.deepcopy(example)
    past_calls["read"]["tracts"][0]["calls"]["end"] = 99
    unordered = copy.deepcopy(example)
    unordered["read"]["tracts"][0]["candidates"].reverse()
    shares = copy.deepcopy(example)
    shares["read"]["tracts"][0]["first"]["shares"][0] = 0.5
    gap = copy.deepcopy(example)
    gap["read"]["tracts"][0]["total"]["lengths"] = [1, 2, 4]
    clear_without_margin = copy.deepcopy(example)
    clear_without_margin["read"]["tracts"][0]["first"]["z_adjusted"] = None
    return [
        ("evidence tract past the calls", past_calls),
        ("evidence tract candidates out of order", unordered),
        ("evidence tract shares not summing to one", shares),
        ("evidence tract total lengths with a gap", gap),
        ("evidence tract clear without a margin", clear_without_margin),
    ]


def validate_signal(document: Any, label: str, errors: list[str]) -> None:
    """Check a run signal's call total (three channels per call), key order, and share order."""
    if not isinstance(document, dict):
        return
    keys = document.get("keys", [])
    total = sum(count for key in keys for _, count in key.get("shares", []))
    calls = document.get("calls")
    if not isinstance(calls, int) or total != 3 * calls:
        errors.append(
            f"{label}: the keys' total {total} is not three times calls {calls}"
        )
    order = [(key.get("context"), key.get("base"), key.get("second")) for key in keys]
    if order != sorted(order) or len(set(order)) != len(order):
        errors.append(f"{label}: keys are not sorted and distinct")
    for key in keys:
        shares = [share for share, _ in key.get("shares", [])]
        if shares != sorted(set(shares)):
            errors.append(
                f"{label}: shares of {key.get('context')} are not sorted and distinct"
            )
        if key.get("context", "  ")[2:3] != key.get("base"):
            errors.append(
                f"{label}: base of {key.get('context')} is not the context's centre"
            )
    traces = document.get("traces", [])
    if traces != sorted(traces):
        errors.append(f"{label}: traces are not sorted")


def validate_run_context(document: Any, label: str, errors: list[str]) -> None:
    """Check a run context's read order, window order, edge order, and each run's histogram."""
    if not isinstance(document, dict):
        return
    reads = document.get("reads", [])
    if reads != sorted(set(reads)):
        errors.append(f"{label}: reads are not sorted and distinct")
    length = document.get("length", {})
    windows = length.get("windows", [])
    if any(first > last for first, last in windows) or any(
        left[1] >= right[0] for left, right in itertools.pairwise(windows)
    ):
        errors.append(f"{label}: length windows are not ordered and disjoint")
    edges = length.get("edges", [])
    if any(left >= right for left, right in itertools.pairwise(edges)):
        errors.append(f"{label}: length edges are not strictly ascending")
    runs = length.get("runs", [])
    if [run.get("run") for run in runs] != [1, 2]:
        errors.append(f"{label}: length runs are not run 1 then run 2")
    for run in runs:
        histogram = run.get("histogram", [])
        if len(histogram) != len(edges) + 1:
            errors.append(
                f"{label}: run {run.get('run')} histogram has {len(histogram)} bins, "
                f"expected {len(edges) + 1}"
            )
        if sum(histogram) != run.get("pairs"):
            errors.append(
                f"{label}: run {run.get('run')} histogram does not sum to its pairs"
            )


def rejected_run_context_shapes(
    example: dict[str, Any],
) -> list[tuple[str, dict[str, Any]]]:
    other_method = copy.deepcopy(example)
    other_method["length"]["method"] = 2
    noisy = copy.deepcopy(example)
    noisy["length"]["noise_forward"] = [0.6, 0.1]
    one_run = copy.deepcopy(example)
    one_run["length"]["runs"].pop()
    no_reads = copy.deepcopy(example)
    no_reads.pop("reads")
    bad_use = copy.deepcopy(example)
    bad_use["extra"] = True
    return [
        ("run context of another length method", other_method),
        ("run context with noise above one half", noisy),
        ("run context with one run", one_run),
        ("run context without its reads", no_reads),
        ("run context with an unknown field", bad_use),
    ]


def read_table(text: str) -> dict[str, Any]:
    """Read one `dna-eval` TSV table as `dna.eval_run/v1` describes it.

    Raises `ValueError` when the text is not a table: no header, a line without
    its newline, a carriage return, or a row whose cells do not match the header.
    """
    if "\r" in text:
        raise ValueError("holds a carriage return")
    if not text.endswith("\n"):
        raise ValueError("does not end in a newline")
    lines = text[:-1].split("\n")
    comments: list[list[str]] = []
    while lines and lines[0].startswith("#"):
        line = lines.pop(0)
        if not line.startswith("# "):
            raise ValueError("a comment line does not start with '# '")
        comments.append(line[2:].split("\t"))
    if not lines or not lines[0]:
        raise ValueError("has no header")
    columns = lines[0].split("\t")
    if len(set(columns)) != len(columns):
        raise ValueError("repeats a column")
    rows = []
    for number, line in enumerate(lines[1:], start=len(comments) + 2):
        cells = line.split("\t")
        if len(cells) != len(columns):
            raise ValueError(
                f"line {number} has {len(cells)} cells, the header {len(columns)}"
            )
        rows.append(dict(zip(columns, cells, strict=True)))
    return {"comments": comments, "columns": columns, "rows": rows}


def read_run(directory: Path) -> dict[str, Any]:
    """Read a `dna-eval run` output directory as one `dna.eval_run/v1` instance."""
    document: dict[str, Any] = {}
    for path in sorted(directory.iterdir()):
        if not path.is_file():
            raise ValueError(f"{path.name} is not a file")
        text = path.read_text(encoding="utf-8")
        if path.suffix == ".json":
            document[path.name] = json.loads(text)
        elif path.suffix == ".tsv":
            try:
                document[path.name] = read_table(text)
            except ValueError as exc:
                raise ValueError(f"{path.name} {exc}") from exc
        else:
            document[path.name] = text
    return document


def check_eval_tables(schema: dict[str, Any], errors: list[str]) -> None:
    """Check that every table of the evaluation run schema lists its columns in row order."""
    for name, definition in schema["$defs"].items():
        properties = definition.get("properties", {})
        if "columns" not in properties:
            continue
        columns = properties["columns"]
        fixed = columns.get("const") or [
            item["const"] for item in columns.get("prefixItems", [])
        ]
        row = properties["rows"]["items"]
        if list(row["properties"]) != fixed or row["required"] != fixed:
            errors.append(f"eval run table {name}: columns disagree with its rows")
        if ("items" in columns) != isinstance(row["additionalProperties"], dict):
            errors.append(
                f"eval run table {name}: open columns disagree with its open rows"
            )


# The measures of `run.json` `samples`, by their status column in `scores.tsv`.
EVAL_MEASURES = {
    "primary": "status",
    "strict": "strict_status",
    "iupac_compatible": "iupac_compatible_status",
    "haplotype": "haplotype_status",
    "haplotype_iupac_compatible": "haplotype_iupac_compatible_status",
    "haplotype_lod": "haplotype_lod_status",
}


def validate_eval_run(document: Any, label: str, errors: list[str]) -> None:
    """Check an evaluation run's sample order and that `run.json` counts what the tables hold."""
    if not isinstance(document, dict):
        return
    header = document.get("run.json", {})
    scores = document.get("scores.tsv", {}).get("rows", [])
    samples = [row.get("sample_id") for row in scores]
    if samples != sorted(set(samples)):
        errors.append(f"{label}: scores.tsv samples are not sorted and distinct")
    for measure, column in EVAL_MEASURES.items():
        counts = header.get("samples", {}).get(measure, {})
        statuses = [row.get(column) for row in scores]
        found = {status: statuses.count(status) for status in set(statuses)}
        expected = {key: value for key, value in counts.items() if key != "total"}
        if counts.get("total") != len(scores) or any(
            found.get(status, 0) != value for status, value in expected.items()
        ):
            errors.append(
                f"{label}: run.json samples.{measure} disagrees with scores.tsv"
            )
    changes = document.get("truth_changes.tsv")
    if changes is not None and len(changes["rows"]) != header.get("truth", {}).get(
        "changed_samples"
    ):
        errors.append(f"{label}: truth_changes.tsv disagrees with changed_samples")
    validate_eval_passes(
        header.get("workload", {}), "recheck" in header, scores, label, errors
    )
    triage = header.get("triage")
    if isinstance(triage, dict):
        policies = triage.get("policies")
        columns = document.get("triage_policies.tsv", {}).get("columns", [])
        sweep = document.get("triage_sweep.tsv", {}).get("rows", [])
        if len(columns) - 3 != policies or len(sweep) != policies:
            errors.append(
                f"{label}: triage policies disagree with the sweep and policy tables"
            )
        names = [row.get("policy") for row in sweep]
        if names != columns[3:] or not triage.get("policy", "").startswith(
            f"{names[0] if names else ''}: "
        ):
            errors.append(
                f"{label}: the sweep's policies are not the policy table's, reported first"
            )


# The classes that split the passed units with an unreviewed block.
EVAL_ESCAPED_SPLIT = ("non_het", "truth_n_only", "het_only", "recheck_only")


def validate_eval_passes(
    workload: Any,
    recheck: bool,
    scores: list[dict[str, Any]],
    label: str,
    errors: list[str],
) -> None:
    """Check that `run.json` `workload` counts what `scores.tsv`'s review columns hold.

    Units split into `EVAL_ESCAPED_SPLIT`, which sum to `total`; without a
    re-check list (`recheck` false) no unit is re-check only.
    """
    if not isinstance(workload, dict):
        return
    escaped_keys = ("total", *EVAL_ESCAPED_SPLIT)
    for name, passes in workload.items():
        escaped = [passes.get("auto_pass", {}).get("unreviewed", {})]
        for key, units_key in (("regions", "by_region"), ("zones", "by_zone")):
            group = passes.get(key, {})
            units = group.get(units_key, [])
            escaped += [group.get("unreviewed", {})]
            escaped += [unit.get("unreviewed", {}) for unit in units]
            summed = [
                sum(unit.get(field, 0) for unit in units)
                for field in ("total", "checked")
            ]
            summed += [
                sum(unit.get("unreviewed", {}).get(field, 0) for unit in units)
                for field in escaped_keys
            ]
            recorded = [group.get("total"), group.get("checked")]
            recorded += [
                group.get("unreviewed", {}).get(field) for field in escaped_keys
            ]
            if summed != recorded:
                errors.append(
                    f"{label}: run.json workload.{name} {key} disagree with {units_key}"
                )
        if any(
            sum(item.get(key, 0) for key in EVAL_ESCAPED_SPLIT) != item.get("total", 0)
            for item in escaped
        ):
            errors.append(
                f"{label}: run.json workload.{name} splits unreviewed units inconsistently"
            )
        if not recheck and any(item.get("recheck_only", 0) for item in escaped):
            errors.append(
                f"{label}: run.json workload.{name} counts re-check units without a list"
            )
    passes = workload.get("all")
    if not isinstance(passes, dict):
        return
    counted = [
        row
        for row in scores
        if row.get("haplotype_lod_status") in ("concordant", "discordant")
        and row.get("review_clear")
    ]
    clear = [row for row in counted if row.get("review_clear") == "Y"]
    zones = [set(row.get("review_zones", "").split()) for row in counted]
    regions = [{zone.split("@", 1)[1] for zone in found} for found in zones]
    by_region = passes.get("regions", {}).get("by_region", [])
    by_zone = passes.get("zones", {}).get("by_zone", [])
    expected = (
        len(counted),
        len(clear),
        sum(1 for row in clear if row.get("unreviewed")),
        sum(len(found) for found in regions),
        sum(int(row.get("review_zone_count") or 0) for row in counted),
        [
            sum(1 for found in regions if unit.get("name") in found)
            for unit in by_region
        ],
        [sum(1 for found in zones if unit.get("name") in found) for unit in by_zone],
    )
    found = (
        passes.get("samples"),
        passes.get("auto_pass", {}).get("samples"),
        passes.get("auto_pass", {}).get("unreviewed", {}).get("total"),
        passes.get("regions", {}).get("checked"),
        passes.get("zones", {}).get("checked"),
        [unit.get("checked") for unit in by_region],
        [unit.get("checked") for unit in by_zone],
    )
    if expected != found:
        errors.append(f"{label}: run.json workload.all disagrees with scores.tsv")


def rejected_eval_run_shapes(
    example: dict[str, Any],
) -> list[tuple[str, dict[str, Any]]]:
    without_scores = copy.deepcopy(example)
    without_scores.pop("scores.tsv")
    unknown_file = copy.deepcopy(example)
    unknown_file["notes.tsv"] = copy.deepcopy(example["checks.tsv"])
    renamed_column = copy.deepcopy(example)
    renamed_column["scores.tsv"]["columns"][14] = "rules"
    unknown_status = copy.deepcopy(example)
    unknown_status["scores.tsv"]["rows"][0]["status"] = "exact"
    row_without_cell = copy.deepcopy(example)
    row_without_cell["windows.tsv"]["rows"][0].pop("unreviewed_upper95")
    unformatted_rate = copy.deepcopy(example)
    unformatted_rate["estimates.tsv"]["rows"][0]["rate"] = "0.5"
    diff_without_baseline = copy.deepcopy(example)
    diff_without_baseline["run.json"]["baseline"] = None
    baseline_without_diff = copy.deepcopy(example)
    baseline_without_diff.pop("diff.tsv")
    unchanged_with_changes = copy.deepcopy(example)
    unchanged_with_changes["run.json"]["truth"]["state"] = "unchanged"
    triage_without_file = copy.deepcopy(example)
    triage_without_file.pop("triage_leaks.tsv")
    triage_files_without_triage = copy.deepcopy(example)
    triage_files_without_triage["run.json"].pop("triage")
    signals_without_run_signal = copy.deepcopy(example)
    signals_without_run_signal["run.json"].pop("run_signal")
    context_table_without_context = copy.deepcopy(example)
    context_table_without_context["run_context.tsv"] = {
        "comments": [],
        "columns": [
            "month",
            "samples",
            "analysed",
            "with_context",
            "run1_calibrated",
            "run2_calibrated",
            "batches",
            "run1_pairs",
            "run2_pairs",
        ],
        "rows": [
            {
                "month": "all",
                "samples": "1",
                "analysed": "1",
                "with_context": "0",
                "run1_calibrated": "0",
                "run2_calibrated": "0",
                "batches": "0",
                "run1_pairs": "0",
                "run2_pairs": "0",
            }
        ],
    }
    run_without_constants = copy.deepcopy(example)
    run_without_constants["run.json"].pop("constants")
    run_without_workload = copy.deepcopy(example)
    run_without_workload["run.json"].pop("workload")
    current_without_passes = copy.deepcopy(example)
    current_without_passes["run.json"]["workload"].pop("current")
    zones_without_rows = copy.deepcopy(example)
    zones_without_rows["run.json"]["workload"]["all"]["zones"].pop("by_zone")
    zone_without_unreviewed = copy.deepcopy(example)
    zone_without_unreviewed["run.json"]["workload"]["all"]["zones"]["by_zone"][0].pop(
        "unreviewed"
    )
    split_without_truth_n = copy.deepcopy(example)
    split_without_truth_n["run.json"]["workload"]["all"]["auto_pass"]["unreviewed"].pop(
        "truth_n_only"
    )
    recheck_without_sha256 = copy.deepcopy(example)
    recheck_without_sha256["run.json"]["recheck"].pop("sha256")
    recheck_with_bad_rows = copy.deepcopy(example)
    recheck_with_bad_rows["run.json"]["recheck"]["rows"] = -1
    forms_without_database = copy.deepcopy(example)
    forms_without_database["forms.tsv"]["comments"] = []
    commented_scores = copy.deepcopy(example)
    commented_scores["scores.tsv"]["comments"] = [["note"]]
    unknown_policy_class = copy.deepcopy(example)
    unknown_policy_class["triage_policies.tsv"]["rows"][0]["draft"] = "approved"
    return [
        ("eval run without scores.tsv", without_scores),
        ("eval run with a file the contract does not name", unknown_file),
        ("eval scores with a renamed column", renamed_column),
        ("eval scores with a status the harness does not write", unknown_status),
        ("eval windows row without a cell", row_without_cell),
        ("eval estimate rate without four decimals", unformatted_rate),
        ("eval diff.tsv without a baseline", diff_without_baseline),
        ("eval baseline without diff.tsv", baseline_without_diff),
        ("eval truth_changes.tsv with an unchanged truth", unchanged_with_changes),
        ("eval triage without one of its files", triage_without_file),
        ("eval triage files without the triage header", triage_files_without_triage),
        ("eval signals.tsv without the run signal", signals_without_run_signal),
        ("eval run_context.tsv without run contexts", context_table_without_context),
        ("eval run.json without its constants", run_without_constants),
        ("eval run.json without its workload", run_without_workload),
        ("eval current data without its passes", current_without_passes),
        ("eval zone passes without their zones", zones_without_rows),
        ("eval zone row without its unreviewed counts", zone_without_unreviewed),
        ("eval unreviewed split without truth_n_only", split_without_truth_n),
        ("eval re-check list without its sha256", recheck_without_sha256),
        ("eval re-check list with negative rows", recheck_with_bad_rows),
        ("eval forms.tsv without its database line", forms_without_database),
        ("eval scores.tsv with a comment line", commented_scores),
        ("eval policy table with an unknown class", unknown_policy_class),
    ]


def validate_eval_run_reader(example: dict[str, Any], errors: list[str]) -> None:
    """Write the example run as files, read it back, and reject malformed tables."""
    with tempfile.TemporaryDirectory() as temporary:
        directory = Path(temporary)
        for name, value in example.items():
            path = directory / name
            if name.endswith(".json"):
                path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
            elif name.endswith(".tsv"):
                lines = ["# " + "\t".join(comment) for comment in value["comments"]] + [
                    "\t".join(value["columns"])
                ]
                lines += [
                    "\t".join(row[column] for column in value["columns"])
                    for row in value["rows"]
                ]
                path.write_text("\n".join(lines) + "\n", encoding="utf-8")
            else:
                path.write_text(value, encoding="utf-8")
        if read_run(directory) != example:
            errors.append("the example eval run does not read back from its files")
    for text, description in [
        ("a\tb\n1\n", "a row without a cell"),
        ("a\tb\r\n1\t2\r\n", "carriage returns"),
        ("a\tb\n1\t2", "a last line without its newline"),
        ("a\ta\n1\t2\n", "a repeated column"),
        ("#note\na\n1\n", "a comment without '# '"),
        ("", "no header"),
    ]:
        try:
            read_table(text)
        except ValueError:
            continue
        errors.append(f"expected an eval table with {description} to be rejected")


# SRS-OUT-002
def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    errors: list[str] = []
    analysis_validator = validator(ANALYSIS_SCHEMA, errors)
    basecall_validator = validator(BASECALL_SCHEMA, errors)
    sample_validator = validator(SAMPLE_SCHEMA, errors)
    variants_validator = validator(VARIANTS_SCHEMA, errors)
    notation_validator = validator(NOTATION_SCHEMA, errors)
    consensus_validator = validator(CONSENSUS_SCHEMA, errors)
    review_validator = validator(REVIEW_SCHEMA, errors)
    signal_validator = validator(SIGNAL_SCHEMA, errors)
    run_context_validator = validator(RUN_CONTEXT_SCHEMA, errors)
    evidence_validator = validator(EVIDENCE_SCHEMA, errors)
    eval_run_validator = validator(EVAL_RUN_SCHEMA, errors)
    placed_validator = validator(PLACED_SCHEMA, errors)
    analysis_paths = args.analysis or [ANALYSIS_EXAMPLE]
    basecall_paths = args.basecalls or [BASECALL_EXAMPLE]
    sample_paths = args.sample_evidence or [SAMPLE_EXAMPLE]
    variants_paths = args.variants or [VARIANTS_EXAMPLE]
    notation_paths = args.notation or [NOTATION_EXAMPLE]
    consensus_paths = args.consensus or [CONSENSUS_EXAMPLE]
    review_paths = args.review or [REVIEW_EXAMPLE]
    signal_paths = args.signal or [SIGNAL_EXAMPLE]
    run_context_paths = args.run_context or [RUN_CONTEXT_EXAMPLE]
    evidence_paths = args.evidence or [EVIDENCE_EXAMPLE, REJECTED_EVIDENCE_EXAMPLE]
    validate_documents(analysis_validator, analysis_paths, errors)
    validate_documents(basecall_validator, basecall_paths, errors)
    validate_documents(sample_validator, sample_paths, errors)
    validate_documents(variants_validator, variants_paths, errors)
    validate_documents(notation_validator, notation_paths, errors)
    validate_documents(consensus_validator, consensus_paths, errors)
    validate_documents(review_validator, review_paths, errors)
    validate_documents(signal_validator, signal_paths, errors)
    validate_documents(run_context_validator, run_context_paths, errors)
    validate_documents(evidence_validator, evidence_paths, errors)
    validate_documents(placed_validator, [PLACED_EXAMPLE], errors)
    for path in evidence_paths:
        validate_evidence(load_json(path), str(path), errors)
    eval_runs: list[tuple[str, Any]] = []
    if not args.eval_run:
        eval_runs.append((str(EVAL_RUN_EXAMPLE), load_json(EVAL_RUN_EXAMPLE)))
    for directory in args.eval_run:
        try:
            eval_runs.append((str(directory), read_run(directory)))
        except (OSError, ValueError) as exc:
            errors.append(f"{directory}: not an evaluation run: {exc}")
    for label, document in eval_runs:
        try:
            eval_run_validator.validate(document)
        except ValidationError as exc:
            errors.append(f"{label}: invalid: {exc.message} at {exc.json_path}")
        validate_eval_run(document, label, errors)
    for path in signal_paths:
        validate_signal(load_json(path), str(path), errors)
    for path in run_context_paths:
        validate_run_context(load_json(path), str(path), errors)
    for path in review_paths:
        validate_review(load_json(path), str(path), errors)
    errors.extend(review_rules_errors())
    for path in consensus_paths:
        document = load_json(path)
        if isinstance(document, dict):
            validate_consensus_rejections(document, str(path), errors)
    validate_sample_support_topology(sample_paths, errors)
    validate_read_callability(analysis_paths + basecall_paths + sample_paths, errors)
    for path in analysis_paths:
        validate_analysis_calls(load_json(path), str(path), errors)

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

    wrong_opposition = copy.deepcopy(sample_example)
    wrong_opposition["variants"][0]["opposition"] = {
        "reads": [],
        "forward_reads": 0,
        "reverse_reads": 0,
    }
    opposition_errors: list[str] = []
    validate_sample_support_topology_document(
        wrong_opposition, "synthetic opposition", opposition_errors
    )
    if not opposition_errors:
        errors.append(
            "expected an opposition that omits a callable read to be rejected"
        )

    rejected_index = {
        read["rejection"]["reason"]: index
        for index, read in enumerate(sample_example["rejected_reads"])
    }
    for description, reason, mutate in [
        (
            "admitted name",
            "callable_calls_below_minimum",
            lambda read: read.update(name=sample_example["reads"][0]["name"]),
        ),
        (
            "callable count",
            "callable_calls_below_minimum",
            lambda read: read["rejection"].update(callable_calls=3),
        ),
        (
            "minimum",
            "callable_calls_below_minimum",
            lambda read: read["rejection"].update(minimum_callable_calls=0),
        ),
        (
            "too-short call count",
            "too_short",
            lambda read: read["rejection"].update(calls=10),
        ),
    ]:
        inconsistent_rejection = copy.deepcopy(sample_example)
        mutate(inconsistent_rejection["rejected_reads"][rejected_index[reason]])
        rejection_errors: list[str] = []
        validate_sample_support_topology_document(
            inconsistent_rejection, "synthetic rejected read", rejection_errors
        )
        if not rejection_errors:
            errors.append(
                f"expected a rejected read with inconsistent {description} to fail"
            )

    short_consensus = load_json(CONSENSUS_EXAMPLE)
    for read in short_consensus["rejected_reads"]:
        if read["rejection"]["reason"] == "too_short":
            read["rejection"]["calls"] = read["rejection"]["minimum_calls"]
    short_consensus_errors: list[str] = []
    validate_consensus_rejections(
        short_consensus, "synthetic consensus too-short read", short_consensus_errors
    )
    if not short_consensus_errors:
        errors.append(
            "expected a consensus too-short read reaching its minimum to fail"
        )

    inconsistent_callability = copy.deepcopy(analysis_example)
    inconsistent_callability["read"]["callability"]["masked_calls"] += 1
    callability_errors: list[str] = []
    validate_callability_view(
        inconsistent_callability["read"]["callability"],
        inconsistent_callability["read"]["call_count"],
        "synthetic inconsistent callability",
        callability_errors,
    )
    if not callability_errors:
        errors.append("expected inconsistent callability counts to be rejected")

    for description, mutate in [
        ("unordered", lambda calls: calls.reverse()),
        ("repeated", lambda calls: calls.append(dict(calls[-1]))),
        ("beyond the read", lambda calls: calls[-1].update(index=10_000)),
    ]:
        mutated = copy.deepcopy(analysis_example)
        mutate(mutated["calls"])
        call_errors: list[str] = []
        validate_analysis_calls(mutated, f"synthetic {description} calls", call_errors)
        if not call_errors:
            errors.append(f"expected {description} aligned calls to be rejected")

    outside_callable = copy.deepcopy(analysis_example)
    outside_callable["alignment"]["callable_reference_segments"] = [
        {"start": 0, "end": 1_000_000}
    ]
    outside_errors: list[str] = []
    validate_callable_segments(
        outside_callable["alignment"], "synthetic callable segment", outside_errors
    )
    if not outside_errors:
        errors.append("expected a callable segment outside the mapping to be rejected")

    unsorted_offsets = copy.deepcopy(analysis_example)
    for segment in unsorted_offsets["read"]["callability"]["segments"]:
        if segment["state"] == "dephased":
            segment["shadow_offsets"] = [1, -1]
    offset_errors: list[str] = []
    validate_callability_view(
        unsorted_offsets["read"]["callability"],
        unsorted_offsets["read"]["call_count"],
        "synthetic unsorted shadow offsets",
        offset_errors,
    )
    if not offset_errors:
        errors.append("expected unsorted shadow offsets to be rejected")

    assert_rejected(analysis_validator, rejected_analysis, errors)
    assert_rejected(basecall_validator, rejected_basecalls, errors)
    assert_rejected(sample_validator, rejected_samples, errors)
    rejected_variants = rejected_variants_shapes(load_json(VARIANTS_EXAMPLE))
    assert_rejected(variants_validator, rejected_variants, errors)
    rejected_notation = rejected_notation_shapes(load_json(NOTATION_EXAMPLE))
    assert_rejected(notation_validator, rejected_notation, errors)
    rejected_consensus = rejected_consensus_shapes(load_json(CONSENSUS_EXAMPLE))
    assert_rejected(consensus_validator, rejected_consensus, errors)
    rejected_run_context = rejected_run_context_shapes(load_json(RUN_CONTEXT_EXAMPLE))
    assert_rejected(run_context_validator, rejected_run_context, errors)
    review_example = load_json(REVIEW_EXAMPLE)
    rejected_review = rejected_review_shapes(review_example)
    assert_rejected(review_validator, rejected_review, errors)
    evidence_example = load_json(EVIDENCE_EXAMPLE)
    rejected_evidence = rejected_evidence_shapes(evidence_example)
    assert_rejected(evidence_validator, rejected_evidence, errors)
    for label, document in inconsistent_tract_shapes(evidence_example):
        tract_errors: list[str] = []
        validate_evidence(document, label, tract_errors)
        if not tract_errors:
            errors.append(f"expected {label} to be rejected, but it validated")
    check_eval_tables(load_json(EVAL_RUN_SCHEMA), errors)
    eval_run_example = load_json(EVAL_RUN_EXAMPLE)
    rejected_eval_runs = rejected_eval_run_shapes(eval_run_example)
    assert_rejected(eval_run_validator, rejected_eval_runs, errors)
    validate_eval_run_reader(eval_run_example, errors)
    miscounted_run = copy.deepcopy(eval_run_example)
    miscounted_run["run.json"]["samples"]["primary"]["failed"] += 1
    run_errors: list[str] = []
    validate_eval_run(miscounted_run, "synthetic miscounted run", run_errors)
    if not run_errors:
        errors.append(
            "expected run.json counts that disagree with scores.tsv to be rejected"
        )
    miscounted_passes = copy.deepcopy(eval_run_example)
    miscounted_passes["run.json"]["workload"]["all"]["auto_pass"]["samples"] += 1
    run_errors = []
    validate_eval_run(miscounted_passes, "synthetic miscounted passes", run_errors)
    if not run_errors:
        errors.append(
            "expected run.json workload counts that disagree with scores.tsv to be rejected"
        )
    miscounted_zone = copy.deepcopy(eval_run_example)
    zone = miscounted_zone["run.json"]["workload"]["all"]["zones"]["by_zone"][0]
    zone["unreviewed"] = {
        "total": 1,
        "non_het": 1,
        "het_only": 0,
        "truth_n_only": 0,
        "recheck_only": 0,
    }
    run_errors = []
    validate_eval_run(miscounted_zone, "synthetic miscounted zone", run_errors)
    if not run_errors:
        errors.append(
            "expected zone rows that disagree with their totals to be rejected"
        )
    for item, description in [
        ({"total": 1, "non_het": 1, "het_only": 1}, "a split that exceeds its total"),
        (
            {"total": 1, "truth_n_only": 0, "recheck_only": 0},
            "a split short of its total",
        ),
    ]:
        missplit = copy.deepcopy(eval_run_example)
        missplit["run.json"]["workload"]["all"]["auto_pass"]["unreviewed"] = item
        run_errors = []
        validate_eval_run(missplit, f"synthetic {description}", run_errors)
        if not any("splits unreviewed units" in error for error in run_errors):
            errors.append(f"expected passes with {description} to be rejected")
    unlisted_recheck = copy.deepcopy(eval_run_example)
    unlisted_recheck["run.json"].pop("recheck")
    for passes in unlisted_recheck["run.json"]["workload"].values():
        unreviewed = passes["auto_pass"]["unreviewed"]
        if unreviewed["non_het"]:
            unreviewed["non_het"] -= 1
            unreviewed["recheck_only"] += 1
    run_errors = []
    validate_eval_run(unlisted_recheck, "synthetic unlisted re-check", run_errors)
    if not any("re-check units without a list" in error for error in run_errors):
        errors.append("expected re-check units without a re-check list to be rejected")
    unsorted_run = copy.deepcopy(eval_run_example)
    unsorted_run["scores.tsv"]["rows"].reverse()
    run_errors = []
    validate_eval_run(unsorted_run, "synthetic unsorted scores", run_errors)
    if not run_errors:
        errors.append("expected scores.tsv out of sample order to be rejected")
    unlisted_reason = copy.deepcopy(evidence_example)
    unlisted_reason["calls"][3]["identity"]["withheld"].append("unlisted_reason")
    evidence_errors: list[str] = []
    validate_evidence(unlisted_reason, "synthetic unlisted reason", evidence_errors)
    if not evidence_errors:
        errors.append("expected a reason read.reasons does not list to be rejected")
    miscounted_warnings = copy.deepcopy(evidence_example)
    miscounted_warnings["warnings"]["unresolved_primary_calls"] += 1
    warning_errors: list[str] = []
    validate_evidence(
        miscounted_warnings, "synthetic miscounted warnings", warning_errors
    )
    if not warning_errors:
        errors.append(
            "expected evidence warnings that disagree with the calls to be rejected"
        )
    uncounted = load_json(REJECTED_EVIDENCE_EXAMPLE)
    uncounted["rejection"]["callable_calls"] += 1
    uncounted_errors: list[str] = []
    validate_evidence(uncounted, "synthetic uncounted rejection", uncounted_errors)
    if not uncounted_errors:
        errors.append(
            "expected rejected evidence whose calls do not add up to be rejected"
        )
    unordered_review = copy.deepcopy(review_example)
    unordered_review["items"].reverse()
    review_errors: list[str] = []
    validate_review(unordered_review, "synthetic unordered review", review_errors)
    if not review_errors:
        errors.append("expected review items out of priority order to be rejected")
    miscounted = copy.deepcopy(review_example)
    miscounted["reference_only"].append(miscounted["reference_only"][-1] + 1)
    count_errors: list[str] = []
    validate_review(miscounted, "synthetic miscounted reference_only", count_errors)
    if not count_errors:
        errors.append("expected a miscounted reference_only to be rejected")
    minor_first = copy.deepcopy(review_example)
    for read in minor_first["reads"]:
        for mixture in read.get("mixtures", []):
            mixture["components"].reverse()
    mixture_errors: list[str] = []
    validate_review(minor_first, "synthetic minor-first mixture", mixture_errors)
    if not mixture_errors:
        errors.append("expected a mixture led by a minor sequence to be rejected")
    # An explained possible heteroplasmy without a difference is a note at
    # the context rank, and leaves the summary's checks.
    noted = copy.deepcopy(review_example)
    note = next(
        item
        for item in noted["items"]
        if [entry["label"] for entry in item["labels"]] == ["possible_heteroplasmy"]
        and not item["calls"]
    )
    note["labels"][0]["explained"] = True
    note["kind"] = "note"
    noted["items"].remove(note)
    table = review_label_table(REVIEW_SCHEMA)
    context = max(rank for rank, _ in table.values())
    at = next(
        index
        for index, item in enumerate(noted["items"])
        if min(
            (table[entry["label"]][0] for entry in item["labels"]),
            default=context + 1,
        )
        >= context
        and item["first"]["position"] > note["first"]["position"]
    )
    noted["items"].insert(at, note)
    noted["summary"]["checks"] -= 1
    noted["summary"]["positions"]["check"] -= (
        note["last"]["position"] - note["first"]["position"] + 1
    )
    for zone in noted["zones"]:
        if zone["first"] <= note["first"]["position"] <= zone["last"]:
            zone["checks"] -= 1
            noted["summary"]["zones"]["check"] -= zone["checks"] == 0
    note_errors: list[str] = []
    validate_review(noted, "synthetic note", note_errors)
    if note_errors:
        errors.extend(note_errors)
    note["kind"] = "check"
    still_checked: list[str] = []
    validate_review(noted, "synthetic explained check", still_checked)
    if not any("kind disagrees" in error for error in still_checked):
        errors.append(
            "expected an explained heteroplasmy listed as a check to be rejected"
        )
    gapped_zones = copy.deepcopy(review_example)
    gapped_zones["zones"].pop()
    zone_errors: list[str] = []
    validate_review(gapped_zones, "synthetic gapped zones", zone_errors)
    if not zone_errors:
        errors.append("expected zones that leave positions out to be rejected")

    for error in errors:
        print(f"FAIL: {error}", file=sys.stderr)
    if errors:
        print(f"{len(errors)} check(s) failed", file=sys.stderr)
        return 1
    print(
        f"OK: validated {len(analysis_paths)} analysis, {len(basecall_paths)} basecall, "
        f"{len(sample_paths)} sample-evidence, {len(variants_paths)} variants, "
        f"{len(notation_paths)} notation, {len(consensus_paths)} consensus, "
        f"{len(review_paths)} review, {len(signal_paths)} run signal, "
        f"{len(run_context_paths)} run context, "
        f"{len(evidence_paths)} read evidence document(s), 1 placed read, "
        f"{len(eval_runs)} evaluation run(s), "
        f"{len(valid_shapes)} analysis "
        "call shapes; rejected "
        f"{len(rejected_analysis) + len(rejected_basecalls) + len(rejected_samples) + len(rejected_variants) + len(rejected_notation) + len(rejected_consensus) + len(rejected_run_context) + len(rejected_review) + len(rejected_evidence) + len(rejected_eval_runs)} "
        "invalid shape(s)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
