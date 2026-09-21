# `src/model/read_observation.rs`

## Purpose

Defines the complete scientific observation produced from one independently
processed read after reference placement.

## Responsibilities

- Group input, reference, and configuration identities with signal-derived calls, basecall-independent locus/signal observations, quality-control result, selected alignment, and read-level variants for one read.
- Establish the explicit boundary between one-read processing and sample-level reconciliation.
- Preserve evidence-derived orientation and mapped reference segments from the
  selected alignment.
- Retain the UTF-8 source basename as reviewer-facing provenance while keeping assay names, primer labels, filename semantics, declared direction, and canonical F/R pairing outside the scientific placement model.

## Non-responsibilities

No sample aggregation, consensus, cross-read conflict resolution, amplicon
classification, metadata-driven placement, filesystem access, or reporting
policy.

## Key type

- `ReadObservation`: immutable ownership bundle for one independently processed read's scientific products, stable input identity, and non-authoritative source basename provenance.

## Invariants

- A read is processed independently before any sample-level reconciliation.
- Covered reference region and orientation are consequences of alignment
  evidence.
- The source basename may be displayed for reviewer traceability but does not alter scientific placement, overlap, or reconciliation. Optional assay metadata may be used later for provenance or post-mapping QC under the same rule.
- Missing canonical F/R partners do not invalidate an otherwise usable read.

## Dependencies

Validated model types only: base calls, signal analysis, quality control,
alignment and variant-calling result.

## Traceability

ADR-0023, ADR-0028, and ADR-0029; `INV-EVID-004` through `INV-EVID-006`; `INV-READ-001` through `INV-READ-004`; `SRS-ALN-007` through `SRS-ALN-011`.

## Status

Implemented as the one-read boundary. Sample-level aggregation is intentionally
outside this module.
