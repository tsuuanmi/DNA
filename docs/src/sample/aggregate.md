# `src/sample/aggregate.rs`

## Purpose

Owns sample-level validation, deterministic read ordering, and assembly of compact `SampleEvidence`.

## Responsibilities

- Require at least one independently completed read.
- Require identical reference/configuration identities.
- Reject duplicate trace content by SHA-256 even when filenames differ.
- Sort reads by SHA-256 so internal registry indexes are deterministic and independent of CLI order; expose the same validated ordering internally for dense validation-locus aggregation.
- Build the one top-level read registry with source basename, stable identity, trace-integrity evidence, and selected post-trim alignment summary.
- Consume the validated `SampleReconciliationConfig`.
- Derive run-length total/forward/reverse coverage topology from the read registry.
- Delegate pairwise overlap admission, sparse production locus selection through the shared locus builder, and normalized variant aggregation/support-topology/profile derivation.

## Non-responsibilities

No filename-driven placement, pair-first merging, JSON filename-stem projection, consensus, or sample-level conflict adjudication.

## Invariants

Internal indexes are deterministic implementation references only. Public reviewer-facing read names are produced later by the report layer and do not influence scientific ordering or placement.

## Tests

Unit tests cover deterministic SHA ordering with coverage/overlap evidence bound to that registry, forward/reverse profile projection, missing-profile preservation, misindexed-locus rejection, sparse all-reference overlap, explicit reference support at differential loci, filtered variant retention with factorized support topology, duplicate locus/variant support rejection, incompatible identities, and renamed duplicate content.

## Status

Implemented.
