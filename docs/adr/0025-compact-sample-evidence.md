# ADR-0025: Factor sample evidence into a read registry and sparse differences

- **Status:** Accepted
- **Date:** 2026-09-19

## Context

ADR-0024 introduced the first multi-read evidence contract. Its dense locus model
proved scientifically explicit but operationally noisy: every covered reference
position repeated read name, SHA-256, orientation, state, call index, and relative
quality even when all reads simply matched the reference.

Exploratory review of a local four-read mtDNA sample output made the imbalance concrete. The document contained 1,027 locus records and 1,814 locus observations; 1,793 observations were reference matches, while only 12 loci contained any non-reference state. The resulting sample JSON was roughly 703 KB even though the scientifically interesting difference set was small. No local sample identifier, filename, hash, or derived sample document is committed as validation evidence; synthetic tests enforce the resulting representation invariants.

Read identity, orientation, and mapped reference segments are already known once
the independent read observations have been placed. Repeating those facts at every
locus does not add evidence. At the same time, simply deleting all reference
observations would lose useful reference-support quality at positions where another
read differs.

The output contract therefore needs to factor shared read metadata once, make
coverage semantics explicit, and retain detailed locus evidence only where a
difference exists.

## Decision

DNA replaces `dna.sample_evidence/v1` with `dna.sample_evidence/v2`. At
adoption, DNA emitted v2 only; no v1 compatibility result or alias was retained.

### Read registry

`reads[]` is the single deterministic registry of contributing reads and is sorted
by input SHA-256. Each entry owns reviewer-facing basename provenance, stable
SHA-256 identity, derived orientation, alignment summary, and mapped reference
segments.

Internal aggregation may use deterministic indexes into the SHA-sorted registry, but the public contract refers to reads by unique reviewer-facing filename stem. Numeric indexes are not exposed as reviewer identifiers. SHA-256 remains the scientific content identity and filename semantics never become placement or merge keys.

### Sparse differential loci

The dense `loci[]` table is replaced by `locus_differences[]`.

A locus is retained only when at least one covering read is alternate, unresolved,
or deletion. Once a locus is retained, every read covering that position is
included, including canonical reference observations with observed base and quality.

Aggregation uses two passes:

1. discover reference positions with at least one non-reference state;
2. collect all covering-read observations only for those positions.

Mapped reference segments define coverage. For a given read:

- a 1-based position inside its mapped segments but absent from `locus_differences[]`
  is a canonical reference match;
- a position outside its mapped segments is uncovered;
- a retained differential locus uses its explicit observations.

This preserves the distinction between reference support and missing coverage
without materializing routine all-reference positions.

### Normalized variants

`variants[]` remains the normalized variant layer. Per-read support keeps configured eligibility and exclusion reasons, refers to the read by unique reviewer-facing name, and exposes reference-oriented called base plus co-located A/C/G/T peak heights and quality. Original call index, PLOC, and mapped call coordinate remain internal.

Read-level filtering changes eligibility, not whether the normalized observation
exists in sample evidence.

### Versioning

The shape change is incompatible with v1, so it is a new schema version rather
than a mutation of the accepted v1 schema. The v1 schema/example and implementation are removed instead of maintained as a compatibility layer. The v2 repository example is synthetic and non-identifying.

## Consequences

- Sample output size scales primarily with reads plus differences rather than with
  every covered reference base.
- Reviewer-facing filenames remain easy to find without being repeated or used as
  scientific keys.
- Reference-support call quality remains explicit at loci where another read
  differs.
- Routine per-base reference quality is intentionally absent from the public sample
  document. Full one-read evidence still exists upstream before sample aggregation;
  future interpretation needing additional focused evidence should derive it there
  rather than restore a dense whole-coverage table.
- Variant support becomes more reviewable because scientific identity/orientation are factored to the read registry while each support uses the human-readable read name and direct peak/quality evidence.
- The internal sample implementation is split into validation/read ordering,
  differential-locus aggregation, and normalized-variant aggregation modules.

## Non-goals

This decision does not add consensus, majority voting, conflict adjudication,
genotype, heteroplasmy estimation, haplogroup inference, metadata-driven
placement, or canonical F/R pairing.

## Follow-up

ADR-0027 later advanced the sample contract to `dna.sample_evidence/v3` solely to add the closed-enum `mixed_supporting_dna` eligibility reason. ADR-0030 subsequently advances the current contract to `dna.sample_evidence/v4` by adding pairwise overlap/admission evidence. The compact read registry and sparse-difference architecture defined here remain authoritative.
