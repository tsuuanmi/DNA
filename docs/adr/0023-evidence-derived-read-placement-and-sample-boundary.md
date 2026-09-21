# ADR-0023: Evidence-derived read placement and sample boundary

- **Status:** Accepted
- **Date:** 2026-09-18

## Context

Tracy research showed that the useful abstraction is not "merge forward and reverse" or "classify this trace as HV1/HV2/HV3." For short references, Tracy lets a trimmed read locate itself by semi-global alignment: the query is consumed, unused reference flanks are free, both orientations are evaluated, and traceback determines the covered interval.

DNA already implements a stricter version of that primitive. `align_best()` evaluates forward and reverse orientations, rejects unresolved orientation ties and equal placements, supports circular projection, and returns mapped reference segments.

The same research also showed why pair-first merging is lossy for tiled Sanger samples. Reads such as HV2F and HV3R can overlap even though they are not a canonical pair. A usable read may also lack its nominal partner.

## Decision

DNA adopts the following production architecture:

1. Every trace is processed independently through signal-derived calling, signal analysis, quality control, evidence-driven reference alignment, and read-level variant extraction.
2. The resulting one-read scientific product is represented by `ReadObservation`.
3. Read orientation and covered reference segments are derived from alignment evidence, not filename, HV label, primer name, declared F/R direction, or an expected-region constraint.
4. Future sample reconciliation will consume `ReadObservation[]` in normalized reference-coordinate/variant space.
5. Canonical F/R pairing, amplicon, primer, and replicate identity are optional provenance/support dimensions, not exclusive merge keys.
6. A consensus sequence, when introduced, is a downstream projection. Sample variants and discordance must not be defined by diffing a flattened consensus string.

## Consequences

### Positive

- Unlabeled reads can be placed scientifically.
- Renaming a file cannot change placement.
- Cross-amplicon overlaps remain available for reconciliation.
- Missing canonical partners do not discard useful evidence.
- Circular origin-spanning reads retain explicit segment topology.
- The one-read and future sample-level models have a clean typed boundary.

### Costs

- Sample aggregation must model coverage, disagreement, support topology, and local contribution eligibility explicitly.
- Assay metadata needs a separate provenance/QC representation if introduced.
- Pairwise F/R consensus cannot be implemented as a special scientific shortcut.

## Non-goals

This ADR does not introduce multi-read CLI input, a sample JSON schema, evidence-profile alignment, FM indexing, or sample consensus behavior. Those require separate implementation and validation.

## Validation

The implementation must preserve existing single-read output while making the read observation boundary explicit. Future sample work must include filename/metadata invariance, reverse-orientation inference, circular-origin placement, ambiguous-placement failure, cross-amplicon overlap, missing-partner, and discordant-overlap fixtures.

## Follow-up

ADR-0024 implements the first sample consumer of this boundary: the multi-read
CLI, `SampleEvidence`, and `dna.sample_evidence/v2`. The non-goals above
describe the scope of ADR-0023 itself, not the repository state after ADR-0024.
