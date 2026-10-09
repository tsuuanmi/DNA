# Sample Callability Requirements

**Requirement namespace:** `SRS-SAMPLE-024..SRS-SAMPLE-027`

These requirements are part of the canonical [sample evidence SRS](README.md) and cover how read callability ([SRS-CALL](../callability.md)) enters sample evidence.

- **SRS-SAMPLE-024:** A read's call that is masked (SRS-CALL-005) MUST observe a retained differential locus with state `masked`, keeping its reference-oriented primary base, relative quality, profile, and noisy-region context for review. A masked observation MUST NOT retain a locus by itself, MUST NOT enter pairwise comparable denominators, and MUST NOT contribute nucleotide mass (SRS-SAMPLE-018).
- **SRS-SAMPLE-025:** Locus support topology MUST count masked observations as `masked_reads`, and the state counts MUST still sum to total reads.
- **SRS-SAMPLE-026:** Every read MUST publish `callable_reference_segments` (SRS-ALN-015) beside its mapped reference segments, so that a covered but masked position can be told apart from a canonical reference match (SRS-SAMPLE-005).
- **SRS-SAMPLE-027:** A read with fewer than `callability.minimum_callable_calls` unmasked calls MUST NOT fail the `sample` operation: it MUST be recorded in `rejected_reads[]`, ordered by SHA-256, with its name, SHA-256, integrity, callability view, and reason `callable_calls_below_minimum` with the callable and minimum counts, and MUST contribute no alignment, coverage, overlap, locus, variant, or notation evidence. Read names MUST be unique across `reads[]` and `rejected_reads[]`. A sample whose every read is rejected MUST fail with a typed sample error and publish nothing.
