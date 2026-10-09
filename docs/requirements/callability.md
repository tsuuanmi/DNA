# Signal-derived Read Callability Requirements

**Requirement namespace:** `SRS-CALL-*`

These requirements are part of the canonical [DNA SRS](SRS.md).

- **SRS-CALL-001:** DNA MUST derive, for every read and without reference, target-profile, primer, or file-name knowledge, one per-position callability record in trace order from the immutable analyzed channels, canonical loci, basecall-independent `LocusEvidence`, and primary calls. The derivation MUST NOT modify channels, loci, calls, selected peaks, or locus evidence.
- **SRS-CALL-002:** Per-position features MUST follow the documented `dna.read_callability/v1` formulas — primary/secondary channel ranking, dominance, secondary ratio, slippage shift offset searched nearest first up to three calls, spacing deviation from the local median spacing, and weakness relative to the read's median primary amplitude — with six-decimal rounding before storage or comparison. Every threshold MUST be configuration (SRS-CFG-008).
- **SRS-CALL-003:** Phase segmentation MUST assign every call position exactly one state from the closed set `in_phase`, `dephased`, `mixed`, `weak`, `irregular`, as ordered 0-based half-open call intervals that partition `[0, call_count)`. A masked stretch MAY return to `in_phase` after it; a callable island shorter than one window MUST be absorbed by its masked neighbours.
- **SRS-CALL-004:** Repeat runs of at least `variant_calling.homopolymer_min_length` identical or alternating canonical primary calls MAY lower the onset threshold in the window that follows the run and MUST attribute a segment that starts in that window (`after_repeat`), but MUST NOT mask by themselves or lower any per-position threshold inside the run.
- **SRS-CALL-005:** The mask MUST mark every position of a non-`in_phase` segment; the callable span MUST be the interval from the first to the last unmasked position and MUST be empty when every position is masked. Until the mask is allowed to act, callability MUST NOT alter trim bounds, alignment, warning totals, or variant eligibility, and public contracts MUST publish it as observation.
- **SRS-CALL-006:** Operational logs MUST record callability as aggregate metrics (repeat and segment counts per state, masked-call count, callable span, one bounded token per segment) and MUST NOT emit per-position arrays.
- **SRS-CALL-007:** Phase states, mask reasons, and callability features MUST NOT be described as Phred, error probability, heteroplasmy, mixture fraction, or artifact classification.
- **SRS-CALL-008:** The callability core MUST operate on plain per-position records; decoding Sanger types into those records MUST live in one adapter that owns no thresholds, and no modality trait is introduced until a second adapter exists.
