# ADR-0062: Read callability from the read's own calls

- **Status:** Accepted
- **Date:** 2026-10-08
- **Related decisions:** [ADR-0027](0027-mixed-supporting-signal-snv-eligibility.md),
  [ADR-0060](0060-separate-variant-canonicalization-nomenclature.md)

## Context

Against the reviewed Sequencher calls for 44 real mtDNA samples, almost every
false positive after per-read nomenclature was evidence, not naming:

- calls within about ten bases of a read's aligned end, where an insertion is
  cheaper to align as flank mismatches (HVS-III reverse reads giving `301C 302C`
  instead of `309.1C 309.2C`);
- calls right after a long homopolymer in the sequencing direction, where
  polymerase slippage shifts phase (forward reads giving `310C` after the
  303-309 C run, `16193.1C`/`16194C` after the C tract created by `16189C`).

Tracy's quality trimming, which DNA's `quality_control` already mirrors, keeps
these regions. The legacy `mtdna_raw` pipeline removes them with per-primer
callable ranges derived from file names, which DNA's sample contract forbids.

## Decision

Variant eligibility gains two read-callability rules that use only the read's
own calls, in trace order (the sequencing direction for both strands):

- `read_end`: a mapped call within `read_end_margin` calls of either end of the
  retained interval;
- `post_homopolymer`: a mapped call within the `post_homopolymer_window` calls
  that start at the last call of a run of at least `homopolymer_min_length`
  identical canonical primary calls.

A variant with any untrusted mapped call is ineligible and stays in observed
evidence with the reason. Thresholds are configuration (schema 6, defaults
`read_end_margin = 8`, `homopolymer_min_length = 8`,
`post_homopolymer_window = 7`); `0` disables a rule.

## Alternatives

- Hard-trimming reads before alignment: loses evidence, and a mid-read window
  cannot be expressed as one retained interval.
- Primer callable ranges from file names: contradicts the contract that
  placement never depends on filenames.
- Removing everything after a long homopolymer: measured recall fell from 0.979
  to 0.83.

## Consequences

Measured with the production binary on the 44-sample Sequencher set
(notation calls inside each sample's analyzed range):

| `read_end_margin` / window | Precision | Recall | Exact samples |
| --- | --- | --- | --- |
| 0 / 0 (rules off) | 0.885 | 0.979 | 21 |
| 0 / 7 | 0.914 | 0.973 | 22 |
| 5 / 7 | 0.934 | 0.970 | 25 |
| **8 / 7 (default)** | **0.949** | **0.970** | **28** |
| 10 / 7 | 0.953 | 0.963 | 29 |
| 12 / 7 | 0.952 | 0.947 | 25 |

Margins above 8 start removing true calls that quality trimming already leaves
near a read end, such as `16189C` where both HVS-I reads are trimmed at the C
tract and `73G` near the start of HVS-II forward reads. The default therefore
favours recall over the last few false positives. Windows of 4-10 calls give
nearly identical results. The thresholds are empirical: revisit them as
validation data grows. The rules are
modality-specific Sanger evidence, applied before canonicalization as
ADR-0060 §5 requires.
