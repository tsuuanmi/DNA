# ADR-0046: Anchor nucleotide evidence to the nearest locus event

- **Status:** Accepted
- **Date:** 2026-09-19

## Context

`LocusEvidence` previously selected the sample with maximum total baseline-corrected
A/C/G/T amplitude anywhere inside the neighboring-midpoint PLOC window.

Real-trace validation showed that a stronger neighboring chromatogram event can occur
inside the same window and displace the selected evidence sample by roughly one local
base spacing. This produced loci where the primary/reference call was coherent at its
own peak but the exported `EvidenceProfile` described the neighboring event instead.

Representative validation observations showed event-to-primary-peak offsets of -5 to -6
samples together with extreme profile disagreement, including forward/reverse Total
Variation near 1.0. The diagnostic evidence supported event mis-association rather than
ordinary low-level background as the primary mechanism.

DNA still needs a nucleotide-agnostic event locator. Selecting the called nucleotide's
peak directly would couple `EvidenceProfile` composition to the basecall verdict.

## Decision

For each validated PLOC locus, DNA derives total non-negative baseline-corrected
A+C+G+T amplitude per chromatogram sample and considers only positive local maxima of
that total-signal series inside the locus window.

The event-selection order is:

1. smallest absolute sample distance to the locus PLOC;
2. greater total corrected signal when two local maxima are equally distant;
3. lower sample coordinate for an exact remaining tie.

If the locus window contains no positive total-signal local maximum, DNA uses the
validated PLOC sample directly.

The local-maximum plateau rule is the same deterministic asymmetric rule used for
channel-local peak detection:

```text
(left <= current > right) OR (left < current >= right)
```

`EvidenceProfile`, corrected amplitudes, SNR, and trace-integrity event-signal summaries
continue to derive from the selected co-located event sample.

## Consequences

- A stronger neighboring base event can no longer steal a locus profile merely because
  it lies inside the same midpoint window.
- Event placement remains independent of primary/ambiguity calls, selected channel
  identity, qualifying-channel membership, and alignment/reference sequence.
- Small chromatogram-sample offsets from PLOC remain allowed when the nearest actual
  total-signal event is slightly shifted.
- Existing validation event diagnostics remain sufficient to measure event/PLOC and
  event/primary-peak offsets after this change.
- Alignment and downstream profile geometry may change because they consume corrected
  locus evidence; this is an intentional scientific correction, not a compatibility path.

## Non-goals

No basecall rescue, called-base weighting, profile threshold, heteroplasmy inference,
artifact classifier, consensus rule, or public schema field is introduced.
