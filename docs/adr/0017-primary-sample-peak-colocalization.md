# ADR-0017: Gate Secondary Calls at the Primary Peak Sample

- **Status:** Accepted
- **Date:** 2026-09-07

## Context

`dna.peak_recall/v2` selects the strongest positive local maximum independently in each A/C/G/T channel inside a PLOC midpoint window. A secondary channel qualifies when its selected peak reaches `secondary_peak_ratio` relative to the uniquely strongest peak. Because the channel maxima may occur at different samples, a strong neighboring event elsewhere in the same window can create a false ambiguity at the current locus.

PLOC positions, midpoint windows, per-channel peak selection, and deterministic primary ranking are otherwise established behavior. The compact result schemas intentionally omit full signal arrays, and one trace does not justify genotype or heteroplasmy claims.

## Options

1. Keep the window-wide selected-peak ratio unchanged.
2. Add an absolute or spacing-derived peak-distance threshold and new configuration fields.
3. Preserve selected peaks and primary ranking, but require each qualifying channel to pass the existing ratio both at its selected peak and at the selected primary peak sample.

## Decision

Choose option 3 and identify the method as `dna.peak_recall/v3`. For a uniquely strongest positive primary peak with height `top` and sample position `p`, a channel qualifies only when both conditions hold:

```text
selected_peak_height > 0
selected_peak_height / top >= secondary_peak_ratio

channel_signal[p] > 0
channel_signal[p] / top >= secondary_peak_ratio
```

The primary channel necessarily passes both checks because its selected height and signal at `p` are the same value. Exact strongest ties and non-positive strongest peaks remain unresolved. PLOC windows, local-maximum selection, PLOC fallback, primary selection, and the existing one/two/three/four-channel call semantics remain unchanged.

No configuration key or output field is added. The existing configuration checksum continues to identify the effective threshold, and the closed `dna.analysis/v5` and `dna.basecalls/v1` schemas remain unchanged.

## Consequences

Qualifying channels under v3 are a subset of those under v2 for every uniquely strongest peak. The change can remove a remote secondary channel but cannot introduce one or change the uniquely strongest primary channel. It may therefore convert prior ambiguous calls to canonical calls and affect downstream ambiguity penalties, trimming, alignments, and variants; the method version and changelog make that scientific behavior change explicit.

The one-sample gate is deliberately narrow. It does not establish complete peak overlap, calibrated mixed-signal evidence, genotype, or heteroplasmy. A future spacing-aware locus-evidence model may use peak width, prominence, local spacing, and a validated distance rule under a new method version.
