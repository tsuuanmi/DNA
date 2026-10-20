---
implements: [SRS-CALL, SRS-SIG]
decided-by: [ADR-0088]
---

# Electrophoretic Spikes

Part of the canonical [signal-processing method](README.md); decided by
[ADR-0088](../../decisions/adr/0088-electrophoretic-spikes.md), required by
SRS-SIG-012 and SRS-CALL-027.

## Signature

A spike is a narrow event every dye channel shows at once (a bubble, a
particle, an electrical transient at the detector). The raw filter channels
show it at one scan; the analysed channels DNA reads show it staggered by
each dye's mobility offset. On development panel v5 the canonical form is:

| Feature | Typical value |
|---|---|
| apex order | T first, A and C a few scans later, G last |
| apex stagger, first to last | 8–12 scans (one call spacing is about 12) |
| width at half height, per channel | 5–8 scans |
| height against the local primary peak | 3 to over 300 times |

Sample 2538054, HV2F: apexes T 1589, A 1594, C 1595, G 1601, at 8.3, 17.2,
12.6, and 14.0 times the local level of 547; in the raw channels all four
peak together at raw scan 2824, about 80 times the raw median.

## Method

1. **Level.** A call's level is the median primary peak, the highest of its
   selected peaks, of the `spike_level_calls` calls on each side, the call
   left out, and at least one. A channel sample takes the level of the call
   whose locus window holds it.
2. **Apexes.** A channel's apex is a local maximum (higher than the sample
   before, at least the sample after) at least `spike_height` times its
   level.
3. **Grouping.** In sample order, an apex not yet used opens a span of
   `spike_span` median call spacings (rounded up to whole samples). When
   every channel has an apex in the span, the highest apex of each channel
   forms a spike, and the next spike starts after the latest of the four.
4. **Extent.** The spike's samples run from the earliest of its four apexes'
   half-height starts to the latest of their half-height ends.
5. **Calls.** The spike's calls are those whose primary peak lies among its
   samples. A spike's leading edge often holds the primary peak of the call
   before its centre; the call after its trailing edge keeps its own peak
   outside and is not marked.

Grouping requires all four channels: a template shows one or two channels at
a call, and a dye blob one. The level is local because a read's signal falls
several-fold along its length.

## Use

The Sanger evidence withholds the calls of every spike (SRS-CALL-027,
reason `signal_spike`); see [read callability](../callability.md#consumers).
Signal processing itself only records the spikes and logs their count
(`spikes`, `spike_calls`) with `signal_processing_completed`.
