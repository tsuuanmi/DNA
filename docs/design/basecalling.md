# Basecalling Method

Part of the canonical [DNA pipeline](pipeline.md).

Re-calls every vendor-defined locus from the channel signals. Vendor base
strings are retained as evidence but never replace signal-derived re-calling.

### Substep 2.1 — Call windows

For each basecall position `p[i]`, a half-open sample window is built from the
midpoints of neighboring positions:

- first window: `[p[0] - (p[1]-p[0])/2, midpoint(p[0], p[1]))`;
- interior window `i`: `[midpoint(p[i-1], p[i]), midpoint(p[i], p[i+1]))`;
- last window: `[midpoint(p[n-2], p[n-1]), p[n-1] + (p[n-1]-p[n-2]+1)/2)`,
  clamped to the sample count.

`midpoint(a, b) = a + (b - a) / 2`. At least two basecall positions are
required.

### Substep 2.2 — Per-channel peak selection

Within each window, each channel is searched for a positive local maximum. A
sample `v` at position `j` is a local maximum when
`(v[j-1] <= v && v > v[j+1]) || (v[j-1] < v && v >= v[j+1])`. The highest such
sample is the channel peak. If no positive local maximum exists, the channel
value at the basecall position is used as a fallback. Each channel peak records
its base, height, position, and source (`local_maximum` or `ploc_fallback`).

### Substep 2.3 — Call decision

The four channel peaks are ranked by height (ties broken by channel order
A < C < G < T). Let `top` be the highest height. If `top <= 0` or the second
peak ties `top`, the call is unresolved (`N`). Otherwise, let
`primary_peak_position` be the selected peak position of the uniquely strongest
channel:

- **Qualifying channels** have a positive selected peak satisfying
  `selected_height / top >= secondary_peak_ratio` and positive channel signal at
  `primary_peak_position` satisfying
  `signal_at_primary_peak / top >= secondary_peak_ratio`.
- The two conditions are an intersection with the v2 selected-peak rule, so a
  remote maximum elsewhere in the same call window cannot create ambiguity and
  no new secondary channel can qualify.
- **Primary** is the strongest base when one to three channels qualify; four qualifying channels produce `N`.
- **Ambiguity** depends on the number of qualifying channels: one base maps to
  itself (canonical); two bases map to the standard two-base IUPAC symbol; three
  bases are unresolved `N` (primary is still the strongest); four bases are
  unresolved `N` for both primary and ambiguity.
The primary and ambiguity sequences are the concatenation of the per-call
primary and ambiguity symbols.
