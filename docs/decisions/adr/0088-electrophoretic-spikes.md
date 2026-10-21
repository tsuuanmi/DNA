---
id: ADR-0088
status: accepted
date: 2026-10-21
family: [signal-evidence, sanger-eligibility]
related: [ADR-0027, ADR-0031, ADR-0069, ADR-0080]
---

# ADR-0088: A call inside an electrophoretic spike shows no base

## Context

An electrophoretic spike is a narrow event that every dye channel shows at
once: an air bubble, a particle, or an electrical transient passing the
detector. In the raw filter channels (`DATA.1`–`DATA.4`) its four peaks
coincide; the instrument's mobility correction then moves each dye by its own
offset, so in the analysed channels (`DATA.9`–`DATA.12`) that DNA reads the
four apexes lie within about one call spacing, T first and G last on the
development panel's instrument. Each channel is many times the read's local
signal and about half a call spacing wide at half height.

A call whose primary peak lies inside a spike reads the spike's channel, not
the template. In sample 2538054 (HV2F, scans 1586–1605, apexes 8–17 times the
local level) the spike's leading T edge read cleanly, with a second share of
0.04, at 186, where the truth has the reference C. The HV3R read covering 186
had been set aside as weak, so DNA passed `186T` unflagged in an interior
zone. No existing rule saw it: the call has one qualifying channel, an
unremarkable second share, regular spacing, and lies in phase.

[ADR-0031](0031-trace-integrity-evidence.md) keeps trace-integrity evidence
unthresholded until an artifact classifier is separately specified and
validated. This decision specifies one, for one artefact, from its own
shape; it does not threshold the event-signal ratio.

## Decision

Signal processing finds a read's spikes (SRS-SIG-012) from the analysed
channels and the read's own calls: within `signal_processing.spike_span`
(1.5) median call spacings, every channel has a local maximum at least
`signal_processing.spike_height` (3.0) times the local level, the median
primary peak of the `signal_processing.spike_level_calls` (6) calls on each
side. A spike's samples run from the earliest apex's half-height start to the
latest apex's half-height end; its calls are those whose primary peak lies
among them. All four channels are required: a real template shows one or
two channels at a call, and a dye blob shows one.

The Sanger evidence withholds every call of a spike for its base and for
indels with the reason `signal_spike` (SRS-CALL-027): it shows no run's end,
confirms no reference (ADR-0080). No rule
that vouches for a call or takes a verdict off an edge lifts it. The core
sees only the withheld verdicts (ADR-0069), so a spiked position is read by
the sample's other reads, or stays undecided and becomes a check.

## Evidence

A census of development panel v5 (66,211 reads, `main` 85cab63) screened
every read for four-channel events and ran the rule on the 103 reads that
showed one: it withholds 124 calls in 60 reads of 59 samples, about one read
in 1,100.
Of the 24 calls it newly withholds (the others were already withheld or
masked), 12 have a truth base: 9 read another base than the truth (one the
unflagged `186T` above, seven already in a strand-disagreement check) and 3
read the truth. The other 12 lie outside the truth's analysed range.

Against `main` 85cab63 on the whole panel (`dna-eval diff`, gate passed):
the unflagged `186T` becomes an undecided check; samples with an unreviewed
block 618 → 617; passed non-het units unchanged in every other zone and
period; five samples whose spiked read had disagreed with the other strand
become clear (current 7,646 → 7,651); positions to review 164,199 → 164,139;
no coverage change; one sample concordant under every measure.

## Consequences

- A spike can no longer pass a base its own channel invented. Where another
  read covers the position the sample is unchanged; where none does, the
  position is undecided and reviewed.
- A spike over a true base costs that read's call. On the panel this was
  3 calls with a known truth.
- `signal_spike` is a new Sanger reason label (`dna.read_evidence/v1`,
  `dna.sample_evidence/v10` exclusion reasons), and three keys join
  `[signal_processing]` (configuration schema 21).

## Alternatives

- **Threshold the event-signal ratio** (ADR-0031): one per read, it cannot
  say which calls a spike covers, and a strong but real read start reaches
  it too.
- **Mark calls whose own window shows all four channels high:** the dyes'
  mobility offsets spread a spike over two call windows, so no single call
  shows all four channels; the spike's leading-edge call, the harmful one,
  shows three.
- **Mark every call whose window overlaps the spike:** marks the true call
  just after the spike's trailing edge too.
