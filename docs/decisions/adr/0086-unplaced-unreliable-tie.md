---
id: ADR-0086
status: proposed
date: 2026-10-20
family: [placement]
related: [ADR-0029, ADR-0069, ADR-0085, PROP-0012]
---

# ADR-0086: A placement tie confined to a read's unreliable end leaves those calls unplaced

## Context

[[ADR-0029]] makes two genuinely distinct, equally scoring placements of
the selected orientation an error, and the sample commands reject the read
as `ambiguous_placement` (SRS-SAMPLE-029). On development panel v5 this
rejected 31 reads, 23 of them on current data.

In every one of them the two placements agree on the read's start and on
every call except the last 3 to 15 query calls at one end of the read. There
the two explanations, an insertion or deletion against substitutions, end
at different reference positions and score exactly the same (fixed-point
integers).

- **Read start (29 reads).** 25 HV1R reads at 16385-16404 and 4 HV3R reads
  at 588-605. Every disputed call lies inside the read's unreliable span
  (`poor_start_signal`, SRS-CALL-024). In 21 of them every disputed call
  lies before the quality-trim start, where SRS-CALL-014 extends the
  informative interval over calls vouched for by their base only.
- **Read end (2 reads).** An HV1F read and an HV2F read whose disputed
  calls are all masked as anchoring: the dephased calls at the end of a
  C run.

Rejecting the read throws away 200 to 750 good calls because a few calls
at its poorest end admit two explanations. The region is then read by one
strand only, or not at all.

## Decision

When the selected orientation's two best placements tie, DNA takes the
calls they place at different reference positions (an inserted call at
none). It looks at the stretch from one end of the query to the farthest
of those calls. When every call of that stretch is unreliable, the stretch
is left unplaced and the rest of the query is aligned again under the same
contract (SRS-ALN-006). A call is unreliable when it lies inside one of the
read's unreliable spans (`ReadEvidence::unreliable`) or is masked.

- **Unplaced, not masked.** The calls stay in the read's evidence and keep
  their verdicts. They have no alignment column, so they observe no position
  and support or oppose nothing. Masking them as `N` would still place them,
  at a position chosen by gap cost alone.
- **Aligned again.** The remaining calls get the authoritative optimum of
  their own query: both orientations, canonicalization, the admission gates,
  and the uniqueness check. A tie there rejects the read. The rule does not
  repeat.
- **Genuine ambiguity stays rejected.** A tie that disputes a reliable call
  still rejects the read, and so does a tie between whole placements.
- **No new number.** The rule reads only the kernel's assertions: the
  read's spans and its masks, never a reason label, a position, or a
  primer. It holds no threshold.

### Owner

The rule belongs to `dna-align`, which owns placement and its ambiguity.

- **Dependency direction.** The tie exists only once the read is aligned.
  `dna-sanger` builds the evidence before any alignment and cannot see it,
  so narrowing SRS-CALL-014 there would have to give up the informative
  extension on every read, including the many that align uniquely with it.
  It would also leave the two read-end cases rejected.
- **Modality-neutral.** `dna-align` depends only on `dna-kernel`. It reads
  `ReadEvidence::unreliable` and `Placement`, which any modality may assert.
- **One more use of the span.** SRS-CALL-024 lets the core decide nothing
  from the unreliable span except where a region's analysed range ends
  ([[ADR-0085]]). This decision adds the second use: which calls of a tied
  placement stay unplaced.

### Variants considered

- **Narrow SRS-CALL-014 in the modality.** Rejected for the reasons under
  Owner.
- **Calls not trusted to place an indel.** No disputed call in the 31 reads
  is trusted for an indel: each is masked, at an edge, withheld, or vouched
  for by its base only. A rule on that verdict admits the same reads. It was
  not chosen. That verdict speaks to an indel a call would place, not to
  where the read lies. It also holds next to any masked call, deep inside a
  read.

## Consequences

Measured on development panel v5 against a fresh run of `main` 85cab63,
with the re-check list `--recheck` (`dna-eval diff`, GATE PASS). A period is
a range of batch dates.

| Period | Samples | Clear | Positions to review | Zones needing a check |
|---|---:|---|---|---|
| current (2026-01 and later) | 10,695 | 7,646 → 7,654 | 82,989 → 81,539 | 29 fewer |
| 2025-12 | 781 | 476 → 476 | −10 | start@HV2 +1, window@HV2 −1 |
| 2025-08..11 | 2,694 | 1,460 → 1,462 | −60 | 4 fewer |
| older | 2,240 | 1,608 → 1,609 | 34,314 → 33,970 | 6 fewer |

- **Admitted reads.** All 31 reads are now placed, one per sample. No other
  read changes.
- **Non-het escapes.** Non-het passed zones are unchanged in every zone and
  every period. No sample gains an unreviewed block.
- **Against the truth.** 20 truth alleles in 9 samples are now covered,
  and DNA calls 19 of them as the truth holds. Four samples become
  strict-concordant. The twentieth makes LN_25_AA1052 primary-discordant:
  the truth holds `16192Y`, and with the admitted HV1R read DNA calls
  `16192T`, inside a check (`16184-16193`). No admitted read adds any other
  allele the truth lacks.
- **2025-12.** In LN_25_AA9646 the admitted HV3R read reads the reference A
  at 111 with a second C, as the HV2F read does. The site becomes an
  `undecided` check, and the sample's 311-316 check goes.
- **Plugins.** `core` 12.

## Links

- Rule: SRS-ALN-006. SRS-CALL-024 allows the core this use of the span.
- Spans: [[PROP-0012]] increment 5. Placement: [[ADR-0029]].
