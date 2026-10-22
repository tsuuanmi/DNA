---
id: ADR-0080
status: accepted
date: 2026-10-19
family: [callability, sanger-eligibility]
related: [ADR-0027, ADR-0069, ADR-0071, ADR-0074, ADR-0078]
---

# ADR-0080: A modality may assert a call confirms only the reference

<!-- docs:begin lifecycle -->
- **Status:** Accepted
- **Date:** 2026-10-19
- **Related:** [ADR-0027](0027-mixed-supporting-signal-snv-eligibility.md), [ADR-0069](0069-plugin-first-modality-core-post-calling.md), [ADR-0071](0071-run-length-edits-need-resolved-run-ends.md), [ADR-0074](0074-run-structure-consensus.md), [ADR-0078](0078-call-the-consensus-decisions.md)
<!-- docs:end lifecycle -->

<!-- docs:begin backlinks -->
- **Governs:** [design/consensus](../../design/consensus.md), [design/second-peak-cause](../../design/second-peak-cause.md), [engineering/improvement-method](../../engineering/improvement-method.md), [reference/consensus](../../reference/consensus.md), [reference/read-evidence](../../reference/read-evidence.md), [reference/review](../../reference/review.md), [requirements/api](../../requirements/api.md), [requirements/callability](../../requirements/callability.md), [requirements/consensus](../../requirements/consensus.md)
<!-- docs:end backlinks -->

## Context

The Sanger plugin withholds the base of a call whose second channel
qualifies (`mixed_supporting_dna`, ADR-0027). Such a position, read by no
other read cleanly, stays undecided, and the review lists it as a check. The
second-peak cause (SRS-CALL-021) shows that many of these second peaks are
the signal's own shape: the tail of a neighbouring call of that base, or a
wide dye blob. They are not a second sequence.

Lifting the veto for those causes was measured on development panel v5
against `main` 65c221e, in steps. Lifting it outright cleared about 930
samples, but samples with an unreviewed difference rose in every period. The
harm came through four paths:

- one-strand calls deciding the reference where a variant hid;
- decided reference neighbours unlocking a weak one-read variant (`530G`)
  that a run end left uncallable;
- an undecided neighbour no longer making a variant's item a check;
- a base-only lift removing a call's run count that had opposed a length.

## Decision

A modality may assert, per call, that the call confirms only the reference.
The core reads that assertion through the kernel contract and never the
modality's reasons.

**Contract** (`dna-kernel`, `CallEvidence`):

- `reference_only`: the modality withholds the call's base, yet asserts that
  its signal shows that base where it is the reference's. The call may
  confirm the reference at its own position, and nothing else.
- `second_unexplained`: the call shows a second peak standing out that the
  modality does not explain by the signal's own shape, so possibly a second
  sequence.
- `SecondBase.stands_out`: the second base's share exceeds the modality's own
  floor, so no threshold lives in the core.

**Sanger** sets `reference_only` on a mixed call withheld for nothing else,
in phase between in-phase single-channel calls, whose second-peak cause is
`neighbour_tail` or `dye_blob` (SRS-CALL-022). The mixed reason stays on the
call verbatim, and its run role is unchanged.

**Consensus** (SRS-CONS-010):

1. A reference-only call whose base is the reference's makes its position's
   base observation clean, as a confirmation of the reference. Its run role
   (counting, run ends) stays its call's.
2. A confirmation is withdrawn where any read shows an unexplained second
   peak at the position, and where a read of the other strand shows the
   confirming call's own second base standing out (narrowed twice below:
   one-strand unexplained seconds, and explained seconds on both strands).
3. A confirmation within or beside a decided site that differs from the
   reference, including a site written partly as `N`, is withdrawn. The
   positions are then decided again, until no confirmation is withdrawn.
4. The positions only confirmations observe cleanly are published as
   `reference_only` in `dna.consensus/v1`. When the decisions are called
   (ADR-0078), the call at such a position shows its base only: it bounds no
   neighbouring run.

### Read-start calls

A second class of calls confirms the reference too. The class was measured
on development panel v5 against `main` b49ff87. Reference-reading calls that
only the peak floor or relative quality withholds, or that sit at a benign
signal edge, were wrong in 0 of about 30,000 cases (per-call research dump,
14 runs), while their differences were wrong 15–99 % of the time.

The Sanger plugin marks such calls reference-only (SRS-CALL-022). It owns
which edges are benign:
- the untrimmed start of a read (not its untrimmed 3' end, where the
  signal decays);
- a weak-signal mask clear of every repeat run (beside one, a weak segment
  is the dephased tail).

The core waives only the modality's own edge verdict for these calls. The
read's unvouched ends still stand.

The core withdraws a confirmation in three more cases:
- at or beside any read's observed difference, clean or not, since a read
  shifted in a repeat shows the reference beside its own difference;
- in a stretch, by its run structure or position by position, since a
  confirmation shows a base, not a length;
- in a reference run of at least `consensus.seeded_run` bases across a
  junction no read observes cleanly, where a length may change unseen.
  Applied to every run, this cost half the gain without removing a real
  error.

Measured result: samples with an unreviewed difference unchanged in every
period; clear samples 68.7 → 70.2 % on current data, 47.4 → 47.7 % in
2025-08..11, and 67.2 → 68.5 % older; no strict or `haplotype_lod` changes.
Plugins `sanger` 12 and `consensus` 12.

### Start-tail confirmations

The two Sanger rules were then widened together, measured on development
panel v5 against `main` 20fc7ca:
- a single-channel call may confirm the reference with a second base
  standing out when its cause is `neighbour_tail` or `dye_blob`;
- an isolated mixed tail or blob call keeps its confirmation at a signal edge
  whose every cause is benign, by the same test as a single-channel call;
- an isolated mixed tail or blob call may also be withheld for the peak floor
  or its relative quality, and a neighbour masked as weak before the trim
  interval counts as single-channel.

Every core guard is unchanged: an unexplained second peak on any read, the
other strand's same second base, beside-difference withdrawal, the run
junction, and stretches. An `off_centre`, `second_template`, or `phase_shift`
second never lets a call confirm, and a read's physical end is never benign.

Measured result: samples with an unreviewed difference and with a non-het
one unchanged in every period, no status of any measure changed; clear
samples 71.00 → 72.05 % on current data (+113), 48.59 → 49.29 % in 2025-08..11
(+19), and 69.41 → 70.84 % older (+43). Plugin `sanger` 14.
### One-strand unexplained seconds

Amended 2026-10-19. Rule 2 withdrew a confirmation wherever any read showed
an unexplained second peak. It now applies to a read of the confirming
call's strand as before. A read of the other strand withdraws it only when
its second base does not stand out, or when a read of the confirming strand
shows that second base standing out too. The other-strand guard on the
confirming call's own second base is unchanged. This narrows rule 2 and
the last of the consequences below.

An unexplained second on one strand alone is not the two-strand signature
of heteroplasmy the guard exists for. The strand that confirms reads its
base cleanly, apart from an explained artefact. A second on the confirming
strand still hides the confirmation: at 2563003 152 (a `Y` heteroplasmy) the
HV1F tail confirms C while HV2F, of the same strand, shows an unexplained C,
and the position stays undecided.

Measured on development panel v5 against `main` 20fc7ca:

| Period | Unreviewed | Non-het | Clear |
|---|---|---|---|
| current | 440 → 440 | 94 → 94 | 71.00 → 71.39 % (+42) |
| 2025-08..11 | 151 → 151 | 51 → 51 | 48.59 → 52.15 % (+96) |
| older | 174 → 174 | 84 → 84 | 69.41 → 70.01 % (+18) |

No status changes; 761 positions decided more, none fewer. One sample whose
only check was the undecided position becomes clear with its existing
unreviewed `309` length difference. Plugin `consensus` 14.

### Vouched calls within the end margin

Amended 2026-10-20. The core trusts no call within
`variant_calling.read_end_margin` calls of a read end the modality does not
vouch for (SRS-VAR-013), and Sanger never vouches. At the start of HVS-I the
first read calls fall inside that margin. The Sanger read-start rule
(SRS-CALL-018) already clears many of them of the modality's own edge from
the read's own signal, but the core margin kept them from deciding, so
16024-16035 stayed undecided wherever no read of the other strand reached
it.

A call within that margin now confirms the reference at its position, as a
reference-only call does, when all of these hold:
- the modality vouches for its base: scored, clear of its edges, and
  withheld for nothing;
- it shows no second base standing out and no second peak left
  unexplained;
- it reads the reference's base, and its alignment column is not shifted.

Every guard of a confirmation applies to it unchanged. The margin exists so
that no indel or run end rests on a read end's calls, and a confirmation
shows neither: its run role stays its call's, and stretches do not count it.
A reference-only call keeps the margin, since its base is withheld.

Measured on development panel v5 against the stack of #169, #171, #170, #173,
and #172 (integration 3498c4e):

| Period | Unreviewed | Non-het | Clear |
|---|---|---|---|
| current | 440 → 440 | 94 → 94 | 72.80 → 73.16 % (+38) |
| 2025-08..11 | 151 → 151 | 51 → 51 | 53.38 → 53.38 % |
| older | 174 → 174 | 84 → 84 | 71.63 → 71.96 % (+10) |

No status of any measure changes. All 48 samples were held only by
read-start checks at 16024-16045. Plugin `consensus` 15.

### Blob-wide neighbours

Amended 2026-10-20. A mixed tail or blob call confirmed the reference only
between in-phase neighbours with one qualifying channel. The guard kept a
call beside another mixed call from confirming, for fear of a second
sequence. A dye blob, however, is wider than one call, so the calls beside a
blob call are mixed too, and their own second peaks are the same blob or a
neighbour's tail. The guard then left whole blob stretches undecided on one
strand.

A neighbour now also counts as showing one base when its own second peak is
a neighbour's tail or a dye blob (SRS-CALL-021). It must still be in phase,
or masked as weak before the trim interval. A neighbour whose second peak is
`second_template`, `phase_shift`, `off_centre`, `noise`, `context`, or a
mixed stretch still stops the call. The modality cannot ask whether the
neighbours read the reference; every core guard of a confirmation holds
that line unchanged. A confirmation is withdrawn beside any read's observed
difference and beside a decided site that differs from the reference. Each
neighbour confirms only its own position, and only where it reads the
reference.

Measured on development panel v5 with `dna-eval diff` against the previous
commit (the edge-deletion zone and the missing-peak rule, PR stack on
`feat/eval-diff` ae98063). GATE PASS:

| Period | Clear | Unreviewed | Non-het passed zones |
|---|---|---|---|
| current (2026-01 and later) | 76.19 → 76.51 % (+35) | 387 → 387 | unchanged |
| earlier | 64.78 → 65.63 % (+49) | 280 → 280 | unchanged |

- **Workload.** Current `interior@HVS-I` checks fall 303 → 268, and
  `interior@HVS-II` 219 → 206. Over all periods, `end@HVS-I` falls
  403 → 338 and `after@HVS-III` 750 → 684. Positions to review fall
  141,588 → 136,579.
- **Unchanged.** No status of any measure changes, and no coverage changes.
- **Het.** Three newly clear samples carry an unreviewed heteroplasmy in a
  zone no check touched before either: 2538274 `146Y`, 2566906 `227R`, and
  LN_25_AA3929 `16182M`. In LN_25_AA2476 the zone of an unreviewed `16266Y`
  loses its check. The gate counts heteroplasmy as optional.
- **Plugin.** `sanger` 19.

### Explained seconds on both strands

Amended 2026-10-22. Rule 2 also withdrew a confirmation wherever a read of
the other strand showed the confirming call's own second base standing out,
whatever its cause. A reference-only call's own second is explained by
construction (a neighbour's tail or a dye blob, SRS-CALL-022). When the other
strand's read explains its second too, neither strand shows an unexplained
second there, and the review's two-strand heteroplasmy rule calls such a
second `explained` (rules version 19). The guard now withdraws only when the
other strand leaves the second unexplained, which the unexplained-second
part of rule 2 already covers; the separate guard is removed.

Census on development panel v5 at `main` af0823d3 (current blocked samples):
20 undecided positions in 19 samples were held only by this guard, all of
them reference in the truth, for example 111 between two C neighbours
(each strand sees the downstream C's tail) and 530. The heteroplasmies the
guard was kept for (`152Y` in 2563003, LN_26_AA3106, LN_26_AE0209) show an
unexplained `phase_shift` second on the other strand and stay undecided.

Measured with `dna-eval diff` against `main` af0823d3, GATE PASS on both
panels:

| Panel | Period | Clear | Unreviewed | Non-het passed zones |
|---|---|---|---|---|
| dev v5 | current | 8,813 → 8,827 (+14) | unchanged | unchanged |
| dev v5 | earlier | 4,049 → 4,065 (+16) | unchanged | unchanged |
| readiness-v2 | current | 15,158 → 15,165 (+7) | unchanged | unchanged |
| readiness-v2 | earlier | 2,011 → 2,024 (+13) | unchanged | unchanged |

- **Unchanged.** No unreviewed block, status, or coverage changes on either
  panel.
- **Het.** Two newly clear dev v5 samples carry an unreviewed heteroplasmy
  in an interior zone no check touched before either (LN_26_AB5943 `16301Y`,
  LN_25_AA9640 `16318R`).
- **Omission.** No newly clear sample omits a truth allele, except one
  earlier readiness-v2 sample (2541327) whose `73G` lay in a `not_analysed`
  stretch before and after; its only check, at 83, is now decided.
- **Plugin.** `consensus` 22.

## Alternatives

- **Lift the veto by cause in the modality.** This gave the largest gain and
  failed the review gate through the four paths above.
- **Require a clean read of the other strand.** Nothing was left to gain:
  where the other strand reads the base cleanly, the position was already
  decided.
- **Unresolve every differing site beside a confirmed position.** This held
  the gate, but turned about 100 correct calls into `N`.
- **Let a modality-vouched call within the end margin also bound runs or
  support variants.** The margin guards exactly those uses, and the
  read-start residue needs only the reference confirmed.
- **Let the core read the modality's reason labels.** This was rejected:
  reasons are reported verbatim, and the core interprets none of them
  (ADR-0069).

## Consequences

Development panel v5 against `main` 65c221e, before the version bumps:

| Period | Unreviewed | Clear |
|---|---|---|
| current | unchanged | 66.1 → 68.7 % |
| 2025-08..11 | unchanged | 42.3 → 47.4 % |
| older | +1 (`538`, truth `N`) | 60.0 → 67.2 % |

- No strict status changes.
- About 630 samples become clear.
- Heteroplasmy stays the review's two-strand rule. A confirmation never
  stands where another read shows an unexplained second peak.

`dna.consensus/v1` gains `reference_only`, revised in place before its first
release. Plugins `sanger` 11, `consensus` 11, and `decided` 3.
