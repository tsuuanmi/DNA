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
- **Governs:** [design/consensus](../../design/consensus.md), [design/second-peak-cause](../../design/second-peak-cause.md), [reference/consensus](../../reference/consensus.md), [reference/read-evidence](../../reference/read-evidence.md), [reference/review](../../reference/review.md), [requirements/api](../../requirements/api.md), [requirements/callability](../../requirements/callability.md), [requirements/consensus](../../requirements/consensus.md)
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
   confirming call's own second base standing out.
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

## Alternatives

- **Lift the veto by cause in the modality.** This gave the largest gain and
  failed the review gate through the four paths above.
- **Require a clean read of the other strand.** Nothing was left to gain:
  where the other strand reads the base cleanly, the position was already
  decided.
- **Unresolve every differing site beside a confirmed position.** This held
  the gate, but turned about 100 correct calls into `N`.
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
