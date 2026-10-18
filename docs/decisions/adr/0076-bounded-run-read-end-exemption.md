# ADR-0076: Bounded runs before a post-repeat phase loss are exempt from `read_end`

- **Status:** Accepted
- **Date:** 2026-10-16
- **Related decisions:** [ADR-0062](0062-read-callability.md),
  [ADR-0067](0067-signal-derived-read-callability.md),
  [ADR-0071](0071-run-length-edits-need-resolved-run-ends.md),
  [ADR-0074](0074-run-structure-consensus.md),
  [ADR-0075](0075-evaluation-harness-in-workspace.md)

## Context

`read_end` (ADR-0062, ADR-0067) distrusts every call within
`variant_calling.read_end_margin` informative calls of an uninformative call.
The margin exists for alignment-edge artifacts: at a read's start or end, and
where callability masks a segment as unresolved, a difference close to the edge
may come from the alignment, not from the template.

The evaluation harness (ADR-0075) measured the HVS-II poly-C window on the
development panel. Where reviewers record `315.1C`, the reverse read often
reads 316 → 303 in phase, with both run ends (`T310`, `G316`) as in-phase
calls, and then loses phase right after the `A302` run exit. Its calls over
303–316 then lie within the margin of that post-repeat segment, so the
consensus has no clean read there and writes `N`. The evidence that decides
the run length is in the trace; the margin hides it.

Lowering the margin for every edit recovered many of these cases but brought
back the artifacts the margin stops (measured on the development panel, about
2,600 samples):

- substitutions at a run exit (`310C`, `16194C`);
- deletions at a read start;
- isolated substitutions near read ends.

The Researcher session's analysis of those regressions gave the rule below.

## Decision

A run of one base that the read bounds is exempt from `read_end` when the
margin's edge is the read losing phase right after a repeat run, not a read
end. All of these must hold:

1. **Run-length evidence only.** The exemption covers the calls of the run:
   - for an insertion or deletion that changes the run's length, the variant
     is not `read_end`;
   - in the consensus, the run's calls are clean where they match the
     reference or are inserted: the run's length may change, its composition
     may not.

   A substitution is never exempt: in `dna sample` eligibility, and in the
   consensus, where a run call aligned to a reference base of another base
   keeps the margin. In the consensus, a run's two bounds are clean where
   they match the reference (revision of 2026-10-16 below); a bound that
   reads another base keeps its own trust, so a substitution at a run exit
   still needs a read that is clean there without the exemption.
2. **Bounded run.** The run (two or more calls of one canonical base) has a
   resolving call on each side (ADR-0071). The bounds and every run call are
   unmasked.
3. **Post-repeat edge.** On each side, either the nearest uninformative call
   is beyond the margin, or it is a masked call whose segment the modality
   places right after a repeat run.
   A read's physical end, the trim interval, or any other masked segment
   never exempts.
4. **Guard.** At least two informative calls separate each bound from such an
   edge. This is a structural minimum, not a tuned value: guards of 0, 2, and
   3 gave identical results on the development panel.

Since `ReadEvidence` v2 (PROP-0008 step 4), the modality applies rules 2 to 4
on the read's own calls and marks each such run and its bounds in
`ReadEvidence` (`post_repeat`, SRS-CALL-016); the Sanger plugin finds the
post-repeat edge where callability marks a segment `after_repeat`. The core
keeps the alignment-dependent half: which of the marked calls it trusts (the
run's matching or inserted calls, a bound that reads the reference), and that
the edit lies on the run and its bounds. It does not interpret reason
labels. The exemption claims nothing
biological (ADR-0009); it only stops distrusting run calls the read shows in
phase.

## Alternatives

- **A smaller global margin (12 → 8).** Rejected. It fixed about 77 samples on
  the primary measure but broke 22, mostly through the artifacts above.
- **Also exempting segments that start inside a run.** Rejected. It added two
  samples but introduced a false `308DEL` and two haplotype regressions.
- **Anchoring short mixed blips after a repeat.** Not needed. Those blips are
  `after_repeat` segments, so this rule already covers them.

## Consequences

- **Measured on the development panel** (dna-eval, against main before this
  change):

  | Measure | Fixes | Regressions |
  |---|---:|---:|
  | Primary | +40 | 0 |
  | Strict | +37 | 0 |
  | Haplotype IUPAC-compatible | +50 | 0 |

  - The strict measure, which counts uncovered truth alleles, rises from
    1,820 to 1,857 concordant samples, so this is an accuracy gain and not a
    reclassification.
  - No new DNA allele appears.
  - `315.1C` is now matched in 55 more samples.
  - Keeping substituted calls under the margin costs 9 samples: there the
    aligner writes two length changes (`309.1C 315.1C`) as two mismatches at
    the run ends instead of two insertions, so the exempt run's calls sit on
    a reference base of another base. Recognising such shifted run ends is
    left to a follow-up decision.
  - A cross-check with the Researcher session's separate consensus → call
    scorer on two 320-sample subsets of the panel gained only true calls
    (`309.1C`, `315.1C`), with no worse sample. An earlier build also let the
    consensus compose a substitution (`16183C`) from an exempt run call; the
    rule now keeps substituted calls under the margin.
- **`dna sample`.** Its eligibility changes the same way: an insertion or
  deletion of such a run is no longer `read_end`.
- **Versions.** `ReadEvidence` masks carry an onset (`MaskOnset`). The `core`
  plugin version rises to 4 and the `sanger` plugin version to 3.
- **Not covered.** Reads that fade without a post-repeat onset (weak tails),
  and runs a read loses for good, are still distrusted. They remain the
  largest HVS-II residual.

## Revision 2026-10-16: reference-matching bounds

The consensus decides a run's length only from a read whose run ends are
clean (ADR-0074). With the bounds outside the exemption, a read that the
exemption trusts across the run still could not show the run's length when a
bound sat within the margin, typically `T310` after a post-repeat phase loss
at the `A302` exit. The consensus now also trusts an exempt run's bound where
it matches the reference. A bound that reads another base keeps its own
trust, as does every substituted run call. `dna sample` eligibility is
unchanged: the bound adds no variant there.

**Measured on the development panel** (dna-eval, panel v4, against the
previous rule):

| Measure | Fixes | Regressions |
|---|---:|---:|
| Primary | +13 | 0 |
| Strict | +11 | 0 |
| Haplotype IUPAC-compatible | +14 | 0 |

- No new DNA allele and no new uncovered truth allele appear; `315.1C` is
  matched in 14 more samples. The draft triage outcomes do not change.
- The rule trusts 621 bound calls in 228 samples. Most lie on HVS-II run
  exits: `T310` (80 samples), `C303` (56), `A300` (54), `C317` (33), `C299`
  (32), `T321` (26). In HVS-I, `C16184` (17), `T16172` (12), and `A16170` (9)
  are the most frequent. No bound lies at `16182`, `16183`, or `16194`.

**Residual risk.** A bound that matches the reference can never create a
variant. If the read is wrong at the bound, it can hide a real alternate
there as a reference call. The risk is largest where an alternate at the bound
is common in the reviewed corpus (about 112,000 samples): `16172C` 22 %,
`16184` 3.4 %, `332` 1.4 %. For comparison, the corpus has `16183C` at 31 %,
`16182C` at 15 %, and `316` at 1.3 %. None of those three positions is a
trusted bound. On the development panel, one sample has a truth alternate at a
trusted bound: `16184Y`, a heteroplasmy that DNA did not report under either
rule. Draft triage already sends that sample to review (`minor_signal`). The
validation panel's haplotype gate, which counts a missed heteroplasmy too, and
the adjudication of its leaks are the checks on this risk.

A trusted bound can also oppose another read's clean alternate. One read
against one read gives the reference and `contested` (SRS-CONS-005), and any
other tie leaves the position undecided and `contested`. The draft triage
sends a contested region to review, so such a tie never passes. On the
development panel, 17 newly trusted observations in 13 samples oppose a clean
alternate at the same position: 12 at `310` (`T` against one `C`) and 5 at
`523` (`A` against deletions). None of them changes a sample's result or its
triage outcome, and no sample's contested count changes. The stretch decision
over run structure (ADR-0074) reads both reads as the same sequence aligned
differently.

**Versions.** The `consensus` plugin version rises to 3; the `core` plugin is
unchanged.
