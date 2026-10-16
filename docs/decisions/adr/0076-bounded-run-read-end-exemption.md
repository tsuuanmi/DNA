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
   - in the consensus, the run's calls are clean for its run observations.

   A substitution is never exempt, and the run's two bounds keep their own
   trust. A bound shows where the run ends but never a clean base, so a
   substitution at a run exit still needs a read that is clean there without
   the exemption.
2. **Bounded run.** The run (two or more calls of one canonical base) has a
   resolving call on each side (ADR-0071). The bounds and every run call are
   unmasked.
3. **Post-repeat edge.** On each side, either the nearest uninformative call
   is beyond the margin, or it is a masked call whose segment the modality
   places right after a repeat run (`MaskOnset::AfterRun` in `ReadEvidence`).
   A read's physical end, the trim interval, or any other masked segment
   never exempts.
4. **Guard.** At least two informative calls separate each bound from such an
   edge. This is a structural minimum, not a tuned value: guards of 0, 2, and
   3 gave identical results on the development panel.

The modality reports the onset as part of the kernel `ReadEvidence` contract
between modality and core (ADR-0069): each mask carries a `MaskOnset`,
`AfterRun` or `Other`. The Sanger plugin sets `AfterRun` exactly where
callability marks a segment `after_repeat` (the `post_homopolymer` reason).
The core does not interpret reason labels. The exemption claims nothing
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
  | Primary | +49 | 0 |
  | Strict | +46 | 0 |
  | Haplotype IUPAC-compatible | +60 | 0 |

  - The strict measure, which counts uncovered truth alleles, rises from
    1,820 to 1,866 concordant samples, so this is an accuracy gain and not a
    reclassification.
  - No new DNA allele appears.
  - `315.1C` is now matched in 67 more samples.
  - A cross-check with the Researcher session's separate consensus → call
    scorer on two 320-sample subsets of the panel gained only true calls
    (`309.1C`, `315.1C`, and one `16183C` that the consensus composes from the
    run's calls), with no worse sample.
- **`dna sample`.** Its eligibility changes the same way: an insertion or
  deletion of such a run is no longer `read_end`.
- **Versions.** `ReadEvidence` masks carry an onset (`MaskOnset`). The `core`
  plugin version rises to 4 and the `sanger` plugin version to 3.
- **Not covered.** Reads that fade without a post-repeat onset (weak tails),
  and runs a read loses for good, are still distrusted. They remain the
  largest HVS-II residual.
