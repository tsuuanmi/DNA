---
id: ADR-0091
status: proposed
date: 2026-10-21
family: [placement]
related: [ADR-0029, ADR-0047, ADR-0072, ADR-0090]
---

# ADR-0091: A re-expressed repeat stretch keeps the read's own anchor base

## Context

After the per-run `reported` rule of the review (SRS-API-020, #226), 312
`window@HV2` items of development panel v5's current period stayed checks
because DNA did not show HVS-II's anchor `T310`. In every one of them a
reverse read reads run 2 (`311-315`) in phase at DNA's length, and also
reads `T310` cleanly, yet its alignment places a C at 310 and the T on an
inserted column.

The read is `C9 T C6` against the reference's `C7 T C5`. Two causes stack.

1. **The affine optimum merges the two length changes.** One gap opening
   (−10) costs more than a mismatch and a lost match (−5 − 3), so the path
   `310C` plus an inserted `TCC` scores 2 more than `309.1C 309.2C 315.1C`.
   [[ADR-0072]] already re-expresses such a stretch, but only when it holds a
   masked call. The reverse read's calls there are all clean.
2. **The read's last in-phase call reads `302` as C.** It is the bleed of
   the long G run in the read's orientation (`phase_shift`, its profile
   favouring the reference's A). The stretch then grows over the `A3`
   run: `A2 C10 T C6` against `A3 C7 T C5`. One run is shorter and another
   longer, a "mixed" change that ADR-0072 leaves as aligned. The masked
   calls before it, the read's dephased end, are a separate disturbance;
   their stretch shares only anchor columns with this one, but the two were
   merged.

The review then reads the decoded run lengths of the corrected reverse read
on the wrong column. A run's inserted bases sit at its 3' end
(SRS-ALN-012), which is the start of a reverse read's run. The decode of
run 2 starts on `315.1` and was credited to the next aligned call in the
read's order, `316`. That lies past the window, so the decode was dropped.

## Decision

Run-structure re-expression (SRS-ALN-016) changes in three ways.

1. **Gate.** A stretch also qualifies by a substitution, where every run
   whose length changes is two reference bases or more. In a repeat run, a
   substitution beside a gap is a length change that the gap-open cost
   merged into the gap. A single base the read repeats is no slipped run,
   and the substitution may be the read's own error. An AC-repeat read that
   reads an A as C, beside a true `524.1A 524.2C`, keeps that substitution.
   A stretch without a masked call or a qualifying substitution is never
   changed.
2. **Run boundaries.** Where one run is shorter and its neighbour longer,
   the boundary between them moves by substitution, from 5' to 3', as far
   as both changes allow (the `16182C 16183C` reading). When the remaining
   changes are all lengthening or all shortening, the stretch is rewritten
   as before: one length edit at each changed run's 3' end. Otherwise it is
   kept.
3. **Adjacent stretches.** Two stretches whose anchors are each other's
   first and last columns share no column between their anchors. They are
   re-expressed separately.

The read's calls and the reference positions are never changed, only how
they pair. The placement, its span and score, and every column outside a
rewritten stretch are kept, as under ADR-0072.

The review credits a decode whose first call sits on an inserted column to
the reference position the insertion follows, on either strand
(SRS-API-020 *Credit*).

No new parameter is added, and nothing is tied to a position. The rule
applies to any anchored repeat window and to each read on its own, so it
holds over any number of reads per strand ([[ADR-0090]]).

## Alternatives

- **Re-express every unmasked stretch.** This was measured on 2026-10-16
  (roadmap): 14 regressions and new false alleles. On today's panel, the
  ungated substitution arm hid one truth `524.1A 524.2C`, in a current
  sample, behind a forward read's own A→C error. The repeat-run condition
  keeps that alignment.
- **Change the gap scores.** A cheaper gap opening would move every
  placement and every score. The merge is a correct optimum of the scoring;
  the representation is what is wrong.
- **Normalise in the consensus or the review.** Both would repair one read
  representation downstream of placement. Placement is the step that owns
  sequence-equivalent edits ([[ADR-0047]]), and every later step reads the
  same columns.
- **Mask the bleed call at `302`.** A masked call keeps its base, so the
  stretch would still be mixed. It would also drop a call that does
  support the read's reading of `302`.

## Consequences

`dna-eval diff` against `main` 95f2edf2 used the evidence cache, the frozen
re-check list, and the default truth. The period is current
(≥ 2026-01). **GATE PASS** on both panels. Unreviewed non-het samples are
unchanged.

| Panel | Clear | window@HV2 checks | Unreviewed non-het | Passed het only | Passed length note only | Passed truth N |
|---|---|---|---|---|---|---|
| dev v5 | 8,358 → 8,813 | 1,054 → 473 | 4 → 4 | 276 → 288 | 14 → 24 | 1 → 7 |
| readiness v2 | 14,408 → 15,159 | 1,659 → 847 | 1 → 1 | 622 → 662 | 33 → 57 | 26 → 34 |

- **Split.** On dev v5 the credit fix alone gives clear 8,620 and 719
  window checks. The alignment change gives the rest.
- **The anchor.** The census found 337 dev-v5 samples, in all periods,
  whose window check had run 2 `undecided_call` and 310 not callable. Of
  these, 205 now show `T310`, and 191 lose their window check.
- **Passed truth N.** The new truth-N passes are `573` blocks at the end
  of HVS-III, where the truth writes `N` or stops.
- **New passed blocks at 303-316** are length notes (never gated): the
  truth reads the convention's length or DNA's length less one, and
  `other` is unchanged. Four more are on the re-check list (DNA `C9`
  against the truth's `C8`).
- **Heteroplasmy.** No current HVS-II heteroplasmy is newly hidden. The
  rise in passed het-only samples is samples that now clear while carrying
  a truth heteroplasmy elsewhere, which `main` already left unflagged. In
  the earlier period, two truth `309Y` per panel are newly passed.
- **Earlier period.** On dev v5, passed non-het rises from 3 to 4:
  2534275 now clears. Its non-het block is the artefact bench's
  `truth_error` `16183C`.
- **Artefact bench.** 53 pass and 2 known, as before.
- **`dna call`.** Consensus records split the same way are now
  re-expressed as well. Their variants already matched through the
  nomenclature window, and their alignment summaries may change.
- **Not fixed.** 132 samples still lack the anchor. In 126 of them the
  reverse read's masked end at 290-301 is read one base shifted (`ccccaaac`
  over `CACCAAA`). Its runs then differ from the reference's, so the
  stretch is kept and 310 stays not shown. This is a limit of placing a
  read's dephased end, not a representation choice; those items stay
  checks.
