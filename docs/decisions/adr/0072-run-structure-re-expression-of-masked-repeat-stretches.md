# ADR-0072: Run-structure re-expression of masked repeat stretches

- **Status:** Accepted
- **Date:** 2026-10-10
- **Related decisions:** [ADR-0029](0029-profile-aware-gotoh.md),
  [ADR-0047](0047-canonical-right-aligned-mtdna-gaps.md),
  [ADR-0065](0065-result-comparison-downstream.md),
  [ADR-0067](0067-signal-derived-read-callability.md),
  [ADR-0069](0069-plugin-first-modality-core-post-calling.md),
  [ADR-0071](0071-run-length-edits-need-resolved-run-ends.md)

## Context

After ADR-0071, the HVS-II poly-C window (rCRS 303–315, `CCCCCCC T CCCCC`) is
the largest remaining Sanger error. On 320 held-out samples, 18 of 316
windows differed from reviewer truth (ADR-0065 comparison).

In many of those reads the calls spell the correct run lengths. A typical
forward read reads `C×8` in phase, then a dephased `T`, then `C×6`. Dephased
calls are anchoring masks (ADR-0067): they keep their base and their mixed
locus profile. The optimal Gotoh path then puts the 8th C on T310 and inserts
the T, or merges both length changes into one two-base gap next to a
substitution. Merging is cheaper: one gap opening plus one mismatch costs less
than two gap openings.

Variant extraction then reports the in-phase substitution `310T>C` as
eligible. It excludes the insertion because its inserted call is masked. The
read's eligible subset names a haplotype the read does not show (`310C`), or
nothing. On reverse reads the same split hides the length change inside a
mixed insertion, which `run_boundary` (ADR-0071) does not recognize.

Clean consensus reads are split the same way. Every edit there is eligible,
though, and the nomenclature window rebuilds the haplotype, so `dna call` is
correct.

## Decision

After orientation and placement are selected and the read passes the
callable-base and identity gates, the core re-expresses masked repeat
stretches of the selected alignment.

1. **Stretch.** Each maximal block of columns that are not exact matches is
   grown to whole runs of one reference base. A run may be a single base. The
   stretch grows run by run until each side ends at an anchor: an exact
   canonical match whose base the read's adjacent call inside the stretch
   does not repeat. A stretch never crosses the circular origin and never
   spans more than 64 reference bases. A stretch with no anchor on either
   side is left unchanged.
2. **Condition.** The stretch is rewritten only when all of these hold:
   - every read call in it is canonical;
   - at least one of them is masked;
   - the read's run bases, in order, equal the reference's;
   - every changed run is longer, or every changed run is shorter, than the
     reference run.

   Mixed lengthening and shortening are also explained by substitutions at
   run boundaries, such as `16182C 16183C`, so they are kept.
3. **Rewrite.** Each run becomes its paired matches followed by one length
   edit at its 3' end: insertion columns when the read run is longer,
   deletion columns when it is shorter. This is the right-most placement of
   SRS-ALN-012.
4. **What is kept.** The placement, its span, its score (the optimum of the
   selected path), and every column outside the stretch. Metrics are
   recomputed from the rewritten columns.

Variant extraction, `read_end`, `run_boundary`, and mask reasons then work on
one length edit per run. The rule needs only `ReadEvidence` and the reference.
It is modality-neutral, target-neutral, and has no parameter. It never changes
a read without masked calls, so consensus calls from `dna call` stay
identical.

## Alternatives

- **Aligning anchoring calls by a one-hot profile of their base.** This keeps
  the dynamic-programming optimum but was much worse on the 320-sample set:
  - precision fell from 0.9894 to 0.9720;
  - 5 samples failed with ambiguous placements;
  - about 50 false `310C` calls appeared.

  A dephased call's primary is often a shifted ladder, and its mixed profile
  is the better anchor.
- **Gating on any uncertain profile instead of a masked call.** This added
  two false insertions outside the poly-C windows and fixed nothing more.
- **Ungated re-expression.** It would also rewrite consensus reads, whose
  split edits are all eligible and already named correctly.
- **Profile-window knowledge in the core, or a post-calling repair.** The
  first breaks ADR-0069's boundary. The second lacks the mask evidence that
  decides whether a read reads a run at all.
- **Extending `run_boundary` to mixed edits.** The edit set would still
  describe the wrong haplotype.

## Consequences

- **320 held-out samples:**
  - 6 more true HVS-II calls (`309.1C` ×5, `309.2C`) and one false `310C`
    fewer;
  - no true call is lost;
  - precision rises from 0.9894 to 0.9897 and recall from 0.9800 to 0.9816.
- **160-sample consensus ceiling:**
  - `dna call` outputs are byte-identical;
  - Sanger precision against consensus rises from 0.9914 to 0.9919, and
    recall from 0.9871 to 0.9882.
- **Local Sequencher set:** unchanged.
- **Published alignment summaries** (`identity`, `gap_opens`, and rarely
  `callable_reference_segments`) change for reads with a rewritten stretch.
  So do their observed variants, locus differences, and overlaps. On the
  320-sample set, 84 documents changed. Outside the repeat windows, only
  observed variants outside the target regions changed.
- **INV-ALN-001 and INV-ALN-003** gain an explicit exception (INV-ALN-007,
  SRS-ALN-016). The selected alignment's columns are no longer always a
  maximum-score traceback, though its placement and score are.
- **Unreadable length heteroplasmy** (a read that loses T310 among several
  length populations) still gives no call. Estimating a dominant length needs
  its own decision under ADR-0009.
