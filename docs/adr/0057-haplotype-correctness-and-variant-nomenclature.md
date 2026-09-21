# ADR-0057: Separate haplotype correctness from variant nomenclature

- **Status:** Accepted
- **Date:** 2026-09-20
- **Related:** ADR-0009, ADR-0019, ADR-0047

## Context

DNA's end-to-end objective is to recover the correct sample variant profile from the
available AB1 traces. A variant profile is both a biological claim and a representation of
that claim relative to a reference.

Those two concerns are not identical.

In homopolymers, tandem repeats, and other sequence-equivalent contexts, the same resulting
haplotype can be described by different variant coordinates or by a different decomposition
into SNVs and indels. One system may report one right-aligned deletion while another reports
an SNV plus a deletion elsewhere in the same repeat. If applying both descriptions to the
same reference produces exactly the same resolved sequence, treating the descriptions as
different biological outcomes creates artificial false positives and false negatives.

ADR-0047 already defines DNA's internal alignment convention for repeat-equivalent
indels: retain an optimal alignment and choose one deterministic 3'/right-most mtDNA gap
placement. That decision makes DNA internally consistent, but it does not imply that an
external reviewer, legacy application, or future nomenclature standard must encode the same
haplotype with the same event boundaries.

Therefore validation, internal canonicalization, and external nomenclature must be separate
contracts.

## Decision

DNA adopts a **haplotype-first correctness model** for representation-ambiguous variant
comparisons.

### 1. Biological equivalence precedes textual representation

When two unambiguous variant descriptions are applied to the same exact reference and
produce the same resulting sequence, they represent the same resolved haplotype for
validation purposes.

Different raw coordinates, anchors, or event decompositions MUST NOT by themselves create a
biological false positive or false negative.

For example, conceptually:

~~~text
reference:
    ... repeat context ...

description A:
    SNV + downstream repeat deletion

description B:
    one canonical right-aligned deletion

if apply(A, reference) == apply(B, reference):
    biological outcome = equivalent
    representation = different
~~~

This rule is generic. It MUST NOT depend on a hard-coded mtDNA position, named repeat,
poly-C tract, sample identifier, or observed validation hotspot.

### 2. DNA still needs one deterministic internal representation

Haplotype-first validation does not make representation arbitrary inside DNA.

For identical reference, evidence, configuration, and scientific method, DNA SHOULD emit
one deterministic canonical internal representation. ADR-0047 remains authoritative for
repeat-equivalent alignment gap placement and currently establishes the right-most/3'
mtDNA convention.

Internal consistency is important for:

- stable sample aggregation;
- reproducible JSON;
- deterministic testing;
- provenance and diffability;
- downstream evidence joins.

A representation difference between DNA runs under the same declared method remains a
regression even when both representations reconstruct the same sequence.

### 3. Nomenclature is a separate boundary

External nomenclature rules may eventually require a representation that differs from
DNA's internal canonical event structure.

Such conversion MUST be an explicit, versioned nomenclature/reporting layer. It MUST NOT
silently alter:

- read evidence;
- alignment evidence;
- the reconstructed sample haplotype;
- variant eligibility;
- phase evidence;
- sample reconciliation.

A future EMPOP-, forensic-, HGVS-, VCF-, or other convention therefore belongs at a named
serialization or interpretation boundary rather than inside signal measurement or
alignment.

### 4. Validation compares representation in layers

Reviewer/profile comparison MUST proceed from stricter to more semantic checks:

~~~text
exact event identity
        ↓
single-event sequence equivalence
        ↓
unambiguous multi-event haplotype equivalence
        ↓
unmatched biological disagreement
~~~

Representation-equivalent matches MUST remain visible as representation disagreements, but
they MUST NOT inflate biological FP/FN counts.

The validation artifact MUST preserve the original source events on both sides so that an
equivalence decision is auditable.

### 5. Multi-event equivalence must be conservative

A multi-event comparison MAY collapse N reviewer events against M DNA events only when:

- all compared reviewer events resolve to one unambiguous alternate;
- compared events can be applied to the reference without contradictory overlapping edits;
- applying the complete reviewer group and complete DNA group to the same reference
  produces exactly the same resulting sequence;
- no proper subgroup already explains the same equivalence;
- no event participates in multiple competing minimal equivalence groups.

If those conditions are not met, validation MUST leave the events unmatched rather than
guessing an equivalence.

This prevents haplotype equivalence from becoming a permissive fuzzy-matching rule.

### 6. Ambiguous biological mixtures remain distinct

A single resolved final sequence is not sufficient to represent every biological state.

IUPAC mixtures, heteroplasmy, multiple haplotypes, unresolved phase, and other ambiguous
states MUST retain their explicit evidence semantics. Haplotype-equivalence collapsing is
only valid where the compared event group resolves to a definite sequence outcome.

DNA MUST NOT force an ambiguous reviewer or DNA observation into one sequence merely
to obtain a match.

### 7. Final-sequence validation is a correctness fallback, not a production shortcut

When nomenclature differs and DNA cannot directly reproduce an external representation,
validation MAY use exact reconstructed sequence/haplotype equality to establish biological
equivalence.

This does not mean production may emit arbitrary variant decompositions. Production still
needs deterministic internal canonicalization, and a future nomenclature layer SHOULD
translate that canonical result into the chosen external convention where required.

## Consequences

- Repeat-region representation differences no longer create artificial biological FP/FN
  when both sides encode exactly the same resolved haplotype.
- DNA's internal right-aligned representation remains deterministic and authoritative
  for current production behavior.
- Reviewer or external nomenclature differences remain visible and auditable instead of
  being rewritten into the ground-truth source.
- Variant-profile metrics can distinguish biological calling errors from representation
  disagreements.
- Future nomenclature support can evolve independently from signal processing, alignment,
  phase interpretation, and variant inference.
- Validation code must compare source events conservatively and preserve provenance rather
  than normalizing one source destructively into the other.
- A stable biological result can remain correct even when reporting conventions change.

## Relationship to phase-aware calling

Phase-aware work is judged by whether it improves the final sample biological variant
profile.

This ADR prevents repeat/nomenclature representation differences from being mistaken for
phase-related false positives or false negatives. Only disagreements that remain after
representation-equivalence handling should enter biological error attribution and later
phase-context analysis.

Thus the intended validation dependency is:

~~~text
reviewer/source variant descriptions
        +
DNA canonical variant descriptions
        ↓
representation-aware haplotype comparison
        ↓
biological missing / extra disagreements
        ↓
phase and other error attribution
~~~

## Validation requirements

Validation tooling implementing this ADR MUST cover at least:

- exact event identity;
- repeat-shifted one-event indel equivalence;
- one-event versus multi-event descriptions that reconstruct the same sequence;
- multi-event versus multi-event equivalence;
- a proper equivalent subgroup that prevents over-grouping;
- competing minimal groups that remain unmatched;
- an ambiguous/IUPAC event that MUST NOT be collapsed into a resolved haplotype;
- a genuinely different final sequence that remains a biological disagreement.

## Non-goals

This ADR does not define a new production variant caller, phylogenetic realignment,
heteroplasmy genotype model, EMPOP/HGVS/VCF renderer, phase-state classifier, or universal
variant-normalization standard.

It also does not claim that final sequence alone captures every mixed biological state.
