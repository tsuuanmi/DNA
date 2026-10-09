# ADR-0068: Variant opposition evidence

- **Status:** Accepted
- **Date:** 2026-10-09
- **Related decisions:** [ADR-0009](0009-biological-semantics.md),
  [ADR-0033](0033-variant-support-topology.md),
  [ADR-0065](0065-result-comparison-downstream.md),
  [ADR-0067](0067-signal-derived-read-callability.md)

## Context

Sample evidence publishes, for every normalized variant, the reads that
support it and their support topology (ADR-0033). It does not say which other
reads observed the same place and saw something else. ADR-0033 deferred
"opposition" to a separate decision, because before ADR-0067 a covered position
was not necessarily an informative one.

Read callability now publishes callable reference segments per read
(ADR-0067), which tell a covered and callable position apart from a masked one.
Most false calls left on the local corpora come from one read whose call
another read contradicts at the same callable position. A downstream consumer
can only see this today by joining `locus_differences[]`, which has no records
for insertions, with each read's callable segments.

DNA stays an evidence producer. Sample consensus, majority votes, and
sample-level verdicts remain out of scope, and nothing here is specific to a
target or organism.

## Decision

1. Every sample variant publishes `opposition`: the admitted reads, in
   read-registry order, whose callable reference segments cover the variant's
   evidence span and that have no support record for this variant, with their
   forward and reverse counts.
2. The evidence span is the anchored reference allele: the position for an SNV
   and the anchor plus the deleted bases for a deletion. For an insertion it is
   the anchor and the base after it, so that both sides of the inserted bases
   are observed.
3. Opposition is evidence, not a vote, verdict, weight, or confidence. A read
   that supports a different event at the same place opposes this one, a
   forward and a reverse read can share an assay artifact, and a read can
   oppose a true minor allele. DNA does not combine support and opposition into
   a sample decision.
4. Opposition is derived only from published evidence (support records,
   callable reference segments, orientation), so contract validation can
   recompute it.
5. `dna.sample_evidence/v10` is revised in place while unreleased (versioning
   policy); the analysis and basecall contracts are unchanged.

## Alternatives

- **Exclude contradicted variants from notation:** a sample-level verdict that
  ADR-0009 and the roadmap keep out of DNA. On held-out data one in five
  opposed calls is true.
- **Count reference-observing reads only:** needs a reference/competing state
  per read and per variant. Insertions have no locus record, and a competing
  event is already visible as that read's own support elsewhere.
- **Leave the join to consumers:** every consumer would re-derive spans and
  callable coverage, and get insertions differently.

## Consequences

- **Breaking (unreleased, revised in place):** each `variants[]` entry of
  `dna.sample_evidence/v10` requires `opposition`.
- Reviewers and downstream adjudication see contradicting reads directly.
- Measured outside the repository (ADR-0065), on top of revision `e470564`,
  in aggregate:
  - On the local reviewed corpus, 4 eligible SNVs had an opposing read, and
    all 4 were false calls; none of 316 correct SNVs was opposed.
  - On a held-out set of 320 samples from a separate reviewed corpus that no
    threshold was tuned on, eligible SNVs with an opposing read were 12 false
    and 3 true calls, against 2 534 correct unopposed SNVs.
- Variant calling and notation results are unchanged.

## Supersession

This record takes the opposition decision that ADR-0033 deferred, for
callable coverage only. Gap handling, contributor admission, and calibration
for a sample-level policy remain open.
