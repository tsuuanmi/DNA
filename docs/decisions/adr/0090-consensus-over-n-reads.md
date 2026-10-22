---
id: ADR-0090
status: accepted
date: 2026-10-21
family: [consensus]
related: [ADR-0069, ADR-0071, ADR-0074, ADR-0080, ADR-0083, ADR-0089]
---

# ADR-0090: Consensus over N reads per site and strand

<!-- docs:begin lifecycle -->
- **Status:** Accepted
- **Date:** 2026-10-21
- **Related:** [ADR-0069](0069-plugin-first-modality-core-post-calling.md), [ADR-0071](0071-run-length-edits-need-resolved-run-ends.md), [ADR-0074](0074-run-structure-consensus.md), [ADR-0080](0080-reference-only-calls.md), [ADR-0083](0083-set-aside-dissent.md), [ADR-0089](0089-insertions-need-both-strands.md)
<!-- docs:end lifecycle -->

<!-- docs:begin backlinks -->
- **Governs:** [design/consensus](../../design/consensus.md), [design/sample-evidence](../../design/sample-evidence/README.md), [design/sample-evidence/read-registry](../../design/sample-evidence/read-registry.md), [engineering/improvement-method](../../engineering/improvement-method.md), [proposals/roadmap](../../proposals/roadmap.md), [reference/consensus](../../reference/consensus.md), [reference/review](../../reference/review.md), [reference/sample-evidence/reads-coverage](../../reference/sample-evidence/reads-coverage.md), [requirements/consensus](../../requirements/consensus.md), [requirements/sample-evidence/core](../../requirements/sample-evidence/core.md), [validation/evaluation-harness](../../validation/evaluation-harness.md), [validation/regression-set](../../validation/regression-set.md)
<!-- docs:end backlinks -->

## Context

The user decided on 2026-10-21 that DNA must never assume exactly two
traces at a position. The consensus and every crate must handle any number
of reads per site and per strand. Re-sequenced reads exist today, and NGS
will bring many reads per strand later.

An audit of `main` (`ba0dc0a2`) found that the core already loops over
every read: no rule holds a fixed forward–reverse pair or one optional read
per strand. The assumptions are narrower:

- **Silent same-strand dissent.** The review's `strands_disagree` fires
  only when a supporting and an opposing read lie on different strands, and
  `read_disagrees` covers only supporting and set-aside reads. A site that
  two forward reads decide against a third forward read got no check, only
  `one_strand` as context.
- **Read identity by file name.** The consensus keyed its reads by the
  input's file name, the review joined consensus sites to reads by it, and
  the reports failed with `DuplicateReadName` when two inputs shared a file
  stem. A re-sequenced trace saved under the same name in another folder
  failed the sample, or two reads collapsed into one map entry. SRS-CONS-002
  already forbade the consensus to use file names.
- **One read standing in for one strand.** The review's `run_deletion`,
  `run_indel`, `edge_deletion` and `reference_second` gate on exactly one
  supporting read, and `long_run` counts reads, where the question is
  whether the support lies on one strand. Two same-strand reads with the
  same artefact remove the check.
- **Counts that ignore strands.** The decision rule (decision rule version
  1) is a read-count majority, so same-strand reruns vote as independent
  evidence. `set_aside_dissent` (ADR-0083) lets one clean supporting read
  set aside every agreeing opposing read, whatever their number.
- **Pairs of one view per strand.** The noise calibration of τ
  (`noise.rs::cross_strand`, PROP-0014) uses exactly one view per strand,
  so a re-sequenced sample adds nothing. Rule A (ADR-0074, SRS-CONS-013)
  takes the largest spread over all forward–reverse pairs, which grows
  stricter with more reads. A strand's length is clear only when every read
  of it is clear and equal, so one outlier breaks the strand.
- **Evaluation limits.** `dna-eval` stores per-position strand counts as
  `u8`, and reads a primer only from `HV<n><F|R>`, so `HV3-350F` gets none.

On the development panels, a region and strand read twice or more is still
rare: 62 of 10,694 current samples of development panel v5 and 18 of 17,307
current samples of readiness-v2 (by the traces' primer names), 57 of 5,716
and 23 of 2,833 earlier ones.

## Decision

**Principle.** Every rule over a site's reads works over sets of reads
grouped by strand, of any size. No rule names "the other read", pairs reads
as forward and reverse, or tests `len() == 1` or `len() == 2` where it asks
a question about strands. A read is known by its input's content SHA-256;
its file name is for display only.

**Read and strand.** A read is one admitted input trace. Its identity is
the SHA-256 of its content; its strand is its placement orientation. Two
inputs of one content are one trace given twice: the sample fails typed
(`DuplicateTrace`), whether or not a read is admitted, rather than counting
it twice or dropping a copy silently. A re-sequenced trace is a new file
with new content, so it is a new read even under the same file name. Its
reviewer-facing name is its file stem; reads that share a stem each gain
`@` and the shortest prefix of their SHA-256, at least eight characters,
that tells them apart (SRS-SAMPLE-004).

**One read and one strand.**
- *One read*: exactly one read supports the decision. It describes
  evidence for the reviewer (`one_read`).
- *One strand*: every supporting read lies on one strand, however many.
  Independence is a question about strands: reads of one strand share its
  systematic errors (dephasing after a C tract, a dye blob at a primer's
  start), so a rule that asks whether the support is independent asks
  whether it is on one strand, never whether there is one read.

**Vote unit (target of increment 3).** Clean observations vote per strand
first, then across strands:
1. Within a strand, a sequence that a strict majority of the strand's clean
   reads shows is the strand's vote. A strand whose clean reads tie (k
   against k) casts no vote, and its dissent is reported.
2. Across strands, two strands that vote the same sequence decide it. One
   strand voting alone decides it on one strand. Two strands that vote
   different sequences are contested: the reference when one of them is
   the reference, otherwise no sequence. This is decision rule version 1's
   one-against-one rule, with strands in place of reads.

Until increment 3 ships, decision rule version 1 stays: a read-count
majority across all reads, one read against one read to the reference.

**Same-strand dissent (step 1).** A site whose supporting and opposing
reads all lie on one strand is a review check, `same_strand_dissent`
(rank 0): clean reads of one strand disagree, and their count alone does
not settle the site. It lists every supporting, then every opposing read.
`strands_disagree` keeps the case where a supporting read lies on another
strand than an opposing one. Both labels come from one rule over all the
site's reads, and a clipped site recomputes its label from the opposing
reads that remain.

**Plan.** The work ships in increments, each measured against `main` with
the `dna-eval diff` gate before it merges:
1. Same-strand dissent check and read identity (this decision's step 1).
2. One-strand gates by strand: `run_deletion`, `run_indel`,
   `edge_deletion` and `reference_second` gate on all support lying on one
   strand rather than on one read, and their labels carry `reads: Vec`;
   `long_run` counts strands.
3. Per-strand voting as above, and a count-aware `set_aside_dissent`: a
   dissent is set aside only by at least as many clean reads of the other
   strand.
4. Pooled per-strand length statistics for the two-strand rule and for τ
   (`noise.rs::cross_strand`), and a per-strand dominant length with a
   dissent count in place of "every read clear and equal".
5. `dna-eval`: `u32` strand counts, and primers from `HV<n>-<pos><F|R>`.

**Roadmap, not now.** The kernel's read evidence and the consensus still
carry Sanger vocabulary (in-phase bounds, phase mixtures, second bases,
anchored ends, run-length evidence; the review reads `PhaseState`). An NGS
modality needs it moved behind the modality boundary of [[ADR-0069]]. That
is future work on the roadmap; this decision does not start it.

The rules read only the core's trust levels, the reads' strands and their
identities, never a modality's reasons or a file name ([[ADR-0069]],
[[ADR-0080]]). The run-structure rules ([[ADR-0071]], [[ADR-0074]]) and the
two-strand insertion rule ([[ADR-0089]]) keep their meaning; increments 3
and 4 restate them over sets of reads.

## Consequences

Step 1: plugins `consensus` 21 and `review` 30. Measured against `main`
`ba0dc0a2` with the evidence cache and the frozen re-check list:

- **Development panel v5, current period:** the label adds no check.
  It marks 7 items, all contested same-strand ties already `undecided`
  checks, whose rank moves from 2 to 0. None lies over a truth difference.
  Clear 8,066 → 8,066, non-het unreviewed blocks 4 → 4.
- **Development panel v5, earlier period:** 21 such items, the same kind;
  no check added, non-het unreviewed 12 → 12.
- **Readiness-v2, current period:** 4 items, the same kind; clear
  14,117 → 14,117, non-het unreviewed 1 → 1. Earlier period: 8 items,
  non-het unreviewed 11 → 11.

The `dna-eval diff` gate passed on both panels with no sample changed; the
scores and unreviewed blocks are byte-identical, and the artefact bench is
unchanged (53 pass, 2 known). No current sample has a decided site whose
clean support and dissent lie on one strand, so the label is a guard for
re-sequenced and NGS samples. No panel sample has two traces of one file
stem, so read identity changes no output today; it removes the failure.

Step 2: plugin `review` 32. `run_deletion`, `run_indel`, `edge_deletion`
and `reference_second` gate on all support lying on one strand, the first
three labels carry `{strand, reads}`, and `long_run` counts strands.
Measured against `main` `95f2edf2` with the evidence cache and the frozen
re-check list, the `dna-eval diff` gate passed on both panels and the
artefact bench is unchanged (53 pass, 2 known). Every added check is a
`reference_second`; none in the current period lies over a difference
between DNA and the truth, and the non-het unreviewed blocks are
unchanged.

- **Development panel v5, current period:** 8 new checks in 8 samples,
  all re-sequenced, 5 of them in the HVS-II C tract where the forward read
  loses phase and a reverse rerun stands beside the first reverse read;
  one existing check gains the label. Clear 8,358 → 8,353; non-het
  unreviewed 4 → 4.
- **Readiness-v2, current period:** 2 new checks at samples without a
  rerun: an HV1 forward read reads through the origin into HVS-II, so its
  tail and the HV2 forward read are the support of one strand, and the
  tail shows the reference second (`249DEL`, `151T`, DNA right). Clear
  14,408 → 14,406; non-het unreviewed 1 → 1.
- **Earlier periods:** 7 new checks on each panel, 4 (development panel
  v5) and 5 (readiness-v2) at re-sequenced samples and the others by HV1
  read-through. One now covers a truth heteroplasmy at 204 that was
  unreviewed (het-only unreviewed 194 → 193), and one a difference where
  the truth leaves HVS-II unanalysed. Clear −4 and −1.

So a second read of a strand is not only a rerun: an amplicon that reads
into the next region adds one too, and the rule treats it alike.

## Alternatives

- **Generalise `read_disagrees` to opposing reads.** `read_disagrees`
  shows a read only when the review itself reads it cleanly and bounded
  over the whole item. An opposing read that the consensus counts as clean
  but the review cannot bound would stay silent, the gap this decision
  closes. It also mixes reads the consensus explained away with reads it
  outvoted, which blurs each label's precision.
- **Widen `strands_disagree` to any opposing read.** The label would lose
  its cross-strand meaning and its measured precision, and reviewers and
  `dna-eval` group checks by label.
- **Deduplicate identical inputs.** Counting a repeated input twice fakes
  support, and dropping a copy silently hides an operator error; the typed
  error says what happened.
- **Key reads by their input index.** An index depends on argument order.
  The SHA-256 is stable, already orders the sample's reads, and keys the
  evidence cache.
- **Keep file-name identity and refuse duplicates.** It blocks
  re-sequencing folders that reuse a trace's file name, and breaks
  SRS-CONS-002.
- **Pool all reads without strands.** Reads of one strand are not
  independent evidence, so a pure read count lets reruns of one strand
  outvote the other strand.
