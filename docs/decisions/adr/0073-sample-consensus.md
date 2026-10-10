# ADR-0073: Sample consensus as the adjudication layer

- **Status:** Accepted
- **Date:** 2026-10-10
- **Related proposal:** [PROP-0003](../../proposals/0003-sample-consensus.md)
- **Related decisions:** [ADR-0023](0023-evidence-derived-read-placement-and-sample-boundary.md),
  [ADR-0033](0033-variant-support-topology.md),
  [ADR-0035](0035-locus-support-topology.md),
  [ADR-0041](0041-mean-nucleotide-evidence-profiles.md),
  [ADR-0065](0065-result-comparison-downstream.md),
  [ADR-0068](0068-variant-opposition-evidence.md),
  [ADR-0069](0069-plugin-first-modality-core-post-calling.md),
  [ADR-0071](0071-run-length-edits-need-resolved-run-ends.md)

## Context

Sample evidence stops before adjudication. ADR-0023, ADR-0033, ADR-0035 and
ADR-0041 call it *pre-consensus* evidence and leave the consensus to a later
layer. Its notation is the union of every read's eligible calls.

Reviewers build a consensus sequence and call variants from it, and the core
reproduces their calls from their consensus sequences. Downstream
measurements (ADR-0065) found two problems with the union:
- a call that one read supports and another callable read opposes is the
  largest tractable source of false calls;
- one read that cannot be placed used to fail a whole sample.

## Decision

A modality-neutral **consensus** module in the core (`dna-core::consensus`,
plugin `consensus`, family `core`) adjudicates between a sample's admitted
reads. Sample evidence stays as it is, as the evidence the consensus
consumes.

**Observations.**
- Every read observes each reference position it covers, and the inserted
  bases at each junction between two covered positions.
- An observation is *clean* when its calls are unmasked, trusted
  (`read_end`), and raise no support veto.
- A base shows where a neighbouring run ends (ADR-0071) when it is canonical,
  informative, and unmasked or anchoring.
- Only clean observations decide.

**Stretches.** Each place where a read cleanly differs from the reference is
grown to whole reference runs on both sides. Overlapping or touching
stretches merge, and a stretch never crosses the origin. Unclean differences
decide nothing and do not widen a stretch. A stretch is decided in this order
of preference:

1. **As a whole.** A read observes the whole stretch cleanly with a run
   structure other than the reference's, for example a substitution inside a
   run. Each such read contributes its sequence over the stretch.
2. **Run by run.** Reads with the reference run structure each contribute the
   length of a run, when the run's calls are clean and the bases on both
   sides show where it ends. An unresolved call (`N`) removes that evidence
   from the bases next to it. One read can decide one run and another read
   the next.
3. **Position by position.** This applies only when no read shares the
   reference run structure, or when run-by-run decides nothing. Reads without
   a clean indel inside the stretch contribute the base their alignment
   places on each position.

Positions outside stretches are decided one by one.

**Rule, version 1.** Over clean observations:
1. all agree: that sequence;
2. a strict majority wins;
3. one read against one read: the reference sequence, marked `contested`;
4. any other tie: no sequence, marked `contested`;
5. no clean observation: `unresolved`.

**Assembly.**
- Undecided positions are written as `N`.
- Segments are contiguous decided positions. A gap of at most
  `alignment.minimum_callable_bases` undecided positions stays inside a
  segment as `N`, while a longer gap splits it.
- A segment with fewer resolved bases than that minimum is dropped, because
  the core could not place it.

**Output.**
- **`dna consensus`** writes `dna.consensus/v1` with the segments, the
  notable sites and the reads behind them, plus a FASTA of the segments.
- **Variants** come from the existing `dna call` over that FASTA. There is
  one variant caller, and it reads the consensus as reviewed sequences.

The rule is versioned through the plugin's method version, and any change to
it raises that version.

## Alternatives

- **Rename `sample` to `consensus`.** Rejected, because sample evidence makes
  no decisions.
- **Vote per position and junction.** Rejected. Reads place one repeat-length
  change at different junctions, or split it into a substitution next to an
  insertion. Per-site votes then broke HVS-II haplotypes.
- **Normalize each read's indels to their 3' place before voting.** Rejected.
  It does not handle a substitution split from an insertion.
- **Vote whole stretches only.** Rejected. No read observes a whole HVS-II
  stretch cleanly, and recall fell sharply.
- **Seed stretches from any difference.** Rejected. Masked reads widened
  stretches across neighbouring runs until no read matched their structure.
- **Emit variants directly from `dna consensus`.** Rejected. It would
  duplicate the call path and the output name of `dna call`.

## Consequences

**Measured downstream.** These figures are on the fair scale, with IUPAC
truth scored separately.

| Set | Rule | Exact samples | Precision | Recall |
|---|---|---|---|---|
| 320 evaluation | union (current) | 288 | 0.9924 | 0.9924 |
| 320 evaluation | consensus | 290 | 0.9967 | 0.9866 |
| disjoint 320 | union (current) | 290 | 0.9934 | 0.9926 |
| disjoint 320 | consensus | 292 | 0.9964 | 0.9901 |

- **Against reviewer consensus calls** (160 samples): 142 identical for the
  consensus, against 140 for the union.
- **Local Sequencher set:** 26 exact samples for the consensus, against 25.
- **Readability:** `N` covers about 0.9 % of decided positions.

**Remaining limits.**
- A run whose only reads meet it near their ends (`read_end`), or beside an
  unresolved call, stays `N`. Most remaining misses are in HVS-II.
- Point heteroplasmy is not reported (SRS-VAR-009).

**Unchanged contracts.** `dna sample`, `dna call` and their contracts are
unchanged, and no configuration is added.
