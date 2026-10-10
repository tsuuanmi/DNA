---
id: PROP-0003
type: proposal
status: implemented
owners: []
created: 2026-10-10
related-requirements: [SRS-SAMPLE-027, SRS-VAR-009, SRS-VAR-015]
related-decisions: [ADR-0023, ADR-0033, ADR-0035, ADR-0041, ADR-0065, ADR-0068, ADR-0069, ADR-0071, ADR-0072]
implementation: []
---

# Proposal: Sample consensus

## Problem

`dna sample` stops at evidence. It aggregates independently called reads
into coverage, overlaps, locus differences, support and opposition topology,
and a per-read notation. That notation is the union of every read's eligible
calls. DNA never adjudicates between reads. ADR-0023, ADR-0033, ADR-0035 and
ADR-0041 describe this evidence explicitly as *pre-consensus* and leave the
consensus for later.

Reviewers work the other way round. In Sequencher they build a consensus
sequence from the traces and call variants from it. Three measurements show
what the missing layer costs.

- **The core is already good enough.** On 160 held-out samples, `dna call` on
  the reviewers' consensus FASTA reproduces reviewer calls exactly for 158.
  The other two differ only in naming.
- **A read that one other read opposes is the largest tractable error.** On
  the fair scale (IUPAC truth scored separately), the 320-sample evaluation
  subset has 28 false and 28 missed calls. Outside the repeat windows, calls
  supported by exactly one read and opposed by another callable read were
  11 false against 3 true. The union cannot see this.
- **One bad read fails a whole sample.** On a second, disjoint 320-sample
  subset, 5 samples fail because a single read cannot be placed: low
  identity, or equally scoring placements. A consensus that tolerates a
  missing read would still cover the sample from the other strands.

## Goals

- Add a **consensus** layer that adjudicates between reads and produces one
  consensus sequence per sample, with evidence for every decision.
- **Keep it modality-neutral.** It consumes `CalledRead` records and sample
  evidence only, so Sanger, consensus input, and a future NGS modality can
  feed it (ADR-0069).
- **Call variants from the consensus through the existing core path.** Treat
  the consensus as a reviewed sequence (the `sequence` modality with vouched
  ends) and call it with `read_call`, as `dna call` does. There is then one
  variant-calling implementation.
- **Make every rule deterministic, versioned, and explained.** Each position
  records which reads agreed and which opposed, and why the rule decided as
  it did.
- **Improve on the union on the fair scale on both disjoint subsets.**
  Precision must be at least the union's, the count of exact samples must not
  fall, and no sample may fail because of a single read.

## Non-goals

- Point heteroplasmy as IUPAC codes. This is future work
  ([research](../research/point-heteroplasmy/README.md)); the consensus model
  only reserves a place for it.
- Dominant length of unreadable length heteroplasmy (M3, roadmap).
- Quantitative heteroplasmy, genotype, haplogroup inference, and clinical
  interpretation (SRS-VAR-009).
- Trace-level signal merging (Tracy-style assembly of chromatograms). The
  consensus works on per-read evidence, not on raw signal.
- Renaming or changing `sample`. Sample evidence stays the pre-consensus
  contract that reviewers and consumers can check.
- An NGS pileup implementation. Only the contract has to allow one.

## Evidence and research

Downstream comparisons only (ADR-0065); aggregates only. All figures use the
fair scale.

| Rule on notation calls | 320 evaluation subset, precision / recall / exact | 320 disjoint subset |
|---|---|---|
| union (current) | 0.9924 / 0.9924 / 288 | 0.9933 / 0.9930 / 286 |
| drop a single-read call opposed by a clean observation (relative quality above the veto threshold) | 0.9956 / 0.9910 / 290 | 0.9955 / 0.9925 / 289 |
| drop a single-read call opposed by any callable read | 0.9970 / 0.9888 / 293 | 0.9983 / 0.9919 / 294 |

- **Guarding repeat-length edits hurt.** Exempting edits that change a run
  length from the rule kept false deletions in two-base runs (`CC`, `TT`), so
  the rule needs no repeat exemption.
- **Unmasked primary calls are 99.98 % correct** (66 wrong of about 405,000).
  Most errors left are therefore disagreements *between* reads, not miscalled
  bases. This is consensus work.
- **The main false-call mechanisms** are single-read deletions, one misaligned
  short reverse-read island, and early-read insertions. In most of them
  another read covers the position and disagrees.
- See also [point-heteroplasmy research](../research/point-heteroplasmy/README.md)
  and the ADR-0071/ADR-0072 consequences.

## Proposed design

### Placement and boundaries

- **Location.** A new module `crates/dna-core/src/consensus/` with its own
  plugin descriptor (family `core`, id `consensus`). It owns a
  `[consensus]` configuration section only if a rule needs a parameter; the
  first increment needs none.
- **Inputs.** The sample's `CalledRead` records: evidence, alignment columns,
  and observed variants with eligibility. No Sanger type is used (ADR-0069).
- **Not here.** Profile windows and notation stay in post-calling. Consensus
  emits a sequence; nomenclature names the variants called from it.

### Per-position model

Over the union of the reads' mapped reference segments, each read contributes
one observation per reference position:

- **States:** `base`, `deletion`, `masked` (covered by a masked call),
  `untrusted` (within `read_end_margin` of an uninformative call), or
  `uncovered`.
- **Insertions:** each junction between two reference positions also gets one
  observation per read: the inserted string, or *none* when the read spans the
  junction with trusted calls.
- **Clean observations:** an observation is clean when its calls are
  unmasked, trusted, and carry no support veto. These are the same
  modality-neutral notions eligibility uses today.

### Decision rule (increment 1)

The rule is applied per position and per junction, over clean observations:

1. **All agree:** call that value.
2. **They disagree:** the value with more reads wins.
3. **A one-read tie** (one read against one read, the case in the evidence
   above): the reference-concordant value wins and the position is marked
   `contested`.
4. **Any other tie:** the position is unresolved (`N`, or no insertion at a
   junction) and marked `contested`.
5. **No clean observation:** the position is `unresolved` if some read
   covers it, and `uncovered` otherwise.

Masked observations are recorded as evidence but never decide.

Run-length decisions use the ADR-0071 notion: a read whose run end is not
resolved gives no clean observation for that run's length.

### Outputs and contracts

- **`dna.consensus/v1`** (new). It holds:
  - the sample identity, provenance (plugins, including `consensus` and its
    method version), and the reference identity;
  - the consensus segments in reference coordinates;
  - a sparse list of non-trivial positions and junctions, each with its
    decision (`called`, `contested`, `unresolved`) and the agreeing and
    opposing reads;
  - the reads rejected from the sample, with reasons.
- **A FASTA record per contiguous consensus segment.** It is reusable as
  `dna call` input.
- **Variants:** a `dna.variants/v1`-shaped result from running the consensus
  through `read_call`, plus notation from post-calling.
- **CLI.** `dna consensus <sample-id> <trace>... --reference <fasta>` writes
  `results/<sample-id>.consensus.json`, `.consensus.fasta`, and
  `.variants.json`. `dna sample` stays as it is.

### Read-failure tolerance (increment 0)

A read that cannot be placed should move to `rejected_reads` with a reason,
just as SRS-SAMPLE-027 already does for reads with too few callable calls. It
should not fail the operation. Both `sample` and `consensus` benefit.

## Alternatives considered

- **Rename `sample` to `consensus`.** Rejected. Sample evidence makes no
  decisions, and renaming it would erase the evidence/verdict boundary that
  ADR-0023/0033/0035/0041 keep.
- **A `contested` label inside sample evidence only.** Superseded. It would
  change the sample contract and then be re-derived by the consensus anyway.
- **Leave adjudication to downstream pipelines.** Rejected as the only answer.
  The core can already call a good consensus well, and producing that
  consensus is the gap between DNA and reviewer output. The evidence stays
  published for consumers who want their own rule.
- **Trace-level consensus**, merging chromatograms before calling. Rejected.
  It needs strand-specific signal registration and loses the per-read
  evidence model.
- **Variant-level adjudication without a sequence.** Rejected. A sequence
  gives the base-level benchmark against reviewer consensus FASTA, reuses the
  core caller, and is what reviewers produce.

## Validation plan

- **Base-level.** Compare the consensus with the 160 reviewer consensus FASTA
  records: identical records, base differences, `N` rate, and contested
  positions, with IUPAC sites reported separately.
- **Variant-level, fair scale, on both disjoint 320-sample subsets.** These
  must hold for the union and for `consensus`:
  - precision at least the union's;
  - the count of exact samples does not fall;
  - every contested decision that changed a call is listed in the downstream
    report.
- **Sequencher set.** No regression.
- **Robustness.** No sample fails because of one unplaceable read.
- **Tests.**
  - Unit tests for every decision branch, insertion junctions, untrusted and
    masked observations, and ties.
  - Integration tests for `dna consensus`.
  - A schema and example for `dna.consensus/v1`.
  - Byte-identity of `dna call` and `dna sample` outputs, apart from rejected
    reads in increment 0.

## Rollout and rollback

1. **Increment 0.** An unplaceable read becomes a rejected read in `sample`.
2. **Increment 1.** Consensus for bases and deletions under the rule above,
   with the `dna consensus` command and contract.
3. **Increment 2.** Insertions and run lengths at junctions.
4. **Increment 3** (future, separate decision). IUPAC point heteroplasmy
   under ADR-0009, once its prerequisites hold.

`dna consensus` is a new command, so rolling back an increment removes it
without touching `sample` or `call`.

## Observability and operations

Consensus records per-sample counts of called, contested, unresolved, and
uncovered positions as operational log events, alongside the existing stage
spans. Contested positions are listed in the result, so a reviewer can audit
every decision.

## Security and privacy

No new input formats or dependencies. Consensus FASTA identifies a sample
like the existing outputs do. It goes to `results/`, which is gitignored for
real data, under the same publication rules.

## Compatibility and migration

- New command and new contract.
- `dna.sample_evidence/v10` is unchanged apart from increment 0 (rejected-read
  reasons), revised in place under the unreleased-cycle rule.
- No migration is required.

## Decision

Accepted on 2026-10-10. Before the consensus contract lands (increment 1), an
ADR has to record three things:
- consensus as the adjudication layer above sample evidence;
- which earlier "pre-consensus" statements it fulfils;
- the decision rule and its versioning.

Increment 0 needs no such ADR. It extends SRS-SAMPLE-027.

## Implementation

- **Increment 0** (read-failure tolerance): an unplaceable read becomes a rejected
  read in `sample` (SRS-SAMPLE-029; `dna_core::read_call::PlacementRejection`).
- **Increment 1** (`dna consensus`, ADR-0073, SRS-CONS):
  - **Scope change:** increment 1 also covers insertions and run lengths,
    because the HVS-II calls depend on them. It decides them by stretches of
    whole reference runs rather than per junction.
  - **Variants:** they come from `dna call` over the consensus FASTA, not from
    the consensus command itself.
- **Increments 2 and 3: deferred** (closed 2026-10-10 by review).
  - **Run lengths:** after increment 1, the HVS-II misses that remain are
    mostly runs that no read reads. The union misses them too. Downstream
    probes (aggregates only, ADR-0065) found little more to gain:
    - dropping the `read_end` trust requirement for consensus bases gained 8
      calls on one 320-sample subset and none on the other;
    - restricting that change to runs with both ends visible gained 2.
  - **IUPAC point heteroplasmy** stays future work under ADR-0009
    ([research](../research/point-heteroplasmy/README.md)).
  - Both are listed in the [roadmap](roadmap.md).
