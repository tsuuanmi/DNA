---
id: ADR-0087
status: accepted
date: 2026-10-21
family: [alignment, evaluation]
related: [ADR-0010, ADR-0029, ADR-0069, ADR-0075, PROP-0011]
---

# ADR-0087: Placed reads as a step contract that dna-align writes and reads

## Context

[[PROP-0011]] makes each step's output a contract the next step can start
from. Placement (`dna-align`) has only a partial one: `dna.analysis/v10`
publishes the alignment's summary and each call's column, not a placement
that a later step can read back. Every analysis therefore searches every
read again.

On development panel v5 with a warm read evidence cache, a probe that
changes only the consensus or the review still spends most of its CPU
re-placing reads: 3,709 CPU-s per run, about 66,000 placements. The
evaluation loop is the main user, but the contract is the same one any step
needs to start from placed reads.

`dna-align` is the only producer of alignments (`align_best`, and
`Alignment::placed` for a consensus record), so the alignment's invariants
hold wherever it is read. A placement read back from a file must not weaken
that.

## Decision

Proposed:

- **Contract.** `dna.placed_read/v1` holds one read's placement as
  `align_best` selected it: orientation, score, origin wrap, reference and
  callable segments, metrics, and every column, with the SHA-256 of the read
  evidence document it places, the reference, and the configuration.
- **Owner.** `dna-align` alone writes (`placed::write`) and reads it back
  (`placed::read`), so it stays the only producer of alignments.
- **Read-back checks.** Reading re-checks the document against the read's
  evidence and the reference:
  - each column's base against the call (`N` for an unresolved call) and
    the reference;
  - the columns' order for the orientation, and the origin wrap;
  - every column inside the reference segments;
  - the callable segments and masked bases, recomputed from the evidence's
    masks;
  - the metrics' consistency.

  The score is recorded as selected. Whoever stores the document keys it by
  the placement code, so a change to the search is a new key, not a stale
  read.
- **Core.** `dna_call::read_call::call_read_placed` starts from a given
  placement. Everything after placement is unchanged.
- **Facade.** `dna::recorded::placed_read_document` writes the document for
  one read evidence document. `analyze_placed_outcome` analyses a sample
  from read evidence documents with their placed-read documents.
- **Evaluation.** `dna-eval`'s `--evidence-cache` directory also holds a
  placement cache. Entries are keyed by the evidence document's SHA-256,
  under a generation named for the configuration, the placement stage's
  code identity (git object ids of the kernel, `dna-align`, and the facade's
  configuration, evidence reading, and writer), and the reference. The cache
  is gzip-compressed, writes are atomic, a corrupt entry is a miss, and the
  run's provenance lists the cache in `run.json`.
- **Status.** Accepted. It touches production crates (`dna-align`,
  `dna-call`, the facade) for a speed-up the evaluation loop uses first. No
  CLI flag reads placed-read files yet.

## Alternatives

- **A restore hook in the core, kept by the harness.** Rejected. `dna-eval`
  would serialize alignments itself and inject them through a production
  hook only it uses. `dna-align` would no longer be the only writer, and no
  contract would describe the bytes.
- **Faster search.** Still open, as separate work. Most reads after a
  homopolymer tract fall back from the proof tier to full Gotoh; making that
  faster must keep the score optimum and the ambiguity outcomes
  (INV-ALN-005). It would help production too.
- **A stratified quick subset.** Kept as a separate pre-filter
  (`dna-eval subset`). It cuts samples, not per-read work, and is never the
  gate.

## Consequences

Measured on development panel v5, from the same commit (0ff4e1c):

- **Byte identity.** All nine tables (`scores`, `checks`, `unreviewed`,
  `estimates`, `components`, `tracts`, `anchors`, `forms`, `windows`) are
  byte-identical across three runs: main (85cab63) without placements, a
  cold placement cache (66,211 placements made and stored, about 200 MB),
  and a warm one (66,211 read back).
- **Speed.** On a machine loaded to about 180-235 of 256 cores:
  - main with a warm evidence cache: 3,709 CPU-s, 42.9 s;
  - warm evidence and placement caches: 2,623 CPU-s, 46.3 s;
  - a consensus-only probe, both caches warm: 2,356 CPU-s, 28.7 s.

  CPU falls 29-36 %. Wall time follows the machine's load more than the
  cache. The gain is smaller than the single-sample trace suggested: the
  remaining time is in parsing documents, the consensus, the review, and the
  harness's comparison.
- **Keys.** A change to `dna-align` or the placement stage misses, as
  intended. The first run after a `Cargo.lock` change also misses the read
  evidence cache, since both stages key on it.
