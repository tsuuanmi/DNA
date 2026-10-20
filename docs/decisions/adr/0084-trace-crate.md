---
id: ADR-0084
status: proposed
date: 2026-10-20
family: [callability]
related: [ADR-0069, ADR-0070, ADR-0079, ADR-0082, ADR-0083, PROP-0012]
---

# ADR-0084: Trace signal engineering is a library crate the modality calls

<!-- docs:begin lifecycle -->
- **Status:** Proposed
- **Date:** 2026-10-20
- **Related:** [ADR-0069](0069-plugin-first-modality-core-post-calling.md), [ADR-0070](0070-workspace-split-by-plugin-family.md), [ADR-0079](0079-core-crates-by-step.md), [ADR-0082](0082-phase-mixtures.md), [ADR-0083](0083-set-aside-dissent.md), [PROP-0012](../../proposals/0012-trace-crate.md)
<!-- docs:end lifecycle -->

<!-- docs:begin backlinks -->
- **Governs:** [design/length-mixture](../../design/length-mixture.md), [engineering/ci-cd](../../engineering/ci-cd.md), [reference/configuration](../../reference/configuration.md), [reference/read-evidence](../../reference/read-evidence.md), [reference/review](../../reference/review.md), [requirements/callability](../../requirements/callability.md)
<!-- docs:end backlinks -->

## Context

Length heteroplasmy after a homopolymer is the largest block of review work
left. Resolving it is numerical work on the trace: fitting shifted copies of
the read's own signal, estimating per-read base signatures and leakage, and
choosing the number of components. [[ADR-0082]] put the first version inside
`dna-sanger`, as the modality's assertion.

The user asked for this work in a crate of its own (2026-10-20), and chose
the name `dna-trace`. "Signal" stays the run-level product (`dna signal`,
`dna.signal/v1`).

The workspace rule ([[ADR-0070]], [[ADR-0079]]) is that plugin crates never
depend on one another. A crate the Sanger plugin depends on needs a stated
place in that graph.

## Decision

`dna-trace` is a **library crate, not a plugin**.

- It has no plugin descriptor, no configuration section of its own, and no
  output document.
- It depends on `dna-kernel` only, and `dna-sanger` depends on it. The
  plugin-to-plugin rule is unchanged: `dna-trace` is not a plugin family, and
  no other plugin depends on it today.
- Its input is plain data: the analysed channels, the call positions, the
  primary calls, and the call ranges the modality flags. It reads no ABIF,
  reference, alignment, or truth.
- Its output is typed results. The modality turns them into kernel
  assertions (`PhaseMixture`, [[ADR-0082]]), and the core reads only those.
- Its parameters belong to the modality's configuration, validated and
  recorded in provenance through the `sanger` plugin. Levels come from the
  read's or the run's own data wherever possible.
- Stage 1 re-processes only the windows the modality flags after a repeat
  run. Cleaning the whole read before basecalling is a separate decision.
- A run's length is decided from decoded components only when the decode is
  clear and both strands agree. Otherwise the review flags it with the
  components ([[PROP-0012]]).

### Two length assertions, for now

Increment 2 adds a second assertion over related windows rather than
extending `PhaseMixture`:

- `PhaseMixture` ([[ADR-0082]]) is the offset model. It fits each dephased
  segment as the read's own calls shifted by whole calls, and its explained
  gate feeds the consensus's set-aside rules ([[ADR-0083]]).
- `TractLengths` is the length decode. Before every call where callability
  sees the read leave phase, not only at dephased segments, it locates the
  tract from the read's primaries, chooses its structure (a single run, or
  two runs around a mark), and reports the first run's and the total's
  length marginals.

The two locate their run differently. Matching them by position would be a
rule of its own, so neither is fitted into the other. The core reads nothing
from `tracts` in increment 2. Increment 3 or 4 decides, with data, whether
one replaces the other, so the duplication stays temporary and visible. The
decode's constants are the `[trace]` section, which the `sanger` plugin
declares and validates.

## Alternatives

- **Keep the decoder in `dna-sanger`.** Simpler, but the user asked for a
  crate. A crate also tests on synthetic traces and serves another
  dye-terminator modality.
- **Make `dna-trace` a plugin family with its own contract.** That adds a
  plugin descriptor and an output document with no consumer besides the
  modality. The kernel assertion already is the contract.
- **Let `dna-trace` read ABIF.** That would duplicate `dna-sanger`'s decoder
  and tie the crate to one instrument format.

## Consequences

- `crates/README.md`, the module-layer validator, and CI's per-crate tests
  gain `dna-trace` as a dependency of `dna-sanger` alone.
- `dna-sanger/src/decompose` moves into `dna-trace`
  ([[PROP-0012]] increment 1, byte-identical).
- The `sanger` plugin version covers every output change that `dna-trace`
  causes.
- Until increment 3 or 4 settles it, a read can carry both a `PhaseMixture`
  and a `TractLengths` for one run.
