---
id: ADR-0079
status: accepted
date: 2026-10-18
related: [ADR-0069, ADR-0070, ADR-0073, ADR-0078, PROP-0008]
---

# ADR-0079: One crate per core step

<!-- docs:begin lifecycle -->
- **Status:** Accepted
- **Date:** 2026-10-18
- **Related:** [ADR-0069](0069-plugin-first-modality-core-post-calling.md), [ADR-0070](0070-workspace-split-by-plugin-family.md), [ADR-0073](0073-sample-consensus.md), [ADR-0078](0078-call-the-consensus-decisions.md), [PROP-0008](../../proposals/0008-pipeline-responsibilities.md)
<!-- docs:end lifecycle -->

<!-- docs:begin backlinks -->
- **Governs:** [architecture/overview](../../architecture/overview.md), [engineering/ci-cd](../../engineering/ci-cd.md), [reference/rust-api](../../reference/rust-api.md)
<!-- docs:end backlinks -->

## Context

ADR-0070 split the workspace by plugin family, and the core caller became
one crate, `dna-core`. It holds four steps, each answering one question:
- placement: where does a read lie on the reference?
- the read assertion: what does a placed read assert?
- sample evidence: what do a sample's reads show together?
- the consensus: which sequence do they agree on?

Inside one crate, nothing kept those steps apart. The consensus imported the
private trust rule and the whole core configuration, and the sample and the
consensus read alignment columns through crate-private fields. PROP-0008 gave
each step one owner (steps 1 to 4) and planned the split as step 5, a move of
code with outputs unchanged.

## Decision

1. The core is four crates in a chain, each depending only on the kernel and
   the steps before it:
   - `dna-align`: alignment, `Alignment`, and the `[alignment]` section;
   - `dna-call` (on `dna-align`): variant calling, `Eligibility`, the
     one-read path `call_read`, `CallConfig`, the `[variant_calling]`
     section, and the read records (`CalledRead`, variants, `CallTrust`);
   - `dna-sample` (on both): sample aggregation, its evidence, and the
     `[sample_reconciliation]` section;
   - `dna-consensus` (on `dna-align` and `dna-call`, not `dna-sample`): the
     consensus and its decided calls.

   `dna-core` is removed.
2. Only the types the next steps read cross a boundary, and each is public
   with the narrowest surface that serves them:
   - an `Alignment` gives a read-only view of its columns;
   - `Alignment::placed` is the one other way to build one, for a source that
     places a read itself (the consensus's decided reads). It checks the
     invariants `align_best` guarantees (calls and reference positions
     increase, positions lie in the covered segments, callable segments lie
     in them) and computes the metrics from the columns;
   - `Eligibility` is the one trust rule, with the queries the sample and the
     consensus use;
   - `CallConfig` holds the sections a called read needs (alignment and
     variant calling), so the consensus no longer takes the whole core
     configuration.

   Unchecked constructors that tests of downstream crates need sit behind
   each crate's `test-support` feature.
3. The kernel names the contract the consensus reads, `CalledReads`, which
   the `core` plugin provides. Contracts are validated at compile time and
   are not published.
4. Plugin identities do not change. `core` (version 11) stays one provenance
   identity across `dna-align`, `dna-call`, and `dna-sample`; its descriptor
   lives in `dna-call` next to `call_read` and still names the
   `sample_reconciliation` section, which `dna-sample` owns. Splitting that
   identity would change result provenance, so it is a later, versioned
   decision. `consensus` and `decided` live in `dna-consensus`.
5. The evaluation harness depends only on the crates whose records it reads
   (`dna-align`, `dna-sample`, `dna-consensus`).
   `validate_module_layers.py` enforces the whole map in sources and
   manifests.

## Consequences

- Every machine-readable output is unchanged. Operational log records of the
  core name the new crates in their target (`dna_call::read_call` instead of
  `dna_core::read_call`).
- The facade's configuration holds the call sections and the sample section
  side by side. The configuration file does not change.
- CI tests each core crate alone and checks that none depends on
  `dna-sanger` or `dna-eval`.
- Rust callers of the unpublished core crates change their paths; the `dna`
  facade API does not change.

## Alternatives considered

- **Keep `dna-core` with module-level rules.** Rejected: privacy inside one
  crate cannot keep a step from reading another step's internals, which is
  how the consensus came to import the private trust rule.
- **Move the decided builder into `dna-align`.** Rejected: it applies
  consensus rules (estimated run ends, gaps moved within runs). A validating
  constructor keeps `dna-align` the owner of what a valid alignment is.
- **Precompute trust as data on each called read.** Rejected: three of the
  consensus's queries judge a candidate site (its bases, orientation, and
  flanks), not a call, so they stay rules on `Eligibility`.
