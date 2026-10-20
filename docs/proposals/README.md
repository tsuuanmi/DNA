# Proposals

Proposals describe **changes under consideration**. They are evolutionary artifacts, not current production truth.

## Lifecycle

```text
Draft -> Proposed -> Accepted -> Implementing -> Implemented
                    \-> Rejected
Accepted/Implemented -> Superseded
```

- Research provides evidence.
- A proposal frames a concrete change, alternatives, design, validation, rollout, and risks.
- Acceptance does not itself change production behavior.
- Architecturally significant accepted choices may require an ADR.
- Current truth is updated only when requirements/design/reference/source/tests/validation are changed.

## Active proposals

<!-- docs:begin index:proposal:open -->
- [PROP-0001 — Modular DNA Analysis Platform](0001-modular-dna-analysis-platform.md) — **Implementing**; public Rust API, canonical contracts, multiple input modalities, replaceable scientific implementations and data providers, and downstream analysis modules.
- [PROP-0004 — Review triage](0004-review-triage.md) — **Accepted**; DNA's canonical per-sample review list (`dna.review/v1`): the positions a reviewer should look at, with evidence-based labels and no conclusion, for the laboratory's review tool. The review list is implemented (`dna review`, version 1); watched positions are pending.
- [PROP-0005 — Evaluation method and fast feedback loop](0005-evaluation-feedback-loop.md) — **Draft**; labels, metrics, fixed test panels drawn from about 100,000 reviewed samples, and a cached evidence → rules → compare loop for validating triage rules. Harness v0 (`dna-eval`, ADR-0075) is implemented.
- [PROP-0008 — One responsibility per pipeline step](0008-pipeline-responsibilities.md) — **Accepted**; the modality decides trust as per-call verdicts, the core asserts only what a placed read shows whole, the consensus votes stretch sequences on bounded reads, and the core splits into align, call, sample, and consensus crates.
- [PROP-0009 — Run signal for second-channel judgement](0009-run-signal.md) — **Draft**; the Sanger mixed verdict and `possible_heteroplasmy` judge a second channel against the level its own run shows, keyed by the read's own 5-mer context, given as a recorded `dna.signal/v1` input from the Sanger signal processing; fixed thresholds stay the fallback.
- [PROP-0010 — Decompose a length mixture into its sequences](0010-length-mixture-decomposition.md) — **Implementing**; after a repeat run with length heteroplasmy, the Sanger plugin decomposes the dephased stretch into the read's own sequence at offsets with fractions, each component is compared with the reference, and the consensus decides a run's dominant length where both strands agree.
- [PROP-0011 — A contract at every step](0011-step-contracts.md) — **Implementing**; every crate's output is a versioned, schema-validated contract producible on its own (`dna.read_evidence/v1`, `dna.analysis/v10` with per-call trust and component calls, schemas for the evaluation outputs), with a crate → contract → command map, so a wrong result is traceable to the step that made it.
- [PROP-0012 — Trace signal engineering in its own crate](0012-trace-crate.md) — **Draft**; a library crate `dna-trace` re-processes the trace windows the Sanger plugin flags after a repeat run: it models the read's own base signatures and leakage, decodes the length mixture, and gives each run's components with fractions and a margin; the consensus decides a run's length only when the decode is clear and both strands agree, and otherwise flags it with the components, so no length is passed unflagged on a weak signal.
<!-- docs:end index:proposal:open -->

## Closed proposals

<!-- docs:begin index:proposal:done -->
- [PROP-0002 — Plugin-first architecture](0002-plugin-first-architecture.md) — **Implemented**; modality plugins, a standalone core over per-read evidence, and post-calling plugins, delivered in phases with exit criteria.
- [PROP-0003 — Sample consensus](0003-sample-consensus.md) — **Implemented**; `dna consensus` adjudicates between a sample's reads into a consensus sequence that `dna call` calls (ADR-0073). Run-length refinement and IUPAC heteroplasmy are deferred.
- [PROP-0006 — Call the consensus's decisions without re-aligning them](0006-call-consensus-decisions.md) — **Implemented**; `call` reports the edits the consensus decided instead of aligning its FASTA again, which loses where `N`s sit (placement, run bounds, AC deletions), and every read places an equivalent indel at one canonical position before the consensus votes, with the laboratory's conventions applied after it.
- [PROP-0007 — EMPOP export and labelled report conventions](0007-empop-export-and-report-conventions.md) — **Implemented**; an offline EMPOP query (ranges and mitotype) in `notation/v1` whose ranges cover only positions DNA decided, and profile conventions that label calls the reviewers do not record instead of deleting them.
<!-- docs:end index:proposal:done -->

Use [0000-template.md](0000-template.md) for new proposals. The [roadmap](roadmap.md) is a non-normative queue of future directions, not a substitute for a reviewed proposal.

An implemented, rejected, or withdrawn proposal stays live while it retains rationale no current document holds; otherwise it is retired to the [ledger](../governance/retired.md) with a one-sentence outcome. Lists above are generated from each proposal's front matter (`status`, `summary`).
