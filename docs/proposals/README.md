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

- [PROP-0001 — Modular DNA Analysis Platform](0001-modular-dna-analysis-platform.md) — **Implementing**; public Rust API, canonical contracts, multiple input modalities, replaceable scientific implementations and data providers, and downstream analysis modules.
- [PROP-0004 — Review triage](0004-review-triage.md) — **Draft**; catalogue errors and sort regions, traces, and samples into pass, review, and rerun with typed reasons, as DNA's [purpose](../requirements/SRS.md#purpose) states.
- [PROP-0005 — Evaluation method and fast feedback loop](0005-evaluation-feedback-loop.md) — **Draft**; labels, metrics, fixed test panels drawn from about 100,000 reviewed samples, and a cached evidence → rules → compare loop for validating triage rules. Harness v0 (`dna-eval`, ADR-0075) is implemented.
- [PROP-0007 — EMPOP export and labelled report conventions](0007-empop-export-and-report-conventions.md) — **Draft**; an offline EMPOP query (ranges and mitotype) in `notation/v1` whose ranges cover only decided positions, and profile conventions that label calls the reviewers do not record instead of deleting them.

## Implemented proposals

- [PROP-0002 — Plugin-first architecture](0002-plugin-first-architecture.md) — **Implemented**; modality plugins, a standalone core over per-read evidence, and post-calling plugins, delivered in phases with exit criteria.
- [PROP-0003 — Sample consensus](0003-sample-consensus.md) — **Implemented**; `dna consensus` adjudicates between a sample's reads into a consensus sequence that `dna call` calls (ADR-0073). Run-length refinement and IUPAC heteroplasmy are deferred.

Use [0000-template.md](0000-template.md) for new proposals. The [roadmap](roadmap.md) is a non-normative queue of future directions, not a substitute for a reviewed proposal.

Implemented/rejected proposals remain as historical evolution records when they retain useful rationale. Temporary planning notes that contain no durable knowledge should be deleted.
