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
- [PROP-0003 — Sample consensus](0003-sample-consensus.md) — **Implementing**; a modality-neutral consensus layer above sample evidence that adjudicates between reads, emits a consensus sequence, and calls it through the core.

## Implemented proposals

- [PROP-0002 — Plugin-first architecture](0002-plugin-first-architecture.md) — **Implemented**; modality plugins, a standalone core over per-read evidence, and post-calling plugins, delivered in phases with exit criteria.

Use [0000-template.md](0000-template.md) for new proposals. The [roadmap](roadmap.md) is a non-normative queue of future directions, not a substitute for a reviewed proposal.

Implemented/rejected proposals remain as historical evolution records when they retain useful rationale. Temporary planning notes that contain no durable knowledge should be deleted.
