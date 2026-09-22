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

Use [0000-template.md](0000-template.md) for new proposals. The [roadmap](roadmap.md) is a non-normative queue of future directions, not a substitute for a reviewed proposal.

Implemented/rejected proposals remain as historical evolution records when they retain useful rationale. Temporary planning notes that contain no durable knowledge should be deleted.
