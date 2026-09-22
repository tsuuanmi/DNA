# Documentation Lifecycle

DNA classifies documentation by authority and durability.

## Canonical / living

Current truth that must change with the system:

- requirements;
- architecture and invariants;
- design;
- reference/contracts;
- validation policy;
- engineering/operations/security/governance policy.

Typical lifecycle: `draft -> active -> superseded/replaced`.

## Historical / durable

Evidence of why or what happened:

- ADRs;
- postmortems when incidents exist;
- implemented/rejected proposals that retain useful rationale;
- durable benchmark/validation reports when retained.

Historical artifacts are not rewritten to look current.

## Evolutionary

Proposals describe changes under review:

`draft -> proposed -> accepted/rejected -> implementing -> implemented -> superseded`.

Acceptance alone does not make a proposal current production truth.

## Exploratory

Research provides evidence and candidate ideas. It is non-normative.

Typical lifecycle: `active -> concluded -> promoted or deleted`.

## Promotion

```text
research -> proposal -> decision (when needed)
        -> requirements/architecture/design/reference
        -> source + tests -> validation -> release
```

Promoted content moves to its canonical home. Delete duplicated temporary research/planning once it contains no unique evidence.

## Metadata

Structured lifecycle artifacts should use lightweight front matter where it helps tooling:

```yaml
---
id: PROP-0017
type: proposal
status: implementing
owners: [signal-core]
created: 2026-09-20
related-requirements: [SRS-042]
related-decisions: [ADR-0054]
implementation: [PR-48]
---
```

Do not add metadata mechanically to documents where it provides no lifecycle or ownership value.
