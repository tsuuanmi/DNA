# Research

Research is **exploratory evidence**, not current production truth.

Use `docs/research/<topic>/` for external-system analysis, experiments, algorithm investigation, benchmarks under exploration, and evidence that has not yet been accepted into the product.

## Promotion

```text
research evidence
      ↓
proposal / review
      ↓
decision when architecturally significant
      ↓
requirements + architecture/design/reference
      ↓
source + tests
      ↓
validation + release
```

A research note does not change production behavior merely by existing. Agents must not implement directly from research when it conflicts with or has not been promoted into canonical current-state documentation.

When a conclusion is fully promoted, retire the note to the [ledger](../governance/retired.md) with a one-sentence outcome unless it still contains unique evidence worth keeping live.

See [proposal lifecycle](../proposals/README.md) and [documentation governance](../governance/documentation.md).

## Topics

Open:

<!-- docs:begin index:research:open -->
- [unused ABIF tags](abif-metadata/README.md) — **Active**; Which ABIF tags beyond the five DNA reads (DATA.9-12, FWO_.1, PLOC.2, PBAS.2, PCON.2) could add evidence, and how many carry an instrument-side failure or run context.
- [signal denoising and baseline correction](denoising/README.md) — **Active**; Does a conditioned projection of the instrument-analyzed channels improve re-calling or read callability on real traces?
- [C-tract length heteroplasmy](length-heteroplasmy/README.md) — **Active**; How forensic practice records C-tract length heteroplasmy, how precise DNA's tract-length review items are per evidence state, and whether decode noise should be calibrated per laboratory batch.
- [when may a region read by one strand pass?](one-strand/README.md) — **Active**; Evidence for the open question on whether a single strand can pass a region.
- [how other systems decide pass, review, and rerun](review-triage/README.md) — **Active**; Evidence supporting the review triage proposal.
<!-- docs:end index:research:open -->

Concluded:

<!-- docs:begin index:research:done -->
- [phase-shift classification and recovery behind homopolymer runs](phase-recovery/README.md) — **Concluded**; Recovered calls as variant evidence after phase loss; concluded, not promoted.
- [point-heteroplasmy evidence in Sanger reads](point-heteroplasmy/README.md) — **Concluded**; Point heteroplasmy as future work; concluded, not promoted.
<!-- docs:end index:research:done -->
