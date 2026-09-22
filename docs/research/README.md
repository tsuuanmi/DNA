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

When a conclusion is fully promoted, delete duplicated temporary research text unless it still contains unique evidence worth retaining.

See [proposal lifecycle](../proposals/README.md), [documentation lifecycle](../governance/lifecycle.md), and [documentation governance](../governance/documentation.md).
