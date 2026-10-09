# ADR-0064: Crate-ready module layering

- **Status:** Accepted
- **Date:** 2026-10-08
- **Related decisions:** [ADR-0002](0002-single-crate-layering.md),
  [ADR-0058](0058-canonical-contracts-and-modular-analysis-composition.md),
  [ADR-0063](0063-target-profiles.md)

## Context

ADR-0002 keeps one crate until a split is justified by reuse, build isolation,
or independent release, and ADR-0058 adds the precondition of a cycle-free
dependency graph. None of the justifications holds yet: DNA is about 18k lines,
one unpublished binary, with no external library consumer and one sequencing
modality.

The precondition did not hold either. The module graph had a cycle through
`profile`, `variant_normalization`, `variant_analysis`, `input`, and
`variant_representation`. The canonical called-variant types (`Variant`,
`VariantKind`, `ReferenceIdentity`, `CalledVariantSet`) lived in the
`variant_analysis` capability, so the generic representation modules depended
upward on a Sanger capability, against ADR-0058's rule that contracts must not
depend on implementations.

## Decision

1. DNA stays one crate (ADR-0002 is not superseded).
2. Every top-level module belongs to one layer. Production code may depend only
   on the same or a lower layer, and the module graph must be acyclic:

   | Layer | Modules | Possible future crate |
   |---|---|---|
   | 0 core | `error`, `checksum`, `model`, `locus`, `variant` | `dna-core` |
   | 1 target data | `config`, `reference`, `profile` | `dna-core` or `dna-profile` |
   | 2 science | `basecalling`, `signal_processing`, `callability`, `quality_control`, `read_processing`, `alignment`, `variant_calling`, `sample`, `variant_representation`, `variant_normalization`, `variant_nomenclature` | `dna-sanger`, `dna-analysis`, `dna-representation` |
   | 3 adapters and capabilities | `input`, `variant_analysis` | `dna-input-sanger`, `dna-variant-analysis` |
   | 4 delivery | `report`, `operation_log`, `pipeline`, `cli` | `dna` (facade and CLI) |

3. The canonical called-variant contracts live in the public core module
   `dna::variant`; capabilities produce or consume them.
4. Profiles name their own vocabulary (`IndelPlacement`) instead of importing a
   capability's policy type; `variant_normalization` maps it to
   `NormalizationPolicy`.
5. `tools/python/scripts/validate_module_layers.py` enforces the layer map and
   acyclicity in CI. A new top-level module must be assigned a layer there and
   in this table.

A split into crates follows when one of these holds:

- a second sequencing modality (for example an NGS adapter) needs the core and
  representation layers without the Sanger stack;
- an external consumer uses the library;
- build times make crate-level compilation isolation worthwhile;
- a component needs an independent release cadence.

## Alternatives

- **Split into a workspace now.** It adds versioning, CI, and release overhead
  and changes public paths with no consumer to benefit.
- **Only break the cycle.** The cycle would return unnoticed; the CI check is
  what keeps the graph crate-ready.

## Consequences

- A future split moves whole layers into crates without redesigning contracts.
- **Breaking (Rust API):** `dna::variant_analysis::{Variant, VariantKind,
  ReferenceIdentity, CalledVariantSet}` move to `dna::variant`.
- Scientific output is unchanged.
