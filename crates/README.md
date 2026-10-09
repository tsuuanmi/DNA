# Workspace crates

DNA is a Cargo workspace of plugin-family crates
([ADR-0070](../docs/decisions/adr/0070-workspace-split-by-plugin-family.md),
[ADR-0069](../docs/decisions/adr/0069-plugin-first-modality-core-post-calling.md)).
The `dna` facade at the repository root ([src](../src/README.md)) composes
them. It owns the CLI, the pipelines, the configuration envelope, the plugin
registry, and the reports.

| Crate | Owns | Depends on |
|---|---|---|
| [dna-kernel](dna-kernel/src/README.md) | shared contracts: errors, read evidence, called variants, target profiles, reference loading, plugin descriptors | — |
| [dna-core](dna-core/src/README.md) | the core caller: alignment, variant calling, sample aggregation | `dna-kernel` |
| [dna-sanger](dna-sanger/src/README.md) | the Sanger modality: ABIF decoding to `ReadEvidence` | `dna-kernel` |
| [dna-post](dna-post/src/README.md) | post-calling plugins: representation, normalization, nomenclature, conformance | `dna-kernel` |

The plugin crates never depend on one another, and only the facade composes
them. `tools/python/scripts/validate_module_layers.py` enforces this in source
paths and manifests, and CI tests each crate on its own.
