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
| [dna-normalize](dna-normalize/src/README.md) | post-calling representation, normalization, and nomenclature in the lab's reviewer convention ([ADR-0077](../docs/decisions/adr/0077-normalize-and-empop-crates.md)) | `dna-kernel` |
| [dna-empop](dna-empop/src/README.md) | post-calling conformance to forensic (EMPOP-style) notation conventions, offline | `dna-kernel` |
| [dna-eval](dna-eval/src/README.md) | the evaluation harness: panel manifests, reviewer truth, comparison, and reports ([ADR-0075](../docs/decisions/adr/0075-evaluation-harness-in-workspace.md)) | `dna` |

The plugin crates never depend on one another, and only the facade composes
them. The evaluation harness may depend on the facade and the plugin crates,
and no crate depends on it; it is the only place for truth formats and
comparison (SRS-COMPAT-003). It is tested with the workspace but not built for
or shipped in a release. `tools/python/scripts/validate_module_layers.py`
enforces this graph in source paths and manifests, and CI tests each plugin
crate on its own.
