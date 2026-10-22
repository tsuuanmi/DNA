# Workspace crates

DNA is a Cargo workspace of plugin-family crates, with the core split by
step ([ADR-0070](../docs/decisions/adr/0070-workspace-split-by-plugin-family.md),
[ADR-0069](../docs/decisions/adr/0069-plugin-first-modality-core-post-calling.md),
[ADR-0079](../docs/decisions/adr/0079-core-crates-by-step.md)). The `dna` facade at the repository root ([src](../src/README.md)) composes
them. It owns the CLI, the pipelines, the configuration envelope, the plugin
registry, and the reports.

| Crate | Owns | Depends on |
|---|---|---|
| [dna-kernel](dna-kernel/src/README.md) | shared contracts: errors, read evidence, called variants, target profiles, reference loading, plugin descriptors | — |
| [dna-align](dna-align/src/README.md) | placement: evidence-profile alignment of one read | `dna-kernel` |
| [dna-call](dna-call/src/README.md) | the read assertion: one placed read's variants, trust, and callable positions; the `core` plugin descriptor | `dna-kernel`, `dna-align` |
| [dna-sample](dna-sample/src/README.md) | sample evidence: what a sample's reads show together, for review | `dna-kernel`, `dna-align`, `dna-call` |
| [dna-consensus](dna-consensus/src/README.md) | the consensus and its decided calls (ADR-0073, ADR-0078) | `dna-kernel`, `dna-align`, `dna-call` |
| [dna-sanger](dna-sanger/src/README.md) | the Sanger modality: ABIF decoding to `ReadEvidence` | `dna-kernel`, `dna-trace` |
| [dna-trace](dna-trace/src/README.md) | trace signal engineering from plain data: phase-mixture decomposition of the windows a modality flags ([ADR-0084](../docs/decisions/adr/0084-trace-crate.md)); a library, not a plugin | `dna-kernel` |
| [dna-normalize](dna-normalize/src/README.md) | post-calling representation, normalization, and nomenclature in the lab's reviewer convention ([ADR-0077](../docs/decisions/adr/0077-normalize-and-empop-crates.md)) | `dna-kernel` |
| [dna-empop](dna-empop/src/README.md) | post-calling conformance to forensic (EMPOP-style) notation conventions, offline | `dna-kernel` |
| [dna-eval-lib](dna-eval-lib/src/README.md) | the evaluation harness: panel manifests, reviewer truth, comparison, reports, the artefact and simulated benches, and the evidence, placement, and forms caches ([ADR-0075](../docs/decisions/adr/0075-evaluation-harness-in-workspace.md)) | `dna`, `dna-kernel`, `dna-align`, `dna-sample`, `dna-consensus`, `dna-sanger`, `dna-trace` (its `simulate` feature; tests: `dna-normalize`) |
| [dna-eval](dna-eval/src/README.md) | the `dna-eval` binary: what it was built from (revision, and each cached stage's code identity from `stages.txt`), over the harness | `dna-eval-lib` |

The core is four crates, one per step, in a chain: each depends only on the
kernel and the steps before it ([ADR-0079](../docs/decisions/adr/0079-core-crates-by-step.md)).
The other plugin crates never depend on one another, and only the facade
composes them. `dna-trace` is a library, not a plugin: of the plugin crates only
`dna-sanger` depends on it ([ADR-0084](../docs/decisions/adr/0084-trace-crate.md)), and the evaluation harness
renders the simulated bench's reads with its `simulate` feature. The evaluation harness may depend on the facade and the plugin crates,
and no crate depends on it but its thin binary, which keeps the build
stamps that change with every commit out of the harness's compilation; it
is the only place for truth formats and comparison (SRS-COMPAT-003). It is tested with the workspace but not built for
or shipped in a release. `tools/python/scripts/validate_module_layers.py`
enforces this graph in source paths and manifests, and CI tests each plugin
crate on its own.
