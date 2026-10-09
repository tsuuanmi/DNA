# ADR-0070: Workspace split by plugin family

- **Status:** Accepted
- **Date:** 2026-10-09
- **Related proposal:** [PROP-0002](../../proposals/0002-plugin-first-architecture.md) (phase 5)
- **Related decisions:** [ADR-0002](0002-single-crate-layering.md),
  [ADR-0058](0058-canonical-contracts-and-modular-analysis-composition.md),
  [ADR-0064](0064-crate-ready-module-layering.md),
  [ADR-0069](0069-plugin-first-modality-core-post-calling.md)

## Context

[ADR-0002](0002-single-crate-layering.md) keeps one crate until reuse, build
isolation, or independent release justifies a split.
[ADR-0064](0064-crate-ready-module-layering.md) names the trigger: "a second
sequencing modality needs the core and representation layers without the
Sanger stack". That trigger now holds. The reviewed-consensus `sequence`
modality, `dna call`, and `dna notation` run the core and the post-calling
plugins without Sanger. On 160 held-out samples, the core alone reproduces
reviewer calls from consensus sequences (PROP-0002 phase 3).

ADR-0069 already composes three plugin families through data contracts.
PROP-0002 phase 5a untangled the sources along those families in place:
- each plugin owns its configuration sections and its descriptor;
- the facade composes the configuration envelope and the registry;
- a validator enforced the target crate map.

## Decision

1. DNA is a Cargo workspace. The `dna` facade crate stays at the repository
   root, and the plugin crates live in `crates/`:
   - `dna-kernel` holds the shared contracts:
     - the one `Error` of every stage;
     - `ReadEvidence`;
     - the public called-variant contracts;
     - target profiles and reference loading;
     - the plugin descriptor types with their compile-time validation.
   - `dna-core` is the core caller: alignment, variant calling, and sample
     aggregation.
   - `dna-sanger` is the Sanger modality, from ABIF decoding to
     `ReadEvidence`.
   - `dna-post` holds the post-calling plugins: representation,
     normalization, nomenclature, and conformance.
   - `dna` holds the CLI, the pipelines, the configuration envelope, input
     orchestration, the consensus-sequence and variants-document adapters,
     the plugin registry and workflow compositions, `analyze_sanger`, and the
     reports.
2. The plugin crates depend only on `dna-kernel`, and only `dna` composes them.
   The module validator enforces this graph in source paths and manifests.
   CI also tests each plugin crate on its own.
3. One `Error` stays in the kernel, so that `dna::error` is unchanged. The
   facade re-exports the existing public paths (`dna::error`, `dna::variant`,
   `dna::profile`, `dna::variant_normalization`, `dna::variant_nomenclature`)
   with the same items.
4. Items used across crates are `pub` and documented. Everything else keeps
   crate visibility. All crates are `publish = false` and share one version,
   which the release tag matches through the `dna` package.
5. Test fixtures that cross crates sit behind the kernel's `test-support`
   feature. The fuzz harness stays outside the workspace and reaches the ABIF
   parser through `dna::fuzzing` and the Sanger crate's `fuzzing` feature.

## Alternatives

- **Facade moved to `crates/dna` under a virtual root:** more symmetric, but it
  moves the integration tests, fixtures, and many documentation links for no
  consumer benefit.
- **One error type per crate:** fully separate, but it changes the public
  `dna::error` API and about 300 signatures. The kernel's error module has no
  dependencies, so keeping it whole costs only shared vocabulary.
- **Make every item public:** fewer edits, but it publishes about 400
  undocumented internals. Only the items used across crates are public.

## Consequences

- Each crate builds and tests without the crates it does not depend on. The
  core and the post-calling plugins can be reused without the Sanger stack.
- Cargo commands use `--workspace`, and binaries are built with `-p dna`.
- Result documents and the public Rust API are unchanged. Operation-log
  payloads are unchanged; the module path of events moved to the plugin crates
  now names that crate (for example `dna_sanger::read_processing`).

## Supersession

- [ADR-0002](0002-single-crate-layering.md) is superseded.
- [ADR-0064](0064-crate-ready-module-layering.md) is superseded in part:
  - its layer table and neutral-module rule give way to the crate graph
    above;
  - its acyclic-module rule still holds within each crate.
- [ADR-0058](0058-canonical-contracts-and-modular-analysis-composition.md) §9
  (single-crate layering) is superseded.
