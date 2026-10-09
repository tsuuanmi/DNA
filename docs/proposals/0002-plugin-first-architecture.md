---
id: PROP-0002
type: proposal
status: implementing
owners: []
created: 2026-10-09
related-requirements: [SRS-VAR-012, SRS-VAR-013]
related-decisions: [ADR-0058, ADR-0060, ADR-0064, ADR-0069]
implementation: []
---

# Proposal: Plugin-first architecture

## Problem

DNA should be usable as three independent parts:

1. modality plugins that turn raw data into clean per-read evidence;
2. a core that calls variants correctly from clean evidence and works without
   any modality or post-calling plugin;
3. post-calling plugins that check and analyze the called variants
   (normalization, nomenclature, EMPOP-style conformance, haplogroups).

Today the core is not usable alone:
- alignment and variant calling take Sanger types;
- one variant-calling filter evaluates chromatogram peaks;
- configuration is one monolithic document with every section required;
- post-calling capabilities exist only as library functions, with neither a
  command nor a JSON contract.

## Goals

- A modality-neutral core with a data contract at its input (`ReadEvidence`)
  and at its output (`CalledVariantSet`).
- Clean input, correct output: the core alone reproduces reviewer haplotypes
  from clean sequences.
- Plugins composed by a small kernel: contracts, a static registry, startup
  validation, per-plugin configuration sections, and provenance.
- Process-level extension through versioned JSON contracts.
- Every phase preserves scientific output unless it says otherwise.

## Non-goals

The non-goals of [PROP-0001](0001-modular-dna-analysis-platform.md) remain:
- no dynamic library loading;
- not every module becomes a trait or a crate;
- no claim of scientific equivalence between implementations.

In addition:
- no runtime workflow scheduler;
- no short-read NGS caller: per-read profile alignment and pairwise overlap do
  not scale to short-read depth, so an NGS modality would bring its own
  placement or join at `CalledVariantSet`;
- no plugin marketplace or third-party ABI.

## Evidence and research

- Plugin runtimes such as the DeepSeek agent harness describe "everything is a
  plugin". What they actually ship is a small kernel that loads plugins and
  resolves the services each one provides and requires. The kernel itself is
  not a plugin.
- Exploration of the current crate (2026-10-09) found that the alignment
  kernels, difference extraction, anchoring, sample aggregation, called-variant
  contracts, and representation capabilities are already modality-neutral. No
  target constants exist in production code; profiles hold them. Sanger types
  leak only through the alignment and variant-calling entry points, the
  support filter, the mask-to-reason mapping, the mixed exclusion vocabulary,
  sample read observations, and the monolithic configuration.

## Proposed design

The decision is recorded in
[ADR-0069](../decisions/adr/0069-plugin-first-modality-core-post-calling.md).
The design lands in phases. Each phase is its own change with its own exit
criterion.

| Phase | Scope | Exit criterion |
|---|---|---|
| 1 | `ReadEvidence` contract inside the crate. The Sanger adapter (`read_processing::evidence`) owns vetoes and mask reasons. Alignment and variant calling consume only `ReadEvidence`. The module validator enforces neutrality. | Every basecall, analysis, and sample result document is byte-identical before and after. |
| 2 | Static plugin registry: identity, family, method version, provided and required contracts, owned configuration sections. Per-plugin configuration sections. Plugin identities in provenance. | Outputs identical apart from provenance; registry validation rejects missing requirements and shared or unknown sections. |
| 2b | Sample attachments: core sample records hold neutral facts, and Sanger per-call and per-read extras move to an attachment that the report joins. | Byte-identical sample documents; `sample` imports no Sanger type. |
| 3 | A sequence adapter (consensus FASTA/FASTQ → `ReadEvidence`) and a core-only `dna call` command. | Reviewer consensus sequences give haplotype-identical calls after normalization and nomenclature. |
| 4 | `dna.called_variants/v1` JSON, emitted by the core and consumed by `dna normalize` and `dna nomenclature`. A conformance plugin (EMPOP-style notation rules). | JSON round trip gives the same result as in-process composition. |
| 5 | Workspace split (kernel, core, Sanger, post-calling, CLI crates), with a new ADR superseding ADR-0002. | Each crate builds and tests without the crates it does not depend on. |

Phase 1 is implemented as follows.
- `ReadEvidence` carries:
  - per call: the base, an optional evidence profile, an optional mask
    (unresolved or anchoring, with a reason label), and veto bits into an
    ordered per-read vocabulary;
  - per read: the informative interval.
- The core emits reasons in this order:
  1. region;
  2. support vetoes, in vocabulary order;
  3. `read_end`;
  4. mask reasons, in mapping order.

  It never interprets a modality label.

## Alternatives considered

See ADR-0069. In short:
- every module a plugin, or one universal trait: no coupling removed;
- dynamic or WASM plugins: no stable ABI, and needless sandbox cost;
- callback traits at the core boundary: they cannot cross a process and they
  hide ordering;
- one universal raw-alignment type: rejected by ADR-0060.

## Validation plan

- Phases 1, 2, and 2b:
  - Every result document is hashed over the local real-trace corpus and over
    a held-out set of 320 samples, before and after the change. Zero
    differences are required, except declared provenance fields.
  - Concordance is measured outside the repository (ADR-0065) and must be
    unchanged.
- Phase 3: reviewer consensus sequences, taken from the same held-out corpus,
  must reproduce reviewer haplotypes through the core alone.
- Phase 4: the round trip is tested against in-process composition.
- Every phase: all repository gates, including the module-layer validator
  with its neutrality rule.

## Rollout and rollback

Each phase is a separate change, revertible on its own. Phases 1 and 2b change
no output. Phase 2 changes the configuration schema, and phases 3 and 4 add
commands and contracts; each follows the versioning policy.

## Observability and operations

Phase 2 records plugin identities and versions in provenance. No other
operational change.

## Security and privacy

Process plugins exchange JSON documents through files or pipes. No code is
loaded at run time. Evaluation data stays outside the repository, and only
aggregates are quoted.

## Compatibility and migration

- Phase 1: none (internal). The Rust `Error` gains an `Evidence` variant.
- Phase 2: the unreleased configuration schema 7 is revised in place, and the
  shipped configuration is updated in the same change.
- Later contracts are new versioned documents.

## Decision

Accepted on 2026-10-09 and recorded by
[ADR-0069](../decisions/adr/0069-plugin-first-modality-core-post-calling.md).

## Implementation

- Phase 1: the `read_evidence` module, the Sanger adapter in
  `read_processing::evidence`, and the neutrality rule in
  `tools/python/scripts/validate_module_layers.py`.
- Phase 2: the compile-time registry and workflow compositions in
  `src/plugin.rs`, the Sanger-owned `[sanger_evidence]` configuration section
  (configuration schema 7, revised in place while unreleased), and
  `provenance.plugins` in the analysis, basecall, and sample documents. Apart
  from `configuration_sha256` and `provenance.plugins`, every result document
  of the local corpus and the 320-sample held-out set is unchanged.
- Phase 2b: `CalledRead` and `SangerAttachment` in `model`, neutral sample
  aggregation over `CalledRead`, report-side joins of Sanger evidence, and the
  neutrality rule for `sample` and the neutral children of `model`. Every
  result document and the sample aggregation metrics are byte-identical.
