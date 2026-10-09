# ADR-0069: Plugin-first composition of modality, core, and post-calling plugins

- **Status:** Accepted
- **Date:** 2026-10-09
- **Related proposal:** [PROP-0002](../../proposals/0002-plugin-first-architecture.md)
- **Related decisions:** [ADR-0002](0002-single-crate-layering.md),
  [ADR-0058](0058-canonical-contracts-and-modular-analysis-composition.md),
  [ADR-0059](0059-reuse-ecosystem-machinery-behind-dna-contracts.md),
  [ADR-0060](0060-separate-variant-canonicalization-nomenclature.md),
  [ADR-0064](0064-crate-ready-module-layering.md),
  [ADR-0067](0067-signal-derived-read-callability.md)

## Context

DNA is meant to grow beyond one Sanger path: other sequencing modalities,
other callers, and analyses on called variants (normalization, nomenclature,
EMPOP-style conformance, haplogroups). The intended shape is plugin-first, in
the spirit of "everything is a plugin" agent runtimes: a small kernel that only
knows contracts and composition, with every capability supplied by a plugin.

Three independent parts follow from that:

1. **Modality plugins** turn raw data (Sanger traces today) into clean per-read
   evidence.
2. **The core** turns clean evidence into correct variant calls and must be
   usable on its own: good evidence in, correct calls out, with no dependency
   on any modality or post-calling plugin.
3. **Post-calling plugins** consume called variants.

Part 3 already works without the others through `CalledVariantSet` (ADR-0060).
Part 2 does not:
- alignment takes Sanger quality-control, callability, and signal types;
- variant calling evaluates chromatogram peak heights, relative quality, and
  qualifying channels;
- the exclusion-reason vocabulary mixes core and Sanger reasons.

The alignment kernels, difference extraction, anchoring, and sample-level
aggregation are already modality-neutral. The coupling is limited to the
stage entry points and one filter.

A literal "everything is a plugin" would make every internal stage a plugin.
ADR-0058 rejected that, together with a universal plugin trait, because it
adds indirection without removing coupling. Rust also has no stable ABI for
loading plugins at run time.

## Decision

1. **Three coarse plugin families.** DNA composes three families:
   - modality plugins: Sanger now;
   - the core caller: evidence-profile alignment, per-read calling with the
     core gates (`outside_target_region`, `indel_length_exceeded`,
     `non_canonical_allele`, `read_end`), and sample-level aggregation;
   - post-calling plugins: normalization, nomenclature, conformance, analyses.

   Stages inside a plugin stay ordinary modules. Basecalling, signal
   processing, callability, quality control, and the alignment tiers are
   modules, not plugins.
2. **Contracts are data, not traits.** Each family boundary is a versioned
   data contract, so a plugin can live in the same build or in another
   process:
   - modality → core: per-read evidence (`ReadEvidence`);
   - core → post-calling: called variants (`CalledVariantSet`);
   - post-calling results.

   A trait may wrap a contract later, but the contract is the record, not
   the trait.
3. **The modality → core contract.** Per read, `ReadEvidence` carries:
   - the called bases;
   - an optional basecall-independent A/C/G/T profile per call;
   - per call, an optional mask with a reason label and how the call aligns:
     *unresolved* (aligns as `N` without a profile) or *anchoring* (keeps its
     base and profile but supports no variant);
   - per call, the modality's support vetoes, as bits into an ordered
     per-read vocabulary of labelled reasons scoped to substitutions or to
     substitutions and insertions;
   - the informative call interval.

   The core reports modality reasons verbatim, in vocabulary order. It never
   interprets them, and a modality may not reuse a core reason label.
4. **Plugin mechanism.**
   - **In-process plugins** are compile-time Rust code: modules now, crates
     after the split. They are composed by a small kernel that owns the
     contract types, a static registry of plugin descriptors (identity,
     family, provided and required contracts, owned configuration section),
     startup validation of that registry, configuration routing, and
     provenance.
   - **Out-of-process plugins** exchange the versioned JSON form of the
     contracts.
   - No dynamic library loading.
   - Workflows stay fixed compositions chosen by the delivery layer; the
     kernel has no runtime scheduler.
5. **Enforced neutrality.** The module validator rejects any dependency from
   a modality-neutral module on a Sanger module or on a Sanger child of
   `model`. The neutral modules are `read_evidence`, `variant`, `alignment`,
   `variant_calling`, and the post-calling modules. `sample` joins them after
   the sample-attachment phase of PROP-0002.
6. **Phasing.** The architecture lands in phases, each with its own exit
   criterion ([PROP-0002](../../proposals/0002-plugin-first-architecture.md)).
   Phase 1 introduces the evidence contract inside the current crate and must
   leave every result document byte-identical.

## Alternatives

- **Every module a plugin, or one universal plugin trait:** rejected again for
  the reasons in ADR-0058. Coarse families keep the seams where variation is
  real.
- **Dynamic libraries or WASM plugins:** no stable Rust ABI; a WASM runtime
  adds a sandbox, a dependency, and a security surface that no current plugin
  needs. Process plugins over JSON contracts cover cross-language extension.
- **Behaviour traits at the core boundary** (e.g. a support-gate trait the
  core calls back): such a trait cannot cross a process boundary, and it hides
  reason ordering in code. Precomputed vetoes keep the contract inspectable.
- **One universal raw-alignment type for all modalities:** still rejected
  (ADR-0060 §6). `ReadEvidence` is per-read, pre-alignment evidence for
  callers that place reads themselves. A short-read NGS modality would supply
  its own placement or join at `CalledVariantSet`.

## Consequences

- **Phase 1:**
  - The core no longer imports Sanger types: alignment and variant calling
    consume `ReadEvidence`.
  - The Sanger adapter owns the peak, relative-quality, and mixed-signal
    vetoes and the phase-state mask reasons.
  - Published reason labels and every result document are unchanged.
- **Later phases** add a second evidence adapter and a core-only command,
  JSON forms of the contracts with post-calling commands, per-plugin
  configuration sections, and finally a workspace split.
- **Scaling limit:** per-read profile alignment and pairwise overlap admission
  do not scale to short-read depth. "Core usable alone" applies to
  low-read-count evidence: Sanger reads, consensus sequences, long reads.

## Supersession

- [ADR-0058](0058-canonical-contracts-and-modular-analysis-composition.md) §4
  is superseded in part. The three families and their data contracts are
  declared variation seams. "A module is not automatically a plugin" still
  holds for stages inside a plugin.
- [ADR-0060](0060-separate-variant-canonicalization-nomenclature.md) §1 and §6
  are superseded in part. The caller is a modality-neutral core over
  `ReadEvidence`, and modality eligibility arrives as data. `CalledVariantSet`
  remains the convergence boundary for imported calls.
- [ADR-0064](0064-crate-ready-module-layering.md) gains the `read_evidence`
  module and the neutrality rule.
- [ADR-0002](0002-single-crate-layering.md) stays in force until the
  workspace-split phase, which needs its own record.

## Revision 2026-10-09 (phase 2)

The registry is a set of compile-time constants in `src/plugin.rs`, so its
validation runs at compile time instead of at start-up. Constant evaluation
rejects a duplicate identity, a configuration section owned by two plugins, and
a required contract that no earlier plugin of a workflow composition provides.
A test checks that the registry's sections are exactly the shipped
configuration's sections. Results record their workflow's plugins in
`provenance.plugins` with a per-plugin method version
([versioning](../../governance/versioning.md#plugins)). The Sanger support
thresholds moved from `[variant_calling]` to the Sanger-owned
`[sanger_evidence]` section.

## Revision 2026-10-09 (phase 2b)

`sample` joined the neutral modules. Each read after the core is a `CalledRead`
(identity, `ReadEvidence`, alignment, variants) plus a `SangerAttachment`
(calls, signal, callability, quality). Sample aggregation consumes only
`CalledRead`. Locus observations and variant calls keep their source call
index, and the sample report joins Sanger quality, peaks, noisy-region context,
integrity, and callability by read identity and call index. Variant call
resolution has one implementation: the neutral part picks the calls and their
bases, and the report adds the Sanger peaks and quality. The validator now also
applies the neutrality rule to the neutral children of `model`.

## Revision 2026-10-09 (phase 3)

A second modality plugin, `sequence`, turns reviewed consensus FASTA records
into `ReadEvidence`, and `dna call` runs it with the core alone into
`dna.variants/v1`. The evidence contract gains one per-read fact: whether
the modality vouches for the read's calls up to its physical ends. Without that,
the core's read-end margin would treat the ends of a reviewed consensus like the
raw ends of a Sanger read and suppress real differences next to them (for
example the first base of a region). Sanger reads never vouch for their ends,
so their results are unchanged. Sanger observation and `dna call` share one
modality-neutral path, `variant_analysis::read_call`, which places the evidence
and calls its variants.
