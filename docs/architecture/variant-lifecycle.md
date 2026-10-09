# Variant Lifecycle

This document defines the architectural boundaries from source-specific evidence
to variant representations. It distinguishes current production behavior from
planned capabilities.

## Core rule

DNA does not require every input type to pass through one universal alignment or
one universal variant-calling implementation.

The first cross-modality convergence boundary is an evidence-backed called
variant set:

~~~text
source-specific evidence
        |
        v
source-specific caller or importer
        |
        v
    CalledVariantSet
        |
        +--------------------------> direct downstream consumer
        |
        +--> optional normalization/canonicalization policy
        |         |
        |         v
        |   NormalizationResult
        |         |
        |         +---------------> downstream consumer
        |         |
        |         +--> optional nomenclature policy
        |                   |
        |                   v
        |          represented variant set
        |
        +--> optional nomenclature policy when that policy accepts
             the unnormalized representation context
~~~

A workflow selects only the capabilities it needs.

## 1. Calling variants from an alignment

An alignment is evidence geometry, not a final variant representation.

For the current Sanger pairwise path, the selected alignment provides
reference-oriented columns plus mappings back to trace calls. The caller:

1. walks the selected alignment columns;
2. groups contiguous query gaps into deletion events;
3. groups contiguous reference gaps into insertion events;
4. converts canonical reference/query mismatches into SNV events;
5. attaches the supporting or flanking source-evidence mappings;
6. applies caller-owned evidence eligibility;
7. constructs a valid anchored REF/ALT representation without moving the event
   through a repeat.

Conceptually:

~~~text
reference + selected alignment + source evidence
                    |
                    v
             difference extraction
                    |
                    v
        evidence-supported source events
                    |
                    v
              allele anchoring
                    |
                    v
              CalledVariantSet
~~~

The caller answers **what biological difference is supported by this evidence**.
It does not answer **which sequence-equivalent coordinate representation should a
later workflow prefer**.

Right/left shifting, repeat representation policy, mtDNA naming conventions, and
external-format normalization therefore do not belong to raw difference
extraction.

## 2. Input-specific paths

Different input types reach the called-variant boundary through different
scientific paths.

| Input/source | Source-specific evidence path | Variant boundary |
|---|---|---|
| Sanger ABIF | ABIF -> chromatogram -> base calling / signal / callability / QC -> selected pairwise alignment | Sanger alignment caller -> `CalledVariantSet` |
| assembled/consensus FASTA sequence | sequence -> reference alignment | sequence-difference caller -> `CalledVariantSet` |
| FASTQ / NGS reads | reads -> QC/preprocessing -> mapping -> read/depth evidence | NGS caller -> `CalledVariantSet` |
| BAM / CRAM | validated aligned-read input -> read/depth/mapping evidence | NGS caller -> `CalledVariantSet` |
| VCF / BCF | validated external called-variant representation | variant importer/adapter -> `CalledVariantSet`; raw calling is bypassed |
| future evidence source | source-specific adapter/analysis | source-specific caller or importer -> `CalledVariantSet` |

These rows are architectural paths, not claims that all sources are currently
implemented. Current production support is Sanger ABIF plus FASTA reference.

### Why there is no universal raw-alignment contract

A Sanger pairwise alignment and an NGS read-alignment collection do not carry
the same biological evidence:

- Sanger has chromatogram-derived per-locus signal and one selected read
  placement;
- NGS may have many reads, base qualities, mapping qualities, CIGAR operations,
  strand support, depth, duplicate state, and caller-specific statistical
  evidence.

DNA therefore does not flatten them into a fake universal `RawAlignment`
object. Each caller consumes the smallest source-specific evidence contract it
needs and emits a common called-variant contract.

## 3. Called variants are the convergence boundary

A called variant is an evidence-backed biological difference relative to a
specific reference identity.

The common boundary must preserve enough information to support later
representation transforms and audit. Conceptually it includes:

~~~text
CalledVariant
├── reference identity
├── source event / edit semantics
├── source representation
├── evidence provenance
└── caller / method provenance
~~~

The current exact Rust shape is owned by
[Rust public API](../reference/rust-api.md); architecture owns only its semantic
boundary.

A `CalledVariantSet` is not required to be right-aligned, left-aligned,
HGVS-formatted, VCF-normalized, or mtDNA-nomenclature-normalized.

## 4. Normalization/canonicalization is optional policy

Normalization is a haplotype-preserving representation transform.

A policy receives called variants plus the reference and may choose another
sequence-equivalent edit representation:

~~~text
CalledVariantSet + Reference + NormalizationPolicy
                         |
                         v
                NormalizationResult
~~~

The invariant is:

~~~text
apply(reference, source_edits)
    ==
apply(reference, normalized_edits)
~~~

The result preserves the source edits and provenance rather than destructively
overwriting them.

Examples of explicit policies may include:

- human-mtDNA 3'/right-most sequence-equivalent indel placement;
- a future VCF interchange left-normalization policy;
- another target-specific normalization policy;
- no normalization at all.

Right alignment is therefore reusable sequence-equivalence machinery selected by
the current mtDNA policy, not a universal rule for every DNA workflow.

## 5. Nomenclature is an optional representation layer

Nomenclature is also optional.

It selects a representation according to a named target/domain convention while
preserving the biological haplotype. It must not manufacture evidence or
reinterpret source signal.

For the current mtDNA direction:

~~~text
CalledVariantSet
      |
      +--> mtDNA right-alignment policy        optional
      |          |
      |          v
      |   NormalizationResult
      |          |
      |          +--> mtDNA nomenclature       optional
      |                    |
      |                    v
      |             represented mtDNA variants
      |
      +--> direct biological/sample analysis   allowed
~~~

The legacy mtDNA Tracy path demonstrates why nomenclature may need the preserved
source edit set and reconstructed haplotype in addition to the generic normalized
edits. Validated repeat-window policy can choose a preferred mtDNA
representation while proving that the represented sequence is unchanged.

Nomenclature is therefore not synonymous with variant calling and is not a
mandatory property of the core called-variant contract.

## 6. Downstream consumers choose the representation they require

Consumers should request the smallest contract they need.

Examples:

- evidence review may consume called variants directly;
- mtDNA sample reconciliation may choose normalized/nomenclature-aware variants
  to avoid representation-only discordance;
- haplogroup analysis may consume a canonical biological representation;
- VCF export may apply its own explicit interchange representation policy;
- debugging may retain source and normalized representations side by side.

No serializer or nomenclature system silently redefines the underlying
biological call.

## Current production truth

Today DNA implements:

~~~text
Sanger ABIF
  -> Sanger evidence
  -> selected pairwise alignment
  -> variant_calling
  -> current Variant Analysis result
~~~

The current Sanger aligner also selects deterministic right-most
repeat-equivalent traceback topology for evidence placement under ADR-0047, and
the current variant builder preserves that selected placement.

DNA now also exposes `VariantAnalysisResult::called_variants()` as the typed
`CalledVariantSet` boundary and an optional standalone
`variant_normalization` capability. The `RightAligned` policy preserves source
calls and the complete reconstructed haplotype while selecting
sequence-equivalent 3'/right-most indel placement without crossing the FASTA
seam.

DNA also exposes the immutable `variant_nomenclature` input seam from a
`VariantNormalizationResult`, carrying reference identity, exact source calls,
the reconstructed alternate haplotype, and normalized variants together.

Target-specific nomenclature is data: a target profile (ADR-0063) declares
reference windows with ordered rules, and the generic engine applies the first
rule that reproduces each window haplotype, verifying that the complete
represented haplotype is unchanged. The shipped human-mtDNA profile declares the
HVS-II 303-315, HVS-III 513-524, and HVS-I 16181-16193 windows.
