# ADR-0060: Separate variant calling, canonicalization, and nomenclature

- **Status:** Accepted
- **Date:** 2026-10-03
- **Related decisions:** [ADR-0047](0047-canonical-right-aligned-mtdna-gaps.md), [ADR-0057](0057-haplotype-correctness-and-variant-nomenclature.md), [ADR-0058](0058-canonical-contracts-and-modular-analysis-composition.md)

## Context

DNA already owns Sanger variant calling. The current production caller extracts
primary-sequence differences from the selected alignment, attaches direct
trace-call provenance, applies eligibility gates, and constructs anchored
reference/alternate alleles.

The legacy mtDNA Tracy path demonstrates a second concern after calling:
sequence-equivalent edit sets may need deterministic canonicalization and then
mtDNA-specific nomenclature. In that path, Tracy's own `variants.rows` are
diagnostic only; final calls are reconstructed from chromatogram/alignment
evidence, normalized by reconstructed haplotype, and then passed through mtDNA
representation policy before multi-read consensus.

These concerns should not remain coupled to one Sanger aligner. Future NGS
callers, imported variant sources, or another validated modality may reach the
same biological variant boundary without using DNA's current pairwise alignment
implementation.

At the same time, Sanger alignment still needs a deterministic traceback for
evidence placement and provenance. ADR-0047 therefore remains authoritative for
alignment topology, but its earlier assumption that final variant representation
must always inherit that exact gap placement is too restrictive for a
modality-independent canonical-variant layer.

## Decision

DNA separates three responsibilities and does not require every workflow to
execute all three:

~~~text
variant calling / importing
          |
          v
   called variants
      /        \
     v          v
direct use   optional normalization/canonicalization
                    |
                    v
             optional nomenclature
~~~

The called-variant boundary is the first cross-modality convergence point.
Normalization/canonicalization and nomenclature are explicitly selected
capabilities, not mandatory hidden stages.

### 1. Variant calling owns evidence-supported biological differences

A variant caller consumes modality-appropriate interpreted evidence and produces
called differences with enough provenance to audit the biological observation.

For the current Sanger implementation this includes:

- the selected reference-oriented alignment;
- SNV/indel extraction;
- direct supporting/flanking trace-call mappings;
- configured evidence eligibility;
- explicit excluded-candidate diagnostics.

Variant calling MUST NOT own mtDNA-specific nomenclature rules.

### 2. Alignment canonicalization and variant canonicalization are separate contracts

The current Sanger aligner continues to canonicalize repeat-equivalent optimal
tracebacks under ADR-0047. That decision is required for deterministic evidence
placement, call/reference mapping, and reproducible Sanger alignment output.

Post-calling variant canonicalization is a separate boundary. It operates on the
called biological changes and reconstructed haplotype, not on raw traceback
state.

The two boundaries may use the same sequence-equivalence machinery, but neither
is defined as an implementation detail of the other.

### 3. Canonicalization is haplotype-preserving

Post-calling canonicalization MUST reconstruct the represented alternate
haplotype and MUST preserve it exactly.

For a reference `R`, source edit set `S`, and canonical edit set `C`:

~~~text
apply(R, S) == apply(R, C)
~~~

Canonicalization may change event coordinates or decomposition only when the
complete represented haplotype remains identical and source provenance remains
recoverable.

The canonicalization result MUST preserve enough immutable context for later
nomenclature and audit, conceptually including the source edit set, reconstructed
alternate haplotype, selected canonical edit set, and source-to-canonical
provenance. Canonical movement must not destructively overwrite the source
description.

This is the production counterpart of the haplotype-first comparison principle
in ADR-0057.

### 4. Right alignment is reusable machinery, not a universal DNA policy

For the current human-mtDNA target, sequence-equivalent indels use the
established 3'/right-most rCRS light-strand convention, with the canonical rCRS
origin seam preserved.

The algorithm for shifting sequence-equivalent indels can be generic and
reusable. Selecting right-most placement is target/method policy, not a universal
requirement for all future nuclear, NGS, VCF-interchange, or other workflows.

A future target may select another explicit canonicalization policy without
changing variant-calling semantics.

### 5. Nomenclature is an optional target-specific representation layer

Nomenclature is not required for every workflow. A composition layer may omit
it when the consumer needs biological calls rather than a named target-specific
representation.

When a nomenclature policy is composed after generic haplotype-preserving
canonicalization, it consumes the canonicalization result rather than merely a
destructively rewritten final edit list, so target rules may inspect the
preserved source description and reconstructed haplotype when scientifically
required. A policy MAY also accept an unnormalized representation context when
its own contract explicitly defines and validates that path.

For mtDNA, nomenclature may select validated representations for unstable or
repeat-rich regions such as the current 309/315, 513-524, and 16189/16193
families, provided the represented haplotype is unchanged. A nomenclature rule
may therefore supersede the generic right-aligned representation inside its
explicit validated window while preserving the same alternate haplotype.

Nomenclature MUST NOT reinterpret chromatogram signal, change evidence
eligibility, manufacture phase, or alter the biological haplotype merely to
obtain a preferred name.

Any Sanger-specific artifact interpretation that can change the inferred
biological sequence belongs before canonicalization, in the modality-specific
evidence/read-interpretation path.

### 6. Called variants are the convergence boundary

The intended architecture is:

~~~text
Sanger evidence --> Sanger caller --------------------+
                                                      |
future NGS evidence --> NGS caller -------------------+--> CalledVariantSet
                                                      |        |       |
VCF/BCF or future called source --> importer/adapter -+        |       |
                                                               |       +--> direct downstream
                                                               |
                                                               +--> optional normalization
                                                                        |
                                                                        +--> optional nomenclature
~~~

A future imported VCF/BCF source therefore bypasses raw calling and adapts into
the called-variant boundary without making VCF records DNA's canonical domain
model.

The input to a caller remains source-specific. DNA does not require a Sanger
pairwise alignment and an NGS aligned-read collection to share one artificial
universal raw-alignment type.

### 7. Workflows choose where representation convergence is required

Where multiple reads/traces describe one sample, a workflow SHOULD prevent
sequence-equivalent representations from appearing biologically discordant.

For the current mtDNA direction, applying the selected mtDNA normalization and
nomenclature policy before representation-based sample reconciliation is one
valid composition. A future sample algorithm may instead compare reconstructed
haplotypes directly and therefore need less pre-reconciliation representation
policy.

True biological or evidence disagreement remains explicit and MUST NOT be hidden
by normalization or nomenclature.

### 8. Notation and serialization remain outer representation concerns

Mechanical rendering into mtDNA position notation, VCF/BCF, JSON, or another
external representation is separate from biological canonicalization and
nomenclature policy.

An external serializer MAY require its own representational convention, such as
VCF-compatible anchoring or left normalization, but that convention MUST NOT
silently redefine DNA's internal biological result.

## Current implementation status

This ADR defines the architectural direction; it does not claim all stages are
implemented.

Current production truth now includes:

- Sanger variant calling in `src/variant_calling`;
- deterministic Sanger alignment right-most gap placement under ADR-0047;
- caller construction that preserves the alignment-selected placement;
- a typed public `CalledVariantSet` projection from Variant Analysis;
- an optional standalone `variant_normalization` capability with the initial
  human-mtDNA 3'/right-most policy, source-variant preservation, and complete
  haplotype-equivalence checks;
- a separate public `variant_nomenclature` capability with the initial
  human-mtDNA HVS-II 309/315 poly-C representation rule, preserving source,
  normalized, and represented variants plus the complete alternate haplotype.

Additional target-specific windows such as HVS-III 513-524 and HVS-I
16189/16193 remain unimplemented and require their own tests and validation.

**Revision (2026-10-08):** the `sample` workflow now composes normalization and
the HVS-II rule per read against the rCRS, as §7 permits, and publishes
per-base notation from the outer report layer, as §8 requires
(`dna.sample_evidence/v9`, SRS-NOM-010 to SRS-NOM-012). Analysis output and the
public API are unchanged.

## Consequences

- DNA's variant caller remains a real scientific module rather than a thin
  serializer around an external tool.
- Right alignment can become reusable sequence-equivalence machinery instead of
  being reachable only through the Sanger aligner.
- Future NGS callers can converge at the called-variant boundary and reuse the
  same mtDNA canonicalization/nomenclature path.
- mtDNA special-region naming policy remains outside generic variant calling.
- Provenance must survive movement from source calls to canonical calls.
- Implementation will require explicit intermediate typed contracts rather than
  treating one final `Variant` struct as every stage of the lifecycle.

## Validation requirements

Implementation MUST cover at least:

- exact SNV calls that require no movement;
- homopolymer insertion/deletion right alignment preserving the complete
  haplotype;
- tandem-repeat movement preserving the complete haplotype;
- multiple nearby edits where independent local shifting would change phase;
- source-call provenance retained after canonical movement;
- no movement across the canonical rCRS origin seam;
- an mtDNA nomenclature-window case where generic canonicalization and final
  nomenclature differ but reconstruct the same haplotype;
- forward/reverse Sanger descriptions converging before sample reconciliation;
- a genuinely different haplotype that MUST remain different.

## Supersession relationship

ADR-0047 remains authoritative for deterministic Sanger alignment topology.

This ADR supersedes ADR-0047 only in the narrower assumption that final
canonical variant representation must always inherit the alignment-selected gap
placement and that positional variant canonicalization has one authoritative
home inside `alignment::canonical`.
