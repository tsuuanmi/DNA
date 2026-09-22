# ADR-0019: Separate signal evidence, read interpretation, and biological claims

## Status

Accepted

## Context

A Sanger chromatogram contains richer evidence than a single primary base string. DNA already preserves selected per-call peak evidence and local SNR internally, while intentionally refusing to call genotype, heteroplasmy fraction, pathogenicity, or clinical significance.

Future mtDNA work will add richer peak geometry, multiple reads per sample, bidirectional evidence, mixed-signal detection, consensus, haplogroup QC, and mtDNA-facing variant notation. Without an explicit evidence hierarchy, those features could accidentally turn an observation into a stronger biological claim than the data support.

This is especially important for Sanger data:

- fluorescence peak heights depend on dye chemistry, local sequence context, and background and therefore are not direct molecule fractions;
- low-level secondary signal is context- and noise-dependent;
- repeat and poly-C regions can create alignment and length-mixture ambiguity;
- a difference observed in one chromatogram is not automatically a sample-level variant.

## Decision

DNA will model scientific information in explicit layers. A downstream layer may summarize an upstream layer, but it MUST NOT erase the distinction between observation and interpretation.

### Layer 0: source evidence

Authoritative source evidence consists of the decoded analyzed A/C/G/T channel arrays plus validated ABIF metadata required to locate calls.

Rules:

- decoded channel arrays are immutable;
- baseline correction, smoothing, denoising, normalization, or locus refinement creates a separate versioned derived representation;
- raw evidence is never overwritten to make a downstream interpretation easier to express.

For the current method, `PLOC.2` remains the locus authority. A future locus-refinement method must be versioned, independently validated, and retain the original PLOC coordinate.

### Layer 1: per-locus observations

A call locus may contain:

- channel peak heights and positions;
- offsets/spacing/width/prominence when implemented;
- corrected amplitude and local SNR when implemented;
- primary and ambiguity calls;
- vendor agreement;
- relative or calibrated confidence.

These are observations about the chromatogram.

A peak-height ratio or local SNR MUST NOT be described as an mtDNA allele fraction or heteroplasmy fraction unless a separately validated quantitative method establishes that interpretation for a defined acquisition domain.

### Layer 2: read interpretation

A single read may produce:

- retained/trimmed interval;
- selected alignment and orientation;
- read coverage;
- normalized primary-sequence differences;
- mixed-signal candidates;
- mapping and quality warnings.

The preferred term for a disagreement derived from one chromatogram is `primary_sequence_difference` or `read_variant_observation`, not a sample-level genotype assertion.

A two-channel or secondary-signal event remains `mixed_base_candidate` or `secondary_signal` until stronger evidence exists.

### Layer 3: sample evidence

Sample-level conclusions require explicit aggregation of independently acquired read evidence.

The future sample model should retain:

- every contributing read identity;
- sequencing direction;
- expected amplicon/region when known;
- per-position coverage and callability;
- concordance/discordance;
- independent support level;
- consensus state.

The system MUST distinguish:

```text
uncovered
covered but no reliable call
single-read support
single-direction support
bidirectional support
discordant support
```

Absence of a reported difference MUST NOT be interpreted as evidence for the reference allele at positions that were not demonstrably covered and callable.

Conflicting high-quality reads should produce an explicit discordant/review state rather than an arbitrary majority answer.

### Layer 4: mtDNA-specific projection

Internal canonical variant representation and mtDNA-facing nomenclature are separate concerns.

For human mtDNA projection:

- the reference identity and version MUST be explicit;
- rCRS / RefSeq `NC_012920.1` is the preferred human mitochondrial reference for HGVS-style `m.` descriptions;
- circular coordinate handling and canonical internal normalization remain algorithmic concerns independent of display notation;
- haplogroup knowledge may be used for QC or review support, but MUST NOT silently rewrite chromatogram evidence.

### Mixed signal and heteroplasmy

Until DNA has a validated method with an explicit study design, domain, limit of detection, repeatability assessment, and false-positive characterization:

- mixed peaks are candidates, not quantitative heteroplasmy calls;
- peak ratios are descriptive features, not molecule fractions;
- one-strand-only evidence is visibly distinct from bidirectional confirmation;
- persistent post-indel phase-shift evidence is distinguished from a simple point mixed base;
- contamination, mixture, artifact, and heteroplasmy are alternative explanations rather than labels inferred from one trace.

A future heteroplasmy method requires its own ADR, SRS changes, validation protocol, and output-contract version.

### Repeat and indel context

Homopolymers, mtDNA poly-C tracts, and repeat-associated indels are high-risk contexts.

DNA should retain sequence-context annotations and observed read mappings separately from the normalized reported allele. Normalization MUST NOT rewrite or fabricate the evidence locations from which the event was observed.

### Validation implications

Scientific validation should be stratified by evidence regime rather than reporting only one aggregate accuracy number.

The approved real-trace corpus should eventually cover, where available:

- forward and reverse reads;
- multiple instruments/runs or explicitly documented acquisition domains;
- HVI/HVII/HVIII and representative coding-region reads;
- strong and weak signal;
- read ends;
- homopolymers and poly-C tracts;
- known SNVs and indels;
- mixed-base examples and clean negative controls;
- origin-spanning/circular-coordinate cases when relevant.

Every truth item must record how truth was established. DNA output itself is not an independent truth source for training or validation.

## Consequences

### Positive

- Future ML, consensus, and haplogroup features can consume richer evidence without overstating biology.
- Read-level and sample-level semantics become explicit.
- Quantitative heteroplasmy remains blocked until the project has evidence to justify it.
- Review tooling can expose why a call exists instead of presenting only a final base.
- mtDNA nomenclature can evolve without contaminating generic alignment and normalization logic.

### Cost

- Internal models become richer than the compact public JSON.
- Sample-level processing needs explicit manifests and evidence aggregation rather than filename heuristics.
- Validation corpus design becomes a first-class engineering task.
- Some user-facing terms must remain conservative even when a stronger interpretation appears plausible.

## References

- Tracy: basecalling, alignment, assembly and deconvolution of Sanger chromatogram trace files (BMC Genomics, 2020).
- HGVS Sequence Variant Nomenclature, mitochondrial reference-sequence guidance.
- ClinGen/ACMG specifications for mitochondrial DNA variant interpretation.
- Reviews and collaborative studies of mtDNA heteroplasmy detection in Sanger electropherograms.
