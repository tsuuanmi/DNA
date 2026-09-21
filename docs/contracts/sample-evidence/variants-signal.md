# Variants and Differential-locus Signal Evidence

Part of the canonical [sample evidence contract](../sample-evidence.md).

## Variants

`variants[]` is the **normalized variant layer**. It answers a different
question:

> Which normalized biological variant was observed, by which reads, and how
> strong was the supporting trace evidence?

A locus difference is not automatically the same thing as a normalized variant.
For example:

- an unresolved aligned base can exist in `locus_differences[]` without becoming
  a canonical variant;
- an insertion exists in `variants[]` even though inserted bases have no
  reference-coordinate locus;
- indel normalization can move the reported normalized allele representation away
  from the exact alignment gap;
- a normalized variant remains in sample evidence even if a read-level reporting
  filter marks that read's support ineligible.

Variants aggregate by `(position, reference, alternate, kind)`. Each variant also contains `support_topology`, a deterministic summary of
the existing per-read support records:

- `reads`: reads that observed this exact normalized variant, including
  ineligible observations;
- `eligible_reads`: the subset passing existing single-read variant
  eligibility;
- `forward_reads` / `reverse_reads`: observing reads grouped by selected
  evidence-derived orientation;
- `eligible_forward_reads` / `eligible_reverse_reads`: the eligible subset
  within each selected orientation.

The counts are recomputed from `support[]` and the read registry by contract
validation. They do not include covering reads that support the reference, are
unresolved, or observe another event. Those local denominator/opposition states
remain in `coverage[]` and `locus_differences[]`.

Orientation support is not a claim of assay independence, and eligible-read
count is not a probability, confidence score, vote weight, genotype, or
heteroplasmy fraction. DNA has no authoritative amplicon/replicate input
contract yet, so those dimensions are not inferred from filenames.

Each support
record contains:

- `read`: the human-readable read name;
- `eligible`: whether that read's variant observation passes configured
  reporting eligibility;
- `exclusion_reasons`: exact failed configured rules when ineligible;
- `calls[]`: reviewer-facing trace evidence.

Each call intentionally omits original call index, aligned call position, and ABIF
PLOC because those implementation coordinates do not help routine variant review.
Instead it contains:

- `role`: `supporting` or `flanking`;
- `base`: the called base projected onto the reference strand;
- `peaks`: raw analyzed A/C/G/T channel heights sampled together at the
  uniquely strongest primary-event coordinate, also projected to reference
  orientation;
- `quality`: the uncalibrated relative quality score for that call.

For reverse reads, both `base` and the A/C/G/T peak labels are
reference-oriented. A reviewer can therefore compare the variant allele directly
with the strongest channel without mentally reverse-complementing the trace.

Supporting calls carry the observed alternate or inserted base evidence. Indels
can also carry flanking calls because a deletion has no signal at the deleted
reference base and an insertion is bounded by aligned reference bases.

An eligible support has an empty exclusion list. An ineligible support retains one
or more reasons such as `outside_configured_region`, `peak_below_minimum`,
`relative_quality_not_above_threshold`, or `mixed_supporting_dna`.
`mixed_supporting_dna` means an SNV's supporting call retained more than one
co-localized qualifying channel under the authoritative basecalling rule; the
normalized observation remains evidence, but it is not presented as a clean SNV. The latter name remains explicit because
the configured gate still operates on the internal relative-quality method even
though the public numeric field is simply `quality`.

## Differential-locus signal evidence in v8

v8 promotes the smallest reviewer-useful part of the internal sample signal
evidence into `locus_differences[]`: the normalized A/C/G/T
`EvidenceProfile`, existing noisy-region membership, and factorized support
topology.

The richer quantitative evidence remains internal: corrected amplitudes,
per-channel SNR, aggregate mean profiles, profile heterogeneity decomposition,
and directional profile distance are not part of the compact public contract.
That keeps routine sample JSON focused while preserving enough chromatogram
shape to distinguish a clean cross-read disagreement from mixed signal within a
read.
