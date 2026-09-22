# ADR-0033: Factor normalized-variant support by eligibility and orientation

- **Status:** Accepted
- **Date:** 2026-09-19

## Context

DNA sample evidence already preserves every independently placed read, a
reference-coordinate coverage map, pairwise overlap/admission evidence, sparse
differential loci, and each read that observes a normalized variant.

The Tracy source audit highlights that a single ordinal support label loses
important topology. For example, two reads can represent two selected
orientations or only one, and those patterns should not be collapsed into one
"support level".

DNA does not yet have authoritative amplicon or replicate metadata, so those
dimensions cannot be inferred from filenames or invented by the sample layer.

## Decision

Each normalized sample variant gains a deterministic `support_topology`
summary derived from its existing per-read `support[]` records and the selected
orientation stored in the sample read registry.

The summary contains:

~~~text
reads
eligible_reads
forward_reads
reverse_reads
eligible_forward_reads
eligible_reverse_reads
~~~

The following invariants hold:

~~~text
reads = forward_reads + reverse_reads
eligible_reads = eligible_forward_reads + eligible_reverse_reads
eligible_reads <= reads
eligible_forward_reads <= forward_reads
eligible_reverse_reads <= reverse_reads
~~~

Every read may contribute at most once to one normalized variant identity.

### What is counted

`reads` counts reads that **observed this normalized variant**, including
observations that failed configured single-read reporting eligibility.

The eligible counts are the subset whose existing
`VariantSupport::eligible` value is true.

Orientation counts use the selected evidence-derived read orientation already
stored in the read registry. They do not re-run placement and do not infer
direction from filenames.

### What is not counted

The summary does not count reads that cover the locus but support the reference,
are unresolved, or observe a competing event unless they independently have a
`support[]` record for this exact normalized variant.

Those local denominator/opposition states remain available through
`coverage[]` and `locus_differences[]`. Insertion/deletion competition remains
a distinct evidence problem.

### Interpretation boundary

Forward/reverse counts are support topology, not proof of biological strand
independence. A forward and reverse chromatogram can share upstream assay
artifacts.

Likewise, `eligible_reads` is not a calibrated confidence score, vote weight,
or sample-level variant verdict. It only summarizes the existing configured
single-read eligibility state.

No amplicon count, replicate count, or independence label is emitted until such
metadata has an explicit authoritative input contract.

### Public contract

`dna.sample_evidence/v6` is replaced by
`dna.sample_evidence/v8`.

Every item in `variants[]` now requires `support_topology` alongside the
authoritative per-read `support[]` records.

The v6 schema/example is removed; no compatibility alias or duplicate output is
retained.

Because JSON Schema cannot express all cross-field/read-registry relationships,
the result-contract validator independently recomputes topology from
`reads[] + variants[].support[]` and rejects inconsistent documents.

## Consequences

- Reviewers and future sample interpretation can distinguish read-count,
  eligibility, and selected-orientation dimensions without reconstructing them.
- Filtered observations remain visible in both the per-read evidence and the
  topology.
- Support topology is deterministic and derived from existing evidence; it adds
  no new threshold.
- Future consensus/sample-variant policy can consume this evidence but must
  define opposition, gap handling, contributor admission, and calibration in a
  separate decision.
