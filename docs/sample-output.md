# DNA Sample Evidence JSON

`signal sample <sample-id> <trace.ab1>... --reference <reference.fasta>`
writes one deterministic `results/<sample-id>.sample.json` document identified as
`dna.sample_evidence/v8`. The authoritative schema is
[`schemas/sample-evidence-v8.schema.json`](schemas/sample-evidence-v8.schema.json)
and the example is
[`examples/sample-evidence-v8.example.json`](examples/sample-evidence-v8.example.json).

The sample identifier and read names are reviewer-facing provenance. They never
constrain scientific placement, orientation, overlap discovery, or variant
reconciliation.

## Reads and post-trim coverage

`reads[]` is the one registry of contributing reads. Records are sorted by
SHA-256, so CLI argument order does not change the scientific document.

Each record contains:

- `name`: the AB1 filename stem, for example
  `D11_20260404_LN_26_AB0442_HV1F_11`;
- `sha256`: stable content identity;
- `integrity`: the same concise PLOC/vendor cardinality, PLOC-spacing,
  exact-clipping, and event-signal-scale observations retained by the one-read
  pipeline; this evidence remains read-local and does not by itself admit/reject
  a read;
- `alignment`: the evidence-derived orientation, callable-base count and
  identity, unresolved-base count, gap-open count, mapped reference segments, and
  origin-wrap state.

The scientific pipeline trims each read before alignment. `reference_segments`
therefore describe where the **retained post-trim sequence** aligned on the
reference, not the untrimmed raw call span.

`reference_segments` are 0-based half-open. For a segment
`{"start": S, "end": E}`, the covered 1-based biological positions are
`S + 1` through `E`, inclusive.

For a circular reference, a read can cross the reference origin. In that case
`wraps_origin` is `true` and `reference_segments` contains two segments:
one from the mapped start to the end of the reference and one from reference
position 0 to the mapped end. A normal non-crossing read has
`wraps_origin: false`.

Read names are unique within one emitted sample document because they are used as
human-readable references from overlap, locus, and variant evidence. SHA-256 remains the
scientific content identity. Filename semantics are never used as placement or
merge keys.

## Coverage topology

`coverage[]` is a run-length encoded summary of selected post-trim reference
coverage. Each item contains a 0-based half-open `reference` interval plus:

- `read_depth`: number of independently placed reads covering every coordinate
  in the interval;
- `forward_depth`: covering reads whose selected orientation is forward;
- `reverse_depth`: covering reads whose selected orientation is reverse.

For every item, `read_depth = forward_depth + reverse_depth`. Adjacent intervals
with identical depth tuples are merged even when the identity of the covering
read changes at the boundary; exact read identities and segments remain in
`reads[]`.

Coverage counts all independently placed reads. It does **not** remove a read
because a pairwise overlap is ineligible, and it does not mean canonical-base
agreement, nucleotide comparability, consensus confidence, or biological strand
independence. Deletion columns remain reference-coordinate coverage; insertions
do not create extra reference coordinates.

Circular origin-spanning reads contribute through their two explicit
`reference_segments`, so the linearized JSON coverage map can contain runs near
both reference ends without an implicit wrapped interval.

## Pairwise overlap admission

`overlaps[]` is the **pre-consensus reconciliation layer** learned from Tracy's
explicit minimum-overlap and minimum-match admission checks, adapted to DNA's
N-read reference-coordinate model.

DNA evaluates every unordered pair only after both reads have independently
completed placement. A pair is present when the two selected alignments share at
least one reference coordinate. Non-overlapping reads have no edge; this does not
reject either read and no canonical F/R partner is required.

Each overlap record contains:

- `left` / `right`: reviewer-facing read names from the SHA-sorted registry;
- `shared_positions`: reference coordinates covered by both reads;
- `comparable_bases`: shared coordinates where both aligned query symbols are
  canonical A/C/G/T;
- `agreements` / `conflicts`: equal versus unequal canonical base/base
  observations;
- optional `agreement = agreements / comparable_bases`;
- `eligible`: whether the edge meets the configured pre-consensus gates;
- `exclusion_reasons`: exact failed rules.

The configured defaults are Tracy-derived:

~~~text
minimum_comparable_bases = 25
minimum_overlap_agreement = 0.50
~~~

DNA intentionally does **not** copy a base-vs-gap scalar match rule.
Deletions and unresolved symbols remain part of shared coverage but do not enter
the nucleotide agreement denominator. Insertions and deletions remain explicit in
the existing locus/normalized-variant evidence so future consensus can treat gap
support as a separate evidence problem.

An eligible edge has no exclusion reasons. An ineligible edge records
`comparable_bases_below_minimum`, `agreement_below_minimum`, or both. This decision
does not erase a read, mutate placement, or change read-level variant eligibility.

## Sparse locus differences

`locus_differences[]` is the **alignment-observation layer**. It contains only
reference positions where at least one covering read is not a canonical reference
match. It answers:

> What did each covering read actually observe at this reference coordinate?

Each retained locus also exposes `support_topology` derived from the same
authoritative observations:

- total and forward/reverse read counts;
- reference/alternate/unresolved/deletion read counts;
- profile-bearing read count and its forward/reverse partition.

These counts are evidence topology only. They are not a vote, confidence score,
or claim of biological independence.

At a retained locus, `observations[]` contains every read that covers that locus,
including reads that agree with the reference. A call-backed observation contains:

- `read`: human-readable read name;
- `state`: `reference`, `alternate`, or `unresolved`;
- `base`: reference-oriented observed base;
- `quality`: the existing uncalibrated relative quality score;
- optional `profile`: normalized basecall-independent A/C/G/T evidence already
  projected to reference orientation;
- `in_noisy_region`: whether the source call falls inside the existing merged
  candidate-noisy call context.

A valid zero-signal call has no `profile`; DNA does not invent a called-base,
uniform, or reference-derived replacement. A deletion observation contains only
`read` and `state: "deletion"`; DNA does not fabricate deleted-base
quality, profile, or noisy-call context.

Dense all-reference positions are omitted. The compact default is explicit:

- inside a read's mapped `reference_segments`, absence from
  `locus_differences[]` means that read is a canonical reference match there;
- outside the read's mapped segments, the position is uncovered;
- at a retained differential locus, the explicit observations are authoritative.

This preserves reference support versus missing coverage without serializing
routine reference loci one by one.

Inserted query bases have no reference-coordinate locus and therefore do not
create a synthetic `locus_differences` entry.

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

## Why coverage, overlaps, locus_differences, and variants are separate

The four arrays intentionally preserve different evidence layers:

```text
selected per-read alignments
      ↓
coverage[]              local mapped-read denominator/orientation topology
      ↓
overlaps[]              which mapped read pairs are eligible for later reconciliation
      ↓
locus_differences[]     what each read observed at differential reference coordinates
      ↓
normalization/filtering
      ↓
variants[]              normalized alleles + eligibility + trace evidence
```

`locus_differences[]` is useful for disagreement and coverage reasoning,
including reference-vs-alternate or unresolved evidence. `variants[]` is the
reviewer-facing normalized biological call layer and is where per-read peak
evidence belongs.

None of these arrays is a consensus result.

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

## Contract boundary

v8 remains compact and difference-focused. It does not serialize per-base
evidence for loci where every covering read agrees with the reference. Pairwise
overlap records summarize only admission-relevant counts rather than dense
per-coordinate comparisons. The scientific pipeline still processes each read
independently before sample aggregation.

The current implementation emits v8 only. Earlier sample-evidence contracts are
not emitted as aliases or compatibility output.

## Non-goals

The v8 contract contains no consensus sequence, sample-level adjudicated variant
verdict, majority-vote result, genotype, heteroplasmy estimate, haplogroup
interpretation, F/R pair object, primer/HV placement rule, or filename-derived
placement. `overlaps[]` is an evidence/admission graph, not a pair-first merge
structure.
