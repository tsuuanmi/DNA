# DNA Pipeline

This document describes the shared read-processing stages used by `signal analyze`, `signal basecall`, and `signal sample`. `basecall` stops after quality control; `analyze` produces one read analysis; `sample` independently produces one `ReadObservation` per trace and then aggregates sample evidence. It is the authoritative description of the current Rust behavior. Every stage,
substep, and formula below is derived from the source under `src/`; where this
document and the source disagree, the source is ground truth and this document
should be corrected.

The scientific pipeline is deterministic: the same required inputs, effective
configuration, command, and DNA version always produce the same result. Each
command creates one JSON document with its command-specific derived suffix. A
separate nondeterministic append-only operational log records aggregate stage
progress and failures without entering the JSON contract.

## Overview

```text
AB1 + TOML ──► decode ──► basecalling ──► signal_processing ──► quality_control
                                                                     ├─► basecalls/v2
FASTA reference ─────────────────────────────────────────────────────┴─► alignment ─► variant_calling
                                                                                              │
                                                                                       ReadObservation
                                                                                         ├─► analysis/v7
                                                                                         └─► sample aggregation
```

`analyze` and `basecall` consume exactly one AB1 trace. `sample` consumes one or more AB1 traces and processes each independently through the same reference-guided observation path. `analyze` and `sample` additionally consume one single-record FASTA reference; `basecall` performs no reference I/O. Each
stage consumes the validated output of the previous stage and produces a new
typed result; no stage mutates shared state.

## Inputs

- **Trace:** one regular ABIF/AB1 file, decoded into four A/C/G/T signal
  channels, basecall positions (`PLOC.2`), and optional vendor evidence
  (`PBAS.2`, `PCON.2`). `P2BA.1` is ignored. Vendor base strings retain uppercase
  IUPAC symbols, and PCON accepts the ABIF one-byte byte or char representation.
- **Reference (`analyze` and `sample`):** one plain FASTA record of A/C/G/T/N bases, up to 50,000 bases, interpreted as linear or circular per configuration. `basecall` does not accept or load a reference.
- **Configuration:** one strict TOML file selected by `DNA_CONFIG` or
  `config/dna.toml`. Unknown keys, missing sections, and out-of-range values
  are errors.

## Stage 1 — Decode (`dna.abif_decode/v1`)

Parses the ABIF container and validates every directory entry, offset, element
size, and element count before access. It extracts:

- the four `DATA.9`–`DATA.12` channels as signed 16-bit samples, reordered into
  canonical A/C/G/T order using the `FWO_.1` channel-order string;
- the `PLOC.2` basecall positions (strictly increasing, within the sample
  range);
- optional vendor base strings and quality values. Their decoded cardinality may
  differ from PLOC and is retained as trace-integrity evidence rather than
  changing the PLOC-defined call series.

The decoded chromatogram records the source file name and SHA-256, the canonical
four channel arrays, the basecall positions, and optional vendor evidence. ABIF
version, channel order, and sample count are validated during decode but are not
duplicated as retained metadata.

## Stage 2 — Basecalling (`dna.peak_recall/v3`)

Re-calls every vendor-defined locus from the channel signals. Vendor base
strings are retained as evidence but never replace signal-derived re-calling.

### Substep 2.1 — Call windows

For each basecall position `p[i]`, a half-open sample window is built from the
midpoints of neighboring positions:

- first window: `[p[0] - (p[1]-p[0])/2, midpoint(p[0], p[1]))`;
- interior window `i`: `[midpoint(p[i-1], p[i]), midpoint(p[i], p[i+1]))`;
- last window: `[midpoint(p[n-2], p[n-1]), p[n-1] + (p[n-1]-p[n-2]+1)/2)`,
  clamped to the sample count.

`midpoint(a, b) = a + (b - a) / 2`. At least two basecall positions are
required.

### Substep 2.2 — Per-channel peak selection

Within each window, each channel is searched for a positive local maximum. A
sample `v` at position `j` is a local maximum when
`(v[j-1] <= v && v > v[j+1]) || (v[j-1] < v && v >= v[j+1])`. The highest such
sample is the channel peak. If no positive local maximum exists, the channel
value at the basecall position is used as a fallback. Each channel peak records
its base, height, position, and source (`local_maximum` or `ploc_fallback`).

### Substep 2.3 — Call decision

The four channel peaks are ranked by height (ties broken by channel order
A < C < G < T). Let `top` be the highest height. If `top <= 0` or the second
peak ties `top`, the call is unresolved (`N`). Otherwise, let
`primary_peak_position` be the selected peak position of the uniquely strongest
channel:

- **Qualifying channels** have a positive selected peak satisfying
  `selected_height / top >= secondary_peak_ratio` and positive channel signal at
  `primary_peak_position` satisfying
  `signal_at_primary_peak / top >= secondary_peak_ratio`.
- The two conditions are an intersection with the v2 selected-peak rule, so a
  remote maximum elsewhere in the same call window cannot create ambiguity and
  no new secondary channel can qualify.
- **Primary** is the strongest base when one to three channels qualify; four qualifying channels produce `N`.
- **Ambiguity** depends on the number of qualifying channels: one base maps to
  itself (canonical); two bases map to the standard two-base IUPAC symbol; three
  bases are unresolved `N` (primary is still the strongest); four bases are
  unresolved `N` for both primary and ambiguity.
The primary and ambiguity sequences are the concatenation of the per-call
primary and ambiguity symbols.

## Stage 3 — DNA processing (`dna.windowed_snr/v1`)

Calculates observation-only signal-quality features from the immutable analyzed channels and basecalling evidence. It uses full-width, stride-one windows of configured size `5..=10` calls. Each base call retains the sample interval used for peak selection, so a rolling call interval maps to one exact channel-sample span.

For each channel, the local baseline is the median sample and noise sigma is the median absolute deviation of first differences divided by `0.67448975 × sqrt(2)`, with a one-channel-unit floor. Selected peak heights are baseline-corrected and divided by channel noise. Every internal window records its minimum primary SNR, maximum secondary SNR, and whether the minimum is strictly below `minimum_primary_snr`. Values are rounded to six decimal places before comparison; compact analysis v7 serializes only each merged region's minimum primary SNR.

DNA processing also derives one authoritative internal `LocusEvidence` record per PLOC-defined locus. The configured rolling-window width selects a deterministic nearby sample context for baseline/noise estimation. Inside the shared locus window, event refinement selects the sample with maximum total non-negative baseline-corrected A/C/G/T signal, breaking ties by nearest PLOC and then lower sample coordinate. At that event, DNA retains raw A/C/G/T values, corrected amplitudes, SNR, and an optional normalized `EvidenceProfile`. The profile is absent when corrected signal mass is zero and is independent of primary/ambiguity calls, selected basecall peaks, qualifying-channel membership, and `secondary_peak_ratio`.

Overlapping or adjacent candidate-noisy windows are unioned into 0-based half-open call and sample intervals only when a consecutive run contains at least `minimum_noisy_windows` windows (default 2). Isolated candidate windows do not form a noisy interval. Clean gaps are never filled. Windows and per-locus evidence remain internal; compact analysis v7 emits only merged regions. These annotations do not alter calls, candidate-noisy classification, quality, trimming, alignment, warning totals, variant eligibility, or the compact JSON contracts. See [`signal-processing.md`](signal-processing.md) for formulas, evidence, and limitations.

## Stage 4 — Quality control (`dna.apollo_relative_quality/v1`,
`dna.apollo_end_trim/v1`)

Computes one bounded, uncalibrated quality value per call and selects one
retained interval. It never removes internal sequence regions.

### Substep 4.1 — Per-call penalty

For each call `i`, a window of `trim_window_size` calls centered on `i` is
examined. The penalty is the sum of two components:

- **Ambiguity penalty:** the count of calls in the window whose ambiguity symbol
  is not a canonical A/C/G/T.
- **Spacing penalty:** with `mean_spacing` the average distance between adjacent
  basecall positions across the whole read, and `min`/`max` the minimum and
  maximum adjacent spacing inside the window, the spacing penalty is
  `floor((|max - mean| + |min - mean|) / 2)`.

The penalty is `ambiguity + spacing_penalty`.

### Substep 4.2 — Best section

The best contiguous section is the window of length
`max(1, floor(call_count * best_section_fraction))` with the minimum summed
penalty. Its average penalty is recorded.

### Substep 4.3 — Relative quality score

Scores are uncalibrated and bounded. Let `max_penalty` be the largest penalty in
the read. If `max_penalty <= 0`, every call receives
`max_relative_quality_score`. Otherwise each call receives

```text
floor(max_relative_quality_score * (1 - penalty / max_penalty))
```

clamped to `[0, max_relative_quality_score]`. These scores are **not** Phred
calibrated; `phred_calibrated` is always `false`.

### Substep 4.4 — End trimming

The trim threshold is `trim_stringency * best_average * trim_window_size`.
Starting from the best section, the algorithm walks outward and stops when a
window's summed penalty exceeds the threshold, producing `trim_start` and
`trim_end`. The retained interval must contain at least
`minimum_retained_bases` calls, otherwise analysis fails. The retained sequence
is `primary_sequence[trim_start..trim_end]`.

### Substep 4.5 — Per-call record

Each call records its penalty, relative quality score, and optional vendor
quality. `vendor_quality_applies` is true only when a vendor quality exists and
the vendor primary agrees with the signal primary. Retention is represented once
by the global trim interval rather than duplicated per call.

## Stage 5 — Alignment (`dna.profile_gotoh/v1`)

Aligns the retained basecall-independent evidence-profile sequence to the reference with affine-gap Gotoh dynamic programming. The retained primary sequence stays attached to traceback coordinates for admission metrics and downstream primary-sequence variant extraction. Alignment is semi-global: the retained query is fully consumed while unaligned reference flanks are allowed.

### Substep 5.1 — Orientation candidates

The retained query evidence is aligned in both orientations:

- **forward:** retained profile order and retained primary sequence as-is;
- **reverse:** reverse profile order with A↔T/C↔G profile complementation, plus the reverse-complemented retained primary sequence for traceback character/provenance mapping.

For a circular reference, the reference is duplicated (concatenated with itself)
so the query may wrap across the origin; the working reference length is the
modulo length. A traceback may consume at most one reference length, so a query
whose required reference span is longer than the circle is unsupported.

### Substep 5.2 — Gotoh scoring

Three dynamic-programming matrices track match, insertion, and deletion states. All score deltas use fixed-point scale 1024. For canonical reference base `r`, profile support is quantized as `u = round(weight[r] × 1024)` and substitution score is `u × match_score + (1024-u) × mismatch_score`. A missing profile or non-canonical reference base receives `1024 × ambiguous_score`. Gap open/extension deltas use the same scale, preserving `open + k × extension` semantics and preserving the clean one-hot method ordering. Endpoint candidates are ranked from the last query row, allowing free reference flanks. For a circular reference, DNA selects the highest-scoring traceback whose consumed reference span is at most one circle rather than letting an invalid unbounded candidate mask a valid placement. Allocation is bounded by a compiled cell cap.

### Substep 5.3 — Traceback

The traceback internally decodes the selected path into equal-length gapped query and
gapped reference strings, an operation-run string (e.g. `5M`, `3M1I1M`), and
alignment metrics. Compact analysis v7 emits only the selected alignment summary and
reference segments, not the rows, operation runs, or score. When multiple paths tie, a documented state order
(match > deletion > insertion) makes the result deterministic. Metrics are:

- `exact_matches`, `mismatches`, `gap_opens`;
- `callable_columns` (columns where both bases are canonical);
- `callable_identity` = `exact_matches / callable_columns` (0 when no callable
  columns);
- `unresolved_query_bases` (query `N` columns).

### Substep 5.4 — Orientation selection

The forward and reverse candidates are compared by fixed-point profile score only. The strictly better orientation is selected; an exact score tie is an error and primary-sequence exact/mismatch metrics do not break it. After placement, `callable_columns`, `callable_identity`, exact/mismatch counts, and unresolved-query count are still computed from the retained primary sequence on the selected traceback. The selected orientation must meet the existing `minimum_callable_bases` and `minimum_identity` primary-sequence admission gates, otherwise analysis fails.

### Substep 5.5 — Reference segments

For a linear reference, the alignment maps to one half-open reference segment.
For a circular reference, the aligned span is projected back onto the reference;
if it crosses the origin it is split into two segments and `wraps_origin` is
`true`.

## Stage 6 — Variant calling (`dna.primary_difference/v4`)

Extracts normalized primary-sequence differences from the selected alignment.
Only differences in the primary sequence are considered; no allele-frequency,
genotype, or heteroplasmy inference is performed.

### Substep 6.1 — Difference extraction

Walking the alignment columns:

- a column where the query is `-` is a **deletion** of the reference bases;
- a column where the reference is `-` is an **insertion** of the query bases;
- a column with unequal canonical query and reference bases is an **SNV**.

Differences whose allele contains a non-canonical base, or whose indel length
exceeds `max_indel_length`, increment the excluded-candidate warning count rather
than being reported.

### Substep 6.2 — Normalization

Reported variants are normalized:

- **linear references:** indels are left-normalized against the reference where
  an equivalent placement exists (`linear_left`). When no aligned left flank
  exists, the actual reference predecessor is derived from the event position; a
  true linear origin insertion/deletion right-anchors to the next reference base.
- **circular references:** indels are placed at the canonical rotation
  (`circular_canonical`). Repeat normalization walks the whole circle, so the
  resulting representation is anchor-independent.

Internally each variant retains its contig, 1-based position, reference/alternate
alleles, kind, normalization, and direct call mappings. Compact v7 emits only
`position`, `reference`, `alternate`, `kind`, and `calls`. Every public call
contains only its supporting/flanking `role`, reference-oriented called `base`,
co-located reference-oriented A/C/G/T primary-event channel heights in `peaks`,
and uncalibrated `quality`. Original call index, PLOC, mapped call position,
trace-strand symbols, selected-peak positions/sources, penalties, and vendor
evidence remain internal. Deletions carry real aligned flanks only and never
fabricate deleted-base dna. The emitted reference allele is validated against
the supplied reference. Normalization may move the allele representation without
changing the underlying observed evidence.

### Substep 6.3 — Configured eligibility

A normalized candidate is retained only when its 1-based anchor `position` lies
inside at least one configured inclusive region. SNV supporting calls and every
inserted-base supporting call must each have a highest A/C/G/T peak greater than
or equal to `minimum_peak_height` and an uncalibrated relative score strictly
greater than `relative_quality_threshold`. Insertion flanks are not evaluated.
Deletions have no supporting trace base, so their flanks are not subjected to
peak or quality thresholds; their normalized anchor must still be in a region.
Vendor PCON is not used by this filter. For SNVs, a supporting call with more than one co-localized qualifying channel is retained as a normalized observation but is ineligible for clean-SNV reporting with `mixed_supporting_dna`. Insertions and deletions are not subjected to this point-mixed-signal gate; persistent mixed-length evidence is a separate method boundary.

Each removed candidate increments `excluded_variant_candidates` once, even when
it fails more than one eligibility condition. The pure variant stage also returns
a concise exclusion diagnostic containing kind, contig, normalized position when
available, and all failed rules. Pipeline orchestration writes one WARN record per
diagnostic without reference/alternate alleles. Sample aggregation logs aggregate
counts of differential-locus observations and variant-associated calls that
retain a basecall-independent profile; those operational counts do not alter the
scientific result.

### Substep 6.4 — Ordering

Reported variants are sorted by `(contig, position, reference, alternate)` and
deduplicated.

## One-read observation boundary

After selected alignment and variant calling, DNA materializes a `ReadObservation` that owns the input identity, base calls, basecall-independent locus/signal observations, quality-control result, selected alignment, and read-level variant result for exactly one trace.

The read has already located itself at this boundary. Its orientation and covered reference segments come from evidence-driven semi-global alignment and circular projection; filenames or nominal HV/F/R labels are not placement inputs. This same one-read product feeds both the current analysis report and implemented sample-level reconciliation.

## Sample evidence aggregation

`signal sample` processes every trace through the one-read observation path before aggregation. `sample::aggregate` requires identical reference/configuration identities, rejects duplicate input SHA-256 values, and sorts reads by SHA-256 independently of CLI trace order. The top-level read registry retains reviewer-facing filename stem, stable SHA-256, and the concise selected-alignment summary (orientation, callable bases/identity, gap opens, unresolved bases, mapped segments, and origin-wrap state).

Before pairwise/locus aggregation, DNA derives a run-length reference coverage topology from every selected post-trim read segment. Each maximal interval records total read depth plus forward/reverse orientation depth. This counts all independently placed reads regardless of later pairwise eligibility and does not imply nucleotide agreement or consensus admission.

DNA then builds a deterministic pairwise overlap graph from the SHA-sorted read registry. Every unordered pair is compared only at shared selected-alignment reference coordinates. `shared_positions` counts all shared coordinates, while the agreement denominator includes only positions where both query observations are canonical A/C/G/T. Equal canonical observations are agreements; unequal canonical observations are conflicts. Unresolved symbols and deletions do not enter that nucleotide denominator, so gap/indel evidence remains separate. A pair is eligible for later consensus reconciliation only when the comparable-base count reaches `sample_reconciliation.minimum_comparable_bases` and the agreement fraction reaches `sample_reconciliation.minimum_overlap_agreement`. Non-overlapping reads produce no pair edge and remain valid sample evidence.

Before sparse locus/variant projection, sample reconciliation resolves every call-backed sample observation back to the authoritative `LocusEvidence` record by original call index. One internal `CallDNAEvidence` preserves baseline-corrected A/C/G/T amplitudes, per-channel SNR, the optional basecall-independent `EvidenceProfile`, and existing merged candidate-noisy-region membership. All channel-valued evidence is projected into reference orientation: forward reads keep A/C/G/T order and reverse reads reorder T/G/C/A into reference-oriented A/C/G/T, including profile complementation. A zero-signal locus remains profile-less and deletion observations have no call signal object. This evidence does not alter placement, overlap admission, callability, variant eligibility, or consensus weighting. The v8 sample report projects only the normalized reference-oriented profile and existing noisy-region membership at retained differential loci; corrected amplitudes and per-channel SNR remain internal.

Sparse locus aggregation then runs in two passes. The first pass identifies reference positions where at least one covering read is alternate, unresolved, or deleted. The second pass retains every covering read only at those positions, including canonical reference support with observed base and quality. Each retained locus also derives factorized support topology over total reads, forward/reverse orientation, reference/alternate/unresolved/deletion states, and the subset of observations with an available basecall-independent profile split by forward/reverse orientation. v8 projects these counts directly into each sparse differential-locus record. Profile availability is only a pre-consensus denominator: it does not imply agreement, contributor eligibility, or weight. Each retained locus observation separately receives structural nucleotide-contribution eligibility: call-backed observations with a real profile are `Eligible`, missing-profile calls are `MissingProfile`, and deletions are `DeletionEvent`. Unresolved symbols and candidate-noisy membership do not suppress real profile evidence in this first policy. Eligible profiles are then accumulated with unit read mass into total and forward/reverse A/C/G/T support vectors. Each non-empty partition also derives its arithmetic mean A/C/G/T evidence profile from support divided by contributor count. DNA then decomposes profile heterogeneity into mean within-read impurity and between-read dispersion, retaining the same geometry for total/forward/reverse partitions. When both orientation means exist, their Total Variation distance is retained as `directional_profile_distance`. These remain descriptive evidence geometry only: amplitude, SNR, quality, and noisy-region state do not alter the unit weight, and no consensus base, biological-mixture label, discordance class, or confidence is selected. Call-backed locus observations additionally preserve whether their source call falls inside an existing merged candidate-noisy region; deletions preserve that call-specific state as absent. Both partitions must sum to the same total read count. The topology is evidence summary only: it introduces no weighting, confidence, consensus, or biological-independence claim. Positions inside a read's mapped segments but absent from `locus_differences[]` are therefore canonical reference matches; positions outside the mapped segments are uncovered. Routine all-reference loci are never materialized in sample evidence.

Canonical normalized variant observations are separately grouped by `(position, reference, alternate, kind)`. Each support publishes the unique reviewer-facing read name plus configured eligibility, exclusion reasons, and reference-oriented base/peak/quality evidence. For every normalized variant, DNA also derives factorized support topology from those same observations: total/eligible observing reads and forward/reverse plus eligible-forward/eligible-reverse counts. The summary never adds reference-supporting or unresolved coverage as variant support and is not a confidence or vote. Internal aggregation remains SHA-ordered and index-based, but numeric indexes do not leak into the reviewer contract. A read-level filter can remove a candidate from `analysis/v7` reporting without erasing the observation from `SampleEvidence`. Insertions are normalized variant evidence rather than fabricated reference-locus observations. No consensus or sample-level conflict verdict is produced in v8; coverage topology and overlap eligibility are pre-consensus evidence only.

## Output

`analyze` publishes `dna.analysis/v7` at `results/<trace-stem>.json`.
`sample` publishes `dna.sample_evidence/v8` at
`results/<sample-id>.sample.json`; its detailed semantics are defined in
[`sample-output.md`](sample-output.md). Both use the same atomic no-overwrite
publisher and keep operational logs outside deterministic JSON.

The completed `dna.analysis/v7` result contains compact provenance, read count
and trim bounds, merged candidate-noisy regions, the selected post-trim alignment
summary, normalized variants, and reviewer-facing reference-oriented call evidence
(`base`, four co-located A/C/G/T `peaks`, and `quality`). It omits filenames,
full sequences, individual rolling windows, gapped rows, operation runs, method
constants, call indexes/PLOC coordinates, selected-peak positions/sources, vendor
data, and redundant fields. The strict configuration remains schema version 5. No compatibility result is emitted. The document is
published atomically to `results/<trace-stem>.json` without overwriting. Operational records are appended
separately to `$DNA_LOG_DIR/<trace-stem>.log` (default `logs/`) and are not part
of deterministic JSON. One run-correlated record summarizes input/decode,
basecalling, signal processing, quality control, alignment, variant calling, and publication
readiness with aggregate metrics and elapsed milliseconds. WARN records identify
removed candidates by kind/contig/position/reasons and summarize final warning
categories; ERROR records identify the active failed stage. Records omit complete
sequences, alleles, region contents, per-call peaks, alignment strings, and JSON
bodies. The JSON shape is defined in
[`json-output.md`](json-output.md), validated by
[`schemas/analysis-v7.schema.json`](schemas/analysis-v7.schema.json), and shown in
[`examples/analysis-v7.example.json`](examples/analysis-v7.example.json).

## Biological limitations

DNA is a primary-sequence analysis tool, not a clinical or population
diagnostic. The following limitations are intentional and documented:

- **No genotype or heteroplasmy calls.** The pipeline reports differences in a
  single primary sequence. It does not estimate allele fractions, genotype
  likelihoods, or heteroplasmy levels, and it does not decompose mixed or
  two-allele signals.
- **Observational SNR.** Rolling SNR values are robust local features, not Phred scores or error probabilities. Candidate-noisy intervals do not suppress calls or variants.
- **Uncalibrated quality.** Quality values are relative, bounded scores derived
  from ambiguity and peak spacing. They are not Phred-calibrated and are not
  error probabilities.
- **Single reference, single orientation.** The query is aligned to one
  reference record in one of two orientations. Multi-contig references,
  alternative references, and reference search/indexing are out of scope.
- **Primary-sequence variants only.** Differences are derived from the conservative
  signal-derived primary sequence. A mixed two- or three-channel call may still
  produce a normalized strongest-base SNV observation, but it is retained only as
  evidence and is not eligible for ordinary clean-SNV reporting. Unresolved N
  differences, indels longer than `max_indel_length`, out-of-region candidates,
  mixed-supporting SNVs, and SNV/insertion candidates below configured supporting-
  signal thresholds are excluded from the single-read report. This does not infer
  genotype or heteroplasmy.
- **Sanger trace limitations.** Basecalling depends on the quality of the
  four-channel signal and the vendor-defined basecall positions. Poor signal,
  mixed templates, and sequencing artifacts can produce unresolved (`N`) calls
  or excluded candidates rather than confident differences.
- **No clinical interpretation.** The output is a technical analysis document.
  It does not assign pathogenicity, disease association, or clinical
  significance to any observed difference.
