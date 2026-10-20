---
implements: [SRS-CALL, SRS-CONS, SRS-VAR]
decided-by: [ADR-0069, ADR-0078]
---

# DNA Variants JSON

`dna call <sample-id> <sequences.fasta>... --reference <reference.fasta>`
runs the core caller alone over reviewed consensus sequences
([[ADR-0069]]).
`dna call <sample-id> <sample-id>.consensus.json --reference
<reference.fasta>`, given the sample's consensus document alone, calls the
edits the consensus decided without placing its segments again
([[ADR-0078]],
SRS-VAR-021). Either writes one deterministic `results/<sample-id>.variants.json` document
identified as `dna.variants/v1`. The authoritative schema is
[`schemas/variants-v1.schema.json`](schemas/variants-v1.schema.json),
and the synthetic example is
[`examples/variants-v1.example.json`](examples/variants-v1.example.json).

The document holds only what the core derives from modality-neutral evidence.
It carries no modality evidence: no peaks, quality, signal, or callability.

## Fields

- `provenance`: the reference `name`, `topology`, and sequence `sha256`;
  `configuration_sha256`; the target `profile`; and `plugins`, the workflow's
  plugins in execution order (`sequence`, then `core`, then `normalization` and
  `nomenclature` when the profile declares notation; for a consensus document,
  `decided` instead of `sequence` and `core`).
- `reads[]`: one entry per FASTA record, in input order, or per consensus
  segment, named as the consensus FASTA names it. A segment's `alignment`
  describes its own columns: the segment's interval as its one reference
  segment, forward, and the counts of its decided bases; its callable segments
  leave out every `N` and every excluded variant's positions.
  - `name`: the record identifier, unique within the call;
  - `sha256`: the SHA-256 of the record's upper-case sequence;
  - `alignment`: the selected-alignment summary of the
    [[reference/sample-evidence/reads-coverage|sample read registry]]. A consensus
    masks no call, so `masked_bases` is `0`. Its `N` and IUPAC codes count
    in `unresolved_bases` and observe no base, so the callable reference
    segments are the mapped ones less those positions, and less every
    position of an excluded variant (SRS-VAR-019): a substitution's
    position, a deletion's deleted positions, and both positions beside an
    insertion's junction. A position DNA excluded a call at is not read as
    the reference.
  - `variants[]`: every observed variant of the read, sorted by `position`,
    `reference`, and `alternate`. Each has the anchored `reference` and
    `alternate` alleles, its `kind` (`SNV`, `INS`, or `DEL`), `eligible`, and
    its `exclusion_reasons`. These are core reasons (`outside_target_region`,
    `indel_length_exceeded`, `non_canonical_allele`, `read_end`,
    `base_only`, `run_boundary`, `unresolved_flank`, `placement_conflict`),
    or modality labels reported verbatim. `base_only` excludes an insertion
    or deletion that maps a call its modality vouches for by its base only
    (SRS-CALL-014). `unresolved_flank` excludes a variant with a
    mapped call read as `N`, such as a deletion beside an unresolved base. `placement_conflict` excludes every variant of a record whose
    name declares a reference interval (`<id>_<first>-<last>`, 1-based, as
    consensus segments are named) that its placement does not overlap, and
    such a record has no callable reference segments.
- `not_analysed`: the reference intervals (0-based, half-open, as the
  segments) the call does not analyse: a reportable region's end that only
  reads' unreliable spans observe, up to the last position the consensus
  does not settle as the reference (SRS-CONS-012). Empty unless the
  document calls a consensus. `notation` lists these positions in the EMPOP
  query's `cut` as `not_analysed`.
- `notation` (only when the profile declares notation): the same per-read view
  as in [[reference/sample-evidence/notation|sample evidence]]. Each rendered call
  lists the reads whose represented eligible calls contain it.

Reads are not reconciled with each other. Reads that disagree contribute
different calls, and DNA makes no sample-level verdict.
