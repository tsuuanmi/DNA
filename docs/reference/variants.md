# DNA Variants JSON

`dna call <sample-id> <sequences.fasta>... --reference <reference.fasta>`
runs the core caller alone over reviewed consensus sequences
([ADR-0069](../decisions/adr/0069-plugin-first-modality-core-post-calling.md))
and writes one deterministic `results/<sample-id>.variants.json` document
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
  `nomenclature` when the profile declares notation).
- `reads[]`: one entry per FASTA record, in input order:
  - `name`: the record identifier, unique within the call;
  - `sha256`: the SHA-256 of the record's upper-case sequence;
  - `alignment`: the selected-alignment summary of the
    [sample read registry](sample-evidence/reads-coverage.md). A consensus
    masks no call, so `masked_bases` is `0`, and the callable reference
    segments equal the mapped ones.
  - `variants[]`: every observed variant of the read, sorted by `position`,
    `reference`, and `alternate`. Each has the anchored `reference` and
    `alternate` alleles, its `kind` (`SNV`, `INS`, or `DEL`), `eligible`, and
    its `exclusion_reasons`. These are core reasons (`outside_target_region`,
    `indel_length_exceeded`, `non_canonical_allele`, `read_end`,
    `run_boundary`), or modality labels reported verbatim.
- `notation` (only when the profile declares notation): the same per-read view
  as in [sample evidence](sample-evidence/notation.md). Each rendered call
  lists the reads whose represented eligible calls contain it.

Reads are not reconciled with each other. Reads that disagree contribute
different calls, and DNA makes no sample-level verdict.
