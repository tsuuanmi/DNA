# Data Flow

## Current Sanger production path

```text
Sanger ABIF
    |
    v
dna_sanger::abif
    |
    v
model::sanger::Chromatogram
    |
    v
read_processing
    |
    +--> basecalling
    +--> signal_processing
    +--> callability
    +--> quality_control
              |
              v
      variant_analysis::observation
              |
              v
   read_processing::evidence (Sanger adapter)
              |
              v
         ReadEvidence
              |
FASTA ------> alignment
              |
              v
        variant_calling <-- ReadEvidence
              |
              v
        ReadObservation = CalledRead + SangerAttachment
           /                 \
          v                   v
   analysis result     sample aggregation (CalledRead only)
                              |
                              v
                       sample evidence --> report joins SangerAttachment
```

A sample's consensus document is called from its decided edits, without
placing a segment again (ADR-0078):

```text
dna.consensus/v1 --> input::consensus --> dna_consensus::decided --> CalledRead --> called variants
```

Reviewed consensus sequences take a shorter path into the same core:

```text
consensus FASTA --> input::sequence --> ReadEvidence (vouched ends)
                                             |
FASTA ---------------------------------> dna_call::read_call
                                             |
                                         CalledRead --> called variants
```

A variants document reaches the post-calling plugins in another process:

```text
variants document --> input::variants --> normalization --> nomenclature
                                                                 |
                                                   conformance --+--> notation document
```

`variant_analysis::observation` owns the authoritative reference-guided
single-read scientific path. The public Rust Variant Analysis capability, CLI
analysis, and sample processing reuse that path rather than owning duplicate
implementations. Sample processing is one pass, `sample_analysis::reads`, that
places each trace once; the sample evidence and the consensus both start from
it, in the `sample` and `consensus` commands and in the public Sample Analysis
capability alike.

Alignment and variant calling see the read only as `ReadEvidence`
([ADR-0069](../decisions/adr/0069-plugin-first-modality-core-post-calling.md)):
the Sanger adapter turns calls, signal, callability, and quality into bases,
profiles, labelled masks, support vetoes, and the informative interval.

Completed typed scientific state is projected by the report layer and serialized
before atomic no-overwrite CLI publication. Logging remains operational side
evidence rather than part of deterministic scientific result contracts.

## Platform direction

Future modalities keep source-specific evidence until a biological difference
has actually been called or imported. The common convergence point is
`CalledVariantSet`, not a universal raw-alignment object.

```text
Sanger evidence --> ReadEvidence --> core caller -----+
consensus FASTA --> ReadEvidence --> core caller -----+
                                                      |
NGS evidence ----> NGS caller ------------------------+--> CalledVariantSet
                                                      |       /     |      \
VCF/BCF --------> validated importer -----------------+      v      v       v
                                                         direct  optional  optional
                                                          use   normalize  nomenclature
                                                                   |
                                                                   +--> optional nomenclature
```

The source-specific paths are intentionally different:

| Source | Path to called variants |
|---|---|
| Sanger ABIF | chromatogram -> base calling/signal/callability/QC -> `ReadEvidence` -> core caller (implemented) |
| consensus document (`dna.consensus/v1`) | the consensus's decided edits (implemented, `call`, `SampleAnalysis::call`) |
| assembled/consensus FASTA | `ReadEvidence` -> core caller (implemented, `call`) |
| FASTQ / NGS reads | read QC/preprocessing -> mapping -> NGS caller |
| BAM / CRAM | validated aligned-read evidence -> NGS caller |
| VCF / BCF | validated variant importer; raw calling is bypassed |

These rows are architecture direction, not claims that all sources are currently
implemented.

For current human mtDNA, one optional normalization policy is
haplotype-preserving 3'/right-most indel placement. Nomenclature is a separate
optional target-specific representation layer. A workflow may consume called
variants directly when it does not need either representation policy.

Sanger alignment still has its own deterministic right-most traceback
canonicalization for evidence placement and provenance. That alignment rule is
separate from any post-calling normalization selected by another workflow.

The detailed boundary is owned by
[variant lifecycle](variant-lifecycle.md) and ADR-0060.

Coordinate/strand semantics are canonical in
[reference/coordinates](../reference/coordinates.md); pipeline invariants are
canonical in [invariants](invariants/README.md).
