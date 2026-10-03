# Data Flow

## Current Sanger production path

```text
Sanger ABIF
    |
    v
input::sanger::abif
    |
    v
model::sanger::Chromatogram
    |
    v
read_processing
    |
    +--> basecalling
    +--> signal_processing
    +--> quality_control
              |
              v
      variant_analysis::observation
              |
FASTA ------> alignment
              |
              v
        variant_calling
              |
              v
        ReadObservation
           /       \
          v         v
   analysis result  sample aggregation
                        |
                        v
                  sample evidence
```

`variant_analysis::observation` owns the authoritative reference-guided
single-read scientific path. The public Rust Variant Analysis capability, CLI
analysis, and sample processing reuse that path rather than owning duplicate
implementations.

Completed typed scientific state is projected by the report layer and serialized
before atomic no-overwrite CLI publication. Logging remains operational side
evidence rather than part of deterministic scientific result contracts.

## Platform direction

Future modalities keep source-specific evidence until a biological difference
has actually been called or imported. The common convergence point is
`CalledVariantSet`, not a universal raw-alignment object.

```text
Sanger evidence --> Sanger caller --------------------+
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
| Sanger ABIF | chromatogram -> base calling/signal/QC -> selected pairwise alignment -> Sanger caller |
| assembled/consensus FASTA | reference alignment -> sequence-difference caller |
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
