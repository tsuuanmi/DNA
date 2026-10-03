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

Future modalities do not need to imitate Sanger evidence. They keep
modality-specific evidence until a downstream semantic contract genuinely
converges.

```text
Sanger evidence -------------------> Sanger analysis path ---+
                                                             |
NGS read/alignment evidence -------> NGS analysis path ------+--> canonical variants
                                                             |
other validated evidence ----------> future path ------------+
                                                                    |
                                                  +-----------------+----------------+
                                                  v                 v                v
                                             haplogroup        nomenclature    targeted/SNP analysis
```

This is architectural direction, not a claim that NGS or those downstream
capabilities are currently implemented.

Coordinate/strand semantics are canonical in
[reference/coordinates](../reference/coordinates.md); pipeline invariants are
canonical in [invariants](invariants/README.md).
