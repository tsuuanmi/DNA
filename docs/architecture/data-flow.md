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
modality-specific evidence until biological differences genuinely converge.
Variant calling, canonicalization, and nomenclature are separate boundaries.

```text
Sanger evidence --> Sanger analysis --> Sanger variant caller ----+
                                                                  |
NGS evidence ----> NGS analysis ----> NGS variant caller ---------+--> called variants
                                                                  |         |
other validated source --> source-specific adapter ---------------+         v
                                                                  canonicalization policy
                                                                            |
                                                                            v
                                                                    target nomenclature
                                                                            |
                                                                            v
                                                                     canonical variants
                                                                      /      |       \
                                                                     v       v        v
                                                                sample   haplogroup  targeted/
                                                              reconcile               SNP analysis
```

For current human mtDNA, the planned canonicalization policy includes
haplotype-preserving sequence-equivalent 3'/right-most indel placement. The
right-shifting algorithm can be reusable machinery, but selecting right-most
placement is not a universal policy for every future biological target.

Sanger alignment still has its own deterministic right-most traceback
canonicalization for evidence placement and provenance. Post-calling variant
canonicalization is a separate contract so future NGS or imported variant
sources do not need to pass through the Sanger aligner.

This is architectural direction, not a claim that NGS, post-calling
canonicalization, or target nomenclature are currently implemented.

Coordinate/strand semantics are canonical in
[reference/coordinates](../reference/coordinates.md); pipeline invariants are
canonical in [invariants](invariants/README.md).
