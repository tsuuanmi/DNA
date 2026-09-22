# Data Flow

```text
AB1 -> decode -> basecalling -> signal_processing -> quality_control
                                                   |
                                                   +-> basecall result
                                                   |
FASTA -----------------------------------------> alignment -> variant_calling -> ReadObservation
                                                                                   |
                                                            +----------------------+
                                                            |                      |
                                                            v                      v
                                                    analysis result          sample aggregation
                                                                                   |
                                                                                   v
                                                                         sample evidence
```

`pipeline::observation` is the authoritative reference-guided read path. Single-read analysis consumes one observation; sample processing independently creates one observation per trace before aggregation.

Completed typed scientific state is projected by the report layer and serialized before atomic no-overwrite publication. Logging remains operational side evidence rather than part of deterministic scientific result contracts.

Coordinate/strand semantics are canonical in [reference/coordinates](../reference/coordinates.md); pipeline invariants are canonical in [invariants](invariants/README.md).
