# DNA Software Requirements Specification

This document is the canonical index for DNA's normative Software Requirements Specification (SRS).

DNA is a deterministic Rust CLI for reference-free base re-calling and research analysis of one Sanger ABIF/AB1 chromatogram against one short reference. Normative terms **MUST**, **SHOULD**, and **MAY** apply to every `SRS-*` item.

## Requirement families

| Namespace | Scope | Specification |
|---|---|---|
| `SRS-IN-*` | inputs and process boundary | [input.md](input.md) |
| `SRS-CFG-*` | configuration | [configuration.md](configuration.md) |
| `SRS-BC-*` | DNA-derived base re-calling | [basecalling.md](basecalling.md) |
| `SRS-SIG-*` | observational signal processing | [signal-processing.md](signal-processing.md) |
| `SRS-QC-*` | quality control and trimming | [quality-control.md](quality-control.md) |
| `SRS-ALN-*` | reference alignment and placement | [alignment.md](alignment.md) |
| `SRS-SAMPLE-*` | multi-read sample evidence | [sample-evidence](sample-evidence/README.md) |
| `SRS-VAR-*` | primary-sequence differences | [variants.md](variants.md) |
| `SRS-OUT-*` | result publication and operational logging | [output.md](output.md) |
| `SRS-BAT-*` | external batch orchestration | [batch.md](batch.md) |
| `SRS-COMPAT-*` | external-reference comparison constraints | [constraints.md](constraints.md) |
| `SRS-NFR-*` | non-functional quality attributes | [quality-attributes.md](quality-attributes.md) |
| `SRS-VAL-*` | validation acceptance criteria | [validation/acceptance criteria](../validation/acceptance-criteria.md) |

## Organization rules

- The SRS is one logical specification split by stable requirement family, not a collection of independent mini-specifications.
- Requirement IDs are stable references. Moving a requirement between files does not change its ID or semantics.
- Add a requirement to the file that owns its system responsibility; do not mirror the `src/` tree here.
- Do not create one file per requirement. Split a family only when it has a durable domain boundary and update this index when doing so.
- Cross-family requirements should reference the authoritative requirement ID and family document rather than restating the requirement.
- Detailed mechanisms belong in design documentation, rationale belongs in decision records, implementation ownership belongs in source-directory `README.md` files, and exact serialized shapes belong in reference schemas/contracts.
