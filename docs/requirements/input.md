# Input and Process Boundary Requirements

**Requirement namespace:** `SRS-IN-*`

These requirements are part of the canonical [DNA SRS](SRS.md).

- **SRS-IN-001:** `analyze` and `basecall` MUST accept exactly one regular non-empty AB1 path. `sample` MUST accept one sample identifier and one or more regular non-empty AB1 paths. Reference-guided commands MUST require one regular non-empty FASTA path; `basecall` MUST NOT accept a reference. Directories, manifest/list compatibility inputs, and globs MUST NOT be accepted by the core CLI.
- **SRS-IN-002:** AB1 bytes MUST begin with `ABIF`; every directory count, size, product, offset, inline payload, and allocation MUST be checked before use.
- **SRS-IN-003:** Canonical decode MUST require `DATA.9-12`, `FWO_.1`, and `PLOC.2`; FWO MUST be an exact A/C/G/T permutation; channels MUST be equally sized; PLOC MUST be strictly increasing and in range.
- **SRS-IN-004:** `PLOC.2` is required; `PBAS.2` and `PCON.2` MAY be consumed as optional vendor evidence. `P2BA.1` MUST be ignored. Uppercase IUPAC vendor bases and the ABIF one-byte byte/char PCON representations MUST be accepted. No alternate-tag fallback is permitted.
- **SRS-IN-008:** Optional PBAS/PCON cardinality MUST NOT define or extend the scientific locus series. A present vendor series whose decoded length differs from the valid PLOC count MUST remain non-fatal trace-integrity evidence; DNA MUST process exactly the PLOC-defined loci and boundedly ignore unavailable/excess vendor positions.
- **SRS-IN-005:** A reference-guided `analyze` or `sample` FASTA MUST contain exactly one non-empty A/C/G/T/N record no longer than 50,000 bases. `basecall` MUST perform no FASTA I/O.
- **SRS-IN-006:** Empty, missing, malformed, unsupported, over-limit, or unreadable input MUST return a typed error without panic or result file.
- **SRS-IN-007:** Generic ABIF directory parsing MUST follow the container specification rather than imposing scientific-tag assumptions globally. A directory item MUST contain at least its logical `element_size × element_count` payload; a larger `data_size` that is permitted by ABIF padding/reserved-space behavior MUST NOT by itself make the file invalid. Decoders for required scientific tags MAY impose stricter tag-specific type, element-size, cardinality, and semantic constraints.
- **SRS-IN-010:** The single-read reference-guided command MUST be `signal analyze <trace.ab1> --reference <reference.fasta>` and derive `results/<trace-stem>.json`. The reference-free command MUST be `signal basecall <trace.ab1>` and derive `results/<trace-stem>.basecalls.json`. The sample command MUST be `signal sample <sample-id> <trace.ab1>... --reference <reference.fasta>` and derive `results/<sample-id>.sample.json`. No output-path or format compatibility option is permitted.
- **SRS-IN-011:** Help/version MUST succeed without reading analysis inputs.
- **SRS-IN-012:** Invalid CLI arguments MUST produce concise stderr and nonzero status.
