# DNA Software Requirements Specification

This document is the canonical index for DNA's normative Software Requirements Specification (SRS).

DNA is a deterministic Rust library and CLI for DNA analysis. The current normative production requirements cover the implemented Sanger ABIF path, reviewed consensus sequences (`call`), sample consensus (`consensus`), and variants documents (`notation`) against short references; future mtDNA, nuclear/genomic DNA, targeted/SNP, NGS, or downstream capabilities become normative only when their requirement families are added. Normative terms **MUST**, **SHOULD**, and **MAY** apply to every `SRS-*` item.

## Purpose

DNA exists to shrink the manual review of sequencing results without
lowering their accuracy. No system is perfect, so DNA sorts the work instead
of hiding it.

The sorting starts from a catalogue of the errors that occur in the data.
Each error class has one outcome:

- **pass** — DNA resolves the error from the evidence, or the error is
  absent, and nobody needs to review the result;
- **review** — DNA cannot resolve the error by itself, but a reviewer can
  settle the result from the traces; DNA says where and why;
- **rerun** — the trace is too poor or unreliable to settle anything; DNA
  names the trace and the cause, so the lab can rerun it quickly.

DNA makes these decisions at four levels, worked on from the bottom up:

- **region** — a stretch of the target;
- **trace** — one read;
- **sample** — the whole result;
- **batch** — a run or plate, aggregated by downstream pipelines (ADR-0065).

The levels are not a tree. A trace covers several regions, and several traces
cover a region. A sample passes only when all its regions pass, and a region
normally needs both strands to pass.

Success is measured on the share of work in each outcome. If the input is
100 %, A % passes, B % goes to review, and C % goes to rerun. The lab's
chemistry sets C; DNA decides it quickly and names the cause. DNA aims to:

- raise A as far as the evidence allows;
- lower B step by step, class by class;
- keep passed results free of errors. A wrong pass costs more than a review.
  The error rate of passed results is measured against ground truth.

Today DNA reports the evidence for these decisions (callability, rejected
reads, eligibility reasons, consensus `N` and `contested` marks), but not the
decisions themselves. [PROP-0004](../proposals/0004-review-triage.md) proposes
how DNA could make them.

## Requirement families

| Namespace | Scope | Specification |
|---|---|---|
| `SRS-IN-*` | inputs and process boundary | [input.md](input.md) |
| `SRS-API-*` | public Rust library boundary | [api.md](api.md) |
| `SRS-CFG-*` | configuration | [configuration.md](configuration.md) |
| `SRS-PRF-*` | target profiles | [profiles.md](profiles.md) |
| `SRS-BC-*` | DNA-derived base re-calling | [basecalling.md](basecalling.md) |
| `SRS-SIG-*` | observational signal processing | [signal-processing.md](signal-processing.md) |
| `SRS-CALL-*` | signal-derived read callability: phase state, mask, callable span | [callability.md](callability.md) |
| `SRS-QC-*` | quality control and trimming | [quality-control.md](quality-control.md) |
| `SRS-ALN-*` | reference alignment and placement | [alignment.md](alignment.md) |
| `SRS-SAMPLE-*` | multi-read sample evidence | [sample-evidence](sample-evidence/README.md) |
| `SRS-CONS-*` | sample consensus | [consensus.md](consensus.md) |
| `SRS-VAR-*` | primary-sequence differences / calling | [variants.md](variants.md) |
| `SRS-VN-*` | optional post-calling variant normalization | [variant-normalization.md](variant-normalization.md) |
| `SRS-NOM-*` | optional target-specific variant nomenclature | [variant-nomenclature.md](variant-nomenclature.md) |
| `SRS-OUT-*` | result publication and operational logging | [output.md](output.md) |
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
