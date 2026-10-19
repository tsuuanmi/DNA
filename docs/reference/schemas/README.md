# JSON Schemas

This directory contains the authoritative machine-readable schemas for current public result contracts.

- [analysis-v9.schema.json](analysis-v9.schema.json) — `dna.analysis/v9`
- [basecalls-v3.schema.json](basecalls-v3.schema.json) — `dna.basecalls/v3`
- [read-evidence-v1.schema.json](read-evidence-v1.schema.json) — `dna.read_evidence/v1`; it references definitions of the basecall schema by `$id`
- [sample-evidence-v10.schema.json](sample-evidence-v10.schema.json) — `dna.sample_evidence/v10`
- [variants-v1.schema.json](variants-v1.schema.json) — `dna.variants/v1`
- [notation-v1.schema.json](notation-v1.schema.json) — `dna.notation/v1`
- [consensus-v1.schema.json](consensus-v1.schema.json) — `dna.consensus/v1`
- [review-v1.schema.json](review-v1.schema.json) — `dna.review/v1`; it references definitions of the notation and consensus schemas by `$id`
- [signal-v1.schema.json](signal-v1.schema.json) — `dna.signal/v1`

Schemas define exact serialized shape. Human-readable semantics remain in the parent [contracts index](../README.md) and the corresponding contract folders.
