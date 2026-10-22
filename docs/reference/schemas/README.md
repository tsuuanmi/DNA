# JSON Schemas

This directory contains the authoritative machine-readable schemas for current public result contracts.

- [analysis-v10.schema.json](analysis-v10.schema.json) — `dna.analysis/v10`; it references definitions of the read evidence schema by `$id`
- [basecalls-v3.schema.json](basecalls-v3.schema.json) — `dna.basecalls/v3`
- [read-evidence-v1.schema.json](read-evidence-v1.schema.json) — `dna.read_evidence/v1`; it references definitions of the basecall schema by `$id`
- [sample-evidence-v10.schema.json](sample-evidence-v10.schema.json) — `dna.sample_evidence/v10`
- [variants-v1.schema.json](variants-v1.schema.json) — `dna.variants/v1`
- [notation-v1.schema.json](notation-v1.schema.json) — `dna.notation/v1`
- [consensus-v1.schema.json](consensus-v1.schema.json) — `dna.consensus/v1`
- [review-v1.schema.json](review-v1.schema.json) — `dna.review/v1`; it references definitions of the notation and consensus schemas by `$id`
- [signal-v1.schema.json](signal-v1.schema.json) — `dna.signal/v1`
- [run-context-v1.schema.json](run-context-v1.schema.json) — `dna.run_context/v1`; it references definitions of the consensus schema by `$id`, and the consensus and review schemas reference its `use`
- [eval-run-v1.schema.json](eval-run-v1.schema.json) — `dna.eval_run/v1`, the output directory of the evaluation harness `dna-eval run`: validation evidence, not a DNA result; its semantics are in the [evaluation harness](../../validation/evaluation-harness.md#outputs)
- [eval-regress-v1.schema.json](eval-regress-v1.schema.json) — `dna.eval_regress/v1`, the output directory of `dna-eval regress check`: validation evidence, not a DNA result; its semantics are in the [regression set](../../validation/regression-set.md#outputs)

Schemas define exact serialized shape. Human-readable semantics remain in the parent [contracts index](../README.md) and the corresponding contract folders.
