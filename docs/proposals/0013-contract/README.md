# PROP-0013 draft contract

Draft, non-normative files of [PROP-0013](../0013-review-contract-v2.md). They move to
`docs/reference/schemas/` and `docs/reference/examples/` with increment 1, and stay here until then.

- [review-v2.schema.json](review-v2.schema.json) — the JSON Schema of `dna.review/v2`; it references definitions of the notation and consensus schemas by `$id`
- [review-v2.labels.json](review-v2.labels.json) — the label registry: ID, category, scopes, title, description, v1 name, and policy 1 class and rank
- [review-v2.example.json](review-v2.example.json) — the `dna.review/v1` synthetic example rewritten as v2; it validates against the draft schema
