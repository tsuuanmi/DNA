# Data and Fixture Governance

This document governs local biological data, derived results, validation truth, and approved real-trace evidence. Operational batch mechanics live in [operations/batch.md](../operations/batch.md).

## Scope

The repository workspace may contain a local AB1 corpus under `data/` for scientific validation. Local data is development evidence, not part of the DNA product interface and not a prerequisite for normal builds or CI.

The contents and filenames of local corpora may change independently of the repository. Documentation, tests, and CI MUST NOT depend on a particular local file count, sample identifier, or filename.

## Privacy and repository policy

Treat raw traces, manifests, filenames, sample identifiers, reviewer-derived truth, hashes, and biological differences as potentially identifying data.

- `data/`, `results/`, `logs/`, and `validation-results/` remain ignored.
- Do not force-add AB1 files, manifests, sample identifiers, reviewer truth, or derived outputs.
- Do not copy local traces into committed test fixtures without explicit approval.
- Repository examples MUST use synthetic names, identities, hashes, and scientific values.
- CI and normal unit tests MUST work when local biological data is absent.

## Derived result handling

Compact analysis results omit many raw trace details but still contain input/reference/configuration identities and biological differences. Sample evidence additionally exposes the sample identifier and read names. Reference-free basecall results contain complete sequence strings.

All of these artifacts inherit the approval, storage, retention, and redistribution restrictions of their source AB1 data.

Operational logs intentionally omit complete sequences, peak arrays, alleles, and JSON bodies, but can still contain trace/reference names, paths, hashes, coordinates, aggregate metrics, thresholds, and failures. Treat them under the same local-data policy.

## Reviewer-derived truth

Reviewer comparison tables and extracted variant truth are validation data, not ordinary repository fixtures. Keep them in approved ignored storage such as:

```text
data/validation/ground-truth/
```

Generated evaluation artifacts under `validation-results/` inherit the same policy because extra/missing/representation rows disclose biological differences.

Repository tests may exercise the same parsing and comparison logic only with synthetic reviewer notation and synthetic identities.

## Approval record for real validation evidence

Before using a real AB1 trace as compatibility, validation, or release evidence, retain an approval record containing:

1. approval and intended use;
2. source run and instrument context;
3. primer/region and expected orientation;
4. AB1 SHA-256 checksum;
5. reference identity and checksum;
6. DNA configuration checksum or exact values;
7. expected stage outputs or biological truth and how they were established;
8. redistribution permission or restriction.

A trace without this record may be used for exploratory local debugging, but not to support compatibility or release claims.

## Operational boundary

Selecting, grouping, cleaning, and rerunning a local corpus is external orchestration. It does not expand the core CLI input contract. See [local batch orchestration](../operations/batch.md).
