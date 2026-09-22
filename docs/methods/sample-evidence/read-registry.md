# Sample Read Registry

Part of the canonical [sample evidence aggregation method](../sample-evidence.md).

Part of the canonical [DNA pipeline](pipeline.md). Normative requirements live in [SRS-SAMPLE-*](../srs/sample-evidence.md), and public serialization semantics live in the [sample evidence contract](../contracts/sample-evidence.md).

`signal sample` processes every trace through the one-read observation path before aggregation. `sample::aggregate` requires identical reference/configuration identities, rejects duplicate input SHA-256 values, and sorts reads by SHA-256 independently of CLI trace order. The top-level read registry retains reviewer-facing filename stem, stable SHA-256, and the concise selected-alignment summary (orientation, callable bases/identity, gap opens, unresolved bases, mapped segments, and origin-wrap state).
