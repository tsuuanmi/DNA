---
decided-by: [ADR-0090]
---

# Sample Read Registry

Part of the canonical [sample evidence aggregation method](README.md).


`dna sample` processes every trace through the one-read observation path before aggregation. `sample::aggregate` requires identical reference/configuration identities, rejects duplicate input SHA-256 values, and sorts reads by SHA-256 independently of CLI trace order. The top-level read registry retains the reviewer-facing filename stem, the stable SHA-256 that identifies the read everywhere downstream ([[ADR-0090]]), and the concise selected-alignment summary (orientation, callable bases/identity, gap opens, unresolved and masked bases, mapped and callable reference segments, and origin-wrap state).

A read with fewer calls than `signal_processing.window_size_bases` is set aside before any stage runs and is recorded the same way, with reason `too_short` and without integrity or callability. A read with fewer than `callability.minimum_callable_calls` unmasked calls stops after its callability stage: it is logged as `sample_read_rejected` and recorded in a separate rejected-read registry, sorted by SHA-256, with its name, SHA-256, integrity, callability view, and rejection counts. It has no alignment and contributes no coverage, overlap, locus, variant, or notation evidence; the other reads continue. Reviewer names are unique across both registries: reads that share a stem (one file name in several folders) each gain `@` and the shortest prefix of their SHA-256, at least eight characters, that tells them apart. A trace content may appear only once, in one of them (`DuplicateTrace`, also when every read is rejected), and a sample whose every read is rejected fails typed.
