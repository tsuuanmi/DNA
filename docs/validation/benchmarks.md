# Performance Evidence

The current non-functional requirement sets a representative runtime/memory target for 500–1,000 base analysis; see [quality attributes](../requirements/quality-attributes.md).

A release record must document the host/toolchain and measured runtime/peak memory used to support that claim. Benchmark results are evidence and must not be copied into the requirement itself.

No benchmark number should be treated as portable across hosts or toolchains without an explicitly comparable setup.


## Corpus throughput evidence

High-throughput validation MUST report enough context to make throughput
measurements interpretable:

- DNA revision and Rust toolchain;
- host CPU and available logical/physical cores when known;
- worker count;
- storage location/type relevant to ABIF reads;
- trace count and representative trace-size distribution;
- total wall-clock time, traces/second, and failure count;
- whether reference/configuration were prepared once or reloaded per trace.

Corpus throughput is measured at the outer trace/sample orchestration boundary.
The scientific kernel remains deterministic per trace; benchmark code must not
change scientific thresholds or semantics to improve throughput.
