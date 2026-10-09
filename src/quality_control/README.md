# Quality Control

Owns relative per-call quality analysis and the trim interval.

Key children: `penalty.rs` (per-call ambiguity and spacing penalties),
`quality.rs` (relative scores), and `trim.rs` (the trim interval: the read's
callable span from [callability](../callability/README.md),
`dna.callable_span_trim/v1`).

This module does not decide where a read is callable, perform Phred
calibration, reference alignment, or variant filtering.

See [quality-control requirements](../../docs/requirements/quality-control.md) and
[quality-control method](../../docs/design/quality-control.md).
