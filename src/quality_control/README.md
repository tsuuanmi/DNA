# Quality Control

Owns relative per-call quality analysis and deterministic end trimming.

Key children: `penalty.rs`, `quality.rs`, and `trim.rs`.

This module does not perform Phred calibration, reference alignment, or variant
filtering.

See [quality-control requirements](../../docs/srs/quality-control.md) and
[quality-control method](../../docs/methods/quality-control.md).
