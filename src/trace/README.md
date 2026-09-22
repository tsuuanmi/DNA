# Trace

Owns bounded, checked decoding of untrusted ABIF/AB1 input into a validated
chromatogram.

Key children: `reader.rs`, `abif.rs`, and `decode.rs`.

This module validates offsets, lengths, records, channels, positions, and available
vendor evidence; it does not perform re-calling, trimming, alignment, or reporting.

See [input requirements](../../docs/requirements/input.md) and
[ABIF decoding method](../../docs/design/abif-decoding.md).
