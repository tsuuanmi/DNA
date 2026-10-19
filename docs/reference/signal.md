# DNA Run Signal JSON

`dna signal <run-id> <trace>...` writes one deterministic document,
`results/<run-id>.signal.json`, identified as `dna.signal/v1` (SRS-API-016,
SRS-CALL-019, [PROP-0009](../proposals/0009-run-signal.md)). It records
the second-channel signal a run's reads share, so a later change can judge a
call against its own run rather than a fixed share. It reads the traces and
the configuration only: no reference, alignment, or truth.

The authoritative schema is
[`schemas/signal-v1.schema.json`](schemas/signal-v1.schema.json); the
synthetic example is
[`examples/signal-v1.example.json`](examples/signal-v1.example.json).

## Command

Every trace of the run is processed as `dna sample` processes it up to its
core evidence. A read the sample commands reject (too short, too few callable
calls) is counted in `rejected_reads` and contributes nothing. The output does
not depend on the order of the traces. The operation log is
`logs/<run-id>.log`.

## Fields

- `schema_version`, `run_id`.
- `provenance`: the configuration's SHA-256, the plugins, `selection`, the
  version of the call selection, and the method parameters the keys and
  shares were made with: `flank` (`signal_processing.run_signal_flank`) and
  `resolution` (`signal_processing.run_signal_resolution`). A command that
  reads a run signal refuses one whose `flank` or `resolution` differs from
  its configuration's.
- `traces`: the content SHA-256 of every trace of the run, sorted.
- `reads`, `rejected_reads`, `calls`: reads that contributed, reads rejected,
  and selected calls. Each selected call counts under its three other
  channels, so the keys' counts add up to three times `calls`.
- `keys`: one entry per key, sorted by `context`, `base`, `second`:
  - `context`: the read's primary calls from `flank` before to `flank` after
    the call, on the read strand (a 5-mer at `flank` 2);
  - `base`: the call's base (the context's centre); `second`: one of its
    other channels;
  - `shares`: `[share, calls]` pairs, by share: that channel's peak against
    the call's primary peak and itself, in `resolution` steps (0 to
    `resolution`; thousandths at `1000`; for the
    strongest other channel it is the call's second-base share,
    SRS-CALL-015), and the selected calls at it.

## Selection

A call is selected when it is informative and scored, its verdict on its base
is clear of a signal edge (including edge calls whose own signal took them
off it, SRS-CALL-017 and SRS-CALL-018) and free of the peak-floor and
relative-quality reasons, and its base and the `flank` calls on each side are
canonical. Calls judged mixed today are kept: their shadow is what a run
signal measures. A read's first and last `flank` calls have no context. The same
selection takes a read's own contributions out of a run signal, so a sample
can be judged without itself.

`dna review --signal` reads a run signal ([review](review.md)).
