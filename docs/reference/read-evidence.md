# DNA Read Evidence Output

`dna evidence <trace.ab1>` reads the strict configuration selected by
`DNA_CONFIG` or `config/dna.toml`, runs the same reference-free Sanger read
processing as `dna basecall`, then the Sanger evidence adapter, and atomically
creates `results/<trace-stem>.evidence.json`. It writes the read's per-call
evidence exactly as the core receives it (`ReadEvidence`,
[ADR-0069](../decisions/adr/0069-plugin-first-modality-core-post-calling.md)),
so the modality → core boundary can be checked on its own
([PROP-0011](../proposals/0011-step-contracts.md)). It does not load a
reference, align, or call variants. Existing targets are never overwritten.

The authoritative contract is
[`schemas/read-evidence-v1.schema.json`](schemas/read-evidence-v1.schema.json); a
synthetic example is
[`examples/read-evidence-v1.example.json`](examples/read-evidence-v1.example.json).
Every object is closed and `schema_version` is `dna.read_evidence/v1`.

## Fields

- `provenance`: as in [`dna.basecalls/v3`](basecalls.md): the input AB1 SHA-256,
  the configuration SHA-256, and `plugins` (`evidence` runs only the `sanger`
  plugin).
- `read`:
  - `call_count`: the number of calls;
  - `informative`: the 0-based half-open call interval the modality retains;
  - `ends`: `unvouched` (the read's ends bound its evidence) or `vouched`;
  - `reasons`: the modality's reason labels in the order the core reports
    them; every reason a verdict holds is listed;
  - `callability`: the read callability view, as in `dna.basecalls/v3`.
  - `mixtures`: the stretches after a run that the modality decomposes into
    the read's own sequence at whole-call offsets (ADR-0082): the stretch
    `calls`, the `run_calls` of the run whose length varies, whether the
    components are `explained` at the read's own noise level, the `residual`,
    the `components` (offset and fraction, the dominant sequence first), and
    the `unexplained` calls. How each component compares with the reference is
    a read assertion of the core, in `dna review` (`reads[].mixtures`).
- `calls`, one per call in call order. First the kernel contract the core reads:
  - `index`, `base` (read strand; `N` when unresolved), and `profile`, the
    basecall-independent A/C/G/T weights when the call has signal;
  - `placement`: `kind` (`scored`, `anchoring`, `unresolved`) and, for a masked
    call, the modality's `reason`;
  - `identity` and `indel`: the modality's verdict on each use, `edge` and the
    `withheld` reasons, verbatim;
  - `bound`: whether the call shows where a run of another base ends
    (`in_phase`, `anchored`, `withheld`), and `post_repeat` (ADR-0076);
  - `second`: the second base of the call's signal, its share in `[0, 0.5]`,
    and `stands_out` (the share exceeds `signal_processing.second_floor`);
  - `second_unexplained` and `reference_only` (ADR-0080).

  Then `signal`, the Sanger evidence the verdicts rest on: the trace sample
  `scan` of the call's locus, the `height` of its strongest channel peak, its
  `relative_quality`, and the `cause` of its second peak when it has one
  (SRS-CALL-021).

Shares and weights are rounded to three decimals. The document is the same
for the same trace and configuration, and its bases and callability equal the
`dna.basecalls/v3` document's.
