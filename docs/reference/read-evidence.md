---
implements: [SRS-API, SRS-CALL]
decided-by: [ADR-0069, ADR-0076, ADR-0080, ADR-0082]
---

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
  - `trim`: the relative-quality trim interval, as in `dna.basecalls/v3`;
    the informative interval can start before it (SRS-CALL-014);
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
- `signal_quality` and `warnings`: the trace integrity, merged noisy regions,
  and the modality's warning counts, as in `dna.basecalls/v3`.
- `calls`, one per call in call order. First the kernel contract the core reads:
  - `index`, `base` (read strand; `N` when unresolved), and `profile`, the
    basecall-independent A/C/G/T weights when the call has signal, at full
    precision;
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
  `scan` of the call's locus, the `height` of its strongest channel peak, the
  A/C/G/T `peaks` at its primary event on the read strand when it has one,
  its `relative_quality`, and the `cause` of its second peak when it has one
  (SRS-CALL-021).

Second-base shares, mixture fractions, and residuals are rounded to three
decimals; the core reads none of them. Profile weights are written at full
precision: the aligner quantizes them to 1/1024, so three decimals change the
quantized weight of about a quarter of all calls, and placed two reads in a
thousand differently. The document is the same for the same trace and
configuration, and its bases, trim, callability, signal quality, and warnings
equal the `dna.basecalls/v3` document's.

## As an input

The document is also the input of `dna analyze --evidence <file> --reference
<fasta>` (SRS-API-018): the core places and calls the read from it alone,
without the trace, and writes the [analysis document](analysis/README.md)
that `dna analyze <trace>` writes for the same trace and configuration, byte
for byte, under `results/<trace-stem>.json` (the file stem without
`.evidence`). The reader:

- accepts only `dna.read_evidence/v1`, each object closed, at most 16 MiB;
- requires the configuration SHA-256 of the loaded configuration and this
  build's evidence plugins in provenance, so the analysis records one
  configuration and the plugins that made it;
- requires the calls in call order, every reason label one the Sanger plugin
  raises (a core label or an unknown one is refused), each verdict's reasons
  listed once, and the trim, callability, noisy regions, and warnings
  consistent with the calls;
- rebuilds the kernel's `ReadEvidence` through its own constructors (reason
  order, informative interval, ends, mixtures), so every kernel invariant is
  checked again.

Numbers are parsed correctly rounded, so every value written is read back
bit for bit. Each failure is a typed `invalid read evidence` error and writes
no output.
The `scan`, `height`, and `cause` of each call are checked for shape only:
the analysis reports none of them.
