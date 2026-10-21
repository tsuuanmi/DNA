---
implements: [SRS-ALN, SRS-API, SRS-CALL, SRS-CONS]
decided-by: [ADR-0069, ADR-0076, ADR-0080, ADR-0082, ADR-0084]
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

A read the sample commands set aside before the core still gets its document:
one too short for a signal window or with too few callable calls is
recorded with its `rejection` instead of its calls (below). A read the core
cannot place has an ordinary document: the core rejects it again when it
reads it.

The authoritative contract is
[`schemas/read-evidence-v1.schema.json`](schemas/read-evidence-v1.schema.json); a
synthetic example is
[`examples/read-evidence-v1.example.json`](examples/read-evidence-v1.example.json).
Every object is closed and `schema_version` is `dna.read_evidence/v1`.

## Fields

- `provenance`: as in [`dna.basecalls/v3`](basecalls.md), the input AB1
  SHA-256, the configuration SHA-256, and `plugins` (`evidence` runs only the
  `sanger` plugin), and the trace's file `name` beside its SHA-256. The
  sample documents name a read by the stem of that name, so a read from its
  document has the name its trace gives it, whatever the document is called.
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
  - `unreliable`: the spans of `calls` whose differences the modality
    asserts are not to stand unchecked, each with its `reason`:
    `poor_start_signal`, the read's start before its first run of sound calls
    (SRS-CALL-024). The calls keep their verdicts. The core reads the spans
    only to end a region's analysed range (SRS-CONS-012) and to leave the
    calls of a placement tie unplaced (SRS-ALN-006).
  - `tracts`: the homopolymer tracts the modality decodes from the read's own
    signal (SRS-CALL-025, [[PROP-0012]]), one per distinct first call, in
    onset order. Evidence only: the core reads none of it. Each tract has:
    - `onset`, the call where the read leaves phase;
    - `calls`, the calls the decode fits from the tract's first call;
    - `extent`, the calls that read the run or the mark;
    - `run` and `mark`, the run base and the chosen mark on the read strand,
      with `mark` null for a single run;
    - `candidates`, every structure fitted with its `residual`: the single
      run first, then each other base in A, C, G, T order;
    - `residual`, `sigma`, and `autocorrelation` of the chosen fit;
    - `first` and `total`, the first run's and the total's length marginals,
      the only identifiable ones. Each lists `lengths` with their `shares` of
      the signal and standard `errors`, the `dominant` and `next` lengths, the
      margin between them `z` and `z_adjusted` (null when a single length
      shows), and the `margin`: `clear` (a single length, or a margin at
      least `trace.near_equal_z`), `near_equal` (below it), or `undefined`
      (two lengths without a finite margin, never clear);
    - `cleaned`, the chosen structure refitted with the read's own base
      signatures and neighbour `leakage`, with its `residual` and `first`
      marginal; null when the refit does not fit.

    Shares are signal shares, never heteroplasmy levels.
- `signal_quality` and `warnings`: the trace integrity, merged noisy regions,
  and the modality's warning counts, as in `dna.basecalls/v3`.
- `calls`, one per call in call order. First the kernel contract the core reads:
  - `index`, `base` (read strand; `N` when unresolved), and `profile`, the
    basecall-independent A/C/G/T weights when the call has signal, at full
    precision;
  - `placement`: `kind` (`scored`, `anchoring`, `unresolved`) and, for a masked
    call, the modality's `reason`;
  - `identity` and `indel`: the modality's verdict on each use, `edge` and the
    `withheld` reasons, verbatim. A Sanger call inside an electrophoretic
    spike is withheld for both with `signal_spike` (SRS-CALL-027);
  - `bound`: whether the call shows where a run of another base ends
    (`in_phase`, `anchored`, `withheld`), and `post_repeat` (ADR-0076). A
    Sanger call beside a junction whose spacing marks a missing or an extra
    peak is `withheld` (SRS-CALL-026), and so is a call inside a spike
    (SRS-CALL-027);
  - `second`: the second base of the call's signal, its share in `[0, 0.5]`,
    and `stands_out` (the share exceeds `signal_processing.second_floor`);
  - `second_unexplained` and `reference_only` (ADR-0080).

  Then `signal`, the Sanger evidence the verdicts rest on: the trace sample
  `scan` of the call's locus, the selected A/C/G/T peak of each channel in
  its window (`channels`, read strand), the A/C/G/T `peaks` at its primary
  event on the read strand when it has one, its `relative_quality`, and the
  `cause` of its second peak when it has one (SRS-CALL-021). A call has a
  second peak when its second base stands out, or when it is masked as mixed
  and its own second peak is small (`mixed_stretch`, its second base then not
  standing out). The review reads `scan`, `relative_quality`, the second peak
  with its cause, and `channels`, from which the run signal counts a read's
  calls.

Every number is written at full precision, so the documents a later step
writes from it are those it writes in process. Profile weights matter most:
the aligner quantizes them to 1/1024, so three decimals change the quantized
weight of about a quarter of all calls, and placed two reads in a thousand
differently. Second-base shares are judged against their threshold by the
review, and the review and the core read mixture fractions and residuals.
The document is the same for the same trace and configuration, and its
bases, trim, callability, signal quality, and warnings equal the
`dna.basecalls/v3` document's.

### A read set aside

A read the Sanger plugin sets aside has `provenance` and, in place of `read`,
`signal_quality`, `warnings`, and `calls`, its `rejection` as the
[[reference/sample-evidence/reads-coverage|sample read registry]] reports it:

- `too_short`, with its `calls` and `minimum_calls`
  (`signal_processing.window_size_bases`): no stage ran, so the document has
  nothing else;
- `callable_calls_below_minimum`, with its `callable_calls` and
  `minimum_callable_calls` (`callability.minimum_callable_calls`), and the
  read's `integrity` and `callability` as in `dna.basecalls/v3`.

A synthetic example is
[`examples/read-evidence-v1.rejected.example.json`](examples/read-evidence-v1.rejected.example.json).

## As an input

The document is also the input of `dna analyze --evidence <file> --reference
<fasta>` (SRS-API-018): the core places and calls the read from it alone,
without the trace, and writes the [analysis document](analysis/README.md)
that `dna analyze <trace>` writes for the same trace and configuration, byte
for byte, under `results/<trace-stem>.json` (the file stem without
`.evidence`). The document of a read set aside fails the command as the
analysis of its trace does, with the same error.

`dna sample`, `dna consensus`, and `dna review` take a sample's documents
with `--evidence <file>...` in place of its traces (SRS-API-019): the core
places every read from its document, a read set aside is rejected with its
recorded cause, and the commands write the documents the trace-based
commands write for the same traces and configuration, byte for byte,
`--signal` included. Traces and `--evidence` exclude each other.

The reader:

- accepts only `dna.read_evidence/v1`, each object closed, at most 16 MiB;
- requires the configuration SHA-256 of the loaded configuration and this
  build's evidence plugins in provenance, so the analysis records one
  configuration and the plugins that made it;
- requires a trace file name without a directory, the calls in call order,
  every reason label one the Sanger plugin raises (a core label or an unknown
  one is refused), each verdict's reasons listed once, callability segments
  that partition the calls, the trim, noisy regions, and warnings consistent
  with the calls, and a second-peak cause exactly when the call has a second
  peak;
- requires of a read set aside either a read or a rejection, never both, the
  loaded configuration's minimum, a count below it, and callable and masked
  calls that add up to the read's calls;
- rebuilds the kernel's `ReadEvidence` through its own constructors (reason
  order, informative interval, ends, mixtures, unreliable spans, tracts), so
  every kernel invariant is checked again.

Numbers are parsed correctly rounded, so every value written is read back
bit for bit. Each failure is a typed `invalid read evidence` error and writes
no output.
