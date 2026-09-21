# DNA Reference-Free Basecall Output

`signal basecall <trace.ab1>` reads the strict configuration selected by
`DNA_CONFIG` or `config/dna.toml`, runs canonical ABIF decode, signal-derived
base re-calling, observational signal analysis, and relative quality trimming,
then atomically creates `results/<trace-stem>.basecalls.json`. It does not load a
reference, align, or call variants. Existing targets are never overwritten.

The authoritative contract is
[`schemas/basecalls-v2.schema.json`](schemas/basecalls-v2.schema.json); a synthetic
example is [`examples/basecalls-v2.example.json`](examples/basecalls-v2.example.json).
Every object is closed and `schema_version` is `dna.basecalls/v2`.

## Fields

- `provenance`: input AB1 SHA-256 and complete strict configuration SHA-256.
  Software/build identity, the trace filename, local paths, timestamps, and host
  data are omitted; software/build provenance is deferred until a stable
  versioning strategy is defined.
- `read.call_count`: number of decoded PLOC call loci.
- `read.primary`: strongest conservative signal-derived base at each locus.
- `read.ambiguity`: canonical/IUPAC ambiguity symbol at each locus.
- `read.retained`: the primary sequence inside `read.trim` after end trimming.
- `read.trim`: 0-based half-open call interval `[start, end)`.
- `signal_quality.integrity`: PLOC/vendor-series cardinality evidence, adjacent
  PLOC spacing summary, exact signed-16-bit clipping count, and optional
  maximum-to-median corrected event-signal ratio. These observations do not
  reclassify artifacts or alter the read.
- `signal_quality.noisy_regions`: merged observation-only call/sample intervals
  and their minimum primary SNR. Individual rolling windows are omitted.
- `warnings`: unresolved-primary and multi-channel-unresolved counts, vendor
  disagreement count when optional vendor calls are available, vendor/PLOC
  cardinality mismatch count, and exact clipped-channel-sample count.

The primary and ambiguity sequence lengths equal `call_count`; trim bounds lie
within that count; and `retained` equals the primary sequence slice selected by
the trim interval. These cross-field invariants are enforced by typed Rust
construction and integration tests because JSON Schema cannot express them all.

## Interpretation and privacy

The result contains complete sequence strings and can identify a sample. It must
follow the same approval, storage, and redistribution policy as its source AB1.
The relative score used during trimming and the rolling SNR method are not
Phred-calibrated error probabilities. This output makes no genotype,
heteroplasmy, phase, pathogenicity, or clinical claim.

Operational records append to `$DNA_LOG_DIR/<trace-stem>.log` (default
`logs/`). They include stage metrics, timings, warning counts, and failures, but
never sequence strings, peak arrays, or JSON bodies.
