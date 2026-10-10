# DNA Reference-Free Basecall Output

`dna basecall <trace.ab1>` reads the strict configuration selected by
`DNA_CONFIG` or `config/dna.toml`, runs canonical ABIF decode, signal-derived
base re-calling, observational signal analysis, read callability, relative quality, and the trim interval,
then atomically creates `results/<trace-stem>.basecalls.json`. It does not load a
reference, align, or call variants. Existing targets are never overwritten.

The authoritative contract is
[`schemas/basecalls-v3.schema.json`](schemas/basecalls-v3.schema.json); a synthetic
example is [`examples/basecalls-v3.example.json`](examples/basecalls-v3.example.json).
Every object is closed and `schema_version` is `dna.basecalls/v3`; v3 adds the
read callability view to `read`.

The versioned v3 schema retains `ploc_*` field names for ABIF-origin integrity
evidence. Internally, DNA projects `PLOC.2` into canonical Sanger locus positions
at the input boundary and scientific stages use locus terminology.

## Fields

- `provenance`: input AB1 SHA-256, complete strict configuration SHA-256, and
  `plugins`: the plugins that produced the document ([ADR-0069](../decisions/adr/0069-plugin-first-modality-core-post-calling.md)), in execution order, each with its `id`, `family` (`modality`, `core`, or `post_calling`), and method `version` ([versioning](../governance/versioning.md#plugins)). `basecall` runs only the `sanger` plugin.
  Software/build identity, the trace filename, local paths, timestamps, and host
  data are omitted; software/build provenance is deferred until a stable
  versioning strategy is defined.
- `read.call_count`: number of canonical Sanger call loci decoded from `PLOC.2`.
- `read.primary`: strongest conservative signal-derived base at each locus.
- `read.ambiguity`: canonical/IUPAC ambiguity symbol at each locus.
- `read.retained`: the primary sequence inside `read.trim`.
- `read.trim`: 0-based half-open call interval `[start, end)`: the callable span
  plus at most `variant_calling.read_end_margin` adjacent dephased calls per side.
- `read.callability`: the `dna.read_callability/v1` view ([method](../design/callability.md)): `callable_span` (0-based half-open calls from the first to the last unmasked call, empty when everything is masked), `segments` (ordered 0-based half-open call intervals partitioning the read, each with `state` `in_phase`/`dephased`/`mixed`/`weak`/`irregular`, `after_repeat`, and `shadow_offsets` on dephased segments), and `masked_calls`. It sets `trim` and therefore `retained`.
- `signal_quality.integrity`: ABIF-origin locus/vendor-series cardinality evidence, merged repeated-position count, adjacent
  locus-spacing summary, exact signed-16-bit clipping count, and optional
  maximum-to-median corrected event-signal ratio. These observations do not
  reclassify artifacts or alter the read.
- `signal_quality.noisy_regions`: merged observation-only call/sample intervals
  and their minimum primary SNR. Individual rolling windows are omitted.
- `warnings`: unresolved-primary and multi-channel-unresolved counts, vendor
  disagreement count when optional vendor calls are available, vendor/locus
  cardinality mismatch count, and exact clipped-channel-sample count.

The primary and ambiguity sequence lengths equal `call_count`; trim bounds lie
within that count; and `retained` equals the primary sequence slice selected by
the trim interval. These cross-field invariants are enforced by typed Rust
construction and integration tests because JSON Schema cannot express them all.

## Interpretation and privacy

The result contains complete sequence strings and can identify a sample. It must
follow the same approval, storage, and redistribution policy as its source AB1.
The relative quality score, the rolling SNR method, and the
callability phase states are not Phred-calibrated error probabilities, mixture
fractions, or artifact classes. This output makes no genotype,
heteroplasmy, phase, pathogenicity, or clinical claim.

Operational records append to `$DNA_LOG_DIR/<trace-stem>.log` (default
`logs/`). They include stage metrics, timings, warning counts, and failures, but
never sequence strings, peak arrays, or JSON bodies.
