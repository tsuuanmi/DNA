# DNA JSON-Only MVP Delivery Record

## Objective

Deliver one auditable AB1-to-primary-difference vertical slice in Rust. The implementation accepts one canonical analyzed ABIF/AB1, one strict config, and one single-record FASTA; then decodes, re-calls, trims, aligns, normalizes variants, and atomically publishes one JSON result.

## Completed implementation phases

1. Repository governance, one-to-one source manuals, strict configuration, versioned output schema, and reference identity.
2. Validated config/domain models, checked ABIF directory parsing, canonical tag decode, and FASTA loading.
3. PLOC-window signal re-calling, ambiguity evidence, relative quality, and end trimming.
4. Bounded affine-gap semi-global alignment, deterministic traceback, strand selection, and circular reference projection.
5. Primary SNV/small-indel extraction, linear/circular normalization, typed JSON assembly, and atomic no-overwrite output.
6. Observation-only rolling signal analysis and merged candidate-noisy regions.
7. Compact `dna.analysis/v7` projection retaining provenance, read/trim, trace-integrity and merged noisy-region evidence, alignment, normalized variants, and reference-oriented base/peak/quality evidence without implementation-coordinate noise.
8. Reference-free `signal basecall` using the same validated decode/re-calling/signal/QC path and publishing one typed `dna.basecalls/v2` JSON result with trace-integrity evidence but without reference, alignment, or variants.
9. Clean external batch reruns with complete preflight/build before selected-only
   destructive cleanup.
10. Multi-read `signal sample` analysis over independently placed reads, publishing compact `dna.sample_evidence/v8` with per-read trace integrity, a deterministic read registry, run-length total/forward/reverse coverage topology, Tracy-derived pairwise overlap/admission evidence, sparse differential loci with factorized support topology plus per-read A/C/G/T profile/noisy context, and normalized variants with factorized read/eligibility/orientation support topology.
11. Synthetic malformed/unit/end-to-end tests, deterministic output checks, and
   CI contract gates.

## Command and output

```text
signal basecall <trace.ab1>
signal analyze <trace.ab1> --reference <reference.fasta>
signal sample <sample-id> <trace.ab1>... --reference <reference.fasta>
```

A successful core CLI invocation creates the command-specific basecall, analysis, or `results/<sample-id>.sample.json` result and appends a separate operational log. Core operation failure creates no JSON result, and an existing target is never overwritten. The external batch wrapper has separate selected-cleanup semantics documented in `data.md`; a later batch failure may leave partial new outputs.

## Acceptance

- canonical ABIF records are bounds-checked and inconsistent arrays are rejected;
- compact analysis v7 JSON exposes input/reference/configuration hashes, call count/trim, merged noisy regions, alignment summary, normalized variants with reviewer-facing reference-oriented peak/quality evidence, and warning counts; effective parameters remain in strict configuration schema v5;
- internal coordinates and external 1-based variants are explicit;
- circular rCRS origin-spanning reads are representable;
- analysis v7 omits filenames, full sequences/windows/gapped rows, method constants, full peaks, and vendor data; basecalls v2 includes full sequences but omits reference/alignment/variant evidence; sample evidence v7 factors read identity/orientation/coverage once, emits coverage topology plus concise overlap/admission evidence, differential loci, and normalized variant support, and omits consensus/interpretation; no current contract emits compatibility output, genotype, heteroplasmy fraction, clinical meaning, VCF/BCF, or hidden regional correction;
- format/check/Clippy/tests/rustdoc/schema/TOML/reference/docs-mirror gates pass;
- batch cleanup is limited to fully preflighted selected sample directories and matching logs, rejects ambiguity/collisions/symlinks, and preserves unselected artifacts;
- approved real-trace validation is recorded before a scientific release claim.

See `pipeline.md` for formulas and `compatibility.md` for intentional Apollo divergences.
