# DNA Roadmap

## Active objective

The current goal is **not** to maximize feature count.

DNA first needs one small, trustworthy Sanger analysis path whose behavior is scientifically understandable, deterministic, and validated on approved real traces.

The active MVP is:

```text
AB1
 ↓
validated chromatogram decode
 ↓
primary base calling
 ↓
basic read quality / end trimming
 ↓
forward-or-reverse alignment to one short reference
 ↓
primary-sequence variant calling
 ↓
versioned JSON result
```

See [ADR-0021](adr/0021-scientific-core-confidence-floor.md).

## MVP acceptance priorities

### 1. Decode the chromatogram correctly

Before improving variant logic, prove that DNA reads the relevant ABIF structures correctly:

- analyzed A/C/G/T channels;
- FWO channel ordering;
- PLOC loci;
- optional vendor evidence where documented;
- malformed and inconsistent input rejection.

The decoded analyzed signal is the evidence foundation. Every later stage depends on it.

The ABIF container parser must also remain specification-conformant: generic directory parsing should reject truncated/overflowing payloads but must not reject a valid entry merely because `data_size` contains permitted trailing/padded bytes beyond `element_size × element_count`. Scientific tags such as `DATA.9-12`, `FWO_.1`, and `PLOC.2` are validated more strictly at their decoder boundary.

### 2. Base-call correctly enough to trust the read

The first scientific target is a conservative primary sequence plus explicit ambiguity. For the MVP this is deliberately **re-calling at `PLOC.2` loci**, not de novo discovery of base-event locations.

The basecaller must:

- derive calls from signal, not silently inherit vendor PBAS;
- retain enough signal evidence to explain a call;
- preserve unresolved states instead of guessing;
- behave deterministically at ties and thresholds.

Validation should compare against approved traces and independently established expected sequence, not only synthetic fixtures.

### 3. Trim only what can be justified

End trimming should remain simple and deterministic.

Current relative quality is not Phred-calibrated. MVP acceptance therefore requires evidence that trimming removes obviously poor tails without silently deleting trustworthy internal sequence.

More sophisticated confidence calibration is deferred.

### 4. Align the read correctly

For one short reference, DNA must:

- evaluate both orientations;
- select deterministically;
- preserve mapping back to original call indexes/PLOC;
- fail explicitly when placement is genuinely ambiguous.

Alignment correctness is more important than adding search/indexing features.

### 5. Validate simple SNVs first, then preserve known-good variant support

The first validation anchor is straightforward primary-sequence A/C/G/T substitution because it is the easiest variant class to reason about end-to-end.

This does **not** mean removing or disabling existing indel support. Current small insertion/deletion handling remains part of the working baseline when its normalization, call mapping, tests, and evidence remain sound.

MVP validation should answer:

- was the read placed correctly?
- does the mapped trace call support the alternate base?
- is the reference coordinate correct?
- is the allele reported on the reference strand?
- are unresolved calls excluded rather than guessed?

### 6. Keep output clear and stable

The JSON result should be easy to inspect and hard to misinterpret.

The schema should make explicit:

- input/config/reference identity;
- read and trim summary;
- alignment orientation and mapped region;
- simple variants and their supporting trace calls;
- warnings;
- coordinate conventions.

Schema stability and interpretability are more important than exposing every internal feature.

## Current supported baseline beyond the confidence floor

The repository already contains or explores functionality beyond this baseline, including:

- rolling SNR annotations;
- circular-reference handling;
- insertion/deletion extraction and normalization;
- batch orchestration;
- independently placed multi-read sample evidence with Tracy-derived pairwise overlap/admission, sparse differential loci, and normalized variant support;
- observation-only trace-integrity evidence for PLOC/vendor cardinality, PLOC spacing, exact clipping, and event-signal scale;
- pre-consensus run-length sample coverage topology with total/forward/reverse read depth;
- factorized normalized-variant support topology across observed/eligible and selected forward/reverse reads;
- richer per-call evidence.

These capabilities are not treated as disposable experiments merely because they are beyond the simplest validation path. Where their current behavior is coherent and tested, future work builds on them.

The distinction is:

```text
core confidence floor = what we validate first
current supported baseline = what the product already does intentionally and keeps
future scope = capabilities not yet implemented or not yet justified
```

A current capability may remain in the product contract while its real-trace evidence continues to strengthen.

## Production-readiness track

Engineering maturity is developed in parallel with scientific confidence. The repository now targets a production-grade delivery baseline with:

- pinned Rust release/MSRV toolchains;
- warnings-clean Rust CI and locked builds;
- Rust-only production source boundaries;
- cargo-deny, RustSec audit, and pull-request dependency review;
- immutable GitHub Action pins maintained by Dependabot;
- CodeQL and OpenSSF Scorecard;
- ABIF fuzz smoke plus scheduled adversarial campaigns;
- protected-branch/tag governance contract;
- SBOM, SHA-256 checksums, and cryptographic release attestations;
- explicit platform/versioning/security policies.

These controls make the software build and supply chain defensible. They do not replace the remaining scientific release evidence. The production-ready milestone is reached only when an exact tagged revision also has approved real-AB1 ground-truth regression, difficult-locus analysis, and runtime/peak-memory evidence recorded with the release template.

## Deferred until after MVP confidence

The following are intentionally deferred:

1. new or more complex indel models beyond the current supported behavior;
2. production mtDNA poly-C detector/state/recovery/weighting;
3. sample-level multi-read consensus;
5. quantitative heteroplasmy;
6. mixed-template or length-mixture decomposition;
7. haplogroup-based QC or inference;
8. calibrated quality/error probabilities;
9. ML-based calling, correction, or training export;
10. SCF, VCF/BCF, multi-contig references, and genome-scale indexing.

Research notes live under `docs/research/<topic>/` and are not active delivery commitments until promoted through the root ADR/SRS/contract process.

## Validation before expansion

Before promoting harder features, complete an approved real-trace regression baseline for the core path.

For each approved trace, retain:

- AB1 checksum;
- source/instrument context when available;
- primer or expected region;
- orientation;
- reference checksum;
- exact configuration;
- expected primary sequence or independently established truth;
- expected simple SNVs;
- explanation of disagreements;
- runtime and peak memory.

The first milestone is not "DNA supports everything Tracy supports".

The first milestone is:

> DNA can take an ordinary Sanger AB1, produce a trustworthy primary read, align it correctly, report straightforward variants correctly, preserve currently validated harder behavior such as small indels where supported, and explain the result through a stable schema.
