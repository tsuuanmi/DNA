# Apollo Compatibility and Scientific Corrections

## Policy

Apollo C++ is evidence for intended behavior, not authority for known unsafe, inconsistent, or biologically misleading behavior. DNA preserves compatible behavior where it is well-defined and records intentional divergences.

## Mapping

| Apollo area | DNA | Treatment |
|---|---|---|
| `preprocessing/abif.h` | `trace` | strict bounds and canonical DATA.9-12/FWO/PLOC decode |
| ABIF `basecall` | `basecalling`; `signal basecall` | PLOC-window signal re-calling with corrected ties/ambiguity and one versioned JSON result; no Apollo FASTA/FASTQ/TSV compatibility surface |
| no Apollo equivalent | `signal_processing` | observation-only rolling SNR and candidate-noisy regions |
| quality helpers/`trim.h` | `quality_control` | safe penalty and end-trim behavior; score named relative |
| `alignment/gotoh.h` | `alignment` | bounded deterministic fixed-point profile-aware semi-global affine DP |
| primary subset of `variant.h` | `variant_calling` | SNV/small-indel extraction, normalization, and configured eligibility |
| `logger.rs` | `logger` | Apollo-style timestamp/level/source records written to per-trace files |
| `report/json.h` | `report` | versioned, nested, schema-governed JSON |

## Exact evidence targets

For an approved fixture, raw decoded bytes, channel remapping, PLOC positions, unique positive local maxima, unaffected penalty arithmetic, affine scores, and unambiguous normalized variants should compare exactly when preconditions match.

## Intentional divergences

- malformed or unequal ABIF records are rejected, never truncated;
- only canonical analyzed tags are accepted; there is no alternate-tag fallback;
- base locations come from PLOC, so DNA does not claim fully de novo location discovery;
- `P2BA.1` is ignored; only optional `PBAS.2` and `PCON.2` vendor evidence is consumed;
- exact strongest-channel ties become N instead of favoring a channel by iteration order;
- three qualifying channels keep the strongest primary with unresolved N ambiguity; four qualifying channels yield unresolved N for both primary and ambiguity;
- `dna.peak_recall/v3` requires each qualifying channel to pass the selected-peak ratio and the same ratio at the primary peak sample, preventing a remote channel maximum in the call window from creating ambiguity;
- per-channel PLOC fallback is kept and Apollo's collective midpoint rescue is not ported;
- the Rust Gotoh traceback breaks ties by state preference Match > Deletion > Insertion, which may differ from C++ tie behavior on equally scoring paths;
- alignments below `minimum_callable_bases` or `minimum_identity` fail rather than being silently accepted;
- rolling SNR annotations are a DNA method and do not alter Apollo-comparable peak calls or variants;
- Apollo's called/other peak-dominance ratio and Phred-like signal quality are not ported; DNA instead requires configured maximum-channel peak and uncalibrated relative-quality thresholds for SNV/inserted-base supporting calls;
- indels are capped at the configured 50 bp changed length; longer candidates are excluded rather than emitted;
- HV-region eligibility is expressed as strict TOML list-of-lists rather than hardcoded sample logic; there is no poly-C correction, rescue coordinate, or sample patch;
- zero maximum penalty yields maximum relative score rather than division by zero;
- the relative score is not labeled Phred; PCON remains distinct;
- N uses an explicit alignment score and cannot produce a reportable SNV;
- rCRS is treated as circular and origin-spanning coordinates are preserved;
- indels are normalized deterministically, including circular repeats;
- compact `dna.analysis/v6` JSON is the sole reference-guided single-read analysis result; it retains provenance, read/trim, merged noisy-region, alignment, normalized-variant, and warning summaries while omitting filenames, sequences, rolling windows, gapped rows, method constants, call indexes/PLOC coordinates, selected-peak positions/sources, vendor data, and redundant fields; append-only operational logs remain separate and no VCF/BCF or legacy JSON compatibility layer exists;
- no ConfirmFilter, PHFinder, genotype, allelic fraction, or two-allele decomposition.

None of these divergences keeps a legacy path: DNA emits a single
`dna.analysis/v6` result, consumes only the documented tags, and exposes no
compatibility switch, alias, or removed behavior.

The previous Apollo Rust port is not ground truth where it uses PBAS as the final result, sample-specific early-region logic, loose goldens, ignored write errors, or inconsistent coordinates/quality fields.

## Differential evidence

Local ignored AB1 files are candidates, not goldens. A compatibility claim requires approval, source/checksums, region/orientation, config identity, generating revision, expected fields, and an exact comparison rule. Missing/extra variants are never ignored.
