# Implemented Rust Source Layout

Every listed `.rs` file is implemented and has an exact same-path `.md` manual under `docs/src/`.

```text
src/
├── lib.rs
├── main.rs
├── cli/mod.rs
├── error/mod.rs
├── logger.rs
├── checksum.rs
├── locus.rs
├── config/{mod,defaults,types,load}.rs
├── model/{mod,coordinate,nucleotide,trace,basecalls,locus_evidence,signal,quality,reference,alignment,variant,result,basecall_result,read_observation,sample_evidence,sample_result}.rs
├── trace/{mod,reader,abif,decode}.rs
├── reference/{mod,fasta}.rs
├── basecalling/{mod,iupac,peak,call}.rs
├── signal_processing/{mod,features,integrity,locus_evidence,regions,statistics}.rs
├── quality_control/{mod,penalty,quality,trim}.rs
├── alignment/{mod,scoring,gotoh,traceback,canonical,orient}.rs
├── variant_calling/{mod,mapping,extract,normalize,filter}.rs
├── sample/{mod,aggregate,call_evidence,contribution,coverage,loci,nucleotide_support,overlap,profile_geometry,variants}.rs
├── report/{mod,json,basecall,sample,signal,variant,atomic}.rs
└── pipeline/{mod,input,read,observation,analyze,basecall,sample,sample_metrics,sample_reads}.rs
```

## Dependency direction

```text
main -> cli -> lib dispatcher -> pipeline
pipeline shared read -> config + trace + basecalling + signal_processing + quality_control
pipeline observation -> reference + alignment + variant_calling; analyze/sample reuse observation; basecall -> no reference
pipeline commands -> report + logger
all stages -> model + error
model -> no filesystem, CLI, or algorithm module
report -> completed models; no scientific computation
```

## File ownership

- `reader`/`abif` parse generic binary structure; `decode` knows required scientific tags.
- `checksum` provides the shared stable SHA-256 identity used by config, trace, and reference loading.
- `locus` owns the shared PLOC-defined sample-window geometry used by evidence stages.
- `peak` selects basecalling peaks inside shared locus windows; `call` classifies them; `iupac` maps ambiguity.
- `features` estimates rolling sample-domain SNR; `integrity` derives whole-trace PLOC/vendor/clipping/event-scale observations; `locus_evidence` derives basecall-independent per-locus A/C/G/T evidence profiles; `statistics` owns shared robust local statistics; `regions` merges candidate-noisy intervals without changing calls.
- `penalty`, `quality`, and `trim` keep distinct QC responsibilities.
- `gotoh` computes DP; `traceback` reconstructs rows and owns reusable alignment metrics; `canonical` score-verifies repeat-equivalent 3'/right-most gap topology; `orient` applies strand/topology policy; `scoring` centralizes arithmetic.
- `extract` finds primary-sequence differences; `mapping` binds original calls to aligned reference positions; `normalize` builds minimal anchored alleles while preserving canonical alignment placement; `filter` applies configured region and supporting-signal eligibility.
- `logger` appends timestamped per-operation operational records without entering scientific stages or JSON.
- `json` assembles analysis v7 and owns shared serialization; `basecall` assembles basecalls v2; `sample` projects sample-evidence v8; `signal` is the shared integrity/noisy-region projection; `variant` projects mapped analysis calls; `atomic` is the one no-overwrite publisher.
- `input` loads command-specific resources and keeps reusable sample science inputs separate from publication targets; `read` sequences reference-independent stages; `observation` owns one authoritative reference-guided read path; `sample_reads` reuses that path across sample operations; `sample/aggregate` validates and orders reads; `sample/coverage` builds run-length total/forward/reverse reference coverage topology; `sample/overlap` builds the Tracy-derived pairwise overlap/admission graph; `sample/loci` builds sparse production reference-coordinate locus evidence; `sample/call_evidence`, `sample/contribution`, `sample/nucleotide_support`, and `sample/profile_geometry` own call projection and threshold-free nucleotide evidence geometry; `sample/variants` aggregates normalized variant support; `sample_metrics` owns operational sample summaries; command modules own orchestration/publication.

VCF has no source file or compatibility path in the MVP.
