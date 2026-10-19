# Data Flow

## Current Sanger production path

```text
Sanger ABIF
    |
    v
dna_sanger::abif
    |
    v
model::sanger::Chromatogram
    |
    v
read_processing
    |
    +--> basecalling
    +--> signal_processing
    +--> callability
    +--> quality_control
              |
              v
      variant_analysis::observation
              |
              v
   read_processing::evidence (Sanger adapter)
              |
              v
         ReadEvidence
              |
FASTA ------> alignment
              |
              v
        variant_calling <-- ReadEvidence
              |
              v
        ReadObservation = CalledRead + SangerAttachment
           /                 \
          v                   v
   analysis result     sample aggregation (CalledRead only)
                              |
                              v
                       sample evidence --> report joins SangerAttachment
```

A sample's consensus document is called from its decided edits, without
placing a segment again (ADR-0078):

```text
dna.consensus/v1 --> input::consensus --> dna_consensus::decided --> CalledRead --> called variants
```

Reviewed consensus sequences take a shorter path into the same core:

```text
consensus FASTA --> input::sequence --> ReadEvidence (vouched ends)
                                             |
FASTA ---------------------------------> dna_call::read_call
                                             |
                                         CalledRead --> called variants
```

A variants document reaches the post-calling plugins in another process:

```text
variants document --> input::variants --> normalization --> nomenclature
                                                                 |
                                                   conformance --+--> notation document
```

`variant_analysis::observation` owns the authoritative reference-guided
single-read scientific path. The public Rust Variant Analysis capability, CLI
analysis, and sample processing reuse that path rather than owning duplicate
implementations. Sample processing is one pass, `sample_analysis::reads`, that
places each trace once; the sample evidence and the consensus both start from
it, in the `sample` and `consensus` commands and in the public Sample Analysis
capability alike.

Alignment and variant calling see the read only as `ReadEvidence`
([ADR-0069](../decisions/adr/0069-plugin-first-modality-core-post-calling.md)):
the Sanger adapter turns calls, signal, callability, and quality into bases,
profiles, labelled masks, support vetoes, and the informative interval.

Completed typed scientific state is projected by the report layer and serialized
before atomic no-overwrite CLI publication. Logging remains operational side
evidence rather than part of deterministic scientific result contracts.

## Step contracts

Each step's output is a versioned contract that can be checked on its own
([PROP-0011](../proposals/0011-step-contracts.md)). The facade serializes every
contract; the crate computes what it holds. Each crate's source README names
its output contract, and `validate_docs_structure.py` checks that the named
schema exists.

| Step (crate) | Input contract(s) | Output contract | Command | Schema |
|---|---|---|---|---|
| shared contracts ([dna-kernel](../../crates/dna-kernel/src/README.md)) | — | none: no step; `dna.read_evidence/v1` serializes its `ReadEvidence` | — | — |
| read processing ([dna-sanger](../../crates/dna-sanger/src/README.md)) | ABIF trace | `dna.basecalls/v3` | `dna basecall` | [basecalls-v3](../reference/schemas/basecalls-v3.schema.json) |
| per-call evidence ([dna-sanger](../../crates/dna-sanger/src/README.md) → `ReadEvidence`) | ABIF trace | `dna.read_evidence/v1` | `dna evidence` | [read-evidence-v1](../reference/schemas/read-evidence-v1.schema.json) |
| run signal ([dna-sanger](../../crates/dna-sanger/src/README.md) `run_signal`) | a run's ABIF traces | `dna.signal/v1` | `dna signal` | [signal-v1](../reference/schemas/signal-v1.schema.json) |
| placement ([dna-align](../../crates/dna-align/src/README.md)) | `dna.read_evidence/v1`, reference FASTA | `dna.analysis/v10` `alignment` and each call's `column` | `dna analyze` | [analysis-v10](../reference/schemas/analysis-v10.schema.json) |
| read assertion ([dna-call](../../crates/dna-call/src/README.md)) | `dna.read_evidence/v1`, `dna.analysis/v10` `alignment` | `dna.analysis/v10` `variants` and each call's `trust`; component calls: none yet ([PROP-0010](../proposals/0010-length-mixture-decomposition.md)) | `dna analyze` | [analysis-v10](../reference/schemas/analysis-v10.schema.json) |
| read assertion over reviewed sequences ([dna-call](../../crates/dna-call/src/README.md), facade `input::sequence`) | consensus FASTA | `dna.variants/v1` | `dna call <FASTA>...` | [variants-v1](../reference/schemas/variants-v1.schema.json) |
| sample evidence ([dna-sample](../../crates/dna-sample/src/README.md)) | the sample's placed reads (`CalledRead`, in process) | `dna.sample_evidence/v10` | `dna sample` | [sample-evidence-v10](../reference/schemas/sample-evidence-v10.schema.json) |
| consensus ([dna-consensus](../../crates/dna-consensus/src/README.md) `build`) | the sample's placed reads (`CalledRead`, in process) | `dna.consensus/v1` and its FASTA | `dna consensus` | [consensus-v1](../reference/schemas/consensus-v1.schema.json) |
| decided calls ([dna-consensus](../../crates/dna-consensus/src/README.md) `decided`) | `dna.consensus/v1` | `dna.variants/v1` | `dna call <consensus.json>` | [variants-v1](../reference/schemas/variants-v1.schema.json) |
| normalization and nomenclature ([dna-normalize](../../crates/dna-normalize/src/README.md)) | `dna.variants/v1` | `dna.notation/v1` `notation` | `dna notation` | [notation-v1](../reference/schemas/notation-v1.schema.json) |
| conformance and EMPOP query ([dna-empop](../../crates/dna-empop/src/README.md)) | `dna.variants/v1`, as represented | `dna.notation/v1` `conformance`, `empop` | `dna notation` | [notation-v1](../reference/schemas/notation-v1.schema.json) |
| review (facade [report](../../src/report/README.md) `review.rs`) | the sample's analysis and call (in process), optional `dna.signal/v1` | `dna.review/v1` | `dna review` | [review-v1](../reference/schemas/review-v1.schema.json) |
| evaluation ([dna-eval](../../crates/dna-eval/src/README.md)) | panel manifest, reviewer truth, ABIF traces | none yet: `run.json` and TSV tables by convention ([increment 3](../proposals/0011-step-contracts.md#increment-3-evaluation-contracts)) | `dna-eval run` | — |

When the profile declares notation, `dna sample` and `dna call` also carry
dna-normalize's `notation` in their documents.

## Platform direction

Future modalities keep source-specific evidence until a biological difference
has actually been called or imported. The common convergence point is
`CalledVariantSet`, not a universal raw-alignment object.

```text
Sanger evidence --> ReadEvidence --> core caller -----+
consensus FASTA --> ReadEvidence --> core caller -----+
                                                      |
NGS evidence ----> NGS caller ------------------------+--> CalledVariantSet
                                                      |       /     |      \
VCF/BCF --------> validated importer -----------------+      v      v       v
                                                         direct  optional  optional
                                                          use   normalize  nomenclature
                                                                   |
                                                                   +--> optional nomenclature
```

The source-specific paths are intentionally different:

| Source | Path to called variants |
|---|---|
| Sanger ABIF | chromatogram -> base calling/signal/callability/QC -> `ReadEvidence` -> core caller (implemented) |
| consensus document (`dna.consensus/v1`) | the consensus's decided edits (implemented, `call`, `SampleAnalysis::call`) |
| assembled/consensus FASTA | `ReadEvidence` -> core caller (implemented, `call`) |
| FASTQ / NGS reads | read QC/preprocessing -> mapping -> NGS caller |
| BAM / CRAM | validated aligned-read evidence -> NGS caller |
| VCF / BCF | validated variant importer; raw calling is bypassed |

These rows are architecture direction, not claims that all sources are currently
implemented.

For current human mtDNA, one optional normalization policy is
haplotype-preserving 3'/right-most indel placement. Nomenclature is a separate
optional target-specific representation layer. A workflow may consume called
variants directly when it does not need either representation policy.

Sanger alignment still has its own deterministic right-most traceback
canonicalization for evidence placement and provenance. That alignment rule is
separate from any post-calling normalization selected by another workflow.

The detailed boundary is owned by
[variant lifecycle](variant-lifecycle.md) and ADR-0060.

Coordinate/strand semantics are canonical in
[reference/coordinates](../reference/coordinates.md); pipeline invariants are
canonical in [invariants](invariants/README.md).
