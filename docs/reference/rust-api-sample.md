---
implements: [SRS-API, SRS-CALL, SRS-NOM, SRS-SAMPLE, SRS-VAR, SRS-VN]
decided-by: [ADR-0070, ADR-0076, ADR-0077, ADR-0078, ADR-0079]
covers:
  - crates/dna-kernel/src/variant.rs
  - crates/dna-normalize/src/
  - src/lib.rs
  - src/model/review_result.rs
  - src/pipeline/review.rs
  - src/pipeline/sample_call.rs
  - src/report/review.rs
  - src/sample_analysis/
  - src/variant_analysis/
---

# Rust Public API: Sample Analysis

Part of the [[reference/rust-api|Rust public API]]: the sample-level capability (`SampleAnalysis`) and its consensus, call, and review entry points.

## Sample Analysis

One placement pass over a sample's Sanger traces yields both the sample
evidence and the consensus:

```rust
pub fn analyze_sanger(
    traces: &[PathBuf],
    reference: &Path,
    config: &Path,
    workers: NonZeroUsize,
) -> dna::error::Result<SampleAnalysis>
```

through:

```rust
dna::sample_analysis::analyze_sanger(...)
```

It loads the traces, reference, configuration, and target profile as the
`sample` and `consensus` commands do, processes and places every trace once,
sets aside rejected reads (SRS-SAMPLE-027, SRS-SAMPLE-029, SRS-SAMPLE-030), aggregates the
admitted reads, and builds their consensus. The two commands run this same
pass, so the result describes exactly what they publish. A caller that needs
both no longer places every trace twice.

Reads are processed on up to `workers` scoped threads; with one worker they run
one after another on the calling thread. A caller that already runs many
samples in parallel passes `NonZeroUsize::MIN`; the CLI passes the available
parallelism. The result does not depend on the number of workers or on
scheduling. When a read fails other than by rejection, the failure of the
first such read in trace order is reported, as a sequential pass would: its
error is returned, or its panic (a defect) resumed. Workers stop taking reads
after the lowest failed one. Like `analyze_sanger` for one read, it opens no
log, derives no output path, serializes nothing, and publishes nothing.

```rust
#[non_exhaustive]
pub struct SampleAnalysis {
    pub reference: ReferenceIdentity,
    pub profile: ProfileIdentity,
    pub evidence: SampleEvidence,
    pub consensus: Consensus,
    pub reads: BTreeMap<String, SangerAttachment>,
    pub rejected: BTreeMap<String, RejectedRead>,
}

impl SampleAnalysis {
    pub fn regions(&self) -> &[[usize; 2]];
    pub fn observations(
        &self,
        first_1based: usize,
        last_1based: usize,
    ) -> dna::error::Result<Vec<LocusObservations>>;
    pub fn trust(&self, input_sha256: &str) -> Option<Vec<CallTrust>>;
    pub fn components(&self) -> Vec<ReadMixtures>;
    pub fn tracts(&self) -> Vec<ReadTracts>;
}

#[non_exhaustive]
pub struct ReadMixtures {
    pub read: String,
    pub sha256: String,
    pub forward: bool,
    pub mixtures: Vec<MixtureCalls>,
}

pub struct ReadTracts {
    pub read: String,
    pub sha256: String,
    pub forward: bool,
    pub tracts: Vec<PlacedTract>,
}

pub struct PlacedTract {
    pub reference: Option<usize>,
    pub onset_reason: Option<&'static str>,
    pub to_end: bool,
    pub tract: TractLengths,
}

#[non_exhaustive]
pub struct RejectedRead {
    pub input_name: String,
    pub input_sha256: String,
    pub cause: RejectionCause,
}

#[non_exhaustive]
pub enum RejectionCause {
    TooShort(ShortRead),
    Callability { sanger: SangerRejection, rejection: ReadRejection },
    Placement { sanger: SangerRejection, rejection: PlacementRejection },
}

impl RejectionCause {
    pub const fn sanger(&self) -> Option<&SangerRejection>;
}
```

- `evidence` is the core's sample evidence (`dna_sample`): the read registry in
  content-SHA-256 order, the rejected reads, coverage, overlaps, locus
  differences, and normalized variants with support and opposition. It carries
  the reference and configuration SHA-256.
- `consensus` is the core's consensus (`dna_consensus`): segments, sites, and
  summary; each site names its reads in registry order.
- `reads` holds each admitted read's Sanger evidence (`dna_sanger`): calls,
  signal analysis, callability, and quality control.
- `regions()` returns the reportable regions of the target profile the
  analysis ran with (inclusive, 1-based), so a caller reads them from the
  profile rather than repeating them.
- `trust(input_sha256)` returns, for the admitted read with that content
  SHA-256, the core's trust in each of its calls (`CallTrust`, `dna_call`):
  its mask (reason and onset, as the modality gave them), whether it is
  unresolved, the support vetoes raised on it by label, whether it is trusted
  under `read_end` and whether an ADR-0076 exemption applies, and what
  withholds trust (`Withholder`: the read's own end, its trim, or a masked
  call, before then after the call), empty when trusted, what the modality
  vouches for (`TrustScope`: `Full`, or `Base` for a leading call vouched
  for by its base only, SRS-CALL-014), and the call's second base and share
  (SRS-CALL-015). It is the verdict
  variant calling and the consensus reach, computed on demand; diagnostic
  evidence for evaluation and review, not a calling contract.
- `observations(first_1based, last_1based)` returns, for each position of
  the inclusive interval that an admitted read covers, a `LocusObservations`
  (`dna_sample`): the position, its reference base, and every covering read's
  observation exactly as `evidence.locus_differences` records it (state,
  reference-strand base, source call index, and reference-oriented A/C/G/T
  profile, with read indices into `evidence.reads`), whether or not a read
  differs there, but without the support aggregates of a differential locus
  (SRS-API-013). It is computed on demand from the placed reads the analysis
  keeps and changes nothing; uncovered positions have no record. The sample
  evidence keeps only differential loci; a caller that needs the minor
  channel at a position where every read calls the reference reads it here.
- `components()` returns, per admitted read in placement order, its phase
  mixtures (`MixtureCalls`, `dna_call::components`, re-exported here with
  `ComponentCall`): the modality's assertion (`PhaseMixture`: stretch, run,
  components with offset and fraction, residual, unexplained calls, and
  whether it is explained, SRS-CALL-023) and, for an explained mixture, each
  component compared with the reference on its own: its span (first and last
  positions in read order on the reference strand, the first exceeding the
  last across the origin), its bases there, and its observed variants
  (SRS-VAR-022). It is computed on demand, changes nothing, and is evidence
  for evaluation and review, not a calling contract.
- `read_evidence_document(trace, config)` returns a trace's
  `dna.read_evidence/v1` document in memory, as `dna evidence` writes it;
  `analyze_evidence_outcome(documents, reference, config, workers)` analyses
  a sample from such documents in memory, as `dna sample --evidence` does from
  files (SRS-API-019), and returns `EvidenceOutcome`: the analysis, or how
  many distinct reads were set aside. Reads from recorded evidence carry no
  Sanger products, so that analysis's `reads` and `rejected` are empty; its
  sample evidence lists the rejected reads. `dna-eval run --evidence-cache`
  uses both.
- `tracts()` returns, per admitted read in placement order, its tract
  decodes (`TractLengths` and `LengthMarginal`, `dna_kernel::read_evidence`,
  re-exported here; SRS-CALL-025), each with the reference position (1-based)
  its first call is aligned to, or `None`; the read's reason for masking the
  onset call, its own assertion of why it left phase there; and whether no
  later call is scored (`to_end`). It is computed on demand, changes
  nothing, and is evidence for evaluation, not a calling contract.
- `rejected` holds each rejected read's name, identity, and cause: fewer
  calls than one signal window, too few callable calls, or a placement the
  core rejected. The last two carry the read's integrity and callability
  (`SangerRejection`); a read too short to process has none.

`reads` and `rejected` are keyed by read content SHA-256, the identity the
core records use. Every coordinate field names its base (`_0based`,
`_1based`, `_0based_exclusive`). The published projections of these records
are `dna.sample_evidence/v10` and `dna.consensus/v1`. Notation is not part of
this result; the consensus call below derives it.

### Sample outcome

`analyze_sanger` fails with `SampleError::NoAdmittedReads` when every read is
rejected, so the caller learns how many but not which. A caller that must name
those traces, such as the evaluation harness asking for a rerun, uses the same
pass through:

```rust
pub fn analyze_sanger_outcome(
    traces: &[PathBuf],
    reference: &Path,
    config: &Path,
    workers: NonZeroUsize,
) -> dna::error::Result<SampleOutcome>

pub enum SampleOutcome {
    Analyzed(Box<SampleAnalysis>),
    Rejected(BTreeMap<String, RejectedRead>),
}

pub fn read_evidence_document(trace: &Path, config: &Path) -> dna::error::Result<Vec<u8>>

pub fn analyze_evidence_outcome(
    documents: &[(String, Vec<u8>)],
    reference: &Path,
    config: &Path,
    workers: NonZeroUsize,
) -> dna::error::Result<EvidenceOutcome>

pub enum EvidenceOutcome {
    Analyzed(Box<SampleAnalysis>),
    Rejected { reads: usize },
}
```

It returns `Analyzed` with exactly the analysis `analyze_sanger` returns, or
`Rejected` with every read, its cause, and its Sanger evidence, keyed by read
content SHA-256, when none is admitted. Every other failure is the error
`analyze_sanger` returns; `analyze_sanger` is this call with `Rejected` turned
into that error. The two cases are complete, so the enum is not marked
non-exhaustive.

### Consensus call

`SampleAnalysis::call` completes the recommended per-sample workflow,
`consensus` then `call` ([[reference/consensus#workflow|workflow]]), in process:

```rust
impl SampleAnalysis {
    pub fn call(&self, sample_id: &str) -> dna::error::Result<ConsensusCall>;
    pub fn review_checks(&self, call: &ConsensusCall) -> dna::error::Result<Vec<ReviewCheck>>;
}

#[non_exhaustive]
pub struct ReviewCheck {
    pub first: usize,
    pub last: usize,
    pub labels: Vec<&'static str>,
    pub fast: bool,
}
```

`review_checks` returns the `check` items of the review list `dna review`
would write for the analysis and that call (SRS-API-015), built by the same
rules: each item's first and last reference positions (1-based, inclusive)
its label names, sorted, and whether it is a fast check, one with a
`confirm` answer ([[reference/review]]). It writes nothing. The evaluation harness reads
it to measure the review workload.

It validates the sample identifier and calls the consensus's decided edits as
`dna call` calls the consensus document `dna consensus` would publish for that
sample (SRS-VAR-021, ADR-0078): each segment, named as in the consensus FASTA,
without being placed again, with the configuration, profile, and reference of
the analysis. It fails where `call` would fail, for example when the consensus
has no segment. Like `analyze_sanger`, it opens no log, derives no output path,
serializes nothing, and publishes nothing. A caller that runs many samples in
parallel calls `analyze_sanger` with `NonZeroUsize::MIN` workers, then `call`.

```rust
#[non_exhaustive]
pub struct ConsensusCall {
    pub callable_segments: Vec<ReferenceSegment>,
    pub notation: Option<Vec<NotationCall>>,
    pub empop: Option<EmpopExport>,
}

#[non_exhaustive]
pub struct NotationCall {
    pub position: usize,
    pub ordinal: usize,
    pub reference: Option<char>,
    pub alternate: Option<char>,
    pub text: String,
    pub convention: Option<Convention>,
}
```

- `callable_segments` holds the reference intervals the consensus sequences
  cover with callable bases (zero-based, half-open), in consensus read order:
  the `callable_reference_segments` of the `dna.variants/v1` reads.
  `ReferenceSegment` is the type of `variant_analysis`, re-exported from
  `sample_analysis`.
- `notation` holds the sample's calls in the profile notation, distinct and in
  reference order: the `notation.calls[].call` values of the `dna.variants/v1`
  document `dna call` would publish. It is `None` when the profile declares no
  notation. Each call gives the 1-based `position` of the changed base (or of
  the base an insertion follows), the insertion `ordinal` (`0` for a
  substitution or deletion), the `reference` base (`None` for an inserted
  base), the called `alternate` base (`None` for a deleted base), and the
  rendered `text` (`16189C`, `309.1C`, `249DEL`). `convention` is the
  profile report convention that labels the call (SRS-NOM-017), the
  `notation.calls[].convention` of the document; the call is kept either way.
  `Convention`, `ConventionKind`, and `ConventionLabel` are the profile types
  of `dna_kernel::profile`, re-exported from `sample_analysis`.
- `empop` is the EMPOP query `dna notation` would export from the
  `dna.variants/v1` document (SRS-NOM-018), or `None` when the profile
  declares no `[empop]`. `EmpopExport` (`ranges`, `mitotype`, `cut`, with
  `ranges_text()` and `mitotype_text()`), `Cut`, and `CutReason` are the types
  of `dna_empop::export`, re-exported from `sample_analysis`.
