# Evaluation harness

`dna-eval` runs DNA over a fixed test panel and compares the calls with
reviewer truth ([PROP-0005](../proposals/0005-evaluation-feedback-loop.md),
[ADR-0075](../decisions/adr/0075-evaluation-harness-in-workspace.md)). It is
the only place in the repository for truth formats and comparison
(SRS-COMPAT-003). Its outputs are validation evidence about a DNA revision,
not DNA contracts. Source: [crates/dna-eval](../../crates/dna-eval/src/README.md).

## Commands

```bash
cargo run --release -p dna-eval -- run \
  --manifest <panel.tsv> --truth <statistics.sqlite> \
  --config config/dna.toml --reference references/rCRS.fasta \
  --output <new-run-directory> [--baseline <earlier-run-directory>] \
  --threads <1-256> [--triage] [--forms <earlier-run-directory>/forms.tsv]

cargo run --release -p dna-eval -- panel \
  --manifest <parent-panel.tsv> --truth <statistics.sqlite> \
  --version <v2> --created <YYYY-MM-DD> --output <new-panel.tsv>

cargo run --release -p dna-eval -- power \
  --dev-run <development-or-readiness-run-directory> --manifest <hold-out-panel.tsv> \
  [--base <sweep-policy> --candidate <sweep-policy> \
   [--true-rate <rate> --bound <rate>...]]

cargo run --release -p dna-eval -- approve \
  --prereg <committed-prereg.toml> --output results/<new-directory> \
  --threads <1-256>

cargo run --release -p dna-eval -- convention \
  --truth <statistics.sqlite> --profile config/profiles/human-mtdna-rcrs.toml \
  --reference references/rCRS.fasta --output <new-directory>
```

Every path is an argument. `run` needs a new or empty output directory, and
`panel` a new output file; nothing is overwritten. `panel` prints the new
body's checksum line (`<sha256>  <file> (body)`). The harness is built from
the workspace but not shipped in a release.

`run` refuses a manifest whose header names a panel other than `dev` or
`readiness`: a hold-out panel is read only by `approve`, once, under a
committed pre-registration ([approval](#approval)). `power` reads a
development or readiness run.

## Approval

`dna-eval approve` is the single approval run of the [approval
protocol](../proposals/0005-evaluation-feedback-loop.md). Its input is a
pre-registration (TOML) committed to the repository:

- **`[inputs]`** — the manifest and its body SHA-256, the panel truth digest,
  a pinned read-only copy of the truth database and its SHA-256, the
  configuration and its SHA-256, the reference, and the DNA revision.
- **`[analysis]`** — `bootstrap_draws` and `bootstrap_seed`.
- **`[utility]`** — `form` (`point` or `lower_bound`) and
  `minimum_pass_rate`.
- **`[[gates]]`** — a measure and an `upper_bound` each; every tested
  hypothesis must meet every gate.
- **`[[hypotheses]]`** — the fixed sequence in order: a name, `requires`
  (`first`, then `previous_pass`: the sequence stops at the first failure, so
  the family-wise alpha holds), the full policy
  (`short_calls`, `causes`, `one_strand` as `everywhere`, `off`,
  `outside_repeats`, or `outside_repeats_h<k>`, `proxy`, `minor` as one or
  more `<kind>@<t>` joined by `+` (any of which marks) or `off`, and
  `single_minor` as `<kind>@<t>` or `off`), and an optional `marginal`
  (`against` an earlier hypothesis, a measure, an `upper_bound`, the
  `projected_n` that `dna-eval power` gave at its lower bound, and
  `projected_from`: the development `run`, the SHA-256 of its
  `triage_policies.tsv`, and the `base` and `candidate` sweep policies
  projected).

The run stops before DNA when any value is still `"pending"`, a policy lacks
the `n` cause, a marginal `upper_bound` is below the count-based bound at 0
leaks among its `projected_n` samples (it cannot be met), the file is not
committed unmodified, the harness build is not the frozen revision or has
local changes, or the manifest, configuration,
truth copy, or panel truth digest differs from the frozen value. The output
directory must be new and lie under `results/`.

Per tested hypothesis, a sample passes when DNA called it and the policy
decides `pass`; failed samples never pass, and retired ones are left out.
Each gate bounds the leak rate (passed samples not concordant under the
gate's measure) three ways, at one-sided alpha 0.05:

- the count-based Clopper–Pearson upper bound;
- the Korn–Graubard upper bound of the corpus-weighted rate: the variance is
  the stratified variance of `estimate.rs` over the passed samples, each
  stratum's passed samples its sample, which gives the effective sample size
  (a zero rate uses the passed count); the degrees-of-freedom correction is
  omitted, and the bootstrap is what covers clustering by batch;
- the 95th percentile of a bootstrap that resamples the panel's batches with
  replacement, each draw the combined ratio `Σ w·leak / Σ w·pass` with fixed
  corpus weights, seeded by `SplitMix64(seed ^ draw · 0x9E3779B97F4A7C15)`; a
  draw without a passed sample counts as rate 1.

A gate passes only when all three are within its bound; at zero leaks the
bootstrap bound is 0 and proves nothing alone. The marginal endpoint bounds
the leaks under its measure among the samples a hypothesis passes and
`against` does not (count-based bound), and reports apart the leaks with a discordant position that the one-strand exception lets one strand read. The utility floor applies to the
corpus-weighted pass rate or its one-sided lower bound. A hypothesis passes
when its gates, marginal, and utility pass; the next one is tested only when
the previous one passed.

The outputs are `approval.json` (identities, every endpoint and verdict,
with each marginal's `projected_n` beside the realised newly passed count),
`approval.md` (aggregates only), `leaks.tsv` (every leak of every gate with
an empty `verdict` column for the adjudicators; the only file with sample
identifiers), and `prereg.toml`. Adjudication never changes the verdict.

`dna-eval power` projects an approval run before the freeze. It reads the
development run's `triage.tsv` and the hold-out manifest's stratum sizes,
never the hold-out truth. Per stratum it reports the development pass rate
and the passes it projects. For the whole panel it reports:

- the stratified rate: the development rates weighted by the panel's
  stratum sizes, with the variance `Σ W_h² p_h(1 − p_h) / n_h` over the
  development samples `n_h`;
- its Korn–Graubard one-sided 95 % lower bound, the gates' method: the
  Clopper–Pearson lower bound at the effective sample size;
- the expected passes at each, with the count-based bound attainable at
  0–10 leaks.

A panel stratum without development samples counts at rate zero. The
development rates come from the panel the rules were tuned on, so even the
lower bound may be optimistic.

With `--base` and `--candidate`, two policy names of the development run's
sweep, it projects the samples the candidate passes and the base does not,
the samples a marginal endpoint against the base counts. It reads them from
`triage_policies.tsv` and reports the same projection, with the marginal
bound attainable at 0–3 leaks. With `--true-rate` and one or more `--bound`,
it also reports, per bound, the most leaks among the lower-bound samples
that still meet it and the binomial chance of at most that many at the true
rate. The pre-registration records the lower-bound projection as
`projected_n`.

## Notation convention

`dna-eval convention` measures how far DNA's notation agrees with the way the
lab's reviewers write their own calls, without traces. Every truth sample
whose calls are plain A/C/G/T substitutions, deletions, and insertions is
rebuilt as anchored variants (runs of deleted positions as one deletion,
insertion ordinals `1, 2, …` without a gap as one insertion) and goes through
DNA's right-aligned [normalization](../design/variant-normalization.md) and
the profile's nomenclature windows, as DNA's own calls do, with a deletion
trade bounded by the windows and every base taken as read. A sample whose
written calls come back unchanged is a fixed point of DNA's notation; a
changed one shows a convention DNA writes differently. A sample with an
IUPAC call, or calls DNA cannot apply, is skipped.

The output directory, new or empty, holds `convention.tsv`: per change
pattern, most samples first and at most 200, its sample count, the
truth-only calls, and the DNA-only calls. It names no sample. The command
prints the samples checked, unchanged (with their share), changed, the
patterns, and the samples skipped. The truth is only read; it is never
rewritten or scored as DNA output.

## Inputs

**Panel manifest.** `#` comment lines with `key=value` fields come first, then
one tab-separated column header, then one row per sample.

- `body_sha256` is required: the SHA-256 of every byte after the column header
  line. A mismatch stops the run. `rows`, when present, must equal the row
  count. `panel`, `version`, `salt`, and from version 3 `truth_digest` (the
  panel truth digest) are recorded.
- The columns `sample_id`, `batch`, `stratum`, `weight`, `weight_corpus`,
  `anchor_candidate` (`0` or `1`), and `trace_paths` (`;`-separated) are read.
  The optional `retired` column (`0` or `1`, from panel version 2) marks
  samples kept for continuity but left out of the estimates. The optional
  `truth_digest` column (from version 3, 64 lowercase hexadecimal digits)
  records each sample's truth digest when the panel was built. Other columns
  are ignored. Sample identifiers are unique.
- Bounds: 256 MiB, 1,000,000 rows, and 64 traces per sample.

**Truth.** The lab's statistics database (`statistics.sqlite`, at most
4 GiB), opened read-only. It is the indexed copy of the lab's merged reviewer
results (`statistics.json`) and the single truth source of the harness. Its
`info` table must declare `schema_version` 1. Per panel sample the harness
reads the batch (`samples`), the reviewed alleles (`variants`: `pos`, `ref`,
`seq`, with `-` for an inserted or deleted base; at most 10,000 per sample),
the analysed intervals (`intervals`), and the sample's status in its batch
(`sample_status`). Sample flags and counts are not used.

**Baseline.** The `scores.tsv` of an earlier run. `sample_id` and `status`
are required; `strict_status` and `iupac_compatible_status` are compared when
present.

## Panel versions

`panel` derives a new version from a parent manifest and the database:

- membership is frozen: every parent row is kept, in order;
- a sample the database has no record of is retired (`retired = 1`, both
  weights `0`), and the weights of the other samples of its stratum are scaled
  in proportion, `w' = w · total / Σ staying w` (written with four decimals),
  so each stratum still sums to its population and keeps its differences
  within; a zero weight stays zero. When every staying weight of a stratum is
  zero, its weights stay zero and `panel` reports the stratum as having lost
  weight. Strata without a newly retired sample keep their weights as written;
- each row records its sample's truth digest (`truth_digest`), and the header
  the panel truth digest;
- the header records the panel, the new version, the salt, the creation date,
  the row and retired counts, and the parent version and body hash; the
  database path, file SHA-256 (`truth_database_sha256`), `created_at`, and
  sample count are recorded as provenance only. Header values escape `%`, `=`,
  and white space, so any path reads back.

Earlier versions stay runnable; report results on the rows both versions
share.

## Panel truth digest

The lab rebuilds the database in place, so its file hash changes with every
refresh even when no panel sample changes. A panel therefore records the
truth it was built against as digests of exactly what the harness uses:

- a sample's digest is the SHA-256 of a canonical text: `batch\t<batch>`,
  `status\t<status>` (empty when unrecorded), one `allele\t<pos>\t<ref>\t<seq>`
  line per allele with the position normalised, and one
  `interval\t<start>\t<end>` line per analysed interval, the alleles and
  intervals sorted, each line ending in a newline; a sample without a record
  digests `missing\n`;
- the panel's digest is the SHA-256 of one `<sample_id>\t<digest>` line per
  sample, in panel order.

`run` computes the digests of the current database. When the panel digest
equals the manifest's, the truth is `unchanged` and the database's file hash
and `created_at` are recorded in `run.json` for provenance only. Otherwise the
run is marked `truth changed`: the summary starts with a notice, the samples
whose digest differs, or that record none under a header that does, are
listed in `truth_changes.tsv` (written whenever the truth changed, with only
its header when no row differs; the notice then says so), and a new panel
version is due. A manifest without a digest (versions 1 and 2) reports
`not_recorded`.

## Running DNA

For each sample the harness removes repeated trace content (keeping the first
path of each content), then runs `dna::sample_analysis::analyze_sanger` with
one worker and `SampleAnalysis::call`: the `consensus` then `call` workflow in
process, with the configuration and reference given. Its notation calls become alleles (`309.1C` →
`309.1`, `-`, `C`; `249DEL` → `249`, `A`, `-`), and the consensus reads'
callable reference segments become DNA's intervals. Samples run on a pool of
scoped threads (`--threads`). A sample whose traces cannot be read, that DNA
rejects, or on which DNA panics gets a failure row; the others continue. The
truth of the whole panel and the reference load before DNA runs, so a
database that cannot be opened, has another schema, or fails a query, or a
reference that cannot be read, stops the run first. The haplotype measures
compare on that reference sequence; a sample that DNA analyses on another
sequence gets a failure row (kind `reference`).

`--triage` also takes each sample's triage evidence from the same analysis
and runs the [draft triage](#draft-triage); without it the run and its
outputs are unchanged.

## Measures

The truth side is compared as stored in the database; DNA's side uses its
notation. Every sample is compared under six measures: primary, strict,
IUPAC-compatible, the two haplotype measures, and haplotype LoD.

**Primary.** The lab's concordance (`mtdna_raw`, `src/core/comparison.py`,
`SampleComparator.is_concordant`), with DNA as tool A (the automated caller)
and the truth as tool B:

- an allele is `(position, reference, alternate)`. Positions are normalised
  (`16193.0` → `16193`, thousands separators removed), and alleles whose
  reference equals the alternate are dropped;
- `N/A` when either side has no allele;
- insertions and deletions at 309, 455, 463, and 573, and at their `.n`
  insertion positions, are left out ("special" positions); substitutions
  there count;
- `N` when a position both sides call holds different allele sets
  (`comparison.py` also skips 309.1, 309.2, and 573.1, which hold insertions
  and are special already);
- `N` when DNA alone calls an allele inside HV1 (16024–16365), HV2 (73–340),
  or HV3 (438–576);
- `N` when the truth alone calls an allele inside an HV region and inside
  DNA's callable intervals, or DNA reports none;
- `Y` otherwise. Interval and region tests use an allele's anchor base, so
  `16193.1` lies in `[16024, 16193]`.

A truth allele where DNA has no callable base is therefore not counted, and
a truth IUPAC allele differs from any definite base DNA calls.

**Strict.** The primary rules, and a truth-only allele in an HV region also
counts inside the truth's analysed intervals: where DNA has no callable base
it is missed (rule `uncovered`). This makes DNA's coverage gaps visible.

**IUPAC-compatible.** The primary rules, and DNA's base matches a truth IUPAC
code that admits it, like `iupac_bases_compatible`: two alleles with the same
position and reference match when their called bases share a possible base
(`16189C` and `16189Y`), and where DNA calls nothing inside its callable
intervals its base is the reference, which matches a truth IUPAC allele whose
code admits the reference (`16189Y` over a reference `T`).

**Haplotype and haplotype IUPAC-compatible.** The primary and the
IUPAC-compatible rules, applied only to alleles that do not describe the same
sequence on both sides, so a haplotype written or placed differently is not a
difference (`16190.1C` and `16193.1C` in the C run 16190–16193, `248DEL` and
`249DEL` over `AA`, `514DEL 515DEL` and `523DEL 524DEL` in the AC repeat,
`16182.1C 16189C 16193DEL` and `16183C 16184A 16189C`):

- the comparable intervals are the HV regions within DNA's callable
  intervals (merged where they touch), or the HV regions when DNA reports
  none, and within the truth's analysed intervals when it records any: a
  stretch the reviewers did not analyse holds no truth to compare. There is
  no circular seam inside them;
- a DNA allele outside the truth's analysed intervals raises no `dna_only`
  (the primary rules, the lab's own comparator, still count it);
- in each, both sides' alleles are applied to the reference: an allele
  replaces as many bases as its reference has (`514`, `CA`, `-` deletes 514
  and 515; equally many alternate bases substitute them one by one),
  insertions follow their anchor base in ordinal order, and IUPAC codes stay
  as written. Only the reference's length is used, not its bases;
- the interval is walked base by base and cut wherever both local sequences
  have the same length, so the sequences are equal exactly when every piece
  is. A piece that differs is a *block*. A block that reaches the end of its
  interval with unequal lengths is *open*: the length difference may
  continue where DNA has no callable base;
- the haplotype measure compares symbols exactly. The haplotype
  IUPAC-compatible measure also lets a symbol that a truth allele wrote match
  a DNA symbol that shares a base with it (`16191Y` against DNA's `T`, or
  against the reference base where DNA calls nothing), as two alleles do
  under the IUPAC-compatible rules. A symbol that comes from the reference on
  the truth side matches only itself, so a DNA code (`Y`, `N`) where the
  truth calls nothing stays a difference;
- alleles whose bases lie in the interval outside every block describe the
  same sequence and leave the comparison; the rules decide the rest, special
  positions still left out of them, and `N/A` is decided on all alleles.
  An interval where one side has two alleles at one base or insertion
  ordinal, or an allele that cannot be applied (an insertion without an
  ordinal, an ordinal on a reference base, an inserted `-`, a reference or
  alternate that is empty or not ASCII, or bases past the interval), cannot
  be built and stays under the rules.

Special-position indels are left out of the rules, as in `comparison.py`, so
their lengths must not decide a difference in the sequences either, yet they
can write a placement (`455.1T 459DEL` equals `456T`; `573.1C` equals
`570.1C`). Each interval is therefore built three ways when the sides'
special-position indels differ: each side with its own, both with DNA's, and
both with the truth's. One way decides the whole interval: the first that
leaves no block over any allele the rules count, or else the first with the
fewest such blocks (a count of blocks, not of unexplained alleles, so a tie
falls to the earlier way; this errs toward discordance). An allele leaves the comparison when the deciding way
finds no block over its bases, so two alleles of one interval are never
explained by two different ways. With their own indels, `455.1T 459DEL` against `456T`
is one sequence; with the same indels, `309.1C 315.1C` against
`309.1C 309.2C 310.1C` is (only the excluded 309.2C length differs). A
difference that remains under every way stays a difference.

Only representation therefore changes a status: a sample is concordant on a
haplotype measure exactly when, apart from special-position indel lengths,
every allele the base rules count describes the same sequence as the other
side's, or the base rules find it concordant anyway. A sample concordant
under the base rules stays concordant.

**Haplotype LoD.** The haplotype measure, except at a truth heteroplasmy
the traces cannot show (`lod.rs`). A truth substitution to a two-base IUPAC
code is *observed* when the covering admitted reads show both of its bases
there, each as some read's call or as some read's second channel (a channel
other than its own call) at 0.10 or more of its normalised A/C/G/T profile
(rounded to thousandths, so 0.0995 counts), or when one of its bases shows
as such a second channel; otherwise it is *below the limit of detection*.
Every observation of an admitted read counts, masked ones included, so the
limit errs toward counting a heteroplasmy as observed. Below the limit the allele compares as under the haplotype
IUPAC-compatible measure: DNA's base there, or the reference where DNA calls
nothing but reads the position, matches when it is one of the code's bases,
in the alleles and in the sequences. At or above the limit the code must
match exactly. A minor reference base where DNA calls the alternate shows the
heteroplasmy as well. The limit is 10 % (user, 2026-10-16; see the
[approval protocol](../proposals/0005-evaluation-feedback-loop.md)); each
sample's count of truth IUPAC positions below it is reported (`below_lod`).

**DNA's callable intervals.** The comparable span is DNA's callable
intervals, which leave out unresolved consensus bases (`N`; SRS-ALN-015), so
neither DNA's sequence nor the IUPAC-compatible rules' "where DNA calls
nothing" ever stands for a position DNA left unresolved.

For each sample the harness also reports, on the primary measure, the alleles
unique to each side (special positions left out), every rule that applies,
the positions those rules name, whether every such position holds a truth
IUPAC allele (`iupac_only`), and their windows: `hv2_polyc` 302–316,
`hv1_polyc` 16180–16195, `hv3_polyc` 560–576, `ac_repeat` 513–526, `other`.
It lists the `uncovered` truth alleles as well.

## Discrepancy classes

Each sample takes one class, which names the part of DNA (or the lab) that
owns its difference from the truth (ADR-0077). A sample that failed under a
measure is `failed`, and one not applicable under a measure (a side without an
allele) is `not_comparable`; any other takes the first of these whose measure
is discordant:

| Class | Measure | Difference | Owner |
|---|---|---|---|
| `failed` | any | DNA failed, or the truth has no record | core or the input |
| `not_comparable` | any | DNA or the truth has no allele | the input, or coverage |
| `sequence` | haplotype IUPAC-compatible | the sequences differ | core, or the truth (adjudication) |
| `heteroplasmy` | haplotype | only at a truth IUPAC code | core minor signal (deferred) |
| `unanalysed` | primary | every primary difference is a DNA call outside the truth's analysed intervals | a review where the truth has no analysis |
| `representation` | primary | the same sequence written differently | `dna-normalize`, or a truth minority form |
| `coverage` | strict | a truth allele where DNA has no callable base | core callability, or the lab's assay |
| `none` | — | every measure concordant | — |

`class_regions` lists the windows of the discordant positions of the class's
measure. Whether a representation is DNA's to fix or a minority form in the
truth is measured on the corpus with `dna-eval convention`, not per sample.
The class follows from the statuses alone, so a baseline run's classes come
from its `scores.tsv` when it reports all four measures.

## Discrepancy causes

Each sample with a discrepancy also has a cause: who accounts for it, so the
difference reaches the one who can act on it. Every discordant position of
the sample's class gives a finding, `cause:evidence@position` (or `@window`)
in `cause_evidence`; the sample's `cause` is the first of `dna`, `lab`,
`technology`, and `undetermined` among them. Each finding also has a scope
(`cause_scope` for the sample): `defect`, or `deferred` for an item left
undone by design, such as the IUPAC codes DNA does not write; a defect comes
before a deferred item. The causes come from the
admitted reads' evidence at each region position (with the core's own trust
verdict on each call, `SampleAnalysis::trust`), the limit-of-detection
signals, and how the reviewers write each sequence over the whole truth
database; nothing is listed by hand, and the truth is only read.

| Class | Evidence | Cause |
|---|---|---|
| `coverage` (each uncovered truth allele) | no admitted read covers the position (`no_read`) | `lab`: a missing or short trace |
| | a truth IUPAC code of two bases that a read shows, judged as under `heteroplasmy` (`iupac_not_reported`, `below_lod`) | `dna` deferred, or `technology` |
| | a read covers it with a call masked for phase loss: dephased or mixed signal (a mixed segment counts as phase loss), or a segment after a repeat run (`phase_loss`) | `technology` |
| | a read covers it with a call masked for weak signal or irregular spacing (`weak_signal`) | `lab`: the run's quality |
| | clean reads tie and the consensus writes `N` (`contested_reads`) | `undetermined`: for the reviewer |
| | it lies in a stretch of runs the consensus left undecided, no read showing both run ends (`run_unread`) | `technology` |
| | a read's call there is `N` (`unresolved_call`) | `dna`: its callability policy |
| | a clean call the core does not trust, withheld by a call masked for phase loss: dephased or mixed signal, or a segment DNA starts after a repeat run (`beside_phase_loss`) | `technology` |
| | a weak-signal veto denies a call that reads the reference (`veto_out_of_scope`) | `dna`: the veto guards variant support |
| | a call there shows more than one qualifying channel (`mixed_call`) | `undetermined` |
| | every untrusted call is withheld by its read's own end, its trim, or a call masked for weak signal (`read_edge`), or carries a weak-signal veto (`weak_call`); an ADR-0076 exemption counts as trust | `lab`: read quality |
| | otherwise (`untrusted_call`) | `dna`: its callability policy |
| `unanalysed` (each window) | DNA calls where the reviewers did not analyse (`outside_truth_intervals`) | `undetermined`: only a review there verifies them |
| `representation` (each window) | DNA writes the reviewers' majority form of the window sequence (`truth_minority_form`) | `lab`: a minority form in the truth |
| | the truth writes it (`dna_minority_form`) | `dna`: its notation |
| | outside the windows, the majority of truth samples carrying each side's calls (`…_proxy`) | as above |
| `heteroplasmy` (each truth IUPAC code) | a read shows the code's second base at the limit of detection (`iupac_not_reported`) | `dna`, deferred: DNA writes no IUPAC code by design |
| | a read shows it only at 3 % or more (`below_lod`) | `technology` |
| | clean reads on both strands show no second base at 3 % (`no_second_base`) | `lab`: the truth, for a re-check |
| | one strand only (`one_strand`), or a code of more than two bases (`iupac_code`) | `undetermined` |
| `sequence` (each SNV or deletion) | a DNA-only call where clean reads on both strands show only the reference (`phantom`) | `dna` |
| | a truth-only call there (`reads_show_reference`) | `lab`: the truth, for a re-check |
| `failed` | no record in the truth, an unreadable trace, or every trace rejected | `lab` |
| | a DNA error | `dna` |

Any other discrepancy is `undetermined` (`adjudicate`, `not_comparable`).
The written forms group the truth's IUPAC-free samples by the sequence of
each repeat window, so DNA's form and the truth's are compared for the same
sequence; outside the windows no sequence groups them, and the count of
samples carrying each side's calls is a proxy. Each run writes this index as
`forms.tsv`, keyed by the truth database's SHA-256; `--forms` reads an
earlier run's instead of indexing the whole database, and a file written for
another database stops the run. A sample can need more than
one party: the summary counts the samples per class and cause, deferred DNA
items apart, and per cause, scope, and evidence the samples with such a
finding. Insertions are not judged under `sequence` yet: reads showing the
reference at an anchor do not exclude one.

## Estimates

For each weighting (`corpus` = `weight_corpus`, `pool` = `weight`), each
measure with the metrics `not_concordant` (anything but `Y`, failures
included) and `discordant`, and `failed` (measure `all`), over samples of
positive weight that are not retired:

- rate `p = Σ w·y / Σ w`, per stratum and over all;
- stratum size `N_h = Σ_h w`, `W_h = N_h / N`;
- variance `Σ_h W_h² · p_h(1 − p_h) / n_h · max(0, 1 − n_h / N_h)`;
- 95 % interval `p ± 1.96·√variance`, and the unweighted rate.

## Anchors and baseline

Statuses fall into classes so that runs of other conventions stay comparable:
`exact` and `concordant` are good; `error` and `discordant` are bad;
`not_applicable` and `failed` are their own classes.

- Anchor candidates are the manifest's `anchor_candidate = 1` samples, listed
  under each measure. With a baseline, the anchor set is the candidates whose
  baseline `status` was good.
- Changes against the baseline, per measure the baseline reports:
  `regression` (good before, not now), `fix` (the reverse), `changed` (another
  change of class), `added`, `removed`. The summary also tabulates the
  baseline `status` against the primary status.

## Draft triage

`run --triage` measures whether the triage of
[PROP-0004](../proposals/0004-review-triage.md) is worth building (phase 0).
Its detectors are research drafts over the evidence DNA already produces:
nothing in DNA reads them, they decide no call, and none is approved.

**Evidence.** Per sample, from the analysis the run already made
(`analyze_sanger_outcome`), kept even when the call then fails:

- the reportable regions of the target profile (`SampleAnalysis::regions`);
- per region position, its clean and masked reads from one source: at a
  consensus site, the site's supporting and opposing reads observe it cleanly
  and its uninformative reads are masked; elsewhere, a read observes it
  cleanly inside its callable reference segments and masked in the rest of its
  reference segments. Clean reads are counted per strand. Sites name reads by
  trace name; a name that two admitted reads share counts a strand only when
  the site names every read of that name or all of them share the strand,
  otherwise as masked;
- the marks of the consensus sites over each position;
- per region position, each read's *minor channel* where that read calls a
  canonical base without a mask (every covered position, through
  `SampleAnalysis::observations`, read once over the span of the regions):
  the strongest channel of its reference-oriented
  A/C/G/T profile other than its own call, kept when its weight is at least
  0.05; the pair (call, minor base); whether the minor base is the base the
  read calls just before or just after in sequencing order (a neighbour's
  shadow); how systematic each of those two shadows is in the read here,
  the share of its calls within 10 on either side (the call itself left
  out) whose minor channel, at 0.05 or more, is the base of the call before
  (or after) them; and the read's local noise, the median minor weight of
  its calls within 10 calls on either side. Weights and shares are kept in
  thousandths, rounded down for shares;
- per input trace: its name, the primer in the name (a `HV<n><F|R>` token),
  whether it was admitted (with its strand) or rejected (too short, too few
  callable calls, or unplaced), and its callable span in calls.

When every read is rejected, the evidence is the traces alone, with no region
evidence. When DNA cannot give an analysed sample's observations, the
stage-1 evidence is kept without minor channels: the sample's call and scores
stand, a policy with a minor-signal rule sends its would-pass regions to
review because the gate could not run, and `triage.tsv` gives the cause
(`minor_unavailable=…`).

**Detectors.**

| Level | Finding | Fires when |
|---|---|---|
| trace | `trace.too_short`, `trace.weak`, `trace.unplaced` | the read was rejected for that cause (SRS-SAMPLE-030, SRS-SAMPLE-027, SRS-SAMPLE-029) |
| trace | `trace.short` | the callable span has fewer than `short_calls` calls (default 100; `0` turns it off) |
| region | `region.unreliable` | a position is `n` (a site without a call or with `N`), `contested`, `masked` (a site with an uninformative read), or `phase_loss` or `reference_frame` (a run of the site has that length evidence); `masked_read` (a covering read reads the position only through a masked call) is measured but not applied by default |
| region | `region.one_strand` | a position is observed cleanly by one strand only, where the policy requires both strands: everywhere (the draft), nowhere (`off`), or only inside the repeat windows ±3 (`outside_repeats`, R6m) and also in or next to a reference homopolymer of at least `k` bases (`outside_repeats_h<k>`, R6m-h4 with 4); the windows are the measures' windows (`compare::Region::WINDOWS`) and the runs are read from the reference |
| region | `region.minor_signal_single_strand` | a position one strand alone reads cleanly shows a minor signal under the policy's single-strand rule (off by default); it applies only where the policy lets one strand pass |
| region | `region.gap` | a position is observed cleanly by no read |
| region | `region.minor_signal` | a position's reads show a second allele under one of the policy's minor rules (default `upstream@100`; off turns it off) |

**Minor rules.** A rule is a kind and a threshold `t` in thousandths; a read
counts when its minor weight (its second channel's share of the read's
normalised A/C/G/T profile, which sums to 1, rounded to thousandths, so 0.0995
counts at `t` = 100) is at least `t` (`>=`) and it passes the kind's shadow
test. The shadow tests compare the minor base with the read's *adjacent
calls* in reference order, not with the reference bases. The cross-strand
kinds mark a position when a forward and a reverse read that count show the
same unordered pair:

| Kind | A read counts when | Marks when |
|---|---|---|
| `strands` | always | a forward and a reverse read show the same pair |
| `unshadowed` | its minor base is neither neighbour's base | as `strands` |
| `upstream` | its minor base is not the base of its 3′-side adjacent call in reference order (the call after it on a forward read, the call before it in sequencing order on a reverse read) | as `strands` |
| `noise` | its minor weight is at least 3 times its local noise | as `strands` |
| `local_shadow` | its minor base is not a neighbour's base whose shadow at least 0.05 of the read's nearby calls show | as `strands` |
| `any` | its minor base is neither neighbour's base | any read counts |
| `lone` | as `unshadowed` | as `unshadowed`, or any read counts where one strand alone reads the position cleanly |

A minor signal never asks for a rerun.

**Blind spot by construction.** A cross-strand kind marks a position only when
a forward and a reverse read count, so none of them can fire where one strand
alone reads cleanly. `upstream` can never mark a heteroplasmy whose alternate
base equals the 3′ neighbour (next reference base), and `unshadowed` and `any`
also miss one equal to the 5′ neighbour. In the merged reviewer truth of about
103,000 samples, 9,368 reference-plus-one-alternate IUPAC sites split as
alternate = previous reference base only 15.5 %, next only 11.3 %, both
10.4 %, neither 62.8 %. `upstream` is therefore blind to about 21.7 % of
corpus heteroplasmy sites (common examples: 146Y, 16189Y, 16183M) and
`unshadowed` to about 37.2 %. The leak bound of a policy with these rules
holds only for the classes they can see; [PROP-0004](../proposals/0004-review-triage.md)
gives the dev-panel split.

A trace with any trace finding is *poor*; it stays a diagnostic unless a
region asks for it. A sample without an admitted read is `rerun` with every
poor trace, the clearest rerun case. A sample DNA cannot analyse (a trace that
cannot be read, an invalid input) or whose consensus cannot be called is its
own class, `failed`; after a call failure its trace findings are still
reported.

**Outcomes.** A region is `rerun` when it has a gap, or a one-strand position
while the policy applies `region.one_strand`, and a poor trace could have read
it; `review` when any applied finding remains; `pass` otherwise. A sample takes
its worst region (rerun > review > pass). Which poor trace could have read a
region is a proxy:

- `primer` (default): the poor trace's primer reads the region on a strand
  the region lacks. A primer's strand is its `F` or `R`; its regions are those
  that at least half of the admitted traces with that primer in the run's
  active samples observe cleanly at 50 positions or more (its footprint). A
  sample's own traces are not left out: a poor trace never contributes to a
  footprint, and a sample's admitted traces are a negligible share of a
  primer's;
- `any`: any poor trace of the sample.

**Measurement.** Retired samples are left out; samples without a comparison
count in A/B/C but not in leaks or catches.

- **A/B/C:** the share of samples in each outcome, unweighted and
  corpus-weighted with the stratified estimator (with its interval); and the
  share of regions of decided samples with region evidence, unweighted and
  weighted by their sample's `weight_corpus` (a point estimate).
- **Leaks:** passed samples that are not concordant under each measure, with
  the one-sided 95 % Clopper–Pearson upper bound on the leak rate among passed
  samples (`1 − 0.05^(1/n)`, about `3/n`, at zero leaks) and the
  corpus-weighted rate. Region leaks are discordant positions inside passed
  regions, each with its kind (`dna_unique`, `truth_unique`, `both`,
  `uncovered`), allele class (`iupac`, `insertion`, `deletion`,
  `substitution`), window, and the detectors that mark the position itself.
- **Catches** per region detector: not-concordant samples it flags (catches),
  concordant samples it flags (burden), and primary discordant positions it
  marks itself (hits). For a detector the policy applies, also the catches no
  other applied detector makes (unique) and the concordant samples only it
  holds back (sole burden: clearing it would pass them). For a detector the
  policy does not apply, instead what applying it would add: the
  not-concordant samples it would newly catch and the concordant samples it
  would newly hold back.
- **Sweep:** the draft and policies that each change one thing
  (`short_calls` 0, 50, 150, 200, 300; the `any` proxy; `region.one_strand`
  off; `masked` off; `masked_read` on; `n` and `contested` only; no unreliable
  cause; one-strand and `masked` off; the gap alone; the minor signal off),
  then the draft with each minor rule but `lone` (which equals `unshadowed`
  while `region.one_strand` applies) at thresholds 60, 80, 100, 120, 150, 200,
  and 250, leaving out the draft's own rule, `lone` at the same
  thresholds with `region.one_strand` off, and the draft's rule joined by
  `local_shadow` at the same thresholds;
  each with A/B/C and its leaks under every measure: a risk–coverage table.
  Every policy runs on the evidence already taken, without DNA.
- **Minor signal:** for the draft and each policy with another minor rule,
  against the same policy without the rule, decided like any policy:
  corpus-weighted A without and
  with the rule and its cost; the compared samples the rule alone moves from
  pass to review and the concordant ones among them (burden); passed samples
  and leaks under every measure without and with the rule; and the recall of
  truth IUPAC alleles (primary discordant positions of class `iupac`) that
  the rule marks itself, among those inside a region that passes without the
  rule and among those inside any region.

Phase-0 results are in
[PROP-0004](../proposals/0004-review-triage.md#phase-0-results).

## Outputs

| File | Content |
|---|---|
| `scores.tsv` | per sample, sorted by `sample_id`: panel and truth batch, stratum, weights, anchor and retired flags; primary status (`concordant`, `discordant`, `not_applicable`, `failed`) and code (`Y`, `N`, `N/A`), allele counts, unique alleles of each side, rules, positions, `iupac_only`, windows; strict status, code, and `uncovered` alleles; IUPAC-compatible status and code; haplotype, haplotype IUPAC-compatible, and haplotype LoD status, code, and differing blocks of the way that decides each interval that hold a discordant position (`<start>-<end>[+]:<DNA>/<truth>`, `+` when open, `-` for an empty sequence, a sequence over 24 symbols as its first and last 10 and its length, `CCCTCCCCCC..TTCTGGCCAC(40)`); `below_lod`, the truth IUPAC positions below the limit of detection; the [discrepancy class](#discrepancy-classes) and its windows (`class`, `class_regions`); its [cause](#discrepancy-causes), scope, and findings (`cause`, `cause_scope`, `cause_evidence`); distinct and rejected traces; failure (`<kind>: <message>`) |
| `estimates.tsv` | per weighting, measure, and metric: each stratum, then `all` with the interval |
| `anchors.tsv` | anchor candidates with their status under each measure, their baseline status, and the primary change |
| `diff.tsv` | with a baseline: per measure the baseline reports, the samples whose status class changed, including added and removed ones |
| `truth_changes.tsv` | whenever the panel's truth changed: each sample whose truth digest differs, with its stratum, retired flag, and recorded (`-` when none) and current digests |
| `summary.md` | the run in aggregate: the panel checked against the truth (truth digest state, changed, retired samples, samples without a record, batch differences, status not `OK`), outcomes and estimates under each measure, primary discordance, discrepancy classes (per class and window, retired samples left out, the causes per class and the samples per cause and evidence, and with a baseline the samples whose class changed), coverage gaps, IUPAC alleles, representation (discordant samples a haplotype measure finds concordant, by window, and the blocks of the rest), anchors, and the baseline |
| `triage.tsv` | with `--triage`, per sample under the draft policy: outcome, region outcomes, findings with their position counts (and `minor_unavailable=…` when the minor channels could not be read), poor traces, the traces a rerun asks for, and the status under each measure |
| `triage_regions.tsv` | with `--triage`, per sample and region: outcome, positions per detector, rerun traces, and the discordant positions under each measure |
| `triage_traces.tsv` | with `--triage`, per input trace: primer, state, callable span, findings, and the regions whose rerun asks for it (`no_admitted_read` for a sample without an admitted read) |
| `triage_leaks.tsv` | with `--triage`, every discordant position inside a passed region under each measure: the sample's outcome, region, position, window, kind, class, the alleles of each side, and the detectors marking it |
| `triage_detectors.tsv` | with `--triage`, per region detector: whether the draft applies it, samples, regions, and positions flagged, catches per measure, burden, and hits; unique catches and sole burden for an applied detector, `if_applied_*` catches and burden for one not applied, `-` in the other set |
| `triage_positions.tsv` | with `--triage`, per region position: the analysed samples each detector marks it in |
| `triage_policies.tsv` | with `--triage`, per sample (retired ones included, with `retired`): its stratum and its class under each policy of the sweep, one column per policy named as in `triage_sweep.tsv` (`pass`, `review`, `rerun`, or `failed`) |
| `triage_minor.tsv` | with `--triage`, per policy of the minor-signal measurement: its rule, A without and with it and the cost, moved samples and burden, IUPAC sites and those caught (in regions that pass without the rule, and in any region), and per measure the passed samples and leaks without and with the rule, the upper bound, and the corpus rate |
| `triage_sweep.tsv` | with `--triage`, per policy: its parameters, A/B/C counts (with the `failed` samples) and corpus-weighted shares, and per measure the passed samples, leaks, rate, upper bound, and corpus rate |
| `triage.md` | with `--triage`, the triage in aggregate: A/B/C overall and by stratum, leaks and what the leaked samples differ in, region outcomes (counts and weighted shares) and region leaks, leak classes, detectors, traces and primer footprints, the sweep, and the minor signal |
| `forms.tsv` | how the reviewers write each sequence over the whole truth database, keyed by its SHA-256: per repeat window and window sequence the samples per written form, and the samples per set of calls outside the windows (`-` for an empty one); no sample identifiers |
| `run.json` | the run header: harness version, the source revision and whether tracked files differed from it, recorded when the harness was built, threads, manifest fields and hash, the panel truth digest, the recorded one, the truth state and the number of changed samples, the database's file hash and `info` (provenance only), paths, DNA's configuration, reference, and profile identities, counts per measure, timings, and with `--triage` the draft policy, the number of policies, and the triage time |

Every file except `run.json` depends only on the panel, the truth, DNA's
calls (and with `--triage` its evidence), and the configuration, reference, profile, and baseline the summary
names, so it is byte-identical whatever the thread count and whatever revision
or database file it ran on. Run outputs and
panel manifests hold sample identifiers and stay outside the repository under
the [data policy](../governance/data.md); DNA documents quote aggregates only.
