---
implements: [SRS-CALL]
decided-by: [ADR-0065, ADR-0067, ADR-0069, ADR-0076]
covers:
  - crates/dna-sanger/src/callability/
  - crates/dna-sanger/src/callability/sanger.rs
  - crates/dna-sanger/src/model/callability.rs
---

# Callability Method

Part of the canonical [[design/pipeline|DNA pipeline]].

**Method identifier:** `dna.read_callability/v1`

Derives, from the read's own signal in trace order, where the read is still
one ladder and therefore callable. The stage runs after signal processing and
before quality control for every read, reference-, profile-, and filename-free.
It never modifies channels, loci, calls, selected peaks, or locus evidence; a
masked position keeps its original call and evidence. Quality control,
alignment, variant eligibility, and sample aggregation consume the mask
([[ADR-0067]]; see
[Consumers](#consumers)).

## Inputs

One record per call position, in call order:

- the non-negative corrected A/C/G/T amplitudes at the locus event
  (`LocusEvidence.corrected_amplitudes`, [[design/signal-processing/locus-evidence|locus evidence]]);
- the canonical locus sample position, used only for spacing;
- the primary call's channel, `None` for an unresolved `N`.

The core consumes nothing else. The Sanger adapter (`crates/dna-sanger/src/callability/sanger.rs`)
is the only code that reads `Chromatogram`, `BaseCalls`, and `SignalAnalysis`.

## Coordinate domains

Features, segments, the mask, and the callable span use call indexes; sample
indexes enter only through the spacing feature. Both are 0-based and all
intervals are half-open.

### Substep 4.1 — Per-position features

For position `i` with amplitudes `a` and total `T = Σ a`:

- `c1`, `c2`: the channels with the highest and second-highest amplitude; ties
  favour the lower channel index; no primary channel when `T ≤ 0`;
- `dominance = (a[c1] − a[c2]) / T`;
- `secondary_ratio = a[c2] / a[c1]`;
- `spacing_deviation = max(|Δ_{i−1} − M|, |Δ_i − M|) / M`, with `Δ` the sample
  spacing between adjacent calls and `M` the median spacing over the
  `spacing_half_window` spacings on each side of `i` (8 by default; available
  terms only; zero when `M` carries no information);
- `weak`: `T ≤ 0`, or `a[c1]` below `weak_amplitude_fraction` of the read's
  median primary amplitude. The read-wide median is used because the local
  noise estimate of the SNR method is inflated inside homopolymers, which must
  not look weak.

Every value is rounded to six decimals before storage or comparison. Fewer than
two positions is a typed failure.

### Substep 4.2 — Repeat runs

Maximal runs of at least `repeat_min_length` identical
canonical primary calls (homopolymers) or alternating pairs of distinct
canonical calls (dinucleotide repeats) are recorded; an unresolved call breaks a
run. Runs are priors: they never mask by themselves.

### Substep 4.3 — Defects and segmentation

Each position is one defect class, by precedence: `weak`; `double` when
`secondary_ratio ≥ double_peak_ratio`; `spacing` when
`spacing_deviation ≥ spacing_defect`; otherwise none. The `repeat_window_calls`
positions after a repeat run form its prior window; the run's own calls are
never inside one.

The forward statistic at `i` is the defect fraction over `[i, i + window_calls)`.
A hysteresis state machine walks the read:

- it starts masked when the first fraction reaches `onset_defect_fraction`;
- an in-phase stretch ends where the fraction reaches `onset_defect_fraction`,
  or `exit_defect_fraction + 1/window_calls` inside a prior window, localized at
  the first defective position of that window; when that position lies in the
  `repeat_window_calls` positions after a repeat run, the masked stretch starts right
  after the run instead, also masking the calls already walked past, because a
  slipped population can lead there before its second peak becomes visible;
- a masked stretch ends where the fraction drops to `exit_defect_fraction`;
- decisions freeze once fewer than half a window remains, so the tail inherits
  the current state;
- a callable island shorter than one window is absorbed by its neighbours.

Direction-awareness needs no orientation knowledge: "after the run" is later in
trace order on both strands, and only the first long run met in trace order
reads its true length.

### Substep 4.4 — Segment classification

Each masked run is labelled by its dominant defect class (ties resolved in the
order weak, double, spacing): `weak`; for double peaks `dephased` or `mixed`
from the shadow model below; `irregular` for spacing. A segment that starts in
a prior window is marked `after_repeat`. Classification never moves a segment
boundary.

**Shadow model.** Offset `k` means that the shadow at call `i` copies the
primary call at `i + k`. Over the segment's calls, every position that is not
weak and has a resolved primary call contributes its normalized amplitudes
`y(i,c) = a[c] / T` and seven one-hot templates `t_k(i,c) = [primary(i + k) = c]`
for `k ∈ {−3,…,+3}`. A template whose neighbour lies outside the read or is
unresolved is zero; neighbours may lie outside the segment. The fit minimizes
`Σ_i Σ_c (y(i,c) − Σ_k w_k t_k(i,c))²` with `w_k ≥ 0`:

- the 7×7 Gram matrix of the templates (integer counts) and their products with
  `y` are accumulated once;
- every subset of the six shadow templates is solved together with the main
  template `k = 0` by Cholesky factorization, fewest templates first and then
  in the fixed bit order `+1, −1, +2, −2, +3, −3`; a subset is singular when a
  pivot falls to `10⁻⁹` of its largest diagonal entry, and it is rejected when
  any weight is negative;
- the residual `y·y − 2 w·b + wᵀGw` (clamped at zero) decides, and a later
  subset replaces the best one only when it is smaller by more than
  `10⁻⁹ · max(y·y, 1)`; the weights are rounded to six decimals.

A fit needs at least eight usable positions, one more than the parameters, and
a positive weight sum; otherwise the segment is `mixed`. With
`share_k = w_k / Σ w`, the main share is `share_0` and the far share is the sum
over `|k| ≥ 2`. An offset is reported when its share reaches
`minimum_shadow_share`. The segment is `dephased` when the main share reaches
`minimum_main_share`, the far share is at most `maximum_far_share`, and `−1` or
`+1` is reported; otherwise it is `mixed`. Shares are rounded to six decimals
before comparison.

### Substep 4.5 — Mask and callable span

Every position of a non-`in_phase` segment is masked. The callable span runs
from the first to the last unmasked position and is empty when every position
is masked. Segments partition `[0, call_count)`; the module checks this and the
mask once and fails typed on an inconsistency.

## Projections

- `read.callability` in `dna.basecalls/v3` and `dna.analysis/v9`, and
  `reads[].callability` in `dna.sample_evidence/v10`: the callable span, the
  ordered segments with state, `after_repeat`, and the ascending
  `shadow_offsets` of dephased segments, and the masked-call count.
  Per-position features, shadow weights and shares, and the mask are not
  serialized.
- `callability_completed`: aggregate counts, the callable span, and one bounded
  token per segment, for example `282..445:dephased(-1,+1)[0.640,0.012]+repeat`
  or `106..275:mixed[0.855,0.131]`: the reported offsets of a dephased segment
  and, for every fitted double-peak segment, its main and far shares.

## Consumers

The core caller does not read callability directly. The Sanger adapter
(`read_processing::evidence`) turns every masked call into a `ReadEvidence` mask
([[ADR-0069]]):
- the mask is *anchoring* for a `dephased` call and *unresolved* otherwise;
- its reason is `post_homopolymer` when the segment starts in the window after
  a repeat run, and otherwise the state label;
- every call within `sanger_evidence.signal_edge_margin` calls of an
  uninformative call or a read end is at a signal edge (SRS-CALL-013);
- before the callable span, a masked call with one qualifying channel and a
  primary peak at least `leading_height_fraction` (1.0) of the read's median
  primary peak is unmasked and
  vouched for by its base only (SRS-CALL-014). A read's start is noisy, and
  its masked stretch only ends where a whole window is clean, so such calls
  are masked with their noisy neighbours; on development panel v5 they read
  the reviewers' base at 1,650 of 1,650 truth variants. A call masked after a
  repeat run (`post_homopolymer`) is a dephased tail, not start noise, and
  stays masked: no base comes from a dephased tail (on development panel v5
  that changed no output);
- at an end of the informative interval, a scored call at a signal edge whose
  own signal is unambiguous has its verdict on its base taken off the edge
  (SRS-CALL-017), within `edge_identity_margin` calls of that end: one
  qualifying channel, a primary peak at least `edge_height_fraction` (1.0) of
  the read's median, a second-channel share below `edge_second_share` (0.10),
  a peak no broader than `edge_width_factor` (1.5) times the read's median
  width (a dye blob is broader), a spacing within `edge_spacing_minimum` to
  `edge_spacing_maximum` (0.5 to 1.5) times the median, no phase
  loss beyond that end, and not the base of its neighbour towards that end (a
  signal one position off). Its indel verdict keeps the edge and its bound is
  unchanged, so it still counts towards runs: as a base-only call it would
  stop counting, and run lengths the consensus estimated from it (`315.1C`)
  were lost. On development panel v5 against the rule's absence: `315.1C`
  called in 126 more samples and none lost, `73G` in 56, `strict` 33 better
  and 1 worse, samples with an unreviewed difference unchanged (4.83 %), and
  every new false call (`310C` twice, `309T`, `16183DEL`, one `73G`) inside a
  check item. The same-base condition gives up 60 correct `16189C` calls that
  otherwise came with one-position shift errors (`310C`, `460C`, `567C`,
  `574C`);
- at the read's start, an edge call whose own signal is unambiguous against
  the start's own norms has its verdict on its base taken off the edge as well
  (SRS-CALL-018). The start's calls sit behind the start noise the leading
  rescue leaves masked, so they lie at no end of the informative interval,
  and a read's signal rises over its first calls, so they stand below the
  read's median: within `start_identity_margin` calls of the callable span's
  start, the rule takes the median primary peak, width, and spacing of the
  first `start_norm_calls` (30) in-phase calls of the callable span and asks
  for `start_height_fraction` (0.5) of that peak and a second-channel share
  below `start_second_share` (0.15). A call reading the base of the call
  before it qualifies only behind a vouched call, judged in read order, so
  the runs of a read's start (`TT` at 16024-16025, `TTT` at 16027-16029)
  chain from their first call. Masked calls before the callable span may lie
  between the call and the start, but no masked call of the span may. On development panel v5:
  positions to review 3.20 to 3.07 %, clear samples 40.56 to 43.79 % (45.11
  to 48.22 % on batches from 2026), samples with an unreviewed difference
  unchanged, `73G` in 29 more samples, HVS-I start window `16024-16033`
  checked in 17 % fewer samples, `haplotype_lod` none better or worse;
- every call carries its second base and its share (SRS-CALL-015).

- **Quality control** trims to the callable span plus up to
  `quality_control.trim_context_margin` calls of an adjacent dephased segment
  ([[design/quality-control|quality control]], substep 5.3). A read with fewer than
  `minimum_callable_calls` unmasked calls fails `analyze` and `basecall`
  typed after the `callability_completed` event; `sample` and `consensus`
  record it as rejected and continue ([[design/sample-evidence/read-registry|read registry]]).
- **Alignment** keeps dephased calls with their call and evidence profile, which
  still read the main ladder and anchor the alignment, and presents every other
  masked call as unresolved `N` without a profile ([[design/alignment|alignment]]).
- **Variant eligibility** turns a masked evidence call into the reason of its
  segment, and a column on a masked call is no SNV candidate
  ([[design/variant-calling|variant calling]]).
- **Sample aggregation** records a masked call as a `masked` observation that
  never retains a locus by itself and carries no nucleotide mass
  ([sample evidence](sample-evidence/README.md)).

## Configuration

`[callability]` holds `window_calls`, `onset_defect_fraction`,
`exit_defect_fraction`, `minimum_main_share`, `maximum_far_share`,
`minimum_shadow_share`, `weak_amplitude_fraction`, `repeat_min_length`,
`minimum_callable_calls`, `repeat_window_calls`, `double_peak_ratio`,
`spacing_defect`, and `spacing_half_window`
([[reference/configuration|configuration]]). `repeat_min_length` sets the
shortest repeat run, and `minimum_callable_calls` the fewest unmasked calls a
read needs. The shadow model's offsets (±1..±3), its minimum fit size (one more
than its seven parameters), and the solver tolerances are the model's
definition and method constants of `dna.read_callability/v1`.

## Interpretation limits

Phase states are read-level signal states, not Phred values, error
probabilities, mixture fractions, or artifact classes. Shadow weights and shares
are model coefficients, not population fractions. A `dephased` segment says
that shadows of neighbouring calls explain its double peaks; event amplitudes
cannot tell a co-located slippage shadow from a neighbour's peak tail, so a
ladder smeared by one call is `dephased` as well. It does not recover the
sequence behind the shift. The templates come from the read's own calls, so a
primary flipped by a strong shadow becomes a wrong template for its neighbours.
Slippage in a dinucleotide repeat leaves shadows two calls away, which the
model counts as far, so such a segment is `mixed`. A `mixed` segment says only
that the shadow model does not explain its double peaks. A slipped population
above one half flips the primary inside the run with no in-run signature; only
cross-strand disagreement at sample scope can reveal it.

## Validation

Unit tests cover the core on plain records (ranking, spacing, weakness, run
detection, onset/exit hysteresis, prior lowering, island absorption, tail
freezing, shadow-model recovery of synthetic mixtures, non-negativity,
coincident templates, determinism, classification thresholds, unchanged
boundaries under any classification rule, repeat attribution, partition
invariants). Synthetic ABIF fixtures provide one- and two-sided shadow ladders
behind a long homopolymer, incoherent double peaks, and an amplitude-collapsed
tail; the CLI
tests pin the published view and the completion event. Real-trace evidence is
produced outside the repository and quoted in aggregate (ADR-0065,
[[governance/data|data governance]]).

## Non-goals

Phase recovery or shadow-ladder deconvolution, calibrated weighting from
callability features, and signal denoising remain research
([phase recovery](../research/phase-recovery/README.md),
[denoising](../research/denoising/README.md)).
