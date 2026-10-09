# ADR-0067: Signal-derived read callability

- **Status:** Accepted
- **Date:** 2026-10-09
- **Related decisions:** [ADR-0013](0013-observational-signal-quality.md),
  [ADR-0019](0019-scientific-evidence-hierarchy.md),
  [ADR-0021](0021-scientific-core-confidence-floor.md),
  [ADR-0027](0027-mixed-supporting-signal-snv-eligibility.md),
  [ADR-0029](0029-profile-aware-gotoh.md),
  [ADR-0053](0053-public-differential-locus-signal-evidence.md),
  [ADR-0060](0060-separate-variant-canonicalization-nomenclature.md),
  [ADR-0062](0062-read-callability.md),
  [ADR-0064](0064-crate-ready-module-layering.md),
  [ADR-0065](0065-result-comparison-downstream.md)

## Context

The reference-free read path decodes a chromatogram, re-calls bases, annotates
rolling SNR, and trims the two ends with a Tracy-style best-section walk-out.
Nothing in that path measures whether the signal is still one ladder. Four
weaknesses follow:

- the trim threshold is relative to the read's own best section, so a read that
  is bad everywhere keeps most of its calls;
- only the two ends can be removed, although a read can lose and regain phase
  in the middle;
- the sequence-only `post_homopolymer` rule of ADR-0062 excludes calls after
  every long homopolymer whether or not the signal there is still in phase, and
  the `read_end` margin removes calls that trimming has pushed to a read end;
- polymerase slippage in long homopolymers and dinucleotide repeats makes the
  signal downstream a superposition of ladders offset by one or more positions,
  which the caller reads as phantom indels and mixed bases.

Measured outside the repository (ADR-0065) against revision `67a25ae` on the
local reviewed corpus (aggregate numbers only, per [data governance](../../governance/data.md)):

- of 54 reviewed samples, 28 agreed exactly with the reviewer calls;
- seven missed calls sat beside a poly-C run, where both strands had observed
  the call and the eligibility rules had excluded it, including a reverse read
  with no double peak after the run;
- one globally mixed read with 42 % two-channel calls was retained to 86 % by
  trimming, aligned at identity 0.81 with 24 gap opens, and produced 13 false
  indels; four samples of the same read class failed the identity gate and
  aborted;
- directly measured, clean stretches showed 0–3 double peaks per 40 calls
  against 14–34 after a poly-C run; among those double peaks, the secondary
  channel matched a neighbouring primary channel in 94–100 % of cases after a
  single-offset slippage, 50–74 % in a mixed read, and 24–71 % after a 13-base
  run with several offsets; a forward read re-synchronised about 30 calls after
  a 9-base run; a 7-base run and the HVS-III dinucleotide repeat dephased
  neither strand in a sample where the run was short. Only the first long run
  met in trace order, and everything before it, is reliably in phase;
- the current trimming kept the post-run section of one read and discarded 229
  in-phase calls before the run.

Severity therefore depends on the run the read actually carries, not on a
reference window, and the information needed to decide callability is in the
read's own signal.

## Decision

1. A new layer-2 module `callability` derives, per read and reference-,
   profile-, and filename-free, in trace order for both strands, the
   `dna.read_callability/v1` view: per-position features (dominance, secondary
   ratio, spacing deviation, weakness), phase-state segments from a hysteresis
   state machine over rolling defect fractions, the classification of every
   masked segment (`weak`, `irregular`, and for double peaks `dephased` or
   `mixed` from a shadow model of the read's own calls shifted by up to three
   positions; see the calibration amendment), a per-position mask, and the
   callable span.
2. Repeat runs found in the read's own primary calls lower the onset threshold
   right after the run and attribute a segment (`after_repeat`); they never mask
   by themselves. The run itself stays callable while its signal is one ladder.
3. The core works on plain per-position records and never on Sanger types; one
   adapter builds those records from the chromatogram loci, basecall-independent
   locus evidence, and primary calls. No trait is introduced until a second
   adapter exists (ADR-0058 §4).
4. Decoded channels, loci, calls, selected peaks, and locus evidence are never
   modified; a masked position keeps its call and evidence (INV-EVID-010).
5. Increment 1 publishes the view as observation in `dna.basecalls/v3`,
   `dna.analysis/v9`, and `dna.sample_evidence/v10` and as aggregate log
   metrics; it alters no trim bound, alignment, warning total, or eligibility.
   Thresholds are configuration (schema 7, `[callability]`), calibrated on the
   local corpus before the mask is allowed to act.
6. Increment 2 makes the mask act: the trim interval becomes the callable span
   plus a context margin, masked positions enter alignment as unresolved,
   variant eligibility reports the mask state in place of the sequence-only
   `post_homopolymer` window, sample aggregation treats masked coverage as
   covered but not callable, and a `sample` read with too few callable calls is
   recorded as rejected instead of failing the operation. That increment
   supersedes ADR-0062 in part.
7. Signal denoising remains research ([denoising](../../research/denoising/README.md)).

## Alternatives

- Static per-primer callable ranges, as the legacy pipeline uses: contradicts
  the contract that placement never depends on filenames (ADR-0062).
- Hard trimming after every long run: ADR-0062 measured recall falling from
  0.979 to 0.83, and it discards in-phase calls before the run.
- Deconvolving the shifted ladders to recover the sequence behind the shift:
  rewrites evidence and needs truth-labelled validation (ADR-0013, ADR-0019);
  it stays deferred and is studied in
  [phase recovery](../../research/phase-recovery/README.md).
- Declaring difficult regions in the target profile: a reference window cannot
  see the run length a sample actually carries, and method parameters do not
  belong in profiles (ADR-0063).
- Reusing `dna.windowed_snr/v1` noisy regions as the mask: amplitude SNR is
  blind to the difference between an in-phase ladder and a shifted one, and
  ADR-0013 made those regions observation-only.
- Several retained intervals instead of one interval plus a mask: breaks the
  `retained == primary[trim]` contract and the single-query alignment.

## Consequences

- One more sequential stage (`callability`) between signal processing and
  quality control, with a stable span name and one completion event.
- **Breaking:** configuration schema 7 adds `[callability]` and renames
  `quality_control.trim_window_size` to `penalty_window_size`; the three result
  contracts gain the required callability view.
- The thresholds shipped in `config/dna.toml` are initial values; increment 1
  exists to record their effect on the local corpus before increment 2.
- A slipped population above one half flips the primary inside the run with no
  in-run signature; only cross-strand disagreement at sample scope can reveal
  it ([known limitations](../../validation/known-limitations.md)).

## Calibration amendment (2026-10-09)

Increment 1 labelled a double-peak segment `dephased` when at least 75 % of its
double peaks shared one modal neighbour offset. On the local corpus that rule
called most post-run slippage `mixed`, because slippage is usually two-sided:
shorter and longer length populations leave shadows one call before and one
call after. Before increment 2, still inside the unreleased cycle, the rule is
replaced by a shadow model ([callability method](../../design/callability.md),
substep 4.4):

- each double-peak segment's normalized corrected amplitudes are fitted by
  exact non-negative least squares to one-hot templates of the read's own
  primary calls shifted by −3..+3;
- the segment is `dephased` when the main-ladder share reaches
  `minimum_main_share` (0.35), the share of shadows two or three calls away
  stays at or below `maximum_far_share` (0.12), and a one-call shadow reaches
  `minimum_shadow_share` (0.10); otherwise it is `mixed`;
- dephased segments publish their `shadow_offsets`; weights and shares are
  logged in aggregate form only and are never mixture fractions.

Measured outside the repository on the 199 local traces with the Rust
implementation on top of revision `c91b1bf` (aggregate only):

- segment boundaries, masks, callable spans, exit statuses, and every other
  output field were identical to `c91b1bf` across all 430 results;
- `dephased` segments went from 19 to 245 (46 569 calls; median length 129
  calls) and `mixed` segments from 322 to 96; the 81 fitted mixed segments
  cover 2 083 calls (median length 14);
- dephased segments: main share q10/median/q90 0.52/0.68/0.79, far share
  0.000/0.007/0.065; shadows on both sides in 131 segments, one call before
  only in 65, one call after only in 46;
- of the 81 fitted mixed segments, 42 exceeded the far-share limit, 38 had no
  one-call shadow, and 4 fell below the main-share minimum; 15 were too short
  to fit.

With baseline-corrected event amplitudes, far shadows are close to zero even
after long runs and in the globally mixed read class, so `dephased` means
"explained by shadows of neighbouring calls". It covers slippage and a ladder
smeared by one call alike; the event amplitudes cannot tell a co-located
shadow from a neighbour's peak tail. How severe a dephased segment is, and
whether its sequence can be recovered, is left to the
[phase-recovery research](../../research/phase-recovery/README.md).

## Supersession

ADR-0013's deferred decision on signal-driven exclusion is taken here; its
observational SNR windows stay observation-only. ADR-0027 and ADR-0053 keep
phase recovery, weighting, and consensus as research; detection is promoted.
ADR-0062 remains the authority for `read_end` and the sequence-only
`post_homopolymer` rule until increment 2 lands, when this record supersedes
that rule in part.
