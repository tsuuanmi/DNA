---
id: PROP-0015
status: accepted
created: 2026-10-21
summary: "the tract decode reads each call's channel amplitudes through the read's own peak model, renders its design from a forward model of the electropherogram (context heights, dye gains), chooses a tract's structure by a penalised criterion, calibrates its margins by parametric bootstrap, and is gated by a frozen simulated benchmark; HVS-II forward reads give run 1 from the total marginal, and a 'shortest dominant' convention is an optional labelled policy on bias-corrected shares."
evidence: [RES-length-heteroplasmy]
---

# Proposal: A forward-model tract decoder

## Problem

The tract decode ([[PROP-0012]] increment 2, [[ADR-0084]],
[[design/length-mixture]]) fits each tract as a mixture of molecules over
the read's per-call profiles: one sample per channel at the call's position,
divided by the row sum, against one-hot rows. Until now its accuracy had no
outside measure, and no published method estimates a C-tract length
distribution from the phase shift after the tract (RES-length-heteroplasmy).

The trace simulator (RES-trace-simulator, draft PR #219) supplies that
measure. It is a physical model of the analysed channels fitted on 1,169
development-panel v5 reads, and it renders reads with known length mixtures
that run through `dna evidence`. On 5,059 simulated reads (revision
`3cad8903`):

- **The core is accurate.** With the right structure, the dominant length is
  right in 96-100 % of decodes at the reads' own noise, for HVS-I and HVS-II
  in both directions, down to a neighbour ratio
  `ρ = share(L−1) / share(L)` of 0.9.
- **Structure choice is the main failure.**
  - The lowest-residual choice has no complexity penalty. It adds a spurious
    mark to a single run (HVS-I with T16189C) in 6 % of forward decodes at
    the reads' own noise, and in 22 % and 65 % at two and four times that
    noise. It drops HVS-II reverse's `A` mark in 13 %.
  - Real reads do the same at HVS-II: 7 of 87 HV2F decodes and 3 of 65 HV3R
    decodes choose another structure.
  - With a spurious mark, the first-run marginal is spread over lengths 1-10,
    yet a third of these marginals are labelled `clear`. The total marginal
    stays right (34 of 34 at the reads' own noise).
- **Margins are over-confident under noise.** At four times the reads' own
  noise, 14-36 % of decodes are `clear` but wrong.
- **Near ties are biased.** The decoder under-reads the minor neighbour by
  about 5 % (ρ̂/ρ 0.94-0.96).
  - The decoder reads one sample at the call position, which sits on the
    strongest peak. A minor molecule's peak in another dye sits off it, by
    its dye offset and its sequence-dependent mobility: about 1 sample
    against a peak sigma of 2.8. With ideal positions, ρ̂/ρ is 1.00.
  - Context-dependent peak heights then tilt exact ties: to the longer length
    in 80-89 % on HVS-II forward and HVS-I, to the shorter in 75 % on HVS-II
    reverse. A signal share is not a molecule share.
  - Forward and reverse reads of one near-tie at HVS-II therefore disagree by
    construction. This is consistent with the "truth = mode − 1" reviewer
    pattern of the 2026-10-19 decoder spike.
  - The cleaned refit (signature and leakage) does not remove the bias.
- **The first run is weakly identified on HVS-II forward reads.** Its
  marginal is 2-4 times less precise than the total's: ρ̂ spread near a tie
  0.077 vs 0.025.

The simulator is optimistic. Real decode residuals are more autocorrelated
(0.45-0.65 vs 0.33-0.42), and real neighbouring peaks are misplaced together
more, so real error rates are at least those above.

## Goals

- Shares that estimate molecule fractions, with no strand- or
  tract-dependent bias near ties.
- A structure (single run, or a mark) chosen with a penalty, and no
  published first-run marginal for a mark the signal does not support.
- Margins whose `clear` state is calibrated: rarely wrong at any noise the
  reads show.
- A frozen simulated benchmark that gates every decoder change, beside the
  adjudicated artefact bench ([[validation/artefact-bench]]).
- On HVS-II forward reads, a run-1 view from the total marginal where run 2
  is known.
- An optional, labelled "shortest dominant" convention, kept out of the
  decoder and applied to bias-corrected shares.

## Non-goals

- **Deciding runs from one strand.** The two-strand rule
  (SRS-CONS-013, [[ADR-0074]]) stays the only decision.
- **Encoding reviewer habit in the decoder.** The convention is a separate,
  optional policy, and the user decides whether to enable it.
- **Heteroplasmy levels.** Shares stay signal shares (SRS-VAR-009); after
  this proposal they are better estimates of molecule fractions, still not
  levels.
- **Re-calling the whole read** (PROP-0012's stage 2). Only flagged tract
  windows change.
- **HVS-III 568-573.** It lies at a read's start or late end, where the
  simulator is least realistic. It is measured, not decided.

## Evidence and research

- **RES-trace-simulator** (`docs/research/trace-simulator`, draft PR #219):
  - the forward model and its fitted parameters;
  - realism checks of simulated against real reads;
  - the tract-decode benchmark, the ablation that attributes the near-tie
    bias, and the replay of real decodes.
- RES-length-heteroplasmy: field practice ("shortest dominant" when lengths
  are equally dominant, with no numeric cut), and the per-batch noise
  question.
- [[PROP-0012]] and [[design/length-mixture]]: the current decoder and its
  marginals.
- [[ADR-0074]] and SRS-CONS-013: the two-strand rule, its `direct` and
  `derived` views, and its noise τ.

Fitted facts the design relies on (median over reads):

- Peaks in the analysed channels are Gaussian. Sigma is 0.23-0.25 of the
  spacing (2.8 samples at 12 samples per call); exponentially modified
  tails add nothing.
- Dye offsets are A −0.23, C −0.05, G +0.24 and T +0.05 samples. Peak
  centres scatter by 1.0 sample (SD) around a local line.
- A trinucleotide table explains 56 % of the log-height variance, for
  example G after A −0.55 log. Of the rest, 77 % of the height variance and
  85 % of the position variance recur at one reference position across
  samples.
- Instruments and periods are indistinguishable; primers differ.

## Proposed design

### Ownership

| Crate | Owns |
|---|---|
| `dna-trace` | The per-read peak model, amplitude read-out, forward-model design, structure criterion, and bootstrap margins. Also the simulator's renderer, as test support only (`#[cfg(test)]` or a dev-only feature), never a runtime path |
| `dna-sanger` | Hands `dna-trace` the analysed channels, the call positions, and the callability onsets (it hands profiles today). Owns the `[trace]` configuration and the height-context table's provenance |
| `dna-kernel` | `ReadEvidence::tracts` (`TractLengths`), revised in place |
| `dna-consensus` | The new run-1 view for HVS-II forward reads (`length::views`) |
| `dna-normalize` | The optional convention policy |
| `dna-eval` | `tracts.tsv` columns, and reporting the simulated bench beside `dna-eval bench` |

`dna-trace` keeps its boundary ([[ADR-0084]]): plain data in, no instrument
format, reference, alignment, or truth. The Python simulator stays research
(ADR-0066). No required step runs Python.

### 1. Per-read peak model and amplitude read-out

`dna-trace` fits the read's own peak model on its in-phase calls before the
tract (the window the signature reads today):

- sigma as a line in the call index;
- a per-dye position offset;
- the local spacing.

A call's four channel amplitudes then come from a least-squares fit of
Gaussians at the call's expected positions: each dye at the grid position
plus its offset, over the call's window, with the neighbouring calls' peaks
in the design. The decode reads these amplitudes in place of one sample per
channel.

SRS-CALL-025's "analysed channels at its calls' locus samples" becomes
"amplitudes read through the read's peak model". In the ablation, this alone
moves ρ̂/ρ from 0.95 to 1.00.

### 2. Decode design from the forward model

A molecule's column is the signal it is expected to show, not a one-hot row.
At each call, its base's expected amplitude is:

- the read's dye gain;
- times a context factor (previous, own, and next base);
- times the envelope, which cancels in the row normalisation.

The context table is chemistry-level data (BigDye v3.1, POP7, 3500-series).
It is shipped as versioned data with its provenance and the procedure to
refit it, and re-estimated under the same triggers as the consensus noise τ:
chemistry, polymer, instrument, basecaller, or primers.

The downstream dynamic programme uses the same expected profiles. The
existing signature and leakage refit becomes the no-context special case and
is retired once the forward design passes the bench.

### 3. Structure by a penalised criterion

Each candidate structure (single run, or a mark of each other base) is
scored by a penalised criterion. The baseline is BIC on the per-call
residual with the read's own noise and the candidate's column count. A mark
must also be supported by its own call: the mark base's share at the mark's
call exceeds the read's in-phase noise quantile.

- An unsupported mark is never chosen.
- When the chosen structure is a single run, or a mark is not supported, the
  first-run marginal is not published. `first` becomes optional in the
  contract, and the total marginal carries a single run's length.
- Every candidate's score is published, as the residuals are today.

### 4. Margins by parametric bootstrap

For a decode whose deflated margin falls in a configured grey zone:

- simulate `B` replicate windows from the fitted mixture, the read's peak
  model, and its residual, resampled in blocks so that its autocorrelation
  is kept;
- refit each replicate;
- `p_dominant` is the share of replicates whose dominant length is the
  fitted one.

A marginal is `clear` only when `p_dominant ≥ 1 − α`. Decodes outside the
grey zone keep the analytic margin, which bounds the cost. The decoder's
cost is already 90 % of the read stage, so this matters.

`B`, `α` and the zone are configuration with provenance. If the bootstrap is
not adopted, `trace.near_equal_z` is retuned on the bench instead, to the
same acceptance criterion.

### 5. HVS-II forward reads: run 1 from the total

On a forward read whose tract's mark is the anchor, run 1 also gets a
`derived` view: total − 1 − run 2. It is used when run 2's length is known
on that sample, because the consensus reads it in phase or has decided it.
Its standard error is the total marginal's.

The forward strand's run-1 view is the more precise of this view and the
`direct` view. The rule itself (SRS-CONS-013) is unchanged. A view is still
one of each strand's inputs to `(p − ½) − max(d, τ) ≥ k · SE`.

### 6. Optional "shortest dominant" convention

A labelled `dna-normalize` policy, off by default. On a run whose decodes
are `clear` or `near_equal` on both strands with the same top two lengths
`L` and `L − 1`, it reports `L − 1` when the bias-corrected
`share(L−1) ≥ ρ · share(L)` on both strands. The output is labelled with the
policy and ρ.

- ρ is configuration with provenance. It must not be set closer to 1 than
  the single-read resolution the bench measures, about 0.85-0.9 today.
- Before item 2 lands, "bias-corrected" means divided by the bench's
  measured strand- and tract-specific ρ̂/ρ.
- ρ is never fitted to reviewer records. It is validated against an
  orthogonal truth (fragment sizing, RES-length-heteroplasmy) when one
  exists.
- The user decides whether to enable it.

### Contracts touched

- **`ReadEvidence::tracts`** (`TractLengths`, `dna.read_evidence/v1`), revised
  in place before the first release, as for [[ADR-0082]]:
  - `first` becomes optional;
  - each candidate carries its penalised score and whether its mark is
    supported;
  - each marginal gains `p_dominant` (when bootstrapped) beside `z` and
    `z_adjusted`;
  - the amplitude read-out and design method are named in the decode;
  - `cleaned` is retired with item 2.
- `dna-consensus` length views: the new `derived` run-1 view
  (SRS-CONS-013).
- The review (SRS-API-020) reads an absent `first` as no direct view.
- `dna-eval`'s `tracts.tsv`: the new columns.
- `[trace]` configuration: the peak-model window, the structure penalty and
  mark-support quantile, the bootstrap's `B`, `α`, zone and seed, and the
  context table's path.
- `dna-normalize` configuration: the policy's switch and ρ.
- The `sanger` plugin version rises with every output-changing increment;
  the `consensus` version with increment 6.

### Increments

Each increment is one PR. Each is measured on the frozen simulated bench, on
`dna-eval bench` (PASS), and on development panel v5 per period (current
batches first) with `dna-eval diff`. Real-data acceptance uses aggregates
only.

1. **Bench.**
   - Port the simulator's renderer to Rust as `dna-trace` test support,
     with the pooled parameter table as committed data: aggregates only,
     no sample data.
   - Freeze a case set: tract × direction × ρ × noise × seed, about 400
     reads.
   - Record the current decoder's numbers as the baseline.
   - *Acceptance:*
     - the Rust renderer matches the research simulator on the realism
       statistics (sigma, spacing, second-share quantiles within 5 %);
     - the current decoder's bench numbers are within 3 points of the
       research benchmark;
     - no DNA output changes.
2. **Structure.** Penalised criterion, mark support, optional `first`.
   - *Acceptance (bench):* wrong structure ≤ 1 % at the reads' own noise and
     ≤ 5 % at twice it (from 6-13 % and 11-22 %).
   - *Acceptance (real reads):* on development panel v5, the anchor mark is
     chosen in ≥ 98 % of HVS-II tract decodes (from 92-95 %).
   - *Acceptance (both):* no first-run marginal is published with an
     unsupported mark; unflagged non-het blocks do not rise in any period.
3. **Amplitude read-out.** The per-read peak model.
   - *Acceptance (bench):* ρ̂/ρ within 0.98-1.02 on HVS-I and HVS-II in both
     directions at the reads' own noise; dominant accuracy given the right
     structure is not lower.
   - *Acceptance (panel):* unflagged non-het blocks do not rise.
4. **Forward-model design.** Context heights and dye gains; retire the
   cleaned refit.
   - *Acceptance (bench):* exact ties go to the shorter length in 40-60 %
     on every tract and direction; ρ̂/ρ within 0.97-1.03.
   - *Acceptance (panel):* forward/reverse disagreement on development-panel
     near ties falls, reported by period.
5. **Margins.** Bootstrap in the grey zone, or retuned `near_equal_z`.
   - *Acceptance (bench):* `clear` but wrong ≤ 1 % at the reads' own noise
     and at twice it, ≤ 5 % at four times it; the `clear` coverage is
     reported.
   - *Acceptance (real reads):* two-strand agreement among `clear` decodes
     is ≥ 99 % on development panel v5.
   - *Acceptance (run context):* where the run context of [[PROP-0014]] is
     available, the per-batch noise does not change the verdict on more
     than a reported share of runs.
   - *Cost:* read-stage time rises by at most 25 %.
6. **HVS-II forward run 1 from the total.**
   - *Acceptance:* on development panel v5, decided 303-309 runs rise, with
     accuracy per bucket not below the existing two-strand bucket.
     Unflagged non-het blocks do not rise in any period. The decision
     buckets are reported as in [[PROP-0012]].
7. **Convention policy.** Off by default.
   - *Acceptance:* its effect is reported against the "mode − 1" reviewer
     cases, never fitted to them, and ρ is bounded as in item 6 of the
     design. Enabling it is the user's decision.

## Alternatives considered

- **Keep one-hot profiles and divide by measured bias factors.** It is
  cheap, but the factors depend on the strand, the tract, and the
  surrounding sequence. It corrects the average and not the case, and it
  hides the cause. It is only the interim correction for the convention
  (item 6) until increment 4.
- **Re-call the whole read through the forward model** (full trace
  deconvolution before basecalling). Every call would move at once. This
  stays PROP-0012's stage 2.
- **Decode the raw channels (DATA.1-4).** The instrument's mobility and
  spectral correction would have to be re-modelled. The analysed channels
  are what the fitted model describes.
- **Commit simulated traces, or Python-generated vectors, as test data.**
  This would need megabytes of channel data, and a CI step would depend on
  Python (ADR-0066). A Rust renderer with an aggregate parameter table is
  smaller and deterministic.
- **Retune `near_equal_z` only.** It is cheap, and it is the fallback of
  item 5. But the replay shows real misfit is structured, so one global
  threshold either over-flags clean reads or under-flags noisy ones.
- **A full Bayesian posterior over lengths (MCMC).** It is better
  calibrated in principle, but the cost per tract is far above the
  bootstrap restricted to a grey zone.
- **Put the convention in the decoder.** Rejected: the decoder must not
  encode a reporting habit ([[PROP-0012]]).

## Validation plan

- **The simulated bench** (increment 1) gates every increment on:
  - structure accuracy;
  - dominant accuracy given the structure;
  - `clear` but wrong;
  - ρ̂/ρ;
  - exact-tie direction;

  each at one, two, and four times the reads' own noise, per tract and
  direction. The bench is synthetic: it bounds the decoder's behaviour, it
  does not replace real-data gates.
- **The artefact bench** (`dna-eval bench`, [[validation/artefact-bench]])
  passes before any panel run.
- **Development panel v5:** `dna-eval diff` per period, current batches
  first.
  - The primary gate is unflagged non-het error blocks per zone kind. They
    must never rise on any period.
  - Decided runs are reported by bucket, as in [[PROP-0012]].
- **Run context.** Where [[PROP-0014]] (`dna.run_context/v1`) is
  merged, margins and the two-strand rule's τ are checked with per-batch
  noise.
- **The current-period test set** is used for aggregates only, at logged
  uses, once per accepted increment. A result that drives a change turns it
  into development data.
- **Unit tests** in `dna-trace` on rendered windows:
  - one component;
  - a near-equal pair;
  - a single run with no mark;
  - a mark with no support;
  - a leakage tail with no mixture;
  - each artefact class next to a tract.
- **Realism.** Before increment 5 is accepted, the simulator's known gaps
  (correlated mobility, residual autocorrelation) are closed or bounded, and
  the bench thresholds are re-checked.

## Rollout and rollback

One PR per increment, each raising the `sanger` plugin version (and the
`consensus` version for increment 6). Increments 1-5 change evidence only,
and no decision, except through the existing two-strand rule's inputs. A
revert rolls an increment back. The convention is off by default.

## Observability and operations

- The `tracts_decoded` log event gains counts of:
  - structures chosen by the penalty over the residual;
  - unsupported marks;
  - bootstrapped decodes, and time spent.
- `dna-eval` reports the bench beside `bench`. The replay script stays in
  research.

## Security and privacy

No new input format or runtime dependency. The committed parameter table
holds aggregates over reads, not sample data. Simulated traces carry no
sample information. Research outputs stay in scratch and `results/`.

## Compatibility and migration

`dna.read_evidence/v1`, `dna.consensus/v1`, and `dna.review/v1` are revised
in place before the first release. SRS-CALL-025, SRS-CONS-013, and
SRS-API-020 are updated with their increments. [[design/length-mixture]]
records the forward model when increment 3 lands. This proposal extends
[[PROP-0012]]'s increments 2 and 4; PROP-0012 stays their record.

## Decision

Pending.

## Implementation

None yet.
