---
id: RES-trace-simulator
status: active
summary: A physical forward model of the electropherogram, fitted on real reads, and a simulator that benchmarks the tract decoder against known length mixtures.
---

# Research: trace simulator and tract-decode benchmark

**Status:** active. This is exploratory evidence, not production authority. The code is in
[`tools/python/research/trace_simulator`](../../../tools/python/research/trace_simulator/README.md)
(research Python, [[ADR-0066]]). Nothing in DNA depends on it, and the comparison with simulated
truth stays there ([[ADR-0065]]).

## Question

How accurate is the tract decode ([[design/length-mixture]], [[ADR-0084]]) when the truth is known?
No published method estimates a C-tract length distribution from the phase shift after the tract,
so the decoder has no outside benchmark. Reviewer records are a convention at near-equal mixtures:
"shortest of equally dominant".

## Method

- **Forward model.** Every template's fragment of a given length migrates to one grid point, so a
  length mixture shares the peak grid. Each peak has:
  - a Gaussian shape;
  - a centre from the grid, a dye offset, a trinucleotide context, and a sequence-determined
    pattern;
  - a height from an envelope, a dye gain, a context, and a pattern, times the template's fraction.

  The background is minor peaks on the grid, cross-talk, a slow floor, white noise, and a noisy read
  start. Calls come from an oracle-grid basecaller. The output is a minimal ABIF that `dna` reads
  like an instrument file.
- **Fit.** 1,169 traces of 300 development-panel v5 samples (HV1F, HV1R, HV2F, HV3R), clean calls
  only.
- **Checks and benchmark.**
  - Realism: real reads against simulated twins, measured with the same code.
  - Artefacts: spike, dye blob, noise burst, point heteroplasmy.
  - Benchmark: 5,059 simulated length mixtures through `dna evidence` (revision `3cad8903`).
  - Ablation: the same mixtures with model parts switched off.
  - Replay: real decodes rendered back as simulated truth.

## Findings (aggregates)

- **Peaks.** Peaks in the analysed data are Gaussian: sigma is 0.23-0.25 of the spacing, and
  exponentially modified tails add nothing. Spacing is 12.0 samples per call.
- **Sequence determines most of the variation.**
  - The trinucleotide context explains 56 % of the peak-height variance (G after A: -0.55 log).
  - 77 % of the remaining height variance and 85 % of the remaining position variance recur at one
    reference position across samples.
- **Groups.** Instruments and periods are indistinguishable; primers differ.
- **Realism.**
  - Simulated reads match real ones on widths, spacing, height spread, floor, second-peak
    quantiles, and withheld calls. On simulated clean templates, `dna` reproduces the template
    (5e-4 calls differ).
  - Gaps: real neighbouring peaks are misplaced together (position residual lag-1 0.46 vs 0.22),
    and real decode residuals are more autocorrelated (0.45-0.65 vs 0.33-0.42). The benchmark is
    therefore optimistic.
- **Accuracy.** At real noise, with the right structure, the decoder gets the dominant length right
  in 96-100 % of decodes for HV1 and HV2 in both directions, down to rho = share(L-1)/share(L) =
  0.9.
- **Structure choice is the main failure.**
  - The lowest-residual choice adds a spurious mark to a single run (HV1, T16189C) in 6 % of
    forward decodes, and drops HV2 reverse's mark in 13 %. Real reads do the same at HV2 in 5-8 %.
  - A spurious mark makes the first-run marginal meaningless, yet a third of those marginals are
    still labelled `clear`.
- **Margins** are over-confident under noise: 14-36 % of decodes are `clear` but wrong at four times
  the real noise.
- **Bias near ties.**
  - Decoded shares under-read the minor neighbour by about 5 %. The decoder reads one sample at the
    call, and a minor molecule's peak in another dye sits off it (ablation: 0.95 -> 1.00 with ideal
    positions).
  - Context-dependent heights tilt exact ties to the longer length on HV2 forward and HV1 (80-89 %),
    but to the shorter on HV2 reverse (75 %).
  - Strands therefore disagree on near-ties, and DNA's forward mode sits at the length after the
    "shortest dominant" convention's choice.

## Candidate next steps

None is adopted; each needs a proposal.

1. Read each call's channel amplitudes by fitting the read's own peak shape, not one sample.
2. Render the decode's design from the forward model (context heights, dye gains).
3. Choose the structure by a penalised criterion, and publish no first-run marginal for an
   unsupported mark.
4. Derive HV2 run 1 from the total marginal on forward reads too.
5. Calibrate margins by parametric bootstrap from the fitted mixture.
6. Keep a frozen simulated case set as the decoder's regression benchmark.
7. Close the realism gaps, and validate shares against fragment sizing, before any rho policy is set.
