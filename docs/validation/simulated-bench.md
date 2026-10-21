---
implements: [SRS-CALL]
decided-by: [ADR-0074, ADR-0075, ADR-0084]
covers: [crates/dna-eval/src/simbench.rs, crates/dna-eval/baselines/, crates/dna-trace/src/simulate/, crates/dna-trace/tables/]
---

# Simulated tract-decode bench

Part of the [[validation/evaluation-harness|evaluation harness]]. It measures
the tract decoder ([[design/length-mixture]], SRS-CALL-025) against a known
truth, and gates every decoder change of [[PROP-0015]] beside the
[[validation/artefact-bench|artefact bench]].

```bash
cargo run --release -p dna-eval -- sim-bench \
  --config config/dna.toml --reference references/rCRS.fasta \
  --threads <1-256> --output <new-directory> [--baseline <bench.tsv>] \
  [--reps <n>] [--keep-traces] [--tolerance <points>]
```

## What it measures

Each case renders one read through a C-tract as a mixture of two molecules:
the dominant length and its shorter neighbour at `rho` times its share. The
mixture is the truth.

- **Rendering.** `dna-trace`'s forward-model renderer (feature `simulate`,
  test support only, [[ADR-0084]]) does the rendering, with read parameters
  drawn from an aggregate table.
- **Decoding.** The case is written as ABIF and runs through DNA's read
  evidence stage in process, exactly as `dna evidence` runs it.
- **Scoring.** The harness scores the decode of the tract against the truth
  on the marginal that carries the varying run: the total for HVS-II, the
  first run for HVS-I.

The scoring lives in the harness, never in a DNA crate.

| Tract | Primer | Molecules | Expected structure |
|---|---|---|---|
| HV2-F | HV2F | `C^a T C^6` (315.1C), `a` 9 or 8 | mark `T` |
| HV2-R | HV3R | the same, read in reverse | mark `A` |
| HV1-F | HV1F | `C^L` with T16189C, `L` 11 or 10 | one run, no mark |
| HV1-R | HV1R | the same, read in reverse | one run, no mark |

The frozen case set (version 1) has four tracts, three noise levels (×1, ×2
and ×4 the read's own), and `rho` 0.5, 0.8, 0.9, 1 and 1.25, plus pure reads.
It comes in two sizes:

- **The gating set** (the default): 60 reads per cell and 30 pure reads per
  tract and noise level, 3,960 reads in about 15 s on 100 threads. Decoder
  increments decide on it; [[PROP-0015]] states the resolution each
  threshold gets.
- **The smoke set** (`--reps 6`): 396 reads in a few seconds, a quick
  check.

Each case's seed is a hash of its name, so a case renders the same in every
run, and the smoke set's cases are the first of the gating set's.

Measures, per tract and noise level (percentages):

| Measure | Meaning |
|---|---|
| `decoded` | reads given a decode at the tract |
| `structure_right` | the chosen mark is the tract's |
| `dominant_right` | the dominant length is right (`rho` 1 excluded) |
| `right_given_structure` | the same, among decodes with the right structure |
| `clear_wrong` | a `clear` margin on a wrong dominant |
| `rho_ratio` | median decoded `rho` over the true one, at `rho` 0.5-0.9 |
| `tie_shorter` | decodes at `rho` 1 whose dominant is the shorter length |
| `pure_decoded` | pure reads (no neighbour) given a decode |

`bench.md` reports the measures with their change from the committed
baseline of the set's size (`crates/dna-eval/baselines/sim-bench.tsv` for the
gating set, `sim-bench-smoke.tsv` for the smoke set) or from `--baseline`.
The baseline reads strictly: a cell is a finite number or `-`, and a tract
and noise level appear once. Against it the bench gives a verdict: `SIM-BENCH
FAIL`, with exit status 2, when `decoded`, `structure_right`,
`dominant_right`, or `right_given_structure` of a cell falls, or
`clear_wrong` rises, by more than `--tolerance` points (default 0: the
frozen cases render the same reads every run), when a cell is missing on
either side, or when a read has no read evidence; otherwise `SIM-BENCH
PASS`. A validation set without a committed baseline is reported, not
gated.
`bench.tsv` holds the measures and `cases.tsv` every case. Any other
`--reps` is a validation run without a committed baseline.

## Baseline (increment 1)

The current decoder on the gating set, at the commit that introduced the
bench (`crates/dna-eval/baselines/sim-bench.tsv`). These are the numbers
later increments must beat.

| Tract | Noise | Structure right | Dominant right | Right given structure | Clear but wrong | `rho_ratio` | Tie → shorter |
|---|---:|---:|---:|---:|---:|---:|---:|
| HV2-F | 1 | 97.7 | 96.2 | 99.1 | 3.8 | 0.953 | 16.7 |
| HV2-F | 2 | 95.3 | 95.4 | 98.7 | 4.6 | 0.954 | 15.0 |
| HV2-F | 4 | 76.1 | 65.5 | 83.7 | 32.3 | 0.943 | 20.7 |
| HV2-R | 1 | 86.4 | 85.9 | 99.5 | 14.1 | 1.012 | 69.5 |
| HV2-R | 2 | 84.3 | 80.9 | 96.9 | 19.1 | 1.008 | 54.4 |
| HV2-R | 4 | 91.6 | 81.8 | 88.5 | 15.8 | 0.994 | 63.5 |
| HV1-F | 1 | 95.3 | 95.0 | 100.0 | 3.3 | 0.950 | 13.3 |
| HV1-F | 2 | 76.5 | 75.0 | 100.0 | 11.9 | 0.955 | 24.6 |
| HV1-F | 4 | 36.2 | 35.5 | 98.6 | 16.7 | 0.940 | 17.6 |
| HV1-R | 1 | 97.3 | 93.8 | 96.2 | 5.4 | 0.964 | 21.7 |
| HV1-R | 2 | 89.9 | 87.3 | 97.2 | 8.4 | 0.965 | 21.7 |
| HV1-R | 4 | 57.3 | 54.3 | 96.9 | 16.5 | 0.975 | 20.7 |

- Structure choice is the main failure. HVS-I gains a spurious mark as the
  noise rises, and HVS-II reverse loses its `A` mark in 8-16 %.
- With the right structure, the dominant length is right in 84-100 %.
- Clear but wrong reaches 14-19 % on HVS-II reverse at the reads' own noise
  and twice it, and up to 32 % at four times it.
- Near ties the decoder under-reads the shorter neighbour by about 5 %
  (`rho_ratio` 0.94-0.98) on all but HVS-II reverse (0.99-1.01), and exact
  ties go to the longer length except on HVS-II reverse.

These are the research benchmark's findings.

## Agreement with the research simulator

The renderer ports the research simulator
(`tools/python/research/trace_simulator`), including its BLAKE2b-hashed
sequence-determined patterns. Those patterns are fixed per sequence context,
so each tract sees one realisation of them. Measured with
`tools/python/research/trace_simulator/port_check.py`:

- **Realism** (pure reads at the reads' own noise, 240 per simulator, the
  research fitter and `dna evidence` on both). Medians, Rust against the
  research simulator: sigma at call 100 −0.4 %, spacing −0.1 %, clean-call
  second/primary q50 +1.3 % and q90 −2.3 %, `dna` second share q50 +2.3 %,
  q90 +0.2 % and q99 +0.2 %. Per primer the second-share medians differ by
  up to 12 % on 60 reads each.
- **Decoder** (the bench's grid at 60 and 120 reads per cell; pooled over
  noise, the largest difference over tracts in points). The comparison is
  against a second research run with other seeds:

| Measure | Rust vs research | Research vs research |
|---|---:|---:|
| decoded | 2.0 | 2.7 |
| structure right | 2.5 | 2.7 |
| dominant right | 1.4 | 2.6 |
| right given structure | 1.0 | 1.0 |
| clear but wrong | 1.0 | 2.6 |
| `rho_ratio` (×100) | 0.5 | 0.5 |
| tie → shorter | 6.7 | 9.9 |
| pure decoded | 3.3 | 3.3 |

Every measure agrees within three points, or within the research
simulator's own run-to-run spread for the two measures read on few
decodes.

## Data

The renderer reads read parameters from
`crates/dna-trace/tables/simulator-params.json`. It holds aggregates only:

- per primer, quantiles over the fitted reads (levels 0.025 to 0.975) of
  each read parameter;
- a Gaussian copula: the correlation of the parameters' normal scores, so
  that, for example, a read with a low envelope also has a low floor;
- the pooled trinucleotide context tables and the calibration constants.

It holds no read's own values and no trace. It is rebuilt from the research
fit (`fit.py`, development panel v5, 300 samples, seed 20261021) by
`table.py`. Rebuild it, re-run the agreement check, and raise the case set
version when the chemistry, polymer, instrument, basecaller, or primers
change, as for the consensus noise τ.

## Limits

- The bench decodes one read per case. The decoder works per read, and any
  aggregation over a sample's reads (any number per strand) is the
  consensus's ([[ADR-0074]]), measured on real panels.
- The simulated signal is cleaner than real dephased signal: real decode
  residuals are more autocorrelated (research note). Real error rates are
  therefore at least the bench's.
- HVS-III 568-573 is not in the bench (a non-goal of [[PROP-0015]]).
