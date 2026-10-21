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
  [--reps <n>] [--keep-traces]
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
and ×4 the read's own), and `rho` 0.5, 0.8, 0.9, 1 and 1.25, six reads per
cell. It also has three pure reads per tract and noise level: 396 reads,
which run in a few seconds on 100 threads. Each case's seed is a hash of its
name, so a case renders the same in every run.

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
baseline (`crates/dna-eval/baselines/sim-bench.tsv`) or from
`--baseline`. `bench.tsv` holds the measures and `cases.tsv` every case.
`--reps` runs a larger set for validation, where a cell's 30 reads are too
few to read a percentage near 0 or 100.

## Baseline (increment 1)

The current decoder on the frozen set, at the commit that introduced the
bench (`crates/dna-eval/baselines/sim-bench.tsv`). These are the numbers later increments must
beat.

| Tract | Noise | Structure right | Dominant right | Clear but wrong | `rho_ratio` | Tie → shorter |
|---|---:|---:|---:|---:|---:|---:|
| HV2-F | 1 | 96.7 | 95.8 | 4.2 | 0.963 | 0.0 |
| HV2-F | 2 | 90.0 | 95.8 | 4.2 | 0.961 | 0.0 |
| HV2-F | 4 | 89.3 | 82.6 | 13.0 | 0.907 | 40.0 |
| HV2-R | 1 | 79.3 | 78.3 | 21.7 | 1.012 | 50.0 |
| HV2-R | 2 | 80.0 | 70.8 | 29.2 | 1.027 | 33.3 |
| HV2-R | 4 | 88.9 | 76.2 | 23.8 | 0.988 | 66.7 |
| HV1-F | 1 | 96.7 | 95.8 | 0.0 | 0.965 | 16.7 |
| HV1-F | 2 | 76.7 | 70.8 | 4.2 | 0.971 | 0.0 |
| HV1-F | 4 | 33.3 | 31.6 | 21.1 | 0.910 | 20.0 |
| HV1-R | 1 | 96.7 | 83.3 | 16.7 | 0.965 | 0.0 |
| HV1-R | 2 | 90.0 | 87.5 | 8.3 | 0.986 | 33.3 |
| HV1-R | 4 | 63.3 | 54.2 | 20.8 | 0.977 | 16.7 |

On 7,920 reads (`--reps 120`), pooled over noise:

- Structure is right in 91 % (HV2-F), 86 % (HV2-R), 69 % (HV1-F) and 81 %
  (HV1-R) of decodes.
- With the right structure, the dominant length is right in 93-100 %.
- Near ties the decoder under-reads the shorter neighbour by about 5 %
  (`rho_ratio` 0.95-0.97) on all but HV2-R (1.00).

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
