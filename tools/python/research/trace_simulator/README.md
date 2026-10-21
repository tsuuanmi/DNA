# Trace simulator (research)

A physical forward model of the Sanger electropherogram, fitted on real reads, and a
simulator that renders synthetic traces with known truth: one template, or a mixture
of templates such as a C-tract length heteroplasmy, with optional artefacts. It gives
DNA's tract decoder (`dna-trace`) an accuracy benchmark that no published method
provides ([research note](../../../../docs/research/trace-simulator/README.md)).

Research only (ADR-0066): nothing in DNA depends on it, it is not run by CI, and the
comparison with simulated truth lives here, never in a DNA crate (ADR-0065, ADR-0075).
Each script is a self-contained `uv` script (PEP 723) that brings its own `numpy` and
`scipy`; the tooling project's lockfile is untouched.

## Model

Every template's fragment ending `n` bases after the primer migrates to one grid
point, so a length mixture shares the peak grid and differs only in the base each
template shows there. Each peak is a Gaussian (fitted exponentially modified
Gaussian tails add nothing) with:

- centre: the grid (spacing and its drift), a per-dye offset, a trinucleotide
  context offset, a sequence-determined pattern smoothed along the template, and a
  per-read residual;
- height: a log envelope along the read, a per-dye gain, a trinucleotide context
  effect, a sequence-determined pattern, and a per-read residual, times the
  template's fraction;
- width: a line in the call index.

The analysed channels add minor peaks on the grid, spectral cross-talk, a slowly
varying floor, white noise and a noisy read start, and are clipped at zero.
Artefacts: `spike` (four narrow channels), `blob` (one wide channel), `noise`
(a burst). A point heteroplasmy is a second template with one base changed.

Calls come from an oracle-grid basecaller: one call per grid point at the strongest
apex within a third of a spacing, its base the strongest channel there, its vendor
quality the real reads' median quality for that primary share. It never merges or
splits calls, unlike the instrument's basecaller in mixed signal.

## Files

| File | Role |
|---|---|
| `abif.py` | Minimal ABIF reader (real traces) and writer (DATA.9-12, FWO_, PLOC, PBAS, PCON) |
| `reference.py` | rCRS, read placement by seed votes, read templates with edits |
| `fit.py` | Fits the model per real read (manifest rows, clean calls only) |
| `aggregate.py` | Parameter distributions per primer, period and instrument; pooled context tables; cross-sample reproducibility |
| `model.py` | The forward model, the renderer and the oracle-grid basecaller; calibration constants |
| `simulate.py` | Renders one trace and its truth from a JSON spec |
| `dnarun.py` | Runs `dna evidence` in a scratch directory and keeps what the checks read |
| `realism.py` | Simulated vs real reads on the same statistics, and the fit's round trip |
| `artefacts.py` | What `dna evidence` makes of each artefact |
| `bench.py` | The tract-decode benchmark grids (`main`, `ablation`, `smoke`, explicit cases) |
| `replay.py` | Real decodes as simulated truth: residual structure, recovery |
| `report.py` | Benchmark tables |
| `table.py` | The aggregate parameter table of the Rust port (`dna-trace` feature `simulate`): per-primer quantiles and a Gaussian copula over fitted reads, no read's own values |
| `port_check.py` | The Rust port against this simulator: the `sim-bench` case grid for `bench.py --cases`, realism on pure reads, and the bench measures side by side ([simulated bench](../../../../docs/validation/simulated-bench.md)) |

## Run

From this directory, with a release `dna` (`DNA_BIN`, default `target/release/dna`):

```bash
uv run --script fit.py ../../../../results/test-panel/dev-v5.tsv OUT/fit.jsonl --samples 300
uv run --script aggregate.py OUT/fit.jsonl OUT/agg
uv run --script realism.py OUT/agg/pooled.json OUT/fit.jsonl OUT/realism
uv run --script bench.py OUT/agg/pooled.json OUT/bench --grid main --reps 10
uv run --script report.py OUT/bench
```

Outputs go outside the repository; quote only aggregates ([data governance](../../../../docs/governance/data.md)).
