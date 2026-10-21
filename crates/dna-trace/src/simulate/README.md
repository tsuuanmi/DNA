# simulate

Owns the forward model of the Sanger electropherogram and its renderer
([PROP-0015](../../../../docs/proposals/0015-forward-model-tract-decoder.md)
increment 1). It sits behind the `simulate` feature as test support for the
[simulated bench](../../../../docs/validation/simulated-bench.md). Only the
evaluation harness enables the feature, and no runtime path reaches it.

It ports the research simulator
(`tools/python/research/trace_simulator`, `model.py`):

- `mod.rs`: `Params`, one read's physical parameters. `render` draws a
  mixture of templates (read-order sequences and fractions) as analysed
  channels: Gaussian peaks on one grid, dye offsets and gains, trinucleotide
  context effects, sequence-determined patterns, minor peaks, cross-talk, a
  slowly varying floor, white noise, a noisy read start, and optional spike,
  blob, and noise artefacts. It then calls the result with the oracle-grid
  basecaller (`basecall`).
- `table.rs`: `Table`, the committed parameter table `params.json`. It holds
  per-primer quantiles over fitted reads of each read parameter, and a
  Gaussian copula (the parameters' normal-score correlation), which draws a
  read's parameters together. It holds aggregates only and is built by the
  research script `table.py`.
- `random.rs`: a deterministic generator (xoshiro256** seeded by
  `SplitMix64`, Box-Muller normals) and `hash_normal`. The latter is the
  research simulator's BLAKE2b-hashed standard normal of a sequence
  context, so both simulators carry the same sequence-determined patterns.

It reads no instrument format, reference, or truth: the harness builds the
templates and writes the ABIF.
