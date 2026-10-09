# Research: signal denoising and baseline correction

**Status:** active — exploratory evidence, not production authority.

## Question

Does a separate conditioned projection of the instrument-analyzed `DATA.9`–`DATA.12`
channels improve DNA's re-calling or read callability on real traces, without
removing genuine secondary peaks?

## Why this is research

The analyzed channels are already baseline-corrected, mobility-shifted, and
spectrally deconvolved by the instrument. On the local reviewed corpus the
dominant failure modes were dephased and globally mixed signal, which no filter
removes; the production answer to those is
[read callability](../../design/callability.md). Any transform that could change
a call must preserve the decoded trace and produce a separate processed
projection ([INV-EVID-001/002](../../architecture/invariants/evidence.md),
[ADR-0013](../../decisions/adr/0013-observational-signal-quality.md),
[ADR-0019](../../decisions/adr/0019-scientific-evidence-hierarchy.md)).

## Candidate methods

- peak-preserving Savitzky–Golay smoothing ([Savitzky and Golay 1964](https://doi.org/10.1021/ac60214a047));
- asymmetric baseline correction ([Eilers 2003](https://doi.org/10.1021/ac034173t); [airPLS](https://doi.org/10.1039/b922045c));
- wavelet soft-thresholding ([Donoho 1995](https://doi.org/10.1109/18.382009)).

## Protocol

1. Apply the candidate to a copy of the channels of each local trace with a
   scratch script outside the repository (Python is limited to research by
   [ADR-0066](../../decisions/adr/0066-python-limited-to-repository-tooling.md)).
2. Re-run peak selection and re-calling on the copy; measure vendor agreement,
   two-channel IUPAC count, and the callability segment map against the
   unconditioned run.
3. Measure secondary-peak retention on known length-variant runs and on
   synthetic major/secondary peaks, baseline drift, impulses, compressed peaks,
   homopolymers, and read ends.
4. Compare downstream concordance with reviewer calls (comparison runs outside
   DNA, [ADR-0065](../../decisions/adr/0065-result-comparison-downstream.md)),
   quoting aggregates only ([data governance](../../governance/data.md)).

## Success criterion

No loss of genuine secondary peaks and a measurable concordance gain on the
reviewed corpus. Promotion then needs an ADR, a `signal_processing` consumer of
a versioned `ConditionedTrace` projection, requirements, contracts, tests, and
validation; the callability core needs no change because it consumes locus
evidence.
