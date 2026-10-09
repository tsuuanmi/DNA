# Research: phase-shift classification and recovery behind homopolymer runs

**Status:** active. This note is exploratory evidence and has no production authority.

## Question

After a long poly-C run (HVS-I 16184–16193, HVS-II 303–315) or the HVS-III AC repeat,
a Sanger read becomes a superposition of the same sequence shifted by one or more
calls. These are the length populations created by slippage. This note asks two
questions:

1. Can the read's own signal tell such a *dephased* stretch apart from a genuinely
   *mixed* one (mixed template, failed read)?
2. Can the sequence behind the shift be recovered without reference knowledge?

[Read callability](../../design/callability.md) already detects and masks these
stretches, but it does not recover them ([ADR-0067](../../decisions/adr/0067-signal-derived-read-callability.md)).

## Prior art

### Tracy `decompose`

Tracy's `decompose` command is described in [Rausch et al. 2020](https://doi.org/10.1186/s12864-020-6635-8),
with source in `src/decompose.h` and `src/indigo.h`. It targets heterozygous indels:

- **Reference.** It always needs a reference: an indexed genome, a FASTA of at most
  50 kbp, or a wild-type trace. There is no reference-free mode.
- **Breakpoint.** It finds a single breakpoint. Per call it takes "best − second"
  profile probability, then looks for the largest difference between the means of two
  25-call windows. A difference below 0.25 means no shift. Both the window and the
  threshold are hard-coded.
- **Offset.** It finds a single constant offset. For each candidate offset it counts
  the positions after the breakpoint where the shifted reference base is *not* among
  the called peaks. Candidates are minima under a median − 5·MAD threshold, and the
  smallest deletion candidate wins.
- **Alleles.** Allele 1 is the reference-consistent peak and allele 2 is "the other
  called peak" at each position. Peak heights are never used to separate the alleles.
- **Allele fraction.** It reports one global allele fraction per read, from a
  least-squares fit of normalised channel intensities to one-hot templates of
  allele 1 and allele 2. The fit is a brute-force grid with step 0.01.
- **Thresholds.** A minor peak must reach `pratio` (default 0.33) of the major peak to
  be called. This makes a population below about 25 % invisible.

These assumptions do not match mtDNA length heteroplasmy:

- **Several populations.** A tract usually holds several length populations at once
  (n, n+1, n+2 cytosines), which means several offsets. Tracy collapses them into one
  mixed "allele 2".
- **More than one breakpoint.** A read may re-synchronise after a short dephased
  stretch, or meet several runs.
- **No repeat handling.** No Tracy code mentions homopolymers or repeats.
- **No reference allowed.** DNA's read-level stages are reference-free.

The downstream `mtdna_raw` pipeline runs `decompose` on every read but consumes only
allele 1 (`ref1align`/`alt1align`). In practice it uses Tracy for basecalling and
alignment, not for decomposition.

Two ideas carry over:

- The global allele-fraction fit generalises to **K shifted templates of the read's
  own calls**.
- The test "is the secondary peak explained by a shifted copy?" can be made
  **reference-free** by shifting the read's own calls instead of the reference.

### Literature

- **Indelligent** ([Dmitriev & Rakitov 2008](https://pmc.ncbi.nlm.nih.gov/articles/PMC2429969/))
  decodes the two allelic sequences of a heterozygous-indel trace from peak-call letters
  by dynamic programming, with no reference. On 104 traces it reconstructed 99.1 % of
  bases. It requires similar alleles, a fragment much longer than the indel, and
  well-spaced indels.
- **Other tools.** Poly Peak Parser and its sangerseqR implementation, CHILD (rare-variant
  indels with a ratio estimate), MSR, and ShiftDetector are compared in
  [a published tool comparison](https://pmc.ncbi.nlm.nih.gov/articles/PMC3385616/). All of them
  target one indel between two alleles.
- **Length heteroplasmy is usually multi-population.**
  - Blood cells carry mixtures in 303–315 such as 7C+8C, 8C+9C, and up to three
    populations ([a blood-cell study](https://www.sciencedirect.com/science/article/pii/S0006497118677134)).
  - Sequence 3′ of the tract is ambiguous in Sanger reads.
  - Forensic practice avoids reading through the tract, for example with
    junction-placed primers ([HV1/HV2 sequencing with length heteroplasmy](https://www.sciencedirect.com/science/article/abs/pii/S0531513102005460)),
    fragment analysis, or massively parallel sequencing.
  - Length calls remain the least reproducible part of mtDNA typing across technologies
    and software ([EDNAP collaborative study](https://www.sciencedirect.com/science/article/pii/S1872497326001742)).
- **No published Sanger deconvolution for poly-C length mixtures was found.**

## Shadow model

The model uses the [callability](../../design/callability.md) offset convention: offset
`k` means that the shadow at call `i` copies the primary call at `i + k`.

**Fit.** For a stretch of calls:

- normalise the corrected amplitudes per call: `y(i,c) = a_c(i) / Σ_c a_c(i)`;
- fit `y(i,c) ≈ Σ_k w_k · [primary(i + k) = c]` for `k ∈ {−3,…,+3}` with `w_k ≥ 0`.

`w_0` is the main ladder. `|k| = 1` are the near shadows of the n±1 populations.
`|k| ≥ 2` are far shadows, which mostly absorb whatever no near shadow explains.

**Peeling.** Recovery iterates
`S(i) = argmax_c [ y(i,c) − Σ_{k≠0} w_k · [S(i + k) = c] ]`,
starting from the primary calls.

## Measurements

These measurements were made outside the repository against revision `c91b1bf`, on the
local exploratory corpus of 199 traces from 52 samples, and are quoted in aggregate only
([data governance](../../governance/data.md)). Method: raw analysed amplitudes at the
basecaller peak position ±2 samples, argmax primaries, and the segments of at least 16
calls published by increment 1. Because the method differs from production, the numbers
are indicative and are not production calibration.

Weights per increment-1 segment label (q10 / median / q90):

| increment-1 label | segments | main `w_0` | near (`w_−1 + w_+1`) | far (`Σ_{\|k\|≥2}`) |
|---|---:|---|---|---|
| `in_phase` | 252 | 0.58 / 0.92 / 0.94 | 0.03 / 0.05 / 0.31 | 0.02 / 0.02 / 0.09 |
| `dephased` | 5 | 0.58 / 0.64 / 0.65 | 0.17 / 0.28 / 0.32 | 0.01 / 0.03 / 0.13 |
| `dephased`, after a run | 1 | 0.60 | 0.33 | 0.04 |
| `mixed` | 184 | 0.28 / 0.41 / 0.58 | 0.16 / 0.30 / 0.36 | 0.07 / 0.14 / 0.21 |
| `mixed`, after a run | 19 | 0.30 / 0.44 / 0.60 | 0.29 / 0.31 / 0.35 | 0.07 / 0.12 / 0.17 |
| `weak` | 49 | 0.27 / 0.44 / 0.53 | 0.08 / 0.23 / 0.29 | 0.08 / 0.18 / 0.30 |
| `irregular` | 8 | 0.49 / 0.56 / 0.59 | 0.09 / 0.13 / 0.22 | 0.12 / 0.17 / 0.30 |

Observations:

- **Labels.** The increment-1 rule labels a stretch `dephased` when at least 75 % of its
  double peaks share one modal offset. It labels most post-run slippage `mixed`, because
  that slippage is usually **two-sided**: n−1 and n+1 populations give shadows at both
  −1 and +1.
- **Moderate slippage after HVS-II and HVS-III runs** has main 0.40–0.60, near ≈ 0.34,
  and far ≤ 0.09.
- **Long runs** (≥ 13 identical calls) have main ≈ 0.24. Several populations of similar
  size remain.
- **A globally mixed read** has main ≈ 0.39 with far ≈ 0.23. Its shadows spread over
  every offset.
- **Recoverable share.** About 99 of 266 masked segments, roughly 12 000 of 50 000 masked
  calls, satisfy main ≥ 0.30 and far ≤ 0.12.
- **Primaries survive moderate slippage.** After a moderate run the primary call is
  usually still right. The damage is the double peaks, which become IUPAC calls, add to
  the trim penalty, and trigger exclusions for mixed support.
- **Peeling** reduced mismatches against rCRS in moderate segments of 40–160 calls.
  Five such segments went 23 → 3, 19 → 5, 16 → 5, 11 → 4, and 10 → 2. Segments that were
  already correct stayed correct. Peeling gave no gain on long-run or globally mixed
  segments. True sample differences are included in both counts, and the result has not
  yet been checked against reviewer truth.

### Production amplitudes

The callability shadow model (ADR-0067 calibration amendment) fits the
baseline-corrected event amplitudes of locus evidence and the basecaller's
primary calls instead of the raw values above. With those inputs the far
shadows almost vanish:

- dephased segments: main share q10/median/q90 0.52/0.68/0.79, far share
  0.000/0.007/0.065;
- the long-run segments and the globally mixed read class are explained by
  one-call shadows as well (main share about 0.4–0.7, far share about 0.03).

The far shares in the table above therefore mostly measured background, which
the corrected amplitudes remove. The production label `dephased` means
"explained by shadows of neighbouring calls". It does not measure severity,
and event amplitudes cannot tell a co-located slippage shadow from a
neighbour's peak tail. Severity, and whether a segment can be recovered, are
what the spike below must establish; the main share is the first candidate.

## Limitations

- **Templates come from the read's own calls.** Where the near shadows land on the main
  channel, the primary can flip. That wrong call then serves as a wrong neighbour
  template, which lowers the near weights. Peeling is in effect iterative re-estimation
  of the templates.
- **AC repeats.** Slippage in a dinucleotide repeat produces ±2 shadows, which a model
  that treats only `|k| = 1` as near will call far.
- **Majority slip.** A slipped population above one half flips the primary inside the
  run, with no in-run signature
  ([known limitations](../../validation/known-limitations.md)).
- **Not mixture fractions.** Weights are model coefficients, not mixture or heteroplasmy
  fractions ([ADR-0009](../../decisions/adr/0009-biological-semantics.md), SRS-CALL-007).

## Protocol for the recovery spike

1. Outside the repository, fit and peel every segment that callability labels
   `dephased` ([ADR-0065](../../decisions/adr/0065-result-comparison-downstream.md),
   [ADR-0066](../../decisions/adr/0066-python-limited-to-repository-tooling.md)).
2. Compare primary and peeled calls with the reviewer consensus and with the opposite
   strand. Measure per-position accuracy, created and destroyed differences, behaviour
   at known 309/315/16189 length variants, and ±2 shadows in the AC repeat.
3. Quote aggregates against a named revision.

## Success criterion and promotion

Peeled calls must agree with reviewer consensus clearly better than primary calls in
dephased segments, without creating differences that the opposite strand contradicts.

Promotion requires a new ADR. Recovered calls would be a separate derived
representation, because decoded calls stay immutable
([INV-EVID-001/002](../../architecture/invariants/evidence.md)). They would describe the
dominant ladder's primary sequence, never a mixture fraction. Requirements, contracts,
tests, and validation would follow as a later callability increment.
