---
id: RES-one-strand
status: active
summary: Evidence for the open question on whether a single strand can pass a region.
---

# Research: when may a region read by one strand pass?

**Status:** exploratory evidence for a decision. This note has no production
authority. It supports open question 6 of
[PROP-0004](../../proposals/0004-review-triage.md#open-questions); the policy
remains that a pass needs both strands until the user decides otherwise.

## Question

In phase 0 the draft detector `region.one_strand` flags 1,920 of the 2,618
active development samples and is the largest lever on A. Which kinds of
one-strand stretches exist, how often does each kind hide an error, and which
exception, if any, would raise A without letting errors pass?

## Data and method

- **Panel and truth.** Development panel v3 only (2,618 active samples), the
  lab's statistics database of 2026-10-15, `config/dna.toml`, rCRS, revision
  a244644. The validation panel was not used.
- **Run.** `dna-eval run --triage` ([draft triage](../../validation/evaluation-harness.md#draft-triage)).
  Its draft outcomes are the baseline: 515 passed samples, A = 0.178
  corpus-weighted.
- **Extra evidence, outside the repository.** A throwaway copy of the harness
  wrote, per sample, every region position that is not read cleanly by both
  strands, with the strand counts and consensus marks the triage reads, plus,
  from each read's alignment columns, what the missing strand's reads do there
  and the relative quality and minor signal of the present strand's calls. An
  analysis script re-implemented the draft decision on that dump; it
  reproduces the harness exactly (0 sample outcomes differ under the draft,
  642 passes with `region.one_strand` off, A = 0.1775).
- **Stretch.** A maximal run of consecutive positions of one region read
  cleanly by the same single strand. There are 9,258 stretches in 1,921
  samples.
- **Classes.**
  - *Location:* `edge` touches the first or last 8 positions of a region
    (16024–16031, 16358–16365, 73–80, 333–340, 438–445, 569–576);
    `repeat` reaches within 3 positions of a repeat window of the
    [measures](../../validation/evaluation-harness.md#measures) (`hv1_polyc`
    16180–16195, `hv2_polyc` 302–316, `hv3_polyc` 560–576, `ac_repeat`
    513–526); `edge_to_repeat` does both, typically the half of a region
    between a primer end and a poly-C; `other` does neither. HVS-III's last
    positions lie in `hv3_polyc`, so its edge counts as `edge_to_repeat`.
  - *Length:* fewer than 5, 5–19, or at least 20 positions.
  - *Cause* (majority over the stretch's positions), from the reads of the
    missing strand whose primer reads the region:
    `dephased` — the read covers the position with a dephased masked call, or
    its alignment stops before the position and its remaining calls, enough to
    reach it, are mostly dephased (lost after a homopolymer);
    `poor_signal` — the same with mixed, weak, or irregular calls;
    `read_end` — the read physically ends before the position;
    `read_start` — the position lies before the read's first aligned call
    (its primer end);
    `rejected` — no such read, but a rejected or poor trace of that strand
    could have read the region;
    `no_read` — no trace of that strand reads the region (a missing trace);
    `site` — the read is in phase and callable at the position but the
    consensus site counts it as uninformative (an untrusted call, a support
    veto, a neighbouring `N`, or an indel inside the site).
  - *Single-strand quality:* per position, the best read of the present
    strand; over the stretch, its lowest relative quality score (`rq_min`,
    uncalibrated, 0–60) and its highest minor signal (`sr_max`, the second
    over the first channel of the call's normalised profile). The present
    strand's calls are in phase by construction: a clean observation is an
    unmasked call.
- **Allowing a class.** Its stretches no longer count for
  `region.one_strand`, and the `masked` mark of their positions is cleared:
  a dephased read of the other strand is what makes those sites masked. Every
  other finding still applies. Rows marked "masked kept" clear only the
  one-strand flag.
- **Measures.** Primary, strict, and IUPAC-compatible from the harness. The
  haplotype measures (representation-normalised, with and without IUPAC
  compatibility) are not on `main`; they come from a run of the unmerged
  haplotype comparison on the same panel and DNA calls, whose three other
  measures equal this run's for every sample. Discordant positions count
  inside a stretch when their anchor base lies in it; for the haplotype
  measure, a difference block overlapping the stretch counts.
- **Rates.** For each class or rule: the samples it newly passes (for a
  class, those where it is the only thing blocking pass), the corpus-weighted
  A gain, the discordant positions inside the allowed stretches of those
  samples, and the newly passed samples that are not concordant (leaks), with
  the one-sided 95 % Clopper–Pearson upper bound (UB). Under the draft,
  50 / 515 passed samples leak on the primary measure (UB 0.121), 9 / 515 on
  IUPAC-compatible (UB 0.030), and 1 / 515 on haplotype IUPAC-compatible
  (UB 0.009).

## Results

### Which stretches exist and how often they hold a difference

Every stretch, whatever blocks its sample otherwise. "Discordant inside" is
primary / strict / IUPAC-compatible positions; "blocks" are haplotype
IUPAC-compatible difference blocks overlapping the stretch.

| Cause | Stretches | Samples | Positions | Discordant inside | Blocks | Stretches with a primary difference |
|---|---:|---:|---:|---:|---:|---:|
| `dephased` | 5,240 | 1,374 | 455,483 | 195 / 221 / 148 | 378 | 3.3 % |
| `poor_signal` | 1,064 | 617 | 68,261 | 35 / 38 / 24 | 31 | 2.3 % |
| `read_start` | 367 | 190 | 2,297 | 5 / 18 / 5 | 6 | 1.4 % |
| `read_end` | 582 | 343 | 74,605 | 17 / 25 / 5 | 9 | 2.6 % |
| `rejected` | 158 | 33 | 9,177 | 8 / 13 / 6 | 15 | 3.8 % |
| `no_read` | 1,148 | 201 | 44,937 | 48 / 74 / 46 | 110 | 3.1 % |
| `site` | 699 | 564 | 3,932 | 38 / 38 / 31 | 64 | 3.9 % |

| Location, length | Stretches | Samples | Positions | Discordant inside | Blocks | Stretches with a primary difference |
|---|---:|---:|---:|---:|---:|---:|
| `edge`, < 5 | 474 | 346 | 929 | 2 / 7 / 2 | 13 | 0.4 % |
| `edge`, 5–19 | 585 | 520 | 5,212 | 5 / 14 / 4 | 39 | 0.7 % |
| `edge`, ≥ 20 | 693 | 572 | 43,728 | 8 / 15 / 5 | 47 | 1.2 % |
| `edge_to_repeat`, < 20 | 217 | 205 | 1,844 | 7 / 7 / 7 | 12 | 3.2 % |
| `edge_to_repeat`, ≥ 20 | 3,298 | 1,306 | 476,481 | 171 / 203 / 113 | 117 | 4.2 % |
| `repeat`, < 5 | 561 | 384 | 1,439 | 16 / 16 / 16 | 117 | 2.9 % |
| `repeat`, 5–19 | 586 | 444 | 7,548 | 89 / 92 / 83 | 112 | 11.9 % |
| `repeat`, ≥ 20 | 894 | 696 | 92,358 | 35 / 53 / 25 | 33 | 2.8 % |
| `other`, < 5 | 1,030 | 513 | 2,123 | 5 / 5 / 4 | 82 | 0.5 % |
| `other`, 5–19 | 544 | 314 | 5,847 | 3 / 8 / 3 | 36 | 0.6 % |
| `other`, ≥ 20 | 376 | 277 | 21,183 | 5 / 7 / 3 | 5 | 1.3 % |

| Single-strand quality | Stretches | Samples | Positions | Discordant inside | Blocks | Stretches with a primary difference |
|---|---:|---:|---:|---:|---:|---:|
| `rq_min` ≥ 40, `sr_max` < 0.3 | 7,202 | 1,861 | 434,655 | 243 / 277 / 206 | 518 | 2.7 % |
| `rq_min` ≥ 40, `sr_max` ≥ 0.3 | 1,481 | 1,019 | 179,336 | 74 / 98 / 38 | 54 | 4.6 % |
| `rq_min` < 40 | 575 | 330 | 44,701 | 29 / 52 / 21 | 41 | 4.3 % |

- **Location decides more than cause.** Stretches that stay away from the
  repeat windows hold a primary difference in 0.4–1.3 % of cases; stretches
  that reach a window, 2.8–11.9 %. Within one location the causes differ
  little.
- **Most one-strand positions come from dephasing at a poly-C.** 73 % of the
  stretches' positions lie in `edge_to_repeat` stretches, three quarters of
  them dephased: one primer reads up to a poly-C and the opposite one dephases
  there, so a whole half of HVS-I or HVS-II is read by one strand.
- **What the differences are.** Of the 346 primary discordant positions
  inside stretches, 247 lie in a repeat window or within 3 positions of one,
  mostly HVS-II poly-C insertions (116), HVS-I poly-C substitutions and truth
  IUPAC alleles (49), and AC-repeat deletions (23). The 99 outside are split
  below.
- **Quality hardly separates.** Clean single strands (`rq_min` ≥ 40, no minor
  signal) hold a difference in 2.7 % of stretches against 4.3–4.6 % for the
  others; most differences sit in clean stretches.
- **`site` is not an absence.** The other strand is present but untrusted
  at a consensus site, and the site is `masked` anyway. Its rare edge and
  edge-to-repeat stretches hold a difference in 6 % and 50 % of cases.

Discordant positions inside stretches but outside the windows ± 3, by what
differs; a reference run is a run of equal bases in rCRS, and "next to" means
one position before or after it. The last column counts those in the samples
that R6m (below) newly passes.

| What differs | Primary | Strict | IUPAC | In R6m newly passed (P / S / I) |
|---|---:|---:|---:|---:|
| truth IUPAC (heteroplasmy) | 53 | 53 | 0 | 1 / 1 / 0 |
| substitution | 18 | 69 | 18 | 0 / 0 / 0 |
| indel in or next to a reference run ≥ 3 | 12 | 12 | 12 | 0 / 0 / 0 |
| indel in or next to a reference run of 2 | 10 | 10 | 10 | 0 / 0 / 0 |
| other indel | 6 | 6 | 6 | 0 / 0 / 0 |

Homopolymer length variants also occur outside the four windows (in the
corpus, among others `249DEL`, `106`–`111DEL`, `459DEL`, `337DEL`,
`16298DEL`, `16166DEL`, `455.1T`, `291.1A`, and `194DEL`): 22 of the 99
outside-window differences are such indels, none of them in a sample R6m
newly passes. On this panel the exception would not have passed one, but
nothing in R6m protects against it; R6m-h below measures the cost of doing
so.

### Each class allowed alone

"Only blocker" counts the samples that pass once the class is allowed and
nothing else changes.

| Class | Only blocker | ΔA | Discordant inside | Leaks primary (UB) | Strict (UB) | IUPAC (UB) | Haplotype (UB) | Haplotype IUPAC (UB) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| every stretch | 434 | 0.178 | 14 / 15 / 2 | 55 (0.156) | 57 (0.161) | 19 (0.064) | 44 (0.129) | 4 (0.021) |
| `dephased` | 133 | 0.040 | 8 / 9 / 1 | 19 (0.204) | 20 (0.213) | 8 (0.107) | 15 (0.170) | 1 (0.035) |
| `poor_signal` | 35 | 0.018 | 0 / 0 / 0 | 3 (0.207) | 3 (0.207) | 3 (0.207) | 0 (0.082) | 0 (0.082) |
| `read_start` | 10 | 0.002 | 0 / 0 / 0 | 1 (0.394) | 2 (0.507) | 0 (0.259) | 1 (0.394) | 0 (0.259) |
| `read_end` | 6 | 0.004 | 0 / 0 / 0 | 1 (0.582) | 1 (0.582) | 0 (0.393) | 1 (0.582) | 0 (0.393) |
| `rejected` | 4 | 0.001 | 0 / 0 / 0 | 0 (0.527) | 0 (0.527) | 0 (0.527) | 0 (0.527) | 0 (0.527) |
| `no_read` | 25 | 0.018 | 0 / 0 / 0 | 3 (0.282) | 3 (0.282) | 0 (0.113) | 3 (0.282) | 0 (0.113) |
| `site` | 62 | 0.022 | 1 / 1 / 0 | 9 (0.240) | 9 (0.240) | 2 (0.098) | 7 (0.202) | 0 (0.047) |
| `edge` | 63 | 0.039 | 1 / 1 / 0 | 4 (0.139) | 5 (0.160) | 0 (0.046) | 4 (0.139) | 0 (0.046) |
| `edge_to_repeat` | 100 | 0.030 | 5 / 6 / 0 | 12 (0.189) | 13 (0.201) | 3 (0.076) | 9 (0.153) | 0 (0.030) |
| `repeat` | 29 | 0.013 | 1 / 1 / 0 | 5 (0.329) | 5 (0.329) | 1 (0.153) | 4 (0.288) | 0 (0.098) |
| `other` | 48 | 0.010 | 0 / 0 / 0 | 4 (0.181) | 4 (0.181) | 1 (0.095) | 3 (0.154) | 0 (0.061) |
| < 5 positions | 90 | 0.041 | 0 / 0 / 0 | 8 (0.155) | 8 (0.155) | 2 (0.068) | 6 (0.127) | 0 (0.033) |
| 5–19 positions | 46 | 0.025 | 2 / 2 / 0 | 5 (0.215) | 6 (0.241) | 0 (0.063) | 5 (0.215) | 0 (0.063) |
| ≥ 20 positions | 123 | 0.040 | 7 / 8 / 1 | 17 (0.202) | 18 (0.211) | 4 (0.073) | 16 (0.192) | 2 (0.051) |
| `edge`, < 5 | 26 | 0.016 | 0 / 0 / 0 | 0 (0.109) | 0 (0.109) | 0 (0.109) | 0 (0.109) | 0 (0.109) |
| `edge`, 5–19 | 23 | 0.016 | 1 / 1 / 0 | 3 (0.304) | 4 (0.355) | 0 (0.122) | 3 (0.304) | 0 (0.122) |
| `edge`, ≥ 20 | 8 | 0.005 | 0 / 0 / 0 | 1 (0.471) | 1 (0.471) | 0 (0.312) | 1 (0.471) | 0 (0.312) |
| `rq_min` ≥ 40, `sr_max` < 0.3 | 324 | 0.137 | 3 / 4 / 0 | 37 (0.147) | 39 (0.154) | 13 (0.063) | 27 (0.113) | 1 (0.015) |

- **One class rarely frees a sample.** Allowing every stretch frees 434
  samples (ΔA 0.178), but the classes alone free far fewer, because a sample
  usually has several kinds of stretch: the two dephased halves of a region,
  a primer-end stretch, and short stretches at sites. The most common sets
  are `dephased` alone (133), `site` alone (62), `dephased` with
  `poor_signal` (37), and `poor_signal` alone (35).
- **The newly passed leaks are the known stage-1 leaks.** Of the 55 primary
  leaks of every stretch, 41 differ only outside the allowed stretches, at
  truth IUPAC alleles both strands read or in how an HVS-I poly-C length is
  written; the haplotype IUPAC-compatible measure removes all but 4.
- **What one strand hid.** Inside allowed stretches, 9 of the 14 samples
  differ at the same heteroplasmy inside the HVS-I poly-C window (truth
  `16192Y`, DNA `16192T`); the others at truth `Y` or `N` alleles elsewhere
  (3) and at two poly-C lengths (`315.1C` called by DNA only, `16193.1C` in
  the truth only). One of the 4 haplotype IUPAC-compatible leaks is a
  truth `16193.1C` that DNA missed inside a dephased single-strand stretch of
  the HVS-I poly-C, read with `rq_min` 52 and a minor signal of 0.02: a single
  clean strand can get a poly-C length wrong.

### Candidate exception rules

Each rule is applied to the draft alone. The headline is the **marginal
risk**: the differences inside the stretches a rule allows, in the samples it
newly passes. Sample-level leaks follow separately; most of them are
heteroplasmy at positions both strands read, which the draft already passes
(50 / 515 primary leaks) and which the minor-signal gate is meant to stop.
The per-stretch bounds treat stretches as independent trials, although
several come from one sample, so they are approximate; the sample is the
unit to quote beside them (R6m: 1 of 60 newly passed samples has a
difference inside an allowed stretch, UB 0.077; R6m-h4: 0 of 40, UB 0.072).

**R6m is defined exactly as the research harness policy `one_strand_windows`:**
`region.one_strand` applies only at a position inside a repeat window
(`hv1_polyc` 16180–16195, `hv2_polyc` 302–316, `hv3_polyc` 560–576,
`ac_repeat` 513–526) or within 3 positions of one (so 16177–16198, 299–319,
510–529, and 557–579); elsewhere a position read cleanly by one strand does
not hold its region back. Every other finding, `masked` included, still
applies. The margin covers the positions a read reaches just after leaving a
repeat, where it may still be dephased. On this panel it changes no pass:
with or without it, R6m passes the same 60 samples. R6m-h3 and R6m-h4 also
hold back one-strand positions in or next to a reference run of at least 3
or 4 equal bases.

| Rule | New pass | ΔA | A | Allowed stretches | With a primary difference (UB) | Strict (UB) | IUPAC (UB) | Positions | Discordant positions P / S / I | Haplotype IUPAC blocks |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| R1 `edge`, < 5 positions, not `site` | 26 | 0.016 | 0.193 | 26 | 0 (0.109) | 0 (0.109) | 0 (0.109) | 55 | 0 / 0 / 0 | 0 |
| R4 `edge`, any length, not `site` | 63 | 0.039 | 0.216 | 72 | 1 (0.064) | 1 (0.064) | 0 (0.041) | 1,765 | 1 / 1 / 0 | 0 |
| R6 outside windows ± 3, not `site`, masked cleared | 70 | 0.042 | 0.219 | 82 | 1 (0.057) | 1 (0.057) | 0 (0.036) | 1,990 | 1 / 1 / 0 | 0 |
| **R6m** outside windows ± 3 (`one_strand_windows`) | 60 | 0.036 | 0.214 | 69 | 1 (0.067) | 1 (0.067) | 0 (0.042) | 1,885 | 1 / 1 / 0 | 0 |
| R6m-h3 R6m, not in or next to a reference run ≥ 3 | 18 | 0.013 | 0.190 | 19 | 0 (0.146) | 0 (0.146) | 0 (0.146) | 48 | 0 / 0 / 0 | 0 |
| R6m-h4 R6m, not in or next to a reference run ≥ 4 | 40 | 0.024 | 0.201 | 43 | 0 (0.067) | 0 (0.067) | 0 (0.067) | 1,323 | 0 / 0 / 0 | 0 |
| R8 other strand dephased/poor/ended, windows included | 252 | 0.091 | 0.268 | 678 | 11 (0.027) | 12 (0.029) | 1 (0.007) | 59,455 | 11 / 12 / 1 | 2 |
| R9 R8 with `rq_min` ≥ 40, `sr_max` < 0.3 | 162 | 0.059 | 0.236 | 422 | 2 (0.015) | 3 (0.018) | 0 (0.007) | 29,955 | 2 / 3 / 0 | 0 |
| every stretch (reference) | 434 | 0.178 | 0.356 | 1,105 | 14 (0.020) | 15 (0.021) | 2 (0.006) | 76,725 | 14 / 15 / 2 | 3 |

Sample-level leaks of the same rules (newly passed samples that are not
concordant, and all passed samples):

| Rule | New pass | Primary (UB) | Strict (UB) | IUPAC (UB) | Haplotype (UB) | Haplotype IUPAC (UB) | All passed primary (UB) | All passed haplotype IUPAC (UB) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| R1 | 26 | 0 (0.109) | 0 (0.109) | 0 (0.109) | 0 (0.109) | 0 (0.109) | 50/541 (0.116) | 1/541 (0.009) |
| R4 | 63 | 4 (0.139) | 5 (0.160) | 0 (0.046) | 4 (0.139) | 0 (0.046) | 54/578 (0.116) | 1/578 (0.008) |
| R6 | 70 | 5 (0.144) | 6 (0.162) | 0 (0.042) | 5 (0.144) | 0 (0.042) | 55/585 (0.116) | 1/585 (0.008) |
| **R6m** | 60 | 5 (0.167) | 6 (0.188) | 0 (0.049) | 5 (0.167) | 0 (0.049) | 55/575 (0.118) | 1/575 (0.008) |
| R6m-h3 | 18 | 0 (0.153) | 0 (0.153) | 0 (0.153) | 0 (0.153) | 0 (0.153) | 50/533 (0.117) | 1/533 (0.009) |
| R6m-h4 | 40 | 3 (0.183) | 3 (0.183) | 0 (0.072) | 3 (0.183) | 0 (0.072) | 53/555 (0.119) | 1/555 (0.009) |
| R8 | 252 | 33 (0.172) | 35 (0.181) | 14 (0.086) | 25 (0.136) | 2 (0.025) | 83/766 (0.129) | 3/766 (0.010) |
| R9 | 162 | 18 (0.160) | 20 (0.174) | 9 (0.095) | 11 (0.110) | 0 (0.018) | 68/677 (0.122) | 1/677 (0.007) |
| every stretch | 434 | 55 (0.156) | 57 (0.161) | 19 (0.064) | 44 (0.129) | 4 (0.021) | 105/948 (0.129) | 5/948 (0.011) |

Other rules explored, on the same windows ± 3 (sample level):

| Rule | New pass | ΔA | Primary (UB) | IUPAC (UB) | Haplotype IUPAC (UB) |
|---|---:|---:|---:|---:|---:|
| R2 `edge`, < 10 positions, not `site` | 37 | 0.023 | 1 (0.122) | 0 (0.078) | 0 (0.078) |
| R3 `edge`, < 20 positions, not `site` | 54 | 0.033 | 3 (0.137) | 0 (0.054) | 0 (0.054) |
| R4m R4, masked kept | 54 | 0.034 | 4 (0.162) | 0 (0.054) | 0 (0.054) |
| R5 `edge` or `other`, not `site` | 70 | 0.042 | 5 (0.144) | 0 (0.042) | 0 (0.042) |
| R7 R6 with `rq_min` ≥ 40, `sr_max` < 0.3 | 68 | 0.040 | 5 (0.148) | 0 (0.043) | 0 (0.043) |
| R10 R9, outside windows ± 3 | 40 | 0.020 | 3 (0.183) | 0 (0.072) | 0 (0.072) |
| R11 any stretch but `site`, windows included | 292 | 0.116 | 36 (0.160) | 14 (0.074) | 2 (0.021) |
| R13 R6 outside windows ± 3, R9 inside | 201 | 0.087 | 21 (0.147) | 10 (0.083) | 0 (0.015) |

- **Outside the windows (R1–R7, R6m)** the allowed stretches of newly passed
  samples hold at most one difference (R4, R6, R6m: 1 in 69–82 stretches, a
  truth `Y`, heteroplasmy read by one strand; R1 and the homopolymer variants:
  none), and none under the IUPAC-compatible or haplotype IUPAC-compatible
  measures. Their sample-level primary leaks are all truth
  IUPAC alleles, the others at positions both strands read. The quality gate
  changes almost nothing (R7 against R6), and edges carry most of the gain
  (R4 against R6).
- **Holding back reference homopolymers costs much of the gain.** Positions
  in or next to a run of at least 3 bases are 22–50 % of each region; R6m-h3
  keeps 0.013 of R6m's 0.036 (cost 0.023) and R6m-h4 keeps 0.024 (cost
  0.012). On this panel neither removes a difference R6m let through.
- **Inside the windows (R8, R9, R11, R13)** A rises to 0.24–0.29. The
  per-stretch rate looks low (0.5–1.6 % of allowed stretches) because the
  allowed stretches are long and many, but the newly passed samples bring
  9–14 IUPAC-compatible leaks, all HVS-I poly-C length variants written at
  other positions (representation: the haplotype measure finds them
  concordant), and without a quality gate a true single-strand poly-C length
  miss (R8, R11). R9 and R13 avoid that sample only because another of its
  stretches fails the gate; the stretch with the miss passes it.
- **Cause-based rules stay out of the decision.** R8, R9, R10, and R13 rest
  on this note's cause heuristic, which has not been validated against
  inspected traces; they are reported for context only.
- **The bounds are wide.** With 18–70 new passes, zero events still bound the
  rate at 0.04–0.15. A 1 % bound at zero leaks needs about 300 newly passed
  samples.

A research build of the harness measured R6m as the sweep policy
`one_strand_windows` and reproduced it: 575 passed samples, A 0.214, primary
leaks 55 / 575, IUPAC-compatible 9 / 575. That sweep row is not on main; it
lands with the rule's implementation. R6m-h3 and R6m-h4 stay in the analysis
script: the triage evidence carries no reference sequence, and adding it only
for a research row would need a new facade accessor or a second FASTA reader
in the harness. If the decision needs R6m-h, the evidence can take the
reference runs then.

### External norms

- The ISFG DNA Commission's revised mtDNA guidelines (Parson et al. 2014,
  *Forensic Sci. Int. Genet.* 13:134–142,
  <https://pubmed.ncbi.nlm.nih.gov/25117402/>) recommend that every reported
  position be covered by forward and reverse sequence, with a same-strand
  replicate from a different primer as the fallback; the motivation is
  phantom mutations, notably behind the HVS-I and HVS-II C-tracts.
- The New York City OCME forensic biology interpretation guidelines
  (<https://www.nyc.gov/assets/ocme/downloads/pdf/technical-manuals/forensic-biology-technical-manuals/interpretation_guidelines_091217.pdf>)
  accept single-strand sequence only under special circumstances.
- Both are summarised, not quoted: neither text could be retrieved verbatim
  for this note. Neither allows an unconditional single-strand pass. R6m is
  narrower than a general single-strand exception, but it is still an
  exception to that practice and needs the user's decision.

## Recommendation for the decision

- **Interaction with the minor-signal gate (PR #56).** The cross-strand minor
  kinds of PR #56 (`strands`, `unshadowed`, `upstream`, `noise`) need a
  forward and a reverse read showing the same pair, so inside an allowed
  one-strand stretch they can never fire; only the single-read kinds (`lone`,
  `any`) can. Under R6m a heteroplasmy in an allowed stretch is therefore
  invisible to the default `upstream` gate at 100, and that is exactly the
  one difference R6m let through on this panel. A one-strand exception
  therefore needs its own single-strand minor rule.
- **Pre-registered rule: a pairing.** The exception to test is "R6m-h4 (or
  R6m) plus a single-strand minor rule applied only on the allowed one-strand
  positions": `lone`-style, possibly with a stricter threshold than on
  positions both strands read.
- **First choice: R6m-h4; second: R6m.** On the development panel:

  | Rule | ΔA | Allowed stretches with a difference | Newly passed samples with a difference inside | Sample leaks: primary / IUPAC / haplotype IUPAC (UB) |
  |---|---:|---:|---:|---:|
  | R6m-h4 | +0.024 (40 samples) | 0 of 43 | 0 of 40 (UB 0.072) | 3 / 0 / 0 (0.072) |
  | R6m | +0.036 (60 samples) | 1 of 69 | 1 of 60 (UB 0.077) | 5 / 0 / 0 (0.049) |

  Neither homopolymer hold removes a difference found on this panel. But 22
  of the 99 outside-window differences inside stretches are homopolymer-type
  indels, such length variants are common outside the windows in the corpus
  (`249DEL` alone in about 28 % of samples), and a single clean strand can
  miscount a run (the `16193.1C` miss above). R6m-h4 keeps about two thirds
  of R6m's gain and closes that failure mode outside the windows. **The
  trade-off is the user's decision:** about 0.012 of A against protection
  from single-strand run miscounts that this panel did not happen to show.
- **Keep `masked` as it is** in either rule. Clearing the masked mark where
  the other strand's read is dephased (R6, +0.042) needs the strand of each
  masked read, which the triage evidence does not record yet.
- **Decision basis: on top of the minor-signal gate.** Every number above was
  measured before the gate (PR #56), the unresolved-call fix (PR #59), and
  the bounded-run exemption (ADR-0076). Re-measured on dev panel v4 with all
  three, the harness policies `outside_repeats` (R6m) and
  `outside_repeats_h4` (R6m-h4):

  | Policy | Passed | ΔA | Primary leaks | Haplotype IUPAC-compatible leaks (UB) |
  |---|---:|---:|---:|---:|
  | draft (both strands) | 458 | — | 28 | 1 (0.0103) |
  | R6m | 513 | +0.034 | 30 | 1 (0.0092) |
  | R6m-h4 | 494 | +0.022 | 29 | 1 (0.0096) |
  | R6m-h4 + single-strand `any@150` | 494 | +0.022 | 29 | 1 (0.0096) |
  | R6m-h4 + single-strand `any@100` | 493 | +0.021 | 29 | 1 (0.0096) |

  R6m-h4 now passes 36 more samples than the draft instead of 40, because the
  gate and the exemption change which samples the draft passes and which
  one-strand stretches block them. Of the 36, one leaks under the primary
  measure, at a heteroplasmic position both strands read, so not inside the
  exception: 0 leaks inside the allowed one-strand positions, as before.
  The single-strand rule at 150 changes no outcome on this panel. R1 is not
  a hypothesis; the fixed sequence stops at the first failure. The choice of
  R6m-h4 is recorded with the approval decisions in
  [PROP-0005](../../proposals/0005-evaluation-feedback-loop.md).
- **Selection effect.** About 17 rules were tried on the development panel,
  so the best of them looks better there than it will elsewhere. Validation
  will pre-register the pairing above (R6m-h4 first, R6m second) with R1 as
  the fallback and test them once on the validation panel; the validation
  plan is written separately.
- **Do not allow yet:** one strand inside a repeat window or its margin,
  whatever its quality or cause (R8–R13). It holds the largest gain (up to
  +0.116), but also single-strand poly-C length errors and representation
  leaks, and its cause rules are unvalidated. Revisit after notation
  normalisation lands and with a dedicated length check.
- **Never allow** the `site` class: the other strand is present but
  untrusted there.
- **Narrowest fallback:** R1 (an edge stretch of fewer than 5 positions)
  gains +0.016 with no difference and no leak under any measure (UB 0.109).

## Limitations

- Development panel only, and small counts per class: the bounds above are
  for discussion, not for approval.
- The cause is a heuristic over the alignment and callability: the nearest
  read of the missing strand, within 400 positions, whose primer reads the
  region; a read whose remaining calls could reach the position counts as
  lost rather than ended.
- Reference runs come from rCRS; a homopolymer length variant can also
  sit where the sample, not the reference, has a run.
- The haplotype measures come from an unmerged harness; per-stretch haplotype
  differences use whole difference blocks, which can be long.
- The analysis script and the per-position dump stay outside the repository
  under the [data policy](../../governance/data.md).
