# Research: point-heteroplasmy evidence in Sanger reads

**Status:** concluded as future work (not promoted). This note is exploratory
evidence and has no production authority. SRS-VAR-009 and ADR-0009 keep
heteroplasmy out of DNA output. A change would need its own decision.

## Question

Reviewers write point heteroplasmy as an IUPAC base, for example `16192Y`.
DNA reports one allele at such a site, or nothing. This note asks whether the
evidence DNA already derives can separate a co-located second allele from
noise well enough to report it.

## Data and method

- **Data:**
  - the fixed 320-sample held-out subset of the science corpus;
  - merged reviewer truth;
  - 40 IUPAC truth sites inside the reviewed intervals.
- **Downstream measurement, outside DNA (ADR-0065, ADR-0066):**
  - for every aligned, unmasked call inside the trim interval, the
    reference-oriented A/C/G/T evidence profile from signal processing;
  - its second-largest weight (the *minor channel*);
  - whether the minor base equals the haplotype base read just before or just
    after it in the sequencing direction.
- Only aggregates are quoted.

## Results

**The current SNV gate does not detect heteroplasmy.** `mixed_supporting_dna`
(two qualifying channels under `secondary_peak_ratio`) flagged at least one
read at only 6 of the 40 sites. Of the 127 flagged read observations, 9 were
at IUPAC sites; the rest were at sites reviewers call homoplasmic.

**The minor channel exists at true sites, but per-read rules drown in noise.**
- In 62 of 66 reads covering an IUPAC site, the two top channels are the IUPAC
  pair.
- A per-read rule "minor ≥ 0.20" fires at 18 true sites and at more than 6,000
  homoplasmic site positions.

**Most noise is a neighbour's signal at the wrong position.** For read
observations with a minor channel ≥ 0.15:

| Site | Minor = previous base | Minor = next base | Minor matches no neighbour |
|---|---|---|---|
| homoplasmic, primary call right | 9,247 | 2,465 | 62 |
| homoplasmic, primary call wrong | 6 | 134 | 294 |
| IUPAC (pair matches) | 1 | 6 | 31 |

About 99 % of the minor channels at homoplasmic sites are shadows of an
adjacent base (dephasing or carry-over), not a second allele. Separately, 434
observations had a primary call that is not the reviewer's base. At such a
call, any heteroplasmy reading would be wrong from the start.

**Cross-strand confirmation helps, but not enough.**

| Site-level rule | IUPAC sites found | Homoplasmic sites flagged |
|---|---|---|
| every callable read shows the same pair, minor ≥ 0.20 | 12 / 40 | 1 |
| both strands show the same pair, minor ≥ 0.15 | 14 / 40 | 121 |
| as above, and no strand's minor is a neighbour's base | 13 / 40 | 6 |

Only one read covers 14 of the 40 sites, so they cannot be confirmed across
strands at all.

## Conclusion

Even the best rule found reaches only about two-thirds precision and recovers
about a third of the sites. That is not enough to report heteroplasmy, and the
gain is small next to the remaining homoplasmic errors.

Reliable point-heteroplasmy evidence depends first on the call at the site
being right. That means:
- the locus is not shifted;
- the call is not masked, dephased, clipped, or a dye blob;
- the minor channel is co-located and not a shadow of a neighbouring call.

Improving the primary calls on the current data therefore comes first. Once
they are right, heteroplasmy evidence can be revisited.

## Prerequisites for a future decision

- Primary-call accuracy at the site.
- A shadow-corrected minor channel. The per-segment shadow fit of read
  callability (ADR-0067) could subtract a neighbour's expected contribution
  before the minor channel is measured.
- Cross-strand confirmation. This is a sample-level rule, and DNA has no
  sample consensus today.
- A labelled output without fractions (SRS-VAR-009), decided under ADR-0009.

**Architecture.** Detection would live in the modality plugin:
- a Sanger read measures a population, so it detects the mixture in its own
  signal;
- NGS would detect it from allele counts across reads.

A modality-neutral record of a mixed site would carry the result through the
core to post-calling notation (`16192Y`). The sequence modality could map
consensus IUPAC codes to the same record.
