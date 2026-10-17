# Discrepancy Catalogue

This catalogue sorts every difference between DNA's calls and reviewer truth on
the development panel by the layer that owns it. Each class says whether the
core caller must fix it, normalization can fix it, or it is not handled yet. A
class that is not handled yet is a TODO here, with its reason.

The layers and their boundary are those of
[ADR-0077](../decisions/adr/0077-normalize-and-empop-crates.md) §5:

- the core caller decides **what** the sample's sequence is;
- `dna-normalize` decides **how** a decided sequence is written, in the lab's
  convention;
- `dna-empop` checks and exports the **format**;
- the [evaluation harness](evaluation-harness.md) classifies.

## Classification rule

A sample is classified from the harness measures. Each row applies only when
the rows above it do not.

| Measures | Class | Owner |
|---|---|---|
| DNA failed on the sample | failed | core or input |
| strict discordant, primary concordant | coverage: a truth allele where DNA has no callable base | core callability, or the laboratory's assay |
| primary discordant, haplotype concordant | representation: the same sequence written another way | `dna-normalize` when DNA's form differs from the reviewers' majority form for that sequence; otherwise a minority form in the truth, recorded only |
| haplotype discordant, haplotype IUPAC-compatible concordant | heteroplasmy: a truth IUPAC code where DNA reports a base | core minor-signal evidence |
| haplotype IUPAC-compatible discordant | sequence: DNA and the truth describe different sequences | core logic, unless adjudication finds a truth error |

The reviewers' majority form per window sequence is measured on the reviewed
corpus. The [design](../design/variant-nomenclature.md) describes the window
rules that follow it.

A convention that changes content is not normalization. For example, the
reviewers do not record the HVS-I C-tract length. A policy for such a
convention classifies a call and never deletes it, and it is a TODO
([roadmap](../proposals/roadmap.md)).

## Development panel v5 (revision 7931e9e)

Of 16,410 active samples, 3,201 differ from the truth on at least one measure.
Region attributions use the first differing allele, so they are approximate.

| Class | Samples | Owner | Status |
|---|---:|---|---|
| Coverage | 1,336 | core callability or assay | TODO; see below |
| Heteroplasmy | 948 | core | TODO, deferred: threshold optimization is not current work |
| Representation | 31 | `dna-normalize` or truth minority form | in progress; see below |
| Sequence | about 870 | core, or adjudication | partly fixed; see below |
| Failed | 9 | core or input | TODO |

**Coverage.**
- About 1,087 samples have the uncovered allele outside the repeat windows,
  mostly `73G` at the HVS-II forward read's start. That case is a matter of
  assay and run quality ([known limitations](known-limitations.md)).
- The rest lie in the C-tracts (HVS-I 118, HVS-II 97) and the AC repeat (31),
  where DNA abstains (`N`).

**Representation.**
- The HVS-I and HVS-II C-tract and AC-repeat window rules are done.
- An HVS-I deletion that trades 3′ across substitutions to `16298` is in
  progress.
- An insertion shift that `16172C` blocks is in progress.
- About two thirds of the remaining 31 are minority forms in the truth. DNA
  already writes the majority form there, for example `523DEL 524DEL` against
  the truth's `514DEL 515DEL`.

**Sequence.**

| Topic | Samples | Owner | Status |
|---|---:|---|---|
| Consensus record re-aligned by `call` and placed elsewhere | 24 | core | fixed (SRS-VAR-017); not counted above |
| Consensus to `call`: the AC-repeat deletion the consensus decides at `523–524` is re-placed beside unresolved `513–516` | 31 | core | open |
| Consensus to `call`: an `N` at a decided run's bound (`310`) | 8 | core | open |
| HVS-II run lengths that no read shows | about 250 | core | TODO: phase recovery |
| Calls at run exits in HVS-II (`310C`, `310DEL 311DEL`) and HVS-III `568–573` | about 125 | core | open |
| HVS-I C-tract lengths that the reviewers do not record (`16193.1C`) | 47 | report convention | TODO |
| AC-repeat insertions DNA misses (`524.1A 524.2C`) | 51 | core | open |
| Read-start and other differences outside the windows (e.g. `16031DEL`) | about 240 | core, or adjudication | open |

## Updating the catalogue

Counts are recomputed from a harness run on the development panel. A change
that fixes a class updates its row and status. A class closed as "not
handled" keeps its reason here, so the next iteration starts from it.
