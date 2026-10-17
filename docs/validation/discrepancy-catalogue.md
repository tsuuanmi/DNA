# Discrepancy Catalogue

This catalogue sorts every difference between DNA's calls and reviewer truth on
the development panel by the layer that owns it. Each class says whether the
core caller must fix it, normalization can fix it, the truth needs a reviewer's
re-check, or it is not handled yet. A class that is not handled yet is a TODO
here, with its reason.

Reviewer truth is the reference, not an infallible one. A difference becomes a
DNA defect only when the read evidence shows that DNA is wrong. When the
evidence supports DNA, DNA is not changed: the case goes to the reviewers for
a re-check (see *Truth re-check*).

The layers and their boundary are those of
[ADR-0077](../decisions/adr/0077-normalize-and-empop-crates.md) §5:

- the core caller decides **what** the sample's sequence is;
- `dna-normalize` decides **how** a decided sequence is written, in the lab's
  convention;
- `dna-empop` checks and exports the **format**;
- the [evaluation harness](evaluation-harness.md) classifies.

## Working hypothesis

Current work improves the core caller. It is judged on the hypothesis that a
read's signal is either good, or clearly bad because of the laboratory, so
every difference is attributed to one of three owners:

- **the laboratory did poorly**: no read, weak signal, a short read, or the
  truth itself (`lab`);
- **the laboratory did well, but the read carries dephasing, noise, or a
  length mixture** (`technology`): the core cannot recover it. Cleaning the
  signal, correcting noise, and recovering phase belong to the Sanger
  modality (`dna-sanger`) and are later work, listed among the
  [roadmap's deferred directions](../proposals/roadmap.md#deferred-product-directions):
  phase recovery, signal denoising, and length-mixture decomposition;
- **the core** (`dna`): the current focus.

## Classification rule

The harness gives each sample one [class](evaluation-harness.md#discrepancy-classes)
and one [cause](evaluation-harness.md#discrepancy-causes), with the evidence
for each finding. Those two sections define the rules; this catalogue only
reads their counts and records what is done about each part.

- The **class** names the layer whose measure first disagrees: `failed`,
  `not_comparable`, `sequence`, `heteroplasmy`, `unanalysed`,
  `representation`, `coverage`.
- The **cause** names who can act: `dna` (a defect, or `deferred` where DNA
  leaves something undone by design), `lab` (the assay, a run, or the truth),
  `technology` (a limit of Sanger reads), or `undetermined` (a reviewer
  decides).

The haplotype measures compare only where the reviewers analysed: a DNA call
outside the truth's analysed intervals is `unanalysed`, information for a
review rather than a defect.

A convention that changes content is not normalization. For example, the
reviewers do not record the HVS-I C-tract length. It is handled by a labelled
convention that classifies a call and never deletes it
([PROP-0007](../proposals/0007-empop-export-and-report-conventions.md)).

## Development panel v5 (revision b075635)

Of 16,410 active samples, 13,302 agree with the truth on every measure. The
harness classes and causes the rest. A sample's cause is its first defect, in
the order `dna`, `lab`, `undetermined`; a sample with deferred findings only
takes the first of them:

| Class | Samples | dna | lab | undetermined | dna (deferred) | lab (deferred) | technology (deferred) |
|---|---:|---:|---:|---:|---:|---:|---:|
| coverage | 1,593 | 91 | 656 | 25 | 415 | 0 | 406 |
| heteroplasmy | 992 | 0 | 84 | 86 | 656 | 0 | 166 |
| sequence | 346 | 12 | 2 | 301 | 0 | 31 | 0 |
| unanalysed | 133 | 0 | 0 | 133 | 0 | 0 | 0 |
| representation | 33 | 11 | 22 | 0 | 0 | 0 | 0 |
| failed | 9 | 9 | 0 | 0 | 0 | 0 | 0 |
| not_comparable | 2 | 0 | 0 | 2 | 0 | 0 | 0 |

Under the working hypothesis this splits the work:

- **the core now**: 123 samples (`dna`);
- **the laboratory**: 764 samples (`lab`), most of them no read or weak
  signal at `73`, plus the truth re-check below;
- **a reviewer**: 547 samples (`undetermined`), of which 301 sequence
  differences to adjudicate and 133 calls where the reviewers did not
  analyse;
- **deferred**: 1,674 samples, being the IUPAC output (1,071), the
  reviewers' HVS-I length convention (31), and signal processing (572).

Since revision 7931e9e the `sequence` class fell from 879 to 346 samples and
`coverage` rose from 1,336 to 1,593: calls DNA could not stand behind became
abstentions (see *Fixed*).

**Coverage** (a truth allele where DNA has no callable base). By finding:
- `lab`: no read covers it (`no_read`, mostly `73G` at the HVS-II forward
  read's start, an assay and run-quality matter, see
  [known limitations](known-limitations.md)), a weak-signal veto
  (`weak_call`), the read's own end (`read_edge`), or masked weak signal;
- `technology`, deferred to signal processing: phase loss, a run no read
  shows whole (`run_unread`), or the margin beside a phase-loss mask;
- `dna`, deferred: a truth IUPAC code shown as mixed calls;
- `dna`: an `N` call or a call DNA does not trust (`untrusted_call`, 89
  samples), most at the AC repeat (`523`/`524`) and a few at `315.1C`.

**Heteroplasmy** (a truth IUPAC code where DNA reports a base). Deferred by
design: DNA writes no IUPAC code yet, and threshold work is not current. The
other causes are a second base below the limit of detection (`technology`,
deferred), no
second base in any clean read (`lab`, for a truth re-check), and one strand or
a code of three bases (`undetermined`). Most sit in the HVS-I C-tract and the
regions `101-301`, `16051-16179`, and `16196-16365`.

**Sequence** (the sequences differ). By the region of the differing blocks:
the HVS-II C-tract (110 samples, half with a length difference), the HVS-I
C-tract (77), the HVS-I read start `16024-16050` (76), and fewer elsewhere.
DNA's defects here are phantom calls inside the HVS-I C-tract (`16185-16193`, 12 samples); 2 samples
have clean reads on both strands showing the reference where the truth has an
allele (`lab`, re-check); in 31 every DNA-only call is an HVS-I C-run length
the reviewers do not record (`lab`, deferred, `convention_not_recorded`); the
other 301 go to adjudication.

**Unanalysed** (DNA calls only where the reviewers did not analyse). Mostly
the HVS-II run exit (`310C`, `310DEL`, `311DEL`, `312-318DEL`) and `574C`.
Earlier counted as DNA phantoms, they cannot be judged against this truth;
only a review there verifies them.

**Representation** (the same sequence written differently). DNA writes a
minority form in 11 samples (`dna_minority_form`); the truth does in 22,
which go to the reviewers for notation consistency.

## Fixed

| Topic | Change |
|---|---|
| A consensus record re-aligned by `call` and placed elsewhere | SRS-VAR-017; a misplaced record has no callable bases (SRS-VAR-019) |
| An AC-repeat deletion re-placed beside unresolved `513-516` | SRS-ALN-017, SRS-VAR-018 |
| An excluded call read as the reference | SRS-VAR-019: its positions are not callable |
| HVS-I C-run lengths the reviewers do not record (`16193.1C`, `16193DEL`) counted as open differences | SRS-NOM-017: the `hvs1_c_length` convention labels them, and the harness gives 31 samples to `lab`, deferred. A length in a run shorter than 10 (with `T16189`, or split by `16185T`) stays evidence |
| An AC-repeat insertion placed differently by two reads, lost in the consensus (`524.1A 524.2C`) | ADR-0047 §1 amendment: one 3'-most placement whatever the profile score; contested AC sites 118 to 4 |
| An insertion into a repeat vetoed on whichever call alignment took as inserted | SRS-VAR-020: judged on the repeat's placements |
| DNA calls where the reviewers did not analyse counted as defects | the `unanalysed` class |
| HVS-I and HVS-II C-tract and AC-repeat notation, the deletion trade to `16298`, the insertion shift `16172C` blocked | `dna-normalize` window rules |

## Open and TODO

| Topic | Owner | Status |
|---|---|---|
| Consensus to `call`: an `N` at a decided run's bound (`310`) drops a decided length | core | [PROP-0006](../proposals/0006-call-consensus-decisions.md) (accepted); the decided length is now not callable rather than reported as the reference |
| Calls DNA does not trust at the AC repeat (`523`/`524`) | core | open |
| Phantom calls inside the HVS-I C-tract (`16185-16193`) | core | open |
| DNA writes the HVS-II C run `303-309` at the reference length without a label in about 600 samples where the truth has `309.1C` or `309.2C`: the forward read shows the run's 5' end and is masked from `310`, and the only read reaching `310` is untrusted, so no read bounds the run, yet the consensus forms no run site and decides it position by position. The measures cannot see it, since the rules leave indel lengths at `309`, `455`, `463`, and `573` out | core | planned: seed a run site where no trusted read bounds a run while some read disagrees on it, so the length is labelled (`phase_loss`, `reference_frame`) or left undecided; an informational count of special-position length differences in the harness first |
| `call` re-places a consensus `N`, so the variants document asserts a position the consensus left undecided (16 positions, 13 samples on dev v5; `empop_undecided` in the harness) | core | PROP-0006 follow-up: the consensus-to-`call` interface loss, visible through the EMPOP query |
| A read whose signal collapses inside a C-tract (`16189C`) is classed `read_edge` (lab), since DNA marks a phase loss only after a run; it is a technology limit | Sanger callability | TODO, signal work: attribute signal loss within a run |
| `73` lost after the 3'-most placement of a nearby insertion merges stretches (5 samples) | core | known, coverage only |
| Weak calls that read the reference leave `N` | core | parked: trusting them decides run lengths (measured: more regressions than fixes); a reviewer flag instead of a call is an option |
| Consensus stretches grown over tandem-repeat units | core | parked: fewer strict concordances than it gains |
| The limit-of-detection check counts dephased reads' mixtures after a C-tract as a second base | evaluation harness | watch: leaks clustered at `316-320` or `16194-16199` would point to it |
| Heteroplasmy output (IUPAC) | core | deferred |

## Truth re-check

A case where the read evidence supports DNA goes on a re-check list for the
reviewers. DNA is not changed for it. The list names samples, so it is kept
with the run outputs under `results/` on the lab server and never in
repository documents. A corrected truth reaches the harness through the
statistics database, whose changes the panel's truth digests detect. Until
then the case still counts as a leak in an approval run, as pre-registered
([PROP-0005](../proposals/0005-evaluation-feedback-loop.md)).

Candidates on development panel v5 (revision 4e6f8d1):

| Candidate | Samples | Evidence |
|---|---:|---|
| Truth heteroplasmy with no second base in any read | 84 | Clean reads on both strands show one base, and the code's other base stays below 3 % in every covering read. Most are C-tract positions such as `16189Y`. |
| Truth allele where clean reads on both strands show the reference | 2 | Every read of the position reads the reference cleanly. |
| Truth heteroplasmy that matches a neighbour's shadow | 1 | The second base is the neighbouring `C`, on one strand, right after the HVS-II C-tract, where the other read is a dephased mixture. |
| Truth substitutions all one position early | 1 | DNA writes a known haplogroup motif. |
| Truth notation minority form | 22 | Same sequence; for example `514DEL 515DEL` is written by 24 of 58,807 corpus samples with that sequence, against `523DEL 524DEL`. |

## Updating the catalogue

Counts are recomputed from a harness run on the development panel. A change
that fixes a class updates its row and status. A class closed as "not
handled" keeps its reason here, so the next iteration starts from it.
