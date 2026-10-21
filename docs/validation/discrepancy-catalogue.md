---
implements: [SRS-ALN, SRS-CALL, SRS-NOM, SRS-VAR]
decided-by: [ADR-0047, ADR-0077]
---

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
[[ADR-0077]] §5:

- the core caller decides **what** the sample's sequence is;
- `dna-normalize` decides **how** a decided sequence is written, in the lab's
  convention;
- `dna-empop` checks and exports the **format**;
- the [[validation/evaluation-harness|evaluation harness]] classifies.

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
  [[proposals/roadmap#deferred-product-directions|roadmap's deferred directions]]:
  phase recovery, signal denoising, and length-mixture decomposition;
- **the core** (`dna`): the current focus.

## Classification rule

The harness gives each sample one [[validation/evaluation-measures#discrepancy-classes|class]]
and one [[validation/evaluation-measures#discrepancy-causes|cause]], with the evidence
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
([[PROP-0007]]).

## Development panel v5 (revision f6b979e)

Of 16,410 active samples, 12,877 agree with the truth on every measure. Since
revision f6b979e every measure except the primary judges insertions and
deletions at `309`, `455`, `463`, and `573` like any other position; only the
primary keeps the laboratory's exception. A sample's cause is its first
defect, in the order `dna`, `lab`, `undetermined`; a sample with deferred
findings only takes the first of them:

| Class | Samples | dna | lab | undetermined | dna (deferred) | lab (deferred) | technology (deferred) |
|---|---:|---:|---:|---:|---:|---:|---:|
| sequence | 1,125 | 13 | 2 | 1,079 | 0 | 31 | 0 |
| coverage | 1,301 | 84 | 579 | 37 | 396 | 0 | 205 |
| heteroplasmy | 941 | 0 | 84 | 76 | 623 | 0 | 158 |
| unanalysed | 126 | 0 | 0 | 126 | 0 | 0 | 0 |
| representation | 29 | 8 | 21 | 0 | 0 | 0 | 0 |
| failed | 9 | 9 | 0 | 0 | 0 | 0 | 0 |
| not_comparable | 2 | 0 | 0 | 2 | 0 | 0 | 0 |

Under the working hypothesis this splits the work:

- **the core now**: 114 samples (`dna`);
- **the laboratory**: 686 samples (`lab`), most of them no read or weak
  signal at `73`, plus the truth re-check below;
- **a reviewer**: 1,320 samples (`undetermined`), most of them sequence
  differences at the HVS-II C-tract (see *Sequence*);
- **deferred**: 1,413 samples, being the IUPAC output (1,019), the
  reviewers' HVS-I length convention (31), and signal processing (363).

Judging the special positions made 738 more samples discordant on the
`haplotype_lod` measure (corpus rate 0.065 to 0.111), 695 of them at `309`:
DNA's unlabelled HVS-II C-run length (see *Open and TODO*) is now visible.

**Coverage** (a truth allele where DNA has no callable base). By finding:
- `lab`: no read covers it (`no_read`, mostly `73G` at the HVS-II forward
  read's start, an assay and run-quality matter, see
  [[validation/known-limitations|known limitations]]), a weak-signal veto
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

**Sequence** (the sequences differ). 873 of the 1,125 samples differ in the
HVS-II C-tract, most of them the `303-309` run length DNA writes at the
reference length where the truth has `309.1C` or `309.2C`; the rest are in
the HVS-I C-tract (77), the HVS-I read start `16024-16050` (76), `568-576`
(52), and elsewhere. The harness does not yet judge an insertion under this
class, so these go to adjudication (`undetermined`) rather than to a cause.
DNA's defects here are phantom calls inside the HVS-I C-tract (13 samples); 2
samples have clean reads on both strands showing the reference where the
truth has an allele (`lab`, re-check); in 31 every DNA-only call is an HVS-I
C-run length the reviewers do not record (`lab`, deferred,
`convention_not_recorded`).

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
| Indel lengths at `309`, `455`, `463`, and `573` never judged | every measure but the primary judges them |
| HVS-I and HVS-II C-tract and AC-repeat notation, the deletion trade to `16298`, the insertion shift `16172C` blocked | `dna-normalize` window rules |
| A substitution one read decides over a dye blob, the reference's own peak its second (`527A`, `530G` at HVS-III reverse read starts; `16326G`, `16328G` at HVS-I reverse read starts), auto-passed | review rules 23: `reference_second` makes it a check; the decision is unchanged |

## Open and TODO

| Topic | Owner | Status |
|---|---|---|
| Consensus to `call`: an `N` at a decided run's bound (`310`) drops a decided length | core | [PROP-0006](../proposals/0006-call-consensus-decisions.md) (accepted); the decided length is now not callable rather than reported as the reference |
| Calls DNA does not trust at the AC repeat (`523`/`524`) | core | open |
| Phantom calls inside the HVS-I C-tract (`16185-16193`) | core | open |
| DNA writes the HVS-II C run `303-309` at the reference length without a label in about 600 samples where the truth has `309.1C` or `309.2C`: the forward read shows the run's 5' end and is masked from `310`, and the only read reaching `310` is untrusted, so no read bounds the run, yet the consensus forms no run site and decides it position by position. Visible since the haplotype measures judge the special positions (695 samples at `309`) | core | planned after PROP-0006: seed a run site where no trusted read bounds a run while some read disagrees on it; a length the consensus only estimates (`phase_loss`, `reference_frame`) is then not callable in `call` |
| A run length decided by one read while a read on the other strand shows another length (`309.x`, `573`) | core | closed after measurement: no such strand contradiction. The signal came from a review defect: `read_disagrees` left out masked inserted calls, so phase-shifted reads that agree looked one base short (fixed in review rules 3). A core rule that turns such a length into an estimate changed 9 samples (`strict` 0 better, 2 worse; `haplotype_lod` 1 better), or none when the contradicting read had to be in phase. The leaks it targeted are truth lengths both strands disagree with (truth re-check) or a run one strand only lower-bounds |
| `call` re-places a consensus `N`, so the variants document asserts a position the consensus left undecided (16 positions, 13 samples on dev v5; `empop_undecided` in the harness) | core | PROP-0006 follow-up: the consensus-to-`call` interface loss, visible through the EMPOP query |
| A read whose signal collapses inside a C-tract (`16189C`) is classed `read_edge` (lab), since DNA marks a phase loss only after a run; it is a technology limit | Sanger callability | TODO, signal work: attribute signal loss within a run |
| `73` lost after the 3'-most placement of a nearby insertion merges stretches (5 samples) | core | known, coverage only |
| Weak calls that read the reference leave `N` | core | parked: trusting them decides run lengths (measured: more regressions than fixes); a reviewer flag instead of a call is an option |
| Consensus stretches grown over tandem-repeat units | core | parked: fewer strict concordances than it gains |
| An edge call whose own signal takes its base off the edge (SRS-CALL-017) can still read a base one position off at `310` (`310C`, 2 samples on dev v5, both inside check items) | Sanger evidence | watch: a recurring false call outside check items, or a rise in samples with an unreviewed difference, would call for a stricter rule |
| A forward HVS-II read that ends in the `568-573` C run, with its calls before the run vouched for against their end's norms (SRS-CALL-017), lets the consensus decide the run's site from the reverse HVS-III read, which reads it longer: `573.1C`-`573.3C` (2 samples on dev v5, inside check items, both strands showing the longer run past the truth's HVS-III interval end at `573`). Behind an ambiguous bound at `567` this also gave `567C` in 5 samples, which the bound guard removed | Sanger evidence, or the truth | re-check the two samples; a false call outside check items would call for a core rule on a run a read loses phase right after |
| The limit-of-detection check counts dephased reads' mixtures after a C-tract as a second base | evaluation harness | watch: leaks clustered at `316-320` or `16194-16199` would point to it |
| A one-read substitution under a four-channel signal spike, roughly eight times its neighbours' peaks (`186T` in one earlier sample: the spike's `T` reads cleanly and the reference `C` vanishes, so no second peak flags it; the other strand's `C` is set aside as weak) | Sanger evidence | TODO, signal work: mark a call inside such a spike as not a base |
| Heteroplasmy output (IUPAC) | core | deferred; a minimal two-strand caller was measured and gains little (most truth heteroplasmy is not shown on both strands at 10 %) |
| Sequence-class causes for insertions | evaluation harness | open: insertions go to adjudication |

## Truth re-check

A case where the read evidence supports DNA goes on a re-check list for the
reviewers. DNA is not changed for it. The list names samples, so it is kept
with the run outputs under `results/` on the lab server and never in
repository documents. A corrected truth reaches the harness through the
statistics database, whose changes the panel's truth digests detect. Until
then the case still counts as a leak in an approval run, as pre-registered
([[PROP-0005]]).

Candidates on development panel v5 (revision 4e6f8d1):

| Candidate | Samples | Evidence |
|---|---:|---|
| Truth heteroplasmy with no second base in any read | 84 | Clean reads on both strands show one base, and the code's other base stays below 3 % in every covering read. Most are C-tract positions such as `16189Y`. |
| Truth allele where clean reads on both strands show the reference | 2 | Every read of the position reads the reference cleanly. |
| Truth heteroplasmy that matches a neighbour's shadow | 1 | The second base is the neighbouring `C`, on one strand, right after the HVS-II C-tract, where the other read is a dephased mixture. |
| Truth substitutions all one position early | 1 | DNA writes a known haplogroup motif. |
| Truth run length both strands contradict | 2 | `2510457`: both reads show `309.2C` and `315.1C` (the forward read's masked calls included), the truth writes `309.1C`. `2510987`: the reverse read shows nine C then the A at `574` in phase (`573.3C`), the forward read loses phase inside the run and only lower-bounds it; the truth writes none. |
| Truth notation minority form | 22 | Same sequence; for example `514DEL 515DEL` is written by 24 of 58,807 corpus samples with that sequence, against `523DEL 524DEL`. |

## Updating the catalogue

Counts are recomputed from a harness run on the development panel. A change
that fixes a class updates its row and status. A class closed as "not
handled" keeps its reason here, so the next iteration starts from it.
