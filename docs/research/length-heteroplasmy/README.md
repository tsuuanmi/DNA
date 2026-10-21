---
id: RES-length-heteroplasmy
status: active
summary: How forensic practice records C-tract length heteroplasmy, how precise DNA's tract-length review items are per evidence state, and whether decode noise should be calibrated per laboratory batch.
---

# Research: C-tract length heteroplasmy

**Status:** active — exploratory evidence, not production authority. Nothing
here changes a requirement, a decision, or DNA's output. It informs
[[PROP-0012]] (the trace crate and its length decision) and the draft
run-context proposal, `PROP-0014`, which is under review in
[PR #212](https://github.com/tsuuanmi/DNA/pull/212) and not yet on `main`.

## Question

The C-tracts at HVS-I 16184–16193 (16189C), HVS-II 303–315 (309.1C, 309.2C,
315.1C) and HVS-III 568–573 carry length heteroplasmy: a population of
molecules whose runs differ in length. After the tract a Sanger read dephases.
DNA decodes the length mixture ([[design/length-mixture]]) and the consensus
decides a run's length only when both strands agree
([[design/length-mixture#decision-increment-4|rule A]], [[ADR-0074]] step 5).
Every other tract-length run becomes a review item.

Three questions follow:

1. **Practice.** How do forensic and clinical laboratories record, resolve,
   and validate these lengths, and how firm is the ground truth?
   ([Field practice and methods](#field-practice-and-methods))
2. **Tiers.** When DNA's evidence for a tract length is of a given kind, how
   often is the length wrong? Could some review items become non-blocking
   notes? ([Tract-length review tiers](#tract-length-review-tiers))
3. **Run calibration.** Does the decode noise that rule A allows for vary per
   laboratory batch, enough to change a decision?
   ([Run-calibrated length noise](#run-calibrated-length-noise))

## Summary

- **The truth is a convention near a tie.** For a clear dominant length, the
  reviewer's length is a measurement. For near-equal mixtures it is a
  convention ("shortest dominant"), and it also moves with primer, platform,
  and PCR. Every guideline drops these indels from comparisons and searches.
- **One kind of evidence is precise; most are not.** On development panel v5,
  a reported HV2 length that one strand reads in phase, with every view on
  both strands agreeing (class A), was never wrong (0 of 458 items; 95 %
  upper bound 0.65 %). Decode-only, broad-mixture, and `long_run` items are
  wrong in 1–50 % of runs, almost always one C too long. No class can yet be
  shown at the auto-pass level (about 0.01 %) on development data alone.
- **Per-batch noise is real but inert.** HVS-II run 1's decode noise varies
  between laboratory batches (between-batch sd 0.31 of the corpus
  dispersion). Calibrating τ per batch moves it by at most about ±35 % and
  changed no output on development panel v5: rule A's decisions sit far from
  the τ boundary.

## Field practice and methods

Literature only, 2026-10-21; no code and no measurement. Primary sources were
fetched where possible. Items marked *(abstract only)* were not read in full;
claims marked *unverified* were not confirmed against a primary source.

### Recording conventions

**ISFG, Parson et al. 2014** (FSI: Genetics 13:134–142;
[PDF](https://next.isfg.org/documents/36/Parson-2014-mtDNA.pdf)), section 3.3.2
and recommendations 10, 11, and 16:

- Length heteroplasmy (LHP) appears when a homopolymer exceeds about 8
  identical bases.
  - 16189C gives an uninterrupted C10, with variants up to C14 and
    "generally uninterpretable" sequence downstream. It occurs in 14 % of West
    Eurasian EMPOP sequences and 36 % in haplogroup B / East Asia.
  - HVS-II insertions after 302–310 occur in 74 % of 1,887 sequences.
  - T310C joins C7 and C5 into C13, but the tract "usually tends to be
    contracted down to nine Cs".
  - The 568–573 tract is LHP-prone and often phylogenetic (haplogroup I;
    0.3 %).
- **Casework:** there is "no routine practice". Some labs skip the 10 bases
  from 16184; others report all variants. Indels in these tracts "are ignored
  in direct forensic comparisons and database searches". LHP and point
  heteroplasmy differences never exclude (recommendation 10).
- **Databases** report the **dominant variant (major molecule)**
  (recommendation 11). It is read from the highest representation of a single
  non-repetitive base downstream, a **signpost**:
  - forward: G16196 (HVS-I), T310 and G316 (HVS-II), G577 (HVS-III C-tract),
    G527 (AC repeat);
  - reverse: G16156, G275, G564, G513.
- **Ties:** "where two or more LHP variants are equally dominant, it has
  become practice to determine the **shortest dominant** LHP variant (e.g.
  EMPOP)".
- Searches disregard length variants in homopolymers (recommendation 16).
  EMPOP does this by default for the C-tracts around 16189, 310, and 568–574,
  and for the T insertions near 452 (section 5.4.2).

**SWGDAM 2024** (effective 2024-09-09, supersedes 2019;
[guidelines](https://www.swgdam.org/_files/ugd/4344b0_9784320e1188491fbdd0dbf040ca69ae.pdf),
[supplement](https://swgdam.org/_files/ugd/4344b0_d373f9cdb74341798799ccf6d9cb2813.pdf)):

- §3.7: labs must write their own rules to confirm and report length
  variants. §3.8: validation must characterise how the software handles
  repeats.
- §5.1.1: "SWGDAM does not recommend the use of such length variants for
  comparisons or searches", meaning those after 16193, 309, and 573.
- FAQ-8: "In most cases laboratories do not report the presence of length
  heteroplasmy in these regions or use it for searching and comparison. By
  default, the CODIS and EMPOP databases ignore indels following positions
  16193 and 309."
- FAQ-7: report heteroplasmy only when it is seen in both the forward and the
  reverse reaction. LHP is found in 63.6 % of individuals' control regions
  (Just 2015) and is attributed to polymerase slippage ("akin to stutter").
- FAQ-4: "additional Sanger sequencing primers" are recommended for the
  16183–16194 and 302–310 tracts.
- Alignment rule 4.2.3: 310C when T310 is absent; 311T when T310 is followed
  by a second T.

**EMPOP query engine**
([documentation §4.1.6](https://empop-use.readthedocs.io/en/latest/using_empop/)):

- Standard queries disregard the length-hotspot windows: 16193 (16189–16194),
  309 (302–316), 573 (567–574), plus 455, 463, 524, and coding-region sites.
- **315.1C is "a phylogenetically stable insertion", not LHP.** It must still
  be called.
- 309.2C cannot be written without 309.1C.

**Laboratory procedures**, which show the range of practice:

- **FBI legacy procedure**
  ([BIO-575](https://fbilabqsd.fbi.gov/file-repository/dna/bio-575-00-interpretation-of-legacy-mitochondrial-dna-data.pdf)):
  - HV1 length variants are "not recorded", and the sequence is truncated to
    the C5TC4 format.
  - HV2: "the major length variant as called by Sequencher" is listed.
  - HV3 length variants are "not noted" and are usually single-stranded.
  - A length difference "is not considered as a base difference" in
    comparisons.
- **NYC OCME**
  ([Sanger](https://www.nyc.gov/assets/ocme/downloads/pdf/technical-manuals/forensic-biology-technical-manuals/sequence_nomenclature_and_alignment_080621.pdf),
  [MPS 2025](https://www.nyc.gov/assets/ocme/downloads/pdf/technical-manuals/forensic-biology-technical-manuals/interpretation_guidelines_mtdna_mps_092625.pdf)):
  - LHP is noted as present; the HVI C count is not recorded.
  - The HVII "major type as determined by the analyst" is used.
  - T310C is truncated to C7TC5.
  - In MPS, "the specific number of C nucleotides will not be considered".
  - Exclusion means at least two differences "excluding length heteroplasmy".

### Methods

| Method | Resolves | Needs | Accuracy and limits | Fit for DNA |
|---|---|---|---|---|
| Direct Sanger, signpost reading ([Parson 2014](https://next.isfg.org/documents/36/Parson-2014-mtDNA.pdf)) | dominant length per strand | existing traces | clear in most samples; near-equal mixtures are ambiguous and primer-dependent: GEDNAP labs split 22:7 between 309.2C and 309.1C ([Berger 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3064993)) | what DNA already does; signposts give a cheap independent cross-check |
| Opposite strand or extra primers (SWGDAM FAQ-4; Parson 2014 §2) | sequence downstream of the tract, in phase from the other side; single-strand artefacts after C-tracts are fixable "with data from either the complementary strand or a different primer" | reverse read; an internal primer for HV1 (for example H16167, or A4-B4 at OCME) | the tract length itself is still read at a mixed exit from both sides | strands already in hand; extra primers would be a laboratory change |
| Cloning and sequencing (Bendall & Sykes 1995 *(abstract only)*; [Stewart 2001](https://ojp.gov/ncjrs/virtual-library/abstracts/length-variation-hv2-human-mitochondrial-dna-control-region)) | per-molecule length distribution | cloning, many colonies | laborious ("not of practical use in casework"); PCR slippage still present | validation only |
| Gel densitometry (Bendall & Sykes 1995 *(abstract only)*) | variant proportions | slab sequencing gels | historical; maternal relatives share proportions | none |
| Fluorescent amplicon sizing ([Prieto 2004](https://investigacion.usc.gal/documentos/5da9762629995222b484fdbe?lang=en) *(abstract only)*; [Berger 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3064993)) | ratio of all length variants, 1-bp resolution | labelled primers, proofreading polymerase (KOD), CE sizing | finds the dominant variant in mixtures up to 1:0.95; confirmed every clear sequence-derived dominant; found at least 6 variants at 573 and favoured 309.2C at 1:0.75; still ties at 1:1:0.6 (buccal) | **best orthogonal truth** for a re-check subset; runs on the same 3500 CE |
| RFLP (MfeI, Lutz-Bonengel 2004; [BsaXI, Aral 2006](https://pmc.ncbi.nlm.nih.gov/articles/PMC1388229/)) | C7 against not-C7 at 303–309 | digest | binary; no length or ratio | no |
| DHPLC | point heteroplasmy (0.5 % limit of detection) | HPLC | no length-tract validation found | no |
| Pyrosequencing / 454 ([Holland 2011](https://www.nist.gov/mml/bmd/genetics/upload/CroatMedJ_52_0299.pdf)) | mixtures, point heteroplasmy | 454 | "inability ... to fully resolve homopolymeric stretches and length heteroplasmy" (more than 4 C) | no |
| MPS, Illumina or Ion ([Lin 2017](https://researchnow.flinders.edu.au/en/publications/investigation-of-length-heteroplasmy-in-mitochondrial-dna-control/) *(abstract only)*; [Just 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4550493); [Sturk-Andreaggi 2020](https://ojp.gov/ncjrs/virtual-library/abstracts/impact-sequencing-method-detection-and-interpretation-mitochondrial) *(abstract only)*) | read-length distribution per tract; LHP in tracts longer than 7 | library and sequencing; count only reads spanning the whole tract and flanks | the major molecule differs from Sanger in 37 % (Ion) and 20 % (Illumina) of 30 regions; Sanger majors are longest, Ion shortest; tool-dependent | not available; would be a second validation source |
| MPS tools: Converge ([Thermo](https://www.thermofisher.com/order/catalog/product/A35131)), GeneMarker HTS, AQME (CLC), MitoSAVE, [mtDNA-Server 2](https://pmc.ncbi.nlm.nih.gov/articles/PMC11223783/) | Converge quantifies point and length heteroplasmy to about 10 % at more than 100×; GeneMarker HTS aligns by motif; AQME applies forensic conventions; MitoSAVE converts VCF to forensic notation; mtDNA-Server 2 has no LHP-specific handling (Mutect2 indels) | as MPS | none publishes a Sanger-comparable length rule | naming and convention reference only |
| Sanger trace deconvolution: TIDE ([Brinkman 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC4267669/)), ICE, DECODR, Tracy/Indigo `decompose` ([Rausch 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7071639)); Poly Peak Parser and CRISP-ID *(unverified: from memory, not re-checked)* | indel spectrum after a breakpoint (shifted copies of one sequence, non-negative fit) | a reference or control trace | net indel size is estimated well; sequence-level deconvolution varies across tools ([Aoki 2024](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10854981/)); fixed allele limits (Tracy: at most two alleles, grid fraction) | the same model family as [[PROP-0010]] and [[PROP-0012]] |

**Gap.** No published method was found that estimates an mtDNA homopolymer
length distribution from the post-tract phase shift of Sanger traces, as a
survival curve or a mixture decode. The nearest work is CRISPR trace
decomposition and the 1995 densitometry. The forensic field treats peak
heights as semi-quantitative: dye and context effects mean point-heteroplasmy
ratios "typically cannot be reliably estimated" from peak heights. That
statement comes from a search snippet attributed to EDNAP and is
*unverified*. DNA's decoder therefore has no external accuracy benchmark.

### How firm the ground truth is

- **Inter-laboratory exercises.**
  - GEDNAP split 22:7 on 309.1C against 309.2C for one sample, depending on
    primer ([Berger 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3064993)).
  - EDNAP 2026 (Kiesler, Romsos, Vallone, Parson, FSI: Genetics 86;
    [NIST](https://www.nist.gov/publications/ednap-collaborative-study-mtdna-heteroplasmy-detection-using-sanger-and-massively)):
    LHP interpretation varied "substantially" more than point heteroplasmy.
    Sanger LHP could be harmonised by review under the ISFG guidelines, so it
    is convention-driven. MPS LHP differed between and within technologies.
- **Platform.** The major molecule changes with platform in 20–37 % of
  regions (Sturk-Andreaggi 2020).
- **In vitro.** C-stretch patterns were inconsistent on repeat PCR under
  identical conditions. More cycles or more template gave more LHP, and
  monoclonal templates became heteroplasmic
  ([Seo 2010](https://ojp.gov/ncjrs/virtual-library/abstracts/alterations-length-heteroplasmy-mitochondrial-dna-under-various)
  *(abstract only)*). Stutter and true LHP therefore cannot be separated in
  Sanger signal.
- **In vivo.** Hairs from one person differ in C count and proportions
  (Stewart 2001; an EDNAP hair study, 302–309). The C-tract mutation rate is
  up to 6 % per generation, with an instability threshold of 8 C
  ([Forster 2010](https://research-portal.st-andrews.ac.uk/en/publications/evaluating-length-heteroplasmy-in-the-human-mitochondrial-dna-con/)
  *(abstract only)*).

**Conclusion.** For clear dominants the truth is a measurement. For
near-equal mixtures it is a convention (signpost reading plus "shortest
dominant"), and it varies with primer, platform, and PCR.

**DNA's own data agree.** A two-strand decoder spike (2026-10-19) found 53
confident disagreements on HVS-II run 1. In all of them the truth is the
mode − 1, in the zone where the shorter neighbour has at least 0.6 × the
mode's share. That matches the published "shortest dominant when (nearly)
equally dominant" practice: these look like the convention, not reviewer or
decoder error.

### Implications for DNA

- **A note, not a decisive difference.** Every guideline drops LHP indels at
  16193, 309, and 573 from comparisons and searches; LHP presence is
  note-level information.
  - **315.1C is not in that class.** It is a stable call and stays a normal
    call.
  - **T310C** gives a contracted (about C9) mixture, not C13.
- **Dominant or shortest.** Database reporting uses the dominant variant;
  when variants are equally dominant, EMPOP takes the shorter. "Equally" has
  no published numeric cut. Parson calls these cases rare; DNA's mode − 1
  zone suggests the reviewers here apply the rule more widely.
- **HVS-I 16189C.** Major labs do not record a count: the FBI and OCME
  truncate to C5TC4 with 16189C, and EMPOP and CODIS ignore 16193.x. The
  laboratory truth here leaves the run unrecorded too. Asserting
  16193.1C/.2C adds risk for little value.
- **Both strands.** SWGDAM requires heteroplasmy to be seen in both reactions
  before it is reported.

**Candidate directions** (none adopted):

1. **Keep the convention out of the decoder**, as [[PROP-0012]]'s non-goals
   already say. An optional, labelled reporting policy in `dna-normalize`
   ([[design/variant-nomenclature#report-conventions|report conventions]])
   could report the dominant length, or the next-shorter one when its share
   is at least ρ × the dominant's (ISFG 2014 §3.3.2). ρ would come from
   validation with provenance, never from fitting reviewers. Without a
   policy, the output lists the components.
2. **LHP as a structured note** for 16193, 309, and 573: the tract, the
   components (length and signal share), strand agreement, and an
   ignored-in-search marker for the EMPOP section. Escalate to a check only
   when the dominant length itself is undecided (a near-tie without a
   configured policy, or strands that disagree), never on a minor share
   alone. Shares stay signal shares, not heteroplasmy levels (SRS-VAR-009).
3. **Signpost cross-checks** as decoder-independent evidence of the dominant
   length (G16196, T310/G316, G577 forward; G16156, G275, G564 reverse), with
   a flag when the signpost disagrees with the decoder's dominant component.
4. **HVS-I 16189C** following FBI and OCME practice: note the LHP, record
   16189C, and assert no 16193.x count unless both strands decode it
   unambiguously. Low return: the truth is unrecorded and searches ignore it.
5. **An orthogonal truth subset** by fluorescent amplicon sizing (Berger
   2011: labelled primers around 16189, 310, and 573, KOD polymerase, the
   existing 3500 CE), about 100–200 samples enriched for HVS-II run-1
   near-ties and the mode − 1 zone. It would calibrate the decoder's shares,
   set ρ, and settle whether the reviewer pattern is convention or signal
   without trusting reviewer records.

## Tract-length review tiers

2026-10-21. Data for one user decision: when the evidence for a tract length
is of a given kind, may its review item become a non-blocking **note**, or
must a reviewer confirm it?

### Data and method

- **Run.** Development panel v5
  ([[validation/evaluation-harness#panel-versions|panel versions]]), DNA at
  `main` `60252579` (Sanger plugin 22, review 25), with the evidence cache.
  Truth: the laboratory's statistics database and the frozen re-check list.
  Periods split at `batch_date` 2026-01: current is 10,695 samples, earlier
  5,715. Retired and failed samples are left out.
- **Evidence.** The review's own length judgement over every anchored window,
  confirmed runs included, was dumped through a research-only side channel.
  The dump run's review and evaluation tables are byte-identical to the clean
  run, and a reimplementation of the unreviewed blocks reproduces the
  harness's 715 blocks with the same classes.
- **Errors.** An error is a haplotype difference block, the unit the gate
  uses ([[validation/evaluation-measures#measures|measures]]). Its class is
  re-check, then truth `N`, then heteroplasmy, then non-heteroplasmy. A block
  counts against an item when it overlaps it.
- **Comparison.** DNA's reported run lengths come from its notation calls
  when the window is callable, and from the review's `called` otherwise. The
  truth's run lengths come from its alleles where its intervals cover the
  window.
- **Bounds.** "95 % upper" is the one-sided Clopper–Pearson bound.
- **"Near zero"** is today's current-period auto-pass rate: 1 non-het block
  in 7,728 clear samples; 3 in 155,946 passed zones; 0 in 8,721 passed
  window@HV2 zones; 1 run-1 length differing from the truth in 8,716
  confirmed without a check (0.011 %).

### The tract items

A tract check has one of two shapes:

- **DNA reports a length.** The window is callable, so a wrong length would
  be an unflagged error.
- **DNA reports none.** The run is `estimated_length` (`phase_loss` or
  `reference_frame`). No silent error is possible, but the result lacks a
  length.

| Tract | Period | Check items | Length labels only¹ | DNA reports | DNA reports none | Non-het blocks (length-only) |
|---|---|---:|---:|---:|---:|---:|
| HV2 303–315 | current | 2,023 | 1,790 | 895 | 895 | 18 (all where DNA reports) |
| HV2 303–315 | earlier | 1,103 | 987 | 509 | 478 | 23 + 1 |
| HV1 16184–16193 | current | 442 | 145 (`estimated_length`) | 0² | 145 | 19 |
| HV3 568–573 | current | 235 | 91 (`estimated_length`) | 0 | 91 | 0 |

¹ `length_unconfirmed`, `estimated_length`, or both. `long_run` is its own
check and is kept apart; in HV2 it holds 19 of the 29 current non-het items.
² No HV1 window is callable throughout, yet in 19 of them DNA still reports a
C-run length ([HV1 and HV3](#hv1-and-hv3)).

### HV2 runs by evidence state

Run-level counts over all samples, checked or not. The states are research
sub-states of the review's length `agreement`:

- `one_strand_clear_call_in_phase`: one strand gives DNA's length clearly,
  at least one of its reads in phase;
- `one_strand_clear_call_decode`: the same, from decoded or derived views
  only;
- `near_equal_pair_with_call`: a near-equal pair (#200) of which DNA's call
  is a member;
- `not_clear_tops_call`: every read tops at DNA's length, not all clearly;
- `not_clear_mixed`: the reads top at different lengths.

"= convention" counts truths equal to the #200 suggestion, the shorter of a
near-equal pair. "Estimate ≠ truth" is, for unreported runs, the review's
`called` against the truth.

**Current period:**

| Run | State | Runs | Reported | Reported ≠ truth | = convention | 95 % upper (rest) | Not reported | Estimate ≠ truth |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 1 | confirmed, no check | 8,717 | 8,716 | 1 | 0 | 0.05 % | 0 | – |
| 1 | `one_strand_clear_call_in_phase` | 300 | 298 | **0** | 0 | 1.0 % | 2 | 0 |
| 1 | `near_equal_pair_with_call` | 525 | 189 | 8 | 4 | 4.8 % | 334 | 17 |
| 1 | `one_strand_clear_call_decode` | 310 | 305 | 9 | 0 | 5.1 % | 5 | 0 |
| 1 | `not_clear_tops_call` | 79 | 66 | 27 | 0 | 52 % | 8 | 1 |
| 1 | `not_clear_mixed` | 479 | 2 | 0 | 0 | – | 406 | 365 |
| 1 | `structure_differs` | 115 | 41 | 1 | 0 | 11 % | 51 | 24 |
| 1 | `one_strand_clear_other` | 103 | 1 | 1 | 0 | – | 87 | 83 |
| 1 | confirmed (in a check item) | 34 | 15 | 1 | 0 | 28 % | 18 | 0 |
| 2 | confirmed, no check | 8,715 | 8,715 | 0 | 0 | 0.03 % | 0 | – |
| 2 | `one_strand_clear_call_in_phase` | 1,300 | 519 | **0** | 0 | 0.58 % | 781 | 0 |
| 2 | confirmed (in a check item) | 469 | 390 | 0 | 0 | 0.77 % | 79 | 0 |
| 2 | `one_strand_clear_other` | 52 | 2 | 2 | 0 | – | 36 | 36 |

**Earlier period**, the same shape: run 1 in phase 1 wrong of 135 (a case on
the re-check list); run 1 decode-only 18 of 151; `not_clear_tops_call` 10 of
35; near-equal pair 4 of 133 (1 convention); run 2 in phase 0 of 323.

- **Run 2** is never wrong when one strand reads it in phase: 0 of 2,030 runs
  over both periods (95 % upper 0.15 %). Every residual error is at run 1.
- **The residual errors have one shape:** DNA is one C longer than the truth.
  - *Decode-only, one strand (HV2F).* Where DNA called the longer of the top
    two (`next` = call − 1, shares 0.59–0.69), 28 of 256 are wrong; where it
    called the shorter (`next` = call + 1), 1 of 201.
  - *Broad mixtures.* When the top share is below 0.5, three or more lengths
    are present (for example C9 at share 0.36–0.38 where the truth is C8). A
    decode margin of `clear` does not protect here.
  - *The near-equal pair.* With DNA's call the shorter member (the
    convention) and no `long_run`, the truth is that member in 310 of 312
    reported runs; with DNA's call the longer member, the truth is the
    shorter in 5 of 7. The convention is right where it fires; it cannot
    handle broad mixtures, where the truth lies below the pair.
  - *`long_run` at C10.* A forward decode and a reverse derived run
    "confirm" C10 at shares 0.40–0.46 where the truth is C9. `long_run`
    catches these.

### HV2 item classes

An item is in a class when every run of its window satisfies it. The classes
use only states and fields the review already computes.

- **A, in phase, all agree.** Every run DNA calls has at least one `in_phase`
  contribution equal to DNA's length; every contribution on either strand (in
  phase, decoded, derived) tops at that length; no structure differs. This is
  the rule of #204.
- **B, near-equal pair at the floor.** Every run is A-like or a near-equal
  pair whose shorter member is DNA's call, and at least one run is a pair.
- **C, floor reading.** For every read, its top and `next` are both in
  {c, c + 1}, where c is DNA's call, and at least one read gives c.
- **D, the rest:** one-strand decode-only reads whose call is the longer of
  the top two, and broad mixtures.

**Items with `length_unconfirmed` only, where DNA reports the length:**

| Class | Current items | Non-het | Length ≠ truth | Earlier items | Non-het | Length ≠ truth |
|---|---:|---:|---:|---:|---:|---:|
| A in phase, all agree | 308 | 1¹ | **0** | 159 | 0 | 1 (re-check, C6 against C7) |
| B near-equal pair at the floor | 59 | 0 | 0 | 33 | 1 | 1 (−1) |
| C floor reading | 150 | 2 | 2 (−1) | 100 | 1 | 1 (−1) |
| D the rest | 326 | 14 + 1 re-check | 15 (12 × −1, 1 convention) | 179 | 20 | 20 (−1) |
| `structure_differs` | 42 | 1 | 1 | 37 | 1 | 1 |

¹ A shift block from outside the tract; the item's runs equal the truth and
other checks cover the block.

95 % upper bounds on the length-error rate, both periods together:

- A: 0 in 458 truth-covered items (0.65 %), or 1 in 459 if the re-check case
  counts (1.0 %);
- B: 1 in 92 (5.1 %);
- A + B + C: 4 in about 800 (1.1 %).

None of these can be shown near the auto-pass level (about 0.01 %) on
development panel v5 alone. **Readiness:** the identical rule (#204, one
strand in phase confirms) failed the gate on the readiness-v2 hold-out panel,
with two window@HV2 zone escapes in 2024-11 batches and 8 newly passed non-het
samples. A note that removes the same checks moves the same blocks out of
review, so that result applies unchanged.

**Items DNA does not report** (`estimated_length`, with or without
`length_unconfirmed`; 895 current, 478 earlier) have no blocks by
construction. The question is what a note could carry.

- **DNA's estimate is the wrong reading.** The review's `called` is the
  phase-loss estimate, and it stops short: it equals the truth in only 41 of
  406 current unreported `not_clear_mixed` run 1s (10 %). Over the D items,
  the truth is +1 in 95 and +2 in 310 of 663.
- **The reads' floor is the right reading.** The minimum over every read of
  its top and `next` matches the truth in 354 of those 406 run 1s (87 %).
- **The near-equal pair works.** Its shorter member matches in 112 of 117 B
  items.
- **The laboratory records HV2 lengths.** The truth covers the window in at
  least 97 % of these items, so passing them drops a length the lab reports
  today. That is a reporting change, not an accuracy question.

### HV1 and HV3

**HV1 16184–16193.** The window's run lengths are already context
(`length_not_recorded`); its checks come from other labels.

| Labels | Current | Non-het | Het | Earlier | Non-het |
|---|---:|---:|---:|---:|---:|
| `estimated_length` | 145 | 19 | 24 | 76 | 2 |
| `undecided` | 208 | 8 | 3 | 143 | 2 |
| `long_run` (± `reference_second`) | 44 | 3 | 32 | 30 | 0 |

- All 19 current non-het `estimated_length` cases are length-only: DNA
  reports C10–C12 where the truth is C9 (16183C 16184A 16189C), from decodes
  that read the longest length at a share near 0.33. DNA is too long.
- The convention does not cover them: the truth's C9 is below the profile's
  `minimum_run` of 10 for `hvs1_c_length`, so the laboratory does record it.
- The heteroplasmy cases are mostly 16192 T/Y.
- **Verdict:** no note candidate in HV1. These are DNA length errors in broad
  mixtures, plus heteroplasmy; they stay checks.

**HV3 568–573.** There is no anchored window, so the review gives no length
agreement.

- **`estimated_length`** (91 current, 61 earlier; DNA reports no length): no
  blocks. The truth's HV3 range ends at or before 573 in 149 of 152 items,
  against 2.5 % of all samples: the laboratory cuts its range exactly where
  the tract is unreadable, and records no length. The forward decodes read
  long, broad mixtures (C10–C11 at shares 0.31–0.37). This is HV1's
  `not_recorded` situation in all but name.
- **`long_run`** (10 current, 11 earlier): every one differs from the truth.
  DNA reports C9; the truth gives the reference C6 with a range cut at 573.
  That is a truth-interval edge problem, which the evaluation now counts as
  truth `N` (#211), not a note question.

### Simulated policies

Each policy is simulated by removing its note items from the checks. Current
baseline: 7,728 clear samples of 10,695 (72.26 %) and 1,974 window@HV2
checks.

| Policy | Note items (current) | Newly clear (current / earlier) | Clear, current | New unreviewed blocks (current / earlier) | Truth differences moved into notes (current) |
|---|---:|---:|---:|---|---|
| P1 HV2 class A | 308 | +96 / +43 | 73.16 % | none / 1 re-check | 0 |
| P2 P1 + B | 367 | +149 / +63 | 73.65 % | none / 1 non-het | 0 |
| P3 P2 + C | 517 | +255 / +112 | 74.64 % | 2 non-het / 2 non-het, 1 re-check | 2 (−1, not convention) |
| P4 P1 + unreported A/B (coverage) | 447 | +208 / +87 | 74.20 % | none / 1 re-check | 0 reported; the reading differs in 5 (truth +1) |
| P5 every length-only HV2 item | 1,790 | +1,209 / +506 | 83.56 % | 17 non-het / 24 non-het, 3 het | 17 reported (1 convention); 439 unreported estimates differ |
| HV3 `estimated_length` only | 91 | +11 / +3 | 72.36 % | none | 0 |
| P1 + HV3 | 399 | +141 / +62 | 73.58 % | none / 1 re-check | 0 |

Some newly clear samples already carry an unreviewed block in another zone:
under P1, 3 re-check and 1 het (current) and 4 het and 2 re-check (earlier);
P3 and P5 add non-het ones (3 and 14, current). Convention differences hardly
move: 4 of 5 current ones sit in `long_run` items, which stay checks.

### Candidate policy

For the user's decision; not adopted.

- **Keep as checks** every category with DNA-length errors of 1–50 %, almost
  all "DNA = truth + 1": one-strand decode-only runs whose call is the longer
  of the top two, `not_clear` runs (broad mixtures, top share below 0.5),
  `structure_differs`, `disagree`, `long_run`, and every HV1 tract item.
- **Note tier, HV2 reported: class A.** The condition uses no positions and
  no new thresholds: an anchored window whose lengths the profile records;
  DNA reports it (callable, no `N`); for every run DNA calls, at least one
  `in_phase` contribution equals DNA's length and every contribution on
  either strand tops at it; no structure differs; no other check label. On
  development panel v5: 0 length errors in 458, +96 current clear samples,
  −308 window@HV2 checks. It is **not adoptable on development evidence
  alone**, because the same rule failed readiness-v2. It needs either the
  user opening readiness-v2 to read the two escapes (making it development
  data), or validation on the planned hold-out set of about 20,000 samples.
  Until then, A is a *fast check*: the review shows the in-phase evidence
  and the reviewer confirms in one step.
- **Convention-consistent fast check: class B** (1 error in 92, too high for
  a note). Show the pair and the convention reading for a one-step accept.
- **Coverage note, DNA reports no length.** A reporting decision, not an
  accuracy one.
  - *HV3 568–573:* the laboratory records no length (it cuts the range at 573
    in 98 % of these items). A profile convention like HV1's
    (`not_recorded`, run length over 568–573, base C, a `minimum_run` chosen
    from data) would make `estimated_length` there context: +11 current
    clear samples, no blocks.
  - *HV2:* the laboratory does record the length, so a coverage note drops a
    reported field. P4 gives +112 current clear samples over P1 with no
    blocks. If accepted, the note carries the reads' floor reading, never
    DNA's phase-loss estimate (wrong for 90 % of unreported
    `not_clear_mixed` run 1s).
- **What a note carries.** The `length_unconfirmed` `runs` payload
  ([[reference/review#labels-rules-version-25|review labels]]): per run,
  `run`, `base`, `called`, `agreement`, `convention`, and the reads (`read`,
  `strand`, `source`, `length`, `share`, `next`, `margin`). Added: each
  decode contribution's `first_z` and `total_z`, DNA's run
  `length_evidence`, the reads' floor reading for unreported runs, and the
  tier reason (for example `assured: "in_phase"`). The pattern is
  `possible_heteroplasmy` with `explained`: the label becomes context and the
  item kind becomes `note`, so no new label is needed.
- **How the evaluation would count note-tier differences.** A note is not a
  review: a block overlapped only by note items stays unreviewed and is gated
  by its class, non-het included, as a passed zone's block is. A `notes.tsv`
  would list the note items, and a `tier` column in `unreviewed.tsv` would
  separate `note` from `pass`. One exception: a block whose truth equals the
  note's `convention` reading is classed `convention` (reported, not gated).
  Exempting note blocks wholesale would hide exactly the errors this tier
  must not make.

## Run-calibrated length noise

2026-10-21. Rule A allows for a replicate noise τ per strand and run,
estimated once from 743 re-sequenced samples
([[design/length-mixture#decision-increment-4|Noise τ]], SRS-CONS-013).
#204 (one-strand in-phase confirmation) was safe on development panel v5 but
added two window@HV2 escapes in two 2024-11 batches of the hold-out panel.
**Hypothesis:** decode noise and in-phase reliability vary per run, so a
fixed τ is too optimistic for some runs.

### Data

Measured from clean committed trees only (`main` `60252579`), development
panel v5 and the frozen re-check list. The hold-out panel was not read or
run.

- **Run identity.** ABIF tags of all 66,211 development traces (`RunN1`,
  `RUND1`/`RUNT1`, the `HCFG4` serial, `MODF1`): 2,391 instrument runs over
  1,000 batches and 10 instruments, one mobility file
  (`KB_3500_POP7_BDTv3`), and `FastSeq50_POP7xl` on all but 35 traces.
  HV2F and HV3R of one sample sit on different instrument runs for 45 % of
  samples (runs are per primer); a development batch holds a median of 3
  development samples per instrument run.
- **Context, truth-free.** Every sample of every development batch folder
  (87,607 samples, 175,483 HV2F/HV3R traces; no hold-out batch) run through
  read, placement, consensus, and review; the rule's views of the HVS-II
  window (303–315) were dumped per sample by a research-only build.
- **Statistic.** A window run gives a cross-strand pair when each strand has
  exactly one decoded view (direct or derived, not in phase) with the same
  top length and a defined `SE`: `d = p_F − p_R`, `v = SE_F² + SE_R²`, and
  `z² = d² / (v + τ_F² + τ_R²)` under the configured τ. A batch's dispersion
  φ is the mean `z²` of its pairs; the corpus reference is the mean over all
  pairs (run 1: 0.4643 over 13,382 pairs; run 2: 0.0563 over 14,384).

### Variation across batches

Run 1 has a systematic strand bias, mean `d` = −0.070 (forward direct below
reverse derived) in every period and on every instrument (−0.063 to −0.080);
rms `d` is 0.096 against rms `SE` 0.056. Run 2's mean `d` is −0.008 and its
rms `d` (0.042) is below rms `SE` (0.050): the decode's own `SE` explains
run 2's strand spread, and run 2's τ is far above it.

Per batch (at least 10 pairs; 732 batches for run 1, 760 for run 2),
dispersion relative to the corpus:

| Run | q10 | q50 | q90 | q99 | Observed var | Sampling var | Between-batch var | Prior pseudo-pairs n0 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.50 | 0.96 | 1.47 | 1.81 | 0.132 | 0.037 | 0.096 (sd 0.31) | 5.8 |
| 2 | 0.27 | 0.74 | 1.92 | 4.34 | 0.723 | 0.732 | ≈ 0 | not identifiable by batch (64 by HV3R instrument run) |

- Grouping by the HV2F or HV3R instrument run gives the same quantiles
  (run 1 between-run sd 0.31–0.32).
- Run 2's `z²` is heavy-tailed (per-pair relative variance 11.8); its batch
  differences are sampling noise.
- The empirical-Bayes shrunk scale per batch,
  `√((n0 + n·φ/ref) / (n0 + n))` with n0 = 6 (run 1) and 64 (run 2):
  run 1 q10 0.79, q50 0.99, q90 1.15, max 1.34 (above 1 in 46 % of batches,
  above 1.2 in 3–5 % per period); run 2 q10 0.92, q50 0.97, q90 1.10, max
  1.57. A run-calibrated τ moves by at most about ±35 %.

By period (pooled dispersion over the corpus). An in-phase check is a read
in phase through the run where the other strand gives a clear length; it
disagrees when the lengths differ.

| Period | Run 1 pairs | Run 1 φ/ref | Run 2 φ/ref | In-phase checks (runs 1 + 2) | Disagreeing | Rate |
|---|---:|---:|---:|---:|---:|---:|
| 2024 | 91 | 1.11 | 1.24 | 1,321 | 1 | 0.08 % |
| 2025-01..07 | 2,395 | 1.04 | 0.95 | 36,144 | 20 | 0.06 % |
| 2025-08..12 | 2,615 | 1.02 | 1.45 | 48,096 | 22 | 0.05 % |
| 2026 | 8,281 | 0.98 | 0.87 | 158,535 | 37 | 0.02 % |

Older runs are slightly noisier, but the trend is small next to the
batch-to-batch spread, and the 2024 development batches give only 91 pairs.
By instrument (HV2F read), run 1 ranges 0.78–1.11 and run 2 0.77–1.10. Of the
two 2024-11 development batches, one is in the noisy tail (run 2 φ/ref 2.28
over 22 pairs, run 1 1.48 over 18; scales 1.15 and 1.17) and one is ordinary.

### Risk against truth

- **One-strand in-phase confirmations** (#204's rule replayed on the review
  records): 2,491 runs, 1 wrong (a run-1 C6 read in phase on the reverse
  strand where the truth has C7, and the forward read shows a different
  structure), 2 with an IUPAC base at 302 or 316, the rest right or
  uncovered. The wrong case's batch is quiet (scale 0.91). The in-phase
  disagreement count is 0 in almost every batch, so a per-run in-phase
  reliability cannot be estimated; the one wrong case is sample-level, not
  run-level.
- **Rule A** (`decide` replayed on every run): 26,958 runs decided, 6 wrong.
  Four are run-1 C9 decoded on both strands where the truth has C8
  (`p_F` ≈ 0.80, `p_R` ≈ 0.70–0.75: a C8/C9 mix the laboratory reads as its
  shorter length); two are both-strand in-phase (C9 against C8, C6 against
  C5). Their batch scales are 0.89–1.19; even a 35 % larger τ keeps each
  decided.
- **Error rate by batch scale.** DNA's window run-1 length error rate is
  2.9 % at scale ≤ 1.0, 4.2 % at 1.0–1.2, and 3.8 % above 1.2;
  review-confirmed and wrong: 4, 4, and 0. Batches with fewer than 10 pairs
  carry the highest error rate (9.7 %), all flagged.

**Answer.** On development panel v5 a per-batch noise level does not separate
the risky runs: the wrong confirmations sit in ordinary batches, and rule A's
decisions sit far from the τ boundary. Only one 2024-11 batch fits the
hypothesis, and development data hold no truth-wrong confirmation there to
test it.

### Design sketch

The draft run-context proposal (`PROP-0014`,
[PR #212](https://github.com/tsuuanmi/DNA/pull/212)) takes this sketch up.

- **Unit: the laboratory batch** (what is submitted together), not the
  instrument run. The statistic is cross-strand, and a sample's two strands
  sit on different instrument runs 45 % of the time; the instrument run gives
  the same spread.
- **Owner.** Not the Sanger run signal, which is reference-free
  ([[design/signal-processing|signal processing]]): the length views
  need placement and the consensus's own length views. The pure functions
  belong with the consensus length rule, and a pipeline command beside
  `dna signal` builds the context once per batch.
- **Contract.** Sufficient statistics only: per recorded anchored-window run,
  `{run, pairs, dispersion}`, with provenance (configuration digest, plugins,
  reference and profile identities, the τ the pairs were normalised with, a
  selection version) and the samples' read digests. The consumer applies its
  own validated calibration and subtracts the sample's own pairs
  (leave-one-out). It refuses a context whose τ or selection differs from its
  configuration.
- **Calibration** (configuration, not code): a minimum pair count (10), the
  prior pseudo-pairs (6 and 64), the corpus references (0.4643 and 0.0563), a
  minimum scale of 1.0 (τ never below the replicate estimate), and a
  one-strand scale of 1.0 (one-strand in-phase confirmation only in a batch
  no noisier than the corpus). `τ_run = τ · scale`, with
  `scale = max(minimum, √(φ̃ / ref))` and `φ̃ = (n0·ref + n·φ) / (n0 + n)`.
  Without a context, or below the minimum, τ is the configured one and no
  one-strand confirmation is made: fail closed, today's behaviour.

### Prototype

An unmerged prototype on `main` `60252579`, evaluated with `dna-eval diff`
against the clean `main` run; contexts built truth-free from the development
batches' full folders. Every diff passed the gate.

| Variant | window@HV2 checks | 2024 | 2025-01..07 | 2025-08..12 | 2026 | Clear samples | Unchecked non-het blocks |
|---|---:|---:|---:|---:|---:|---:|---|
| `main` | 3,068 | 11 | 385 | 699 | 1,975 | 11,340 | baseline |
| V0 prototype, no context | 3,068 | 11 | 385 | 699 | 1,975 | 11,340 | identical |
| V3 calibrated τ only | 3,068 | 11 | 385 | 699 | 1,975 | 11,340 | identical |
| V1 τ + one-strand where scale ≤ 1 | 2,939 (−129) | 10 | 356 | 678 | 1,897 | 11,370 | unchanged in every period × zone |
| V2 τ + one-strand wherever a context exists | 2,775 (−293) | 10 | 333 | 639 | 1,795 | 11,408 | unchanged in every period × zone |

The period columns count samples with a check over 303–315. V1 and V2 add
only re-check-only passes (24 → 27 and 29) and het-only passes (347 → 349
and 350); passed non-het stays 12 in every period and source, and the
unchecked non-het blocks are identical per period and zone. For reference,
#204 alone removed 475 window@HV2 checks on its base; requiring a context
(at least 10 pairs) costs about 630 would-be confirmations in the 270
batches with few pairs.

### Conclusions

- Calibrating τ per batch is safe and inert on development panel v5: rule
  A's decisions are not near the τ boundary.
- Gating #204 on the batch scale is not supported by development truth: the
  gate halves the yield (−129 against −293) and blocks nothing that is
  wrong. Whether it blocks the two hold-out escapes needs a run on that
  panel by its owner (V1 against V2, contexts built truth-free from the
  hold-out batches' full folders).
- The one wrong one-strand case is sample-level: the other strand's raw runs
  show a different structure. A sample-level guard (no contradicting
  structure on the other strand, run 1) would block it at the cost of one
  right case. For run 2 the forward read's structure differs in 1,395 right
  cases, so the guard does not transfer. One case is too few to adopt it.
- Run 2's batch differences are not identifiable with mean `z²`; a robust
  dispersion (trimmed or median) is the next step if run-2 calibration is
  pursued.

## Open points

- An orthogonal truth subset (amplicon sizing) to calibrate shares and any
  reporting ratio ρ.
- Class A as a note: blocked on the hold-out evidence (readiness-v2 or the
  planned hold-out set).
- The HV2 coverage note and the HV3 `not_recorded` convention: reporting
  decisions for the user.
- Run calibration: the hold-out run of V1 against V2; a robust run-2
  dispersion.

## Related

- [[PROP-0012]]: the trace crate, its decode, and the length decision;
  [[ADR-0074]]: how a run length is known; [[ADR-0082]] and [[ADR-0084]]:
  phase mixtures and the trace crate.
- [[design/length-mixture]], [[design/variant-nomenclature]],
  [[reference/review]].
- [Phase recovery](../phase-recovery/README.md): the earlier, concluded study
  of dephased stretches behind homopolymer runs.
