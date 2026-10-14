# Research: how other systems decide pass, review, and rerun

**Status:** exploratory. This note has no production authority. It supports
[PROP-0004](../../proposals/0004-review-triage.md).

## Question

DNA's [purpose](../../requirements/SRS.md#purpose) is to sort work into pass,
review, and rerun. Other fields already release results without a human, hand
the rest to people, and send failures back. What can DNA reuse from them for
the *system*: rules, thresholds, approval, and determinism? The biology is out
of scope here.

## Clinical laboratory autoverification

Clinical chemistry labs release most results without a human. CLSI AUTO10-A
(2006) and AUTO15 (2019, reaffirmed 2025) guide this practice.

- **Rules, not scores.**
  - Each rule checks one condition: instrument flags, QC status, measurement
    range, critical values, a delta check against the patient's last result,
    and sample indices.
  - A result is released only when every rule passes; any failing rule holds
    it for review.
- **Each lab owns its thresholds.** The guidelines are templates. Each lab sets
  and validates its own limits.
- **Validation before release.** Labs test the rule set on simulated and real
  cases and compare it with manual review. The goal stated in the literature
  is zero defects in auto-released results.
- **Lifecycle.** The scope runs from design through testing, validation, and
  go-live to ongoing support. Rules are signed off before go-live and
  monitored after it.

Lessons for DNA:

- pass is an allow-list of conditions, each checked by a named rule;
- thresholds belong to a versioned policy;
- a rule acts only after it has been validated and signed off.

## Selective classification

Machine learning calls the same idea *selective classification* or the
*reject option*: a classifier may abstain when unsure.

- **Risk–coverage trade-off.**
  - *Coverage* is the share of inputs that get answered.
  - *Selective risk* is the error rate on the answered inputs.
  - Abstaining more lowers the risk.
- **Bounded-improvement model.** Fix a target risk, then maximise coverage
  subject to it. This avoids pricing a reject against an error.
- **Guaranteed risk.** Geifman and El-Yaniv choose the rejection threshold on
  held-out data. The target risk then holds with a stated probability.

Lessons for DNA:

- A is the coverage and the error rate of passed results is the selective
  risk;
- choose thresholds to maximise A subject to a target error rate with a
  stated confidence;
- tune on one set and confirm on a disjoint one.

## Sanger trace quality metrics

Vendor software summarises a trace with a few numbers:

- **Trace (sample) score** — the mean base quality in the clear range;
- **QV20+** — the number of bases with quality 20 or more;
- **Q20 CRL** — the longest stretch of such bases.

The thresholds are user-set flag settings. No published pass/fail cut-offs
were found; they depend on the lab and the instrument.

Lessons for DNA:

- rerun thresholds must come from this lab's data;
- DNA's relative quality is not calibrated (SRS-QC-002), so a vendor-style
  QV20 cut-off cannot simply be copied;
- callability (callable span, masked fraction) already gives a
  signal-derived counterpart.

## Lab practice for failed traces

- **The cause decides the action.** Sequencing core facilities choose the
  rerun action by cause:
  - off-scale signal → re-inject;
  - dye blobs → clean up again;
  - a failed or weak reaction → fix the template;
  - a mixed template → PCR again or extract again.

  Repeating an identical reaction rarely fixes a failure, so a rerun request
  must carry its cause.
- **Both strands.** Forensic mtDNA practice following SWGDAM asks for
  concordant data from both strands. It accepts single-strand data only when
  it is confirmed.

## Summary for the design

| Need | Practice to reuse |
|---|---|
| pass is safe | allow-list of validated rules; zero known defects on validation data |
| thresholds are justified | maximise coverage under a target error rate on held-out data |
| behaviour is stable | versioned policy; rule sign-off; monitoring after go-live |
| rerun is useful | a cause code that maps to a lab action |
| review is fast | flags that say which rule held the result and why |

## Sources

- CLSI AUTO15, *Autoverification of Medical Laboratory Results for Specific Disciplines*: <https://clsi.org/shop/standards/auto15/>
- LIS-based validation system for autoverification toward zero defects: <https://LIMSwiki.org/index.php/Journal:Development_and_implementation_of_an_LIS-based_validation_system_for_autoverification_toward_zero_defects_in_the_automated_reporting_of_laboratory_test_results>
- Autoverification rule design and validation in practice: <https://www.biochemia-medica.com/en/journal/32/3/10.11613/BM.2022.030704>
- El-Yaniv and Wiener, *On the Foundations of Noise-free Selective Classification*, JMLR 2010: <https://jmlr2020.csail.mit.edu/papers/volume11/el-yaniv10a/el-yaniv10a.pdf>
- Geifman and El-Yaniv, *Selective Classification for Deep Neural Networks*, NeurIPS 2017: <https://papers.NeurIPS.cc/paper_files/paper/2017/file/4a8423d5e91fda00bb7e46540e2b0cf1-Paper.pdf>
- Franc, Prusa and Voracek, *Optimal strategies for reject option classifiers*: <https://arxiv.org/abs/2101.12523>
- Thermo Fisher, Sequencing Analysis Software user guide (trace score, QV20+, CRL): <https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_064539.pdf>
- Thermo Fisher, *Troubleshooting Sanger sequencing data*: <https://assets.thermofisher.com/TFS-Assets/LSG/manuals/MAN0014435_Trbleshoot_Sanger_seq_data_UB.pdf>
- SWGDAM mtDNA interpretation guidelines (2024): <https://swgdam.org/_files/ugd/4344b0_f61de6abf3b94c52b28139bff600ae98.pdf>
