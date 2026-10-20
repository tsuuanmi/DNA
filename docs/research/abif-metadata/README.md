---
id: RES-abif-metadata
status: active
summary: Which ABIF tags beyond the five DNA reads (DATA.9-12, FWO_.1, PLOC.2, PBAS.2, PCON.2) could add evidence, and how many carry an instrument-side failure or run context.
---

# Research: unused ABIF tags

**Status:** active — exploratory evidence, not production authority.

## Question

DNA reads five things from an AB1: the analyzed channels (`DATA.9`–`DATA.12`),
the channel order (`FWO_.1`), the call positions (`PLOC.2`), and the vendor
call and quality (`PBAS.2`, `PCON.2`). A trace holds roughly 200 more tags. Which
could add evidence, and for whom: the core (what a read shows), the modality
(how far to trust a call), or the laboratory (which batch, run, or instrument
caused a loss)?

## Survey

300 raw AB1 files of the science corpus (3500-series instruments, BDTv3.1
chemistry, KB basecaller 1.4.x), tags listed from the ABIF directory only.
Counts below are over those 300 files; they are a survey of what exists, not a
measured effect on DNA's results.

| Tags | Content | Seen in the sample |
|---|---|---|
| `BcRs.1`, `QcRs.1` | the basecaller's failure reason, for example *Failed estimation of average peak width* or *Analysis stopped due to low …* | 99 of 300 files non-empty (63 average peak width, 18 stopped early, 7 width curve, rest mixed) |
| `PuSc`, `TrSc`, `QV20`, `CRLn` (`.1` value, `.2` pass/check/fail) | instrument sample, trim, QV20, and contiguous-read scores with a verdict | PuSc: 6 check, 6 fail; QV20: 60 check; CRLn: 88 check, 1 fail |
| `CLRG.1/.2`, `CRLn.1` | instrument clear range and contiguous read length | present in all |
| `P1AM`, `P1WD`, `P1RL`, `P2AM`, `P2RL`, `P2BA` | the basecaller's per-base peak amplitude, width, resolution, and its own second peak (base and amplitude) | present in all, one value per called base |
| `NOIS.1`, `S/N%.1`, `SPAC.1`, `Scal.1` | per-channel noise, per-dye signal-to-noise, average peak spacing, scale | present in all |
| `DATA.1`–`DATA.4` | raw (unanalyzed) channels, before baseline and mobility correction | present in all |
| `DATA.5`–`DATA.8`, `AmbT.1` | run voltage, current, power, temperature, and ambient temperature over the run | present in all |
| `HCFG.4`, `MCHN`, `SVER`, `InjN`, `RunN`, `RUND/RUNT`, `Tmpr`, `EPVt`, `InVt`, `InSc`, `LsrP` | instrument serial, software, run name, run time, temperature, voltages, injection, laser power | 3 instrument serials, 3 runs |
| `ABLt`, `CBLt`, `SMLt` and their expiry dates | anode buffer, cathode buffer, and sequencing-mix lot | buffers one lot each, sequencing mix 2 lots (225 and 75 files) |
| `SMPL`, `TUBE`, `LANE`, `CTNM`, `RGNm`, `LIMS` | sample name, well, lane, container, region name, LIMS id | present in all |

## Candidate uses

None is adopted; each needs a proposal, and none may change placement or a
call (SRS-ALN-008, SRS-CALL-001).

1. **Laboratory causes, not DNA output.** `BcRs`/`QcRs`, the pass/check/fail
   scores, the instrument serial, run, lots, temperature, and voltages let the
   triage name a *cause* (laboratory, instrument, chemistry lot) beside a
   region's review label and group failures by run or lot. This serves the
   stated goal of routing each discrepancy to its owner. It belongs in the
   evaluation harness and the review list as labels, not in the core.
2. **A second opinion on the second peak.** `P2BA`/`P2AM` are the instrument's
   own second base and amplitude per call. Comparing them with DNA's own
   second-peak cause label and run-signal level would test both. The comparison
   belongs in `dna-eval`, not in a plugin crate.
3. **Instrument noise as a prior.** `NOIS` and `S/N%` give a per-channel noise
   floor that DNA estimates today only from the trace itself; they could seed or
   sanity-check the windowed noise estimate.
4. **Raw channels.** `DATA.1`–`DATA.4` are the input of the instrument's own
   conditioning. The denoising research reads the analyzed channels only; the
   raw ones would show what the instrument removed (baseline, mobility shift).
5. **Run-level health.** `DATA.5`–`DATA.8` and `AmbT` show voltage or current
   excursions during a run, a candidate explanation for a whole run's early
   signal loss.
6. **Label consistency.** `SMPL`, `TUBE`, `LANE`, and the file name should agree
   with the primer direction implied by placement; a disagreement points to a
   sample or plate mix-up. It can be reported as provenance, never used to place
   a read.

## Not yet established

- Whether `BcRs` flags predict DNA's own rejected, masked, or review-listed
  reads. The survey counts flags; it does not join them to results.
- Whether the lot and run groupings in the sample are large enough to
  separate causes; it holds 3 runs and 2 sequencing-mix lots.
- Privacy: tag values carry sample names, container names, user names, and
  instrument serials; any adoption must follow [[governance/data]].

## Next step

Join `BcRs`/`PuSc`/`CRLn` flags and run/lot to the current slice's per-sample
outcomes (batch date from 2026-01) in `dna-eval` and report how many rejected,
masked, and review-listed reads each flag explains, before proposing anything.
