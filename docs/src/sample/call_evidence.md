# `src/sample/call_evidence.rs`

## Purpose

Resolves one source call's authoritative basecall-independent signal evidence and
projects all A/C/G/T channels into the selected reference orientation for sample
reconciliation.

## Responsibilities

- Look up the matching `LocusEvidence` by original call index and require the stored index to match.
- Project corrected A/C/G/T amplitudes and per-channel SNR into reference orientation.
- Preserve the optional normalized `EvidenceProfile` and complement reverse-read A/T and C/G channels.
- Preserve existing merged candidate-noisy-region membership for the same call.
- Return one coherent `CallDNAEvidence` value from a single authoritative lookup.

## Non-responsibilities

No signal recomputation, profile construction, SNR thresholding, read admission,
local contribution eligibility, consensus weighting, variant filtering, artifact
classification, or public JSON projection.

## Invariants

Forward evidence is unchanged. Reverse evidence reorders A/C/G/T as T/G/C/A so
corrected amplitudes, SNRs, and profile weights all share the same
reference-oriented channel frame. No value is reconstructed from the primary
call, ambiguity code, relative quality, or reference base.

Candidate-noisy membership remains observation-only and cannot by itself reject
or down-weight the call.

## Status

Implemented.
