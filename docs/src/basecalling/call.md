# `src/basecalling/call.rs`

## Purpose

Re-calls every vendor-defined PLOC locus from the analyzed channel signals,
producing a primary sequence plus per-call ambiguity evidence.

## Responsibilities

- Build call windows and per-channel peaks for each PLOC locus.
- Rank channels by selected peak height and select the primary call, applying
  the configured secondary-peak ratio both to each selected channel peak and to
  that channel's signal at the uniquely strongest primary peak position.
- Emit the primary and IUPAC ambiguity calls with their qualifying channels.
- Preserve the unique primary event's selected sample position and all four
  analyzed-channel values at that shared sample.
- Retain each call's validated sample-window bounds for downstream
  observation-only signal analysis.
- Record whether the vendor primary agrees with the signal-derived call, without
  retaining duplicate vendor/call fields or using it as the result.

## Non-responsibilities

No ABIF parsing, end trimming, reference alignment, or variant calling.

## Key types and functions

- `call(trace, config) -> Result<BaseCalls>`: the module entry point, re-exported
  from `mod.rs`.

## Invariants and errors

- Calls derive from the four signal channels at validated PLOC loci.
- Every selected peak position must remain inside its call window; a violated
  invariant returns `Error::Basecalling`.
- A non-positive or exactly tied strongest peak yields an unresolved `N` call
  with no primary-event evidence.
- For a uniquely strongest positive peak, primary-event evidence stores exactly
  the A/C/G/T trace values used by the co-localization eligibility check. A
  channel qualifies only when its selected peak is positive and `selected_height
  / primary_height >= secondary_peak_ratio`, and its stored primary-event signal
  is positive and `colocated_height / primary_height >= secondary_peak_ratio`.
- This intersection can only remove remote secondary channels from the prior
  selected-peak rule; it cannot add a new qualifying channel or change the
  uniquely strongest primary channel.
- One qualifying channel gives a canonical call, two give a two-base IUPAC code,
  three keep the strongest primary with unresolved `N` ambiguity, and four are
  unresolved `N` for both primary and ambiguity.
- The call vector and retained `primary_sequence` have equal lengths; ambiguity is
  retained per call rather than duplicated in a second aggregate string.
- Vendor PBAS is evidence only and never replaces the signal-derived call.

## Dependencies

- `iupac` and `peak`.
- `config` for `BasecallingConfig`.
- `model::basecalls` and `model::trace`.
- `error` for `Result`.

## Biological semantics

The secondary-peak ratio distinguishes a clean single-channel call from a mixed
signal position without inferring genotype. A secondary channel must be strong
both at its selected peak and at the primary event's sample, preventing a remote
maximum elsewhere in the same PLOC window from creating ambiguity. A qualifying
second channel produces a two-base IUPAC ambiguity; a tie or a non-positive
strongest peak is left unresolved rather than guessed.

## Tests

- `calls_unambiguous_strongest_channel`: a single strong channel yields the
  canonical call and stores A/C/G/T values at its primary-event sample.
- `exact_strongest_tie_is_unresolved`: an exact tie between two channels yields
  `N` with no primary-event evidence.
- `non_positive_channels_have_no_primary_evidence`: an all-non-positive locus
  yields no primary-event evidence.
- `off_locus_secondary_peak_does_not_create_ambiguity`: a remote secondary
  maximum fails the primary-sample signal gate while remaining distinct from its
  lower primary-event dna.
- `overlapping_offset_secondary_peak_still_qualifies`: an offset secondary peak
  remains eligible when its signal overlaps the primary sample; its selected and
  primary-event values remain distinct.
- `offset_secondary_signal_at_exact_ratio_threshold_qualifies`: an offset peak
  qualifies when its signal at the primary sample is exactly on the inclusive
  ratio boundary.
- `ploc_fallback_primary_always_qualifies`: a uniquely strongest PLOC fallback
  remains its own qualifying primary channel and stores PLOC-sampled evidence.
- `ploc_fallback_secondary_must_be_colocated`: a strong remote PLOC fallback does
  not create ambiguity.
- `three_qualifying_channels_keep_primary_but_not_ambiguity`: three qualifying
  channels keep the strongest primary while the ambiguity is unresolved `N`.
- `four_qualifying_channels_are_fully_unresolved`: four qualifying channels yield
  unresolved `N` for both primary and ambiguity.

## Status

Implemented.
