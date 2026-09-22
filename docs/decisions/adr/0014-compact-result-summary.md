# ADR-0014: Project a Compact Analysis Result Summary

- **Status:** Superseded by ADR-0026
- **Date:** 2026-08-25

## Context

`dna.analysis/v4` exposed bounded rolling signal windows, complete called and retained sequences, gapped alignment rows, method identifiers, and rich per-variant peak/vendor evidence. Those fields were useful while establishing the scientific stages, but they duplicated internal state and produced a result larger and more identifying than downstream review requires.

The stable review needs are narrower: deterministic input identities, read/trim extent, merged noisy regions, alignment quality/placement, normalized differences with enough mapped supporting evidence to audit configured eligibility, and warning counts. Full scientific evidence remains available from the identified AB1, reference, and configuration; it need not be duplicated in every JSON result. Software/build provenance is deferred during pre-release development until a stable versioning strategy is defined.

Configuration and result contracts have independent version domains. The signal-processing settings introduced with strict configuration schema version 4 remain scientifically current even though the JSON projection changes.

## Options considered

1. Keep v4 unchanged.
2. Add a second summary file while retaining v4 output.
3. Replace v4 with one compact v5 result and no compatibility output.

Option 1 preserves unnecessary payload and privacy exposure. Option 2 creates two authorities, complicates atomic publication, and violates DNA's one-result boundary. Option 3 keeps one deterministic contract and makes every retained field intentional.

## Decision

Adopt `dna.analysis/v5` as the only scientific result.

V5 retains:

- input AB1 SHA-256, reference name/topology/SHA-256, and configuration SHA-256;
- original read call count and 0-based half-open trim interval;
- merged candidate-noisy call/sample regions with minimum primary SNR;
- selected alignment orientation, callable bases/identity, unresolved bases, gap opens, reference segments, and origin-wrap flag;
- normalized variants with position, reference/alternate alleles, kind, and concise original-call mappings;
- for supporting calls only, `maximum_peak_height` and uncalibrated `relative_quality`, matching the configured variant evidence gates;
- counts of unresolved primary calls, multi-channel unresolved calls, and excluded variant candidates.

V5 removes:

- trace filenames and full primary, ambiguity, and retained sequences;
- individual rolling signal windows and secondary-SNR details;
- gapped query/reference rows, operation runs, alignment score, traceback columns, and redundant alignment counts;
- program/determinism/method constants;
- full A/C/G/T peak height/position/source objects, penalties, Phred flags, and vendor PBAS/PCON data;
- variant contig/classification/normalization labels when the contract already fixes or implies them;
- warning totals, vendor disagreement counts, and duplicated origin-wrap warning state.

Supporting calls preserve role, original call index, PLOC, trace-strand primary/ambiguity, optional biological position, maximum peak height, and relative quality. Flanking calls preserve mapping context only. Inserted supporting calls continue to omit a fabricated biological position; deletions continue to use flanks without fabricated deleted-base evidence.

Strict configuration remains schema version 4. DNA emits no v4 compatibility result, aliases, duplicate fields, VCF/BCF, or alternative output path.

## Consequences

Results are smaller, less identifying, and focused on review decisions while retaining deterministic input identities and the exact compact evidence used for supporting-call eligibility. Consumers must use the schema rather than rely on removed draft fields. Detailed reconstruction requires the source AB1, reference, and configuration; software/build identity is intentionally unavailable in the current pre-release JSON contract.

Rolling windows, complete sequences, full peaks, vendor evidence, and traceback artifacts remain internal algorithm inputs or diagnostics; removing them from JSON does not change scientific behavior. Candidate-noisy regions remain observation-only. The core CLI continues atomic no-overwrite publication and separate append-only operational logging.

## Relationship to prior decisions

This ADR supersedes only the output projection selected by ADR-0013. ADR-0013 remains accepted for the observation-only `dna.windowed_snr/v1` method, merged-region semantics, configuration schema version 4, and prohibition on noise-driven scientific behavior without validation. ADR-0012's concise mapped-call principle remains in force, with v5 reducing the mapped evidence to supporting maximum peak height and relative quality.

## Follow-up

ADR-0026 replaces the v5 public call-pointer shape with analysis v6 reference-oriented base/peak/quality evidence. The v5 decision remains historical context for the compact top-level summary.
