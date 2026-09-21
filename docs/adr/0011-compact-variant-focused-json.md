# ADR-0011: Emit Compact Variant-Focused JSON

- **Status:** Superseded in part by ADR-0012
- **Date:** 2026-08-23

## Context

The complete `dna.analysis/v1` document duplicated large decoded channel arrays, every call and quality record, all intermediate stage details, and the full reference/config expansion. Typical result files were approximately one megabyte each even when users primarily needed the called sequence, selected alignment, and variants. Variant review still requires local chromatogram peaks and quality.

The command also accepted an output prefix, while operational use expects results consistently under the repository `results/` directory.

## Options

1. Keep the complete v1 document.
2. Remove bulk trace/intermediate records and retain peak/quality evidence only around variants.
3. Keep both schemas or add a compatibility switch.

## Decision

Choose option 2 without a compatibility path. `signal analyze <trace.ab1> --reference <reference.fasta>` writes `results/<trace-stem>.json` using `dna.analysis/v2`.

The document retains deterministic input/reference/config identities, called and retained sequences, the selected alignment including gapped rows, normalized variants, and compact warning counts. Every variant-associated call carries its PLOC position, all four selected channel peaks, relative quality, and applicable vendor quality. Insertions include inserted calls and flanks; deletions include flanks only because no deleted-base trace call exists.

Complete channel arrays, non-variant call/quality tables, decoded-tag inventories, losing orientation candidates, excluded-candidate details, full reference/config expansion, and verbose warning records are omitted. Existing result targets remain no-overwrite.

## Consequences

Results are substantially smaller while retaining the local evidence needed to inspect each reported difference. The schema change is intentionally breaking. Consumers must migrate to `dna.analysis/v2`; DNA does not emit duplicate v1 fields or accept the removed `--out-prefix` option.

The output location depends on the process working directory. Production invocations must run from the intended project/work directory so `results/` resolves correctly.

## Supersession

This ADR supersedes ADR-0008's complete `dna.analysis/v1` shape and prefix-derived filename. ADR-0008's JSON-only and no-VCF decisions remain in force.
