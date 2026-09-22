# ADR-0005: Use Versioned JSON and Standards-Compliant VCF

- **Status:** Superseded by ADR-0008
- **Date:** 2026-08-22

## Context

Apollo's flat JSON mixes trace, alignment, and variant fields without a clean
version boundary. The previous Rust port also exposed incomplete BCF behavior.

## Options

1. Reproduce Apollo JSON and all file names exactly.
2. Define a DNA JSON contract and valid VCF text, deferring BCF.
3. Emit JSON only.

## Decision

Choose option 2. JSON identifies `dna.analysis/v1` and is governed by a checked
schema. VCF is valid VCF 4.2 text with explicit coordinates and normalized
variants. DNA never writes VCF text with a `.bcf` extension.

## Consequences

Consumers receive explicit, evolvable contracts and standards-compliant variants.
Apollo consumers need an adapter if exact legacy shape is required. Schema
changes require review, compatibility notes, and either backward compatibility
or a new schema version.

## Supersession

BCF support requires a compliant writer, an independent parser test, and a new
ADR or amendment through supersession.
