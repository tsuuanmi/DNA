# ADR-0001: Deliver an End-to-End MVP First

- **Status:** Accepted
- **Date:** 2026-08-22

## Context

Apollo contains many commands and specialized algorithms. Porting every module
before producing a working result would delay feedback and make parity failures
hard to isolate.

## Options

1. Port modules in Apollo file order and integrate at the end.
2. Copy the existing Rust port wholesale.
3. Deliver one narrow AB1-to-variant vertical slice, then add features by value.

## Decision

Choose option 3. The MVP covers one AB1, one rCRS FASTA, signal-derived
re-calling at PLOC loci, end-quality control, direct alignment, and basic
primary-sequence variants.

## Consequences

The project obtains usable end-to-end feedback early and can validate each stage
against one coherent flow. Many Apollo commands remain explicitly deferred, and
the MVP is not a drop-in Apollo replacement.

## Supersession

A broader release scope may supersede this ADR only after the MVP acceptance
gates pass.
