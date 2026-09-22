# Production Release Evidence

Copy this template for every release candidate that is intended to carry the
**production-ready** label. A GitHub Release artifact alone is only
**build verified**.

## Identity

- DNA version:
- source revision:
- release tag:
- supported target:
- Rust toolchain:
- Cargo toolchain:
- Cargo.lock SHA-256:
- release artifact SHA-256:
- SPDX SBOM:
- GitHub artifact attestation:

## Engineering gates

- Rust quality CI:
- MSRV:
- dependency policy (`cargo-deny`):
- RustSec lockfile audit:
- packaged auditable binary `.dep-v0` present: yes/no
- auditable binary dependency scan:
- dependency review:
- CodeQL:
- OpenSSF Scorecard reviewed:
- ABIF fuzz smoke:
- extended fuzz campaign:
- schema/config/reference validation:

Every item above must identify the exact CI run, artifact, or review record used.

## Scientific validation

- approved corpus identifier/version:
- number of samples/traces:
- ground-truth authority:
- evaluated loci/regions:
- SNV results:
- insertion results:
- deletion results:
- noisy/mixed-signal results:
- reverse-strand results:
- circular-boundary results:
- known difficult loci:
- unresolved disagreements:
- exclusions and rationale:

A disagreement must be analyzed and classified; it must not be hidden by changing
the denominator or silently suppressing an observation.

## Resource and failure validation

- malformed/truncated ABIF corpus:
- maximum accepted AB1 size:
- maximum accepted reference size:
- alignment allocation bound:
- peak memory measurement:
- runtime measurement:
- panic/crash count:
- partial-publication checks:

## Release decision

- build verified: yes/no
- scientifically validated for documented corpus/domain: yes/no
- production-ready contract satisfied: yes/no
- approver:
- date:
- remaining limitations:

The production-ready label is valid only if all mandatory requirements in
ADR-0018 are satisfied for the exact release revision.
