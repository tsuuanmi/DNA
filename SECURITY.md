# Security Policy

DNA processes untrusted sequencing files and potentially identifying biological
data. Security reports should therefore be handled privately.

## Supported versions

Until the first evidence-backed production release, only the latest revision on
`main` is supported for security fixes. After a production release exists, the
latest production release line and `main` are the supported security surfaces
unless a release note states otherwise.

## Reporting a vulnerability

Use GitHub private vulnerability reporting from the repository **Security**
page. Do not open a public issue containing exploit details, private sequencing
data, sample identifiers, complete sequences, or chromatogram payloads.

A useful report includes:

- the affected commit or release;
- the smallest non-sensitive reproducer possible;
- the expected and observed behavior;
- impact and preconditions;
- whether malformed ABIF/FASTA/configuration input is involved;
- whether the issue can cause panic, resource exhaustion, path escape,
  overwrite, information disclosure, or incorrect result publication.

If a biological sample is needed to reproduce the issue, replace it with a
synthetic/minimized fixture whenever possible.

## Security boundaries

Particularly sensitive surfaces include:

- ABIF/AB1 binary offsets, lengths, cardinalities, and allocations;
- FASTA and TOML parsing;
- filesystem publication and cleanup paths;
- dependency and build supply chain;
- release provenance and artifact integrity;
- logs and output containing sample-linked biological evidence.

The detailed engineering trust model is documented in
[docs/security/README.md](docs/security/README.md).

## Coordinated disclosure

Please allow maintainers to reproduce, fix, test, and prepare a release before
public disclosure. Security fixes should include a regression test or minimized
fuzz corpus entry when practical.
