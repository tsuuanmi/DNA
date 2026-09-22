# System Context

DNA is a local scientific CLI that processes Sanger ABIF/AB1 traces.

## External inputs

- ABIF/AB1 chromatogram files;
- one short FASTA reference for reference-guided operations;
- strict TOML configuration;
- explicit sample identifiers and trace paths for sample aggregation.

## External outputs

- one versioned command-specific JSON result per successful core invocation;
- separate append-only operational logs;
- batch orchestration may group per-trace and aggregate outputs without changing core scientific semantics.

## Trust boundary

Input bytes, paths, filenames, manifests, and configuration are untrusted. Sensitive biological data is governed by [data policy](../governance/data.md) and [security](../security/README.md).

The core is not a network service and currently has no deployment topology. A deployment document should be added only if a real deployed runtime boundary is introduced.
